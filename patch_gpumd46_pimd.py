#!/usr/bin/env python3
"""
Patch clean GPUMD-4.6 source to support three PIMD internal-mode spectra:

  trotter   : original GPUMD ring-polymer frequencies
  matsubara : omega_k / omega_n = 2*pi*min(k,P-k)/P
  eco       : optimized Eco-PIMD frequencies following eco.f strictly

This patcher is intentionally source-code only. It does not edit the makefile and
it does not add machine-specific LAPACK/OpenBLAS linker flags. The patched source
calls LAPACK dsyev_ for the Eco optimizer, exactly as eco.f does, so the user must
link LAPACK/OpenBLAS when compiling GPUMD.

Tested against the GPUMD-4.6 PIMD source layout:
  src/integrate/ensemble_pimd.cuh
  src/integrate/ensemble_pimd.cu
  src/integrate/integrate.cuh
  src/integrate/integrate.cu
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

PATCH_MARK = "ECO_PIMD_PUBLIC_STRICT_FORTRAN_V3"

HELPER_BLOCK = r'''
// === BEGIN ECO_PIMD_PUBLIC_STRICT_FORTRAN_V3 ===
// Three PIMD internal-mode frequency choices for GPUMD:
//   PIMD_MODE_TROTTER   : original GPUMD/Trotter frequencies
//   PIMD_MODE_MATSUBARA : y_k = 2*pi*min(k,P-k)
//   PIMD_MODE_ECO       : optimized y_k from the Fortran eco.f in the SI
//
// For all schemes, this helper returns lambda_k = omega_k/omega_n.
// The free ring-polymer propagator uses omega_k = omega_n * lambda_k.
// The TRPMD/PIMD internal-mode Langevin damping uses gamma_k = omega_k/2.
// This preserves the original GPUMD damping exactly for Trotter, because
// lambda_k = 2*sin(k*pi/P) gives gamma_k = omega_n*sin(k*pi/P).
//
// The Eco optimizer below follows eco.f line by line:
//   l = n/2
//   m = nint(10*xmax)
//   x(j) = (j-0.5)*dx
//   f(j) = x(j)^2 / [ (x(j)/2)/tanh(x(j)/2) - 1 ]
//   w(k) = 2 except the even-P Nyquist mode w(l)=1
//   initial y(k) = 2*pi*k (Matsubara)
//   minimize s=(1/2m)*sum_j r(j)^2 using analytic gradient/Hessian
//   shift: diagonalize Hessian with LAPACK DSYEV, then
//          delta=max(1e-16*d(l),-2*d(1))
//          dy=-(H+delta I)^(-1)g
//   line search: c=1,1/2,1/4,..., enforcing positivity, ordering, and downhill s

extern "C" {
void dsyev_(
  char* jobz,
  char* uplo,
  int* n,
  double* a,
  int* lda,
  double* w,
  double* work,
  int* lwork,
  int* info);
}

static const double ECO_CM_TO_K = 1.4387768775039338; // h*c/k_B in K cm

struct EcoObjFortranStrict
{
  double s;
  std::vector<double> g;
  std::vector<double> h; // column-major l*l; upper triangle is used by DSYEV
};

struct PIMDModeBuildResultFortranStrict
{
  std::vector<double> factor; // lambda_k = omega_k / omega_n, k=0,...,P-1
  double rmst;
  double rmsm;
  double rmse;
  int newton_iterations;
};

static inline int eco_idx(const int row, const int col, const int l)
{
  return row + col * l;
}

static EcoObjFortranStrict eco_objfun_fortran_strict(
  const std::vector<double>& f,
  const std::vector<double>& x,
  const std::vector<double>& w,
  const std::vector<double>& y)
{
  const int m = static_cast<int>(x.size());
  const int l = static_cast<int>(y.size());

  EcoObjFortranStrict obj;
  obj.s = 0.0;
  obj.g.assign(l, 0.0);
  obj.h.assign(l * l, 0.0);

  std::vector<double> dg(l, 0.0);
  std::vector<double> dh(l, 0.0);

  for (int j = 0; j < m; ++j) {
    const double x2 = x[j] * x[j];
    double r = -1.0;

    for (int k = 0; k < l; ++k) {
      const double d = 1.0 / (y[k] * y[k] + x2);
      const double e = f[j] * w[k] * d;
      r += e;
      dg[k] = -2.0 * d * e * y[k];
      dh[k] = 2.0 * d * d * e * (3.0 * y[k] * y[k] - x2);
    }

    obj.s += 0.5 * r * r;

    for (int k = 0; k < l; ++k) {
      obj.g[k] += r * dg[k];
      for (int i = 0; i <= k; ++i) {
        obj.h[eco_idx(i, k, l)] += dg[i] * dg[k];
      }
      obj.h[eco_idx(k, k, l)] += r * dh[k];
    }
  }

  obj.s /= static_cast<double>(m);
  for (int k = 0; k < l; ++k) {
    obj.g[k] /= static_cast<double>(m);
  }
  for (int k = 0; k < l; ++k) {
    for (int i = 0; i <= k; ++i) {
      obj.h[eco_idx(i, k, l)] /= static_cast<double>(m);
    }
  }

  return obj;
}

static std::vector<double> eco_shift_fortran_strict(
  const std::vector<double>& g,
  const std::vector<double>& h_in,
  const int l)
{
  std::vector<double> h = h_in; // DSYEV overwrites h with eigenvectors
  std::vector<double> d(l, 0.0);

  char jobz = 'V';
  char uplo = 'U';
  int n = l;
  int lda = l;
  int lwork = 34 * l;
  int info = 0;
  std::vector<double> work(lwork, 0.0);

  dsyev_(&jobz, &uplo, &n, h.data(), &lda, d.data(), work.data(), &lwork, &info);
  if (info != 0) {
    PRINT_INPUT_ERROR("LAPACK DSYEV failed in Eco-PIMD Hessian diagonalization.");
  }

  // dy = transpose(h) * g  in the Fortran code. Here this is the gradient in
  // the Hessian eigenvector basis.
  std::vector<double> dy_eig(l, 0.0);
  for (int col = 0; col < l; ++col) {
    double sum = 0.0;
    for (int row = 0; row < l; ++row) {
      sum += h[eco_idx(row, col, l)] * g[row];
    }
    dy_eig[col] = sum;
  }

  const double delta = std::max(1.0e-16 * d[l - 1], -2.0 * d[0]);

  std::vector<double> f(l, 0.0);
  for (int j = 0; j < l; ++j) {
    f[j] = -dy_eig[j] / (d[j] + delta);
  }

  // dy = h * f
  std::vector<double> dy(l, 0.0);
  for (int row = 0; row < l; ++row) {
    double sum = 0.0;
    for (int col = 0; col < l; ++col) {
      sum += h[eco_idx(row, col, l)] * f[col];
    }
    dy[row] = sum;
  }

  return dy;
}

static void eco_newton_fortran_strict(
  const std::vector<double>& f,
  const std::vector<double>& x,
  const std::vector<double>& w,
  std::vector<double>& y,
  EcoObjFortranStrict& obj)
{
  const int l = static_cast<int>(y.size());
  const std::vector<double> dy = eco_shift_fortran_strict(obj.g, obj.h, l);

  double g0 = 0.0;
  for (int j = 0; j < l; ++j) {
    g0 += obj.g[j] * dy[j];
  }
  if (g0 > 0.0) {
    PRINT_INPUT_ERROR("dy is not a descent direction in Eco-PIMD Newton step.");
  }

  const double s0 = obj.s;
  double c = 2.0;

  for (int iter = 1; iter <= 60; ++iter) {
    c *= 0.5;

    std::vector<double> z(l, 0.0);
    for (int j = 0; j < l; ++j) {
      z[j] = y[j] + c * dy[j];
    }

    if (z[0] < 0.0) {
      continue;
    }

    bool ordered = true;
    for (int j = 1; j < l; ++j) {
      if (z[j] < z[j - 1]) {
        ordered = false;
        break;
      }
    }
    if (!ordered) {
      continue;
    }

    EcoObjFortranStrict trial = eco_objfun_fortran_strict(f, x, w, z);
    if (trial.s <= s0) {
      y.swap(z);
      obj = trial;
      return;
    }
  }

  PRINT_INPUT_ERROR("failed to go downhill in Eco-PIMD Newton step.");
}

static PIMDModeBuildResultFortranStrict pimd_make_mode_factors_fortran_strict(
  const int n,
  const int mode,
  const double xmax)
{
  const double pi = acos(-1.0);

  if (n < 2) {
    PRINT_INPUT_ERROR("PIMD number of beads should be >= 2.");
  }

  PIMDModeBuildResultFortranStrict result;
  result.factor.assign(n, 0.0);
  result.rmst = 0.0;
  result.rmsm = 0.0;
  result.rmse = 0.0;
  result.newton_iterations = 0;

  if (mode == PIMD_MODE_TROTTER) {
    for (int k = 1; k < n; ++k) {
      result.factor[k] = 2.0 * sin(static_cast<double>(k) * pi / static_cast<double>(n));
    }
    return result;
  }

  if (mode == PIMD_MODE_MATSUBARA) {
    for (int k = 1; k < n; ++k) {
      const int kk = std::min(k, n - k);
      result.factor[k] = 2.0 * pi * static_cast<double>(kk) / static_cast<double>(n);
    }
    return result;
  }

  if (mode != PIMD_MODE_ECO) {
    PRINT_INPUT_ERROR("Unknown PIMD internal-mode scheme.");
  }

  if (!(xmax > 0.0)) {
    PRINT_INPUT_ERROR("Eco-PIMD requires xmax > 0.");
  }

  const int l = n / 2;
  const int m = std::max(1, static_cast<int>(floor(10.0 * xmax + 0.5)));
  const double dx = xmax / static_cast<double>(m);

  std::vector<double> f(m, 0.0);
  std::vector<double> x(m, 0.0);
  for (int j = 0; j < m; ++j) {
    x[j] = (static_cast<double>(j) + 0.5) * dx;
    const double halfx = 0.5 * x[j];
    f[j] = x[j] * x[j] / (halfx / tanh(halfx) - 1.0);
  }

  std::vector<double> w(l, 2.0);
  if (2 * l == n) {
    w[l - 1] = 1.0;
  }

  std::vector<double> y(l, 0.0);

  // Trotter initial guess and RMS error: y(k)=2*n*sin(k*pi/n)
  for (int k = 0; k < l; ++k) {
    y[k] = 2.0 * static_cast<double>(n) *
      sin(static_cast<double>(k + 1) * pi / static_cast<double>(n));
  }
  EcoObjFortranStrict obj_trotter = eco_objfun_fortran_strict(f, x, w, y);
  result.rmst = sqrt(2.0 * obj_trotter.s);

  // Matsubara initial guess and RMS error: y(k)=2*k*pi
  for (int k = 0; k < l; ++k) {
    y[k] = 2.0 * static_cast<double>(k + 1) * pi;
  }
  EcoObjFortranStrict obj = eco_objfun_fortran_strict(f, x, w, y);
  result.rmsm = sqrt(2.0 * obj.s);
  result.rmse = result.rmsm;

  for (int iter = 1; iter <= 10000; ++iter) {
    const double sp = obj.s;
    eco_newton_fortran_strict(f, x, w, y, obj);
    result.rmse = sqrt(2.0 * obj.s);
    result.newton_iterations = iter;
    if (obj.s >= sp) {
      break;
    }
  }

  // Return y(k)/n for k=1,...,n-1 with y(n-k)=y(k).
  result.factor[0] = 0.0;
  for (int k = 1; k <= l; ++k) {
    result.factor[k] = y[k - 1] / static_cast<double>(n);
  }
  for (int k = 1; k <= (n - 1) / 2; ++k) {
    result.factor[n - k] = result.factor[k];
  }

  return result;
}
// === END ECO_PIMD_PUBLIC_STRICT_FORTRAN_V3 ===
'''


def fail(msg: str) -> None:
    raise RuntimeError(msg)


def backup_once(path: Path) -> None:
    bak = path.with_name(path.name + ".before_eco_pimd_strict_fortran_v3")
    if not bak.exists():
        shutil.copy2(path, bak)


def replace_exact(text: str, old: str, new: str, what: str) -> str:
    count = text.count(old)
    if count != 1:
        fail(f"Expected exactly one occurrence for {what}, found {count}.")
    return text.replace(old, new)


def insert_after_exact(text: str, anchor: str, insertion: str, what: str) -> str:
    count = text.count(anchor)
    if count != 1:
        fail(f"Expected exactly one anchor for {what}, found {count}.")
    return text.replace(anchor, anchor + insertion, 1)


def insert_before_exact(text: str, anchor: str, insertion: str, what: str) -> str:
    count = text.count(anchor)
    if count != 1:
        fail(f"Expected exactly one anchor for {what}, found {count}.")
    return text.replace(anchor, insertion + anchor, 1)


def add_include(text: str, include: str) -> str:
    if include in text:
        return text
    # Insert after the last existing #include line in the initial include block.
    matches = list(re.finditer(r"^#include .*?$", text, flags=re.MULTILINE))
    if not matches:
        fail("Cannot find include block.")
    pos = matches[-1].end()
    return text[:pos] + "\n" + include + text[pos:]


def patch_ensemble_pimd_cuh(path: Path) -> None:
    text = path.read_text()
    if PATCH_MARK in text:
        print(f"Already patched {path}")
        return

    old1 = """  Ensemble_PIMD(\n    int number_of_atoms_input, int number_of_beads_input, bool thermostat_internal, Atom& atom);"""
    new1 = """  Ensemble_PIMD(\n    int number_of_atoms_input,\n    int number_of_beads_input,\n    bool thermostat_internal,\n    Atom& atom,\n    int internal_mode_scheme_input = 0,\n    double eco_omega_max_cm1_input = 0.0);"""
    text = replace_exact(text, old1, new1, "RPMD/TRPMD Ensemble_PIMD constructor declaration")

    old2 = """  Ensemble_PIMD(\n    int number_of_atoms_input, int number_of_beads_input, double temperature_coupling, Atom& atom);"""
    new2 = """  Ensemble_PIMD(\n    int number_of_atoms_input,\n    int number_of_beads_input,\n    double temperature_coupling,\n    Atom& atom,\n    int internal_mode_scheme_input = 0,\n    double eco_omega_max_cm1_input = 0.0);"""
    text = replace_exact(text, old2, new2, "NVT-PIMD Ensemble_PIMD constructor declaration")

    old3 = """    double pressure_coupling[6],\n    Atom& atom);"""
    new3 = """    double pressure_coupling[6],\n    Atom& atom,\n    int internal_mode_scheme_input = 0,\n    double eco_omega_max_cm1_input = 0.0);"""
    text = replace_exact(text, old3, new3, "NPT-PIMD Ensemble_PIMD constructor declaration tail")

    text = insert_after_exact(
        text,
        "  double omega_n;\n",
        "  int internal_mode_scheme = 0; // 0=trotter, 1=matsubara, 2=eco\n"
        "  double eco_omega_max_cm1 = 0.0;\n"
        "  double last_internal_mode_temperature = -1.0;\n"
        "  GPU_Vector<double> internal_mode_factors;\n",
        "Ensemble_PIMD internal-mode members",
    )

    text = insert_after_exact(
        text,
        "  void initialize(Atom& atom);\n",
        "  void update_internal_mode_factors();\n",
        "Ensemble_PIMD update_internal_mode_factors declaration",
    )

    # Header-local constants so integrate.cu can use the same symbolic values.
    text = insert_after_exact(
        text,
        "#include <vector>\n",
        "\n// === BEGIN ECO_PIMD_PUBLIC_STRICT_FORTRAN_V3_DECL ===\n"
        "static const int PIMD_MODE_TROTTER = 0;\n"
        "static const int PIMD_MODE_MATSUBARA = 1;\n"
        "static const int PIMD_MODE_ECO = 2;\n"
        "static inline const char* pimd_mode_name(const int mode)\n"
        "{\n"
        "  if (mode == PIMD_MODE_TROTTER) { return \"trotter\"; }\n"
        "  if (mode == PIMD_MODE_MATSUBARA) { return \"matsubara\"; }\n"
        "  if (mode == PIMD_MODE_ECO) { return \"eco\"; }\n"
        "  return \"unknown\";\n"
        "}\n"
        "// === END ECO_PIMD_PUBLIC_STRICT_FORTRAN_V3_DECL ===\n",
        "PIMD mode constants in ensemble_pimd.cuh",
    )

    backup_once(path)
    path.write_text(text)
    print(f"Patched {path}")


def patch_integrate_cuh(path: Path) -> None:
    text = path.read_text()
    if PATCH_MARK in text:
        print(f"Already patched {path}")
        return

    text = insert_after_exact(
        text,
        "  int number_of_beads;\n",
        "  int pimd_internal_mode_scheme = 0; // 0=trotter, 1=matsubara, 2=eco\n"
        "  double pimd_eco_omega_max_cm1 = 0.0;\n"
        "  int pimd_base_num_param = 0;\n"
        "  // ECO_PIMD_PUBLIC_STRICT_FORTRAN_V3\n",
        "Integrate PIMD internal-mode members",
    )

    backup_once(path)
    path.write_text(text)
    print(f"Patched {path}")


def patch_ensemble_pimd_cu(path: Path) -> None:
    text = path.read_text()
    if PATCH_MARK in text:
        print(f"Already patched {path}")
        return

    for inc in [
        "#include <algorithm>",
        "#include <cmath>",
        "#include <cstdio>",
        "#include <vector>",
    ]:
        text = add_include(text, inc)

    text = insert_before_exact(
        text,
        "void Ensemble_PIMD::initialize_rng()\n",
        HELPER_BLOCK + "\n\n",
        "strict Fortran Eco helper block",
    )

    old_sig1 = """Ensemble_PIMD::Ensemble_PIMD(\n  int number_of_atoms_input, int number_of_beads_input, bool thermostat_internal_input, Atom& atom)\n{\n  number_of_atoms = number_of_atoms_input;\n  number_of_beads = number_of_beads_input;"""
    new_sig1 = """Ensemble_PIMD::Ensemble_PIMD(\n  int number_of_atoms_input,\n  int number_of_beads_input,\n  bool thermostat_internal_input,\n  Atom& atom,\n  int internal_mode_scheme_input,\n  double eco_omega_max_cm1_input)\n{\n  number_of_atoms = number_of_atoms_input;\n  number_of_beads = number_of_beads_input;\n  internal_mode_scheme = internal_mode_scheme_input;\n  eco_omega_max_cm1 = eco_omega_max_cm1_input;"""
    text = replace_exact(text, old_sig1, new_sig1, "RPMD/TRPMD constructor definition")

    old_sig2 = """Ensemble_PIMD::Ensemble_PIMD(\n  int number_of_atoms_input,\n  int number_of_beads_input,\n  double temperature_coupling_input,\n  Atom& atom)\n{\n  number_of_atoms = number_of_atoms_input;\n  number_of_beads = number_of_beads_input;"""
    new_sig2 = """Ensemble_PIMD::Ensemble_PIMD(\n  int number_of_atoms_input,\n  int number_of_beads_input,\n  double temperature_coupling_input,\n  Atom& atom,\n  int internal_mode_scheme_input,\n  double eco_omega_max_cm1_input)\n{\n  number_of_atoms = number_of_atoms_input;\n  number_of_beads = number_of_beads_input;\n  internal_mode_scheme = internal_mode_scheme_input;\n  eco_omega_max_cm1 = eco_omega_max_cm1_input;"""
    text = replace_exact(text, old_sig2, new_sig2, "NVT-PIMD constructor definition")

    old_sig3 = """  double target_pressure_input[6],\n  double pressure_coupling_input[6],\n  Atom& atom)\n{\n  number_of_atoms = number_of_atoms_input;\n  number_of_beads = number_of_beads_input;"""
    new_sig3 = """  double target_pressure_input[6],\n  double pressure_coupling_input[6],\n  Atom& atom,\n  int internal_mode_scheme_input,\n  double eco_omega_max_cm1_input)\n{\n  number_of_atoms = number_of_atoms_input;\n  number_of_beads = number_of_beads_input;\n  internal_mode_scheme = internal_mode_scheme_input;\n  eco_omega_max_cm1 = eco_omega_max_cm1_input;"""
    text = replace_exact(text, old_sig3, new_sig3, "NPT-PIMD constructor definition tail")

    # Allocate and initialize mode factors after the transformation matrix is constructed.
    text = insert_after_exact(
        text,
        "  transformation_matrix.copy_from_host(transformation_matrix_cpu.data());\n",
        "\n  internal_mode_factors.resize(number_of_beads);\n"
        "  std::vector<double> initial_mode_factors =\n"
        "    pimd_make_mode_factors_fortran_strict(number_of_beads, PIMD_MODE_TROTTER, 1.0).factor;\n"
        "  internal_mode_factors.copy_from_host(initial_mode_factors.data());\n",
        "internal_mode_factors initialization",
    )

    update_method = r'''
void Ensemble_PIMD::update_internal_mode_factors()
{
  bool need_update = false;

  if (last_internal_mode_temperature < 0.0) {
    need_update = true;
  }

  if (internal_mode_scheme == PIMD_MODE_ECO) {
    if (temperature <= 0.0) {
      PRINT_INPUT_ERROR("Eco-PIMD requires a positive temperature.");
    }
    if (eco_omega_max_cm1 <= 0.0) {
      PRINT_INPUT_ERROR("Eco-PIMD requires a positive omega_max in cm^-1.");
    }
    if (fabs(temperature - last_internal_mode_temperature) > 1.0e-12) {
      need_update = true;
    }
  }

  if (!need_update) {
    return;
  }

  double xmax = 1.0;
  if (internal_mode_scheme == PIMD_MODE_ECO) {
    xmax = ECO_CM_TO_K * eco_omega_max_cm1 / temperature;
  }

  PIMDModeBuildResultFortranStrict result =
    pimd_make_mode_factors_fortran_strict(number_of_beads, internal_mode_scheme, xmax);

  internal_mode_factors.copy_from_host(result.factor.data());
  last_internal_mode_temperature =
    (internal_mode_scheme == PIMD_MODE_ECO) ? temperature : 0.0;

  if (internal_mode_scheme == PIMD_MODE_ECO) {
    printf(
      "    Eco-PIMD internal modes: omega_max=%g cm^-1, T=%g K, xmax=%g, "
      "RMSE(trotter)=%g, RMSE(matsubara)=%g, RMSE(eco)=%g, Newton iterations=%d.\n",
      eco_omega_max_cm1,
      temperature,
      xmax,
      result.rmst,
      result.rmsm,
      result.rmse,
      result.newton_iterations);
  }
}

'''
    text = insert_before_exact(
        text,
        "Ensemble_PIMD::~Ensemble_PIMD(void)\n",
        update_method,
        "update_internal_mode_factors method implementation",
    )

    # gpu_nve_1 signature and omega_k.
    text = replace_exact(
        text,
        """  const int number_of_beads,\n  const double omega_n,\n  const double time_step,""",
        """  const int number_of_beads,\n  const double omega_n,\n  const double* internal_mode_factors,\n  const double time_step,""",
        "gpu_nve_1 signature add internal_mode_factors",
    )
    text = replace_exact(
        text,
        "      double omega_k = 2.0 * omega_n * sin(k * PI / number_of_beads);",
        "      double omega_k = omega_n * internal_mode_factors[k];",
        "gpu_nve_1 omega_k replacement",
    )

    # gpu_langevin signature and damping.
    text = replace_exact(
        text,
        """  const double temperature_coupling,\n  const double omega_n,\n  const double time_step,""",
        """  const double temperature_coupling,\n  const double omega_n,\n  const double* internal_mode_factors,\n  const double time_step,""",
        "gpu_langevin signature add internal_mode_factors",
    )
    text = replace_exact(
        text,
        """      double c1 = (k == 0) ? exp(-0.5 / temperature_coupling)\n                           : exp(-time_step * omega_n * sin(k * PI / number_of_beads));""",
        """      double c1 = (k == 0) ? exp(-0.5 / temperature_coupling)\n                           : exp(-0.5 * time_step * omega_n * internal_mode_factors[k]);""",
        "gpu_langevin damping replacement",
    )

    # Kernel call sites and update calls.
    # There are two identical omega_n assignments, one in compute1 and one in compute2.
    anchor = "  omega_n = number_of_beads * K_B * temperature / HBAR;\n\n"
    if text.count(anchor) != 2:
        fail(f"Expected exactly two omega_n anchors in compute1/compute2, found {text.count(anchor)}.")
    text = text.replace(anchor, anchor + "  update_internal_mode_factors();\n\n", 2)

    text = replace_exact(
        text,
        """      omega_n,\n      time_step,\n      transformation_matrix.data(),""",
        """      omega_n,\n      internal_mode_factors.data(),\n      time_step,\n      transformation_matrix.data(),""",
        "gpu_langevin call add internal_mode_factors",
    )
    text = replace_exact(
        text,
        """    omega_n,\n    time_step,\n    transformation_matrix.data(),""",
        """    omega_n,\n    internal_mode_factors.data(),\n    time_step,\n    transformation_matrix.data(),""",
        "gpu_nve_1 call add internal_mode_factors",
    )

    backup_once(path)
    path.write_text(text)
    print(f"Patched {path}")


def patch_integrate_cu(path: Path, max_beads: int | None) -> None:
    text = path.read_text()
    if PATCH_MARK in text:
        print(f"Already patched {path}")
        return

    # Pass parsed PIMD mode settings into Ensemble_PIMD constructors.
    text = replace_exact(
        text,
        """        ensemble.reset(\n          new Ensemble_PIMD(number_of_atoms, number_of_beads, temperature_coupling, atom));""",
        """        ensemble.reset(new Ensemble_PIMD(\n          number_of_atoms,\n          number_of_beads,\n          temperature_coupling,\n          atom,\n          pimd_internal_mode_scheme,\n          pimd_eco_omega_max_cm1));""",
        "NVT-PIMD constructor call",
    )

    text = replace_exact(
        text,
        """          target_pressure,\n          pressure_coupling,\n          atom));""",
        """          target_pressure,\n          pressure_coupling,\n          atom,\n          pimd_internal_mode_scheme,\n          pimd_eco_omega_max_cm1));""",
        "NPT-PIMD constructor call",
    )

    old_pimd_branch = """  } else if (strcmp(param[1], "pimd") == 0) {\n    type = 33;\n    if (num_param != 6 && num_param != 9 && num_param != 13 && num_param != 19) {\n      PRINT_INPUT_ERROR("ensemble pimd should have 4 or 7 or 11 or 17 parameters.");\n    }"""
    new_pimd_branch = """  } else if (strcmp(param[1], "pimd") == 0) {\n    type = 33;\n    pimd_internal_mode_scheme = PIMD_MODE_TROTTER;\n    pimd_eco_omega_max_cm1 = 0.0;\n    pimd_base_num_param = num_param;\n\n    // Optional final arguments for the PIMD internal-mode spectrum:\n    //   trotter\n    //   matsubara\n    //   eco omega_max_cm1\n    if (num_param >= 4 && strcmp(param[num_param - 2], "eco") == 0) {\n      pimd_internal_mode_scheme = PIMD_MODE_ECO;\n      pimd_base_num_param = num_param - 2;\n      if (!is_valid_real(param[num_param - 1], &pimd_eco_omega_max_cm1)) {\n        PRINT_INPUT_ERROR("Eco-PIMD omega_max should be a number in cm^-1.");\n      }\n      if (pimd_eco_omega_max_cm1 <= 0.0) {\n        PRINT_INPUT_ERROR("Eco-PIMD omega_max should > 0.");\n      }\n    } else if (num_param >= 3 && strcmp(param[num_param - 1], "trotter") == 0) {\n      pimd_internal_mode_scheme = PIMD_MODE_TROTTER;\n      pimd_base_num_param = num_param - 1;\n    } else if (num_param >= 3 && strcmp(param[num_param - 1], "matsubara") == 0) {\n      pimd_internal_mode_scheme = PIMD_MODE_MATSUBARA;\n      pimd_base_num_param = num_param - 1;\n    }\n\n    if (\n      pimd_base_num_param != 6 && pimd_base_num_param != 9 && pimd_base_num_param != 13 &&\n      pimd_base_num_param != 19) {\n      PRINT_INPUT_ERROR(\n        "ensemble pimd should have 4, 7, 11, or 17 parameters, optionally followed by "\n        "trotter, matsubara, or eco omega_max_cm1.");\n    }\n    // ECO_PIMD_PUBLIC_STRICT_FORTRAN_V3"""
    text = replace_exact(text, old_pimd_branch, new_pimd_branch, "PIMD parser branch")

    # Inside the PIMD pressure parser use the base count excluding mode arguments.
    text = replace_exact(
        text,
        "      if (num_param >= 9) {\n        if (num_param == 13) {",
        "      if (pimd_base_num_param >= 9) {\n        if (pimd_base_num_param == 13) {",
        "PIMD pressure parser base-count start",
    )
    text = replace_exact(
        text,
        "        } else if (num_param == 9) { // isotropic",
        "        } else if (pimd_base_num_param == 9) { // isotropic",
        "PIMD isotropic pressure base-count",
    )

    # Printed information should use base PIMD argument count.
    text = replace_exact(
        text,
        "      if (num_param >= 9) {\n        printf(\"Use NPT-PIMD for this run.\\n\");",
        "      if (pimd_base_num_param >= 9) {\n        printf(\"Use NPT-PIMD for this run.\\n\");",
        "PIMD print NPT/NVT base-count",
    )
    text = replace_exact(
        text,
        "      if (num_param >= 9) {\n        if (num_target_pressure_components == 1) {",
        "      if (pimd_base_num_param >= 9) {\n        if (num_target_pressure_components == 1) {",
        "PIMD print pressure base-count",
    )

    old_print_block = """      printf(\"    number of beads is %d.\\n\", number_of_beads);
      printf(\"    initial temperature is %g K.\\n\", temperature1);
      printf(\"    final temperature is %g K.\\n\", temperature2);
      printf(\"    tau_T is %g time_step.\\n\", temperature_coupling);"""
    new_print_block = old_print_block + """
      printf(\"    PIMD internal mode scheme is %s.\\n\", pimd_mode_name(pimd_internal_mode_scheme));
      if (pimd_internal_mode_scheme == PIMD_MODE_ECO) {
        printf(\"    Eco-PIMD omega_max is %g cm^-1.\\n\", pimd_eco_omega_max_cm1);
      }"""
    text = replace_exact(text, old_print_block, new_print_block, "PIMD mode print block")

    if max_beads is not None:
        text = text.replace(
            'PRINT_INPUT_ERROR("number of beads should <= 128.");',
            f'PRINT_INPUT_ERROR("number of beads should <= {max_beads}.");',
        )

    backup_once(path)
    path.write_text(text)
    print(f"Patched {path}")


def patch_max_num_beads(root: Path, max_beads: int | None) -> None:
    if max_beads is None:
        return

    found = []
    for path in (root / "src").rglob("*"):
        if not path.is_file():
            continue
        if path.suffix not in {".cuh", ".cu", ".h", ".hpp"}:
            continue
        text = path.read_text(errors="ignore")
        new_text, n = re.subn(
            r"(#\s*define\s+MAX_NUM_BEADS\s+)(\d+)",
            rf"\g<1>{max_beads}",
            text,
        )
        if n > 0:
            backup_once(path)
            path.write_text(new_text)
            found.append(str(path))

    if found:
        print("Updated MAX_NUM_BEADS in:")
        for p in found:
            print("  " + p)
    else:
        print("WARNING: did not find '#define MAX_NUM_BEADS <number>' under src/.")
        print("         If compilation still says MAX_NUM_BEADS=128, edit that definition manually.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Patch clean GPUMD-4.6 source for Trotter/Matsubara/Eco PIMD internal modes."
    )
    parser.add_argument("gpumd_root", help="Path to clean GPUMD source root, e.g. /path/to/GPUMD-4.6")
    parser.add_argument(
        "--max-beads",
        type=int,
        default=None,
        help="Optionally replace MAX_NUM_BEADS and the parser error message, e.g. 256 or 512.",
    )
    args = parser.parse_args()

    root = Path(args.gpumd_root).resolve()
    if not root.exists():
        fail(f"GPUMD root does not exist: {root}")

    ensemble_cuh = root / "src" / "integrate" / "ensemble_pimd.cuh"
    ensemble_cu = root / "src" / "integrate" / "ensemble_pimd.cu"
    integrate_cuh = root / "src" / "integrate" / "integrate.cuh"
    integrate_cu = root / "src" / "integrate" / "integrate.cu"

    for path in [ensemble_cuh, ensemble_cu, integrate_cuh, integrate_cu]:
        if not path.exists():
            fail(f"Required file not found: {path}")

    print("Patching GPUMD source root:", root)
    patch_ensemble_pimd_cuh(ensemble_cuh)
    patch_integrate_cuh(integrate_cuh)
    patch_ensemble_pimd_cu(ensemble_cu)
    patch_integrate_cu(integrate_cu, args.max_beads)
    patch_max_num_beads(root, args.max_beads)

    print("\nPatch completed.")
    print("New PIMD syntax examples:")
    print("  ensemble pimd P T1 T2 tau_T")
    print("  ensemble pimd P T1 T2 tau_T trotter")
    print("  ensemble pimd P T1 T2 tau_T matsubara")
    print("  ensemble pimd P T1 T2 tau_T eco 3500")
    print("For NPT-PIMD, append the same optional mode at the end of the original command.")
    print("\nImportant: this patch does not edit src/makefile.")
    print("The Eco optimizer calls LAPACK dsyev_, so you must link LAPACK/OpenBLAS when compiling.")


if __name__ == "__main__":
    main()

