# Eco-PIMD patch for GPUMD

This repository provides a lightweight source-code patch for adding **Eco-PIMD internal-mode frequencies** to the GPUMD implementation of path-integral molecular dynamics.

The main goal is to accelerate the convergence of PIMD simulations by replacing the primitive Trotter ring-polymer internal-mode frequencies with optimized Eco frequencies. In practice, this can reduce the number of beads required to converge quantum nuclear effects, especially for systems containing high-frequency light-atom vibrations such as H-containing materials, molecular crystals, water, ice, perovskites, and metal-organic frameworks.

<p align="center">
  <img src="docs/eco_pimd_concept.png" width="850">
</p>

## 1. Motivation

In conventional primitive PIMD, each quantum particle is represented by a ring polymer with \(P\) beads. The internal normal-mode frequencies are

\[
\Omega_k^{\rm Trotter}
=
2\omega_P \sin\left(\frac{k\pi}{P}\right),
\qquad
k=1,\ldots,P-1,
\]

where

\[
\omega_P = \frac{Pk_{\rm B}T}{\hbar}.
\]

Large bead numbers are often required when the system contains high-frequency vibrational modes. For example, hydrogen-containing systems may require many beads at low temperature because the dimensionless quantum frequency

\[
x=\beta\hbar\omega
\]

can be large.

Eco-PIMD changes the internal-mode frequency spectrum while leaving the physical potential, bead force evaluation, normal-mode transformation, and GPUMD PIMD structure unchanged. The Eco internal frequencies are optimized to reproduce the quantum harmonic fluctuation over a chosen frequency range with fewer beads.

The patched GPUMD supports three internal-mode choices:

\[
\lambda_k =
\frac{\Omega_k}{\omega_P}
=
\begin{cases}
2\sin(k\pi/P), & \text{Trotter}, \\
2\pi\min(k,P-k)/P, & \text{Matsubara}, \\
y_k/P, & \text{Eco}.
\end{cases}
\]

Here \(y_k\) are optimized dimensionless Eco frequencies.

## 2. What this patch changes

The patch modifies only the PIMD internal-mode frequency generation and the associated internal-mode Langevin thermostat factors.

It does **not** modify:

- NEP or other interatomic potentials;
- bead coordinate, velocity, force, or virial storage;
- physical force evaluation;
- normal-mode transformation;
- pressure control;
- kinetic-energy or virial estimators;
- thermodynamic output format.

The original GPUMD PIMD behavior is preserved as the default.

## 3. User-facing syntax

After patching, the original GPUMD syntax remains valid:

```text
ensemble pimd P T1 T2 tau_T
````

This defaults to primitive Trotter PIMD.

The patched version also accepts explicit internal-mode choices:

```text
ensemble pimd P T1 T2 tau_T trotter
ensemble pimd P T1 T2 tau_T matsubara
ensemble pimd P T1 T2 tau_T eco omega_max_cm1
```

For NPT-PIMD, append the internal-mode option at the end of the original pressure-control arguments. For example:

```text
ensemble pimd 160 60 60 100 0.0 100.0 1000 eco 3500
```

This means:

* bead number: (P=160);
* temperature: (T=60) K;
* internal-mode scheme: Eco;
* physical maximum frequency cutoff: (\omega_{\max}=3500~{\rm cm}^{-1}).

## 4. Meaning of `omega_max_cm1`

The Eco input parameter is a **wavenumber** in ({\rm cm}^{-1}), not an angular frequency.

The dimensionless fitting range is

[
x_{\max}
========

# \beta\hbar\omega_{\max}

\frac{hc}{k_{\rm B}}
\frac{\tilde{\nu}_{\max}}{T},
]

where (\tilde{\nu}_{\max}) is the maximum vibrational wavenumber in ({\rm cm}^{-1}). Numerically,

[
x_{\max}
========

1.438776877
\frac{\tilde{\nu}_{\max}^{\rm cm^{-1}}}{T}.
]

For example,

[
3500~{\rm cm}^{-1}
\approx 105~{\rm THz}
]

as a linear frequency.

For physical simulations, use a fixed physical cutoff (\omega_{\max}^{\rm cm^{-1}}) for the material across all temperatures. The corresponding (x_{\max}) naturally changes with temperature.

## 5. Eco frequency optimization

The Eco frequency optimization follows the supplied Fortran reference implementation.

For a chosen (P) and (x_{\max}), a grid is constructed as

[
m = {\rm nint}(10x_{\max}),
]

[
x_j = \left(j-\frac{1}{2}\right)\Delta x,
\qquad
\Delta x = \frac{x_{\max}}{m}.
]

The objective minimizes the fractional error in the harmonic quantum fluctuation function:

[
r_j
===

f_j
\sum_{k=1}^{l}
\frac{w_k}{x_j^2+y_k^2}
-1,
]

where

[
f_j
===

\frac{x_j^2}
{
(x_j/2)/\tanh(x_j/2)-1
}.
]

The objective is

[
s=
\frac{1}{2m}
\sum_{j=1}^{m} r_j^2,
]

and the reported RMSE is

[
{\rm RMSE} = \sqrt{2s}.
]

Only the independent half of the internal modes is optimized:

[
l = \left\lfloor \frac{P}{2}\right\rfloor.
]

The mode weights are

[
w_k=2,
]

except for the Nyquist mode when (P) is even:

[
w_{P/2}=1.
]

The full spectrum is reconstructed using

[
y_{P-k}=y_k.
]

## 6. Numerical optimizer

The Eco optimizer uses the same shifted-Hessian Newton procedure as the Fortran code.

The Hessian is diagonalized as

[
H = UDU^{\rm T}.
]

The Hessian shift is

[
\delta
======

\max\left(10^{-16}d_{\max}, -2d_{\min}\right).
]

The Newton step is

[
\Delta y
========

-U(D+\delta I)^{-1}U^{\rm T}g.
]

A backtracking line search is then used. A trial step is accepted only if:

[
z_1 \ge 0,
]

[
z_k \ge z_{k-1},
]

and

[
s(z) \le s(y).
]

Because the optimizer diagonalizes the Hessian using LAPACK `DSYEV`, the patched GPUMD must be linked with a LAPACK/BLAS backend, such as OpenBLAS.

## 7. Requirements

### Python patcher

The patcher itself only requires standard Python:

```bash
python3
```

No special Python packages are required for the patching step.

### GPUMD compilation

The patched GPUMD source requires:

* CUDA;
* a C++14-capable compiler;
* LAPACK/BLAS for `dsyev_`.

On many HPC systems, OpenBLAS provides the needed LAPACK symbol.

Example module setup:

```bash
module purge
module load CUDA/10.1.243-GCC-8.3.0
module load OpenBLAS/0.3.7-GCC-8.3.0
```

Loading the OpenBLAS module makes the library available, but the final GPUMD link command must still explicitly link it, for example:

```make
-lopenblas
```

If compilation fails with

```text
undefined reference to `dsyev_'
```

then LAPACK/OpenBLAS is not linked in the final GPUMD executable.

## 8. How to apply the patch

Start from a clean GPUMD 4.6 source tree.

```bash
cd /path/to/software
cp -r GPUMD-4.6 GPUMD-4.6-ecopimd
cd GPUMD-4.6-ecopimd
```

Copy the patcher into the GPUMD source folder:

```bash
cp /path/to/patch_gpumd46_pimd_modes_strict_fortran.py .
```

Run the patcher:

```bash
python3 patch_gpumd46_pimd_modes_strict_fortran.py . --max-beads 512
```

The `--max-beads` option is optional. It increases the compile-time maximum bead number in GPUMD. For example, `--max-beads 512` allows PIMD runs with up to 512 beads.

The patcher modifies:

```text
src/integrate/ensemble_pimd.cu
src/integrate/ensemble_pimd.cuh
src/integrate/integrate.cu
src/integrate/integrate.cuh
```

Backup files are written before modification.

## 9. Compilation

After patching, compile GPUMD as usual.

```bash
cd src
make clean
make -j
```

If the final link step fails with

```text
undefined reference to `dsyev_'
```

edit the GPUMD `src/makefile` and add OpenBLAS or LAPACK/BLAS to the final link line. For example:

```make
LIBS = -lcublas -lcusolver -lcufft -lopenblas
```

or, depending on the system,

```make
LIBS = -lcublas -lcusolver -lcufft -llapack -lblas
```

Then recompile:

```bash
make clean
make -j
```

## 10. Quick test

A minimal Eco-PIMD test can be run with:

```text
potential        nep.txt
velocity         100
ensemble         pimd 16 100 100 100 eco 3500
time_step        0.5
dump_thermo      1
run              1
```

For a successful Eco setup, GPUMD should print a line similar to:

```text
Eco-PIMD internal modes: omega_max=3500 cm^-1, T=100 K, xmax=50.3572,
RMSE(trotter)=..., RMSE(matsubara)=..., RMSE(eco)=..., Newton iterations=...
```

The RMSE values depend only on (P) and (x_{\max}), not on the material or system size.

## 11. RMSE check

For a fixed (x_{\max}), Eco should give a much smaller RMSE than primitive Trotter PIMD at the same bead number.

A typical trend is

[
{\rm RMSE}*{\rm Eco}
\ll
{\rm RMSE}*{\rm Trotter}
<
{\rm RMSE}_{\rm Matsubara}.
]

The finite Matsubara curve can be worse than Trotter because a finite Matsubara representation is a hard truncation of the infinite Matsubara sum, whereas primitive Trotter PIMD is a finite imaginary-time discretization of the path-integral partition function.

When plotting RMSE versus bead number, make sure all data points use the same (x_{\max}). Mixing different temperatures or different (\omega_{\max}) values can create artificial non-monotonic behavior.

## 12. Example input

For MOF-5 at low temperature, using a maximum frequency around (3500~{\rm cm}^{-1}):

```text
potential        nep.txt
velocity         60

ensemble         pimd 160 60 60 100 0.0 100.0 1000 eco 3500

time_step        0.5
dump_thermo      200
run              100000
```

For primitive Trotter PIMD with the same bead number:

```text
ensemble         pimd 160 60 60 100 0.0 100.0 1000 trotter
```

or simply omit the final keyword:

```text
ensemble         pimd 160 60 60 100 0.0 100.0 1000
```

## 13. Notes on physical consistency

For physical simulations over a temperature range, keep the physical cutoff (\omega_{\max}^{\rm cm^{-1}}) fixed.

For example, use:

```text
eco 3500
```

at 60 K, 70 K, 80 K, 90 K, and 100 K.

Then

[
x_{\max}=\beta\hbar\omega_{\max}
]

naturally changes with temperature. This is physically correct because the same vibrational mode becomes more quantum mechanical at lower temperature.

Do not force the same (x_{\max}) across different temperatures for production calculations, because that would require changing the physical frequency cutoff with temperature.

## 14. Citation

If you use this patch in published work, please cite the relevant GPUMD and Eco-PIMD references, and cite this repository.

A suggested software citation is:

```text
Zezhu Zeng, eco-pimd-patch-to-GPUMD, GitHub repository,
https://github.com/ZengZezhu/eco-pimd-patch-to-GPUMD
```

## 15. Status

This patch is under active development. Current features include:

* Trotter internal-mode frequencies;
* Matsubara internal-mode frequencies;
* Eco optimized internal-mode frequencies;
* Fortran-consistent Eco optimizer;
* RMSE reporting for Trotter, Matsubara, and Eco;
* optional compile-time maximum bead-number update.

Planned improvements include:

* cleaner makefile integration for LAPACK/OpenBLAS;
* more examples;
* automated RMSE validation scripts;
* benchmark data for representative materials.
