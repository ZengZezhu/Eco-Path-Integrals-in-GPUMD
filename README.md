# Eco-PIMD patch for GPUMD

This repository provides a source-code patch for adding **Eco-PIMD internal-mode frequencies** to the GPUMD implementation of path-integral molecular dynamics.

Eco-PIMD accelerates PIMD convergence by replacing the primitive Trotter ring-polymer internal-mode frequencies with optimized frequencies. The goal is to reproduce quantum harmonic fluctuations with fewer beads, thereby reducing the computational cost of PIMD simulations, especially for systems containing high-frequency light-atom vibrations.

<p align="center">
  <img src="docs/eco_pimd_concept.png" width="800">
</p>

## 1. Basic idea

In conventional PIMD, the primitive internal-mode frequencies are

$
\Omega_k^{\mathrm{Trotter}}
===========================

2\omega_P\sin\left(\frac{k\pi}{P}\right),
\qquad
\omega_P=\frac{Pk_{\mathrm{B}}T}{\hbar}.
$

Eco-PIMD replaces these primitive frequencies by optimized frequencies,

$$
\Omega_k^{\mathrm{Eco}}
=======================

\omega_P\lambda_k^{\mathrm{Eco}}.
$$

The optimized factors $\lambda_k^{\mathrm{Eco}}$ are chosen to reproduce the exact quantum harmonic fluctuation over a selected frequency range. Therefore, fewer beads can often give the same converged result.

The patched GPUMD supports three internal-mode choices:

```text
trotter
matsubara
eco
```

The default is still `trotter`, so the original GPUMD behavior is preserved unless another option is explicitly used.

## 2. Requirements

### Python

The patcher only needs standard Python:

```bash
python3
```

No additional Python package is required for applying the patch.

### GPUMD compilation

The patched GPUMD requires:

```text
CUDA
C++14-compatible compiler
LAPACK/BLAS
```

LAPACK/BLAS is needed because the Eco optimizer diagonalizes a Hessian using `dsyev_`.

On many HPC systems, OpenBLAS is sufficient. For example:

```bash
module load CUDA
module load OpenBLAS
```

Loading OpenBLAS makes the library available, but the final GPUMD link command must still explicitly link it. If compilation fails with

```text
undefined reference to `dsyev_'
```

then LAPACK/OpenBLAS is not linked. Add the proper LAPACK/BLAS library to the GPUMD `src/makefile`, for example:

```make
-lopenblas
```

or, depending on your system,

```make
-llapack -lblas
```

## 3. How to apply the patch

Start from a clean GPUMD 4.6 source tree:

```bash
cp -r GPUMD-4.6 GPUMD-4.6-ecopimd
cd GPUMD-4.6-ecopimd
```

Copy the patcher into this folder:

```bash
cp /path/to/patch_gpumd46_pimd_modes_strict_fortran.py .
```

Apply the patch:

```bash
python3 patch_gpumd46_pimd_modes_strict_fortran.py . --max-beads 512
```

The `--max-beads` option is optional. It changes the compile-time maximum bead number allowed by GPUMD.

The patch modifies the GPUMD PIMD source files under:

```text
src/integrate/
```

Backups of the original files are automatically generated.

## 4. Compile GPUMD

After patching, compile GPUMD as usual:

```bash
cd src
make clean
make -j
```

If LAPACK/OpenBLAS is linked correctly, the executable `gpumd` should be generated.

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

## 5. How to run Eco-PIMD

The original GPUMD PIMD syntax still works:

```text
ensemble pimd P T1 T2 tau_T
```

This uses the original Trotter internal-mode frequencies.

To explicitly choose the internal-mode scheme:

```text
ensemble pimd P T1 T2 tau_T trotter
ensemble pimd P T1 T2 tau_T matsubara
ensemble pimd P T1 T2 tau_T eco omega_max_cm1
```

For example:

```text
ensemble pimd 160 100 100 100 eco 3500
```

For NPT-PIMD, append the option after the pressure-control parameters. For example:

```text
ensemble pimd 160 60 60 100 0.0 100.0 1000 eco 3500
```

Here `3500` means

$$
\omega_{\max}=3500~\mathrm{cm}^{-1}.
$$

This is a wavenumber cutoff, not an angular frequency.

The dimensionless fitting range is

$$
x_{\max}
========

# \beta\hbar\omega_{\max}

1.438776877
\frac{\omega_{\max}^{\mathrm{cm}^{-1}}}{T}.
$$

For physical simulations over a temperature range, keep the same physical cutoff, for example `eco 3500`, at all temperatures.

## 6. Example input

A minimal Eco-PIMD test:

```text
potential        nep.txt
velocity         100

ensemble         pimd 16 100 100 100 eco 3500

time_step        0.5
dump_thermo      1
run              1
```

A production-style NPT Eco-PIMD input:

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

## 7. RMSE output

When using the `eco` option, GPUMD prints the fitting errors:

```text
RMSE(trotter)=...
RMSE(matsubara)=...
RMSE(eco)=...
```

These RMSE values depend only on

$$
P
$$

and

$$
x_{\max}=\beta\hbar\omega_{\max}.
$$

They do not depend on the material, potential, cell size, or atom number.

For a fixed $x_{\max}$, Eco should give a much smaller RMSE than primitive Trotter PIMD at the same bead number.

When plotting RMSE versus bead number, make sure all data points use the same $x_{\max}$. Mixing different temperatures or different $\omega_{\max}$ values can create artificial non-monotonic behavior.

## 8. Notes on `omega_max_cm1`

The Eco input cutoff is a wavenumber in $\mathrm{cm}^{-1}$.

For example,

$$
3500~\mathrm{cm}^{-1}
\approx
105~\mathrm{THz}
$$

as a linear frequency.

There is no extra $2\pi$ factor in the conversion because

$$
\hbar\omega = hf = hc\tilde{\nu}.
$$

For a material with maximum vibrational frequency around $100~\mathrm{THz}$, a cutoff around

$$
3300\text{--}3500~\mathrm{cm}^{-1}
$$

is a reasonable starting choice.

## 9. Citation

If you use this patch, please cite this repository:

```text
Zezhu Zeng, eco-pimd-patch-to-GPUMD,
https://github.com/ZengZezhu/eco-pimd-patch-to-GPUMD
```

## 10. Status

This patch is under active development. Current features include:

```text
Trotter internal-mode frequencies
Matsubara internal-mode frequencies
Eco optimized internal-mode frequencies
Fortran-consistent Eco optimizer
RMSE reporting for Trotter, Matsubara, and Eco
optional compile-time maximum bead-number update
```
