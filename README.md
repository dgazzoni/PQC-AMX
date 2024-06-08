This is the source code accompanying the paper "PQC-AMX: Accelerating Saber and FrodoKEM on the Apple M1 and M3 SoCs". This is a list of folders and their contents:

- `aes`: Implementations of routines related to FrodoKEM-AES "A" matrix generation, using ARMv8's Cryptographic Extensions instructions.
- `amx`: [AMX macros from Peter Cawley](https://github.com/corsix/amx) (`aarch64.h`), plus helper macros of our own (`amx.h`) and implementations of routines described in our paper;
- `googletest`: a copy of the [Google Test](https://github.com/google/googletest/) library;
- `neon-ntt`: relevant files from the paper ["Neon NTT: Faster Dilithium, Kyber, and Saber on Cortex-A72 and Apple M1"](https://eprint.iacr.org/2021/986), obtained from the associated [GitHub repository](https://github.com/neon-ntt/neon-ntt), with some modifications as described in our paper. Also includes our AMX implementation of Saber;
- `PQCrypto-LWEKE`: relevant files from the reference and optimized implementations of [FrodoKEM](https://frodokem.org), obtained from the associated [GitHub repository](https://github.com/microsoft/PQCrypto-LWEKE), with some modifications as described in our paper. Also includes our NEON and AMX implementations of FrodoKEM;
- `rng_opt`: an optimized implementation of the NIST `randombytes` routines, based on AES-256 CTR-DRBG, implemented using AES instructions available in the ARMv8-A Cryptographic Extensions;
- `speed`: benchmarking harnesses for our implementations, constant-time experiment for the "genlut" instruction and our AES-ECB routines for FrodoKEM.
- `speed_results_M1`: raw benchmark results in the M1, and an Excel spreadsheet compiling them, generated using a Python script discussed [below](#Helper-scripts-for-benchmarking-and-constant-time-experiments);
- `speed_results_M3`: raw benchmark results in the M3, and an Excel spreadsheet compiling them, generated using a Python script discussed [below](#Helper-scripts-for-benchmarking-and-constant-time-experiments);
- `test`: tests (using the [Google Test](https://github.com/google/googletest/) library) to validate various aspects of the implementation;

The root folder also includes some files of note:

- `run_benchmarks.sh`: a helper script to run benchmarks (see instructions [below](#Helper-scripts-for-benchmarking-and-constant-time-experiments));
- `consolidate_benchmarks.py`: a Python 3 script to consolidate benchmark results from different systems into a single Microsoft Excel file.

# Building the code

We use [CMake](https://cmake.org) as our build system. It can be installed using [Homebrew](https://brew.sh) with the command `brew install cmake`. A typical sequence of commands to build the code, starting from the root folder of the repository, is:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

**NOTE**: for the tested compilers, there is a register allocation issue when the optimized `randombytes` routine is compiled in Debug mode (i.e. passing `-DCMAKE_BUILD_TYPE=Debug` to CMake), and the build fails. However, in RelWithDebInfo and Release mode, there is no issue.

# Running tests

Compilation produces many test binaries in the build folder (`build/test_*` if using the directions in [Building the code](#building-the-code) above). While it is possible to run each binary directly, we recommend using the `ctest` utility from CMake to run all available tests with a single invocation. `ctest` also runs additional tests that automate the process of comparing KATs using the `PQCgenKAT_kem_*` binaries.

# Running benchmarks

Compilation produces many benchmarking binaries in the build folder (`build/speed_*` if using the directions in [Building the code](#building-the-code) above). Each binary may be run directly, or a full benchmark set can be automatically run using the helper scripts described in [Benchmarking helper scripts](#benchmarking-helper-scripts) below.

Note that binaries must be run with `sudo` to allow access the cycle counter.

# Helper scripts for benchmarking and constant-time experiment for the "genlut" instruction

We provide a helper script to automatically run all available benchmarks (except for those related to the RNG), in the form of `run_benchmarks.sh`. It must be run from the root folder of the repository, and places their results in a folder called `speed_results_Mx`, where `Mx` will be replaced by the CPU name in the machine where the script is run, e.g. `M1`, `M2` or `M3`; this is obtained from `sysctl -n machdep.cpu.brand_string`. Each executable file that is run creates an associated text file containing the benchmark results, with a self-explanatory naming scheme.

A Python 3 script, `consolidate_benchmarks.py`, can be run afterwards (also from the root folder of the repository). It requires the [`pandas`](https://pandas.pydata.org) and [`xlsxwriter`](https://pypi.org/project/XlsxWriter/) packages, which can be installed using `pip`.

This script reads all results from the `speed_results_Mx` folder and generates a Microsoft Excel file displaying them in a tabular form, in a format suitable for comparison with the results presented in our paper.

# License

Our work builds upon many other libraries and implementations, with different licenses for each. Any modifications that we make to an existing work is released under the same original license as that work. As for our original code, we release it under the [Creative Commons CC0 1.0 Universal (CC0 1.0)
Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/).
