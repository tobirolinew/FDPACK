# Installation

## Requirements

- Fortran compiler: **Intel ifort**, **ifx**, or **gfortran**
- GNU Make
- *(optional)* LAPACK — required only for `diag_num_hess`

## Tested Compilers

| Compiler   | Version                          |
|------------|----------------------------------|
| `gfortran` | GCC 11.5.0 (Red Hat 11.5.0-5)   |
| `ifort`    | Intel oneAPI 2024.2.1            |
| `ifx`      | Intel oneAPI 2024.2.1            |

> **Note:** `ifort` emits deprecation warning 10448.
>
> The `Makefile` suppresses it automatically with
> 
> `-diag-disable=10448`.

## Building

### 1. Select a compiler

Uncomment the desired compiler in the `Makefile`:
```makefile
FC = ifort    # Intel Fortran (classic)
#FC = ifx     # Intel Fortran (LLVM-based)
#FC = gfortran
```

### 2. Select a build mode
```makefile
MODE = release   # optimized (-O3 / -O2)
#MODE = debug    # runtime checks and traceback
```

### 3. Build
```bash
make
```

Compiles all modules and test programs.<br>
The build system uses stamp files (`.cmp_stamp`, `.dist_stamp`) <br>
to avoid redundant work on subsequent invocations.

On success, the distributable library is written to `dist/$(TC)`:
```
dist/
├── libfdpack.a
├── mod_fdbase.mod
└── mod_fdtune.mod
```

## Running the Tests
```bash
make test
```

Presents an interactive menu to run individual test <br> 
programs or all of them at once. <br>
Executables are in `bin/`:
```
bin/
├── test_demo.x
├── test_func_noise.x
├── test_grad_offsets.x
└── test_hess_offsets.x
```

## Optional: LAPACK Support

LAPACK is required for `diag_num_hess`. <br>
Tested with **LAPACK 3.12.0**. <br>
Set the following in the `Makefile` before building:
```makefile
USE_LAPACK  = 1
LAPACK_PATH = /path/to/lapack-3.12.0/lib/$(TC)
```

`$(TC)` is set automatically by the `Makefile` based on <br> 
the selected compiler: `gnu` for gfortran, `intel` for ifort/ifx. <br>
The build expects the libraries at:
```
$(LAPACK_PATH)/liblapack.a
$(LAPACK_PATH)/librefblas.a
```

## Directory Structure

| Directory    | Contents                               |
|--------------|----------------------------------------|
| `src/fdlib/` | core library source files              |
| `src/iface/` | interface include files                |
| `src/test/`  | test program source files              |
| `mk/`        | Makefile includes                      |
| `build/`     | object files, `.mod` files, stamps     |
| `bin/`       | compiled test executables (`.x`)       |
| `dist/`      | distributable library and `.mod` files |
| `doc/`       | documentation                          |

## Cleaning
```bash
make clean
```

Removes all files under `build/`, `bin/`, and `dist/`.
