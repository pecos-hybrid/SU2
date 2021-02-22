Using the Meson Build System
============================

This is a modified version of the build instructions available [on the SU2 website](https://su2code.github.io/docs_v7/Build-SU2-Linux-MacOS/). Most of the content is the same.  There are two main types of changes:

1. A better description of the meson build system, and how to properly do an out-of-source build.
2. Description of options that differ between the two branches.

Note that the following guide works only on Linux/MacOS and on Windows using Cygwin or the [Linux Subsystem](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

---

- [Quick Compilation Guide](#quick-compilation-guide)
- [Requirements](#requirements)
  - [Compilers](#compilers)
  - [MPI](#mpi)
  - [Python](#python)
  - [Optional: swig and mpi4py](#optional-swig-and-mpi4py)
- [Automatically installed dependencies](#automatically-installed-dependencies)
  - [Meson and Ninja](#meson-and-ninja)
  - [CoDiPack and MeDiPack](#codipack-and-medipack)
- [Configuration and Compilation](#configuration-and-compilation)
  - [Basic Configuration](#basic-configuration)
  - [Advanced Configuration](#advanced-configuration)
    - [Build Type](#build-type)
    - [Compiler optimizations](#compiler-optimizations)
    - [Warning level](#warning-level)
    - [Linear algebra options](#linear-algebra-options)
    - [Setting the compiler](#setting-the-compiler)
  - [Compilation](#compilation)
  - [Running Unit Tests](#running-unit-tests)
  - [Setting environment variables](#setting-environment-variables)
- [Troubleshooting](#troubleshooting)
  - [MPI installation is not found](#mpi-installation-is-not-found)
  - [mpi4py library is not found](#mpi4py-library-is-not-found)
  - [MKL is not found](#mkl-is-not-found)
  - [OpenMPI is used in place of MPICH](#openmpi-is-used-in-place-of-mpich)
  - [The wrong version of OpenMPI is used](#the-wrong-version-of-openmpi-is-used)

---

## Quick Compilation Guide ##

This is a quick guide to compile and install a *basic version* of SU2. For more information on the requirements and a more detailed description of the build system **continue reading** the rest of this page.

Short summary of the minimal requirements:

- C/C++ compiler
- Python 3

**Note:** all other necessary build tools and dependencies are shipped with the source code or are downloaded automatically.

If you have these tools installed, you can create a configuration using the `meson.py` found in the root source code folder:
```
./meson.py builddir
```
Use `ninja` to compile and install the code

```
./ninja -C builddir install
```
---


## Requirements ##

### Compilers ###
Installing SU2 from source requires a C++ compiler. The GNU compilers (gcc/g++) are open-source, widely used, and reliable for building SU2. The Intel compiler set has been optimized to run on Intel hardware and has also been used successfully by the development team to build the source code, though it is commercially licensed. The Apple LLVM compiler (Clang) is also commonly used by the developers.

- GNU gcc / g++ 
- Intel icc / icpc 
- Apple LLVM (clang)
  
**Note**: SU2 uses some C++11 features, that means at least GCC >= v4.7, Clang >= v3.0 or Intel C++ >= v12.0 is necessary.

### MPI ###
In order to build SU2 with parallel support, you need a suitable MPI installation on your machine. During the configuration the build tool does a check and enables MPI support. If no installation is found, a serial version of SU2 will be compiled.

### Python ###

SU2 requires Python 3 for compilation and for running the python scripts. Make sure that you have properly set up the `python3` executables in your environment. 

### Optional: swig and mpi4py ###
If you want to use the python wrapper capabilities, also `swig` and `mpi4py` are required. On **Linux** `swig` should be available in the package manager of your distribution and `mpi4py` can be installed using [pip](https://pip.pypa.io/en/stable/).

On **Mac OS X**, you can use the [Homebrew](http://brew.sh/) package manager. Once it is installed on your system, you can install Swig by running:

    $ sudo brew install swig

Install mpi4py with Python pip using easy install:

    $ sudo easy_install pip
    $ sudo pip install mpi4py
    
---

## Automatically installed dependencies ##

The following dependencies are automatically downloaded (or initialized if source code was cloned using `git`) during the [configuration](#configuration-and-compilation). 

### Meson and Ninja ###
The build system of SU2 is based on a combination of [meson](http://mesonbuild.com/) (as the front-end) and [ninja](https://ninja-build.org/) (as the back-end). Meson is an open source build system meant to be both extremely fast, and, even more importantly, as user friendly as possible. Ninja is a small low-level build system with a focus on speed. 

### CoDiPack and MeDiPack ###
In order to use the discrete adjoint solver the compilation requires two additional (header-only) libraries. [CoDi](https://github.com/SciCompKL/CoDiPack) provides the AD datatype and [MeDi](https://github.com/SciCompKL/MeDiPack) provides the infrastructure for the MPI communication when the reverse mode of AD is used. 

--- 
## Configuration and Compilation ##

Like mentioned above, SU2 uses meson and ninja for configuration and compilation, respectively. A configuration using meson is generated first and then an invocation of ninja is used to compile SU2 with this configuration. 

### Basic Configuration ###

The following steps assume that you have preinstalled version of meson and ninja available on your machine.  If you do not, you can use the versions distributed with SU2.  Just replace `meson` and `ninja` calls with `./meson.py` and `./ninja`.   

The only required argument for `meson` is a name of a directory where it should store the configuration. Meson is designed to perform out-of-source builds.  This means that the configuration-specific build commands and/or the build files should be placed in their own folder, independent from the source directory. You can have multiple configurations in different folders next to each other. To generate a basic configuration that will be stored in the folder `builddir` use

```
meson builddir
```

Options can be passed to the script to enable or disable different features of SU2.  Below you find a list of project options and their default values:
 
| Option | Default value | Description |
|---| --- | --- |
| `-Denable-autodiff`  | `false`   |   enable AD (reverse) support (needed for discrete adjoint solver)  |
| `-Denable-directdiff` | `false`     |  enable AD (forward) support |
| `-Denable-pywrapper` | `false`      |    enable Python wrapper support|
| `-Dwith-mpi`       | `auto` |   Set dependency mode for MPI (`auto`,`enabled`,`disabled`)  |
| `-Denable-cgns`     | `true`    |       enable CGNS support           |        
| `-Denable-tecio`    |  `true`       |    enable TECIO support         |
| `-Denable-mkl`      |  `false`      |    enable Intel MKL support     |
| `-Denable-boost-utf` | `false` | enable unit tests that require the Boost UTF |

Openblas and Pastix support have not been integrated into the hybrid branch. As an example, to enable AD support pass the option to the `meson` script along with a value:
```
meson builddir -Denable-autodiff=true
```
To set a installation directory for the binaries and python scripts, use the `--prefix` option, e.g.:

```
meson builddir -Denable-autodiff=true --prefix=/home/username/SU2
```
If you are not interested in setting custom compiler flags and other options you can now go directly to the [Compilation](#compilation) section, otherwise continue reading the next section.

### Advanced Configuration ###
In general meson appends flags set with the environment variable `CXX_FLAGS`. It is however recommended to use 
mesons built-in options to set debug mode, warning levels and optimizations. All options can be found [here](https://mesonbuild.com/Builtin-options.html) or by using `meson configure`. An already created configuration can be modified by using the `--reconfigure` flag, e.g.:
```
meson builddir --reconfigure --buildtype=debug
```
Note that it is only possible to change one option at once.

#### Build Type ####

The debug mode can be enabled by using the `--buildtype=debug` option. This adds `-g` flag and disables all compiler optimizations. If you still want to have optimizations, use `--buildtype=debugoptimized`. The default build type is `release`.

Many debugging checks are guarded by a `NDEBUG` preprocessor flag.  If `NDEBUG` is defined, these checks will not be compiled or executed.  By default, the `NDEBUG` flag is defined for release builds (when `buildtype=release`).  This is different than the upstream branch of SU2 and the defaults for meson. If you wish to change this behavior, use the compiler option `-Db_ndebug`. For example, `-Db_ndebug=false`.

#### Compiler optimizations ####

The optimization level can be set with `--optimization=level`, where `level` corresponds to a number between 0 (no optimization) and 3 (highest level of optimizations). The default level is 3.

#### Warning level ####

The warning level can be set with `--warnlevel=level`, where  `level` corresponds to a number between 0 (no warnings) and 3 (highest level of warning output). Level 1 corresponds to `-Wall`, level 2 to `-Wall -Wextra` and level 3 to `-Wall -Wextra -Wpedantic`. The default level is 0.

**Note:** The warning flags `-Wno-unused-parameter`, `-Wno-empty-body` and `-Wno-format-security` are always added by default.

#### Linear algebra options ####

Compiling with support for a BLAS library (`-Denable-mkl` or `-Denable-openblas`) is highly recommended if you use the high order finite element solver, or radial basis function interpolation in fluid structure interaction problems.
`-Denable-mkl` takes precedence over `-Denable-openblas`, by default the build system looks for MKL in `/opt/intel/mkl`, this can be changed via option `-Dmkl_root`.
When OpenBLAS support is requested the build system uses [pkg-config](https://en.wikipedia.org/wiki/Pkg-config) to search the system for package `openblas`, option `-Dblas-name`, if the library was built from source it may be necessary to set the environment variable PKG_CONFIG_PATH.

For large structural FEA problems on highly anisotropic grids iterative linear solvers might fail. Version 7 introduces experimental support for the direct sparse solver [PaStiX](https://gforge.inria.fr/projects/pastix/) (`-Denable-pastix`) see detailed instructions in `TestCases/pastix_support/readme.txt`.

**Note:** The BLAS library needs to provide support for LAPACK functions.

#### Setting the compiler ####

The preferred way to specify the compiler in meson is to use the `CC` and `CXX` environmental variables.  For example, if `'ccache mpicc'` is the C compiler and `'ccache mpicxx'` is the C++ compiler, then the meson command `meson builddir` would be replaced by:

```
export CC='ccache mpicc' && export CXX='ccache mpicxx' && meson builddir
```

### Compilation ###

Finally to compile and install SU2 use 
```
ninja -C builddir install
```
where `build` is again a folder with a configuration created using a call to `meson.py` described in the previous section. By default ninja uses all available cores in your system for the compilation. You can set the number of cores manually by using the `-jN` flag, where `N` is the number of cores you want to use.

Alternately, you can invoke the separate steps of the build using:
```
ninja -C builddir
ninja -C builddir test
ninja -C builddir install
```
This will compile, run unit tests, and then install SU2 in successive commands.


### Running Unit Tests ###

Meson and ninja together offer a simple unit testing framework, very similar to that offered by autotools.  To compile and run the unit tests, just build `test` as a target:

```
ninja -C builddir test
```

While autotools has the target `check` (e.g. `make check`), the target `check` is not a valid meson build target.


### Setting environment variables ###
Set the environment variables to use the executables from any directory without explicity specifying the path as described in the [installation section](/docs_v7/SU2-Linux-MacOS).

---

## Troubleshooting ##

### MPI installation is not found ###
Meson looks for an MPI installation using [pkg-config](https://en.wikipedia.org/wiki/Pkg-config). But if your MPI implementation does not provide them, it will search for the standard wrapper executables, `mpic`, `mpicxx`, `mpic++`. If these are not in your path, they can be specified by setting the standard environment variables `MPICC`, `MPICXX` during configuration.

### mpi4py library is not found ###
Meson imports the mpi4py module and searches for the include path. If it is installed in a custom location, make sure to add this path to the `PYTHONPATH` environment variable prior calling `meson.py`.

### MKL is not found ###

Meson looks for an MKL installation using [pkg-config](https://en.wikipedia.org/wiki/Pkg-config). You can check if `pkg-config` can find MKL by using the command:
```
pkg-config --libs mkl-static-lp64-seq
```
If pkg-config can't find MKL, it will explicitly say so.  If it can find MKL, it will print some linker flags and linker paths.  If pkg-config can't find MKL and it is indeed installed on your system, then you can tell pkg-config where to find MKL by adding its location to the `PKG_CONFIG_PATH` environmental variable.  For example, if the pkg-config script for MKL is located at `$MKLROOT/bin/pkgconfig`, you can set the environmental variable as:

```
export PKG_CONFIG_PATH=$MKLROOT/bin/pkgconfig
```

Note that the pkg-config path is one of the variables cached by meson.  That means you'll need to clear your cache every time you change the pkg-config path.  Otherwise, pkg-config will fail to recognize the updated `PKG_CONFIG_PATH`. You can clear the meson build cache by typing `meson configure --clearcache builddir` or deleting the entire build directory.

### OpenMPI is used in place of MPICH ###


This problem can be diagnosed in two ways:

1. You try to compile with MPICH, but when you run the executable it runs
   with several serial jobs instead of one parallel job.
2. When you compile the program with `ninja ... -v`, you see that the
   headers and linked files for MPI are *not* the MPICH headers and linked
   libraries.

This problem occurs when you have both an OpenMPI implementation and an
MPICH implementation on your system. If you try to compile and run with
an MPICH implementation, it will compile your program with OpenMPI no
matter what. As of Meson 0.52.1, Meson first uses pkg-config to look for
OpenMPI, then looks at the `mpicc` and `mpicxx` wrappers.  It *does* not
look for MPICH using pkg-config.  If it finds OpenMPI, it will halt
before it ever looks at `mpicc` or `mpicxx`, and then assume that you want
the OpenMPI implementation it found.

This is a problem with Meson, not with SU2.  Unfortunately, you cannot
trick meson by settings `CC=mpicc`.  SU2 uses a `WITH_MPI` preprocessor
flag that is only enabled if MPI is found.

The only known workarounds are to manually correct the meson files to look
for MPICH or to switch to OpenMPI.

### The wrong version of OpenMPI is used ###

If you've got multiple versions of OpenMPI, then you must make sure that
the path to your desired OpenMPI directory is at the start of the
pkg-config search path. For example, you could set:

```
export PKG_CONFIG_PATH=$MPI_DIR/lib/pkgconfig:$PKG_CONFIG_DIR
```

You can check that Meson will find the correct version of OpenMPI using
pkg-config by running the following command:
```
pkg-config --libs ompi-c
```
The library include folders shown should match the OpenMPI version you
are trying to use.
