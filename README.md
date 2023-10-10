# 1. WHAT IS LIBAMTRACK?

* libamtrack provides computational routines for the prediction of detector response and radiobiological efficiency in heavy charged particle beams.
* libamtrack is designed for research in proton and ion dosimetry and radiotherapy.
* libamtrack also provides many auxiliary routines for the work with proton and ion beams.
* libamtrack is a program library.
* libamtrack works both under Linux, Windows, and Mac OS.
* libamtrack is written in ANSI C.
* libamtrack is free, open-source software under GNU GPL licence.
* libamtrack is intended to facilitate the comparison of and the communication on amorphous track models for particle beam research.


# 2. HOW CAN I USE LIBAMTRACK?

libamtrack is a program library and cannot be run as a single executable. A number of interfaces is provided with different complexities depending on your needs and experience. They are given below in order of complexity. Please be aware that for option (iv)-(vi) the GNU Scientific Library (GSL) has to be installed on your system.
Refer to README_DEVELOPERS on how to do that.

## i. Web-interface
The web interface covers a subset of libamtrack functions for easy access and can be found on libamtrack's website: https://libamtrack.github.io/web/. The sources can be found at (/distributions/JavaScript). The interface development was started by Christoph Kolb within his Bachelor thesis and consolidated by Leszek Grzanka.

## ii. R
To access (almost any) function of libamtrack, we recommend to use the R environment. libamtrack is a contributed package (http://cran.r-project.org/web/packages/libamtrack) on the Comprehensive R Archive Network (CRAN, http://cran.r-project.org) and can easily be installed and used from within R. It comes with documentation.

## iii. Python, Matlab wrappers
libamtrack comes with a collection of wrappers for Python and Matlab (/distributions/Python, /distributions/Matlab).

## iv. Use precompiled binaries
If you want to use libamtrack in your own code, please try to use the precompiled binaries (incl. headers) for your OS (Win, Mac OS X) are found on libamtrack's website. In case you are running Linux (or some exotic OS) you will have to compile and install libamtrack using autotools. But in the future, rpm/deb packages will be provided. 

## vi. Work with full sources / development

If you want to compile it on your own latest source code of libamtrack, please use this to build and install according to your OS.

### Linux

Requirements:
- git
- gcc
- libtool
- gfortran
- libgsl-dev
- cmake

Get the latest source code of the libamtrack:

```bash
git clone https://github.com/libamtrack/library.git
```

Go to `library` directory and then create `build` directory:

```bash
cd library && mkdir build
```

Generate build files for the libamtrack:
```bash
cmake -S . -B build
```

Build the libamtrack:
```bash
cmake -S . -B build
```

Install the libamtrack:
```bash
cmake --install build
```

### Windows(MSYS2)

Get first the [MSYS2](https://www.msys2.org/) and please follow the installation guide.

Requirements:
- git
- gcc
- libtool
- mingw-w64-x86_64-gcc-libgfortran
- mingw-w64-x86_64-gsl
- cmake

Get the latest source code of the libamtrack:

```bash
git clone https://github.com/libamtrack/library.git
```

Go to `library` directory and then create `build` directory:

```bash
cd library && mkdir build
```

Generate build files for the libamtrack:
```bash
cmake -S . -B build
```

Build the libamtrack:
```bash
cmake -S . -B build
```

Install the libamtrack:
```bash
cmake --install build
```
### MacOS

Get first the [Homebrew](https://brew.sh/) to install the required packages.

Requirements:
- git
- gcc
- libtool
- gsl
- cmake

Get the latest source code of the libamtrack:

```bash
git clone https://github.com/libamtrack/library.git
```

Go to `library` directory and then create `build` directory:

```bash
cd library && mkdir build
```

Generate build files for the libamtrack:
```bash
cmake -S . -B build
```

Build the libamtrack:
```bash
cmake -S . -B build
```

Install the libamtrack:
```bash
cmake --install build

# 3. CAN I USE LIBAMTRACK IN MY RESEARCH AND/OR MODIFY THE CODE?

Everybody is welcome to read, use and modify (preferably to improve) the code according to GNU GPL 3.


# 4. WHERE DO I FIND DOCUMENTATION ON LIBAMTRACK?

- The libamtrack manual is found here: /docs/libamtrackManual.pdf