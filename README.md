# What is libamtrack?

* libamtrack provides computational routines for the prediction of detector response and radiobiological efficiency in heavy charged particle beams.
* libamtrack is designed for research in proton and ion dosimetry and radiotherapy.
* libamtrack also provides many auxiliary routines for the work with proton and ion beams.
* libamtrack is a program library.
* libamtrack works both under Linux, Windows, and Mac OS.
* libamtrack is written in ANSI C.
* libamtrack is free, open-source software under GNU GPL licence.
* libamtrack is intended to facilitate the comparison of and the communication on amorphous track models for particle beam research.


# 2. How can I use libamtrack?

libamtrack is a program library and cannot be run as a single executable. A number of interfaces is provided with different complexities depending on your needs and experience. They are given below in order of complexity. Please be aware that for option (iv)-(vi) the GNU Scientific Library (GSL) has to be installed on your system.
Refer to README_DEVELOPERS on how to do that.

## i. Web-interface
The web interface covers a subset of libamtrack functions for easy access and can be found on libamtrack's website: https://libamtrack.github.io/web/. The sources can be found at (/distributions/JavaScript). The interface development was started by Christoph Kolb within his Bachelor thesis and consolidated by Leszek Grzanka.

## ii. R

To access most functions of **libamtrack**, we originally recommended using the R environment. **libamtrack** was a contributed package on the [Comprehensive R Archive Network (CRAN)](https://CRAN.R-project.org/package=libamtrack) and could be easily installed and used within R, along with documentation.

**However, please note that this package is no longer maintained.** It was removed from CRAN on *September 3, 2019*. Formerly available versions can still be accessed from the [CRAN archive](https://cran.r-project.org/src/contrib/Archive/libamtrack/), but no further updates or support are provided.

## iii. Python, Matlab wrappers
libamtrack comes with a collection of wrappers for Python and Matlab (/distributions/Python, /distributions/Matlab).

## iv. Use precompiled binaries
If you want to use **libamtrack** in your own code, precompiled binaries (including headers) were previously available for Windows and macOS on the **libamtrack** website. However, **these are no longer available at the moment**.

For Linux (or other less common operating systems), the library must be compiled and installed manually using **CMake**. While rpm/deb packages were planned, **they are currently not available**.


## vi. Work with full sources / development

If you want to compile it on your own latest source code of libamtrack, please use this to build and install according to your OS.

Here is the updated documentation with **proper header levels**, maintaining clarity and structure.

---

### **Linux Installation**

#### **Requirements**
Ensure the following dependencies are installed before building **libamtrack**:

- `git`
- `gcc`
- `libtool`
- `gfortran`
- `libgsl-dev`
- `cmake`

#### **Getting the Source Code**
Clone the latest version of **libamtrack** from GitHub:

```bash
git clone https://github.com/libamtrack/library.git
```

#### **Building the Library**
Navigate to the source directory and create a `build/` directory:

```bash
cd library && mkdir build && cd build
```

Then, configure the build and compile the library:

```bash
cmake ..
cmake --build . --parallel
```

#### **Installing the Library**
The **preferred installation method** is to use a **custom directory** rather than installing system-wide. This prevents conflicts with system libraries and allows easy removal.

##### **Recommended Installation Locations**
Instead of installing to `/usr/local/`, consider these options:
- **Inside your home directory (recommended for personal use)**
  ```bash
  $HOME/.local
  ```
- **A dedicated user directory for local installations**
  ```bash
  $HOME/usr
  ```

To install **libamtrack** into a custom directory (e.g., `$HOME/.local`), run:

```bash
cmake --install . --prefix $HOME/.local
```

This will install:
- **Libraries** → `$HOME/.local/lib`
- **Executables** → `$HOME/.local/bin`
- **Headers** → `$HOME/.local/include` (if `BUILD_DEV=ON`)
- **CMake package files** → `$HOME/.local/lib/cmake/libamtrack` (if `BUILD_DEV=ON`)

If you prefer installing **inside the build directory** for an isolated setup:

```bash
cmake --install . --prefix ./install
```

This will place all installed files inside `build/install/`.

#### **Installing Development Files (Headers & CMake Config)**
If you need headers and CMake package files for development, enable `BUILD_DEV`:

```bash
cmake .. -DBUILD_DEV=ON
cmake --build . --parallel
cmake --install . --prefix $HOME/.local
```

---

#### **Uninstalling libamtrack**
To **uninstall** the library, you need to remove the installed files manually, as CMake does not provide a built-in uninstall command.

##### **If installed in `$HOME/.local` or another custom directory:**
Simply delete the relevant files:
```bash
rm -rf $HOME/.local/lib/libamtrack*
rm -rf $HOME/.local/include/libamtrack
rm -rf $HOME/.local/bin/amtrack*
rm -rf $HOME/.local/lib/cmake/libamtrack
```

If installed in `$HOME/usr`, adjust the paths accordingly:
```bash
rm -rf $HOME/usr/lib/libamtrack*
rm -rf $HOME/usr/include/libamtrack
rm -rf $HOME/usr/bin/amtrack*
rm -rf $HOME/usr/lib/cmake/libamtrack
```

##### **If installed in the system-wide `/usr/local` (not recommended)**
You will need `sudo` to remove files:
```bash
sudo rm -rf /usr/local/lib/libamtrack*
sudo rm -rf /usr/local/include/libamtrack
sudo rm -rf /usr/local/bin/amtrack*
sudo rm -rf /usr/local/lib/cmake/libamtrack
```

##### **Alternative: Using CMake’s Install Manifest**
If you originally installed using CMake, it may have generated an **install_manifest.txt** in your `build` directory. You can use it to remove installed files:
```bash
xargs rm < install_manifest.txt
```

---

#### **Final Notes**
- **Avoid system-wide installation** (`/usr/local/`) unless necessary.
- Using `$HOME/.local` or `$HOME/usr` is **preferred** for user-level installations.
- If using a **custom directory**, you may need to add it to your environment:
  ```bash
  export PATH="$HOME/.local/bin:$PATH"
  export LD_LIBRARY_PATH="$HOME/.local/lib:$LD_LIBRARY_PATH"
  export CMAKE_PREFIX_PATH="$HOME/.local:$CMAKE_PREFIX_PATH"
  ```
  Add these lines to your `.bashrc` or `.zshrc` for persistent configuration.

Would you like an **automated uninstall script** for easier removal? 🚀

### Windows (MSVC)

Start 64 bit "x64 Native Tools Command Prompt for VS 2022" and then run the following commands:

```cmd
vcpkg install --triplet x64-windows
```

That will install the required dependencies for libamtrack.

Then create the build directory and run the following commands:

```cmd
mkdir build
```

Then go to the `build` directory and run the following commands:

```cmd
cd build
```

Then run the following command:

```cmd
cmake .. -DGSL_INCLUDE_DIR="../vcpkg_installed/x64-windows/include" -DGSL_LIBRARY="../vcpkg_installed/x64-windows/lib/gsl.lib" -DGSL_CBLAS_LIBRARY="../vcpkg_installed/x64-windows/lib/gslcblas.lib" -DGETOPT_LIBRARY="../vcpkg_installed/x64-windows/lib/getopt.lib"  -DGETOPT_INCLUDE_DIR="../vcpkg_installed/x64-windows/include" -G "Ninja"
```

Then run the following command:

```cmd
cmake --build . --parallel
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
cd library && mkdir build && cd build
cmake -S .. -B .
cmake --build . --parallel
cmake --install .
```

### MacOS

First, install [Homebrew](https://brew.sh/) to manage the required packages.

Requirements:
- git
- gcc
- libtool
- gsl
- cmake

```bash
brew install git gcc libtool gsl cmake
```

Get the latest source code of the libamtrack:

```bash
git clone https://github.com/libamtrack/library.git
```

Go to `library` directory and then create `build` directory:

```bash
cd library && mkdir build && cd build
cmake -S .. -B .
cmake --build . --parallel
sudo cmake --install .
```

# 3. Can I use libamtrack in my research and/or modify the code?

Everybody is welcome to read, use and modify (preferably to improve) the code according to GNU GPL 3.


# 4. Where do I find documentation on libamtrack?

- The **libamtrack** manual is available here: `/docs/libamtrackManual.pdf`.

⚠ **Note:** This documentation may not be up to date. A new and improved version is currently being developed and will be released soon. Please check the repository or project website for updates.