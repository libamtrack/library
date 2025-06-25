# libamtrack Distributions

This directory contains language-specific and platform-specific distributions of the libamtrack library. The libamtrack library provides computational routines for the prediction of detector response and radiobiological efficiency in heavy charged particle beams.

## Directory Structure

### Java

- **Status:** Unmaintained
- **Description:** Contains Java GUI for libamtrack with compilation scripts for both Windows and Linux
- **Features:** 
  - Web-based interface with Java support
  - Compilation scripts (make_jar.bat for Windows, make_jar.sh for Linux)
  - SWIG interface for C/Java integration
- **Requirements:** Java JDK, SWIG, GSL, GCC compiler
- **Note:** Not actively maintained, contributors welcome

### JavaScript

- **Status:** active development moved to separate repository
- **Description:** Contains only an README.md file
- **Note:** Active development moved to https://github.com/libamtrack/web/

### Linux

- **Status:** Unknown, likely not fully working
- **Description:** Contains scripts and files for Debian packaging of libamtrack
- **Features:** Creates .deb files for installation on Ubuntu systems
- **Note:** Contributions and improvements welcome

### Matlab

- **Status:** Unknown, potentially incomplete
- **Description:** Contains Matlab wrapper for library's distribution functionalities
- **Features:** 
  - Wrapper generation scripts
  - NAMESPACE definitions
- **Note:** May require updates, contributions welcome

### Python (pyamtrack)

- **Status:** Legacy version 0.14.0 (2022-10-03), active development moved to separate repository
- **Description:** Python wrapper for libamtrack library
- **Features:** 
  - Binary wheel package for Linux only
  - Scripts for wheel package generation and testing
- **Limitations:** 
  - No Windows or macOS support
  - Limited documentation
- **Note:** Active development moved to https://github.com/libamtrack/pyamtrack/

### R

- **Status:** Work in progress, no stable version available
- **Description:** R package for libamtrack aimed at physicists and researchers
- **History:** Previously existed in CRAN repository, removed in 2019
- **Note:** Contributions and feedback welcome

## Overall Status

The libamtrack distributions are in varying states of maintenance. The Python distribution has been actively moved to a separate repository, while others (Java, Matlab, R) are either unmaintained or work in progress. Contributions to any of these distributions are welcome.

For the latest information on libamtrack, refer to the main project repository and documentation.
