# Converter for polynomials with multi-precision interger coefficients
Convert polynomials with multi-precision interger coefficients to RNS form.

## Build MPpoly-converter

### Dependencies
We have tested MPpoly-converter on Ubuntu.

- CMake >= 3.14
- gcc >= 9.0
- [GMP](https://gmplib.org/) library

Some example codes require Intel [HEXL](https://github.com/intel/hexl).

We assume that external libraries are installed in `$HOME/.local`. If you install them in other directories, please specify the directories in CMakeLists.txt.

### Options
#### ENABLE_TEST
Set ON to enable building of unit test. Default is ON.

#### ENABLE_EXAMPLE
Set ON to enable building of example codes. Default is OFF.

## Testing
The unit-test executable itself is located at `build/test/unit-test`, or just do `ctest` command in `build/`.

## Examples
You have to set ENABLE_EXAMPLE ON to build examples.

### speed_test
Measure times of convertion and polynomial multiplication. The executable is located at `build/example/speed_test`.