# simpledsp

Simple DSP functions.

## Overview

Functions for DSP, mostly targeting real time application.

## Implemented
- FFT (radix2 / radix4)

## TODO
- IIR filter
- FIR filter
- Others

## Required to build tests
- Conan2
- CMake

## Getting started
Commands for Ubuntu 22.04.  Other platforms may vary.
```
git clone https://github.com/mike919192/simpledsp.git
cd simpledsp
conan install . --build=missing
cd build/Release/generators/
cmake ../../.. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release
cmake --build .
./testFFT
```
