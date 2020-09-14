# Acoustic

Wide-anlge mode parabolic equation (WAMPE) solver for Helmholtz equation.

# Installing dependencies

Both [fftw3](http://www.fftw.org/) and [Boost](https://www.boost.org/) (namely [Boost.Program_options](https://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html)) are required.

## Linux

Just use your favorite (or not) package manager. Nice isn't it?

## Windows

[vcpkg](https://github.com/microsoft/vcpkg) is highly recommended. You can use something like
```powershell
PS> ./vcpkg install --triplet x64-windows boost-program-options fftw3 nlohmann-json
```

# Building 

## Linux

```bash
$ cd build/cmake
$ mkdir ACOUSTIC && cd ACOUSTIC
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
```

## Windows

```powershell
PS> cd build/cmake
PS> mkdir ACOUSTIC
PS> cd ACOUSTIC
PS> cmake -DCMAKE_TOOLCHAIN_FILE=[vcpkg root]\scripts\buildsystems\vcpkg.cmake -DUSE_VCPKG=true -DCMAKE_BUILD_TYPE=Release ..
```
