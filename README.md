# ViennaHRLE

ViennaHRLE is a header-only C++ library for storing sparse spatial data efficiently. In the worst case, traversing the whole data structure is achieved in O(N), where N is the number of data points stored in the structure. Random access is achieved in O(log(N)).

## Support
For help with getting started, have a look at the [Examples](https://github.com/ViennaTools/viennahrle/tree/master/Examples).

Bug reports and suggestions should be filed on GitHub.

## Contributing
If you want to contribute to ViennaHRLE, make sure to follow the [LLVM Coding guidelines](https://llvm.org/docs/CodingStandards.html). Before creating a pull request, make sure ALL files have been formatted by clang-format, which can be done using the format-project.sh script in the root directory.

## Releases
Releases are tagged on the maser branch and available in the [releases section](https://github.com/ViennaTools/viennahrle/releases).

## Building

### System Requirements

* C++ Compiler with OpenMP support

### Installing

Since this is a header only project, it does not require any installation. However, we recommend the following procedure:

```
git clone github.com/ViennaTools/viennahrle.git
cd viennahrle
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
make install
```

This will install the necessary header and cmake files to the specified path. If DCMAKE_INSTALL_PREFIX is not specified, it will be installed to the standard path for your system, usually /usr/local/ .

### Building examples

The examples can be built using CMake (recommended) or GNU make.\
Building for Windows is only supported using CMake with Visual Studio.

Build using CMake (recommended):

```
mkdir build && cd build
cmake .. -DVIENNAHRLE_BUILD_EXAMPLES=ON
make
```

Build using GNU make:

```
cd Examples
make
```

### Integration in CMake projects

In order to use this library in your CMake project, add the following lines to the CMakeLists.txt of your project

```
set(ViennaHRLE_DIR "/path/to/your/custom/install/")
find_package(ViennaHRLE REQUIRED)
target_link_libraries(${PROJECT_NAME} ViennaHRLE)
```

ViennaHRLE also sets the variable VIENNAHRLE_INCLUDE_DIRS, which might be more convenient for certain projects. Therefore, it is also possible to integrate ViennaHRLE as

```
set(ViennaHRLE_DIR "/path/to/your/custom/install/")
find_package(ViennaHRLE REQUIRED)
target_include_directories(${PROJECT_NAME} ${VIENNAHRLE_INCLUDE_DIRS})
```

## Authors

Current contributors: Lado Filipovic, Paul Manstetten, Xaver Klemenschits and Josef Weinbub

Founder and initial developer: Otmar Ertl

Contact us via: viennats@iue.tuwien.ac.at

ViennaHRLE was developed under the aegis of the 'Institute for Microelectronics' at the 'TU Wien'.
http://www.iue.tuwien.ac.at/

License
--------------------------
See file LICENSE in the base directory.
