<div align="center">

![](https://raw.githubusercontent.com/ViennaTools/ViennaLS/master/assets/logo.png)

<h1>ViennaHRLE</h1>

</div>

ViennaHRLE is a header-only C++ library for storing sparse spatial data efficiently. In the worst case, traversing the whole data structure is achieved in O(N), where N is the number of data points stored in the structure. Random access is achieved in O(log(N)).

## Support
For help with getting started, have a look at the [examples](Examples/).

Bug reports and suggestions should be filed on [GitHub](https://github.com/ViennaTools/ViennaHRLE/issues/new).

## Contributing
If you want to contribute to ViennaHRLE, make sure to follow the [LLVM Coding guidelines](https://llvm.org/docs/CodingStandards.html). Before creating a pull request, make sure ALL files have been formatted by clang-format, which can be done using the format-project.sh script in the root directory.

## Releases
Releases are tagged on the maser branch and available in the [releases section](https://github.com/ViennaTools/viennahrle/releases).

## Building

### System Requirements

* C++11 Compiler with OpenMP support

### Installing

Since this is a header only project, it does not require any installation. However, we recommend the following procedure:

```bash
git clone https://github.com/ViennaTools/ViennaHRLE.git
cd ViennaHRLE

cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
cmake --build build
```

This will install the necessary headers and CMake files to the specified path. If `-DCMAKE_INSTALL_PREFIX` is not specified, it will be installed to the standard path for your system, usually `/usr/local`.

## Building examples

The examples can be built using CMake:

```bash
cmake -B build -DVIENNAHRLE_BUILD_EXAMPLES=ON
cmake --build build
```

## Integration in CMake projects

We recommend using [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake) to consume this library.

* Installation with CPM
  ```cmake
  CPMAddPackage("gh:viennatools/viennahrle@0.4.0")
  ```

* With a local installation
    > In case you have ViennaHRLE installed in a custom directory, make sure to properly specify the `CMAKE_MODULE_PATH`

    
    ```cmake
    find_package(ViennaHRLE REQUIRED)
    target_link_libraries(${PROJECT_NAME} ViennaTools::ViennaHRLE)
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
