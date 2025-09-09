<div align="center">

![](https://raw.githubusercontent.com/ViennaTools/ViennaLS/master/assets/logo.png)

<h1>ViennaHRLE</h1>

</div>

ViennaHRLE is a header-only C++ library for storing sparse spatial data efficiently. In the worst case, traversing the
whole data structure is achieved in O(N), where N is the number of data points stored in the structure. Random access is
achieved in O(log(N)).

## Support

For help with getting started, have a look at the [examples](examples/).

Bug reports and suggestions should be filed on [GitHub](https://github.com/ViennaTools/ViennaHRLE/issues/new).

## Releases

Releases are tagged on the maser branch and available in
the [releases section](https://github.com/ViennaTools/viennahrle/releases).

## Building

### System Requirements

* C++11 Compiler with OpenMP support

### Installing

Since this is a header only project, it does not require any installation. However, we recommend the following
procedure:

```bash
git clone https://github.com/ViennaTools/ViennaHRLE.git
cd ViennaHRLE

cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
cmake --build build
```

This will install the necessary headers and CMake files to the specified path. If `-DCMAKE_INSTALL_PREFIX` is not
specified, it will be installed to the standard path for your system, usually `/usr/local`.

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
  CPMAddPackage("gh:viennatools/viennahrle@0.7.0")
  ```

* With a local installation
  > In case you have ViennaHRLE installed in a custom directory, make sure to properly specify the `CMAKE_MODULE_PATH`

    ```cmake
    find_package(ViennaHRLE REQUIRED)
    target_link_libraries(${PROJECT_NAME} ViennaTools::ViennaHRLE)
    ```

## Authors

Contact us via: viennats@iue.tuwien.ac.at

ViennaHRLE was developed under the aegis of the 'Institute for Microelectronics' at the 'TU Wien'.
http://www.iue.tuwien.ac.at/

License
--------------------------
See file LICENSE in the base directory.
