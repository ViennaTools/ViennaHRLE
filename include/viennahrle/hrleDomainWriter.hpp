#ifndef HRLE_DOMAIN_WRITER_HPP
#define HRLE_DOMAIN_WRITER_HPP

#define HRLE_FILE_VERSION_NUMBER 1

#include <fstream>
#include <iostream>

/*
*********************************************************************
*********************** THE HRLE FILE FORMAT ***********************
*********************************************************************
*    File Header: 8 Bytes     *
********************************
* 4 Bytes   Identification Bytes (HRLE)
* 1 Byte    File Version Number
* 1 Byte    Endianess - Little Endian (0) or Big Endian (1)
* 1 Byte    Dimension (2 or 3)
* 1 Byte    Bytes Per Distance
*********************************************************************
*    Grid: 15+ Bytes    *
*************************
* 1 Byte    This byte contains the number of bytes
            used for the grid min and grid max of each dimension.
* NOTE: The following block is repeated for each dimension.
* x Bytes   Grid minimum
* x Bytes   Grid maximum
* 1 Byte    Boundary Condition
* NOTE: END
* 8 Bytes   GridDelta
*********************************************************************
*    H-RLE Block Header: 15 Bytes    *
**************************************
* NOTE: The following is repeated for each dimension
* 1 Byte    This byte contains the number of bytes
            used for each start index
* 1 Byte    This byte contains the number of bits(!!!!!!)
            used for each run type
* 1 Byte    This byte contains the number of bytes
            used for each runbreak
* 4 Bytes   Number of saved Start Indices
* 4 Bytes   Number of saved Runtypes
* 4 Bytes   Number of saved Runbreaks
*********************************************************************
*    H-RLE Block Data    *
**************************
* Start Indices   using adaptive number of bytes(delta encoded)
* Runtypes        using n bits per runtype
* Indices of runtypes  using adaptive number of bytes(delta encoded)
* Runbreaks       using adaptive number of bytes
* NOTE: END
*********************************************************************
*    Values Header: 4 Bytes    *
***********************************
* 4 Bytes  Number of defined values
* 4 Bytes  Number of undefined values
*********************************************************************
*    Values Data    *
************************
* Defined Values - using sizeOf Bytes per value or
                   serialize/deserialize of underlying data structure
* Undefined Values - same as Defined Values
*********************************************************************
*/

/// Class which handles the output of an hrleDomain
/// to the binary .hrle format
template <class hrleDomain> class hrleDomainWriter {
  typedef typename hrleDomain::hrleValueType hrleValueType;

  static constexpr int D = hrleDomain::dimension;

  union {
    uint16_t shortVar;  // binary  number of length 16 Bits
    uint8_t charVar[2]; // 2 binary numbers, each 8 Bits
  } test_endianness;

  hrleDomain *domain;
  std::string filePath;
  int valueTypeByteSize =
      sizeof(hrleValueType); // if hrleValueType is float or double, reduce to
                             // this bytesize before saving

  bool bigEndian() {
    test_endianness.shortVar = 0x8000; // MSB of 16
    return test_endianness.charVar[0] != 0;
  }

public:
  hrleDomainWriter() = default;

  explicit hrleDomainWriter(hrleDomain *domainPointer)
      : domain(domainPointer) {}

  explicit hrleDomainWriter(hrleDomain &passedDomain) : domain(&passedDomain) {}

  // setters and getters
  void setDomain(hrleDomain *domainPointer) { domain = domainPointer; }
  void setDomain(hrleDomain &passedDomain) { domain = &passedDomain; }
  hrleDomain &getDomain() { return *domain; }
  void setFilePath(std::string path) {
    if (path.find(".hrle") != path.size() - 5) {
      std::cout << "File path not ending in '.hrle'!" << std::endl;
      path.append(".hrle");
      std::cout << "Using '" << path << "' instead." << std::endl;
    }
    filePath = path;
  }
  const std::string &getFilePath() const { return filePath; }

  void apply() {
    std::ofstream fout(filePath, std::ofstream::binary);
    if (!fout.is_open()) {
      std::cout << "ERROR: Could not open the file: " << filePath << std::endl;
      return;
    }

    // write file header
    {
      fout << "HRLE";
      fout << HRLE_FILE_VERSION_NUMBER;
      fout << (bigEndian() ? 1 : 0);
      fout << D;
      fout << char(valueTypeByteSize);
    }

    // now write grid serialization
    domain->getGrid().serialize(fout);

    // now write hrle structure
    domain->serialize(fout);
  }
};

#endif // HRLE_DOMAIN_WRITER_HPP
