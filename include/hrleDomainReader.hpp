#ifndef HRLE_DOMAIN_READER_HPP
#define HRLE_DOMAIN_READER_HPP

#define HRLE_FILE_READ_VERSION_NUMBER 1

#include <cmath>
#include <fstream>
#include <iostream>

#include <hrleGrid.hpp>
#include <hrleIndexType.hpp>
#include <hrleRunTypeValues.hpp>

/// Class which handles the input of an hrleDomain
/// from a binary .hrle file
template <class hrleDomain> class hrleDomainReader {
  typedef typename hrleDomain::hrleValueType hrleValueType;
  static constexpr int D = hrleDomain::dimension;

  hrleDomain *domain = nullptr;
  std::string filePath;
  int valueTypeByteSize;

  union {
    uint16_t shortVar;  // binary  number of length 16 Bits
    uint8_t charVar[2]; // 2 binary numbers, each 8 Bits
  } test_endianness;

  bool bigEndian() {
    test_endianness.shortVar = 0x8000; // MSB of 16
    return test_endianness.charVar[0] != 0;
  }

public:
  hrleDomainReader() {}
  hrleDomainReader(hrleDomain *domainPointer) : domain(domainPointer) {}
  hrleDomainReader(hrleDomain &passedDomain) : domain(&passedDomain) {}
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
    if (domain == nullptr) {
      std::cout
          << "ERROR: In order to read an hrleDomain, you first have to set the "
             "object to read it into. Use hrleDomainReader.setDomain()"
          << std::endl;
      return;
    }

    if (filePath.find(".hrle") == std::string::npos) {
      std::cout
          << "ERROR: File name does not have the correct file ending: '.hrle'"
          << std::endl;
      return;
    }

    std::ifstream fin(filePath, std::ifstream::binary);
    if (!fin.is_open()) {
      std::cout << "ERROR: Could not open the file: " << filePath << std::endl;
      return;
    }

    // FILE HEADER
    char buff[8] = {};
    fin.read(buff, 8);
    // Comparing Identification Bytes
    if (std::string(buff).compare(0, 4, "HRLE")) {
      std::cout << "ERROR: File is not an HRLE file." << std::endl;
      return;
    }
    if (HRLE_FILE_READ_VERSION_NUMBER != buff[4] - 48) {
      std::cout << "WARNING: File of version " << buff[4] - 48
                << " is read by this reader(Version "
                << HRLE_FILE_READ_VERSION_NUMBER << ")!" << std::endl;
      if (HRLE_FILE_READ_VERSION_NUMBER < buff[4] - 48)
        return;
    }
    if (bigEndian() != bool(buff[5] - 48)) {
      std::cout << "WARNING: File was written in a different byte order than "
                   "it is being read. Results may be incorrect!"
                << std::endl;
    }
    const int dimension = buff[6] - 48;
    if (dimension != D) {
      std::cout << "ERROR: Domain in file has " << dimension
                << " dimensions, but trying to read domain with " << D
                << " dimensions." << std::endl;
      return;
    }

    valueTypeByteSize = int(buff[7]);

    // read grid
    domain->getGrid().deserialize(fin);

    // create new hrleDomain with new grid
    hrleDomain newDomain(domain->getGrid());

    // deserialize into new domain
    newDomain.deserialize(fin);

    // now copy into old domain
    domain->deepCopy(newDomain);
  }
};

#endif // HRLE_DOMAIN_READER_HPP
