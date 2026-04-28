#ifndef HRLE_DOMAIN_READER_HPP
#define HRLE_DOMAIN_READER_HPP

#define HRLE_FILE_READ_VERSION_NUMBER 1

#include <fstream>
#include <iostream>

#include <hrleGrid.hpp>
#include <vcLogger.hpp>

namespace viennahrle {
using namespace viennacore;
/// Class which handles the input of an hrleDomain
/// from a binary .hrle file
template <class hrleDomain> class DomainReader {
  typedef typename hrleDomain::ValueType ValueType;
  static constexpr int D = hrleDomain::dimension;

  hrleDomain *domain = nullptr;
  std::string filePath;
  int valueTypeByteSize = 0;

  bool bigEndian() {
    uint16_t number = 0x1;
    char *numPtr = reinterpret_cast<char *>(&number);
    return (numPtr[0] != 1);
  }

public:
  DomainReader() = default;
  explicit DomainReader(hrleDomain *domainPointer) : domain(domainPointer) {}
  explicit DomainReader(hrleDomain &passedDomain) : domain(&passedDomain) {}
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
      VIENNACORE_LOG_ERROR(
          "In order to read an hrleDomain, you first have to set the "
          "object to read it into. Use DomainReader.setDomain()");
      return;
    }

    if (filePath.find(".hrle") == std::string::npos) {
      VIENNACORE_LOG_ERROR(
          "File name does not have the correct file ending: '.hrle'");
      return;
    }

    std::ifstream fin(filePath, std::ifstream::binary);
    if (!fin.is_open()) {
      VIENNACORE_LOG_ERROR("Could not open the file: " + filePath);
      return;
    }

    // FILE HEADER
    char buff[9] = {}; // 1 extra byte for string constructor
    fin.read(buff, 8);
    // Comparing Identification Bytes
    if (std::string(buff).compare(0, 4, "HRLE")) {
      VIENNACORE_LOG_ERROR("File is not an HRLE file.");
      return;
    }
    if (HRLE_FILE_READ_VERSION_NUMBER != buff[4] - 48) {
      VIENNACORE_LOG_WARNING("File of version " + std::to_string(buff[4] - 48) +
                             " is read by this reader(Version " +
                             std::to_string(HRLE_FILE_READ_VERSION_NUMBER) +
                             ")!");
      if (HRLE_FILE_READ_VERSION_NUMBER < buff[4] - 48)
        return;
    }
    if (bigEndian() != bool(buff[5] - 48)) {
      VIENNACORE_LOG_WARNING("File was written in a different byte order than "
                             "it is being read. Results may be incorrect!");
    }
    const int dimension = buff[6] - 48;
    if (dimension != D) {
      VIENNACORE_LOG_ERROR("ERROR: Domain in file has " +
                           std::to_string(dimension) +
                           " dimensions, but trying to read domain with " +
                           std::to_string(D) + " dimensions.");
      return;
    }

    valueTypeByteSize = int(buff[7]);

    // read grid
    domain->getGrid().deserialize(fin);
    // read hrleDomain
    domain->deserialize(fin);
  }
};
} // namespace viennahrle

#endif // HRLE_DOMAIN_READER_HPP
