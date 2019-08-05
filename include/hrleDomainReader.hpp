#ifndef HRLE_DOMAIN_READER_HPP
#define HRLE_DOMAIN_READER_HPP

#define HRLE_FILE_READ_VERSION_NUMBER 0

#include <cmath>
#include <fstream>
#include <iostream>

#include <hrleGrid.hpp>
#include <hrleIndexType.hpp>
#include <hrleRunTypeValues.hpp>

#include <bitset> // TODO remove

/// Class which handles the input of an hrleDomain
/// from a binary .hrle file
template <class hrleDomain> class hrleDomainReader {
  typedef typename hrleDomain::hrleValueType hrleValueType;
  static constexpr int D = hrleDomain::dimension;

  hrleDomain *domain = NULL;
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

  bool read() {
    if (domain == NULL) {
      std::cout
          << "ERROR: In order to read an hrleDomain, you first have to set the "
             "object to read it into. Try using hrleDomainReader.setDomain()"
          << std::endl;
      return false;
    }

    if (filePath.find(".hrle") == std::string::npos) {
      std::cout
          << "ERROR: File name does not have the correct file ending: '.hrle'"
          << std::endl;
      return false;
    }

    std::ifstream fin(filePath);
    if (!fin.is_open()) {
      std::cout << "ERROR: Could not open the file: " << filePath << std::endl;
      return false;
    }

    // FILE HEADER
    char buff[9] = {};
    unsigned char byte;
    fin.read(buff, 9);
    // Comparing Identification Bytes
    if (std::string(buff).compare(0, 4, "HRLE")) {
      std::cout << "ERROR: File is not an HRLE file." << std::endl;
    }
    if (HRLE_FILE_READ_VERSION_NUMBER != buff[4] - 48)
      std::cout << "ERROR: File version is deprecated and cannot be read!"
                << std::endl;
    if (bigEndian() != buff[5] - 48)
      std::cout << "WARNING: File was written in a different byte order than "
                   "it is being read. Results may be incorrect!"
                << std::endl;

    const int dimension = buff[6] - 48;
    if (dimension != D) {
      std::cout << "ERROR: Domain in file has " << dimension
                << " dimensions, but trying to read domain with " << D
                << " dimensions." << std::endl;
    }
    valueTypeByteSize = int(buff[7]);
    const int gridBoundaryBytes = int(buff[8]);

    // READ GRID
    hrleIndexType gridMin[D], gridMax[D];
    typename hrleGrid<D>::boundaryType boundaryConditions[D];
    double gridDelta;
    for (int i = D - 1; i >= 0; --i) {
      unsigned condition = 0;
      char min, max;
      fin.read(&min, gridBoundaryBytes);
      fin.read(&max, gridBoundaryBytes);
      fin.read((char *)&condition, 1);
      boundaryConditions[i] = typename hrleGrid<D>::boundaryType(condition);
      gridMin[i] = hrleIndexType(min);
      gridMax[i] = hrleIndexType(max);
    }
    fin.read((char *)&gridDelta, sizeof(double));

    // create hrleDomain with grid properties
    hrleGrid<D> grid(gridMin, gridMax, gridDelta, boundaryConditions);
    hrleDomain newDomain(grid);

    // READ HRLE PROPERTIES
    for (int dim = D - 1; dim >= 0; --dim) {
      // get the start indices, runtypes and runbreaks vectors
      std::vector<hrleSizeType> &startIndices =
          newDomain.domainSegments[0].startIndices[dim];
      std::vector<hrleSizeType> &runTypes =
          newDomain.domainSegments[0].runTypes[dim];
      std::vector<hrleIndexType> &runBreaks =
          newDomain.domainSegments[0].runBreaks[dim];

      int32_t numberOfStartIndices, numberOfRunTypes, numberOfRunBreaks;
      char bytesPerIndex, bitsPerRunType, bytesPerRunBreak;
      // reading in the 15 byte H-RLE header
      fin.read(&bytesPerIndex, 1);
      fin.read(&bitsPerRunType, 1);
      fin.read(&bytesPerRunBreak, 1);
      fin.read((char *)&numberOfStartIndices, 4);
      fin.read((char *)&numberOfRunTypes, 4);
      fin.read((char *)&numberOfRunBreaks, 4);

      // READ START INDICES
      {
        startIndices.clear();
        unsigned long long sum = 0;
        // push the 0, it was not written to the file
        startIndices.push_back(0);
        for (int i = 0; i < numberOfStartIndices - 1; ++i) {
          unsigned long long current = 0;
          fin.read((char *)&current, bytesPerIndex);
          sum += current;
          startIndices.push_back(sum);
        }
      }

      // READ RUN TYPES
      {
        runTypes.clear();
        int numberOfValuesPerByte = 8 / bitsPerRunType;
        int numberOfBytes = (numberOfRunTypes - 1) / numberOfValuesPerByte + 1;
        // Read defined run IDs with second file stream
        std::ifstream runIdFin(filePath);
        runIdFin.seekg((long)fin.tellg() + (long)numberOfBytes);

        unsigned long long definedId = 0;
        if (bitsPerRunType > 4) {
          for (unsigned i = 0; i < numberOfRunTypes; ++i) {
            hrleSizeType current = 0;
            fin.read((char *)&current, (bitsPerRunType - 1) / 8 + 1);
            if (current == 0) { // defined run
              unsigned long long relativeId = 0;
              runIdFin.read((char *)&relativeId, bytesPerIndex);
              definedId += relativeId;
              runTypes.push_back(definedId);
            } else { // undefined run
              runTypes.push_back(current - 1 + hrleRunTypeValues::UNDEF_PT);
            }
          }
        } else { // if there are several values in each byte
          unsigned readValues = 0;
          unsigned long long bitMask = ((1 << bitsPerRunType) - 1)
                                       << (8 - bitsPerRunType);
          for (unsigned i = 0; i < numberOfBytes; ++i) {
            char byte;
            fin.read(&byte, 1);
            for (int j = 0; j < numberOfValuesPerByte; ++j) {
              if (readValues == numberOfRunTypes)
                break;
              hrleSizeType current = (byte & bitMask) >> (8 - bitsPerRunType);
              byte <<= bitsPerRunType;
              if (current == 0) { // defined run
                unsigned long long relativeId = 0;
                runIdFin.read((char *)&relativeId, bytesPerIndex);
                definedId += relativeId;
                runTypes.push_back(definedId);
              } else { // undefined run
                runTypes.push_back(current - 1 + hrleRunTypeValues::UNDEF_PT);
              }
              ++readValues;
            }
          }
        }
        fin.seekg(runIdFin.tellg());
        runIdFin.close();
      }

      // READ RUNBREAKS
      {
        runBreaks.clear();
        // bitmask sets all higher bytes to FF for negative numbers
        long long bitMask = 0;
        --bitMask <<= bytesPerRunBreak * 8;
        for (unsigned i = 0; i < numberOfRunBreaks; ++i) {
          long long runBreak = 0;
          fin.read((char *)&runBreak, bytesPerRunBreak);
          if (runBreak >> (bytesPerRunBreak * 8 - 1))
            runBreak |= bitMask;
          runBreaks.push_back(runBreak);
        }
      }
    }

    // READ DEFINED VALUES
    {
      // HEADER
      unsigned numberOfDefinedValues = 0, numberOfUndefinedValues = 0;
      fin.read((char *)&numberOfDefinedValues, 4);
      fin.read((char *)&numberOfUndefinedValues, 4);

      // DEFINED VALUES
      for (unsigned i = 0; i < numberOfDefinedValues; ++i) {
        hrleValueType definedValue;
        fin.read((char *)&definedValue, valueTypeByteSize);
        newDomain.domainSegments[0].definedValues.push_back(definedValue);
      }

      // UNDEFINED VALUES
      for (unsigned i = 0; i < numberOfUndefinedValues; ++i) {
        hrleValueType undefinedValue;
        fin.read((char *)&undefinedValue, valueTypeByteSize) >> undefinedValue;
        newDomain.domainSegments[0].undefinedValues.push_back(undefinedValue);
      }
    }

    fin.close();

    // now swap with passed domain
    domain->deepCopy(newDomain);

    return fin.good();
  }
};

#endif // HRLE_DOMAIN_READER_HPP
