#ifndef HRLE_DOMAIN_WRITER_HPP
#define HRLE_DOMAIN_WRITER_HPP

#define HRLE_FILE_VERSION_NUMBER 0

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include <hrleIndexType.hpp>
#include <hrleRunTypeValues.hpp>

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
* 4 Bytes  Number of distances
*********************************************************************
*    Values Data    *
************************
* Distances - using sizeOf Bytes per value
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

  char getBitSizeOfNumber(int number) {
    number = std::abs(number);
    int numBits = 0;
    while (number != 0) {
      number /= 2;
      ++numBits;
    }
    return numBits;
  }

  char getByteSizeOfNumber(int number) {
    int numBits = getBitSizeOfNumber(number);
    return numBits / 8 + (numBits % 8) ? 1 : 0;
  }

public:
  hrleDomainWriter() {}

  hrleDomainWriter(hrleDomain *domainPointer) : domain(domainPointer) {}

  hrleDomainWriter(hrleDomain &passedDomain) : domain(&passedDomain) {}

  // setters and getters
  void setDomain(hrleDomain *domainPointer) { domain = domainPointer; }
  void setDomain(hrleDomain &passedDomain) { domain = &passedDomain; }
  hrleDomain &getDomain() { return *domain; }
  // void setValueTypeByteSize(int size) { valueTypeByteSize = size; }
  const int &getValueTypeByteSize() const { return valueTypeByteSize; }
  void setFilePath(std::string path) {
    if (path.find(".hrle") != path.size() - 5) {
      std::cout << "File path not ending in '.hrle'!" << std::endl;
      path.append(".hrle");
      std::cout << "Using '" << path << "' instead." << std::endl;
    }
    filePath = path;
  }
  const std::string &getFilePath() const { return filePath; }

  bool write() {
    std::ofstream fout(filePath);
    if (!fout.is_open()) {
      std::cout << "ERROR: Could not open the file: " << filePath << std::endl;
      return false;
    }

    // write file header
    {
      fout << "HRLE";
      fout << HRLE_FILE_VERSION_NUMBER;
      fout << (bigEndian() ? 1 : 0);
      fout << D;
      fout << char(valueTypeByteSize);
    }

    // GRID PROPERTIES
    hrleIndexType bounds[2 * D];
    domain->getDomainBounds(bounds);
    {
      // find number of bytes needed to represent the highest grid extent
      char gridBoundaryBytes;
      for (unsigned i = 0; i < 2 * D; ++i) {
        if (bounds[i] != 0) {
          gridBoundaryBytes =
              std::max(getByteSizeOfNumber(bounds[i]), gridBoundaryBytes);
        }
      }
      gridBoundaryBytes =
          std::min(gridBoundaryBytes, char(8)); // maximum of 8 Bytes

      // grid properties
      fout << gridBoundaryBytes;
      for (int dim = D - 1; dim >= 0; --dim) {
        fout.write((char *)&bounds[2 * dim], gridBoundaryBytes);
        fout.write((char *)&bounds[2 * dim + 1], gridBoundaryBytes);
        auto boundaryCondition = domain->getGrid().getBoundaryConditions(dim);
        fout.write((char *)&boundaryCondition, 1);
      }
      double gridDelta = domain->getGrid().getGridDelta();
      // std::cout << gridDelta << std::endl;
      fout.write((char *)&gridDelta, sizeof(double));
      // fout << gridDelta;
    }

    bool structureIsSerial = true;
    if (domain->getNumberOfSegments() > 1) {
      domain->serialize();
      structureIsSerial = false;
    }

    // START INDICES, RUN TYPES, RUN BREAKS FOR EACH DIMENSION
    for (int dim = D - 1; dim >= 0; --dim) {
      // get start indices, runbreaks and runtypes
      const std::vector<hrleSizeType> &startIndices =
          domain->domainSegments[0].startIndices[dim];
      const std::vector<hrleSizeType> &runTypes =
          domain->domainSegments[0].runTypes[dim];
      const std::vector<hrleIndexType> &runBreaks =
          domain->domainSegments[0].runBreaks[dim];

      const char bytesPerIndex =
          getByteSizeOfNumber(bounds[2 * dim + 1] - bounds[2 * dim]);

      char bitsPerRunType =
          getBitSizeOfNumber(domain->getNumberOfUndefinedValues());
      if (bitsPerRunType == 3) {
        ++bitsPerRunType;
      } else if (bitsPerRunType > 4 && bitsPerRunType < 8) {
        bitsPerRunType = 8;
      } else if (bitsPerRunType > 8 && (bitsPerRunType % 8)) {
        bitsPerRunType += 8 - bitsPerRunType % 8;
      }

      const char bytesPerRunBreak =
          std::max(getBitSizeOfNumber(bounds[2 * dim]),
                   getBitSizeOfNumber(bounds[2 * dim + 1])) /
              8 +
          1;

      // HRLE BLOCK HEADER
      fout.write((char *)&bytesPerIndex, 1);
      fout.write((char *)&bitsPerRunType, 1);
      fout.write((char *)&bytesPerRunBreak, 1);
      {
        uint32_t numberOfValues = startIndices.size();
        fout.write((char *)&numberOfValues, 4);
        numberOfValues = runTypes.size();
        fout.write((char *)&numberOfValues, 4);
        numberOfValues = runBreaks.size();
        fout.write((char *)&numberOfValues, 4);
      }

      // uint32_t values_written = 0;
      // Write start indices; only save the difference to the next start index
      // (delta encoding)

      // First index is always 0, no need to write explicitly
      for (unsigned int i = 1; i < startIndices.size(); i++) {
        unsigned long diff = startIndices[i] - startIndices[i - 1];
        fout.write((char *)&diff, bytesPerIndex);
        // values_written++;
      }

      // write all runtypes to the file, skipping all segments and indices
      int count = 8 / bitsPerRunType - 1;
      unsigned char byte = 0;
      std::vector<hrleSizeType>
          definedRunIndices; // store all indices of defined runtypes

      // each runType needs at least one byte
      if (bitsPerRunType > 4) {
        for (typename std::vector<hrleSizeType>::const_iterator it =
                 runTypes.begin();
             it != runTypes.end(); ++it) {
          hrleSizeType PtId = 0;
          // if undefined point, need to shift id
          if (*it >= hrleRunTypeValues::UNDEF_PT)
            PtId = (*it) - hrleRunTypeValues::UNDEF_PT + 1;
          else
            definedRunIndices.push_back(*it);

          fout.write((char *)&PtId, (bitsPerRunType - 1) / 8 + 1);
        }
      } else { // can fit more than one value in a byte
        for (typename std::vector<hrleSizeType>::const_iterator it =
                 runTypes.begin();
             it != runTypes.end(); ++it) {
          hrleSizeType PtId = 0;
          if (*it >= hrleRunTypeValues::UNDEF_PT)
            PtId = (*it) - hrleRunTypeValues::UNDEF_PT + 1;
          else
            definedRunIndices.push_back(*it);

          byte |= (PtId << (count * bitsPerRunType));
          --count;

          if (count < 0) { // push byte to stream and start again
            fout << byte;
            count = 8 / bitsPerRunType - 1;
            byte = 0;
          }
        }
        // if last byte is not completely filled, just push it
        if (count >= 0 && count != 8 / bitsPerRunType - 1) {
          fout << byte;
        }
      }

      // Write indices of defined runtypes; only save the difference to the next
      // defined runtype write the first runtype(always 0) explicitly; makes
      // reading easier
      fout.write((char *)&definedRunIndices[0], bytesPerIndex);
      for (unsigned int i = 0; i < definedRunIndices.size() - 1; i++) {
        unsigned long diff = definedRunIndices[i + 1] - definedRunIndices[i];
        fout.write((char *)&diff, bytesPerIndex);
      }

      // Write runbreaks
      for (typename std::vector<hrleIndexType>::const_iterator it =
               runBreaks.begin();
           it != runBreaks.end(); ++it) {
        fout.write((char *)&(*it), bytesPerRunBreak);
      }
    }

    // DEFINED VALUES, UNDEFINED VALUES
    // IMPORTANT NOTE:
    // in order for this to work, you have to make sure that the type
    // you store in the HRLE structure is serializable, i.e.: has an
    // operator overload of the form:
    // ostream& operator<<(ostream&, hrleValueType)
    // HEADER
    {
      uint32_t numberOfDefinedValues = domain->getNumberOfPoints();
      uint32_t numberOfUndefinedValues = domain->getNumberOfUndefinedValues();
      fout.write((char *)&numberOfDefinedValues, 4);
      fout.write((char *)&numberOfUndefinedValues, 4);
    }
    // DATA
    {
      std::vector<hrleValueType> &definedValues =
          domain->domainSegments[0].definedValues;
      for (typename std::vector<hrleValueType>::const_iterator it =
               definedValues.begin();
           it != definedValues.end(); ++it) {
        fout << *it;
      }

      std::vector<hrleValueType> &undefinedValues =
          domain->domainSegments[0].undefinedValues;
      for (typename std::vector<hrleValueType>::const_iterator it =
               undefinedValues.begin();
           it != undefinedValues.end(); ++it) {
        fout << *it;
      }
    }

    // CLEANUP
    fout.close();

    // if the hrleDomain was segmented before, segment it again
    if (!structureIsSerial) {
      domain->segment();
    }

    return fout.good();
  }
};

#endif // HRLE_DOMAIN_WRITER_HPP
