#ifndef HRLE_UTIL
#define HRLE_UTIL

namespace hrleUtil {
inline unsigned getBitSizeOfNumber(long number) {
  number = std::abs(number);
  char bitSize = 0;
  while (number != 0) {
    number >>= 1;
    ++bitSize;
  }
  // one additional bit for the sign
  return bitSize + 1;
}

inline char getByteSizeOfNumber(const long number) {
  return (getBitSizeOfNumber(number) + 7) / 8;
}
} // namespace hrleUtil
#endif