/**
 * \file RandomType.hpp
 * \brief Class to hold bit-width and unsigned type
 *
 * This provides a simple class to couple a bit-width and an unsigned type
 * capable of holding all the bits.  In addition is offers static methods for
 * I/O and checksumming.
 *
 * Copyright (c) Charles Karney (2006-2011) <charles@karney.com> and licensed
 * under the LGPL.  For more information, see http://randomlib.sourceforge.net/
 **********************************************************************/
#if !defined(RANDOMLIB_RANDOMTYPE_HPP)
#define RANDOMLIB_RANDOMTYPE_HPP "$Id: 109954f78e4584c0d87a5f99bdf9947f4fdbe651 $"

#include <limits>
#include <string>
#include <iostream>

namespace RandomLib {
  /**
   * \brief Class to hold bit-width and unsigned type
   *
   * This provides a simple class to couple a bit-width and an unsigned type
   * capable of holding all the bits.  In addition is offers static methods for
   * I/O and checksumming.
   **********************************************************************/
  template<int bits, typename UIntType>
  class RANDOM_EXPORT RandomType {
  public:
    /**
     * The unsigned C++ type
     **********************************************************************/
    typedef UIntType type;
    /**
     * The number of significant bits
     **********************************************************************/
    static const unsigned width = bits;
    /**
     * A mask for the significant bits.
     **********************************************************************/
    static const type mask =
      ~type(0) >> (std::numeric_limits<type>::digits - width);
    /**
     * The minimum representable value
     **********************************************************************/
    static const type min = type(0);
    /**
     * The maximum representable value
     **********************************************************************/
    static const type max = mask;
    /**
     * A combined masking and casting operation
     **********************************************************************/
    template<typename IntType> static type cast(IntType x) throw()
    { return type(x) & mask; }
    /**
     * Read a data value from a stream of 32-bit quantities (binary or text)
     **********************************************************************/
    static void Read32(std::istream& is, bool bin, type& x);
    /**
     * Read the data value to a stream of 32-bit quantities (binary or text)
     **********************************************************************/
    static void Write32(std::ostream& os, bool bin, int& cnt, type x);
    /**
     * Accumulate a checksum of a integer into a 32-bit check.  This implements
     * a very simple checksum and is intended to avoid accidental corruption
     * only.
     **********************************************************************/
    static void CheckSum(type n, uint32_t& check) throw();
  };

  /**
   * The standard unit for 32-bit quantities
   **********************************************************************/
  typedef RandomType<32, uint32_t> Random_u32;
  /**
   * The standard unit for 64-bit quantities
   **********************************************************************/
  typedef RandomType<64, uint64_t> Random_u64;

  /**
   * The integer type of constructing bitsets.  This used to be unsigned long.
   * C++0X has made this unsigned long long.
   **********************************************************************/
#if defined(_MSC_VER) && _MSC_VER >= 1600
  typedef unsigned long long bitset_uint_t;
#else
  typedef unsigned long bitset_uint_t;
#endif

  /// \cond SKIP

  // Accumulate a checksum of a 32-bit quantity into check
  template<>
  inline void Random_u32::CheckSum(Random_u32::type n, Random_u32::type& check)
    throw() {
    // Circular shift left by one bit and add new word.
    check = (check << 1 | (check >> 31 & Random_u32::type(1))) + n;
    check &= Random_u32::mask;
  }

  // Accumulate a checksum of a 64-bit quantity into check
  template<>
  inline void Random_u64::CheckSum(Random_u64::type n, Random_u32::type& check)
    throw() {
    Random_u32::CheckSum(Random_u32::cast(n >> 32), check);
    Random_u32::CheckSum(Random_u32::cast(n      ), check);
  }
  /// \endcond

} // namespace RandomLib

#endif  // RANDOMLIB_RANDOMTYPE_HPP
