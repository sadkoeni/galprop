/**
 * \file ExactPower.hpp
 * \brief Header for ExactPower
 *
 * Sample exactly from a power distribution.
 *
 * Copyright (c) Charles Karney (2006-2011) <charles@karney.com> and licensed
 * under the LGPL.  For more information, see http://randomlib.sourceforge.net/
 **********************************************************************/

#if !defined(RANDOMLIB_EXACTPOWER_HPP)
#define RANDOMLIB_EXACTPOWER_HPP "$Id: 7fd110f02deb7128b3427362c500fd894f09a869 $"

#include <RandomLib/RandomNumber.hpp>

namespace RandomLib {
  /**
   * \brief Sample exactly from a power distribution.
   *
   * Sample exactly from power distribution (<i>n</i> + 1)
   * <i>x</i><sup><i>n</i></sup> for \e x in (0,1) and integer \e n >= 0 using
   * infinite precision.  The template parameter \a bits specifies the number
   * of bits in the base used for RandomNumber (i.e., base =
   * 2<sup><i>bits</i></sup>).
   **********************************************************************/
  template<int bits = 1> class ExactPower {
  public:
    /**
     * Return the random deviate with a power distribution, (<i>n</i> + 1)
     * <i>x</i><sup><i>n</i></sup> for \e x in (0,1) and integer \e n >= 0.
     * Returned result is a RandomNumber with base 2<sup><i>bits</i></sup>.
     * For \e bits = 1, the number of random bits in the result and consumed
     * are as follows:\n
     \verbatim
      n    random bits
          result  consumed
      0   0        0
      1   2        4
      2   2.33     6.67
      3   2.67     9.24
      4   2.96    11.71
      5   3.20    14.11
      6   3.41    16.45
      7   3.59    18.75
      8   3.75    21.01
      9   3.89    23.25
     10   4.02    25.47
     \endverbatim
     * The relative frequency of the results with \a bits = 1 can be shown via
     * a histogram\n <img src="powerhist.png" width=580 height=750 alt="exact
     * binary sampling of power distribution">\n The base of each rectangle
     * gives the range represented by the corresponding binary number and the
     * area is proportional to its frequency.  A PDF version of this figure
     * <a href="powerhist.pdf">here</a>.  This allows the figure to be
     * magnified to show the rectangles for all binary numbers up to 9 bits.
     **********************************************************************/
    template<class Random>
    RandomNumber<bits> operator()(Random& r, unsigned n) const;
  };

  template<int bits> template<class Random> RandomNumber<bits>
  ExactPower<bits>::operator()(Random& r, unsigned n) const {
    // Return max(u_0, u_1, u_2, ..., u_n).  Equivalent to taking the
    // (n+1)th root of u_0.
    RandomNumber<bits> x;
    for (RandomNumber<bits> y; n--;) {
      y.Init();
      if (x.LessThan(r, y))
        x = y;
    }
    return x;
  }
} // namespace RandomLib

#endif  // RANDOMLIB_EXACTPOWER_HPP
