#ifndef _utl_Interpolation_h
#define _utl_Interpolation_h

#include <valarray>
#include <vector>


namespace utl {
	/** \brief Assuming array is sorted ascendingly, find the index where the array value is <= the searched for value.
	 *
	 * \param array the ascendingly sorted array
	 * \param value the value to search for
	 *
	 * \note Returns the lowest index in the array if value is lower than all values in the array
	 */
	size_t indexBelow(const std::valarray< double >& array, double value);
	/** \overload */
        size_t indexBelow(const std::vector< double >& array, double value);

	/** \brief Assuming array is sorted ascendingly, find the index where the array value is >= the searched for value.
	 *
	 * \param array the ascendingly sorted array
	 * \param value the value to search for
	 * 
	 * \note Returns the highest index in the array if value is higher than all values in the array
	 */
	size_t indexAbove(const std::valarray< double >& array, double value);
	/** \overload */
        size_t indexAbove(const std::vector< double >& array, double value);

	/** \brief Find the index range which brackets the range from lowValue to highValue
	 *
	 * \param array the ascendingly sorted array to search
	 * \param lowValue the lower limit of the range to bracket
	 * \param highValue the upper limit of the range to bracket
	 * \param lowIndex a reference to a integer, which will contain the lower index in the array
	 * \param highIndex a reference to a integer, which will contain the higher index in the array
	 *
	 * If the range is outside the array span, the index will be in the end closer to the range.  The
	 * index will only be the same if the array size is less than 2.  Throws an exception if lowValue 
         * is higher than highValue.
	 */
	void findIndexRange(const std::vector< double >& array, double lowValue, double highValue, size_t& lowIndex, size_t& highIndex);
	/** \overload */
	void findIndexRange(const std::valarray< double >& array, double lowValue, double highValue, size_t& lowIndex, size_t& highIndex);

        /** \brief Base interpolation class for 1D interpolation
         *
         * xValues should be sorted ascendingly. Behaviour is undefined if not fulfilled.
         *
         * The interpolators should be lazy, such that evalutating a single value
         * should be fast but they should also cache results so that evaluating 
         * many values is also fast.
         *
         * The x and y values are only stored as references in the class, so they should not be temporary variables
         *
         * Extrapolation happens automatically if possible.
         */
        class Interpolation1D {
           public:
              ///Initialize the interpolator
              Interpolation1D(const std::vector<double> &xValues, const std::vector<double> &yValues);

              ///Return interpolated value at x.  Do extrapolation as well.
              virtual double operator () (double x) const = 0;
           protected:
              const std::vector<double>& fxValues, fyValues;
        };
           

	/// Linear interpolation of values 
        class LinearInterpolation : public Interpolation1D {
           public:
              LinearInterpolation(const std::vector<double> &xValues, const std::vector<double> &yValues);

              virtual double operator () (double x) const;
           private:
              mutable std::vector<double> slope; //Mutable so operator() can be const.
        };

	/// Power law interpolation of values 
        class PowerLawInterpolation : public Interpolation1D {
           public:
              PowerLawInterpolation(const std::vector<double> &xValues, const std::vector<double> &yValues);

              virtual double operator () (double x) const;
           private:
              mutable std::vector<double> index; //Mutable so operator() can be const.
        };

	/** \brief Calculate the n-1th order polynomial coefficients, given n numbers of x and y values.
	 *
	 * \param x the x values to use
	 * \param y the y values to use
	 * \param coeff on return, stores the calculated coefficients.  The coefficients are stored in order, so the total
	 * function will be \f$ f = \sum^i coeff_i x^i \f$
	 */
	//void polynomialApproximation(const std::valarray<double> & x,const std::valarray<double> & y, std::valarray<double> & coeff);

	/** \brief Calculate the 3rd order Lagrange Coefficients, given x and y values.
	 *
	 * \param x the x values to use
	 * \param y the y values to use
	 * \param coeff on return, stores the calculated coefficients.  The coefficients are stored in order, so the total
	 * function will be \f$ f = \sum^i coeff_i x^i \f$
	 */
	void lagrangeCoeff4(const std::valarray<double> & x,const std::valarray<double> & y, std::valarray<double> & coeff);

}
#endif
