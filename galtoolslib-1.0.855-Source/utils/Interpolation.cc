#include "Interpolation.h"

#include <valarray>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <limits>
/* The polynomial approximation is not being used
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
*/

size_t utl::indexBelow(const std::valarray<double> &array, double value) {
   
   //Use standard algorithms
   size_t index = std::upper_bound(&array[0], &array[0]+array.size(), value) - &array[0];

   //Since we want the index that is below the value (or at the value), we reduce it by 1
   if(index > 0) --index;

   return index;
}

size_t utl::indexBelow(const std::vector<double> &array, double value) {
   
   //Use standard algorithms
   size_t index = std::upper_bound(array.begin(), array.end(), value) - array.begin();

   //Since we want the index that is below the value (or at the value), we reduce it by 1
   if(index > 0) --index;

   return index;
}

size_t utl::indexAbove(const std::valarray<double> &array, double value){

   //Use standard algorithms
   size_t index = std::lower_bound(&array[0], &array[0]+array.size(), value) - &array[0];

   //Make sure we return an index to a value within the array
   if (index >= array.size()) index = array.size()-1;

   return index;
}

size_t utl::indexAbove(const std::vector<double> &array, double value){

   //Use standard algorithms
   size_t index = std::lower_bound(array.begin(), array.end(), value) - array.begin();

   //Make sure we return an index to a value within the array
   if (index >= array.size()) index = array.size()-1;

   return index;
}


void utl::findIndexRange(const std::vector<double> &array, double lowValue, double highValue, size_t &lowIndex, size_t &highIndex){

   //This has to be the case
   if ( lowValue > highValue )
      throw std::runtime_error("Lower value should be larger than higher value");

   //Nothing to do
   if ( array.size() == 0 )
      return;

   if ( array.size() == 1 ) {
      lowIndex = highIndex = 0;
      return;
   }

   lowIndex = indexBelow(array,lowValue);
   highIndex = indexAbove(array,highValue);

   //Make sure they are not the same
   if ( lowIndex == highIndex ) {

      if ( lowIndex == 0 )
         ++highIndex;
      else
         --lowIndex;
   }

}

void utl::findIndexRange(const std::valarray<double> &array, double lowValue, double highValue, size_t &lowIndex, size_t &highIndex){

   //This has to be the case
   if ( lowValue > highValue )
      throw std::runtime_error("Lower value should be larger than higher value");

   //Nothing to do
   if ( array.size() == 0 )
      return;

   if ( array.size() == 1 ) {
      lowIndex = highIndex = 0;
      return;
   }

   lowIndex = indexBelow(array,lowValue);
   highIndex = indexAbove(array,highValue);

   //Make sure they are not the same
   if ( lowIndex == highIndex ) {

      if ( lowIndex == 0 )
         ++highIndex;
      else
         --lowIndex;
   }

}

utl::Interpolation1D::Interpolation1D (
      const std::vector<double> &xValues, 
      const std::vector<double> &yValues
      ) :
   fxValues(xValues),
   fyValues(yValues)
{
   if (fxValues.size() != fyValues.size())
      throw std::length_error("x and y arrays should have sime size in interpolation");
   
   if (fxValues.size() < 2)
      throw std::length_error("need at least 2 values for x and y arrays in interpolation");
}

utl::LinearInterpolation::LinearInterpolation (
      const std::vector<double> &xValues,
      const std::vector<double> &yValues
      ) :
   Interpolation1D(xValues, yValues),
   slope(fxValues.size()-1, std::numeric_limits<double>::min())
{
}

double utl::LinearInterpolation::operator () ( double x ) const {

   //Find the i to the nearest point
   size_t i = indexBelow(fxValues, x);

   if ( i == slope.size() )
      --i;

   if ( slope[i] == std::numeric_limits<double>::min() )
      slope[i] = (fyValues[i+1]-fyValues[i])/(fxValues[i+1]-fxValues[i]);

   return fyValues[i] + (x-fxValues[i])*slope[i];

}

utl::PowerLawInterpolation::PowerLawInterpolation (
      const std::vector<double> &xValues,
      const std::vector<double> &yValues
      ) :
   Interpolation1D(xValues, yValues),
   index(fxValues.size()-1, std::numeric_limits<double>::min())
{
}

double utl::PowerLawInterpolation::operator () ( double x ) const {
   
   //Make sure we have positive values
   if ( x <= 0 )
      throw (std::domain_error("x should be positive when evaluating power law interpolation"));

   //Find the index to the nearest point
   size_t i = indexBelow(fxValues, x);

   if ( i == index.size() )
      --i;

   if ( index[i] == std::numeric_limits<double>::min() ) {

      //Make sure we have strictly positive values
      if ( fyValues[i] <= 0 || fyValues[i+1] <=0 || fxValues[i] <=0 || fxValues[i+1] <= 0 )
         throw (std::domain_error("x and y interpolation points should be strictly positive"));

      index[i] = log(fyValues[i+1]/fyValues[i])/log(fxValues[i+1]/fxValues[i]);

   }

   return fyValues[i] * pow (x/fxValues[i], index[i]);

}

/*
//Compute the coefficients, c_j, for the polynomial p(x) \sum_{j=0}^{n-1} c_j*x^j
//fulfilling the conditions p(x_j) = y_j for all j=0,n-1 where n is the number
//of points in the input arrays
void utl::polynomialApproximation(const std::valarray<double> & x, const std::valarray<double> &y, std::valarray<double> &coeff){
   if (x.size() != y.size()) {
      std::cout<<"Sizes of arrays not the same in polynomial approximation"<<std::endl;
      return;
   }

   //Solve the equations p(x_j) = y_j for all j
   //create the power of x array
   std::vector<double> data_x(x.size()*x.size());
   for (size_t i = 0; i < x.size(); ++i){
      for (size_t j = 0; j < x.size(); ++j){
	 data_x[x.size()*i+j] = pow(x[i],j);
      }
   }
   std::vector<double> data_y(y.size());
   for (size_t i = 0; i < y.size(); ++i){
      data_y[i] = y[i];
   }

   //Create gsl matrixes
   gsl_matrix_view m = gsl_matrix_view_array (&data_x[0], x.size(), x.size());
   gsl_vector_view b = gsl_vector_view_array (&data_y[0], x.size());

   //The output vector
   gsl_vector *r = gsl_vector_alloc (x.size());

   int s;
   gsl_permutation * p = gsl_permutation_alloc (x.size());
   gsl_linalg_LU_decomp (&m.matrix, p, &s);
   gsl_linalg_LU_solve (&m.matrix, p, &b.vector, r);
     
   //Put the output in the coeff array
   coeff.resize(x.size());
   for (size_t i = 0; i < x.size(); ++i){
      coeff[i] = gsl_vector_get(r,i);
   }
}
*/

void utl::lagrangeCoeff4(const std::valarray<double> & x,const std::valarray<double> & y, std::valarray<double> & coeff){
	if (x.size() != 4 || y.size() != 4) 
           throw(std::runtime_error("Cannot compute Lagrange coefficient unless sizes are == 4"));

	coeff.resize(4);
	double tmp = y[0]/((x[0]-x[1])*(x[0]-x[2])*(x[0]-x[3]));
	coeff[0] = tmp * (- ( x[1] * x[2] * x[3] ));
	coeff[1] = tmp * (x[1]*x[2] + x[1]*x[3] + x[2]*x[3]);
	coeff[2] = tmp * (-x[1] - x[2] - x[3]);
	coeff[3] = tmp;
	tmp = y[1]/((x[1]-x[0])*(x[1]-x[2])*(x[1]-x[3]));
	coeff[0] += tmp * (- ( x[0] * x[2] * x[3]));
	coeff[1] += tmp * (x[0]*x[2] + x[0]*x[3] + x[2]*x[3]);
	coeff[2] += tmp * (-x[0] - x[2] - x[3]);
	coeff[3] += tmp;
	tmp = y[2]/((x[2]-x[0])*(x[2]-x[1])*(x[2]-x[3]));
	coeff[0] += tmp * (- ( x[0] * x[1] * x[3]));
	coeff[1] += tmp * (x[0]*x[1] + x[0]*x[3] + x[1]*x[3]);
	coeff[2] += tmp * (-x[0] - x[1] - x[3]);
	coeff[3] += tmp;
	tmp = y[3]/((x[3]-x[0])*(x[3]-x[1])*(x[3]-x[2]));
	coeff[0] += tmp * (- ( x[0] * x[1] * x[2]));
	coeff[1] += tmp * (x[0]*x[1] + x[0]*x[2] + x[1]*x[2]);
	coeff[2] += tmp * (-x[0] - x[1] - x[2]);
	coeff[3] += tmp;
}

