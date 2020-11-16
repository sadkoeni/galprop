#ifndef _integ_h_
#define _integ_h_

#include<vector>
#include<valarray>
#include<cmath>

namespace integ {
template <typename T>
class vecFun {
public:
   virtual ~vecFun() {}
   virtual std::vector<T> operator () ( const double x ) const = 0;
};

/** \brief Simpson integration method for function of single variable returning vectors of valarrays
 * 
 * \param lbound gives lower bound and \param ubound gives upper bound for integration
 * \param step is the absolute step size
 * \param relEps gives relative precision and \param absEps gives absolute precision
 * if \param relEps and absEps are both 0, the program uses a fixed step size.
*/
template <typename T>
std::vector<T> simvec ( const double lbound, const double ubound, double step, double relEps, double absEps, const vecFun<T> &fu, bool useCache=true );

inline double max ( double value ) { return value; }
inline double max ( const std::valarray<double> &value ) { return value.max(); }

inline double fabs ( const double value ) { return std::abs(value); }
inline std::valarray<double> fabs ( const std::valarray<double> &value ) {return abs(value); }

template <typename T>
std::vector<T> simvec ( const double lbound, const double ubound, double step, const double minStep, double relEps, double absEps, const vecFun<T> &fu, bool useCache ) {
    step = ubound >= lbound ? fabs ( step ) : -fabs ( step );
    const double sign = step >= 0 ? 1. : -1.;
    double x ( lbound );
    bool C ( false );

    //Set up the cache
    std::vector<bool> validCache(7,false);
    std::vector<T> f0( fu(x) );
    std::vector<T> integ( f0 ), absInteg( f0 );
    std::vector< std::vector<T> > fvalue(7, f0);
    std::vector<T> DI1 ( f0 );
    std::vector<T> DI2 ( f0 );
    std::vector<T> DI3 ( f0 );

    //Trick needed to properly initialize if T is a valarray
    for ( size_t i = 0; i < f0.size(); ++i ) {
        integ[i] = 0;
        absInteg[i] = 0;
    }

    std::vector<double> P(5);
    P[1]=P[3]=4.;
    P[2]=2.;
    P[4]=1.;

    if ( ubound==lbound )
        return integ;

    relEps = fabs ( relEps );
    absEps = fabs ( absEps );
    while ( ! C ) {
        const double x0 = x;
        
        //Check if we have reached the boundary and fix the step size in that
        //case.
        if ( ( x0 + 4.*step - ubound ) *sign > 0 ) {
            step = ( ubound - x0 ) /4.;
            if ( step == 0 )
                return integ;
            //All cached values are invalid in this case, delete them
            for ( int k = 1; k < 7; ++k ) {
               validCache[k] = false;
            }
            C=true;
        }

        //The first value is always pre-calculated
        for ( size_t i = 0; i < DI3.size(); ++i ) {
            DI2[i] = fvalue[0][i];
            DI3[i] = fabs ( fvalue[0][i] );
        }

        // Calculate the function value and the integral with finer step size
        // needed for the error calculation
        for ( int k = 1; k < 5; ++k ) {
            x += step;
            if ( ( x-ubound ) *sign >= 0 )
                x = ubound;
            //Only calculate non-cached values
            if ( ! validCache[k] ) {
                fvalue[k] = fu ( x );
                validCache[k] = useCache;
            }
            for ( size_t i = 0; i < DI2.size(); ++i ) {
                DI2[i] = DI2[i] + P[k]* fvalue[k][i];
                DI3[i] = DI3[i] + P[k]*fabs ( fvalue[k][i] );
            }
        }

        // Doing the actual integral
        for ( size_t i = 0; i < DI1.size(); ++i ) {
            DI1[i] = ( fvalue[0][i] + 4.* fvalue[2][i] + fvalue[4][i] ) *2./3.*step;
            DI2[i] *= step / 3.;
            DI3[i] *= step / 3.;
        }

        // If any value is not up to accuracy we don't accept step
        bool stepValid (true);
        bool stepIncrease (true);
        std::vector<T> DELTA (DI3), EPS (DI3), DIFF (DI3);
        for ( size_t i = 0; i < DI1.size() && stepValid; ++i ) {
           DELTA[i] = fabs ( DI2[i] - DI1[i] );
           EPS[i] = absInteg[i] + DI3[i];
           EPS[i] *= relEps;
           //Make sure DELTA is less then EPS in all cases
           DIFF[i] = DELTA[i]-EPS[i];
           if ( max (DIFF[i]) > 0 ) {
              // make sure we are also above the absolute error
              DIFF[i] = DELTA[i] - absEps;
              if ( max (DIFF[i]) > 0 )
                 stepValid = false;
           }
           if (stepIncrease) {
              DIFF[i] = DELTA[i] - EPS[i]/8.;
              if ( max (DIFF[i]) > 0 ) {
                 // make sure we are also above the absolute error
                 DIFF[i] = DELTA[i] - absEps/8.;
                 if ( max (DIFF[i]) > 0 )
                    stepIncrease = false;
              }
           }
        }

        //The error is within range and we have a successful step.  If errors are both 0 we have a fixed step size
        if ( stepValid || ( relEps == 0 && absEps == 0 ) ) {
            //The step is good enough to increase the step size
            if ( stepIncrease && ( relEps != 0 || absEps != 0) ) {
                step = step*2.;
                for ( int k = 0; k < 3; ++k ) {
                   fvalue[k] = fvalue[k+4];
                   validCache[k] = validCache[k+4];
                }
                for ( int k = 3; k < 7; ++k ) {
                    validCache[k] = false;
                }
            } else {
                fvalue[0] = fvalue[4];
                validCache[0] = validCache[4];
                fvalue[2] = fvalue[5];
                validCache[2] = validCache[5];
                fvalue[4] = fvalue[6];
                validCache[4] = validCache[6];
                validCache[1] = false;
                validCache[3] = false;
                validCache[5] = false;
                validCache[6] = false;
   	    }
            for ( size_t i = 0; i < DI1.size(); ++i ) {
                integ[i] += DI2[i] + ( DI2[i]-DI1[i] ) /15.;
                absInteg[i] += DI3[i];
            }
        //Step wasn't successful, reduce the step size
        } else {
            //Have a minimum step size
            if ( step/2.0 > minStep ) {
               step = step/2.0;
               fvalue[6] = fvalue[4];
               validCache[6] = validCache[4];
               fvalue[5] = fvalue[3];
               validCache[5] = validCache[3];
               fvalue[4] = fvalue[2];
               validCache[4] = validCache[2];
               fvalue[2] = fvalue[1];
               validCache[2] = validCache[1];
               validCache[1] = false;
               validCache[3] = false;
               x = x0;
               C = false;
            } else {
               //Continue with an error message
               std::cerr<<"Minimum step size reached in simvec, accuracy may not be reached for all elements"<<std::endl;
               fvalue[0] = fvalue[4];
               validCache[0] = validCache[4];
               fvalue[2] = fvalue[5];
               validCache[2] = validCache[5];
               fvalue[4] = fvalue[6];
               validCache[4] = validCache[6];
               validCache[1] = false;
               validCache[3] = false;
               validCache[5] = false;
               validCache[6] = false;
            }
        }
    }
    return integ;
}
}

#endif
