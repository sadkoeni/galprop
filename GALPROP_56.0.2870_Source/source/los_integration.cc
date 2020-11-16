#include "los_integration.h"
#include <valarray>
#include <iostream>

using namespace std;

template <> void SM::LOSintegrator<valarray<double> >::output (valarray<double> t) const{
   for (int i = 0; i < t.size(); ++i) {
      cout<<t[i]<<", "<<endl;
   }
}


