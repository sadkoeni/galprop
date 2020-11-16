/* 
 * Base class for propagation of particles.
 * This defines the interface for derived classes.
 * It is actually fairly simple, we only need a creator and an operator
 *
 * Currently this is designed around the old propel routine, but it should be flexible enough for most
 */

#include "Particle.h"
#include "Galdef.h"

class PropelBase {
 protected:
  const Galdef &galdef; //and the Galdef class
  
 public:
  //The particle is just to get the size of the grid
  //We should have a better way of handling the grid
  PropelBase ( const Particle &particle,  const Galdef &_galdef ) :
  galdef(_galdef)
  { }
  virtual ~PropelBase(){}
  
  //This is non-constant because we will likely need to change some of the temporaries
  virtual void operator() (Particle &particle) = 0;
  
};
