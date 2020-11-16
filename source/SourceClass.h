#ifndef SOURCE_CLASS_H
#define SOURCE_CLASS_H

#include <Registry.h>
#include <Parameters.h>

#include "Particle.h"

#include <memory>

/** \brief Interface class for source classes in galprop.
 *
 * Each source class provides injection spectrum for arbitrary types
 * of CR particles.  There is in general no limitations to the 
 * flexibility of source classes.
 *
 * Galprop expects the units for the injection to be density per unit total momentum 
 * per unit time [cm^-3 MeV^-1 s^-1].  If specified, these are normalized such that the flux of protons at 
 * a specified energy at the solar location is a given number.  All primary fluxes, apart 
 * from primary electrons and positrons, are adjusted with the same number.  There is a 
 * special treatement for electrons and positrons, they are normalized such that the total
 * electrons at a specified energy at the solar location is a given number.
 *
 * Note that it is not possible to mix relative and absolute source classes.  That would 
 * need the memory and CPU to store and propagate separately contributions from fixed and 
 * relative source classes.
 *
 * Uses the Registry.h class to create instances.  There is a reserved parameter called
 * "SourceClassType" that gives the name of the class in the registry.  This allows for 
 * easy extension of the capabilities of GALPROP, e.g. to hook it up with DM sources.
 */
class SourceClass {

   public:
      //! Should be the only initializer.  Set everything up from the parameters object.
      SourceClass( utl::Parameters &&pars) : fpars(std::move(pars)) {}
      virtual ~SourceClass() {}


      //! Create a source class from the parameters object.  The type of class is set with the SourceClassType parameter.
      static std::unique_ptr<SourceClass> create( utl::Parameters && pars );

      //! Get the current parameters
      const utl::Parameters& getPars() const { return fpars; }

      /** \brief Set the parameters of the object
       *
       * The parameters object should be updated to add the new parameters.
       */
      virtual void setPars( const utl::Parameters &pars ) = 0;

      /** \brief Add the steady state source density to the particle.
       *
       * The routine should add to the primary_source_function only.  It is
       * expected that this routine is called using the same particle for many
       * different source classes.
       */
      virtual void addSource( Particle &particle ) const = 0;

      /** \brief Add the time varying source to the particle
       *
       * The routine should add to the primary_source_function only.  It is
       * expected that this routine is called using the same particle for many
       * different source classes.
       *
       * This should be in addition to any specified steady state source.
       * Currently only called for timestep == 2 and 3D.
       * time is in years.
       */
      virtual void addSource( Particle &particle, double time ) const = 0;

   protected:
      utl::Parameters fpars;

      //! Convenience function to throw away special particles
      bool isSecondary(const Particle &particle) const;

};

#endif
