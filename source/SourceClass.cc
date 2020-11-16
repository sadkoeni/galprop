#include "SourceClass.h"
#include "SourceClass_Compatibility.h"
#include "SourceClass_UniformSpectra.h"

#include <string>

#include <exception>

#include <Registry.h>

const static utl::Registry1<SourceClass,utl::Parameters&&>::Registrar<SourceClass_Compatibility> registrarCompatibility("Compatibility");
const static utl::Registry1<SourceClass,utl::Parameters&&>::Registrar<SourceClass_UniformSpectra> registrarUS("UniformSpectra");

std::unique_ptr<SourceClass> SourceClass::create( utl::Parameters && pars ) 
{
   std::string name;
   try {
      pars.getParameter("SourceClassType", name);
   } catch (utl::Parameters::ParameterError) {
      ERROR("SourceClassType parameter missing.");
      utl::Registry1<SourceClass,utl::Parameters&&>::QueryRegistered();
      throw(std::runtime_error("Fix the parameter file"));
   }
   return utl::Registry1<SourceClass,utl::Parameters&&>::create(name, std::move(pars));
}

bool SourceClass::isSecondary(const Particle &particle) const
{
   //All the special secondary particles
   if ( 
         particle.name == "secondary_positrons" ||
         particle.name == "knock_on_electrons" ||
         particle.name == "secondary_electrons" ||
         particle.name == "tertiary_antiprotons" ||
         particle.name == "secondary_antiprotons" ||
         particle.name == "secondary_protons" 
         )
      return true;

   //K_capture nuclei are never primary
   if ( particle.K_electron > 0 ) 
      return true;

   return false;
}
