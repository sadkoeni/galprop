#include "SourceClass_UniformSpectra.h"
#include <Registry.h>
#include "PowerLawSD.h"
#include "SmoothPowerLawSD.h"
#include "GalstructXMLDistribution.h"
#include <constants.h>

#include <cmath>
#include <iomanip>

//Registered spectral distributions
static const utl::Registry1<SpectralDistribution,const utl::Parameters&>::Registrar<PowerLawSD> registrarPowerLawSD("PowerLaw");
static const utl::Registry1<SpectralDistribution,const utl::Parameters&>::Registrar<SmoothPowerLawSD> registrarSmoothPowerLawSD("SmoothPowerLaw");

//Registered spatial distributions
static const utl::Registry1<SpatialDistribution,const utl::Parameters&>::Registrar<GalstructXMLDistribution> registrarGalstruct("GalstructXML");

SourceClass_UniformSpectra::SourceClass_UniformSpectra ( utl::Parameters && pars) :
   SourceClass(std::move(pars))
{

   setPars(fpars);

}

void SourceClass_UniformSpectra::setPars(const utl::Parameters &pars)
{

   INFO("Entry");

   //Set all the new parameters
   fpars.setParameters(pars);

   //Spectral distribution
   std::string spectrumName;
   fpars.getParameter("spectrum_type", spectrumName);
   std::ostringstream buf;
   buf<<"Using spectrum of type: "<<spectrumName;
   INFO(buf.str());
   spectralDist = utl::Registry1<SpectralDistribution,const utl::Parameters&>::create(spectrumName, fpars);

   //Spatial distribution
   std::string spatialName;
   fpars.getParameter("spatial_type", spatialName);
   buf.str("");
   buf<<"Using spatial distribution of type: "<<spatialName;
   INFO(buf.str());
   spatialDist = utl::Registry1<SpatialDistribution,const utl::Parameters&>::create(spatialName, fpars);

   //Normalization
   fpars.getParameter("source_normalization", source_normalization);
   buf.str("");
   buf<<"Source normalization: "<<source_normalization;
   INFO(buf.str());
   
   if (source_normalization <= 0) {
      ERROR("source_normalization must be positive");
      throw(std::runtime_error("Fix your configuration"));
   }

   //Abundances
   fpars.getParameter("Elow_abundance", ELowAbundance);
   fpars.getParameter("Ehigh_abundance", EHighAbundance);

   if (ELowAbundance > EHighAbundance) {
      buf.str("");
      buf<<"Elow_abundance cannot be higher than Ehigh_abundance: "<<ELowAbundance<<" > "<<EHighAbundance;
      ERROR(buf.str());
      throw(std::runtime_error("Fix your configuration"));
   }

   //This is very excessive, but isn't called often so we don't have to optimize it.
   for (int iZ(-1); iZ<90; ++iZ) {
      for (int iA(0); iA<3*abs(iZ); ++iA) {
         try {
            double abundance;
            std::ostringstream parname;
            parname << "iso_abundance_" << std::setw(2) << std::setfill('0') << iZ << '_' << std::setw(3) << std::setfill('0') << iA;
            fpars.getParameter(parname.str(), abundance);

            //Must have non-negative abundances, there are no sinks.
            if (abundance >= 0) {
               abundances[std::make_pair(iZ,iA)] = abundance;
            } else {
               WARNING("Abundance cannot be negative");
            }

            //Re-use the parname string stream
            parname.str("");
            parname<<"Abundance for ("<< iZ <<", "<< iA <<") is "<<abundance;
            INFO(parname.str());
         } catch (utl::Parameters::ParameterError) {}
      }
   }

   //Print a warning message if nothing is added
   if (abundances.size() == 0)
      WARNING("No element added, this SourceClass will do nothing.");

   //Get spectral parameters, make sure the defaults are valid
   fpars.getParameter("spectral_pars", defaultSpectralPars);
   if (! spectralDist->validateParameters(&defaultSpectralPars) ) {
      ERROR("Invalid default spectral parameters");
      throw(std::runtime_error("Fix your configuration"));
   }

   buf.str("");
   buf<<"Using default spectral parameters :";
   for (const auto val : defaultSpectralPars)
      buf<<"  "<<val;
   INFO(buf.str());

   //Get specialized parameters and assert validity
   for (auto it = abundances.begin(); it != abundances.end(); ++it) {
      try {
         std::vector<double> tmpPars(defaultSpectralPars);

         const int iZ = it->first.first;
         const int iA = it->first.second;

         std::ostringstream parname;
         parname << "spectral_pars_" << std::setw(2) << std::setfill('0') << iZ << '_' << std::setw(3) << std::setfill('0') << iA;
         fpars.getParameter(parname.str(), tmpPars);

         if (! spectralDist->validateParameters(&tmpPars) ) {
            buf.str("");
            buf<<"Invalid spectral parameters for isotope ("<<iZ<<", "<<iA<<")";
            ERROR(buf.str());
            throw(std::runtime_error("Fix your configuration"));
         }

         spectralPars[it->first] = tmpPars;

         buf.str("");
         buf<<"Using separate spectral parameters for ("<<iZ<<", "<<iA<<") :";
         for (const auto val : tmpPars)
            buf<<"  "<<val;
         INFO(buf.str());

      } catch (utl::Parameters::ParameterError) {}
   }
}

void SourceClass_UniformSpectra::addSource( Particle &particle ) const 
{
   // Ignore secondary particles
   if (isSecondary(particle))
      return;

   // Ignore DM particles
   if (particle.name.find("DM") != string::npos)
      return;

   // Make sure we have non-zero abundance
   const auto type = std::make_pair(particle.Z,particle.A);
   const auto abIt = abundances.find(type);

   if (abIt == abundances.end())
      return;

   //Find the spectral parameters, use a pointer for convenience
   const std::vector<double> *spPars = &defaultSpectralPars;

   const auto spIt = spectralPars.find(type);
   if (spIt != spectralPars.end())
      spPars = &spIt->second;

   std::ostringstream buf;
   buf<<"Using spectral parameters :";
   for (const auto val : *spPars)
      buf<<"  "<<val;
   INFO(buf.str());

   //Normalize the spectrum
   double spNorm(0);
   if (ELowAbundance < EHighAbundance)
      spNorm = spectralDist->integralEkin(spPars, particle, ELowAbundance, EHighAbundance);
   else
      spNorm = spectralDist->valueEkin(spPars, particle, ELowAbundance);

   if (spNorm <= 0) {
      std::ostringstream buf;
      buf<<"Cannot normalize spectrum for particle "<<particle.name;
      ERROR(buf.str());
      throw(std::runtime_error("Fix your configuration"));
   }

   //Apply the abundance and source normalization
   spNorm = source_normalization * abIt->second / spNorm;

   //Convert from momentum per nucleon to total momentum
   if (particle.A > 1)
      spNorm /= particle.A;

   buf.str("");
   buf<<"Spectral normalization for "<<particle.name<<" is "<<spNorm;
   INFO(buf.str());
      
   //Calculate the spectrum
   auto spectra = spectralDist->evaluate(spPars, particle);
   for (size_t i(0); i < spectra.size(); ++i)
      spectra[i] *= spNorm;

   //Fill the primary_source_function.
   if (2 == particle.n_spatial_dimensions) {

#pragma omp parallel for default(shared) schedule(static)
      for (int ir = 0; ir < particle.n_rgrid; ++ir) {
         for (int iz(0); iz < particle.n_zgrid; ++iz) {

            const double sd = (*spatialDist)(particle.r[ir], particle.z[iz]);

            for (int ip(0); ip < particle.n_pgrid; ++ip)
               particle.primary_source_function.d2[ir][iz].s[ip] += sd*spectra[ip];
         }
      }

   } else if (3 == particle.n_spatial_dimensions) {

#pragma omp parallel for default(shared) schedule(static)
      for (int ix = 0; ix < particle.n_xgrid; ++ix) {
         for (int iy(0); iy < particle.n_ygrid; ++iy) {
            for (int iz(0); iz < particle.n_zgrid; ++iz) {

               const double sd = (*spatialDist)(particle.x[ix], particle.y[iy], particle.z[iz]);

               for (int ip(0); ip < particle.n_pgrid; ++ip)
                  particle.primary_source_function.d3[ix][iy][iz].s[ip] += sd*spectra[ip];
            }
         }
      }

   }

}
