#ifndef SOURCE_CLASS_UNIFORM_SPECTRA_H
#define SOURCE_CLASS_UNIFORM_SPECTRA_H

#include "SourceClass.h"
#include "SpatialDistribution.h"
#include "SpectralDistribution.h"
#include <map>

/** \brief Source class with uniform spectra throughout the galaxy and 
 * only the normalization varies.
 *
 * Each species has the same spectra throughout the galaxy and the spatial distribution
 * of the normalization of all species is the same.  The relative abundance can be adjusted
 * and the spectral parameters are not identical.  This should accomodate most homogeneous 
 * source classes where we don't expect their properties to vary throughout the Galaxy.
 *
 * Built up of two main components, the spatial distribution and the spectral form.
 * Only one spatial distribution with a single set of parameters is allowed, but each 
 * isotope can have its own spectral parameters and abundance.  Abundance is defined 
 * over a user defined range in kinetic energy per nucleon.  This should be more stable
 * than the old method of using a single normalization.
 *
 * Unlike the classic GALPROP method primary electrons and positrons are treated on 
 * equal grounds with the nuclei and have their abundance and spectral parametes set
 * using the same mechanism.  Note though that the automatic post-process galprop 
 * normalization handles electrons and positrons differently from protons and nuclei.
 */
class SourceClass_UniformSpectra : public SourceClass
{
   public:
      /** \brief Constructs the class from a parameters object.
       *
       * Uses the following parameters:
       * * spectrum_type: Name of the spectrum class to use
       * * spatial_type: Name of the spatial class to use
       * * spectral_pars: Default spectral parameters to use.  Vector of doubles.
       * * source_normalization: Normalization parameter for the source class.
       * * iso_abundance_ZZ_AAA: Isotopic abundance of isotope with charge Z and mass number A.  Use A=0 for electrons and positrons.
       * * spectral_pars_ZZ_AAA: Spectral parameters to use for isotope with charge Z and mass number A.  Vector of doubles.
       * * Elow_abundance: Lower kinetic energy per nucleon in MeV/nuc for the range at which we determine the abundance ratios
       * * Ehigh_abundance: Higher kinetic energy per nucleon in MeV/nuc for the range at which we determine the abundance ratios
       * 
       * The meaning of spectral_pars depends on the spectrum_type selected.
       * The parameters object is also used to initialize the spectrum and spatial
       * distributions which may require more parameters, consult the documentation
       * for each class.
       *
       * There are no default values, they must all be specified.
       *
       * The source is normalized such that
       * iso_abundance_ZZ_AAA = int_Elow^Ehigh Ekin^2 Q(Sun,Ekin) dEkin / source_normalization
       * where Q(sun,Ekin) is the value of the source distribution at the location 
       * of the sun in units of (momentum per nucleon)-1.
       *
       * Only isotopes with positive iso_abundance value are used.
       * Options in galprop can turn off various isotopes globally, i.e.
       * use_Z, primary_electrons, primary_positrons.
       */
      SourceClass_UniformSpectra(utl::Parameters &&pars);

      /** \brief Sets the parameters from the new parameters object
       *
       * Overrides only the parameters given. Silently ignores 
       * irrelevant parameters.
       */
      virtual void setPars( const utl::Parameters &pars ) ;

      //! Adds the specified distribution to particle.
      virtual void addSource( Particle &particle ) const ;

      //! No option for time dependent source distribution.
      virtual void addSource( Particle &particle, double time ) const {}

   protected:
      std::unique_ptr<SpatialDistribution> spatialDist;  
      std::unique_ptr<SpectralDistribution> spectralDist;

      double ELowAbundance, EHighAbundance;
      double source_normalization;
      std::vector<double> defaultSpectralPars;
      std::map<std::pair<int, int>, double> abundances;
      std::map<std::pair<int, int>, std::vector<double> > spectralPars;
};

#endif
