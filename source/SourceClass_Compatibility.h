#ifndef SOURCE_CLASS___COMPATIBILITY_H
#define SOURCE_CLASS___COMPATIBILITY_H

#include "SourceClass.h"

#include <vector>
#include <map>
#include <string>

class SourceClass_Compatibility : public SourceClass
{
   public: 
      //! Storage for the individual isotopic injection spectra
      struct specProperties {
         double rigid_br0, rigid_br;
         double g_0, g_1, g_2;
      };

      /** \brief Backwards compatibility source class
       *
       * This should only be used for backwards compatibility.  It is controlled
       * with the following parameters:
       *
       * * inj_spectrum_type:  rigidity | beta_rig | Etot
       * * nuc_rigid_br0: First break of spectrum
       * * nuc_rigid_br: Second break of spectrum
       * * nuc_g_0: First index (below br0)
       * * nuc_g_1: Second index (between br0 and br)
       * * nuc_g_2: Third index (above br)
       * * electron_<...>: Same as above for electrons.  Defaults to nuc values if not provided.
       * * nuc_<...>_ZZ_AAA: Same as above but for different species.  Defaults to nuc values if not provided.
       * * source_normalization: Normalization values for protons at proton_norm_Ekin.  Note this value is also multiplied with the abundance and source model.
       * * electron_source_normalization:  Normalization values for electrons at proton_norm_Ekin.
       * * proton_norm_Ekin: Kinetic energy/nuclei at which to normalize the spectrum.  Also used for electrons
       * * iso_abundance_ZZ_AAA: The isotopic source abundance specified as the normalization relative to source_normalization at proton_norm_Ekin.  Defaults to 0.
       * * source_model:  Integer specifying the source model, see the source.  Defaults to 1.
       * * source_model_electron: Source model for electrons.  Defaults to nuclei sources.
       * * source_parameters: Vector with parameter values for the source model.  Meaning depends on the source model.
       * * source_parameters_electron: Vector with parameter values for the electron source model.  Meaning depends on the source model.
       * * source_values: The values for linear interpolation if source model == 8
       * * source_radius: Radius for the interpolation values if source model == 8
       * * source_values_electron: Same as above for electrons
       * * source_radius_electron: Same as above for electrons
       * * source_xmlFile: File containing the XML specifications of libgalstruct.  Only used if source model == 15.  Same for electrons and nuclei.
       * * source_specification:  0: Use full distribution, 1: Only add at the GC, 2: Only add plane.
       * * n_cr_sources: Number of cr point sources (steady state).  Same for electrons and nuclei and 3D only.
       * * cr_source_x_ii: x coordinate of source number ii (number starts at 0) in kpc
       * * cr_source_y_ii: y coordinate of source number ii (number starts at 0) in kpc
       * * cr_source_z_ii: z coordinate of source number ii (number starts at 0) in kpc
       * * cr_source_w_ii: Gaussian width of source number ii (number starts at 0) in kpc
       * * cr_source_L_ii: Luminosity of source number ii (number starts at 0) relative to source_normalization and abundance.
       * * SNR_events: Set to true for time variable SNR sources.  Only in timestep=2 and 3D.
       * * SNR_interval: time in years between SNRs in 1 kpc^3
       * * SNR_livetime: CR-producing live-time in years of an SNR
       * * SNR_electron_sdg: delta electron source index Gaussian sigma
       * * SNR_nuc_sdg: delta nucleus  source index Gaussian sigma
       * * SNR_electron_dgpivot: delta electron source index pivot rigidity (MeV)   AWS20010410
       * * SNR_nuc_dgpivot: delta nucleus  source index pivot rigidity (MeV)   AWS20010410 
       */
      SourceClass_Compatibility ( utl::Parameters &&pars );

      virtual void setPars( const utl::Parameters &pars ) ;

      virtual void addSource( Particle &particle ) const ;

      virtual void addSource( Particle &particle, double time ) const ;

   private:
      void setDefaultParameters();

      double source_distribution( const double x,
            const double y, 
            const double z, 
            int srcModel, 
            int n_spatial_dimension,
            const std::vector<double> &parameters,
            const std::vector<double> &source_values,
            const std::vector<double> &source_radius) const;

      void create_SNR(const Particle &particle) const;
      void source_SNR_event(Particle &particle, double t) const;

      //This is set when steady state is added because we must remove it for stocastic SNRs.
      mutable bool removeSteadyState;

      //These are changed in const calls because they depend on knowledge of the grid size
      mutable bool createSNRDistributions;
      mutable Distribution SNR_cell_time;           // time between SNR for each cell
      mutable Distribution SNR_cell_phase;          // phase of SNR for each cell
      mutable Distribution SNR_electron_dg;         // electron injection spectral index delta (Gaussian distributed) AWS20010410 
      mutable Distribution SNR_nuc_dg;              // nucleus  injection spectral index delta (Gaussian distributed) AWS20010410

      std::vector<double> createSpecShape(const Particle &particle) const;

      specProperties spDefault;
      specProperties spElectron;

      std::map<std::pair<int, int>, specProperties> iso_inj_spectra;

      std::map<std::pair<int, int>, double> isotopic_abundance;          // isotopic abundances
  
      std::string inj_spectrum_type;             // "rigidity"||"beta_rig"||"Etot": choose the nucleon injection spectrum IMOS20000613.1

      int source_specification;             // 0,1,2
      double source_normalization;          // 1.                                         IMOS20030129
      double electron_source_normalization; // 1.                                         IMOS20031016
  
      double   proton_norm_Ekin;            // kinetic energy for normalization (MeV)
      double   spectra_norm_Rigidity;       // rigidity for spectral normalization (GV)

      int source_model;                     //1= 2= 3= 5= 6=S&Mattox with cutoff
      int source_model_electron;
      std::vector<double> source_parameters; // for source_model=1
      std::vector<double> source_parameters_electron;
      std::vector<double> source_radius, source_values;// for source_model=8
      std::vector<double> source_radius_electron, source_values_electron;// for source_model=8

      std::string source_xmlFile;   // for source_model=15

      int   n_cr_sources;                   // number of pointlike cosmic-ray sources
      std::vector<double> cr_source_x;      // source x positions
      std::vector<double> cr_source_y;      // source y positions
      std::vector<double> cr_source_z;      // source z positions
      std::vector<double> cr_source_w;      // source width sigma in kpc
      std::vector<double> cr_source_L;      // source luminosity in TBD units

      int    SNR_events;                    // handle stochastic SNR events
      double SNR_interval;                  // time in years between SNRs in 1 kpc^3
      double SNR_livetime;                  // CR-producing live-time in years of an SNR
      double SNR_electron_sdg;              // delta electron source index Gaussian sigma
      double SNR_nuc_sdg;                   // delta nucleus  source index Gaussian sigma
      double SNR_electron_dgpivot;          // delta electron source index pivot rigidity (MeV)   AWS20010410
      double SNR_nuc_dgpivot;               // delta nucleus  source index pivot rigidity (MeV)   AWS20010410 

};

#endif
