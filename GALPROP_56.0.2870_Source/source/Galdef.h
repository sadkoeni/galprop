
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Galdef.h *                                    galprop package * 10/12/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef Galdef_h
#define Galdef_h

#include "ErrorLogger.h"
#include <Parameters.h>
#include "SourceClass.h"

#include <cstdio>     //AWS20050624
#include <string>     //AWS20050624
#include <valarray>
#include <vector>
#include <map>
#include <sstream>
#include <ostream>
#include <algorithm>
#include <iterator>
#include <memory>

//The name of the Parameters object
#define PARAMETERS_NAME pars
//Set a variable parameterName to the value from a parameter with the same name.  Parameters is the name of the Parameters object.
//This will throw an exception if not defined
#define SETPAR(parameterName) (PARAMETERS_NAME).getParameter(#parameterName,parameterName) 
//Optionally do the same as a above and issue a WARNING if parameter not found
#define SETPAROPT(parameterName) try{SETPAR(parameterName);}catch(...){ostringstream buf; buf << "Using default value of '"<<parameterName<<"' for '"<<#parameterName<<"'"; WARNING(buf.str());}
//Try to read parameter, otherwise set to default value and issue a WARNING
#define SETPARDEF(parameterName,value) try{SETPAR(parameterName);}catch(...){ostringstream buf; buf << "Using default value of '"<<(value)<<"' for '"<<#parameterName<<"'"; WARNING(buf.str()); parameterName = (value);}
//Try to read a parameter but fall back to another name if first not found and
//throw an exception if neither are found
#define SETPARNAME(parameterName,fallback) try{SETPAR(parameterName);}catch(...){(PARAMETERS_NAME).getParameter(#fallback,parameterName);}
//No throwing, just a warning
#define SETPARNAMEOPT(parameterName,fallback) try{SETPAR(parameterName);}catch(...){try{(PARAMETERS_NAME).getParameter(#fallback,parameterName);INFO("The parameter name '"#fallback"'is deprecated, use'"#parameterName"' in stead");}catch(...){ostringstream buf; buf << "Using default value of '"<<parameterName<<"' for '"<<#parameterName<<"'"; WARNING(buf.str());}}
//Try to read a parameter but use contents of string for the galdef name
//throw an exception if not found
#define SETPARSTR(parameterName,string) (PARAMETERS_NAME).getParameter(string,parameterName)
//No throwing, just a warning
#define SETPARSTROPT(parameterName,string) try{SETPARSTR(parameterName,string);}catch(...){ostringstream buf; buf << "Using default value of '"<<parameterName<<"' for '"<<string<<"'"; WARNING(buf.str());}
//Set a parameter with a different name, should probably not be used
//#define SETPARDIFF(parameterName,galdefName) (PARAMETERS_NAME).getParameter(#galdefName,parameterName)
//#define SETPARDIFFOPT(parameterName,galdefName) try{SETPARDIFF(parameterName,galdefName);}catch(...){ostringstream buf; buf << "Using default value of '"<<parameterName<<"' for '"<<#galdefName<<"'"; WARNING(buf.str());}

//Define a << operator for a vector
template <typename T>
std::ostream& operator<< (std::ostream& out, std::vector<T> &v ){
  std::copy(v.begin(),v.end(),std::ostream_iterator<T>(out," "));
   return out;
}

class Galdef {
  
  
 public:
  
  std::string version;    // galprop version number
  std::string run_no;    // identifier of galprop run
  std::string galdef_ID; // full identifier e.g. 01_123456
  
  int n_spatial_dimensions;             // 1,2 or 3
  double z_min,z_max,dz;                // for 1,2,3D    
  double r_min,r_max,dr;                // for 2D 
  double x_min,x_max,dx,y_min,y_max,dy; // for 3D 
  double p_min,p_max,p_factor;          // momentum start, end, factor
  double Ekin_min,Ekin_max,Ekin_factor; // kinetic energy per nucleon start, end, factor
  
  std::string p_Ekin_grid;                    // "p"||"Ekin": construct grid in p or Ekin
  
  double E_gamma_min,E_gamma_max,E_gamma_factor; // gamma-ray energy (MeV) start, end, factor
  double long_min,long_max;                      // gamma-ray skymap longitude min,max (degrees)
  double  lat_min, lat_max;                      // gamma-ray skymap latitude  min,max (degrees)
  double d_long,d_lat;                           // gamma-ray skymap longitude,latitude binsize (degrees)       
  int integration_mode;                          // integr.over the particle spectrum: =1-old E*logE; !=1-power-law analyt.
  int healpix_order;                             // The order for healpix maps.
  int lat_substep_number;                        // latitude bin splitting (0,1=no split, 2=split in 2...) IMOS20080114
  double LoS_step;                               // kpc, Line of Sight (LoS) integration step              IMOS20080114
  double LoS_minStep;                            // kpc, minimum Line of Sight (LoS) integration step when integration mode > 0
  int LoS_substep_number;                        // number of substeps per LoS integration step for gas number density averaging (0,1=no substeps) IMOS20080114
  int los_integration_mode;                      // 1 gives new los integration, != 1 gives old one
  double los_integration_accuracy;               // The relative accuracy for the LOS integration.  Put to 0 for fixed step size
  int anisoHealpixOrder;  // Healpix binning for calculating anisotropic IC

  double nu_synch_min,nu_synch_max,nu_synch_factor;// synchrotron frequency min,max (Hz), factor
  
  double D0_xx;                         // diffusion coefficient at reference rigidity
  double D_rigid_ref;                   // reference rigidity
  double D_rigid_br;                    // break     rigidity for diffusion coefficient in MV
  double D_g_1;                         // diffusion coefficient index below break  rigidity
  double D_g_2;                         // diffusion coefficient index above break  rigidity 
  double D_eta;                         // The power of beta (D_xx \propto beta^eta).  Defaults to 1.  Not applied if diff_reacc < 0

  double Dxx_plane_scale;               // Scale parameter for adjusting the diffusion coefficient in the plane
  double Dxx_plane_scale_height;        // The scale height for the diffusion coefficient adjustment in the plane, Gaussian.

  int B_dep_diffusion;                  // Rescale diffusion coefficient and Alfven velocity with B-field
  double D_xx_max;                      // Maximum value of the diffusion coefficient at 10 GV, used when linked with B-field
  double D_xx_min;                      // Minimum value of the diffusion coefficient at 10 GV, used when linked with B-field
  
  int diff_reacc;                       // 1,2=incl.diff.reacc.; 11=Kolmogorov+damping; 12=Kraichnan+damping
  double v_Alfven;                      // Alfven speed in km s-1
  double damping_p0;                    // ~ 1.e5 MV -characteristic rigidity         IMOS20030129
  double damping_max_path_L;            // ~ 1.*kpc2cm; Lmax~1 kpc, max free path     IMOS20030129
  double damping_const_G;               // a const derived from fitting B/C           IMOS20030129
  double damping_const_K;               // a const derived from fitting B/C, used in damping approximation (diff_reacc = 2)           
  
  int convection;                       // 1=include convection                       AWS20010323
  double   v0_conv;                     // v=v0_conv+dvdz_conv*abs(z)                 AWS20010323
  double dvdz_conv;                     // v=v0_conv+dvdz_conv*abs(z)                 AWS20010323
  
  double He_H_ratio;                    // He/H of ISM, by number
  int     n_X_CO;                       // option for X_CO values IMOS20080114
  double  X_CO;                         // conversion factor CO integrated temperature -> H2 column density IMOS20080114
  std::vector<double> X_CO_radius;           // for n_X_CO=1, radial values for the interpolation
  std::vector<double> X_CO_values;           // for n_X_CO=1, XCO values at corresponding radius
  std::vector<double> X_CO_parameters;       // for n_X_CO=2
  int    propagation_X_CO;              // option to control H2 from CO for propagation  AWS20090623
  int    nHI_model;                     // selection of model for HI  gas density        AWS20090814
  int    nH2_model;                     // selection of model for H2  gas density        AWS20090814
  int    nHII_model;                    // selection of model for HII gas density        AWS20090814
  std::string HI_xmlFilename;           // xml description for HI distribution if HI_model == 9
  std::string H2_xmlFilename;           // xml description for H2 distribution if H2_model == 9
  std::string HII_xmlFilename;          // xml description for HII distribution if HII_model == 9

  double  HII_Te;                       // free electron temperature                     AWS20110701
  double  HII_clumping_factor;          // free electron clumping factor                 AWS20110701

  std::string  COR_filename;                 // CO -molecular gas file IMOS20080114
  std::string  HIR_filename;                 // HI -atomic    gas file IMOS20080114

  std::vector<std::string> source_class_files;   // List of filenames pointing to source class files.
  std::vector<std::unique_ptr<SourceClass> > source_classes;  // Pointers to source classes.

  int fragmentation;                    // 1=include fragmentation
  int momentum_losses;                  // 1=include momentum losses
  int radioactive_decay;                // 1=include radioactive_decay
  int K_capture;                        // 1=include K_capture                        AWS20010731
  int ionization_rate;                  // 1=compute ionization rate          IMOS20060420

  bool ionization_losses;               // False=turn off ionization losses
  bool coulomb_losses;                  // False=turn off coulomb losses
  bool bremss_losses;                   // False=turn off bremss losses
  bool IC_losses;                       // False=turn off IC losses
  bool sync_losses;                     // False=turn off sync losses
  
  double start_timestep;                // start time step in years
  double   end_timestep;                //   end time step in years
  double       timestep_factor;         //   factor to multiply timestep
  int          timestep_repeat;         //   number of times to repeat for factor
  int          timestep_repeat2;        //   number of  times to repeat in timestep_mode=2
  int          timestep_print ;         //   number of timesteps between printings
  int          timestep_diagnostics;    //   number of timesteps between propel diagnostics
  int           control_diagnostics;    //   control details of propel diagnostics
 
  int  solution_convergence;            //   control use of convergence diagnostic in solution of propagation AWS20110118
  int  solution_method;                 //   method for solving propagation equation                          AWS20110118     
  double solution_rel_accuracy;         //   Stopping criteria (relative change) for iterations in CN solver
  
  int  network_iterations;              //   number of iterations of the protons (Needed for damping)
  int  network_iter_compl;              //   number of iterations for the entire network, 1 should be enough
                                        //   network_iterations >= network_iter_compl
  
  int prop_r;                           // for 2D: 1=propagate in r;
  int prop_x,prop_y,prop_z;             // for 2D: 1=propagate in z;for 3D 1=propagate in x,y,z
  int prop_p;                           // propagate in momentum
  int use_symmetry;                     // xyz symmetry (3D)
  
  int vectorized;                       // 0=unvectorized   code, 1=vectorized   code
  
  int  max_Z;                           // maximum   nucleus  Z in galdef file
  std::vector<int> use_Z;               // 1=use this nucleus Z
  
  int total_cross_section;              // total inel. cross-sections option AWS20010620
  int cross_section_option;             // controls which cross-sections used
  
  double t_half_limit;                  // lower limit on radioactive half-life for explicit inclusion  AWS20010731
  
  int primary_electrons;                // 1=propagate primary electrons
  int primary_positrons;                // 1=propagate primary positrons
  int secondary_positrons;              // 1=propagate secondary positrons
  int secondary_electrons;              // 1=propagate secondary electrons
  int knock_on_electrons;               // 1=propagate knock-on electrons     IMOS20060504
  int tertiary_antiprotons;             // 1=propagate tertiary antiprotons   IMOS20000605.13
  int secondary_antiprotons;            // 1=propagate secondary antiprotons
  int secondary_protons;                // 1=propagate secondary protons      IMOS20000605.14
  
  int gamma_rays;                       // 1=compute gamma-ray emission
  int pi0_decay;                        // 1 - Dermer 1986 formalism, 2 - Blattnig et al. 2000,PRD 62,094030  IMOS20050216
  int IC_isotropic;                     // 1=compute isotropic inverse Compton IMOS20060420
  int IC_anisotropic;                   // 1=compute anisotropic inverse Compton
  int bremss;                           // 1=compute bremsstrahlung            IMOS20060420
  int synchrotron;                      // 1=compute synchrotron emission
  int free_free_absorption;             // >=1 free-free absorption for synchrotron     AWS20110701

  int pair_production;                   // 1=compute secondary pairs from gamma absorption on ISRF

  int globalLuminosities;               // 1=compute global luminosities
  // DM parameters IMOS20050912
  int DM_positrons;                     // 1=compute positrons from DM
  int DM_electrons;                     // 1=compute electrons from DM
  int DM_antiprotons;                   // 1=compute antiprotons from DM
  int DM_gammas;                        // 1=compute gamma rays from DM
  double                                // user-defined params of DM (double)
    DM_double0, DM_double1, DM_double2, DM_double3, DM_double4,
    DM_double5, DM_double6, DM_double7, DM_double8, DM_double9;
  int                                   // user-defined params of DM (int)
    DM_int0, DM_int1, DM_int2, DM_int3, DM_int4,
    DM_int5, DM_int6, DM_int7, DM_int8, DM_int9;
  
  double local_bubble_radius;           // If >0 and 3D, specifies a radius around the sun affected by the local bubble
  double local_bubble_source_fraction;  // Multiply the source distribution with this number within the local bubble
  double local_bubble_gas_fraction;     // Multiply the total gas distribution with this number within the local bubble
  
//  int HI_survey;                        // HI survey : 8=original 8 rings+high-latitudes, 9= 9 rings all-sky  AWS20050913
//  int CO_survey;                        // CO survey : 8=original 8 rings                 9= 9 rings all-sky  AWS20050913
  
  int     B_field_model;                // >1000=parameterized model
  std::string  B_field_name;                 // 3D B-field model name ("galprop_original" uses B_field_model)      AWS20080313 
  std::vector<double> B_field_parameters;    // parameters for 3D B-field models                                   AWS20080313

  std::string ISRF_file;   // ISRF input file AWS20050301
  int ISRF_filetype;        // 0 for CMB, 1 for old (<Jan07), 2 for new FITS, 3 for new HEALPix
  int ISRF_healpixOrder; // Healpix binning for ISRF skymaps
  std::vector<double> ISRF_factors;  // ISRF factors for inverse Compton calculation       AWS20050301          
  
  double rigid_min,inj_Ekin_min;            // minimum rigidity/Ekin of injection spectra
  double rigid_max,inj_Ekin_max;            // maximum rigidity/Ekin of injection spectra
  
  double   proton_norm_Ekin;            // proton   kinetic energy for normalization (MeV)
  double   proton_norm_flux;            // flux of protons   at normalization energy (cm^-2 sr^-1 s^-1 MeV^-1)
  int      proton_norm_type;            // 0=no normalization, 1= to flux at sun, 2=to luminosity in particles s-1, 3= to luminosity in erg s-1   AWS20101202
  double electron_norm_Ekin;            // electron kinetic energy for normalization (MeV)
  double electron_norm_flux;            // flux of electrons at normalization energy (cm^-2 sr^-1 s^-1 MeV^-1)
  int    electron_norm_type;            // 0=no normalization, 1= to flux at sun, 2=to luminosity in particles s-1, 3= to luminosity in erg s-1   AWS20101202

  std::vector<double> fCameraLocation; // Location of camera for sky map generation
  int skymap_format;                    // Select the output format for the skymap
  
  int output_gcr_full;                  // output full 2D or 3D gcr array
  int warm_start;                       // read nuclei file and continue run   AWS20010121
  
  int verbose;                          // verbosity: 0=min 10=max 
  int test_suite;                       // run test suit instead of normal run 
  
  //interface functions prototypes
  Galdef();
  ~Galdef();
  int read(const std::string& version, 
	   const std::string& runNumber,
	   const std::string& galdefDirectory);
	   //char *version_,char  *run_no_, const std::string& galdef_directory);//char *galdef_directory);
  int AssignParameters( const utl::Parameters &pars );
  void AssignDefaultParameters();
  void print();

 private:
  std::string gdDirectory;
};

#endif










