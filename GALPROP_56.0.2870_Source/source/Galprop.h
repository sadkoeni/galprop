#ifndef GALPROP_H
#define GALPROP_H

#include "Distribution.h"
#include "Particle.h"
#include "Spectrum.h"
#include "Galaxy.h"
#include "Configure.h"
#include "Galdef.h"
#include "los_integration.h"
#include "DistributionFunction.h"

#include <string>
#include <valarray>
#include <vector>

//using namespace std;

class Galprop {

 public:
  
  Galprop();
  ~Galprop();
  
  int Run ( std::string galdefPath,
	    std::string fitsPath,
	    std::string outputPath,
	    std::string outputPrefix,
	    std::string runNumber );//int argc, char** argv);
  
  int AvXCO ( std::string galdefPath,
	      std::string fitsPath,
	      std::string outputPath,
	      std::string outputPrefix,
	      std::string runNumber );//int argc, char** argv);
  
  int create_gcr();
  int fill_transport_arrays ( Particle& );
  int create_galaxy();
  //int create_SNR();
  int cr_luminosity();
  int propagate_particles();
  int e_KN_loss ( Particle& );
  
  double D_pp ( double,double,double,double );
  double alfvenVelocity(int,int,int,int);
  int    D_xx ( Particle&,int,int,int,int,int,int );
  double Dxx_from_Bfield_model(double, int,int,int,int,int);
  static double fu ( double );
  
  // DM routines IMOS20050912
  int gen_DM_source ( Particle& );
  int gen_DM_emiss();
  double DM_profile ( double, double, double );
  double DM_profile_av ( double,double,double,double,double );
  double DM_profile_av ( double,double,double,double,double,double,double );
  int store_DM_emiss();
  int store_DM_skymap();
  
  double IC_anisotropy_factor ( double,double,double,double,int );//IMOS20060420
  void decayed_cross_sections ( int iz,int ia,int jz,int ja, double *Ekin,int np,double *sigma );
  
  int nuclei_normalize();
  int electrons_normalize();
  
  int propel ( Particle& );

  int propel_diagnostics (); //AWS20110113
  int propel_diagnostics ( Particle& particle,
			   Distribution& alpha1_r,
			   Distribution& alpha1_z,
			   Distribution& alpha1_p,
			   Distribution& alpha2_r,
			   Distribution& alpha2_z,
			   Distribution& alpha2_p,
			   Distribution& alpha3_r,
			   Distribution& alpha3_z,
			   Distribution& alpha3_p,
			   Distribution& total_source_function,
			   double dt );
  int propel_diagnostics ( Particle& particle,
			   Distribution& alpha1_x,
			   Distribution& alpha1_y,
			   Distribution& alpha1_z,
			   Distribution& alpha1_p,
			   Distribution& alpha2_x,
			   Distribution& alpha2_y,
			   Distribution& alpha2_z,
			   Distribution& alpha2_p,
			   Distribution& alpha3_x,
			   Distribution& alpha3_y,
			   Distribution& alpha3_z,
			   Distribution& alpha3_p,
			   Distribution& total_source_function,
			   double dt );
  
  /*
  void protri ( Particle& particle,
		Distribution& alpha1_x,
		Distribution& alpha1_y,
		Distribution& alpha1_z,
		Distribution& alpha1_p,
		Distribution& alpha2_x,
		Distribution& alpha2_y,
		Distribution& alpha2_z,
		Distribution& alpha2_p,
		Distribution& alpha3_x,
		Distribution& alpha3_y,
		Distribution& alpha3_z,
		Distribution& alpha3_p,
		Distribution& Nx1_,
		Distribution& Ny1_,
		Distribution& Nz1_,
		Distribution& Np1_,
		Distribution& Nx2_,
		Distribution& Ny2_,
		Distribution& Nz2_,
		Distribution& Np2_,
		Distribution& Nx3_,
		Distribution& Ny3_,
		Distribution& Nz3_,
		Distribution& Np3_,
		Distribution& total_source_function,
		double dt,
		int nrept_outer,
		double f_use );
  */
  
  //double source_distribution(const double x, const double y, const double z, int srcModel, const std::vector<double> &parameters, const std::vector<double> &source_values, const std::vector<double> &source_radius);
  //int source_SNR_event ( Particle &particle,double t );
  //int source_SNR_event_vec ( Particle &particle,double t,
//			     float *total_source_function_x );
  
  int test_Particle();
  //int test_source_SNR_event();
  int test_isotope_cs();
  int test_cfactor();
  
  int print_BC();
  float isrf_energy_density ( float rr, float zz );
  int HIR_iRing ( double RR );
  
  int read_gcr();
  int read_isrf ( const int version );
  
  int read_gas_maps (const std::string &type); //IMOS20080114
  int gas_iRing ( double ); //IMOS20080114
  double fX_CO ( double ) const;  //IMOS20080114
  
  //Classes for use with LOSintegrators
  class GasFunction : public SM::LOSfunction<double> {
  private:
    GasFunction();
    enum TYPE { HI, H2, CO, HII } ftype;
    const double frInd;
    const Galprop &fgp; //To have access to fX_CO
  public:
    // The functions can be optionally multiplied with radius to the index
    // rPLindex
    // Type is one from the type enum above (represented as a string)
    GasFunction(const std::string& type, double rPLindex, const Galprop& gp);
    virtual double operator () ( const double x, const double y, const double z, const vec3 &dir ) const;
  };

  class GasEmissFunction : public SM::LOSfunction<std::valarray<double> > {
  private:
    GasEmissFunction();
    Distribution test;
    enum TYPE { BREMSS, PION, TEST } ftype;
    const GasFunction fgf;
    DistributionFunction *fdf;
    const Galprop &fgp;
  public:
    //Type is one of the above, TEST has unit emissivity.  For gas_type see GasFunction
    GasEmissFunction(const std::string &type, const std::string &gas_type, const Galprop& gp);
    ~GasEmissFunction();
    virtual std::valarray<double> operator () ( const double x, const double y, const double z, const vec3 &dir ) const;
  };

  class AnisoICFunction : public SM::LOSfunction< std::valarray<double> > {
  private:
    AnisoICFunction();
    std::valarray<double> isoCrossSection, anisoCrossSection;//, cosZetaTest;
    const double cameraX, cameraY, cameraZ;
    std::vector< Skymap<double> > optAngDist, irAngDist;
    std::valarray<double> targetE, gammaE, electronE, cmbNumberDensity, rGrid, xGrid, yGrid;
    std::valarray< std::valarray<double> > optISRF, irISRF;
    Skymap<double> cmbSkymap;
    Distribution electrons;
    //Healpix_Base hpTest;
    const Galprop &fgp;
    const size_t cosThetaBins;
    size_t rBins, nEGammaBins, targetBins, electronBins;
    double rMax;
    std::valarray<vec3> dirTarget;
  public:
    AnisoICFunction(const Galprop &gp);
    //Returns the emissivity for iso and aniso IC for x,y,z and dir.
    //Outputs emissivity arrays, iso opt, iso ir, iso cmb, then aniso in same order
    //The arrays are joined in the output so iso ir starts at index 1*emiss.size() and so forth
    virtual std::valarray<double> operator () ( const double x, const double y, const double z, const vec3 &dir ) const;

    const std::valarray<double> & get_isoCrossSection() const { return isoCrossSection; }
    const std::valarray<double> & get_anisoCrossSection() const { return anisoCrossSection; }

    const std::vector< Skymap<double> > & get_optAngDist() const { return optAngDist; }
    const std::vector< Skymap<double> > & get_irAngDist() const { return irAngDist; }

    const std::valarray< std::valarray<double> > & get_optISRF() const { return optISRF; }
    const std::valarray< std::valarray<double> > & get_irISRF() const { return irISRF; }

    const std::valarray<double> & get_targetE() const { return targetE; }
    const std::valarray<double> & get_gammaE() const { return gammaE; }
    const std::valarray<double> & get_electronE() const { return electronE; }
    const std::valarray<double> & get_rGrid() const { return rGrid; }
    const std::valarray<double> & get_cmbNumberDensity() const { return cmbNumberDensity; }

    const Skymap<double> & get_cmbSkymap() const { return cmbSkymap; }

    const Distribution & get_electrons() const { return electrons; }

    int get_cosThetaBins() const { return cosThetaBins; }
    unsigned int get_rBins() const { return rBins; }
    double get_rMax() const { return rMax; }
  };
    

  int store_gcr();
  int store_gcr_full();
  int store_gcr_source_functions ( Particle &particle );//AWS2010031
  
  void store_IC_emiss(); // TAP20090312
  int store_IC_skymap ( const std::string& type );  //AWS20090415
  int store_IC_skymap_comp ( const std::string& type );
  int store_bremss_emiss();
  int store_bremss_ionized_skymap();
  int store_bremss_skymap();
  int store_pi0_decay_emiss();
  int store_pi0_decay_skymap();
  int store_pi0_decay_H2R_skymap(); //AWS20041215
  int store_pi0_decay_HIR_skymap(); //AWS20041215
  int store_pi0_decay_HII_skymap(); //IMOS20080114*
  int store_bremss_H2R_skymap();    //AWS20041215
  int store_bremss_HIR_skymap();    //AWS20041215
  int store_bremss_HII_skymap();    //IMOS20080114*
  int store_synch_emiss();          //AWS20080314
  int store_synch_skymap();
  int store_ionization_rate();
  int store_skymap ( float array[], long naxes[4], const std::string name, double crval[4], double cdelt[4] );
  int store_mapcube_skymap ( float array[],
			     double energy[],
			     const int nComponents,
			     const int nEnergies,
			     const std::string& name,
			     const bool MeV ) const;
  
  int gen_secondary_source ( Particle& );
  int gen_isrf_energy_density();
  bool gen_IC_emiss(); // TAP20090311

  void gen_skymaps();
  int gen_bremss_emiss();
  int gen_pi0_decay_emiss();
  int gen_synch_emiss();
  int gen_secondary_positron_source ( Particle &particle );
  void gen_secondary_pair_source(Particle& particle);
  int gen_tertiary_antiproton_source ( Particle &particle );
  int gen_secondary_antiproton_source ( Particle &particle );
  int gen_secondary_proton_source ( Particle &particle );
  int gen_knock_on_electron_source ( Particle &particle ); //IMOS20060504
  double knock_on_cross_section ( double,double,int );     //IMOS20060504
  int gen_ionization_rate();
  int gen_luminosity(); //AWS20100121
  void store_luminosity();
  int test_suite ();
  
  //////////////////////////////////////////////////
  //Spectrum S;
  //Distribution Dist;
  //Particle P;
  
  int n_species;
  int isrf_energy_density_i_comp; // required for cfactor
  
  ///////////////////////////////////
  gp::Configure configure;
  Galdef galdef;
  
  Particle* gcr; // all species
  Galaxy galaxy;

};
#endif
