
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Particle.h *                                  galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef _Particle_h_
#define _Particle_h_

#include <iostream>   //AWS20050624
#include <cmath>      //AWS20050624
#include <string>
#include <valarray>
#include <map>

#include "Distribution.h"

using namespace std;
  
class Particle {

 public: 
 
  //char name[100];
  std::string name;
  int Z,A;
  int K_electron;           // number of K-electrons (0 or 1)      AWS20010731
  double mass;              // in MeV
  double t_half;            // half-life in years
  double primary_abundance; // primary isotopic abundance
  
  double z_min, z_max, dz;                // for 1,2,3D    
  double r_min, r_max, dr;                // for 2D 
  double x_min, x_max, dx, y_min, y_max, dy; // for 3D 
  
  double p_min, p_max, p_factor; // momentum start, end, factor
  double Ekin_min, Ekin_max, Ekin_factor; // kinetic energy/nucleon start,end,factor
  
  int n_spatial_dimensions;// 1,2,3D
  int n_pgrid;             // number of points in momentum
  int n_rgrid;             // number of points in radius (2D)   
  int n_zgrid;             // number of points in z (1D,2D)  
  int n_xgrid;             // number of points in x (3D)
  int n_ygrid;             // number of points in y (3D)    
  
  std::valarray<double> x; // x grid
  std::valarray<double> y; // y grid
  std::valarray<double> z; // z grid 
  std::valarray<double> r; // r grid 
  
  std::valarray<double> p; // total momentum of particle in MV
  std::valarray<double> Etot; // total energy   of particle in MeV
  std::valarray<double> Ekin; // kinetic energy per nucleon in MeV
  std::valarray<double> beta; // velocity/c
  std::valarray<double> gamma; // Etot/mass
  std::valarray<double> rigidity; // rigidity p/Z in MV  
  
  std::string species;
  
  bool arrays_assigned;
  
  double normalization_factor;// to normalize protons or electrons AWS20010121

  std::map<Particle*, std::valarray<double> > dependencies; //Store pointers to particles this one depends on and its cross section
  
  Distribution cr_density;
  Distribution primary_source_function;
  Distribution secondary_source_function;
  
  Distribution fragment;// fragmentation destruction rate
  Distribution decay   ;// radioactive decay rate        
  Distribution dpdt;    // momentum change rate       
  
  Distribution Dxx ;    // spatial  diffusion coefficient
  Distribution Dyy ;    // spatial  diffusion coefficient AWS20090625
  Distribution Dzz ;    // spatial  diffusion coefficient AWS20090625
  
  Distribution Dpp ;    // momentum diffusion coefficient
  Distribution v_conv;  // convection velocity in z-direction
  
  //interface functions prototypes
  Particle();
  Particle(const Particle& old);
  ~Particle();

  void init(const std::string &name_, 
	    int Z_, int A_, 
	    int K_electron_, 
	    double t_half_, //AWS20010731
	    double r_min_, double r_max_, double dr_, 
	    double z_min_, double z_max_, double dz_,
	    double p_min_, double p_max_, double p_factor_,
	    double Ekin_min_, double Ekin_max_, double Ekin_factor_,
	    const std::string& p_Ekin_grid);

  void init(const std::string &name_,
	    int Z_, int A_,
	    int K_electron_, 
	    double t_half_, //AWS20010731
	    double x_min_, double x_max_, double dx_, 
	    double y_min_, double y_max_, double dy_, 
	    double z_min_, double z_max_, double dz_,
	    double p_min_, double p_max_, double p_factor_,
	    double Ekin_min_, double Ekin_max_, double Ekin_factor_,
	    const std::string& p_Ekin_grid);
  int delete_transport_arrays();
  int create_transport_arrays();
  int delete_arrays();  //Gulli20070810
  void print();
  Particle& operator=(const Particle& particle);
};
  
#endif
