
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Particle.cc *                                 galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include"Particle.h"
#include"galprop_internal.h"
#include <cstring>
#include "constants.h"

//Mass of elements in units of atomic mass.  Table retreived from
//http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
//Uses information on mass number (A) only because it is identical up to
//four significant digits in all cases with multiple Z values but same A.
static vector<double> massNumber;
const static double mu = 931.494; //Atomic mass unit in MeV

static void setMassNumber() {
   if (massNumber.size() == 0) {
      massNumber.resize(100);

      //Electrons and positrons
      massNumber[0] = 5.485799e-4;
      massNumber[1] = 1.007825;
      massNumber[2] = 2.014101;
      massNumber[3] = 3.016;
      massNumber[4] = 4.0026;
      massNumber[5] = 5; //No data available
      massNumber[6] = 6.015;
      massNumber[7] = 7.016;
      massNumber[8] = 8; //No data available
      massNumber[9] = 9.012;
      massNumber[10] = 10.012;
      massNumber[11] = 11.009;
      massNumber[12] = 12.0;  //Definition of u
      massNumber[13] = 13.003;
      massNumber[14] = 14.003;
      massNumber[15] = 15.0;
      massNumber[16] = 15.995;
      massNumber[17] = 16.999;
      massNumber[18] = 17.999;
      massNumber[19] = 18.998;
      massNumber[20] = 19.992;
      massNumber[21] = 20.994;
      massNumber[22] = 21.991;
      massNumber[23] = 22.990;
      massNumber[24] = 23.985;
      massNumber[25] = 24.986;
      massNumber[26] = 25.982;
      massNumber[27] = 26.982;
      massNumber[28] = 27.977;
      massNumber[29] = 28.976;
      massNumber[30] = 29.973;
      massNumber[31] = 30.974;
      massNumber[32] = 31.972;
      massNumber[33] = 32.971;
      massNumber[34] = 33.968;
      massNumber[35] = 34.969;
      massNumber[36] = 35.967;
      massNumber[37] = 36.965;
      massNumber[38] = 37.963;
      massNumber[39] = 38.963;
      massNumber[40] = 39.963;
      massNumber[41] = 40.962;
      massNumber[42] = 41.959;
      massNumber[43] = 42.959;
      massNumber[44] = 43.956;
      massNumber[45] = 44.956;
      massNumber[46] = 45.953;
      massNumber[47] = 46.952;
      massNumber[48] = 47.950;
      massNumber[49] = 48.948;
      massNumber[50] = 49.946;
      massNumber[51] = 50.944;
      massNumber[52] = 51.941;
      massNumber[53] = 52.941;
      massNumber[54] = 53.939;
      massNumber[55] = 54.938;
      massNumber[56] = 55.935;
      massNumber[57] = 56.935;
      massNumber[58] = 57.934;
      massNumber[59] = 58.933;
      massNumber[60] = 59.931;
      massNumber[61] = 60.931;
      massNumber[62] = 61.928;
      massNumber[63] = 62.930;
      massNumber[64] = 63.929;
      massNumber[65] = 64.928;
      massNumber[66] = 65.926;
      massNumber[67] = 66.927;
      massNumber[68] = 67.925;
      massNumber[69] = 68.926;
      massNumber[70] = 69.924;
      massNumber[71] = 70.925;
      massNumber[72] = 71.922;
      massNumber[73] = 72.923;
      massNumber[74] = 73.922;
      massNumber[75] = 74.922;
      massNumber[76] = 75.920;
      massNumber[77] = 76.920;
      massNumber[78] = 77.918;
      massNumber[79] = 78.918;
      massNumber[80] = 79.916;
      massNumber[81] = 80.916;
      massNumber[82] = 81.915;
      massNumber[83] = 82.914;
      massNumber[84] = 83.912;
      massNumber[85] = 84.912;
      massNumber[86] = 85.910;
      massNumber[87] = 86.909;
      massNumber[88] = 87.906;
      massNumber[89] = 88.906;
      massNumber[90] = 89.904;
      massNumber[91] = 90.906;
      massNumber[92] = 91.906;
      massNumber[93] = 92.906;
      massNumber[94] = 93.906;
      massNumber[95] = 94.906;
      massNumber[96] = 95.906;
      massNumber[97] = 96.906;
      massNumber[98] = 97.906;
      massNumber[99] = 98.906;
   }
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
void Particle::init(const std::string &name_,int Z_,int A_,int K_electron_,double t_half_,  //AWS20010731
                    double r_min_, double r_max_, double dr_, 
                    double z_min_, double z_max_, double dz_,
                    double    p_min_,double    p_max_,double    p_factor_,
                    double Ekin_min_,double Ekin_max_,double Ekin_factor_,
                    const std::string& p_Ekin_grid) {
  
  int ir,iz,ip;
  
  delete_arrays();
  //cout<<" Particle: initializing 2D\n";
  n_spatial_dimensions=2;     // 2D    
  name = name_;
  //strcpy(name,name_);
  Z=Z_;
  A=A_; 
  mass=massNumber[A]*mu;
  K_electron=K_electron_;                                                     //AWS20010731                        
  t_half=t_half_;
  r_min=r_min_;
  r_max=r_max_;
  dr   =   dr_;
  z_min=z_min_;
  z_max=z_max_;
  dz   =   dz_;
  
  p_min=p_min_;
  p_max=p_max_;
  p_factor=p_factor_;
  
  Ekin_min=Ekin_min_;
  Ekin_max=Ekin_max_;
  Ekin_factor=Ekin_factor_;
  
  n_rgrid=(int) ((r_max-r_min)/dr + 1.5);
  n_zgrid=(int) ((z_max-z_min)/dz + 1.5);
  
  if(p_Ekin_grid == "p")
    n_pgrid=(int) (log(p_max/p_min)/log(p_factor) + 1.9);
  
  if(p_Ekin_grid == "Ekin")
    n_pgrid=(int) (log(Ekin_max/Ekin_min)/log(Ekin_factor) + 1.9);

  r.resize(n_rgrid);
  z.resize(n_zgrid);
  p.resize(n_pgrid);
  Ekin.resize(n_pgrid);
  Etot.resize(n_pgrid);
  beta.resize(n_pgrid);
  gamma.resize(n_pgrid);
  rigidity.resize(n_pgrid);
  
  //r=    new double[n_rgrid];
  //z=    new double[n_zgrid];
  //p=    new double[n_pgrid];
  //Ekin= new double[n_pgrid];
  //Etot= new double[n_pgrid];
  //beta= new double[n_pgrid];
  //gamma=new double[n_pgrid];
  //rigidity=new double[n_pgrid];
  for(ir=0; ir<n_rgrid; ir++) r[ir]=r_min+ir*dr; 
  for(iz=0; iz<n_zgrid; iz++) z[iz]=z_min+iz*dz; 
  
  species = "nucleus";
  if(A==0) species = "electron";
  
  if(p_Ekin_grid == "p")
    for(ip=0; ip<n_pgrid; ip++)
      {
	p[ip]   =exp(log(   p_min)+ip*log(   p_factor));
	Ekin[ip]=-1; 
	kinematic(Z,A,mass,p[ip],Ekin[ip],Etot[ip],beta[ip],gamma[ip],rigidity[ip],0);
      }
  
  if(p_Ekin_grid == "Ekin")
    for(ip=0; ip<n_pgrid; ip++)
      {
	Ekin[ip]=exp(log(Ekin_min)+ip*log(Ekin_factor));
	p[ip]=-1;
	kinematic(Z,A,mass,p[ip],Ekin[ip],Etot[ip],beta[ip],gamma[ip],rigidity[ip],0);
      }
  
  cr_density.init(n_rgrid, n_zgrid, n_pgrid);
  //   cout<<"cr"<<cr_density.d2[0][0].s[0]<<endl;
  
  arrays_assigned = true;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Particle::init(const std::string &name_,int Z_,int A_,int K_electron_,double t_half_,  //AWS20010731   
                    double x_min_, double x_max_, double dx_, 
                    double y_min_, double y_max_, double dy_, 
                    double z_min_, double z_max_, double dz_,
                    double p_min_, double p_max_, double p_factor_,
                    double Ekin_min_,double Ekin_max_,double Ekin_factor_,
                    const std::string& p_Ekin_grid) { 

  int ix,iy,iz,ip;    
  
  delete_arrays();
  //cout<<" Particle: initializing 3D\n";
  n_spatial_dimensions=3;     // 3D    
  name = name_;
  //strcpy(name,name_);
  Z=Z_;
  A=A_;
  mass=massNumber[A]*mu;
  K_electron=K_electron_;                                                     //AWS20010731
  t_half=t_half_;
  x_min=x_min_;
  x_max=x_max_;
  dx   =   dx_;
  y_min=y_min_;
  y_max=y_max_;
  dy   =   dy_;
  z_min=z_min_;
  z_max=z_max_;
  dz   =   dz_;
  
  p_min=p_min_;
  p_max=p_max_;
  p_factor=p_factor_;
  
  Ekin_min=Ekin_min_;
  Ekin_max=Ekin_max_;
  Ekin_factor=Ekin_factor_;
  
  n_xgrid=(int)((x_max-x_min)/dx + 1.5);
  n_ygrid=(int)((y_max-y_min)/dy + 1.5);
  n_zgrid=(int)((z_max-z_min)/dz + 1.5);
  if(p_Ekin_grid == "p")
    n_pgrid=(int) (log(p_max/p_min)/log(p_factor) + 1.9);
  
  if(p_Ekin_grid == "Ekin")
    n_pgrid=(int) (log(Ekin_max/Ekin_min)/log(Ekin_factor) + 1.9);
  //cout << "n_pgrid=" << n_pgrid << endl;
  
  x.resize(n_xgrid);
  y.resize(n_ygrid);
  z.resize(n_zgrid);
  p.resize(n_pgrid);
  Ekin.resize(n_pgrid);
  Etot.resize(n_pgrid);
  beta.resize(n_pgrid);
  gamma.resize(n_pgrid);
  rigidity.resize(n_pgrid);

  //x=    new double[n_xgrid];
  //y=    new double[n_ygrid];
  //z=    new double[n_zgrid];
  //p=    new double[n_pgrid];
  //Ekin =new double[n_pgrid];
  //Etot =new double[n_pgrid];
  //beta =new double[n_pgrid];
  //gamma=new double[n_pgrid];
  //rigidity=new double[n_pgrid];
  
  for(ix=0; ix<n_xgrid; ix++) x[ix]=x_min+ix*dx; 
  for(iy=0; iy<n_ygrid; iy++) y[iy]=y_min+iy*dy; 
  for(iz=0; iz<n_zgrid; iz++) z[iz]=z_min+iz*dz; 
  
  species = "nucleus";
  if(A==0) species = "electron";
  
  if(p_Ekin_grid == "p")
    for(ip=0; ip<n_pgrid; ip++)
      {
	p[ip]   =exp(log(   p_min)+ip*log(   p_factor));
         Ekin[ip]=-1; 
         kinematic(Z,A,mass,p[ip],Ekin[ip],Etot[ip],beta[ip],gamma[ip],rigidity[ip],0);
      }
  
  if(p_Ekin_grid == "Ekin")
    for(ip=0; ip<n_pgrid; ip++)
      {
	Ekin[ip]=exp(log(Ekin_min)+ip*log(Ekin_factor));
	p[ip]=-1;
	kinematic(Z,A,mass,p[ip],Ekin[ip],Etot[ip],beta[ip],gamma[ip],rigidity[ip],0);
      }
  
  cr_density.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
  
  arrays_assigned = true;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// flags object as having arrays not yet assigned

Particle::Particle() : 
  name(""), 
  Z(0), A(0), 
  K_electron(0), 
  mass(0), 
  t_half(0), 
  primary_abundance(0), 
  z_min(0), z_max(0), dz(0),
  r_min(0), r_max(0), dr(0), 
  x_min(0), x_max(0), dx(0),
  y_min(0), y_max(0), dy(0), 
  p_min(0), p_max(0), p_factor(0),
  Ekin_min(0), Ekin_max(0), Ekin_factor(0), 
  n_spatial_dimensions(-1), 
  n_pgrid(0), n_rgrid(0), n_zgrid(0), n_xgrid(0), n_ygrid(0),
  arrays_assigned(false), 
  normalization_factor(1) {

  //Initalize the mass numbers
  setMassNumber();
  
}

Particle::~Particle(){

  delete_arrays();

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

//Gulli20070810
int Particle::delete_arrays() {

  if (arrays_assigned){
    delete_transport_arrays();
    cr_density.delete_array();
    //if (n_spatial_dimensions == 3){
    //delete[] x;
    //delete[] y;
    //}
    //if (n_spatial_dimensions == 2){
    //delete[] r;
    //}
    //delete[] z;
    //delete[] p;
    //delete[] Ekin;
    //delete[] Etot;
    //delete[] beta;
    //delete[] gamma;
    //delete[] rigidity;
    arrays_assigned = false;
  }
  return 0;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Particle::delete_transport_arrays() {

  //cout<<"Particle delete_transport_arrays "<<name<<endl;
  
  primary_source_function.delete_array();
  secondary_source_function.delete_array();
  fragment.delete_array();
  //if(t_half!=0.0)      //Gulli20070810
  decay.delete_array();
  dpdt.delete_array();
  Dxx.delete_array();
  Dpp.delete_array();
  v_conv.delete_array();
  return 0;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Particle::create_transport_arrays() {
  
  if (2 == n_spatial_dimensions) {
    //cr_density                 .init(n_rgrid, n_zgrid, n_pgrid); //Gulli20070810
    primary_source_function.init(n_rgrid, n_zgrid, n_pgrid);
    secondary_source_function.init(n_rgrid, n_zgrid, n_pgrid);
    fragment.init(n_rgrid, n_zgrid, n_pgrid);
    decay.init(n_rgrid, n_zgrid, n_pgrid);
    dpdt.init(n_rgrid, n_zgrid, n_pgrid);
    Dxx.init(n_rgrid, n_zgrid, n_pgrid);
    Dpp.init(n_rgrid, n_zgrid, n_pgrid);
    v_conv.init(n_rgrid, n_zgrid, n_pgrid);
  
  } // 2D
  
  if (3 == n_spatial_dimensions) {

    //cr_density                 .init(n_xgrid, n_ygrid, n_zgrid, n_pgrid); //Gulli20070810
    primary_source_function.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    secondary_source_function.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    fragment.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    decay.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    dpdt.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    Dxx.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    Dpp.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    v_conv.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);

  } // 3D
  
  return 0;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Particle::print() {

  cout<<"name= " <<name <<endl;
  cout<<"species = " <<species <<endl;
  cout<<"Z   = " <<Z    <<endl;
  cout<<"A   = " <<A    <<endl;
  cout<<"K_electron = " <<K_electron    <<endl;                               //AWS20010731       
  cout<<"t_half="<<t_half<<endl;
  cout<<"primary abundance    ="<< primary_abundance   <<endl;
  cout<<"n_spatial_dimensions ="<< n_spatial_dimensions<<endl;
  
  if(n_spatial_dimensions==2)
    {
      cout<<"r_min="<<r_min<<endl;
      cout<<"r_max="<<r_max<<endl;
      cout<<"dr   ="<<dr   <<endl;
      cout<<"z_min="<<z_min<<endl;
      cout<<"z_max="<<z_max<<endl;
      cout<<"dz   ="<<dz   <<endl;
      cout<<"p_min="<<p_min<<endl;
      cout<<"p_max="<<p_max<<endl;
      cout<<"p_factor="<<p_factor<<endl;
      cout<<"Ekin_min="<<Ekin_min<<endl;
      cout<<"Ekin_max="<<Ekin_max<<endl;
      cout<<"Ekin_factor="<<Ekin_factor<<endl;

      cout<<"n_rgrid="<<n_rgrid<<endl;
      cout<<"n_zgrid="<<n_zgrid<<endl;
      cout<<"n_pgrid="<<n_pgrid<<endl;

      cout<<"r grid:"<<endl;
      for(int ir=0; ir<n_rgrid; ir++) cout<<r[ir]<<" ";  cout<<endl;
      cout<<"z grid:"<<endl;
      for(int iz=0; iz<n_zgrid; iz++) cout<<z[iz]<<" ";  cout<<endl;
      cout<<"p grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<p[ip]<<" ";  cout<<endl;
      cout<<"Etot grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<Etot[ip]<<" ";  cout<<endl;
      cout<<"Ekin grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<Ekin[ip]<<" ";  cout<<endl;
      cout<<"beta grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<beta[ip]<<" ";  cout<<endl;
      cout<<"gamma grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<gamma[ip]<<" ";  cout<<endl;
      cout<<"rigidity grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<rigidity[ip]<<" ";  cout<<endl;

      cout<<"cr_density[0][0].s[0]="<<cr_density.d2[0][0].s[0]<<endl;
//cout<<"primary_source_function[0][0].s[0]="<<primary_source_function.d2[0][0].s[0]<<endl;
//cout<<"secondary_source_function[0][0].s[0]="<<secondary_source_function.d2[0][0].s[0]<<endl;

   } // (n_spatial_dimensions==2

   if(n_spatial_dimensions==3)
   {
      cout<<"x_min="<<x_min<<endl;
      cout<<"x_max="<<x_max<<endl;
      cout<<"dx   ="<<dx   <<endl;
      cout<<"y_min="<<y_min<<endl;
      cout<<"y_max="<<y_max<<endl;
      cout<<"dy   ="<<dy   <<endl;
      cout<<"z_min="<<z_min<<endl;
      cout<<"z_max="<<z_max<<endl;
      cout<<"dz   ="<<dz   <<endl;
      cout<<"p_min="<<p_min<<endl;
      cout<<"p_max="<<p_max<<endl;
      cout<<"p_factor="<<p_factor<<endl;
      cout<<"Ekin_min="<<Ekin_min<<endl;
      cout<<"Ekin_max="<<Ekin_max<<endl;
      cout<<"Ekin_factor="<<Ekin_factor<<endl;

      cout<<"n_xgrid="<<n_xgrid<<endl;
      cout<<"n_ygrid="<<n_ygrid<<endl;
      cout<<"n_zgrid="<<n_zgrid<<endl;
      cout<<"n_pgrid="<<n_pgrid<<endl;

      cout<<"x grid:"<<endl;
      for(int ix=0; ix<n_xgrid; ix++) cout<<x[ix]<<" ";  cout<<endl;
      cout<<"y grid:"<<endl;
      for(int iy=0; iy<n_ygrid; iy++) cout<<y[iy]<<" ";  cout<<endl;
      cout<<"z grid:"<<endl;
      for(int iz=0; iz<n_zgrid; iz++) cout<<z[iz]<<" ";  cout<<endl;
      cout<<"p grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<p[ip]<<" ";  cout<<endl;
      cout<<"Etot grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<Etot[ip]<<" ";  cout<<endl;
      cout<<"Ekin grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<Ekin[ip]<<" ";  cout<<endl;
      cout<<"beta grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<beta[ip]<<" ";  cout<<endl;
      cout<<"gamma grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<gamma[ip]<<" ";  cout<<endl;
      cout<<"rigidity grid:"<<endl;
      for(int ip=0; ip<n_pgrid; ip++) cout<<rigidity[ip]<<" ";  cout<<endl;

      cout<<"cr_density[0][0][0].s[0]="<<cr_density.d3[0][0][0].s[0]<<endl;
//cout<<"primary_source_function[0][0][0].s[0]="<<primary_source_function.d3[0][0][0].s[0]<<endl;
//cout<<"secondary_source_function[0][0][0].s[0]="<<secondary_source_function.d3[0][0][0].s[0]<<endl;

   } // (n_spatial_dimensions==3
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Particle::Particle(const Particle& particle) {

  name = particle.name;
  //   strcpy( name,particle.name);
   
  Z=particle.Z;
  A=particle.A;
  K_electron=particle.K_electron;                                             //AWS20010731
  mass=particle.mass ;
  t_half=particle.t_half;        
  primary_abundance=particle.primary_abundance; 
  
  species = particle.species;

  dependencies = particle.dependencies;
  
  n_spatial_dimensions=particle.n_spatial_dimensions;
  n_pgrid=particle.n_pgrid;             // number of points in momentum
  n_zgrid=particle.n_zgrid;             // number of points in z (1D,2D)  
  // number of points in radius (2D)   
  if(n_spatial_dimensions==2) n_rgrid=particle.n_rgrid; 
  if(n_spatial_dimensions==3)
    {
      n_xgrid=particle.n_xgrid;          // number of points in x (3D)
      n_ygrid=particle.n_ygrid;          // number of points in y (3D)    
    }
  
  z_min=particle.z_min;
  z_max=particle.z_max;
  dz=particle.dz;                   // for 1,2,3D    
  
  r_min=particle.r_min;
  r_max=particle.r_max;
  dr=particle.dr;                // for 2D 
  
  x_min=particle.x_min;
  x_max=particle.x_max;
  dx=particle.dx;
  y_min=particle.y_min;
  y_max=particle.y_max;
  dy=particle.dy;                // for 3D 
  
  p_min   =particle.p_min;  
  p_max   =particle.p_max;
  p_factor   =particle.p_factor;
  Ekin_min   =particle.Ekin_min; 
  Ekin_max   =particle.Ekin_max;
  Ekin_factor=particle.Ekin_factor;
  
  normalization_factor = particle.normalization_factor;
  
  p.resize(n_pgrid);
  p = particle.p;

  Etot.resize(n_pgrid);
  Etot = particle.Etot;

  Ekin.resize(n_pgrid);
  Ekin = particle.Ekin;

  beta.resize(n_pgrid);
  beta = particle.beta;

  gamma.resize(n_pgrid);
  gamma = particle.gamma;

  rigidity.resize(n_pgrid);
  rigidity = particle.rigidity;

  //p=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) p[ip]=particle.p[ip];
  //Etot=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) Etot[ip]=particle.Etot[ip];
  //Ekin=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) Ekin[ip]=particle.Ekin[ip];
  //beta=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) beta[ip]=particle.beta[ip];
  //gamma=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) gamma[ip]=particle.gamma[ip];
  //rigidity=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) rigidity[ip]=particle.rigidity[ip];
  
  z.resize(n_zgrid);
  z = particle.z;

  r.resize(n_rgrid);
  r = particle.r;

  x.resize(n_xgrid);
  x = particle.x;

  y.resize(n_ygrid);
  y = particle.y;

  //z=new double[n_zgrid];
  //for(int iz=0; iz<n_zgrid; iz++) z[iz]=particle.z[iz];
  
  //if(n_spatial_dimensions==2)
  //{
  //  r=new double[n_rgrid];
  //  for(int ir=0; ir<n_rgrid; ir++) r[ir]=particle.r[ir];
  //}
  
  //if(n_spatial_dimensions==3)
  //{
  //  x=new double[n_xgrid];
  //  for(int ix=0; ix<n_xgrid; ix++) x[ix]=particle.x[ix];
  //  y=new double[n_ygrid];
  //  for(int iy=0; iy<n_ygrid; iy++) y[iy]=particle.y[iy];
  //}
  
  arrays_assigned=1;
  cr_density =particle. cr_density;
  primary_source_function = particle.primary_source_function;
  secondary_source_function = particle.secondary_source_function;
  fragment = particle.fragment;
  decay    = particle.decay;
  dpdt     = particle.dpdt;
  Dxx      = particle.Dxx;
  Dpp      = particle.Dpp;
  v_conv   = particle.v_conv;
}

Particle& Particle::operator=(const Particle& particle) {
  
  //strcpy( name,particle.name);
  name = particle.name;
  
  Z=particle.Z;
  A=particle.A;
  K_electron=particle.K_electron;                                             //AWS20010731
  mass=particle.mass ;
  t_half=particle.t_half;        
  primary_abundance=particle.primary_abundance; 
  
  species = particle.species;

  dependencies = particle.dependencies;
  
  n_spatial_dimensions=particle.n_spatial_dimensions;
  n_pgrid=particle.n_pgrid;             // number of points in momentum
  n_zgrid=particle.n_zgrid;             // number of points in z (1D,2D)  
  // number of points in radius (2D)   
  if(n_spatial_dimensions==2) n_rgrid=particle.n_rgrid; 
  if(n_spatial_dimensions==3)
    {
      n_xgrid=particle.n_xgrid;          // number of points in x (3D)
      n_ygrid=particle.n_ygrid;          // number of points in y (3D)    
    }
  
  z_min=particle.z_min;
  z_max=particle.z_max;
  dz=particle.dz;                   // for 1,2,3D    
  
  if(n_spatial_dimensions==2)
    {
      r_min=particle.r_min;
      r_max=particle.r_max;
      dr=particle.dr;                // for 2D 
    }
  
  if(n_spatial_dimensions==3)
    {
      x_min=particle.x_min;
      x_max=particle.x_max;
      dx=particle.dx;
      y_min=particle.y_min;
      y_max=particle.y_max;
      dy=particle.dy;                // for 3D 
    }
  
  p_min   =particle.p_min;  
  p_max   =particle.p_max;
  p_factor   =particle.p_factor;
  Ekin_min   =particle.Ekin_min; 
  Ekin_max   =particle.Ekin_max;
  Ekin_factor=particle.Ekin_factor;
  normalization_factor = particle.normalization_factor;
  
  //if (arrays_assigned){
  //delete[] p;
  //delete[] Etot;
  //delete[] Ekin;
  //delete[] beta;
  //delete[] gamma;
  //delete[] rigidity;
  //delete[] z;
  //}
  
  p.resize(n_pgrid);
  p = particle.p;

  Etot.resize(n_pgrid);
  Etot = particle.Etot;

  Ekin.resize(n_pgrid);
  Ekin = particle.Ekin;
  
  beta.resize(n_pgrid);
  beta = particle.beta;

  gamma.resize(n_pgrid);
  gamma = particle.gamma;

  rigidity.resize(n_pgrid);
  rigidity = particle.rigidity;

  //p=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) p[ip]=particle.p[ip];
  //Etot=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) Etot[ip]=particle.Etot[ip];
  //Ekin=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) Ekin[ip]=particle.Ekin[ip];
  //beta=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) beta[ip]=particle.beta[ip];
  //gamma=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) gamma[ip]=particle.gamma[ip];
  //rigidity=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) rigidity[ip]=particle.rigidity[ip];
  
  z.resize(n_zgrid);
  z = particle.z;

  r.resize(n_rgrid);
  r = particle.r;

  x.resize(n_xgrid);
  x = particle.x;

  y.resize(n_ygrid);
  y = particle.y;

  //z=new double[n_zgrid];
  //for(int iz=0; iz<n_zgrid; iz++) z[iz]=particle.z[iz];
  
  //if(n_spatial_dimensions==2)
  //{
  //      if (arrays_assigned){
  //	delete[] r;
  //  }
  //  r=new double[n_rgrid];
  //  for(int ir=0; ir<n_rgrid; ir++) r[ir]=particle.r[ir];
  //}
  
  //if(n_spatial_dimensions==3)
  //{
  ///  if (arrays_assigned){
  //	delete[] x;
  //	delete[] y;
  //  }
  //  x=new double[n_xgrid];
  //  for(int ix=0; ix<n_xgrid; ix++) x[ix]=particle.x[ix];
  //  y=new double[n_ygrid];
  //  for(int iy=0; iy<n_ygrid; iy++) y[iy]=particle.y[iy];
  //}
  
  arrays_assigned=1;
  cr_density =particle. cr_density;
  return *this;
}

