
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_DM_source.cc *                             galprop package * 9/09/2005 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|


//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// The routine gen_DM_source calculates the source functions of the products of the
// dark matter (DM) particle annihilation [cm^-3 s^-1 MeV^-1].
// The routine can be used to calculate source function of positrons, electrons,
// and antiprotons.
// Use gen_DM_emiss to define gamma-ray emissivity (cm^-3 s^-1 MeV^-1)
// in terms (dn/dEdt *c/4pi), where n is the number density, c is speed of light.
// The user must use the parameters DM_double0-9 and DM_int0-9 (galdef-file) to 
// specify the Galactic DM profile, branching, decay channels, and spectra (see 
// the template below). The DM profile is defined in the DM_profile routine.
// The profile is then averaged over the grid step (dR,dz) or (dx,dy,dz) with 
// a smaller step: normally 1/10 of the grid size.          IMOS20050912
//
// See example in Moskalenko I.V., Strong A.W. 1999, Phys. Rev. D 60, 063003
// and realization below.
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!
using namespace std;
#include"galprop_classes.h"
#include"galprop_internal.h"

//#include <fort_interface.h>

#include <cstring>

//extern "C" void RHO_DARKSUSY_F77(double*,double*,double*,double*); //IMOS20060901

#include <ErrorLogger.h>
void printProfile(Particle& particle);
double dbar_injection_spectrum(double Ekin);

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int Galprop::gen_DM_source(Particle& particle) {
   
  INFO("Entry");
  
  ostringstream os;
  os<<"generating "<<particle.name<<" source function for n_spatial_dimensions="
      <<gcr[0].n_spatial_dimensions<<endl;
  INFO(os.str());

  //injection spectrum for WW channel with DM mass 100GeV, from Ibarra's paper
  
  double DMwidth,DMbranching,   // annihilation product distribution
    DMsecondary_spectrum,       // spectrum of secondaries from DM annihilation
    DME0,                       // delta function energy used for Green's function
    DMmass  =galdef.DM_double2, // DM particle mass, in GeV
    DMcs_v  =galdef.DM_double9, // DM <cross_sec*V> -thermally overaged, cm3 s-1 
    dzz=0.01;                   // kpc, gas averaging step
  int stat=0;

  // define the spectra of annihilation products: positrons, electrons, antiprotons
  
  if("DM_positrons" == particle.name) //strcmp(particle.name,"DM_positrons")==0)
    {
       DMwidth     =galdef.DM_double3;
       DMbranching =galdef.DM_double4;
     }

  if("DM_electrons" == particle.name) //strcmp(particle.name,"DM_electrons")==0)
     {
       DMwidth     =galdef.DM_double5;
       DMbranching =galdef.DM_double6;
     }

  if("DM_antiprotons" == particle.name) //strcmp(particle.name,"DM_antiprotons")==0)
     {
       DMwidth     =galdef.DM_double7;
       DMbranching =galdef.DM_double8;
     }
 		cout<<"#######DM distribution #######"<<endl;
        for(int ir=0; ir<gcr[0].n_rgrid; ir++){
				cout<<DM_profile(galaxy.r[ir], 0,0)<<", ";
		}
		cout<<endl<<"###### end of DM distribution #####"<<endl;

// assign the source function (2D)

   if(galaxy.n_spatial_dimensions==2)
     {
       for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	 {
	   for(int iz=0; iz<gcr[0].n_zgrid; iz++)
	     {
	       for(int ip=0; ip<particle.n_pgrid; ip++)
		 {


				 if(particle.name == "DM_antiprotons"){
						 /* if (2*massDM>=Etot) {
						  * 	source.add = 1/2 * DM_profile(r)**2 /massDM*DMcsv * (spectra_ww + spectra_bb);
						  * 	}
						  * else {source.add = 0;}
						  */
						 //if(DMmass<particle.Etot[ip]*1e-3) {INFO("");continue;}//Was commented out since its included in the injection spectra. 
						 //double sec_source = 1/2 * DMcs_v * pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz)/DMmass, 2)*dbar_inj_spectrum[ip];//todo: fix p values
						 double sec_source = 0.5 * DMcs_v * pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz)/DMmass, 2) * dbar_injection_spectrum(particle.Ekin[ip]*1.0e-3);//done: fix p values

						 particle.secondary_source_function.d2[ir][iz].s[ip] += sec_source;
						// =0.5* (rho(r)/m_DM)**2 * <csv> * dN/dT (T) 
						cout<<endl<<"#########DM filling pbar #######"<<endl;
						cout<<ip<<" "<<particle.Ekin[ip]<<" "<<DMcs_v<<" "<<pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz)/DMmass, 2)<<" "<< dbar_injection_spectrum(particle.Ekin[ip]*1.0e-3)<<endl;
						cout<<sec_source<<endl;






				 }
				 /*
// test of electron propagation vs analytical calculations IMOS20061030
// to run test, assign galdef.DM_int0=99, other parameters:
// galdef.DM_double6 - the half thickness of the disk source distribution (e.g. 0.1 kpc), the source 
//                     distribution is uniform within the disk; normalization =1 at 1 GeV
// galdef.DM_double7 - the photon field energy density (e.g. 1 eV/cc)
// galdef.DM_double8 - the injection spectral index of electrons (e.g. 2.4)
		   if(abs(galdef.DM_int0)==99 && particle.A==0) 
		     {
		       if("DM_electrons" == particle.name) //strcmp(particle.name,"DM_electrons")==0) //numerical  solution "DM_electrons"
                          if (particle.Ekin[ip] > galdef.inj_Ekin_min && particle.Ekin[ip] < galdef.inj_Ekin_max) {
			 particle.secondary_source_function.d2[ir][iz].s[ip]= 
			   (galdef.DM_double6 <= fabs(galaxy.z[iz])) ? 0.:
			   C/4./Pi*pow(particle.Ekin[ip]*1e-3,-galdef.DM_double8); //The spectrum is normalized at 1 GeV
                          }
			 if("DM_positrons" == particle.name) //strcmp(particle.name,"DM_positrons")==0) //analytical solution "DM_positrons"
			 particle.secondary_source_function.d2[ir][iz].s[ip]=0.;
		       continue;		     
		     }
// end of the test area
		   if(galdef.DM_int1==9) // Green's function to work with DarkSUSY IMOS20060901
		     {
		       if(DME0<particle.Ekin[ip] || DME0/DMwidth>particle.Ekin[ip]) continue;
		       particle.secondary_source_function.d2[ir][iz].s[ip]
			 +=pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz),2)
			 *DMsecondary_spectrum*DMbranching/4./Pi*C;
		       continue;
		     }
		   if(particle.Etot[ip]*1.e-3<=DMmass) 
		     particle.secondary_source_function.d2[ir][iz].s[ip]+= DMcs_v*
		       pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz)/DMmass,2)
		       *C/4./Pi*DMbranching*exp(-pow((DMmass-particle.Etot[ip]*1.e-3)/DMwidth,2))/DMmass*1.e-3;
			   */
		 } // ip
	     }  //  iz
	 }  //  ir
     }  //  particle.n_spatial_dimensions==2
   
// assign the source function (3D)

   if(galaxy.n_spatial_dimensions==3)
     {
       for(int ix=0; ix<gcr[0].n_xgrid; ix++)
	 {
	   for(int iy=0; iy<gcr[0].n_ygrid; iy++)
	     {
	       for(int iz=0; iz<gcr[0].n_zgrid; iz++)
		 {
		   for(int ip=0; ip<particle.n_pgrid; ip++)
		     {
		       if(galdef.DM_int1==9) // Green's function to work with DarkSUSY IMOS20060901
			 {
			   if(DME0<particle.Ekin[ip] || DME0/DMwidth>particle.Ekin[ip]) continue;
			   particle.secondary_source_function.d3[ix][iy][iz].s[ip]
			     +=pow(DM_profile_av(galaxy.r[ix], galaxy.r[iy], galaxy.z[iz], galaxy.dx, galaxy.dy, galaxy.dz, dzz),2)
			     *DMsecondary_spectrum*DMbranching/4./Pi*C;
			   continue;
			 }
		       if(particle.Etot[ip]*1.e-3<=DMmass) 
			 particle.secondary_source_function.d3[ix][iy][iz].s[ip]+= DMcs_v*
			   pow(DM_profile_av(galaxy.x[ix], galaxy.y[ix], galaxy.z[iz], galaxy.dx, galaxy.dy, galaxy.dz, dzz)/DMmass,2)
		       *C/4./Pi*DMbranching*exp(-pow((DMmass-particle.Etot[ip]*1.e-3)/DMwidth,2))/DMmass*1.e-3;
		     } //ip
		 }  //  iz
	     }  //  iy
	 }  //  ix
     }  //  particle.n_spatial_dimensions==3
 
 // test printout

   if(galdef.verbose>=2)
     {
       cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
       particle.secondary_source_function.print();
     }
   INFO("Exit");
   return stat;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int Galprop::gen_DM_emiss()
{
   INFO("Entry");
   double 
     DMmass      =galdef.DM_double2, // DM particle mass
     DMcs_v      =galdef.DM_double9, // DM <cross_sec*V> -thermally overaged, cm3 s-1 , this should be of the order 1e-26, need to check exact number
     DMbranching =0.1,
     dzz=0.01;                       // kpc, gas averaging step
   int stat=0;

   galaxy.DM_emiss=0.;

// define the spectra of annihilation products: gammas
   
   if(galdef.n_spatial_dimensions==2)
     {
       ostringstream os;
       os<<"generating DM emissivity for n_spatial_dimensions="<<galdef.n_spatial_dimensions;
       INFO(os.str());
       for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	 {
	   for(int iz=0; iz<gcr[0].n_zgrid; iz++)
	     {
               for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
		 {
		   if(galaxy.E_gamma[iEgamma]*1.e-3>DMmass) 
		     {
		       galaxy.DM_emiss.d2[ir][iz].s[iEgamma]=0;
		       continue;
		     }
		   galaxy.DM_emiss.d2[ir][iz].s[iEgamma]= DMcs_v *DMbranching/(4.*Pi)// sr^-1 IMOS20060420, Need to understand where this formula comes from. 
		     *pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz)/DMmass,2)
		     /galaxy.E_gamma[iEgamma];
		 }
	     }
	 }
     }
   if(galdef.n_spatial_dimensions==3)
     {
       ostringstream os;
       os<<"generating DM emissivity for n_spatial_dimensions="<<galdef.n_spatial_dimensions;
       INFO(os.str());
       for(int ix=0; ix<gcr[0].n_rgrid; ix++)
	 {
	   for(int iy=0; iy<gcr[0].n_rgrid; iy++)
	     {
	       for(int iz=0; iz<gcr[0].n_zgrid; iz++)
		 {
		   for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
		     {
		       if(galaxy.E_gamma[iEgamma]*1.e-3>DMmass) 
			 {
			   galaxy.DM_emiss.d3[ix][iy][iz].s[iEgamma]=0;
			   continue;
			 }
		       galaxy.DM_emiss.d3[ix][iy][iz].s[iEgamma]=  DMcs_v *DMbranching/(4.*Pi) // sr^-1 IMOS20060420
			 *pow(DM_profile_av(galaxy.x[ix], galaxy.y[ix], galaxy.z[iz], galaxy.dx, galaxy.dy, galaxy.dz, dzz)/DMmass,2)
			 /galaxy.E_gamma[iEgamma];
		     }
		 }
	     }
	 }
     }
   INFO("Exit");
   return(stat);
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

double Galprop::DM_profile(double Xkpc, double Ykpc, double Zkpc)
{
  double R=sqrt(Xkpc*Xkpc+Ykpc*Ykpc+Zkpc*Zkpc),
    Rsun =8.5,                     //kpc, galactocentric distance of the solar system 
    Rc         =galdef.DM_double0, //core radius
    rho0       =galdef.DM_double1; //local DM mass density| As its used, this is actually the normalization factor and needs to reproduce a local DM density at the sun of about 0.39GeVcm-3. 
  int profile_key =galdef.DM_int0; //profile type
	  double alpha = 0.17;
	  double rho=rho0;
  
  switch(profile_key)
    {//important ones are 0 (NFW), 1 (isothermal), and 11 (Einasto)
    case 0:   //NFW profile
	  Rc = 24.42; //As stated in Ibarra@s paper
	  rho = rho0 * (8.5/Rc*pow(1.+8.5/Rc,2));//This makes sure we use the correct normalization.
      return(rho/(R/Rc*pow(1.+R/Rc,2)));
      
    case 1:   //isothermal profile
	  Rc = 4.38;
	  rho = rho0*(pow(8.5, 2) + pow(Rc, 2));
	  return(rho /(pow(R, 2) +pow(Rc, 2)));//according to ibarra paper arxiv.1209.5539v2
      //return(rho0*(pow(Rc,2)+pow(Rsun,2))/(pow(Rc,2)+pow(R,2)));
      
    case 2:   //Evans profile
      return(rho0*pow(pow(Rc,2)+pow(Rsun,2),2)/(3.*pow(Rc,2)+pow(Rsun,2))
	     *(3.*pow(Rc,2)+pow(R,2))/pow(pow(Rc,2)+pow(R,2),2));
      
    case 3:   //alternative profile
      return(rho0*pow(Rc+Rsun,2)/pow(Rc+R,2));
      
    case 9:   //DarkSUSY profile (use only if the DarkSUSY and GALPROP combined) IMOS20060901
      rho0 = rho_darksusy_cc(Xkpc,Ykpc,Zkpc);

      if(rho0<0.)
	{
	  FATAL("gen_DM_source: rho_darksusy() function is not defined");
	  exit(0);
	}
      return(rho0);
	case 11: 	//Einasto
	  Rc = 28.44;//according to ibarra paper arxiv.1209.5539v2
	  rho =rho0 / pow(2.71828, -2/alpha * (pow(8.5/Rc, alpha)-1));
	  return(rho*pow(2.71828, -2/alpha * (pow(R/Rc, alpha)-1)));
	case 10: 	//For inserting new profiles
	  return(rho0);

    default:
      return(rho0);
    }
}
  

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

 double Galprop::DM_profile_av(double r,double z,double dr,double dz,double dzz)
   {  
     double DM_profile_av_=0.0;
     int nuse=0;
     
     for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
       for (double rr=r-dr/2.; rr<=r+dr/2.; rr+=dr/10.)
	 { 
	   if (rr<0.) continue;
	   DM_profile_av_+=DM_profile(rr,0,zz);
	   nuse++; 
	 }
     return (DM_profile_av_/nuse);
   }
 
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
 
 double Galprop::DM_profile_av(double x,double y,double z,double dx,double dy,double dz,double dzz)
   {  
     double DM_profile_av_=0.0;
     int nuse=0;
     
     for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
       for (double xx=x-dx/2.; xx<=x+dx/2.; xx+=dx/10.)
	 for (double yy=y-dy/2.; yy<=y+dy/2.; yy+=dy/10.)
	   {
	     DM_profile_av_+=DM_profile(xx,yy,zz);
	     nuse++;
	   }
     return DM_profile_av_/nuse;
   }
 
double dbar_injection_spectrum(double Ekin){
		//This function atm is for WW channel at DM mass of 100GeV
		if(Ekin>100){return 0;}
		double Ekin_values[80] = {0.0010865359752156843, 0.0012547121698960636, 0.001448919009766593, 0.0016731855697525507, 0.0019321645529926125, 0.0022312288172514285, 0.002576582841881917, 0.002975391447865205, 0.003435928441393822, 0.003967748264803714, 0.004581884215978187, 0.005291077348544593, 0.006110040801697761, 0.007055765043518109, 0.008147870360456343, 0.009409013905840187, 0.010865359752156842, 0.012547121698960637, 0.01448919009766593, 0.016731855697525506, 0.019321645529926122, 0.022312288172514287, 0.02576582841881917, 0.029753914478652048, 0.03435928441393822, 0.039677482648037145, 0.04581884215978187, 0.052910773485445935, 0.061100408016977616, 0.07055765043518109, 0.08147870360456343, 0.09409013905840186, 0.10865359752156842, 0.12547121698960637, 0.1448919009766593, 0.16731855697525508, 0.19321645529926124, 0.22312288172514286, 0.25765828418819176, 0.2975391447865205, 0.3435928441393822, 0.39677482648037143, 0.4581884215978187, 0.5291077348544593, 0.611004080169776, 0.7055765043518109, 0.8147870360456343, 0.9409013905840187, 1.086535975215684, 1.2547121698960637, 1.448919009766593, 1.6731855697525506, 1.9321645529926126, 2.231228817251429, 2.576582841881917, 2.9753914478652046, 3.4359284413938216, 3.967748264803714, 4.581884215978187, 5.291077348544593, 6.110040801697761, 7.055765043518108, 8.147870360456343, 9.409013905840187, 10.86535975215684, 12.547121698960636, 14.489190097665931, 16.731855697525507, 19.321645529926123, 22.31228817251429, 25.765828418819176, 29.753914478652046, 34.35928441393822, 39.67748264803714, 45.818842159781866, 52.91077348544593, 61.10040801697761, 70.55765043518109, 81.47870360456344, 94.09013905840187};
		double dN_dEkin_values[80] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00020445265711969667, 8.85243534409128e-05, 0.00045995359097008394, 0.0005974551002721614, 0.0050587762694545855, 0.0006471517289325627, 0.0010777127794284202, 0.0015678777931574786, 0.0018103017238342454, 0.0017636143469068272, 0.004630163088415025, 0.00424047244812508, 0.004326533010243058, 0.003337329956471753, 0.0033398680289233092, 0.004072700306898642, 0.005336221883631368, 0.004957369004565674, 0.005772424451609159, 0.00661848248791319, 0.005478429882650507, 0.006944444457013567, 0.00776384196321989, 0.0074474221912903, 0.008330623022377569, 0.007653525760642393, 0.007822795107481543, 0.008919690248459178, 0.008582167498143044, 0.00920441254822519, 0.00909125110605348, 0.0086896010638846, 0.008538970152555161, 0.008296507656192885, 0.008139645332663824, 0.007185398390261336, 0.007357071108643559, 0.006626835024086447, 0.006132185749703947, 0.005440164167534501, 0.005090182894285346, 0.00452857061893763, 0.003932731544073425, 0.0035953761571021698, 0.0029367710579046107, 0.002590521255222876, 0.0020816230665938193, 0.0018065067722537547, 0.0013530685534869208, 0.0011862273436266594, 0.0008511446606796828, 0.0007317504446149783, 0.0005561784824980283, 0.0004081042999680341, 0.000265763998258263, 0.0002092744042523811, 0.00013538471649493318, 9.30516403308249e-05, 6.765856952922518e-05, 3.5708167592017884e-05, 2.090553282938425e-05, 1.3160422818966862e-05, 6.447991150100834e-06, 2.290155324315829e-06, 6.951391875465582e-07, 2.567206913719536e-07, 1.5331786365669463e-08, 0.0, 0.0, 0.0};

		//#####Pseudocode#######
		//loop over Ekin in while loop
		//	if Ekin_value >Ekin: break;
		//diff = Ekin_values[i]-Ekin_values[i-1];
		//Ekin_value = Ekin_values[i];
		//interpolation = (Ekin-Ekin_value) /diff;
		//
		//dN_dEkin_value = dN_dEkin_values[i-1] + (interpolation * (dN_dEkin_values[i]-dN_dEkin_values[i-1]));
		//return dN_dEkin_value;

		int index = 0;
		int N = 80;//length of arrays
		while(index<N and Ekin_values[index]<Ekin){
				if(Ekin_values[index]>=Ekin) break;
				index++;
		}
		cout<<endl<<"index: "<<index<<endl;
		
		double difference = Ekin_values[index]-Ekin_values[index-1];
		double Ekin_value = Ekin_values[index-1];
		double interpolation = (Ekin-Ekin_value)/difference;

		double dN_dEkin_value = dN_dEkin_values[index-1]+(interpolation* ( dN_dEkin_values[index]-dN_dEkin_values[index-1]));//The interpolation should make it better, but is switched off for now
		cout<<endl<<"#### dN_dEkin value: "<<dN_dEkin_value<<" ######"<<endl;
		return dN_dEkin_value;
}
/*
void printProfile(Particle& particle){
		cout<<"#######DM distribution #######"<<endl;
        for(int ir=0; ir<gcr[0].n_rgrid; ir++){
				cout<<DM_profile(galaxy.r[ir], 0,0)<<", ";
		}
		cout<<endl<<"###### end of DM distribution #####"<<endl;
		//return 0;
}

*/
