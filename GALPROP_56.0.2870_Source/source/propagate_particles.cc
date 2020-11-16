//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * propagate_particles.cc *                      galprop package * 1/29/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include "galprop_classes.h"
#include "galprop_internal.h"

#include <PropelOperatorSplitting.h>

#include <cstring>
#include <algorithm>
#include <sstream>
#include <set>
#include <stdexcept>

using namespace std;

#include <ErrorLogger.h>

int Galprop::propagate_particles() { //AWS20050816

  //cout<<">>>>propagate_particles"<<endl;

  INFO("Entry");
  
  int i,net_iter;           //IMOS20032901
  Particle particle;
  
  if (0 == galdef.warm_start) 
    for (i = 0; i < n_species; ++i) 
      gcr[i].cr_density = 0.; //AWS20010121
  else if (1 == galdef.warm_start) 
    read_gcr(); //AWS20010121  
  else if (2 == galdef.warm_start) {
    read_gcr();
    return 0; // This is for simply loading and then generating skymaps -- a `hot' start. Useful for testing with 3D distributions that take a long time to calculate.
  } else {

    std::ostringstream buf;
    buf << "Bad value " << galdef.warm_start << " for galdef.warm_start passed to propagate_particles";

    INFO(buf.str());
    exit(-1);

  }
  
  particle = gcr[0];
  particle.create_transport_arrays();

  //galdef.network_iterations gives the total number of iterations for protons
  //(is automatically increased if the following is larger)
  //galdef.network_iter_compl specifies the total number of iterations for
  //nuclei.
  //This was implemented to speed up damping calculations
  
  int network_iterations = max(galdef.network_iterations, galdef.network_iter_compl);

  Particle* protons = nullptr;
  std::vector<Particle*> particlesToPropagate;
  particlesToPropagate.reserve(n_species);  //Maximum number we have to worry about

  //Use the new propel classes.  
  PropelBase *propelClass = nullptr;
  if (galdef.solution_method == 3) {
     propelClass = new PropelCrankNicolsonD(particle, galdef);
  } else if (galdef.solution_method == 4) {
     propelClass = new PropelCrankNicolsonVector<double>(particle, galdef);
  }

  //Count backwards since we test how many iterations are left.
  for (int r_net_iter=network_iterations; r_net_iter>=1; --r_net_iter) { //IMOS20030129
    
    //For pretty printing, turn that into number of iterations
    net_iter = network_iterations-r_net_iter+1;
    //cout<<"    Network iteration "<<net_iter<<endl;   //IMOS20030129
    
    std::set<Particle*> particlesDone;
    bool loop = true;
    
    while (loop) {
      
      particlesToPropagate.clear();    
      
      // Do protons only ?
      if ( r_net_iter > galdef.network_iter_compl ) {
	
	if ( ! protons ) {
	  
	  for ( i = n_species-1; i >= 0; --i ) {
	    
	    if ( gcr[i].Z == 1 && gcr[i].A == 1 ) {
	      protons = &gcr[i];
	      break;
	    }
	    
	  }
	  
	}
	
	particlesToPropagate.push_back(protons);
	
	loop = false;
	
      } else {
	
	//Loop over the entire species array and add particles that aren't done and have all dependencies fulfilled.
	for (size_t i = 0; i < size_t(n_species); ++i) {
	  
	  std::set<Particle*>::iterator sit = particlesDone.find(&gcr[i]);

	  if ( sit == particlesDone.end() ) {
	    
	    std::map<Particle*, std::valarray<double> >::iterator mit;
	    
	    bool add = true;
	    
	    for (mit = gcr[i].dependencies.begin(); mit != gcr[i].dependencies.end(); ++mit) {
	      
	      if ( particlesDone.find(mit->first) == particlesDone.end() ) {
		// Do not exclude if the only difference is the K_electron and 
		// we are currently dealing with non-K_electron particle
		if ( ! ( gcr[i].K_electron == 0 && 
			 mit->first->K_electron == 1 && 
			 gcr[i].A == mit->first->A && 
			 gcr[i].Z == mit->first->Z ) ) {
		  //std::cerr<<" Particle "<<mit->first->name<<" not yet propagated, dependency for "<< gcr[i].name<<std::endl;
		  add = false;
		  break;
		}
	      }
	      
	    }
	    
	    if (add) 
	      particlesToPropagate.push_back(&gcr[i]);
	    
	  }
	  
	}
        
	if (particlesToPropagate.size() == 0) {
	  
	  if (particlesDone.size() != size_t(n_species)) {
	    ERROR("Nothing found to propagate before all particles have been propagated.");
	    throw(std::logic_error("Error in dependency calculation during propagation"));
	  } else {
	    loop = false;
	    break;
	  }
	  
	}
	
      }
      
      //This may be parallelized in the future
      std::ostringstream buf;
      buf<<"Number of Particles to propagate this loop: "<< particlesToPropagate.size();
      INFO(buf.str());
      for (size_t i = 0; i < particlesToPropagate.size(); ++i) {
	
	particle = *particlesToPropagate[i];
	
	//Reset Dxx once iterations for wave damping have finished
	//This should be made more elegant
	if ( r_net_iter <= galdef.network_iter_compl ) 
	  particle.Dxx = 0;
	
	fill_transport_arrays(particle);
	
	if (0 != gen_secondary_source(particle)) 
	  return 1;
	
	if(galdef.verbose==10) particle.secondary_source_function.print();
	
	if(galdef.verbose>= 0) {
	  
	  std::ostringstream lvl0Buf;
	  lvl0Buf << "\n Network iteration "<<net_iter   //IMOS20030129
		  <<" species "<<i<<" "<<particle.name<<" (Z,A) = ("<<particle.Z
		  <<","<<particle.A<<")"; 
	  INFO(lvl0Buf.str());
	  
	}
	
	if (galdef.solution_method >= 3 && galdef.solution_method <= 9) {
	  (*propelClass)(particle);
	} else {
	  if (propel(particle) != 0) 
	    return 1;
	}
	
	*particlesToPropagate[i] = particle;

	
	particlesDone.insert(particlesToPropagate[i]);
	
	if (galdef.verbose >= 1) {
	  
	  std::ostringstream lvl1Buf;
	  lvl1Buf <<"Network iteration "<<net_iter     //IMOS20030129
		  <<" species "<<i<<"  "<<particle.name;
	  INFO(lvl1Buf.str());
	  particlesToPropagate[i]->cr_density.print();
	  particle.cr_density.print();
	  particle.print();
	}
      }
      
    }
    
    // test of electron propagation vs analytical calculations (only runs for galdef.DM_int0>0)  IMOS20061030
    if (99 == galdef.DM_int0) {

      INFO(" ***** Analytical test of electron propagation  *****"
	   " ***** results are stored in DM_positrons array *****");
	
	int iDM_positrons=-1;// identify test DM_positrons array
	for(i=0; i<n_species; i++)  
	  if ("DM_positrons" == gcr[i].name) //strcmp(gcr[i].name,"DM_positrons")==0)
	    {
	      iDM_positrons=i;
	      cout<<"  DM_positrons found as species #"<<iDM_positrons<<endl;
	    }
	if(iDM_positrons==-1) { 
           INFO("  DM_positrons not found!");  
        } else {
	
	// assigning analytical solution to DM_positrons array
	const double elossconst=32./9.*Pi*pow(Rele/Mele*1.e3,2)*C*galdef.DM_double7*1.e-9;// 1/(GeV s)
	int iz1=0, iz2=gcr[0].n_zgrid-1;
	if(!galdef.output_gcr_full) iz1=iz2=(int)(1.e-6-galdef.z_min/galdef.dz);//z=0,Galactic plane
#pragma omp parallel for schedule(dynamic) default(shared)
	for(int ir=0; ir<gcr[0].n_rgrid;  ir++)
	  {
             cout<<">> r= "<<galaxy.r[ir]<<endl;
	    for(int iz=iz1; iz<iz2+1; cout<<" "<<galaxy.z[iz]<<":", iz++)
	      {
		for(int ip=0; ip<gcr[iDM_positrons].n_pgrid; ip++)
		  {
		    gcr[iDM_positrons].cr_density.d2[ir][iz].s[ip]= C/(4.*Pi)/1.e3
		      *eprop(gcr[iDM_positrons].Ekin[ip]/1.e3,galaxy.r[ir],galaxy.z[iz],
			     galdef.r_max,galdef.DM_double6,galdef.z_max,elossconst,
			     galdef.DM_double8,galdef.D0_xx*pow(kpc2cm,-2),galdef.D_g_1,configure.fGlobalDataPath);
		  }
	      }
	  }
        }
	INFO(" *****   End of test of electron propagation    *****");
      }
    
    // convert density per momentum to flux per (KE/nucleon) for output NOT done in store_gcr
    // only for nuclei
    // for(i=0; i<n_species; i++) if(gcr[i].A!=0) gcr[i].cr_density*=gcr[i].A;
    
    if (galdef.proton_norm_flux > 0) 
      nuclei_normalize(); //IMOS20030129
    
    if (1 == galdef.primary_electrons && 
	r_net_iter <= galdef.network_iter_compl && 
	galdef.electron_norm_flux > 0) 
      electrons_normalize();             //IMOS20030129
    
  } //nuc_iter IMOS20030129
  
  //for (int i = 0; i < n_species; ++i) {

  //cout << gcr[i].name << endl;

  //if ("secondary_positrons" == gcr[i].name) {

  //  for (int iR = 0; iR < 1; ++iR)
  //	for (int iZ = galaxy.n_zgrid/2; iZ < galaxy.n_zgrid/2+1; ++iZ)
  //	  for (int iE = 0; iE < galaxy.n_pgrid; ++iE)
  //	    cout << iR << " " << iZ << " " << iE << " " << galaxy.p[iE] << " " << gcr[i].cr_density.d2[iR][iZ].s[iE] << endl;//" " << gcr[i].secondary_source_function.d2[iR][iZ].s[iE] << endl;

  //}

  //}

  //Gulli20070810 clean up the gcr transport arrays
  //This is unnecessary
  /*
  ostringstream cleanBuf;
  cleanBuf << "Cleaning up for: "<<n_species<<" number of species.";
  INFO(cleanBuf.str());
  for(i=0; i<n_species; i++)
   cout<< "Deleting transport arrays:" << i << endl;
   gcr[i].delete_transport_arrays();
  particle.delete_arrays();
  */

  delete propelClass;
  
  INFO("Exit");
  return 0;
}


