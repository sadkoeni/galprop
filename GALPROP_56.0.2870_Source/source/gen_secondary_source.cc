//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_secondary_source.cc *                     galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>

using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"

#include <ErrorLogger.h>

#include <Nuclei_Interface.h>

int Galprop::gen_secondary_source(Particle& particle) {

  INFO("Entry");

  //AWS20091009 fixed the changes introduced by TAP in r338: particle.name is a char array not a string

  //string particle_name;                                             //AWS20091009
  //particle_name = particle.name; // transfer from char to string    //AWS20091009

  const string particle_name = particle.name; 

  if (galdef.verbose >= 1) {

    ostringstream buf;
    buf << "Generating for (Z, A, K_electron) = (" 
	<< particle.Z << ", "
	<< particle.A << ", "
	<< particle.K_electron << ")";
    INFO(buf.str());

  }

  particle.secondary_source_function = 0;

  if (particle_name == "primary_electrons") { // Primary

    INFO("Exit");
    return 0;

  }

  if (particle_name == "secondary_electrons" && galdef.secondary_electrons || 
      particle_name == "secondary_positrons" && galdef.secondary_positrons) {

    gen_secondary_positron_source(particle);

    INFO("Particle name: " + particle_name);
    ostringstream buf;
    buf << "secondary_source_function.max() = " << particle.secondary_source_function.max();
    INFO(buf.str());
    if ( particle.secondary_source_function.max() == 0 ){
      INFO("Source function is 0");
    } else {
      INFO("Source function is non-zero");
    }

    INFO("Exit");
    return 0;

  }

  if (particle_name == "knock_on_electrons") {

    gen_knock_on_electron_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "secondary_protons") {

    gen_secondary_proton_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "secondary_antiprotons") {

    gen_secondary_antiproton_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "tertiary_antiprotons") {

    gen_tertiary_antiproton_source(particle);
    INFO("Exit");
    return 0;

  }

  /*if(strcmp(particle.name,"primary_electrons"    )==0)       // since primary
    {if(galdef.verbose>=1) cout<<" <<<<gen_secondary_source"<<endl  ;  return 0;}

   if(strcmp(particle.name,"secondary_electrons"  )==0)       // secondary e-
      {gen_secondary_positron_source  (particle);  return 0;}

   if(strcmp(particle.name,"knock_on_electrons"  )==0)        // knock-on e-
      {gen_knock_on_electron_source  (particle);  return 0;}

   if(strcmp(particle.name,"secondary_positrons"  )==0)       // secondary e+
      {gen_secondary_positron_source  (particle);  return 0;}

   if(strcmp(particle.name,"secondary_antiprotons")==0)       // secondary p-
      {gen_secondary_antiproton_source(particle);  return 0;}

   if(strcmp(particle.name,"tertiary_antiprotons")==0)        // tertiary p- IMOS20000605
      {gen_tertiary_antiproton_source(particle);  return 0;}  //             IMOS20000605

   if(strcmp(particle.name,"secondary_protons")==0)           // secondary protons IMOS20000605
      {gen_secondary_proton_source(particle);  return 0;}     //                   IMOS20000605
  */

  // DM source

  if (particle_name == "DM_positrons" ||
      particle_name == "DM_electrons" ||
      particle_name == "DM_antiprotons") {

    gen_DM_source(particle);
    INFO("Exit");
    return 0;

  }

  /*if(strcmp(particle.name,"DM_positrons"  )==0)              // DM e+:  IMOS20050912
    {gen_DM_source(particle);  return 0;}
  
  if(strcmp(particle.name,"DM_electrons"  )==0)              // DM e-:  IMOS20050912
    {gen_DM_source(particle);  return 0;}

  if(strcmp(particle.name,"DM_antiprotons")==0)              // DM p-:  IMOS20050912
    {gen_DM_source(particle);  return 0;}
  */

  // Use the pre-calculated cross-sections and dependencies
  std::map<Particle*, std::valarray<double> >::iterator it;

  ostringstream buf;
  buf<<"Generating secondary source for (Z, A, K_electron) = (" 
     << particle.Z << ", "
     << particle.A << ", "
     << particle.K_electron << ")";
  INFO(buf.str());


  for (it = particle.dependencies.begin(); it != particle.dependencies.end(); ++it) {

    valarray<double> cross_section = it->second[slice(0,particle.n_pgrid,1)];

    buf.str("");
    buf<<"Adding contributions from "<<it->first->name<<".  Maximum cross section: "<<cross_section.max();
    INFO(buf.str());

    //Some pre-calculations to speed things up a bit
    cross_section *= 1e-27; // mb -> cm^2
    cross_section *= C; 
    cross_section *= it->first->beta;
    
    if (2 == galdef.n_spatial_dimensions)
      for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	for(int iz=0; iz<gcr[0].n_zgrid; iz++) {
          const double gas = (galaxy.n_HI.d2[ir][iz].s[0] + 2.*galaxy.n_H2.d2[ir][iz].s[0] + galaxy.n_HII.d2[ir][iz].s[0]);
	  for(int ip=0; ip<gcr[0].n_pgrid; ip++)
	    particle.secondary_source_function.d2[ir][iz].s[ip] += it->first->cr_density.d2[ir][iz].s[ip]*cross_section[ip]*gas;
        }
    
    if (3 == galdef.n_spatial_dimensions)
      for(int ix=0; ix<gcr[0].n_xgrid; ix++)
	for(int iy=0; iy<gcr[0].n_ygrid; iy++)
	  for(int iz=0; iz<gcr[0].n_zgrid; iz++) {
            const double gas = (galaxy.n_HI.d3[ix][iy][iz].s[0] + 2.*galaxy.n_H2.d3[ix][iy][iz].s[0] + galaxy.n_HII.d3[ix][iy][iz].s[0]);
	    for(int ip=0; ip<gcr[0].n_pgrid; ip++)
	      particle.secondary_source_function.d3[ix][iy][iz].s[ip] += it->first->cr_density.d3[ix][iy][iz].s[ip] * cross_section[ip]*gas;
          }
    
      
    if (galdef.verbose == -601 || galdef.verbose == 10) { // selectable debug
       cout << "Progenitor density of " << it->first->name << " for " << particle.name << ":\n";
       it->first->cr_density.print();
    }
      
    // RADIOACTIVE DECAY SOURCE (beta+/-, EC) AWS20010806 IMOS20010816

    const double branching_ratio = it->second[particle.n_pgrid];
    const double t_half = it->second[particle.n_pgrid+1];
    
    if (branching_ratio > 0 && t_half > 0.) {

       //Some precalculations
       std::valarray<double> mult(1., particle.gamma.size());
       mult /= it->first->gamma;
       mult *= branching_ratio/t_half * log(2.0);

       if (2 == galdef.n_spatial_dimensions)
          for(int ir=0; ir<gcr[0].n_rgrid; ir++)
             for(int iz=0; iz<gcr[0].n_zgrid; iz++)
                for(int ip=0; ip<gcr[0].n_pgrid; ip++)
                   particle.secondary_source_function.d2[ir] [iz].s[ip] += mult[ip] * it->first->cr_density.d2[ir][iz].s[ip];

       if (3 == galdef.n_spatial_dimensions)
          for(int ix=0; ix<gcr[0].n_xgrid; ix++)
             for(int iy=0; iy<gcr[0].n_ygrid; iy++)
                for(int iz=0; iz<gcr[0].n_zgrid; iz++)
                   for(int ip=0; ip<gcr[0].n_pgrid; ip++)
                      particle.secondary_source_function.d3[ix][iy][iz].s[ip] += mult[ip] * it->first->cr_density.d3[ix][iy][iz].s[ip];
    }
	  
  }  //  dependencies
     
  if(galdef.verbose==-603 ||galdef.verbose==10) { // selectable debug
      
    cout<<" secondary source for "<<particle.name<<":\n";
    particle.secondary_source_function.print();
    
  }
  
  if(galdef.verbose==10) 
    particle.print();
  
  INFO("Exit");
  
  return 0;
    
}
