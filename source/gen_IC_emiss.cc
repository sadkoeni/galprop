
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_IC_emiss.cc *                             galprop package * 4/25/2002 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate inverse Compton emissivity
/*
  ISRF is in Hz eV cm-3 Hz-1
  -> N(nu)*nu = nu*density in Hz cm-3 Hz-1 = ISRF * eV_to_erg / (h_planck * nu)
  integral  density Hz-1  d(nu) = integral (nu*density Hz-1) d(log nu)
  d(log nu) is constant in this ISRF

Etarget= (h_planck*nu)* erg_to_eV ! target photon energy in eV

CR density  gcr.cr_density is in c/4pi * n(E) [ cm s^-1 sr^-1 cm^-3 MeV^-1]

emissivity (cm^-3 s^-1 sr^-1 MeV^-1)=
(c/4pi)*integral[integral(sigma{Egamma,Eelectron,Etarget} N(nu) *nu dlog nu ) n(E)E dlog(E)]
factor=(eV_to_erg / h_planck)  *  log(nu(2)/nu(1)) * log(Ekin_factor)
*/

#include <cassert>
#include <sstream>
#include <valarray>

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

#include <ErrorLogger.h>

bool Galprop::gen_IC_emiss() {

  // Rewrite to clean up code -- TAP 20090311

  INFO("Entry: gen_IC_emiss");
  
  Distribution electrons;

  if (2 == galdef.n_spatial_dimensions) 
    electrons.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);

  if (3 == galdef.n_spatial_dimensions) 
    electrons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
  
  valarray<double> electronGamma(0., electrons.n_pgrid), electronEkin(0., electrons.n_pgrid);
  
  //electronGamma.resize(electrons.n_pgrid);
  //electronEkin.resize(electrons.n_pgrid);

  electrons = 0.;

  bool foundElectrons = false;

  for (auto i = 0; i < n_species; ++i) {

    if (100 == 100*abs(gcr[i].Z) + gcr[i].A) {

      electrons += gcr[i].cr_density;
      foundElectrons = true;

      for (int iP = 0; iP < electrons.n_pgrid; ++iP) {

	electronGamma[iP] = gcr[i].gamma[iP];
	electronEkin[iP] = gcr[i].Ekin[iP];

      }

      ostringstream buf;
      buf << "CR " << gcr[i].name << " found as species #" << i;
      INFO(buf.str());

    }

  }

  if (foundElectrons) {

    // Zero the emissivity arrays

    for (int iComp = 0; iComp < galaxy.n_ISRF_components; ++iComp) 
      galaxy.IC_iso_emiss[iComp] = 0.;

    //for (int iE = 0; iE < galaxy.n_pgrid; ++iE)
    //cout << "ic: " << iE << " " << galaxy.p[iE] << " " << electrons.d2[0][galaxy.n_zgrid/2].s[iE] << endl;

    std::valarray<double> targetEnergy(0., galaxy.ISRF[0].n_pgrid), oneOnNu(0., galaxy.ISRF[0].n_pgrid);

    //targetEnergy.resize(galaxy.ISRF[0].n_pgrid);
    
    // Target photon energy in eV
    
    for (auto i = 0; i < galaxy.ISRF[0].n_pgrid; ++i) {
      targetEnergy[i] = h_planck*galaxy.nu_ISRF[i]*erg_to_eV;
      oneOnNu[i] = 1./galaxy.nu_ISRF[i];
    }
      
    assert(targetEnergy.size() > 0); // We're dead if this is not the case
  
    const double factor = (eV_to_erg/h_planck)*
      std::log(galaxy.nu_ISRF[1]/galaxy.nu_ISRF[0])*std::log(galdef.Ekin_factor);

    std::valarray<double> icCrossSection;

    icCrossSection.resize(galaxy.n_E_gammagrid*electrons.n_pgrid*galaxy.ISRF[0].n_pgrid);

    for (auto iE = 0; iE < galaxy.n_E_gammagrid; ++iE)
      for (auto iP = 0; iP < electrons.n_pgrid; ++iP)
	for (auto iT = 0; iT < galaxy.ISRF[0].n_pgrid; ++iT) {

	  const size_t index = iE*electrons.n_pgrid*galaxy.ISRF[0].n_pgrid +
	    iP*galaxy.ISRF[0].n_pgrid + iT; // explicit encoding -- use it! The `inner loop' is over the target photon field bins so we want these loaded into the cache-line when we do the emissivity calculation

	  const auto epsG = galaxy.E_gamma[iE]/Mele;

	  const auto epsT = targetEnergy[iT]*1e-6/Mele;

	  icCrossSection[index] = fjones_cc(electronGamma[iP], epsT, epsG);

	}

    assert(icCrossSection.size() > 0); // Again, death will ensure ... 

    // 2D case

    if (2 == galdef.n_spatial_dimensions) {

#pragma omp parallel for schedule(dynamic)    
      for (auto iR = 0; iR < electrons.n_rgrid; ++iR) {
	
	for (auto iZ = 0; iZ < electrons.n_zgrid; ++iZ) {
	  
	  for (auto iComp = 0; iComp < galaxy.n_ISRF_components; ++iComp) {
	    
	    valarray<double> isrf_over_nu(0., galaxy.ISRF[0].n_pgrid);
	    
	    for (auto i = 0; i < galaxy.ISRF[0].n_pgrid; ++i) 
	      isrf_over_nu[i] = 
		galaxy.ISRF[iComp].d2[iR][iZ].s[i]*oneOnNu[i];///galaxy.nu_ISRF[i];
	   	  
	    for (auto iE = 0; iE < galaxy.n_E_gammagrid; ++iE) {
	      
	      //const double epsG = galaxy.E_gamma[iE]/Mele;
	      
	      for (auto iP(0); iP < electrons.n_pgrid; ++iP) {
		
		auto sum(0.);
		
		for (size_t iT(0); iT < isrf_over_nu.size(); ++iT) {
	
		  const auto epsT = targetEnergy[iT]*1e-6/Mele;

		  const auto fac = 4.*epsT*electronGamma[iP];
		  
		  if (galaxy.E_gamma[iE] > 
		      fac*Mele*electronGamma[iP]/(1. + fac)) 
		    continue;

		  const auto index = 
		    iE*electrons.n_pgrid*galaxy.ISRF[0].n_pgrid +
		    iP*galaxy.ISRF[0].n_pgrid + iT; 
		    
		  sum += isrf_over_nu[iT]*icCrossSection[index];
		  
		}

		sum *= factor;

		galaxy.IC_iso_emiss[iComp].d2[iR][iZ].s[iE] += 
		  sum*electrons.d2[iR][iZ].s[iP]*electronEkin[iP];
		
	      }

	    }

	  }

	}
	
      }

    }

    // 3D case

    if (3 == galdef.n_spatial_dimensions) {

#pragma omp parallel for schedule(dynamic) 
      for (auto iXY = 0; iXY < electrons.n_xgrid*electrons.n_ygrid; ++iXY) {

	const auto iX = iXY/electrons.n_ygrid;
	const auto iY = iXY%electrons.n_ygrid;

	//#pragma omp critical
	//{
	//std::cout << iXY << " " << iX << " " << iY << std::endl;
	//}
	
	//for (int iX = 0; iX < electrons.n_xgrid; ++iX) {
	
	//for (int iY = 0; iY < electrons.n_ygrid; ++iY) {

	for (auto iZ(0); iZ < electrons.n_zgrid; ++iZ) {
	  
	  for (auto iComp(0); iComp < galaxy.n_ISRF_components; ++iComp) {
	    
	    valarray<double> isrf_over_nu(0., galaxy.ISRF[0].n_pgrid);
	    
	    for (auto i(0); i < galaxy.ISRF[0].n_pgrid; ++i) {
	      
	      isrf_over_nu[i] = galaxy.ISRF[iComp].d3[iX][iY][iZ].s[i]*oneOnNu[i];///galaxy.nu_ISRF[i];
	      
	      //cout << iX << " " << iY << " " << iZ << " " << iComp << " " << galaxy.nu_ISRF[i] << " " << galaxy.ISRF[iComp].d3[iX][iY][iZ].s[i] << endl;
	      
	    }
	    
	    for (auto iE(0); iE < galaxy.n_E_gammagrid; ++iE) {
	      
	      //const double epsG = galaxy.E_gamma[iE]/Mele;
	      
	      for (auto iP(0); iP < electrons.n_pgrid; ++iP) {
		
		auto sum(0.);
		
		for (auto iT(0); iT < isrf_over_nu.size(); ++iT) {
		  
		  const auto epsT = targetEnergy[iT]*1e-6/Mele;
		  
		  const auto fac = 4.*epsT*electronGamma[iP];
		  
		  if (galaxy.E_gamma[iE] > 
		      fac*Mele*electronGamma[iP]/(1. + fac)) 
		    continue;
		  
		  const auto index = 
		    iE*electrons.n_pgrid*galaxy.ISRF[0].n_pgrid +
		    iP*galaxy.ISRF[0].n_pgrid + iT; 
		  
		  sum += isrf_over_nu[iT]*icCrossSection[index];
		  
		}
		
		sum *= factor;
		
		galaxy.IC_iso_emiss[iComp].d3[iX][iY][iZ].s[iE] += 
		  sum*electrons.d3[iX][iY][iZ].s[iP]*electronEkin[iP];
		
	      }
	      
	    }
	    
	  }
	  
	}
	
      }
      
    }

    //}

    electrons.delete_array();

  } else {

    INFO("CR electrons not found. No IC gamma rays will be calculated");

  }

  INFO("Exit: gen_IC_emiss");

  //exit(0);

  return foundElectrons;

}
