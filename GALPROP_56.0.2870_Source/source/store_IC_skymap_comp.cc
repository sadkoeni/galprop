//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_IC_skymap_comp.cc *                     galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <string>
#include <cstring>
#include <valarray>
#include <sstream>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"
#include "fitsio.h"

#include <ErrorLogger.h>
#include <BaseSkyFitsIO.h>

int Galprop::store_IC_skymap_comp(const std::string& IC_type) {

  INFO("Entry");

  int status = 0;

  if (3 == galdef.skymap_format) {
  
    string fileprefix = configure.fOutputDirectory + configure.fOutputPrefix;
    
    if (IC_type == "isotropic") fileprefix += "ics_isotropic_comp_";
    if (IC_type == "anisotropic") fileprefix += "ics_anisotropic_comp_";

    //if(strcmp(IC_type,  "isotropic")==0) fileprefix += "ics_isotropic_comp_";
    //if(strcmp(IC_type,"anisotropic")==0) fileprefix += "ics_anisotropic_comp_";
    
    if (IC_type == "isotropic") { //strcmp(IC_type,  "isotropic")==0){
	
      for (int i_comp=0; i_comp<galaxy.n_ISRF_components; i_comp++){
        ostringstream os;
        os << i_comp+1;
	string filename = fileprefix + os.str();
	filename += "_healpix_";
	filename += galdef.galdef_ID;
	filename += ".gz";
        SM::writeToFits(*galaxy.IC_iso_hp_skymap[i_comp], filename, true, true, "Energy", "MeV");

        if (2 == galdef.gamma_rays && galaxy.IC_iso_rings_hp_skymap.size() > 0) {

           for (int iR(0); iR < galaxy.n_Ring; ++iR) {
              filename = fileprefix + os.str();
              filename += "_ring_";
              ostringstream osr;
              osr << iR+1;
              filename += osr.str();
              filename += "_healpix_";
              filename += galdef.galdef_ID;
              filename += ".gz";
              SM::writeToFits(*galaxy.IC_iso_rings_hp_skymap[i_comp][iR], filename, true, true, "Energy", "MeV");
           }

        }

      }
     
    }

    if (IC_type == "anisotropic") { //strcmp(IC_type,"anisotropic")==0){
	
      for (int i_comp=0; i_comp<galaxy.n_ISRF_components; i_comp++){
        ostringstream osc;
        osc << i_comp+1;
	string filename = fileprefix + osc.str();
	filename += "_healpix_";
	filename += galdef.galdef_ID;
	filename += ".gz";
        SM::writeToFits(*galaxy.IC_aniso_hp_skymap[i_comp], filename, true, true, "Energy", "MeV");

        if (2 == galdef.gamma_rays && galaxy.IC_aniso_rings_hp_skymap.size() > 0) {

           for (int iR(0); iR < galaxy.n_Ring; ++iR) {
              filename = fileprefix + osc.str();
              filename += "_ring_";
              ostringstream osr;
              osr << iR+1;
              filename += osr.str();
              filename += "_healpix_";
              filename += galdef.galdef_ID;
              filename += ".gz";
              SM::writeToFits(*galaxy.IC_aniso_rings_hp_skymap[i_comp][iR], filename, true, true, "Energy", "MeV");
           }

        }
      }
    
    }
  
  } else {
  
    //fitsfile *fptr;       // pointer to the FITS file; defined in fitsio.h
    long nelements;
    long naxes[4];
    string filename;
    double crval[4], cdelt[4];
    
    naxes[0]=galaxy.n_long;
    naxes[1]=galaxy.n_lat;             
    naxes[2]=galaxy.n_E_gammagrid;
    naxes[3]=galaxy.n_ISRF_components;
    nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];

    valarray<float> array(0., nelements);
    
    assert (IC_type == "isotropic" || IC_type == "anisotropic");

    if (IC_type == "isotropic") filename = "ics_skymap_comp_";
    if (IC_type == "anisotropic") filename = "ics_anisotropic_comp_";
    
    if (galdef.skymap_format != 1) { //!< Old output
      
      int i=0;
      
      for (int i_comp=0; i_comp<naxes[3]; i_comp++)
	for (int ip=0; ip<naxes[2]; ip++)
	  for (int ib=0; ib<naxes[1]; ib++)
	    for (int il=0; il<naxes[0]; il++) {

	      if (IC_type == "isotropic")
		array[i] = galaxy.IC_iso_skymap[i_comp].d2[il][ib].s[ip];

	      if (IC_type == "anisotropic")
		array[i] = galaxy.IC_aniso_skymap[i_comp].d2[il][ib].s[ip];

	      array[i]*=pow(galaxy.E_gamma[ip],2);
	      i++;
	      
	    }
      
      crval[0] = galaxy.long_min;
      crval[1] = galaxy. lat_min;
      crval[2] = log10(galaxy.E_gamma_min);
      crval[3] = 0; // IMOS20080114
      
      cdelt[0] = galaxy.d_long;
      cdelt[1] = galaxy.d_lat;
      cdelt[2] = log10(galaxy.E_gamma_factor);
      cdelt[3] = 1;
      
      //Use the standard method to store the skymap
      status = store_skymap(&array[0], naxes, filename, crval, cdelt);

      //Store the rings if needed
      if (2 == galdef.gamma_rays && galaxy.IC_iso_rings_skymap) {

         for (int iR(0); iR < galaxy.n_Ring; ++iR) {

            ostringstream os;
            os << filename << "ring_" << iR+1 << "_";

            i = 0;

            for (int i_comp=0; i_comp<naxes[3]; i_comp++)
               for (int ip=0; ip<naxes[2]; ip++)
                  for (int ib=0; ib<naxes[1]; ib++)
                     for (int il=0; il<naxes[0]; il++) {

                        if (IC_type == "isotropic")
                           array[i] = galaxy.IC_iso_rings_skymap[i_comp].d3[il][ib][iR].s[ip];

                        if (IC_type == "anisotropic")
                           array[i] = galaxy.IC_aniso_rings_skymap[i_comp].d3[il][ib][iR].s[ip];

                        array[i]*=pow(galaxy.E_gamma[ip],2);
                        i++;

                     }

            status = store_skymap(&array[0], naxes, os.str(), crval, cdelt);
         
         }

      }
      
    }
    
    if (galdef.skymap_format == 1 || galdef.skymap_format == 2) { //!< Mapcube output compatable with Glast science tools
      
      int i=0;
      
      for (int i_comp=0; i_comp<naxes[3]; i_comp++) 
	for (int ip=0; ip<naxes[2]; ip++)
	  for (int ib=0; ib<naxes[1]; ib++)
	    for (int il=0; il<naxes[0]; il++) {

	      if (IC_type == "isotropic")
		array[i] = galaxy.IC_iso_skymap[i_comp].d2[il][ib].s[ip];

	      if (IC_type == "anisotropic")
		array[i] = galaxy.IC_aniso_skymap[i_comp].d2[il][ib].s[ip];

	      ++i;
	      
	    }
      
      //filename = filename.replace(filename.find("skymap"), 6, "mapcube");
   
      string file = filename + "mapcube_";
   
      status = store_mapcube_skymap(&array[0], &galaxy.E_gamma[0], galaxy.n_ISRF_components, galaxy.n_E_gammagrid, file, true);
    
      //Store the rings if needed
      if (2 == galdef.gamma_rays && galaxy.IC_iso_rings_skymap) {

         for (int iR(0); iR < galaxy.n_Ring; ++iR) {

            ostringstream os;
            os << filename << "ring_" << iR+1 << "_mapcube_";

            int i=0;

            for (int i_comp=0; i_comp<naxes[3]; i_comp++) 
               for (int ip=0; ip<naxes[2]; ip++)
                  for (int ib=0; ib<naxes[1]; ib++)
                     for (int il=0; il<naxes[0]; il++) {

                        if (IC_type == "isotropic")
                           array[i] = galaxy.IC_iso_rings_skymap[i_comp].d3[il][ib][iR].s[ip];

                        if (IC_type == "anisotropic")
                           array[i] = galaxy.IC_aniso_rings_skymap[i_comp].d3[il][ib][iR].s[ip];

                        ++i;

                     }

            status = store_mapcube_skymap(&array[0], &galaxy.E_gamma[0], galaxy.n_ISRF_components, galaxy.n_E_gammagrid, os.str(), true);

         }
      }
    }

  }

  INFO("Exit");

  return status;

}
