#include <cassert>
#include <string>
#include <cstring>
#include <valarray>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"
#include "fitsio.h"

#include <ErrorLogger.h>
#include <BaseSkyFitsIO.h>

int Galprop::store_bremss_HIR_skymap() {

  INFO("Entry");

  int status = 0;

  if (3 == galdef.skymap_format) {
    
    char index[3];
    string fileprefix = configure.fOutputDirectory + configure.fOutputPrefix;
    fileprefix += "bremss_HIR_ring_";

    for (int i_Ring=0; i_Ring<galaxy.n_Ring; i_Ring++){
      sprintf(index, "%d", i_Ring+1);
      string filename = fileprefix + index;
      filename += "_healpix_";
      filename += galdef.galdef_ID;
      filename += ".gz";
      SM::writeToFits(*galaxy.bremss_HIR_hp_skymap[i_Ring], filename, true, true, "Energy", "MeV");
    }

  }else{
    
    //fitsfile *fptr;       // pointer to the FITS file; defined in fitsio.h
    int status, ii, jj;
    long nelements;
    long naxes[4]; 
    double crval[4], cdelt[4];
    
    naxes[0]=galaxy.n_long;
    naxes[1]=galaxy.n_lat;   
    naxes[2]=galaxy.bremss_HIR_skymap.n_zgrid; // number of Galactocentric rings
    naxes[3]=galaxy.n_E_gammagrid;
    
    nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];

    valarray<float> array(0., nelements);
        
    if (galdef.skymap_format != 1){ //!< Old output
      int i=0; 
      for     (int ip       =0;            ip <naxes[3];        ip++)
	for    (int i_Ring   =0;        i_Ring <naxes[2];    i_Ring++)
	  for   (int ib       =0;            ib <naxes[1];        ib++)
	    for (int il       =0;            il <naxes[0];        il++)
	      {
		array[i]=0.0;
		array[i]+=galaxy.bremss_HIR_skymap .d3[il][ib][i_Ring].s[ip] *pow(galaxy.E_gamma[ip],2);
		i++;
	      }
      
      crval[0]=galaxy.long_min;
      crval[1]=galaxy. lat_min;
      crval[2]=0; //IMOS20080114
      crval[3]=log10(galaxy.E_gamma_min);
      
      cdelt[0]=galaxy.d_long;
      cdelt[1]=galaxy.d_lat;
      cdelt[2]=1;
      cdelt[3]=log10(galaxy.E_gamma_factor);
      
      //Use the standard method to store the skymap
      status = store_skymap(&array[0], naxes, "bremss_HIR_skymap_", crval, cdelt);
      
    }

    if(galdef.skymap_format == 1 || galdef.skymap_format == 2){ //!< Mapcube output compatable with Glast science tools
      int i=0; 
      for       (int i_Ring   =0;        i_Ring <naxes[2];    i_Ring++)
        for     (int ip       =0;            ip <naxes[3];        ip++) // IMOS20080114
	  for   (int ib       =0;            ib <naxes[1];        ib++)
	    for (int il       =0;            il <naxes[0];        il++)
	      {
		array[i]=galaxy.bremss_HIR_skymap .d3[il][ib][i_Ring].s[ip];
		i++;
	      }
      
      status = store_mapcube_skymap(&array[0], &galaxy.E_gamma[0], galaxy.bremss_HIR_skymap.n_zgrid, galaxy.n_E_gammagrid, "bremss_HIR_mapcube_", true);
    }

  }

  INFO("Exit");

  return status;

}
