//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_synch_skymap.cc *                       galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <cstring>
#include <string>
#include <sstream>
#include <valarray>

#include "galprop_classes.h"
#include "galprop_internal.h"

using namespace std;

#include "fitsio.h" 

#include <ErrorLogger.h>
#include <BaseSkyFitsIO.h>

int Galprop::store_synch_skymap() { //AWS20050817

  INFO("Entry");

  int status = 0;

  if (3 == galdef.skymap_format) {

    string filename;                                                                                                                                //AWS20100709
    filename =       configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_healpix_" + galdef.galdef_ID + ".gz";                      //AWS20100107
    SM::writeToFits(*galaxy.synchrotron_hp_skymap, filename, true, true, "Frequency", "Hz");

    filename =       configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_Q_healpix_" + galdef.galdef_ID + ".gz";                    //AWS20100709
    SM::writeToFits(*galaxy.synchrotron_Q_hp_skymap, filename, true, true, "Frequency", "Hz");
    
    filename =       configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_U_healpix_" + galdef.galdef_ID + ".gz";             //AWS20100709
    SM::writeToFits(*galaxy.synchrotron_U_hp_skymap, filename, true, true, "Frequency", "Hz");

    filename =       configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_P_healpix_" + galdef.galdef_ID + ".gz";             //AWS20110328
    SM::writeToFits(*galaxy.synchrotron_P_hp_skymap, filename, true, true, "Frequency", "Hz");

    filename =       configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_polang_healpix_" + galdef.galdef_ID + ".gz";             //AWS20110919
    SM::writeToFits(*galaxy.synchrotron_polang_hp_skymap, filename, true, true, "Frequency", "Hz");

    filename =       configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_polfra_healpix_" + galdef.galdef_ID + ".gz";             //AWS20110922
    SM::writeToFits(*galaxy.synchrotron_polfra_hp_skymap, filename, true, true, "Frequency", "Hz");

    filename =       configure.fOutputDirectory + configure.fOutputPrefix + "free_free_healpix_" + galdef.galdef_ID + ".gz";             //AWS20110906
    if ( galdef.free_free_absorption >= 1 ) 
       SM::writeToFits(*galaxy.free_free_hp_skymap, filename, true, true, "Frequency", "Hz");


  } else {

    long naxes[4]; 
    double crval[4], cdelt[4];
        
    naxes[0] = galaxy.n_long;
    naxes[1] = galaxy.n_lat;             
    naxes[2] = galaxy.n_nu_synchgrid;
    naxes[3] = 1;
    
    const long nElements = naxes[0]*naxes[1]*naxes[2]*naxes[3];
    
    valarray<float> array      (0., nElements);
    valarray<float> arrayQ     (0., nElements);//AWS20100709
    valarray<float> arrayU     (0., nElements);//AWS20100709
    valarray<float> arrayP     (0., nElements);//AWS20110328
    valarray<float> arraypolang(0., nElements);//AWS20110919
    valarray<float> arraypolfra(0., nElements);//AWS20110922

    valarray<float> array_free_free(0., nElements);//AWS20110906

    //float *array;          
    //array=new float[nelements];
  
    int i = 0;
  
    for (int ip = 0; ip < naxes[2]; ++ip) {

      for (int ib = 0; ib < naxes[1]; ++ib) {
	
	for (int il = 0; il < naxes[0]; ++il) {
	  
	  array [i]      = galaxy.synchrotron_skymap       .d2[il][ib].s[ip];
	  arrayQ[i]      = galaxy.synchrotron_Q_skymap     .d2[il][ib].s[ip];//AWS20100709
	  arrayU[i]      = galaxy.synchrotron_U_skymap     .d2[il][ib].s[ip];//AWS20100709
	  arrayP[i]      = galaxy.synchrotron_P_skymap     .d2[il][ib].s[ip];//AWS20110328
	  arraypolang[i] = galaxy.synchrotron_polang_skymap.d2[il][ib].s[ip];//AWS20110919
	  arraypolfra[i] = galaxy.synchrotron_polfra_skymap.d2[il][ib].s[ip];//AWS20110922

	  array_free_free[i] = galaxy.free_free_skymap.d2[il][ib].s[ip];//AWS20110906

	  ++i;
	  
	}
	
      }
      
    }

    if (galdef.skymap_format != 1) {

      crval[0] = galaxy.long_min;
      crval[1] = galaxy. lat_min;
      crval[2] = log10(galaxy.nu_synch_min);
      crval[3] = 1;

      cdelt[0] = galaxy.d_long;
      cdelt[1] = galaxy.d_lat;
      cdelt[2] = log10(galaxy.nu_synch_factor);
      cdelt[3] = 1;

      //Use the standard method to store the skymap
      status = store_skymap(&array      [0], naxes, "synchrotron_skymap_",        crval, cdelt);
      status = store_skymap(&arrayQ     [0], naxes, "synchrotron_Q_skymap_",      crval, cdelt);//AWS20100709
      status = store_skymap(&arrayU     [0], naxes, "synchrotron_U_skymap_",      crval, cdelt);//AWS20100709
      status = store_skymap(&arrayP     [0], naxes, "synchrotron_P_skymap_",      crval, cdelt);//AWS20110328
      status = store_skymap(&arraypolang[0], naxes, "synchrotron_polang_skymap_", crval, cdelt);//AWS20110919
      status = store_skymap(&arraypolfra[0], naxes, "synchrotron_polfra_skymap_", crval, cdelt);//AWS20110922

      status = store_skymap(&array_free_free[0], naxes, "free_free_skymap_",      crval, cdelt);//AWS20110906

    }

    if ((1 == galdef.skymap_format) || 
	(2 == galdef.skymap_format)) { //!< Mapcube output compatable with Glast science tools

      status = store_mapcube_skymap(&array [0], &galaxy.nu_synch[0], 1, galaxy.n_nu_synchgrid, "synchrotron_mapcube_",   false);
      status = store_mapcube_skymap(&arrayQ[0], &galaxy.nu_synch[0], 1, galaxy.n_nu_synchgrid, "synchrotron_Q_mapcube_", false);
      status = store_mapcube_skymap(&arrayU[0], &galaxy.nu_synch[0], 1, galaxy.n_nu_synchgrid, "synchrotron_U_mapcube_", false);
      status = store_mapcube_skymap(&arrayP[0], &galaxy.nu_synch[0], 1, galaxy.n_nu_synchgrid, "synchrotron_P_mapcube_", false);//AWS20110328
      status = store_mapcube_skymap(&arraypolang[0], &galaxy.nu_synch[0], 1, galaxy.n_nu_synchgrid, "synchrotron_polang_mapcube_", false);//AWS20110919

      status = store_mapcube_skymap(&arraypolfra[0], &galaxy.nu_synch[0], 1, galaxy.n_nu_synchgrid, "synchrotron_polfra_mapcube_", false);//AWS20110919

      status = store_mapcube_skymap(&array_free_free[0], &galaxy.nu_synch[0], 1, galaxy.n_nu_synchgrid, "free_free_mapcube_", false);//AWS20110906

    }

  }

  INFO("Exit");

  return status;

}
