
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * read_isrf.cc *                                galprop package * 5/22/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include"galprop_classes.h"
#include"galprop_internal.h"
#include"fitsio.h"

#include <Units.h>
#include <PhysicalConstants.h>
#include <GalacticRadiationField.h>
#include <RadiationField.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace utl;
using namespace rf;

#include <ErrorLogger.h>

namespace {

  double BlackBodyNumberDensity(const double energy, const double kT) {   
    
    const double energyOnPi = energy/utl::kPi,
      constant = 1./(utl::kPlanckReduced/utl::s/utl::eV*utl::kSpeedOfLight_SI/utl::cm);
    
    return constant*constant*constant*
      energyOnPi*energyOnPi/(exp(energy/kT) - 1.); // eV^-1 cm^-3
    
  }

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//indexing of input ISRF array as read into 1D array

namespace {

  int isrf_index(int ir,
		 int iz,
		 int inu,
		 int icomp,
		 int nr,
		 int nz,
		 int nnu,
		 int ncomp) {
    
    return icomp*nnu*nz*nr + inu*nz*nr + iz*nr + ir;
    
  }

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
/* 
ISRF is in Hz eV cm-3 Hz-1 (or micron eV cm^-3 micron^-1)
integral  energy density Hz-1  d(nu) = integral (nu* energy density Hz-1) d(log nu)
d(log nu) is constant in this ISRF
factor=  LOG(nu(2)/nu(1)) 
*/

using namespace rf;

//void ReadISRFFormatV0(Galprop& g); // V0 format -- CMB only
//void ReadISRFFormatV1(Galprop& g); // V1 format -- FITS-based SM2000
//void ReadISRFFormatV2(Galprop& g); // V2 format -- FITS-based e.g. Porter2008
//void ReadISRFFormatV3(Galprop& g); // V3 format -- Healpix-based

static void ReadISRFFormatV0(Galprop& g) {

  INFO("Entry");

  Galaxy& galaxy = g.galaxy;
  Galdef& galdef = g.galdef;

  ostringstream buf1;
  buf1 << "Reading ISRF from no file: CMB only";
  INFO(buf1.str());

  double targetEMin = 13.6e-5, targetEMax = 13.6; // eV

  const int targetBins = 51;

  valarray<double> targetE(0., targetBins), eps1(0., targetBins), targetFreq(0., targetBins);

  for (int i = 0; i < targetBins; ++i) {

    eps1[i] = targetEMin/kElectronMass*pow(10., i*log10(targetEMax/targetEMin)/(targetBins-1));

  }

  targetE = eps1*utl::kElectronMass; // eV

  targetFreq = targetE*1./(utl::kPlanck_SI/utl::e_SI);
  
  galaxy.nu_ISRF.resize(targetFreq.size());
  galaxy.nu_ISRF = targetFreq;
  
  // Only three components: stellar + scattered, infrared, CMB
  
  const unsigned long components = 3;
  
  galaxy.n_ISRF_components = components; 
  delete[] galaxy.ISRF;
  galaxy.ISRF = new Distribution[components];
 
  ostringstream fBuf;
  fBuf << "ISRF frequency grid (Hz): ";
  
  for (int inu = 0; inu < targetBins; ++inu) 
    fBuf << galaxy.nu_ISRF[inu] << " ";

  fBuf << ends;

  INFO(fBuf.str());
    
  std::valarray<double> cmbFlux(0., targetBins);

  for (unsigned int i = 0; i < cmbFlux.size(); ++i)
    cmbFlux[i] = BlackBodyNumberDensity(targetE[i], 2.735*kBoltzmann_SI/e_SI);

  if (2 == galdef.n_spatial_dimensions) { // 2D
    
    galaxy.ISRF[0].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[1].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[2].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);

    galaxy.fTotalISRFNumberDensity.init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    
    std::ostringstream iBuf;
    iBuf << "galaxy.ISRF initialised with frequency axis dimension = " << targetBins;
    INFO(iBuf.str());
    
    for (int i = 0; i < galaxy.n_rgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_zgrid; ++j) {
	
	for (int k = 0; k < targetBins; ++k) {

	  const double factor = targetE[k]*targetE[k];

	  galaxy.ISRF[0].d2[i][j].s[k] = 0;
	  
	  galaxy.ISRF[1].d2[i][j].s[k] = 0;
	  
	  galaxy.ISRF[2].d2[i][j].s[k] = cmbFlux[k]*factor;

	  galaxy.fTotalISRFNumberDensity.d2[i][j].s[k] = cmbFlux[k];
	
	  //cout << i << " " << j << " " << r << " " << z << " " << k << " " 
	  //   << setprecision(5) << targetFreq[k] << " " 
	    //<< total*targetE[k]*targetE[k] << " " 
	    //   << stellar*targetE[k]*targetE[k] << " " 
	  // << scattered*targetE[k]*targetE[k] << " " 
	  //   << transient*targetE[k]*targetE[k] << " " 
	  //   << thermal*targetE[k]*targetE[k] << " " 
	  //   << (stellar + scattered + transient + thermal)*targetE[k]*targetE[k] << " " 
	  //   << cmbFlux[k]*targetE[k]*targetE[k] << endl;

	} // target energy bins

      } // z
      
    } // r
    
  } // 2D
    
  if (3 == galdef.n_spatial_dimensions) { // 3D
    
    galaxy.ISRF[0].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[1].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[2].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);

    galaxy.fTotalISRFNumberDensity.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    
    std::ostringstream iBuf;
    iBuf << "galaxy.ISRF initialised with frequency axis dimension = " << targetBins;
    INFO(iBuf.str());
     
#pragma omp parallel for schedule(static) default(shared)
    for (int i = 0; i < galaxy.n_xgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_ygrid; ++j) {
		
	for (int k = 0; k < galaxy.n_zgrid; ++k) {
	  
	  // Recall when retrieving the data from the skymaps it is returned
	  // in units eV cm^-2 s^-1 sr^-1
	  
	  for (int l = 0; l < targetBins; ++l) {

	    const double factor = targetE[l]*targetE[l];
	    
	    galaxy.ISRF[0].d3[i][j][k].s[l] = 0;
	    
	    galaxy.ISRF[1].d3[i][j][k].s[l] = 0;
	    
	    galaxy.ISRF[2].d3[i][j][k].s[l] = cmbFlux[l]*factor;

	    galaxy.fTotalISRFNumberDensity.d3[i][j][k].s[l] = cmbFlux[l];
	    
	  } // target energy bins

	} // z
	
      } // y
      
    } // x
    
  } // 3D
 
  INFO("Exit");

  //exit(0);

}

static void ReadISRFFormatV1(Galprop& g) {

  INFO("Entry");

  Galaxy& galaxy = g.galaxy;
  Galdef& galdef = g.galdef;

  int status = 0;

  fitsfile* fptr;
  //char ISRF_filename[200];
  int NAXIS,NAXIS1,NAXIS2,NAXIS3,NAXIS4;
  float CRVAL1,CRVAL2,CRVAL3;
  float CDELT1,CDELT2,CDELT3;
  char comment[100];
  
  //strcpy(ISRF_filename,configure.fits_directory);
  //strcat(ISRF_filename,"isrf_interp_04_000015");
  //strcat(ISRF_filename,"porter_ISRF.fits");           //AWS20050225
  //strcat(ISRF_filename,"porter_RFScattering10kpc.fits");//AWS20050301
  //strcat(ISRF_filename,galdef.ISRF_file);//AWS20050301
  
  const std::string fitsDirectory = g.configure.fFITSDataDirectory;
  const std::string isrfFilename = g.galdef.ISRF_file;
  const std::string filename = fitsDirectory + isrfFilename;
  
  //if (g.galdef.verbose>=1)cout<<" reading ISRF from "<< filename <<endl;

  std::ostringstream buf;
  buf << "Reading ISRF from " << filename;
  INFO(buf.str());

  if( fits_open_file(&fptr,filename.c_str(),READONLY,&status) ) {
    buf.str("");
    buf<<"read isrf open status= "<<status;
    WARNING(buf.str());
  }
  if (0 == fptr){
    buf.str("");
    buf<<"Cannot open file "<<filename;
    FATAL(buf.str());
    exit(-2);
  }
  
  if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) {
    buf.str("");
    buf<<"read NAXIS open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) {
    buf.str("");
    buf<<"read NAXIS1 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) {
    buf.str("");
    buf<<"read NAXIS2 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) {
    buf.str("");
    buf<<"read NAXIS3 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TINT,"NAXIS4",&NAXIS4,comment,&status) ) {
    buf.str("");
    buf<<"read NAXIS4 open status= "<<status;
    WARNING(buf.str());
  }
  
  if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) {
    buf.str("");
    buf<<"read CRVAL1 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) {
    buf.str("");
    buf<<"read CRVAL2 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TFLOAT,"CRVAL3",&CRVAL3,comment,&status) ) {
    buf.str("");
    buf<<"read CRVAL3 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) {
    buf.str("");
    buf<<"read CDELT1 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) {
    buf.str("");
    buf<<"read CDELT2 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TFLOAT,"CDELT3",&CDELT3,comment,&status) ) {
    buf.str("");
    buf<<"read CDELT3 open status= "<<status;
    WARNING(buf.str());
  }
  
  buf.str("");
  buf<<" NAXIS = "<<NAXIS <<endl;
  DEBUGLOG(buf.str());
  buf.str("");
  buf<<" NAXIS1,2,3,4 = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<" "<<NAXIS4<<endl;
  DEBUGLOG(buf.str());
  buf.str("");
  buf<<" CRVAL1,2,3 = "<<CRVAL1<<" "<<CRVAL2<<" "<<CRVAL3<<endl;
  DEBUGLOG(buf.str());
  buf.str("");
  buf<<" CDELT1,2,3 = "<<CDELT1<<" "<<CDELT2<<" "<<CDELT3<<endl;
  DEBUGLOG(buf.str());
  
  long nelements=NAXIS1*NAXIS2*NAXIS3*NAXIS4, felement=1;
  //float *isrf_in=new float[nelements];
  std::valarray<float> isrf_in(0., nelements);
  float nulval=0;
  int anynul;
    
  if (fits_read_img(fptr, TFLOAT, felement, nelements, &nulval, &isrf_in[0], 
		    &anynul, &status)){
    buf.str("");
    buf<<"read isrf open status= "<<status;
    WARNING(buf.str());
  }
  
  // for(int i=0; i<nelements; i++) cout<<isrf_in[i]<<" ";
  fits_close_file(fptr,&status);
  
  INFO("generating galaxy.ISRF:");
  
  galaxy.n_ISRF_components = NAXIS4;
  delete[] galaxy.ISRF;
  galaxy.ISRF = new Distribution[NAXIS4];
  
  if (2 == galdef.n_spatial_dimensions) {

    galaxy.fTotalISRFNumberDensity.init(galaxy.n_rgrid, galaxy.n_zgrid, NAXIS3);
    galaxy.fTotalISRFNumberDensity = 0;

  }

  if (3 == galdef.n_spatial_dimensions) {

    galaxy.fTotalISRFNumberDensity.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, NAXIS3);
    galaxy.fTotalISRFNumberDensity = 0;

  }

  // Create the array of ISRF frequencies
  // using wavelength in microns for axis 3 of input ISRF on log10 scale.
  // Reverse scale so that frequency increases.
  
  galaxy.nu_ISRF.resize(galaxy.ISRF[0].n_pgrid);
  
  // microns -> cm; nu=c/lambda
  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
    galaxy.nu_ISRF[galaxy.ISRF[0].n_pgrid-1-inu]= 
      c/(pow(10.,1.*CRVAL3+inu*CDELT3)*1.0e-4);

  std::valarray<double> targetE(0., galaxy.ISRF[0].n_pgrid);

  targetE = (kPlanck_SI/e_SI)*galaxy.nu_ISRF;
  
  DEBUGLOG(" ISRF frequency grid (Hz):");
  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++) {
     buf.str("");
    buf<<"inu galaxy.nu_ISRF[inu] "<<inu<<" "<< galaxy.nu_ISRF[inu];
    DEBUGLOG(buf.str());
  }

  for (int i = 0; i < NAXIS4; ++i) {
    
    if (2 == galdef.n_spatial_dimensions) {             // ==== 2D ====
      
      galaxy.ISRF[i].init(galaxy.n_rgrid, galaxy.n_zgrid, NAXIS3);
      buf.str("");
      buf<<" galaxy.ISRF initialized with frequency axis dimension="
	  <<galaxy.ISRF[i].n_pgrid;
      INFO(buf.str());
      
      for(int ir=0; ir<galaxy.n_rgrid; ir++) { 
	
	for(int iz=0; iz<galaxy.n_zgrid; iz++) {
	  
	  int irr=(int)((galaxy.r[ir] -CRVAL1)/CDELT1+0.5);//IMOS20060420
	  int izz=(int)((fabs(galaxy.z[iz])-CRVAL2)/CDELT2+0.5);//IMOS20060420
	  if(irr>NAXIS1-2) irr=NAXIS1-2;
	  if(izz>NAXIS2-2) izz=NAXIS2-2;
	  float rr=CRVAL1+irr*CDELT1;
	  float zz=CRVAL2+izz*CDELT2;
	  
	  // cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"     "<<irr<<" "<<izz<<"     "<<rr<<" "<<zz<<endl;
	  
	  for(int inu=0; inu<galaxy.ISRF[i].n_pgrid; inu++) {
	    float v1=isrf_in[isrf_index(irr  ,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	    float v2=isrf_in[isrf_index(irr+1,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	    float v3=isrf_in[isrf_index(irr  ,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	    float v4=isrf_in[isrf_index(irr+1,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	    float v5=v1+(v2-v1)*(galaxy.r[ir]-rr)/CDELT1;
	    float v6=v3+(v4-v3)*(galaxy.r[ir]-rr)/CDELT1;
	    float value=v5+(v6-v5)*(fabs(galaxy.z[iz])-zz)/CDELT2;
	    if(value<0.0) value=0.0;
	    // reverse scale from wavelength to frequency
	    const int nuIndex = galaxy.ISRF[i].n_pgrid-1-inu;
	    galaxy.ISRF[i].d2[ir][iz].s[nuIndex] = value;

	    galaxy.fTotalISRFNumberDensity.d2[ir][iz].s[nuIndex] += value/targetE[nuIndex]/targetE[nuIndex];
	    
	     // cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"   "<<irr<<" "<<izz<<"   "<<rr<<" "<<zz<<" "<< v5+(v6-v5)*(galaxy.z[iz]-zz)/CDELT2<<endl;
	  }  //  inu
	}  //  iz
      }  //  ir
    }  //  2D
    
    if (3 == galdef.n_spatial_dimensions) {              // ==== 3D ====
      
      galaxy.ISRF[i].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid,NAXIS3);
#pragma omp parallel for schedule(static) default(shared)
      for(int ix=0; ix<galaxy.n_xgrid; ix++) {  
	
	for(int iy=0; iy<galaxy.n_ygrid; iy++) {
	  
	  for(int iz=0; iz<galaxy.n_zgrid; iz++) {
	    
	    float r=sqrt(pow(galaxy.x[ix],2)+pow(galaxy.y[iy],2));
	    int irr=(int)((            r     -CRVAL1) /CDELT1+0.5);//IMOS20060420
	    int izz=(int)((fabs(galaxy.z[iz])-CRVAL2) /CDELT2+0.5);//IMOS20060420
	    if(irr>NAXIS1-2) irr=NAXIS1-2;
	    if(izz>NAXIS2-2) izz=NAXIS2-2;
	    
	    float rr=CRVAL1+irr*CDELT1;
	    float zz=CRVAL2+izz*CDELT2;
	    
	    // cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"     "<<irr<<" "<<izz<<"     "<<rr<<" "<<zz<<endl;
	    
	    for(int inu=0; inu<galaxy.ISRF[i].n_pgrid; inu++) {
	      
	      float v1=isrf_in[isrf_index(irr  ,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	      float v2=isrf_in[isrf_index(irr+1,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	      float v3=isrf_in[isrf_index(irr  ,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	      float v4=isrf_in[isrf_index(irr+1,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	      float v5=v1+(v2-v1)*(       r    -rr)/CDELT1;
	      float v6=v3+(v4-v3)*(       r    -rr)/CDELT1;
	      float value=v5+(v6-v5)*(fabs(galaxy.z[iz])-zz)/CDELT2;
	      if(value<0.0) value=0.0;
	      // reverse scale from wavelength to frequency

	      const int nuIndex = galaxy.ISRF[i].n_pgrid-1-inu;
	      galaxy.ISRF[i].d3[ix][iy][iz].s[nuIndex]= value;

	      galaxy.fTotalISRFNumberDensity.d3[ix][iy][iz].s[nuIndex] += value/targetE[nuIndex]/targetE[nuIndex];

	    }  //  inu   
	  }  //  ix
	}  //  iy
      }  //  iz
    }  //  3D
    
    if(galdef.verbose>=10) {
      
      buf.str("");
      buf<<"ISRF component "<<i+1;
      INFO(buf.str());
      galaxy.ISRF[i].print();
      
    }
    
  }  // isrf component i
  
  //delete[] isrf_in;//AWS20010216
  
  INFO("Exit");

}

static void ReadISRFFormatV2(Galprop& g) {

  INFO("Entry");

  Galaxy& galaxy = g.galaxy;
  Galdef& galdef = g.galdef;
  gp::Configure& configure = g.configure;

  // Read in using new format
  
  const std::string fitsDirectory = configure.fFITSDataDirectory;
  const std::string isrfFilename = galdef.ISRF_file;
  const std::string filename = fitsDirectory + isrfFilename;

  std::ostringstream buf1;
  buf1 << "Reading ISRF from " << filename;
  INFO(buf1.str());

  //galaxy.fISRFModel = new rf::GalacticRadiationField(filename, true); // For now, no angular information is read in -- TAP 09082007

  rf::GalacticRadiationField rf(filename, true); //(galaxy.fISRFModel);

  INFO("Generating galaxy.ISRF");
  std::ostringstream buf;

  // Buffer should be initialized before their str method is used to delete the string
  buf << 'a';

  // Only three components: stellar + scattered, infrared, CMB
  
  const unsigned long components = 3;
  
  galaxy.n_ISRF_components = components; 
  delete[] galaxy.ISRF;
  galaxy.ISRF = new Distribution[components];
  
  const std::valarray<double>& wl = rf.GetWavelengthData();
  
  // Cheat a bit. Get what the wavelength log(delta) is and extend the 
  // range up to 10000 microns (if it already doesn't go that far). This 
  // ensures we get the full CMB as well.
  
  const unsigned long rawWlBins = wl.size();

  std::valarray<double> wavelengthData;
  
  if (wl[rawWlBins-1] < 1e4) {
    
    const double wlDelta = log10(wl[1]/wl[0]);
    
    const unsigned long bins = long(log10(1e4/wl[0])/wlDelta) + 1;
    
    wavelengthData.resize(bins);
    
    for (unsigned long i = 0; i < bins; ++i)
      wavelengthData[i] = 
	(i < rawWlBins ? wl[i] : wl[0]*pow(10.0, i*wlDelta));
    
  } else {
    
    wavelengthData.resize(rawWlBins);
    wavelengthData = wl;
    
  }
  
  const int wlBins = wavelengthData.size();
  
  // Create the array of ISRF frequencies using wavelength in microns. 
  // Reverse scale so that frequency increases.
  
  galaxy.nu_ISRF.resize(wlBins);
  
  for (int i = 0; i < wlBins; ++i)
    galaxy.nu_ISRF[wlBins - 1 - i] = 
      utl::kSpeedOfLight_SI/(wavelengthData[i]*utl::micron/utl::m);

  std::valarray<double> targetE(0., wlBins);

  targetE = (utl::kPlanck_SI/utl::e_SI)*galaxy.nu_ISRF;
  
  DEBUGLOG(" ISRF frequency grid (Hz):");
  for(int inu=0; inu<wlBins; inu++) {
    buf.str("");
    buf<<"inu galaxy.nu_ISRF[inu] "<<inu<<" "<< galaxy.nu_ISRF[inu];
    DEBUGLOG(buf.str());
  }
  
  // Some of this is not particularly good. However, I can't fix it until
  // Galprop undergoes a bigger re-write. I can't see that happening in the
  // near future -- TAP20072301
  
  if (2 == galdef.n_spatial_dimensions) { // 2D
    
    galaxy.ISRF[0].init(galaxy.n_rgrid, galaxy.n_zgrid, wlBins);
    galaxy.ISRF[1].init(galaxy.n_rgrid, galaxy.n_zgrid, wlBins);
    galaxy.ISRF[2].init(galaxy.n_rgrid, galaxy.n_zgrid, wlBins);

    galaxy.fTotalISRFNumberDensity.init(galaxy.n_rgrid, galaxy.n_zgrid, wlBins);
    
    buf.str("");
    buf << " galaxy.ISRF initialized with frequency axis dimension = "
	 << wlBins << endl;
    INFO(buf.str());
    
    for (int i = 0; i < galaxy.n_rgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_zgrid; ++j) {
	
	const double r = galaxy.r[i], z = fabs(galaxy.z[j]);
	
	for (int k = 0; k < wlBins; ++k) {
	  
	  const double wl = wavelengthData[k];
	  
	  // Scale is reversed from wavelength to frequency
	  
	  const double stellar = 
	    rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::STELLAR);
	  
	  const double scattered = 
	    rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::SCATTERED);
	  
	  const double infrared = 
	    rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::INFRARED);
	  
	  const double cmb = 
	    rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::CMB);
	  
	  const int nuIndex = galaxy.ISRF[0].n_pgrid - 1 - k;

	  galaxy.ISRF[0].d2[i][j].s[nuIndex] = stellar + scattered;
	  
	  galaxy.ISRF[1].d2[i][j].s[nuIndex] = infrared;
	  
	  galaxy.ISRF[2].d2[i][j].s[nuIndex] = cmb;

	  galaxy.fTotalISRFNumberDensity.d2[i][j].s[nuIndex] = (stellar + scattered + infrared + cmb)/targetE[nuIndex]/targetE[nuIndex];
	    
	  //cout << i << " " << j << " " << r << " " << z << " " << galaxy.ISRF[0].n_pgrid - 1 - k << " " 
	  //   << setprecision(5) << galaxy.nu_ISRF[galaxy.ISRF[0].n_pgrid - 1 - k] << " "
	  //   << stellar << " "
	  //   << scattered << " "
	  //   << infrared << " "
	  //   << (stellar + scattered + infrared) << " " 
	  //   << cmb << endl;

	} // wl
	
      } // z
      
    } // r
    
  } // 2D
    
  if (3 == galdef.n_spatial_dimensions) { // 3D
    
    galaxy.ISRF[0].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, wlBins);
    galaxy.ISRF[1].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, wlBins);
    galaxy.ISRF[2].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, wlBins);

    galaxy.fTotalISRFNumberDensity.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, wlBins);
    
    std::ostringstream buf;
    buf<< " galaxy.ISRF initialized with frequency axis dimension = " << wlBins;
    INFO(buf.str());
    
#pragma omp parallel for schedule(static) default(shared)
    for (int i = 0; i < galaxy.n_xgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_ygrid; ++j) {
	
	const double r = sqrt(galaxy.x[i]*galaxy.x[i] + galaxy.y[j]*galaxy.y[j]);
	
	for (int k = 0; k < galaxy.n_zgrid; ++k) {
	  
	  const double z = fabs(galaxy.z[k]);
	  
	  for (int l = 0; l < wlBins; ++l) {
	    
	    const double wl = wavelengthData[l];
	    
	    // Scale is reversed from wavelength to frequency
	    
	    const double stellar = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::STELLAR);
	    
	    const double scattered = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::SCATTERED);
	    
	    const double infrared = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::INFRARED);
	    
	    const double cmb = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::CMB);
	    
	    const int nuIndex = galaxy.ISRF[0].n_pgrid - 1 - l;

	    galaxy.ISRF[0].d3[i][j][k].s[nuIndex] = stellar + scattered;
	    
	    galaxy.ISRF[1].d3[i][j][k].s[nuIndex] = infrared;
	    
	    galaxy.ISRF[2].d3[i][j][k].s[nuIndex] = cmb;

	    galaxy.fTotalISRFNumberDensity.d3[i][j][k].s[nuIndex] = (stellar + scattered + infrared + cmb)/targetE[nuIndex]/targetE[nuIndex];
	    
	  } // wl
	  
	} // z
	
      } // y
      
    } // x
    
  } // 3D
 
  INFO("Exit");

}

static void ReadISRFFormatV3(Galprop& g) {

  INFO("Entry");

  Galaxy& galaxy = g.galaxy;
  Galdef& galdef = g.galdef;
  gp::Configure& configure = g.configure;

  const std::string fitsDirectory = configure.fFITSDataDirectory;
  const std::string isrfFilename = galdef.ISRF_file;
  const std::string filename = fitsDirectory + isrfFilename;

  std::ostringstream buf1;
  buf1 << "Reading ISRF from " << filename;
  INFO(buf1.str());

  if (g.galdef.verbose >= 1) {

    std::ostringstream addInfoBuf;
    addInfoBuf << "FITS directory: " << fitsDirectory << endl 
	       << "ISRF filename: " << isrfFilename << endl 
	       << "ISRF healpix order: " << galdef.ISRF_healpixOrder;
    INFO(addInfoBuf.str());

  }

  const double targetEMin = 13.6e-5, targetEMax = 13.6; // eV

  const size_t targetBins = 101;

  std::valarray<double> targetE(0., targetBins), eps1(0., targetBins), targetFreq(0., targetBins);

  for (int i = 0; i < targetBins; ++i) {

    eps1[i] = targetEMin/utl::kElectronMass*pow(10., i*log10(targetEMax/targetEMin)/(targetBins-1));

  }

  targetE = eps1*utl::kElectronMass; // eV

  targetFreq = targetE*1./(utl::kPlanck_SI/utl::e_SI);

  if ( filename != galaxy.fISRFloadedfile ) {
     delete galaxy.fISRFrf;
     galaxy.fISRFrf = new rf::RadiationField(filename, targetFreq, galdef.ISRF_healpixOrder);
     galaxy.fISRFloadedfile = filename;
  }

  rf::RadiationField& rf = *galaxy.fISRFrf;
  
  galaxy.nu_ISRF.resize(targetFreq.size());
  galaxy.nu_ISRF = targetFreq;
  
  // Only three components: stellar + scattered, infrared, CMB
  
  const unsigned int components = 3;
  
  galaxy.n_ISRF_components = components; 
  delete[] galaxy.ISRF;
  galaxy.ISRF = new Distribution[components];
 
  std::ostringstream fBuf;
  fBuf << "ISRF frequency grid (Hz): ";
  
  for (int inu = 0; inu < targetBins; ++inu) 
    fBuf << galaxy.nu_ISRF[inu] << " ";

  fBuf << ends;

  INFO(fBuf.str());
    
  std::valarray<double> cmbFlux(0., targetBins);

  for (unsigned long i = 0; i < cmbFlux.size(); ++i)
    cmbFlux[i] = BlackBodyNumberDensity(targetE[i], 2.735*utl::kBoltzmann_SI/utl::e_SI);

  if (2 == galdef.n_spatial_dimensions) { // 2D
    
    galaxy.ISRF[0].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[1].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[2].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);

    galaxy.fTotalISRFNumberDensity.init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    
    std::ostringstream iBuf;
    iBuf << "galaxy.ISRF initialised with frequency axis dimension = " << targetBins;
    INFO(iBuf.str());
    
    for (int i = 0; i < galaxy.n_rgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_zgrid; ++j) {
	
	const double r = galaxy.r[i], z = galaxy.z[j];

	rf::RadiationField::ThreeVector pos(r, 0, z);

	// Recall when retrieving the data from the skymaps it is returned
	// in units eV^-1 cm^-3 sr^-1 -- simply sum over all pixels and 
	// multiply by solid angle to get number density

	//const Skymap<double> skymapTotal = rf.GetSkymap(pos, RadiationField::TOTAL, galdef.ISRF_healpixOrder);
	
	//const Skymap<double> skymapDirect = rf.GetSkymap(pos, RadiationField::DIRECT, galdef.ISRF_healpixOrder);

	//const Skymap<double> skymapScattered = rf.GetSkymap(pos, RadiationField::SCATTERED, galdef.ISRF_healpixOrder);

	//const Skymap<double> skymapTransient = rf.GetSkymap(pos, RadiationField::TRANSIENT, galdef.ISRF_healpixOrder);

	//const Skymap<double> skymapThermal = rf.GetSkymap(pos, RadiationField::THERMAL, galdef.ISRF_healpixOrder);

	valarray<double> numberDensityDirect(0., targetBins);
	valarray<double> numberDensityScattered(0., targetBins);
	valarray<double> numberDensityTransient(0., targetBins);
	valarray<double> numberDensityThermal(0., targetBins);

	numberDensityDirect = rf.GetNumberDensity(pos, RadiationField::DIRECT);
	numberDensityScattered = rf.GetNumberDensity(pos, RadiationField::SCATTERED);
	numberDensityTransient = rf.GetNumberDensity(pos, RadiationField::TRANSIENT);
	numberDensityThermal = rf.GetNumberDensity(pos, RadiationField::THERMAL);
	  
	for (int k = 0; k < targetBins; ++k) {

	  const double factor = targetE[k]*targetE[k];

	  //const double total = skymapTotal.sum(k)*skymapTotal.solidAngle();

	  const double stellar = numberDensityDirect[k];//skymapDirect.sum(k)*skymapDirect.solidAngle();
	  const double scattered = numberDensityScattered[k];//skymapScattered.sum(k)*skymapScattered.solidAngle();
	  const double transient = numberDensityTransient[k];//skymapTransient.sum(k)*skymapTransient.solidAngle();
	  const double thermal = numberDensityThermal[k];//skymapThermal.sum(k)*skymapThermal.solidAngle();

	  galaxy.ISRF[0].d2[i][j].s[k] = (stellar + scattered)*factor;
	  
	  galaxy.ISRF[1].d2[i][j].s[k] = (transient + thermal)*factor;
	  
	  galaxy.ISRF[2].d2[i][j].s[k] = cmbFlux[k]*factor;

	  galaxy.fTotalISRFNumberDensity.d2[i][j].s[k] = stellar + scattered + transient + thermal + cmbFlux[k];

	  //cout << i << " " << j << " " << r << " " << z << " " << k << " " << galaxy.ISRF[0].d2[i][j].s[k] << " " << galaxy.ISRF[1].d2[i][j].s[k] << " " << galaxy.ISRF[2].d2[i][j].s[k] << endl;


	} // target energy bins

	//exit(0);

      } // z
      
    } // r
    
  } // 2D

  //exit(0);
    
  if (3 == galdef.n_spatial_dimensions) { // 3D
    
    galaxy.ISRF[0].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[1].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[2].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    
    galaxy.fTotalISRFNumberDensity.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);

    std::ostringstream iBuf;
    iBuf << "galaxy.ISRF initialised with frequency axis dimension = " << targetBins;
    INFO(iBuf.str());

#pragma omp parallel for schedule (static)     
    for (int ij = 0; ij < galaxy.n_xgrid*galaxy.n_ygrid; ++ij) {

      const int i = ij/galaxy.n_ygrid;
      const int j = ij%galaxy.n_ygrid;

      const double x = galaxy.x[i], y = galaxy.y[j];
		
      //std::cout << ij << " " << i << " " << j << " " << x << " " << y << std::endl;

      for (int k = 0; k < galaxy.n_zgrid; ++k) {
	  
	const double z = galaxy.z[k];
	
	rf::RadiationField::ThreeVector pos(x, y, z);
	
	// Recall when retrieving the data from the skymaps it is returned
	// in units eV cm^-2 s^-1 sr^-1
	
	//const Skymap<double> skymapTotal = rf.GetSkymap(pos, RadiationField::TOTAL, galdef.ISRF_healpixOrder);
	
	//const Skymap<double> skymapDirect = rf.GetSkymap(pos, RadiationField::DIRECT, galdef.ISRF_healpixOrder);
	
	//const Skymap<double> skymapScattered = rf.GetSkymap(pos, RadiationField::SCATTERED, galdef.ISRF_healpixOrder);
	
	//const Skymap<double> skymapTransient = rf.GetSkymap(pos, RadiationField::TRANSIENT, galdef.ISRF_healpixOrder);
	
	//const Skymap<double> skymapThermal = rf.GetSkymap(pos, RadiationField::THERMAL, galdef.ISRF_healpixOrder);
	
	/*std::valarray<double> numberDensityDirect(0., targetBins);
	std::valarray<double> numberDensityScattered(0., targetBins);
	std::valarray<double> numberDensityTransient(0., targetBins);
	std::valarray<double> numberDensityThermal(0., targetBins);
	
	numberDensityDirect = rf.GetNumberDensity(pos, RadiationField::DIRECT);
	numberDensityScattered = rf.GetNumberDensity(pos, RadiationField::SCATTERED);
	numberDensityTransient = rf.GetNumberDensity(pos, RadiationField::TRANSIENT);
	numberDensityThermal = rf.GetNumberDensity(pos, RadiationField::THERMAL);
	*/

	std::valarray<double> numberdensityoptical(0., targetBins), numberdensityinfrared(0., targetBins);

	if (rf.Is3D()) {

	  numberdensityoptical = rf.GetNumberDensity(pos, RadiationField::OPTICAL);
	  numberdensityinfrared = rf.GetNumberDensity(pos, RadiationField::INFRARED);

	} else {

	  std::valarray<double> numberdensitydirect(0., targetBins), numberdensityscattered(0., targetBins), numberdensitytransient(0., targetBins), numberdensitythermal(0., targetBins);

	  numberdensitydirect = rf.GetNumberDensity(pos, RadiationField::DIRECT);
	  numberdensityscattered = rf.GetNumberDensity(pos, RadiationField::SCATTERED);
	  numberdensitytransient = rf.GetNumberDensity(pos, RadiationField::TRANSIENT);
	  numberdensitythermal = rf.GetNumberDensity(pos, RadiationField::THERMAL);
	  numberdensityoptical = numberdensitydirect + numberdensityscattered;
	  numberdensityinfrared = numberdensitytransient + numberdensitythermal;

	}

	for (int l = 0; l < targetBins; ++l) {
	  
	  const double factor = targetE[l]*targetE[l];
	  
	  //const double total = skymapTotal.sum(l)*skymapTotal.solidAngle();
	  
	  //const double stellar = numberDensityDirect[l];//skymapDirect.sum(l)*skymapDirect.solidAngle();
	  //const double scattered = numberDensityScattered[l];//skymapScattered.sum(l)*skymapScattered.solidAngle();
	  //const double transient = numberDensityTransient[l];//skymapTransient.sum(l)*skymapTransient.solidAngle();
	  //const double thermal = numberDensityThermal[l];//skymapThermal.sum(l)*skymapThermal.solidAngle();
          //const double optical = numberDensityDirect[l] + numberDensityScattered[l];
          //const double infrared = numberDensityTransient[l] + numberDensityThermal[l];
	  
#pragma omp atomic write
	  galaxy.ISRF[0].d3[i][j][k].s[l] = numberdensityoptical[l]*factor;
	  
#pragma omp atomic write
	  galaxy.ISRF[1].d3[i][j][k].s[l] = numberdensityinfrared[l]*factor;
	  
#pragma omp atomic write
	  galaxy.ISRF[2].d3[i][j][k].s[l] = cmbFlux[l]*factor;
	  
#pragma omp atomic write
	  galaxy.fTotalISRFNumberDensity.d3[i][j][k].s[l] = numberdensityoptical[l] + numberdensityinfrared[l] + cmbFlux[l];

	  /*#pragma omp critical (printoutisrf)
	  {	  
	    if ((x == 0 && y == 0 && z == 0) || (x==15 && y==15 && z==-6))
	      std::cout << i << " " << j << " " << k << " " << x << " " << y << " " << z << " " << kSpeedOfLight_SI/targetFreq[l]*(m/micron) << " " << galaxy.ISRF[0].d3[i][j][k].s[l] << " " << galaxy.ISRF[1].d3[i][j][k].s[l] << " " << galaxy.fTotalISRFNumberDensity.d3[i][j][k].s[l]*factor << std::endl;
	  }
	  */
	} // target energy bins
	
      } // z
        
    } // xy
    
  } // 3D
 
  //exit(0);

  INFO("Exit");

  /*for (auto ij = 0; ij < galaxy.n_xgrid*galaxy.n_ygrid; ++ij) {

    auto i = ij/galaxy.n_ygrid;
    auto j = ij%galaxy.n_ygrid;
    
    auto x = galaxy.x[i], y = galaxy.y[j];

    auto sumOptical(0.), sumIR(0.);

    // ISRF is in eV cm^-3 micron^-1 micron

    for (auto k = 0; k < targetBins; ++k) {

      sumOptical += galaxy.ISRF[0].d3[i][j][2].s[k];
      sumIR += galaxy.ISRF[1].d3[i][j][2].s[k];
	  	  
    }

    std::cout << x << " " << y << " " << sumOptical*log(targetE[1]/targetE[0]) << " " << sumIR*log(targetE[1]/targetE[0]) << " " << log(targetE[1]/targetE[0]) << std::endl;

  }

  exit(0);
  */
}

int Galprop::read_isrf(const int version) {

  INFO("Entry");

  assert (version >= 0 && version <= 3); // Update later for other versions -- make sure you have the right form of the ISRF file structure before running!

  int status = 0;

  if (0 == version)
    ReadISRFFormatV0(*this);

  if (1 == version)
    ReadISRFFormatV1(*this);
 
  if (2 == version)
    ReadISRFFormatV2(*this);

  if (3 == version)
    ReadISRFFormatV3(*this);

  for (int i = 0; i < galaxy.n_ISRF_components; ++i) {

    if (2 == galdef.n_spatial_dimensions) {

      for (int iR = 0; iR < galaxy.ISRF->n_rgrid; ++iR)
	for (int iZ = 0; iZ < galaxy.ISRF->n_zgrid; ++iZ) 	  
	  for (int iP = 0; iP < galaxy.ISRF->n_pgrid; ++iP)
	    galaxy.ISRF[i].d2[iR][iZ].s[iP] *= galaxy.fISRFFactors[i]; 

    }

    if (3 == galdef.n_spatial_dimensions) {

      for (int iX = 0; iX < galaxy.ISRF->n_xgrid; ++iX)
	for (int iY = 0; iY < galaxy.ISRF->n_ygrid; ++iY)
	  for (int iZ = 0; iZ < galaxy.ISRF->n_zgrid; ++iZ)
	    for (int iP = 0; iP < galaxy.ISRF->n_pgrid; ++iP)
	      galaxy.ISRF[i].d3[iX][iY][iZ].s[iP] *= galaxy.fISRFFactors[i]; 

    }

  }
      
  INFO("Exit");

  return status;

}
