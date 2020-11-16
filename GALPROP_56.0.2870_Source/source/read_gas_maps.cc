
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * read_gas_maps.cc *                             galprop package * 1/14/2008 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
//written by IMOS20080114

using namespace std;
#include"galprop_classes.h"
#include"galprop_internal.h"
#include"fitsio.h"
#include"Skymap.h"
#include <vector>
#include <algorithm>
#include <cstring>

#include <string>
#include <constants.h>
#include <cassert>

#include <ErrorLogger.h>
#include <FullSky.h>
#include <SparseSky.h>
#include <BaseSkyFitsIO.h>

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Calculation of the ring number
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
int Galprop::gas_iRing(double R)
{
  int i_Ring;
  for(i_Ring=0; i_Ring<galaxy.n_Ring-1; i_Ring++) 
    if( R <= galaxy.R_bins[galaxy.n_Ring+i_Ring] ) break;
  return i_Ring;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// functional dependence of X_CO
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double Galprop::fX_CO(double R) const //IMOS20080114
{
  double fX_CO_=galdef.X_CO;
  double a,b,Rcut;

  switch(galdef.n_X_CO)
    {
    case 0:
      return fX_CO_;

			//Tabulated values
		case 1:
      if ( R < galdef.X_CO_radius[0] ) {
				fX_CO_ = galdef.X_CO_values[0];
			} else if ( R >= galdef.X_CO_radius[galdef.X_CO_values.size()-1]) {
				fX_CO_ = galdef.X_CO_values[galdef.X_CO_values.size()-1];
			} else {
				int i = 0;
				while (R > galdef.X_CO_radius[i]) 
					i++;
				fX_CO_ = galdef.X_CO_values[i-1] + (galdef.X_CO_values[i]-galdef.X_CO_values[i-1]) /
					                                 (galdef.X_CO_radius[i]-galdef.X_CO_radius[i-1]) *
																				   (R-galdef.X_CO_radius[i-1]);
			}
			return fX_CO_;
			//Exponential function with a constant and a linear term
		case 2:
			{
			double X0 = galdef.X_CO_parameters[0];
			double A = galdef.X_CO_parameters[1];
			double B = galdef.X_CO_parameters[2];
			double C = galdef.X_CO_parameters[3];
			return X0 + A*R + B*pow(10.0, C*R);
			}

			//Tabulated values with power law interpolation if possible
			//No extrapolation is performed, last value used as a constant
			//Linear interpolation is done if one of the values or radius is 0
		case 3:
      if ( R < galdef.X_CO_radius[0] ) {
				fX_CO_ = galdef.X_CO_values[0];
			} else if ( R >= galdef.X_CO_radius[galdef.X_CO_values.size()-1]) {
				fX_CO_ = galdef.X_CO_values[galdef.X_CO_values.size()-1];
			} else {
				int i = 0;
				while (R > galdef.X_CO_radius[i]) 
					i++;
				if ( galdef.X_CO_values[i-1] != 0 && galdef.X_CO_values[i] != 0 && galdef.X_CO_radius[i] != 0 && galdef.X_CO_radius[i-1] != 0) {
					double index = log(galdef.X_CO_values[i]/galdef.X_CO_values[i-1])/log(galdef.X_CO_radius[i]/galdef.X_CO_radius[i-1]);
					fX_CO_ = galdef.X_CO_values[i]*pow(R/galdef.X_CO_radius[i],index);
				} else {
					fX_CO_ = galdef.X_CO_values[i-1] + (galdef.X_CO_values[i]-galdef.X_CO_values[i-1]) /
					                                   (galdef.X_CO_radius[i]-galdef.X_CO_radius[i-1]) *
									  											   (R-galdef.X_CO_radius[i-1]);
				}
			}
			return fX_CO_;

			//Strong et al. 2004
    case 9: 
      fX_CO_=0.4E20;
      if (R> 3.5) fX_CO_= 0.6E20;
      if (R> 5.5) fX_CO_= 0.8E20;
      if (R> 7.5) fX_CO_= 1.5E20;
      if (R> 9.5) fX_CO_=10.0E20;
      return fX_CO_;

                               
      // an exponential function: AWS20090615
    case 10:
      a = -0.4;  b = 0.066; Rcut = 15.0;
      if(R< Rcut)  fX_CO_ = 1.0E20 * pow(10.0, a + b * R   );
      if(R>=Rcut)  fX_CO_ = 1.0E20 * pow(10.0, a + b * Rcut);
      return fX_CO_;
      
    default: return fX_CO_;
    }  
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// read_gas_maps
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
int Galprop::read_gas_maps(const std::string &type) {

  ostringstream buf;
  buf << "Entry: " << type;
  INFO(buf.str());
  int status=0;

  //Make the re-normalization with mapcubes optional
  if ( galdef.HIR_filename == "" || galdef.COR_filename == "" ) {
     INFO("HIR or COR filename is empty, gasmaps will not be used to renormalize the gamma-ray skymaps");
     galaxy.renormGas = false;

     //Have a single ring that encompasses the entire Galaxy
     //It could be possible in the future to specify binning in this case as well.
     galaxy.n_Ring = 1;
     galaxy.R_bins.resize(2);
     galaxy.R_bins[0] = 0;
     if (galdef.n_spatial_dimensions == 2)
        galaxy.R_bins[1] = galdef.r_max;
     else
        //Symmetry is not guaranteed
        galaxy.R_bins[1] = sqrt(galdef.x_max*galdef.x_max + galdef.y_max*galdef.y_max);
        galaxy.R_bins[1] = galaxy.R_bins[1] < sqrt(galdef.x_max*galdef.x_max + galdef.y_min*galdef.y_min) ? sqrt(galdef.x_max*galdef.x_max + galdef.y_min*galdef.y_min) : galaxy.R_bins[1];
        galaxy.R_bins[1] = galaxy.R_bins[1] < sqrt(galdef.x_min*galdef.x_min + galdef.y_min*galdef.y_min) ? sqrt(galdef.x_min*galdef.x_min + galdef.y_min*galdef.y_min) : galaxy.R_bins[1];
        galaxy.R_bins[1] = galaxy.R_bins[1] < sqrt(galdef.x_min*galdef.x_min + galdef.y_max*galdef.y_max) ? sqrt(galdef.x_min*galdef.x_min + galdef.y_max*galdef.y_max) : galaxy.R_bins[1];
  
  } else {

  galaxy.renormGas = true;

  fitsfile *fptr;
  
  int NAXIS,NAXIS1,NAXIS2,NAXIS3,NAXIS4;
  float CRVAL1,CRVAL2,CRVAL3;
  float CDELT1,CDELT2,CDELT3;
  char comment[100];//, filename[200];
  
  //strcpy(filename,configure.fFITSDataDirectory);

  const std::string fitsDirectory = configure.fFITSDataDirectory;

  std::string filename = fitsDirectory;
  
  if( type == "HIR" ) filename += galdef.HIR_filename;//strcat(filename,galdef.HIR_filename);
  if( type == "COR" ) filename += galdef.COR_filename;//strcat(filename,galdef.COR_filename);
  
  Distribution nuse;      // number of cells used in rebinned map
  
  //   cout<<"galaxy.n_long,n_lat  "<<galaxy.n_long<<" "<<galaxy.n_lat<<endl;
  
  ostringstream ost;
  ost<<" reading "<<type<<" from "<<filename;
  INFO(ost.str());

  
  fits_open_file(&fptr,filename.c_str(),READONLY,&status);
  if( status ) {
     ost.str("");
     ost<<"read "<<type<<" open status= "<<status;
     ERROR(ost.str());
     throw std::runtime_error("Failed to open skymap \"" + filename + "\".  Check your paths and galdef file.");
  }
  
  
  // for(int i=0; i<nelements; i++) cout<<image[i]<<" ";
 
  ost.str("");
  ost<<"generating galaxy."<<type<<":";
  INFO(ost.str());
  
  // HEALPIX format 
  if (galdef.skymap_format == 3 || 4 == galdef.skymap_format)
    {
       std::vector<double> Rmin, Rmax;

       // Read the old style files into a FullSky map and then convert to SparseSky when needed
       std::unique_ptr< SM::BaseSky<double> > input;

       try {

          INFO("Trying skymap format");

          //In case we have an healpix input 
          input.reset( new SM::FullSky<double>(SM::SpectralBinning({1}, {2}), galdef.healpix_order, RING, SM::CoordSys::GAL, 0) );
          std::map<std::string, std::string> keywords;
          SM::readFromFits(input, filename, keywords);

          //We need to make sure the input is of the right order
          input = input->ResetOrder( galdef.healpix_order, false );

          //Replace the ring boundaries and assert that it is binned
          Rmin = input->GetBinning().GetBinsMin();
          Rmax = input->GetBinning().GetBinsMax();

          assert(Rmin.size() > 0);

          galaxy.n_Ring = Rmin.size();
          galaxy.R_bins.resize(2*galaxy.n_Ring);

          for (size_t i = 0; i < Rmin.size(); ++i) {
             galaxy.R_bins[i] = Rmin[i];
             galaxy.R_bins[i+Rmin.size()] = Rmax[i];
          }

          INFO("Success");

       } catch ( ... ) {

          INFO("Reverting to WCS mapcube format");

          //Get the ring information from the file
          if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) cerr<<"3read "<<type<<" status= "<<status<<endl;
          galaxy.n_Ring=NAXIS3;

          ost.str("");
          ost<<"Reading "<<galaxy.n_Ring<<" rings";
          INFO(ost.str());

          //read an extension (2nd header): a table of ring radii
          int HDU_TYPE;
          float nulval=1e-15;
          int anynul;
          if( fits_movabs_hdu(fptr,2,&HDU_TYPE,&status) ) cerr<<"8read "<<type<<" status= "<<status<<endl;

          //galaxy.R_bins=new float[2*galaxy.n_Ring]; // two columns: Rmin, Rmax
          galaxy.R_bins.resize(2*galaxy.n_Ring);

          long long int FIRSTROW = 1, FIRSTELEM = 1, NELEMENTS = galaxy.n_Ring;

          float* rBins = new float[2*galaxy.n_Ring];

          for(int i=0;i<2;i++) // read each column
             if( fits_read_col(fptr,TFLOAT,i+1,FIRSTROW,FIRSTELEM,NELEMENTS,&nulval,&rBins[i*galaxy.n_Ring],&anynul,&status) ) 
                cout<<"12read "<<type<<" status= "<<status<<endl;

          for (size_t i = 0; i < galaxy.R_bins.size(); ++i) 
             galaxy.R_bins[i] = double(rBins[i]);

          delete[] rBins;

          Rmin.resize(galaxy.n_Ring), Rmax.resize(galaxy.n_Ring);
          for ( int i = 0; i < galaxy.n_Ring; ++i ) {
             Rmin[i] = galaxy.R_bins[i];
             Rmax[i] = galaxy.R_bins[i+galaxy.n_Ring];
             ost.str("");
             ost<<"Ring boundary "<<i<<": "<<Rmin[i]<<" - "<<Rmax[i];
             INFO(ost.str());
          }

          input.reset( new SM::FullSky<double>(SM::SpectralBinning(Rmin, Rmax), galdef.healpix_order, RING, SM::CoordSys::GAL, 0));
          //Assume it is WCS format
          SM::readFromFitsWcs(input, filename, false);
       }

       //Now create the pointers
       if( type == "HIR" )
       {
          galaxy.hpHIR.resize( galaxy.n_Ring );
          
#pragma omp parallel for schedule(dynamic) default(shared)
          for ( int i = 0; i < galaxy.n_Ring; ++i ) {

             //Sparse sky except local ring
             if ( Rmin[i] < Rsun && Rmax[i] > Rsun ) {
                galaxy.hpHIR[i].reset( new SM::FullSky<double>(SM::SpectralBinning(std::vector<double>({Rmin[i]}), std::vector<double>({Rmax[i]})), galdef.healpix_order, RING, SM::CoordSys::GAL, 0));
             } else {
                galaxy.hpHIR[i].reset( new SM::SparseSky<double>(SM::SpectralBinning(std::vector<double>({Rmin[i]}), std::vector<double>({Rmax[i]})), galdef.healpix_order, RING, SM::CoordSys::GAL, 0));
             }

             //Loop over the input and insert if value is not 0
             for ( int ihp = 0; ihp < input->Npix(); ++ihp ) 
                if (input->GetReference(ihp, i) != 0) 
                   galaxy.hpHIR[i]->SetValue(ihp,0,input->GetReference(ihp, i)*1e20);

          }

	  if(galdef.verbose==-1003) // selectable debug
	    {
	      filename = configure.fFITSDataDirectory; // + "HIR_healpix.fits";

              for ( int i = 0; i < galaxy.n_Ring; ++i ) {
                 galaxy.hpHIR[i]->ApplyFunction([](double &a){ a*= 1e-20; });

                 ostringstream ost;
                 ost << filename << "HIR_healpix_ring_"<<i+1<<".fits";
                 SM::writeToFits(*galaxy.hpHIR[i], ost.str(), true, true, "Radius", "kpc", "1e20 cm2");

                 galaxy.hpHIR[i]->ApplyFunction([](double &a){ a*= 1e20; });
              }
	    }
	  
	}

      if( type == "COR" )
	{
          galaxy.hpCOR.resize( galaxy.n_Ring );
          
#pragma omp parallel for schedule(dynamic) default(shared)
          for ( int i = 0; i < galaxy.n_Ring; ++i ) {

             //Sparse sky 
             galaxy.hpCOR[i].reset( new SM::SparseSky<double>(SM::SpectralBinning(std::vector<double>({Rmin[i]}), std::vector<double>({Rmax[i]})), galdef.healpix_order, RING, SM::CoordSys::GAL, 0));

             //Loop over the input and insert if value is not 0
             for ( int ihp = 0; ihp < input->Npix(); ++ihp ) 
                if (input->GetReference(ihp, i) != 0) 
                   galaxy.hpCOR[i]->SetValue(ihp,0,input->GetReference(ihp, i));

          }

	  
	  if(galdef.verbose==-1003) // selectable debug
	    {
	      filename = configure.fFITSDataDirectory; 
              for ( int i = 0; i < galaxy.n_Ring; ++i ) {
                 ostringstream ost;
                 ost << filename << "COR_healpix_ring_"<<i+1<<".fits";
                 SM::writeToFits(*galaxy.hpCOR[i], ost.str(), true, true, "Radius", "kpc", "K (km/s)");
              }
	    }
	  
	}
      
      
    }
  else // (l,b) format IMOS20080114
    {   
       int HDU_TYPE;
  if( fits_movabs_hdu(fptr,1,&HDU_TYPE,&status) ) cout<<"8read "<<type<<" status= "<<status<<endl;

  if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) cerr<<"0read "<<type<<" status= "<<status<<endl;
  if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) cerr<<"1read "<<type<<" status= "<<status<<endl;
  if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) cerr<<"2read "<<type<<" status= "<<status<<endl;
  if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) cerr<<"3read "<<type<<" status= "<<status<<endl;
  
  if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) cerr<<"4read "<<type<<" status= "<<status<<endl;
  if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) cerr<<"5read "<<type<<" status= "<<status<<endl;
  
  if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) cerr<<"6read "<<type<<" status= "<<status<<endl;
  if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) cerr<<"7read "<<type<<" status= "<<status<<endl;
  
  ost.str("");
  ost<<" NAXIS      = "<<NAXIS;
  INFO(ost.str());
  ost.str("");
  ost<<" NAXIS1,2,3 = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3;
  INFO(ost.str());
  ost.str("");
  ost<<" CRVAL1,2   = "<<CRVAL1<<" "<<CRVAL2;
  INFO(ost.str());
  ost.str("");
  ost<<" CDELT1,2   = "<<CDELT1<<" "<<CDELT2;
  INFO(ost.str());
  
  long nelements=NAXIS1*NAXIS2*NAXIS3, felement=1;
  float *image=new float[nelements];
  float nulval=1e-15;
  int anynul;

  int i_long,i_lat,i_Ring;
  int i_long_in,i_lat_in;
  int n_long_in=NAXIS1;
  int n_lat_in =NAXIS2;
  galaxy.n_Ring=NAXIS3;

  if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,image,&anynul,&status) )
    cout<<"#read "<<type<<" status= "<<status<<endl;
  
  //Read the ring boundaries
  if( fits_movabs_hdu(fptr,2,&HDU_TYPE,&status) ) cerr<<"8read "<<type<<" status= "<<status<<endl;

  //galaxy.R_bins=new float[2*galaxy.n_Ring]; // two columns: Rmin, Rmax
  galaxy.R_bins.resize(2*galaxy.n_Ring);

  long long int FIRSTROW = 1, FIRSTELEM = 1, NELEMENTS = galaxy.n_Ring;

  float* rBins = new float[2*galaxy.n_Ring];

  for(int i=0;i<2;i++) // read each column
     if( fits_read_col(fptr,TFLOAT,i+1,FIRSTROW,FIRSTELEM,NELEMENTS,&nulval,&rBins[i*galaxy.n_Ring],&anynul,&status) ) 
        cout<<"12read "<<type<<" status= "<<status<<endl;

  for (size_t i = 0; i < galaxy.R_bins.size(); ++i)
     galaxy.R_bins[i] = double(rBins[i]);

  delete[] rBins;

  //The gas rings store the values of the center of the pixel
  
  Distribution HIR_input;
  Distribution COR_input;
  if( type == "HIR" ) HIR_input.init(n_long_in,n_lat_in,galaxy.n_Ring,1); //IMOS20080114
  if( type == "COR" ) COR_input.init(n_long_in,n_lat_in,galaxy.n_Ring,1); //IMOS20080114
  cout<<type<<"_input initialized"<<endl;
  
  for(i_Ring=0; i_Ring<galaxy.n_Ring; i_Ring++) //IMOS20080114
    for(i_lat_in=0; i_lat_in<n_lat_in; i_lat_in++)
      for(i_long_in=0; i_long_in<n_long_in; i_long_in++)
	{
	  long i=i_Ring*n_lat_in*n_long_in +i_lat_in*n_long_in +i_long_in;
	  
	  if( type == "HIR" ) HIR_input.d3[i_long_in][i_lat_in][i_Ring].s[0] = image[i]; 
	  if( type == "COR" ) COR_input.d3[i_long_in][i_lat_in][i_Ring].s[0] = image[i]; 
	}
  
  if( type == "HIR" )
    {
      HIR_input*= 1.0e20; // assigning correct units
      if(galdef.verbose== -201) HIR_input.print();  // selectable debug
    }
  
  if( type == "COR" ) // selectable debug to test vs gas maps IMOS20080114
    {
      if(galdef.verbose== -301) COR_input.print();
    }  

      //redefining the longitude and latitude ranges after reading gas maps in IMOS20080114

  //This does not work for resolution above 2.0 in latitude.  90 has to be divisible by that number.
  //This is anyway a horrible algorithm and should be abandoned.
      //The resolution must be in the form 2^i, with i an integer.  The pixel size will be
      //decreased automatically to accomodate for that.
      galaxy.d_long = pow(2,floor(log(galaxy.d_long)/log(2.)));
      galaxy.d_lat = pow(2,floor(log(galaxy.d_lat)/log(2.)));
   
      //The grid points must align with the gas grid points.  This could be a
      //little bit smarter, but it should be good enough for most purposes. 
      galaxy.long_min = galaxy.d_long/2.+ galaxy.d_long* floor(galaxy.long_min/galaxy.d_long);
      galaxy.long_max =-galaxy.d_long/2.+ galaxy.d_long* ceil (galaxy.long_max/galaxy.d_long);
      galaxy.lat_min  = galaxy.d_lat/2. + galaxy.d_lat * floor(galaxy.lat_min /galaxy.d_lat);
      galaxy.lat_max  =-galaxy.d_lat/2. + galaxy.d_lat * ceil (galaxy.lat_max /galaxy.d_lat);

      //We must stay within bounds, no wrapping of l
      if( galaxy.long_min<=galaxy.d_long/2.     ) galaxy.long_min = galaxy.d_long/2.; 
      if( galaxy.long_max>=360.-galaxy.d_long/2.) galaxy.long_max = 360.-galaxy.d_long/2.;
      if( galaxy.lat_min <=-90.+galaxy.d_lat/2. ) galaxy.lat_min  = -90.+galaxy.d_lat/2.;
      if( galaxy.lat_max >= 90.-galaxy.d_lat/2. ) galaxy.lat_max  = +90.-galaxy.d_lat/2.;

      galaxy.n_long=(int)((galaxy.long_max-galaxy.long_min)/galaxy.d_long + 1.001);
      galaxy.n_lat =(int)((galaxy. lat_max-galaxy. lat_min)/galaxy.d_lat  + 1.001);
      
      //create column density map corresponding to required gamma skymaps by rebinning the input map
      
      double l,b;
      
      cout<<"galaxy."<<type<<".init( "<<galaxy.n_long<<" "<<galaxy.n_lat<<" "<<galaxy.n_Ring<<",1)"<<endl;
      
      if( type == "HIR" ) galaxy.HIR.init(galaxy.n_long,galaxy.n_lat,galaxy.n_Ring,1);
      if( type == "COR" ) galaxy.COR.init(galaxy.n_long,galaxy.n_lat,galaxy.n_Ring,1);
      
      nuse.init(galaxy.n_long,galaxy.n_lat,galaxy.n_Ring,1);

      // choose which grid (input or output gas map) is finer (output grid can be finer!)

      double delt1, delt2, val1, val2, n_long, n_lat;

      delt1 = CDELT1; 
      val1  = CRVAL1; 
      n_long= n_long_in;
 
      if(CDELT1 > galaxy.d_long) 
	{ 
	  delt1 = galaxy.d_long; 
	  val1  = galaxy.long_min; 
	  n_long= galaxy.n_long; 
	}
      
      delt2 = CDELT2; 
      val2  = CRVAL2; 
      n_lat = n_lat_in;
 
      if(CDELT2 > galaxy.d_lat) 
	{ 
	  delt2 = galaxy.d_lat; 
	  val2  = galaxy.lat_min; 
	  n_lat = galaxy.n_lat; 
	}

      for(i_Ring=0; i_Ring<galaxy.n_Ring; i_Ring++)
	{
	  for(int i_b=0; i_b<n_lat; i_b++)
	    {
	      for(int i_l=0; i_l<n_long;i_l++)
		{
		  
		  l=val1 +i_l*delt1;
		  b=val2 +i_b*delt2;

		  i_long=    (int) ((l-galaxy.long_min)/galaxy.d_long+0.5);
		  i_lat =    (int) ((b-galaxy.lat_min )/galaxy.d_lat +0.5);
		  
		  if(i_long<0 || i_long>=galaxy.n_long || i_lat<0 || i_lat>=galaxy.n_lat) continue;
		  
		  i_long_in= (int) ((l-CRVAL1         )/CDELT1       +0.5);
		  i_lat_in = (int) ((b-CRVAL2         )/CDELT2       +0.5);

		  if( type == "HIR" ) galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]+= HIR_input.d3[i_long_in][i_lat_in][i_Ring].s[0];
		  if( type == "COR" ) galaxy.COR.d3[i_long][i_lat][i_Ring].s[0]+= COR_input.d3[i_long_in][i_lat_in][i_Ring].s[0];
  
		  nuse.d3[i_long][i_lat][i_Ring].s[0]+=1;
		  
		  
		  // cout<<"l b i_long_in i_lat_in i_long i_lat i_Ring "<<l<<" "<< b<<" "<< i_long_in<<" "<< i_lat_in<<" "<< i_long<<" "<< i_lat<<" "<< i_Ring<<endl;
		  
		}  //  i_long
	    }   //  i_lat
	}     //  i_Ring
      
      // normalize by number of cells used in rebinning
      for(i_Ring=0; i_Ring<galaxy.n_Ring; i_Ring++)
	for(i_lat=0; i_lat<galaxy.n_lat; i_lat++)
	  for(i_long=0; i_long<galaxy.n_long; i_long++)
	    if(nuse.d3[i_long][i_lat][i_Ring].s[0]>0) 
	      {
		if( type == "HIR" ) galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0] /=nuse.d3[i_long][i_lat][i_Ring].s[0];	      
		if( type == "COR" ) galaxy.COR.d3[i_long][i_lat][i_Ring].s[0] /=nuse.d3[i_long][i_lat][i_Ring].s[0];
	      }
      
      nuse.delete_array();
      
  
  if(galdef.verbose== -203)// selectable debug
    {
      cout<<"read_gas_maps "<<type<<": galaxy."<<type<<":"<<endl;
      if( type == "HIR" ) galaxy.HIR.print();
    }
  
  if(galdef.verbose== -303)// selectable debug
    {
      cout<<"read_gas_maps "<<type<<": galaxy."<<type<<":"<<endl;
      if( type == "COR" ) galaxy.COR.print();
    }

  delete[] image;
    }//HEALPIX OR NOT
  
  fits_close_file(fptr, &status);

  INFO("Exit");

  } //End check for empty COR and HIR names

  return status;
}

