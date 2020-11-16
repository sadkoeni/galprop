#include "BaseSkyFitsIO.h"

#include <fitsio.h>
#include <iostream>
#include <string>
#include <vector>

#ifdef HAVE_WCS
extern "C" {
#include <wcslib/wcs.h>
#include <wcslib/wcshdr.h>
#include <wcslib/wcserr.h>
#include <wcslib/wcsfix.h>
}
#endif


namespace SM {

template <typename T>
void writeToFits(const BaseSky<T> &map,
      const std::string &fileName,
      bool overwrite,
      bool compress,
      const std::string &energyType,
      const std::string &energyUnit,
      const std::string &mapUnit,
      const std::map<std::string, std::string> &additionalFitsKeywords
      ) 
{

   //Type mapping from c++ to fits
   static std::map<const char*, std::string> formatMap;
   formatMap[typeid(char).name()] = "I";
   formatMap[typeid(short).name()] = "I";
   formatMap[typeid(int).name()] = "J";
   formatMap[typeid(long).name()] = "K";
   formatMap[typeid(float).name()] = "E";
   formatMap[typeid(double).name()] = "D";
   formatMap[typeid(unsigned char).name()] = "B";
   formatMap[typeid(unsigned short).name()] = "U";
   formatMap[typeid(unsigned int).name()] = "V";
   formatMap[typeid(unsigned long).name()] = "K";

   auto typeIt = formatMap.find(typeid(T).name());

   if (typeIt == formatMap.end()) {
      ERROR("Incorrect type supplied.  We only support basic types");
      throw(std::invalid_argument("Only basic types supported"));
   }

   //Create the FITS object
   //Prepend ! so the file gets overwritten
   std::string fname = (fileName[0] != '!' && overwrite) ? "!" + fileName : fileName;
   //Append .gz to allow compression
   if ( compress && fname.substr(fname.size()-3).compare(".gz") )
      fname += ".gz";

   //Create a CCFITS object with the filename
   CCfits::FITS fits(fname, CCfits::Write );
         
   //Create the data format string
   std::ostringstream strForm;
   strForm << 1 << typeIt->second;

   //Need to create the columns before creating the table, the number of columns depends on the type
   std::vector<size_t> activePixels;
   std::string indexScheme, object;
   std::vector<std::string> colNames, colForm, colUnit;
   if (map.IsFullSky()) {
      indexScheme = "IMPLICIT";
      object = "FULLSKY";

      colNames.resize(map.GetBinning().GetSize());
      colForm.resize(map.GetBinning().GetSize(), strForm.str());
      colUnit.resize(map.GetBinning().GetSize(), mapUnit);
      
      for (size_t i(0); i < colNames.size(); ++i) {
         std::ostringstream oss;
         oss<<"Bin"<<i;
         colNames[i] = oss.str();
      }

   } else {

      indexScheme = "EXPLICIT";
      object = "PARTIAL";

      //Retrieve the active pixels for the map
      auto pixSet = map.GetActivePixels();
      activePixels.reserve(pixSet.size());
      activePixels.insert(activePixels.begin(), pixSet.begin(), pixSet.end());

      colNames.resize(map.GetBinning().GetSize()+1);
      colForm.resize(map.GetBinning().GetSize()+1, strForm.str());
      colUnit.resize(map.GetBinning().GetSize()+1, mapUnit);

      colNames[0] = "HPXINDEX";
      colForm[0] = "1K";
      colUnit[0] = "";
      for (size_t i(1); i < colNames.size(); ++i) {
         std::ostringstream oss;
         oss<<"Bin"<<i-1;
         colNames[i] = oss.str();
      }
   }

   //Get the number of rows to write
   const size_t totalRows = activePixels.size() > 0 ? activePixels.size() : map.Npix();

   CCfits::Table *imgTable = fits.addTable("SKYMAP2", totalRows, colNames, colForm, colUnit);

   //Write keywords
   imgTable->addKey("PIXTYPE", "HEALPIX", "Healpix pixel scheme");
   std::string keyvalue = "RING";
   if (NEST == map.Scheme()) keyvalue = "NESTED";
   imgTable->addKey("ORDERING", keyvalue, "Ring or nested ordering of pixels");
   imgTable->addKey("NSIDE", map.Nside(), "Number of sides in a base pixel");

   imgTable->addKey("INDXSCHM", indexScheme, "Indexing scheme");
   imgTable->addKey("OBJECT", object, "Type of object, partial or fullsky");
   if (activePixels.size() > 0) {
      imgTable->addKey("FIRSTPIX", activePixels.front(), "Number of first pixel");
      imgTable->addKey("LASTPIX", activePixels.back(), "Number of last pixel");
   } else {
      imgTable->addKey("FIRSTPIX", 0, "Number of first pixel");
      imgTable->addKey("LASTPIX", map.Npix()-1, "Number of last pixel");
   }
   imgTable->addKey("MAPTYPE", map.name(), "Type of map");

   std::ostringstream osst;
   osst << "Number of "<<energyType<<" bins in spectra";
   imgTable->addKey("NBRBINS", map.GetBinning().GetSize(), osst.str());
   std::string coord;
   switch (map.GetCoordinateSystem()) {
      case CoordSys::EQ:
         coord = "EQ";
         break;
      case CoordSys::GAL:
         coord = "GAL";
         break;
      case CoordSys::ECL:
         coord = "ECL";
         break;
   }
   imgTable->addKey("COORDSYS", coord, "Coordinate system");

   imgTable->addKey("EMPTYVAL", map.GetEmptyValue(), "Value for non-assigned pixels");


   //Write additional keywords
   std::map<std::string, std::string>::const_iterator it;
   for (it = additionalFitsKeywords.begin(); it != additionalFitsKeywords.end(); ++it)
      imgTable->addKey(it->first, it->second, "");


   //Get the most efficient number of rows to write out.
   long nrows(totalRows);
   
   //Only important for non-compressed files
   if ( fname.substr(fname.size()-3).compare(".gz") ) {
      long nrows(10);
      int status(0);
      fits_get_rowsize(fits.fitsPointer(),&nrows,&status);
      //default to writing the whole map at a time if this fails
      if (nrows < 1)
         nrows = totalRows;
   }

   //Get the values of the rows into a temporary storage
   //Consider adding an interface to get a unique pointer to continuous
   //memory for a list of pixels.
   std::vector< std::vector<T> > values(map.GetBinning().GetSize());

   //Temporary array for the pixels, if needed
   std::vector<size_t> apWrite;

   for (size_t irow(0); irow < totalRows; irow += nrows) {

      const size_t stoprow = std::min(irow+nrows,totalRows);

      if (activePixels.size() > 0) {
         apWrite.resize(stoprow-irow);
         for (size_t j = irow; j < stoprow; ++j) 
            apWrite[j-irow] = activePixels[j];

         //Write the pixels to the file
         imgTable->column("HPXINDEX").write(apWrite, irow+1);

         //Get the pixels into the values vectors
         for (size_t i = 0; i < values.size(); ++i) {
            values[i].resize(stoprow-irow);
            for (size_t j = irow; j < stoprow; ++j) {
               values[i][j-irow] = map.GetValue(activePixels[j], i);
            }
         }

      } else {

         //Get the pixels into the values vectors
         for (size_t i = 0; i < values.size(); ++i) {
            values[i].resize(stoprow-irow);
            for (size_t j = irow; j < stoprow; ++j) {
               values[i][j-irow] = map.GetValue(j, i);
            }
         }

      }

      //Write the values
      for (size_t i = 0; i < values.size(); ++i) {
         std::ostringstream oss;
         oss<<"Bin"<<i;
         imgTable->column(oss.str()).write(values[i], irow+1);
      }

   }

   //Create another table to store the energy
   if (map.GetBinning().IsIntegrated()){
      colNames.resize(3); colUnit.resize(3); colForm.resize(3);
      colNames[0] = "CHANNEL"; colNames[1] = energyType + "_MIN"; colNames[2] = energyType + "_MAX";
      colUnit[0] = ""; colUnit[1] = energyUnit; colUnit[2] = energyUnit;
      colForm[0] = "I"; colForm[1] = "D"; colForm[2] = "D";
      CCfits::Table *eTable = fits.addTable("EBOUNDS", values.size(), colNames, colForm, colUnit);

      //Create the channel array
      std::valarray<int> channel(values.size());
      for (int i = 0; i < int(channel.size()); ++i){
         channel[i] = i+1;
      }

      eTable->column("CHANNEL").write(channel,1);
      eTable->column(energyType + "_MIN").write(map.GetBinning().GetBinsMin(),1);
      eTable->column(energyType + "_MAX").write(map.GetBinning().GetBinsMax(),1);

   } else {

      colNames.resize(1); colUnit.resize(1); colForm.resize(1);
      colNames[0] = energyType;
      colUnit[0] = energyUnit;
      colForm[0] = "D";
      CCfits::Table *eTable = fits.addTable("ENERGIES", values.size(), colNames, colForm, colUnit);

      eTable->column(energyType).write(map.GetBinning().GetBins(),1);
   }
         
}




template <typename T>
void readFromFits(std::unique_ptr< BaseSky<T> > &map_ptr,
      const std::string &fileName,
      std::map<std::string, std::string> &keywords
      ) 
{

   //Open the first hdu that has the pixtype keyword equal to healpix.
   std::vector< std::string > lookupKeyword(1,"");
   std::vector< std::string > values(1,"");
   lookupKeyword[0] = "PIXTYPE";
   values[0] = "HEALPIX";
   CCfits::FITS fits(fileName, CCfits::Read, lookupKeyword, values);

   //A reference to the table containing the skymap data
   CCfits::ExtHDU &skymapTable = fits.currentExtension(); 

   //Get the properties of the binning, these are required
   int nSide(1);
   std::string ordering;
   skymapTable.readKey("NSIDE", nSide);
   skymapTable.readKey("ORDERING", ordering);
   
   //Calculate the order
   const size_t hporder = size_t(log(double(nSide))/log(2.0)+0.1);

   //Try to get the indexing scheme, we only support EXCPLICIT with index in the first row
   std::string indexScheme = "IMPLICIT";
   std::string object;
   try {
      skymapTable.readKey("INDXSCHM", indexScheme);
   } catch (CCfits::HDU::NoSuchKeyword) {
      try {
         skymapTable.readKey("OBJECT", object);
         if (object == "PARTIAL") {
            indexScheme = "EXPLICIT";
         }
      } catch (CCfits::HDU::NoSuchKeyword) {}
   }

   //Extract the coordinate system, defaults to Galactic
   CoordSys co = CoordSys::GAL;
   try {
      std::string coordsys;
      skymapTable.readKey("COORDSYS", coordsys);
      if (coordsys == "EQ")
         co = CoordSys::EQ;
      else if (coordsys == "ECL")
         co = CoordSys::ECL;
   } catch (CCfits::HDU::NoSuchKeyword) {}

   //Get the maptype, needed if map_ptr is null.  Default to FullSky
   std::string maptype("FullSky");
   try {
      skymapTable.readKey("MAPTYPE", maptype);
   } catch (CCfits::HDU::NoSuchKeyword) {}

   //The value for empty pixels
   //Read it in as double and do casting later
   double empty(0);
   try {
      skymapTable.readKey("EMPTYVAL", empty);
   } catch (CCfits::HDU::NoSuchKeyword) {}

   //Now read all keywords into the string structure
   skymapTable.readAllKeys();
   for (auto kit : skymapTable.keyWord()) {
      //There is no automatic conversion from numerical values to strings so we need the switch statement.
      //Find the value type
      CCfits::ValueType vtype = kit.second->keytype();
      std::string value;
      std::ostringstream os;
      switch(vtype) {
         case CCfits::Tstring:
            kit.second->value(value);
            break;

         case CCfits::Tlogical:
            bool btmp;
            kit.second->value(btmp);
            os << btmp;
            value = os.str();
            break;

         case CCfits::Tbyte:
            char ctmp;
            kit.second->value(ctmp);
            os << ctmp;
            value = os.str();
            break;

         case CCfits::Tshort:
            short stmp;
            kit.second->value(stmp);
            os << stmp;
            value = os.str();
            break;

         case CCfits::Tushort:
            unsigned short ustmp;
            kit.second->value(ustmp);
            os << ustmp;
            value = os.str();
            break;

         case CCfits::Tint:
            int itmp;
            kit.second->value(itmp);
            os << itmp;
            value = os.str();
            break;

         case CCfits::Tuint:
            unsigned int uitmp;
            kit.second->value(uitmp);
            os << uitmp;
            value = os.str();
            break;

         case CCfits::Tlong:
            long ltmp;
            kit.second->value(ltmp);
            os << ltmp;
            value = os.str();
            break;

         case CCfits::Tulong:
            unsigned long ultmp;
            kit.second->value(ultmp);
            os << ultmp;
            value = os.str();
            break;

         case CCfits::Tlonglong:
            long long lltmp;
            kit.second->value(lltmp);
            os << lltmp;
            value = os.str();
            break;

         case CCfits::Tfloat:
            float ftmp;
            kit.second->value(ftmp);
            os << ftmp;
            value = os.str();
            break;
         case CCfits::Tdouble:
            double dtmp;
            kit.second->value(dtmp);
            os << dtmp;
            value = os.str();
            break;

         default:
            value = "";
      }

      keywords[kit.first] = value;
   }


   //Try to get the number of energy planes, assume 1 as the default
   int nPlanes(1);
   try {
      skymapTable.readKey("NBRBINS", nPlanes);
   } catch (CCfits::HDU::NoSuchKeyword) {
      //Assume one plane per column
      nPlanes = (indexScheme == "EXPLICIT") ? skymapTable.numCols() - 1 : skymapTable.numCols();
   }

   //Look for ENERGIES or EBOUNDS extension to get the energy planes
   std::unique_ptr<SpectralBinning> sp_pointer;

   try {

      //Try the ENERGIES table
      fits.read(std::vector<std::string> (1,"ENERGIES"));
      CCfits::ExtHDU &energyTable = fits.extension("ENERGIES");

      //Extract the data
      int nSpectra;
      energyTable.readKey("NAXIS2", nSpectra);
      std::vector<double> energies(nSpectra);
      energyTable.column(1).read(energies, 1, nSpectra);

      //Create the binning
      sp_pointer.reset(new SpectralBinning(energies));

   } catch (CCfits::FITS::NoSuchHDU) {

      try { 

         //Then EBOUNDS table
         fits.read(std::vector<std::string> (1,"EBOUNDS"));
         CCfits::ExtHDU &energyTable = fits.extension("EBOUNDS");

         //Extract the data
         int nSpectra;
         energyTable.readKey("NAXIS2", nSpectra);
         std::vector<double> eMin(nSpectra), eMax(nSpectra);
         energyTable.column(2).read(eMin, 1, nSpectra);
         energyTable.column(3).read(eMax, 1, nSpectra);

         //Create the binning
         sp_pointer.reset(new SpectralBinning(eMin, eMax));

      } catch (CCfits::FITS::NoSuchHDU) { }

   }

   //Create a default binning with indexes if neither table is found
   if (sp_pointer.get() == nullptr) {
      WARNING("Could not read energy extensions (ENERGIES or EBOUNDS), defaulting to indexing for the planes.");
      std::vector<double> energies(nPlanes);
      for (size_t i(0); i < energies.size(); ++i)
         energies[i] = i;
      sp_pointer.reset(new SpectralBinning(energies));
   }

   //Make sure the number of planes is identical to the size of the binning
   if (size_t(nPlanes) != sp_pointer->GetSize()) {
      throw(std::runtime_error("Number of planes does not match energies extension"));
   }

   const Healpix_Ordering_Scheme sch = (ordering == "RING") ? RING : NEST;

   //Create the map, or reset it depending on input ptr
   if (map_ptr.get() == nullptr) {
      map_ptr = BaseSky<T>::create(maptype, *sp_pointer, hporder, sch, co, T(empty));
   } else {
      map_ptr->Reset(*sp_pointer, hporder, sch, co, T(empty));
   }

   //Get the best number of rows to read from the file
   skymapTable.makeThisCurrent();
   int totalRows;
   skymapTable.readKey("NAXIS2", totalRows);

   long nrows(10);
   int status(0);
   fits_get_rowsize(fits.fitsPointer(),&nrows,&status);
   //default to reading the whole map at a time if this fails
   if (nrows < 1)
      nrows = totalRows;

   //Now we read the data, what happens now depends on the indexing scheme
   if (indexScheme == "EXPLICIT") {

      //Make sure we have the correct number of columns
      if ( size_t(skymapTable.numCols()) != 1+sp_pointer->GetSize() ) {
         throw(std::runtime_error("Incorrect number of columns in the file"));
      }

      //Storage for reading in the maps
      std::vector<size_t> hpxIndex;
      std::vector<T> values;

      //We simply assume the first column is Healpix index and the planes are stored after that
      for (size_t irow(0); irow < size_t(totalRows); irow += nrows) {

         const size_t stop = std::min(irow+nrows, size_t(totalRows));
         skymapTable.column(1).read(hpxIndex, irow+1, stop);

         //Throw an error if the index is out of bounds
         for (size_t j(0); j < hpxIndex.size(); ++j)
            if ( hpxIndex[j] >= size_t(map_ptr->Npix()) ) 
               throw(std::runtime_error("Healpix index out of bounds"));

         for (size_t isp(0); isp < sp_pointer->GetSize(); ++isp) {
            skymapTable.column(isp+2).read(values, irow+1, stop);

            for (size_t j(0); j < values.size(); ++j) 
               if (values[j] != T(empty))
                  map_ptr->SetValue(hpxIndex[j], isp, values[j]);
         }
      }

   } else {

      //IMPLICIT scheme, no index column

      //Make sure the number of rows is correct
      if ( totalRows != map_ptr->Npix() ) {
         throw(std::runtime_error("Incorrect number of rows in the file"));
      }

      if (skymapTable.numCols() == 1 && sp_pointer->GetSize() > 1) {

         //Old format with a single vector column.
         std::vector< std::valarray<T> > values;

         for (size_t irow(0); irow < size_t(totalRows); irow += nrows) {
            const size_t stop = std::min(irow+nrows, size_t(totalRows));

            skymapTable.column(1).readArrays(values, irow+1, stop);

            //Assert that the spectral size is correct
            if (values[0].size() != sp_pointer->GetSize())
               throw(std::runtime_error("Incorrect number of elements in vector row"));

            for (size_t j(0); j < values.size(); ++j) 
               for (size_t isp(0); isp < values[j].size(); ++isp)
                  map_ptr->SetValue(irow+j, isp, values[j][isp]);
         }
         
      } else {

         //Make sure we have the correct number of columns
         if ( size_t(skymapTable.numCols()) != sp_pointer->GetSize() ) {
            throw(std::runtime_error("Incorrect number of columns in the file"));
         }

         std::vector<T> values;

         for (size_t irow(0); irow < size_t(totalRows); irow += nrows) {

            const size_t stop = std::min(irow+nrows, size_t(totalRows));

            for (size_t isp(0); isp < sp_pointer->GetSize(); ++isp) {
               skymapTable.column(isp+1).read(values, irow+1, stop);

               for (size_t j(0); j < values.size(); ++j) 
                  map_ptr->SetValue(irow+j, isp, values[j]);
            }
         }
      }

   }
}


/** \brief  This routine reads a WCS map into a BaseSky.
 *
 * Any data in map is erased and the spectral information is read 
 * from the fits file.  Spatial information from the input map
 * is kept as is and coordinate transform is performed as needed.
 *
 * If solidAngleCorrection is false then map is assumed to be in units 
 * sr^-1 (intensity), otherwise we take into account the different solid
 * angle of the input and output.
 */
template <typename T>
void readFromFitsWcs( std::unique_ptr< BaseSky<T> > &map_ptr, const std::string &fileName, bool solidAngleCorrection)
{
#ifdef HAVE_WCS
   //Read from the primary image
   CCfits::FITS fits(fileName);
   CCfits::PHDU &mapCube = fits.pHDU();

   wcserr_enable(1);

   std::ostringstream oss;

   DEBUGLOG("Using WCS to read the fits file");
   //Read the file using the wcslib
   char *header;
   int nkeys, status(0);

   mapCube.makeThisCurrent();

   fits_convert_hdr2str(mapCube.fitsPointer(), 1, NULL, 0, &header, &nkeys, &status);
   if (status != 0) {
      ERROR("Could not read header of file");
      throw(std::runtime_error("Broken file"));
   }

   //Read the WCS keywords
   wcsprm *wcsData;
   int ctrl=0, nreject, nwcs;
   status = wcspih(header, nkeys, 0, ctrl, &nreject, &nwcs, &wcsData);
   if (status != 0) {
      ERROR("Could not read WCS information");
      throw(std::runtime_error("Broken file"));
   }

   fits_free_memory(header, &status);
   if (status != 0) {
      ERROR("Could not free memory");
      throw(std::runtime_error("Major error"));
   }

   //Assert that we have at least 1 projection and the first has 2 axis, nothing complicated here
   if (nwcs < 1 || wcsData[0].naxis < 2) {
      ERROR("Number of axis less than 2, cannot continue");
      throw(std::runtime_error("Number of axis less than 2, cannot continue"));
   }

   //Make sure the files do not have too many axis
   if (wcsData[0].naxis > 3) {
      ERROR("More than 3 axis in the image, cannot continue");
      throw(std::runtime_error("Can at most handle a single spectrum"));
   }

   wcsprm *wcs = &wcsData[0];

   //Apply fixes to the header, e.g. units may be the wrong type
   int fixstat[NWCSFIX];
   status = wcsfix(0, NULL, wcs, fixstat);
   if (status != 0) {
      ERROR("WCSFIX returned an error");
      throw(std::runtime_error("Broken file"));
   }
   status = wcsset( wcs );
   if (status != 0) {
      ERROR("WCSSET returned an error");
      throw(std::runtime_error("Broken file"));
   }

   //Just for fun
   if (mapCube.axes() != wcs->naxis) {
      ERROR("Dimensions for wcs and FITS do not match");
      throw(std::runtime_error("This just shouldn't happen"));
   }

   //Need two celestial axis in the map
   if (wcs->lat < 0 || wcs->lng < 0) {
      ERROR("Need two celestial axis");
      throw(std::runtime_error("Cannot convert this map"));
   }

   //Find the longitude and latitude range of the input map
   //Loop over all pixels and find the range
   double lonmin = std::numeric_limits<double>::max();
   double lonmax = -std::numeric_limits<double>::max();
   double latmin = std::numeric_limits<double>::max();
   double latmax = -std::numeric_limits<double>::max();

   size_t mapPix = mapCube.axis(wcs->lat);
   double worldTmp[mapPix][3];
   double phiTmp[mapPix], thetaTmp[mapPix];
   double imgcrdTmp[mapPix][3], pixcrdTmp[mapPix][3];
   int statTmp[mapPix];

   for (int i0 = 0; i0 < mapCube.axis(wcs->lng); ++i0) {

      size_t count(0);
      for (int i1 = 0; i1 < mapCube.axis(wcs->lat); ++i1) {
         pixcrdTmp[count][wcs->lng] = i0+1;
         pixcrdTmp[count][wcs->lat] = i1+1;
         ++count;
      }

      wcsp2s(wcs, count, 3, pixcrdTmp[0], imgcrdTmp[0], phiTmp, thetaTmp, worldTmp[0], statTmp);

      for (size_t i = 0; i < count; ++i) {
         lonmin = std::min(lonmin, worldTmp[i][wcs->lng] - fabs(wcs->cdelt[wcs->lng]));
         lonmax = std::max(lonmax, worldTmp[i][wcs->lng] + fabs(wcs->cdelt[wcs->lng]));
         latmin = std::min(latmin, worldTmp[i][wcs->lat] - fabs(wcs->cdelt[wcs->lat]));
         latmax = std::max(latmax, worldTmp[i][wcs->lat] + fabs(wcs->cdelt[wcs->lat]));
      }

   }

   oss.str("");
   oss<<"Axis for longitude: "<<wcs->lng+1;
   DEBUGLOG(oss.str());
   oss.str("");
   oss<<"Axis for latitude: "<<wcs->lat+1;
   DEBUGLOG(oss.str());

   int spec = wcs->spec;
   if (spec < 0) {
      int i = 0;
      while (spec < 0 && i < mapCube.axes()) {
         if (i != wcs->lat && i != wcs->lng) 
            spec = i;
         ++i;
      }
   }

   //We use the Coordinate class to handle any projection
   std::string Coord = wcs->lngtyp;
   //Try to use ctype if Coord is empty, handles old gasmaps for galprop
   if (Coord == "    ")
      Coord = wcs->ctype[wcs->lng];

   CoordSys cSys = CoordSys::GAL;
   oss.str("");
   oss<<"Coordinates of input map is";
   if (Coord == "RA") {
      cSys = CoordSys::EQ;
      oss<<" Equatorial";
   } else if (Coord == "ELON") {
      cSys = CoordSys::ECL;
      oss<<" Ecliptic";
   } else if (Coord == "GLON") {
      cSys = CoordSys::GAL;
      oss<<" Galactic";
   } else {
      ERROR("Unknown coordinate system \""+Coord+"\", assuming Galactic");
      oss<<" Galactic";
   }
   DEBUGLOG(oss.str());
   oss.str("");

   oss<<"Maximum dimensions of input map (longitude) : (latitude) : ("<<lonmin<<", "<<lonmax<<") : ("<<latmin<<", "<<latmax<<")";
   DEBUGLOG(oss.str());
   oss.str("");

   oss<<"Projection name: "<<wcs->cel.prj.name;
   DEBUGLOG(oss.str());

   for ( int i(0); i < wcs->naxis; ++i) {
      oss.str("");
      oss<<"CRVAL"<<i+1<<": "<<wcs->crval[i];
      DEBUGLOG(oss.str());
      oss.str("");
      oss<<"CRPIX"<<i+1<<": "<<wcs->crpix[i];
      DEBUGLOG(oss.str());
      oss.str("");
      oss<<"CDELT"<<i+1<<": "<<wcs->cdelt[i];
      DEBUGLOG(oss.str());
      oss.str("");
      oss<<"CUNIT"<<i+1<<": "<<wcs->cunit[i];
      DEBUGLOG(oss.str());
      oss.str("");
      oss<<"CTYPE"<<i+1<<": "<<wcs->ctype[i];
      DEBUGLOG(oss.str());
      for (int j(0); j < wcs->naxis; ++j) {
         oss.str("");
         oss<<"PC"<<i+1<<"_"<<j+1<<": "<<wcs->pc[wcs->naxis*i+j];
         DEBUGLOG(oss.str());
      }
      oss.str("");
      oss<<"CROTA"<<i+1<<": "<<wcs->crota[i];
      DEBUGLOG(oss.str());
      for (int j(0); j < wcs->naxis; ++j) {
         oss.str("");
         oss<<"CD"<<i+1<<"_"<<j+1<<": "<<wcs->cd[wcs->naxis*i+j];
         DEBUGLOG(oss.str());
      }
   }

   oss.str("");
   oss<<"cel.flag: "<<wcs->cel.flag;
   DEBUGLOG(oss.str());
   oss.str("");
   oss<<"cel.offset: "<<wcs->cel.offset;
   DEBUGLOG(oss.str());
   oss.str("");
   oss<<"cel.phi0: "<<wcs->cel.phi0;
   DEBUGLOG(oss.str());
   oss.str("");
   oss<<"cel.theta0: "<<wcs->cel.theta0;
   DEBUGLOG(oss.str());
   oss.str("");
   oss<<"cel.ref[0]: "<<wcs->cel.ref[0];
   DEBUGLOG(oss.str());
   oss.str("");
   oss<<"cel.ref[1]: "<<wcs->cel.ref[1];
   DEBUGLOG(oss.str());
   oss.str("");
   oss<<"cel.ref[2]: "<<wcs->cel.ref[2];
   DEBUGLOG(oss.str());
   oss.str("");
   oss<<"cel.ref[3]: "<<wcs->cel.ref[3];
   DEBUGLOG(oss.str());

   //If this is a proper mapcube, read in the energies value
   std::valarray<double> energies, emin, emax;
   try {
      //Try energies extension first
      CCfits::ExtHDU & energyTable = fits.extension("ENERGIES");
      int nSpectra;
      energyTable.readKey("NAXIS2", nSpectra);
      energyTable.column(1).read(energies, 1, nSpectra);
      DEBUGLOG("Energies from ENERGIES extension");
   } catch (CCfits::FITS::NoSuchHDU) {
      try{
         //Then ebounds
         CCfits::ExtHDU & energyTable = fits.extension("EBOUNDS");
         int nSpectra;
         energyTable.readKey("NAXIS2", nSpectra);
         energyTable.column(2).read(emin, 1, nSpectra);
         energyTable.column(3).read(emax, 1, nSpectra);
         energies.resize(emin.size());
         DEBUGLOG("Energies from EBOUNDS extension");
      } catch (CCfits::FITS::NoSuchHDU) {}
   }

   //If no energies extension is found, create the values from CRVAL and CDELT values
   //Assume the GALPROP convention for storing the spectral information axis
   if (energies.size() == 0) {
      if (wcs->naxis < 3) {
         DEBUGLOG("No energy axis in map");
         energies.resize(1);
         energies[0] = 1;
      } else {
         energies.resize(mapCube.axis(spec));
         //Assume CDELT represents logarithmic values
         for (size_t i = 0; i < energies.size(); ++i) {
            energies[i] = pow(10,wcs->crval[spec] + i*wcs->cdelt[spec]);
         }
         oss.str("");
         oss<<"Energy axes created from axis "<<spec<<" from input map assuming CDELT is log10";
         DEBUGLOG(oss.str());
      }
   }

   //Now we can set up the map
   std::unique_ptr<SpectralBinning> spbinning;
   if (emin.size() == 0) {
      spbinning.reset(new SpectralBinning(std::vector<double>(&energies[0], &energies[0]+energies.size())));
   } else {
      spbinning.reset(new SpectralBinning(std::vector<double>(&emin[0], &emin[0]+emin.size()), std::vector<double>(&emax[0], &emax[0]+emax.size())));
   }

   //Create the map, or reset it depending on input ptr
   if (map_ptr.get() == nullptr) {
      //Calculate the order
      const double res = std::min(fabs(wcs->cdelt[wcs->lng]), fabs(wcs->cdelt[wcs->lat]));
      size_t hporder = size_t( log( sqrt(3./utl::kPi)*60/res ) / log(2.0) ) + 1;
      if (hporder > 13)
         hporder = 13;

      //Use fullSky only if both angles cover more than half of the sky
      std::string mapType = "SparseSky";
      if ( latmax-latmin > 90 && lonmax-lonmin > 180 )
         mapType = "FullSky";

      map_ptr = BaseSky<T>::create(mapType, *spbinning, hporder, RING, cSys, 0);
   } else {
      map_ptr->Reset(*spbinning, map_ptr->Order(), map_ptr->Scheme(), map_ptr->GetCoordinateSystem(), map_ptr->GetEmptyValue());
   }

   //Read the image data and fill the skymap
   std::valarray<T> image;
   long npix(1);
   for (int i = 0; i < wcs->naxis; ++i){
      npix *= mapCube.axis(i);
   }
   mapCube.read(image,1,npix);

   //Correct the map for solid angle if necessary
   if (solidAngleCorrection) {
      //One more than number of pixels
      size_t nvert = (mapCube.axis(wcs->lng)+1)*(mapCube.axis(wcs->lat)+1);

      //Arrays for conversion
      double world[nvert][3];
      double phi[nvert];  //Also used for solid angle
      double theta[nvert];
      double imgcrd[nvert][3];
      double pixcrd[nvert][3];
      int stat[nvert];

      size_t ind(0);
      for (int ib = 0; ib < mapCube.axis(wcs->lat)+1; ++ib) {
         for (int il = 0; il < mapCube.axis(wcs->lng)+1; ++il) {
            //Reset the values, in case spec is not defined
            pixcrd[ind][0] = pixcrd[ind][1]= pixcrd[ind][2] = 1;
            pixcrd[ind][wcs->lng] = il+0.5;
            pixcrd[ind][wcs->lat] = ib+0.5;
            ++ind;
         }
      }
      wcsp2s(wcs, ind, 3, pixcrd[0], imgcrd[0], phi, theta, world[0], stat);

      //For conversion between il, ib, ispec to map index
      size_t ilmul = 1;
      if (wcs->lng == 1)
         ilmul *= mapCube.axis(0);
      else if (wcs->lng == 2)
         ilmul *= mapCube.axis(0)*mapCube.axis(1);
      size_t ibmul = 1;
      if (wcs->lat == 1)
         ibmul *= mapCube.axis(0);
      else if (wcs->lat == 2)
         ibmul *= mapCube.axis(0)*mapCube.axis(1);
      size_t ispecmul = 1;
      if (spec == 1)
         ispecmul *= mapCube.axis(0);
      else if (spec == 2)
         ispecmul *= mapCube.axis(0)*mapCube.axis(1);

      //Loop through the pixels and apply the solid angle correction
#pragma omp parallel for default(shared) schedule(static)
      for (int ib = 0; ib < mapCube.axis(wcs->lat); ++ib) {
         for (int il = 0; il < mapCube.axis(wcs->lng); ++il) {
            double solidAngle(0);

            //The indices to the pixel corners
            const size_t i1 = ib*(mapCube.axis(wcs->lng)+1)+il;
            const size_t i2 = i1+1;
            const size_t i4 = i2+mapCube.axis(wcs->lng);
            const size_t i3 = i4+1;

            //Need the distances between the corners, cache sin and cos of the lat
            const double sinlat1 = std::sin(world[i1][wcs->lat]*utl::kConvertDegreesToRadians);
            const double coslat1 = std::cos(world[i1][wcs->lat]*utl::kConvertDegreesToRadians);
            const double sinlat2 = std::sin(world[i2][wcs->lat]*utl::kConvertDegreesToRadians);
            const double coslat2 = std::cos(world[i2][wcs->lat]*utl::kConvertDegreesToRadians);
            const double sinlat3 = std::sin(world[i3][wcs->lat]*utl::kConvertDegreesToRadians);
            const double coslat3 = std::cos(world[i3][wcs->lat]*utl::kConvertDegreesToRadians);
            const double sinlat4 = std::sin(world[i4][wcs->lat]*utl::kConvertDegreesToRadians);
            const double coslat4 = std::cos(world[i4][wcs->lat]*utl::kConvertDegreesToRadians);

            //Distances between adjacent corners
            const double d12 = chacos(sinlat1*sinlat2+coslat1*coslat2*std::cos((world[i1][wcs->lng]-world[i2][wcs->lng])*utl::kConvertDegreesToRadians));
            const double d14 = chacos(sinlat1*sinlat4+coslat1*coslat4*std::cos((world[i1][wcs->lng]-world[i4][wcs->lng])*utl::kConvertDegreesToRadians));
            const double d32 = chacos(sinlat3*sinlat2+coslat3*coslat2*std::cos((world[i3][wcs->lng]-world[i2][wcs->lng])*utl::kConvertDegreesToRadians));
            const double d34 = chacos(sinlat3*sinlat4+coslat3*coslat4*std::cos((world[i3][wcs->lng]-world[i4][wcs->lng])*utl::kConvertDegreesToRadians));

            //In case of triangles
            if (d12 <= 0 || d14 <= 0) {
               //Diagonal
               const double d24 = chacos(sinlat2*sinlat4+coslat2*coslat4*std::cos((world[i2][wcs->lng]-world[i4][wcs->lng])*utl::kConvertDegreesToRadians));

               //Cache for sines and cosines
               const double sind32 = std::sin(d32);
               const double sind24 = std::sin(d24);
               const double sind34 = std::sin(d34);
               const double cosd32 = std::cos(d32);
               const double cosd24 = std::cos(d24);
               const double cosd34 = std::cos(d34);

               //The angles
               const double a2 = chacos((cosd34-cosd24*cosd32)/(sind24*sind32));
               const double a3 = chacos((cosd24-cosd34*cosd32)/(sind34*sind32));
               const double a4 = chacos((cosd32-cosd24*cosd34)/(sind24*sind34));

               solidAngle = (a2+a3+a4) - utl::kPi;

            } else if ( d32 <= 0 || d34 <= 0) {
               //Diagonal
               const double d24 = chacos(sinlat2*sinlat4+coslat2*coslat4*std::cos((world[i2][wcs->lng]-world[i4][wcs->lng])*utl::kConvertDegreesToRadians));

               //Cache for sines and cosines
               const double sind12 = std::sin(d12);
               const double sind24 = std::sin(d24);
               const double sind14 = std::sin(d14);
               const double cosd12 = std::cos(d12);
               const double cosd24 = std::cos(d24);
               const double cosd14 = std::cos(d14);

               //The angles
               const double a1 = chacos((cosd24-cosd12*cosd14)/(sind12*sind14));
               const double a2 = chacos((cosd14-cosd12*cosd24)/(sind12*sind24));
               const double a4 = chacos((cosd12-cosd14*cosd24)/(sind14*sind24));

               solidAngle = (a1+a2+a4) - utl::kPi;

            } else {
               //Complete polygon
               const double d13 = chacos(sinlat1*sinlat3+coslat1*coslat3*std::cos((world[i1][wcs->lng]-world[i3][wcs->lng])*utl::kConvertDegreesToRadians));
               const double d24 = chacos(sinlat2*sinlat4+coslat2*coslat4*std::cos((world[i2][wcs->lng]-world[i4][wcs->lng])*utl::kConvertDegreesToRadians));

               //Cache for sines and cosines
               const double sind12 = std::sin(d12);
               const double sind14 = std::sin(d14);
               const double sind32 = std::sin(d32);
               const double sind34 = std::sin(d34);
               const double cosd12 = std::cos(d12);
               const double cosd13 = std::cos(d13);
               const double cosd14 = std::cos(d14);
               const double cosd24 = std::cos(d24);
               const double cosd32 = std::cos(d32);
               const double cosd34 = std::cos(d34);

               //The angles
               const double a1 = chacos((cosd13-cosd34*cosd14)/(sind34*sind14));
               const double a2 = chacos((cosd24-cosd32*cosd34)/(sind32*sind34));
               const double a3 = chacos((cosd13-cosd12*cosd32)/(sind12*sind32));
               const double a4 = chacos((cosd24-cosd14*cosd12)/(sind14*sind12));

               solidAngle = (a1+a2+a3+a4) - utl::kTwoPi;
            }

            if (solidAngle > 0) {
               if (spec < 0) 
                  image[il*ilmul+ib*ibmul] /= solidAngle;
               else {
                  const size_t imind = il*ilmul+ib*ibmul;
                  for (int ispec = 0; ispec < mapCube.axis(spec); ++ispec)
                     image[imind+ispec*ispecmul] /= solidAngle;
               }
            }

         }
      }

   }

   //Divide each HEALPix pixel into 256 smaller ones and fill those with the
   //value from the map pixel directly beneath its center.  The HEALPix
   //pixel value is the average of those 256 smaller ones.
   const int dOrder = std::min(5,13-map_ptr->Order()); //Order can not exceed 13
   const int nPix = (1<<dOrder)*(1<<dOrder);

   //Create a nested healpix base for us to work with
   const Healpix_Base finer(map_ptr->Order()+dOrder, NEST);

   //The solid angle
   const double SA = solidAngleCorrection ? map_ptr->solidAngle() : 1.0;

#pragma omp parallel for default(shared) schedule(static) private(status)
   for (int p = 0; p<map_ptr->Npix(); ++p){

      //Use pixel numbering from the nested scheme
      const int pp = map_ptr->Scheme() == NEST ? p : map_ptr->ring2nest(p);

      //Set up arrays for all the pixels, use only the first index in the spectral coordinate
      double world[nPix][3];
      double phi[nPix];
      double theta[nPix];
      double imgcrd[nPix][3];
      double pixcrd[nPix][3];
      int stat[nPix];

      //Loop over all of the subpixels
      size_t count(0);
      for (int sp = pp*nPix; sp<(pp+1)*nPix; ++sp){
         const pointing pnt = finer.pix2ang(sp);

         Coordinate co(pnt, map_ptr->GetCoordinateSystem());


         double tl, tb;
         co.getCoordinates(tl, tb, cSys);

         //Convert to degrees
         tb *= 180./utl::kPi;
         tl *= 180./utl::kPi;
         if (tb > latmin && tb < latmax) {

            if ( tl > lonmax )
               tl -= 360;

            if ( tl < lonmin )
               tl += 360;

            if ( tl > lonmin && tl < lonmax ) {
               //Reset the values, in case spec is not defined
               world[count][0] = world[count][1]= world[count][2] = 1;
               world[count][wcs->lng] = tl;
               world[count][wcs->lat] = tb;
               ++count;
            } 
         }

      }

      //std::cout<<"Pixel number "<<p<<" has "<<count<<" valid sub-pixels"<<std::endl;

      //No valid pixels
      if (count == 0)
         continue;

      //Do the conversion
      status = wcss2p(wcs, count, 3, world[0], phi, theta, imgcrd[0], pixcrd[0], stat);
      if (status != 0) {
         ERROR("wcss2p failed, cannot convert map to healpix");
         throw(std::runtime_error("WCS failure"));
      }

      //To store the total value for all the pixels
      std::vector<T> totSpectra(spbinning->GetSize(), T(0));

      //Now we must loop over the pixels and take the average
      //What happens next depends on spec
      size_t numValid(0);
      if ( spec < 0 ) {
         for (size_t i(0); i < count; ++i) {
            //No spectral information, only a single value
            int i0 = int(round(pixcrd[i][0]-1));
            int i1 = int(round(pixcrd[i][1]-1));

            //They must both be within bounds to do something useful
            if (i0 >= 0 && i0 < mapCube.axis(0) && i1 >= 0 && i1 < mapCube.axis(1)) {
               ++numValid;
               const int ind1 = i0 + i1*mapCube.axis(0);
               //There is only one value in the spectra
               totSpectra[0] += image[ind1];
            } 
         }
      } else if (spec == 0) {
         for (size_t i(0); i < count; ++i) {
            //Spectral information in first axis
            int i1 = int(round(pixcrd[i][1]-1));
            int i2 = int(round(pixcrd[i][2]-1));
            //std::cout<<i1<<", "<<i2<<" are the pixel coordinates"<<std::endl;

            //They must both be within bounds to do something useful
            if (i1 >= 0 && i1 < mapCube.axis(1) && i2 >= 0 && i2 < mapCube.axis(2)) {
               ++numValid;
               const int ind1 = i1*mapCube.axis(0) + i2*mapCube.axis(0)*mapCube.axis(1);

#if (_OPENMP >= 201307)
//#pragma omp simd
#endif
               for (size_t is=0; is<spbinning->GetSize(); ++is){
                  totSpectra[is] += image[ind1+is];
               }
            }
         }
      } else if (spec == 1) {
         for (size_t i(0); i < count; ++i) {
            //Spectral information in second axis
            int i0 = int(round(pixcrd[i][0]-1));
            int i2 = int(round(pixcrd[i][2]-1));
            //std::cout<<i0<<", "<<i2<<" are the pixel coordinates"<<std::endl;

            //They must both be within bounds to do something useful
            if (i0 >= 0 && i0 < mapCube.axis(0) && i2 >= 0 && i2 < mapCube.axis(2)) {
               ++numValid;
               const int ind1 = i0 + i2*mapCube.axis(0)*mapCube.axis(1);

#if (_OPENMP >= 201307)
//#pragma omp simd
#endif
               for (size_t is=0; is<spbinning->GetSize(); ++is){
                  totSpectra[is] += image[ind1+is*mapCube.axis(0)];
               }
            }
         }
      } else if (spec == 2) {
         for (size_t i(0); i < count; ++i) {
            //Spectral information in third axis
            int i0 = int(round(pixcrd[i][0]-1));
            int i1 = int(round(pixcrd[i][1]-1));
            //std::cout<<i0<<", "<<i1<<" are the pixel coordinates"<<std::endl;

            //They must both be within bounds to do something useful
            if (i0 >= 0 && i0 < mapCube.axis(0) && i1 >= 0 && i1 < mapCube.axis(1)) {
               ++numValid;
               const int ind1 = i0 + i1*mapCube.axis(0);

#if (_OPENMP >= 201307)
//#pragma omp simd
#endif
               for (size_t is=0; is<spbinning->GetSize(); ++is){
                  totSpectra[is] += image[ind1+is*mapCube.axis(0)*mapCube.axis(1)];
               }
            }
         }
      }

      if (numValid > 0) {
         //std::cout<<"Doing pixel number "<<p<<" with values ";
         for (size_t l = 0; l < spbinning->GetSize(); ++l){
            map_ptr->SetValue(p, l, T(totSpectra[l]/double(numValid)*SA));
            //std::cout<<T(totSpectra[l]/double(numValid)*SA)<<", ";
         }
         //std::cout<<std::endl;
      }
   }


   wcsvfree(&nwcs, &wcsData);



#else
   ERROR("This functionality requires the WCS library");
   throw(std::runtime_error("Recompile with WCS"));
#endif
}



//Instantiate them
template void writeToFits(const BaseSky<char> &map, const std::string &, bool, bool, const std::string &, const std::string &, const std::string &, const std::map<std::string, std::string> &);
template void writeToFits(const BaseSky<short> &map, const std::string &, bool, bool, const std::string &, const std::string &, const std::string &, const std::map<std::string, std::string> &);
template void writeToFits(const BaseSky<int> &map, const std::string &, bool, bool, const std::string &, const std::string &, const std::string &, const std::map<std::string, std::string> &);
template void writeToFits(const BaseSky<long> &map, const std::string &, bool, bool, const std::string &, const std::string &, const std::string &, const std::map<std::string, std::string> &);
template void writeToFits(const BaseSky<float> &map, const std::string &, bool, bool, const std::string &, const std::string &, const std::string &, const std::map<std::string, std::string> &);
template void writeToFits(const BaseSky<double> &map, const std::string &, bool, bool, const std::string &, const std::string &, const std::string &, const std::map<std::string, std::string> &);
template void writeToFits(const BaseSky<unsigned char> &map, const std::string &, bool, bool, const std::string &, const std::string &, const std::string &, const std::map<std::string, std::string> &);
template void writeToFits(const BaseSky<unsigned short> &map, const std::string &, bool, bool, const std::string &, const std::string &, const std::string &, const std::map<std::string, std::string> &);
template void writeToFits(const BaseSky<unsigned int> &map, const std::string &, bool, bool, const std::string &, const std::string &, const std::string &, const std::map<std::string, std::string> &);
template void writeToFits(const BaseSky<unsigned long> &map, const std::string &, bool, bool, const std::string &, const std::string &, const std::string &, const std::map<std::string, std::string> &);

template void readFromFits(std::unique_ptr< BaseSky<char> > &, const std::string &, std::map<std::string, std::string> &);
template void readFromFits(std::unique_ptr< BaseSky<short> > &, const std::string &, std::map<std::string, std::string> &); 
template void readFromFits(std::unique_ptr< BaseSky<int> > &, const std::string &, std::map<std::string, std::string> &);
template void readFromFits(std::unique_ptr< BaseSky<long> > &, const std::string &, std::map<std::string, std::string> &); 
template void readFromFits(std::unique_ptr< BaseSky<float> > &, const std::string &, std::map<std::string, std::string> &); 
template void readFromFits(std::unique_ptr< BaseSky<double> > &, const std::string &, std::map<std::string, std::string> &); 
template void readFromFits(std::unique_ptr< BaseSky<unsigned char> > &, const std::string &, std::map<std::string, std::string> &);
template void readFromFits(std::unique_ptr< BaseSky<unsigned short> > &, const std::string &, std::map<std::string, std::string> &); 
template void readFromFits(std::unique_ptr< BaseSky<unsigned int> > &, const std::string &, std::map<std::string, std::string> &);
template void readFromFits(std::unique_ptr< BaseSky<unsigned long> > &, const std::string &, std::map<std::string, std::string> &); 

template void readFromFitsWcs(std::unique_ptr< BaseSky<char> > &, const std::string &, bool);
template void readFromFitsWcs(std::unique_ptr< BaseSky<short> > &, const std::string &, bool);
template void readFromFitsWcs(std::unique_ptr< BaseSky<int> > &, const std::string &, bool);
template void readFromFitsWcs(std::unique_ptr< BaseSky<long> > &, const std::string &, bool);
template void readFromFitsWcs(std::unique_ptr< BaseSky<float> > &, const std::string &, bool);
template void readFromFitsWcs(std::unique_ptr< BaseSky<double> > &, const std::string &, bool);
template void readFromFitsWcs(std::unique_ptr< BaseSky<unsigned char> > &, const std::string &, bool);
template void readFromFitsWcs(std::unique_ptr< BaseSky<unsigned short> > &, const std::string &, bool);
template void readFromFitsWcs(std::unique_ptr< BaseSky<unsigned int> > &, const std::string &, bool);
template void readFromFitsWcs(std::unique_ptr< BaseSky<unsigned long> > &, const std::string &, bool);
}//Namespace
