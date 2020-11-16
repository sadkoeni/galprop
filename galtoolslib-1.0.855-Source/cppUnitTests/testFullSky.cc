#include "testFullSky.h"
#include <cmath>
#include "BaseSkyFitsIO.h"

CPPUNIT_TEST_SUITE_REGISTRATION( testFullsky );

using namespace SM;

void testFullsky::setUp() {
}

void testFullsky::tearDown() {
}

//Test functions that should be independent of the storage
void testFullsky::basicFunctions() {
   FullSky<int> mapInt(SpectralBinning(std::vector<double>(4,1.0)), 1, RING, CoordSys::GAL, 0.0);

   SpectralBinning sp(mapInt.GetBinning());

   //Now lets see if they are identical
   CPPUNIT_ASSERT(mapInt.GetBinning().IsIdentical(sp));

   //Create a new map that should be identical to the other
   FullSky<int> mapInt2(std::move(sp), 1, RING, CoordSys::GAL, 0.0);
   CPPUNIT_ASSERT(mapInt.Equivalent(mapInt2));

   //Now lets try to do the same with another type
   FullSky<char> mapChar(mapInt.GetBinning(), 1, RING, CoordSys::GAL, 0.0);
   CPPUNIT_ASSERT(mapInt.Equivalent(mapChar));

   //Modify one item at a time, first the spectrum
   FullSky<int> mapInt3(SpectralBinning(std::vector<double>(3,1.0)), 1, RING, CoordSys::GAL, 0.0);
   CPPUNIT_ASSERT(!mapInt.Equivalent(mapInt3));
   FullSky<int> mapInt4(SpectralBinning(std::vector<double>(4,2.0)), 1, RING, CoordSys::GAL, 0.0);
   CPPUNIT_ASSERT(!mapInt.Equivalent(mapInt4));
   FullSky<int> mapInt5(SpectralBinning(std::vector<double>(4,1.0), std::vector<double>(4,1.0)), 1, RING, CoordSys::GAL, 0.0);
   CPPUNIT_ASSERT(!mapInt.Equivalent(mapInt5));

   //Now the order
   FullSky<int> mapInt6(SpectralBinning(std::vector<double>(4,1.0)), 0, RING, CoordSys::GAL, 0.0);
   CPPUNIT_ASSERT(!mapInt.Equivalent(mapInt6));
   //And the scheme
   FullSky<int> mapInt7(SpectralBinning(std::vector<double>(4,1.0)), 1, NEST, CoordSys::GAL, 0.0);
   CPPUNIT_ASSERT(!mapInt.Equivalent(mapInt7));
   //And coordinate system
   FullSky<int> mapInt8(SpectralBinning(std::vector<double>(4,1.0)), 1, RING, CoordSys::EQ, 0.0);
   CPPUNIT_ASSERT(!mapInt.Equivalent(mapInt8));

   //Make sure the order is correct
   CPPUNIT_ASSERT_EQUAL(mapInt.GetHealpixOrder(), size_t(1));

   //Reset one of the earlier map to make it agree
   mapInt3.Reset(mapInt.GetBinning(), mapInt.Order(), mapInt.Scheme(), mapInt.GetCoordinateSystem(), 0.0);
   CPPUNIT_ASSERT(mapInt.Equivalent(mapInt3));

   //Clone 
   auto mapIntClone = mapInt.clone();
   CPPUNIT_ASSERT(mapIntClone->Equivalent(mapInt));

   //Full sky
   CPPUNIT_ASSERT(mapInt.IsFullSky());

   //Active pixels
   auto pixels = mapInt.GetActivePixels();
   CPPUNIT_ASSERT_EQUAL(pixels.size(), size_t(mapInt.Npix()));

   size_t j = 0;
   for (auto pix : pixels) {
      CPPUNIT_ASSERT(j == pix);
      ++j;
   }
}

void testFullsky::assignment() {

   //Set up a map with 48 pixels and 5 energy bins
   std::vector<double> energies(5);
   for (size_t i = 0; i < energies.size(); ++i)
      energies[i] = 0.5*(i+1);

   FullSky<double> map(SpectralBinning(std::move(energies)), 1, RING, CoordSys::GAL, 0.0);

   //Do simple assignment and retrieval
   CPPUNIT_ASSERT_DOUBLES_EQUAL(map.GetValue(0,0), 0.0, 1e-8);

   double val = 1.3;
   map.SetValue(0,0, val);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(map.GetValue(0,0), val, 1e-8);

   map.GetReference(0,1) = val;
   CPPUNIT_ASSERT_DOUBLES_EQUAL(map.GetValue(0,1), val, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(map.GetReference(0,0), val, 1e-8);

   //Use the coordinate access to set a pixel
   SM::Coordinate co(20,0);
   size_t hpind = map.ang2pix(pointing(M_PI/2.,M_PI/9.));
   const double energy = 2.1;
   const size_t eind = 3;
   
   map.SetValue(co, energy, val);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(map.GetValue(hpind, eind), val, 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(map.GetValue(co, energy), val, 1e-8);

   //Now use a slightly different coordinate to get the same pixel
   SM::Coordinate co2(22,0);
   val = 2.3;
   map.GetReference(co2, energy) =  val;
   CPPUNIT_ASSERT_DOUBLES_EQUAL(map.GetValue(hpind, eind), val, 1e-8);

   //Lets assign values to the entire map
   //This is not the way to do things efficiently but should work for us
   for (size_t hpi=0; hpi < size_t(map.Npix()); ++hpi)
      for (size_t ei=0; ei < 5; ++ei)
         map.SetValue(hpi, ei, 40*hpi+2*ei);

   //Now extract spectra for a specific index
   std::vector<double> sp(map[20]);
   CPPUNIT_ASSERT_EQUAL(size_t(5), sp.size());
   for (size_t i=0; i < sp.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(40*20+2*i, sp[i], 1e-8);

   //And through a coordinate
   sp = map[co];
   for (size_t i=0; i < sp.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(40*hpind+2*i, sp[i], 1e-8);

   //And with the GetSpectrum interface
   sp.resize(0);
   map.GetSpectrum(20, sp);
   CPPUNIT_ASSERT_EQUAL(size_t(5), sp.size());
   for (size_t i=0; i < sp.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(40*20+2*i, sp[i], 1e-8);

   map.GetSpectrum(co, sp);
   for (size_t i=0; i < sp.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(40*hpind+2*i, sp[i], 1e-8);

}

void testFullsky::iterators()
{
   //Set up a map with 48 pixels and 5 energy bins
   std::vector<double> energies(5);
   for (size_t i = 0; i < energies.size(); ++i)
      energies[i] = 0.5*(i+1);

   FullSky<double> map(SpectralBinning(std::move(energies)), 1, RING, CoordSys::GAL, 0.0);

   //Lets assign values to the entire map
   //This is not the way to do things efficiently but should work for us
   for (size_t hpi=0; hpi < size_t(map.Npix()); ++hpi)
      for (size_t ei=0; ei < 5; ++ei)
         map.SetValue(hpi, ei, 40*hpi+2*ei);

   //Make sure we loop correctly over the entire map and get correct indices.
   auto it = map.begin();
   for (size_t j = 0; j < size_t(map.Npix()); ++j)
     for (size_t i = 0; i < 5; ++i,++it) {
       CPPUNIT_ASSERT_EQUAL(it.energyIndex(), i);
       CPPUNIT_ASSERT_EQUAL(it.healpixIndex(), j);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(*it, double(40*j+2*i), 1e-8);
     }

   CPPUNIT_ASSERT(it == map.end());

   //Now lets try to modify the map through the iterators
   const auto endIt = map.end();
   for ( it = map.begin(); it != endIt; ++it )
      *it = 3.4*it.energyIndex() + 54.2*it.healpixIndex();

   //And make sure it happened
   for ( it = map.begin(); it != endIt; ++it )
       CPPUNIT_ASSERT_DOUBLES_EQUAL(*it, double(54.2*it.healpixIndex()+3.4*it.energyIndex()), 1e-8);

   //Now do a single energy plane
   const auto mapEndIt3 = map.mapEnd(3);
   for ( auto mit = map.mapBegin(3); mit != mapEndIt3; ++mit )
      *mit = 5.4*mit.energyIndex() + 5.2*mit.healpixIndex();

   for ( it = map.begin(); it != endIt; ++it )
      if (it.energyIndex() != 3)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(*it, double(54.2*it.healpixIndex()+3.4*it.energyIndex()), 1e-8);
      else
         CPPUNIT_ASSERT_DOUBLES_EQUAL(*it, double(5.2*it.healpixIndex()+5.4*it.energyIndex()), 1e-8);

   //Revert the last change
   for ( auto mit = map.mapBegin(3); mit != mapEndIt3; ++mit )
      *mit = 3.4*mit.energyIndex() + 54.2*mit.healpixIndex();

   //Change a single spectrum
   const auto spEndIt20 = map.spectraEnd(20);
   for ( auto sit = map.spectraBegin(20); sit != spEndIt20; ++sit)
      *sit = 5.4*sit.energyIndex() + 5.2*sit.healpixIndex();

   for ( it = map.begin(); it != endIt; ++it )
      if (it.healpixIndex() != 20)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(*it, double(54.2*it.healpixIndex()+3.4*it.energyIndex()), 1e-8);
      else
         CPPUNIT_ASSERT_DOUBLES_EQUAL(*it, double(5.2*it.healpixIndex()+5.4*it.energyIndex()), 1e-8);

}

void testFullsky::interpolation()
{
   //Set up a map with 48 pixels and 5 energy bins
   std::vector<double> energies(5);
   for (size_t i = 0; i < energies.size(); ++i)
      energies[i] = 0.5*(i+1);

   FullSky<double> map(SpectralBinning(energies), 3, RING, CoordSys::GAL, 0.0);

   //Now lets modify the map through the iterators
   for ( auto it = map.begin(); it != map.end(); ++it )
      *it = 3.4*it.energyIndex() + 54.2*it.healpixIndex();

   //Get the interpolated vector near the GC
   Coordinate co(0.0003,0.002, CoordSys::GAL);
   std::vector<double> sp(map.GetSkyInterpolatedSpectra(co));

   //Need to get the weights using the healpix method for comparison.
   fix_arr<int, 4> pix;
   fix_arr<double, 4> weight;
   map.get_interpol(co.healpixAng(), pix, weight);
   for (size_t i=0; i < sp.size(); ++i) {
      double val = 3.4*i;
      for (size_t j=0; j < 4; ++j)
         val += 54.2*pix[j]*weight[j];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(val, sp[i], 1e-8);
   }

   // Do some tests with the interpolation routine, start by filling the map by ring, linear in theta, power-law in energy
   // NOTE: rings are numbered from 1
   for (int i = 1; i < 4*map.Nside(); ++i) {
      int start, npix;
      double theta;
      bool shifted;
      map.get_ring_info2(i, start, npix, theta, shifted);
      for (int j = start; j < start+npix; ++j) {
         for (size_t l = 0; l < map.GetBinning().GetSize(); ++l) {
            map.SetValue(j, l, theta * pow(energies[l], -2.4) );
         }
      }
   }

   // Create a new map with one higher order, but same binning
   FullSky<double> map2(map.GetBinning(), 4, RING, CoordSys::GAL, 0.0);

   map.Interpolate(map2, false);

   // Loop over the rings in the second map
   for (int i = 2; i < 4*map2.Nside()-1; ++i) {
      int start, npix;
      double theta;
      bool shifted;
      map2.get_ring_info2(i, start, npix, theta, shifted);

      for (int j = start; j < start+npix; ++j) {
         for (size_t l = 0; l < map.GetBinning().GetSize(); ++l) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL( theta*pow(energies[l], -2.4), map2.GetValue(j,l), 1e-8);
         }
      }
   }

   // Do interpolation and extrapolation in energy
   for (size_t i=0; i < energies.size(); ++i) {
      energies[i] = 0.2*(3*i+1);
   }
   FullSky<double> map3(SpectralBinning(energies), 3, RING, CoordSys::GAL, 0.0);

   map.Interpolate(map3, true);

   for (int i = 2; i < 4*map3.Nside()-1; ++i) {
      int start, npix;
      double theta;
      bool shifted;
      map3.get_ring_info2(i, start, npix, theta, shifted);

      for (int j = start; j < start+npix; ++j) {
         for (size_t l = 0; l < map.GetBinning().GetSize(); ++l) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL( theta*pow(energies[l], -2.4), map3.GetValue(j,l), 1e-8);
         }
      }
   }

   // Now make sure we reset the map properly and do not extrapolate
   map.Interpolate(map3, false);
   for (int i = 2; i < 4*map3.Nside()-1; ++i) {
      int start, npix;
      double theta;
      bool shifted;
      map3.get_ring_info2(i, start, npix, theta, shifted);

      for (int j = start; j < start+npix; ++j) {
         for (size_t l = 0; l < map.GetBinning().GetSize(); ++l) {
            if (energies[l] >= 0.5 && energies[l] <= 2.5) {
               CPPUNIT_ASSERT_DOUBLES_EQUAL( theta*pow(energies[l], -2.4), map3.GetValue(j,l), 1e-8);
            } else {
               CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, map3.GetValue(j,l), 1e-8);
            }
         }
      }
   }

   // Try with linear interpolation, once if energies are <= 0, another with values <=0
   // Need the map to be linear in energy
   for (int i = 1; i < 4*map.Nside(); ++i) {
      int start, npix;
      double theta;
      bool shifted;
      map.get_ring_info2(i, start, npix, theta, shifted);
      for (int j = start; j < start+npix; ++j) {
         for (size_t l = 0; l < map.GetBinning().GetSize(); ++l) {
            map.SetValue(j, l, (i-4)*theta * map.GetBinning().GetBins()[l]*0.1 );
         }
      }
   }

   // Make sure we have linear interpolation in the first 4 rings
   // and power-law interpolation for the rest
   map.Interpolate(map3, true);

   for (int i = 1; i < 4*map3.Nside(); ++i) {
      int start, npix;
      double theta;
      bool shifted;
      map3.get_ring_info2(i, start, npix, theta, shifted);

      if ( i < 5 ) {
         for (int j = start; j < start+npix; ++j) {
            for (size_t l = 0; l < map.GetBinning().GetSize(); ++l) {
               CPPUNIT_ASSERT_DOUBLES_EQUAL( (i-4)*theta*energies[l]*0.1, map3.GetValue(j,l), 1e-8);
            }
         }
      } else {
         auto interp = utl::PowerLawInterpolation(map.GetBinning().GetBins(), map[start]);
         for (int j = start; j < start+npix; ++j) {
            for (size_t l = 0; l < energies.size(); ++l) {
               CPPUNIT_ASSERT_DOUBLES_EQUAL( interp(energies[l]), map3.GetValue(j,l), 1e-8 );
            }
         }
      }
   }
   
   for (size_t i=0; i < energies.size(); ++i) {
      energies[i] = 0.2*(3*i);
   }
   map3.Reset(SpectralBinning(energies), 3, RING, CoordSys::GAL, 0.0);
   
   map.Interpolate(map3, true);

   for (int i = 2; i < 4*map3.Nside()-1; ++i) {
      int start, npix;
      double theta;
      bool shifted;
      map3.get_ring_info2(i, start, npix, theta, shifted);

      for (int j = start; j < start+npix; ++j) {
         for (size_t l = 0; l < map3.GetBinning().GetSize(); ++l) {
            CPPUNIT_ASSERT_DOUBLES_EQUAL( (i-4)*theta*energies[l]*0.1, map3.GetValue(j,l), 1e-8);
         }
      }
   }
}

static int pixelCounter;
static inline void count(const double &a) {
#pragma omp critical
   ++pixelCounter;
}

static inline void countExtra(const double &a, size_t i) {
#pragma omp critical
   ++pixelCounter;
}

void testFullsky::applyFunctions()
{
   //Set up a map with 48 pixels and 5 energy bins
   std::vector<double> energies(5);
   for (size_t i = 0; i < energies.size(); ++i)
      energies[i] = 0.5*(i+1);

   FullSky<double> map(SpectralBinning(std::move(energies)), 1, RING, CoordSys::GAL, 0.0);

   //Perform the operations through a BaseSky reference.
   BaseSky<double> &mapRef = map;

   //Assign values of val to all pixels 
   const double expected=2.5;
   mapRef.ApplyFunction( [=](double& pv) { pv = expected; } );

   //And make sure it happened
   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, val, 1e-8);

   //This does nothing to the map but should count the pixels
   pixelCounter = 0;
   mapRef.ApplyFunction( &count );
   CPPUNIT_ASSERT_EQUAL(48*5, pixelCounter);

   //Accumulate the map
   double sum = mapRef.Accumulate();
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expected*48*5, sum, 1e-8);

   //Now lets accumulate the map with a function
   sum = mapRef.ApplyFunctionAccumulate( [](const double &pv) { return exp(pv); } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected)*48*5, sum, 1e-8);

   //Again, this time modifying the map value at the same time
   sum = mapRef.ApplyFunctionAccumulate( [](double &pv) -> double { pv = exp(pv); return log(pv); } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL(sum, expected*48*5, 1e-8);
   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected), val, 1e-8);

   //Now lets try two maps
   auto map2ptr = map.clone();

   //Assign expected to map2 
   map2ptr->ApplyFunction( [=](double& pv) { pv = expected; } );

   //Add map2 to map
   ApplyFunction<double,double>( map, *map2ptr, [](double& pv1, const double& pv2) { pv1 += pv2; } );

   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected)+expected, val, 1e-8);

   //Add it again through the internal apply function
   map.ApplyFunction( [](double& pv1, const double& pv2) { pv1 += pv2; }, *map2ptr);

   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected)+2*expected, val, 1e-8);

   //Delete it using the BaseSky reference
   ApplyFunction<double,double>( mapRef, *map2ptr, [](double& pv1, const double& pv2) { pv1 -= 2*pv2; } );

   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected), val, 1e-8);

   //Now lets see if this works for maps of a different type
   FullSky<int> mapInt(map.GetBinning(), 1, RING, CoordSys::GAL, 0.0);
   mapInt.ApplyFunction( [=](int& pv) { pv = expected; } );

   ApplyFunction<double,int>( map, mapInt, [](double& pv1, const int& pv2) { pv1 += pv2; } );
   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected)+int(expected), val, 1e-8);

   //Again with BaseSky reference
   ApplyFunction<double,int>( mapRef, mapInt, [](double& pv1, const int& pv2) { pv1 -= pv2; } );
   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected), val, 1e-8);

   //A simple accumulate, just sum up the maps
   sum = ApplyFunctionAccumulate<double,int>( *map2ptr, mapInt, [](double &pv1, const int& pv2) { return pv1 + pv2; } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL((expected+int(expected))*48*5, sum, 1e-8);

   BaseSky<double> &map2Ref = *map2ptr;
   const BaseSky<int> &mapIntRef = mapInt;
   sum = ApplyFunctionAccumulate<double,int>( map2Ref, mapIntRef, [](double &pv1, const int& pv2) { return pv1 + pv2; } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL(sum, (expected+int(expected))*48*5, 1e-8);

   const BaseSky<double> &cmap2Ref = *map2ptr;
   sum = ApplyFunctionAccumulate<double,int>( cmap2Ref, mapIntRef, [](const double &pv1, const int& pv2) { return pv1 + pv2; } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL(sum, (expected+int(expected))*48*5, 1e-8);


   //Now test the spectral version of the functions
   map.ApplyFunction( [=](double &a, size_t i){ a = expected*i; } );

   for (auto it = map.begin(); it != map.end(); ++it)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected*it.energyIndex(), *it, 1e-8);

   //This does nothing to the map but should count the pixels
   pixelCounter = 0;
   mapRef.ApplyFunction( &countExtra );
   CPPUNIT_ASSERT_EQUAL(48*5, pixelCounter);

   //Accumulate with a twist
   sum = map.ApplyFunctionAccumulate( [](const double &a, size_t i) { return (i > 0) ? a/i : 0; } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL( expected*48*4, sum, 1e-8);

   //Now modify it as well
   sum = map.ApplyFunctionAccumulate( [=](double &a, size_t i)->double { a = (i > 0) ? a/i : expected; return a; } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL( expected*48*5, sum, 1e-8);

   map.ApplyFunction( [=](const double &a) { CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, a, 1e-8); } );

   //And with another map
   map2ptr->ApplyFunction( [=](double &a){ a = expected; } );
   map.ApplyFunction( [](double &a, const double &b, size_t i){ a = b*(i+1); }, *map2ptr );
   for (auto it = map.begin(); it != map.end(); ++it)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected*(it.energyIndex()+1), *it, 1e-8);

   //Now accumulate
   sum = map.ApplyFunctionAccumulate( [](const double &a, const double &b, size_t i) { return a/(i+1)+b; }, *map2ptr );
   CPPUNIT_ASSERT_DOUBLES_EQUAL( 2*expected*48*5, sum, 1e-8);

   sum = map.ApplyFunctionAccumulate( [](double &a, const double &b, size_t i)->double { a = a/(i+1)+b; return a;}, *map2ptr );
   CPPUNIT_ASSERT_DOUBLES_EQUAL( 2*expected*48*5, sum, 1e-8);
   map.ApplyFunction( [=](const double &a) { CPPUNIT_ASSERT_DOUBLES_EQUAL(2*expected, a, 1e-8); } );

}

void testFullsky::coordinateConversion()
{
   //Set up a map with 48 pixels and 5 energy bins
   std::vector<double> energies(5);
   for (size_t i = 0; i < energies.size(); ++i)
      energies[i] = 0.5*(i+1);

   FullSky<double> map(SpectralBinning(std::move(energies)), 4, RING, CoordSys::GAL, 0.0);

   //Now lets modify the map through the iterators
   const auto endIt = map.end();
   for ( auto it = map.begin(); it != endIt; ++it )
      *it = 1e-3 + 1e-5*it.energyIndex() * 5.2*it.healpixIndex();

   //DEBUG: write out the maps and investigate manually
   writeToFits(map, "OriginalConversionTmpFile.fits.gz");

   //Convert the map to ECL and back
   auto eclMapPtr = map.CoordinateConversion(CoordSys::ECL);
   writeToFits(*eclMapPtr, "ECLConversionTmpFile.fits.gz");
   auto eclGalMapPtr = eclMapPtr->CoordinateConversion(CoordSys::GAL);
   writeToFits(*eclGalMapPtr, "ECLGALConversionTmpFile.fits.gz");

   //Assert the maps are same
   CPPUNIT_ASSERT(eclGalMapPtr->Equivalent(map));

   //There is a rather large fluctuation involved in double conversion at low resolution
   map.ApplyFunction([](const double &a, const double &b){ CPPUNIT_ASSERT_DOUBLES_EQUAL(a, b, 0.8); }, *eclGalMapPtr);

   //Clean up the memory
   delete eclGalMapPtr.release();

   //Now convert twice, to eq and then to ecl and compare to ecl map
   auto eqMapPtr = map.CoordinateConversion(CoordSys::EQ);
   auto eqEclMapPtr = eqMapPtr->CoordinateConversion(CoordSys::ECL);

   //Assert they are the same
   CPPUNIT_ASSERT(eqEclMapPtr->Equivalent(*eclMapPtr));

   eqEclMapPtr->ApplyFunction([](const double &a, const double &b)->void { CPPUNIT_ASSERT_DOUBLES_EQUAL(a, b, 0.8); }, *eclMapPtr);

}
