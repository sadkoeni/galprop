#include "testSparseSky.h"
#include <cmath>
#include "BaseSkyFitsIO.h"
#include "FullSky.h"

CPPUNIT_TEST_SUITE_REGISTRATION( testSparsesky );

using namespace SM;

void testSparsesky::setUp() {
}

void testSparsesky::tearDown() {
}

//Test functions that should be independent of the storage
void testSparsesky::basicFunctions() {
   SparseSky<int> mapInt(SpectralBinning(std::vector<double>(4,1.0)), 1, RING, CoordSys::GAL, 0.0);

   SpectralBinning sp(mapInt.GetBinning());

   //Create a new map that should be identical to the other
   FullSky<int> mapInt2(std::move(sp), 1, RING, CoordSys::GAL, 0.0);
   CPPUNIT_ASSERT(mapInt.Equivalent(mapInt2));

}

void testSparsesky::assignment() {

   //Set up a map with 48 pixels and 5 energy bins
   std::vector<double> energies(5);
   for (size_t i = 0; i < energies.size(); ++i)
      energies[i] = 0.5*(i+1);

   SparseSky<double> map(SpectralBinning(std::move(energies)), 1, RING, CoordSys::GAL, 0.0);

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

   //Now use a slightly different coordinate to get the same pixel
   SM::Coordinate co2(22,0);
   val = 2.3;
   map.GetReference(co2, energy) =  val;
   CPPUNIT_ASSERT_DOUBLES_EQUAL(map.GetValue(hpind, eind), val, 1e-8);

   //Get the spectrum for a single pixel
   std::vector<double> sp(map[hpind]);
   CPPUNIT_ASSERT_EQUAL(size_t(5), sp.size());
   for (size_t i = 0; i < sp.size(); ++i) {
      if ( i == eind )
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val, sp[i], 1e-8);
      else
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, sp[i], 1e-8);
   }

   //Lets assign values to the entire map
   for (size_t hpi=0; hpi < size_t(map.Npix()); ++hpi)
      for (size_t ei=0; ei < 5; ++ei)
         map.SetValue(hpi, ei, 40*hpi+2*ei);

   //Now extract spectra for a specific index
   sp = map[20];
   for (size_t i=0; i < sp.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sp[i], 40*20+2*i, 1e-8);

   //And through a coordinate
   sp = map[co];
   for (size_t i=0; i < sp.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sp[i], 40*hpind+2*i, 1e-8);

   //Now extract spectra for a specific index
   sp.resize(0);
   map.GetSpectrum(20, sp);
   CPPUNIT_ASSERT_EQUAL(size_t(5), sp.size());
   for (size_t i=0; i < sp.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sp[i], 40*20+2*i, 1e-8);

   //And through a coordinate
   map.GetSpectrum(co, sp);
   for (size_t i=0; i < sp.size(); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sp[i], 40*hpind+2*i, 1e-8);


}

void testSparsesky::iterators()
{
   //Set up a map with 48 pixels and 5 energy bins
   std::vector<double> energies(5);
   for (size_t i = 0; i < energies.size(); ++i)
      energies[i] = 0.5*(i+1);

   SparseSky<double> map(SpectralBinning(std::move(energies)), 1, RING, CoordSys::GAL, 0.0);

   //Lets make sure that all start and end iterators are the same
   CPPUNIT_ASSERT(map.begin() == map.end());
   for (size_t j = 0; j < map.GetBinning().GetSize(); ++j)
      CPPUNIT_ASSERT(map.mapBegin(j) == map.mapEnd(j));
   for (int i = 0; i < map.Npix(); ++i)
      CPPUNIT_ASSERT(map.spectraBegin(i) == map.spectraEnd(i));

   //Lets assign values to the entire map
   for (size_t hpi=0; hpi < size_t(map.Npix()); ++hpi)
      for (size_t ei=0; ei < 5; ++ei)
         map.SetValue(hpi, ei, 40*hpi+2*ei);

   //Make sure we loop correctly over the entire map and get correct indices.
   auto it = map.begin();
   for (size_t i = 0; i < 5; ++i) 
      for (size_t j = 0; j < size_t(map.Npix()); ++j) {
       CPPUNIT_ASSERT_EQUAL(i, it.energyIndex());
       CPPUNIT_ASSERT_EQUAL(j, it.healpixIndex());
       CPPUNIT_ASSERT_DOUBLES_EQUAL(double(40*j+2*i), *it, 1e-8);
       ++it;
     }

   CPPUNIT_ASSERT(it == map.end());

   //Now lets try to modify the map through the iterators
   const auto itEnd = map.end();
   for ( it = map.begin(); it != itEnd; ++it )
      *it = 3.4*it.energyIndex() + 54.2*it.healpixIndex();

   //And make sure it happened
   for ( it = map.begin(); it != itEnd; ++it )
       CPPUNIT_ASSERT_DOUBLES_EQUAL(double(54.2*it.healpixIndex()+3.4*it.energyIndex()), *it, 1e-8);

   //Clear the map 
   map.Reset(map.GetBinning(), map.Order(), map.Scheme(), map.GetCoordinateSystem(), 0.0);

   //Make sure beginning and end is the same here
   CPPUNIT_ASSERT(map.begin() == map.end());

   //Assign only a few values near the North pole
   for (size_t hpi = 0; hpi < 4; ++hpi)
      for (size_t ei = 0; ei < 5; ++ei)
         map.SetValue(hpi, ei, 40*hpi+2*ei);

   //Loop over those values and assert they are fine
   size_t count(0);
   for ( it = map.begin(); it != map.end(); ++it ) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(double(40*it.healpixIndex()+2*it.energyIndex()), *it, 1e-8);
      ++count;
   }
   CPPUNIT_ASSERT_EQUAL(size_t(20), count);

   //Get the spectrum at the south pole and assert it is 0
   auto spSp = map[Coordinate(0,-utl::kPiOnTwo, CoordSys::GAL)];
   CPPUNIT_ASSERT_EQUAL(spSp.size(), size_t(5));
   for (auto val : spSp)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0, val, 1e-8);

   //Make sure the number of active pixels is correct
   auto activePixels = map.GetActivePixels();
   CPPUNIT_ASSERT_EQUAL(activePixels.size(), size_t(4));
   auto pit = activePixels.begin();
   for (size_t i = 0; i < activePixels.size(); ++i, ++pit)
      CPPUNIT_ASSERT_EQUAL(i, *pit);

   //Get the healpix map at middle energies
   auto hpMap = map.GetHealpixMap(size_t(2));
   CPPUNIT_ASSERT_EQUAL(hpMap.Npix(), map.Npix());
   for (size_t i = 0; i < size_t(hpMap.Npix()); ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(map.GetValue(i,2), hpMap[i], 1e-8);

}

void testSparsesky::interpolation()
{
   //Set up a map with 48 pixels and 5 energy bins
   std::vector<double> energies(5);
   for (size_t i = 0; i < energies.size(); ++i)
      energies[i] = 0.5*(i+1);

   SparseSky<double> map(SpectralBinning(std::move(energies)), 1, RING, CoordSys::GAL, 1.0);

   //Now lets modify the map 
   for (size_t hpi=10; hpi < size_t(map.Npix()-10); ++hpi)
      for (size_t ei=0; ei < 5; ++ei)
         map.SetValue(hpi, ei, 3.4*ei + 54.2*hpi);

   //Get the interpolated vector at the GC
   Coordinate co(0,0);
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

   //Make sure the Interpolation routine works the same here as for FullSky
   FullSky<double> fsmap(map.GetBinning(), 1, RING, CoordSys::GAL, 1.0);

   ApplyFunction<double, double>(fsmap, map, [](double &a, const double &b){ a = b; });

   SparseSky<double> intMap(map.GetBinning(), 3, RING, CoordSys::ECL, 1.0);
   FullSky<double> fsintMap(map.GetBinning(), 3, RING, CoordSys::ECL, 1.0);

   map.Interpolate(intMap, false);
   fsmap.Interpolate(fsintMap, false);

   for ( auto it = fsintMap.begin(); it != fsintMap.end(); ++it ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(*it, intMap.GetValue(it.healpixIndex(), it.energyIndex()), 1e-8);

}

void testSparsesky::applyFunctions()
{
   //Set up a map with 48 pixels and 5 energy bins
   std::vector<double> energies(5);
   for (size_t i = 0; i < energies.size(); ++i)
      energies[i] = 0.5*(i+1);

   SparseSky<double> map(SpectralBinning(std::move(energies)), 1, RING, CoordSys::GAL, 0.0);

   //Perform the operations through a BaseSky reference.
   BaseSky<double> &mapRef = map;

   //Set a few values of the map
   for (size_t hpi=10; hpi < 20; ++hpi) 
      for (size_t ei=1; ei < 4; ++ei)
         mapRef.SetValue(hpi, ei, 2.4*(ei+0.3) + 4.3*(hpi-0.2));

   //Assign values of val to all active pixels 
   const double expected=2.5;
   mapRef.ApplyFunction( [=](double& pv) { pv = expected; } );

   //And make sure it happened
   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, val, 1e-8);

   //Make sure the first pixel and the empty value was not affected
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, mapRef.GetValue(0,0), 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, mapRef.GetEmptyValue(), 1e-8);

   //Accumulate the map
   double sum = mapRef.Accumulate();
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expected*10*3, sum, 1e-8);

   //Now lets accumulate the map with a function
   sum = mapRef.ApplyFunctionAccumulate( [](const double &pv) { return exp(pv); } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected)*10*3 + 38*5+10*2, sum, 1e-8);

   //Again, this time modifying the map value at the same time
   sum = mapRef.ApplyFunctionAccumulate( [=](double &pv) -> double { pv = exp(expected); return log(pv); } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL(expected*48*5, sum, 1e-8);
   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected), val, 1e-8);

   //Make sure the empty value is also affected
   CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected), mapRef.GetValue(0,0), 1e-8);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected), mapRef.GetEmptyValue(), 1e-8);

   //Now lets try two maps, clone the map
   auto map2ptr = map.clone();

   //Assign expected to map2  
   map2ptr->ApplyFunction( [=](double& pv) { pv = expected; } );
   map2ptr->SetEmptyValue(expected);

   //Add map2 to map
   ApplyFunction<double,double>( map, *map2ptr, [](double& pv1, const double& pv2) { pv1 += pv2; } );
   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected)+expected, val, 1e-8);

   //Make sure the empty value is also affected
   CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected)+expected, mapRef.GetEmptyValue(), 1e-8);
   //The active pixels should only include the currently active pixels
   CPPUNIT_ASSERT_EQUAL(size_t(10), map.GetActivePixels().size());
   //The iterators should only iterate over the 3 set planes
   size_t count(0);
   for (auto val : *map2ptr )
      ++count;
   CPPUNIT_ASSERT_EQUAL(size_t(30), count);


   //Make sure map2 was not affected 
   CPPUNIT_ASSERT_EQUAL(size_t(10), map2ptr->GetActivePixels().size());

   //Do the same with reference
   ApplyFunction<double,double>( mapRef, *map2ptr, [](double& pv1, const double& pv2) { pv1 -= pv2; } );
   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected), val, 1e-8);

   //Make sure the empty value is also affected
   CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected), mapRef.GetEmptyValue(), 1e-8);
   //The active pixels should only include the currently active pixels
   CPPUNIT_ASSERT_EQUAL(size_t(10), map.GetActivePixels().size());
   //The iterators should now iterate over all 5 planes
   count = 0;
   for (auto val : map )
      ++count;
   CPPUNIT_ASSERT_EQUAL(size_t(50), count);

   //Make sure map2 was not affected 
   CPPUNIT_ASSERT_EQUAL(size_t(10), map2ptr->GetActivePixels().size());


   //Now lets see if this works for maps of a different type
   SparseSky<float> mapFloat(map.GetBinning(), 1, RING, CoordSys::GAL, expected);

   CPPUNIT_ASSERT_EQUAL(float(expected), mapFloat.GetEmptyValue());

   ApplyFunction<double,float>( map, mapFloat, [](double& pv1, const float& pv2) { pv1 += pv2; } );
   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected)+float(expected), val, 1e-8);

   //A simple accumulate, just sum up the maps
   sum = ApplyFunctionAccumulate<double,float>( *map2ptr, mapFloat, [](const double &pv1, const float& pv2) { return pv1 + pv2; } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL((expected+float(expected))*48*5, sum, 1e-8);

   //Now lets see if this works for maps of a different type
   SparseSky<long> mapInt(map.GetBinning(), 1, RING, CoordSys::GAL, expected);

   CPPUNIT_ASSERT_EQUAL(long(expected), mapInt.GetEmptyValue());

   ApplyFunction<double,long>( map, mapInt, [](double& pv1, const long& pv2) { pv1 += pv2; } );
   for ( auto val : map ) 
      CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(expected)+expected+long(expected), val, 1e-8);

   //A simple accumulate, just sum up the maps
   sum = ApplyFunctionAccumulate<double,long>( *map2ptr, mapInt, [](const double &pv1, const long& pv2) { return pv1 + pv2; } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL((expected+long(expected))*48*5, sum, 1e-8);

   //Test the BaseSky implementations by using a reference
   BaseSky<double> &map2Ref = *map2ptr;
   const BaseSky<long> &mapIntRef = mapInt;
   sum = ApplyFunctionAccumulate<double,long>( map2Ref, mapIntRef, [](double &pv1, const long& pv2) { return pv1 + pv2; } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL((expected+long(expected))*48*5, sum, 1e-8);

   const BaseSky<double> &cmap2Ref = *map2ptr;
   sum = ApplyFunctionAccumulate<double,long>( cmap2Ref, mapIntRef, [](const double &pv1, const long& pv2) { return pv1 + pv2; } );
   CPPUNIT_ASSERT_DOUBLES_EQUAL((expected+long(expected))*48*5, sum, 1e-8);

   //Test the internal implementation, reset the maps
   map.Reset(map.GetBinning(), map.Order(), map.Scheme(), map.GetCoordinateSystem(), 0.0);
   map2ptr->Reset(map.GetBinning(), map.Order(), map.Scheme(), map.GetCoordinateSystem(), expected);

   //Assign a few values to maps
   map.SetValue(5,1,expected);
   map.SetValue(8,2,expected);
   map.SetValue(40,4,expected);

   map2ptr->SetValue(1,0, expected);
   map2ptr->SetValue(8,2, expected);

   //Accumulate the maps using a constant function that adds their values
   sum = map.ApplyFunctionAccumulate( [](const double &a, const double &b) { return a + b; }, *map2ptr );
   CPPUNIT_ASSERT_DOUBLES_EQUAL( expected*(48*5+3), sum, 1e-8);

   //Do the same with an index
   sum = map.ApplyFunctionAccumulate( [](const double &a, const double &b, size_t i) { return a*i + b; }, *map2ptr );
   CPPUNIT_ASSERT_DOUBLES_EQUAL( expected*(48*5+7), sum, 1e-8);

   //Now add map2 to the map using the apply function, it only affects active pixels
   map.ApplyFunction( [](double &a, const double &b) { a += b; }, *map2ptr );

   count = 0;
   for (auto val : map) {
      if (count == 0)
         CPPUNIT_ASSERT_DOUBLES_EQUAL( 1*expected, val, 1e-8);
      else
         CPPUNIT_ASSERT_DOUBLES_EQUAL( 2*expected, val, 1e-8);
      ++count;
   }
   CPPUNIT_ASSERT_EQUAL( size_t(4), count);
   
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, map.GetEmptyValue(), 1e-8);

   //Now accumulate with the non-const version, but without assigning anything
   sum = map.ApplyFunctionAccumulate( [](double &a, const double &b) { return a + b; }, *map2ptr );
   CPPUNIT_ASSERT_DOUBLES_EQUAL( expected*(48*5+7), sum, 1e-8);

   //Do the same with an index
   sum = map.ApplyFunctionAccumulate( [](double &a, const double &b, size_t i) { return a*i + b; }, *map2ptr );
   CPPUNIT_ASSERT_DOUBLES_EQUAL( expected*(48*5+14), sum, 1e-8);


}

void testSparsesky::coordinateConversion()
{
   //Set up a map with 48 pixels and 5 energy bins
   std::vector<double> energies(5);
   for (size_t i = 0; i < energies.size(); ++i)
      energies[i] = 0.5*(i+1);

   SparseSky<double> map(SpectralBinning(std::move(energies)), 4, RING, CoordSys::GAL, 1e-10);

   //Assign values to a few pixels only
   for ( size_t hpi = 10; hpi < 20; ++hpi)
      for (size_t ei = 1; ei < 4; ++ei)
         map.SetValue(hpi, ei, 1e-3 + 1e-5*ei + 5.2*hpi);

   //Create the same map with FullSky
   FullSky<double> fsMap(map.GetBinning(), map.Order(), map.Scheme(), map.GetCoordinateSystem(), 1e-10);

   auto mapPtr = map.clone();
   fsMap.ApplyFunction([](double &a, const double &b){ a = b; }, *mapPtr);

   //Convert both maps to ECL
   auto eclMapPtr = map.CoordinateConversion(CoordSys::ECL);
   auto eclFsMapPtr = fsMap.CoordinateConversion(CoordSys::ECL);

   //Make sure the values of the fsMap are the same as map
   eclFsMapPtr->ApplyFunction([](const double &a, const double &b){ CPPUNIT_ASSERT_DOUBLES_EQUAL(a, b, 1e-8); }, *eclMapPtr);

}
