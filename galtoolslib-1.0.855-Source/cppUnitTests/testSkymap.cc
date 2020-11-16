#include "testSkymap.h"

CPPUNIT_TEST_SUITE_REGISTRATION( testSkymap );

void testSkymap::setUp() {

   std::valarray<double> specMin1(1), specMin4(4), specMax1(1), specMax4(4), spec1(1), spec4(4), spec10(10);

   specMin1[0] = 10;
   specMax1[0] = 20;
   spec1[0] = 15;

   for (size_t i(0); i < specMin4.size(); ++i) {
      specMin4[i] = 10 + i*10;
      specMax4[i] = 20 + i*10;
      spec4[i] = 15 + i*10;
   }

   for (size_t i(0); i < spec10.size(); ++i)
      spec10[i] = 100+i*10;

   NestMapBinned1.Resize(4,specMin1,specMax1,NEST);
   NestMap1.Resize(4,spec1,NEST);
   NestMapBinned4.Resize(4,specMin4,specMax4,NEST);
   NestMap4.Resize(4,spec4,NEST);
   NestMap10.Resize(4,spec10,NEST);
   RingMap4.Resize(4,spec4,RING);
   
}

void testSkymap::tearDown() {
}

void testSkymap::assignment() {

   const double val = 2.3;

   NestMap4 = val;

   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val, NestMap4[i][j], 1e-8);

   NestMap4 = 0;

}

void testSkymap::addition() {

   //First test simple addition between equal files
   const double val = 1.2;
   NestMap4 = val;
   Skymap<double> map(NestMap4);

   NestMap4 += map;

   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2*val, NestMap4[i][j], 1e-8);

   //Now make sure adding maps of wrong type is not allowed.
   NestMapBinned1 = val;
   NestMap1 = val;
   NestMapBinned4 = val;
   NestMap10 = val;
   map.setSpectra(map.getSpectra()+1.0);

   std::cerr<<"Expect 5 warning messages now:"<<std::endl;
   NestMap4 += NestMapBinned1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2*val, NestMap4[i][j], 1e-8);

   NestMap4 += NestMap1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2*val, NestMap4[i][j], 1e-8);

   NestMap4 += NestMapBinned4;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2*val, NestMap4[i][j], 1e-8);

   NestMap4 += NestMap10;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2*val, NestMap4[i][j], 1e-8);

   NestMap4 += map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2*val, NestMap4[i][j], 1e-8);

   //Make sure adding RingMap to NestMap works as expected
   for (int i(0); i < RingMap4.Npix(); ++i)
      RingMap4[i] = 0.1*i;

   NestMap4 += RingMap4;
   for (int i(0); i < NestMap4.Npix(); ++i) {
      int ir = NestMap4.nest2ring(i);
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2*val+RingMap4[ir][j], NestMap4[i][j], 1e-8);
   }
   NestMap4 = 0;

   //Now test that changing the mask works as expected.
   const unsigned int originalMask = Skymap<double>::MathMask();

   //First test spectral values
   Skymap<double>::SetMathMask(SkymapDim::ORDERBIT | SkymapDim::SPSIZEBIT | SkymapDim::SPBINNEDBIT);
   NestMap4 += map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val, NestMap4[i][j], 1e-8);

   //Then binning
   Skymap<double>::SetMathMask(SkymapDim::SPSIZEBIT );
   NestMap4 += NestMapBinned4;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(2*val, NestMap4[i][j], 1e-8);

   //Finally spectral size
   Skymap<double>::SetMathMask(0);
   NestMap4 += NestMap1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3*val, NestMap4[i][j], 1e-8);

   //Test wrong order and wrong spectral shapes
   std::cerr<<"Expect 2 warning messages now:"<<std::endl;
   NestMap4 += NestMap10;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3*val, NestMap4[i][j], 1e-8);

   map.Resize(1, map.getSpectra(), NEST, val);
   NestMap4 += map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(3*val, NestMap4[i][j], 1e-8);

   NestMap1 = 0;
   NestMapBinned1 = 0;
   NestMap4 = 0;
   NestMapBinned4 = 0;
   NestMap10 = 0;
   RingMap4 = 0;

   Skymap<double>::SetMathMask(originalMask);
}

void testSkymap::subtraction() {

   //First test simple addition between equal files
   const double val = 1.2;
   NestMap4 = val;
   Skymap<double> map(NestMap4);

   NestMap4 -= map;

   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0*val, NestMap4[i][j], 1e-8);

   //Now make sure adding maps of wrong type is not allowed.
   NestMapBinned1 = val;
   NestMap1 = val;
   NestMapBinned4 = val;
   NestMap10 = val;
   map.setSpectra(map.getSpectra()+1.0);

   std::cerr<<"Expect 5 warning messages now:"<<std::endl;
   NestMap4 -= NestMapBinned1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0*val, NestMap4[i][j], 1e-8);

   NestMap4 -= NestMap1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0*val, NestMap4[i][j], 1e-8);

   NestMap4 -= NestMapBinned4;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0*val, NestMap4[i][j], 1e-8);

   NestMap4 -= NestMap10;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0*val, NestMap4[i][j], 1e-8);

   NestMap4 -= map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(0*val, NestMap4[i][j], 1e-8);

   //Make sure adding RingMap to NestMap works as expected
   for (int i(0); i < RingMap4.Npix(); ++i)
      RingMap4[i] = 0.1*i;

   NestMap4 -= RingMap4;
   for (int i(0); i < NestMap4.Npix(); ++i) {
      int ir = NestMap4.nest2ring(i);
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-RingMap4[ir][j], NestMap4[i][j], 1e-8);
   }
   NestMap4 = 0;

   //Now test that changing the mask works as expected.
   const unsigned int originalMask = Skymap<double>::MathMask();

   //First test spectral values
   Skymap<double>::SetMathMask(SkymapDim::ORDERBIT | SkymapDim::SPSIZEBIT | SkymapDim::SPBINNEDBIT);
   NestMap4 -= map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-val, NestMap4[i][j], 1e-8);

   //Then binning
   Skymap<double>::SetMathMask(SkymapDim::SPSIZEBIT );
   NestMap4 -= NestMapBinned4;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-2*val, NestMap4[i][j], 1e-8);

   //Finally spectral size
   Skymap<double>::SetMathMask(0);
   NestMap4 -= NestMap1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-3*val, NestMap4[i][j], 1e-8);

   //Test wrong order and wrong spectral shapes
   std::cerr<<"Expect 2 warning messages now:"<<std::endl;
   NestMap4 -= NestMap10;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-3*val, NestMap4[i][j], 1e-8);

   map.Resize(1, map.getSpectra(), NEST, val);
   NestMap4 -= map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(-3*val, NestMap4[i][j], 1e-8);

   NestMap1 = 0;
   NestMapBinned1 = 0;
   NestMap4 = 0;
   NestMapBinned4 = 0;
   NestMap10 = 0;
   RingMap4 = 0;

   Skymap<double>::SetMathMask(originalMask);
}

void testSkymap::multiplication() {

   //First test simple addition between equal files
   const double val = 1.2;
   NestMap4 = val;
   Skymap<double> map(NestMap4);

   NestMap4 *= map;

   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val, NestMap4[i][j], 1e-8);

   //Now make sure adding maps of wrong type is not allowed.
   NestMapBinned1 = val;
   NestMap1 = val;
   NestMapBinned4 = val;
   NestMap10 = val;
   map.setSpectra(map.getSpectra()+1.0);

   std::cerr<<"Expect 5 warning messages now:"<<std::endl;
   NestMap4 *= NestMapBinned1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val, NestMap4[i][j], 1e-8);

   NestMap4 *= NestMap1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val, NestMap4[i][j], 1e-8);

   NestMap4 *= NestMapBinned4;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val, NestMap4[i][j], 1e-8);

   NestMap4 *= NestMap10;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val, NestMap4[i][j], 1e-8);

   NestMap4 *= map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val, NestMap4[i][j], 1e-8);

   //Make sure adding RingMap to NestMap works as expected
   for (int i(0); i < RingMap4.Npix(); ++i)
      RingMap4[i] = 0.1*i;

   NestMap4 *= RingMap4;
   for (int i(0); i < NestMap4.Npix(); ++i) {
      int ir = NestMap4.nest2ring(i);
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val*RingMap4[ir][j], NestMap4[i][j], 1e-8);
   }
   NestMap4 = 1;

   //Now test that changing the mask works as expected.
   const unsigned int originalMask = Skymap<double>::MathMask();

   //First test spectral values
   Skymap<double>::SetMathMask(SkymapDim::ORDERBIT | SkymapDim::SPSIZEBIT | SkymapDim::SPBINNEDBIT);
   NestMap4 *= map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val, NestMap4[i][j], 1e-8);

   //Then binning
   Skymap<double>::SetMathMask(SkymapDim::SPSIZEBIT );
   NestMap4 *= NestMapBinned4;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val, NestMap4[i][j], 1e-8);

   //Finally spectral size
   Skymap<double>::SetMathMask(0);
   NestMap4 *= NestMap1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val*val, NestMap4[i][j], 1e-8);

   //Test wrong order and wrong spectral shapes
   std::cerr<<"Expect 2 warning messages now:"<<std::endl;
   NestMap4 *= NestMap10;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val*val, NestMap4[i][j], 1e-8);

   map.Resize(1, map.getSpectra(), NEST, val);
   NestMap4 *= map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val*val, NestMap4[i][j], 1e-8);

   NestMap1 = 0;
   NestMapBinned1 = 0;
   NestMap4 = 0;
   NestMapBinned4 = 0;
   NestMap10 = 0;
   RingMap4 = 0;

   Skymap<double>::SetMathMask(originalMask);
}

void testSkymap::division() {

   //First test simple addition between equal files
   const double val = 1.2;
   NestMap4 = val;
   Skymap<double> map(NestMap4);

   NestMap4 /= map;

   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, NestMap4[i][j], 1e-8);

   //Now make sure adding maps of wrong type is not allowed.
   NestMapBinned1 = val;
   NestMap1 = val;
   NestMapBinned4 = val;
   NestMap10 = val;
   map.setSpectra(map.getSpectra()+1.0);

   std::cerr<<"Expect 5 warning messages now:"<<std::endl;
   NestMap4 /= NestMapBinned1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, NestMap4[i][j], 1e-8);

   NestMap4 /= NestMap1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, NestMap4[i][j], 1e-8);

   NestMap4 /= NestMapBinned4;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, NestMap4[i][j], 1e-8);

   NestMap4 /= NestMap10;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, NestMap4[i][j], 1e-8);

   NestMap4 /= map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, NestMap4[i][j], 1e-8);

   //Make sure adding RingMap to NestMap works as expected
   for (int i(0); i < RingMap4.Npix(); ++i)
      RingMap4[i] = 2.3*0.1*i;

   NestMap4 /= RingMap4;
   for (int i(0); i < NestMap4.Npix(); ++i) {
      int ir = NestMap4.nest2ring(i);
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/RingMap4[ir][j], NestMap4[i][j], 1e-8);
   }
   NestMap4 = val;

   //Now test that changing the mask works as expected.
   const unsigned int originalMask = Skymap<double>::MathMask();

   //First test spectral values
   Skymap<double>::SetMathMask(SkymapDim::ORDERBIT | SkymapDim::SPSIZEBIT | SkymapDim::SPBINNEDBIT);
   NestMap4 /= map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, NestMap4[i][j], 1e-8);

   //Then binning
   Skymap<double>::SetMathMask(SkymapDim::SPSIZEBIT );
   NestMap4 /= NestMapBinned4;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/val, NestMap4[i][j], 1e-8);

   //Finally spectral size
   Skymap<double>::SetMathMask(0);
   NestMap4 /= NestMap1;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/val/val, NestMap4[i][j], 1e-8);

   //Test wrong order and wrong spectral shapes
   std::cerr<<"Expect 2 warning messages now:"<<std::endl;
   NestMap4 /= NestMap10;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/val/val, NestMap4[i][j], 1e-8);

   map.Resize(1, map.getSpectra(), NEST, val);
   NestMap4 /= map;
   for (int i(0); i < NestMap4.Npix(); ++i)
      for (int j(0); j < NestMap4.nSpectra(); ++j)
         CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/val/val, NestMap4[i][j], 1e-8);

   NestMap1 = 0;
   NestMapBinned1 = 0;
   NestMap4 = 0;
   NestMapBinned4 = 0;
   NestMap10 = 0;
   RingMap4 = 0;

   Skymap<double>::SetMathMask(originalMask);
}
