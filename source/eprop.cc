using namespace std;
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <gsl/gsl_sf_hyperg.h>
#define sign(a,b) (((b) > 0.) ? (a) : (-a))

double bessj(int, double);

//** see a sample main routine at the end **
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Propagation of electrons in the Galaxy cylindrical geometry
// INPUT:
// Ee         = electron kinetic energy, GeV
// R, Z       = (R,Z) point coordinates, kpc
// Rg, hg     = radius and 1/2 thickness of the area filled with sources, kpc
// hh         = halo size, kpc
// elossconst = energy loss constant, 1/(GeV s); can not be =0
// elossconst=32./9.*Pi*pow(Rele/Mele,2)*C*Ugevcc; Ugevcc = ISRF energy density
// gamma_e    = electron spectral index 
// Dxx, g1    = diffusion coefficient normalization @ 1 GeV (kpc^2/s), and index
// requires file j0zero500.dat with zeros of Bessel J0 function
// OUTPUT:
// Electron density @ (R,Z), 1/(cc GeV) 
// -calculated for a uniform distribution of sources (within R<Rg, |Z|<hg), 
// the source normalization is 1/(cc GeV) at 1 GeV; 
// the energy losses dE/dt= -elossconst*E^2
// REFERENCE: Bulanov, S.V., Dogiel, V.A. 1974, Astrophys. Spa. Sci. 29,305,
//            their eq.(8) has to be divided by 2pi (-perhaps an error)
//                                                  I.V.Moskalenko  11/01/2006
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double eprop(double Ee, double R, double Z, 
	     double Rg, double hg, double hh, double elossconst,
	     double gamma_e, double Dxx, double g1, const string &path)
{
  static int readarray=0;
  static double J0zero[500];
  int i, j, m=100, n=200;
  double x,S1,S2,S3,Pi,eps=1.e-16;

  Pi=acos(-1.);
#pragma omp single
  if(readarray==0)
    {
      ifstream data;
      const string fullFilename = path + "/j0zero500.dat";
      data.open(fullFilename.c_str());                    // open file if exists
      if(data.fail())
	{
	  cerr<<" >>eprop>> Error opening file "<<fullFilename<<endl;
	  exit(1);
	}
      for(i=0; i<500; data >> J0zero[i++]);
      data.close();
      readarray=1;
    }

  for(S2=1.,S3=0., i=0; i<m || S2>eps*S3; i++)
    {
      for(S1=0., j=0;j<n;j++)
	{
	  x=-(pow(Pi*(i+0.5),2)+pow(J0zero[j]*hh/Rg,2))
	    *Dxx*pow(Ee,g1-1)/(hh*hh*elossconst)/(1.-g1);
	  S1+= 
	    bessj(0,J0zero[j]*R/Rg) / bessj(1,J0zero[j]) / J0zero[j]*
            gsl_sf_hyperg_1F1(1., (gamma_e-g1)/(1.-g1), x);
	}
      S2=S1*sin(Pi*hg/hh*(i+0.5))*cos(Pi*Z/hh*(i+0.5))/(i+0.5);
      S3+=S2;
    }
  S3*=4.*pow(Ee,-gamma_e-1)/(Pi*(gamma_e-1.)*elossconst);
  return S3*1000.; // units 1/(cm3 GeV)
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Bessel J0, J1 functions: adaptations of CERNLIB C312 fortran routines 
//                                                  I.V.Moskalenko  11/01/2006
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double bessj(int order, double x)
{
  double A,B,F,P,Q,V;

  switch(order)
    {
    case 0: // bessj0
      V=fabs(x);
      if(V<8.)
	{
	  F=0.0625*x*x-2.;
	  A =           - 0.0000000000000008;
	  B = F * A     + 0.0000000000000413;
	  A = F * B - A - 0.0000000000019438;
	  B = F * A - B + 0.0000000000784870;
	  A = F * B - A - 0.0000000026792535;
	  B = F * A - B + 0.0000000760816359;
	  A = F * B - A - 0.0000017619469078;
	  B = F * A - B + 0.0000324603288210;
	  A = F * B - A - 0.0004606261662063;
	  B = F * A - B + 0.0048191800694676;
	  A = F * B - A - 0.0348937694114089;
	  B = F * A - B + 0.1580671023320973;
	  A = F * B - A - 0.3700949938726498;
	  B = F * A - B + 0.2651786132033368;
	  A = F * B - A - 0.0087234423528522;
	  A = F * A - B + 0.3154559429497802;
	  return 0.5*(A-B);
	}
      else
	{
	  F=256./(x*x)-2.;
	  B =           + 0.0000000000000007;
	  A = F * B     - 0.0000000000000051;
	  B = F * A - B + 0.0000000000000433;
	  A = F * B - A - 0.0000000000004305;
	  B = F * A - B + 0.0000000000051683;
	  A = F * B - A - 0.0000000000786409;
	  B = F * A - B + 0.0000000016306465;
	  A = F * B - A - 0.0000000517059454;
	  B = F * A - B + 0.0000030751847875;
	  A = F * B - A - 0.0005365220468132;
	  A = F * A - B + 1.9989206986950373;
	  P=A-B;
	  B =           - 0.0000000000000006;
	  A = F * B     + 0.0000000000000043;
	  B = F * A - B - 0.0000000000000334;
	  A = F * B - A + 0.0000000000003006;
	  B = F * A - B - 0.0000000000032067;
	  A = F * B - A + 0.0000000000422012;
	  B = F * A - B - 0.0000000007271916;
	  A = F * B - A + 0.0000000179724572;
	  B = F * A - B - 0.0000007414498411;
	  A = F * B - A + 0.0000683851994261;
	  A = F * A - B - 0.0311117092106740;
	  Q=8.0*(A-B)/V;
	  F=V-0.785398163397448;
	  A=cos(F);
	  B=sin(F);
	  F=0.398942280401432/sqrt(V);
	  return F*(P*A-Q*B);
	}

    case 1: // bessj1
      V=fabs(x);
      if(V<8.)
	{
	  F=0.0625*x*x-2.;
	  B =           + 0.0000000000000114;
	  A = F * B     - 0.0000000000005777;
	  B = F * A - B + 0.0000000000252812;
	  A = F * B - A - 0.0000000009424213;
	  B = F * A - B + 0.0000000294970701;
	  A = F * B - A - 0.0000007617587805;
	  B = F * A - B + 0.0000158870192399;
	  A = F * B - A - 0.0002604443893486;
	  B = F * A - B + 0.0032402701826839;
	  A = F * B - A - 0.0291755248061542;
	  B = F * A - B + 0.1777091172397283;
	  A = F * B - A - 0.6614439341345433;
	  B = F * A - B + 1.2879940988576776;
	  A = F * B - A - 1.1918011605412169;
	  A = F * A - B + 1.2967175412105298;
	  return 0.0625*(A-B)*x;
	}
      else
	{
	  F=256./(x*x)-2.;
	  B =           - 0.0000000000000007;
	  A = F * B     + 0.0000000000000055;
	  B = F * A - B - 0.0000000000000468;
	  A = F * B - A + 0.0000000000004699;
	  B = F * A - B - 0.0000000000057049;
	  A = F * B - A + 0.0000000000881690;
	  B = F * A - B - 0.0000000018718907;
	  A = F * B - A + 0.0000000617763396;
	  B = F * A - B - 0.0000039872843005;
	  A = F * B - A + 0.0008989898330859;
	  A = F * A - B + 2.0018060817200274;
	  P=A-B;
	  B =           + 0.0000000000000007;
	  A = F * B     - 0.0000000000000046;
	  B = F * A - B + 0.0000000000000360;
	  A = F * B - A - 0.0000000000003264;
	  B = F * A - B + 0.0000000000035152;
	  A = F * B - A - 0.0000000000468636;
	  B = F * A - B + 0.0000000008229193;
	  A = F * B - A - 0.0000000209597814;
	  B = F * A - B + 0.0000009138615258;
	  A = F * B - A - 0.0000962772354916;
	  A = F * A - B + 0.0935555741390707;
	  Q=8.*(A-B)/V;
	  F=V-2.356194490192345;
	  A=cos(F);
	  B=sin(F);
	  F=0.398942280401432/sqrt(V);
	  return( x>0. ? F*(P*A-Q*B): -F*(P*A-Q*B));
	}
    default:
      cout<<" >>bessj>> routine has been called with order="<<order<<endl
	  <<" >>bessj>> routine calculates Bessel function of orders 0 and 1"<<endl;
      exit(1);
    }
}

