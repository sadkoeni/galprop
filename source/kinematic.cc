
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * kinematic.cc *                                galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|


// kinematic computation:
// Z        input         nuclear charge
// A        input         nuclear mass number
// mass     input         mass in MeV
// p        input/output  total momentum in MV
// Ekin     input/output  kinetic energy PER NUCLEON in MeV
// Etot     output        total energy of nucleus in MeV
// beta     output        velocity/c 
// gamma    output        total energy/ total mass
// rigidity output        p/Z  in MV 
//
// On input, if p   <0 then Ekin is used to compute all quantities
//          if Ekin<0 then p    is used to compute all quantities
// If both > 0 , error condition -1 returned.
// If both < 0 , error condition -2 returned.

using namespace std;//AWS20050624
#include"constants.h"
#include<cmath>
#include<iostream>
#include<string> //IMOS20020112
//#include <cstring>

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int kinematic(int Z, int A, const double m0, 
	      double &p, double &Ekin,double &Etot, 
	      double &beta, double &gamma, double &rigidity,int test) {

  Etot =beta =gamma =rigidity =-1.;  // preset values in case of error condition
  
  if(p>=0.0 && Ekin>=0.0)
    {  cout<<"kinematic: error:p or Ekin to be evaluated ?"<<endl; return -1; } // one of them has to be <0 
  if(p<=0.0 && Ekin<=0.0)
    {  cout<<"kinematic: error:p or Ekin has to be >0     "<<endl; return -2; } // one of them has to be >0 
  
  if (A == 0) A = 1;

// P input, compute Ekin etc.
   if(Ekin<0.)
   {
      Etot=sqrt(p*p+m0*m0);
      Ekin=(Etot-m0)/A;
   }

// Ekin input, compute P etc.
   if(p<0.)
   {
      Etot=Ekin*A+m0;
      p=sqrt(Etot*Etot-m0*m0);
   }

   gamma=Etot/m0;
   beta=sqrt((1.-1./gamma)*(1.+1./gamma));
   rigidity=p/fabs(1.*Z);
 
   if(test)
   {
      cout<<1.-1./(gamma*gamma)-beta*beta<<"   test: = 0?"<<endl;
      cout<<Etot*Etot-p*p-m0*m0          <<"   test: = 0?"<<endl;
   }
   return 0;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int test_kinematic()
{
   cout<<">>>>test_kinematic"<<endl;

   double p, Ekin,Etot, beta, gamma, rigidity;
   int Z=6, A=12; 
   const double mu=931.5; //Atomic mass unit in MeV
   double m0=A*mu;
   int test=0;

   double p0=1.e6;

   //strcpy(species,"nucleus");
   for(int i=0; i<10; i++)
   {
      p=p0;
      Ekin=-1;
      kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test); 
      cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
         <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<" p-p0="<<p-p0<<endl;
      p=-1; // reverse transformation
      kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test);
      cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
         <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<" p-p0="<<p-p0<<endl;
      Z++; A++;
      m0+=mu;
   }

   //strcpy(species,"electron");
   Z=1;
   A=0;
   p0=1.e0;
   m0=m_electron;
   for(int i=0; i<10; i++)
   {
      p=p0;
      Ekin=-1;
      kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test);
      cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
         <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<" p-p0="<<p-p0<<endl;
      p=-1; // reverse transformation
      kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test);
      cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
         <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<" p-p0="<<p-p0<<endl;
      Z--;
   }

   p=-1;
   test=1;
   kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test);
   cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
      <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<" p-p0="<<p-p0<<endl;
   test=0;
   cout<<"\n------------------\n";
   //strcpy(species,"nucleus");
   
   Z=6;
   A=12;
   m0=A*mu;
   for (p=.001; p<1e6; p*=2)
   {
      Ekin=-1;
      kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test);
      cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
         <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<endl;
   }
   cout<<"\n------------------\n";
   //strcpy(species,"nucleus");
   for (Ekin=.001; Ekin<1e6; Ekin*=2)
   {
      p=-1;
      kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test);
      cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
         <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<endl;
   }
   cout<<"\n------------------\n";
   //strcpy(species,"electron");
   Z=-1;
   A=0;
   m0=m_electron;
   for (p=.001; p<1e6; p*=2)
   {
      Ekin=-1;
      kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test);
      cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
         <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<endl;
   }
   cout<<"\n------------------\n";
   //strcpy(species,"electron");
   for (Ekin=.001;Ekin<1e6;Ekin*=2)
   {
      p=-1;
      kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test);
      cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
         <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<endl;
   }
   cout<<"                             ------- Testing error conditions ---------"<<endl;
   Ekin=1; p=1;
   kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test);
   cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
      <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<" p-p0="<<p-p0<<endl;
   cout<<"<<<<test_kinematic"<<endl;
   //strcpy(species,"nucleus");
   Ekin=-1; p=-1;
   kinematic(Z, A, m0, p, Ekin,Etot, beta, gamma, rigidity,test);
   cout<<"Z="<<Z<<" A="<< A <<" m0="<<m0 <<" p="<<p<<" Ekin="<< Ekin <<" Etot="<<Etot
      <<" beta="<< beta<<" gamma="<< gamma<<" rigidity="<< rigidity<<" p-p0="<<p-p0<<endl;
   cout<<"<<<<test_kinematic"<<endl;
   return 0;
}
