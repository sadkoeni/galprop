
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * create_gcr.cc *                               galprop package * 08/16/2001 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>
#include <string>

#include <ErrorLogger.h>

#include <Nuclei_Interface.h>
//using namespace std;


int Galprop::create_gcr()//AWS20050816
{
  //cout<<" >>>> create_gcr"<<endl;
  INFO("Entry"); 
  if(galdef.verbose >= 1)galdef.print();

   int i=0,j,Z,A,Z1,A1, stat=0;
   std::string name;
   //string name;
   static const char *element[]=
      {
        "Hydrogen",    "Helium",   "Lithium", "Beryllium",     "Boron",
          "Carbon",  "Nitrogen",    "Oxygen",  "Fluorine",      "Neon",
	  "Sodium", "Magnesium",  "Aluminium",   "Silicon","Phosphorus",
          "Sulphur",  "Chlorine",     "Argon", "Potassium",   "Calcium",
        "Scandium",  "Titanium",  "Vanadium",  "Chromium", "Manganese",
            "Iron",    "Cobalt",    "Nickel",    "Copper",      "Zinc"
      };
   int galdef_network_par=0;         // imos network, use everywhere
   int K_electron;                                                     // AWS20010731

// calculate the number of species
   if(!galdef.secondary_antiprotons) galdef.tertiary_antiprotons=0;    // IMOS20000802

   n_species=0;
   if(galdef.secondary_positrons || galdef.pair_production)   n_species++;
   if(galdef.knock_on_electrons)    n_species++; //IMOS20060504
   if(galdef.secondary_electrons || galdef.pair_production)   n_species++;
   if(galdef.primary_electrons)     n_species++;
   if(galdef.primary_positrons)     n_species++;

// DM decay species: IMOS20050912
   if(galdef.DM_positrons)          n_species++; 
   if(galdef.DM_electrons)          n_species++;
   if(galdef.DM_antiprotons)        n_species++;

   if(galdef.tertiary_antiprotons)  n_species++;                       // IMOS20000605
   if(galdef.secondary_antiprotons) n_species++;
   if(galdef.secondary_protons)     n_species++;                       // IMOS20000605
   for(Z=1; Z<=galdef.max_Z; Z++)
   {
      if(!galdef.use_Z[Z]) continue;
      for(A=2*Z-2; A<2.5*Z+4.2; A++)                                   // IMOS20010816 whole loop
      {
         double t_half[2];
 	 if(!nucdata(galdef_network_par,Z,A,0,Z,A,&Z1,&A1,&t_half[0])) continue;
 	 if(!nucdata(galdef_network_par,Z,A,1,Z,A,&Z1,&A1,&t_half[1])) continue;
	 for(K_electron=0;K_electron<=galdef.K_capture;K_electron++)
	 {
            if(t_half[K_electron]>0. && t_half[K_electron]/year2sec<galdef.t_half_limit) continue;
            if(K_electron) if(t_half[0] == t_half[1]) continue;
            n_species++;
         }
      }
   }

   if(galdef.verbose >= 1)cout<<endl<<"Number of species to create: n_species= "<<n_species<<endl<<endl;
   if(!n_species) {cout<<"create_gcr.cc: No particles specified -exit"<<endl; exit(1);}

// create a Particle array
   delete[] gcr;
   gcr=new Particle[n_species];
   K_electron=0; // for all secondaries and non-nuclei                    AWS20010731

// DM POSITRONS   IMOS20050912

   if(galdef.DM_positrons)
   {
      name = "DM_positrons";
      Z=1;  A=0;  
      const double t_half(0);

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half,                               // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half,                               // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print(); 
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// DM ELECTRONS   IMOS20050912

   if(galdef.DM_electrons)
   {
      name = "DM_electrons";
      Z=-1;  A=0;
      const double t_half(0);

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half,                               // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half,                               // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print(); 
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// SECONDARY POSITRONS

   if (galdef.secondary_positrons || galdef.pair_production) {

     name = "secondary_positrons";
     Z = 1;  
     A = 0;  
     const double t_half = 0;                                                  // IMOS20010816
     
     gcr[i].primary_abundance = 0;
     
     if (2 == galdef.n_spatial_dimensions)
       gcr[i].init(name, 
		   Z, A, 
		   K_electron, t_half, 
		   galdef.r_min, galdef.r_max, galdef.dr,  
		   galdef.z_min, galdef.z_max, galdef.dz,
		   galdef.p_min, galdef.p_max, galdef.p_factor,
		   galdef.Ekin_min, galdef.Ekin_max, galdef.Ekin_factor,
		   galdef.p_Ekin_grid);  
     
      if (3 == galdef.n_spatial_dimensions)
	gcr[i].init(name,
		    Z, A, 
		    K_electron, t_half,
		    galdef.x_min, galdef.x_max, galdef.dx,  
		    galdef.y_min, galdef.y_max, galdef.dy,
		    galdef.z_min, galdef.z_max, galdef.dz,
		    galdef.p_min, galdef.p_max, galdef.p_factor,
		    galdef.Ekin_min, galdef.Ekin_max, galdef.Ekin_factor,
		    galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print(); 
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// KNOCK-ON ELECTRONS

   if(galdef.knock_on_electrons)                                         // IMOS20060504
   {
      name = "knock_on_electrons";
      Z=-1;  A=0;
      const double t_half(0);

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// SECONDARY ELECTRONS

   if (galdef.secondary_electrons || galdef.pair_production) {
     
     name = "secondary_electrons";
     Z = -1;  
     A = 0;  
     const double t_half = 0;                                         // IMOS20010816

     gcr[i].primary_abundance = 0;

     if (2 == galdef.n_spatial_dimensions)
       gcr[i].init(name,
		   Z, A, 
		   K_electron, t_half,
		   galdef.r_min, galdef.r_max, galdef.dr,  
		   galdef.z_min, galdef.z_max, galdef.dz,
		   galdef.p_min, galdef.p_max, galdef.p_factor,
		   galdef.Ekin_min, galdef.Ekin_max, galdef.Ekin_factor,
		   galdef.p_Ekin_grid);  

     if (3 == galdef.n_spatial_dimensions)
       gcr[i].init(name,
		   Z, A, 
		   K_electron, t_half,
		   galdef.x_min, galdef.x_max, galdef.dx,  
		   galdef.y_min, galdef.y_max, galdef.dy,
		   galdef.z_min, galdef.z_max, galdef.dz,
		   galdef.p_min, galdef.p_max, galdef.p_factor,
		   galdef.Ekin_min, galdef.Ekin_max, galdef.Ekin_factor,
		   galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// PRIMARY ELECTRONS

   if(galdef.primary_electrons)
   {
      name = "primary_electrons";
      Z=-1;  A=0;
      const double t_half(0);

      gcr[i].primary_abundance=1.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx, 
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print(); 
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// PRIMARY Positrons

   if(galdef.primary_positrons)
   {
      name = "primary_positrons";
      Z=1;  A=0;
      const double t_half(0);

      gcr[i].primary_abundance=1.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx, 
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print(); 
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// DM ANTIPROTONS   IMOS20050912

   if(galdef.DM_antiprotons)
   {
      name = "DM_antiprotons";
      Z=-1;  A=1;
      const double t_half(0);

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// TERTIARY ANTIPROTONS

   if(galdef.tertiary_antiprotons)
   {
      name = "tertiary_antiprotons";
      Z=-1;  A=1;
      const double t_half(0);

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// SECONDARY ANTIPROTONS

   if(galdef.secondary_antiprotons)
   {
      name = "secondary_antiprotons";
      Z=-1;  A=1;
      const double t_half(0);

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// NUCLEONS: SECONDARY PROTONS

   if(galdef.secondary_protons)
   {
      name = "secondary_protons";
      Z=1;  A=1;
      const double t_half(0);

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half,                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// OTHER NUCLEONS

   name = "particle";

   for(Z=1; Z<=galdef.max_Z; Z++)
   {
      if(!galdef.use_Z[Z]) continue;
      for(A=2*Z-2; A<2.5*Z+4.2; A++)
      {
         double t_half[2];
 	 if(!nucdata(galdef_network_par,Z,A,0,Z,A,&Z1,&A1,&t_half[0])) continue;          // IMOS20010816
 	 if(!nucdata(galdef_network_par,Z,A,1,Z,A,&Z1,&A1,&t_half[1])) continue;          // IMOS20010816
         for(K_electron=galdef.K_capture;K_electron>=0;K_electron--)                      // AWS20010731
	 {
            if(galdef.verbose >= 1)cout<<"nucleus being tested Z A K_electron "<<Z<<" "<<A<<" "<<K_electron<<endl;//AWS20010731
	                                                                                  // IMOS20010816 next 2 lines
            if(t_half[K_electron]>0. && t_half[K_electron]/year2sec<galdef.t_half_limit) continue;
            if(K_electron) if(t_half[0] == t_half[1]) continue;

            t_half[K_electron]/= year2sec;
            if(galdef.verbose >= 1)cout<<"nucleus used Z A t_half "<<Z<<" "<<A<<" "<<t_half[K_electron]<<endl;
            std::ostringstream oss;
            if (Z > 0 && Z < 31) oss << element[Z-1] <<"_"<<A;
            name = oss.str();

            gcr[i].primary_abundance=0.0; //galdef.isotopic_abundance[Z][A];  This is now handled by the source classes
            //if(K_electron >0) gcr[i].primary_abundance= 0.0;   // no prim. K-capture nuclei AWS20010731
     
            if(K_electron==1) name += "*";

            if(galdef.n_spatial_dimensions==2)
               gcr[i].init(name,Z,A,K_electron,t_half[K_electron],                      // IMOS20010816
                  galdef.r_min,  galdef.r_max, galdef.dr,  
                  galdef.z_min,  galdef.z_max, galdef.dz,
                  galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
                  galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
                  galdef.p_Ekin_grid);  

            if(galdef.n_spatial_dimensions==3)
               gcr[i].init(name,Z,A,K_electron,t_half[K_electron],                      // IMOS20010816
                  galdef.x_min,  galdef.x_max, galdef.dx,  
                  galdef.y_min,  galdef.y_max, galdef.dy,
                  galdef.z_min,  galdef.z_max, galdef.dz,
                  galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
                  galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
                  galdef.p_Ekin_grid); 

            if(galdef.verbose >= 1)gcr[i].print();
            if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
            i++;
	 } // K_electron
      } //A
   } //Z


   //Loop over the species and calculate dependencies, use the propel loop and the gen_secondary_source algorithm.
//#pragma omp parallel for default(shared) schedule(dynamic) private(i,j,Z,A,Z1,A1,K_electron)
   for (int i = n_species-1; i >= 0; --i) {

      //Check for special particles, i.e. everything that is not a nuclei
      if (gcr[i].name == "primary_electrons" ||
          gcr[i].name == "primary_positrons" ||
          gcr[i].name == "DM_positrons" ||
          gcr[i].name == "DM_electrons" ||
          gcr[i].name == "DM_antiprotons" )
         continue;

      if ((gcr[i].name == "secondary_electrons" && galdef.secondary_electrons) ||
          (gcr[i].name == "secondary_positrons" && galdef.secondary_positrons) ||
           gcr[i].name == "knock_on_electrons" ||
           gcr[i].name == "secondary_antiprotons") 
      {

         //Look for proton and helium particles
         for ( int j = 0; j < n_species; ++j ) {
            
            if ( (gcr[j].A == 1 && gcr[j].Z == 1) || (gcr[j].A == 4 && gcr[j].Z == 2) ) 
               gcr[i].dependencies[&gcr[j]] = std::valarray<double>(0); //Not the cross section, just for dependency calculation

         }

         //In the special case of secondary_antiprotons with option 3
         if (galdef.secondary_antiprotons == 3) {
            //Add all nuclei
            for ( int j = 0; j < n_species; ++j ) {
               if ( (gcr[j].A > 0 && gcr[j].Z > 0 && gcr[j].A < 61) ) 
                  gcr[i].dependencies[&gcr[j]] = std::valarray<double>(0); //Not the cross section, just for dependency calculation
            }
         }
      }

         if ( gcr[i].name == "secondary_protons" ) {
         
            //Look for proton and helium particles
            for ( int j = 0; j < n_species; ++j ) {
               
               if ( gcr[j].A == 1 && gcr[j].Z == 1 && j != i ) 
                  gcr[i].dependencies[&gcr[j]] = std::valarray<double>(0); //Not the cross section, just for dependency calculation

            }
         }

         if ( gcr[i].name == "tertiary_antiprotons" ) {

            //Look for secondary antiprotons
            for ( int j = 0; j < n_species; ++j ) {

               if (j == i)
                  continue;
               
               if ( gcr[j].A == 1 && gcr[j].Z == -1 ) 
                  gcr[i].dependencies[&gcr[j]] = std::valarray<double>(0); //Not the cross section, just for dependency calculation

            }
         }

         //Now the nuclei
         std::valarray<double> cross_section(gcr[i].n_pgrid+2);
         int Z1, A1;
         double t_half;
         
         for (j = 0; j < n_species; ++j) {

            if (gcr[j].A < 1)
               continue;

            cross_section = 0;

            if (gcr[j].A > gcr[i].A && !gcr[i].K_electron && galdef.fragmentation ) {

               decayed_cross_sections(gcr[j].Z,gcr[j].A, gcr[i].Z,gcr[i].A, &gcr[i].Ekin[0], gcr[i].n_pgrid, &cross_section[0]);

               if (galdef.verbose==-602 ||galdef.verbose>=1) { // selectable debug                  //AWS20010828

               ostringstream buf;
               buf << "Cross_section for Z, A = " <<gcr[j].Z << ", " <<gcr[j].A << "-> " << gcr[i].Z << ", " <<gcr[i].A << " :";
               for (int ip = 0; ip < gcr[0].n_pgrid; buf << cross_section[ip++] << " "); 
               INFO(buf.str());

            }

            // apply factors as shown in galprop notes

            cross_section *= double(gcr[j].A)/gcr[i].A;

            if (galdef.verbose==-602 ||galdef.verbose>=1) { // selectable debug                  //AWS20010828

               ostringstream buf;
               buf << "Cross_section with Aprim/Asec factors: ";
               for (int ip = 0; ip < gcr[0].n_pgrid; buf << cross_section[ip++] << " ");
               INFO(buf.str());

            }

            //Add to dependencies if cross section is > 0 for any energy
            for (int ip = 0; ip < gcr[i].n_pgrid; ++ip) {
               
               if (cross_section[ip] > 0) {

                  gcr[i].dependencies[&gcr[j]] = std::valarray<double>(0);
                  gcr[i].dependencies[&gcr[j]].resize(cross_section.size());
                  gcr[i].dependencies[&gcr[j]] = cross_section;
                  break;

               }

            }
               

         }  //  gcr.A>particle.A

         if (1 == galdef.K_capture) { //AWS20010828

            // ELECTRON ATTACHMENT or STRIPPING SOURCE   IMOS20010816

            double attach_H, strip_H, attach_He, strip_He;
            int zH = 1, zHe = 2;

            // Source from isotopes quickly decaying after EC 
            if (gcr[j].A == gcr[i].A && gcr[j].Z != gcr[i].Z)
               if (nucdata(galdef_network_par,
                        gcr[j].Z,
                        gcr[j].A,
                        1,
                        gcr[i].Z,
                        gcr[i].A, 
                        &Z1,
                        &A1,
                        &t_half)) {

                  if (t_half > 0. && t_half/year2sec < galdef.t_half_limit) {

                     for (int ip = 0; ip<gcr[i].n_pgrid; ++ip) {

                        nuclei::Kcapture_cs(gcr[j].Ekin[ip],gcr[j].Z,zH, &attach_H ,&strip_H );
                        nuclei::Kcapture_cs(gcr[j].Ekin[ip],gcr[j].Z,zHe,&attach_He,&strip_He);
                        cross_section[ip] = attach_H + galdef.He_H_ratio*attach_He;

                     }

                     if (galdef.verbose==-602 ||galdef.verbose>=1) { // selectable debug                  //AWS20010828

                        ostringstream buf;
                        buf << "gen_secondary_source: electron attachment and fast decay cross_section for Z,A ="
                           <<gcr[j].Z<<","<<gcr[j].A<<" K_electron="<<gcr[j].K_electron
                           <<"->>"
                           <<gcr[i].Z<<","<<gcr[i].A<<"  K_electron="<<gcr[i].K_electron
                           <<" t_half="<<t_half/year2sec<<" yr : ";
                        for(int ip=0; ip<gcr[0].n_pgrid; buf << cross_section[ip++]<<" "); 
                        INFO(buf.str());

                     }

                     //Add to dependencies if cross section is > 0 for any energy
                     for (int ip = 0; ip < gcr[i].n_pgrid; ++ip) {

                        if (cross_section[ip] > 0) {

                           gcr[i].dependencies[&gcr[j]] = std::valarray<double>(cross_section.size());
                           gcr[i].dependencies[&gcr[j]].resize(cross_section.size());
                           gcr[i].dependencies[&gcr[j]] = cross_section;
                           break;

                        }

                     }

                  }

               }

            // Source from electron attachment or stripping (for slow EC decays) 
            if (100*gcr[j].Z+gcr[j].A == 100*gcr[i].Z+gcr[i].A)
               if (gcr[j].K_electron != gcr[i].K_electron) {

                  for (int ip = 0; ip < gcr[0].n_pgrid; ++ip) {

                     nuclei::Kcapture_cs(gcr[j].Ekin[ip],gcr[j].Z,zH, &attach_H ,&strip_H );
                     nuclei::Kcapture_cs(gcr[j].Ekin[ip],gcr[j].Z,zHe,&attach_He,&strip_He);

                     if(gcr[i].K_electron) 
                        cross_section[ip] = attach_H + galdef.He_H_ratio*attach_He;
                     else                    
                        cross_section[ip] = strip_H + galdef.He_H_ratio*strip_He ;

                  }

                  if (galdef.verbose==-602||galdef.verbose>=1) { // selectable debug                  //AWS20010828

                     ostringstream buf;

                     if (gcr[i].K_electron) 
                        buf<<"gen_secondary_source: electron attachment for Z,A =";
                     else
                        buf<<"gen_secondary_source: stripping cross_section for Z,A =";

                     buf <<gcr[j].Z<<","<<gcr[j].A<<" K_electron="<<gcr[j].K_electron<<"->"
                        <<gcr[i].Z<<","<<gcr[i].A<<"  K_electron="<<gcr[i].K_electron<<" :";
                     for (int ip=0; ip<gcr[0].n_pgrid; buf<<cross_section[ip++]<<" "); 

                     INFO(buf.str());

                  }

                  //Add to dependencies if cross section is > 0 for any energy
                  for (int ip = 0; ip < gcr[i].n_pgrid; ++ip) {

                     if (cross_section[ip] > 0) {

                        if ( gcr[i].K_electron == 0 ) {
                           std::ostringstream ost;
                           ost << "Circular dependency for K_electron attachment and stripping for "<<gcr[i].name;
                           ost << ".  A few network iterations will be needed to get an accurate solution.";
                           if (galdef.network_iter_compl < 2) {
                              ost << "  Increase \"network_iter_compl\" to at least 2 for an accurate solution.";
                              WARNING(ost.str());
                           } else {
                              INFO(ost.str());
                           }
                        }

                        gcr[i].dependencies[&gcr[j]] = std::valarray<double>(cross_section.size());
                        gcr[i].dependencies[&gcr[j]].resize(cross_section.size());
                        gcr[i].dependencies[&gcr[j]] = cross_section;
                        break;

                     }

                  }

               } // gcr[i].K_electron != particle.K_electron

         } //galdef.K_capture==1

         // RADIOACTIVE DECAY SOURCE (beta+/-, EC) AWS20010806 IMOS20010816

         //double branching_ratio;

         if (gcr[j].A >= gcr[i].A && galdef.radioactive_decay) {

            if (100*gcr[j].Z+gcr[j].A != 100*gcr[i].Z+gcr[i].A) {

               const double branching_ratio = nucdata(galdef_network_par,gcr[j].Z,  
                     gcr[j].A, gcr[j].K_electron,
                     gcr[i].Z,gcr[i].A, 
                     &Z1, &A1, &t_half);

               if (branching_ratio != 0 && t_half != 0. && 100*Z1+A1 == 0) {

                  if(galdef.verbose==-602||galdef.verbose>=1) { // selectable debug                  //AWS20010828

                     ostringstream buf;
                     buf <<"gen_secondary_source: (Z,A,K_electron) "
                        << gcr[j].Z<<" "<<  gcr[j].A<<" "<<  gcr[j].K_electron<<"->"
                        <<gcr[i].Z<<" "<<gcr[i].A
                        <<" t_half="<<t_half/year2sec<<" yr  branching ratio="<<branching_ratio;
                     INFO(buf.str());

                  }

                  // Make sure the particle is not already in dependencies
                  std::map<Particle*, std::valarray<double> >::iterator it = gcr[i].dependencies.find(&gcr[j]);
                  
                  if ( it == gcr[i].dependencies.end() ) {

                     //Add it with 0 cross sections but branching ratio and half time
                     cross_section[gcr[i].n_pgrid] = branching_ratio;
                     cross_section[gcr[i].n_pgrid+1] = t_half;

                     gcr[i].dependencies[&gcr[j]] = std::valarray<double>(cross_section.size());
                     gcr[i].dependencies[&gcr[j]].resize(cross_section.size());
                     gcr[i].dependencies[&gcr[j]] = cross_section;

                  } else {

                     it->second[gcr[i].n_pgrid] = branching_ratio;
                     it->second[gcr[i].n_pgrid+1] = t_half;

                  }


               } // branching_ratio != 0 && t_half != 0. && 100*Z1+A1 == 0

            } // 100*gcr[i].Z+gcr[i].A != 100*particle.Z+particle.A

         } // gcr[i].A >= particle.A

      }  //  species

   }

   //Check for circular dependencies
   bool circularDependency(false);
   for (int i = 0; i < n_species; ++i) {

      std::map<Particle*, std::valarray<double> >::iterator it;
      for ( it = gcr[i].dependencies.begin(); it != gcr[i].dependencies.end(); ++it) {

         std::map<Particle*, std::valarray<double> >::iterator it2;
         for ( it2 = it->first->dependencies.begin(); it2 != it->first->dependencies.end(); ++it2) {

            if ( it2->first == &gcr[i] ) {

               ostringstream buf;
               buf << "Circular dependency detected: "<<gcr[i].name<<" depends on itself through "<<it->first->name;

               //Do not abort unless we are not dealing with K_electron
               if ( ! ( gcr[i].A == it->first->A && gcr[i].Z == it->first->Z && gcr[i].K_electron != it->first->K_electron ) ) {
                  ERROR(buf.str());
                  circularDependency = true;
               } else {
                  INFO(buf.str());
               }
               
            }
         }
      }
   }

   if (circularDependency)
      return 1;
      

   //Print the dependencies, make sure they agree with the debugging output from above

   if (galdef.verbose==-602 ||galdef.verbose>=1) { // selectable debug                  //AWS20010828
      
      for (int i = n_species-1; i >= 0; --i) {

         ostringstream buf;

         buf << "Number of dependencies for "<<gcr[i].name<<" is "<<gcr[i].dependencies.size()<<":"<<std::endl;

         //Loop over the dependencies, print out relevant information
         std::map<Particle*, std::valarray<double> >::iterator it;
         for ( it = gcr[i].dependencies.begin(); it != gcr[i].dependencies.end(); ++it) {

            buf <<"  "<<it->first->name<<": ";
            for (size_t j = 0; j < it->second.size(); ++j)
               buf << it->second[j]<< " ";
            buf << std::endl;

         }

         INFO(buf.str());

      }

   }


   INFO("Exit");
   //cout<<" <<<< create_gcr"<<endl;
   return stat;
}




