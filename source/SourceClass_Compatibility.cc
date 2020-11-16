#include "SourceClass_Compatibility.h"
#include "Galdef.h"
#include "galprop_internal.h"
#include <utility>
#include <iomanip>
#include <Registry.h>


SourceClass_Compatibility::SourceClass_Compatibility ( utl::Parameters &&pars ) :
   SourceClass(std::move(pars)),
   removeSteadyState(false),
   createSNRDistributions(false)
{
   //Set default values for parameters
   setDefaultParameters();

   //Read in all the parameters
   setPars(fpars);

}

void SourceClass_Compatibility::setDefaultParameters()
{
   //Spectrum type
   inj_spectrum_type    ="rigidity";

   //Nuclear injection spectrum
   spDefault.rigid_br0 = 0;
   spDefault.rigid_br = 9.0e3;
   spDefault.g_0 = 0;
   spDefault.g_1 = 1.82;
   spDefault.g_2 = 2.36;

   //Electron injection spectrum
   spElectron.rigid_br0 = 4.0e3;
   spElectron.rigid_br = 1.0e9;
   spElectron.g_0 = 1.6;
   spElectron.g_1 = 2.5;
   spElectron.g_2 = 5.0;

   //The model
   source_specification = 0;

   source_model         = 0;
   source_model_electron = 0;

   source_parameters.resize(10);
   source_parameters_electron.resize(10);
   for (int iS = 0; iS < 10; ++iS) {

      source_parameters[iS] = source_parameters_electron[iS] = 0;


   }
   source_parameters[0] = source_parameters_electron[0] = 0.2;
   source_parameters[4] = source_parameters_electron[4] = 100;

   source_xmlFile       = "";

   n_cr_sources         = 0;

   SNR_events           = 0;

   proton_norm_Ekin     = 1.00e+5;
   spectra_norm_Rigidity= -1.0;  //Default to normalization on Ekin

   //Assume we have normalization to flux by default
   source_normalization = 1;
   electron_source_normalization = 1;

}

void SourceClass_Compatibility::setPars( const utl::Parameters &pars )
{

   fpars.setParameters(pars);
   
   //Spectrum type
   try { fpars.getParameter("inj_spectrum_type", inj_spectrum_type); } catch (utl::Parameters::ParameterError) { fpars.setParameter("inj_spectrum_type", inj_spectrum_type); }

   //Nuclear injection spectrum
   try { fpars.getParameter("nuc_rigid_br0", spDefault.rigid_br0); } catch (utl::Parameters::ParameterError) { fpars.setParameter("nuc_rigid_br0", spDefault.rigid_br0); }
   try { fpars.getParameter("nuc_rigid_br", spDefault.rigid_br); } catch (utl::Parameters::ParameterError) { fpars.setParameter("nuc_rigid_br", spDefault.rigid_br); }
   try { fpars.getParameter("nuc_g_0", spDefault.g_0); } catch (utl::Parameters::ParameterError) { fpars.setParameter("nuc_g_0", spDefault.g_0); }
   try { fpars.getParameter("nuc_g_1", spDefault.g_1); } catch (utl::Parameters::ParameterError) { fpars.setParameter("nuc_g_1", spDefault.g_1); }
   try { fpars.getParameter("nuc_g_2", spDefault.g_2); } catch (utl::Parameters::ParameterError) { fpars.setParameter("nuc_g_2", spDefault.g_2); }

   //Electron injection spectrum
   try { fpars.getParameter("electron_rigid_br0", spElectron.rigid_br0); } catch (utl::Parameters::ParameterError) { fpars.setParameter("electron_rigid_br0", spElectron.rigid_br0); }
   try { fpars.getParameter("electron_rigid_br", spElectron.rigid_br); } catch (utl::Parameters::ParameterError) { fpars.setParameter("electron_rigid_br", spElectron.rigid_br); }
   try { fpars.getParameter("electron_g_0", spElectron.g_0); } catch (utl::Parameters::ParameterError) { fpars.setParameter("electron_g_0", spElectron.g_0); }
   try { fpars.getParameter("electron_g_1", spElectron.g_1); } catch (utl::Parameters::ParameterError) { fpars.setParameter("electron_g_1", spElectron.g_1); }
   try { fpars.getParameter("electron_g_2", spElectron.g_2); } catch (utl::Parameters::ParameterError) { fpars.setParameter("electron_g_2", spElectron.g_2); }

   //Read in the isotopic abundances.  Only use those defined in the parameters object
   //max_Z is not accessible, only particles created will be added anyway.
   for ( int iZ = 1; iZ < 90; ++iZ ) {
      for ( int iA=iZ; iA<3*iZ; ++iA ) {
         try {
            double abundance;
            std::ostringstream parname;
            parname << "iso_abundance_" << std::setw(2) << std::setfill('0') << iZ << '_' << std::setw(3) << std::setfill('0') << iA;
            fpars.getParameter(parname.str(), abundance);
            isotopic_abundance[std::make_pair(iZ,iA)] = abundance;

            //Re-use the parname string stream
            parname.str("");
            parname<<"Abundance for ("<< iA <<", "<< iZ <<") is "<<abundance;
            INFO(parname.str());

         } catch (utl::Parameters::ParameterError) {}
      }
   }

   //Now set the spectral parameters of those if found.
   for (auto it = isotopic_abundance.begin(); it != isotopic_abundance.end(); ++it) {
      const int iZ = it->first.first;
      const int iA = it->first.second;

      //Only add it if needed
      bool found(false);

      //Use default parameters
      specProperties spec = spDefault;

      std::ostringstream parname;
      parname << "nuc_rigid_br0_" << std::setw(2) << std::setfill('0') << iZ << '_' << std::setw(3) << std::setfill('0') << iA;
      try { 
         fpars.getParameter(parname.str(), spec.rigid_br0); 
         found = true;
      } catch (utl::Parameters::ParameterError) {}

      parname.str("");
      parname << "nuc_rigid_br_" << std::setw(2) << std::setfill('0') << iZ << '_' << std::setw(3) << std::setfill('0') << iA;
      try { 
         fpars.getParameter(parname.str(), spec.rigid_br); 
         found = true;
      } catch (utl::Parameters::ParameterError) {}
      
      parname.str("");
      parname << "nuc_g_0_" << std::setw(2) << std::setfill('0') << iZ << '_' << std::setw(3) << std::setfill('0') << iA;
      try { 
         fpars.getParameter(parname.str(), spec.g_0); 
         found = true;
      } catch (utl::Parameters::ParameterError) {}
      
      parname.str("");
      parname << "nuc_g_1_" << std::setw(2) << std::setfill('0') << iZ << '_' << std::setw(3) << std::setfill('0') << iA;
      try { 
         fpars.getParameter(parname.str(), spec.g_1); 
         found = true;
      } catch (utl::Parameters::ParameterError) {}
      
      parname.str("");
      parname << "nuc_g_2_" << std::setw(2) << std::setfill('0') << iZ << '_' << std::setw(3) << std::setfill('0') << iA;
      try { 
         fpars.getParameter(parname.str(), spec.g_2); 
         found = true;
      } catch (utl::Parameters::ParameterError) {}

      if (found)
         iso_inj_spectra[it->first] = spec;
   }

   try { fpars.getParameter("source_specification", source_specification); } catch (utl::Parameters::ParameterError) { fpars.setParameter("source_specification", source_specification); }

   try { fpars.getParameter("source_normalization", source_normalization); } catch (utl::Parameters::ParameterError) { fpars.setParameter("source_normalization", source_normalization); }
   try { fpars.getParameter("electron_source_normalization", electron_source_normalization); } catch (utl::Parameters::ParameterError) { fpars.setParameter("electron_source_normalization", electron_source_normalization); }

   try { fpars.getParameter("proton_norm_Ekin", proton_norm_Ekin); } catch (utl::Parameters::ParameterError) { fpars.setParameter("proton_norm_Ekin", proton_norm_Ekin); }
   try { fpars.getParameter("spectra_norm_Rigidity", spectra_norm_Rigidity); } catch (utl::Parameters::ParameterError) {  }

   //CR source parameters
   try { fpars.getParameter("source_model", source_model); } catch (utl::Parameters::ParameterError) { fpars.setParameter("source_model", source_model); }

   //If not defined, electron source model is the same as nuclei
   source_model_electron = source_model;
   try { fpars.getParameter("source_model_electron", source_model_electron); } catch (utl::Parameters::ParameterError) {}

   //Resize the source_prameters vectors
   for (int iS = 0; iS < 10; ++iS) {

      std::ostringstream buf1, buf2;
      buf1 << "source_parameters_" << iS;
      buf2 << "source_pars_elec_" << iS;

      try { fpars.getParameter(buf1.str(), source_parameters[iS]); } catch (utl::Parameters::ParameterError) { fpars.setParameter(buf1.str(), source_parameters[iS]); }

      source_parameters_electron[iS] = source_parameters[iS]; 
      try { fpars.getParameter(buf2.str(), source_parameters_electron[iS]); } catch (utl::Parameters::ParameterError) {}

   }

   //Informative output
   std::ostringstream srcBuf;

   srcBuf << "Source model nuclei: " << source_model << ", electrons: " << source_model_electron;
   INFO(srcBuf.str());

   for (int iS = 0; iS < 10; ++iS) {

      srcBuf.str("");
      srcBuf << iS << " " << source_parameters[iS] << " " << source_parameters_electron[iS];
      INFO(srcBuf.str());

   }

   //Tabulated parameters, only needed if source_model == 8
   try { fpars.getParameter("source_values", source_values); } catch (utl::Parameters::ParameterError) {}
   try { fpars.getParameter("source_radius", source_radius); } catch (utl::Parameters::ParameterError) {}

   if (source_values.size() != source_radius.size()) {
      ERROR("source_values vector not the same size as source_radius vector");
      throw(std::runtime_error("Fix your configuration file."));
   }

   source_values_electron = source_values;
   source_radius_electron = source_radius;
   try { fpars.getParameter("source_values_electron", source_values_electron); } catch (utl::Parameters::ParameterError) {}
   try { fpars.getParameter("source_radius_electron", source_radius_electron); } catch (utl::Parameters::ParameterError) {}
   if (source_values_electron.size() != source_radius_electron.size()) {
      ERROR("source_values_electron vector not the same size as source_radius_electron vector");
      throw(std::runtime_error("Fix your configuration file."));
   }

   // The xml file for source model == 15
   try { fpars.getParameter("source_xmlFile", source_xmlFile); } catch (utl::Parameters::ParameterError) { fpars.setParameter("source_xmlFile", source_xmlFile); }

   //Individual cr sources for 3D
   try { fpars.getParameter("n_cr_sources", n_cr_sources); } catch (utl::Parameters::ParameterError) { fpars.setParameter("n_cr_sources", n_cr_sources); }
   cr_source_x.resize(n_cr_sources,0.0);
   cr_source_y.resize(n_cr_sources,0.0);
   cr_source_z.resize(n_cr_sources,0.0);
   cr_source_w.resize(n_cr_sources,0.0);
   cr_source_L.resize(n_cr_sources,0.0);

   //These are not optional
   for ( int i_cr_source=0; i_cr_source<n_cr_sources; i_cr_source++ ) {
      std::ostringstream buf;
      buf << "cr_source_x_" << std::setw(2) << std::setfill('0') << i_cr_source+1;
      pars.getParameter(buf.str(), cr_source_x[i_cr_source]);

      buf.str("");
      buf << "cr_source_y_" << std::setw(2) << std::setfill('0') << i_cr_source+1;
      pars.getParameter(buf.str(), cr_source_y[i_cr_source]);

      buf.str("");
      buf << "cr_source_z_" << std::setw(2) << std::setfill('0') << i_cr_source+1;
      pars.getParameter(buf.str(), cr_source_z[i_cr_source]);

      buf.str("");
      buf << "cr_source_w_" << std::setw(2) << std::setfill('0') << i_cr_source+1;
      pars.getParameter(buf.str(), cr_source_w[i_cr_source]);

      buf.str("");
      buf << "cr_source_L_" << std::setw(2) << std::setfill('0') << i_cr_source+1;
      pars.getParameter(buf.str(), cr_source_L[i_cr_source]);
   }

   //SNR events
   try { fpars.getParameter("SNR_events", SNR_events); } catch (utl::Parameters::ParameterError) { fpars.setParameter("SNR_events", SNR_events); }
   try { fpars.getParameter("SNR_interval", SNR_interval); } catch (utl::Parameters::ParameterError) {}
   try { fpars.getParameter("SNR_livetime", SNR_livetime); } catch (utl::Parameters::ParameterError) {}

   try { fpars.getParameter("SNR_electron_sdg", SNR_electron_sdg); } catch (utl::Parameters::ParameterError) {}
   try { fpars.getParameter("SNR_nuc_sdg", SNR_nuc_sdg); } catch (utl::Parameters::ParameterError) {}

   try { fpars.getParameter("SNR_electron_dgpivot", SNR_electron_dgpivot); } catch (utl::Parameters::ParameterError) {}
   try { fpars.getParameter("SNR_nuc_dgpivot", SNR_nuc_dgpivot); } catch (utl::Parameters::ParameterError) {}

   if (SNR_events == 1)
      createSNRDistributions = true;
}


void SourceClass_Compatibility::addSource( Particle &particle ) const 
{
   //Make sure the source needs adding
   if (isSecondary(particle))
      return;

   const auto type = std::make_pair(particle.Z, particle.A);

   const auto isoIt = isotopic_abundance.find(type);
   
   // Only add injection for found isotopes and primary electrons.  Ignore K-capture nuclei
   if ( ( isoIt == isotopic_abundance.end() && particle.name != "primary_electrons" ) || particle.K_electron > 0 )
      return;

   //Create the SNR distributions if needed.
   create_SNR(particle);

   std::vector<double> spec_shape = createSpecShape(particle);
  //This we cannot do
  //particle.primary_source_function = 0.; // IMOS20020418 whole 2D/3D particle.primary_source_function loops are changed
  
  //To allow for different source distribution of electrons and nuclei
  const std::vector<double>* parameters;
  const std::vector<double>* values;
  const std::vector<double>* radius;
  int model;
  
  const string priElecStr = "primary_electrons";

  if (priElecStr != particle.name) {
    
    model = source_model;
    parameters = &source_parameters;
    values = &source_values;
    radius = &source_radius;
    
  } else {
    
    model = source_model_electron;
    parameters = &source_parameters_electron;
    values = &source_values_electron;
    radius = &source_radius_electron;
    
  }

  if (2 == particle.n_spatial_dimensions) {


     int ir=0, iz=particle.n_zgrid/2;
     if(source_specification==1)
     {
        for (int ip = 0; ip < particle.n_pgrid; ++ip) {
           particle.primary_source_function.d2[ir][iz].s[ip]+=spec_shape[ip];
        }
     }
     for(ir=0; ir<particle.n_rgrid; ir++)
     {
        if(source_specification==2)
        {
           for (int ip = 0; ip < particle.n_pgrid; ++ip) {
              particle.primary_source_function.d2[ir][iz].s[ip]+=spec_shape[ip];
           }
        }
        for(iz=0; iz<particle.n_zgrid; iz++)
        {
           if (0 == source_specification) {
              const double sd = source_distribution(particle. r[ir], 0, particle.z[iz], model, particle.n_spatial_dimensions, *parameters, *values, *radius);
              for (int ip = 0; ip < particle.n_pgrid; ++ip) {
                 particle.primary_source_function.d2[ir][iz].s[ip] += sd *spec_shape[ip];


              }
           }
        }
     }
  }

  //if (priElecStr == particle.name)
  //exit(0);

  if (3 == particle.n_spatial_dimensions) {

     int ix = particle.n_xgrid/2;
     int iy = particle.n_ygrid/2;
     int iz = particle.n_zgrid/2;

     if (1 == source_specification)
        for (int ip = 0; ip < particle.n_pgrid; ++ip) {
           particle.primary_source_function.d3[ix][iy][iz].s[ip] += spec_shape[ip];
        }

     for(ix=0; ix<particle.n_xgrid; ix++) {

        for(iy=0; iy<particle.n_ygrid; iy++) {

           if (2 == source_specification)
              for (int ip = 0; ip < particle.n_pgrid; ++ip) {
                 particle.primary_source_function.d3[ix][iy][iz].s[ip] += spec_shape[ip];
              }

           for(iz=0; iz<particle.n_zgrid; iz++) {

              const double sd = source_distribution(particle.x[ix],
                    particle.y[iy],
                    particle.z[iz],
                    model,
                    particle.n_spatial_dimensions,
                    *parameters,
                    *values,
                    *radius);

              for (int ip = 0; ip < particle.n_pgrid; ++ip) {
                 // source spectral index dispersion

                 double spec_dg_ratio(1);                             //AWS20010411

                 if (1 == SNR_events) {

                    if ("primary_electrons" == particle.name) { 

                       spec_dg_ratio=
                          pow(particle.rigidity[ip]/SNR_electron_dgpivot,1.*SNR_electron_dg.d3[ix][iy][iz].s[0]);

                       /*if(galdef.verbose==-501) { // selectable debug

                         std::ostringstream ost1;
                         ost1<<"SNR_electron_dg="<<SNR_electron_dg.d3[ix][iy][iz].s[0]
                         <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio;
                         INFO(ost1.str());
                         }
                         */
                    } //electrons
                    else {

                       spec_dg_ratio=
                          pow(particle.rigidity[ip]/SNR_nuc_dgpivot,     1.*SNR_nuc_dg.     d3[ix][iy][iz].s[0]);

                       /*
                          if(verbose==-501) { // selectable debug

                          std::ostringstream ost1;
                          ost1<<"SNR_nuc_dg="<<SNR_nuc_dg.d3[ix][iy][iz].s[0] 
                          <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio;
                          INFO(ost1.str());
                          }
                          */
                    }//nucleons		     
                 }// if SNR_events

                 if (0 == source_specification)
                    particle.primary_source_function.d3[ix][iy][iz].s[ip] += sd*spec_dg_ratio*spec_shape[ip];

              } //ip

           } //iz

        } //iy

     } //ix

  } // 3D
  
  //Mark steady state for removal in stocastic SNR events
  removeSteadyState = true;
}

void SourceClass_Compatibility::addSource(Particle &particle, double t) const
{
   //Only in 3D
   if (particle.n_spatial_dimensions != 3)
      return;

   //Make sure the source needs adding
   const auto type = std::make_pair(particle.Z, particle.A);

   const auto isoIt = isotopic_abundance.find(type);
   
   // Only add injection for found isotopes and primary electrons.  Ignore K-capture nuclei
   if ( ( isoIt == isotopic_abundance.end() && particle.name != "primary_electrons" ) || particle.K_electron > 0 )
      return;

   source_SNR_event(particle, t);
}


std::vector<double> SourceClass_Compatibility::createSpecShape(const Particle &particle) const
{
   //Now set up the spectral parameters, code copied from fill_transport_arrays
   double g_0=spDefault.g_0, 
          rigid_br0=spDefault.rigid_br0,
          rigid_br=spDefault.rigid_br,                           // IMOS20000607
          g_1=spDefault.g_1,
          g_2=spDefault.g_2;

   // Check for specific injection spectra
   const auto type = std::make_pair(particle.Z, particle.A);
   const auto specIt = iso_inj_spectra.find(type);
   if (specIt != iso_inj_spectra.end()) {
      g_0 = specIt->second.g_0;
      g_1 = specIt->second.g_1;
      g_2 = specIt->second.g_2;
      rigid_br0 = specIt->second.rigid_br0;
      rigid_br = specIt->second.rigid_br;
   }


   const string priElecStr = "primary_electrons";

   if (priElecStr == particle.name) {
      g_0=spElectron.g_0;                                    // IMOS20031012
      rigid_br0=spElectron.rigid_br0;
      g_1=spElectron.g_1;
      rigid_br=spElectron.rigid_br;
      g_2=spElectron.g_2;
   }
   std::ostringstream ost;
   ost<<particle.name<<" g_0="<<g_0<<"  rigid_br0= "<<rigid_br0  // IMOS20031012
      <<" g_1="<<g_1<<"  rigid_br= " <<rigid_br <<" g_2="<<g_2;
   INFO(ost.str());

   //Normalization for the spectral shape should be evaluated at proton_norm_Ekin
   //Calculate for convenience the corresponding rigidity and momentum
   double pNorm(-1.), EtotNorm(1), betaNorm(0), gammaNorm(1), rigidityNorm(0);
   double EkinNorm(proton_norm_Ekin);
   kinematic(particle.Z, particle.A, particle.mass, pNorm, EkinNorm,
         EtotNorm, betaNorm, gammaNorm, rigidityNorm, 0);

   //Ignore the primary_abundance because that can differ between classes.
   const auto isoIt = isotopic_abundance.find(type);
   double specNorm = isoIt->second; //particle.primary_abundance;
   if (priElecStr == particle.name) 
      specNorm = 1;

   //Switch to spectra_norm_Rigidity for backwards compatibility if necessary
   if (spectra_norm_Rigidity > 0)
      rigidityNorm = spectra_norm_Rigidity;

   //Use the same stupid swithces as in the spectral shapes loop.  This should be streamlined.
   if (inj_spectrum_type == "Etot") {

      specNorm *= pow(EtotNorm, g_2);

   } else {

      if (rigidityNorm < rigid_br0) {
         specNorm *= pow(rigidityNorm/rigid_br0,g_0)*pow(rigid_br0/rigid_br,g_1);
      } else if (rigidityNorm < rigid_br) {
         specNorm *= pow(rigidityNorm/rigid_br,g_1);
      } else {
         specNorm *= pow(rigidityNorm/rigid_br,g_2);
      }

   }

   if (inj_spectrum_type == "beta_rig") 
      specNorm *= sqrt(1. + pow(2.e3/rigidityNorm, 2.)); 

   if (inj_spectrum_type=="dirac") 
      specNorm = isoIt->second; //particle.primary_abundance;

   if(inj_spectrum_type=="step")
      specNorm = isoIt->second; //particle.primary_abundance;

   if(inj_spectrum_type=="expcutoff")
      specNorm = isoIt->second * pow(EkinNorm,g_2)*exp(EkinNorm/rigid_br);

   if(inj_spectrum_type=="doubleexpcutoff")
      specNorm = isoIt->second * pow(EkinNorm,g_2)*exp(EkinNorm/rigid_br)*exp(rigid_br0/EkinNorm);


   if (particle.A > 1)
      specNorm /= particle.A;


   //Normalize the spectrum
   // CASE: PRIMARY NUCLEI                                                         AWS20000601.1

   if ("primary_electrons" != particle.name) //strcmp(particle.name,"primary_electrons")!=0)
   {
      specNorm *= source_normalization;
      ost.str("");
      ost<<" >>>>>>>>>>>>>>>>>>> norm "<<source_normalization<<" >>>>>>>>>>>>>>>>>>>";
      INFO(ost.str());
   }

   // CASE: PRIMARY ELECTRONS                                                      IMOS20031016

   if ("primary_electrons" == particle.name) //strcmp(particle.name,"primary_electrons")==0)
   {
      specNorm *= electron_source_normalization;
      ost.str("");
      ost<<" >>>>>>>>>> electron_norm "<<electron_source_normalization<<" >>>>>>>>>>>>>>>>>>>";
      INFO(ost.str());
   }


   //Spectral shape is (currently) independent of position.  Move it outside the loop.
   std::vector<double> spec_shape(particle.n_pgrid, 0.0);
   for (int ip = 0; ip < particle.n_pgrid; ++ip) {

      //TODO: This switch statement is rather arcane, should be fixed
      if (inj_spectrum_type == "Etot") {
         spec_shape[ip] = pow(particle.Etot[ip],-g_2); // IMOS20000613
         continue;
      }
      else {

         if(particle.rigidity[ip]< rigid_br0) // IMOS20031012
            spec_shape[ip] = pow(particle.rigidity[ip]/rigid_br0,-g_0)*pow(rigid_br0/rigid_br,-g_1);
         if (rigid_br0 <= particle.rigidity[ip] && particle.rigidity[ip]< rigid_br)
            spec_shape[ip] = pow(particle.rigidity[ip]/rigid_br, -g_1);

         if(rigid_br <= particle.rigidity[ip])
            spec_shape[ip] =pow(particle.rigidity[ip]/rigid_br, -g_2);

      }

      if (inj_spectrum_type == "beta_rig") 
         spec_shape[ip] /= sqrt(1. + pow(2.e3/particle.rigidity[ip], 2.)); // IMOS20011210

      //new spectral shapes  AWS20101202
      if (inj_spectrum_type=="dirac") {
         const double Edirac = rigid_br;

         if (log(Edirac) >= log(particle.Ekin[ip]) - 0.5*log(particle.Ekin_factor) && 
               log(Edirac) < log(particle.Ekin[ip]) + 0.5*log(particle.Ekin_factor)) {
            spec_shape[ip] = 1.;

         } else 
            spec_shape[ip] = 0.;

         continue;

      } 

      if(inj_spectrum_type=="step")
      {
         if (rigid_br >= particle.Ekin[ip]) spec_shape[ip]=1.0;
         else spec_shape[ip]=0.0;
         continue;
      }

      if(inj_spectrum_type=="expcutoff")
      {
         spec_shape[ip]=pow(particle.Ekin[ip],-g_2)*exp(-particle.Ekin[ip]/rigid_br);
         continue;
      }

      if(inj_spectrum_type=="doubleexpcutoff")
      {
         spec_shape[ip]=pow(particle.Ekin[ip],-g_2)*exp(-particle.Ekin[ip]/rigid_br)*exp(-rigid_br0/particle.Ekin[ip]);
         continue;
      }
      //                      AWS20101201
      spec_shape[ip] *= specNorm;
   }


   return spec_shape;
}
