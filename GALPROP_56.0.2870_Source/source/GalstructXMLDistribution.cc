#include "GalstructXMLDistribution.h"

GalstructXMLDistribution::GalstructXMLDistribution(const utl::Parameters &pars) :
   SpatialDistribution(pars)
{
   std::string filename;
   pars.getParameter("galstructXMLfile", filename);
   readXMLfile(filename);
}

void GalstructXMLDistribution::readXMLfile(const std::string &filename)
{
   profiles.clear();

   utl::Reader xmlreader(filename);
   // Add all cylindrical profiles in the top branch to the model
   for ( utl::Branch rp = xmlreader.GetTopBranch().GetFirstChild(); rp; rp = rp.GetNextSibling()) {
      //Make sure it is the correct type before adding
      if ( rp.GetBranchNameString() == "CylindricalProfile" ) {
         profiles.push_back( GalacticStructure::CylindricalProfile::createProfile(rp));
      }
   }

}

double GalstructXMLDistribution::operator() (double x, double y, double z) const
{
   //The profiles expect cylindrical coordinates
   const double theta = atan2(y, x);
   const double r = sqrt(x*x+y*y);

   double results(0);
   for (size_t i(0); i < profiles.size(); ++i)
      results += profiles[i]->operator()(r, theta, z);

   return results;
}
