#ifndef GALSTRUCT_XML_DISTRIBUTION_H
#define GALSTRUCT_XML_DISTRIBUTION_H

#include "SpatialDistribution.h"

#include <cylindricalprofiles.h>
#include <vector>
#include <memory>

class GalstructXMLDistribution : public SpatialDistribution {

   public:
      /** \brief Read in a xml file compatible with galstruct cylindrical distributions
       *
       * Only single parameters
       * * galstructXMLfile: path to the xml file
       */
      GalstructXMLDistribution(const utl::Parameters &pars);

      /** \brief Parse an xml file for cylindrical profiles
       *
       * Assumes there is a single global branch that contains many
       * cylindrical profiles that are each added to make a final model.
       */
      void readXMLfile(const std::string &filename);

      virtual double operator() (double x, double y, double z) const ;

   protected:
      std::vector<std::unique_ptr<GalacticStructure::CylindricalProfile> > profiles;
};
#endif
