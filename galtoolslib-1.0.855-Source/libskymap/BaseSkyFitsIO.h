#ifndef BASESKYFITSIO_H
#define BASESKYFITSIO_H

#include "BaseSky.h"

#include <CCfits/CCfits>

namespace SM {

/** \brief Write map to fits file
 *
 * The format should be readable by aladdin, i.e. one map per spectral bin.
 * * We do two different schemes, implicit and explicit. 
 * * * Implicit assumes full sky.
 * * * Explicit uses partial sky, but identical for all energies.
 * * The map is stored in the most appropriate unit for the type of the map.
 *   We only handle basic types, maps of special types cannot be stored with this method.
 * * If overwrite is true or file does not exist, we store the map in the first
 *   extension in a binary table.  If overwrite is false, we
 *   add a new extension to the file if it exists, otherwise we create it.
 * * If compress is true we add a ".gz" extension for automatic compression.
 * * Information on energy units and map units can be added.
 * * Any additional keywords can also be added.
 */
template <typename T>
void writeToFits(const BaseSky<T> &map,
      const std::string &filename,
      bool overwrite=true,
#ifdef ENABLE_COMPRESSION
      bool compress=true,
#else
      bool compress=false,
#endif
      const std::string &energyType = "energy",
      const std::string &energyUnit = "MeV",
      const std::string &mapUnit = "",
      const std::map<std::string, std::string> &additionalFitsKeywords = (std::map<std::string, std::string>())
      );

/** \brief Read map from fits file
 *
 * Handles HEALPix formats from previous versions transparently.
 * Does not handle other formats at all and will throw an exception.
 * The map_ptr is replaced with a new map containing correct dimensions.
 * If map_ptr is empty, we create a skymap corresponding to the type that was saved.
 * If map_ptr is not empty, we reset the map to correct dimensions.
 *
 * We try to gracefully handle other HEALPix formats, but may fail.
 * We only support single values per row in each column.
 *
 * keywords holds a list of all the keywords found in the HEALPix map header
 *
 * Reading general WCS maps is now handled by another function.
 */
template <typename T>
void readFromFits(std::unique_ptr< BaseSky<T> > &map_ptr,
      const std::string &filename,
      std::map<std::string, std::string> &keywords
      );


/** \brief  This routine reads a WCS map into a BaseSky.
 *
 * If solidAngleCorrection is false then map is assumed to be in units 
 * sr^-1 (intensity), otherwise we take into account the different solid
 * angle of the input and output.
 *
 * If map_ptr is not empty, any data in map is erased and the spectral information is read 
 * from the fits file.  Spatial information from the input map
 * is kept as is and coordinate transform is performed as needed.
 *
 * If map_ptr is empty, we create a skymap corresponding using
 * SparseSky if the map covers less than half the sky but FullSky
 * otherwise.  The order of the map is set to agree roughly with 
 * the spatial dimension of the input.  No coordinate transform
 * is done in this case.
 */
template <typename T>
void readFromFitsWcs( std::unique_ptr< BaseSky<T> > &map_ptr, const std::string &fileName, bool solidAngleCorrection);

/* TODO:
 *
 * Write a routine that converts healpix maps to WCS.
 */

}//Namespace

#endif
