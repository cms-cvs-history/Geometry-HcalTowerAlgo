#ifndef GEOMETRY_HCALTOWERALGO_HCALFLEXIHARDCODEGEOMETRYLOADER_H
#define GEOMETRY_HCALTOWERALGO_HCALFLEXIHARDCODEGEOMETRYLOADER_H 1

/** \class HcalFlexiHardcodeGeometryLoader
 *
 * $Date: 2008/04/21 22:19:36 $
 * $Revision: 1.6 $
 * \author F.Ratnikov, UMd
*/

#include "Geometry/CaloGeometry/interface/CaloVGeometryLoader.h"

class HcalTopology;

class HcalFlexiHardcodeGeometryLoader 
{
   public:

      HcalFlexiHardcodeGeometryLoader();
  
      CaloSubdetectorGeometry* load(const HcalTopology* fTopology);
  
};

#endif
