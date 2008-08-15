#include "Geometry/HcalTowerAlgo/interface/HcalFlexiHardcodeGeometryLoader.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/IdealObliquePrism.h"
#include "Geometry/CaloGeometry/interface/IdealZPrism.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>

namespace {
  const int MAX_HCAL_PHI = 72;
  const double DEGREE2RAD = M_PI / 180.;

  // ----------> HB <-----------
  struct HBHOCellParameters {
    HBHOCellParameters (int f_eta, int f_depth, int f_phiFirst, int f_phiStep, int f_dPhi, float f_rMin, float f_rMax, float f_etaMin, float f_etaMax)
      : eta(f_eta), depth(f_depth), phiFirst(f_phiFirst), phiStep(f_phiStep), dphi(f_dPhi), rMin(f_rMin), rMax(f_rMax), etaMin(f_etaMin), etaMax(f_etaMax)
    {}
 
    int eta;
    int depth;
    int phiFirst;
    int phiStep;
    int dphi;
    float rMin;
    float rMax;
    float etaMin;
    float etaMax;
  };

  const float HBRMIN = 181.225;
  const float HBRMAX = 289.009;

  HBHOCellParameters hbCells [] = {
    HBHOCellParameters ( 1, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*0, 0.087*1),
    HBHOCellParameters ( 2, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*1, 0.087*2),
    HBHOCellParameters ( 3, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*2, 0.087*3),
    HBHOCellParameters ( 4, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*3, 0.087*4),
    HBHOCellParameters ( 5, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*4, 0.087*5),
    HBHOCellParameters ( 6, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*5, 0.087*6),
    HBHOCellParameters ( 7, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*6, 0.087*7),
    HBHOCellParameters ( 8, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*7, 0.087*8),
    HBHOCellParameters ( 9, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*8, 0.087*9),
    HBHOCellParameters (10, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*9, 0.087*10),
    HBHOCellParameters (11, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*10, 0.087*11),
    HBHOCellParameters (12, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*11, 0.087*12),
    HBHOCellParameters (13, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*12, 0.087*13),
    HBHOCellParameters (14, 1, 1, 1, 5, HBRMIN, HBRMAX, 0.087*13, 0.087*14),
    HBHOCellParameters (15, 1, 1, 1, 5, HBRMIN,258.635, 0.087*14, 0.087*15),
    HBHOCellParameters (15, 2, 1, 1, 5,258.635,289.009, 0.087*14, 0.087*15),
    HBHOCellParameters (16, 1, 1, 1, 5, HBRMIN,190.595, 0.087*15, 0.087*16),
    HBHOCellParameters (16, 2, 1, 1, 5,190.595,232.837, 0.087*15, 0.087*16)
  };

  // ----------> HO <-----------
  const float HORMIN = 390.342;
  const float HORMAX = 413.959;

  HBHOCellParameters hoCells [] = {
    HBHOCellParameters ( 1, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*0, 0.087*1),
    HBHOCellParameters ( 2, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*1, 0.087*2),
    HBHOCellParameters ( 3, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*2, 0.087*3),
    HBHOCellParameters ( 4, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*3, 0.087*4),
    HBHOCellParameters ( 5, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*4, 0.087*5),
    HBHOCellParameters ( 6, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*5, 0.087*6),
    HBHOCellParameters ( 7, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*6, 0.087*7),
    HBHOCellParameters ( 8, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*7, 0.087*8),
    HBHOCellParameters ( 9, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*8, 0.087*9),
    HBHOCellParameters (10, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*9, 0.087*10),
    HBHOCellParameters (11, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*10, 0.087*11),
    HBHOCellParameters (12, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*11, 0.087*12),
    HBHOCellParameters (13, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*12, 0.087*13),
    HBHOCellParameters (14, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*13, 0.087*14),
    HBHOCellParameters (15, 4, 1, 1, 5, HORMIN, HORMAX, 0.087*14, 0.087*15)
  };

  // ----------> HE <-----------
  struct HECellParameters {
    HECellParameters (int f_eta, int f_depth, int f_phiFirst, int f_phiStep, int f_dPhi, float f_zMin, float f_zMax, float f_etaMin, float f_etaMax)
      : eta(f_eta), depth(f_depth), phiFirst(f_phiFirst), phiStep(f_phiStep), dphi(f_dPhi), zMin(f_zMin), zMax(f_zMax), etaMin(f_etaMin), etaMax(f_etaMax)
    {}
 
    int eta;
    int depth;
    int phiFirst;
    int phiStep;
    int dphi;
    float zMin;
    float zMax;
    float etaMin;
    float etaMax;
  };

  const float HEZMIN = 400.458;
  const float HEZMID = 436.168;
  const float HEZMAX = 549.268;

  HECellParameters heCells [] = {
    HECellParameters ( 16, 3, 1, 1, 5,391.883,427.468, 0.087*15, 0.087*16),
    HECellParameters ( 17, 1, 1, 1, 5, HEZMIN, HEZMAX, 0.087*16, 0.087*17),
    HECellParameters ( 18, 1, 1, 1, 5, HEZMIN,427.468, 0.087*17, 0.087*18),
    HECellParameters ( 18, 2, 1, 1, 5,427.468, HEZMAX, 0.087*17, 0.087*18),
    HECellParameters ( 19, 1, 1, 1, 5, HEZMIN, HEZMID, 0.087*18, 0.087*19),
    HECellParameters ( 19, 2, 1, 1, 5, HEZMID, HEZMAX, 0.087*18, 0.087*19),
    HECellParameters ( 20, 1, 1, 1, 5, HEZMIN, HEZMID, 0.087*19, 1.74),
    HECellParameters ( 20, 2, 1, 1, 5, HEZMID, HEZMAX, 0.087*19, 1.74),
    HECellParameters ( 21, 1, 1, 2,10, HEZMIN, HEZMID, 1.74, 1.83),
    HECellParameters ( 21, 2, 1, 2,10, HEZMID, HEZMAX, 1.74, 1.83),
    HECellParameters ( 22, 1, 1, 2,10, HEZMIN, HEZMID, 1.83, 1.93),
    HECellParameters ( 22, 2, 1, 2,10, HEZMID, HEZMAX, 1.83, 1.93),
    HECellParameters ( 23, 1, 1, 2,10, HEZMIN, HEZMID, 1.93, 2.043),
    HECellParameters ( 23, 2, 1, 2,10, HEZMID, HEZMAX, 1.93, 2.043),
    HECellParameters ( 24, 1, 1, 2,10, HEZMIN, HEZMID, 2.043, 2.172),
    HECellParameters ( 24, 2, 1, 2,10, HEZMID, HEZMAX, 2.043, 2.172),
    HECellParameters ( 25, 1, 1, 2,10, HEZMIN, HEZMID, 2.172, 2.322),
    HECellParameters ( 25, 2, 1, 2,10, HEZMID, HEZMAX, 2.172, 2.322),
    HECellParameters ( 26, 1, 1, 2,10, HEZMIN, HEZMID, 2.322, 2.500),
    HECellParameters ( 26, 2, 1, 2,10, HEZMID, HEZMAX, 2.322, 2.500),
    HECellParameters ( 27, 1, 1, 2,10, HEZMIN,418.768, 2.500, 2.650),
    HECellParameters ( 27, 2, 1, 2,10,418.768, HEZMID, 2.500, 2.650),
    HECellParameters ( 27, 3, 1, 2,10, HEZMID, HEZMAX, 2.500, 2.650),
    HECellParameters ( 28, 1, 1, 2,10, HEZMIN,418.768, 2.650, 2.868),
    HECellParameters ( 28, 2, 1, 2,10,418.768, HEZMID, 2.650, 2.868),
    HECellParameters ( 28, 3, 1, 2,10, HEZMID, HEZMAX, 2.650, 3.000),
    HECellParameters ( 29, 1, 1, 2,10, HEZMIN, HEZMID, 2.868, 3.000),
    HECellParameters ( 29, 2, 1, 2,10,418.768, HEZMID, 2.868, 3.000)
  };

  // ----------> HF <-----------
  struct HFCellParameters {
    HFCellParameters (int f_eta, int f_depth, int f_phiFirst, int f_phiStep, int f_dPhi, float f_zMin, float f_zMax, float f_rMin, float f_rMax)
      : eta(f_eta), depth(f_depth), phiFirst(f_phiFirst), phiStep(f_phiStep), dphi(f_dPhi), zMin(f_zMin), zMax(f_zMax), rMin(f_rMin), rMax(f_rMax)
    {}
 
    int eta;
    int depth;
    int phiFirst;
    int phiStep;
    int dphi;
    float zMin;
    float zMax;
    float rMin;
    float rMax;
  };

  const float HFZMIN1 = 1115.;
  const float HFZMIN2 = 1137.;
  const float HFZMAX = 1280.1;

  HFCellParameters hfCells [] = {
    HFCellParameters (29, 1, 1, 2, 10, HFZMIN1, HFZMAX,116.2,130.0),
    HFCellParameters (29, 2, 1, 2, 10, HFZMIN2, HFZMAX,116.2,130.0),
    HFCellParameters (30, 1, 1, 2, 10, HFZMIN1, HFZMAX, 97.5,116.2),
    HFCellParameters (30, 2, 1, 2, 10, HFZMIN2, HFZMAX, 97.5,116.2),
    HFCellParameters (31, 1, 1, 2, 10, HFZMIN1, HFZMAX, 81.8, 97.5),
    HFCellParameters (31, 2, 1, 2, 10, HFZMIN2, HFZMAX, 81.8, 97.5),
    HFCellParameters (32, 1, 1, 2, 10, HFZMIN1, HFZMAX, 68.6, 81.8),
    HFCellParameters (32, 2, 1, 2, 10, HFZMIN2, HFZMAX, 68.6, 81.8),
    HFCellParameters (33, 1, 1, 2, 10, HFZMIN1, HFZMAX, 57.6, 68.6),
    HFCellParameters (33, 2, 1, 2, 10, HFZMIN2, HFZMAX, 57.6, 68.6),
    HFCellParameters (34, 1, 1, 2, 10, HFZMIN1, HFZMAX, 48.3, 57.6),
    HFCellParameters (34, 2, 1, 2, 10, HFZMIN2, HFZMAX, 48.3, 57.6),
    HFCellParameters (35, 1, 1, 2, 10, HFZMIN1, HFZMAX, 40.6, 48.3),
    HFCellParameters (35, 2, 1, 2, 10, HFZMIN2, HFZMAX, 40.6, 48.3),
    HFCellParameters (36, 1, 1, 2, 10, HFZMIN1, HFZMAX, 34.0, 40.6),
    HFCellParameters (36, 2, 1, 2, 10, HFZMIN2, HFZMAX, 34.0, 40.6),
    HFCellParameters (37, 1, 1, 2, 10, HFZMIN1, HFZMAX, 28.6, 34.0),
    HFCellParameters (37, 2, 1, 2, 10, HFZMIN2, HFZMAX, 28.6, 34.0),
    HFCellParameters (38, 1, 1, 2, 10, HFZMIN1, HFZMAX, 24.0, 28.6),
    HFCellParameters (38, 2, 1, 2, 10, HFZMIN2, HFZMAX, 24.0, 28.6),
    HFCellParameters (39, 1, 1, 2, 10, HFZMIN1, HFZMAX, 20.1, 24.0),
    HFCellParameters (39, 2, 1, 2, 10, HFZMIN2, HFZMAX, 20.1, 24.0),
    HFCellParameters (40, 1, 3, 4, 20, HFZMIN1, HFZMAX, 16.9, 20.1),
    HFCellParameters (40, 2, 3, 4, 20, HFZMIN2, HFZMAX, 16.9, 20.1),
    HFCellParameters (41, 1, 3, 4, 20, HFZMIN1, HFZMAX, 12.5, 16.9),
    HFCellParameters (41, 2, 3, 4, 20, HFZMIN2, HFZMAX, 12.5, 16.9)
  };

  //
  // Convert constants to appropriate cells
  //
  void fillHBHO (CaloSubdetectorGeometry* fGeometry, const HBHOCellParameters* fCells, int nCells , bool fHB) {
    for (int iCell = 0; iCell < nCells; ++iCell) {
      const HBHOCellParameters& param = fCells[iCell];
      for (int iPhi = param.phiFirst; iPhi <= MAX_HCAL_PHI; iPhi += param.phiStep) {
	for (int iside = -1; iside <= 1; iside += 2) { // both detector sides are identical
	  HcalDetId hid (fHB ? HcalBarrel : HcalOuter, param.eta*iside, iPhi, param.depth);
	  double phiCenter = ((iPhi-1)*360./MAX_HCAL_PHI + 0.5*param.dphi) * DEGREE2RAD; // middle of the cell
	  double etaCenter = 0.5*(param.etaMin + param.etaMax);
	  double x = param.rMin* cos (phiCenter);
	  double y = param.rMin* sin (phiCenter);
	  double z = iside * param.rMin * sinh(etaCenter);
	  // make cell geometry
	  GlobalPoint refPoint (x,y,z); // center of the cell's face
	  std::vector<float> cellParams; cellParams.reserve (3);
	  cellParams.push_back (0.5 * (param.etaMax - param.etaMin)); // deta_half
	  cellParams.push_back (0.5 * param.dphi * DEGREE2RAD);  // dphi_half
	  cellParams.push_back (0.5 * (param.rMax - param.rMin) * cosh (etaCenter)); // dr_half
	  
// 	  std::cout << "HcalFlexiHardcodeGeometryLoader::fillHBHO-> " << hid << hid.ieta() << '/' << hid.iphi() << '/' << hid.depth()
// 		    << refPoint << '/' << cellParams [0] << '/' << cellParams [1] << '/' << cellParams [2] << std::endl;
	  
	  CaloCellGeometry* newcell = 
	    new calogeom::IdealObliquePrism (refPoint, fGeometry->cornersMgr(),
					     CaloCellGeometry::getParmPtr (cellParams, fGeometry->parMgr(), fGeometry->parVecVec()));
	  // ... and store it
	  fGeometry->addCell (hid, newcell);						       
	}
      }
    }
  }
  
  void fillHE (CaloSubdetectorGeometry* fGeometry, const HECellParameters* fCells, int nCells ) {
    for (int iCell = 0; iCell < nCells; ++iCell) {
      const HECellParameters& param = fCells[iCell];
      for (int iPhi = param.phiFirst; iPhi <= MAX_HCAL_PHI; iPhi += param.phiStep) {
	for (int iside = -1; iside <= 1; iside += 2) { // both detector sides are identical
	  HcalDetId hid (HcalEndcap, param.eta*iside, iPhi, param.depth);
	  double phiCenter = ((iPhi-1)*360./MAX_HCAL_PHI + 0.5*param.dphi) * DEGREE2RAD; // middle of the cell
	  double etaCenter = 0.5 * (param.etaMin + param.etaMax);

	  double z = iside * param.zMin;
	  double perp = z / sinh (etaCenter);
	  double x = perp * cos (phiCenter);
	  double y = perp * sin (phiCenter);
	  // make cell geometry
	  GlobalPoint refPoint (x,y,z); // center of the cell's face
	  std::vector<float> cellParams; cellParams.reserve (3);
	  cellParams.push_back (0.5 * (param.etaMax - param.etaMin)); // deta_half
	  cellParams.push_back (0.5 * param.dphi * DEGREE2RAD);  // dphi_half
	  cellParams.push_back (-0.5 * (param.zMax - param.zMin) / tanh (etaCenter)); // dz_half, "-" means edges in Z
	  
// 	  std::cout << "HcalFlexiHardcodeGeometryLoader::fillHE-> " << hid << refPoint << '/' << cellParams [0] << '/' << cellParams [1] << '/' << cellParams [2] << std::endl;
	  
	  CaloCellGeometry* newcell = 
	    new calogeom::IdealObliquePrism (refPoint, fGeometry->cornersMgr(),
				       CaloCellGeometry::getParmPtr (cellParams, fGeometry->parMgr(), fGeometry->parVecVec()));
	  // ... and store it
	  fGeometry->addCell (hid, newcell);						       
	}
      }
    }
  }

  void fillHF (CaloSubdetectorGeometry* fGeometry, const HFCellParameters* fCells, int nCells ) {
    for (int iCell = 0; iCell < nCells; ++iCell) {
      const HFCellParameters& param = fCells[iCell];
      for (int iPhi = param.phiFirst; iPhi <= MAX_HCAL_PHI; iPhi += param.phiStep) {
	for (int iside = -1; iside <= 1; iside += 2) { // both detector sides are identical
	  HcalDetId hid (HcalForward, param.eta*iside, iPhi, param.depth);
	  double phiCenter = ((iPhi-1)*360./MAX_HCAL_PHI + 0.5*param.dphi) * DEGREE2RAD; // middle of the cell
	  GlobalPoint inner (param.rMin, 0, param.zMin);
	  GlobalPoint outer (param.rMax, 0, param.zMin);
	  double etaCenter = 0.5 * (inner.eta() + outer.eta());

	  double z = iside * param.zMin;
	  double perp = z / sinh (etaCenter);
	  double x = perp * cos (phiCenter);
	  double y = perp * sin (phiCenter);
	  // make cell geometry
	  GlobalPoint refPoint (x,y,z); // center of the cell's face
	  std::vector<float> cellParams; cellParams.reserve (3);
	  cellParams.push_back (0.5 * (inner.eta() - outer.eta())); // deta_half
	  cellParams.push_back (0.5 * param.dphi * DEGREE2RAD);  // dphi_half
	  cellParams.push_back (0.5 * (param.zMax - param.zMin)); // dz_half
	  
// 	  std::cout << "HcalFlexiHardcodeGeometryLoader::fillHF-> " << hid << refPoint << '/' << cellParams [0] << '/' << cellParams [1] << '/' << cellParams [2] << std::endl;
	  
	  CaloCellGeometry* newcell = 
	    new calogeom::IdealZPrism (refPoint, fGeometry->cornersMgr(),
				       CaloCellGeometry::getParmPtr (cellParams, fGeometry->parMgr(), fGeometry->parVecVec()));
	  // ... and store it
	  fGeometry->addCell (hid, newcell);						       
	}
      }
    }
  }

}

HcalFlexiHardcodeGeometryLoader::HcalFlexiHardcodeGeometryLoader()
{
}

CaloSubdetectorGeometry* HcalFlexiHardcodeGeometryLoader::load(const HcalTopology* fTopology) {
  CaloSubdetectorGeometry* hcalGeometry = new HcalGeometry (fTopology);
  if (!hcalGeometry->cornersMgr()) hcalGeometry->allocateCorners (9072) ;
  if (!hcalGeometry->parMgr()) hcalGeometry->allocatePar (500, 3);
  fillHBHO (hcalGeometry, hbCells, sizeof(hbCells)/sizeof(HBHOCellParameters), true);
  fillHBHO (hcalGeometry, hoCells, sizeof(hoCells)/sizeof(HBHOCellParameters), false);
  fillHE (hcalGeometry, heCells, sizeof(heCells)/sizeof(HECellParameters));
  fillHF (hcalGeometry, hfCells, sizeof(hfCells)/sizeof(HFCellParameters));
  return hcalGeometry;
}



