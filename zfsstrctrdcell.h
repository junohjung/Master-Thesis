#ifndef ZFSSTRCTRDCELL
#define ZFSSTRCTRDCELL

//#include "zfslist.h"
#include "zfscellproperties.h"

//contains the structs to store cell information and surface information
class ZFSStrctrdCell {
 public:
  ZFSStrctrdCell() = default;
  ~ZFSStrctrdCell() = default;
  ZFSFloat** coordinates;
  ZFSFloat** mgOldCoordinates;
  ZFSFloat** variables;
  ZFSFloat** pvariables;
  ZFSFloat** rightHandSide;
  ZFSFloat**  dxt; //volume fluxes for the three cell surfaces
  ZFSFloat*  flux; //contains the convective flux over the surface
  ZFSFloat** area; //contains 3 of the areas of the volume
  ZFSFloat** oldVariables;
  ZFSFloat*  cellVolume;
  ZFSFloat* uVWT; //array of central derivative of primitive variables.
  ZFSFloat* eFlux;
  ZFSFloat* fFlux;
  ZFSFloat* gFlux;
  ZFSFloat* temperature;
  ZFSFloat* lamvisc;
  ZFSFloat* ds;

  ZFSFloat* cellJac;
  ZFSFloat* oldCellJac;
  ZFSFloat* cornerJac;
  ZFSFloat* viscousFlux;
  ZFSFloat** cellMetrics;
  ZFSFloat*  cornerMetrics;
  ZFSFloat** surfaceMetrics;
  ZFSFloat** surfaceDist;//contains the surface distance
  ZFSFloat** cellLength;
  ZFSFloat* localTimeStep;

  //fq field
  ZFSFloat** fq;
  ZFSFloat** stg_fq;

  //spongeLayer
  ZFSFloat* spongeFactor;

  //auxillary data
  ZFSFloat* cf;
  ZFSFloat* cp;
  ZFSFloat* powerVisc;
  ZFSFloat* powerPres;
  ZFSInt* cfOffsets;
  ZFSInt* cpOffsets;
  ZFSInt* powerOffsets;

  //least squares
  ZFSFloat** reconstructionConstants;
  ZFSInt* nghbr;
  ZFSInt* numOfNghbr;

  //RANS
  ZFSFloat* wallDistance;
  ZFSFloat* saFlux1;
  ZFSFloat* saFlux2;
  ZFSFloat* prodDest;

 private:
  //no private data yet
};



#endif
