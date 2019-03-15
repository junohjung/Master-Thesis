#ifndef ZFSSTRCTRDINTERPOLATION_H
#define ZFSSTRCTRDINTERPOLATION_H

#include "zfsglobals.h"
#include "zfspointbox.h"
#include "zfskdtree.h"
#include <vector>

using namespace std;

template <ZFSInt nDim>
class ZFSStrctrdInterpolation
{
 public:
  ZFSStrctrdInterpolation(const MPI_Comm strctrdCommunicator);
  ZFSStrctrdInterpolation(ZFSInt* noDonorCells, ZFSFloat** donorCoordinates, ZFSFloat** donorVariables, const MPI_Comm strctrdCommunicator);
  ZFSStrctrdInterpolation(ZFSInt* noDonorCellsDir, ZFSFloat** donorCoordinates, const MPI_Comm strctrdCommunicator); //junoh
  ~ZFSStrctrdInterpolation();
  void interpolateAtPoint(ZFSFloat* intPoint);
  void prepareInterpolationField(ZFSInt* noReceiverCells, ZFSFloat** receiverCoordinates) ;
  void interpolateField(ZFSString, ZFSFloat*);
  void loadDonorGrid();
  void loadDonorVariable(ZFSString varName);
  void prepareInterpolation(ZFSInt, ZFSFloat**, ZFSInt*);
  void interpolateVariables(ZFSFloat**);
  ZFSFloat getInterpolatedVariable(ZFSInt, ZFSInt);
  void prepareZonalInterpolation(ZFSInt, ZFSFloat**, ZFSInt*, ZFSBool);  //junoh */
  ZFSFloat interpolateVariableZonal(ZFSFloat*, ZFSInt);//junoh
  /* void interpolateVariablesNoFallback(ZFSFloat*);  //junoh */
  ZFSFloat getInterpolatedVariableZonal(ZFSFloat*, ZFSInt); //junoh

 protected:
  void buildDonorTree();
  inline void crossProduct( ZFSFloat result[3], ZFSFloat vec1[3], ZFSFloat vec2[3]);
  inline ZFSFloat scalarProduct( ZFSFloat vec1[3], ZFSFloat vec2[3]);
  inline ZFSId getCellIdfromCell( ZFSId origin, ZFSInt incI, ZFSInt incJ, ZFSInt incK, ZFSInt blockId);
  inline ZFSBool approx( const ZFSFloat&, const ZFSFloat&, const ZFSFloat); //junoh
  inline ZFSId getBlockId(ZFSInt cellId);
  inline void trilinearInterpolation(ZFSFloat*, ZFSInt, ZFSInt, ZFSFloat*, ZFSId);
  inline void trilinearInterpolation(ZFSFloat*, ZFSInt, ZFSFloat*, ZFSId);

  inline void computeInterpolationCoefficients(ZFSFloat*, ZFSInt);
  inline void transformPoint(ZFSId hexOrigin, ZFSFloat intPoint[3], ZFSFloat transformedPoint[3]);
  inline ZFSId findSurroundingHexahedron(ZFSFloat intPoint[3], ZFSId centerCellId, ZFSId stencil);
  inline void nearestNeighbourInterpolation(ZFSId, ZFSInt, ZFSFloat*) ;
  inline void nearestNeighbourInterpolation(ZFSId, ZFSFloat*) ;
  inline ZFSId cellIndex(ZFSInt i, ZFSInt j, ZFSInt k, ZFSInt blockId);
  inline ZFSId ic(ZFSInt, ZFSInt, ZFSInt);
  void computeCellCentreCoordinates(ZFSInt*, ZFSFloatScratchSpace&, ZFSInt, ZFSInt);

  //index variables
  static const ZFSId xsd = 0;
  static const ZFSId ysd = 1;
  static const ZFSId zsd = 2;

  ZFSInt** m_noDonorCellsDir;
  ZFSInt** m_noDonorPointsDir;
  ZFSInt* m_noDonorCells;
  ZFSInt* m_noDonorPoints;
  ZFSFloat** m_donorCoordinates;
  ZFSFloat** m_donorVariables;
  ZFSFloat* m_donorVar;
  ZFSId m_noDonorDims;
  ZFSId m_totalNoDonorCells;
  ZFSId* m_donorBlockOffsets;
  ZFSInt m_noBlocks;

  const MPI_Comm m_zfsStrctrdComm;
  ZFSInt m_domainId;
  ZFSInt m_noDonorVariables;
  ZFSBool m_donorIsCellCentered;
  ZFSBool m_isFieldInterpolation;
  ZFSFloat m_eps;

  ZFSId* m_donorOriginId;                 //needed if donorIds are saved for multiple interpolation
  ZFSFloat** m_interpolationCoefficients; //needed if coefficients are saved for multiple interpolation
  ZFSFloat* m_donorDistance; //junoh
  ZFSFloat* m_globalDonorDistanceMin; //junoh
  ZFSBool m_hasInterpolationPartnerDomain;//junoh
  ZFSId* m_hasInterpolationPartnersZonal;  //junoh
  ZFSId* m_hasInterpolationPartnersZonalGlobal; //junoh

  ZFSFloat** m_transformedReceiverPoints;
  ZFSFloat** m_receiverVariables;
  ZFSBool* m_hasInterpolationPartners;
  ZFSId m_noReceiverCells;
  ZFSId m_currentReceiverId;

  vector<Point<3>> m_donorPoints;
  KDtree<3>* m_donorTree;

	//pyramidPoints contains the ijk-combinations
//for all possible tetraeders
//inside a hexahedron, use together with function ic(tetraeder,side,dim)
	static constexpr ZFSInt m_pyramidPoints[72] = {0,1,0, 1,0,0, 0,0,0, 0,0,1,
																								 1,0,1, 0,1,1, 1,1,1, 1,1,0,
																								 1,0,1, 0,1,0, 1,0,0, 0,0,1,
																								 1,0,1, 0,1,0, 0,1,1, 0,0,1,
																								 1,0,1, 0,1,0, 0,1,1, 1,1,0,
																								 1,0,1, 0,1,0, 1,0,0, 1,1,0};
};

#endif
