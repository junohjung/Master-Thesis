#include "zfsstrctrdblck.h"
#include "zfsstrctrdinterpolation.h"
#include "zfsiolib.h"

template <ZFSInt nDim>
constexpr ZFSInt ZFSStrctrdInterpolation<nDim>::m_pyramidPoints[72];

/** \brief Constructor for property given donor
 *
 * The constructor for when the grid and the Q file
 * are specified in the property file.
 *
 * \author Marian Albers, Nov 2015
 */
template <ZFSInt nDim>
ZFSStrctrdInterpolation<nDim>::ZFSStrctrdInterpolation(const MPI_Comm strctrdCommunicator)
  :m_zfsStrctrdComm(strctrdCommunicator),
   m_donorIsCellCentered(true),
   m_isFieldInterpolation(false),
   m_eps(std::numeric_limits<ZFSFloat>::epsilon())
  {
    MPI_Comm_rank(m_zfsStrctrdComm, &m_domainId);

    loadDonorGrid();
    buildDonorTree();
  }

/** \brief Constructor for manually given donor
 *
 * The constructor when the grid and the variable field
 * are handed over manually, donorCoordinates are the coordinates
 * at which the donorVariables are located, both have the amount of cells
 * given by noDonorCells
 *
 * \author Marian Albers, Nov 2015
 */
template <ZFSInt nDim>
ZFSStrctrdInterpolation<nDim>::ZFSStrctrdInterpolation(ZFSInt* noDonorCellsDir,
                                                       ZFSFloat** donorCoordinates,
                                                       ZFSFloat** donorVariables,
                                                       const MPI_Comm strctrdCommunicator)
  :m_donorCoordinates(donorCoordinates),
   m_donorVariables(donorVariables),
   m_noBlocks(1),
   m_zfsStrctrdComm(strctrdCommunicator),
   m_noDonorVariables(5),
   m_donorIsCellCentered(true),
   m_isFieldInterpolation(false),
   m_eps(std::numeric_limits<ZFSFloat>::epsilon())
  {
    MPI_Comm_rank(m_zfsStrctrdComm, &m_domainId);
    zfsAlloc(m_noDonorCellsDir, m_noBlocks, nDim, __CALLING_FUNCTION__, 0, "m_noDonorCellsDir");
    zfsAlloc(m_donorBlockOffsets, m_noBlocks, __CALLING_FUNCTION__, 0, "m_donorBlockOffset");

    m_totalNoDonorCells = 1;
    for(ZFSId dim=0; dim<nDim; dim++) {
      m_noDonorCellsDir[0][dim] = noDonorCellsDir[dim];
      m_totalNoDonorCells *= noDonorCellsDir[dim];
    }

    buildDonorTree();
  }

template <ZFSInt nDim> //junoh
ZFSStrctrdInterpolation<nDim>::ZFSStrctrdInterpolation(ZFSInt* noDonorCellsDir,
                                                       ZFSFloat** donorCoordinates,
                                                       const MPI_Comm strctrdCommunicator)
  :m_donorCoordinates(donorCoordinates),
   m_noBlocks(1),
   m_zfsStrctrdComm(strctrdCommunicator),
   m_noDonorVariables(5),
   m_donorIsCellCentered(true),
   m_isFieldInterpolation(false),
   m_eps(std::numeric_limits<ZFSFloat>::epsilon())
  {
    MPI_Comm_rank(m_zfsStrctrdComm, &m_domainId);
    zfsAlloc(m_noDonorCellsDir, m_noBlocks, nDim, __CALLING_FUNCTION__, 0, "m_noDonorCellsDir");
    zfsAlloc(m_donorBlockOffsets, m_noBlocks, __CALLING_FUNCTION__, 0, "m_donorBlockOffset");

    m_totalNoDonorCells = 1;
    for(ZFSId dim=0; dim<nDim; dim++) {
      m_noDonorCellsDir[0][dim] = noDonorCellsDir[dim];
      m_totalNoDonorCells *= noDonorCellsDir[dim];
    }

    buildDonorTree();
  }





template <ZFSInt nDim>
ZFSStrctrdInterpolation<nDim>::~ZFSStrctrdInterpolation()
{
  if(m_isFieldInterpolation) {
    zfsDeallocate(m_donorOriginId);
    zfsDeallocate(m_transformedReceiverPoints);
    zfsDeallocate(m_hasInterpolationPartners);
    zfsDeallocate(m_donorVar);
    zfsDeallocate(m_hasInterpolationPartnersZonal);  //junoh
    zfsDeallocate(m_hasInterpolationPartnersZonalGlobal); //junoh
    zfsDeallocate(m_globalDonorDistanceMin);  //junoh
  } else {
    zfsDeallocate(m_donorOriginId);
    zfsDeallocate(m_interpolationCoefficients);
    zfsDeallocate(m_hasInterpolationPartners);
    zfsDeallocate(m_hasInterpolationPartnersZonal);   //junoh
    zfsDeallocate(m_hasInterpolationPartnersZonalGlobal); //junoh
    zfsDeallocate(m_globalDonorDistanceMin); //junoh
  }

  delete m_donorTree;
}


/** \brief Builds a kd-tree
* Creates a kd-tree from the predefined grid-data in
* m_donorCoordinates
*
* \author Marian Albers, Nov 2015
*/
template <ZFSInt nDim>
void ZFSStrctrdInterpolation<nDim>::buildDonorTree() {
  zfs_log << "Building up kd-tree..." << endl;

  // cout << "totalNoCells: " << m_totalNoDonorCells << endl;

  //first create kd tree from whole grid
  for(ZFSInt globalId=0; globalId<m_totalNoDonorCells; globalId++) {
    Point<3> a(m_donorCoordinates[0][globalId],
               m_donorCoordinates[1][globalId],
               m_donorCoordinates[2][globalId], globalId);
    m_donorPoints.push_back(a);
  }
  
  zfs_log << "Created points for kd-tree" << endl;

  //build up the tree and fill it
  m_donorTree = new KDtree<3> (m_donorPoints);

  zfs_log << "Building up kd-tree... FINISHED!" << endl;
}

/** \brief interpolates variables at point
* For a given 3D coordinate the method interpolates from
* the donor grid
*
* \author Marian Albers, Nov 2015
*/
template <ZFSInt nDim>
void ZFSStrctrdInterpolation<nDim>::interpolateAtPoint(ZFSFloat intPoint[3])
{
  //now go through own, fine cells and look for closest partner neighbour
  ZFSFloat distance = 0;
  Point<3> pt(intPoint[0], intPoint[1], intPoint[2]);
  cout << "Finding nearest point..." << endl;
  //find point on the grid that is closest to intPoint
  ZFSId centerCellId = m_donorTree->nearest(pt, distance);
  cout << "Finding nearest point... FINISHED! Point"
       << " x: " << m_donorCoordinates[0][centerCellId]
       << " y: " <<  m_donorCoordinates[1][centerCellId]
       << " z: " <<  m_donorCoordinates[2][centerCellId] << endl;
  cout << "Finding surrounding hexahedron..." << endl;
  //now eight hexahedron could be candidates for a new home for intPoint,
  //all around centerCellId
  ZFSId hexahedronOriginId = findSurroundingHexahedron(intPoint, centerCellId, 1);
  ZFSFloat interpolatedVariables[5];

  if(hexahedronOriginId != -1) {
    ZFSFloat transformedPoint[3];
    //now hexahedronOriginId is the id of the hexahedron in which intPoint is immersed
    transformPoint(hexahedronOriginId, intPoint, transformedPoint);
    ZFSId currentBlockId = getBlockId(hexahedronOriginId);
    //interpolate variables at transformed coordinate
    trilinearInterpolation(transformedPoint, hexahedronOriginId, interpolatedVariables, currentBlockId);
  } else {
    //fallback to nearest neighbour interpolation
    nearestNeighbourInterpolation(centerCellId, interpolatedVariables);
  }
}

/** \brief interpolates a field
 * interpolates a given varName and
 * varF
 *
 * \author Marian Albers, Nov 2015
 */
template <ZFSInt nDim>
void ZFSStrctrdInterpolation<nDim>::interpolateField(ZFSString varName, ZFSFloat* receiverVar)
{
  //first load the current variable
  loadDonorVariable(varName);
  for(ZFSInt cellId=0; cellId<m_noReceiverCells; cellId++) {
    if(cellId % 50000 == 0 && m_domainId==0) {
      cout << "Variable " << varName << " interpolation progress: " << (ZFSInt)((ZFSFloat)cellId/((ZFSFloat)m_noReceiverCells)*100.0) << " percent"  << endl;
    }

    if(m_hasInterpolationPartners[cellId]) {
      ZFSFloat transformedPoint[3] = {m_transformedReceiverPoints[0][cellId],
                                      m_transformedReceiverPoints[1][cellId],
                                      m_transformedReceiverPoints[2][cellId]};

      //interpolate variables at transformed coordinate
      ZFSId currentBlockId = getBlockId(m_donorOriginId[cellId]);
      trilinearInterpolation(transformedPoint, m_donorOriginId[cellId], cellId, receiverVar, currentBlockId);
    } else {
      //fallback to nearest neighbour interpolation
      nearestNeighbourInterpolation(m_donorOriginId[cellId], cellId, receiverVar);
    }
  }
}

/** \brief Prepares interpolation for field
* For a given 3D coordinate field the method
* computes the transformed Points and stores them
* to perform interpolations later
*
* \author Marian Albers, Nov 2015
*/
template <ZFSInt nDim>
void ZFSStrctrdInterpolation<nDim>::prepareInterpolationField(ZFSInt* noReceiverCells, ZFSFloat** receiverCoordinates)
{
  m_isFieldInterpolation = true;
  m_noReceiverCells = noReceiverCells[0]*noReceiverCells[1]*noReceiverCells[2];
  ZFSInt noTrilinear = 0;
  ZFSInt noFallback = 0;
  zfsAlloc(m_donorOriginId, m_noReceiverCells, __CALLING_FUNCTION__, -1, "m_donorOriginId");
  zfsAlloc(m_transformedReceiverPoints, nDim, m_noReceiverCells, __CALLING_FUNCTION__, F0, "m_interpolationCoefficients");
  zfsAlloc(m_hasInterpolationPartners, m_noReceiverCells, __CALLING_FUNCTION__, "m_hasInterpolationPartners");
  zfsAlloc(m_donorVar, m_totalNoDonorCells, "m_donorVariables", F0, __CALLING_FUNCTION__);

  for(ZFSInt cellId=0; cellId<m_noReceiverCells; cellId++) {
    if(cellId % 50000 == 0) {
      cout << "Interpolation progress: " << (ZFSInt)((ZFSFloat)cellId/((ZFSFloat)m_noReceiverCells )*100.0) << " percent"  << endl;
    }
    m_currentReceiverId = cellId;
    ZFSFloat intPoint[3] = {receiverCoordinates[0][cellId],receiverCoordinates[1][cellId],receiverCoordinates[2][cellId]};
    //now go through own, fine cells and look for closest partner neighbour
    ZFSFloat distance = 0;
    Point<3> pt(intPoint[0], intPoint[1], intPoint[2]);
    //find point on the grid that is closest to intPoint
    ZFSId centerCellId = m_donorTree->nearest(pt, distance);

    //now eight hexahedron could be candidates for a new home for intPoint,
    //all around centerCellId
    ZFSInt hexahedronOriginId = -1;
    const ZFSInt maxNghbrRadius = 4;
    for(ZFSId nghbrRadius = 1; nghbrRadius < maxNghbrRadius; nghbrRadius++) {
      hexahedronOriginId = findSurroundingHexahedron(intPoint, centerCellId, nghbrRadius);
      if(hexahedronOriginId != -1) {
        break;
      }
    }

    if(hexahedronOriginId != -1) {
      ZFSFloat transformedPoint[3] = {F0,F0,F0};
      //now hexahedronOriginId is the id of the hexahedron in which intPoint is immersed
      transformPoint(hexahedronOriginId, intPoint, transformedPoint);
      //interpolate variables at transformed coordinate

      for(ZFSId dim=0; dim<nDim; dim++) {
        m_transformedReceiverPoints[dim][cellId] = transformedPoint[dim];
      }
      m_donorOriginId[cellId] = hexahedronOriginId;
      m_hasInterpolationPartners[cellId] = true;
      noTrilinear++;
    } else {
      //fallback to nearest neighbour interpolation
      m_hasInterpolationPartners[cellId] = false;
      m_donorOriginId[cellId] = centerCellId;
      noFallback++;
    }
  }

  ZFSInt noLocal[3] = {noTrilinear, noFallback, m_noReceiverCells };
  ZFSInt noGlobal[3] = {0,0,0};
  MPI_Allreduce(noLocal, noGlobal, 3, MPI_INT, MPI_SUM, m_zfsStrctrdComm);
  if(m_domainId==0) {
    cout << "Trilinear: " << noGlobal[0] << " (" << 100.0*((ZFSFloat)noGlobal[0])/((ZFSFloat)noGlobal[2]) << "%) "
         << "Fallback: " << noGlobal[1]  << " (" <<  100.0*((ZFSFloat)noGlobal[1])/((ZFSFloat)noGlobal[2]) << "%)" << endl;
  }
}

/** \brief Prepares interpolation neighbours and coefficients
* For a given number of points the methods computes
* the interpolation partners and coefficients for later use
* use together with interpolateVariables()
*
* \author Marian Albers, Jan 2015
*/
template <ZFSInt nDim>
void ZFSStrctrdInterpolation<nDim>::prepareInterpolation(ZFSInt noReceiverCells, ZFSFloat** receiverCellCoordinates, ZFSInt* interpolationPartner)
{
  m_isFieldInterpolation = false;
  ZFSInt noTrilinear = 0;
  ZFSInt noFallback = 0;
  m_noReceiverCells = noReceiverCells;
  zfsAlloc(m_donorOriginId, m_noReceiverCells, __CALLING_FUNCTION__, "m_donorOriginId");
  zfsAlloc(m_interpolationCoefficients, m_noReceiverCells, 8, __CALLING_FUNCTION__, "m_interpolationCoefficients");
  zfsAlloc(m_hasInterpolationPartners, m_noReceiverCells, __CALLING_FUNCTION__, "m_hasInterpolationPartners");

  for(ZFSId receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
    if(m_domainId==0 && receiverId % 10000 == 0) {cout << "receiver no: " << receiverId << endl;}
    //now go through own, fine cells and look for closest partner neighbour
    ZFSFloat intPoint[3] = {receiverCellCoordinates[0][receiverId],receiverCellCoordinates[1][receiverId],receiverCellCoordinates[2][receiverId]};
    Point<3> pt(intPoint[0], intPoint[1], intPoint[2]);
    ZFSFloat dist = 0;
    ZFSId closestCellId = m_donorTree->nearest(pt, dist);

    //now eight hexahedron could be candidates for a new home for intPoint,
    //all around closestCellId
    //ZFSId hexahedronOriginId = findSurroundingHexahedron(intPoint, closestCellId, 1);
    ZFSInt hexahedronOriginId = -1;
    const ZFSInt maxNghbrRadius = 4;
    for(ZFSId nghbrRadius = 1; nghbrRadius < maxNghbrRadius; nghbrRadius++) {
      hexahedronOriginId = findSurroundingHexahedron(intPoint, closestCellId, nghbrRadius);
      if(hexahedronOriginId != -1) {
        break;
      }
    }

    if(hexahedronOriginId != -1) {
      m_hasInterpolationPartners[receiverId] = true;
      m_donorOriginId[receiverId] = hexahedronOriginId;
      ZFSFloat transformedPoint[3];
      //now hexahedronOriginId is the id of the hexahedron in which intPoint is immersed
      transformPoint(hexahedronOriginId, intPoint, transformedPoint);
      //interpolate variables at transformed coordinate
      computeInterpolationCoefficients(transformedPoint, receiverId);
      noTrilinear++;
    } else {
      //fallback to nearest neighbour interpolation
      m_hasInterpolationPartners[receiverId] = false;
      m_donorOriginId[receiverId] = closestCellId;
      noFallback++;
    }
    interpolationPartner[receiverId]=m_hasInterpolationPartners[receiverId];
  }

  ZFSInt noLocal[3] = {noTrilinear, noFallback, noReceiverCells };
  ZFSInt noGlobal[3] = {0,0,0};
  MPI_Allreduce(noLocal, noGlobal, 3, MPI_INT, MPI_SUM, m_zfsStrctrdComm);
  if(m_domainId==0) {
    cout << "Trilinear: " << noGlobal[0] << " (" << 100.0*((ZFSFloat)noGlobal[0])/((ZFSFloat)noGlobal[2]) << "%)"
         << "Fallback: " << noGlobal[1]  << " (" <<  100.0*((ZFSFloat)noGlobal[1])/((ZFSFloat)noGlobal[2]) << "%)" << endl;
  }
}

/** \brief computes the interpolation coefficients
 * Interpolates all variables with precomputed interpolation coefficients
 *
 * \author Marian Albers, Jan 2015
*/

template <ZFSInt nDim>   //junoh
void ZFSStrctrdInterpolation<nDim>::prepareZonalInterpolation(ZFSInt noReceiverCells, ZFSFloat** receiverCellCoordinates, ZFSInt* interpolationPartner, ZFSBool hasInterpolationPartnerDomain)
{
 
  m_isFieldInterpolation = false;
  ZFSInt noTrilinear = 0;
  ZFSInt noFallback = 0;
  m_noReceiverCells = noReceiverCells;
  m_hasInterpolationPartnerDomain = hasInterpolationPartnerDomain;
  
  zfsAlloc(m_donorOriginId, m_noReceiverCells, __CALLING_FUNCTION__, "m_donorOriginId");
  zfsAlloc(m_interpolationCoefficients, m_noReceiverCells, 8, __CALLING_FUNCTION__, "m_interpolationCoefficients");
  zfsAlloc(m_hasInterpolationPartnersZonal, m_noReceiverCells,__CALLING_FUNCTION__,0,"m_hasInterpolationPartnersZonal");
  zfsAlloc(m_hasInterpolationPartnersZonalGlobal, m_noReceiverCells,__CALLING_FUNCTION__,0,"m_hasInterpolationPartnersZonalGlobal");
  zfsAlloc(m_donorDistance, m_noReceiverCells,__CALLING_FUNCTION__,"m_donorDistance");
  zfsAlloc(m_globalDonorDistanceMin,m_noReceiverCells,__CALLING_FUNCTION__,F0,"m_globalDonorDistanceMin");

  if(m_hasInterpolationPartnerDomain==true){
    for(ZFSId receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
      // if(m_domainId==0 && receiverId % 10000 == 0) {cout << "receiver no: " << receiverId << endl;}
      //now go through own, fine cells and look for closest partner neighbour
      ZFSFloat intPoint[3] = {receiverCellCoordinates[0][receiverId],receiverCellCoordinates[1][receiverId],receiverCellCoordinates[2][receiverId]};
      Point<3> pt(intPoint[0], intPoint[1], intPoint[2]);
      ZFSFloat dist = 0;
      ZFSId closestCellId;
      closestCellId = m_donorTree->nearest(pt, dist);
      //now eight hexahedron could be candidates for a new home for intPoint,
      //all around closestCellId
      //ZFSId hexahedronOriginId = findSurroundingHexahedron(intPoint, closestCellId, 1);
      ZFSInt hexahedronOriginId = -1;
      const ZFSInt maxNghbrRadius = 4;
      for(ZFSId nghbrRadius = 1; nghbrRadius < maxNghbrRadius; nghbrRadius++) {
	hexahedronOriginId = findSurroundingHexahedron(intPoint, closestCellId, nghbrRadius);
	if(hexahedronOriginId != -1) {
	  break;
	}
      }
      if(hexahedronOriginId != -1) {
	m_hasInterpolationPartnersZonal[receiverId] = true;
	m_donorOriginId[receiverId] = hexahedronOriginId;
	ZFSFloat transformedPoint[3];
	//now hexahedronOriginId is the id of the hexahedron in which intPoint is immersed
	transformPoint(hexahedronOriginId, intPoint, transformedPoint);
	//interpolate variables at transformed coordinate
	computeInterpolationCoefficients(transformedPoint, receiverId);
	noTrilinear++;
      } else {
	//fallback to nearest neighbour interpolation
	m_hasInterpolationPartnersZonal[receiverId] = false;
	m_donorOriginId[receiverId] = closestCellId;
	m_donorDistance[receiverId] = dist;
	noFallback++;
      }
    }
  }else{
    for(ZFSId receiverId=0; receiverId<m_noReceiverCells; receiverId++){
      m_donorDistance[receiverId] = 1000;    //avoid to pick interpolation cells in receiverDomain using MPI_Allreduce(MIN) below
      m_hasInterpolationPartnersZonal[receiverId] = false;
    }
  }
      if(m_domainId==0){
	cout<<" nearst Distance calculating... "<<endl;
  }
  MPI_Allreduce(m_hasInterpolationPartnersZonal, m_hasInterpolationPartnersZonalGlobal, m_noReceiverCells, MPI_INT, MPI_SUM, m_zfsStrctrdComm);

  MPI_Allreduce(m_donorDistance, m_globalDonorDistanceMin, m_noReceiverCells, MPI_DOUBLE, MPI_MIN,m_zfsStrctrdComm);
   
  if(m_hasInterpolationPartnerDomain==true){
    for(ZFSId receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
      if(m_hasInterpolationPartnersZonalGlobal[receiverId]==0 && approx(m_donorDistance[receiverId],m_globalDonorDistanceMin[receiverId],m_eps)){
	m_hasInterpolationPartnersZonal[receiverId]=true;
      }
    }
  }

  for(ZFSId receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {    
    interpolationPartner[receiverId]=m_hasInterpolationPartnersZonal[receiverId];
  }

     
  // ZFSInt noLocal[3] = {noTrilinear, noFallback, noReceiverCells };
  // ZFSInt noGlobal[3] = {0,0,0};
  // MPI_Allreduce(noLocal, noGlobal, 3, MPI_INT, MPI_SUM, m_zfsStrctrdComm);
  // if(m_domainId==0) {
  //   // cout << "Trilinear: " << noGlobal[0] << " (" << 100.0*((ZFSFloat)noGlobal[0])/((ZFSFloat)noGlobal[2]) << "%)"
  //   //      << "Fallback: " << noGlobal[1]  << " (" <<  100.0*((ZFSFloat)noGlobal[1])/((ZFSFloat)noGlobal[2]) << "%)" << endl;
  // }
}



template <ZFSInt nDim>  //junoh
ZFSFloat ZFSStrctrdInterpolation<nDim>::interpolateVariableZonal(ZFSFloat* donorVars, ZFSInt cellIdBC)
{
  ZFSFloat interpolatedVariable=0.0;

  if(m_hasInterpolationPartnersZonalGlobal[cellIdBC]>0) {
    interpolatedVariable = getInterpolatedVariableZonal(donorVars, cellIdBC);
  } else if(m_hasInterpolationPartnersZonalGlobal[cellIdBC]==0) {
    interpolatedVariable = donorVars[m_donorOriginId[cellIdBC]];
 }
  return interpolatedVariable;
}






template <ZFSInt nDim>
void ZFSStrctrdInterpolation<nDim>::interpolateVariables(ZFSFloat** interpolatedVariables)
{
  for(ZFSId receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
    if(m_hasInterpolationPartners[receiverId]) {
      for(ZFSId var = 0; var<m_noDonorVariables; var++) {
        interpolatedVariables[var][receiverId] = getInterpolatedVariable(receiverId, var);
      }
    } else {
      for(ZFSId var = 0; var<m_noDonorVariables; var++) {
        interpolatedVariables[var][receiverId] = m_donorVariables[var][m_donorOriginId[receiverId]];
      }
    }
  }
}

/** \brief Finds surrounding hexahedron for point
 * For a given cellId it finds out which of
 * the 8 eight surrounding hexahedrons contains
 * intPoint. If point is inside any of the 6
 * possible tetraeders inside the hexahedron it
 * is inside the hexahedron, the origin cellId
 * is returned then.
 *
 * \author Marian Albers, Nov 2015
*/

// template <ZFSInt nDim>   //junoh but not used anywhere
// void ZFSStrctrdInterpolation<nDim>::interpolateVariablesNoFallback(ZFSFloat* interpolatedVariables)
// {
//   for(ZFSId receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
//     if(m_hasInterpolationPartners[receiverId]) {
//       for(ZFSId var = 0; var<m_noDonorVariables; var++) {
//         interpolatedVariables[var*m_noReceiverCells+receiverId] = getInterpolatedVariable(receiverId, var);
//       }
//     }
//   }
// }




template <ZFSInt nDim>
inline ZFSId ZFSStrctrdInterpolation<nDim>::findSurroundingHexahedron(ZFSFloat intPoint[3],
                                                                      ZFSId centerCellId,
                                                                      ZFSId stencil) {
  ZFSFloat v1[3] = {F0,F0,F0};
  ZFSFloat v2[3] = {F0,F0,F0};
  ZFSFloat v3[3] = {F0,F0,F0};
  ZFSFloat vp[3] = {F0,F0,F0};
  ZFSFloat vn[3] = {F0,F0,F0};

  ZFSBool isInside = false;

  for(ZFSId i = -stencil; i<stencil; i++) {
    for(ZFSId j = -stencil; j<stencil; j++) {
      for(ZFSId k = -stencil; k<stencil; k++) {
        //find out current blockId
        const ZFSId currentBlockId = getBlockId(centerCellId);
        const ZFSId currentOffset = m_donorBlockOffsets[currentBlockId];
        const ZFSId IJK = getCellIdfromCell(centerCellId,i,j,k,currentBlockId);

        //compute i,j,k from cellId, but use local cellIds
        ZFSId noLocalCellsDir[3] = {m_noDonorCellsDir[currentBlockId][0]-currentOffset,
                                    m_noDonorCellsDir[currentBlockId][1]-currentOffset,
                                    m_noDonorCellsDir[currentBlockId][2]-currentOffset};
        ZFSInt trueI = (IJK%(noLocalCellsDir[2]*noLocalCellsDir[1]))%noLocalCellsDir[2];
        ZFSInt trueJ = ((IJK - trueI)/noLocalCellsDir[2])%noLocalCellsDir[1];
        ZFSInt trueK = ((IJK - trueI)/noLocalCellsDir[2] - trueJ)/noLocalCellsDir[1];

        //check if inside regular grid
        if(trueI < 0 || trueI >= noLocalCellsDir[2]-1 || 
           trueJ < 0 || trueJ >= noLocalCellsDir[1]-1 || 
           trueK < 0 || trueK >= noLocalCellsDir[0]-1) {
          continue;
        }

        //now loop over all 6 possible tetragonal pyramids
        for(ZFSId tetra=0; tetra<6; tetra++) {
          //and now loop over all 4 sides of the pyramid
          isInside = true;
          for(ZFSId side=0; side<4; side++) {
            ZFSId ixp=(side+1)%4;
            ZFSId ixp2=(ixp+1)%4;
            ZFSId ixp3=(ixp2+1)%4;

            //compute vectors
            for(ZFSId dim=0; dim<nDim; dim++) {
              v1[dim] = m_donorCoordinates[dim][getCellIdfromCell(IJK,ic(tetra,side,0),ic(tetra,side,1),ic(tetra,side,2),currentBlockId)] -
                        m_donorCoordinates[dim][getCellIdfromCell(IJK,ic(tetra,ixp ,0),ic(tetra,ixp ,1),ic(tetra,ixp ,2),currentBlockId)];
              v2[dim] = m_donorCoordinates[dim][getCellIdfromCell(IJK,ic(tetra,ixp2,0),ic(tetra,ixp2,1),ic(tetra,ixp2,2),currentBlockId)] -
                        m_donorCoordinates[dim][getCellIdfromCell(IJK,ic(tetra,ixp ,0),ic(tetra,ixp ,1),ic(tetra,ixp ,2),currentBlockId)];
              v3[dim] = m_donorCoordinates[dim][getCellIdfromCell(IJK,ic(tetra,ixp3,0),ic(tetra,ixp3,1),ic(tetra,ixp3,2),currentBlockId)] -
                        m_donorCoordinates[dim][getCellIdfromCell(IJK,ic(tetra,ixp ,0),ic(tetra,ixp ,1),ic(tetra,ixp ,2),currentBlockId)];
              vp[dim] = intPoint[dim] -
                        m_donorCoordinates[dim][getCellIdfromCell(IJK,ic(tetra,ixp ,0),ic(tetra,ixp ,1),ic(tetra,ixp ,2),currentBlockId)];
            }

            //compute perpendicular vector vn
            crossProduct(vn,v1,v2);

            //check direction of normal vector and correct if necessary
            if(scalarProduct(v3,vn) > m_eps) {
              //invert direction of vn
              for(ZFSId dim=0; dim<nDim; dim++) {
                vn[dim] = - vn[dim];
              }
            }

            //check scalar product, if vp*vn > 0 then point is outside of this pyramid
            if(scalarProduct(vp,vn) > m_eps) {
              //not inside pyramid
              isInside = false;
              break;
            }
          }

          //point is inside the current tetra and
          //thus inside the current hexahedron
          if(isInside) {
            return IJK;
          }
        }
      }
    }
  }

  //surprisingly the point isn't inside any of the 8 hexahedrons, return invalid cellId
  return -1;
}

/** \brief Transforms point from physical to computational space
 *
 * Transforms a given point from physical space
 * to computational space with help of the surrounding
 * eight points (which define a hexahedron). Iterative
 * solution with Newton to solve system of equations.
 * Benek (1987) - Chimera - A grid-embedding technique
 *
 * /author Marian Albers, Nov 23, 2015
 */
template <ZFSInt nDim>
inline void ZFSStrctrdInterpolation<nDim>::transformPoint(ZFSId hexOrigin,
                                                          ZFSFloat intPoint[3],
                                                          ZFSFloat transformedPoint[3]) {

  ZFSFloatScratchSpace a(3,8, __CALLING_FUNCTION__, "a");
  const ZFSId currentBlockId = getBlockId(hexOrigin);

  for(ZFSId dim=0; dim<nDim; dim++) {
    a(dim,0) = m_donorCoordinates[dim][hexOrigin];
    a(dim,1) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin,1,0,0,currentBlockId)] - a(dim,0);
    a(dim,2) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin,0,1,0,currentBlockId)] - a(dim,0);
    a(dim,3) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin,0,0,1,currentBlockId)] - a(dim,0);
    a(dim,4) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin,1,1,0,currentBlockId)] - a(dim,0) - a(dim,1) - a(dim,2);
    a(dim,5) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin,1,0,1,currentBlockId)] - a(dim,0) - a(dim,1) - a(dim,3);
    a(dim,6) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin,0,1,1,currentBlockId)] - a(dim,0) - a(dim,2) - a(dim,3);
    a(dim,7) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin,1,1,1,currentBlockId)] - a(dim,0) - a(dim,1) - a(dim,2) - a(dim,3) - a(dim,4) - a(dim,5) - a(dim,6) - a(dim,7);
  }

  //set initial values 0.5
  ZFSFloatScratchSpace dxezdxyz(3,3, __CALLING_FUNCTION__, "dxezdxyz");
  ZFSFloatScratchSpace rhs(3, __CALLING_FUNCTION__, "rhs");

  ZFSFloat xi = F1B2;
  ZFSFloat eta = F1B2;
  ZFSFloat zeta = F1B2;

  ZFSFloat xicor,etacor,zetacor;

  ZFSInt noIterations = 50;
  for(ZFSId newton=0; newton < noIterations; newton++) {
    for(ZFSId dim=0; dim<nDim; dim++) {
      dxezdxyz(dim,0) = a(dim,1) + a(dim,4)*eta + a(dim,5)*zeta + a(dim,7)*eta*zeta;
      dxezdxyz(dim,1) = a(dim,2) + a(dim,4)*xi  + a(dim,6)*zeta + a(dim,7)*xi*zeta;
      dxezdxyz(dim,2) = a(dim,3) + a(dim,5)*xi  + a(dim,6)*eta  + a(dim,7)*xi*eta;

      rhs[dim] = intPoint[dim] - (a(dim,0) + a(dim,1)*xi + a(dim,2)*eta + a(dim,3)*zeta + a(dim,4)*xi*eta + a(dim,5)*xi*zeta + a(dim,6)*eta*zeta + a(dim,7)*xi*eta*zeta);
    }

    ZFSFloat jac = dxezdxyz(0,0)*(dxezdxyz(1,1)*dxezdxyz(2,2) - dxezdxyz(1,2)*dxezdxyz(2,1)) +
      dxezdxyz(0,1)*(dxezdxyz(1,2)*dxezdxyz(2,0) - dxezdxyz(1,0)*dxezdxyz(2,2)) +
      dxezdxyz(0,2)*(dxezdxyz(1,0)*dxezdxyz(2,1) - dxezdxyz(1,1)*dxezdxyz(2,0));

    if(fabs(jac) < m_eps ) {
      xicor = F0;
      etacor = F0;
      zetacor = F0;
    } else {
      xicor = (rhs[0]*(dxezdxyz(1,1)*dxezdxyz(2,2) - dxezdxyz(1,2)*dxezdxyz(2,1)) +
               dxezdxyz(0,1)*(dxezdxyz(1,2)*rhs[2] - rhs[1]*dxezdxyz(2,2)) +
               dxezdxyz(0,2)*(rhs[1]*dxezdxyz(2,1) - dxezdxyz(1,1)*rhs[2]))/jac;
      etacor = (dxezdxyz(0,0)*(rhs[1]*dxezdxyz(2,2) - dxezdxyz(1,2)*rhs[2]) +
                rhs[0]*(dxezdxyz(1,2)*dxezdxyz(2,0) - dxezdxyz(1,0)*dxezdxyz(2,2)) +
                dxezdxyz(0,2)*(dxezdxyz(1,0)*rhs[2] - rhs[1]*dxezdxyz(2,0)))/jac;
      zetacor = (dxezdxyz(0,0)*(dxezdxyz(1,1)*rhs[2] - rhs[1]*dxezdxyz(2,1)) +
                 dxezdxyz(0,1)*(rhs[1]*dxezdxyz(2,0) - dxezdxyz(1,0)*rhs[2]) +
                 rhs[0]*(dxezdxyz(1,0)*dxezdxyz(2,1) - dxezdxyz(1,1)*dxezdxyz(2,0)))/jac;
    }

    xi = xi + xicor;
    eta = eta + etacor;
    zeta = zeta + zetacor;

    ZFSFloat sumcor = fabs(xicor) + fabs(etacor) + fabs(zetacor);

    if(sumcor < m_eps) {
      break;
    }
  }

  xi = zfsMAX(zfsMIN(xi, F1), F0);
  eta = zfsMAX(zfsMIN(eta, F1), F0);
  zeta = zfsMAX(zfsMIN(zeta, F1), F0);

  transformedPoint[0] = xi;
  transformedPoint[1] = eta;
  transformedPoint[2] = zeta;
}

/** \brief Trilinear Interpolation for field
 *
 * Performs a trilinear interpolation in computational
 * space for a given point dx. The 8 values are specified
 * at the edges of a uniform cube around point dx
 *
 * \author Marian Albers, Nov 2015
 */
template <ZFSInt nDim>
inline void ZFSStrctrdInterpolation<nDim>::trilinearInterpolation(ZFSFloat dx[3],
                                                                  ZFSInt hexOrigin,
                                                                  ZFSId receiverId,
                                                                  ZFSFloat* receiverVar,
                                                                  ZFSId currentBlockId) {

  ZFSFloatScratchSpace v(8, __CALLING_FUNCTION__, "v");

  v[0]=dx[0]*dx[1]*dx[2];
  v[1]=(F1-dx[0])*dx[1]*dx[2];
  v[2]=dx[0]*(F1-dx[1])*dx[2];
  v[3]=dx[0]*dx[1]*(F1-dx[2]);
  v[4]=(F1-dx[0])*(F1-dx[1])*dx[2];
  v[5]=(F1-dx[0])*dx[1]*(F1-dx[2]);
  v[6]=dx[0]*(F1-dx[1])*(F1-dx[2]);
  v[7]=(F1-dx[0])*(F1-dx[1])*(F1-dx[2]);

  //now the actual interpolation
  receiverVar[receiverId] = 
    v[0]*m_donorVar[getCellIdfromCell(hexOrigin,1,1,1,currentBlockId)] +
    v[1]*m_donorVar[getCellIdfromCell(hexOrigin,0,1,1,currentBlockId)] +
    v[2]*m_donorVar[getCellIdfromCell(hexOrigin,1,0,1,currentBlockId)] +
    v[3]*m_donorVar[getCellIdfromCell(hexOrigin,1,1,0,currentBlockId)] +
    v[4]*m_donorVar[getCellIdfromCell(hexOrigin,0,0,1,currentBlockId)] +
    v[5]*m_donorVar[getCellIdfromCell(hexOrigin,0,1,0,currentBlockId)] +
    v[6]*m_donorVar[getCellIdfromCell(hexOrigin,1,0,0,currentBlockId)] +
    v[7]*m_donorVar[getCellIdfromCell(hexOrigin,0,0,0,currentBlockId)];
}

/** \brief Trilinear interpolation for single point
 *
 * Performs a trilinear interpolation in computational
 * space for a given point dx. The 8 values are specified
 * at the edges of a uniform cube around point dx
 *
 * \author Marian Albers, Nov 2015
 */
template <ZFSInt nDim>
inline void ZFSStrctrdInterpolation<nDim>::trilinearInterpolation(ZFSFloat dx[3],
                                                                  ZFSInt hexOrigin,
                                                                  ZFSFloat* receiverVariables,
                                                                  ZFSId currentBlockId) {

  ZFSFloatScratchSpace v(8, __CALLING_FUNCTION__, "v");

  v[0]=dx[0]*dx[1]*dx[2];
  v[1]=(F1-dx[0])*dx[1]*dx[2];
  v[2]=dx[0]*(F1-dx[1])*dx[2];
  v[3]=dx[0]*dx[1]*(F1-dx[2]);
  v[4]=(F1-dx[0])*(F1-dx[1])*dx[2];
  v[5]=(F1-dx[0])*dx[1]*(F1-dx[2]);
  v[6]=dx[0]*(F1-dx[1])*(F1-dx[2]);
  v[7]=(F1-dx[0])*(F1-dx[1])*(F1-dx[2]);

  //now the actual interpolation
  for(ZFSInt var=0; var<m_noDonorVariables; var++) {
    receiverVariables[var] =
      v[0]*m_donorVariables[var][getCellIdfromCell(hexOrigin,1,1,1,currentBlockId)] +
      v[1]*m_donorVariables[var][getCellIdfromCell(hexOrigin,0,1,1,currentBlockId)] +
      v[2]*m_donorVariables[var][getCellIdfromCell(hexOrigin,1,0,1,currentBlockId)] +
      v[3]*m_donorVariables[var][getCellIdfromCell(hexOrigin,1,1,0,currentBlockId)] +
      v[4]*m_donorVariables[var][getCellIdfromCell(hexOrigin,0,0,1,currentBlockId)] +
      v[5]*m_donorVariables[var][getCellIdfromCell(hexOrigin,0,1,0,currentBlockId)] +
      v[6]*m_donorVariables[var][getCellIdfromCell(hexOrigin,1,0,0,currentBlockId)] +
      v[7]*m_donorVariables[var][getCellIdfromCell(hexOrigin,0,0,0,currentBlockId)];

  }
}


/** \brief Computes trilinear interpolation coefficients
 *
 * \author Marian Albers, Jan 2016
 */
template <ZFSInt nDim>
inline void ZFSStrctrdInterpolation<nDim>::computeInterpolationCoefficients(ZFSFloat dx[3],
                                                                            ZFSInt receiverId) {
  m_interpolationCoefficients[receiverId][0]=dx[0]*dx[1]*dx[2];
  m_interpolationCoefficients[receiverId][1]=(F1-dx[0])*dx[1]*dx[2];
  m_interpolationCoefficients[receiverId][2]=dx[0]*(F1-dx[1])*dx[2];
  m_interpolationCoefficients[receiverId][3]=dx[0]*dx[1]*(F1-dx[2]);
  m_interpolationCoefficients[receiverId][4]=(F1-dx[0])*(F1-dx[1])*dx[2];
  m_interpolationCoefficients[receiverId][5]=(F1-dx[0])*dx[1]*(F1-dx[2]);
  m_interpolationCoefficients[receiverId][6]=dx[0]*(F1-dx[1])*(F1-dx[2]);
  m_interpolationCoefficients[receiverId][7]=(F1-dx[0])*(F1-dx[1])*(F1-dx[2]);
}

/** \brief Trilinear interpolation for single point
 *
 * Performs trilinear interpolation from previously
 * computed interpolation coefficients
 * (see computeInterpolationCoefficients)
 *
 * \author Marian Albers, Jan 2016
 */
template <ZFSInt nDim>
ZFSFloat ZFSStrctrdInterpolation<nDim>::getInterpolatedVariable(ZFSInt receiverId,
                                                                ZFSInt var) {

  const ZFSId donorOrigin = m_donorOriginId[receiverId];
  const ZFSId currentBlockId = getBlockId(donorOrigin);
  const ZFSFloat interpolatedVar  =
    m_interpolationCoefficients[receiverId][0]*m_donorVariables[var][getCellIdfromCell(donorOrigin,1,1,1,currentBlockId)] +
    m_interpolationCoefficients[receiverId][1]*m_donorVariables[var][getCellIdfromCell(donorOrigin,0,1,1,currentBlockId)] +
    m_interpolationCoefficients[receiverId][2]*m_donorVariables[var][getCellIdfromCell(donorOrigin,1,0,1,currentBlockId)] +
    m_interpolationCoefficients[receiverId][3]*m_donorVariables[var][getCellIdfromCell(donorOrigin,1,1,0,currentBlockId)] +
    m_interpolationCoefficients[receiverId][4]*m_donorVariables[var][getCellIdfromCell(donorOrigin,0,0,1,currentBlockId)] +
    m_interpolationCoefficients[receiverId][5]*m_donorVariables[var][getCellIdfromCell(donorOrigin,0,1,0,currentBlockId)] +
    m_interpolationCoefficients[receiverId][6]*m_donorVariables[var][getCellIdfromCell(donorOrigin,1,0,0,currentBlockId)] +
    m_interpolationCoefficients[receiverId][7]*m_donorVariables[var][getCellIdfromCell(donorOrigin,0,0,0,currentBlockId)];
  
    return interpolatedVar;
}

template <ZFSInt nDim>//junoh
ZFSFloat ZFSStrctrdInterpolation<nDim>::getInterpolatedVariableZonal(ZFSFloat* donorVars, ZFSInt receiverId) {

  const ZFSId donorOrigin = m_donorOriginId[receiverId];
  const ZFSId currentBlockId = getBlockId(donorOrigin);
  const ZFSFloat interpolatedVar  =
    m_interpolationCoefficients[receiverId][0]*donorVars[getCellIdfromCell(donorOrigin,1,1,1,currentBlockId)] +
    m_interpolationCoefficients[receiverId][1]*donorVars[getCellIdfromCell(donorOrigin,0,1,1,currentBlockId)] +
    m_interpolationCoefficients[receiverId][2]*donorVars[getCellIdfromCell(donorOrigin,1,0,1,currentBlockId)] +
    m_interpolationCoefficients[receiverId][3]*donorVars[getCellIdfromCell(donorOrigin,1,1,0,currentBlockId)] +
    m_interpolationCoefficients[receiverId][4]*donorVars[getCellIdfromCell(donorOrigin,0,0,1,currentBlockId)] +
    m_interpolationCoefficients[receiverId][5]*donorVars[getCellIdfromCell(donorOrigin,0,1,0,currentBlockId)] +
    m_interpolationCoefficients[receiverId][6]*donorVars[getCellIdfromCell(donorOrigin,1,0,0,currentBlockId)] +
    m_interpolationCoefficients[receiverId][7]*donorVars[getCellIdfromCell(donorOrigin,0,0,0,currentBlockId)];
  
    return interpolatedVar;
}



/** \brief Loads a grid file
 *
 * Loads the grid file given in the property file.
 * No domain decompositioning, all domains read
 * whole grid file as it is necessary for neighbour
 * search.
 *
 * \author Marian Albers, Nov 2015
 */
template <ZFSInt nDim>
void ZFSStrctrdInterpolation<nDim>::loadDonorGrid()
{
  if(m_domainId==0) {cout << "Reading in donor grid file..." << endl;}

  /*! \page propertyPage1
    \section donorGridName
    <code>ZFSInt ZFSStrctrdInterpolation::m_donorGridName </code>\n
    default = <code> "" </code>\n \n
    Name of the donor grid file to interpolate from.\n
    Possible values are:\n
    <ul>
    <li>String containing file name</li>
    </ul>
    Keywords: <i>INTERPOLATION, STRCTRD</i>
  */
  ZFSString donorGridName =*(ZFSContext::getProperty("donorGridName", 0, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString(0));

  ZFSFloat translation [3] = {F0,F0,F0};
  ZFSFloat scale[3] = {F1,F1,F1};

  /*! \page propertyPage1
    \section donorTranslation
    <code>ZFSInt ZFSStrctrdInterpolation::m_donorTranslation </code>\n
    default = <code> 0.0, 0.0, 0.0 </code>\n \n
    Translation of the donor grid in 3 space dimensions.\n
    Possible values are:\n
    <ul>
    <li>Float</li>
    </ul>
    Keywords: <i>INTERPOLATION, STRCTRD</i>
  */
  if(ZFSContext::propertyExists("donorTranslation", 0 )){
    for(ZFSInt dim=0; dim<nDim; dim++){
      translation[dim] =*(ZFSContext::getProperty("donorTranslation", 0, __CALLING_FUNCTION__,(ZFSFloat*)  &translation[dim] , nDim )->asFloat(dim));
    }
  }

  /*! \page propertyPage1
    \section donorScale
    <code>ZFSInt ZFSStrctrdInterpolation::m_donorScale </code>\n
    default = <code> 1.0, 1.0, 1.0 </code>\n \n
    Scaling of the donor grid.\n
    Possible values are:\n
    <ul>
    <li>Float > 0.0</li>
    </ul>
    Keywords: <i>INTERPOLATION, STRCTRD</i>
  */
  if(ZFSContext::propertyExists("donorScale", 0 )){
    for(ZFSInt dim=0; dim<nDim; dim++){
      scale[dim] =*(ZFSContext::getProperty("donorScale", 0, __CALLING_FUNCTION__,(ZFSFloat*)  &scale[dim] , nDim )->asFloat(dim));
    }
  }

  if(m_domainId==0) {
    cout << "Translating donorGrid by deltaX: " << translation[0] << " deltaY: " << translation[1] << " deltaZ: " << translation[2] << endl;
    cout << "Scaling donorGrid by scaleX: " << scale[0] << " scaleY: " << scale[1] << " scaleZ: " << scale[2] << endl;
  }

  //file id to access the file
  ZFSId file_id = -1;
  file_id = io_openfile("hdf5", donorGridName.c_str(), "collective", m_zfsStrctrdComm);
  if(file_id<0) zfsTerm(1, __CALLING_FUNCTION__, "ERROR: grid File could not be opened :-(");
  //create the string to contain the datasetname in the file
  ZFSInt inputBlockId = 0;
  ZFSString sBlockName = "/block";
  stringstream dummy1;
  dummy1 << inputBlockId <<"/";
  sBlockName += dummy1.str();
  m_noDonorDims =  io_getDatasetDimension(file_id, sBlockName.c_str(), "x");

  if(m_domainId==0) {cout << "Donor grid has " << m_noDonorDims << " dimensions" << endl;}
  m_noBlocks = 0;
  io_read_iattribute1(file_id, "", "noBlocks", &m_noBlocks);

  zfsAlloc(m_noDonorCellsDir, m_noBlocks, nDim, __CALLING_FUNCTION__, 0, "m_noDonorCellsDir");
  zfsAlloc(m_noDonorPointsDir, m_noBlocks, nDim, __CALLING_FUNCTION__, 0, "m_noDonorPointsDir");
  zfsAlloc(m_noDonorPoints, m_noBlocks, __CALLING_FUNCTION__, 1, "m_noDonorPoints");
  zfsAlloc(m_noDonorCells, m_noBlocks, __CALLING_FUNCTION__, 1, "m_noDonorCells");
  m_totalNoDonorCells = 0;
  for(ZFSId blockId=0; blockId<m_noBlocks; blockId++) {
    ZFSInt ijk_max[3];
    stringstream blockName;
    blockName << "/block" << blockId <<"/";
    ZFSString blockNameStr = blockName.str();
    io_getDatasetSize(file_id , blockNameStr.c_str(), "x", m_noDonorDims , ijk_max );
    for(ZFSId j=0; j<m_noDonorDims; j++) {
      if(m_donorIsCellCentered == true) {
        m_noDonorCellsDir[blockId][j] = ijk_max[j]-1;
      } else {
        m_noDonorCellsDir[blockId][j] = ijk_max[j];
      }
      m_noDonorPointsDir[blockId][j] = ijk_max[j];
      //no of points in each block
      m_noDonorPoints[blockId] *= m_noDonorPointsDir[blockId][j];
      //no of cells in each block
      m_noDonorCells[blockId] *= m_noDonorCellsDir[blockId][j];
    }
    if(m_domainId==0) { cout << "Donor block " << blockId << " has " << m_noDonorCells[blockId] << " cells" << endl; }
    //total number of cells in all blocks together
    m_totalNoDonorCells += m_noDonorCells[blockId];
  }

  if(m_domainId==0) { cout << "totalNoDonorCells: " << m_totalNoDonorCells << endl; }

  if(m_noDonorDims==2) {
    m_totalNoDonorCells = 0;
    for(ZFSId blockId=0; blockId<m_noBlocks; blockId++) {
      m_noDonorPointsDir[blockId][2] = m_noDonorPointsDir[blockId][1];
      m_noDonorPointsDir[blockId][1] = m_noDonorPointsDir[blockId][0];
      m_noDonorPointsDir[blockId][0] = 3;
      m_noDonorCellsDir[blockId][2] = m_noDonorCellsDir[blockId][1];
      m_noDonorCellsDir[blockId][1] = m_noDonorCellsDir[blockId][0];
      m_noDonorCellsDir[blockId][0] = 2;
      m_noDonorCells[blockId] = m_noDonorCellsDir[blockId][0] *
                                m_noDonorCellsDir[blockId][1] *
                                m_noDonorCellsDir[blockId][2];
      m_noDonorPoints[blockId] = m_noDonorPointsDir[blockId][0] *
                                 m_noDonorPointsDir[blockId][1] *
                                 m_noDonorPointsDir[blockId][2];
      m_totalNoDonorCells += m_noDonorCells[blockId];
    }
  }

  zfsAlloc(m_donorBlockOffsets, m_noBlocks, __CALLING_FUNCTION__, 0, "m_donorBlockOffset");
  if(m_noBlocks>1) {
    for(ZFSId blockId=1; blockId<m_noBlocks; blockId++) {
      m_donorBlockOffsets[blockId] = m_donorBlockOffsets[blockId-1] + m_noDonorCells[blockId-1];
    }
  }

  ZFSInt start[3] = {0,0,0};
  zfsAlloc(m_donorCoordinates, nDim, m_totalNoDonorCells, "m_donorGridCoordinates", __CALLING_FUNCTION__);
  for(ZFSId blockId=0; blockId<m_noBlocks; blockId++) {
    ZFSInt offset = m_donorBlockOffsets[blockId];
    stringstream blockName;
    blockName << "/block" << blockId <<"/";
    ZFSString blockNameStr = blockName.str();

    if(m_noDonorDims==3) {
      ZFSFloatScratchSpace gridCoordinates(m_noDonorPoints[blockId], __CALLING_FUNCTION__, "gridCoordinates");
      for(ZFSId dim = 0; dim<m_noDonorDims; dim++) {
        ZFSString coordName[3] = {"x","y","z"};
        io_read_ddataset_part1d1(file_id, blockNameStr.c_str(), coordName[dim].c_str(), m_noDonorDims, start, m_noDonorPointsDir[blockId], gridCoordinates.begin());
        if(!m_donorIsCellCentered) {
          for(ZFSId cellId = 0; cellId < m_noDonorCells[blockId]; cellId++)  {
            m_donorCoordinates[dim][offset + cellId] = gridCoordinates[cellId];
          }
        } else {
          computeCellCentreCoordinates(m_noDonorPointsDir[blockId], gridCoordinates, dim, blockId);
        }
        for(ZFSId cellId = 0; cellId < m_noDonorCells[blockId]; cellId++) {
          m_donorCoordinates[dim][offset + cellId] = translation[dim] + scale[dim]*m_donorCoordinates[dim][offset + cellId];
        }
      }
    } else {
      ZFSInt gridPoints2D[2] = {m_noDonorPointsDir[blockId][1],m_noDonorPointsDir[blockId][2]};
      ZFSFloatScratchSpace gridCoordinates(m_noDonorPoints[blockId], __CALLING_FUNCTION__, "gridCoordinates");

      for(ZFSId dim = 0; dim<nDim; dim++) {
        ZFSString coordName[3] = {"x","y","z"};
        if(dim<2) {
          io_read_ddataset_part1d1(file_id, blockNameStr.c_str(), coordName[dim].c_str(), m_noDonorDims, start, gridPoints2D, gridCoordinates.begin());

          //duplicate the values in z-direction
          for(ZFSInt i=0; i<m_noDonorPointsDir[blockId][2]; i++) {
            for(ZFSInt j=0; j<m_noDonorPointsDir[blockId][1]; j++) {
              for(ZFSInt k=1; k<m_noDonorPointsDir[blockId][0]; k++) {
                ZFSId pointId2D = i+j*m_noDonorPointsDir[blockId][2];
                ZFSId pointId = i+(j+k*m_noDonorPointsDir[blockId][1])*m_noDonorPointsDir[blockId][2];
                gridCoordinates(pointId) = gridCoordinates(pointId2D);
              }
            }
          }
        } else {
          ZFSFloat minZ = F0;
          ZFSFloat maxZ = F1;
          for(ZFSInt i=0; i<m_noDonorPointsDir[blockId][2]; i++) {
            for(ZFSInt j=0; j<m_noDonorPointsDir[blockId][1]; j++) {
              for(ZFSInt k=0; k<m_noDonorPointsDir[blockId][0]; k++) {
                ZFSId pointId = i+(j+k*m_noDonorPointsDir[blockId][1])*m_noDonorPointsDir[blockId][2];
                gridCoordinates(pointId) = (maxZ-minZ)/(m_noDonorPointsDir[blockId][0]-1)*k;
              }
            }
          }
        }

        computeCellCentreCoordinates(m_noDonorPointsDir[blockId], gridCoordinates, dim, blockId);
        for(ZFSId cellId = 0; cellId < m_noDonorCells[blockId]; cellId++) {
          m_donorCoordinates[dim][cellId] = translation[dim] + scale[dim]*m_donorCoordinates[dim][offset + cellId];
        }
      }
    }
  }
  
  io_closefile(file_id);

  ZFSBool write2D3DGrid = false;
  if(write2D3DGrid) {
    const char* fileName = "rescaledGrid.hdf5";
    ZFSInt file = io_openfile("hdf5", fileName, "collective", m_zfsStrctrdComm);
    for(ZFSId blockId=0; blockId<m_noBlocks; blockId++) {
      ZFSInt offset = m_donorBlockOffsets[blockId];
      stringstream blockName;
      blockName << "/block" << blockId <<"/";
      ZFSString blockNameStr = blockName.str();
      io_create_ddataset(file, blockNameStr.c_str() , "x", 3,  m_noDonorCells);
      io_create_ddataset(file, blockNameStr.c_str() , "y", 3,  m_noDonorCells);
      io_create_ddataset(file, blockNameStr.c_str() , "z", 3,  m_noDonorCells);
      ZFSId localOffset[3] = {0,0,0};

      if(m_domainId==0) {
        io_write_ddataset_part(file, blockNameStr.c_str(), "x", nDim, m_noDonorCellsDir[blockId], localOffset, &m_donorCoordinates[0][offset + 0]);
        io_write_ddataset_part(file, blockNameStr.c_str(), "y", nDim, m_noDonorCellsDir[blockId], localOffset, &m_donorCoordinates[1][offset + 0]);
        io_write_ddataset_part(file, blockNameStr.c_str(), "z", nDim, m_noDonorCellsDir[blockId], localOffset, &m_donorCoordinates[2][offset + 0]);
      } else {
        ZFSInt size[3] = {0,0,0};
        io_write_ddataset_part(file, blockNameStr.c_str(), "x", nDim, size, localOffset, NULL);
        io_write_ddataset_part(file, blockNameStr.c_str(), "y", nDim, size, localOffset, NULL);
        io_write_ddataset_part(file, blockNameStr.c_str(), "z", nDim, size, localOffset, NULL);
      }

    }
    io_closefile(file);
  }
  if(m_domainId==0) {cout << "Reading in donor grid file... FINISHED!" << endl;}
}

/** \brief Loads a Q file
 *
 * loads a given variable (with varName) from donorFile
 *
 * \author Marian Albers, Nov 2015
 */
template <ZFSInt nDim>
void ZFSStrctrdInterpolation<nDim>::loadDonorVariable(ZFSString varName)
{
  stringstream donorFileName;
  ZFSInt donorFileId =-1;

  if(m_domainId==0) {cout << "Reading in " << varName << " from donor file..." << endl;}
  //reading in a specified donor file from the properties file

  /*! \page propertyPage1
    \section donorVars
    <code>ZFSInt ZFSStrctrdInterpolation::m_donorVars </code>\n
    default = <code> "" </code>\n \n
    Name of the donor var file to interpolate from.\n
    Possible values are:\n
    <ul>
    <li>String containing file name</li>
    </ul>
    Keywords: <i>INTERPOLATION, STRCTRD</i>
  */
  ZFSString donorFile = *(ZFSContext::getProperty("donorVars", 0, __CALLING_FUNCTION__, (ZFSString*) NULL )->asString(0));

  //open the file
  donorFileId = io_openfile("hdf5",donorFile.c_str(),"collective", m_zfsStrctrdComm);

  ZFSInt start[ 3 ] = {0,0,0};
  for(ZFSId blockId=0; blockId<m_noBlocks; blockId++) {
    stringstream blockName;
    blockName << "/block" << blockId <<"/";
    ZFSString blockNameStr = blockName.str();

    ZFSBool fieldExists = io_checkObj(donorFileId, blockNameStr.c_str(), varName.c_str());
    if(fieldExists) {
      if(m_domainId==0) {cout << "Field " << varName << " exist, loading from donorFile. " << endl;}
      ZFSInt offset = m_donorBlockOffsets[blockId];
      ZFSInt localNoDonorCellsDir[3] = {0,0,0};
      if(m_noDonorDims==3) {
        localNoDonorCellsDir[0] = m_noDonorCellsDir[blockId][0];
        localNoDonorCellsDir[1] = m_noDonorCellsDir[blockId][1];
        localNoDonorCellsDir[2] = m_noDonorCellsDir[blockId][2];
      } else {
        localNoDonorCellsDir[0] = m_noDonorCellsDir[blockId][1];
        localNoDonorCellsDir[1] = m_noDonorCellsDir[blockId][2];
      }

      io_read_ddataset_part1d1(donorFileId, blockNameStr.c_str(), varName.c_str(), m_noDonorDims, start, localNoDonorCellsDir, &m_donorVar[offset] );

      //enlarge 2D field to 3D
      if(m_noDonorDims==2) {
        for(ZFSId k=1; k<m_noDonorCellsDir[blockId][0]; k++) {
          for(ZFSId j=0; j<m_noDonorCellsDir[blockId][1]; j++) {
            for(ZFSId i=0; i<m_noDonorCellsDir[blockId][2];i++) {
              const ZFSId cellId = cellIndex(i,j,k,blockId);
              const ZFSId cellId2D = cellIndex(i,j,0,blockId);
              m_donorVar[cellId] = m_donorVar[cellId2D];
            }
          }
        }
      }
    } else {
      if(m_domainId==0) {cout << "Field " << varName << " doesn't exist, setting zero. " << endl;}
      for(ZFSId k=0; k<m_noDonorCellsDir[blockId][0]; k++) {
        for(ZFSId j=0; j<m_noDonorCellsDir[blockId][1]; j++) {
          for(ZFSId i=0; i<m_noDonorCellsDir[blockId][2];i++) {
            const ZFSId cellId = cellIndex(i,j,k,blockId);
            m_donorVar[cellId] = F0;
          }
        }
      }
    }
  }

  io_closefile(donorFileId);

  if(m_domainId==0) {cout << "Reading in " << varName << " from donor file... FINISHED!" << endl;}
}

/** \brief Computes cell centered coordinates
 *
 * Computes the cell-centered donorCoordinates from
 * node-centered grid
 *
 * \author Marian Albers, Nov 2015
 */
template <ZFSInt nDim>
void ZFSStrctrdInterpolation<nDim>::computeCellCentreCoordinates(ZFSInt* noPoints,
                                                                 ZFSFloatScratchSpace& coordinates,
                                                                 ZFSInt dim,
                                                                 ZFSInt blockId)
{
  //function to compute the coordinates at cell centre
  //do it over I, J, K loop but change to one array

  for(ZFSId k=0; k<m_noDonorCellsDir[blockId][0]; k++) {
    for(ZFSId j=0; j<m_noDonorCellsDir[blockId][1]; j++) {
      for(ZFSId i=0; i<m_noDonorCellsDir[blockId][2];i++) {
        ZFSId pointId = i+(j+k*noPoints[1])*noPoints[2];
        ZFSId IJK = pointId;
        ZFSId IP1JK= pointId+1;
        ZFSId IJP1K= pointId+noPoints[2];
        ZFSId IP1JP1K= IJP1K+1;
        ZFSId IJKP1= pointId+noPoints[2]*noPoints[1];
        ZFSId IP1JKP1= IJKP1+1;
        ZFSId IJP1KP1= pointId+noPoints[2]+noPoints[2]*noPoints[1];
        ZFSId IP1JP1KP1= IJP1KP1+1;
        ZFSId cellId = cellIndex(i,j,k,blockId);

        //average the coordinates for cell centre data
        m_donorCoordinates[dim][cellId]=F1B8*(coordinates[IJK]+
                                              coordinates[IP1JK]+
                                              coordinates[IJP1K]+
                                              coordinates[IP1JP1K]+
                                              coordinates[IJKP1]+
                                              coordinates[IP1JKP1]+
                                              coordinates[IJP1KP1]+
                                              coordinates[IP1JP1KP1]);
      }
    }
  }
}

/** \brief Nearest neighbour interpolation
 *
 * Fallback routine for when no trilinear interpolation
 * is possible. Takes values from nearest neighbour.
 *
 * \author Marian Albers, Nov 2015
 */
template <ZFSInt nDim>
inline void ZFSStrctrdInterpolation<nDim>::nearestNeighbourInterpolation(ZFSId donorId, ZFSId receiverId, ZFSFloat* receiverVar)
{
  receiverVar[receiverId] = m_donorVar[donorId];
}

template <ZFSInt nDim>
inline void ZFSStrctrdInterpolation<nDim>::nearestNeighbourInterpolation(ZFSId donorId, ZFSFloat* receiverVariables)
{
    for(ZFSInt var=0; var<m_noDonorVariables; var++) {
      receiverVariables[var] = m_donorVariables[var][donorId];
    }
}

template <ZFSInt nDim>
inline void ZFSStrctrdInterpolation<nDim>::crossProduct( ZFSFloat result[3], ZFSFloat vec1[3], ZFSFloat vec2[3])
{
  result[xsd] = vec1[ysd] * vec2[zsd] - vec1[zsd] * vec2[ysd];
  result[ysd] = vec1[zsd] * vec2[xsd] - vec1[xsd] * vec2[zsd];
  result[zsd] = vec1[xsd] * vec2[ysd] - vec1[ysd] * vec2[xsd];
}

template <ZFSInt nDim>
inline ZFSFloat ZFSStrctrdInterpolation<nDim>::scalarProduct( ZFSFloat vec1[3], ZFSFloat vec2[3])
{
  return  vec1[xsd] * vec2[xsd] + vec1[ysd] * vec2[ysd] + vec1[zsd] * vec2[zsd];
}

template <ZFSInt nDim>
inline ZFSId ZFSStrctrdInterpolation<nDim>::getCellIdfromCell(ZFSId origin,
                                                              ZFSInt incI,
                                                              ZFSInt incJ,
                                                              ZFSInt incK,
                                                              ZFSInt blockId )
{
  return origin + incI + incJ * m_noDonorCellsDir[blockId][2] + incK * m_noDonorCellsDir[blockId][2] * m_noDonorCellsDir[blockId][1];
}

template <ZFSInt nDim>
inline ZFSId ZFSStrctrdInterpolation<nDim>::cellIndex(ZFSInt i, ZFSInt j, ZFSInt k, ZFSInt blockId)
{
  const ZFSInt offset = m_donorBlockOffsets[blockId];
  return offset + i+(j+k*m_noDonorCellsDir[blockId][1])*m_noDonorCellsDir[blockId][2];
}

//returns the ijk-coordinate for one of the 4 side
//of one of the 6 tetraeders inside a hexahedron
template <ZFSInt nDim>
inline ZFSId ZFSStrctrdInterpolation<nDim>::ic(ZFSInt tetra, ZFSInt side, ZFSInt dim)
{
  return m_pyramidPoints[dim + side*3 + tetra*3*4];
}

template <ZFSInt nDim>
inline ZFSId ZFSStrctrdInterpolation<nDim>::getBlockId(ZFSInt cellId) {
  ZFSId currentBlockId = 0;
  if(m_noBlocks>1) {
    for(ZFSId blockId=0; blockId<m_noBlocks; blockId++) {
      if(cellId < m_donorBlockOffsets[blockId+1]) {
        currentBlockId = blockId;
        break;
      }
    }
  }

  return currentBlockId;
}

template<ZFSInt nDim>   //junoh
inline ZFSBool ZFSStrctrdInterpolation<nDim>:: approx (const ZFSFloat& a, const ZFSFloat& b, const ZFSFloat eps) { 
  return abs(a - b) < eps;
}


// Explicit instantiations for 2D and 3D
template class ZFSStrctrdInterpolation<2>;
template class ZFSStrctrdInterpolation<3>;
