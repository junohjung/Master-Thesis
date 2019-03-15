 
#include <cmath>
#include <iostream>
#include <fstream>

#include "zfsstrctrdblck.h"
#include "zfsstrctrdbndrycnd3d.h"
#include "zfsstrctrdblck3d.h"
#include "zfstypes.h"
#include "zfsiolib.h"
#include "zfspointbox.h"
#include "zfskdtree.h"

using namespace std;

template <ZFSBool isRans>
constexpr ZFSInt ZFSStrctrdBndryCnd3D<isRans>::m_reverseCellIdDim[18];
template <ZFSBool isRans>
constexpr ZFSInt ZFSStrctrdBndryCnd3D<isRans>::m_reverseCellIdGC[18];

/*
 *    random number generator
 *    @author Marian Albers
 *    @date 01.01.1010
 *
 */
template <ZFSBool isRans>
ZFSFloat ZFSStrctrdBndryCnd3D<isRans>::generate_rand()
{
  return rand() / double(RAND_MAX);
}

/**
 *  Generate weighted random numbers, for a distribution with
 *  a higher probability towards the higher bound pass a number < 1.0
 *  for a higher probability towards lower bound pass a number > 1.0
 *  (see tfs)
 * @author Marian Albers
 * @date 01.01.1010
 * 
 */
template <ZFSBool isRans>
ZFSFloat ZFSStrctrdBndryCnd3D<isRans>::generate_rand_weighted()
{
  return pow((rand() / double(RAND_MAX)),m_block->m_stgEddieDistribution);
}


/**
 *   Constructor of the Boundary Block of the Strctrd Block in 3D
 *   @author: Pascal Meysonnat
 *   @date: 01.01.1010
 *   
 */
template <ZFSBool isRans>
ZFSStrctrdBndryCnd3D<isRans>::ZFSStrctrdBndryCnd3D( ZFSStrctrdBlck<3>* block, ZFSId noSpecies)
  : ZFSStrctrdBndryCnd<3>( block ), m_noSpecies(noSpecies)
{
  TRACE();  

  const ZFSLong oldAllocatedBytes = allocatedBytes();

  m_block = static_cast<class ZFSStrctrdBlck3D*> (block);
  m_noStrctrdCells=m_block->m_noStrctrdCells;

  m_startCommPeriodic=m_block->m_cmnctnFlag->noNghbrDomainsNormal+m_block->m_cmnctnFlag->noNghbrDomainsChannel;
  m_endCommPeriodic=m_block->m_cmnctnFlag->noNghbrDomainsNormal+m_block->m_cmnctnFlag->noNghbrDomainsChannel+m_block->m_cmnctnFlag->noNghbrDomainsPeriodic;
  m_periodicS=m_block->m_cmnctnFlag->noNghbrDomainsPeriodicS;

  m_startCommChannel=m_block->m_cmnctnFlag->noNghbrDomainsNormal;
  m_endCommChannel=m_block->m_cmnctnFlag->noNghbrDomainsNormal+m_block->m_cmnctnFlag->noNghbrDomainsChannel;

  // aux map variables
m_bCfCpCoeff = m_block->m_bCfCpCoeff;
//m_bCf = m_block->m_bCf;
//m_bCp = m_block->m_bCp;
  m_bCl = m_block->m_bCl;
  m_bPower = m_block->m_bPower;
  m_bCd = m_block->m_bCd;
  m_bCpLineAveraging = m_block->m_bCpLineAveraging;
  m_cpAveragingDir = m_block->m_cpAveragingDir;

  ZFSId noWalls=0;
  for(ZFSUint j=0; j<m_globalStrctrdBndryMaps.size(); ++j){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_globalStrctrdBndryMaps[j]->BC)/1000.0);
    if(firstDigit==1){
      noWalls++;
    }
  }

  m_noPeriodicConnections=0;
  for(ZFSUint j=0; j<m_globalStrctrdBndryMaps.size(); ++j){
    if (m_globalStrctrdBndryMaps[j]->BC > 4400 && m_globalStrctrdBndryMaps[j]->BC < 4407){
      m_noPeriodicConnections++;
    }
  }
  m_noPeriodicConnections/=2;

  if(m_bCl){
    zfsAlloc(m_cL, nDim*noWalls, "m_cL", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_cLp, nDim*noWalls, "m_cLp", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_cLv, nDim*noWalls, "m_cLv", F0, __CALLING_FUNCTION__);
  }

  if(m_bCd){
    zfsAlloc(m_cD, nDim*noWalls, "m_cD", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_cDp, nDim*noWalls, "m_cDp", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_cDv, nDim*noWalls, "m_cDv", F0, __CALLING_FUNCTION__);
  }

  if(m_bPower && m_bCfCpCoeff){
    zfsAlloc(m_Powerv, nDim*noWalls, "m_Powerv", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_Powerp, nDim*noWalls, "m_Powerp", F0, __CALLING_FUNCTION__);
  }

  zfsAlloc(m_cArea, noWalls, "m_cArea", F0, __CALLING_FUNCTION__);

  //allocate virtual box
  if(m_block->m_stgIsActive)
    {
      zfsAlloc(m_stgVbStart, 3, "m_stgVbStart", 0.0, __CALLING_FUNCTION__);
      zfsAlloc(m_stgVbEnd, 3, "m_stgVbEnd", 0.0, __CALLING_FUNCTION__);

      //srand ( time(NULL) );
    }

  // set rotational bc properties
  setRotationalBCProperties();

  //read rescaling bc properties
  if(*m_block->m_rescalingCommGrRoot >= 0) {

    /*! \page propertyPage1
      \section rescalingBLT
      <code>ZFSInt ZFSStrctrdBlck::m_rescalingBLT </code>\n
      default = <code> 1.0 </code>\n \n
      Delta0 thickness at the inflow to be rescaled to\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>RESCALING, STRCTRD</i>
    */
      m_rescalingBLT = 5.95;
      m_rescalingBLT = *(ZFSContext::getProperty("rescalingBLT", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL)->asFloat(0));
    }

  // compute sponge coefficients for each cell
  if(m_block->m_useSponge&&m_block->m_computeSpongeFactor) {
    RECORD_TIMER_START(m_block->m_tbuildSponge);
    readAndDistributeSpongeCoordinates();
    RECORD_TIMER_STOP(m_block->m_tbuildSponge);
  }

  printAllocatedMemory( oldAllocatedBytes, "ZFSStrctrdBndryCnd3D", m_zfsStrctrdComm );
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::setRotationalBCProperties()
{
  /*! \page propertyPage1
    \section periodicRotationBC
    <code>ZFSInt ZFSStrctrdBlck::m_periodicRotationBC </code>\n
    default = <code> 1.0 </code>\n \n
    Trigger the use of the periodic rotation BC\n
    Possible values are:\n
    <ul>
    <li>true/false</li>
    </ul>
    Keywords: <i>ROTATION, BC, STRCTRD</i>
  */
  m_rotBC.hasRotationPeriodicBC=0;
  m_rotBC.hasRotationPeriodicBC=*(ZFSContext::getProperty("periodicRotationBC", m_blockId, __CALLING_FUNCTION__, &m_rotBC.hasRotationPeriodicBC)->asInt(0));
  if(m_rotBC.hasRotationPeriodicBC){
    zfsAlloc(m_rotBC.perNormals, 3, "rotation periodic axis normals", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_rotBC.prevNormal, 3, "normal buffer for two-side domains", F0, __CALLING_FUNCTION__);
    m_rotBC.initRun = 1;

    m_rotBC.prevNormal[0] = -1;
    m_rotBC.prevNormal[1] = -1;
    m_rotBC.prevNormal[2] = -1;

    /*! \page propertyPage1
      \section rotationAxisNormals
      <code>ZFSInt ZFSStrctrdBlck::m_rotationAxisNormals </code>\n
      default = <code> 0.0 </code>\n \n
      Normal of the rotation axis for the Rotation BC\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>ROTATION, BC, STRCTRD</i>
    */
    for(ZFSId i=0; i<3; ++i ) m_rotBC.perNormals[i] =*(ZFSContext::getProperty("rotationAxisNormals", m_blockId, __CALLING_FUNCTION__, &m_rotBC.perNormals[i] ,3)->asFloat(i));
    gatherPeriodic=&ZFSStrctrdBndryCnd3D<isRans>::gatherPeriodicRotation; //use periodic exchange with rotation

    //rotation matrix is needed for the periodic communication
    zfsAlloc(m_rotBC.rotationMatrix4001, 3,3, "rotMatrix4001", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_rotBC.rotationMatrix4002, 3,3, "rotMatrix4002", F0, __CALLING_FUNCTION__);
    //initialisation of the rotation Matrix
    //containing only values on the matrix diagonal for NO ratition

    m_rotBC.rotationMatrix4001[0][0]=F1;
    m_rotBC.rotationMatrix4001[1][1]=F1;
    m_rotBC.rotationMatrix4001[2][2]=F1;

    m_rotBC.rotationMatrix4002[0][0]=F1;
    m_rotBC.rotationMatrix4002[1][1]=F1;
    m_rotBC.rotationMatrix4002[2][2]=F1;
  }else{
    gatherPeriodic=&ZFSStrctrdBndryCnd3D<isRans>::gather;
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::readAndDistributeSpongeCoordinates(){

  vector<ZFSStrctrdWindowMap*> spongeInfo = m_block->m_windowInfo->m_spongeInfoMap;
  ZFSInt noSpongeInfo = spongeInfo.size(); //contains the size of the maps
  ZFSInt memSize = 0;

  //1)
  //allocate the space for all the coordinates of the sponge in Scratch!
  //memory will not be needed later ==> Scratch!

  //1.1) determine the storage size
  //determine the size and to store the whole data
  for(ZFSId i=0; i<noSpongeInfo; ++i){
    ZFSInt size=1;
    for(ZFSInt dim=0; dim<3; ++dim){
      size*=(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]+1);
    }
    memSize+=size;
  }

  //1.2) allocate the memory
  ZFSFloatScratchSpace coordMem(3*memSize, __CALLING_FUNCTION__, "spongeCoordinates" );
  ZFSFloatPointerScratchSpace spongeCoords(noSpongeInfo,__CALLING_FUNCTION__, "spongeCoordPointer");

  ZFSInt totMemSize=0;
  for(ZFSId i=0; i<noSpongeInfo; ++i){
    ZFSInt size=1;
    for(ZFSId dim=0; dim<3; ++dim){
      size*=(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]+1);
    }
    spongeCoords[i]=&coordMem[totMemSize];
    totMemSize+=3*size;
  }


  //2)
  // we do not need the corner points but the centre coordinates of the face
  // we need to allocate the memory for the faces (==cell size) for the sponge face
  // ->determine the size and allocate the memory (again only scratch)

  //2.1) calculate the number of cells to store.
  ZFSInt cellmemSize = 0;
  //determine the size and to store the whole data
  for(ZFSId i=0; i<noSpongeInfo; ++i){
    ZFSInt size=1;
    for(ZFSInt dim=0; dim<3; ++dim){
      if(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]==0) continue;
      size*=(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]);
    }
    cellmemSize+=size;
  }
  //2.2) allocate the space for all the coordinates in Scratch!
  //memory will not be needed later!
  ZFSFloatScratchSpace coordCellMem(3*cellmemSize, __CALLING_FUNCTION__, "spongeCellCoordinates" );
  ZFSFloatPointerScratchSpace spongeSurfCoords(noSpongeInfo,__CALLING_FUNCTION__, "spongeCellCoordPointer");

  ZFSInt totCellMemSize=0;
  for(ZFSId i=0; i<noSpongeInfo; ++i){
    ZFSInt size=1;
    for(ZFSId dim=0; dim<3; ++dim){
      if(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]==0) continue;
      size*=(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]);
    }
    spongeSurfCoords[i]=&coordCellMem[totCellMemSize];
    totCellMemSize+=3*size;
  }

  //3) read in the coordinates of the grid points
  //open file for reading the grid data
  ZFSString gridFileName = m_block->m_gridInputFileName;

  //file id to access the file
  ZFSId file_id = -1;
  //unique identifier needed to associate grid and solution in case of restart
  //const char* aUID= new char[18];

  zfs_log << "Loading and broadcasting all sponge windows..." << endl;
  //open file and read number of blocks and uid
  file_id = io_openfile("hdf5", gridFileName.c_str(), "collective", m_zfsStrctrdComm);

  //read the data in and distribute ist

  //the split of reading and distributing is done on purpose!
  //this is because we then know what the library is doing
  //and not what the io_library such as hdf or netcdf should do
  //but does not
  //a good library should handle that directly! TEST IT
  memSize =1;
  //read in the data if  processor zero else read nothing!!!
  if(m_block->domainId()==0){
    for(ZFSId i=0; i<noSpongeInfo; ++i){
      ZFSId offset[3]={0,0,0};
      ZFSId size[3]={0,0,0};
      memSize = 1;
      for(ZFSId dim=2; dim>=0; --dim){
        size[dim]=(spongeInfo[i]->end1[2-dim]-spongeInfo[i]->start1[2-dim]+1);
        memSize*=size[dim];
        offset[dim]= spongeInfo[i]->start1[2-dim];
      }
      //read in the data if  processor zero else read nothing!!!
      //determine the Block name
      ZFSString bName="block";
      stringstream number;
      number << spongeInfo[i]->Id1;
      bName += number.str();

      io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, &spongeCoords[i][0]);
      io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, &spongeCoords[i][memSize]);
      io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, &spongeCoords[i][memSize*2]);
    }
  } else {
    for(ZFSId i=0; i<noSpongeInfo; ++i){
      ZFSId offset[3]={0,0,0};
      ZFSId size[3]={0,0,0};
      ZFSString bName="block";
      stringstream number;
      number << spongeInfo[i]->Id1;
      bName += number.str();
      io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, NULL);
      io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, NULL);
      io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, NULL);
    }
  }

  io_closefile(file_id);

  //now broadcast the information to everyone!!!
  MPI_Bcast(&spongeCoords[0][0], totMemSize, MPI_DOUBLE, 0, m_zfsStrctrdComm);
  zfs_log << "Loading and broadcasting all sponge windows... SUCCESSFUL!" << endl;
  zfs_log << "Computing sponge surface center coordinates..." << endl;

  //4) computing the coordinates of surface center from corner points;

  for(ZFSId ii=0; ii<noSpongeInfo; ++ii){
    ZFSInt label,size1,size2,count=0;
    for(label=0; label<3; label++){ //3== dimensions
      if(spongeInfo[ii]->end1[label]-spongeInfo[ii]->start1[label]==0) break;
    }
    switch(label){
    case 0:{
      size1=spongeInfo[ii]->end1[1]-spongeInfo[ii]->start1[1]+1;
      size2=spongeInfo[ii]->end1[2]-spongeInfo[ii]->start1[2]+1;
      for(ZFSInt j=0;j<size2-1;j++){
        for(ZFSInt i=0;i<size1-1;i++){
          ZFSInt IJ=i+j*size1;
          ZFSInt IPJ=i+(j+1)*size1;
          ZFSInt IJP=i+1+j*size1;
          ZFSInt IPJP=i+1+(j+1)*size1;
          spongeSurfCoords[ii][count]=0.25*(spongeCoords[ii][IJ]+spongeCoords[ii][IPJ]+spongeCoords[ii][IJP]+spongeCoords[ii][IPJP]);
          spongeSurfCoords[ii][count+(size1-1)*(size2-1)]=0.25*(spongeCoords[ii][IJ+size1*size2]+spongeCoords[ii][IPJ+size1*size2]+spongeCoords[ii][IJP+size1*size2]+spongeCoords[ii][IPJP+size1*size2]);
          spongeSurfCoords[ii][count+2*(size1-1)*(size2-1)]=0.25*(spongeCoords[ii][IJ+2*size1*size2]+spongeCoords[ii][IPJ+2*size1*size2]+spongeCoords[ii][IJP+2*size1*size2]+spongeCoords[ii][IPJP+2*size1*size2]);
          count++;
        }
      }
      break;
    }
    case 1:{
      size1=spongeInfo[ii]->end1[0]-spongeInfo[ii]->start1[0]+1;
      size2=spongeInfo[ii]->end1[2]-spongeInfo[ii]->start1[2]+1;
      for(ZFSInt j=0;j<size2-1;j++){
        for(ZFSInt i=0;i<size1-1;i++){
          ZFSInt IJ=i+j*size1;
          ZFSInt IPJ=i+(j+1)*size1;
          ZFSInt IJP=i+1+j*size1;
          ZFSInt IPJP=i+1+(j+1)*size1;
          spongeSurfCoords[ii][count]=0.25*(spongeCoords[ii][IJ]+spongeCoords[ii][IPJ]+spongeCoords[ii][IJP]+spongeCoords[ii][IPJP]);
          spongeSurfCoords[ii][count+(size1-1)*(size2-1)]=0.25*(spongeCoords[ii][IJ+size1*size2]+spongeCoords[ii][IPJ+size1*size2]+spongeCoords[ii][IJP+size1*size2]+spongeCoords[ii][IPJP+size1*size2]);
          spongeSurfCoords[ii][count+2*(size1-1)*(size2-1)]=0.25*(spongeCoords[ii][IJ+2*size1*size2]+spongeCoords[ii][IPJ+2*size1*size2]+spongeCoords[ii][IJP+2*size1*size2]+spongeCoords[ii][IPJP+2*size1*size2]);
          count++;
        }
      }
      break;
    }
    case 2:{
      size1=spongeInfo[ii]->end1[0]-spongeInfo[ii]->start1[0]+1;
      size2=spongeInfo[ii]->end1[1]-spongeInfo[ii]->start1[1]+1;
      for(ZFSInt j=0;j<size2-1;j++){
        for(ZFSInt i=0;i<size1-1;i++){
          ZFSInt IJ=i+j*size1;
          ZFSInt IPJ=i+(j+1)*size1;
          ZFSInt IJP=i+1+j*size1;
          ZFSInt IPJP=i+1+(j+1)*size1;
          spongeSurfCoords[ii][count]=0.25*(spongeCoords[ii][IJ]+spongeCoords[ii][IPJ]+spongeCoords[ii][IJP]+spongeCoords[ii][IPJP]);
          spongeSurfCoords[ii][count+(size1-1)*(size2-1)]=0.25*(spongeCoords[ii][IJ+size1*size2]+spongeCoords[ii][IPJ+size1*size2]+spongeCoords[ii][IJP+size1*size2]+spongeCoords[ii][IPJP+size1*size2]);
          spongeSurfCoords[ii][count+2*(size1-1)*(size2-1)]=0.25*(spongeCoords[ii][IJ+2*size1*size2]+spongeCoords[ii][IPJ+2*size1*size2]+spongeCoords[ii][IJP+2*size1*size2]+spongeCoords[ii][IPJP+2*size1*size2]);
          count++;
        }
      }
      break;
    }
    default: zfsTerm(1,__CALLING_FUNCTION__, "sponge direction is messed up");
    }
  }

  zfs_log << "Computing sponge surface center coordinates... SUCCESSFUL!" << endl;
  zfs_log << "Determining shortest distance and sponge factor for each cell..." << endl;

  //cout << "seraching for the nearest points(building sigma sponge)" << endl;
  //build a k-d-tree for a quick search:
  //1) rearrange the coordinates into points;
  const ZFSInt spongeTotalWorkload = noSpongeInfo*m_noStrctrdCells;
  const ZFSInt spongeWorkloadPercentage = (ZFSInt)(spongeTotalWorkload/10);

  for(ZFSId i=0; i<noSpongeInfo; ++i){
    ZFSInt noPoints=1;
    for(ZFSId dim=0; dim<3; ++dim){
      if(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]==0)continue;
      noPoints*=(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]);
    }
    //build up the points (surface centres)
    vector < Point<3> > pts;
    for(ZFSId j=0; j<noPoints; ++j){
      Point<3> a(spongeSurfCoords[i][j],spongeSurfCoords[i][j+noPoints] ,spongeSurfCoords[i][j+2*noPoints]);
      pts.push_back(a);
    }

    //build up the tree
    KDtree<3> tree(pts);
    ZFSFloat distance = -1.0;
    //go through all the cells an determine the closest distance
    ZFSFloat spongeThickness= spongeInfo[i]->spongeThickness;
    for(ZFSInt id=0; id<m_noStrctrdCells; ++id){
      if(m_block->domainId()==0) {
        if((i*m_noStrctrdCells+id)%spongeWorkloadPercentage==0) {
          ZFSInt progress = ceil((ZFSFloat)(i*m_noStrctrdCells+id)/(ZFSFloat)spongeTotalWorkload*100.0);
          cout << "Sponge computation - " << progress << " percent" << endl;
        }
      }
      distance=-1.1111111111111111; //to check
      Point<3> pt(m_cells->coordinates[0][id],m_cells->coordinates[1][id], m_cells->coordinates[2][id]);
      (void) tree.nearest(pt, distance);
      if(distance<=spongeThickness){
        ZFSFloat spongeFactor = spongeInfo[i]->sigma*pow((spongeThickness - distance)/spongeThickness, spongeInfo[i]->beta);
        m_cells->fq[FQ->SPONGE_FACTOR][id] = zfsMAX(m_cells->fq[FQ->SPONGE_FACTOR][id], spongeFactor);
      }
    }
  }

  zfs_log << "Determining shortest distance and sponge factor for each cell... SUCCESSFUL!" << endl;
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::updateSpongeLayer(){
  const ZFSFloat gammaMinusOne = m_block->m_gamma - 1.0;
  switch(m_spongeLayerType){
  case 1: {
    ZFSFloat deltaRhoE =F0, deltaRho=F0;
    for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++){
      for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++){
        for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers;i++){
          // compute the forcing term
          // for the pressure or engery and the velocity
          ZFSId cellId = cellIndex(i, j, k);

          const ZFSFloat rhoE = m_cells->pvariables[PV->P][cellId]/gammaMinusOne+F1B2*m_cells->pvariables[PV->RHO][cellId]*(POW2(m_cells->pvariables[PV->U][cellId])+
                                                                                                                            POW2(m_cells->pvariables[PV->V][cellId])+
                                                                                                                            POW2(m_cells->pvariables[PV->W][cellId]));

          deltaRhoE= rhoE-CV->rhoEInfinity;
          deltaRho = m_cells->pvariables[ PV->RHO ][ cellId ] - (CV->rhoInfinity * m_targetDensityFactor);

          m_cells->rightHandSide[ CV->RHO_E ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRhoE*m_cells->cellJac[cellId];//deltaP * m_cells->cellJac[cellId];
          m_cells->rightHandSide[ CV->RHO ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
        }
      }
    }
    break;
  }
  case 2: {
    // damp to values in the FQ field (set at startup, from RANS etc.)
    ZFSFloat deltaRhoE =F0, deltaRho=F0;
    for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++){
      for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++){
        for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers;i++){
          ZFSId cellId = cellIndex(i, j, k);

          const ZFSFloat rhoE = m_cells->pvariables[PV->P][cellId]/gammaMinusOne+F1B2*m_cells->pvariables[PV->RHO][cellId]*(POW2(m_cells->pvariables[PV->U][cellId])+
                                                                                                                            POW2(m_cells->pvariables[PV->V][cellId])+
                                                                                                                            POW2(m_cells->pvariables[PV->W][cellId]));

          deltaRhoE= rhoE-m_cells->fq[FQ->SPONGE_RHO_E][cellId];
          deltaRho = m_cells->pvariables[ PV->RHO ][ cellId ] - (m_cells->fq[FQ->SPONGE_RHO][cellId] * m_targetDensityFactor);

          m_cells->rightHandSide[ CV->RHO_E ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRhoE*m_cells->cellJac[cellId];
          m_cells->rightHandSide[ CV->RHO ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
        }
      }
    }
    break;
  }
  case 3: {
    //damp to rhoInfinity and pInfinity
    const ZFSFloat FgammaMinusOne = F1/(m_block->m_gamma - F1);
    for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++){
      for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++){
        for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers;i++){
          // compute the forcing term
          // for the pressure or engery and the velocity
          const ZFSId cellId = cellIndex(i, j, k);
          const ZFSFloat deltaP= (m_cells->pvariables[PV->P][cellId] - PV->PInfinity) * FgammaMinusOne;
          const ZFSFloat deltaRho = m_cells->pvariables[PV->RHO][cellId] - (CV->rhoInfinity * m_targetDensityFactor);

          m_cells->rightHandSide[ CV->RHO_E ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaP*m_cells->cellJac[cellId];
          m_cells->rightHandSide[ CV->RHO ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
        }
      }
    }
    break;
  }
  case 4: {
    // damp to rho from FQ and pInfinity
    const ZFSFloat FgammaMinusOne = F1/(m_block->m_gamma - F1);
    for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++){
      for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++){
        for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers;i++){
          const ZFSId cellId = cellIndex(i, j, k);
          const ZFSFloat deltaP = (m_cells->pvariables[PV->P][cellId] - PV->PInfinity) * FgammaMinusOne;
          const ZFSFloat deltaRho = m_cells->pvariables[ PV->RHO][cellId] - (m_cells->fq[FQ->SPONGE_RHO][cellId] * m_targetDensityFactor);

          m_cells->rightHandSide[ CV->RHO_E ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaP*m_cells->cellJac[cellId];
          m_cells->rightHandSide[ CV->RHO ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
        }
      }
    }
    break;
  }
  default: zfsTerm(1,__CALLING_FUNCTION__, "Sponge type doesn't exist");
  }
}

//>marian: new function to compute distance to nearest wall for each cell
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::computeWallDistances(){

  vector<ZFSStrctrdWindowMap*> wallDistInfo = m_block->m_windowInfo->m_wallDistInfoMap;
  ZFSInt noWallDistInfo = wallDistInfo.size(); //contains the size of the maps
  ZFSInt memSize = 0;

  //initialize array with high numbers
  for(ZFSInt id=0; id<m_noStrctrdCells; ++id){
    m_cells->fq[FQ->WALLDISTANCE][id] = 99999;
  }

  //1)
  //memory will not be needed later ==> Scratch!

  //1.1) determine the storage size
  //determine the size and to store the whole data
  for(ZFSId i=0; i<noWallDistInfo; ++i){
    ZFSInt size=1;
    for(ZFSInt dim=0; dim<3; ++dim){
      size*=(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]+1);
    }
    memSize+=size;
  }

  //1.2) allocate the memory
  ZFSFloatScratchSpace coordMem(3*memSize, __CALLING_FUNCTION__, "wallCoordinates" );
  ZFSFloatPointerScratchSpace wallDistCoords(noWallDistInfo,__CALLING_FUNCTION__, "wallCoordsPointer");

  ZFSInt totMemSize=0;
  for(ZFSId i=0; i<noWallDistInfo; ++i){
    ZFSInt size=1;
    for(ZFSId dim=0; dim<3; ++dim){
      size*=(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]+1);
    }
    wallDistCoords[i]=&coordMem[totMemSize];
    totMemSize+=3*size;
  }


  //2)
  // we do not need the corner points but the centre coordinates of the face
  // we need to allocate the memory for the faces (==cell size) for the sponge face
  // ->determine the size and allocate the memory (again only scratch)

  //2.1) calculate the number of cells to store.
  ZFSInt cellmemSize = 0;
  //determine the size and to store the whole data
  for(ZFSId i=0; i<noWallDistInfo; ++i){
    ZFSInt size=1;
    for(ZFSInt dim=0; dim<3; ++dim){
      if(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]==0) continue;
      size*=(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]);
    }
    cellmemSize+=size;
  }
  //2.2) allocate the space for all the coordinates in Scratch!
  //memory will not be needed later!
  ZFSFloatScratchSpace coordCellMem(3*cellmemSize, __CALLING_FUNCTION__, "wallDistCellCoordinates" );
  ZFSFloatPointerScratchSpace wallDistSurfCoords(noWallDistInfo,__CALLING_FUNCTION__, "wallDistCellCoordPointer");

  ZFSInt totCellMemSize=0;
  for(ZFSId i=0; i<noWallDistInfo; ++i){
    ZFSInt size=1;
    for(ZFSId dim=0; dim<3; ++dim){
      if(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]==0) continue;
      size*=(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]);
    }
    wallDistSurfCoords[i]=&coordCellMem[totCellMemSize];
    totCellMemSize+=3*size;
  }




  //3) read in the coordinates of the grid points
  //open file for reading the grid data
  ZFSString gridFileName = m_block->m_gridInputFileName;
  ZFSId file_id = -1;
  //open file and read number of blocks and uid

  //read the data in and distribute ist

  //the split of reading and distributing is done on purpose!
  //this is because we then know what the library is doing
  //and not what the io_library such as hdf or netcdf should do
  //but does not
  //a good library should handle that directly! TEST IT
  memSize =1;
  //read in the data if  processor zero else read nothing!!!
  if(m_block->domainId()==0){
    file_id = io_openfile("hdf5", gridFileName.c_str(), "collective",MPI_COMM_SELF); //junoh cmcter changed
    for(ZFSId i=0; i<noWallDistInfo; ++i){
      ZFSId offset[3]={0,0,0};
      ZFSId size[3]={0,0,0};
      memSize = 1;
      for(ZFSId dim=2; dim>=0; --dim){
        size[dim]=(wallDistInfo[i]->end1[2-dim]-wallDistInfo[i]->start1[2-dim]+1);
        memSize*=size[dim];
        offset[dim]= wallDistInfo[i]->start1[2-dim];
      }
      //read in the data if  processor zero else read nothing!!!
      //determine the Block name
      ZFSString bName="block";
      stringstream number;
      number << wallDistInfo[i]->Id1;
      bName += number.str();

      io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, &wallDistCoords[i][0]);
      io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, &wallDistCoords[i][memSize]);
      io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, &wallDistCoords[i][memSize*2]);
    }
    io_closefile(file_id); //junoh
  }//  else {
  //   for(ZFSId i=0; i<noWallDistInfo; ++i){
  //     ZFSId offset[3]={0,0,0};
  //     ZFSId size[3]={0,0,0};
  //     ZFSString bName="block";
  //     stringstream number;
  //     number << wallDistInfo[i]->Id1;
  //     bName += number.str();
  //     io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, NULL);
  //     io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, NULL);
  //     io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, NULL);
  //   }
  // }



  //cout << "broadcasting wall-bc information" << endl;
  //now broadcast the information to everyone!!!
  MPI_Bcast(&wallDistCoords[0][0], totMemSize, MPI_DOUBLE, 0, m_block->m_zfsStrctrdComm);
  //cout << "broadcast end" << endl;

  if(!m_block->m_rans) {    //junoh
    return;
  }


  //4) computing the coordinates of surface center from corner points;

  for(ZFSId ii=0; ii<noWallDistInfo; ++ii){
    ZFSInt label,size1,size2,count=0;
    for(label=0; label<3; label++){ //3== dimensions
      if(wallDistInfo[ii]->end1[label]-wallDistInfo[ii]->start1[label]==0) break;
    }
    switch(label){
    case 0:{
      size1=wallDistInfo[ii]->end1[1]-wallDistInfo[ii]->start1[1]+1;
      size2=wallDistInfo[ii]->end1[2]-wallDistInfo[ii]->start1[2]+1;
      for(ZFSInt j=0;j<size2-1;j++){
        for(ZFSInt i=0;i<size1-1;i++){
          ZFSInt IJ=i+j*size1;
          ZFSInt IPJ=i+(j+1)*size1;
          ZFSInt IJP=i+1+j*size1;
          ZFSInt IPJP=i+1+(j+1)*size1;
          wallDistSurfCoords[ii][count]=0.25*(wallDistCoords[ii][IJ]+wallDistCoords[ii][IPJ]+wallDistCoords[ii][IJP]+wallDistCoords[ii][IPJP]);
          wallDistSurfCoords[ii][count+(size1-1)*(size2-1)]=0.25*(wallDistCoords[ii][IJ+size1*size2]+wallDistCoords[ii][IPJ+size1*size2]+wallDistCoords[ii][IJP+size1*size2]+wallDistCoords[ii][IPJP+size1*size2]);
          wallDistSurfCoords[ii][count+2*(size1-1)*(size2-1)]=0.25*(wallDistCoords[ii][IJ+2*size1*size2]+wallDistCoords[ii][IPJ+2*size1*size2]+wallDistCoords[ii][IJP+2*size1*size2]+wallDistCoords[ii][IPJP+2*size1*size2]);
          count++;
        }
      }
      break;
    }
    case 1:{
      size1=wallDistInfo[ii]->end1[0]-wallDistInfo[ii]->start1[0]+1;
      size2=wallDistInfo[ii]->end1[2]-wallDistInfo[ii]->start1[2]+1;
      for(ZFSInt j=0;j<size2-1;j++){
        for(ZFSInt i=0;i<size1-1;i++){
          ZFSInt IJ=i+j*size1;
          ZFSInt IPJ=i+(j+1)*size1;
          ZFSInt IJP=i+1+j*size1;
          ZFSInt IPJP=i+1+(j+1)*size1;
          wallDistSurfCoords[ii][count]=0.25*(wallDistCoords[ii][IJ]+wallDistCoords[ii][IPJ]+wallDistCoords[ii][IJP]+wallDistCoords[ii][IPJP]);
          wallDistSurfCoords[ii][count+(size1-1)*(size2-1)]=0.25*(wallDistCoords[ii][IJ+size1*size2]+wallDistCoords[ii][IPJ+size1*size2]+wallDistCoords[ii][IJP+size1*size2]+wallDistCoords[ii][IPJP+size1*size2]);
          wallDistSurfCoords[ii][count+2*(size1-1)*(size2-1)]=0.25*(wallDistCoords[ii][IJ+2*size1*size2]+wallDistCoords[ii][IPJ+2*size1*size2]+wallDistCoords[ii][IJP+2*size1*size2]+wallDistCoords[ii][IPJP+2*size1*size2]);
          count++;
        }
      }
      break;
    }
    case 2:{
      size1=wallDistInfo[ii]->end1[0]-wallDistInfo[ii]->start1[0]+1;
      size2=wallDistInfo[ii]->end1[1]-wallDistInfo[ii]->start1[1]+1;
      for(ZFSInt j=0;j<size2-1;j++){
        for(ZFSInt i=0;i<size1-1;i++){
          ZFSInt IJ=i+j*size1;
          ZFSInt IPJ=i+(j+1)*size1;
          ZFSInt IJP=i+1+j*size1;
          ZFSInt IPJP=i+1+(j+1)*size1;
          wallDistSurfCoords[ii][count]=0.25*(wallDistCoords[ii][IJ]+wallDistCoords[ii][IPJ]+wallDistCoords[ii][IJP]+wallDistCoords[ii][IPJP]);
          wallDistSurfCoords[ii][count+(size1-1)*(size2-1)]=0.25*(wallDistCoords[ii][IJ+size1*size2]+wallDistCoords[ii][IPJ+size1*size2]+wallDistCoords[ii][IJP+size1*size2]+wallDistCoords[ii][IPJP+size1*size2]);
          wallDistSurfCoords[ii][count+2*(size1-1)*(size2-1)]=0.25*(wallDistCoords[ii][IJ+2*size1*size2]+wallDistCoords[ii][IPJ+2*size1*size2]+wallDistCoords[ii][IJP+2*size1*size2]+wallDistCoords[ii][IPJP+2*size1*size2]);
          count++;
        }
      }
      break;
    }
    default: zfsTerm(1,__CALLING_FUNCTION__, "wall direction is messed up");
    }
  }

  //now everyone can determine the shortest distance

  zfs_log << "wall distance computation: searching for the nearest wall" << endl;
  //cout << "seraching for the nearest" << endl;
  //build a k-d-tree for a quick search:
  //1) rearragne the coordinates into points;
  for(ZFSId i=0; i<noWallDistInfo; ++i){
    ZFSInt noPoints=1;
    for(ZFSId dim=0; dim<3; ++dim){
      if(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]==0)continue;
      noPoints*=(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]);
    }
    //build up the points (surface centres)
    vector < Point<3> > pts;
    for(ZFSId j=0; j<noPoints; ++j){
      Point<3> a(wallDistSurfCoords[i][j],wallDistSurfCoords[i][j+noPoints] ,wallDistSurfCoords[i][j+2*noPoints] );
      pts.push_back(a);
    }
    //build up the tree
    KDtree<3> tree(pts);
    ZFSFloat distance = -1.0;

    //go through all the cells an determine the closest distance
    for(ZFSInt id=0; id<m_noStrctrdCells; ++id){
      distance=-1.1111111111111111; //to check
      Point<3> pt(m_cells->coordinates[0][id],m_cells->coordinates[1][id], m_cells->coordinates[2][id]);
      (void) tree.nearest(pt, distance);

      //take minimum because another bc1000 might be further away than the current but would overwrite the actually closer distance
      m_cells->fq[FQ->WALLDISTANCE][id] = zfsMIN(m_cells->fq[FQ->WALLDISTANCE][id], distance);
    }
  }

  //correct the wall distance for the effective solution
  if(m_block->m_bc2601IsActive) {
    cout << "Correcting wall distance with gammEpsilon: " << m_block->m_bc2601GammaEpsilon << endl;
    for(ZFSInt id=0; id<m_noStrctrdCells; ++id){
      m_cells->fq[FQ->WALLDISTANCE][id] += m_block->m_bc2601GammaEpsilon;
    }
  }

  zfs_log << "Wall Distance Computation SUCESSFUL: Saved minimum distance to next wall for all cells " << endl;
}
//<marian

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::computeAuxDataRoot()
{
  if(m_bCfCpCoeff && m_bPower){
    computeFrictionPressureCoef<true>();
  }else{
    computeFrictionPressureCoef<false>();
  }
  if(m_bCl) computeLiftCoefRoot();
  if(m_bCd) computeDragCoefRoot();
  if(m_bPower) computePowerCoefRoot();
}

template <ZFSBool isRans>
ZFSStrctrdBndryCnd3D<isRans>::~ZFSStrctrdBndryCnd3D()
{

}

template <ZFSBool isRans>
inline ZFSId ZFSStrctrdBndryCnd3D<isRans>::cellIndex(ZFSInt i, ZFSInt j, ZFSInt k)
{
  return i+(j+k*m_nCells[1])*m_nCells[2];
}

template <ZFSBool isRans>
inline ZFSId ZFSStrctrdBndryCnd3D<isRans>::cellIndexBC(ZFSInt i, ZFSInt j, ZFSInt k)
{
  return i + (j+k*m_block->m_stgBoxSize[1])*m_block->m_stgBoxSize[2];
}

//function to correct the index values in the map for the different boundary conditions
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::correctBndryCndIndices()
{
  //cout << "in correctBndryCndIndices " << endl;
  //in correcting cell Information
  for(ZFSId bcId=0; bcId < m_noBndryCndIds; bcId++)
    {
      (this->*initBndryCndHandler[bcId])(bcId);
    }
}

template <ZFSBool isRans>
inline ZFSFloat ZFSStrctrdBndryCnd3D<isRans>::dist(ZFSFloat* a, ZFSFloat* b)
{
  ZFSFloat dist1=F0;
  for(ZFSId dim=0; dim<3; dim++)
    {
      dist1+=POW2(a[dim*m_noStrctrdCells]-b[dim*m_noStrctrdCells]);
    }
  return sqrt(dist1);
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc1000(ZFSId bcId)
{
  (void) bcId;
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc1003(ZFSId bcId)
{
  (void) bcId;

  /*! \page propertyPage1
    \section isothermalWallTemperature
    <code>ZFSInt ZFSStrctrdBlck::m_isothermalWallTemperature </code>\n
    default = <code> 1.0 </code>\n \n
    Isothermal wall temperature as a factor of T8\n
    Possible values are:\n
    <ul>
    <li>Float > 0.0</li>
    </ul>
    Keywords: <i>ISOTHERMAL, WALL, BC, STRCTRD</i>
  */
  m_isothermalWallTemperature = F1;
  if(ZFSContext::propertyExists("isothermalWallTemperature", m_blockId )){
    m_isothermalWallTemperature = *(ZFSContext::getProperty("isothermalWallTemperature", m_blockId, __CALLING_FUNCTION__,&m_isothermalWallTemperature)->asFloat(0));
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc1004(ZFSId bcId)
{// moving adiabatic wall
  m_block->initializeRungeKutta();
  (void) bcId;
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc1006(ZFSId bcId)
{// moving adiabatic wall
  m_block->initializeRungeKutta();
  (void) bcId;

  m_isothermalWallTemperature = F1;
  if(ZFSContext::propertyExists("isothermalWallTemperature", m_blockId )){
    m_isothermalWallTemperature = *(ZFSContext::getProperty("isothermalWallTemperature", m_blockId, __CALLING_FUNCTION__,&m_isothermalWallTemperature)->asFloat(0));
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc1007(ZFSId bcId)
{// oscillating non-moving wall
  m_block->initializeRungeKutta();
  (void) bcId;

  m_block->m_waveBeginTransition = 0.0;
  m_block->m_waveEndTransition = 0.0;
  m_block->m_waveAmplitude = 0.0;
  m_block->m_waveAmplitudePlus = 0.0;
  m_block->m_waveTimePlus = 0.0;
  m_block->m_waveTime = 0.0;
  /*
    Thomas Luerkens
    August 2016
    reads in wave parameters from property file
  */
  //time needs to be constant for traveling wave
  m_block->m_constantTimeStep = true;
  m_block->m_waveAmplitudePlus = *(ZFSContext::getProperty("waveAmplitudePlus", m_blockId, __CALLING_FUNCTION__, &m_block->m_waveAmplitudePlus)->asFloat(0));
  m_block->m_waveTimePlus = *(ZFSContext::getProperty("waveTimePlus", m_blockId, __CALLING_FUNCTION__, &m_block->m_waveTimePlus)->asFloat(0));
  m_block->m_waveBeginTransition = *(ZFSContext::getProperty("waveBeginTransition", m_blockId, __CALLING_FUNCTION__, &m_block->m_waveBeginTransition)->asFloat(0));
  m_block->m_waveEndTransition = *(ZFSContext::getProperty("waveEndTransition", m_blockId, __CALLING_FUNCTION__, &m_block->m_waveEndTransition)->asFloat(0));

  //compute Wave parameters
  const ZFSFloat cf = 0.024*pow(m_block->m_Re, -F1B4);
  const ZFSFloat uTau = sqrt(cf*F1B2)*PV->UInfinity;
  const ZFSFloat mu8 = zfsSUTHERLANDLAW(PV->TInfinity);
  m_block->m_waveAmplitude = m_block->m_waveAmplitudePlus*uTau;
  m_block->m_waveTime = m_block->m_waveTimePlus*mu8/(POW2(uTau)*m_block->m_Re0);

  cout << "Oscillation speed amplitude: " << m_block->m_waveAmplitude << " time: " << m_block->m_waveTime << endl;
}


/* Initialize with standard pressure extrapolation
 * at inflow or prescribe p_inf at outflow
 *
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc2004(ZFSId bcId)
{
  //call simple in/outflow bc to initialize ghost-cells
  //because bc2024 uses values in the ghost-cells
  bc2003(bcId);
}

/* Characteristic boundary condition supersonic after/with shock
 * Loosen/Meysonnat 2017
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc2009(ZFSId bcId)
{
  (void) bcId;
  m_sigma = *(ZFSContext::getProperty("shockAngle", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL)->asFloat(0));
  m_sigma = (m_sigma/180.0)*PI;
}


/** Rescaling inflow
 *
 *  put values from field to gc to avoid nan, only matters
 *  for first Runge-Kutta Step at first timeStep
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc2500(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  for(ZFSId var=0; var < PV->noVariables; var++) {
    for(ZFSId i=start[0]; i<end[0]; i++) {
      for(ZFSId j=start[1]; j<end[1]; j++) {
        for(ZFSId k=start[2]; k< end[2]; k++) {
          ZFSId cellId = cellIndex(i,j,k);
          ZFSId cellIdAdj = cellIndex(m_noGhostLayers, j,k);
          m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdAdj];
        }
      }
    }
  }
}


/** Prescribe profile BC
 *  In case of the initialStartup2600 flag is set
 *  (no matter if it was a restart or a initial startup)
 *  we load the values from the field to the GC and
 *  save them in the restart field
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc2600(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  m_block->m_bc2600 = true;
  const ZFSInt maxNoVariables=6; //junoh

  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
      {
        if(m_block->m_bc2600InitialStartup) {
          //First copy values from the field into the ghostcells
          for(ZFSInt i = start[0]; i<end[0]; i++) {
            for(ZFSInt j = start[1]; j<end[1]; j++) {
              for(ZFSInt k = start[2]; k<end[2]; k++) {
                const ZFSId cellId = cellIndex(i,j,k);
                const ZFSId cellIdAdj = cellIndex(m_noGhostLayers,j,k);

                for(ZFSId var=0; var<maxNoVariables; var++) { //junoh
                  m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdAdj];
                }
              }
            }
          }

          //Fix diagonal cells at end of domain
          if(m_block->m_bc2600noOffsetCells[2] == 0 && 
             m_block->m_bc2600noOffsetCells[1] + m_block->m_bc2600noActiveCells[1] == m_block->m_nInputBlockCells[1]) {
            for(ZFSInt i = 0; i<m_block->m_bc2600noCells[2]; i++) {
              for(ZFSInt k = 0; k<m_block->m_bc2600noCells[0]; k++) {
                const ZFSId cellIdA2 = cellIndex(i,m_noGhostLayers+m_block->m_bc2600noActiveCells[1]-2,k);
                const ZFSId cellIdA1 = cellIndex(i,m_noGhostLayers+m_block->m_bc2600noActiveCells[1]-1,k);
                const ZFSId cellIdG1 = cellIndex(i,m_noGhostLayers+m_block->m_bc2600noActiveCells[1],  k);
                for(ZFSId var=0; var<maxNoVariables; var++) { //junoh
                  const ZFSFloat distA1A2 = sqrt(POW2(m_cells->coordinates[0][cellIdA1]-m_cells->coordinates[0][cellIdA2])+
                                                 POW2(m_cells->coordinates[1][cellIdA1]-m_cells->coordinates[1][cellIdA2])+
                                                 POW2(m_cells->coordinates[2][cellIdA1]-m_cells->coordinates[2][cellIdA2]));
                  const ZFSFloat slope = (m_cells->pvariables[var][cellIdA1]-m_cells->pvariables[var][cellIdA2])/distA1A2;
                  const ZFSFloat distG1A1 = sqrt(POW2(m_cells->coordinates[0][cellIdG1]-m_cells->coordinates[0][cellIdA1])+
                                                 POW2(m_cells->coordinates[1][cellIdG1]-m_cells->coordinates[1][cellIdA1])+
                                                 POW2(m_cells->coordinates[2][cellIdG1]-m_cells->coordinates[2][cellIdA1]));
                  m_cells->pvariables[var][cellIdG1] = m_cells->pvariables[var][cellIdA1] + distG1A1*slope;
                }
              }
            }
          }
        }

        //Then copy the values from the ghostcells into restart field
        for(ZFSInt i = start[0]; i<end[0]; i++) {
          for(ZFSInt j = start[1]; j<end[1]; j++) {
            for(ZFSInt k = start[2]; k<end[2]; k++) {
              const ZFSId cellId = cellIndex(i,j,k);
              const ZFSId cellIdBc = i + (j + k*m_block->m_bc2600noCells[1])*m_block->m_bc2600noCells[2];
              for(ZFSId var=0; var<maxNoVariables; var++) { //junoh
                m_block->m_bc2600Variables[var][cellIdBc] = m_cells->pvariables[var][cellId];
              }
            }
          }
        }

        break;
      }
    default:
      {
        zfsTerm(1, __CALLING_FUNCTION__, "Face not implemented");
        break;
      }
    }
}

/** Prescribe profile BC
 *  In case of the initialStartup2601 flag is set
 *  (no matter if it was a restart or a initial startup)
 *  we load the values from the field to the GC and
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc2601(ZFSId bcId)
{
  m_block->initializeRungeKutta();
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  m_block->m_bc2601 = true;

  m_2601wave = true;

  if(m_2601wave) {
    const ZFSFloat timeConversion = 11.59597694e6;
    const ZFSFloat timeShift = 10.624277309487987/0.015795;
    //string fileName = "coefficients_actuated_surface.dat";
    ZFSInt step = 0;
    m_2601noCoeff = 5; //number of effective coefficients
    m_2601noPos = 4; //number of streamwise positions
    const ZFSString fileNames[4] = {"coefficients_actuated_surface_00138.dat",
                                    "coefficients_actuated_surface_0014452.dat",
                                    "coefficients_actuated_surface_00148.dat",
                                    "coefficients_actuated_surface_00153.dat"};
    ifstream infile1(fileNames[0].c_str());
    ZFSInt noLines =count(std::istreambuf_iterator<char>(infile1),
                          std::istreambuf_iterator<char>(), '\n');

    cout << "NoLines: " << noLines << endl;

    zfsAlloc(m_2601streamwisePos, m_2601noPos, __CALLING_FUNCTION__, F0, "m_2601streamwisePos");
    m_2601streamwisePos[0] = 0.0138;
    m_2601streamwisePos[1] = 0.014452;
    m_2601streamwisePos[2] = 0.0148;
    m_2601streamwisePos[3] = 0.0153;

    noLines++;
    zfsAlloc(m_2601effConst, noLines, m_2601noCoeff*m_2601noPos, __CALLING_FUNCTION__, F0, "m_2601effConst");

    cout << "Reading in coefficients from file..." << endl;
    for(ZFSId pos=0; pos<m_2601noPos; pos++) {
      step = 0;
      ifstream infile2(fileNames[pos].c_str());
      while (!infile2.eof() && step<noLines-1) {
        string line[5];
        for(int i=0; i<5; i++) {
          infile2 >> line[i];
        }

        for(int i=0; i<5; i++) {
          m_2601effConst[step][pos*m_2601noCoeff+i] = atof(line[i].c_str());

          // if(m_block->domainId() == 0) {
          //   cout << "step: " << step << " pos: " << pos << " coeff: " << i << endl;
          // }
        }

        step++;
      }
    }
    cout << "Reading in coefficients from file... FINISHED!" << endl;

    const ZFSFloat zeroTime = m_2601effConst[0][0];
    for(ZFSId pos=0; pos<m_2601noPos; pos++) {
      for(ZFSInt iter=0; iter<noLines; iter++) {
        m_2601effConst[iter][pos*m_2601noCoeff+0] = (timeConversion*(m_2601effConst[iter][pos*m_2601noCoeff+0]-zeroTime)) - timeShift; //+m_block->m_time
      }
    }

    m_2601noSteps = step;
  }

  switch(m_physicalBCMap[bcId]->face)
    {
    case 2:
      {
        //Copy into zeroth order solution field
        for(ZFSInt i = start[0]; i<end[0]; i++) {
          for(ZFSInt j = start[1]; j<end[1]; j++) {
            for(ZFSInt k = start[2]; k<end[2]; k++) {
              const ZFSId cellId = cellIndex(i,j,k);
              const ZFSId cellIdBc = i + (j + k*m_noGhostLayers)*m_block->m_nCells[2];
              for(ZFSId var=0; var<PV->noVariables; var++) {
                m_block->m_bc2601ZerothOrderSolution[var][cellIdBc] = m_cells->pvariables[var][cellId];
              }
            }
          }
        }


        //Copy the values from the ghostcells into restart field
        for(ZFSInt i = m_noGhostLayers; i<end[0]-m_noGhostLayers; i++) {
          for(ZFSInt j = 0; j<m_noGhostLayers; j++) {
            for(ZFSInt k = m_noGhostLayers; k<end[2]-m_noGhostLayers; k++) {
              const ZFSId cellId = cellIndex(i,j,k);
              const ZFSId cellIdBc = i-m_noGhostLayers + (j + (k - m_noGhostLayers)*m_noGhostLayers)*m_block->m_nActiveCells[2];


              for(ZFSId var=0; var<PV->noVariables; var++) {
                m_block->m_bc2601Variables[var][cellIdBc] = m_cells->pvariables[var][cellId];
              }
            }
          }
        }
        break;
      }
    default:
      {
        zfsTerm(1, __CALLING_FUNCTION__, "Face not implemented");
        break;
      }
    }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc2300(ZFSId bcId) {
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  const ZFSFloat t0=0;
  const ZFSFloat tN=10;
  const ZFSFloat h=0.05;
  const ZFSInt N=(tN-t0)/h; //no of iteration steps
  const ZFSFloat k=0.334;

  ZFSFloatScratchSpace A(3,3, __CALLING_FUNCTION__, "A");
  ZFSFloatScratchSpace y(3, __CALLING_FUNCTION__, "y");
  ZFSFloatScratchSpace k1(3, __CALLING_FUNCTION__, "k1");
  ZFSFloatScratchSpace k2(3, __CALLING_FUNCTION__, "k2");
  ZFSFloatScratchSpace k3(3, __CALLING_FUNCTION__, "k3");
  ZFSFloatScratchSpace k4(3, __CALLING_FUNCTION__, "k4");
  ZFSFloatScratchSpace result(N+1,3, __CALLING_FUNCTION__, "result");
  ZFSFloatScratchSpace t(N+1, __CALLING_FUNCTION__, "t");

  y(0) = F0;
  y(1) = F0;
  y(2) = k;

  for(ZFSId j=0;j<3;j++) {
    result(0,j) = y(j);
  }

  for(ZFSId n=0; n<N; n++) {
    A(0,0) = 0;
    A(0,1) = 1;
    A(0,2) = 0;
    A(1,0) = 0;
    A(1,1) = 0;
    A(1,2) = 1;
    A(2,0) = -y(2)/2.0;
    A(2,1) = 0;
    A(2,2) = 0;

    //k1=A*y
    for(ZFSId i=0;i<3;i++) {
      k1(i) = F0;
      for(ZFSId j=0;j<3;j++) {
        k1(i) += A(i,j)*y(j);
      }
    }

    //k2=A*(y+(h/2)*k1)
    for(ZFSId i=0;i<3;i++) {
      k2(i) = F0;
      for(ZFSId j=0;j<3;j++) {
        k2(i) += A(i,j)*(y(j)+(h/F2)*k1(j));
      }
    }

    //k3=A*(y+(h/2)*k2)
    for(ZFSId i=0;i<3;i++) {
      k3(i) = F0;
      for(ZFSId j=0;j<3;j++) {
        k3(i) += A(i,j)*(y(j)+(h/F2)*k2(j));
      }
    }

    //k4=A*(y+h*k3)
    for(ZFSId i=0;i<3;i++) {
      k4(i) = F0;
      for(ZFSId j=0;j<3;j++) {
        k4(i) += A(i,j)*(y(j)+h*k3(j));
      }
    }

    //y=y+(h/6)*(k1+2*k2+2*k3+k4)
    for(ZFSId j=0;j<3;j++) {
      y(j) += (h/6.0)*(k1(j)+F2*k2(j)+F2*k3(j)+k4(j));
      result(n+1,j) = y(j);
    }
    t(n+1) = (n+1)*h;
  }

  ZFSFloatScratchSpace bl(6,N+1, __CALLING_FUNCTION__, "bl");

  //position at which ReTheta is defined and ReX is known:
  const ZFSFloat rex = 100000;
  const ZFSFloat rexPos = 420.0;

  //transform to position of the boundary condition
  const ZFSFloat myPos = 420;
  const ZFSFloat myRe = rex/rexPos*myPos;
  const ZFSFloat lengthScaling = myPos/sqrt(myRe);

  for(ZFSId n=0;n<N+1; n++) {
    const ZFSFloat rho = CV->rhoInfinity;
    bl(PV->U,n) = PV->UInfinity*result(n,1);
    bl(PV->V,n) = PV->UInfinity*result(n,0)/sqrt(rex);
    cout << "n: " << n << " bl(u): " << bl(PV->U,n) << endl;
    bl(PV->W,n) = F0;
    bl(PV->P,n) = PV->PInfinity;
    bl(PV->RHO,n) = rho;
    bl(5,n) = lengthScaling*t(n);
  }

  //interpolation part
  cout.precision(10);
  FILE* f_blasius;
  f_blasius = fopen("./blasius.dat", "a+");
  ZFSFloat theta = F0;
  for(ZFSId j = start[1]; j<end[1]; j++) {
    ZFSFloat yCoord = m_cells->coordinates[1][cellIndex(start[0],j,start[2])];
    ZFSFloatScratchSpace pvariables(5, __CALLING_FUNCTION__, "pvariables");
    ZFSBool hasIntStencil = false;
    for(ZFSId n=0;n<N; n++) {
      if(bl(5,n)<yCoord&&bl(5,n+1)>=yCoord) {
        for(ZFSId var=0; var<PV->noVariables; var++) {
          pvariables(var) = bl(var,n) + (yCoord-bl(5,n))/(bl(5,n+1)-bl(5,n))*(bl(var,n+1)-bl(var,n));
          if(var==0)
            cout << "yCoord: " << yCoord << " var: " << var << " ybl(n): " << bl(5,n) << " ybl(n+1): " << bl(5,n+1) << " var bl(n): " << bl(var,n) << " var bl(n+1): " << bl(var,n+1) << " interpolated var: " << pvariables(var) << endl;
        }
        hasIntStencil = true;
        break;
      }
    }

    if(!hasIntStencil) {
      if(yCoord<=bl(5,0)) {
        for(ZFSId var=0; var<PV->noVariables; var++) {
          pvariables(var) = bl(var,0);
        }
      } else {
        for(ZFSId var=0; var<PV->noVariables; var++) {
          pvariables(var) = bl(var,N);
        }
      }
    }

    ZFSId p1 = getPointIdFromCell(start[0],j,start[2]);
    ZFSId p2 = getPointIdfromPoint(p1,0,1,0);
    ZFSFloat yDelta = m_coordinates[1][p2]-m_coordinates[1][p1];
    theta += pvariables(PV->U)*pvariables(PV->RHO)/CV->rhoInfinity/PV->UInfinity*(F1 - pvariables(PV->U)/PV->UInfinity)*yDelta;

    fprintf(f_blasius, "%f", yCoord);
    fprintf(f_blasius, " %f", pvariables(0));
    fprintf(f_blasius, " %f", pvariables(1));
    fprintf(f_blasius, " %f", pvariables(2));
    fprintf(f_blasius, " %f", pvariables(3));
    fprintf(f_blasius, " %f", pvariables(4));
    fprintf(f_blasius, "\n");


    for(ZFSId i = start[0]; i<end[0]; i++) {
      for(ZFSId kk = start[2]; kk<end[2]; kk++) {
        const ZFSId cellId = cellIndex(i,j,kk);
        for(ZFSId var=0; var<PV->noVariables; var++) {
          m_cells->pvariables[var][cellId] = pvariables(var);
        }
      }
    }
  }
  fclose(f_blasius);
  cout << "Theta: " << theta << endl;
}


/** \brief init for the acoustic and entropy waves
 * version: init for cut-off boundary condition
 * author: Thomas Schilden, 12.2.2015
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc2700(ZFSId bcId)
{
  TRACE();
  zfs_log << endl << "initBc2700 for bcId " << bcId << endl;
  //zfs_log << " running for " << m_cutOffBndryCndIds[bcId] << endl;

  const  ZFSFloat time = m_block->m_physicalTime;
  //cout << "restart time is: " << time << endl;
  // mach number
  m_Ma = *(ZFSContext::getProperty("Ma", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));
  //1. check the modes and allocate memory
  // modeSr
  ZFSProperty * property;
  property = ZFSContext::getProperty("modeSr", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL);
  m_modes = property->count();
  ZFSFloatScratchSpace modeSr(m_modes,__CALLING_FUNCTION__,"modeSr");
  for(ZFSInt i = 0; i < m_modes; i++)
    modeSr(i) = *(property->asFloat(i));
  // m_modeAmp
  property = ZFSContext::getProperty("modeAmp", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL);
  if(m_modes != property->count())
    zfsTerm(1, __CALLING_FUNCTION__, "modeAmp does not fit modeSr");
  zfsAlloc(m_modeAmp, m_modes, "m_modeAmp", F0, __CALLING_FUNCTION__);
  for(ZFSInt i = 0; i < m_modes; i++)
    m_modeAmp[i] = *(property->asFloat(i));
  // m_modeType
  property = ZFSContext::getProperty("modeType", m_blockId, __CALLING_FUNCTION__, (ZFSInt*) NULL);
  if(m_modes != property->count())
    zfsTerm(1, __CALLING_FUNCTION__, "modeType does not fit modeSr");
  zfsAlloc(m_modeType, m_modes, "m_modeType", 0, __CALLING_FUNCTION__);
  for(ZFSInt i = 0; i < m_modes; i++)
    m_modeType[i] = *(property->asInt(i));
   // m_modePhi
  property = ZFSContext::getProperty("modePhi", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL);
  if(m_modes != property->count())
    zfsTerm(1, __CALLING_FUNCTION__, "modePhi does not fit modeSr");
  zfsAlloc(m_modePhi, m_modes, "m_modePhi", F0, __CALLING_FUNCTION__);
  for(ZFSInt i = 0; i < m_modes; i++)
    m_modePhi[i] = *(property->asFloat(i)) * PI / 180.0;
  // modeAngle
  property = ZFSContext::getProperty("modeAngle", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL);
  if(m_modes*2 != property->count())
    zfsTerm(1, __CALLING_FUNCTION__, "modeAngle does not fit modeSr");
  ZFSFloatScratchSpace modeAngle(m_modes, 2, __CALLING_FUNCTION__,"modeAngle");
  modeAngle.fill(F0);
  for(ZFSInt i = 0; i < m_modes; i++)
    for(ZFSInt j= 0; j < nDim - 1; j++)
      modeAngle(i,j) = *(property->asFloat(i*2+j));
  // nmbrOfModes
  property = ZFSContext::getProperty("nmbrOfModes", m_blockId, __CALLING_FUNCTION__, (ZFSInt*) NULL);
  if(m_modes != property->count())
    zfsTerm(1, __CALLING_FUNCTION__, "nmbrOfModes does not fit modeSr");
  zfsAlloc(m_nmbrOfModes, m_modes, "m_nmbrOfModes", 0, __CALLING_FUNCTION__);
  for(ZFSInt i = 0; i < m_modes; i++)
    m_nmbrOfModes[i] = *(property->asInt(i));
  // memberVariables
  zfsAlloc(m_modeOmega, m_modes, "m_modeOmega", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_modeEtaMin, m_modes, "m_modeEtaMin", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_modeK, m_modes, nDim, "m_modeK", F0, __CALLING_FUNCTION__);

  // non-mode specific stuff
  ZFSFloat initTime;

  m_block->m_restartBc2800 = F0; // has to be changed, if you add BC2800

  if(m_block->m_restartBc2800)
    initTime = m_block->m_restartTimeBc2800;
  else{
    initTime = time;
    m_block->m_restartTimeBc2800 = time;
  }

  zfs_log <<   "time        = " << time << endl;
  zfs_log <<   "initTime    = " << initTime << endl;

  ZFSFloat UInfinity = F0;
  for(ZFSInt i = 0; i < nDim; i++)
    UInfinity += POW2(PV->VVInfinity[i]);
  UInfinity = sqrt(UInfinity);

  zfs_log << "block: UInfinity " << UInfinity << endl << "block: m_referenceLength " << m_block->m_referenceLength << endl;
  for(ZFSInt mode = 0; mode < m_modes; mode++){
    zfs_log << "-- mode " << mode << " --" << endl;
    zfs_log << "   modeType = " << m_modeType[mode] << endl;
    zfs_log << "   modeSr = " << modeSr(mode) << endl;
    zfs_log << "   modeAmp = " << m_modeAmp[mode] << endl;
    zfs_log << "   modePhi = " << m_modePhi[mode] << endl;
    zfs_log << "   nmbrOfModes= " << m_nmbrOfModes[mode] << endl;
    for( ZFSInt i = 0; i < (nDim - 1); i++ )
      zfs_log << "   modeAngle(" << i << ") = " << modeAngle(mode,i) << " (rad)" << endl;
    //3. calculate nondimensional mode properties
    //3.1 omega
    m_modeOmega[mode] = F2 * PI * modeSr(mode) * UInfinity / m_block->m_referenceLength;
    //3.2 wave vector modeK
    m_modeK[mode][0] =  cos( modeAngle(mode,0) ) * cos( modeAngle(mode,1) );
    m_modeK[mode][1] =  sin( modeAngle(mode,0) ) * cos( modeAngle(mode,1) );
    if(nDim == 3)
      m_modeK[mode][2] =  sin( modeAngle(mode,1) );

    ZFSFloat Uk = F0;
    for(ZFSInt i = 0; i < nDim; i++)
      Uk += PV->VVInfinity[i] * m_modeK[mode][i];

    const ZFSFloat propVel = (ZFSFloat)(m_modeType[mode]) * sqrt( PV->TInfinity );
    const ZFSFloat K = m_modeOmega[mode] / ( Uk + propVel );

    for(ZFSInt i = 0; i < nDim; i++)
      m_modeK[mode][i] *= K;
    // output omega and k
    zfs_log << "   Uk = " << Uk << endl;
    zfs_log << "   modeOmega = " << m_modeOmega[mode] << endl;
    zfs_log << "   K = " << K << endl;
    for(ZFSInt i = 0; i < nDim; i++)
      zfs_log << "   modeK[" << i << "] = " << m_modeK[mode][i] << endl;
    //4. get the min wave phase
      ZFSFloat modeEtaMin = F0;

    // here new loop...
        ZFSInt* start = m_physicalBCMap[bcId]->start1;
        ZFSInt* end = m_physicalBCMap[bcId]->end1;
        ZFSId cellId=-1;

        for(ZFSId k=start[2]; k<end[2] ; k++)
          {
            for(ZFSId j=start[1]; j<end[1] ; j++)
              {
                for(ZFSId i=start[0]; i<end[0]; i++)
                  {
                    cellId=cellIndex(i,j,k);

                    // for( ZFSId id = 0; id < m_sortedCutOffCells[bcId]->size(); id++ ){
                    // ZFSId cellId = m_sortedCutOffCells[bcId]->a[ id ];

                    ZFSFloat eta = F0;
                    for(ZFSInt dim = 0; dim < nDim; dim++){
                      //eta += m_block->a_coordinate( cellId , i) * m_modeK[mode][i];
                      eta += m_cells->coordinates[dim][cellId] * m_modeK[mode][dim];
                    }
                    eta -= m_modeOmega[mode] * initTime;
                    //if(cellId == 0)
                    if((k==start[2]) && (j==start[1]) && (i==start[0]))
                      modeEtaMin = eta;
                    else
                      modeEtaMin = zfsMIN(modeEtaMin, eta);
                  }
              }
          }

    m_modeEtaMin[mode] = modeEtaMin;

    }
  zfs_log << "leaving the initialization" << endl;
}


/**
 *      Synthetic Turbulence Generation
 *
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc7909(ZFSId bcId)
{
  const ZFSId NU_T = 5;
  // const ZFSId SIJSIJ = 6;
  // const ZFSId LENGTH_SCALE = 7;
  // const ZFSId FLUC_UU = 8;
  // const ZFSId FLUC_VV = 9;
  // const ZFSId FLUC_WW = 10;
  // const ZFSId FLUC_UV = 11;
  // const ZFSId FLUC_VW = 12;
  // const ZFSId FLUC_UW = 13;
  // const ZFSId FLUC_U = 14;
  // const ZFSId FLUC_V = 15;
  // const ZFSId FLUC_W = 16;

  ZFSInt *commStgRoot = m_block->m_commStgRoot;
  MPI_Comm commStg = *m_block->m_commStg;
  MPI_Comm_rank(commStg, &m_block->m_stgMyRank);

  //initialize RK Step to zero
  m_block->m_RKStep = 0;

  zfsAlloc(m_stgMaxVel, nDim, "m_stgMaxVel", -99999.9, __CALLING_FUNCTION__);

  const ZFSId noCellsJ = m_block->m_totalGridBlockDim[m_block->m_inputBlockId][1]-1;
  if(m_block->m_stgMyRank == *commStgRoot) { cout << "Init stgGlobalLengthscales with " << noCellsJ << " cells" << endl; }
  zfsAlloc(m_stgGlobalLengthScales, noCellsJ, 6, "m_stgGlobalLengthScales", F0, __CALLING_FUNCTION__);

  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  if(m_block->m_stgMyRank == *commStgRoot) {cout << "Initializing BC 7909..." << endl;}

  //Compute size of the Box, take inflow as middle x-coordinate
  ZFSFloatScratchSpace bcast_vb(6, __CALLING_FUNCTION__, "bcast_vb");

  ZFSId minIndex = getPointIdFromCell(m_noGhostLayers,m_noGhostLayers,m_noGhostLayers);
  ZFSId maxIndex = getPointIdFromCell(m_noGhostLayers,m_noGhostLayers,m_noGhostLayers+m_block->m_nActiveCells[0]);
  ZFSFloat inflowStartLocal[3] = {m_coordinates[0][minIndex],
                                   m_coordinates[1][minIndex],
                                   m_coordinates[2][minIndex]};
  ZFSFloat inflowEndLocal[3] = {m_coordinates[0][maxIndex],
                                m_coordinates[1][maxIndex],
                                m_coordinates[2][maxIndex]};
  ZFSFloat inflowStart[3] = {99999.9,99999.9,99999.9};
  ZFSFloat inflowEnd[3] = {F0,F0,F0};

  MPI_Allreduce(inflowStartLocal,inflowStart,nDim, MPI_DOUBLE, MPI_MIN, commStg);
  MPI_Allreduce(inflowEndLocal,inflowEnd,nDim, MPI_DOUBLE, MPI_MAX, commStg);

  const ZFSFloat vbDepth = (inflowEnd[2]-inflowStart[2])*m_block->m_stgBLT3;
  const ZFSFloat zOffset = F1B2*(inflowEnd[2]-inflowStart[2])*(m_block->m_stgBLT3-F1);

  if(m_block->m_stgMyRank == *commStgRoot) {
    m_block->m_stgRootRank = true;

    //Get the coordinate of the inflow
    bcast_vb[0] = inflowStart[0] - m_block->m_stgBLT1*F1B2;
    bcast_vb[1] = inflowStart[1];
    bcast_vb[2] = inflowStart[2] - zOffset;

    bcast_vb[3] = inflowStart[0] + m_block->m_stgBLT1*F1B2;
    bcast_vb[4] = inflowStart[1] + m_block->m_stgBLT2;
    bcast_vb[5] = inflowStart[2] -zOffset + vbDepth;

    if(m_block->m_stgBLT3 > 1.5) {
      stringstream errorMsg;
      errorMsg << "The factor for the depth of the virtual eddy box is probably too large, BLT3: "  << m_block->m_stgBLT3
               << ". This has been changed in the last commit, instead of setting the absolute depth "
               << "in z-direction of the domain, now only a scaling factor needs to be given. A factor "
               << "of 1.0 will make the virtual box as wide as your domain in z-direction." << endl;
      zfsTerm(1, __CALLING_FUNCTION__, errorMsg.str());
    }
  }

  MPI_Bcast(bcast_vb.begin(), 6, MPI_DOUBLE, *commStgRoot, commStg);

  m_stgVbStart[0] = bcast_vb[0];
  m_stgVbStart[1] = bcast_vb[1];
  m_stgVbStart[2] = bcast_vb[2];
  m_stgVbEnd[0] = bcast_vb[3];
  m_stgVbEnd[1] = bcast_vb[4];
  m_stgVbEnd[2] = bcast_vb[5];

  if(m_block->m_stgMyRank == *commStgRoot) {
    cout << "STG Virtual Box depth factor: " << m_block->m_stgBLT3 << ", total depth: " << vbDepth << endl
         << "STG Virtual Box Information:" << endl
         << "STG Virtual Box Start X: " <<  bcast_vb[0] << " Y: " <<  bcast_vb[1] << " Z: " <<  bcast_vb[2] << endl
         << "STG Virtual Box End X: " <<  bcast_vb[3] << " Y: " <<  bcast_vb[4] << " Z: " <<  bcast_vb[5] << endl;
  }

  //Done computing virtual box
   if(!m_block->m_zonal) {
    //junoh
    if(m_block->m_stgInitialStartup) {
      //Reading in RANS Profile for BC
      //Take 3!!! rows in i-direction because of
      //gradient calculations
      for(ZFSId k=start[2]; k<end[2] ; k++) {
	for(ZFSId j=start[1]; j<end[1] ; j++) {
	  for(ZFSId i=start[0]; i<end[0]+1; i++) {
	    const ZFSId cellId=cellIndex(i,j,k);
	    const ZFSId cellIdadj=cellIndex(m_noGhostLayers,j,k); //constant in K anyway
	    const ZFSId IBC = cellIndexBC(i,j,k);

	    //Set all variables like in the field
	    for(int var = 0; var < PV->noVariables; var++) {
	      m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdadj];
	    }

	    //fill up the STG FQ field with the RANS values
	    for(int var = 0; var < PV->noVariables; var++) {
	      m_cells->stg_fq[var][IBC] = m_cells->pvariables[var][cellId];
	    }

	    m_cells->stg_fq[NU_T][IBC] = m_cells->fq[FQ->NU_T][cellIdadj];
	  }
	}
      }
    } else {
      for(ZFSId k=start[2]; k<end[2] ; k++) {
	for(ZFSId j=start[1]; j<end[1] ; j++) {
	  for(ZFSId i=start[0]; i<end[0]; i++) {
	    const ZFSId cellId=cellIndex(i,j,k);
	    const ZFSId IBC = cellIndexBC(i,j,k);
	    for(int var = 0; var < PV->noVariables; var++) {
	      m_cells->pvariables[var][cellId] = m_cells->stg_fq[var][IBC];
	    }
	  }
	}
      }
    }
   }


  //create new eddies if it's an initial start or the createNewEddie flag is set
  if(!m_block->m_restart || m_block->m_stgCreateNewEddies) {
    ZFSFloat xk1t, xk2t, xk3t, epsik1, epsik2, epsik3;
    ZFSInt nran = m_block->m_stgMaxNoEddies;
    ZFSFloatScratchSpace bcast_eddies(nran*6, __CALLING_FUNCTION__, "bcast_eddies");

    ZFSFloat eps = 1e-7;

    if(m_block->m_stgMyRank == *commStgRoot) {
      // cout << "Creating new Eddies inside Virtual Box"  << endl;
      for(ZFSInt n = 0; n < m_block->m_stgMaxNoEddies; n++) {
        xk1t = m_stgVbStart[0] + generate_rand()*(m_stgVbEnd[0] - m_stgVbStart[0]);
        xk2t = m_stgVbStart[1] + generate_rand_weighted()*(m_stgVbEnd[1] - m_stgVbStart[1]);
        xk3t = m_stgVbStart[2] + generate_rand()*(m_stgVbEnd[2] - m_stgVbStart[2]);

        // cout << "Creating eddie at x: " << xk1t << " , y: " << xk2t << " , z: " << xk3t << endl;

        epsik1 = 2.0*generate_rand() - 1.0;
        epsik1 = epsik1/max(abs(epsik1),eps);
        epsik2 = 2.0*generate_rand() - 1.0;
        epsik2 = epsik2/max(abs(epsik2),eps);
        epsik3 = 2.0*generate_rand() - 1.0;
        epsik3 = epsik3/max(abs(epsik3),eps);

        bcast_eddies[n + nran*0] = xk1t;
        bcast_eddies[n + nran*1] = xk2t;
        bcast_eddies[n + nran*2] = xk3t;
        bcast_eddies[n + nran*3] = epsik1;
        bcast_eddies[n + nran*4] = epsik2;
        bcast_eddies[n + nran*5] = epsik3;
      }
    }

    //Broadcast the new/updated eddies to all relevant processes
    MPI_Bcast(bcast_eddies.begin(), 6*nran, MPI_DOUBLE, *commStgRoot, commStg);

    //Copy data into m_FQeddie vector
    for(ZFSInt n = 0; n < nran; n++)
    {
      m_block->m_stgEddies[n][0] = bcast_eddies[n + nran*0]; //bcast_buffer(n,var);
      m_block->m_stgEddies[n][1] = bcast_eddies[n + nran*1];
      m_block->m_stgEddies[n][2] = bcast_eddies[n + nran*2];
      m_block->m_stgEddies[n][3] = bcast_eddies[n + nran*3];
      m_block->m_stgEddies[n][4] = bcast_eddies[n + nran*4];
      m_block->m_stgEddies[n][5] = bcast_eddies[n + nran*5];
    }
  }
}


/** Channel flow / Pipe Flow
 *  Authors Pascal Meysonnat and Marian Albers
 *  Last changes Pascal Meysonnat 29.05.2017
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc2402(ZFSId bcId)
{
  //for the channel flow we need the surface area of the entering and the outflow domain
  //==>determine Channelflow surface (assuming it does not change)
  ZFSId normal=0;
  ZFSId Id=m_channelSurfacIndexMap[bcId];
  if(Id<0) cerr << "id smaller than zero ==> error" << endl;
  ZFSId* startface = &m_block->m_windowInfo->channelSurfaceIndices[Id]->start1[0];
  ZFSId* endface = &m_block->m_windowInfo->channelSurfaceIndices[Id]->end1[0];
  //find out which face it is
  ZFSId fixedind=-1;


  if(m_block->m_initialCondition == 1233) {
    /////////////////////////////////////////////
    ///////// Turbulent channel flow ////////////
    /////////////////////////////////////////////
    ZFSFloat uTau= m_block->m_ReTau*m_block->m_Ma*sqrt(PV->TInfinity)/m_block->m_Re;
    m_block->m_deltaP = -CV->rhoInfinity*POW2(uTau)*F2*(m_block->m_channelLength)/m_block->m_channelHeight;
    m_block->m_channelPresInlet=PV->PInfinity;
    m_block->m_channelPresOutlet=PV->PInfinity+m_block->m_deltaP;
    zfs_log << "=========== Turb. Channel Flow BC Summary =========== " <<endl;
    zfs_log << "-->Turbulent channel flow deltaP: " << m_block->m_deltaP << endl;
    zfs_log << "-->channel friciton velocity: " << uTau << endl;
    zfs_log << "-->Channel pressure inflow: " << m_block->m_channelPresInlet << endl;
    zfs_log << "-->Channel pressure outflow: " << m_block->m_channelPresOutlet << endl;
     zfs_log << "=========== Turb. Channel Flow BC Summary Finished =========== " <<endl;
  } else if(m_block->m_initialCondition == 1234) {
    cout.precision(10);
    /////////////////////////////////////////////
    ///////// Laminar channel flow //////////////
    /////////////////////////////////////////////

    m_block->m_deltaP = -12.0*PV->UInfinity*zfsSUTHERLANDLAW(PV->TInfinity)*m_block->m_channelLength/(POW2(m_block->m_channelHeight)*m_block->m_Re0);
    zfs_log << "=========== Lam. Channel Flow BC Summary =========== " <<endl;
    zfs_log << "Laminar channel deltaP: " << m_block->m_deltaP << endl;
    zfs_log << "Theoretical cD total (both channel walls): " << m_block->m_deltaP*m_block->m_channelHeight*m_block->m_channelWidth/(0.5*CV->rhoInfinity*POW2(PV->UInfinity)) << endl;
    zfs_log << "Theoretical cD (single wall): " << m_block->m_deltaP*m_block->m_channelHeight*m_block->m_channelWidth/(0.5*CV->rhoInfinity*POW2(PV->UInfinity)*2.0) << endl;
    zfs_log << "Theoretical cF total (both channel walls): " << m_block->m_deltaP*m_block->m_channelHeight/(0.5*CV->rhoInfinity*POW2(PV->UInfinity)*m_block->m_channelLength) << endl;
    zfs_log << "Theoretical cF (single wall): " << m_block->m_deltaP*m_block->m_channelHeight/(0.5*CV->rhoInfinity*POW2(PV->UInfinity)*m_block->m_channelLength*2.0) << endl;
    zfs_log << "=========== Lam. Channel Flow BC Summary Finished =========== " <<endl;
  } else if(m_block->m_initialCondition == 1236) {

    /////////////////////////////////////////////
    ///////// Turbulent pipe flow ///////////////
    /////////////////////////////////////////////
    ZFSFloat uTau= m_block->m_ReTau*m_block->m_Ma*sqrt(PV->TInfinity)/m_block->m_Re;
    m_block->m_deltaP = -4.0*CV->rhoInfinity*POW2(uTau)*(m_block->m_channelLength)/m_block->m_channelHeight;

    m_block->m_channelPresInlet=PV->PInfinity;
    m_block->m_channelPresOutlet=PV->PInfinity+m_block->m_deltaP;

    zfs_log << "=========== Turb. Pipe Flow BC Summary =========== " <<endl;
    zfs_log << "-->Turbulent pipe flow deltaP: " << m_block->m_deltaP << endl;
    zfs_log << "-->pipe friciton velocity: " << uTau << endl;
    zfs_log << "-->pipe pressure inflow: " << m_block->m_channelPresInlet << endl;
    zfs_log << "-->pipe pressure outflow: " << m_block->m_channelPresOutlet << endl;
    zfs_log << "=========== Turb. Pipe Flow BC Summary Finished =========== " <<endl;
  }

  for(ZFSId dim=0; dim<nDim; dim++) {
    if(startface[dim]==endface[dim]) {
      //this is the normal
      fixedind=dim;
      if(startface[dim]==m_noGhostLayers) {
        normal=-1;
      } else {
        normal=1;
      }
    }
  }

  if(m_physicalBCMap[bcId]->face==0) {
    MPI_Comm commChannelIn = *m_block->m_commChannelIn;
    MPI_Comm_rank(commChannelIn, &m_channelInflowRank);
  }

  switch(fixedind)
    {
    case 0:
      {
        ZFSFloat surface=F0;
        ZFSInt ii=0;
        if(normal<0) {
          ii=startface[0];

          ZFSId cellId=0;
          for(ZFSId k=startface[2]; k<endface[2] ; k++) {
            for(ZFSId j=startface[1]; j<endface[1] ; j++) {
              cellId=cellIndex(ii,j,k);
              surface+=sqrt(POW2(m_cells->cellMetrics[cellId][0])+POW2(m_cells->cellMetrics[cellId][1])+POW2(m_cells->cellMetrics[cellId][2]));
            }
          }
        } else {
          ii=endface[0];
          //activate this or the method below with the 4 points for the surface calculation
          //this is less exact (for straight surf.) but works better for curved surfaces

          ZFSId cellId=0;
          for(ZFSId k=startface[2]; k<endface[2] ; k++) {
            for(ZFSId j=startface[1]; j<endface[1] ; j++) {
              cellId=cellIndex(ii-1,j,k);
              surface+=sqrt(POW2(m_cells->cellMetrics[cellId][0])+POW2(m_cells->cellMetrics[cellId][1])+POW2(m_cells->cellMetrics[cellId][2]));
            }
          }
        }

        if(normal<0) {
          MPI_Allreduce(&surface, &m_channelSurfaceIn, 1, MPI_DOUBLE, MPI_SUM, m_block->m_commChannelIn[0]);
          cout << "ChannelInSurface: " << m_channelSurfaceIn << endl;
        } else {
          MPI_Allreduce(&surface, &m_channelSurfaceOut, 1, MPI_DOUBLE, MPI_SUM, m_block->m_commChannelOut[0]);
          cout << "ChannelOutSurface: " << m_channelSurfaceOut << endl;
        }

        break;
      }
    case 1:
      {
        zfsTerm(1, __CALLING_FUNCTION__, "surface calculation for j faces(channel not implemented)");
        break;
      }
    case 2:
      {
        zfsTerm(1, __CALLING_FUNCTION__, "surface calculation for k faces(channel not implemented)");
        break;
      }
    default:
      {
        zfsTerm(1, __CALLING_FUNCTION__, "surface calculation for given faces(channel not implemented)");
      }
    }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::initBc4001(ZFSId bcId){
  //determine the angle for the periodic rotation etc
  //check for consistency

  //1) compute normals one and two
  //2) broadcast to own surface and to rotationWorld
  //3) check rotation axis and compute rotation angle
  //4) compute and broadcast rotation matrix

  if(!(m_rotBC.firstRunBcId == bcId && m_rotBC.initRun == 2))
    {

      ZFSFloat epsilon=0.0000000001;
      ZFSFloat length=F0;
      ZFSFloat noBlades;
      ZFSFloat normalNow[3]={F0,F0,F0};
      ZFSFloat normalLocal[3]={F0,F0,F0};
      ZFSFloat normalOne[3]={F0,F0,F0};
      ZFSFloat normalTwo[3]={F0,F0,F0};

      ZFSFloat rotMatrix[3][3];

      if((m_block->m_commPerRotGroup)==0 || (m_block->m_commPerRotGroup)>3){
        zfsTerm(1, __CALLING_FUNCTION__, "no such periodic rotation group");
      }

      //1)
      ZFSInt* startInfo = m_physicalBCMap[bcId]->start1;
      ZFSInt* endInfo = m_physicalBCMap[bcId]->end1;

      ZFSId point=-1, point1=-1, point2=-1, point3=-1;
      ZFSFloat surfVect1[3]={F0,F0,F0};
      ZFSFloat surfVect2[3]={F0,F0,F0};
      ZFSInt i=0, j=0, k=0;

      switch(m_physicalBCMap[bcId]->face){
      case 0:{
        i=startInfo[0]+m_noGhostLayers;
        for(j=startInfo[1]+m_noGhostLayers;j<=endInfo[1]-m_noGhostLayers;++j){
          for(k=startInfo[2]+m_noGhostLayers;k<=endInfo[2]-m_noGhostLayers;++k){
            point=getPointIdFromCell(i,j,k);
            point1=getPointIdfromPoint(point,0,0,0);
            point2=getPointIdfromPoint(point,0,0,1);
            point3=getPointIdfromPoint(point,0,1,0);

            surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
            surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
            surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
            surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
            surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
            surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

            normalNow[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
            normalNow[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
            normalNow[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];
          }
        }

        length=sqrt(POW2(normalNow[0])+POW2(normalNow[1])+POW2(normalNow[2]));
        normalNow[0]/=length;
        normalNow[1]/=length;
        normalNow[2]/=length;

        for(j=startInfo[0]+m_noGhostLayers;j<=endInfo[0]-m_noGhostLayers;++j){
          for(k=startInfo[1]+m_noGhostLayers;k<=endInfo[1]-m_noGhostLayers;++k){
            point=getPointIdFromCell(i,j,k);
            point1=getPointIdfromPoint(point,0,0,0);
            point2=getPointIdfromPoint(point,0,0,1);
            point3=getPointIdfromPoint(point,0,1,0);

            surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
            surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
            surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
            surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
            surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
            surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

            normalLocal[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
            normalLocal[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
            normalLocal[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];

            length=sqrt(POW2(normalLocal[1]*normalNow[2]-normalLocal[2]*normalNow[1])+POW2(normalLocal[2]*normalNow[0]-normalLocal[0]*normalNow[2])+POW2(normalLocal[0]*normalNow[1]-normalLocal[1]*normalNow[0]));
            if(length>epsilon){exit(1);}
          }
        }
        break;
      }
      case 1:{
        i=endInfo[0]-m_noGhostLayers-1;
        for(j=startInfo[1]+m_noGhostLayers;j<=endInfo[1]-m_noGhostLayers;++j){
          for(k=startInfo[2]+m_noGhostLayers;k<=endInfo[2]-m_noGhostLayers;++k){
            point=getPointIdFromCell(i,j,k);
            point1=getPointIdfromPoint(point,1,0,0);
            point2=getPointIdfromPoint(point,1,1,0);
            point3=getPointIdfromPoint(point,1,0,1);

            surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
            surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
            surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
            surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
            surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
            surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

            normalNow[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
            normalNow[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
            normalNow[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];
          }
        }

        length=sqrt(POW2(normalNow[0])+POW2(normalNow[1])+POW2(normalNow[2]));
        normalNow[0]/=length;
        normalNow[1]/=length;
        normalNow[2]/=length;

        for(j=startInfo[1]+m_noGhostLayers;j<=endInfo[1]-m_noGhostLayers;++j){
          for(k=startInfo[2]+m_noGhostLayers;k<=endInfo[2]-m_noGhostLayers;++k){
            point=getPointIdFromCell(i,j,k);
            point1=getPointIdfromPoint(point,1,0,0);
            point2=getPointIdfromPoint(point,1,1,0);
            point3=getPointIdfromPoint(point,1,0,1);

            surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
            surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
            surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
            surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
            surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
            surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

            normalLocal[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
            normalLocal[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
            normalLocal[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];

            length=sqrt(POW2(normalLocal[1]*normalNow[2]-normalLocal[2]*normalNow[1])+POW2(normalLocal[2]*normalNow[0]-normalLocal[0]*normalNow[2])+POW2(normalLocal[0]*normalNow[1]-normalLocal[1]*normalNow[0]));
            if(length>epsilon){
              exit(1);
            }
          }
        }
        break;
      }
      case 2:{
        j=startInfo[1]+m_noGhostLayers;
        for(i=startInfo[0]+m_noGhostLayers;i<=endInfo[0]-m_noGhostLayers;++i){
          for(k=startInfo[2]+m_noGhostLayers;k<=endInfo[2]-m_noGhostLayers;++k){
            point=getPointIdFromCell(i,j,k);
            point1=getPointIdfromPoint(point,0,0,0);
            point2=getPointIdfromPoint(point,1,0,0);
            point3=getPointIdfromPoint(point,0,0,1);

            surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
            surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
            surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
            surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
            surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
            surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

            normalNow[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
            normalNow[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
            normalNow[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];
          }
        }

        length=sqrt(POW2(normalNow[0])+POW2(normalNow[1])+POW2(normalNow[2]));
        normalNow[0]/=length;
        normalNow[1]/=length;
        normalNow[2]/=length;

        for(i=startInfo[0]+m_noGhostLayers;i<=endInfo[0]-m_noGhostLayers;++i){
          for(k=startInfo[2]+m_noGhostLayers;k<=endInfo[2]-m_noGhostLayers;++k){
            point=getPointIdFromCell(i,j,k);
            point1=getPointIdfromPoint(point,0,0,0);
            point2=getPointIdfromPoint(point,1,0,0);
            point3=getPointIdfromPoint(point,0,0,1);

            surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
            surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
            surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
            surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
            surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
            surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

            normalLocal[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
            normalLocal[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
            normalLocal[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];

            length=sqrt(POW2(normalLocal[1]*normalNow[2]-normalLocal[2]*normalNow[1])+POW2(normalLocal[2]*normalNow[0]-normalLocal[0]*normalNow[2])+POW2(normalLocal[0]*normalNow[1]-normalLocal[1]*normalNow[0]));
            if(length>epsilon){
              exit(1);
            }
          }
        }
        break;
      }
      case 3:{
        j=endInfo[1]-m_noGhostLayers-1;
        for(i=startInfo[0]+m_noGhostLayers;i<=endInfo[0]-m_noGhostLayers;++i){
          for(k=startInfo[2]+m_noGhostLayers;k<=endInfo[2]-m_noGhostLayers;++k){
            point=getPointIdFromCell(i,j,k);
            point1=getPointIdfromPoint(point,0,1,0);
            point2=getPointIdfromPoint(point,0,1,1);
            point3=getPointIdfromPoint(point,1,1,0);

            surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
            surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
            surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
            surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
            surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
            surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

            normalNow[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
            normalNow[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
            normalNow[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];
          }
        }

        length=sqrt(POW2(normalNow[0])+POW2(normalNow[1])+POW2(normalNow[2]));
        normalNow[0]/=length;
        normalNow[1]/=length;
        normalNow[2]/=length;

        for(i=startInfo[0]+m_noGhostLayers;i<=endInfo[0]-m_noGhostLayers;++i){
          for(k=startInfo[2]+m_noGhostLayers;k<=endInfo[2]-m_noGhostLayers;++k){
            point=getPointIdFromCell(i,j,k);
            point1=getPointIdfromPoint(point,0,1,0);
            point2=getPointIdfromPoint(point,0,1,1);
            point3=getPointIdfromPoint(point,1,1,0);

            surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
            surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
            surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
            surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
            surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
            surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

            normalLocal[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
            normalLocal[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
            normalLocal[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];

            length=sqrt(POW2(normalLocal[1]*normalNow[2]-normalLocal[2]*normalNow[1])+POW2(normalLocal[2]*normalNow[0]-normalLocal[0]*normalNow[2])+POW2(normalLocal[0]*normalNow[1]-normalLocal[1]*normalNow[0]));
            if(length>epsilon){
              exit(1);
            }
          }
        }
        break;
      }
      case 4:
        {
          k=startInfo[2]+m_noGhostLayers;
          for(i=startInfo[0]+m_noGhostLayers;i<=endInfo[0]-m_noGhostLayers;++i)
            {
              for(j=startInfo[1]+m_noGhostLayers;j<=endInfo[1]-m_noGhostLayers;++j)
                {
                  point=getPointIdFromCell(i,j,k);
                  point1=getPointIdfromPoint(point,0,0,0);
                  point2=getPointIdfromPoint(point,0,1,0);
                  point3=getPointIdfromPoint(point,1,0,0);

                  surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
                  surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
                  surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
                  surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
                  surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
                  surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

                  normalNow[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
                  normalNow[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
                  normalNow[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];
                }
            }

          length=sqrt(POW2(normalNow[0])+POW2(normalNow[1])+POW2(normalNow[2]));
          normalNow[0]/=length;
          normalNow[1]/=length;
          normalNow[2]/=length;

          for(i=startInfo[0]+m_noGhostLayers;i<=endInfo[0]-m_noGhostLayers;++i)
            {
              for(j=startInfo[1]+m_noGhostLayers;j<=endInfo[1]-m_noGhostLayers;++j)
                {
                  point=getPointIdFromCell(i,j,k);
                  point1=getPointIdfromPoint(point,0,0,0);
                  point2=getPointIdfromPoint(point,0,1,0);
                  point3=getPointIdfromPoint(point,1,0,0);

                  surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
                  surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
                  surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
                  surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
                  surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
                  surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

                  normalLocal[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
                  normalLocal[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
                  normalLocal[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];

                  length=sqrt(POW2(normalLocal[1]*normalNow[2]-normalLocal[2]*normalNow[1])+POW2(normalLocal[2]*normalNow[0]-normalLocal[0]*normalNow[2])+POW2(normalLocal[0]*normalNow[1]-normalLocal[1]*normalNow[0]));
                  if(length>epsilon)
                    {
                      cout << "ERROR: Discrepancy between local and average normal surface vector @ cpu: " << m_block->domainId() << " , (i=" << i << ",j=" << j << ",k=" << k << ")" << endl;
                      exit(1);
                    }
                }
            }
          break;
        }
      case 5:
        {
          k=endInfo[2]-m_noGhostLayers-1;
          for(i=startInfo[0]+m_noGhostLayers;i<=endInfo[0]-m_noGhostLayers;++i)
            {
              for(j=startInfo[1]+m_noGhostLayers;j<=endInfo[1]-m_noGhostLayers;++j)
                {
                  point=getPointIdFromCell(i,j,k);
                  point1=getPointIdfromPoint(point,0,0,1);
                  point2=getPointIdfromPoint(point,1,0,1);
                  point3=getPointIdfromPoint(point,0,1,1);

                  surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
                  surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
                  surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
                  surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
                  surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
                  surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

                  normalNow[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
                  normalNow[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
                  normalNow[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];
                }
            }

          length=sqrt(POW2(normalNow[0])+POW2(normalNow[1])+POW2(normalNow[2]));
          normalNow[0]/=length;
          normalNow[1]/=length;
          normalNow[2]/=length;

          for(i=startInfo[0]+m_noGhostLayers;i<=endInfo[0]-m_noGhostLayers;++i)
            {
              for(j=startInfo[1]+m_noGhostLayers;j<=endInfo[1]-m_noGhostLayers;++j)
                {
                  point=getPointIdFromCell(i,j,k);
                  point1=getPointIdfromPoint(point,0,0,1);
                  point2=getPointIdfromPoint(point,1,0,1);
                  point3=getPointIdfromPoint(point,0,1,1);

                  surfVect1[0]=m_coordinates[0][point2]-m_coordinates[0][point1];
                  surfVect1[1]=m_coordinates[1][point2]-m_coordinates[1][point1];
                  surfVect1[2]=m_coordinates[2][point2]-m_coordinates[2][point1];
                  surfVect2[0]=m_coordinates[0][point3]-m_coordinates[0][point1];
                  surfVect2[1]=m_coordinates[1][point3]-m_coordinates[1][point1];
                  surfVect2[2]=m_coordinates[2][point3]-m_coordinates[2][point1];

                  normalLocal[0]+=surfVect1[1]*surfVect2[2]-surfVect1[2]*surfVect2[1];
                  normalLocal[1]+=surfVect1[2]*surfVect2[0]-surfVect1[0]*surfVect2[2];
                  normalLocal[2]+=surfVect1[0]*surfVect2[1]-surfVect1[1]*surfVect2[0];

                  length=sqrt(POW2(normalLocal[1]*normalNow[2]-normalLocal[2]*normalNow[1])+POW2(normalLocal[2]*normalNow[0]-normalLocal[0]*normalNow[2])+POW2(normalLocal[0]*normalNow[1]-normalLocal[1]*normalNow[0]));
                  if(length>epsilon)
                    {
                      cout << "ERROR: Discrepancy between local and average normal surface vector @ cpu: " << m_block->domainId() << " , (i=" << i << ",j=" << j << ",k=" << k << ")" << endl;
                      exit(1);
                    }
                }
            }
          break;
        }
      default:{
        zfsTerm(1, __CALLING_FUNCTION__, "no such face");
      }

        //cout << m_block->domainId() << " initBc4004 computing normal surface vector for face " << m_physicalBCMap[bcId]->face << ": (" << normalNow[0] << ";" << normalNow[1] << "," << normalNow[2] << ")" << endl;


      }


      //2)

      if(m_block->m_commPerRotGroup==1||m_block->m_commPerRotGroup==2||m_block->m_commPerRotGroup==3)
        {
          normalOne[0]=normalNow[0];
          normalOne[1]=normalNow[1];
          normalOne[2]=normalNow[2];
          normalTwo[0]=normalNow[0];
          normalTwo[1]=normalNow[1];
          normalTwo[2]=normalNow[2];
        }
      switch(m_block->m_commPerRotGroup)
        {
        case 0:
          {
            break;
          }
        default:
          {
            //Check if we are a domain with both ends of rotating boundary
            //in that case we will run through this initBc twice, but we should
            //only do the Bcast once because no one else will participate
            if(m_rotBC.initRun == 1)
              {
                cout << "First run of bcInit4001 bcId: " << m_physicalBCMap[bcId]->BC << " , commPerRotGroup: " << m_block->m_commPerRotGroup << endl;
                MPI_Bcast(normalOne,3,MPI_DOUBLE,m_block->m_commPerRotRoots[2],*(m_block->m_commPerRotWorld));
                MPI_Bcast(normalTwo,3,MPI_DOUBLE,m_block->m_commPerRotRoots[3],*(m_block->m_commPerRotWorld));

                //Save the normal of the current side if
                //it's a domain with both in and out
                //also set the secondRun flag to true
                //so no second Bcast will be done
                if(m_block->m_commPerRotGroup == 3)
                  {
                    if(m_physicalBCMap[bcId]->BC == 4001)
                      {
                        m_rotBC.prevNormal[0] = normalOne[0];
                        m_rotBC.prevNormal[1] = normalOne[1];
                        m_rotBC.prevNormal[2] = normalOne[2];
                      }
                    else if(m_physicalBCMap[bcId]->BC == 4002)
                      {
                        m_rotBC.prevNormal[0] = normalTwo[0];
                        m_rotBC.prevNormal[1] = normalTwo[1];
                        m_rotBC.prevNormal[2] = normalTwo[2];
                      }
                  }

                m_rotBC.firstRunBcId = bcId;
              }
            else
              {
                if(m_rotBC.initRun == 2)
                  cout << "Second run of bcInit4001, must be a domain containing rotationOne AND rotationTwo, bcId: " << m_physicalBCMap[bcId]->BC << endl;
                else if(m_rotBC.initRun == 3)
                  cout << "Third run of bcInit4001, must be a domain containing rotationOne AND rotationTwo, bcId: " << m_physicalBCMap[bcId]->BC << endl;

                //We didn't get the second normal because
                //no Bcast was done, so we take the
                //previously saved one, and the other
                //one for the third run
                if(m_physicalBCMap[bcId]->BC == 4002)
                  {
                    normalOne[0] = m_rotBC.prevNormal[0];
                    normalOne[1] = m_rotBC.prevNormal[1];
                    normalOne[2] = m_rotBC.prevNormal[2];

                    m_rotBC.prevNormal[0] = normalTwo[0];
                    m_rotBC.prevNormal[1] = normalTwo[1];
                    m_rotBC.prevNormal[2] = normalTwo[2];
                  }
                else if(m_physicalBCMap[bcId]->BC == 4001)
                  {
                    normalTwo[0] = m_rotBC.prevNormal[0];
                    normalTwo[1] = m_rotBC.prevNormal[1];
                    normalTwo[2] = m_rotBC.prevNormal[2];

                    m_rotBC.prevNormal[0] = normalOne[0];
                    m_rotBC.prevNormal[1] = normalOne[1];
                    m_rotBC.prevNormal[2] = normalOne[2];
                  }
              }

            // if(m_block->m_commPerRotGroup==1||m_block->m_commPerRotGroup==3){
            //   MPI_Bcast(normalNow,3,MPI_DOUBLE,m_block->m_commPerRotRoots[0],*(m_block->m_commPerRotWorld));
            // }else{
            //   MPI_Bcast(normalNow,3,MPI_DOUBLE,m_block->m_commPerRotRoots[1],*(m_block->m_commPerRotWorld));
            // }
            break;
          }
        }


      //3)

      //computing normal axis from normal vectors
      if(m_block->m_commPerRotGroup==1 || m_block->m_commPerRotGroup==3)
        {
          ZFSFloat compRotAx[3];
          compRotAx[0]=normalOne[1]*normalTwo[2]-normalOne[2]*normalTwo[1];
          compRotAx[1]=normalOne[2]*normalTwo[0]-normalOne[0]*normalTwo[2];
          compRotAx[2]=normalOne[0]*normalTwo[1]-normalOne[1]*normalTwo[0];
          length=sqrt(POW2(compRotAx[0])+POW2(compRotAx[1])+POW2(compRotAx[2]));
          if(length<epsilon){
            cout << "Warning: Could not compute rotation axis from normal vectors. Using rotation axis from properties.cdl, secondRun: " << m_rotBC.initRun << " , gDomainId: " << m_block->domainId() << endl;
          }else{
            length=sqrt(POW2(m_rotBC.perNormals[1]*compRotAx[2]-m_rotBC.perNormals[2]*compRotAx[1])+POW2(m_rotBC.perNormals[2]*compRotAx[0]-m_rotBC.perNormals[0]*compRotAx[2])+POW2(m_rotBC.perNormals[0]*compRotAx[1]-m_rotBC.perNormals[1]*compRotAx[0]));
            if(length<epsilon){
              cout << m_block->domainId() << " " << m_block->m_commPerRotGroup << " No discrepancy between computed rotation axis and rotation axis from properties.cdl" << endl;
            }else{
              zfsTerm(1, __CALLING_FUNCTION__, "ERROR: Discrepancy between computed rotation axis and rotation axis from properties.cdl");
            }
          }

          //computing angle and check if number of blades is natural number
          length=normalOne[0]*normalTwo[0]+normalOne[1]*normalTwo[1]+normalOne[2]*normalTwo[2];
          if(length<-1 && length>-1-epsilon) length=-1;
          if(length>1 && length<1+epsilon) length=1;
          m_rotBC.perAngle=M_PI-acos(length);
          noBlades=2*M_PI/m_rotBC.perAngle;
          if((abs(noBlades-roundf(noBlades))>epsilon)||(noBlades<0.5)){
            cout << "Warning: No natural number of blades";
          }else{
            cout << m_block->domainId() << " rounding..." << endl;
            noBlades=roundf(noBlades);
            m_rotBC.perAngle=2*M_PI/noBlades;

            cout << m_block->domainId() << " rotation angle: " << fabs(m_rotBC.perAngle) << " , rotation angle (deg): " << fabs(m_rotBC.perAngle)*180/M_PI << endl
                 << m_block->domainId() << " number of blades: " << 2*M_PI/fabs(m_rotBC.perAngle) << endl;
          }
        }

      //Only do the Bcasting in the first run
      //because in the higher runs not all domains
      //are participating anymore
      if(m_rotBC.initRun == 1)
        {
          cout << m_block->domainId() << "BCasting rotPerAngle..." << endl;
          MPI_Bcast(&m_rotBC.perAngle,1,MPI_DOUBLE,m_block->m_commPerRotRoots[2],*m_block->m_commPerRotWorld);

          m_rotBC.initRun = 2;
        }
      else
        {
          m_rotBC.initRun++;
        }
      //4)

      //computing rotation matrix
      if(m_block->m_commPerRotGroup==1||m_block->m_commPerRotGroup==2 || m_block->m_commPerRotGroup==3)
        {
          ZFSFloat sinus, cosinus;
          ZFSFloat normalOneRot[3];

          cosinus=cos(m_rotBC.perAngle);
          sinus=sin(m_rotBC.perAngle);
          rotMatrix[0][0]=POW2(m_rotBC.perNormals[0])*(1-cosinus)+cosinus;
          rotMatrix[0][1]=m_rotBC.perNormals[0]*m_rotBC.perNormals[1]*(1-cosinus)-m_rotBC.perNormals[2]*sinus;
          rotMatrix[0][2]=m_rotBC.perNormals[0]*m_rotBC.perNormals[2]*(1-cosinus)+m_rotBC.perNormals[1]*sinus;
          rotMatrix[1][0]=m_rotBC.perNormals[1]*m_rotBC.perNormals[0]*(1-cosinus)+m_rotBC.perNormals[2]*sinus;
          rotMatrix[1][1]=POW2(m_rotBC.perNormals[1])*(1-cosinus)+cosinus;
          rotMatrix[1][2]=m_rotBC.perNormals[1]*m_rotBC.perNormals[2]*(1-cosinus)-m_rotBC.perNormals[0]*sinus;
          rotMatrix[2][0]=m_rotBC.perNormals[2]*m_rotBC.perNormals[0]*(1-cosinus)-m_rotBC.perNormals[1]*sinus;
          rotMatrix[2][1]=m_rotBC.perNormals[2]*m_rotBC.perNormals[1]*(1-cosinus)+m_rotBC.perNormals[0]*sinus;
          rotMatrix[2][2]=POW2(m_rotBC.perNormals[2])*(1-cosinus)+cosinus;

          normalOneRot[0]=rotMatrix[0][0]*normalOne[0]+rotMatrix[0][1]*normalOne[1]+rotMatrix[0][2]*normalOne[2]+normalTwo[0];
          normalOneRot[1]=rotMatrix[1][0]*normalOne[0]+rotMatrix[1][1]*normalOne[1]+rotMatrix[1][2]*normalOne[2]+normalTwo[1];
          normalOneRot[2]=rotMatrix[2][0]*normalOne[0]+rotMatrix[2][1]*normalOne[1]+rotMatrix[2][2]*normalOne[2]+normalTwo[2];
          length=sqrt(POW2(normalOneRot[0])+POW2(normalOneRot[0])+POW2(normalOneRot[0]));

          if(length<epsilon){
            if(m_block->m_commPerRotGroup==2 || (m_block->m_commPerRotGroup==3 && m_physicalBCMap[bcId]->BC == 4002)){
              m_rotBC.perAngle*=-1;
            }
          }
          else if(m_block->m_commPerRotGroup==1 || (m_block->m_commPerRotGroup==3 && m_physicalBCMap[bcId]->BC == 4001)){
            m_rotBC.perAngle*=-1;
          }

          cosinus=cos(m_rotBC.perAngle);
          sinus=sin(m_rotBC.perAngle);
          rotMatrix[0][0]=POW2(m_rotBC.perNormals[0])*(1-cosinus)+cosinus;
          rotMatrix[0][1]=m_rotBC.perNormals[0]*m_rotBC.perNormals[1]*(1-cosinus)-m_rotBC.perNormals[2]*sinus;
          rotMatrix[0][2]=m_rotBC.perNormals[0]*m_rotBC.perNormals[2]*(1-cosinus)+m_rotBC.perNormals[1]*sinus;
          rotMatrix[1][0]=m_rotBC.perNormals[1]*m_rotBC.perNormals[0]*(1-cosinus)+m_rotBC.perNormals[2]*sinus;
          rotMatrix[1][1]=POW2(m_rotBC.perNormals[1])*(1-cosinus)+cosinus;
          rotMatrix[1][2]=m_rotBC.perNormals[1]*m_rotBC.perNormals[2]*(1-cosinus)-m_rotBC.perNormals[0]*sinus;
          rotMatrix[2][0]=m_rotBC.perNormals[2]*m_rotBC.perNormals[0]*(1-cosinus)-m_rotBC.perNormals[1]*sinus;
          rotMatrix[2][1]=m_rotBC.perNormals[2]*m_rotBC.perNormals[1]*(1-cosinus)+m_rotBC.perNormals[0]*sinus;
          rotMatrix[2][2]=POW2(m_rotBC.perNormals[2])*(1-cosinus)+cosinus;
          for(i=0; i<=2; ++i)
            {
              for(j=0; j<=2; ++j)
                {
                  //Distinguish the side for domains
                  //which contain both sides
                  //(two rotational matrices are needed then)
                  if(m_physicalBCMap[bcId]->BC == 4001)
                    m_rotBC.rotationMatrix4001[i][j]=rotMatrix[i][j];
                  else
                    m_rotBC.rotationMatrix4002[i][j]=rotMatrix[i][j];
                }
            }
        }


      cout << m_block->domainId() << " normal 1: (" << normalOne[0] << "," << normalOne[1] << "," << normalOne[2] << ")" << endl
           << m_block->domainId() << " normal 2: (" << normalTwo[0] << "," << normalTwo[1] << "," << normalTwo[2] << ")" << endl
           << m_block->domainId() << " normal now: (" << normalNow[0] << "," << normalNow[1] << "," << normalNow[2] << ")" << endl
           << m_block->domainId() << " rotation axis: (" << m_rotBC.perNormals[0] << "," << m_rotBC.perNormals[1] << "," << m_rotBC.perNormals[2] << ")" << endl
           << m_block->domainId() << " rotation angle: " << fabs(m_rotBC.perAngle) << " , rotation angle (deg): " << fabs(m_rotBC.perAngle)*180/M_PI << endl
           << m_block->domainId() << " number of blades: " << 2*M_PI/fabs(m_rotBC.perAngle) << endl;

      if(m_block->m_commPerRotGroup==1||m_block->m_commPerRotGroup==2||m_block->m_commPerRotGroup==3)
        {
          cout << endl << "rotation matrix " << m_physicalBCMap[bcId]->BC << " of cpu " << m_block->domainId() << " is" << endl;
          for(i=0; i<3; ++i)
            {
              cout << m_block->domainId() << " " << m_block->m_commPerRotGroup << " ||| ";
              for(j=0; j<3; ++j)
                {
                  if(m_physicalBCMap[bcId]->BC == 4001)
                    cout << m_rotBC.rotationMatrix4001[i][j] << " ||| ";
                  else
                    cout << m_rotBC.rotationMatrix4002[i][j] << " ||| ";
                }
              cout << endl;
            }
        }

      //If our domain contains two sides, we have to do a third run
      //with the bc that was handled in the first run because the two
      //normals weren't correct then, only from the second run on we know
      //both correct normals
      if(m_rotBC.initRun == 3)
        {
          cout << m_block->domainId() << "Starting third run of rotInit, bcId: " << m_rotBC.firstRunBcId <<  " , CurrentBC: " <<  m_physicalBCMap[bcId]->BC << endl;
          initBc4001(m_rotBC.firstRunBcId);
        }
    }
}


//solid wall bc
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc1000(ZFSId bcId)
{
  const ZFSInt* start = m_physicalBCMap[bcId]->start1;
  const ZFSInt* end = m_physicalBCMap[bcId]->end1;
  const ZFSId face = m_physicalBCMap[bcId]->face;

  for(ZFSId k=start[2]; k<end[2] ; k++) {
    for(ZFSId j=start[1]; j<end[1] ; j++) {
      for(ZFSId i=start[0]; i<end[0]; i++) {
        ZFSId cellId;
        ZFSId cellIdAdj;
        tie(cellId, cellIdAdj) = getMirrorCellIdPair(i,j,k,face);

        m_cells->pvariables[PV->RHO][cellId]=m_cells->pvariables[PV->RHO][cellIdAdj];
        m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdAdj];
        m_cells->pvariables[PV->U][cellId]=-m_cells->pvariables[PV->U][cellIdAdj];
        m_cells->pvariables[PV->V][cellId]=-m_cells->pvariables[PV->V][cellIdAdj];
        m_cells->pvariables[PV->W][cellId]=-m_cells->pvariables[PV->W][cellIdAdj];
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellId]=-m_cells->pvariables[PV->RANS_VAR[0]][cellIdAdj];
        }
      }
    }
  }
}

//isothermal Wall
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc1003(ZFSId bcId){
   //cout << "applying bc 1000" << endl;
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  ZFSFloat temp = m_isothermalWallTemperature*PV->TInfinity;

  const ZFSFloat gamma = m_block->m_gamma;

  switch(m_physicalBCMap[bcId]->face){
  case 2:{
    const ZFSInt cellShift = 2*m_noGhostLayers-1;
    for(ZFSId k=start[2]; k<end[2] ; k++) {
      for(ZFSId j=start[1]; j<end[1] ; j++){
        for(ZFSId i=start[0]; i<end[0]; i++){
          const ZFSId cellId=cellIndex(i,j,k); //ghost
          const ZFSId cellIdAdj=cellIndex(i,cellShift-j,k); // field
          const ZFSFloat rhoActive = m_cells->pvariables[PV->RHO][cellIdAdj];
          const ZFSFloat u = (-1.0)*(m_cells->pvariables[PV->U][cellIdAdj]);
          const ZFSFloat v = (-1.0)*(m_cells->pvariables[PV->V][cellIdAdj]);
          const ZFSFloat w = (-1.0)*(m_cells->pvariables[PV->W][cellIdAdj]);
          const ZFSFloat pressure1=m_cells->pvariables[PV->P][cellIdAdj];
          const ZFSFloat rhoWall = pressure1*gamma/temp;
          const ZFSFloat rho = F2*rhoWall-rhoActive;

          m_cells->pvariables[PV->RHO][cellId]=rho;
          m_cells->pvariables[PV->U][cellId]=u;
          m_cells->pvariables[PV->V][cellId]=v;
          m_cells->pvariables[PV->W][cellId]=w;
          m_cells->pvariables[PV->P][cellId]=pressure1;

          if(m_block->m_rans) {
            m_cells->pvariables[PV->RANS_VAR[0]][cellId]=-m_cells->pvariables[PV->RANS_VAR[0]][cellIdAdj];
          }
        }
      }
    }
    break;
  }
  case 3:{
    const ZFSInt cellShift = 2*(m_nCells[1]-1)-2*m_noGhostLayers+1;
    for(ZFSId k=start[2]; k<end[2] ; k++){
      for(ZFSId j=start[1]; j<end[1] ; j++){
        for(ZFSId i=start[0]; i<end[0]; i++){
          const ZFSId cellId=cellIndex(i,j,k); //ghost
          const ZFSId cellIdAdj=cellIndex(i,cellShift-j,k); // field
          const ZFSFloat rhoActive = m_cells->pvariables[PV->RHO][cellIdAdj];
          const ZFSFloat u = (-1.0)*(m_cells->pvariables[PV->U][cellIdAdj]);
          const ZFSFloat v = (-1.0)*(m_cells->pvariables[PV->V][cellIdAdj]);
          const ZFSFloat w = (-1.0)*(m_cells->pvariables[PV->W][cellIdAdj]);
          const ZFSFloat pressure1=m_cells->pvariables[PV->P][cellIdAdj];
          const ZFSFloat rhoWall = pressure1*gamma/temp;
          const ZFSFloat rho = F2*rhoWall-rhoActive;

          m_cells->pvariables[PV->RHO][cellId]=rho;
          m_cells->pvariables[PV->U][cellId]=u;
          m_cells->pvariables[PV->V][cellId]=v;
          m_cells->pvariables[PV->W][cellId]=w;
          m_cells->pvariables[PV->P][cellId]=pressure1;

          if(m_block->m_rans) {
            m_cells->pvariables[PV->RANS_VAR[0]][cellId]=-m_cells->pvariables[PV->RANS_VAR[0]][cellIdAdj];
          }
        }
      }
    }
    break;
  }
  default:
    {
      zfsTerm(1, __CALLING_FUNCTION__, "Face direction not implemented)");
    }
  }
}

/** Moving rigid wall functions
 * 
 * /author Marian Albers
 * /see tfs mvex & bound.f (from Pascal Meysonnat)
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc1004(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  ZFSFloat fdt = F0;
  if (m_block->m_RKStep!=0) {
    fdt = F1/(m_block->m_timeStep*m_block->m_RKalpha[m_block->m_RKStep-1]);
  } else {
    fdt = F1/(m_block->m_timeStep*m_block->m_RKalpha[m_block->m_noRKSteps-1]);
  }

  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
      {
        ZFSFloat gridVel[3] = {F0,F0,F0};
        ZFSFloat mP1=F0, mP2=F0;

        for(ZFSId k=start[2]; k<end[2] ; k++) {
          for(ZFSId j=start[1]; j<end[1]; j++) {
            // determine global point ID for local cell IDs
            ZFSId ijk    = getPointIdFromCell( start[0]+1,j, k );
            ZFSId ipjk   = getPointIdfromPoint( ijk, 1, 0, 0 );
            ZFSId ipjkp  = getPointIdfromPoint( ijk, 1, 0, 1 );
            ZFSId ipjpk  = getPointIdfromPoint( ijk, 1, 1, 0 );
            ZFSId ipjpkp = getPointIdfromPoint( ijk, 1, 1, 1 );

            for(ZFSId dim = 0; dim<nDim; dim++) {
              // arbitrarily moving wall
              if(m_block->m_time > F0) {
                mP2=(              m_coordinates[dim][ipjk]+              m_coordinates[dim][ipjpk]+              m_coordinates[dim][ipjkp]+              m_coordinates[dim][ipjpkp])/F4;
                mP1=(m_block->m_mgOldCoordinates[dim][ipjk]+m_block->m_mgOldCoordinates[dim][ipjpk]+m_block->m_mgOldCoordinates[dim][ipjkp]+m_block->m_mgOldCoordinates[dim][ipjpkp])/F4;
                gridVel[dim]=(mP2-mP1)*fdt;
              }
            }

            ZFSId cellIdG2=cellIndex(start[0],j,k); //ghost
            ZFSId cellIdG1=cellIndex(start[0]+1,j,k); //ghost
            ZFSId cellIdA1=cellIndex(start[0]+2,j,k); // field
            ZFSId cellIdA2=cellIndex(start[0]+3,j,k); // field

            ZFSFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
            ZFSFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
            ZFSFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
            ZFSFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
            ZFSFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

            ZFSFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
            ZFSFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
            ZFSFloat w2 = m_cells->pvariables[PV->W][cellIdA2];


            ZFSFloat pG1 = p1;
            ZFSFloat pG2 = p1;
            ZFSFloat rhoG1 = rho1;
            ZFSFloat rhoG2 = rho1;

            ZFSFloat uG1 = F2*gridVel[0]-u1;
            ZFSFloat vG1 = F2*gridVel[1]-v1;
            ZFSFloat wG1 = F2*gridVel[2]-w1;

            const ZFSFloat uG2 = F2*gridVel[0]-u2;
            const ZFSFloat vG2 = F2*gridVel[1]-v2;
            const ZFSFloat wG2 = F2*gridVel[2]-w2;

            m_cells->pvariables[PV->RHO][cellIdG1]  = rhoG1;
            m_cells->pvariables[PV->U][cellIdG1]=uG1;
            m_cells->pvariables[PV->V][cellIdG1]=vG1;
            m_cells->pvariables[PV->W][cellIdG1]=wG1;
            m_cells->pvariables[PV->P][cellIdG1]= pG1;

            m_cells->pvariables[PV->RHO][cellIdG2]  = rhoG2;
            m_cells->pvariables[PV->U][cellIdG2]= uG2;
            m_cells->pvariables[PV->V][cellIdG2]= vG2;
            m_cells->pvariables[PV->W][cellIdG2]= wG2;
            m_cells->pvariables[PV->P][cellIdG2]= pG2;
          }
        }
        break;
      }
    case 1:
      {
        ZFSFloat gridVel[3] = {F0,F0,F0};
        ZFSFloat mP1=F0, mP2=F0;

        for(ZFSId k=start[2]; k<end[2] ; k++) {
          for(ZFSId j=start[1]; j<end[1]; j++) {
            // determine global point ID for local cell IDs
            ZFSId ijk    = getPointIdFromCell( end[0]-3,j, k );
            ZFSId ipjk   = getPointIdfromPoint( ijk, 1, 0, 0 );
            ZFSId ipjkp  = getPointIdfromPoint( ijk, 1, 0, 1 );
            ZFSId ipjpk  = getPointIdfromPoint( ijk, 1, 1, 0 );
            ZFSId ipjpkp = getPointIdfromPoint( ijk, 1, 1, 1 );

            for(ZFSId dim = 0; dim<nDim; dim++) {
              // arbitrarily moving wall
              if(m_block->m_time > F0) {
                mP2=(              m_coordinates[dim][ipjk]+              m_coordinates[dim][ipjpk]+              m_coordinates[dim][ipjkp]+              m_coordinates[dim][ipjpkp])/F4;
                mP1=(m_block->m_mgOldCoordinates[dim][ipjk]+m_block->m_mgOldCoordinates[dim][ipjpk]+m_block->m_mgOldCoordinates[dim][ipjkp]+m_block->m_mgOldCoordinates[dim][ipjpkp])/F4;
                gridVel[dim]=(mP2-mP1)*fdt;
              }
            }

            ZFSId cellIdG2=cellIndex(end[0]-1,j,k); //ghost
            ZFSId cellIdG1=cellIndex(end[0]-2,j,k); //ghost
            ZFSId cellIdA1=cellIndex(end[0]-3,j,k); // field
            ZFSId cellIdA2=cellIndex(end[0]-4,j,k); // field

            ZFSFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
            ZFSFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
            ZFSFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
            ZFSFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
            ZFSFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

            ZFSFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
            ZFSFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
            ZFSFloat w2 = m_cells->pvariables[PV->W][cellIdA2];


            ZFSFloat pG1 = p1;
            ZFSFloat pG2 = p1;
            ZFSFloat rhoG1 = rho1;
            ZFSFloat rhoG2 = rho1;

            ZFSFloat uG1 = F2*gridVel[0]-u1;
            ZFSFloat vG1 = F2*gridVel[1]-v1;
            ZFSFloat wG1 = F2*gridVel[2]-w1;

            const ZFSFloat uG2 = F2*gridVel[0]-u2;
            const ZFSFloat vG2 = F2*gridVel[1]-v2;
            const ZFSFloat wG2 = F2*gridVel[2]-w2;

            m_cells->pvariables[PV->RHO][cellIdG1]  = rhoG1;
            m_cells->pvariables[PV->U][cellIdG1]=uG1;
            m_cells->pvariables[PV->V][cellIdG1]=vG1;
            m_cells->pvariables[PV->W][cellIdG1]=wG1;
            m_cells->pvariables[PV->P][cellIdG1]= pG1;

            m_cells->pvariables[PV->RHO][cellIdG2]  = rhoG2;
            m_cells->pvariables[PV->U][cellIdG2]= uG2;
            m_cells->pvariables[PV->V][cellIdG2]= vG2;
            m_cells->pvariables[PV->W][cellIdG2]= wG2;
            m_cells->pvariables[PV->P][cellIdG2]= pG2;
          }
        }
        break;
      }
    case 2:
      {
        ZFSFloat gridVel[3] = {F0,F0,F0};
        ZFSFloat mP1=F0, mP2=F0;

        //const ZFSInt cellShift = 2*m_noGhostLayers-1;
        for(ZFSId k=start[2]; k<end[2] ; k++) {
          for(ZFSId i=start[0]; i<end[0]; i++) {
            // determine global point ID for local cell IDs
            ZFSId ijk    = getPointIdFromCell( i, 1, k );
            ZFSId ipjpk  = getPointIdfromPoint( ijk, 1, 1, 0 );
            ZFSId ipjpkp = getPointIdfromPoint( ijk, 1, 1, 1 );
            ZFSId ijpk   = getPointIdfromPoint( ijk, 0, 1, 0 );
            ZFSId ijpkp  = getPointIdfromPoint( ijk, 0, 1, 1 );

            for(ZFSId dim = 0; dim<nDim; dim++) {
              // arbitrarily moving wall
              if(m_block->m_time > F0) {
                mP2=(m_coordinates[dim][ijpk]+m_coordinates[dim][ipjpk]+m_coordinates[dim][ijpkp]+m_coordinates[dim][ipjpkp])*F1B4;
                mP1=(m_block->m_mgOldCoordinates[dim][ijpk]+m_block->m_mgOldCoordinates[dim][ipjpk]+m_block->m_mgOldCoordinates[dim][ijpkp]+m_block->m_mgOldCoordinates[dim][ipjpkp])*F1B4;
                gridVel[dim]=(mP2-mP1)*fdt;
              }
            }

            ZFSId cellIdG2=cellIndex(i,0,k); //ghost
            ZFSId cellIdG1=cellIndex(i,1,k); //ghost
            ZFSId cellIdA1=cellIndex(i,2,k); // field
            ZFSId cellIdA2=cellIndex(i,3,k); // field

            ZFSFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
            ZFSFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
            ZFSFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
            ZFSFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
            ZFSFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

            ZFSFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
            ZFSFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
            ZFSFloat w2 = m_cells->pvariables[PV->W][cellIdA2];


            ZFSFloat pG1 = p1;
            ZFSFloat pG2 = p1;
            ZFSFloat rhoG1 = rho1;
            ZFSFloat rhoG2 = rho1;

            ZFSFloat uG1 = F2*gridVel[0]-u1;
            ZFSFloat vG1 = F2*gridVel[1]-v1;
            ZFSFloat wG1 = F2*gridVel[2]-w1;

            const ZFSFloat uG2 = F2*gridVel[0]-u2;
            const ZFSFloat vG2 = F2*gridVel[1]-v2;
            const ZFSFloat wG2 = F2*gridVel[2]-w2;

            m_cells->pvariables[PV->RHO][cellIdG1]  = rhoG1;
            m_cells->pvariables[PV->U][cellIdG1]=uG1;
            m_cells->pvariables[PV->V][cellIdG1]=vG1;
            m_cells->pvariables[PV->W][cellIdG1]=wG1;
            m_cells->pvariables[PV->P][cellIdG1]= pG1;

            m_cells->pvariables[PV->RHO][cellIdG2]  = rhoG2;
            m_cells->pvariables[PV->U][cellIdG2]= uG2;
            m_cells->pvariables[PV->V][cellIdG2]= vG2;
            m_cells->pvariables[PV->W][cellIdG2]= wG2;
            m_cells->pvariables[PV->P][cellIdG2]= pG2;
          }
        }
        break;
      }
    case 3:  // indenting wall
      {
        ZFSFloat gridVel[3] = {F0,F0,F0};
        ZFSFloat mP1=F0, mP2=F0;

        for(ZFSId k=start[2]; k<end[2] ; k++){
          for(ZFSId i=start[0]; i<end[0]; i++) {
            // determine global point ID for local cell IDs
            ZFSId ijk    = getPointIdFromCell( i, end[1]-3, k );
            ZFSId ipjpk  = getPointIdfromPoint( ijk, 1, 1, 0 );
            ZFSId ipjpkp = getPointIdfromPoint( ijk, 1, 1, 1 );
            ZFSId ijpk   = getPointIdfromPoint( ijk, 0, 1, 0 );
            ZFSId ijpkp  = getPointIdfromPoint( ijk, 0, 1, 1 );


            for(ZFSId dim = 0; dim<nDim; dim++) {
              // arbitrarily moving wall
              if(m_block->m_time > F0) {
                mP2=(m_coordinates[dim][ijpk]+m_coordinates[dim][ipjpk]+m_coordinates[dim][ijpkp]+m_coordinates[dim][ipjpkp])/F4;
                mP1=(m_block->m_mgOldCoordinates[dim][ijpk]+m_block->m_mgOldCoordinates[dim][ipjpk]+m_block->m_mgOldCoordinates[dim][ijpkp]+m_block->m_mgOldCoordinates[dim][ipjpkp])/F4;
                gridVel[dim]=(mP2-mP1)*fdt;
              }
            }

            ZFSId cellIdG2=cellIndex(i,end[1]-1,k); //ghost
            ZFSId cellIdG1=cellIndex(i,end[1]-2,k); //ghost
            ZFSId cellIdA1=cellIndex(i,end[1]-3,k); // field
            ZFSId cellIdA2=cellIndex(i,end[1]-4,k); // field

            ZFSFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
            ZFSFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
            ZFSFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
            ZFSFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
            ZFSFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

            ZFSFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
            ZFSFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
            ZFSFloat w2 = m_cells->pvariables[PV->W][cellIdA2];


            ZFSFloat pG1 = p1;
            ZFSFloat pG2 = p1;
            ZFSFloat rhoG1 = rho1;
            ZFSFloat rhoG2 = rho1;

            ZFSFloat uG1 = F2*gridVel[0]-u1;
            ZFSFloat vG1 = F2*gridVel[1]-v1;
            ZFSFloat wG1 = F2*gridVel[2]-w1;

            const ZFSFloat uG2 = F2*gridVel[0]-u2;
            const ZFSFloat vG2 = F2*gridVel[1]-v2;
            const ZFSFloat wG2 = F2*gridVel[2]-w2;

            m_cells->pvariables[PV->RHO][cellIdG1]  = rhoG1;
            m_cells->pvariables[PV->U][cellIdG1]=uG1;
            m_cells->pvariables[PV->V][cellIdG1]=vG1;
            m_cells->pvariables[PV->W][cellIdG1]=wG1;
            m_cells->pvariables[PV->P][cellIdG1]= pG1;

            m_cells->pvariables[PV->RHO][cellIdG2]  = rhoG2;
            m_cells->pvariables[PV->U][cellIdG2]= uG2;
            m_cells->pvariables[PV->V][cellIdG2]= vG2;
            m_cells->pvariables[PV->W][cellIdG2]= wG2;
            m_cells->pvariables[PV->P][cellIdG2]= pG2;
          }
        }
        break;
      }
    default:
      {
        cout << "bc1004: face not implemented" << endl;
      }
    }
  //exit(1);
}

/** Moving rigid wall - isothermal
 *
 * /author Marian Albers
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc1006(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  const ZFSFloat temp = m_isothermalWallTemperature*PV->TInfinity;
  const ZFSFloat gamma = m_block->m_gamma;

  ZFSFloat fdt = F0;
  if (m_block->m_RKStep!=0) {
    fdt = F1/(m_block->m_timeStep*m_block->m_RKalpha[m_block->m_RKStep-1]);
  } else {
    fdt = F1/(m_block->m_timeStep*m_block->m_RKalpha[m_block->m_noRKSteps-1]);
  }

  switch(m_physicalBCMap[bcId]->face)
    {
    case 2:
      {
        ZFSFloat gridVel[3] = {F0,F0,F0};
        ZFSFloat mP1=F0, mP2=F0;

        //const ZFSInt cellShift = 2*m_noGhostLayers-1;
        for(ZFSId k=start[2]; k<end[2] ; k++) {
          for(ZFSId i=start[0]; i<end[0]; i++) {
            // determine global point ID for local cell IDs
            ZFSId ijk    = getPointIdFromCell( i, 1, k );
            ZFSId ipjpk  = getPointIdfromPoint( ijk, 1, 1, 0 );
            ZFSId ipjpkp = getPointIdfromPoint( ijk, 1, 1, 1 );
            ZFSId ijpk   = getPointIdfromPoint( ijk, 0, 1, 0 );
            ZFSId ijpkp  = getPointIdfromPoint( ijk, 0, 1, 1 );

            for(ZFSId dim = 0; dim<nDim; dim++) {
              // arbitrarily moving wall
              if(m_block->m_time > F0) {
                mP2=(m_coordinates[dim][ijpk]+m_coordinates[dim][ipjpk]+m_coordinates[dim][ijpkp]+m_coordinates[dim][ipjpkp])*F1B4;
                mP1=(m_block->m_mgOldCoordinates[dim][ijpk]+m_block->m_mgOldCoordinates[dim][ipjpk]+m_block->m_mgOldCoordinates[dim][ijpkp]+m_block->m_mgOldCoordinates[dim][ipjpkp])*F1B4;
                gridVel[dim]=(mP2-mP1)*fdt;
              }
            }

            const ZFSId cellIdG2=cellIndex(i,0,k); //ghost
            const ZFSId cellIdG1=cellIndex(i,1,k); //ghost
            const ZFSId cellIdA1=cellIndex(i,2,k); // field
            const ZFSId cellIdA2=cellIndex(i,3,k); // field

            const ZFSFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
            const ZFSFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
            const ZFSFloat rho2 = m_cells->pvariables[PV->RHO][cellIdA2];
            const ZFSFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
            const ZFSFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
            const ZFSFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

            const ZFSFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
            const ZFSFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
            const ZFSFloat w2 = m_cells->pvariables[PV->W][cellIdA2];

            const ZFSFloat rhoWall = p1*gamma/temp;
            const ZFSFloat rhoG1 = F2*rhoWall-rho1;
            const ZFSFloat rhoG2 = F2*rhoWall-rho2;

            const ZFSFloat pG1 = p1;
            const ZFSFloat pG2 = p1;

            const ZFSFloat uG1 = F2*gridVel[0]-u1;
            const ZFSFloat vG1 = F2*gridVel[1]-v1;
            const ZFSFloat wG1 = F2*gridVel[2]-w1;

            const ZFSFloat uG2 = F2*gridVel[0]-u2;
            const ZFSFloat vG2 = F2*gridVel[1]-v2;
            const ZFSFloat wG2 = F2*gridVel[2]-w2;

            m_cells->pvariables[PV->RHO][cellIdG1]  = rhoG1;
            m_cells->pvariables[PV->U][cellIdG1]=uG1;
            m_cells->pvariables[PV->V][cellIdG1]=vG1;
            m_cells->pvariables[PV->W][cellIdG1]=wG1;
            m_cells->pvariables[PV->P][cellIdG1]= pG1;

            m_cells->pvariables[PV->RHO][cellIdG2]  = rhoG2;
            m_cells->pvariables[PV->U][cellIdG2]= uG2;
            m_cells->pvariables[PV->V][cellIdG2]= vG2;
            m_cells->pvariables[PV->W][cellIdG2]= wG2;
            m_cells->pvariables[PV->P][cellIdG2]= pG2;
          }
        }
        break;
      }
    default:
      {
        cout << "bc1006: face not implemented" << endl;
      }
    }
}


/** Oscillating wall
 *
 * /author Marian Albers
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc1007(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  const ZFSFloat t = m_block->m_time + m_block->m_timeStep*m_block->m_RKalpha[m_block->m_RKStep];

  switch(m_physicalBCMap[bcId]->face)
    {
    case 2:
      {
        ZFSFloat gridVel[3] = {F0,F0,F0};

        for(ZFSId k=start[2]; k<end[2] ; k++) {
          for(ZFSId i=start[0]; i<end[0]; i++) {
            ZFSId cellIdG2=cellIndex(i,0,k); //ghost
            ZFSId cellIdG1=cellIndex(i,1,k); //ghost
            ZFSId cellIdA1=cellIndex(i,2,k); // field
            ZFSId cellIdA2=cellIndex(i,3,k); // field

            const ZFSFloat transitionLength = m_block->m_waveEndTransition - m_block->m_waveBeginTransition;
            const ZFSFloat xInit = m_cells->coordinates[0][cellIdG1];
            ZFSFloat transitionFactor = F0;
            if(xInit <= m_block->m_waveBeginTransition) {
              transitionFactor = F0;
            } else if (xInit > m_block->m_waveBeginTransition && xInit< m_block->m_waveEndTransition) {
              transitionFactor = (1-cos((xInit-m_block->m_waveBeginTransition)/transitionLength*PI))*F1B2;
            } else {
              transitionFactor = F1;
            }

            gridVel[2] = m_block->m_waveAmplitude*sin(F2*PI/m_block->m_waveTime*t)*transitionFactor;

            ZFSFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
            ZFSFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
            ZFSFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
            ZFSFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
            ZFSFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

            ZFSFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
            ZFSFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
            ZFSFloat w2 = m_cells->pvariables[PV->W][cellIdA2];

            ZFSFloat pG1 = p1;
            ZFSFloat pG2 = p1;
            ZFSFloat rhoG1 = rho1;
            ZFSFloat rhoG2 = rho1;

            ZFSFloat uG1 = F2*gridVel[0]-u1;
            ZFSFloat vG1 = F2*gridVel[1]-v1;
            ZFSFloat wG1 = F2*gridVel[2]-w1;

            const ZFSFloat uG2 = F2*gridVel[0]-u2;
            const ZFSFloat vG2 = F2*gridVel[1]-v2;
            const ZFSFloat wG2 = F2*gridVel[2]-w2;

            m_cells->pvariables[PV->RHO][cellIdG1]  = rhoG1;
            m_cells->pvariables[PV->U][cellIdG1]=uG1;
            m_cells->pvariables[PV->V][cellIdG1]=vG1;
            m_cells->pvariables[PV->W][cellIdG1]=wG1;
            m_cells->pvariables[PV->P][cellIdG1]= pG1;

            m_cells->pvariables[PV->RHO][cellIdG2]  = rhoG2;
            m_cells->pvariables[PV->U][cellIdG2]= uG2;
            m_cells->pvariables[PV->V][cellIdG2]= vG2;
            m_cells->pvariables[PV->W][cellIdG2]= wG2;
            m_cells->pvariables[PV->P][cellIdG2]= pG2;
          }
        }
        break;
      }
    default:
      {
        cout << "bc1007: face not implemented" << endl;
      }
    }
}


/** Subsonic rotational inflow
 *
 *  rho=rho_inf, u=u_inf, v=v_inf, w=w_inf, dp/dn=0
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2014(ZFSId bcId)
{
  //implemented for i-direction only for the moment
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
    {
      ZFSId cellId=-1;
      ZFSId cellIdadj=-1;
      for(ZFSId k=start[2]; k<end[2] ; k++) {
        for(ZFSId j=start[1]; j<end[1] ; j++) {
          for(ZFSId i=start[0]; i<end[0]; i++) {
            cellId=cellIndex(m_noGhostLayers-1-i,j,k);
            cellIdadj=cellIndex(m_noGhostLayers-i,j,k);

            ZFSFloat y= m_cells->coordinates[1][cellId];
            ZFSFloat z= m_cells->coordinates[2][cellId];
            ZFSFloat phi=atan2(y,z);
            ZFSFloat r=sqrt(POW2(y)+POW2(z));
            ZFSFloat rmax = 10.0;

            m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
            m_cells->pvariables[PV->U][cellId]= PV->UInfinity;
            m_cells->pvariables[PV->V][cellId]=-(r/rmax)*cos(phi)*PV->UInfinity;
            m_cells->pvariables[PV->W][cellId]=(r/rmax)*sin(phi)*PV->UInfinity;
            m_cells->pvariables[PV->P][cellId]= m_cells->pvariables[PV->P][cellIdadj];
          }
        }
      }
      break;
    }
    case 2:
      {
        for(ZFSId k=start[2]; k<end[2] ; k++) {
          for(ZFSId j=start[1]; j<end[1] ; j++) {
            for(ZFSId i=start[0]; i<end[0]; i++) {
              const ZFSId cellId=cellIndex(i,m_noGhostLayers-j-1,k);//ghost
              const ZFSId cellIdadj=cellIndex(i,m_noGhostLayers-j,k);//field
              m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
              m_cells->pvariables[PV->U][cellId]= PV->UInfinity;
              m_cells->pvariables[PV->V][cellId]= PV->VInfinity;
              m_cells->pvariables[PV->W][cellId]= PV->WInfinity;

              m_cells->pvariables[PV->P][cellId]= m_cells->pvariables[PV->P][cellIdadj];
            }
          }
        }
        break;
      }
    default:
      {
        zfsTerm(1, __CALLING_FUNCTION__, "Face direction not implemented)");
      }
    }

}



/** Subsonic Inflow
 *
 *  rho=rho_inf, u=u_inf, v=v_inf, w=w_inf, dp/dn=0
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2001(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  ZFSInt face = m_physicalBCMap[bcId]->face;
  const ZFSInt n = face%2;
  const ZFSInt dim = face/2;
  ZFSInt offset[3] = {0,0,0};
  ZFSInt inc[3] = {1,1,1};
  ZFSInt adjInc[3] = {0,0,0};

  //on face 0,2,4 set offset and
  //reverse the incrementor
  if(n==0) {
    offset[dim] = m_noGhostLayers-1;
    inc[dim] = -1;
    adjInc[dim] = 1;
  } else {
    inc[dim] = 1;
    adjInc[dim] = -1;
  }

  for(ZFSId k=start[2]; k<end[2] ; k++) {
    for(ZFSId j=start[1]; j<end[1] ; j++) {
      for(ZFSId i=start[0]; i<end[0]; i++) {
        const ZFSId ii = inc[0]*(i-start[0])+start[0]+offset[0];
        const ZFSId jj = inc[1]*(j-start[1])+start[1]+offset[1];
        const ZFSId kk = inc[2]*(k-start[2])+start[2]+offset[2];
        const ZFSId cellId=cellIndex(ii,jj,kk);
        const ZFSId cellIdadj=cellIndex(ii+adjInc[0],jj+adjInc[1],kk+adjInc[2]);

        m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
        m_cells->pvariables[PV->U][cellId]=PV->UInfinity;
        m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
        m_cells->pvariables[PV->W][cellId]=PV->WInfinity;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellId]= CV->ransInfinity[0]/CV->rhoInfinity;
        }

        //extrapolate pressure from the field
        m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
      }
    }
  }
}



/** Subsonic Inflow with u=(y/delta)^(1/7)
 *
 *  rho=rho_inf, u=u_inf, v=v_inf, w=w_inf, dp/dn=0
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2099(ZFSId bcId)
{
  //implemented for i-direction only for the moment
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  ZFSFloat y=F0;
  switch(m_physicalBCMap[bcId]->face)
  {
  case 0:
  {
    for(ZFSId k=start[2]; k<end[2] ; k++) {
      for(ZFSId j=start[1]; j<end[1] ; j++) {
        for(ZFSId i=start[0]; i<end[0]; i++) {
          const ZFSId cellId=cellIndex(m_noGhostLayers-1-i,j,k);
          const ZFSId cellIdadj=cellIndex(m_noGhostLayers-i,j,k);
          m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
          m_cells->pvariables[PV->U][cellId]=PV->UInfinity;
          y=m_cells->coordinates[1][cellId];
          if(y<1.0) {
            m_cells->pvariables[PV->U][cellId]=PV->UInfinity*pow(y,(1.0/7.0));
          }
          m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
          m_cells->pvariables[PV->W][cellId]=PV->WInfinity;
          m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
        }
      }
    }
    break;
  }
  default:
  {
    zfsTerm(1, __CALLING_FUNCTION__, "Face direction not implemented)");
  }
  }
}

/** Supersonic Inflow
 *
 *  rho=rho_inf, u=u_inf, v=v_inf, w=w_inf, p=p_inf
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2002(ZFSId bcId)
{
  TRACE();
  //implemented for i-direction only for the moment
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  ZFSId cellId=-1;
  for(ZFSId k=start[2]; k<end[2] ; k++) {
    for(ZFSId j=start[1]; j<end[1] ; j++) {
      for(ZFSId i=start[0]; i<end[0]; i++) {
        cellId=cellIndex(i,j,k);
        m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
        m_cells->pvariables[PV->U][cellId]=PV->UInfinity;
        m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
        m_cells->pvariables[PV->W][cellId]=PV->WInfinity;
        m_cells->pvariables[PV->P][cellId]=PV->PInfinity;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellId]= CV->ransInfinity[0]/CV->rhoInfinity;
        }
      }
    }
  }
}

/** Supersonic outflow
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2005(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  ZFSInt face = m_physicalBCMap[bcId]->face;
  const ZFSInt n = face%2;
  const ZFSInt dim = face/2;
  ZFSInt offset[3] = {0,0,0};
  ZFSInt inc[3] = {1,1,1};
  ZFSInt adjInc[3] = {0,0,0};

  //on face 0,2,4 set offset and
  //reverse the incrementor
  if(n==0) {
    offset[dim] = m_noGhostLayers-1;
    inc[dim] = -1;
    adjInc[dim] = 1;
  } else {
    inc[dim] = 1;
    adjInc[dim] = -1;
  }

  for(ZFSId k=start[2]; k<end[2] ; k++) {
    for(ZFSId j=start[1]; j<end[1] ; j++) {
      for(ZFSId i=start[0]; i<end[0]; i++) {
        const ZFSId ii = inc[0]*(i-start[0])+start[0]+offset[0];
        const ZFSId jj = inc[1]*(j-start[1])+start[1]+offset[1];
        const ZFSId kk = inc[2]*(k-start[2])+start[2]+offset[2];
        const ZFSId cellId=cellIndex(ii,jj,kk);
        const ZFSId cellIdadj=cellIndex(ii+adjInc[0],jj+adjInc[1],kk+adjInc[2]);

        for(ZFSId var=0; var<PV->noVariables; var++) {
          m_cells->pvariables[var][cellId]=m_cells->pvariables[var][cellIdadj];
        }
      }
    }
  }
}

/* Subsonic in/outflow simple
 * Simple extrapolation/prescription depending
 * on in/outflow condition
 *
 * /author Marian Albers
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2003(ZFSId bcId)
{
  const ZFSInt IJK[nDim] = {1,m_nCells[2],m_nCells[1]*m_nCells[2]};
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  //Here we find out the normal direction of the
  //boundary and the two tangential directions.
  //This way we can make a general formulation of
  //the boundary condition
  const ZFSInt face = m_physicalBCMap[bcId]->face;
  const ZFSInt normalDir = face/2;
  const ZFSInt firstTangentialDir = (normalDir+1)%nDim;
  const ZFSInt secondTangentialDir = (normalDir+2)%nDim;
  const ZFSInt normalDirStart = start[normalDir];
  const ZFSInt firstTangentialStart = start[firstTangentialDir];
  const ZFSInt firstTangentialEnd = end[firstTangentialDir];
  const ZFSInt secondTangentialStart = start[secondTangentialDir];
  const ZFSInt secondTangentialEnd = end[secondTangentialDir];
  const ZFSInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const ZFSInt n = (face%2)*2-1; //-1,+1
  const ZFSInt g1 = normalDirStart + (ZFSId)(0.5-(0.5*(ZFSFloat)n)); //+1,0
  const ZFSInt g2 = normalDirStart + (ZFSId)(0.5+(0.5*(ZFSFloat)n)); //0,+1
  const ZFSInt a1 = normalDirStart + (ZFSId)(0.5-(1.5*(ZFSFloat)n)); //+2,-1
  const ZFSInt a2 = normalDirStart + (ZFSId)(0.5-(2.5*(ZFSFloat)n)); //+3,-2

  for(ZFSId t1=firstTangentialStart; t1<firstTangentialEnd; t1++) {
    for(ZFSId t2=secondTangentialStart; t2<secondTangentialEnd ; t2++) {
      const ZFSId cellIdG1 = g1*inc[0] + t1*inc[1] + t2*inc[2];
      const ZFSId cellIdG2 = g2*inc[0] + t1*inc[1] + t2*inc[2];
      const ZFSId cellIdA1 = a1*inc[0] + t1*inc[1] + t2*inc[2];
      const ZFSId cellIdA2 = a2*inc[0] + t1*inc[1] + t2*inc[2];

      const ZFSFloat dxidx = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+0];
      const ZFSFloat dxidy = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+1];
      const ZFSFloat dxidz = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+2];
      //leaving domain of integration in positive coordinate direction,
      //therefore multiply with positive F1
      const ZFSFloat gradxi = n*F1 / sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz);

      const ZFSFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
      const ZFSFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
      const ZFSFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
      const ZFSFloat wInner = m_cells->pvariables[PV->W][cellIdA1];
      const ZFSFloat pInner = m_cells->pvariables[PV->P][cellIdA1];

      const ZFSFloat maContravariant = (dxidx*uInner +
                                        dxidy*vInner +
                                        dxidz*wInner - m_cells->dxt[normalDir][cellIdA1])*gradxi;


      if(maContravariant < F0) {
        //inflow
        const ZFSFloat p = pInner;
        const ZFSFloat rho = CV->rhoInfinity;

        m_cells->pvariables[PV->RHO][cellIdG1] = rho;
        m_cells->pvariables[PV->U][cellIdG1] = PV->UInfinity;
        m_cells->pvariables[PV->V][cellIdG1] = PV->VInfinity;
        m_cells->pvariables[PV->W][cellIdG1] = PV->WInfinity;
        m_cells->pvariables[PV->P][cellIdG1] = p;

        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1]= CV->ransInfinity[0]/CV->rhoInfinity;
        }
      } else {
        //outflow
        const ZFSFloat p =PV->PInfinity;
        const ZFSFloat rho = rhoInner;

        m_cells->pvariables[PV->RHO][cellIdG1] = rho;
        m_cells->pvariables[PV->U][cellIdG1] = uInner;
        m_cells->pvariables[PV->V][cellIdG1] = vInner;
        m_cells->pvariables[PV->W][cellIdG1] = wInner;
        m_cells->pvariables[PV->P][cellIdG1] = p;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1]= (F2*m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1]
                                                           -m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2]);
        }
      }

      //extrapolate into second ghost cell
      for(ZFSId var=0; var<PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdG2] = F2*m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
      }
    }
  }
}


/** Characteristic boundary condition, detects
 *  whether it is in- or outflow and prescribes
 *  appropriately (see tfs)
 *  /author Marian Albers
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2004(ZFSId bcId)
{
  const ZFSInt IJK[nDim] = {1,m_nCells[2],m_nCells[1]*m_nCells[2]};
  const ZFSFloat gamma= m_block->m_gamma;
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  //Here we find out the normal direction of the
  //boundary and the two tangential directions.
  //This way we can make a general formulation of
  //the boundary condition
  const ZFSInt face = m_physicalBCMap[bcId]->face;
  const ZFSInt normalDir = face/2;
  const ZFSInt firstTangentialDir = (normalDir+1)%nDim;
  const ZFSInt secondTangentialDir = (normalDir+2)%nDim;
  const ZFSInt normalDirStart = start[normalDir];
  const ZFSInt firstTangentialStart = start[firstTangentialDir];
  const ZFSInt firstTangentialEnd = end[firstTangentialDir];
  const ZFSInt secondTangentialStart = start[secondTangentialDir];
  const ZFSInt secondTangentialEnd = end[secondTangentialDir];
  const ZFSInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const ZFSInt n = (face%2)*2-1; //-1,+1
  const ZFSInt g1 = normalDirStart + (ZFSId)(0.5-(0.5*(ZFSFloat)n)); //+1,0
  const ZFSInt g2 = normalDirStart + (ZFSId)(0.5+(0.5*(ZFSFloat)n)); //0,+1
  const ZFSInt a1 = normalDirStart + (ZFSId)(0.5-(1.5*(ZFSFloat)n)); //+2,-1
  const ZFSInt a2 = normalDirStart + (ZFSId)(0.5-(2.5*(ZFSFloat)n)); //+3,-2

  for(ZFSId t1=firstTangentialStart; t1<firstTangentialEnd; t1++) {
    for(ZFSId t2=secondTangentialStart; t2<secondTangentialEnd ; t2++) {
      const ZFSId cellIdG1 = g1*inc[0] + t1*inc[1] + t2*inc[2];
      const ZFSId cellIdG2 = g2*inc[0] + t1*inc[1] + t2*inc[2];
      const ZFSId cellIdA1 = a1*inc[0] + t1*inc[1] + t2*inc[2];
      const ZFSId cellIdA2 = a2*inc[0] + t1*inc[1] + t2*inc[2];

      const ZFSFloat dxidx = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+0];
      const ZFSFloat dxidy = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+1];
      const ZFSFloat dxidz = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+2];
      //multiply with n, so it will be -1 or +1 depending if we enter
      //or leave the domain of integration in positive direction
      const ZFSFloat gradxi = n*F1 / sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz);

      const ZFSFloat dxHelp = dxidx*gradxi;
      const ZFSFloat dyHelp = dxidy*gradxi;
      const ZFSFloat dzHelp = dxidz*gradxi;


      const ZFSFloat cBC = sqrt(gamma*m_cells->pvariables[PV->P][cellIdG1]/m_cells->pvariables[PV->RHO][cellIdG1]);
      const ZFSFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];

      const ZFSFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
      const ZFSFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
      const ZFSFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
      const ZFSFloat wInner = m_cells->pvariables[PV->W][cellIdA1];
      const ZFSFloat pInner = m_cells->pvariables[PV->P][cellIdA1];

      const ZFSFloat maContravariant = (dxidx*uInner +
                                        dxidy*vInner +
                                        dxidz*wInner - m_cells->dxt[normalDir][cellIdA1])*gradxi;

      if(maContravariant < F0) {
        //inflow
        const ZFSFloat p =F1B2*(pInner+PV->PInfinity+rhoBC*cBC*(dxHelp*(uInner - PV->UInfinity)+
                                                                dyHelp*(vInner - PV->VInfinity)+
                                                                dzHelp*(wInner - PV->WInfinity)));

        const ZFSFloat rho = CV->rhoInfinity + (p-PV->PInfinity)/POW2(cBC);
        const ZFSFloat help = (p-PV->PInfinity)/(rhoBC*cBC);

        m_cells->pvariables[PV->RHO][cellIdG1] = rho;
        m_cells->pvariables[PV->U][cellIdG1] = (PV->UInfinity + help*dxHelp);
        m_cells->pvariables[PV->V][cellIdG1] = (PV->VInfinity + help*dyHelp);
        m_cells->pvariables[PV->W][cellIdG1] = (PV->WInfinity + help*dzHelp);
        m_cells->pvariables[PV->P][cellIdG1] = p;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1]= PV->ransInfinity[0];
        }
      } else {
        //outflow
        const ZFSFloat p =PV->PInfinity;
        const ZFSFloat rho = rhoInner + (p-pInner)/POW2(cBC);
        const ZFSFloat help = (p-pInner)/(rhoBC*cBC);

        m_cells->pvariables[PV->RHO][cellIdG1] = rho;
        m_cells->pvariables[PV->U][cellIdG1] = (uInner - help*dxHelp);
        m_cells->pvariables[PV->V][cellIdG1] = (vInner - help*dyHelp);
        m_cells->pvariables[PV->W][cellIdG1] = (wInner - help*dzHelp);
        m_cells->pvariables[PV->P][cellIdG1] = p;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1]= (F2*m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1]
                                                           -m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2]);
        }
      }

      //extrapolate into second ghost cell
      for(ZFSId var=0; var<PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdG2] = F2*m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
      }
    }
  }
}


/** outflow condition after shock
 *
 *  Uses characteristic subsonic inflow (only correct if the machnumber normal to the Boundary is smaller than 1)
 *  Loosen/Meysonnat 2017
 */
template <ZFSBool isRans>   //junoh
void ZFSStrctrdBndryCnd3D<isRans>::bc2222(ZFSId bcId)
{
  const ZFSInt IJK[nDim] = {1,m_nCells[2],m_nCells[1]*m_nCells[2]};
  const ZFSFloat gamma= m_block->m_gamma;
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  //Here we find out the normal direction of the
  //boundary and the two tangential directions.
  //This way we can make a general formulation of
  //the boundary condition
  const ZFSInt face = m_physicalBCMap[bcId]->face;
  const ZFSInt normalDir = face/2;
  const ZFSInt firstTangentialDir = (normalDir+1)%nDim;
  const ZFSInt secondTangentialDir = (normalDir+2)%nDim;
  const ZFSInt normalDirStart = start[normalDir];
  const ZFSInt firstTangentialStart = start[firstTangentialDir];
  const ZFSInt firstTangentialEnd = end[firstTangentialDir];
  const ZFSInt secondTangentialStart = start[secondTangentialDir];
  const ZFSInt secondTangentialEnd = end[secondTangentialDir];
  const ZFSInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const ZFSInt n = (face%2)*2-1; //-1,+1
  const ZFSInt g1 = normalDirStart + (ZFSId)(0.5-(0.5*(ZFSFloat)n)); //+1,0
  const ZFSInt g2 = normalDirStart + (ZFSId)(0.5+(0.5*(ZFSFloat)n)); //0,+1
  const ZFSInt a1 = normalDirStart + (ZFSId)(0.5-(1.5*(ZFSFloat)n)); //+2,-1
  const ZFSInt a2 = normalDirStart + (ZFSId)(0.5-(2.5*(ZFSFloat)n)); //+3,-2

  for(ZFSId t1=firstTangentialStart; t1<firstTangentialEnd; t1++) {
    for(ZFSId t2=secondTangentialStart; t2<secondTangentialEnd ; t2++) {
      const ZFSId cellIdG1 = g1*inc[0] + t1*inc[1] + t2*inc[2];
      const ZFSId cellIdG2 = g2*inc[0] + t1*inc[1] + t2*inc[2];
      const ZFSId cellIdA1 = a1*inc[0] + t1*inc[1] + t2*inc[2];
      const ZFSId cellIdA2 = a2*inc[0] + t1*inc[1] + t2*inc[2];

      const ZFSFloat dxidx = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+0];
      const ZFSFloat dxidy = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+1];
      const ZFSFloat dxidz = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+2];
      //multiply with n, so it will be -1 or +1 depending if we enter
      //or leave the domain of integration in positive direction
      const ZFSFloat gradxi = n*F1 / sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz);

      const ZFSFloat dxHelp = dxidx*gradxi;
      const ZFSFloat dyHelp = dxidy*gradxi;
      const ZFSFloat dzHelp = dxidz*gradxi;


      const ZFSFloat cBC = sqrt(gamma*m_cells->pvariables[PV->P][cellIdG1]/m_cells->pvariables[PV->RHO][cellIdG1]);
      const ZFSFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];

      const ZFSFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
      const ZFSFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
      const ZFSFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
      const ZFSFloat wInner = m_cells->pvariables[PV->W][cellIdA1];
      const ZFSFloat pInner = m_cells->pvariables[PV->P][cellIdA1];

      const ZFSFloat maContravariant = (dxidx*uInner +
                                        dxidy*vInner +
                                        dxidz*wInner - m_cells->dxt[normalDir][cellIdA1])*gradxi;

      if(maContravariant < F0) {
        //inflow
        const ZFSFloat p =F1B2*(pInner+PV->PInfinity+rhoBC*cBC*(dxHelp*(uInner - PV->UInfinity)+
                                                                dyHelp*(vInner - PV->VInfinity)+
                                                                dzHelp*(wInner - PV->WInfinity)));

        // const ZFSFloat rho = CV->rhoInfinity + (p-PV->PInfinity)/POW2(cBC);
        // const ZFSFloat help = (p-PV->PInfinity)/(rhoBC*cBC);

        m_cells->pvariables[PV->RHO][cellIdG1] = m_cells->fq[FQ->AVG_RHO][cellIdG1];
        m_cells->pvariables[PV->U][cellIdG1] = m_cells->fq[FQ->AVG_U][cellIdG1];
        m_cells->pvariables[PV->V][cellIdG1] = m_cells->fq[FQ->AVG_V][cellIdG1];
        m_cells->pvariables[PV->W][cellIdG1] = m_cells->fq[FQ->AVG_W][cellIdG1];
        m_cells->pvariables[PV->P][cellIdG1] = p;       
        
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1]= PV->ransInfinity[0];
	  // m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1]= m_cells->fq[FQ->RECONST_NUTILDE][cellIdG1];
	  //LES-RANS inflow
        }
      } else {
      //outflow
      const ZFSFloat p = m_cells->fq[FQ->AVG_P][cellIdG1];
      const ZFSFloat rho = rhoInner + (p-pInner)/POW2(cBC);
      const ZFSFloat help = (p-pInner)/(rhoBC*cBC);

      m_cells->pvariables[PV->RHO][cellIdG1] = rho;
      m_cells->pvariables[PV->U][cellIdG1] = (uInner - help*dxHelp);
      m_cells->pvariables[PV->V][cellIdG1] = (vInner - help*dyHelp);
      m_cells->pvariables[PV->W][cellIdG1] = (wInner - help*dzHelp);
      m_cells->pvariables[PV->P][cellIdG1] = p;
      if(isRans) {
	m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1]= (F2*m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1]
							 -m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2]);
      }
      }

      //extrapolate into second ghost cell
      for(ZFSId var=0; var<PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdG2] = F2*m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
      }
      if(isRans) { //junoh
	m_cells->pvariables[PV->RANS_VAR[0]][cellIdG2]= (F2*m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1]-m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1]);
      }
    }
  }
}


template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2009(ZFSId bcId)
{
  const ZFSFloat gamma=m_block->m_gamma;
  const ZFSFloat gammaMinusOne=gamma-1.0;
  const ZFSFloat gammaPlusOne=gamma+1.0;
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  ZFSFloat rho2=F0, u2=F0, v2tmp=F0, p2=F0;     //Conditions behind the shock
  ZFSFloat u1n=F0, v1n=F0, u2n=F0, v2n=F0  ; // Velocities normal to the shock upstream and downstream of the shock
  ZFSFloat rho21=F0,p21=F0;                  //quotient of the values upstream and downstream of the shock
  ZFSFloat beta=F0, M_mean=m_block->m_Ma;       
  
  //pressure and density ratio over the shock
  rho21=gammaPlusOne*POW2(M_mean*sin(m_sigma))/(gammaMinusOne*POW2(M_mean*sin(m_sigma))+2);
  p21=1.0+2*gamma/gammaPlusOne*(POW2(M_mean*sin(m_sigma))-1);
    

  //normal and parallel velocities upstream of the shock
  u1n=M_mean*sqrt(PV->TInfinity)*sin(m_sigma);
  v1n=M_mean*sqrt(PV->TInfinity)*cos(m_sigma);
  //normal and parallel velocity downstream of the shock
  u2n=u1n/rho21;
  v2n=v1n;

  beta=m_sigma-atan(u2n/v2n);  //redirection angle of the flow

  u2 = sqrt(POW2(u2n)+POW2(v2n))*cos(beta);
  v2tmp = sqrt(POW2(u2n)+POW2(v2n))*sin(beta);
 
  rho2=rho21*CV->rhoInfinity;
  p2=p21*PV->PInfinity;

 switch(m_physicalBCMap[bcId]->face){
 case 3:{
   ZFSId cellId=-1;
   ZFSId cellIdadj=-1;
   ZFSId pIJK=0, pIJPK=0, pIJKP=0, pIJPKP=0;
   //fully new bc
   ZFSId j= start[1];
   ZFSFloat pBC=F0, rho=F0, u=F0, v=F0, w=F0; 
   ZFSFloat drho=F0, du=F0, dv=F0, dw=F0, dp=F0;
   ZFSFloat yBC=F0;
   ZFSFloat v2=F0;
   //pBC=PV->PInfinity;
   ZFSFloat pInner=F0, c02=F0, distance=F0;
      
   //Change sign of the v velocity
   v2=-v2tmp;
      
   for(ZFSId k=start[2]; k<end[2] ; k++){
     for(ZFSId i=start[0]; i<end[0] ; i++){
       //cellId=cellIndex(i,j,k);
       cellIdadj=cellIndex(i,j-1,k);
       //to determine the face coordinates!!!!!!
       pIJK=getPointIdFromCell(i,j,k);
       pIJPK=getPointIdfromPoint(pIJK,1,0,0);
       pIJKP=getPointIdfromPoint(pIJK,0,0,1);
       pIJPKP=getPointIdfromPoint(pIJK,1,0,1);
	  
       //values at the inner point
       pInner=m_cells->pvariables[PV->P][cellIdadj];
       c02=sqrt(gamma*pInner/m_cells->pvariables[PV->RHO][cellIdadj]);
       u=m_cells->pvariables[PV->U][cellIdadj];
       v=m_cells->pvariables[PV->V][cellIdadj];
       w=m_cells->pvariables[PV->W][cellIdadj];
	  
       ZFSFloat dxidx=m_cells->surfaceMetrics[cellIdadj][3];
       ZFSFloat dxidy=m_cells->surfaceMetrics[cellIdadj][4];
       ZFSFloat dxidz=m_cells->surfaceMetrics[cellIdadj][5];

       //leaving domain of integration in positive coordinate direction,
       //therefore multiply with positive F1
       ZFSFloat gradxi = F1 / sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz);
		
       ZFSFloat dxHelp = dxidx*gradxi;
       ZFSFloat dyHelp = dxidy*gradxi;
       ZFSFloat dzHelp = dxidz*gradxi;

       //values at the boundary
       pBC=F1B2*(pInner+p2+m_cells->pvariables[PV->RHO][cellIdadj]*c02*(dxHelp*(u - u2)+dyHelp*(v - v2)+dzHelp*(w - PV->WInfinity)));
       rho=rho2+((pBC-p2)/(c02*c02));

       u= u2+dxHelp*(pBC-p2)/(m_cells->pvariables[PV->RHO][cellIdadj]*c02);
       v= v2+dyHelp*(pBC-p2)/(m_cells->pvariables[PV->RHO][cellIdadj]*c02);
       w= PV->WInfinity+dzHelp*(pBC-p2)/(m_cells->pvariables[PV->RHO][cellIdadj]*c02);

       //extrapolate the variables into the ghost cells
       //gradients		
       yBC=F1B4*(m_coordinates[1][pIJK]+m_coordinates[1][pIJPK]+m_coordinates[1][pIJKP]+m_coordinates[1][pIJPKP]);
       distance=(yBC-m_cells->coordinates[1][cellIdadj]);

       drho=(rho-m_cells->pvariables[PV->RHO][cellIdadj])/distance;
       du=(u-m_cells->pvariables[PV->U][cellIdadj])/distance;
       dv=(v-m_cells->pvariables[PV->V][cellIdadj])/distance;
       dw=(w-m_cells->pvariables[PV->W][cellIdadj])/distance;
       dp=(pBC-m_cells->pvariables[PV->P][cellIdadj])/distance;
		
       //extrapolate:
       for(ZFSId jj=start[1]; jj<end[1]; ++jj){
	 cellId=cellIndex(i,jj,k);
	 distance=(m_cells->coordinates[1][cellId]-m_cells->coordinates[1][cellIdadj]);
	 m_cells->pvariables[PV->RHO][cellId]=m_cells->pvariables[PV->RHO][cellIdadj]+drho*distance;
	 m_cells->pvariables[PV->U][cellId]=m_cells->pvariables[PV->U][cellIdadj]+du*distance;
	 m_cells->pvariables[PV->V][cellId]=m_cells->pvariables[PV->V][cellIdadj]+dv*distance;
	 m_cells->pvariables[PV->W][cellId]=m_cells->pvariables[PV->W][cellIdadj]+dw*distance;
	 m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj]+dp*distance;
       }
		
     }
   }
   break;
 }
 default :{zfsTerm(1,__CALLING_FUNCTION__, "BC-face not implemented");}
 }
}


template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc3000(ZFSId bcId)
{
  const ZFSInt* start = m_physicalBCMap[bcId]->start1;
  const ZFSInt* end = m_physicalBCMap[bcId]->end1;
  const ZFSId face = m_physicalBCMap[bcId]->face;

  const ZFSId dim = face/2;
  const ZFSInt n1m1[9] = {-1,1,1, 1,-1,1, 1,1,-1};

  for(ZFSId k=start[2]; k<end[2] ; k++) {
    for(ZFSId j=end[1]-1; j>=start[1] ; j--) {
      for(ZFSId i=start[0]; i<end[0]; i++) {
        ZFSId cellId = -1;
        ZFSId cellIdAdj = -1;
        tie(cellId, cellIdAdj) = getMirrorCellIdPair(i,j,k,face);

        m_cells->pvariables[PV->RHO][cellId]=m_cells->pvariables[PV->RHO][cellIdAdj];
        m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdAdj];
        m_cells->pvariables[PV->U][cellId]=n1m1[dim*3  ]*m_cells->pvariables[PV->U][cellIdAdj];
        m_cells->pvariables[PV->V][cellId]=n1m1[dim*3+1]*m_cells->pvariables[PV->V][cellIdAdj];
        m_cells->pvariables[PV->W][cellId]=n1m1[dim*3+2]*m_cells->pvariables[PV->W][cellIdAdj];
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellId]=m_cells->pvariables[PV->RANS_VAR[0]][cellIdAdj];
        }
      }
    }
  }
}


/** \brief Reformulated Synthetic Turbulence Generation
 * Synthetic Turbulence Generation Method
 * Computes fluctuations from a given averaged
 * velocity and nu_t profile and adds them at
 * the inflow of the domain.
 * Same computations as STG in TFS by Benedikt Roidl
 * /author Marian Albers
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc7909(ZFSId bcId)
{
  const ZFSId NU_T = 5;
  const ZFSId SIJSIJ = 6;
  const ZFSId LENGTH_SCALE = 7;
  const ZFSId FLUC_UU = 8;
  const ZFSId FLUC_VV = 9;
  const ZFSId FLUC_WW = 10;
  const ZFSId FLUC_UV = 11;
  const ZFSId FLUC_VW = 12;
  const ZFSId FLUC_UW = 13;
  const ZFSId FLUC_U = 14;
  const ZFSId FLUC_V = 15;   
  const ZFSId FLUC_W = 16;
  const ZFSId LENGTH_X = 17;
  const ZFSId LENGTH_Y = 18;
  const ZFSId LENGTH_Z = 19;

  const ZFSFloat gamma= m_block->m_gamma;
  const ZFSFloat gammaMinusOne = gamma - 1.0;

  ZFSFloat xlengthmean, ylengthmean, zlengthmean, fcount;

  ZFSInt *commStgRoot = m_block->m_commStgRoot;
  MPI_Comm commStg = *m_block->m_commStg;

  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  //Initialization of variables that will be used
  ZFSFloat Temp, xmu, frho, utau2,
    utaux, xlength, ylength, zlength,
    rss, SijSij, SijSijm,Vb;
  ZFSId I, IBC;
  const ZFSId ii = 1;
  const ZFSId nran = m_block->m_stgMaxNoEddies;

  //Arrays for MPI operations
  ZFSFloatScratchSpace PMAX(4, __CALLING_FUNCTION__, "PMAX");
  ZFSFloatScratchSpace PMAXH(4, __CALLING_FUNCTION__, "PMAXH");
  ZFSFloat PMIN=F0, PMINH=F0;

  ZFSFloatScratchSpace bcast_buffer(nran*m_block->m_stgNoEddieProperties, __CALLING_FUNCTION__, "bcast_buffer");

  PMAX.fill(F0);
  PMAXH.fill(F0);
  bcast_buffer.fill(F0);


  /*********************** B *****************/
  //Value initialization for some variables
  /*********************** B *****************/

  I = cellIndex(start[0],start[1],start[2]);
  const ZFSFloat fre = 1.0 / m_block->m_Re0;

  //Limiters
  const ZFSFloat epss = 1e-34;
  const ZFSFloat eps = 1e-16;
  const ZFSFloat epsl = 1e-13;

  const ZFSFloat F072 = 0.72; //Exponent for laminar viscosity calculation
  const ZFSFloat xlengthfac = 0.1; //Potential factor to scale number of eddies

  //Box volume must be specified manually and is sign sensitive
  const ZFSFloat delta_in = m_block->m_stgDelta99Inflow;

  const ZFSFloat BLT1 = abs(m_stgVbEnd[0] - m_stgVbStart[0]);
  const ZFSFloat BLT2 = abs(m_stgVbEnd[1] - m_stgVbStart[1]);
  const ZFSFloat BLT3 = abs(m_stgVbEnd[2] - m_stgVbStart[2]);

  //Manual parameters
  const ZFSFloat c_mu = 0.09;
  const ZFSFloat a1 = 1 / sqrt(c_mu);

  const ZFSFloat timsm = 0.3; //Smoothing factor for new eddies (in time)
  const ZFSFloat aniso = 1.0; //Anisotropic, clustered eddies
  const ZFSFloat exple = m_block->m_stgExple; //to scale length of eddies

  Vb = 0.0; //Initial box volume

  ZFSId I1 = cellIndexBC(ii,start[1] + m_noGhostLayers,start[2] + m_noGhostLayers);

  ZFSFloat rhoRANSI1 = m_cells->stg_fq[PV->RHO][I1];
  ZFSFloat pressure1 = m_cells->stg_fq[PV->P][I1];
  frho = F1/rhoRANSI1;

  Temp = m_block->m_gamma*pressure1*frho;
  xmu = pow(max(Temp,F0), F072);

  //get delta0 boundary layer thickness:
  if(globalTimeStep % 250 == 0 && m_block->m_RKStep == 0) {
    ZFSInt bndryLayerIndex = 0;
    ZFSFloat u99_0 = F0, u99_1 = F0;
    for(ZFSId j = start[1] + m_noGhostLayers; j < end[1] - m_noGhostLayers; j++) {
      if( m_cells->stg_fq[PV->U][cellIndexBC(0,j,2)]  >= 0.99 * PV->UInfinity &&
          m_cells->stg_fq[PV->U][cellIndexBC(0,j-1,2)] < 0.99 * PV->UInfinity ) {
        u99_0 = m_cells->stg_fq[PV->U][cellIndexBC(0,j-1,2)];
        u99_1 = m_cells->stg_fq[PV->U][cellIndexBC(0,j,2)];
        bndryLayerIndex = j-1;
        break;
      }
    }
    ZFSFloat bndryLayerThicknessLocal = 0.0;
    ZFSFloat bndryLayerThicknessGlobal = 0.0;
    cout<<"bndryLayerIndex"<<bndryLayerIndex<<endl;
    if(bndryLayerIndex > 0) {
      bndryLayerThicknessLocal = m_cells->coordinates[1][cellIndex(0,bndryLayerIndex,2)] +
        (0.99*PV->UInfinity - u99_0)/(u99_1 - u99_0)*(m_cells->coordinates[1][cellIndex(0,bndryLayerIndex+1,2)]-m_cells->coordinates[1][cellIndex(0,bndryLayerIndex,2)]);
    }
    MPI_Allreduce(&bndryLayerThicknessLocal, &bndryLayerThicknessGlobal, 1, MPI_DOUBLE, MPI_MAX, commStg);
    if(m_block->m_stgMyRank == *commStgRoot) {cout << "Boundary Layer thickness delta99: " << bndryLayerThicknessGlobal << endl;}
  }

  ZFSFloat sr1, sr2, sr3, srt, rr1, rr2, rr3;

  //this block is only executed at the beginning as the RANS profile doesn't vary
  if((globalTimeStep <= m_block->m_restartTimeStep) ||          //junoh
     ((m_block->m_restartTimeStep==0 && globalTimeStep < 2) &&
(m_block->m_RKStep == 0)) ||
     ((globalTimeStep%m_block->m_zonalExchangeInterval==0) &&
(m_block->m_RKStep == 0))) {

    //Initialize max values
    ZFSFloat utaumax = F0, umax = F0, vmax = F0, wmax = F0, xminlength = F0;

    if(m_block->m_stgMyRank == *commStgRoot) {
      cout << "STG - Computing Reynolds Tensor Components and Length Scales..." << endl;
    }

    for(ZFSId k = start[2] + m_noGhostLayers + 1; k < end[2] - m_noGhostLayers - 1; k++) {
      for(ZFSId j = start[1] + m_noGhostLayers + 1; j < end[1] - m_noGhostLayers - 1; j++) {
        /*********************** C *****************/
        //This is the metrics / strain tensor / max shear part
        /*********************** C *****************/

        ZFSFloat dxdi, dxdj, dxdk, dydi, dydj, dydk, dzdi, dzdj, dzdk;
        ZFSFloat dxidx, dxidy, dxidz, detadx, detady, detadz, dzetadx, dzetady, dzetadz;
        ZFSFloat dudxi, dudeta, dudzeta, dvdxi, dvdeta, dvdzeta, dwdxi, dwdeta, dwdzeta;
        ZFSFloat dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dxl, dyl, dzl;
        ZFSFloat s11, s12, s13, s21, s22, s23, s31, s32, s33;
        ZFSFloat u, v, w;
        ZFSInt IPJK, IMJK, IJPK, IJMK, IJKP, IJKM;

        I = cellIndex(ii,j,k);
        IPJK = cellIndex(ii+1,j,k);
        IMJK = cellIndex(ii-1,j,k);
        IJPK = cellIndex(ii,j+1,k);
        IJMK = cellIndex(ii,j-1,k);
        IJKP = cellIndex(ii,j,k+1);
        IJKM = cellIndex(ii,j,k-1);


        dxdi = F1B2 * (m_cells->coordinates[0][IPJK] - m_cells->coordinates[0][IMJK]);
        dxdj = F1B2 * (m_cells->coordinates[0][IJPK] - m_cells->coordinates[0][IJMK]);
        dxdk = F1B2 * (m_cells->coordinates[0][IJKP] - m_cells->coordinates[0][IJKM]);
        dydi = F1B2 * (m_cells->coordinates[1][IPJK] - m_cells->coordinates[1][IMJK]);
        dydj = F1B2 * (m_cells->coordinates[1][IJPK] - m_cells->coordinates[1][IJMK]);
        dydk = F1B2 * (m_cells->coordinates[1][IJKP] - m_cells->coordinates[1][IJKM]);
        dzdi = F1B2 * (m_cells->coordinates[2][IPJK] - m_cells->coordinates[2][IMJK]);
        dzdj = F1B2 * (m_cells->coordinates[2][IJPK] - m_cells->coordinates[2][IJMK]);
        dzdk = F1B2 * (m_cells->coordinates[2][IJKP] - m_cells->coordinates[2][IJKM]);

        dxl = sqrt(dxdi*dxdi + dydi*dydi + dzdi*dzdi);
        dyl = sqrt(dxdj*dxdj + dydj*dydj + dzdj*dzdj);
        dzl = sqrt(dxdk*dxdk + dydk*dydk + dzdk*dzdk);

        dxidx = (1./max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][0];
        dxidy = (1./max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][1];
        dxidz = (1./max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][2];

        detadx = (1./max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][3 + 0];
        detady = (1./max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][3 + 1];
        detadz = (1./max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][3 + 2];

        dzetadx = (1./max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][6 + 0];
        dzetady = (1./max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][6 + 1];
        dzetadz = (1./max(m_cells->cellJac[I], epss))*m_cells->cellMetrics[I][6 + 2];

        IBC = cellIndexBC(ii,j,k);
        IPJK = cellIndexBC(ii+1,j,k);
        IMJK = cellIndexBC(ii-1,j,k);
        IJPK = cellIndexBC(ii,j+1,k);
        IJMK = cellIndexBC(ii,j-1,k);
        IJKP = cellIndexBC(ii,j,k+1);
        IJKM = cellIndexBC(ii,j,k-1);

        frho = 1.0 / m_cells->stg_fq[PV->RHO][IBC];

        //dud?
        dudxi = F1B2*(m_cells->stg_fq[PV->U][IPJK] - m_cells->stg_fq[PV->U][IMJK]);
        dudeta = F1B2*(m_cells->stg_fq[PV->U][IJPK] - m_cells->stg_fq[PV->U][IJMK]);
        dudzeta = F1B2*(m_cells->stg_fq[PV->U][IJKP] - m_cells->stg_fq[PV->U][IJKM]);

        //dvd?
        dvdxi = F1B2*(m_cells->stg_fq[PV->V][IPJK] - m_cells->stg_fq[PV->V][IMJK]);
        dvdeta = F1B2*(m_cells->stg_fq[PV->V][IJPK] - m_cells->stg_fq[PV->V][IJMK]);
        dvdzeta = F1B2*(m_cells->stg_fq[PV->V][IJKP] - m_cells->stg_fq[PV->V][IJKM]);

        //dwd?
        dwdxi = F1B2*(m_cells->stg_fq[PV->W][IPJK] - m_cells->stg_fq[PV->W][IMJK]);
        dwdeta = F1B2*(m_cells->stg_fq[PV->W][IJPK] - m_cells->stg_fq[PV->W][IJMK]);
        dwdzeta = F1B2*(m_cells->stg_fq[PV->W][IJKP] - m_cells->stg_fq[PV->W][IJKM]);

        dudx = dudxi*dxidx + dudeta*detadx + dudzeta*dzetadx;
        dudy = dudxi*dxidy + dudeta*detady + dudzeta*dzetady;
        dudz = dudxi*dxidz + dudeta*detadz + dudzeta*dzetadz;

        dvdx = dvdxi*dxidx + dvdeta*detadx + dvdzeta*dzetadx;
        dvdy = dvdxi*dxidy + dvdeta*detady + dvdzeta*dzetady;
        dvdz = dvdxi*dxidz + dvdeta*detadz + dvdzeta*dzetadz;

        dwdx = dwdxi*dxidx + dwdeta*detadx + dwdzeta*dzetadx;
        dwdy = dwdxi*dxidy + dwdeta*detady + dwdzeta*dzetady;
        dwdz = dwdxi*dxidz + dwdeta*detadz + dwdzeta*dzetadz;


        s11 = 2.0*dudx;
        s12 = dvdx+dudy;
        s13 = dwdx+dudz;

        s21 = dudy+dvdx;
        s22 = 2.0*dvdy;
        s23 = dwdy+dvdz;

        s31 = dudz+dwdx;
        s32 = dvdz+dwdy;
        s33 = 2.0*dwdz;

        //Strain tensor
        SijSij = F1B4 * (s11*s11 + s12*s12 + s13*s13 +
                         s21*s21 + s22*s22 + s23*s23 +
                         s31*s31 + s32*s32 + s33*s33);

        if(std::isnan(SijSij)) {
          cout << " dudxi  : " << dudxi
               << " dudeta : " << dudeta
               << " dudzeta: " << dudzeta
               << " dvdxi  : " << dvdxi
               << " dvdeta : " << dvdeta
               << " dvdzeta: " << dvdzeta
               << " dwdxi  : " << dwdxi
               << " dwdeta : " << dwdeta
               << " dwdzeta: " << dwdzeta << endl;
        }


        //>marian: in TFS code this isn't SQRT
        m_cells->stg_fq[SIJSIJ][cellIndexBC(ii,j,k)] = SijSij;
        //Assume a one-, or two equation turbulence model
        //Assume a simplified directivity:

        //Read from RANS profile
        ZFSFloat nu_t = m_cells->stg_fq[NU_T][cellIndexBC(m_noGhostLayers,j,k)];
        sr1 = (s12 + s21)*(s12 + s21);
        sr2 = (s23 + s32)*(s23 + s32);
        sr3 = (s13 + s31)*(s13 + s31);
        srt = max(sqrt(sr1 + sr2 + sr3), epsl);

        rr1 = sqrt(sr1)/srt;
        rr2 = sqrt(sr2)/srt;
        rr3 = sqrt(sr3)/srt;

        const ZFSFloat uv = -sqrt(2.0*SijSij)*rr1*nu_t*fre*frho;
        const ZFSFloat vw = -sqrt(2.0*SijSij)*rr2*nu_t*fre*frho;
        const ZFSFloat uw = -sqrt(2.0*SijSij)*rr3*nu_t*fre*frho;

        //>marian: edited to fit TFS
        //Normal stresses are assumed to
        //uu = abs(uv*2.0d0) + abs(vw*2.0d0) + abs(uw*2.0d0)
        const ZFSFloat uu = a1*abs(uv);// sqrt(2.*SijSij)*nu_t*fre*frho;
        const ZFSFloat vv = uu*0.7; // uu
        const ZFSFloat ww = uu*0.8; // uu

        m_cells->stg_fq[FLUC_UU][cellIndexBC(ii,j,k)] = uu;
        m_cells->stg_fq[FLUC_VV][cellIndexBC(ii,j,k)] = vv;
        m_cells->stg_fq[FLUC_WW][cellIndexBC(ii,j,k)] = ww;
        m_cells->stg_fq[FLUC_UV][cellIndexBC(ii,j,k)] = uv;
        m_cells->stg_fq[FLUC_VW][cellIndexBC(ii,j,k)] = vw;
        m_cells->stg_fq[FLUC_UW][cellIndexBC(ii,j,k)] = uw;

        //Get utau using laminar viscosity
        //                utau2 = sqrt(fre*sqrt(2.0*SijSij)*xmu/(m_FQ[IBC][4]*F1));
        utau2 = sqrt(fre*sqrt(2.0*SijSij)*xmu*frho);

        //Save values if they are the new maximum
        if (utau2 >= utaumax) {
          xminlength = pow((dxl*dyl*dzl),0.33);
          utaumax = max(utau2, utaumax);
        }

        //In which direction aims the maximum averaged velocity?
        u = m_cells->stg_fq[PV->U][IBC];
        v = m_cells->stg_fq[PV->V][IBC];
        w = m_cells->stg_fq[PV->W][IBC];

        umax = max(u, umax);
        vmax = max(v, vmax);
        wmax = max(w, wmax);

        //We need a global length scale to compare
        m_cells->stg_fq[LENGTH_SCALE][cellIndexBC(ii,j,k)] = pow(sqrt(max(2.0*SijSij, epss)), -exple);



        /***************** End of C *****************/

        /*********************** D *****************/
        // Zero gradient extrapolation of boundary values
        /*********************** D *****************/
        ZFSId noSTGVariables = 14;

        if(k == start[2] + m_noGhostLayers + 1) {
          //1st layer
          for(int var = 6; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j, k-1)] = m_cells->stg_fq[var][IBC];
          }
          //2nd layer
          for(int var = 0; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j, k-2)] = m_cells->stg_fq[var][IBC];
          }
          //3rd layer
          for(int var = 0; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j, k-3)] = m_cells->stg_fq[var][IBC];
          }
        } else if(k == end[2] - m_noGhostLayers - 2) {
          //1st layer
          for(int var = 6; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j, k+1)] = m_cells->stg_fq[var][IBC];
          }
          //2nd layer
          for(int var = 0; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j, k+2)] = m_cells->stg_fq[var][IBC];
          }
          //3rd layer
          for(int var = 0; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j, k+3)] = m_cells->stg_fq[var][IBC];
          }
        }

        if(j == start[1] + m_noGhostLayers + 1) {
          //1st layer
          for(int var = 6; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j-1, k)] = m_cells->stg_fq[var][IBC];
          }
          //2nd layer
          for(int var = 0; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j-2, k)] = m_cells->stg_fq[var][IBC];
          }
          //3rd layer
          for(int var = 0; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j-3, k)] = m_cells->stg_fq[var][IBC];
          }
        } else if(j == end[1] - m_noGhostLayers - 2) {
          //1st layer
          for(int var = 6; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j+1, k)] = m_cells->stg_fq[var][IBC];
          }
          //2nd layer
          for(int var = 0; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j+2, k)] = m_cells->stg_fq[var][IBC];
          }
          //3rd layer
          for(int var = 0; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j+3, k)] = m_cells->stg_fq[var][IBC];
          }
        }

        //lower left edge
        if(k == start[2] + m_noGhostLayers + 1 && j == start[1] + m_noGhostLayers + 1) {
          for(int var = 0; var < noSTGVariables; var++) {
            if(var > 5) {
              m_cells->stg_fq[var][cellIndexBC(ii, j-1, k-1)] = m_cells->stg_fq[var][cellIndexBC(ii, j-2, k)];
            }

            m_cells->stg_fq[var][cellIndexBC(ii, j-1, k-2)] = m_cells->stg_fq[var][cellIndexBC(ii, j-2, k-1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j-1, k-3)] = m_cells->stg_fq[var][cellIndexBC(ii, j-2, k-2)];

            m_cells->stg_fq[var][cellIndexBC(ii, j-2, k-1)] = m_cells->stg_fq[var][cellIndexBC(ii, j-2, k)];
            m_cells->stg_fq[var][cellIndexBC(ii, j-2, k-2)]= m_cells->stg_fq[var][cellIndexBC(ii, j-2, k-1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j-2, k-3)] = m_cells->stg_fq[var][cellIndexBC(ii, j-2, k-2)];

            m_cells->stg_fq[var][cellIndexBC(ii, j-3, k-1)] = m_cells->stg_fq[var][cellIndexBC(ii, j-3, k)];
            m_cells->stg_fq[var][cellIndexBC(ii, j-3, k-2)] = m_cells->stg_fq[var][cellIndexBC(ii, j-3, k-1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j-3, k-3)] = m_cells->stg_fq[var][cellIndexBC(ii, j-3, k-2)];
          }
        }

        //lower right edge
        if(k == end[2] - m_noGhostLayers - 2 && j == start[1] + m_noGhostLayers + 1) {
          for(int var = 0; var < noSTGVariables; var++) {
            if(var > 5) {
              m_cells->stg_fq[var][cellIndexBC(ii, j-1, k+1)] = m_cells->stg_fq[var][cellIndexBC(ii, j-1, k)];
            }

            m_cells->stg_fq[var][cellIndexBC(ii, j-1, k+2)] = m_cells->stg_fq[var][cellIndexBC(ii, j-1, k+1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j-1, k+3)] = m_cells->stg_fq[var][cellIndexBC(ii, j-1, k+2)];

            m_cells->stg_fq[var][cellIndexBC(ii, j-2, k+1)] = m_cells->stg_fq[var][cellIndexBC(ii, j-2, k)];
            m_cells->stg_fq[var][cellIndexBC(ii, j-2, k+2)] = m_cells->stg_fq[var][cellIndexBC(ii, j-2, k+1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j-2, k+3)] = m_cells->stg_fq[var][cellIndexBC(ii, j-2, k+2)];

            m_cells->stg_fq[var][cellIndexBC(ii, j-3, k+1)] = m_cells->stg_fq[var][cellIndexBC(ii, j-3, k)];
            m_cells->stg_fq[var][cellIndexBC(ii, j-3, k+2)] = m_cells->stg_fq[var][cellIndexBC(ii, j-3, k+1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j-3, k+3)] = m_cells->stg_fq[var][cellIndexBC(ii, j-3, k+2)];
          }
        }

        //upper left edge
        if(k == start[2] + m_noGhostLayers + 1 && j == end[1] - m_noGhostLayers - 2) {
          for(int var = 0; var < noSTGVariables; var++) {
            if(var > 5) {
              m_cells->stg_fq[var][cellIndexBC(ii, j+1, k-1)] = m_cells->stg_fq[var][cellIndexBC(ii, j+1, k)];
            }

            m_cells->stg_fq[var][cellIndexBC(ii, j+1, k-2)] = m_cells->stg_fq[var][cellIndexBC(ii, j+1, k-1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j+1, k-3)] = m_cells->stg_fq[var][cellIndexBC(ii, j+1, k-2)];

            m_cells->stg_fq[var][cellIndexBC(ii, j+2, k-1)] = m_cells->stg_fq[var][cellIndexBC(ii, j+2, k)];
            m_cells->stg_fq[var][cellIndexBC(ii, j+2, k-2)] = m_cells->stg_fq[var][cellIndexBC(ii, j+2, k-1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j+2, k-3)] = m_cells->stg_fq[var][cellIndexBC(ii, j+2, k-2)];

            m_cells->stg_fq[var][cellIndexBC(ii, j+3, k-1)] = m_cells->stg_fq[var][cellIndexBC(ii, j+3, k)];
            m_cells->stg_fq[var][cellIndexBC(ii, j+3, k-2)] = m_cells->stg_fq[var][cellIndexBC(ii, j+3, k-1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j+3, k-3)] = m_cells->stg_fq[var][cellIndexBC(ii, j+3, k-2)];
          }
        }

        //upper right edge
        if(k == end[2] - m_noGhostLayers - 2 && j == end[1] - m_noGhostLayers - 2) {
          for(int var = 0; var < noSTGVariables; var++) {
            if(var > 5) {
              m_cells->stg_fq[var][cellIndexBC(ii, j+1, k+1)] = m_cells->stg_fq[var][cellIndexBC(ii, j+1, k)];
            }

            m_cells->stg_fq[var][cellIndexBC(ii, j+1, k+2)] = m_cells->stg_fq[var][cellIndexBC(ii, j+1, k+1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j+1, k+3)] = m_cells->stg_fq[var][cellIndexBC(ii, j+1, k+2)];

            m_cells->stg_fq[var][cellIndexBC(ii, j+2, k+1)] = m_cells->stg_fq[var][cellIndexBC(ii, j+2, k)];
            m_cells->stg_fq[var][cellIndexBC(ii, j+2, k+2)] = m_cells->stg_fq[var][cellIndexBC(ii, j+2, k+1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j+2, k+3)] = m_cells->stg_fq[var][cellIndexBC(ii, j+2, k+2)];

            m_cells->stg_fq[var][cellIndexBC(ii, j+3, k+1)] = m_cells->stg_fq[var][cellIndexBC(ii, j+3, k)];
            m_cells->stg_fq[var][cellIndexBC(ii, j+3, k+2)] = m_cells->stg_fq[var][cellIndexBC(ii, j+3, k+1)];
            m_cells->stg_fq[var][cellIndexBC(ii, j+3, k+3)] = m_cells->stg_fq[var][cellIndexBC(ii, j+3, k+2)];
          }
        }

        /***************** End of D *****************/

        /*********************** E *****************/
        // Storage of values
        /*********************** E *****************/

        //Save max direction vector and max tau
        PMAX[0] = umax;
        PMAX[1] = vmax;
        PMAX[2] = wmax;
        PMAX[3] = utaumax;

        //To get the dy (or any other direction) determine the maximum value of the stress matrix
        //Get corner value of virtual box
        //Save corner values of the virtual box
        PMIN = xminlength;

        /***************** End of E *****************/
      }
    }

    /*********************** F *****************/
    // Communication: Exchange min and max values
    /*********************** F *****************/

    MPI_Allreduce(PMAX.begin(), PMAXH.begin(), 4, MPI_DOUBLE, MPI_MAX, commStg);
    MPI_Allreduce(&PMIN, &PMINH, 1, MPI_DOUBLE, MPI_MIN, commStg);

    PMIN = PMINH;
    for (ZFSInt i = 0; i < 4; i++) {
      PMAX[i] = PMAXH[i];
    }

    /***************** End of F *****************/

    /*********************** G *****************/
    //
    /*********************** G *****************/

    //Maximum convection velocities at inflow
    m_stgMaxVel[0] = PMAX[0];
    m_stgMaxVel[1] = PMAX[1];
    m_stgMaxVel[2] = PMAX[2];

    xlengthmean = F0;
    ylengthmean = F0;
    zlengthmean = F0;
    fcount = 1.0;

    for (ZFSInt k = start[2]; k < end[2]; k++) {
      for (ZFSInt j = start[1]; j < end[1]; j++) {
        I = cellIndex(ii,j,k);
        IBC = cellIndexBC(ii,j,k);
        utaux = PMAX[3];
        xminlength = PMIN;

        //Length scale in main flow direction
        xlength = 1.5*delta_in;


        //Length scale in the direction of main shear
        ylength = max(min(0.6
                          *max(m_cells->stg_fq[LENGTH_SCALE][IBC]*delta_in*pow(utaux/delta_in, exple), eps)
                          , delta_in*0.66), xminlength);

        //Length scale in the direction perpendicular of x and y
        zlength = ylength*1.5;

        //>marian
        m_cells->stg_fq[LENGTH_X][IBC]= xlength;
        m_cells->stg_fq[LENGTH_Y][IBC]= ylength;
        m_cells->stg_fq[LENGTH_Z][IBC]= zlength;
        //<marian

        //Get maximum strain and local strain to weight the length scales
        //The higher the strain, the more eddies are swimming there...
        SijSijm = (utaux*utaux/fre*F1B2/xmu*m_cells->stg_fq[PV->RHO][IBC])*
          (utaux*utaux/fre*F1B2/xmu*m_cells->stg_fq[PV->RHO][IBC])*F1B2;
        SijSij = m_cells->stg_fq[SIJSIJ][IBC] * m_cells->stg_fq[SIJSIJ][IBC];

        rss = pow((SijSij/SijSijm),xlengthfac);
        xlengthmean = xlengthmean + xlength*rss;
        ylengthmean = ylengthmean + ylength*rss;
        zlengthmean = zlengthmean + zlength*rss;
        fcount++;
      }
    }

    PMAX[0] = ylengthmean;
    PMAX[1] = zlengthmean;
    PMAX[2] = xlengthmean;
    PMAX[3] = fcount;

    // /***************** End of G *****************/


    // /*********************** H *****************/
    // // Communication to determine number of eddy cores
    // /*********************** H *****************/


    MPI_Allreduce(PMAX.begin(), PMAXH.begin(), 4, MPI_DOUBLE, MPI_SUM, commStg);


    for (ZFSInt i = 0; i < 4; i++) {
      PMAX[i] = PMAXH[i];
    }

    //Averaged length scale in at the inflow
    ylengthmean = PMAX[0] / PMAX[3];
    zlengthmean = PMAX[1] / PMAX[3];
    xlengthmean = PMAX[2] / PMAX[3];

    const ZFSId noCellsJ = m_block->m_totalGridBlockDim[m_block->m_inputBlockId][1]-1;
    const ZFSId noVars = 6;
    ZFSFloatScratchSpace localLengthScales(noCellsJ, noVars, __CALLING_FUNCTION__, "localLengthScales");
    ZFSFloatScratchSpace globalLengthScales(noCellsJ, noVars, __CALLING_FUNCTION__, "localLengthScales");
    localLengthScales.fill(F0);
    globalLengthScales.fill(F0);

    for (ZFSInt k = m_noGhostLayers; k < m_nCells[0]-m_noGhostLayers; k++) {
      for(ZFSId j = m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
        IBC = cellIndexBC(ii,j,k);
        ZFSId cellId = cellIndex(ii,j,k);
        ZFSId localJ = m_block->m_nOffsetCells[1] + j - m_noGhostLayers;
        localLengthScales(localJ,0) += localJ;
        localLengthScales(localJ,1) += m_cells->coordinates[1][cellId];
        localLengthScales(localJ,2) += m_cells->stg_fq[LENGTH_X][IBC];
        localLengthScales(localJ,3) += m_cells->stg_fq[LENGTH_Y][IBC];
        localLengthScales(localJ,4) += m_cells->stg_fq[LENGTH_Z][IBC];
        localLengthScales(localJ,5) += F1;
      }
      //manually set first coordinate to zero
      localLengthScales(0,1) = F0;
    }

    MPI_Allreduce(localLengthScales.begin(), globalLengthScales.begin(), noVars*noCellsJ, MPI_DOUBLE, MPI_SUM, commStg);

    for(ZFSId j = 0; j<noCellsJ; j++) {
      ZFSFloat FNoCells = F1/globalLengthScales(j,5);
      for(ZFSId var=0; var<noVars-1; var++) {
        m_stgGlobalLengthScales[j][var] = FNoCells*globalLengthScales(j,var);
      }
    }

    //Create new eddies if it's an initial start or the createNewEddie flag is set
    if(m_block->m_stgCreateNewEddies&&m_block->m_stgEddieLengthScales) {
      ZFSFloat xk1t, xk2t, xk3t, epsik1, epsik2, epsik3;
      ZFSFloat l1,lmin1,lmax1,l2,lmin2,lmax2,l3,lmin3,lmax3; //each eddy has diff length
      ZFSFloat beta1=F1,beta2=F1,beta3=F1;

      if(m_block->m_stgMyRank == *commStgRoot) {
        for(ZFSInt n = 0; n < m_block->m_stgMaxNoEddies; n++) {
          xk1t = m_stgVbStart[0] + generate_rand()*(m_stgVbEnd[0] - m_stgVbStart[0]);
          xk2t = m_stgVbStart[1] + generate_rand_weighted()*(m_stgVbEnd[1] - m_stgVbStart[1]);
          xk3t = m_stgVbStart[2] + generate_rand()*(m_stgVbEnd[2] - m_stgVbStart[2]);

          epsik1 = 2.0*generate_rand() - 1.0;
          epsik1 = epsik1/max(abs(epsik1),eps);
          epsik2 = 2.0*generate_rand() - 1.0;
          epsik2 = epsik2/max(abs(epsik2),eps);
          epsik3 = 2.0*generate_rand() - 1.0;
          epsik3 = epsik3/max(abs(epsik3),eps);

          bcast_buffer[n + nran*0] = xk1t;
          bcast_buffer[n + nran*1] = xk2t;
          bcast_buffer[n + nran*2] = xk3t;
          bcast_buffer[n + nran*3] = epsik1;
          bcast_buffer[n + nran*4] = epsik2;
          bcast_buffer[n + nran*5] = epsik3;

          ZFSFloat omega=generate_rand();

          //find correct cell
          for(ZFSId j=0; j<noCellsJ-1; j++) {
            if(m_stgGlobalLengthScales[j][1] < xk2t && m_stgGlobalLengthScales[j+1][1] > xk2t) {
              lmin1 = m_stgGlobalLengthScales[j][2]*0.95;
              lmax1 = m_stgGlobalLengthScales[j][2]*1.0;
              lmin2 = m_stgGlobalLengthScales[j][3]*0.95;
              lmax2 = m_stgGlobalLengthScales[j][3]*1.0;
              lmin3 = m_stgGlobalLengthScales[j][4]*0.95;
              lmax3 = m_stgGlobalLengthScales[j][4]*1.0;

              cout << "This eddie has xk2t: " << xk2t << " and located between j: " << j << " and j+1: " << j+1 << " lx: " << lmax1 << " ly: " << lmax2 << endl;
              break;
            }
          }

          l1=lmin1+(lmax1-lmin1)*pow(omega,beta1);
          l2=lmin2+(lmax2-lmin2)*pow(omega,beta2);
          l3=lmin3+(lmax3-lmin3)*pow(omega,beta3);

          bcast_buffer[n + nran*6] = l1;
          bcast_buffer[n + nran*7] = l2;
          bcast_buffer[n + nran*8] = l3;
        }
      }
      MPI_Bcast(bcast_buffer.begin(), m_block->m_stgNoEddieProperties*nran, MPI_DOUBLE, *commStgRoot, commStg);
      for(ZFSInt n = 0; n < nran; n++) {
        for(ZFSInt p = 0; p < m_block->m_stgNoEddieProperties; p++) {
          m_block->m_stgEddies[n][p] = bcast_buffer[n + nran*p]; //bcast_buffer(n,var);
        }
      }
    }
  }

  if(m_block->m_RKStep == 0) {
    /*********************** J *****************/
    // The virtual box part - executed by Master Block at the inflow
    /*********************** J *****************/
    ZFSFloat xk1t, xk2t, xk3t, epsik1, epsik2, epsik3;
    ZFSFloat l1,lmin1,lmax1,l2,lmin2,lmax2,l3,lmin3,lmax3; //each eddy has diff length
    ZFSFloat beta1,beta2,beta3;
    // ZFSFloat direciton;
    if(m_block->m_stgEddieLengthScales) {
      lmin1=0.050;lmax1=0.37582; //average : 0.21291
      lmin2=0.005;lmax2=0.01807;    //  0.0115
      lmin3=0.010;lmax3=0.03614;    //0.0115
      beta1=1.0; beta2=1.0; beta3=1.0;
    }

    if(m_block->m_stgMyRank == *commStgRoot) {
      for(ZFSInt n = 0; n < nran; n++) {
        xk1t = m_block->m_stgEddies[n][0];
        xk2t = m_block->m_stgEddies[n][1];
        xk3t = m_block->m_stgEddies[n][2];

        //Check if the eddie has left the Virtual Box
        if (xk1t > m_stgVbEnd[0] || xk1t < m_stgVbStart[0] ||
            xk2t > m_stgVbEnd[1] || xk2t < m_stgVbStart[1] ||
            xk3t > m_stgVbEnd[2] || xk3t < m_stgVbStart[2] ) {
          //Get coordinates of eddie cores and their signs
          // cout  << "Old eddie with position: " << xk1t << " , " << xk2t << " , " << xk3t << endl;
          xk1t = m_stgVbStart[0];
          xk2t = m_stgVbStart[1] + generate_rand_weighted()*BLT2;
          xk3t = m_stgVbStart[2] + generate_rand()*BLT3;
          // cout  << "Creating new eddie with position: " << xk1t << " , " << xk2t << " , " << xk3t << endl;
          epsik1 = 2.0*generate_rand() - 1.0;
          epsik1 = epsik1/max(abs(epsik1),eps);
          epsik2 = 2.0*generate_rand() - 1.0;
          epsik2 = epsik2/max(abs(epsik2),eps);
          epsik3 = 2.0*generate_rand() - 1.0;
          epsik3 = epsik3/max(abs(epsik3),eps);

          if(m_block->m_stgEddieLengthScales) {
            ZFSFloat omega=generate_rand();
            //find correct cell
            const ZFSId noCellsJ = m_block->m_totalGridBlockDim[m_block->m_inputBlockId][1]-1;
            for(ZFSId j=0; j<noCellsJ-1; j++) {
              if(m_stgGlobalLengthScales[j][1] < xk2t && m_stgGlobalLengthScales[j+1][1] > xk2t) {
                lmin1 = m_stgGlobalLengthScales[j][2]*0.95;
                lmax1 = m_stgGlobalLengthScales[j][2]*1.1;
                lmin2 = m_stgGlobalLengthScales[j][3]*0.95;
                lmax2 = m_stgGlobalLengthScales[j][3]*1.0;
                lmin3 = m_stgGlobalLengthScales[j][4]*0.95;
                lmax3 = m_stgGlobalLengthScales[j][4]*1.0;
              }
            }

            l1=lmin1+(lmax1-lmin1)*pow(omega,beta1);
            l2=lmin2+(lmax2-lmin2)*pow(omega,beta2);
            l3=lmin3+(lmax3-lmin3)*pow(omega,beta3);
          }
        } else {
          //>marian: difference to TFS Code
          //Convect eddy cores with a certain velocity
          xk1t = m_stgMaxVel[0]*m_block->m_timeStep + m_block->m_stgEddies[n][0];
          xk2t = m_stgMaxVel[1]*m_block->m_timeStep + m_block->m_stgEddies[n][1];
          xk3t = m_stgMaxVel[2]*m_block->m_timeStep + m_block->m_stgEddies[n][2];

          epsik1 = m_block->m_stgEddies[n][3];
          epsik2 = m_block->m_stgEddies[n][4];
          epsik3 = m_block->m_stgEddies[n][5];

          if(m_block->m_stgEddieLengthScales) {
            l1 = m_block->m_stgEddies[n][6];
            l2 = m_block->m_stgEddies[n][7];
            l3 = m_block->m_stgEddies[n][8];
          }
        }

        bcast_buffer[n + nran*0] = xk1t;
        bcast_buffer[n + nran*1] = xk2t;
        bcast_buffer[n + nran*2] = xk3t;
        bcast_buffer[n + nran*3] = epsik1;
        bcast_buffer[n + nran*4] = epsik2;
        bcast_buffer[n + nran*5] = epsik3;
        if(m_block->m_stgEddieLengthScales) {
          bcast_buffer[n + nran*6] = l1;
          bcast_buffer[n + nran*7] = l2;
          bcast_buffer[n + nran*8] = l3;
        }
      }
    }

    //Broadcast the new/updated eddies to all relevant processes
    MPI_Bcast(bcast_buffer.begin(), nran*m_block->m_stgNoEddieProperties, MPI_DOUBLE, *commStgRoot, commStg);

    //Copy data into m_FQeddie vector
    for(ZFSInt n = 0; n < nran; n++) {
      for(ZFSInt p = 0; p < m_block->m_stgNoEddieProperties; p++) {
        m_block->m_stgEddies[n][p] = bcast_buffer[n+nran*p];
      }
    }

    Vb = BLT2*BLT3*BLT1;

    //Summary of synth turb parameters
    if (m_block->m_stgMyRank == *commStgRoot && globalTimeStep == m_block->m_restartTimeStep) {
      if(!std::isnan(xlengthmean)) {
        cout << "**************************" << endl
             << "Synthetic turbulence:" << endl
             << "zones: 1" << endl
             << "utau_max: " << utaux << endl
             << "nr. eddies: " << nran << endl
             << "conv. vel: " << sqrt(POW2(m_stgMaxVel[0]) + POW2(m_stgMaxVel[1]) + POW2(m_stgMaxVel[2])) << endl
             << "umax = " << m_stgMaxVel[0] << endl
             << "vmax = " << m_stgMaxVel[1] << endl
             << "wmax = " << m_stgMaxVel[2] << endl
             << "virtual box volume: " << Vb << endl
             << "Vb/nran = " << Vb/nran << endl
             << "xlengthmean = " << xlengthmean << endl
             << "ylengthmean = " << ylengthmean << endl
             << "zlengthmean = " << zlengthmean << endl
             << "**************************" << endl;
      }
    }

    /*********************** L *****************/
    // Calculation of the fluctuation induced by all eddies on each cell
    /*********************** L *****************/

    ZFSId iStart = 1, iEnd = 2;
    if(m_block->m_stgShapeFunction!=4) {
      iStart = 0;
      iEnd = 3;
    }

    //only compute fluctuations for second
    for(ZFSInt k = start[2]; k < end[2]; k++) {
      for(ZFSInt j = start[1]; j < end[1]; j++) {
        for(ZFSId i = iStart; i < iEnd; i++) {
          ZFSFloat xLb1, xLb2, xLb3,
            zacfq1, zacfq2, zacfq3,
            rol1H, rol2H, rol3H, fl1, fl2, fl3, fH1, fH2, fH3, velmax,
            ufluc, vfluc, wfluc;
          ZFSFloat help1=F0, help2=F0, help3=F0, help4=F0, help5=F0, help6=F0;

          const ZFSFloat umax = m_stgMaxVel[0];
          const ZFSFloat vmax = m_stgMaxVel[1];
          const ZFSFloat wmax = m_stgMaxVel[2];

          //the tensor components and xyzlengths are only saved in one row (ii = 1)
          const ZFSId cellIdBC = cellIndexBC(i,j,k);
          const ZFSId cellIdBCFirst = cellIndexBC(ii,j,k);
          const ZFSId cellId = cellIndex(ii,j,k);

          const ZFSFloat uu = m_cells->stg_fq[FLUC_UU][cellIdBCFirst];
          const ZFSFloat vv = m_cells->stg_fq[FLUC_VV][cellIdBCFirst];
          const ZFSFloat ww = m_cells->stg_fq[FLUC_WW][cellIdBCFirst];
          const ZFSFloat uv = m_cells->stg_fq[FLUC_UV][cellIdBCFirst];
          const ZFSFloat vw = m_cells->stg_fq[FLUC_VW][cellIdBCFirst];
          const ZFSFloat uw = m_cells->stg_fq[FLUC_UW][cellIdBCFirst];

          //Cholesky decomposition of the Reynolds stress tensor
          const ZFSFloat a11 = sqrt(max(uu,epsl));
          const ZFSFloat a21 = uv / a11;
          const ZFSFloat a31 = uw / a11;
          const ZFSFloat a22 = sqrt(max((vv - a21*a21),epsl));
          const ZFSFloat a32 = (vw - a21*a31) / a22;
          const ZFSFloat a33 = sqrt(max((ww - a31*a31 - a32*a32),epsl));

          for (ZFSInt n = 0; n < nran; n++) {
            //Tent function to determine the symmetric function to model the
            //decay of the fluctuations
            xk1t = m_block->m_stgEddies[n][0];
            xk2t = m_block->m_stgEddies[n][1];
            xk3t = m_block->m_stgEddies[n][2];

            if(m_block->m_stgEddieLengthScales) {
              xLb1 = m_block->m_stgEddies[n][6];
              xLb2 = m_block->m_stgEddies[n][7];
              xLb3 = m_block->m_stgEddies[n][8];
            } else {
              //Length scales
              xLb1 = m_cells->stg_fq[LENGTH_X][cellIdBCFirst];
              xLb2 = m_cells->stg_fq[LENGTH_Y][cellIdBCFirst];
              xLb3 = m_cells->stg_fq[LENGTH_Z][cellIdBCFirst];
            }

            zacfq1 = (m_cells->coordinates[0][cellId] - xk1t) /
              abs(m_cells->coordinates[0][cellId] - xk1t);
            rol1H = zacfq1*min(abs(m_cells->coordinates[0][cellId] - xk1t)/
                               xLb1,1.0);
            zacfq2 = (m_cells->coordinates[1][cellId] - xk2t) /
              abs(m_cells->coordinates[1][cellId] - xk2t);
            rol2H = zacfq2*min(abs(m_cells->coordinates[1][cellId] - xk2t)/
                               xLb2,1.0);
            zacfq3 = (m_cells->coordinates[2][cellId] - xk3t) /
              abs(m_cells->coordinates[2][cellId] - xk3t);
            rol3H = zacfq3*min(abs(m_cells->coordinates[2][cellId] - xk3t)/
                               xLb3,1.0);

            switch(m_block->m_stgShapeFunction) {
            case 0: {
              fl1 = 2.0/sqrt(xLb1)/pow(6.28318531,0.25)*
                exp(-((m_cells->coordinates[0][cellId] - xk1t)*2/xLb1)*
                    ((m_cells->coordinates[0][cellId] - xk1t)*2/xLb1));
              fl2 = 2.0/sqrt(xLb2)/pow(6.28318531,0.25)*
                exp(-((m_cells->coordinates[1][cellId] - xk2t)*2/xLb2)*
                    ((m_cells->coordinates[1][cellId] - xk2t)*2/xLb2));
              fl3 = 2.0/sqrt(xLb3)/pow(6.28318531,0.25)*
                exp(-((m_cells->coordinates[2][cellId] - xk3t)*2/xLb3)*
                    ((m_cells->coordinates[2][cellId] - xk3t)*2/xLb3));
              break;
            }
            case 1: {
              fl1=sqrt(1.5/xLb1)/pow(cosh((m_cells->coordinates[0][cellId] - xk1t)*2/xLb1),2);
              fl2=sqrt(1.5/xLb2)/pow(cosh((m_cells->coordinates[1][cellId] - xk2t)*2/xLb2),2);
              fl3=sqrt(1.5/xLb3)/pow(cosh((m_cells->coordinates[2][cellId] - xk3t)*2/xLb3),2);
              break;
            }
            case 2: {
              fl1 = 2.0/sqrt(xLb1)/pow(6.28318531,0.25)*
                exp(-((m_cells->coordinates[0][cellId] - xk1t)*2/xLb1)*
                    ((m_cells->coordinates[0][cellId] - xk1t)*2/xLb1));
              fl2 = 2.0/sqrt(xLb2)/pow(6.28318531,0.25)*
                exp(-((m_cells->coordinates[1][cellId] - xk2t)*2/xLb2)*
                    ((m_cells->coordinates[1][cellId] - xk2t)*2/xLb2));
              fl3 = 2.0/sqrt(xLb3)/pow(6.28318531,0.25)*
                exp(-((m_cells->coordinates[2][cellId] - xk3t)*2/xLb3)*
                    ((m_cells->coordinates[2][cellId] - xk3t)*2/xLb3));
              break;
            }
            case 3: {
              fl1 = sqrt(27/pow(PI,1.5)/xLb1/xLb2/xLb3)*
                (-3*(m_cells->coordinates[1][cellId] - xk2t)/xLb2+3*(m_cells->coordinates[2][cellId] - xk3t)/xLb3)*
                exp(-4.5*(pow((m_cells->coordinates[0][cellId] - xk1t)/xLb1,2)+
                          pow((m_cells->coordinates[1][cellId] - xk2t)/xLb2,2)+
                          pow((m_cells->coordinates[2][cellId] - xk3t)/xLb3,2)));
              fl2 = sqrt(27/pow(PI,1.5)/xLb1/xLb2/xLb3)*
                (-3*(m_cells->coordinates[2][cellId] - xk3t)/xLb3+3*(m_cells->coordinates[0][cellId] - xk1t)/xLb1)*
                exp(-4.5*(pow((m_cells->coordinates[0][cellId] - xk1t)/xLb1,2)+
                          pow((m_cells->coordinates[1][cellId] - xk2t)/xLb2,2)+
                          pow((m_cells->coordinates[2][cellId] - xk3t)/xLb3,2)));
              fl3 = sqrt(27/pow(PI,1.5)/xLb1/xLb2/xLb3)*
                (-3*(m_cells->coordinates[0][cellId] - xk1t)/xLb1+3*(m_cells->coordinates[1][cellId] - xk2t)/xLb2)*
                exp(-4.5*(pow((m_cells->coordinates[0][cellId] - xk1t)/xLb1,2)+
                          pow((m_cells->coordinates[1][cellId] - xk2t)/xLb2,2)+
                          pow((m_cells->coordinates[2][cellId] - xk3t)/xLb3,2)));
              break;
            }
            case 4:
            default: {
              fl1 = 2.0/sqrt(3.141*xLb1)*
                exp(-((m_cells->coordinates[0][cellId] - xk1t)*2/xLb1)*
                    ((m_cells->coordinates[0][cellId] - xk1t)*2/xLb1));
              fl2 = 2.0/sqrt(3.141*xLb2)*
                exp(-((m_cells->coordinates[1][cellId] - xk2t)*2/xLb2)*
                    ((m_cells->coordinates[1][cellId] - xk2t)*2/xLb2));
              fl3 = 2.0/sqrt(3.141*xLb3)*
                exp(-((m_cells->coordinates[2][cellId] - xk3t)*2/xLb3)*
                    ((m_cells->coordinates[2][cellId] - xk3t)*2/xLb3));
            }
            }

            //Normalization factor cannot be chosen as Pamies did... or we...
            if(m_block->m_stgShapeFunction==4) {
              fH1 = (F1 - cos(2.0*3.141*rol1H)) / (2.0*3.141*rol1H*0.44);
              fH2 = (F1 - cos(2.0*3.141*rol2H)) / (2.0*3.141*rol2H*0.44);
              fH3 = (F1 - cos(2.0*3.141*rol3H)) / (2.0*3.141*rol3H*0.44);
            } else {
              fH1 = (F1 - cos(2.0*3.141592654*rol1H)) / (2.0*3.141592654*rol1H*0.65409);
              fH2 = (F1 - cos(2.0*3.141592654*rol2H)) / (2.0*3.141592654*rol2H*0.65409);
              fH3 = (F1 - cos(2.0*3.141592654*rol3H)) / (2.0*3.141592654*rol3H*0.65409);
            }

            fH1 = aniso*fH1/sqrt(xLb1) + abs(aniso-1.0)*fl1;
            fH2 = aniso*fH2/sqrt(xLb2) + abs(aniso-1.0)*fl2;
            fH3 = aniso*fH3/sqrt(xLb3) + abs(aniso-1.0)*fl3;

            epsik1 = m_block->m_stgEddies[n][3];
            epsik2 = m_block->m_stgEddies[n][4];
            epsik3 = m_block->m_stgEddies[n][5];

            switch(m_block->m_stgShapeFunction) {
            case 2: {
              help4 += sqrt(Vb/nran)*epsik1*fl1*fl2*fl3;
              help5 += sqrt(Vb/nran)*epsik2*fl1*fl2*fl3;
              help6 += sqrt(Vb/nran)*epsik3*fl1*fl2*fl3;
              break;
            }
            case 3: {
              help4 += sqrt(Vb/nran)*epsik1*fl1;
              help5 += sqrt(Vb/nran)*epsik2*fl2;
              help6 += sqrt(Vb/nran)*epsik3*fl3;
              break;
            }
            case 4:
            default: {
              help4 += sqrt(Vb/nran)*epsik1*fl1*fl2*fH3;
              help5 += sqrt(Vb/nran)*epsik2*fl1*fl2*fH3;
              help6 += sqrt(Vb/nran)*epsik3*fl1*fH2*fl3;
            }
            }
          }

          /***************** End of L *****************/


          /*********************** M *****************/
          // Use Cholesky-Trafo to scale random fluctuations
          /*********************** M *****************/

          help1 = help4*a11; //Fluctuation u'
          help2 = help4*a21 + help5*a22; //Fluctuation v'
          help3 = help4*a31 + help5*a32 + help6*a33; //Fluctuation w'

          velmax = sqrt(umax*umax + vmax*vmax + wmax*wmax);

          ufluc = min(max(help1,-0.3*velmax),0.3*velmax);
          vfluc = min(max(help2,-0.3*velmax),0.3*velmax);
          wfluc = min(max(help3,-0.3*velmax),0.3*velmax);

          if(!std::isnan(abs(ufluc)))
            m_cells->stg_fq[FLUC_U][cellIdBC] += timsm*(ufluc - m_cells->stg_fq[FLUC_U][cellIdBC]);
          else
            cout << "ufluc isNAN, i: " << ii << " , j: " << j << " , k: " << k << " domainId: " << m_block->domainId() <<  endl;

          if(!std::isnan(abs(vfluc)))
            m_cells->stg_fq[FLUC_V][cellIdBC] += timsm*(vfluc - m_cells->stg_fq[FLUC_V][cellIdBC]);
          else
            cout << "vfluc isNAN, i: " << ii << " , j: " << j << " , k: " << k << " domainId: " << m_block->domainId() << endl;

          if(!std::isnan(abs(wfluc)))
            m_cells->stg_fq[FLUC_W][cellIdBC] += timsm*(wfluc - m_cells->stg_fq[FLUC_W][cellIdBC]);
          else
            cout << "wfluc isNAN, i: " << ii << " , j: " << j << " , k: " << k << " domainId: " << m_block->domainId() << endl;

          /***************** End of M *****************/
        }
      }
    }
  } //RKStep end if

  /////////////////////////////////////////////////////////
  ////////////// APPLY TO BC //////////////////////////////
  /////////////////////////////////////////////////////////
  //Now comes the BC stuff that we need to do every RK step
  for(ZFSId j=start[1]; j<end[1] ; j++) {
    for(ZFSId k=start[2]; k<end[2]; k++) {
      const ZFSId cellIdG1=cellIndex(m_noGhostLayers-1,j,k);
      IBC = cellIndexBC(m_noGhostLayers-1,j,k);
      const ZFSId cellIdA1=cellIndex(m_noGhostLayers,j,k);

      ZFSFloat dxidx = m_cells->cellMetrics[cellIdA1][0];
      ZFSFloat dxidy = m_cells->cellMetrics[cellIdA1][1];
      ZFSFloat dxidz = m_cells->cellMetrics[cellIdA1][2];

      ZFSFloat gradxi = -F1 / sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz);

      ZFSFloat dxHelp = dxidx*gradxi;
      ZFSFloat dyHelp = dxidy*gradxi;
      ZFSFloat dzHelp = dxidz*gradxi;

      const ZFSFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];
      const ZFSFloat pBC = m_cells->pvariables[PV->P][cellIdG1];
      const ZFSFloat fRhoBC = F1/rhoBC;
      const ZFSFloat aBC = sqrt(gamma*pBC*fRhoBC);
      const ZFSFloat uBC = m_cells->pvariables[PV->U][cellIdG1];
      const ZFSFloat vBC = m_cells->pvariables[PV->V][cellIdG1];
      const ZFSFloat wBC = m_cells->pvariables[PV->W][cellIdG1];

      const ZFSFloat maBC = (dxHelp*uBC +
                             dyHelp*vBC +
                             dzHelp*wBC)/aBC;

      /**This is where the fluctuations should be added to velocity!**/

      //get mean values from the rans
      const ZFSFloat rhoRANS = m_cells->stg_fq[PV->RHO][IBC];
      const ZFSFloat uRANS = m_cells->stg_fq[PV->U][IBC];
      const ZFSFloat vRANS = m_cells->stg_fq[PV->V][IBC];
      const ZFSFloat wRANS = m_cells->stg_fq[PV->W][IBC];
      const ZFSFloat pRANS = m_cells->stg_fq[PV->P][IBC];

      //fluctuation values from the STG
      const ZFSFloat u_prime = m_cells->stg_fq[FLUC_U][IBC];
      const ZFSFloat v_prime = m_cells->stg_fq[FLUC_V][IBC];
      const ZFSFloat w_prime = m_cells->stg_fq[FLUC_W][IBC];

      //superpose onto mean RANS variables
      const ZFSFloat uSTG = max(uRANS + u_prime, epsl);
      const ZFSFloat vSTG = vRANS + v_prime;
      const ZFSFloat wSTG = wRANS + w_prime;

      //compute correct density
      const ZFSFloat u9a = PV->UInfinity;
      const ZFSFloat u9ff = u_prime;
      const ZFSFloat alok = sqrt(gamma*PV->PInfinity/CV->rhoInfinity);
      const ZFSFloat flucc = u9ff/u9a*POW2((PV->UInfinity/alok))*gammaMinusOne*m_cells->stg_fq[PV->RHO][IBC];
      const ZFSFloat zdir = flucc / max(fabs(flucc), 0.0000001);
      const ZFSFloat rhoSTG = rhoRANS + zdir*min(fabs(flucc),0.1*rhoRANS);

      //get field values inside the integration domain
      const ZFSFloat pField = m_cells->pvariables[PV->P][cellIdA1];
      const ZFSFloat rhoField = m_cells->pvariables[PV->RHO][cellIdA1];
      const ZFSFloat uField = m_cells->pvariables[PV->U][cellIdA1];
      const ZFSFloat vField = m_cells->pvariables[PV->V][cellIdA1];
      const ZFSFloat wField = m_cells->pvariables[PV->W][cellIdA1];
      const ZFSFloat aField = sqrt(gamma*pField/rhoField);

      /////////////////////////////////////////////////
      //////////// SUBSONIC PART //////////////////////
      /////////////////////////////////////////////////
      const ZFSFloat pSub = F1B2 * (pField + pRANS + rhoField * aField *
                                    (  + dxHelp*(uField - uSTG)
                                       + dyHelp*(vField - vSTG)
                                       + dzHelp*(wField - wSTG)));
      const ZFSFloat rhoSub = rhoSTG + (pSub - pRANS)/(POW2(aField));
      const ZFSFloat rhoSubHelp = (pSub - pRANS) / (rhoField * aField);

      //Multiply velocities with density
      const ZFSFloat uSub = uSTG + dxHelp*rhoSubHelp;
      const ZFSFloat vSub = vSTG + dyHelp*rhoSubHelp;
      const ZFSFloat wSub = wSTG + dzHelp*rhoSubHelp;

      /////////////////////////////////////////////////
      //////////// SUPERSONIC PART ////////////////////
      /////////////////////////////////////////////////
      const ZFSFloat rhoSup = rhoSTG;
      const ZFSFloat uSup = uSTG;
      const ZFSFloat vSup = vSTG;
      const ZFSFloat wSup = wSTG;
      const ZFSFloat pSup = pRANS;

      //////////////////////////////////////////////////
      /////////// SUB/SUP INTERPOLATION ////////////////
      //////////////////////////////////////////////////

      //by default the subsonic formulation is used
      //switch on "stgSubSup" to get the mixed formulation
      //or "stgSupersonic" to use the pure supersonic formulation
      ZFSFloat xSub = F1;
      ZFSFloat xSup = F0;

      if(m_block->m_stgSubSup) {
        const ZFSFloat maBCAbs = fabs(maBC);
        const ZFSFloat alpha = 14.0;
        const ZFSFloat b = 0.95;
        const ZFSFloat count = alpha*(maBCAbs-b);
        const ZFSFloat denom = (F1-0.99*b)*maBCAbs+b;
        const ZFSFloat ratio = count/denom;
        const ZFSFloat wfun = F1B2*(F1+tanh(ratio)/tanh(alpha));

        xSub = fabs(wfun-F1);
        xSup = fabs(wfun);
      } else if (m_block->m_stgSupersonic) {
        xSub = F0;
        xSup = F1;
      }

      m_cells->pvariables[PV->RHO][cellIdG1] = rhoSub*xSub + rhoSup*xSup;
      m_cells->pvariables[PV->U][cellIdG1] = uSub*xSub + uSup*xSup;
      m_cells->pvariables[PV->V][cellIdG1] = vSub*xSub + vSup*xSup;
      m_cells->pvariables[PV->W][cellIdG1] = wSub*xSub + wSup*xSup;
      m_cells->pvariables[PV->P][cellIdG1] = pSub*xSub + pSup*xSup;


      //////////////////////////////////////////////////
      //////// EXTRAPOLATE TO SECOND GC ////////////////
      //////////////////////////////////////////////////
      const ZFSId cellIdG2 = cellIndex(0,j,k);
      //extrapolate into second ghost cell
      for(ZFSId var=0; var<PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdG2] = F2*m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
      }
    }
  }
}
//<marian




/** Rescaling Boundary Conditions
 *  BC2500 and BC2501 a combined rescaling boundary condition
 *  and can only be used together. The 2501 is the recycling
 *  station from where the values are taken and 2500 is the
 *  inflow plane where the rescaled values are prescribed.
 *
 */

//Inlet station
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2500(ZFSId bcId){
  (void) bcId;
  // Communication between the Recycling and Inlet station.
  //___COMMGr(oup)______________________
  //|______________|_________________| |
  //|Inlet station |Recycling station| |
  //|              |                 | |
  //|______________|_________________| |
  //|__________________________________|
  //
  // infographic of the Groups

  MPI_Comm rescalingCommGrComm = *m_block->m_rescalingCommGrComm;
  ZFSInt rescalingCommGrRootGlobal = *m_block->m_rescalingCommGrRootGlobal;

  // ZFSFloat delta_inMax = delta_in+2.5*delta_in; //limit the integration height
  const ZFSFloat gamma = m_block->m_gamma;
  const ZFSFloat gammaMinusOne=gamma-F1;

  const ZFSFloat rescalEPS = pow(10,-16.0);
  const ZFSFloat alpha= 4.0;
  const ZFSFloat b=0.2;
  const ZFSFloat rc= pow(m_block->m_Pr,F1B3);
  const ZFSFloat ctema = F1B2*gammaMinusOne*POW2(m_block->m_Ma)*rc;
  const ZFSFloat maxIntegrationHeight = 2.0*m_rescalingBLT;

  //van Driest constant & transformed velocity
  const ZFSFloat b_vd = sqrt(ctema/(F1+ctema));
  const ZFSFloat uvd8 = PV->UInfinity * asin(b_vd)/b_vd;

  //compute the momentum thickness at the inlet (i==x, j==y, k==z)
  const ZFSId i = m_noGhostLayers - 1;


  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE THETA AND EXCHANGE /////////////////
  //////////////////////////////////////////////////////////////

  //allocate space in k direction (all k-Cells)
  ZFSFloatScratchSpace thetaLocal(2, __CALLING_FUNCTION__, "thetaLocalIn");
  thetaLocal.fill(F0); //initialize scratch space

  ZFSFloatScratchSpace thetaGlobal(2, __CALLING_FUNCTION__, "thetaGlobalIn");
  thetaGlobal.fill(F0);

  //the offest position in k-direction is the offset
  const ZFSInt thetaLocalOffset= m_block->m_nOffsetCells[0];

  //compute the local moment thickness j=direction of integration
  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; ++k) {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; ++j) {
      const ZFSId cellId = cellIndex(i,j,k);
      const ZFSId pointIdM1 = getPointIdFromCell(i,j,k);
      const ZFSId pointIdP1 = getPointIdfromPoint(pointIdM1,0,1,0);

      if(m_coordinates[1][pointIdM1] > maxIntegrationHeight) {
        continue;
      }

      const ZFSFloat urat= m_cells->pvariables[PV->U][cellId]/PV->UInfinity;
      const ZFSFloat momThick=(m_cells->pvariables[PV->U][cellId]*m_cells->pvariables[PV->RHO][cellId]*fabs(F1-urat))/(CV->rhoUInfinity);

      //integrate normal to the wall
      const ZFSFloat ydist= m_coordinates[1][pointIdP1] - m_coordinates[1][pointIdM1];
      thetaLocal(0) += momThick*ydist;
    }
  }

  MPI_Allreduce(thetaLocal.begin(),thetaGlobal.begin(),2, MPI_DOUBLE, MPI_SUM,rescalingCommGrComm);

  //determine spanwise average: we assume equally spaced in z-direction
  for(ZFSId ii=0; ii<2; ++ii){
    thetaGlobal(ii)/=(m_block->m_totalGridBlockDim[0][0]-1);
  }

  if(globalTimeStep%50 == 0 && m_block->m_RKStep == 0 && m_block->domainId() == rescalingCommGrRootGlobal) {
    cout << "ThetaInflow: " << thetaGlobal(0) <<" ThetaRecyclingStation: " << thetaGlobal(1) << endl;

    FILE* f_channel;
    f_channel = fopen("./theta_inflow.dat", "a+");
    fprintf(f_channel, "%d", globalTimeStep);
    fprintf(f_channel, " %f", m_block->m_physicalTime);
    fprintf(f_channel, " %f", m_block->m_time);
    fprintf(f_channel, " %f", m_block->m_timeStep);
    fprintf(f_channel, " %f", thetaGlobal[0]);
    fprintf(f_channel, " %f", thetaGlobal[1]);
    fprintf(f_channel, "\n");
    fclose(f_channel);
  }

  //////////////////////////////////////////////////////////////
  ///////////////// EXCHANGE WALL PROPERTIES ///////////////////
  //////////////////////////////////////////////////////////////
  const ZFSId noVar = 2; //for more variables if wanted
  ZFSFloatScratchSpace wallPropertiesLocal(m_block->m_totalGridBlockDim[0][0]-1,noVar, __CALLING_FUNCTION__, "wallPropertiesLocalInlet");
  ZFSFloatScratchSpace wallProperties(m_block->m_totalGridBlockDim[0][0]-1,noVar, __CALLING_FUNCTION__, "wallPropertiesInlet");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  MPI_Allreduce(wallPropertiesLocal.begin(),wallProperties.begin(),(*m_block->m_totalGridBlockDim[0]-1)*noVar, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  //////////////////////////////////////////////////////////////
  ///////////////// GAMS, UTAUIN, INNER OUTER COORD ////////////
  //////////////////////////////////////////////////////////////
  ZFSFloatScratchSpace utauIn(m_nCells[0], __CALLING_FUNCTION__, "u_tauIn");
  ZFSFloatScratchSpace gams(m_nCells[0], __CALLING_FUNCTION__, "gams");

  ZFSId kStart = 0;
  ZFSId kEnd = m_nCells[0];

  if(m_block->m_nOffsetCells[0]==0) {
    kStart = m_noGhostLayers;
  }
  if(m_block->m_nOffsetCells[0]+m_block->m_nActiveCells[0]==m_block->m_totalGridBlockDim[0][0]-1) {
    kEnd = m_nCells[0]-m_noGhostLayers;
  }

  for(ZFSId k=kStart; k<kEnd; ++k) {
    const ZFSFloat utauRe=wallProperties(thetaLocalOffset+(k-m_noGhostLayers),0);

    // estimate the friction velocity at the inlet
    // according to standard power law approximations
    // utau_in = utau_re*(theta_re/theta_in)**(1/2*(n-1))
    // where theta is the momentum thickness
    // see Thomas S. Lund, p241

    //when take into account the variance of wall density
    //utau_in = utau_re*(rho_wall_re/rho_wall_in)**0.5
    //* (theta_re/theta_in)**(1/2*(n-1))
    //here n = 5

    gams(k) =pow(thetaGlobal(1)/fabs(thetaGlobal(0)),F1B8);
    utauIn(k)= utauRe*min(max(gams(k),F1),2.5);
  }

  ZFSFloatScratchSpace coordInInner(m_nCells[0]*m_nCells[1], __CALLING_FUNCTION__, "coordInInner");
  ZFSFloatScratchSpace coordInOuter(m_nCells[0]*m_nCells[1], __CALLING_FUNCTION__, "coordInOuter");

  for(ZFSId k=0; k<m_nCells[0]; ++k){
    for(ZFSId j=0; j<m_nCells[1]; ++j){
      const ZFSId cellId=cellIndex(i,j,k);
      const ZFSId faceId=j+k*m_nCells[1];
      const ZFSFloat rho=m_cells->pvariables[PV->RHO][cellId];
      const ZFSFloat frho=F1/rho;
      const ZFSFloat p = m_cells->pvariables[PV->P][cellId];
      const ZFSFloat temp=p*gamma*frho;
      const ZFSFloat mu=zfsSUTHERLANDLAW(temp);

      coordInInner(faceId)=utauIn(k)*rho*m_cells->coordinates[1][cellId]/(mu*sqrt(m_block->m_Re0));
      coordInOuter(faceId)=m_cells->coordinates[1][cellId]*rho/(m_rescalingBLT*CV->rhoInfinity);
    }
  }

  //////////////////////////////////////////////////////////////
  ///////////////// NOW EXCHANGE VAR SLICE /////////////////////
  //////////////////////////////////////////////////////////////
  const ZFSId noVariables = 6;
  const ZFSId totalCells[2]={m_block->m_totalGridBlockDim[0][0]-1, m_block->m_totalGridBlockDim[0][1]};
  ZFSFloatScratchSpace varSliceLocal(noVariables, totalCells[0]*totalCells[1], __CALLING_FUNCTION__, "varSliceLocal");
  ZFSFloatScratchSpace varSlice(noVariables, totalCells[0]*totalCells[1], __CALLING_FUNCTION__, "varSlice");

  //we are at the inlet, only fill with zeros
  varSlice.fill(F0);
  varSliceLocal.fill(F0);

  MPI_Allreduce(varSliceLocal.begin(),varSlice.begin(),noVariables*totalCells[0]*totalCells[1], MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);


  ///////////////////////////////////////////////////////////////
  ///////////////// RESCALING ///////////////////////////////////
  ///////////////////////////////////////////////////////////////

  ZFSId jStart = 0;
  ZFSId jEnd = m_nCells[1];

  if(m_block->m_nOffsetCells[1]==0) {
    jStart = m_noGhostLayers;
  }
  if(m_block->m_nOffsetCells[1]+m_block->m_nActiveCells[1]==m_block->m_totalGridBlockDim[0][1]-1) {
    jEnd = m_nCells[1]-m_noGhostLayers;
  }

  for(ZFSId k=kStart; k<kEnd; ++k) {
    const ZFSFloat ctem1= (F1+ctema)*(F1-POW2(gams(k)));
    const ZFSFloat ctem2= F2*ctema*gams(k)*(F1-gams(k));
    const ZFSFloat ctem3= (F1-gams(k))*(F1+gams(k)+F2*ctema*gams(k));

    for(ZFSId j=jStart; j<jEnd; ++j) {
      const ZFSId faceId = j+k*m_nCells[1];
      // const ZFSId faceIdM1 = (j-1)+k*m_nCells[1];
      const ZFSId cellId = cellIndex(i,j,k);

      if(coordInOuter(faceId) < 1.05){
        ZFSFloat uInner=F0,vInner=F0,wInner=F0,TInner=F0;
        ZFSFloat uOuter=F0,vOuter=F0,wOuter=F0,TOuter=F0;
        const ZFSFloat count = alpha*(coordInOuter(faceId) -b);
        const ZFSFloat denom = (F1-F2*b)*coordInOuter(faceId) +b;
        const ZFSFloat ratio = count/denom;
        const ZFSFloat wfun  = F1B2*(F1+tanh(ratio)/tanh(alpha));

        for(ZFSId jj=0; jj<totalCells[1]-1; ++jj) {
          const ZFSId localId = jj + (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];
          const ZFSId localIdP1 = (jj+1) + (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];

          const ZFSFloat yInnerRe   =  varSlice(4,localId);
          const ZFSFloat yInnerReP1 =  varSlice(4,localIdP1);

          if((yInnerRe-coordInInner(faceId)) < rescalEPS && yInnerReP1>coordInInner(faceId)) {
            const ZFSFloat dy1 = coordInInner(faceId)- yInnerRe;
            const ZFSFloat dy2 = yInnerReP1 - coordInInner(faceId);
            const ZFSFloat dy  = yInnerReP1 - yInnerRe;

            const ZFSFloat u   = varSlice(0,localId);
            const ZFSFloat uP1 = varSlice(0,localIdP1);
            const ZFSFloat v   = varSlice(1,localId);
            const ZFSFloat vP1 = varSlice(1,localIdP1);
            const ZFSFloat w   = varSlice(2,localId);
            const ZFSFloat wP1 = varSlice(2,localIdP1);
            const ZFSFloat t   = varSlice(3,localId);
            const ZFSFloat tP1 = varSlice(3,localIdP1);
            uInner= (uP1*dy1+u*dy2)/dy;
            vInner= (vP1*dy1+v*dy2)/dy;
            wInner= (wP1*dy1+w*dy2)/dy;
            TInner= (tP1*dy1+t*dy2)/dy;
          }
        }

        //outer region
        for(ZFSId jj=0; jj<totalCells[1]-1; ++jj) {
          const ZFSId localId = jj + (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];
          const ZFSId localIdP1 = (jj+1) + (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];

          const ZFSFloat yOuterRe   =  varSlice(5,localId);
          const ZFSFloat yOuterReP1 =  varSlice(5,localIdP1);

          if((yOuterRe-coordInOuter(faceId))<rescalEPS && yOuterReP1>coordInOuter(faceId)) {
            const ZFSFloat dy1= coordInOuter(faceId) - yOuterRe;
            const ZFSFloat dy2= yOuterReP1 - coordInOuter(faceId);
            const ZFSFloat dy = yOuterReP1 - yOuterRe;

            const ZFSFloat u     = varSlice(0,localId);
            const ZFSFloat uP1   = varSlice(0,localIdP1);
            const ZFSFloat v     = varSlice(1,localId);
            const ZFSFloat vP1   = varSlice(1,localIdP1);
            const ZFSFloat w     = varSlice(2,localId);
            const ZFSFloat wP1   = varSlice(2,localIdP1);
            const ZFSFloat t     = varSlice(3,localId);
            const ZFSFloat tP1   = varSlice(3,localIdP1);
            uOuter= (uP1*dy1+u*dy2)/dy;
            vOuter= (vP1*dy1+v*dy2)/dy;
            wOuter= (wP1*dy1+w*dy2)/dy;
            TOuter= (tP1*dy1+t*dy2)/dy;
          }
        }

        const ZFSFloat TInnerA= POW2(gams(k))*TInner+ctem1*PV->TInfinity;
        const ZFSFloat TOuterA= POW2(gams(k))*TOuter -(ctem2*(uOuter/PV->UInfinity)-ctem3)*PV->TInfinity;

        //van Driest transformation
        const ZFSFloat uvdInner= PV->UInfinity * asin(b_vd *uInner/PV->UInfinity)/b_vd;
        const ZFSFloat uvdOuter= PV->UInfinity * asin(b_vd *uOuter/PV->UInfinity)/b_vd;

        //scaling of transformed inner and outer velocities
        uInner = gams(k)*uvdInner;
        uOuter = gams(k)*uvdOuter+(F1-gams(k))*uvd8;
        uInner= PV->UInfinity*sin(b_vd*uInner/PV->UInfinity)/b_vd;
        uOuter= PV->UInfinity*sin(b_vd*uOuter/PV->UInfinity)/b_vd;

        const ZFSFloat pres = PV->PInfinity;
        const ZFSFloat uMean= uInner*(F1-wfun)+uOuter*wfun;
        const ZFSFloat vMean= vInner*(F1-wfun)+vOuter*wfun;
        const ZFSFloat wMean= (wInner*(F1-wfun)+wOuter*wfun)*gams(k);
        const ZFSFloat tMean= TInnerA*(F1-wfun)+TOuterA*wfun;
        const ZFSFloat rhoIn =gamma*pres/tMean;

        // //clebanoff factor is optional
        // const ZFSFloat clebf = 6.1;
        // const ZFSFloat blt   = m_rescalingBLT;
        // const ZFSFloat cleb  = F1/(F1+pow((m_cells->coordinates[1][cellId]/(clebf*blt)), 6.0));

        m_cells->pvariables[PV->RHO][cellId]=rhoIn;
        m_cells->pvariables[PV->U][cellId]=uMean;
        m_cells->pvariables[PV->V][cellId]=vMean;
        m_cells->pvariables[PV->W][cellId]=wMean;
        m_cells->pvariables[PV->P][cellId]= pres;
      }else{
        // if(!edgePointIsSet(k)) {
        //   edgePointJ(k) = j;
        //   edgePointIsSet(k) = 1;
        // }

        const ZFSFloat pres = PV->PInfinity;
        const ZFSFloat rhoIn = gamma*pres/PV->TInfinity;

        const ZFSFloat uMean = PV->UInfinity;
        const ZFSFloat vMean = PV->VInfinity;
        const ZFSFloat wMean = PV->WInfinity;

        m_cells->pvariables[PV->RHO][cellId]=rhoIn;
        m_cells->pvariables[PV->U][cellId]=uMean;
        m_cells->pvariables[PV->V][cellId]=vMean;
        m_cells->pvariables[PV->W][cellId]=wMean;
        m_cells->pvariables[PV->P][cellId]= pres;
      }
    }
  }

  for(ZFSId k=kStart; k<kEnd; ++k) {
    for(ZFSId j=0; j<m_nCells[1]; ++j) {
      //extrapolation for second GC
      const ZFSId cellId = cellIndex(1,j,k);
      const ZFSId cellIdM1 = cellIndex(0,j,k);
      const ZFSId cellIdadj = cellIndex(2,j,k);

      for(ZFSId var=0; var<PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdM1] = 2.0*m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
      }
    }
  }
}


//Recycling station
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2501(ZFSId bcId){
  ZFSInt* start = m_physicalBCMap[bcId]->start1;

  MPI_Comm rescalingCommGrComm = *m_block->m_rescalingCommGrComm;

  const ZFSFloat gamma = m_block->m_gamma;

  //things to move to init or elsewhere
  ZFSFloat F727=72.0/7.0;
  const ZFSId i= start[0]; //position at which recycle is taken
  const ZFSFloat yWall=F0; //this has been fixed else method does not work
  const ZFSFloat maxIntegrationHeight = 2.0*m_rescalingBLT;

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE THETA AND EXCHANGE /////////////////
  //////////////////////////////////////////////////////////////

  //thetaLocal.fill(F0); //initialize scratch space to zero // only for parallel use
  ZFSFloatScratchSpace thetaLocal(2, __CALLING_FUNCTION__, "thetaLocalRe");
  ZFSFloatScratchSpace thetaGlobal(2, __CALLING_FUNCTION__, "thetaGlobalRe");
  thetaLocal.fill(F0);
  thetaGlobal.fill(F0);

  //the offest position in k-direction is the offset
  const ZFSInt thetaLocalOffset= m_block->m_nOffsetCells[0];

  //compute the local moment thickness j=direction of integration
  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; ++k) {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; ++j) {
      const ZFSId cellId = cellIndex(i,j,k);
      const ZFSId pointIdM1 = getPointIdFromCell(i,j,k);
      const ZFSId pointIdP1 = getPointIdfromPoint(pointIdM1,0,1,0);

      if(m_coordinates[1][pointIdM1] > maxIntegrationHeight) {
        continue;
      }

      const ZFSFloat urat= m_cells->pvariables[PV->U][cellId]/PV->UInfinity;
      const ZFSFloat momThick=(m_cells->pvariables[PV->U][cellId]*m_cells->pvariables[PV->RHO][cellId]*fabs(F1-urat))/(CV->rhoUInfinity);

      //integrate normal to the wall
      const ZFSFloat ydist= m_coordinates[1][pointIdP1] - m_coordinates[1][pointIdM1];
      thetaLocal(1) += momThick*ydist;
    }
  }

 //communicate the Thickness across the plane
  MPI_Allreduce(thetaLocal.begin(), thetaGlobal.begin(), 2, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  //determine spanwise average: we assume equally spaced in z-direction
  for(ZFSId ii=0; ii<2; ++ii){
    thetaGlobal(ii)/=(m_block->m_totalGridBlockDim[0][0]-1);
  }

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE WALL PROPERTIES ///////
  //////////////////////////////////////////////////////////////

  const ZFSFloat delta = F727*thetaGlobal(1);
  const ZFSId noVar = 2; //for more variables if wanted
  const ZFSInt wallLocalOffset=m_block->m_nOffsetCells[1]; //Offset in j-direction

  ZFSFloatScratchSpace wallPropertiesLocal((m_block->m_totalGridBlockDim[0][0]-1),noVar, __CALLING_FUNCTION__, "wallPropertiesLocalRe");
  ZFSFloatScratchSpace wallProperties((m_block->m_totalGridBlockDim[0][0]-1),noVar, __CALLING_FUNCTION__,"wallPropertiesRe");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  //determine the wall stuff if wall is contained whithin the partition
  if(wallLocalOffset==0 && m_block->m_nActiveCells[1]>=m_noGhostLayers) {
    for(ZFSId k=m_noGhostLayers; k<m_block->m_nCells[0] - m_noGhostLayers; ++k ) {
      const ZFSId cellId = cellIndex(i, m_noGhostLayers, k);
      const ZFSId localId = m_block->m_nOffsetCells[0] + (k-m_noGhostLayers);
      const ZFSFloat rho = m_cells->pvariables[PV->RHO][cellId];
      const ZFSFloat p = m_cells->pvariables[PV->P][cellId];
      const ZFSFloat t=p*gamma/rho;
      const ZFSFloat mu=zfsSUTHERLANDLAW(t);
      const ZFSFloat uWall = fabs(m_cells->pvariables[PV->U][cellId]);
      const ZFSFloat ydist=m_cells->coordinates[1][cellId]-yWall;
      const ZFSFloat uTau= sqrt(uWall*mu/(ydist*rho));

      wallPropertiesLocal(localId,0)= uTau;
      wallPropertiesLocal(localId,1)= rho;
    }
  }

  MPI_Allreduce(wallPropertiesLocal.begin(),wallProperties.begin(),noVar*(m_block->m_totalGridBlockDim[0][0]-1), MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE VAR SLICE /////////////
  //////////////////////////////////////////////////////////////

  const ZFSId totalCells[2]={m_block->m_totalGridBlockDim[0][0]-1, m_block->m_totalGridBlockDim[0][1]};

  ZFSFloatScratchSpace varSliceLocal(6, totalCells[0]*totalCells[1], __CALLING_FUNCTION__, "varSliceLocal");
  ZFSFloatScratchSpace varSlice(6, totalCells[0]*totalCells[1], __CALLING_FUNCTION__, "varSlice");

  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; ++k) {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; ++j) {
      const ZFSId cellId=cellIndex(i,j,k);
      const ZFSId localId = (m_block->m_nOffsetCells[1]+(j-m_noGhostLayers) + 1) + (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];
      const ZFSFloat rho=m_cells->pvariables[PV->RHO][cellId];
      const ZFSFloat frho=F1/rho;
      const ZFSFloat p = m_cells->pvariables[PV->P][cellId];
      const ZFSFloat temp=p*gamma*frho;
      const ZFSFloat mu=zfsSUTHERLANDLAW(temp);
      const ZFSFloat uTauRe=wallProperties(thetaLocalOffset+(k-m_noGhostLayers),0);
      const ZFSFloat yIn=(m_cells->coordinates[1][cellId]-yWall)*uTauRe*rho/(mu*sqrt(m_block->m_Re0));
      const ZFSFloat yOut=(m_cells->coordinates[1][cellId]-yWall)*rho/(delta*CV->rhoInfinity);
      const ZFSFloat u=m_cells->pvariables[PV->U][cellId];
      const ZFSFloat v=m_cells->pvariables[PV->V][cellId];
      const ZFSFloat w=m_cells->pvariables[PV->W][cellId];

      //save the variables u,v,w,t,yI,yO
      varSliceLocal(0,localId) = u ;
      varSliceLocal(1,localId) = v ;
      varSliceLocal(2,localId) = w ;
      varSliceLocal(3,localId) = temp ;
      varSliceLocal(4,localId) = yIn ;
      varSliceLocal(5,localId) = yOut ;
    }

    //set first value at the wall manually
    if(m_block->m_nOffsetCells[0]==0) {
      const ZFSId localId = 0 + (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];
      const ZFSId localIdP1 = 1 + (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];
      varSliceLocal(0,localId) = 0.0 ; //u
      varSliceLocal(1,localId) = 0.0 ; //v
      varSliceLocal(2,localId) = 0.0 ; //w
      varSliceLocal(3,localId) = varSliceLocal(3,localIdP1); //t
      varSliceLocal(4,localId) = 0.0 ; //yIn
      varSliceLocal(5,localId) = 0.0 ; //yOut
    }
  }

  //communicate the slice
  MPI_Allreduce(varSliceLocal.begin(),varSlice.begin(),6*totalCells[0]*totalCells[1], MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  //participate in communication but only fill with zeros
  // ZFSFloatScratchSpace blEdgeVValueLocal(m_block->m_totalGridBlockDim[0][0]-1,__CALLING_FUNCTION__, "blEdgeVValueLocal");
  // ZFSFloatScratchSpace blEdgeVValueGlobal(m_block->m_totalGridBlockDim[0][0]-1,__CALLING_FUNCTION__, "blEdgeVValueLocal");
  // blEdgeVValueLocal.fill(F0);
  // blEdgeVValueGlobal.fill(F0);

  // MPI_Allreduce(blEdgeVValueLocal.begin(),blEdgeVValueGlobal.begin(), m_block->m_totalGridBlockDim[0][0]-1, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);
}


//Inlet station for RANS
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2510(ZFSId bcId){
  (void) bcId;

  cout.precision(8);
  // ZFSInt* start = m_physicalBCMap[bcId]->start1;
  // ZFSInt* end = m_physicalBCMap[bcId]->end1;
  // Communication between the Recycling and Inlet station.
  //___COMMGr(oup)______________________
  //|______________|_________________| |
  //|Inlet station |Recycling station| |
  //|              |                 | |
  //|______________|_________________| |
  //|__________________________________|
  //
  // infographic of the Groups

  MPI_Comm rescalingCommGrComm = *m_block->m_rescalingCommGrComm;
  ZFSInt rescalingCommGrRootGlobal = *m_block->m_rescalingCommGrRootGlobal;

  // ZFSFloat delta_inMax = delta_in+2.5*delta_in; //limit the integration height
  const ZFSFloat gamma = m_block->m_gamma;
  const ZFSFloat gammaMinusOne=gamma-F1;

  const ZFSFloat rescalEPS = pow(10,-16.0);
  const ZFSFloat alpha= 4.0;
  const ZFSFloat b=0.2;
  const ZFSFloat rc= pow(m_block->m_Pr,F1B3);
  const ZFSFloat ctema = F1B2*gammaMinusOne*POW2(m_block->m_Ma)*rc;

  //van Driest constant & transformed velocity
  const ZFSFloat b_vd = sqrt(ctema/(F1+ctema));
  const ZFSFloat uvd8 = PV->UInfinity * asin(b_vd)/b_vd;
  const ZFSId i = m_noGhostLayers - 1;

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE THETA AND EXCHANGE /////////////////
  //////////////////////////////////////////////////////////////

  //allocate space in k direction (all k-Cells)
  ZFSFloatScratchSpace thetaLocal(2, __CALLING_FUNCTION__, "thetaLocalIn");
  thetaLocal.fill(F0); //initialize scratch space

  ZFSFloatScratchSpace thetaGlobal(2, __CALLING_FUNCTION__, "thetaGlobalIn");
  thetaGlobal.fill(F0);

  //the offest position in k-direction is the offset
  const ZFSInt thetaLocalOffset= m_block->m_nOffsetCells[0];
  //compute the local moment thickness j=direction of integration
  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; ++k) {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; ++j) {
      const ZFSId cellId = cellIndex(i,j,k);
      const ZFSId pointIdM1 = getPointIdFromCell(i,j,k);
      const ZFSId pointIdP1 = getPointIdfromPoint(pointIdM1,0,1,0);

      const ZFSFloat urat= m_cells->pvariables[PV->U][cellId]/PV->UInfinity;
      const ZFSFloat momThick=(m_cells->pvariables[PV->U][cellId]*m_cells->pvariables[PV->RHO][cellId]*fabs(F1-urat))/(CV->rhoUInfinity);
      //integrate normal to the wall
      const ZFSFloat ydist= m_coordinates[1][pointIdP1] - m_coordinates[1][pointIdM1];
      thetaLocal(0) += momThick*ydist;
    }
  }

  MPI_Allreduce(thetaLocal.begin(),thetaGlobal.begin(),2, MPI_DOUBLE, MPI_SUM,rescalingCommGrComm);

  //determine spanwise average: we assume equally spaced in z-direction
  for(ZFSId ii=0; ii<2; ++ii){
    thetaGlobal(ii)/=(m_block->m_totalGridBlockDim[0][0]-1);
  }

  if(globalTimeStep%5==0&&m_block->m_RKStep==0&&m_block->domainId() == rescalingCommGrRootGlobal) {
    cout << m_block->domainId()<< " ThetaInflow " << thetaGlobal(0) <<" ThetaRecyclingStation " << thetaGlobal(1) << endl;

    FILE* f_channel;
    f_channel = fopen("./theta_inflow.dat", "a+");
    fprintf(f_channel, "%d", globalTimeStep);
    fprintf(f_channel, " %f", m_block->m_physicalTime);
    fprintf(f_channel, " %f", m_block->m_time);
    fprintf(f_channel, " %f", m_block->m_timeStep);
    fprintf(f_channel, " %f", thetaGlobal[0]);
    fprintf(f_channel, " %f", thetaGlobal[1]);
    fprintf(f_channel, "\n");
    fclose(f_channel);
  }

  //////////////////////////////////////////////////////////////
  ///////////////// EXCHANGE WALL PROPERTIES ///////////////////
  //////////////////////////////////////////////////////////////
  const ZFSId noWallProperties = 3;
  ZFSFloatScratchSpace wallPropertiesLocal(*m_block->m_totalGridBlockDim[0]-1,noWallProperties, __CALLING_FUNCTION__, "wallPropertiesLocalInlet");
  ZFSFloatScratchSpace wallProperties(*m_block->m_totalGridBlockDim[0]-1,noWallProperties, __CALLING_FUNCTION__, "wallPropertiesInlet");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  MPI_Allreduce(wallPropertiesLocal.begin(),wallProperties.begin(),(*m_block->m_totalGridBlockDim[0]-1)*noWallProperties, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  //////////////////////////////////////////////////////////////
  ///////////////// GAMS, UTAUIN, INNER OUTER COORD ////////////
  //////////////////////////////////////////////////////////////

  ZFSFloatScratchSpace utauIn(m_block->m_nActiveCells[0],2, __CALLING_FUNCTION__, "u_tauIn");
  ZFSFloatScratchSpace gams(m_block->m_nActiveCells[0],2, __CALLING_FUNCTION__, "gams");

  for(ZFSId k=0; k<m_block->m_nActiveCells[0]; ++k) {
    ZFSFloat utauRe=wallProperties(thetaLocalOffset+(k),0);

    // estimate the friction velocity at the inlet
    // according to standard power law approximations
    // utau_in = utau_re*(theta_re/theta_in)**(1/2*(n-1))
    // where theta is the momentum thichness
    // see Thomas S. Lund, p241

    //when take into account the variance of wall density
    //utau_in = utau_re*(rho_wall_re/rho_wall_in)**0.5
    //* (theta_re/theta_in)**(1/2*(n-1))
    //here n = 5
    const ZFSFloat n = 5.0;
    const ZFSFloat facc = F1/(F2*(n-F1));

    gams(k,i) =pow(thetaGlobal(1)/fabs(thetaGlobal(0)),facc);
    utauIn(k,i)= utauRe*min(max(gams(k,i),F1),2.5);
  }

  ZFSFloatScratchSpace coordInInner(m_nCells[0]*m_nCells[1],2, __CALLING_FUNCTION__, "coordInInner");
  ZFSFloatScratchSpace coordInOuter(m_nCells[0]*m_nCells[1],2, __CALLING_FUNCTION__, "coordInOuter");

  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; ++k){
    for(ZFSId j=0; j<m_nCells[1]; ++j){
      const ZFSId cellId=cellIndex(i,j,k);
      const ZFSId localId=j+(k-m_noGhostLayers)*m_nCells[1];
      const ZFSFloat rho=m_cells->pvariables[PV->RHO][cellId];
      const ZFSFloat frho=F1/rho;
      const ZFSFloat p = m_cells->pvariables[PV->P][cellId];
      const ZFSFloat temp=p*gamma*frho;
      const ZFSFloat mu=zfsSUTHERLANDLAW(temp);

      coordInInner(localId,i)=utauIn(k-m_noGhostLayers,i)*rho*m_cells->coordinates[1][cellId]/(mu*sqrt(m_block->m_Re0));
      coordInOuter(localId,i)=m_cells->coordinates[1][cellId]*rho/(m_rescalingBLT*CV->rhoInfinity);
    }
  }

  const ZFSId wallLocalOffset=m_block->m_nOffsetCells[1]; //Offset in j-direction
  ZFSFloatScratchSpace tempWallInletLocal((*m_block->m_totalGridBlockDim[0]-1), __CALLING_FUNCTION__, "tempWallInletLocal");
  ZFSFloatScratchSpace tempWallInletGlobal((*m_block->m_totalGridBlockDim[0]-1), __CALLING_FUNCTION__, "tempWallInletGlobal");
  tempWallInletLocal.fill(F0);
  tempWallInletGlobal.fill(F0);
  //determine the wall stuff if wall is contained whithin the partition
  if(wallLocalOffset==0 && m_block->m_nActiveCells[1]>=m_noGhostLayers) {
    for(ZFSId k=m_noGhostLayers; k<m_block->m_nCells[0] - m_noGhostLayers; ++k ) {
      const ZFSId cellId=cellIndex(1, m_noGhostLayers, k);
      const ZFSId localId = m_block->m_nOffsetCells[0] + (k-m_noGhostLayers);
      tempWallInletLocal(localId)= temperature(cellId);
    }
  }

  MPI_Allreduce(tempWallInletLocal.begin(),tempWallInletGlobal.begin(),(*m_block->m_totalGridBlockDim[0]-1), MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);


  //////////////////////////////////////////////////////////////
  ///////////////// NOW EXCHANGE VAR SLICE /////////////////////
  //////////////////////////////////////////////////////////////
  const ZFSId noVariables = PV->noVariables+1;
  ZFSId totalCells[2]={m_block->m_totalGridBlockDim[0][0]-1, m_block->m_totalGridBlockDim[0][1]-1};
  ZFSFloatScratchSpace varSliceLocal(noVariables, totalCells[0]*totalCells[1], __CALLING_FUNCTION__, "varSliceLocal");
  ZFSFloatScratchSpace varSlice(noVariables, totalCells[0]*totalCells[1], __CALLING_FUNCTION__, "varSlice");

  //we are at the inlet, only fill with zeros
  varSlice.fill(F0);
  varSliceLocal.fill(F0);

  MPI_Allreduce(varSliceLocal.begin(),varSlice.begin(),noVariables*totalCells[0]*totalCells[1], MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  if(globalTimeStep%5==0&&m_block->m_RKStep==0&&m_block->domainId() == rescalingCommGrRootGlobal) {
    cout << m_block->domainId()<< " ThetaInflow " << thetaGlobal(0) <<" ThetaRecyclingStation " << thetaGlobal(1) << endl;
  }

  ///////////////////////////////////////////////////////////////
  ///////////////// RESCALING ///////////////////////////////////
  ///////////////////////////////////////////////////////////////

  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; ++k) {
    const ZFSFloat ctem1= (F1+ctema)*(F1-POW2(gams(k-m_noGhostLayers,i)));
    const ZFSFloat ctem2= F2*ctema*gams(k-m_noGhostLayers,i)*(F1-gams(k-m_noGhostLayers,i));
    const ZFSFloat ctem3= (F1-gams(k-m_noGhostLayers,i))*(F1+gams(k-m_noGhostLayers,i)+F2*ctema*gams(k-m_noGhostLayers,i));


    ZFSId jStart = 0;
    ZFSId jEnd = m_nCells[1];

    if(m_block->m_nOffsetCells[1]==0) {
      jStart = m_noGhostLayers;
    }
    if(m_block->m_nOffsetCells[1]+m_block->m_nActiveCells[1]==m_block->m_totalGridBlockDim[0][1]-1) {
      jEnd = m_nCells[1]-m_noGhostLayers;
    }

    for(ZFSId j=jStart; j<jEnd; ++j) {
      const ZFSId faceId=j+(k-m_noGhostLayers)*m_nCells[1];
      const ZFSId cellId = cellIndex(i,j,k);

      if( coordInOuter(faceId,i) < 1.05){
        ZFSFloat uInner=F0,vInner=F0,wInner=F0, TInner=F0, mutInner=F0;
        ZFSFloat uOuter=F0,vOuter=F0,wOuter=F0, TOuter=F0, mutOuter=F0;
        const ZFSFloat count = alpha*(coordInOuter(faceId,i) -b);
        const ZFSFloat denom = (F1-F2*b)*coordInOuter(faceId,i) +b;
        const ZFSFloat ratio = count/denom;
        const ZFSFloat wfun  = F1B2*(F1+tanh(ratio)/tanh(alpha));

        for(ZFSId jj=0; jj<totalCells[1]-1; ++jj) {
          const ZFSId localId   =  jj+   (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];
          const ZFSId localIdP1 =  jj+1+   (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];

          const ZFSFloat yInnerRe   =  varSlice(4, localId);
          const ZFSFloat yInnerReP1 =  varSlice(4, localIdP1);

          if((yInnerRe-coordInInner(faceId,i)) < rescalEPS && yInnerReP1>coordInInner(faceId,i)) {
            const ZFSFloat dy1 = coordInInner(faceId,i)- yInnerRe;
            const ZFSFloat dy2 = yInnerReP1 - coordInInner(faceId,i);
            const ZFSFloat dy  = yInnerReP1 - yInnerRe;

            const ZFSFloat u     = varSlice(0,localId);
            const ZFSFloat uP1   = varSlice(0,localIdP1);
            const ZFSFloat v     = varSlice(1,localId);
            const ZFSFloat vP1   = varSlice(1,localIdP1);
            const ZFSFloat w     = varSlice(2,localId);
            const ZFSFloat wP1   = varSlice(2,localIdP1);
            const ZFSFloat t     = varSlice(3,localId);
            const ZFSFloat tP1   = varSlice(3,localIdP1);
            const ZFSFloat mut   = varSlice(6,localId);
            const ZFSFloat mutP1 = varSlice(6,localIdP1);
            uInner= (uP1*dy1+u*dy2)/dy;
            vInner= (vP1*dy1+v*dy2)/dy;
            wInner= (wP1*dy1+w*dy2)/dy;
            TInner= (tP1*dy1+t*dy2)/dy;
            mutInner = (mutP1*dy1+mut*dy2)/dy;
          }
        }

        //>marian: catch those cells that didn't pass the if-clause (where TInner is still zero)
        if(TInner < 0.5) {
          for(ZFSId jj=0; jj<totalCells[1]-1; ++jj){
            const ZFSId localId   =  jj+   (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];
            const ZFSId localIdP1 =  jj+1+   (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];

            const ZFSFloat yInnerRe   =  varSlice(4, localId);
            const ZFSFloat yInnerReP1 =  varSlice(4, localIdP1);
            const ZFSFloat diffPercent = (abs(yInnerRe - coordInInner(faceId,i)) / coordInInner(faceId,i))*100.0;

            if(( diffPercent < 10.0) && yInnerReP1>coordInInner(faceId,i)){
              const ZFSFloat dy1 = coordInInner(faceId,i)- yInnerRe;
              const ZFSFloat dy2 = yInnerReP1 - coordInInner(faceId,i);
              const ZFSFloat dy  = yInnerReP1 - yInnerRe;

              const ZFSFloat u     = varSlice(0,localId);
              const ZFSFloat uP1   = varSlice(0,localIdP1);
              const ZFSFloat v     = varSlice(1,localId);
              const ZFSFloat vP1   = varSlice(1,localIdP1);
              const ZFSFloat w     = varSlice(2,localId);
              const ZFSFloat wP1   = varSlice(2,localIdP1);
              const ZFSFloat t     = varSlice(3,localId);
              const ZFSFloat tP1   = varSlice(3,localIdP1);
              const ZFSFloat mut   = varSlice(6,localId);
              const ZFSFloat mutP1 = varSlice(6,localIdP1);
              uInner= (uP1*dy1+u*dy2)/dy;
              vInner= (vP1*dy1+v*dy2)/dy;
              wInner= (wP1*dy1+w*dy2)/dy;
              TInner= (tP1*dy1+t*dy2)/dy;
              mutInner = (mutP1*dy1+mut*dy2)/dy;
            }
          }
        }
        //<marian

        //outer region
        for(ZFSId jj=0; jj<totalCells[1]-1; ++jj) {
            const ZFSId localId = jj+   (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];
            const ZFSId localIdP1 = jj+1 +   (m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];

            const ZFSFloat yOuterRe   =varSlice(5, localId);
            const ZFSFloat yOuterReP1 =varSlice(5, localIdP1);

          if((yOuterRe-coordInOuter(faceId,i))<rescalEPS && yOuterReP1>coordInOuter(faceId,i)) {
            const ZFSFloat dy1= coordInOuter(faceId,i) - yOuterRe;
            const ZFSFloat dy2= yOuterReP1 - coordInOuter(faceId,i);
            const ZFSFloat dy = yOuterReP1 - yOuterRe;


            const ZFSFloat u     = varSlice(0,localId);
            const ZFSFloat uP1   = varSlice(0,localIdP1);
            const ZFSFloat v     = varSlice(1,localId);
            const ZFSFloat vP1   = varSlice(1,localIdP1);
            const ZFSFloat w     = varSlice(2,localId);
            const ZFSFloat wP1   = varSlice(2,localIdP1);
            const ZFSFloat t     = varSlice(3,localId);
            const ZFSFloat tP1   = varSlice(3,localIdP1);
            const ZFSFloat mut   = varSlice(6,localId);
            const ZFSFloat mutP1 = varSlice(6,localIdP1);
            uOuter= (uP1*dy1+u*dy2)/dy;
            vOuter= (vP1*dy1+v*dy2)/dy;
            wOuter= (wP1*dy1+w*dy2)/dy;
            TOuter= (tP1*dy1+t*dy2)/dy;
            mutOuter = (mutP1*dy1+mut*dy2)/dy;
         }
        }

        const ZFSFloat TInnerA= POW2(gams(k-m_noGhostLayers,i))*TInner+ctem1*PV->TInfinity;
        const ZFSFloat TOuterA= POW2(gams(k-m_noGhostLayers,i))*TOuter -(ctem2*(uOuter/PV->UInfinity)-ctem3)*PV->TInfinity;

        //van Driest transformation
        const ZFSFloat uvdInner= PV->UInfinity * asin(b_vd *uInner/PV->UInfinity)/b_vd;
        const ZFSFloat uvdOuter= PV->UInfinity * asin(b_vd *uOuter/PV->UInfinity)/b_vd;

        //scaling of transformed inner and outer velocities
        uInner = gams(k-m_noGhostLayers,i)*uvdInner;
        uOuter = gams(k-m_noGhostLayers,i)*uvdOuter+(F1-gams(k-m_noGhostLayers,i))*uvd8;

        uInner= PV->UInfinity*sin(b_vd*uInner/PV->UInfinity)/b_vd;
        uOuter= PV->UInfinity*sin(b_vd*uOuter/PV->UInfinity)/b_vd;

        const ZFSFloat pres = PV->PInfinity;

        const ZFSFloat uMean= uInner*(F1-wfun)+uOuter*wfun;
        const ZFSFloat vMean= vInner*(F1-wfun)+vOuter*wfun;
        const ZFSFloat wMean= (wInner*(F1-wfun)+wOuter*wfun)*gams(k-m_noGhostLayers,i);
        const ZFSFloat tMean= TInnerA*(F1-wfun)+TOuterA*wfun;
        const ZFSFloat rhoIn =gamma*pres/tMean;

        //turbulent viscosity
        const ZFSFloat tempWallInlet = tempWallInletGlobal(thetaLocalOffset+(k-m_noGhostLayers));
        const ZFSFloat tempWallRecycling = wallProperties(thetaLocalOffset+(k-m_noGhostLayers), 2);

        const ZFSFloat viscWallInlet =  zfsSUTHERLANDLAW(tempWallInlet);
        const ZFSFloat viscWallRecycling = zfsSUTHERLANDLAW(tempWallRecycling);
        const ZFSFloat thetaInlet = thetaGlobal(0);
        const ZFSFloat thetaRecycling = thetaGlobal(1);
        mutInner = mutInner*(viscWallInlet/viscWallRecycling);
        mutOuter = mutOuter*gams(k-m_noGhostLayers,i)*(thetaInlet/thetaRecycling);

        ZFSFloat mutMean = mutInner*(F1-wfun)+mutOuter*wfun;
        const ZFSFloat clebf = 6.6;
        const ZFSFloat blt   = m_rescalingBLT;
        const ZFSFloat cleb=F1/(F1+pow((m_cells->coordinates[1][cellId]/(clebf*blt)), 6.0));

        m_cells->pvariables[PV->RHO][cellId]=rhoIn*cleb;
        m_cells->pvariables[PV->U][cellId]=uMean*cleb;
        m_cells->pvariables[PV->V][cellId]=vMean*cleb;
        m_cells->pvariables[PV->W][cellId]=wMean*cleb;

        m_cells->pvariables[PV->P][cellId]= pres;
        m_cells->pvariables[PV->RANS_VAR[0]][cellId] = mutMean/rhoIn;

      }else{
        //const ZFSFloat pres = pressure(cellIndex(m_noGhostLayers,j,k)); //for supersonic PV->PInfinity
        const ZFSFloat pres = PV->PInfinity;
        const ZFSFloat rhoIn = gamma*pres/PV->TInfinity;

        const ZFSFloat uMean = PV->UInfinity;
        const ZFSFloat vMean = PV->VInfinity;
        const ZFSFloat wMean = PV->WInfinity;

        m_cells->pvariables[PV->RHO][cellId]= rhoIn;
        m_cells->pvariables[PV->U][cellId]= uMean;
        m_cells->pvariables[PV->V][cellId]= vMean;
        m_cells->pvariables[PV->W][cellId]= wMean;
        m_cells->pvariables[PV->P][cellId]= pres;
        m_cells->pvariables[PV->RANS_VAR[0]][cellId]=PV->ransInfinity[0];
      }

    }
  }

  for(ZFSId k=0; k<m_nCells[0]; ++k) {
    for(ZFSId j=0; j<m_nCells[1]; ++j) {
      //extrapolation for second GC
      const ZFSId cellId = cellIndex(1,j,k);
      const ZFSId cellIdM1 = cellIndex(0,j,k);
      const ZFSId cellIdadj = cellIndex(2,j,k);

      for(ZFSId var=0; var<PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdM1] = 2.0*m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
      }
    }
  }
}


//Recycling station for RANS
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2511(ZFSId bcId){
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  cout.precision(8);

  //ZFSInt *rescalingRoot = m_block->m_rescalingRoot;
  MPI_Comm rescalingCommGrComm = *m_block->m_rescalingCommGrComm;

  //scaling between delta0 and delta2
  const ZFSFloat F727=72.0/8.0;
  const ZFSId i= start[0]; //position at which recycle is taken
  const ZFSFloat yWall=F0; //this has been fixed else method does not works

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE THETA AND EXCHANGE /////////////////
  //////////////////////////////////////////////////////////////

  //thetaLocal.fill(F0); //initialize scratch space to zero // only for parallel use
  ZFSFloatScratchSpace thetaLocal(2, __CALLING_FUNCTION__, "thetaLocalRe");
  ZFSFloatScratchSpace thetaGlobal(2, __CALLING_FUNCTION__, "thetaGlobalRe");
  thetaLocal.fill(F0); //initialize scratch space
  thetaGlobal.fill(F0);

  //the offest position in k-direction is the offset
  const ZFSInt thetaLocalOffset= m_block->m_nOffsetCells[0];
  //compute the local moment thickness j=direction of integration
  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; ++k) {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; ++j) {
      const ZFSId cellId = cellIndex(i,j,k);
      const ZFSId pointIdM1 = getPointIdFromCell(i,j,k);
      const ZFSId pointIdP1 = getPointIdfromPoint(pointIdM1,0,1,0);

      const ZFSFloat urat= m_cells->pvariables[PV->U][cellId]/PV->UInfinity;
      const ZFSFloat momThick=(m_cells->pvariables[PV->U][cellId]*m_cells->pvariables[PV->RHO][cellId]*fabs(F1-urat))/(CV->rhoUInfinity);

      //integrate normal to the wall
      const ZFSFloat ydist= m_coordinates[1][pointIdP1] - m_coordinates[1][pointIdM1];
      thetaLocal(1) += momThick*ydist;
    }
  }


 //communicate the Thickness across the plane
  MPI_Allreduce(thetaLocal.begin(), thetaGlobal.begin(), 2 ,MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  thetaGlobal(1)/=(*m_block->m_totalGridBlockDim[0]-1); //we have now the averaged momentum thi

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE WALL PROPERTIES ///////
  //////////////////////////////////////////////////////////////

  const ZFSFloat delta = F727*thetaGlobal(1);
  const ZFSId noVar =3; //for more variables if wanted
  const ZFSId wallLocalOffset=m_block->m_nOffsetCells[1]; //Offset in j-direction

  ZFSFloatScratchSpace wallPropertiesLocal((*m_block->m_totalGridBlockDim[0]-1),noVar, __CALLING_FUNCTION__, "wallPropertiesLocalRe");
  ZFSFloatScratchSpace wallProperties((*m_block->m_totalGridBlockDim[0]-1),noVar, __CALLING_FUNCTION__,"wallPropertiesRe");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  //determine the wall stuff if wall is contained whithin the partition
  if(wallLocalOffset==0 && m_block->m_nActiveCells[1]>=m_noGhostLayers) {
    for(ZFSId k=m_noGhostLayers; k<m_block->m_nCells[0] - m_noGhostLayers; ++k ) {
      const ZFSId cellId=cellIndex(i, m_noGhostLayers, k);
      const ZFSFloat rho = m_cells->pvariables[PV->RHO][cellId];
      const ZFSFloat t=temperature(cellId);
      const ZFSFloat nu=zfsSUTHERLANDLAW(t);
      const ZFSFloat uWall = fabs(m_cells->pvariables[PV->U][cellId]);
      const ZFSFloat ydist=m_cells->coordinates[1][cellId]-yWall;
      const ZFSFloat uTau= sqrt(uWall*nu/(ydist*rho));
      const ZFSId localId = m_block->m_nOffsetCells[0] + (k-m_noGhostLayers);

      wallPropertiesLocal(localId,0)= uTau;
      wallPropertiesLocal(localId,1)= rho;
      wallPropertiesLocal(localId,2)= t;
    }
  }


  MPI_Allreduce(wallPropertiesLocal.begin(),wallProperties.begin(),noVar*(*m_block->m_totalGridBlockDim[0]-1), MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  ZFSFloatScratchSpace tempWallInletLocal((*m_block->m_totalGridBlockDim[0]-1), __CALLING_FUNCTION__, "tempWallInletLocal");
  ZFSFloatScratchSpace tempWallInletGlobal((*m_block->m_totalGridBlockDim[0]-1), __CALLING_FUNCTION__, "tempWallInletGlobal");
  tempWallInletLocal.fill(F0);
  tempWallInletGlobal.fill(F0);

  MPI_Allreduce(tempWallInletLocal.begin(),tempWallInletGlobal.begin(),(*m_block->m_totalGridBlockDim[0]-1), MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE VAR SLICE /////////////
  //////////////////////////////////////////////////////////////

  ZFSId totalCells[2]={m_block->m_totalGridBlockDim[0][0]-1, m_block->m_totalGridBlockDim[0][1]-1};

  const ZFSId noVariables = PV->noVariables+1;
  ZFSFloatScratchSpace varSliceLocal(noVariables, totalCells[0]*totalCells[1], __CALLING_FUNCTION__, "varSliceLocal");
  ZFSFloatScratchSpace varSlice(noVariables, totalCells[0]*totalCells[1], __CALLING_FUNCTION__, "varSlice");


  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; ++k) {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; ++j) {
      const ZFSId cellId=cellIndex(i,j,k);
      const ZFSFloat rho=m_cells->pvariables[PV->RHO][cellId];
      const ZFSFloat temp=temperature(cellId);
      const ZFSFloat nu=zfsSUTHERLANDLAW(temp);
      const ZFSFloat uTauRe=wallProperties(thetaLocalOffset+(k-m_noGhostLayers),0);
      const ZFSFloat yIn=(m_cells->coordinates[1][cellId]-yWall)*uTauRe*rho/(nu*sqrt(m_block->m_Re0));
      const ZFSFloat yOut=(m_cells->coordinates[1][cellId]-yWall)*rho/(delta*CV->rhoInfinity);
      const ZFSFloat u=m_cells->pvariables[PV->U][cellId];
      const ZFSFloat v=m_cells->pvariables[PV->V][cellId];
      const ZFSFloat w=m_cells->pvariables[PV->W][cellId];
      //>RANS
      const ZFSFloat nut = m_cells->pvariables[PV->RANS_VAR[0]][cellId];
      //<RANS

      const ZFSId localId = m_block->m_nOffsetCells[1]+(j-m_noGhostLayers)+(m_block->m_nOffsetCells[0]+(k-m_noGhostLayers))*totalCells[1];

      //save the variables u,v,w,t,yI,yO
      varSliceLocal(0,localId) = u ;
      varSliceLocal(1,localId) = v ;
      varSliceLocal(2,localId) = w ;
      varSliceLocal(3,localId) = temp ;
      varSliceLocal(4,localId) = yIn ;
      varSliceLocal(5,localId) = yOut ;
      varSliceLocal(6,localId) = nut ;
    }
  }

  //communicate the slice
  MPI_Allreduce(varSliceLocal.begin(),varSlice.begin(),noVariables*totalCells[0]*totalCells[1], MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);
}

/* laminar Poiseuille inflow
 *
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2020(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
      {
        for(ZFSId k=start[2]; k<end[2] ; k++) {
          for(ZFSId j=end[1]-1; j>=start[1] ; j--) {
            const ZFSId cellIdG1 = cellIndex(start[0]+1,j,k);
            const ZFSId cellIdG2 = cellIndex(start[0]+0,j,k);
            const ZFSId cellIdA1 = cellIndex(start[0]+2,j,k);
            const ZFSFloat y_max = F1;  // channel height

            const ZFSFloat x=m_cells->coordinates[0][cellIdG1];
            const ZFSFloat y=m_cells->coordinates[1][cellIdG1];
            const ZFSFloat pG1= PV->PInfinity - F3*(x+15.0)*zfsSUTHERLANDLAW(PV->TInfinity)*PV->UInfinity*POW2(F2/y_max)/m_block->m_Re0;

            m_cells->pvariables[PV->RHO][cellIdG1] = CV->rhoInfinity;
            m_cells->pvariables[PV->U][cellIdG1] = (-(F3/F2)*PV->UInfinity*(POW2(y-y_max/F2)-POW2(y_max/F2))/POW2(y_max/F2));
            m_cells->pvariables[PV->V][cellIdG1] = F0;
            m_cells->pvariables[PV->W][cellIdG1] = F0;
            m_cells->pvariables[PV->P][cellIdG1] = pG1;

            //extrapolate into second ghost cell
            for(ZFSId var=0; var<PV->noVariables; var++) {
              m_cells->pvariables[var][cellIdG2] = F2*m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
            }
          }
        }
        break;
      }
    default:
      {
        cout << "bc2020: face not implemented" << endl;
      }
    }
}

/** Prescribe given profile BC
 *
 *  Precribes a profile from the restart file
 *  extrapolate pressure from computational domain
 *  author: Marian Albers
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2600(ZFSId bcId){
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face)
  {
  case 0:
  {
    for(ZFSInt i = start[0]; i<end[0]; i++) {
      for(ZFSInt j = start[1]; j<end[1]; j++) {
        for(ZFSInt k = start[2]; k<end[2]; k++) {
          const ZFSId cellId=cellIndex(m_noGhostLayers-1-i,j,k);
          const ZFSId cellIdadj=cellIndex(m_noGhostLayers-i,j,k);

          //extrapolate pressure to ghost cells
          m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
        }
      }
    }
    break;
  }
  default:
  {
    zfsTerm(1, __CALLING_FUNCTION__, "Face not implemented");
    break;
  }
  }
}


/** Prescribe given profile BC
 *
 *  Precribes a profile from the restart file
 *  extrapolate pressure from computational domain
 *  author: Marian Albers 2016
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2601(ZFSId bcId){
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  ZFSFloat theta100000 = 0.002065747418;
  ZFSFloat rexPosQuad = 0.014452;
  ZFSFloat gammaEpsilon = m_block->m_bc2601GammaEpsilon;

  if(m_block->m_rans) {
    theta100000 = 0.002171197777811;
    rexPosQuad = 0.14452;
  }

  const ZFSFloat lengthFactor = theta100000*rexPosQuad;

  ZFSFloatScratchSpace effConst(m_2601noPos, m_2601noCoeff, __CALLING_FUNCTION__, "effConst");
  effConst.fill(F0);


  switch(m_physicalBCMap[bcId]->face)
    {
    case 2:
      {
        for(ZFSInt k = start[2]; k<end[2]; k++) {

          ////////////////////////////////////////////////////////////
          /////////// TIME INTERPOLATION ACTUATED SOLUTION ///////////
          ////////////////////////////////////////////////////////////
          //manipulate time with spanwise coordinate
          const ZFSId spanCell = cellIndex(m_noGhostLayers,m_noGhostLayers,k);
          const ZFSFloat spanCoord = m_cells->coordinates[2][spanCell];
          const ZFSFloat waveSpeed = 0.015795;
          const ZFSFloat timeShift = spanCoord/waveSpeed;
          const ZFSFloat t = m_block->m_time + m_block->m_timeStep*m_block->m_RKalpha[m_block->m_RKStep] - timeShift;
        
          //////////////////////////////////////////
          /////////// TIME INTERPOLATION ///////////
          //////////////////////////////////////////

          //interpolate between time values
          for(ZFSInt pos=0; pos<m_2601noPos; pos++) {
            effConst(pos,0) = m_2601streamwisePos[pos];
            for(ZFSInt iter=0; iter<m_2601noSteps-1; iter++) {
              if(m_2601effConst[0][pos*m_2601noCoeff+0] >= t) {
                effConst(pos,1) = m_2601effConst[0][pos*m_2601noCoeff+1];
                effConst(pos,2) = m_2601effConst[0][pos*m_2601noCoeff+2];
                effConst(pos,3) = m_2601effConst[0][pos*m_2601noCoeff+3];
                effConst(pos,4) = m_2601effConst[0][pos*m_2601noCoeff+4];
                break;
              } else if(m_2601effConst[iter][pos*m_2601noCoeff+0] < t &&
                        t <= m_2601effConst[iter+1][pos*m_2601noCoeff+0]) {
                for(ZFSInt var=1; var<5; var++) {
                  //linear interpolation between two time values
                  effConst(pos,var) = m_2601effConst[iter][pos*m_2601noCoeff+var] +
                    (t - m_2601effConst[iter][pos*m_2601noCoeff+0])/
                    (m_2601effConst[iter+1][pos*m_2601noCoeff+0]-m_2601effConst[iter][pos*m_2601noCoeff+0])*(m_2601effConst[iter+1][pos*m_2601noCoeff+var]-m_2601effConst[iter][pos*m_2601noCoeff+var]);
                }
                break;
              } else {
                effConst(pos,1) = m_2601effConst[iter][pos*m_2601noCoeff+1];
                effConst(pos,2) = m_2601effConst[iter][pos*m_2601noCoeff+2];
                effConst(pos,3) = m_2601effConst[iter][pos*m_2601noCoeff+3];
                effConst(pos,4) = m_2601effConst[iter][pos*m_2601noCoeff+4];
              }
            }
          }

          for(ZFSId pos=0; pos<m_2601noPos; pos++) {
            effConst(pos,0) = (effConst(pos,0) - rexPosQuad)/lengthFactor + 36.57;
          }

          //////////////////////////////////////////
          /////////// SPACE INTERPOLATION //////////
          //////////////////////////////////////////

          for(ZFSInt i = start[0]; i<end[0]; i++) {
            const ZFSId cellId=cellIndex(i,1,k);

            ZFSId cellIdG2=cellIndex(i,0,k); //ghost
            ZFSId cellIdG1=cellIndex(i,1,k); //ghost
            ZFSId cellIdA1=cellIndex(i,2,k); // field
            ZFSId cellIdA2=cellIndex(i,3,k); // field

            const ZFSId cellIdBcG1 = i + (1 + k*m_noGhostLayers)*m_block->m_nCells[2];
            const ZFSId cellIdBcG2 = i + (0 + k*m_noGhostLayers)*m_block->m_nCells[2];

            const ZFSFloat rho_zerothG1 = m_block->m_bc2601ZerothOrderSolution[PV->RHO][cellIdBcG1];
            const ZFSFloat u_zerothG1 = m_block->m_bc2601ZerothOrderSolution[PV->U][cellIdBcG1];
            const ZFSFloat v_zerothG1 = m_block->m_bc2601ZerothOrderSolution[PV->V][cellIdBcG1];
            const ZFSFloat w_zerothG1 = m_block->m_bc2601ZerothOrderSolution[PV->W][cellIdBcG1];

            //const ZFSFloat rho_zerothG2 = m_block->m_bc2601ZerothOrderSolution[PV->RHO][cellIdBcG2];
            const ZFSFloat u_zerothG2 = m_block->m_bc2601ZerothOrderSolution[PV->U][cellIdBcG2];
            const ZFSFloat v_zerothG2 = m_block->m_bc2601ZerothOrderSolution[PV->V][cellIdBcG2];
            const ZFSFloat w_zerothG2 = m_block->m_bc2601ZerothOrderSolution[PV->W][cellIdBcG2];

            const ZFSFloat uInfQuad = 1.029342800042823e+02;
            const ZFSFloat rhoInfQuad = 1.209997633989123;
            // const ZFSFloat heightG1 = m_cells->coordinates[1][cellIdG1];
            // const ZFSFloat heightG2 = m_cells->coordinates[1][cellIdG2];
            const ZFSFloat xPos = m_cells->coordinates[0][cellId];

            ZFSFloat uFactor = F0;
            ZFSFloat vFactor = F0;
            ZFSFloat wFactor = F0;
            ZFSFloat rhoFactor = F0;

            //interpolated from the actuated solution
            ZFSFloat uActuated = F0, vActuated = F0, wActuated = F0;

            //space interpolation
            if(xPos < effConst(0,0)) {
              const ZFSFloat relPos = (xPos-effConst(0,0))/(effConst(1,0)-effConst(0,0));
              uFactor = effConst(0,1) + relPos*(effConst(1,1)-effConst(0,1));
              vFactor = effConst(0,2) + relPos*(effConst(1,2)-effConst(0,2));
              wFactor = effConst(0,3) + relPos*(effConst(1,3)-effConst(0,3));
              rhoFactor = effConst(0,4) + relPos*(effConst(1,4)-effConst(0,4));
            } else if(effConst(0,0) <= xPos && xPos < effConst(m_2601noPos-1,0)) {
              for(ZFSId pos=0; pos<m_2601noPos-1; pos++) {
                if(effConst(pos,0) <= xPos && xPos < effConst(pos+1,0)) {
                  const ZFSFloat relPos = (xPos-effConst(pos,0))/(effConst(pos+1,0)-effConst(pos,0));
                  uFactor = effConst(pos,1) + relPos*(effConst(pos+1,1)-effConst(pos,1));
                  vFactor = effConst(pos,2) + relPos*(effConst(pos+1,2)-effConst(pos,2));
                  wFactor = effConst(pos,3) + relPos*(effConst(pos+1,3)-effConst(pos,3));
                  rhoFactor = effConst(pos,4) + relPos*(effConst(pos+1,4)-effConst(pos,4));
                  break;
                }
              }
            } else {
              const ZFSFloat relPos = (xPos-effConst(m_2601noPos-2,0))/(effConst(m_2601noPos-1,0)-effConst(m_2601noPos-2,0));
              uFactor   = effConst(m_2601noPos-2,1) + relPos*(effConst(m_2601noPos-1,1)-effConst(m_2601noPos-2,1));
              vFactor   = effConst(m_2601noPos-2,2) + relPos*(effConst(m_2601noPos-1,2)-effConst(m_2601noPos-2,2));
              wFactor   = effConst(m_2601noPos-2,3) + relPos*(effConst(m_2601noPos-1,3)-effConst(m_2601noPos-2,3));
              rhoFactor = effConst(m_2601noPos-2,4) + relPos*(effConst(m_2601noPos-1,4)-effConst(m_2601noPos-2,4));
            }

            // //for now take uniform distribution
            // uFactor   = effConst(0,1);
            // vFactor   = effConst(0,2);
            // rhoFactor = effConst(0,4);

            // ZFSFloat u_correctedG1 = u_zerothG1 + heightG1*lengthFactor*(uFactor)*(PV->UInfinity/uInfQuad);
            // ZFSFloat u_correctedG2 = u_zerothG2 + heightG2*lengthFactor*(uFactor)*(PV->UInfinity/uInfQuad);


            ZFSFloat u_zerothBC = u_zerothG1 + F1B2*(u_zerothG1-u_zerothG2);
            ZFSFloat u_correctedBC = u_zerothBC + gammaEpsilon*lengthFactor*(uFactor)*(PV->UInfinity/uInfQuad);
            ZFSFloat v_zerothBC = v_zerothG1 + F1B2*(v_zerothG1-v_zerothG2);
            ZFSFloat v_correctedBC = v_zerothBC + gammaEpsilon*lengthFactor*(vFactor)*(PV->UInfinity/uInfQuad);
            ZFSFloat w_zerothBC = w_zerothG1 + F1B2*(w_zerothG1-w_zerothG2);
            ZFSFloat w_correctedBC = w_zerothBC + gammaEpsilon*lengthFactor*(wFactor)*(PV->UInfinity/uInfQuad);// + wActuated;
            ZFSFloat rho_correctedG1 = rho_zerothG1 + gammaEpsilon*lengthFactor*(rhoFactor)*(CV->rhoInfinity/rhoInfQuad);

            if(m_block->m_nOffsetCells[2]+i-2==80&&m_block->m_nOffsetCells[0]+k-2==50) {

              if(globalTimeStep%20 == 0 && m_block->m_RKStep == 0) {

                cout.precision(7);
                cout << "globalTimeStep: " << globalTimeStep << " time: " << t << " gammaEpsilon: " << gammaEpsilon
                     << " uZeroth: " << u_zerothBC << " u_correctedBC: " << u_correctedBC
                     << " vZeroth: " << v_zerothBC << " v_correctedBC: " << v_correctedBC
                     << " rhoZeroth: " << rho_zerothG1 << " rho_correctedG1: " << rho_correctedG1
                     << " vActuated: " << vActuated << endl;


                FILE* f_effective;
                f_effective = fopen("./effective_boundary.dat", "a+");
                fprintf(f_effective, "%d", globalTimeStep);
                fprintf(f_effective, " %f", m_block->m_physicalTime);
                fprintf(f_effective, " %f", m_block->m_time);
                fprintf(f_effective, " %f", m_block->m_timeStep);
                fprintf(f_effective, " %f", u_zerothBC);
                fprintf(f_effective, " %f", u_correctedBC);
                fprintf(f_effective, " %f", v_zerothBC);
                fprintf(f_effective, " %f", v_correctedBC);
                fprintf(f_effective, " %f", w_zerothBC);
                fprintf(f_effective, " %f", w_correctedBC);
                fprintf(f_effective, " %f", rho_zerothG1);
                fprintf(f_effective, " %f", rho_correctedG1);
                fprintf(f_effective, " %f", uActuated);
                fprintf(f_effective, " %f", vActuated);
                fprintf(f_effective, " %f", wActuated);
                fprintf(f_effective, "\n");
                fclose(f_effective);
              }
            }

            const ZFSFloat u_uncorrectedBC = u_zerothG1 + F1B2*(u_zerothG1-u_zerothG2);
            const ZFSFloat v_uncorrectedBC = v_zerothG1 + F1B2*(v_zerothG1-v_zerothG2);
            const ZFSFloat w_uncorrectedBC = w_zerothG1 + F1B2*(w_zerothG1-w_zerothG2);
            const ZFSFloat  rho_uncorrectedG1 = rho_zerothG1;

            ZFSFloat u_appliedBC = F0, v_appliedBC = F0, w_appliedBC = F0, rho_appliedG1 = F0;

            if(m_cells->coordinates[0][cellId] < 12.4) {
              u_appliedBC = u_uncorrectedBC;
              v_appliedBC = v_uncorrectedBC;
              w_appliedBC = w_uncorrectedBC;
              rho_appliedG1 = rho_uncorrectedG1;
            } else if(m_cells->coordinates[0][cellId] >= 12.4 && m_cells->coordinates[0][cellId] <15.4) {
              const ZFSFloat fader = (m_cells->coordinates[0][cellId] - 12.4)/3.0;
              u_appliedBC = u_uncorrectedBC*(1.0-fader) + u_correctedBC*fader;
              v_appliedBC = v_uncorrectedBC*(1.0-fader) + v_correctedBC*fader;
              w_appliedBC = w_uncorrectedBC*(1.0-fader) + w_correctedBC*fader;
              rho_appliedG1 = rho_uncorrectedG1*(1.0-fader) + rho_correctedG1*fader;
            } else if(m_cells->coordinates[0][cellId] >= 15.4 && m_cells->coordinates[0][cellId] < 60.0) {
              u_appliedBC = u_correctedBC;
              v_appliedBC = v_correctedBC;
              w_appliedBC = w_correctedBC;
              rho_appliedG1 = rho_correctedG1;
            } else if(m_cells->coordinates[0][cellId] >= 60.0 && m_cells->coordinates[0][cellId] < 63.0) {
              const ZFSFloat fader = (m_cells->coordinates[0][cellId] - 60.0)/3.0;
              u_appliedBC = u_uncorrectedBC*(fader) + u_correctedBC*(1.0-fader);
              v_appliedBC = v_uncorrectedBC*(fader) + v_correctedBC*(1.0-fader);
              w_appliedBC = w_uncorrectedBC*(fader) + w_correctedBC*(1.0-fader);
              rho_appliedG1 = rho_uncorrectedG1*(fader) + rho_correctedG1*(1.0-fader);
            } else {
              u_appliedBC = u_uncorrectedBC;
              v_appliedBC = v_uncorrectedBC;
              w_appliedBC = w_uncorrectedBC;
              rho_appliedG1 = rho_uncorrectedG1;
            }

            // if(m_cells->coordinates[0][cellId] < 15.4) {
            //   u_appliedBC = u_uncorrectedBC;
            //   v_appliedBC = v_uncorrectedBC;
            //   w_appliedBC = w_uncorrectedBC;
            //   rho_appliedG1 = rho_uncorrectedG1;
            // } else if(m_cells->coordinates[0][cellId] >= 15.4 && m_cells->coordinates[0][cellId] < 60.0) {
            //   u_appliedBC = u_correctedBC;
            //   v_appliedBC = v_correctedBC;
            //   w_appliedBC = w_correctedBC;
            //   rho_appliedG1 = rho_correctedG1;
            // } else {
            //   u_appliedBC = u_uncorrectedBC;
            //   v_appliedBC = v_uncorrectedBC;
            //   w_appliedBC = w_uncorrectedBC;
            //   rho_appliedG1 = rho_uncorrectedG1;
            // }

            // const ZFSFloat rho_correctedG2 = rho_zerothG2 + heightG2*lengthFactor*(rhoFactor)*(CV->rhoInfinity/rhoInfQuad);
            // const ZFSFloat rho_correctedBC = rho_correctedG1 + F1B2*(rho_correctedG1-rho_correctedG2);
            const ZFSFloat pField  = pressure(cellIndex(i,m_noGhostLayers,k));


            const ZFSFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
            const ZFSFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
            const ZFSFloat w1 = m_cells->pvariables[PV->W][cellIdA1];
            const ZFSFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];

            const ZFSFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
            const ZFSFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
            const ZFSFloat w2 = m_cells->pvariables[PV->W][cellIdA2];
            // ZFSFloat rho2 = m_cells->pvariables[PV->RHO][cellIdA2];

            const ZFSFloat uG1 = F2*u_appliedBC-u1;
            const ZFSFloat vG1 = F2*v_appliedBC-v1;
            const ZFSFloat wG1 = F2*w_appliedBC-w1;

            const ZFSFloat uG2 = F2*u_appliedBC-u2;
            const ZFSFloat vG2 = F2*v_appliedBC-v2;
            const ZFSFloat wG2 = F2*w_appliedBC-w2;

            const ZFSFloat rhoG1 = rho_appliedG1;//F2*rho_correctedBC-rho1;
            const ZFSFloat rhoG2 = F2*rho_appliedG1-rho1;

            m_cells->pvariables[PV->RHO][cellIdG1] = rhoG1;
            m_cells->pvariables[PV->U][cellIdG1] = uG1;
            m_cells->pvariables[PV->V][cellIdG1] = vG1;
            m_cells->pvariables[PV->W][cellIdG1] = wG1;
            m_cells->pvariables[PV->P][cellIdG1] = pField;

            m_cells->pvariables[PV->RHO][cellIdG2] = rhoG2;
            m_cells->pvariables[PV->U][cellIdG2] = uG2;
            m_cells->pvariables[PV->V][cellIdG2] = vG2;
            m_cells->pvariables[PV->W][cellIdG2] = wG2;
            m_cells->pvariables[PV->P][cellIdG2] = pField;

            if(m_block->m_rans) {
              const ZFSFloat nutilde_zerothG1 = m_block->m_bc2601ZerothOrderSolution[PV->RANS_VAR[0]][cellIdBcG1];
              const ZFSFloat nutilde_zerothG2 = m_block->m_bc2601ZerothOrderSolution[PV->RANS_VAR[0]][cellIdBcG2];
              m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] = nutilde_zerothG1;
              m_cells->pvariables[PV->RANS_VAR[0]][cellIdG2] = nutilde_zerothG2;
            }
          }
        }
        break;
      }
    default:
      {
        zfsTerm(1, __CALLING_FUNCTION__, "Face not implemented");
        break;
      }
    }
}

/** \brief supersonic inflow with imposed acoustic or entropy waves
 * version: cut-off boundary condition
 * author: Thomas Schilden, 12.2.2015
 * edited by: Leo Hoening
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2700(ZFSId bcId)
{

  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  const ZFSFloat gamma = m_block->m_gamma;
  const ZFSFloat time = m_block->m_time;
  ZFSFloat velocity[3] = {F0,F0,F0};

  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
    case 3:
      {
        for(ZFSId k=start[2]; k<end[2] ; k++) {
          for(ZFSId j=start[1]; j<end[1] ; j++) {
            for(ZFSId i=start[0]; i<end[0]; i++) {
              ZFSId cellId=cellIndex(i,j,k);
              ZFSFloat rho = CV->rhoInfinity;
              for(ZFSInt dim = 0; dim < nDim; dim++){
                velocity[dim] =  PV->VVInfinity[dim];
              }
              ZFSFloat p   = PV->PInfinity;

              // add the modes
              for(ZFSInt mode = 0; mode < m_modes; mode++) {
                //1. pressure
                const ZFSFloat pressure_f = gamma * m_Ma * PV->PInfinity * m_modeAmp[mode] * (ZFSFloat)(m_modeType[mode]);

                //2. density
                ZFSFloat density_f;
                const ZFSFloat a = sqrt( PV->TInfinity );
                if(m_modeType[mode])
                  density_f = pressure_f / POW2(a);
                else
                  density_f = m_modeAmp[mode] * CV->rhoInfinity * m_Ma;

                //3. velocity
                ZFSFloat K = F0;
                for(ZFSInt dim = 0; dim < nDim; dim++){
                  K += pow(m_modeK[mode][dim],2);
                }

                K = sqrt(K);
                const ZFSFloat acImp = a * CV->rhoInfinity;
                const ZFSFloat modeVelocity = (ZFSFloat)(m_modeType[mode]) * pressure_f / acImp;
                ZFSFloat velocity_f[3];
                for(ZFSInt dim = 0; dim < nDim; dim++){
                  velocity_f[dim] =  ( m_modeK[mode][dim] * modeVelocity ) / K;
                }

                //1. calculate cell dependant trigonometry
                ZFSFloat trigTerm = F0;
                for(ZFSInt dim = 0; dim < nDim; dim++) {
                  trigTerm += m_modeK[mode][dim] * m_cells->coordinates[dim][cellId];
                }

                trigTerm -= m_modeOmega[mode] * (time);//time;
                trigTerm += m_modePhi[mode];
                trigTerm = sin(trigTerm);

                rho += trigTerm * density_f;
                for(ZFSInt dim = 0; dim < nDim; dim++) {
                  velocity[dim] += trigTerm * velocity_f[0];
                }
                p += trigTerm * pressure_f;
              }

              // conservatives:
              m_cells->pvariables[PV->RHO][cellId]  =rho;
              m_cells->pvariables[PV->U][cellId]=velocity[PV->U];
              m_cells->pvariables[PV->V][cellId]=velocity[PV->V];
              m_cells->pvariables[PV->W][cellId]=velocity[PV->W];
              m_cells->pvariables[PV->P][cellId]=p;
            }
          }
        }
        break;
      }
    default:
      {
        zfsTerm(1, __CALLING_FUNCTION__, "Face direction not implemented)");
      }
    }
}


//Jet Freund Inlet
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2900(ZFSId bcId){
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
      {
        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;
        ZFSId pIJK=0, pIJPK=0, pIJKP=0, pIJPKP=0;
        //fully new bc
        ZFSId i= end[0]-1;
        ZFSFloat pBC=F0, rho=F0, u=F0, v=F0, w=F0;
        ZFSFloat drho=F0, du=F0, dv=F0, dw=F0, dp=F0;
        ZFSFloat xBC=F0;
        //pBC=PV->PInfinity;
        ZFSFloat pInner=F0, c02=F0, distance=F0;
        for(ZFSId k=start[2]; k<end[2] ; k++)
          {
            for(ZFSId j=start[1]; j<end[1] ; j++)
              {
                cellId=cellIndex(i,j,k);
                cellIdadj=cellIndex(i+1,j,k);
                //to determine the face coordinates!!!!!!
                pIJK=getPointIdFromCell(i+1,j,k);
                pIJPK=getPointIdfromPoint(pIJK,0,1,0);
                pIJKP=getPointIdfromPoint(pIJK,0,0,1);
                pIJPKP=getPointIdfromPoint(pIJK,0,1,1);
                ZFSFloat r= sqrt(POW2(m_cells->coordinates[0][cellId])+POW2(m_cells->coordinates[1][cellId])+POW2(m_cells->coordinates[2][cellId]));
                ZFSFloat uInf = F1B2*(F1-tanh(12.5*(fabs(r/0.5)-fabs(0.5/r))))*PV->VVInfinity[0];
                ZFSFloat vInf =0;
                ZFSFloat wInf =0;

                //values at the inner point
                pInner=m_cells->pvariables[PV->P][cellIdadj];
                c02=sqrt(m_block->m_gamma*pInner/m_cells->pvariables[PV->RHO][cellIdadj]);
                u=m_cells->pvariables[PV->U][cellIdadj];
                v=m_cells->pvariables[PV->V][cellIdadj];
                w=m_cells->pvariables[PV->W][cellIdadj];

                ZFSFloat dxidx=m_cells->surfaceMetrics[cellId][0];
                ZFSFloat dxidy=m_cells->surfaceMetrics[cellId][1];
                ZFSFloat dxidz=m_cells->surfaceMetrics[cellId][2];

                //values at the boundary
                pBC=F1B2*(pInner+PV->PInfinity-m_cells->pvariables[PV->RHO][cellIdadj]*c02*(dxidx*(uInf-u)+dxidy*(vInf-v)+dxidz*(wInf-w)));
                rho=CV->rhoInfinity+((pBC-PV->PInfinity)/(c02*c02));

                u= uInf-dxidx*(PV->PInfinity-pBC)/(m_cells->pvariables[PV->RHO][cellIdadj]*c02);
                v= vInf-dxidy*(PV->PInfinity-pBC)/(m_cells->pvariables[PV->RHO][cellIdadj]*c02);
                w= wInf-dxidz*(PV->PInfinity-pBC)/(m_cells->pvariables[PV->RHO][cellIdadj]*c02);

                //extrapolate the variables into the ghost cells
                //gradients

                xBC=F1B4*(m_coordinates[0][pIJK]+m_coordinates[0][pIJPK]+m_coordinates[0][pIJKP]+m_coordinates[0][pIJPKP]);


                distance=(xBC-m_cells->coordinates[0][cellIdadj]);

                drho=(rho-m_cells->pvariables[PV->RHO][cellIdadj])/distance;
                du=(u-m_cells->pvariables[PV->U][cellIdadj])/distance;
                dv=(v-m_cells->pvariables[PV->V][cellIdadj])/distance;
                dw=(w-m_cells->pvariables[PV->W][cellIdadj])/distance;
                dp=(pBC-m_cells->pvariables[PV->P][cellIdadj])/distance;

                //extrapolate:
                for(ZFSId ii=start[0]; ii<end[0]; ++ii)
                  {
                    cellId=cellIndex(ii,j,k);
                    distance=(m_cells->coordinates[0][cellId]-m_cells->coordinates[0][cellIdadj]);
                    m_cells->pvariables[PV->RHO][cellId]=m_cells->pvariables[PV->RHO][cellIdadj]+drho*distance;
                    m_cells->pvariables[PV->U][cellId]=m_cells->pvariables[PV->U][cellIdadj]+du*distance;
                    m_cells->pvariables[PV->V][cellId]=m_cells->pvariables[PV->V][cellIdadj]+dv*distance;
                    m_cells->pvariables[PV->W][cellId]=m_cells->pvariables[PV->W][cellIdadj]+dw*distance;
                    m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj]+dp*distance;
                  }

              }
          }
        break;
      }
    default:{
      zfsTerm(1, __CALLING_FUNCTION__, "Face not implemented");
      break;
    }
    }


}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc6000(ZFSId bcId)
{
  cout << "applying bc " << bcId << endl;
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::bc2402(ZFSId bcId)
{
  TRACE();

  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  //first we need the average pressure and Temperature
  //==> calculate from the field
  //==> average over the surface
  //==> store in variable

  ZFSId Id=m_channelSurfacIndexMap[bcId];
  ZFSId* startface = &m_block->m_windowInfo->channelSurfaceIndices[Id]->start1[0];
  ZFSId* endface = &m_block->m_windowInfo->channelSurfaceIndices[Id]->end1[0];
  ZFSId cellId=0;
  //global Pressure and Temperature at in-/outlet
  ZFSFloat globalTin[2]={F0,F0};
  ZFSFloat globalTout[2]={F0,F0};
  ZFSFloat globalPin[2]={F0,F0};
  ZFSFloat globalPout[2]={F0,F0};
  cout.precision(10);

  //find out which face it is
  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
      {
        //average the pressure
        //use values from the inside to determine the pressure at the face!!!
        ZFSFloat surface=F0;
        ZFSFloat localPin[2]={F0,F0};
        ZFSFloat localTin[2]={F0,F0};
        ZFSFloat localVel = F0, globalVel = F0;
        ZFSFloat localMassFlux = F0, globalMassFlux = F0;

        for(ZFSId k=startface[2]; k<endface[2] ; k++) {
            for(ZFSId j=startface[1]; j<endface[1] ; j++) {
              ZFSId ii=end[0]-1;
              ZFSId cellIdP1=cellIndex(ii+1,j,k);
              ZFSId cellIdP2=cellIndex(ii+2,j,k);

              surface=sqrt(POW2(m_cells->surfaceMetrics[cellIdP1][0])+POW2(m_cells->surfaceMetrics[cellIdP1][1])+POW2(m_cells->surfaceMetrics[cellIdP1][2]));
              localPin[0]+=surface*m_cells->pvariables[PV->P][cellIdP1];
              localPin[1]+=surface*m_cells->pvariables[PV->P][cellIdP2];
              localTin[0]+=surface*temperature(cellIdP1);
              localTin[1]+=surface*temperature(cellIdP2);
              localVel += surface*(m_cells->pvariables[PV->U][cellIdP1]);
              localMassFlux += surface*(m_cells->pvariables[PV->U][cellIdP1]*m_cells->pvariables[PV->RHO][cellIdP1]);
            }
        }

        //next now that the pressure and temperature are known:
        //make it a global variable for the complete inflow plane
        MPI_Allreduce(localPin, globalPin, 2, MPI_DOUBLE, MPI_SUM, m_block->m_commChannelIn[0]);
        MPI_Allreduce(localTin, globalTin, 2, MPI_DOUBLE, MPI_SUM, m_block->m_commChannelIn[0]);
        MPI_Allreduce(&localMassFlux, &globalMassFlux, 1, MPI_DOUBLE, MPI_SUM, m_block->m_commChannelIn[0]);
        MPI_Allreduce(&localVel, &globalVel, 1, MPI_DOUBLE, MPI_SUM, m_block->m_commChannelIn[0]);

        MPI_Barrier(m_block->m_commChannelIn[0]);
        globalPin[0]/=m_channelSurfaceIn;
        globalTin[0]/=m_channelSurfaceIn;
        globalPin[1]/=m_channelSurfaceIn;
        globalTin[1]/=m_channelSurfaceIn;
        globalMassFlux/=m_channelSurfaceIn;
        globalVel/=m_channelSurfaceIn;
        ZFSFloat currentRe = m_block->m_Re/PV->UInfinity*globalVel;

        //set a Barrier to ensure that all cpus having a channel side are at same time level
        //now make the variables known everywhere in the channel world
        MPI_Barrier(m_block->m_commChannelWorld[0]);
        MPI_Bcast(globalTin, 2, MPI_DOUBLE, m_block->m_channelRoots[2], m_block->m_commChannelWorld[0]);
        MPI_Bcast(globalPin, 2, MPI_DOUBLE, m_block->m_channelRoots[2], m_block->m_commChannelWorld[0]);
        MPI_Bcast(globalTout, 2, MPI_DOUBLE, m_block->m_channelRoots[3], m_block->m_commChannelWorld[0]);
        MPI_Bcast(globalPout, 2, MPI_DOUBLE, m_block->m_channelRoots[3], m_block->m_commChannelWorld[0]);
        //now every value is known and can be used to apply the BC !!!!
        //cannot take values from above at here the window is extended to the ghost layers

        if(globalTimeStep > 1 && globalTimeStep%50==0) {
          if(m_channelInflowRank == 0 && globalTimeStep%m_block->m_residualInterval == 0 && m_block->m_RKStep == 0) {
            cout.precision(6);

            FILE* f_channel;
            f_channel = fopen("./massflow.dat", "a+");
            fprintf(f_channel, "%d", globalTimeStep);
            fprintf(f_channel, " %f", m_block->m_physicalTime);
            fprintf(f_channel, " %f", m_block->m_time);
            fprintf(f_channel, " %f", m_block->m_timeStep);
            fprintf(f_channel, " %f", currentRe);
            fprintf(f_channel, " %f", globalMassFlux);
            fprintf(f_channel, " %f", globalPin[0]);
            fprintf(f_channel, " %f", globalPin[1]);
            fprintf(f_channel, " %f", globalTin[0]);
            fprintf(f_channel, " %f", globalTout[0]);
            fprintf(f_channel, " %f", globalPout[0]);
            fprintf(f_channel, " %f", globalPout[1]);
            fprintf(f_channel, " %f", (globalTin[0] - globalTout[1]));
            fprintf(f_channel, "\n");
            fclose(f_channel);
          }
        }

        for(ZFSId k=start[2]; k<end[2] ; k++){
          for(ZFSId j=start[1]; j<end[1] ; j++){
            ZFSId i = 1;
            cellId=cellIndex(i,j,k);

            ZFSFloat pressureFluctuation = m_cells->pvariables[PV->P][cellId] - globalPout[1];
            ZFSFloat x = m_cells->coordinates[0][cellId];
            ZFSFloat pressureInflowMean = m_block->m_deltaP/m_block->m_channelLength*(x - m_block->m_channelInflowPlaneCoordinate) + PV->PInfinity;
            ZFSFloat pressureNew = pressureInflowMean + pressureFluctuation;
            ZFSFloat temperatureFlucOutflow = temperature(cellId) - globalTout[1];
            ZFSFloat temperatureNew = PV->TInfinity + temperatureFlucOutflow;
            ZFSFloat rhoNew = m_block->m_gamma*pressureNew/temperatureNew;

            m_cells->pvariables[PV->RHO][cellId] = rhoNew;
            m_cells->pvariables[PV->P][cellId]= pressureNew;

            for(ZFSId var=0; var<PV->noVariables; var++) {
              m_cells->pvariables[var][cellIndex(start[0],j,k)] = 2.0*m_cells->pvariables[var][cellIndex(start[0]+1,j,k)] - m_cells->pvariables[var][cellIndex(start[0]+2,j,k)];
            }
          }
        }
        break;
      }
    case 1:
      {
        //average the pressure
        //use values from the inside to determine the pressure at the face!!!
        ZFSFloat surface=F0;
        ZFSFloat localPout[]={F0,F0};
        ZFSFloat localTout[]={F0,F0};
        for(ZFSId k=startface[2]; k<endface[2] ; k++){
          for(ZFSId j=startface[1]; j<endface[1] ; j++){
            ZFSId ii=start[0];
            ZFSId cellIdM1=cellIndex(ii-1,j,k);
            ZFSId cellIdM2=cellIndex(ii-2,j,k);
            surface=sqrt(POW2(m_cells->surfaceMetrics[cellIdM1][0])+POW2(m_cells->surfaceMetrics[cellIdM1][1])+POW2(m_cells->surfaceMetrics[cellIdM1][2]));
            localPout[0]+=surface*m_cells->pvariables[PV->P][cellIdM2];
            localPout[1]+=surface*m_cells->pvariables[PV->P][cellIdM1];
            localTout[0]+=surface*temperature(cellIdM2);
            localTout[1]+=surface*temperature(cellIdM1);
          }
        }

        //next now that the pressure and temperature are known:
        //make it a global variable for the complete outflow plane
        MPI_Allreduce(localPout, globalPout, 2, MPI_DOUBLE, MPI_SUM, m_block->m_commChannelOut[0]);
        MPI_Allreduce(localTout, globalTout, 2, MPI_DOUBLE, MPI_SUM, m_block->m_commChannelOut[0]);

        MPI_Barrier(m_block->m_commChannelOut[0]);
        globalPout[0]/=m_channelSurfaceOut;
        globalTout[0]/=m_channelSurfaceOut;
        globalPout[1]/=m_channelSurfaceOut;
        globalTout[1]/=m_channelSurfaceOut;
        //set a Barrier to ensure that all cpus having a channel side are at same time level
        //now make the variables known everywhere in the channel world
        MPI_Barrier(m_block->m_commChannelWorld[0]);
        MPI_Bcast(globalTin, 2, MPI_DOUBLE, m_block->m_channelRoots[2], m_block->m_commChannelWorld[0]);
        MPI_Bcast(globalPin, 2, MPI_DOUBLE, m_block->m_channelRoots[2], m_block->m_commChannelWorld[0]);
        MPI_Bcast(globalTout, 2, MPI_DOUBLE, m_block->m_channelRoots[3], m_block->m_commChannelWorld[0]);
        MPI_Bcast(globalPout, 2, MPI_DOUBLE, m_block->m_channelRoots[3], m_block->m_commChannelWorld[0]);
        //now every value is known and can be used to apply the BC !!!!
        //cannot take values from above at here the window is extended to the ghost layers

        for(ZFSId k=start[2]; k<end[2] ; k++){
          for(ZFSId j=start[1]; j<end[1] ; j++){
            ZFSId i = start[0];
            cellId=cellIndex(i,j,k);

            ZFSFloat pressureFluctuation = m_cells->pvariables[PV->P][cellId] - globalPin[0];
            ZFSFloat x = m_cells->coordinates[0][cellId];
            ZFSFloat pressureOutflowMean = m_block->m_deltaP/m_block->m_channelLength*(x - m_block->m_channelInflowPlaneCoordinate) + PV->PInfinity;
            ZFSFloat pressureNew = pressureOutflowMean + pressureFluctuation;
            ZFSFloat deltaT = (globalTin[0] - globalTout[1]);
            ZFSFloat temperatureInflow = temperature(cellId);
            ZFSFloat temperatureNew = temperatureInflow - deltaT;
            ZFSFloat rhoNew = m_block->m_gamma*pressureNew/temperatureNew;

            m_cells->pvariables[PV->RHO][cellId] = rhoNew;
            m_cells->pvariables[PV->P][cellId]=pressureNew;

            for(ZFSId var=0; var<PV->noVariables; var++) {
              m_cells->pvariables[var][cellIndex(start[0]+1,j,k)] = 2.0*m_cells->pvariables[var][cellIndex(start[0],j,k)] - m_cells->pvariables[var][cellIndex(start[0]-1,j,k)];
            }
          }
        }
        break;
      }
    default:{
      zfsTerm(1, __CALLING_FUNCTION__, "Face not implemented");
      break;
    }
    }
}


template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::periodicPointsChange(ZFSFloat* pt,ZFSId type)
{
  ZFSFloat angle=0;
  /*
    case 4401 4402: first periodic direction
    case 4403 4404: second periodic direction
    case 4405 4405: third periodic direction
    case 4011 4012:  rotation X axis clockwise and anticlockwise
   */
  ZFSFloatScratchSpace rotationMatrix(3, 3, AT_, "rotation Matrix");
  rotationMatrix.fill(F0);
  ZFSFloat tmp[3];

  // SND Map BC is reversed, therefore we need
  // to subtract instead of adding and vice versa compared
  // to the periodicPointsChange in zfsstrctrdblckwindowinfo.cpp
  switch(type)
    {
    case 4401:
    case 4403:
    case 4405:
      {
        const ZFSId displacementId = (ZFSFloat)(type-4400+1)/2.0 - 1;
        for (ZFSInt dim=0;dim<nDim;++dim)
          pt[dim]=pt[dim]-m_block->m_periodicDisplacements[dim*nDim+displacementId];
        break;
      }

    case 4402:
    case 4404:
    case 4406:
      {
        const ZFSId displacementId = (ZFSFloat)(type-4400+1)/2.0 - 1;
        for (ZFSInt dim=0;dim<nDim;++dim)
          pt[dim]=pt[dim]+m_block->m_periodicDisplacements[dim*nDim+displacementId];
        break;
      }

    case 4012:
      {
        rotationMatrix(1,1)=cos(angle);
        rotationMatrix(1,2)=-sin(angle);
        rotationMatrix(2,1)=sin(angle);
        rotationMatrix(2,2)=cos(angle);
        tmp[0]=pt[0];
        tmp[1]=pt[1];
        tmp[2]=pt[2];
        for (ZFSInt i=0;i<3;++i)
          {
      pt[i]=rotationMatrix(i,0)*tmp[0]+rotationMatrix(i,1)*tmp[1]
                  +rotationMatrix(i,2)*tmp[2];
          }
        break;
      }

    case 4011:
      {
        rotationMatrix(1,1)=cos(-angle);
        rotationMatrix(1,2)=-sin(-angle);
        rotationMatrix(2,1)=sin(-angle);
        rotationMatrix(2,2)=cos(-angle);
        tmp[0]=pt[0];
        tmp[1]=pt[1];
        tmp[2]=pt[2];
        for (ZFSInt i=0;i<3;++i)
          {
      pt[i]=rotationMatrix(i,0)*tmp[0]+rotationMatrix(i,1)*tmp[1]
                  +rotationMatrix(i,2)*tmp[2];
          }
        break;
      }

    default:
      {
        cout<<"ERROR!!! periodic type is wrong!!! in BC call BC: "<<type<<endl;
      }
    }


}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::exchangePointsPeriodic()
{
  gatherPoints();

  sendPoints();

  receivePoints();

  for(ZFSId i=m_startCommPeriodic; i<m_endCommPeriodic-m_periodicS; i++)
    {
      MPI_Wait(&(m_block->m_cmnctnFlag->mpi_sndRequest[i]), &(m_block->m_cmnctnFlag->mpi_sndStatus[i]));
    }
  for(ZFSId i=m_startCommPeriodic; i<m_endCommPeriodic-m_periodicS; i++)
    {
      MPI_Wait(&(m_block->m_cmnctnFlag->mpi_rcvRequest[i]), &(m_block->m_cmnctnFlag->mpi_rcvStatus[i]));
    }

  scatterPoints();
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::gatherPoints()
{
  ZFSId pointId;
  ZFSFloat tmppoints[3];

  for(ZFSInt nghbr=m_startCommPeriodic; nghbr < m_endCommPeriodic-m_periodicS; nghbr++) {
    ZFSInt* startInfo=m_block->m_cmnctnFlag->startInfoSNDpoints[nghbr];
    ZFSInt* endInfo= m_block->m_cmnctnFlag->endInfoSNDpoints[nghbr];
    ZFSFloat* bufferSnd = m_block->m_cmnctnFlag->m_bufferPointsSnd[nghbr];
    ZFSInt pos=0;
    ZFSInt BC=m_block->m_cmnctnFlag->bcId[nghbr];

    for(ZFSId dim=0; dim<nDim; dim++) {
      for(ZFSInt k=startInfo[2]; k<endInfo[2]+1; k++) {
        for(ZFSInt j=startInfo[1]; j<endInfo[1]+1; j++) {
          for(ZFSInt i=startInfo[0]; i<endInfo[0]+1; i++) {
            pointId = i +(j+k*m_block->m_nPoints[1])*m_block->m_nPoints[2];
            tmppoints[0]=m_coordinates[0][pointId];
            tmppoints[1]=m_coordinates[1][pointId];
            tmppoints[2]=m_coordinates[2][pointId];
            periodicPointsChange(tmppoints,BC);
            bufferSnd[pos]=tmppoints[dim];
            pos++;
          }
        }
      }
    }
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::sendPoints()
{
  for(ZFSInt nghbr=m_startCommPeriodic; nghbr < m_endCommPeriodic-m_periodicS; nghbr++) {
    ZFSInt tag= m_block->domainId()+m_block->m_cmnctnFlag->m_tagHelperSND[nghbr]*m_block->noDomains();
    //if(m_block->m_cmnctnFlag->m_sndNghbrId[nghbr]!=m_block->domainId()) {

      ZFSInt err= MPI_Isend((void*)&m_block->m_cmnctnFlag->m_bufferPointsSnd[nghbr][0], m_block->m_cmnctnFlag->m_noNghbrDomainPointBufferSizeSnd[nghbr], MPI_DOUBLE, m_block->m_cmnctnFlag->m_sndNghbrId[nghbr], tag, MPI_COMM_WORLD, &m_block->m_cmnctnFlag->mpi_sndRequest[nghbr]);
      if(err) cout << "rank " << m_block->domainId() << " sending periodic points throws error " << endl;
      //}
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::receivePoints()
{
  for(ZFSInt nghbr=m_startCommPeriodic; nghbr < m_endCommPeriodic-m_periodicS; nghbr++) {
    ZFSInt tag=m_block->m_cmnctnFlag->m_rcvNghbrId[nghbr]+(m_block->m_cmnctnFlag->m_tagHelperRCV[nghbr])*m_block->noDomains();
    // if(m_block->m_cmnctnFlag->m_sndNghbrId[nghbr]!=m_block->domainId()) {
      ZFSInt err=MPI_Irecv((void*)&m_block->m_cmnctnFlag->m_bufferPointsRcv[nghbr][0],m_block->m_cmnctnFlag->m_noNghbrDomainPointBufferSizeRcv[nghbr], MPI_DOUBLE,m_block->m_cmnctnFlag->m_rcvNghbrId[nghbr] ,tag, MPI_COMM_WORLD, &m_block->m_cmnctnFlag->mpi_rcvRequest[nghbr]);
      if(err) cout << "rank " << m_block->domainId() << " recv periodic points throws error " << endl;
    // } else {
    //   ZFSInt i,tag1;
    //   for(i = m_startCommPeriodic; i < m_endCommPeriodic-m_periodicS; i++) {
    //          tag1= m_block->domainId()+m_block->m_cmnctnFlag->m_tagHelperSND[i]*m_block->noDomains();
    //          if( tag1==tag&&m_block->domainId()== m_block->m_cmnctnFlag->m_sndNghbrId[i]) {
    //            break;
    //          }
    //   }
    //   if(m_block->m_cmnctnFlag->m_noNghbrDomainPointBufferSizeRcv[nghbr]!=m_block->m_cmnctnFlag->m_noNghbrDomainPointBufferSizeSnd[i]) {
    //          cout<<"TAG "<<tag<<" TAG1 "<<tag1<<" Id "<<nghbr<<" "<<i<<endl;
    //          cout<<"ERROR!!!!!!!!!!!!!SndpointBuffer "<<m_block->m_cmnctnFlag->m_noNghbrDomainPointBufferSizeSnd[i]<<" does not match RecvpointBuffer "<<m_block->m_cmnctnFlag->m_noNghbrDomainPointBufferSizeRcv[nghbr]<<endl;
    //          zfsTerm(1, __CALLING_FUNCTION__, "receive Points in periodicPointsExchange");
    //   } else {
    //          memcpy(m_block->m_cmnctnFlag->m_bufferPointsRcv[nghbr], m_block->m_cmnctnFlag->m_bufferPointsSnd[i], m_block->m_cmnctnFlag->m_noNghbrDomainPointBufferSizeRcv[nghbr]*sizeof(ZFSFloat));
    //   }
    // }
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::scatterPoints()
{
  ZFSId pointId;
  //  ZFSInt cellId;
  //the ordering of the grid points can be different from
  //sending instance ==> reorder it and copy it to the
  //right place
  for(ZFSInt nghbr=m_startCommPeriodic; nghbr < m_endCommPeriodic-m_periodicS; nghbr++) {
    ZFSInt k2, j2, i2, id2;
    ZFSInt* step1 = m_block->m_cmnctnFlag->stepInfoRCV[nghbr];
    ZFSInt step2[3];

    ZFSInt* order = m_block->m_cmnctnFlag->orderInfo[nghbr];
    ZFSInt start1[3];
    //    ZFSInt end1[3];
    ZFSInt start2[3];
    ZFSInt end2[3];
    ZFSInt len2[3];
    ZFSInt totalPoints=1;
    ZFSInt len1[3];

    for(ZFSInt j=0; j<nDim; j++) {
      len1[j]=m_block->m_cmnctnFlag->endInfoRCVpoints[nghbr][j] - m_block->m_cmnctnFlag->startInfoRCVpoints[nghbr][j]+1;
      totalPoints*=len1[j];
      step2[order[j]]=step1[j];
    }

    for(ZFSInt j=0; j<nDim; j++) {
      start2[j]=0;
      end2[j]=len1[j]-1;
      len2[order[j]]=len1[j];
      if(step2[j]<0) {
        ZFSInt dummy=start2[j];
        start2[j]=end2[j];
        end2[j]=dummy;
      }
    }

    ZFSInt* startInfo=m_block->m_cmnctnFlag->startInfoRCVpoints[nghbr];
    ZFSInt* endInfo= m_block->m_cmnctnFlag->endInfoRCVpoints[nghbr];

    ZFSFloat* bufferRcv = m_block->m_cmnctnFlag->m_bufferPointsRcv[nghbr];

    for(ZFSId dim=0; dim<nDim; dim++) {
      k2=start2[2];
      for(ZFSInt k=startInfo[2]; k<endInfo[2]+1; k++) {
        j2=start2[1];
        for(ZFSInt j=startInfo[1]; j<endInfo[1]+1; j++) {
          i2=start2[0];
          for(ZFSInt i=startInfo[0]; i<endInfo[0]+1; i++) {
            start1[order[0]]=i2;
            start1[order[1]]=j2;
            start1[order[2]]=k2;
            id2=dim*totalPoints+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
            pointId = i +(j+k*m_nPoints[1])*m_nPoints[2];
            m_coordinates[dim][pointId]= bufferRcv[id2];
            i2+=step2[0];
          }
          j2+=step2[1];
        }
        k2+=step2[2];
      }
    }
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::exchangePointsPeriodicS()
{
  receivePointsS();

  gatherPointsS();

  sendPointsS();

  for(ZFSId i=m_endCommPeriodic-m_periodicS; i< m_endCommPeriodic; i++)
    {
      MPI_Wait(&(m_block->m_cmnctnFlag->mpi_sndRequest[i]), &(m_block->m_cmnctnFlag->mpi_sndStatus[i]));
    }
  for(ZFSId i=m_endCommPeriodic-m_periodicS; i< m_endCommPeriodic; i++)
    {
      MPI_Wait(&(m_block->m_cmnctnFlag->mpi_rcvRequest[i]), &(m_block->m_cmnctnFlag->mpi_rcvStatus[i]));
    }

  scatterPointsS();
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::gatherPointsS()
{
  ZFSId cellId;
  ZFSFloat tmppoints[3];

  for(ZFSInt nghbr=m_endCommPeriodic-m_periodicS; nghbr<m_endCommPeriodic; nghbr++) {
    ZFSInt* startInfo=m_block->m_cmnctnFlag->startInfoSNDpoints[nghbr];
    ZFSInt* endInfo= m_block->m_cmnctnFlag->endInfoSNDpoints[nghbr];
    ZFSFloat* bufferSnd = m_block->m_cmnctnFlag->m_bufferPointsSnd[nghbr];
    ZFSInt BC=m_block->m_cmnctnFlag->bcId[nghbr];
    ZFSInt pos=0;

    for(ZFSId dim=0; dim<nDim; dim++) {
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            tmppoints[0]=m_cells->coordinates[dim][cellId];
            tmppoints[1]=m_cells->coordinates[dim][cellId];
            tmppoints[2]=m_cells->coordinates[dim][cellId];
            periodicPointsChange(tmppoints,BC);
            bufferSnd[pos]=tmppoints[dim];
            pos++;
          }
        }
      }
    }
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::sendPointsS()
{
  for(ZFSInt nghbr= m_endCommPeriodic-m_periodicS; nghbr< m_endCommPeriodic; nghbr++) {
    ZFSInt tag= m_block->domainId()+(m_block->m_cmnctnFlag->m_tagHelperSND[nghbr])*m_block->noDomains();
    ZFSInt err= MPI_Isend((void*)&m_block->m_cmnctnFlag->m_bufferPointsSnd[nghbr][0], m_block->m_cmnctnFlag->m_noNghbrDomainPointBufferSizeSnd[nghbr], MPI_DOUBLE, m_block->m_cmnctnFlag->m_sndNghbrId[nghbr], tag, MPI_COMM_WORLD, &m_block->m_cmnctnFlag->mpi_sndRequest[nghbr]);
    if(err) {
      cout << "rank " << m_block->domainId() << " sending throws error " << endl;
    }
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::receivePointsS()
{
  for(ZFSInt nghbr= m_endCommPeriodic-m_periodicS; nghbr< m_endCommPeriodic; nghbr++) {
    ZFSInt tag=m_block->m_cmnctnFlag->m_rcvNghbrId[nghbr]+(m_block->m_cmnctnFlag->m_tagHelperRCV[nghbr])*m_block->noDomains();

    ZFSInt err=MPI_Irecv((void*)&m_block->m_cmnctnFlag->m_bufferPointsRcv[nghbr][0],m_block->m_cmnctnFlag->m_noNghbrDomainPointBufferSizeRcv[nghbr], MPI_DOUBLE,m_block->m_cmnctnFlag->m_rcvNghbrId[nghbr] ,tag, MPI_COMM_WORLD, &m_block->m_cmnctnFlag->mpi_rcvRequest[nghbr]);
    if(err) cout << "rank " << m_block->domainId() << " recv periodic points throws error " << endl;
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::scatterPointsS()
{
  ZFSInt cellId;
  //the ordering of the grid points can be different from
  //sending instance ==> reorder it and copy it to the
  //right place
  for(ZFSInt nghbr= m_endCommPeriodic-m_periodicS; nghbr< m_endCommPeriodic; nghbr++) {
    ZFSInt k2, j2, i2, id2;
    ZFSInt* step1 = m_block->m_cmnctnFlag->stepInfoRCV[nghbr];
    ZFSInt  step2[3];

    ZFSInt* order = m_block->m_cmnctnFlag->orderInfo[nghbr];
    ZFSInt start1[3];
    ZFSInt start2[3];
    ZFSInt end2[3];
    ZFSInt len2[3];
    ZFSInt totalCells=1;
    ZFSInt len1[3];

    for(ZFSInt j=0; j<nDim; j++) {
      len1[j]=m_block->m_cmnctnFlag->endInfoRCVpoints[nghbr][j] - m_block->m_cmnctnFlag->startInfoRCVpoints[nghbr][j];
      if(len1[j]!=0) {
        totalCells*=len1[j];
      }
      //added    check the step for RCV part !!!!!!!!important
      step2[order[j]]=step1[j];
    }

    for(ZFSInt j=0; j<nDim; j++) {
      start2[j]=0;
      end2[j]=len1[j]-1;
      len2[order[j]]=len1[j];

      if(step2[j]<0) {
        ZFSInt dummy=start2[j];
        start2[j]=end2[j];
        end2[j]=dummy;
      }
    }

    ZFSInt* startInfo=m_block->m_cmnctnFlag->startInfoRCVpoints[nghbr];
    ZFSInt* endInfo= m_block->m_cmnctnFlag->endInfoRCVpoints[nghbr];
    ZFSFloat* bufferRcv =m_block->m_cmnctnFlag->m_bufferPointsRcv[nghbr];

    for(ZFSId dim=0; dim<nDim; dim++) {
      k2=start2[2];
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        j2=start2[1];
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          i2=start2[0];
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            start1[order[0]]=i2;
            start1[order[1]]=j2;
            start1[order[2]]=k2;
            id2=dim*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
            cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            m_cells->coordinates[dim][cellId]=bufferRcv[id2];

            i2+=step2[0];
          }
          j2+=step2[1];
        }
        k2+=step2[2];
      }
    }
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::periodicExchange()
{

  //  if(!(m_block->m_nonBlockingComm)|| m_noPeriodicConnections > 1 || m_rotBC.hasRotationPeriodicBC >0){
  for(ZFSId periodicDir=0; periodicDir<nDim; periodicDir++) {
    m_currentPeriodicDirection = periodicDir;
    RECORD_TIMER_START(m_block->m_tgather);
    (this->*gatherPeriodic)();
    RECORD_TIMER_STOP(m_block->m_tgather);

    RECORD_TIMER_START(m_block->m_tsend);
    send();
    RECORD_TIMER_STOP(m_block->m_tsend);

    RECORD_TIMER_START(m_block->m_treceive);
    receive();
    RECORD_TIMER_STOP(m_block->m_treceive);

    RECORD_TIMER_START(m_block->m_tsendWait);
    for(ZFSId nghbr=m_startCommPeriodic; nghbr<m_endCommPeriodic; nghbr++) {
      if(skipPeriodicDirection(nghbr)) {
        continue;
      }
      MPI_Wait(&(m_block->m_cmnctnFlag->mpi_sndRequest[nghbr]), &(m_block->m_cmnctnFlag->mpi_sndStatus[nghbr]));
    }
    RECORD_TIMER_STOP(m_block->m_tsendWait);

    RECORD_TIMER_START(m_block->m_treceiveWait);
    for(ZFSId nghbr=m_startCommPeriodic; nghbr<m_endCommPeriodic; nghbr++) {
      if(skipPeriodicDirection(nghbr)) {
        continue;
      }
      MPI_Wait(&(m_block->m_cmnctnFlag->mpi_rcvRequest[nghbr]), &(m_block->m_cmnctnFlag->mpi_rcvStatus[nghbr]));
    }
    RECORD_TIMER_STOP(m_block->m_treceiveWait);

    RECORD_TIMER_START(m_block->m_tscatter);
    scatter();
    RECORD_TIMER_STOP(m_block->m_tscatter);
  }
  //}else{
  //periodicNonBlockingExchange();
  //}
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::periodicNonBlockingExchange(){
  RECORD_TIMER_START(m_block->m_tgatherAndSend);
  RECORD_TIMER_START(m_block->m_tgatherAndSendWait);
  MPI_Waitall(m_block->m_cmnctnFlag->noNghbrDomainsPeriodic, m_block->m_cmnctnFlag->mpi_sndRequestCells, MPI_STATUSES_IGNORE);
  RECORD_TIMER_STOP(m_block->m_tgatherAndSendWait);
  for(ZFSInt nghbr=m_startCommPeriodic; nghbr < m_endCommPeriodic; nghbr++){
    ZFSInt* startInfo=m_block->m_cmnctnFlag->startInfoSNDcells[nghbr];
    ZFSInt* endInfo= m_block->m_cmnctnFlag->endInfoSNDcells[nghbr];
    ZFSFloat* bufferSnd = m_block->m_cmnctnFlag->m_bufferCellsSnd[nghbr];
    ZFSInt pos=0;
    for(ZFSId var=0; var<PV->noVariables; var++){
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++){
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++){
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++){
            ZFSId cellId = cellIndex(i,j,k); //i +(j+k*m_block->m_nCells[1])*m_block->m_nCells[2];
            bufferSnd[pos]=m_cells->pvariables[var][cellId];
            pos++;
          }
        }
      }
    }
    MPI_Start(&m_block->m_cmnctnFlag->mpi_sndRequestCells[nghbr]);
  }
  RECORD_TIMER_STOP(m_block->m_tgatherAndSend);
  RECORD_TIMER_START(m_block->m_tscatterAndReceive);
  ZFSIntScratchSpace array_of_indices(m_block->m_cmnctnFlag->noNghbrDomainsPeriodic, __CALLING_FUNCTION__, "array_of_indices");
  ZFSInt outcount;
  while(true){ //exit only in case of completed_index == MPI_UNDEFINED
     RECORD_TIMER_START(m_block->m_tscatterWaitSome);
    MPI_Waitsome(m_block->m_cmnctnFlag->noNghbrDomainsPeriodic, &m_block->m_cmnctnFlag->mpi_rcvRequestCells[m_startCommPeriodic], &outcount,array_of_indices.getPointer(), MPI_STATUSES_IGNORE);
    RECORD_TIMER_STOP(m_block->m_tscatterWaitSome);
    if(outcount == MPI_UNDEFINED) { break; }
    for(ZFSInt a = 0; a < outcount; a++) {
      ZFSInt nghbr = array_of_indices[a];
      ZFSInt k2, j2, i2, id2;
      ZFSInt* step1 = m_block->m_cmnctnFlag->stepInfoRCV[nghbr];
      ZFSInt step2[3];
      ZFSInt* order = m_block->m_cmnctnFlag->orderInfo[nghbr];
      ZFSInt start1[3];
      ZFSInt start2[ 3 ];
      ZFSInt end2[ 3 ];
      ZFSInt len2[ 3 ];
      ZFSInt totalCells=1;
      ZFSInt len1[ 3 ];
      for(ZFSInt j=0; j<nDim; j++){
        len1[j]=m_block->m_cmnctnFlag->endInfoRCVcells[nghbr][j] - m_block->m_cmnctnFlag->startInfoRCVcells[nghbr][j];
        if(len1[j]!=0) {
          totalCells*=len1[j];
        }
        //added    check the step for RCV part !!!!!!!!important
        step2[order[j]]=step1[j];
      }
      for(ZFSInt j=0; j<nDim; j++) {
        start2[j]=0;
        end2[j]=len1[j]-1;
        len2[order[j]]=len1[j];
        if(step2[j]<0) {
          ZFSInt dummy=start2[j];
          start2[j]=end2[j];
          end2[j]=dummy;
        }
        //  reversed already in building the bc maps
      }
      ZFSInt* startInfo=m_block->m_cmnctnFlag->startInfoRCVcells[nghbr];
      ZFSInt* endInfo=m_block->m_cmnctnFlag->endInfoRCVcells[nghbr];
      ZFSFloat* bufferRcv =m_block->m_cmnctnFlag->m_bufferCellsRcv[nghbr];
      ZFSInt pos=0;
      for(ZFSId var=0; var<PV->noVariables; var++) {
        k2=start2[2];
        for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
          j2=start2[1];
          for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
            i2=start2[0];
            for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
              start1[order[0]]=i2;
              start1[order[1]]=j2;
              start1[order[2]]=k2;
              id2=var*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
              ZFSId cellId = i +(j+k*m_nCells[1])*m_nCells[2];
              m_cells->pvariables[var][cellId]= bufferRcv[id2];
              i2+=step2[0];
              pos++;
            }
            j2+=step2[1];
          }
          k2+=step2[2];
        }
      }
    }
  }
  RECORD_TIMER_STOP(m_block->m_tscatterAndReceive);
  //open new receive
  RECORD_TIMER_START(m_block->m_treceive);
  MPI_Startall(m_block->m_cmnctnFlag->noNghbrDomainsPeriodic, &(m_block->m_cmnctnFlag->mpi_rcvRequestCells[m_startCommPeriodic]));
  RECORD_TIMER_STOP(m_block->m_treceive);
}

template <ZFSBool isRans>
ZFSBool ZFSStrctrdBndryCnd3D<isRans>::skipPeriodicDirection(ZFSId nghbr) {
  ZFSInt currentDirection = (m_block->m_cmnctnFlag->bcId[nghbr] - 4401)/2;
  return ((ZFSBool)(m_currentPeriodicDirection-currentDirection));
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::gather()
{
  ZFSId cellId;
  for(ZFSInt nghbr=m_startCommPeriodic; nghbr < m_endCommPeriodic; nghbr++) {
    if(skipPeriodicDirection(nghbr)) {
      continue;
    }

    ZFSInt* startInfo=m_block->m_cmnctnFlag->startInfoSNDcells[nghbr];
    ZFSInt* endInfo= m_block->m_cmnctnFlag->endInfoSNDcells[nghbr];
    ZFSFloat* bufferSnd = m_block->m_cmnctnFlag->m_bufferCellsSnd[nghbr];
    ZFSInt pos=0;

    for(ZFSId var=0; var<PV->noVariables; var++) {
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            cellId = i +(j+k*m_block->m_nCells[1])*m_block->m_nCells[2];
            bufferSnd[pos]=m_cells->pvariables[var][cellId];
            pos++;
          }
        }
      }
    }
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::gatherPeriodicRotation()
{
  ZFSId cellId=-1;
  for(ZFSInt nghbr=m_startCommPeriodic; nghbr < m_endCommPeriodic; nghbr++)
    {
      ZFSInt* startInfo=m_block->m_cmnctnFlag->startInfoSNDcells[nghbr];
      ZFSInt* endInfo= m_block->m_cmnctnFlag->endInfoSNDcells[nghbr];
      ZFSFloat* bufferSnd = m_block->m_cmnctnFlag->m_bufferCellsSnd[nghbr];
      ZFSInt pos=0;
      ZFSFloat** correctMatrix;

      //Take correct rotational Matrix
      if(m_block->m_cmnctnFlag->bcId[nghbr] == 4441)
          correctMatrix = m_rotBC.rotationMatrix4001;
      else
          correctMatrix = m_rotBC.rotationMatrix4002;

      for(ZFSId var=0; var<PV->noVariables; var++) {
        int dir=-1;
        if(var==PV->U){dir=0;}
        if(var==PV->V){dir=1;}
        if(var==PV->W){dir=2;}
        switch(var){
        case 0://u
        case 1://v
        case 2:{//w
          for(ZFSInt k=startInfo[2]; k<endInfo[2]; ++k){
            for(ZFSInt j=startInfo[1]; j<endInfo[1]; ++j){
              for(ZFSInt i=startInfo[0]; i<endInfo[0]; ++i){
                cellId = i +(j+k*m_block->m_nCells[1])*m_block->m_nCells[2];
                bufferSnd[pos]=(correctMatrix[dir][0]*m_cells->pvariables[PV->U][cellId]+
                                correctMatrix[dir][1]*m_cells->pvariables[PV->V][cellId]+
                                correctMatrix[dir][2]*m_cells->pvariables[PV->W][cellId]);
                ++pos;
              }
            }
          }
          break;
        }
        default:{
          for(ZFSInt k=startInfo[2]; k<endInfo[2]; ++k){
            for(ZFSInt j=startInfo[1]; j<endInfo[1]; ++j){
              for(ZFSInt i=startInfo[0]; i<endInfo[0]; ++i){
                cellId = i +(j+k*m_block->m_nCells[1])*m_block->m_nCells[2];
                bufferSnd[pos]=m_cells->pvariables[var][cellId];
                ++pos;
              }
            }
          }
          break;
        }
        }
      }
    }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::scatter()
{
  ZFSId cellId;
  //the ordering of the grid points can be different from
  //sending instance ==> reorder it and copy it to the
  //right place
  for(ZFSInt nghbr=m_startCommPeriodic; nghbr < m_endCommPeriodic; nghbr++) {
    if(skipPeriodicDirection(nghbr)) {
      continue;
    }

    ZFSInt k2, j2, i2, id2;
    ZFSInt* step1 = m_block->m_cmnctnFlag->stepInfoRCV[nghbr];
    ZFSInt step2[3];

    ZFSInt* order = m_block->m_cmnctnFlag->orderInfo[nghbr];
    ZFSInt start1[3];
    ZFSInt start2[ 3 ];
    ZFSInt end2[ 3 ];
    ZFSInt len2[ 3 ];
    ZFSInt totalCells=1;
    ZFSInt len1[ 3 ];

    for(ZFSInt j=0; j<nDim; j++) {
      len1[j]=m_block->m_cmnctnFlag->endInfoRCVcells[nghbr][j] - m_block->m_cmnctnFlag->startInfoRCVcells[nghbr][j];

      if(len1[j]!=0) {
        totalCells*=len1[j];
      }

      //added    check the step for RCV part !!!!!!!!important
      step2[order[j]]=step1[j];
    }

    for(ZFSInt j=0; j<nDim; j++) {
      start2[j]=0;
      end2[j]=len1[j]-1;
      len2[order[j]]=len1[j];

      if(step2[j]<0) {
        ZFSInt dummy=start2[j];
        start2[j]=end2[j];
        end2[j]=dummy;
      }
      //  reversed already in building the bc maps
    }
    ZFSInt* startInfo=m_block->m_cmnctnFlag->startInfoRCVcells[nghbr];
    ZFSInt* endInfo=m_block->m_cmnctnFlag->endInfoRCVcells[nghbr];
    ZFSFloat* bufferRcv =m_block->m_cmnctnFlag->m_bufferCellsRcv[nghbr];
    ZFSInt pos=0;

    for(ZFSId var=0; var<PV->noVariables; var++) {
      k2=start2[2];
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        j2=start2[1];
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          i2=start2[0];
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            start1[order[0]]=i2;
            start1[order[1]]=j2;
            start1[order[2]]=k2;
            id2=var*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
            cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            m_cells->pvariables[var][cellId]= bufferRcv[id2];
            i2+=step2[0];
            pos++;
          }
          j2+=step2[1];
        }
        k2+=step2[2];
      }
    }
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::send()
{
  for(ZFSId nghbr=m_startCommPeriodic; nghbr < m_endCommPeriodic; nghbr++) {
    if(skipPeriodicDirection(nghbr)) {
      continue;
    }

      ZFSInt tag = m_block->domainId()+(m_block->m_cmnctnFlag->m_tagHelperSND[nghbr]+6)*m_block->noDomains();
      ZFSInt err = MPI_Isend((void*)&m_block->m_cmnctnFlag->m_bufferCellsSnd[nghbr][0], m_block->m_cmnctnFlag->m_noNghbrDomainCellBufferSizeSnd[nghbr], MPI_DOUBLE, m_block->m_cmnctnFlag->m_sndNghbrId[nghbr], tag, m_zfsStrctrdComm, &m_block->m_cmnctnFlag->mpi_sndRequest[nghbr]);
      if(err) cout << "rank " << m_block->domainId() << " sending throws error(periodic Exchange) " << endl;
  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::receive()
{
  for(ZFSId nghbr=m_startCommPeriodic; nghbr < m_endCommPeriodic; nghbr++) {
    if(skipPeriodicDirection(nghbr)) {
      continue;
    }

      ZFSInt tag = m_block->m_cmnctnFlag->m_rcvNghbrId[nghbr]+(m_block->m_cmnctnFlag->m_tagHelperRCV[nghbr]+6)*m_block->noDomains();

      ZFSInt err = MPI_Irecv((void*)&m_block->m_cmnctnFlag->m_bufferCellsRcv[nghbr][0],m_block->m_cmnctnFlag->m_noNghbrDomainCellBufferSizeRcv[nghbr], MPI_DOUBLE,m_block->m_cmnctnFlag->m_rcvNghbrId[nghbr] ,tag, MPI_COMM_WORLD, &m_block->m_cmnctnFlag->mpi_rcvRequest[nghbr]);
      if(err) cout << "rank " << m_block->domainId() << " sending throws error (periodic exchange) " << endl;
  }
}

template <ZFSBool isRans>
inline ZFSFloat ZFSStrctrdBndryCnd3D<isRans>::pressure(ZFSId cellId){
  return m_cells->pvariables[PV->P][cellId];
}

template <ZFSBool isRans>
inline ZFSFloat ZFSStrctrdBndryCnd3D<isRans>::temperature(ZFSId cellId){
  const ZFSFloat gamma = m_block->m_gamma;
  ZFSFloat t = gamma*m_cells->pvariables[PV->P][cellId]/m_cells->pvariables[PV->RHO][cellId];
  return t;
}

template <ZFSBool isRans>
inline ZFSFloat ZFSStrctrdBndryCnd3D<isRans>::pressure(ZFSId i, ZFSId j, ZFSId k){
  ZFSId cellId = cellIndex(i, j, k);
  return pressure(cellId);
}

template <ZFSBool isRans>
inline ZFSId ZFSStrctrdBndryCnd3D<isRans>::pointIndex(ZFSInt i, ZFSInt j, ZFSInt k)
{
   return i+(j+k*m_nPoints[1])*m_nPoints[2];
}

template <ZFSBool isRans>
inline ZFSId ZFSStrctrdBndryCnd3D<isRans>::getPointIdFromCell( ZFSInt i, ZFSInt j, ZFSInt k )
{
  return i + (j + k * m_nPoints[1] ) * m_nPoints[2];
}

template <ZFSBool isRans>
inline ZFSId ZFSStrctrdBndryCnd3D<isRans>::getPointIdfromPoint( ZFSId origin, ZFSInt incI,
                                                        ZFSInt incJ, ZFSInt incK )
{
  return origin + incI + incJ * m_nPoints[2] + incK * m_nPoints[2] * m_nPoints[1];
}

/**
 *     New function to compute the skin friction and pressure coefficient
 *     and the part for the lift and drag coefficient
 *     works in all directions and offers line averaging in one direction
 *     @author: Pascal Meysonnat, Marian Albers
 *     @date: 01.01.1010  
 */
template <ZFSBool isRans>
template <ZFSBool computePower>
void ZFSStrctrdBndryCnd3D<isRans>::computeFrictionPressureCoef(){
  //cf=nu*du/dy/(rho*u_8*u_8*0.5)
  //cp=2*(p-p_8)/(rho8*u_8*u_8*0.5)
  const ZFSFloat stagnationPressure = PV->UInfinity*CV->rhoUInfinity*F1B2;
  const ZFSFloat stagnationEnergy =  POW2(PV->UInfinity)*CV->rhoUInfinity;
  //get the number of walls
  ZFSId noWalls=0;
  for(ZFSUint j=0; j<m_globalStrctrdBndryMaps.size(); ++j){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_globalStrctrdBndryMaps[j]->BC)/1000.0);
    if(firstDigit==1){ //count all walls (1000,1004,etc)
      noWalls++;
    }
  }

  //for moving grids and therefore for the power computation
  ZFSFloat fdt = F0;
  if(computePower){
     if (m_block->m_RKStep!=0) {
       fdt = F1/(m_block->m_timeStep*m_block->m_RKalpha[m_block->m_RKStep-1]);
     } else {
       fdt = F1/(m_block->m_timeStep*m_block->m_RKalpha[m_block->m_noRKSteps-1]);
     }
  }

  //reset values to zero for the calculation
  for(ZFSId i =0; i< (ZFSInt)3*noWalls; ++i){
    //lift values are set to zero
    if(m_bCl) {
      m_cL[i]=F0;
      m_cLp[i]=F0;
      m_cLv[i]=F0;
    }
    //drag values are set to zero
    if(m_bCd) {
      m_cD[i]=F0;
      m_cDp[i]=F0;
      m_cDv[i]=F0;
    }
    //power values are set to zero
    if(computePower){
      m_Powerv[i]=F0;
      m_Powerp[i]=F0;
    }
  }
  //set surface area to zero
  for(ZFSId i =0; i< (ZFSInt)noWalls; ++i){ m_cArea[i]=F0; }
 
  //loop over all maps 
  for(ZFSId map=0; map<(ZFSId)m_auxDataMap.size();++map){
    const ZFSInt mapOffsetCf = m_cells->cfOffsets[map];
    const ZFSInt mapOffsetCp = m_cells->cpOffsets[map];
    const ZFSInt mapOffsetPower = m_cells->powerOffsets[map];
    ZFSInt* start = m_auxDataMap[map]->start1;
    ZFSInt* end = m_auxDataMap[map]->end1;
    
    ZFSInt n=0;
    ZFSId normal = 0;
    //area
    ZFSFloat area = F0;    
    ZFSFloat tempArea = F0;
    //skin-friction
    ZFSFloat cf[3]={F0,F0,F0};
    ZFSFloat tempCf[3]={F0,F0,F0};
    //pressure coefficient
    ZFSFloat cp[3]={F0,F0,F0};
    ZFSFloat tempCp[3]={F0,F0,F0};
    //power computation
    ZFSFloat powerp[3]={F0,F0,F0};
    ZFSFloat powerv[3]={F0,F0,F0};
    ZFSFloat tempPowerV[3]={F0,F0,F0};
    ZFSFloat tempPowerP[3]={F0,F0,F0};

    //indices
    ZFSId n10[3]={0,0,0};
    ZFSId n01[3]={0,0,0};
    ZFSInt n1m1[3]={0,0,0};
    ZFSId pCoordDir[9]={0,0,0,0,0,0,0,0,0};
    ZFSId firstTangential = -1, secondTangential = -1;
    ZFSId i=0, k=0, j=0;
    ZFSInt jj = 0, kk = 0;
    ZFSId *reali=nullptr, *realj=nullptr, *realk=nullptr;
    ZFSInt *realjj=nullptr, *realkk=nullptr;
    ZFSInt sizeJ = 0, sizeK = 0;

    ZFSFloat dxidx=-F1, dxidy=-F1, dxidz=-F1;
    ZFSFloat detadx=-F1, detady=-F1, detadz=-F1;
    ZFSFloat dzetadx=-F1, dzetady=-F1, dzetadz=-F1;

    ZFSFloat supportVec[3] = {F0,F0,F0};
    ZFSFloat firstVec[3] = {F0,F0,F0};
    ZFSFloat secondVec[3] = {F0,F0,F0};
    ZFSFloat normalVec[3] = {F0,F0,F0};
    ZFSFloat cellVec[3] = {F0,F0,F0};

    //determine the wall 
    if((m_auxDataMap[map]->face)%2==0){
      n=-1; //bottom wall
      normal = m_auxDataMap[map]->face/2;
      i = start[normal];
    } else {
      n=1; //top wall
      normal = (m_auxDataMap[map]->face - 1)/2;
      i = end[normal]-1;
    }
    //determine the normals
    n1m1[normal] = -1*n;
    n10[normal] = (ZFSId)(0.5+(0.5*(ZFSFloat)n));
    n01[normal] = (ZFSId)(-0.5+(0.5*(ZFSFloat)n));

    switch(normal) {
    case 0:
      {
        if(m_cpAveragingDir == 2) {
          firstTangential = 2;
          secondTangential = 1;
          //pointers to the correct direction counter
          reali = &i;
          realj = &k;
          realk = &j;
          realjj = &kk;
          realkk = &jj;
        } else if(m_cpAveragingDir == 1) {
          firstTangential = 1;
          secondTangential = 2;

          reali = &i;
          realj = &j;
          realk = &k;
          realjj = &jj;
          realkk = &kk;
        } else {
          zfsTerm(1, __CALLING_FUNCTION__, "Cant average in normal direction");
        }
        //Saves normal directions for three other points of surface
        pCoordDir[0] = 0; pCoordDir[1] = 1; pCoordDir[2] = 0;
        pCoordDir[3] = 0; pCoordDir[4] = 0; pCoordDir[5] = 1;
        pCoordDir[6] = 0; pCoordDir[7] = 1; pCoordDir[8] = 1;
        //determine the sizes
        sizeJ =end[1]-start[1];
        sizeK =end[2]-start[2];
        break;
      }
    case 1:
      {
        if(m_cpAveragingDir == 2) {
          firstTangential = 2;
          secondTangential = 0;
          reali = &k;
          realj = &i;
          realk = &j;
          realjj = &kk;
          realkk = &jj;
        } else if(m_cpAveragingDir == 0) {
          firstTangential = 0;
          secondTangential = 2;
          reali = &j;
          realj = &i;
          realk = &k;
          realjj = &jj;
          realkk = &kk;
        } else {
          zfsTerm(1, __CALLING_FUNCTION__, "Cant average in normal direction");
        }

        sizeJ =end[0]-start[0];
        sizeK =end[2]-start[2];

        pCoordDir[0] = 1; pCoordDir[1] = 0; pCoordDir[2] = 0;
        pCoordDir[3] = 0; pCoordDir[4] = 0; pCoordDir[5] = 1;
        pCoordDir[6] = 1; pCoordDir[7] = 0; pCoordDir[8] = 1;

        break;
      }
    case 2:
      {
        if(m_cpAveragingDir == 1) {
          firstTangential = 1;
          secondTangential = 0;

          reali = &k;
          realj = &j;
          realk = &i;
          realjj = &kk;
          realkk = &jj;
        } else if(m_cpAveragingDir == 0) {
          firstTangential = 0;
          secondTangential = 1;

          reali = &j;
          realj = &k;
          realk = &i;
          realjj = &jj;
          realkk = &kk;
        } else {
          zfsTerm(1, __CALLING_FUNCTION__, "Cant average in normal direction");
        }

        pCoordDir[0] = 1; pCoordDir[1] = 0; pCoordDir[2] = 0;
        pCoordDir[3] = 0; pCoordDir[4] = 1; pCoordDir[5] = 0;
        pCoordDir[6] = 1; pCoordDir[7] = 1; pCoordDir[8] = 0;

        sizeJ =end[0]-start[0];
        sizeK =end[1]-start[1];
        break;
      }
    default:
      {
        zfsTerm(1, __CALLING_FUNCTION__, "Normal direction not implemented");
      }
    }

    //for power computation
    ZFSFloat gridVel[3] = {F0,F0,F0}; //grid velocity

    //main loop over the two tangential directions of the wall
    for(k=start[secondTangential]; k<end[secondTangential] ; k++) {
      jj= 0;
      //this loop always goes into the line averaging direction
      for(j=start[firstTangential]; j<end[firstTangential]; j++) {
        //get id of boundary cell and cell above that one for extrapolation
        const ZFSId cellId = cellIndex(*reali,*realj,*realk);
        const ZFSId cellIdP1 = cellIndex(*reali+n1m1[0],*realj+n1m1[1],*realk+n1m1[2]);

        //get point id and ids of three other surface corners
        const ZFSId pIJK  = getPointIdFromCell(*reali+n10[0],*realj+n10[1],*realk+n10[2]);
        const ZFSId pIJPK = getPointIdfromPoint(pIJK,pCoordDir[0],pCoordDir[1],pCoordDir[2]);
        const ZFSId pIJKP = getPointIdfromPoint(pIJK,pCoordDir[3],pCoordDir[4],pCoordDir[5]);
        const ZFSId pIJPKP= getPointIdfromPoint(pIJK,pCoordDir[6],pCoordDir[7],pCoordDir[8]);

        //first get the position of the wall (reference!!)
        //x is always used as coordinate in normal direction
        const ZFSFloat xRef=F1B4*(m_coordinates[normal][pIJK]+m_coordinates[normal][pIJPK]+m_coordinates[normal][pIJKP]+m_coordinates[normal][pIJPKP]);
        const ZFSFloat yRef=F1B4*(m_coordinates[firstTangential][pIJK]+m_coordinates[firstTangential][pIJPK]+
                                  m_coordinates[firstTangential][pIJKP]+m_coordinates[firstTangential][pIJPKP]);
        const ZFSFloat zRef=F1B4*(m_coordinates[secondTangential][pIJK]+m_coordinates[secondTangential][pIJPK]+
                                  m_coordinates[secondTangential][pIJKP]+m_coordinates[secondTangential][pIJPKP]);

        //compute the pressure coefficient
        //-->
        //get the distcance
        const ZFSFloat dx2=sqrt(POW2(m_cells->coordinates[normal][cellIdP1]-m_cells->coordinates[normal][cellId])+
                               POW2(m_cells->coordinates[firstTangential][cellIdP1]-m_cells->coordinates[firstTangential][cellId])+
                               POW2(m_cells->coordinates[secondTangential][cellIdP1]-m_cells->coordinates[secondTangential][cellId]));
        const ZFSFloat dx1=sqrt(POW2(m_cells->coordinates[normal][cellId]-xRef)+
                                POW2(m_cells->coordinates[firstTangential][cellId]-yRef)+
                                POW2(m_cells->coordinates[secondTangential][cellId]-zRef));
        //pressures
        const ZFSFloat p1 = m_cells->pvariables[PV->P][cellId];
        const ZFSFloat p2 = m_cells->pvariables[PV->P][cellIdP1];
        //extrapolation to the face
        const ZFSFloat pW=((p1-p2)/dx2)*dx1+p1;
        const ZFSFloat cpn = (pW-PV->PInfinity)/stagnationPressure;
        //save pressure to map
        m_cells->cp[mapOffsetCp+(*realjj) + (*realkk)*sizeJ] = cpn;

        //for the power consumption due to pressre:
        if(computePower){
          if(m_block->m_movingGrid){
            const ZFSFloat xRefOld=F1B4*(m_block->m_mgOldCoordinates[normal][pIJK]+
                                         m_block->m_mgOldCoordinates[normal][pIJPK]+
                                         m_block->m_mgOldCoordinates[normal][pIJKP]+
                                         m_block->m_mgOldCoordinates[normal][pIJPKP]);
            const ZFSFloat yRefOld=F1B4*(m_block->m_mgOldCoordinates[firstTangential][pIJK]+
                                         m_block->m_mgOldCoordinates[firstTangential][pIJPK]+
                                         m_block->m_mgOldCoordinates[firstTangential][pIJKP]+
                                         m_block->m_mgOldCoordinates[firstTangential][pIJPKP]);
            const ZFSFloat zRefOld=F1B4*(m_block->m_mgOldCoordinates[secondTangential][pIJK]+
                                         m_block->m_mgOldCoordinates[secondTangential][pIJPK]+
                                         m_block->m_mgOldCoordinates[secondTangential][pIJKP]+
                                         m_block->m_mgOldCoordinates[secondTangential][pIJPKP]);
            //determine grid velocity in cartesian components
            if(m_block->m_time > F0){
              gridVel[0]=(xRef-xRefOld)*fdt;
              gridVel[1]=(yRef-yRefOld)*fdt;
              gridVel[2]=(zRef-zRefOld)*fdt;
            }
          }else{
            gridVel[0]=PV->UInfinity;
            gridVel[1]=PV->VInfinity;
            gridVel[2]=PV->WInfinity;
          }
          //compute the power (p*v)
          const ZFSFloat dp= (pW-PV->PInfinity);
          
          const ZFSFloat P_px=(-1.0)*dp*m_cells->surfaceMetrics[cellIndex(*reali+n01[0],*realj+n01[1],*realk+n01[2])][nDim*normal + 0]*gridVel[0];
          const ZFSFloat P_py=(-1.0)*dp*m_cells->surfaceMetrics[cellIndex(*reali+n01[0],*realj+n01[1],*realk+n01[2])][nDim*normal + 1]*gridVel[1];
          const ZFSFloat P_pz=(-1.0)*dp*m_cells->surfaceMetrics[cellIndex(*reali+n01[0],*realj+n01[1],*realk+n01[2])][nDim*normal + 2]*gridVel[2];
          //save power due to pressure to map
          m_cells->powerPres[mapOffsetPower+0*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_px/stagnationEnergy;
          m_cells->powerPres[mapOffsetPower+1*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_py/stagnationEnergy;
          m_cells->powerPres[mapOffsetPower+2*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_pz/stagnationEnergy;
        }
        //<--- finished computation of cp

        //compute the skin-friction coefficient
        //-->
        //compute orthogonal distance surface-plane to cell center
        for(ZFSId dim = 0; dim<nDim; dim++) {
          supportVec[dim] = m_coordinates[dim][pIJK];
          firstVec[dim] = m_coordinates[dim][pIJPK] - m_coordinates[dim][pIJK];
          secondVec[dim] = m_coordinates[dim][pIJKP] - m_coordinates[dim][pIJK];
          cellVec[dim] = m_cells->coordinates[dim][cellId];
        }

        crossProduct(normalVec, firstVec, secondVec);
        const ZFSFloat normalLength = sqrt(POW2(normalVec[0]) + POW2(normalVec[1]) + POW2(normalVec[2]));
        ZFSFloat orthDist = F0;

        for(ZFSId dim = 0; dim<nDim; dim++) {
          orthDist += (cellVec[dim] - supportVec[dim])*normalVec[dim];
        }

        orthDist = fabs(orthDist/normalLength);

        //extrapolate the Temperature to the wall
        const ZFSFloat T1 = temperature(cellId);
        const ZFSFloat T2 = temperature(cellIdP1);

        //extrapolation to the face
        const ZFSFloat tBc=((T1-T2)/dx2)*dx1+T1;
        const ZFSFloat nu=zfsSUTHERLANDLAW(tBc);

        //velocities
        const ZFSFloat u1=m_cells->pvariables[PV->U][cellId];
        const ZFSFloat v1=m_cells->pvariables[PV->V][cellId];
        const ZFSFloat w1=m_cells->pvariables[PV->W][cellId];

        //now use cell metrics to transform to contravariant velocities
        //first get metrics
        dxidx = m_cells->cellMetrics[cellId][0];
        dxidy = m_cells->cellMetrics[cellId][1];
        dxidz = m_cells->cellMetrics[cellId][2];

        detadx = m_cells->cellMetrics[cellId][3];
        detady = m_cells->cellMetrics[cellId][4];
        detadz = m_cells->cellMetrics[cellId][5];

        dzetadx = m_cells->cellMetrics[cellId][6];
        dzetady = m_cells->cellMetrics[cellId][7];
        dzetadz = m_cells->cellMetrics[cellId][8];

        const ZFSFloat invJac = (dxidx*(detady*dzetadz - detadz*dzetady)-
                                 detadx*(dxidy*dzetadz - dxidz*dzetady)+
                                 dzetadx*(dxidy*detadz - dxidz*detady));
        const ZFSFloat det = F1/sqrt( invJac );

        //compute contravariant velocities parallel to the wall
        const ZFSFloat du1 = (dxidx*u1 + dxidy*v1 + dxidz*w1) - m_cells->dxt[0][cellId];
        const ZFSFloat dv1 = (detadx*u1 + detady*v1 + detadz*w1) - m_cells->dxt[1][cellId];
        const ZFSFloat dw1 = (dzetadx*u1 + dzetady*v1 + dzetadz*w1) - m_cells->dxt[2][cellId];

        // using taylor expansion a second-order approximation
        // of du/dn can be achieved at the wall
        // const ZFSFloat u2=m_cells->pvariables[PV->U][cellIdP1];
        // const ZFSFloat v2=m_cells->pvariables[PV->V][cellIdP1];
        // const ZFSFloat w2=m_cells->pvariables[PV->W][cellIdP1];
        // const ZFSFloat du2 = (dxidx*u2 + dxidy*v2 + dxidz*w2) - m_cells->dxt[0][cellId];
        // const ZFSFloat dv2 = (detadx*u2 + detady*v2 + detadz*w2) - m_cells->dxt[1][cellId];
        // const ZFSFloat dw2 = (dzetadx*u2 + dzetady*v2 + dzetadz*w2) - m_cells->dxt[2][cellId];
        // const ZFSFloat dx3 = dx1+dx2;
        // const ZFSFloat dudx= (du1*dx3*dx3 - du2*dx1*dx1)/(dx1*dx3*dx3 - dx3*dx1*dx1);
        // const ZFSFloat dvdx= (dv1*dx3*dx3 - dv2*dx1*dx1)/(dx1*dx3*dx3 - dx3*dx1*dx1);
        // const ZFSFloat dwdx= (dw1*dx3*dx3 - dw2*dx1*dx1)/(dx1*dx3*dx3 - dx3*dx1*dx1);

        //compute gradient
        const ZFSFloat dudx=du1/orthDist;
        const ZFSFloat dvdx=dv1/orthDist;
        const ZFSFloat dwdx=dw1/orthDist;

        //friction coefficient in computational space
        ZFSFloat cfuu = nu * dudx /stagnationPressure/(m_block->m_Re0);
        ZFSFloat cfvv = nu * dvdx /stagnationPressure/(m_block->m_Re0);
        ZFSFloat cfww = nu * dwdx /stagnationPressure/(m_block->m_Re0);

        const ZFSFloat dxdxi = (detady*dzetadz - dzetady*detadz) * det;
        const ZFSFloat dxdeta = (dzetady*dxidz - dxidy*dzetadz) * det;
        const ZFSFloat dxdzeta = (dxidy*detadz - detady*dxidz) * det;

        const ZFSFloat dydxi = (dzetadx*detadz - detadx*dzetadz) * det;
        const ZFSFloat dydeta = (dxidx*dzetadz - dzetadx*dxidz) * det;
        const ZFSFloat dydzeta = (detadx*dxidz - dxidx*detadz) * det;

        const ZFSFloat dzdxi = (detadx*dzetady - dzetadx*detady) * det;
        const ZFSFloat dzdeta = (dzetadx*dxidy - dxidx*dzetady) * det;
        const ZFSFloat dzdzeta = (dxidx*detady - detadx*dxidy) * det;


        //now transform back to physical space
        //skin-friction
        const ZFSFloat cfx = (cfuu*dxdxi + cfvv*dxdeta + cfww*dxdzeta) * det;
        const ZFSFloat cfy = (cfuu*dydxi + cfvv*dydeta + cfww*dydzeta) * det;
        const ZFSFloat cfz = (cfuu*dzdxi + cfvv*dzdeta + cfww*dzdzeta) * det;

        //save to map
        //skin-friction
        m_cells->cf[mapOffsetCf+0*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=cfx;
        m_cells->cf[mapOffsetCf+1*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=cfy;
        m_cells->cf[mapOffsetCf+2*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=cfz;

        //Power
        //Power due to viscous forces in the computational space ( tau_w*n*u)
        if(computePower){
          ZFSFloat P_cfuu=nu*dudx/(m_block->m_Re0);
          ZFSFloat P_cfvv=nu*dvdx/(m_block->m_Re0);
          ZFSFloat P_cfww=nu*dwdx/(m_block->m_Re0);
          /*if(m_block->m_movingGrid){
            P_cfuu= nu * dudx * m_cells->dxt[0][cellId]/(m_block->m_Re0);
            P_cfvv= nu * dvdx * m_cells->dxt[1][cellId]/(m_block->m_Re0);
            P_cfww= nu * dwdx * m_cells->dxt[2][cellId]/(m_block->m_Re0);
          }else{
            const ZFSFloat u8 = (dxidx*PV->UInfinity + dxidy*PV->VInfinity + dxidz*PV->WInfinity);
            const ZFSFloat v8 = (detadx*PV->UInfinity + detady*PV->VInfinity + detadz*PV->WInfinity);
            const ZFSFloat w8 = (dzetadx*PV->UInfinity + dzetady*PV->VInfinity + dzetadz*PV->WInfinity);
            cout << "u8 " << u8 << " v8 " << v8 << " w8 " << w8 << " PV->INF " << PV->UInfinity << " " << PV->VInfinity << " " << PV->WInfinity << " det " << det <<  endl; 
            const ZFSFloat backU = (u8*dxdxi + v8*dxdeta + w8*dxdzeta)* det;
            const ZFSFloat backV = (u8*dydxi + v8*dydeta + w8*dydzeta)* det;
            const ZFSFloat backW = (u8*dzdxi + v8*dzdeta + w8*dzdzeta)* det;
            cout << "back: u8 " << backU << " v8 " << backV << " w8 " << backW << endl; 
            P_cfuu= nu * dudx/(m_block->m_Re0);
            P_cfvv= nu * dvdx/(m_block->m_Re0);
            P_cfww= nu * dwdx/(m_block->m_Re0);
          }*/
          //Power due to viscous forces in the computational space ( tau_w*n*u)
          const ZFSFloat P_cfx = (P_cfuu*dxdxi + P_cfvv*dxdeta + P_cfww*dxdzeta) * det*gridVel[0];//*PV->UInfinity;
          const ZFSFloat P_cfy = (P_cfuu*dydxi + P_cfvv*dydeta + P_cfww*dydzeta) * det*gridVel[1];//*PV->VInfinity;
          const ZFSFloat P_cfz = (P_cfuu*dzdxi + P_cfvv*dzdeta + P_cfww*dzdzeta) * det*gridVel[2];//*PV->WInfinity;
          //save power to map
          m_cells->powerVisc[mapOffsetPower+0*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_cfx/stagnationEnergy;
          m_cells->powerVisc[mapOffsetPower+1*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_cfy/stagnationEnergy;
          m_cells->powerVisc[mapOffsetPower+2*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_cfz/stagnationEnergy;
        }
        //<-- finished computation of skin-friction
        
        //beware: although called dxidxyz the surface metrics in normal direction are taken,
        //so it can be detadxyz or dzetadxyz depening on the surface normal
        dxidx = m_cells->surfaceMetrics[cellIndex(*reali+n01[0],*realj+n01[1],*realk+n01[2])][nDim*normal + 0];
        dxidy = m_cells->surfaceMetrics[cellIndex(*reali+n01[0],*realj+n01[1],*realk+n01[2])][nDim*normal + 1];
        dxidz = m_cells->surfaceMetrics[cellIndex(*reali+n01[0],*realj+n01[1],*realk+n01[2])][nDim*normal + 2];
        
        ZFSFloat considerValue = F1;
        if(m_block->m_auxDataCoordinateLimits) {
          if(m_block->m_auxDataLimits[0] <= zRef && zRef <= m_block->m_auxDataLimits[1] &&
             m_block->m_auxDataLimits[2] <= yRef && yRef <= m_block->m_auxDataLimits[3]) {
            considerValue = F1;
          } else {
            considerValue = F0;
          }
        }
        
        //this should be the valid method as in TFS
        //cp:
        //sum up all lengths and all cps with their surface width contribution
        tempCp[0]+=(-1.0)*cpn*m_cells->surfaceMetrics[cellIndex(*reali+n01[0],*realj+n01[1],*realk+n01[2])][nDim*normal + 0]*considerValue;
        tempCp[1]+=(-1.0)*cpn*m_cells->surfaceMetrics[cellIndex(*reali+n01[0],*realj+n01[1],*realk+n01[2])][nDim*normal + 1]*considerValue;
        tempCp[2]+=(-1.0)*cpn*m_cells->surfaceMetrics[cellIndex(*reali+n01[0],*realj+n01[1],*realk+n01[2])][nDim*normal + 2]*considerValue;
        //cf
        tempCf[0]+= cfx*(sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz))*considerValue;
        tempCf[1]+= cfy*(sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz))*considerValue;
        tempCf[2]+= cfz*(sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz))*considerValue;
        //area
	tempArea+= sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz)*considerValue;
        //pressure power
        if(computePower) {
          tempPowerP[0]+=m_cells->powerPres[mapOffsetPower+0*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*considerValue;
          tempPowerP[1]+=m_cells->powerPres[mapOffsetPower+1*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*considerValue;
          tempPowerP[2]+=m_cells->powerPres[mapOffsetPower+2*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*considerValue;
          //viscous power
          //this should be the valid method as in TFS (this is the t_w*n*u*da)
          tempPowerV[0]+=  m_cells->powerVisc[mapOffsetPower+0*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*(sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz))*considerValue;
          tempPowerV[1]+=  m_cells->powerVisc[mapOffsetPower+1*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*(sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz))*considerValue;
          tempPowerV[2]+= m_cells->powerVisc[mapOffsetPower+2*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*(sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz))*considerValue;
        }

        if(m_block->m_detailAuxData) {
          //save surface contributions
          m_cells->cf[mapOffsetCf+3*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]= dxidx;
          m_cells->cf[mapOffsetCf+4*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]= dxidy;
          m_cells->cf[mapOffsetCf+5*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]= dxidz;

          //save surface coordinates
          m_cells->cf[mapOffsetCf+6*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=zRef;
          m_cells->cf[mapOffsetCf+7*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=xRef;
          m_cells->cf[mapOffsetCf+8*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=yRef;
        }

        ++jj;
      }

      cp[0]+=tempCp[0];
      cp[1]+=tempCp[1];
      cp[2]+=tempCp[2];
      cf[0]+=tempCf[0];
      cf[1]+=tempCf[1];
      cf[2]+=tempCf[2];

      if(computePower) {
        powerv[0]+=tempPowerV[0];
        powerv[1]+=tempPowerV[1];
        powerv[2]+=tempPowerV[2];
        powerp[0]+=tempPowerP[0];
        powerp[1]+=tempPowerP[1];
        powerp[2]+=tempPowerP[2];

        for(int m = 0; m<nDim; ++m) {
          tempPowerV[m] =F0;
          tempPowerP[m] =F0;
        }
      }

      area+=tempArea;

      for(int m = 0; m<nDim; ++m) {
        tempCp[m] = F0;
        tempCf[m] = F0;
      }
      tempArea = F0;

      ++kk;
    }

    //now add to lift and drag coefficients
    ZFSId count=0;
    for(ZFSUint id=0; id<m_globalStrctrdBndryMaps.size(); ++id){
      //check if its a wall map
      ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_globalStrctrdBndryMaps[id]->BC)/1000.0);
      if(firstDigit!=1) {
        continue;
      }

      if(m_auxDataMap[map]->Id2==m_globalStrctrdBndryMaps[id]->Id2){
        if(m_bCl) {
          m_cLv[count*3+0]+=F0;
          m_cLv[count*3+1]+=cf[1];
          m_cLv[count*3+2]+=F0;
          m_cLp[count*3+0]+=F0;
          m_cLp[count*3+1]+=cp[1];
          m_cLp[count*3+2]+=F0;
        }
        if(m_bCd) {
          m_cDv[count*3+0]+=cf[0];
          m_cDv[count*3+1]+=F0;
          m_cDv[count*3+2]+=F0;
          m_cDp[count*3+0]+=cp[0];
          m_cDp[count*3+1]+=F0;
          m_cDp[count*3+2]+=F0;
        }
        if(computePower){
          m_Powerp[count*3+0]+=powerp[0];
          m_Powerp[count*3+1]+=powerp[1];
          m_Powerp[count*3+2]+=powerp[2];
          m_Powerv[count*3+0]+=powerv[0];
          m_Powerv[count*3+1]+=powerv[1];
          m_Powerv[count*3+2]+=powerv[2];
        }
	m_cArea[count]+=area;
      }

      count++;
    }
  }
}
template void ZFSStrctrdBndryCnd3D<true>::computeFrictionPressureCoef<true>();
template void ZFSStrctrdBndryCnd3D<false>::computeFrictionPressureCoef<true>();
template void ZFSStrctrdBndryCnd3D<true>::computeFrictionPressureCoef<false>();
template void ZFSStrctrdBndryCnd3D<false>::computeFrictionPressureCoef<false>();



template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::computeAuxData()
{
  if(m_bCfCpCoeff && m_bPower){
    this->computeFrictionPressureCoef<true>();
  }else{
    this->computeFrictionPressureCoef<false>();
  }
//if(m_bCp) computePressureCoef();
//  if(m_bCf) computeSkinFriction();
  if(m_bCl) computeLiftCoef();
  if(m_bCd) computeDragCoef();
}

/**
 *     Function to compute the lift coefficient cl, split
 *     split into the viscous part cLv and the pressure part cLp
 *     @author: Marian Albers/Pascal Meysonnat
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::computeLiftCoef(){
  ZFSInt noWalls=0;
  for(ZFSUint j=0; j<m_globalStrctrdBndryMaps.size(); ++j){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_globalStrctrdBndryMaps[j]->BC)/1000.0);
    if(firstDigit==1){
      noWalls++;
    }
  }

  ZFSFloatScratchSpace test1(3*noWalls, __CALLING_FUNCTION__, "test1");
  ZFSFloatScratchSpace test2(3*noWalls, __CALLING_FUNCTION__, "test2");

  for(ZFSId i=0; i<3*noWalls; ++i) {
    test1[i]=m_cLp[i];
    test2[i]=m_cLv[i];
  }

  MPI_Allreduce(test1.begin(), m_cLp, 3*noWalls, MPI_DOUBLE, MPI_SUM, m_zfsStrctrdComm);
  MPI_Allreduce(test2.begin(), m_cLv, 3*noWalls, MPI_DOUBLE, MPI_SUM, m_zfsStrctrdComm);

  if(m_bCpLineAveraging) {
    for(ZFSId i=0; i<3*noWalls; ++i) {
      m_cLp[i] /= m_block->m_globalDomainWidth;
      m_cLv[i] /= m_block->m_globalDomainWidth;
    }
  }
}

/**
 *     Function to compute the drag coefficient, split
 *     into the vicous part cDv and the pressure part cDp
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::computeDragCoef(){
  ZFSInt noWalls=0;
  for(ZFSUint j=0; j<m_globalStrctrdBndryMaps.size(); ++j){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_globalStrctrdBndryMaps[j]->BC)/1000.0);
    if(firstDigit==1){
      noWalls++;
    }
  }

  ZFSFloatScratchSpace test1(3*noWalls, __CALLING_FUNCTION__, "test1");
  ZFSFloatScratchSpace test2(3*noWalls, __CALLING_FUNCTION__, "test2");

  test1.fill(F0);
  test2.fill(F0);

  for(ZFSId i=0; i<3*noWalls; ++i) {
    test1[i]=m_cDp[i];
    test2[i]=m_cDv[i];
  }

  MPI_Allreduce(test1.begin(), m_cDp, 3*noWalls, MPI_DOUBLE, MPI_SUM, m_zfsStrctrdComm);
  MPI_Allreduce(test2.begin(), m_cDv, 3*noWalls, MPI_DOUBLE, MPI_SUM, m_zfsStrctrdComm);

  if(m_bCpLineAveraging) {
    for(ZFSId i=0; i<3*noWalls; ++i) {
      m_cDp[i] /= m_block->m_globalDomainWidth;
      m_cDv[i] /= m_block->m_globalDomainWidth;
    }
  }
}

/**  compute Power consumption and split it into viscous and pressure part
 *   @author Pascal Meysonnat
 *   @ 01.01.2018
 *
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::computePowerCoefRoot(){
  ZFSInt noWalls=0;
  for(ZFSUint j=0; j<m_globalStrctrdBndryMaps.size(); ++j){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_globalStrctrdBndryMaps[j]->BC)/1000.0);
    if(firstDigit==1){
      noWalls++;
    }
  }

  ZFSFloatScratchSpace localData(6*noWalls, __CALLING_FUNCTION__, "localData");
  ZFSFloatScratchSpace globalData(6*noWalls, __CALLING_FUNCTION__, "globalData");
  for(ZFSId i=0; i<3*noWalls; ++i) {
    localData[i]=m_Powerp[i];
    localData[3*noWalls+i]=m_Powerv[i];
  }

   MPI_Reduce(localData.begin(), globalData.begin(), 6*noWalls, MPI_DOUBLE, MPI_SUM, 0, m_zfsStrctrdComm);
   
   for(ZFSId i=0; i<3*noWalls; ++i) {
    m_Powerp[i] = globalData[i];
    m_Powerv[i] = globalData[3*noWalls+i];
  }

  if(m_bCpLineAveraging) {
    for(ZFSId i=0; i<3*noWalls; ++i) {
      m_Powerp[i] /= m_block->m_globalDomainWidth;
      m_Powerv[i] /= m_block->m_globalDomainWidth;
    }
  }
}

/**
 *     Function to compute the lift coefficient cl, split
 *     split into the viscous part cLv and the pressure part cLp
 *     The ROOT version is faster due to an MPI_Reduce instead
 *     of an MPI_Allreduce, but only root rank has data
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::computeLiftCoefRoot(){
  ZFSInt noWalls=0;
  for(ZFSUint j=0; j<m_globalStrctrdBndryMaps.size(); ++j){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_globalStrctrdBndryMaps[j]->BC)/1000.0);
    if(firstDigit==1){
      noWalls++;
    }
  }

  ZFSFloatScratchSpace localData(6*noWalls, __CALLING_FUNCTION__, "localData");
  ZFSFloatScratchSpace globalData(6*noWalls, __CALLING_FUNCTION__, "globalData");

  for(ZFSId i=0; i<3*noWalls; ++i) {
    localData[i]=m_cLp[i];
    localData[3*noWalls+i]=m_cLv[i];
  }

  MPI_Reduce(localData.begin(), globalData.begin(), 6*noWalls, MPI_DOUBLE, MPI_SUM, 0, m_zfsStrctrdComm);

  for(ZFSId i=0; i<3*noWalls; ++i) {
    m_cLp[i] = globalData[i];
    m_cLv[i] = globalData[3*noWalls+i];
  }


  if(m_bCpLineAveraging) {
    for(ZFSId i=0; i<3*noWalls; ++i) {
      m_cLp[i] /= m_block->m_globalDomainWidth;
      m_cLv[i] /= m_block->m_globalDomainWidth;
    }
  }
}

/**
 *     Function to compute the drag coefficient, split
 *     into the vicous part cDv and the pressure part cDp
 *     The ROOT version is faster due to an MPI_Reduce instead
 *     of an MPI_Allreduce, but only root rank has data
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::computeDragCoefRoot(){
  ZFSInt noWalls=0;
  for(ZFSUint j=0; j<m_globalStrctrdBndryMaps.size(); ++j){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_globalStrctrdBndryMaps[j]->BC)/1000.0);
    if(firstDigit==1){
      noWalls++;
    }
  }

  ZFSFloatScratchSpace localData(7*noWalls, __CALLING_FUNCTION__, "localData");
  ZFSFloatScratchSpace globalData(7*noWalls, __CALLING_FUNCTION__, "globalData");

  for(ZFSId i=0; i<3*noWalls; ++i) {
    localData[i]=m_cDp[i];
    localData[3*noWalls+i]=m_cDv[i];
  }
  for(ZFSId i=0; i<noWalls; ++i) {
    localData[6*noWalls+i] = m_cArea[i];
  }

  MPI_Reduce(localData.begin(), globalData.begin(), 7*noWalls, MPI_DOUBLE, MPI_SUM, 0, m_zfsStrctrdComm);

  for(ZFSId i=0; i<3*noWalls; ++i) {
    m_cDp[i] = globalData[i];
    m_cDv[i] = globalData[3*noWalls+i];
  }
  for(ZFSId i=0; i<noWalls; ++i) {
    m_cArea[i] = globalData[6*noWalls+i];
  }

  if(m_bCpLineAveraging) {
    for(ZFSId i=0; i<3*noWalls; ++i) {
      m_cDp[i] /= m_block->m_globalDomainWidth;
      m_cDv[i] /= m_block->m_globalDomainWidth;
    }
    for(ZFSId i=0; i<noWalls; ++i) {
      m_cArea[i] /= m_block->m_globalDomainWidth;
    }
  }
}

/**
 *     Function to compute the momentum coefficient (LE): cm_le
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd3D<isRans>::computeMomentCoef(){

}

template <ZFSBool isRans>
inline void ZFSStrctrdBndryCnd3D<isRans>::crossProduct( ZFSFloat* result,
                                                        ZFSFloat* vec1,
                                                        ZFSFloat* vec2)
{
  result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
  result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
  result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

template <ZFSBool isRans>
inline ZFSId ZFSStrctrdBndryCnd3D<isRans>::getReverseCellId(ZFSId i, ZFSId j, ZFSId k, ZFSId face) {
  const ZFSId i_new = m_reverseCellIdGC[face*3+0]*(m_noGhostLayers-1)+m_reverseCellIdDim[face*3+0]*i;
  const ZFSId j_new = m_reverseCellIdGC[face*3+1]*(m_noGhostLayers-1)+m_reverseCellIdDim[face*3+1]*j;
  const ZFSId k_new = m_reverseCellIdGC[face*3+2]*(m_noGhostLayers-1)+m_reverseCellIdDim[face*3+2]*k;
  return cellIndex(i_new,j_new,k_new);
}

template <ZFSBool isRans>
inline ZFSId ZFSStrctrdBndryCnd3D<isRans>::getExtrNghbrId(ZFSId cellId, ZFSId face) {
  static const ZFSId cellShift[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};
  const ZFSId n = 1-2*(face%2);
  const ZFSId dim = face/2;
  const ZFSId nghbrId = cellId + n*cellShift[dim];
  return nghbrId;
}

template <ZFSBool isRans>
inline pair<ZFSId, ZFSId> ZFSStrctrdBndryCnd3D<isRans>::getMirrorCellIdPair(ZFSId i, ZFSId j, ZFSId k, ZFSId face) {
  const ZFSId cellShift[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};
  const ZFSId gcPos[6]={m_noGhostLayers, m_nCells[2]-m_noGhostLayers-1, m_noGhostLayers , m_nCells[1]-m_noGhostLayers-1, m_noGhostLayers, m_nCells[0]-m_noGhostLayers-1};
  const ZFSId ijk_new[3] = {m_reverseCellIdGC[face*3+0]*(m_noGhostLayers-1)+m_reverseCellIdDim[face*3+0]*i,
                            m_reverseCellIdGC[face*3+1]*(m_noGhostLayers-1)+m_reverseCellIdDim[face*3+1]*j,
                            m_reverseCellIdGC[face*3+2]*(m_noGhostLayers-1)+m_reverseCellIdDim[face*3+2]*k};
  const ZFSId dim = face/2;
  const ZFSId n = 1-2*(face%2);
  const ZFSId mirror = ((gcPos[face]-ijk_new[dim])*2-n);
  const ZFSId cellId = cellIndex(ijk_new[0],ijk_new[1],ijk_new[2]);
  return make_pair(cellId, cellId + mirror*cellShift[dim]);
}

//instantanisation of the class
template class ZFSStrctrdBndryCnd3D<true>;
template class ZFSStrctrdBndryCnd3D<false>;
