#include <cmath>
#include <iostream>
#include <fstream>

#include "zfsstrctrdblck.h"
#include "zfsstrctrdbndrycnd2d.h"
#include "zfsstrctrdblck2d.h"
#include "zfstypes.h"
#include "zfsiolib.h"
#include "zfspointbox.h"
#include "zfskdtree.h"


using namespace std;

template <ZFSBool isRans>
ZFSStrctrdBndryCnd2D<isRans>::ZFSStrctrdBndryCnd2D( ZFSStrctrdBlck<2>* block, ZFSId noSpecies )
  : ZFSStrctrdBndryCnd<2>( block ), m_noSpecies(noSpecies){

  TRACE();

  const ZFSLong oldAllocatedBytes = allocatedBytes();

  m_block = static_cast<class ZFSStrctrdBlck2D*> (block);
  m_noStrctrdCells=m_block->m_noStrctrdCells;

  m_startCommChannel=m_block->m_cmnctnFlag->noNghbrDomainsNormal;
  m_endCommChannel=m_block->m_cmnctnFlag->noNghbrDomainsNormal+m_block->m_cmnctnFlag->noNghbrDomainsChannel;

  m_startCommPeriodic=m_block->m_cmnctnFlag->noNghbrDomainsNormal+m_block->m_cmnctnFlag->noNghbrDomainsChannel;
  m_endCommPeriodic=m_block->m_cmnctnFlag->noNghbrDomainsNormal+m_block->m_cmnctnFlag->noNghbrDomainsChannel+m_block->m_cmnctnFlag->noNghbrDomainsPeriodic;


  // aux map variables
  m_bCfCpCoeff = m_block->m_bCfCpCoeff;
  m_bCl = m_block->m_bCl;
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

  //read rescaling bc properties
  if(*m_block->m_rescalingCommGrRoot >= 0)
    {
      m_rescalingBLT = 5.95;
      m_rescalingBLT = *(ZFSContext::getProperty("rescalingBLT", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL)->asFloat(0));

      cout << "Rescaling BLT: " << m_rescalingBLT << endl;
    }

  // compute sponge coefficients for each cell
  if(m_block->m_useSponge)
    {
      RECORD_TIMER_START(m_block->m_tbuildSponge);
      readAndDistributeSpongeCoordinates();
      RECORD_TIMER_STOP(m_block->m_tbuildSponge);
    }

  printAllocatedMemory( oldAllocatedBytes, "ZFSStrctrdBndryCnd2D", m_zfsStrctrdComm );
}

template <ZFSBool isRans>
ZFSStrctrdBndryCnd2D<isRans>::~ZFSStrctrdBndryCnd2D()
{
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::readAndDistributeSpongeCoordinates(){
  TRACE();

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
    for(ZFSInt dim=0; dim<nDim; ++dim){
      size*=(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]+1);
    }
    memSize+=size;
  }

  //1.2) allocate the memory
  ZFSFloatScratchSpace coordMem(nDim*memSize, __CALLING_FUNCTION__, "spongeCoordinates" );
  ZFSFloatPointerScratchSpace spongeCoords(noSpongeInfo,__CALLING_FUNCTION__, "spongeCoordPointer");

  ZFSInt totMemSize=0;
  for(ZFSId i=0; i<noSpongeInfo; ++i){
    ZFSInt size=1;
    for(ZFSId dim=0; dim<nDim; ++dim){
      size*=(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]+1);
    }
    spongeCoords[i]=&coordMem[totMemSize];
    totMemSize+=nDim*size;
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
    for(ZFSInt dim=0; dim<nDim; ++dim){
      if(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]==0) continue;
      size*=(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]);
    }
    cellmemSize+=size;
  }
  //2.2) allocate the space for all the coordinates in Scratch!
  //memory will not be needed later!
  ZFSFloatScratchSpace coordCellMem(nDim*cellmemSize, __CALLING_FUNCTION__, "spongeCellCoordinates" );
  ZFSFloatPointerScratchSpace spongeSurfCoords(noSpongeInfo,__CALLING_FUNCTION__, "spongeCellCoordPointer");

  ZFSInt totCellMemSize=0;
  for(ZFSId i=0; i<noSpongeInfo; ++i){
    ZFSInt size=1;
    for(ZFSId dim=0; dim<nDim; ++dim){
      if(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]==0) continue;
      size*=(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]);
    }
    spongeSurfCoords[i]=&coordCellMem[totCellMemSize];
    totCellMemSize+=nDim*size;
  }




  //3) read in the coordinates of the grid points
  //open file for reading the grid data
  ZFSString gridFileName =*(ZFSContext::getProperty("gridInputFileName", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString(0));

  //file id to access the file
  ZFSId file_id = -1;
  //unique identifier needed to associate grid and solution in case of restart
  //const char* aUID= new char[18];


  //open file and read number of blocks and uid
  file_id = io_openfile("", gridFileName.c_str(), "collective", m_zfsStrctrdComm);

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
      ZFSId offset[2]={0,0};
      ZFSId size[2]={0,0};
      memSize = 1;
      for(ZFSId dim=1; dim>=0; --dim){
        size[dim]=(spongeInfo[i]->end1[1-dim]-spongeInfo[i]->start1[1-dim]+1);
        memSize*=size[dim];
        offset[dim]= spongeInfo[i]->start1[1-dim];
      }
      //read in the data if  processor zero else read nothing!!!
      //determine the Block name
      ZFSString bName="block";
      stringstream number;
      number << spongeInfo[i]->Id1;
      bName += number.str();

      io_read_ddataset_part1d1(file_id, bName.c_str(), "x", nDim, offset, size, &spongeCoords[i][0]);
      io_read_ddataset_part1d1(file_id, bName.c_str(), "y", nDim, offset, size, &spongeCoords[i][memSize]);
    }
  } else {
    for(ZFSId i=0; i<noSpongeInfo; ++i){
      ZFSId offset[2]={0,0};
      ZFSId size[2]={0,0};
      ZFSString bName="block";
      stringstream number;
      number << spongeInfo[i]->Id1;
      bName += number.str();
      io_read_ddataset_part1d1(file_id, bName.c_str(), "x", nDim, offset, size, NULL);
      io_read_ddataset_part1d1(file_id, bName.c_str(), "y", nDim, offset, size, NULL);
    }
  }

  io_closefile(file_id);

  //now broadcast the information to everyone!!!
  MPI_Bcast(&spongeCoords[0][0], totMemSize, MPI_DOUBLE, 0, m_zfsStrctrdComm);

  //4) computing the coordinates of surface center from corner points;

  for(ZFSId ii=0; ii<noSpongeInfo; ++ii){
    ZFSInt label,size1,count=0;
    for(label=0; label<nDim; label++){ //2== dimensions
      if(spongeInfo[ii]->end1[label]-spongeInfo[ii]->start1[label]==0) break;
    }
    switch(label){
    case 0:{
      size1=spongeInfo[ii]->end1[1]-spongeInfo[ii]->start1[1]+1;
      for(ZFSInt i=0;i<size1-1;i++){
        ZFSInt I=i;
        ZFSInt IP=i+1;
        spongeSurfCoords[ii][count]=0.5*(spongeCoords[ii][I]+spongeCoords[ii][IP]);
        spongeSurfCoords[ii][count+(size1-1)]=0.5*(spongeCoords[ii][I+size1]+spongeCoords[ii][IP+size1]);
        count++;
      }
      break;
    }
    case 1:{
      size1=spongeInfo[ii]->end1[0]-spongeInfo[ii]->start1[0]+1;
      for(ZFSInt i=0;i<size1-1;i++){
        ZFSInt I=i;
        ZFSInt IP=i+1;
        spongeSurfCoords[ii][count]=0.5*(spongeCoords[ii][I]+spongeCoords[ii][IP]);
        spongeSurfCoords[ii][count+(size1-1)]=0.5*(spongeCoords[ii][I+size1]+spongeCoords[ii][IP+size1]);
        count++;
      }
      break;
    }
    default: zfsTerm(1,__CALLING_FUNCTION__, "sponge direction is messed up");
    }
  }

  //now everyone can determine the shortest distance

  zfs_log << "sponge layer build: searching for the nearest points (building simgma sponge)" << endl;
  //cout << "seraching for the nearest points(building sigma sponge)" << endl;
  //build a k-d-tree for a quick search:
  //1) rearragne the coordinates into points;
  for(ZFSId i=0; i<noSpongeInfo; ++i){
    ZFSInt noPoints=1;
    for(ZFSId dim=0; dim<nDim; ++dim){
      if(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]==0)continue;
      noPoints*=(spongeInfo[i]->end1[dim]-spongeInfo[i]->start1[dim]);
    }
    //build up the points (surface centres)
    vector < Point<2> > pts;
    for(ZFSId j=0; j<noPoints; ++j){
      Point<2> a(spongeSurfCoords[i][j],spongeSurfCoords[i][j+noPoints]);
      pts.push_back(a);
    }

    //build up the tree
    KDtree<2> tree(pts);
    ZFSFloat distance = -1.0;
    //go through all the cells an determine the closest distance
    ZFSFloat spongeThickness= spongeInfo[i]->spongeThickness;
    for(ZFSInt id=0; id<m_noStrctrdCells; ++id){
      distance=-1.1111111111111111; //to check
      Point<2> pt(m_cells->coordinates[0][id],m_cells->coordinates[1][id]);
      (void) tree.nearest(pt, distance);
      if(distance<=spongeThickness){
        ZFSFloat spongeFactor = spongeInfo[i]->sigma*pow((spongeThickness - distance)/spongeThickness, spongeInfo[i]->beta);
        m_cells->fq[FQ->SPONGE_FACTOR][id] = zfsMAX(m_cells->fq[FQ->SPONGE_FACTOR][id], spongeFactor);
      }
    }
  }

  zfs_log << "Spone layer SUCESSFUL: finished building sigma sponge " << endl;
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::updateSpongeLayer(){
  // damp to infinity values
  // for the moment only rho and rhoE are damped
  const ZFSFloat gammaMinusOne = m_block->m_gamma-1.0;
  if(m_spongeLayerType == 1){
    ZFSFloat deltaRhoE =F0, deltaRho=F0;
    for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++){
      for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers;i++){
        // compute the forcing term
        // for the pressure or engery and the velocity
        ZFSId cellId = cellIndex(i,j);
        const ZFSFloat rhoE = m_cells->pvariables[PV->P][cellId]/gammaMinusOne+F1B2*m_cells->pvariables[PV->RHO][cellId]*(POW2(m_cells->pvariables[PV->U][cellId])+
                                                                                                                          POW2(m_cells->pvariables[PV->V][cellId]));

        deltaRhoE= rhoE-CV->rhoEInfinity;
        deltaRho = m_cells->pvariables[ PV->RHO ][ cellId ] - (CV->rhoInfinity * m_targetDensityFactor);

        m_cells->rightHandSide[ CV->RHO_E ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRhoE*m_cells->cellJac[cellId];//deltaP * m_cells->cellJac[cellId];
        m_cells->rightHandSide[ CV->RHO ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
      }
    }
  }else{
    // damp to values in the FQ field (set at startup, from RANS etc.)
    ZFSFloat deltaRhoE =F0, deltaRho=F0;
    for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++){
      for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers;i++){
        ZFSId cellId = cellIndex(i,j);
        const ZFSFloat rhoE = m_cells->pvariables[PV->P][cellId]/gammaMinusOne+F1B2*m_cells->pvariables[PV->RHO][cellId]*(POW2(m_cells->pvariables[PV->U][cellId])+
                                                                                                                          POW2(m_cells->pvariables[PV->V][cellId]));

        deltaRhoE= rhoE-m_cells->fq[FQ->SPONGE_RHO_E][cellId];
        deltaRho = m_cells->pvariables[ PV->RHO ][ cellId ] - (m_cells->fq[FQ->SPONGE_RHO][cellId] * m_targetDensityFactor);
        m_cells->rightHandSide[ CV->RHO_E ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRhoE*m_cells->cellJac[cellId];
        m_cells->rightHandSide[ CV->RHO ][ cellId ] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
      }
    }
  }
}


//>marian: new function to compute distance to nearest wall for each cell
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::computeWallDistances(){

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
    for(ZFSInt dim=0; dim<2; ++dim){
      size*=(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]+1);
    }
    memSize+=size;
  }

  //1.2) allocate the memory
  ZFSFloatScratchSpace coordMem(2*memSize, __CALLING_FUNCTION__, "wallCoordinates" );
  ZFSFloatPointerScratchSpace wallDistCoords(noWallDistInfo,__CALLING_FUNCTION__, "wallCoordsPointer");

  ZFSInt totMemSize=0;
  for(ZFSId i=0; i<noWallDistInfo; ++i){
    ZFSInt size=1;
    for(ZFSId dim=0; dim<2; ++dim){
      size*=(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]+1);
    }
    wallDistCoords[i]=&coordMem[totMemSize];
    totMemSize+=2*size;
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
    for(ZFSInt dim=0; dim<2; ++dim){
      if(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]==0) continue;
      size*=(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]);
    }
    cellmemSize+=size;
  }
  //2.2) allocate the space for all the coordinates in Scratch!
  //memory will not be needed later!
  ZFSFloatScratchSpace coordCellMem(2*cellmemSize, __CALLING_FUNCTION__, "wallDistCellCoordinates" );
  ZFSFloatPointerScratchSpace wallDistSurfCoords(noWallDistInfo,__CALLING_FUNCTION__, "wallDistCellCoordPointer");

  ZFSInt totCellMemSize=0;
  for(ZFSId i=0; i<noWallDistInfo; ++i){
    ZFSInt size=1;
    for(ZFSId dim=0; dim<2; ++dim){
      if(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]==0) continue;
      size*=(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]);
    }
    wallDistSurfCoords[i]=&coordCellMem[totCellMemSize];
    totCellMemSize+=2*size;
  }

  //3) read in the coordinates of the grid points
  //open file for reading the grid data
  ZFSString gridFileName =*(ZFSContext::getProperty("gridInputFileName", 0, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString(0));
  ZFSId file_id = -1;
  //open file and read number of blocks and uid
  file_id = io_openfile("", gridFileName.c_str(), "collective", m_block->m_zfsStrctrdComm);

  //read the data in and distribute ist

  //the split of reading and distributing is done on purpose!
  //this is because we then know what the library is doing
  //and not what the io_library such as hdf or netcdf should do
  //but does not
  //a good library should handle that directly! TEST IT
  memSize =1;
  //read in the data if  processor zero else read nothing!!!
  if(m_block->domainId()==0){
    for(ZFSId i=0; i<noWallDistInfo; ++i){
      ZFSId offset[2]={0,0};
      ZFSId size[2]={0,0};
      memSize = 1;
      for(ZFSId dim=1; dim>=0; --dim){
        size[dim]=(wallDistInfo[i]->end1[1-dim]-wallDistInfo[i]->start1[1-dim]+1);
        memSize*=size[dim];
        offset[dim]= wallDistInfo[i]->start1[1-dim];
      }
      //read in the data if  processor zero else read nothing!!!
      //determine the Block name
      ZFSString bName="block";
      stringstream number;
      number << wallDistInfo[i]->Id1;
      bName += number.str();

      io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 2, offset, size, &wallDistCoords[i][0]);
      io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 2, offset, size, &wallDistCoords[i][memSize]);
    }
  } else {
    for(ZFSId i=0; i<noWallDistInfo; ++i){
      ZFSId offset[2]={0,0};
      ZFSId size[2]={0,0};
      ZFSString bName="block";
      stringstream number;
      number << wallDistInfo[i]->Id1;
      bName += number.str();
      io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 2, offset, size, NULL);
      io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 2, offset, size, NULL);
    }
  }

  io_closefile(file_id);

  //cout << "broadcasting wall-bc information" << endl;
  //now broadcast the information to everyone!!!
  MPI_Bcast(&wallDistCoords[0][0], totMemSize, MPI_DOUBLE, 0, m_block->m_zfsStrctrdComm);
  //cout << "broadcast end" << endl;

  //4) computing the coordinates of surface center from corner points;

  for(ZFSId ii=0; ii<noWallDistInfo; ++ii){
    ZFSInt label,size1,count=0;
    for(label=0; label<2; label++){ //2== dimensions
      if(wallDistInfo[ii]->end1[label]-wallDistInfo[ii]->start1[label]==0) break;
    }
    switch(label){
    case 0:{
      size1=wallDistInfo[ii]->end1[1]-wallDistInfo[ii]->start1[1]+1;
      for(ZFSInt i=0;i<size1-1;i++){
        ZFSInt I=i;
        ZFSInt IP=i+1;
        wallDistSurfCoords[ii][count]=0.5*(wallDistCoords[ii][I]+wallDistCoords[ii][IP]);
        wallDistSurfCoords[ii][count+(size1-1)]=0.5*(wallDistCoords[ii][I+size1]+wallDistCoords[ii][IP+size1]);
        count++;
      }
      break;
    }
    case 1:{
      size1=wallDistInfo[ii]->end1[0]-wallDistInfo[ii]->start1[0]+1;
      for(ZFSInt i=0;i<size1-1;i++){
        ZFSInt I=i;
        ZFSInt IP=i+1;
        wallDistSurfCoords[ii][count]=0.5*(wallDistCoords[ii][I]+wallDistCoords[ii][IP]);
        wallDistSurfCoords[ii][count+(size1-1)]=0.5*(wallDistCoords[ii][I+size1]+wallDistCoords[ii][IP+size1]);
        count++;
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
    for(ZFSId dim=0; dim<2; ++dim){
      if(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]==0)continue;
      noPoints*=(wallDistInfo[i]->end1[dim]-wallDistInfo[i]->start1[dim]);
    }
    //build up the points (surface centres)
    vector < Point<2> > pts;
    for(ZFSId j=0; j<noPoints; ++j){
      Point<2> a(wallDistSurfCoords[i][j],wallDistSurfCoords[i][j+noPoints] );
      pts.push_back(a);
    }
    //build up the tree
    KDtree<2> tree(pts);
    ZFSFloat distance = -1.0;

    //go through all the cells an determine the closest distance
    for(ZFSInt id=0; id<m_noStrctrdCells; ++id){
      distance=-1.1111111111111111; //to check
      Point<2> pt(m_cells->coordinates[0][id],m_cells->coordinates[1][id]);
      (void) tree.nearest(pt, distance);

      //take minimum because another bc1000 might be further away than the current but would overwrite the actually closer distance
      m_cells->fq[FQ->WALLDISTANCE][id] = zfsMIN(m_cells->fq[FQ->WALLDISTANCE][id], distance);
    }
  }

  zfs_log << "Wall Distance Computation SUCESSFUL: Saved minimum distance to next wall for all cells " << endl;
}
//<marian


//function to correct the index values in the map for the different boundary conditions
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::correctBndryCndIndices(){

  TRACE();

  //in correcting cell Information
  for(ZFSId bcId=0; bcId < m_noBndryCndIds; bcId++)
    {
      (this->*initBndryCndHandler[bcId])(bcId);
    }
}


template <ZFSBool isRans>
inline ZFSId ZFSStrctrdBndryCnd2D<isRans>::cellIndex(ZFSInt i, ZFSInt j)
{
  return i + (j * m_nCells[1]);
}

template <ZFSBool isRans>
inline ZFSId ZFSStrctrdBndryCnd2D<isRans>::getPointIdFromCell( ZFSInt i, ZFSInt j )
{
  return i + (j * (m_nCells[1] + 1));
}

template <ZFSBool isRans>
inline ZFSId ZFSStrctrdBndryCnd2D<isRans>::getPointIdFromPoint( ZFSId origin, ZFSInt incI, ZFSInt incJ )
{
  return origin + incI + incJ * m_nPoints[1];
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::initBc1000(ZFSId bcId)
{
  (void) bcId;
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::initBc1003(ZFSId bcId)
{
  (void) bcId;
  m_isothermalWallTemperature = F1;
  if(ZFSContext::propertyExists("isothermalWallTemperature", m_blockId )){
    m_isothermalWallTemperature = *(ZFSContext::getProperty("isothermalWallTemperature", m_blockId, __CALLING_FUNCTION__,&m_isothermalWallTemperature)->asFloat(0));

  }
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::initBc2002(ZFSId bcId)
{
  (void) bcId;
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::initBc2004(ZFSId bcId)
{
  (void) bcId;
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::initBc2005(ZFSId bcId)
{
  (void) bcId;
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::initBc2300(ZFSId bcId) {
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  cout << "In initBC2300" << endl;

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

  ZFSFloatScratchSpace bl(5,N+1, __CALLING_FUNCTION__, "bl");

  //position at which ReTheta is defined and ReX is known:
  const ZFSFloat rex = 100000;
  const ZFSFloat rexPos = 470.59;

  //transform to position of the boundary condition
  const ZFSFloat myPos = 434.0;
  const ZFSFloat myRe = rex/rexPos*myPos;
  const ZFSFloat lengthScaling = myPos/sqrt(myRe);

  for(ZFSId n=0;n<N+1; n++) {
    const ZFSFloat rho = CV->rhoInfinity;
    bl(PV->U,n) = PV->UInfinity*result(n,1);
    bl(PV->V,n) = PV->UInfinity*result(n,0)/sqrt(rex);
    cout << "n: " << n << " bl(rho_u): " << bl(PV->U,n) << endl;
    bl(PV->P,n) = PV->PInfinity;
    bl(PV->RHO,n) = rho;
    bl(4,n) = lengthScaling*t(n);
  }

  //interpolation part
  cout.precision(10);
  FILE* f_blasius;
  f_blasius = fopen("./blasius.dat", "a+");
  ZFSFloat theta = F0;
  for(ZFSId j = start[1]; j<end[1]; j++) {
    ZFSFloat yCoord = m_cells->coordinates[1][cellIndex(start[0],j)];
    ZFSFloatScratchSpace variables(PV->noVariables, __CALLING_FUNCTION__, "variables");
    ZFSBool hasIntStencil = false;
    for(ZFSId n=0;n<N; n++) {
      if(bl(4,n)<yCoord&&bl(4,n+1)>=yCoord) {
        for(ZFSId var=0; var<PV->noVariables; var++) {
          variables(var) = bl(var,n) + (yCoord-bl(4,n))/(bl(4,n+1)-bl(4,n))*(bl(var,n+1)-bl(var,n));
          if(var==0)
            cout << "yCoord: " << yCoord << " var: " << var << " ybl(n): " << bl(4,n) << " ybl(n+1): " << bl(4,n+1) << " var bl(n): " << bl(var,n) << " var bl(n+1): " << bl(var,n+1) << " interpolated var: " << variables(var) << endl;
        }
        hasIntStencil = true;
        break;
      }
    }

    if(!hasIntStencil) {
      if(yCoord<=bl(4,0)) {
        for(ZFSId var=0; var<PV->noVariables; var++) {
          variables(var) = bl(var,0);
        }
      } else {
        for(ZFSId var=0; var<PV->noVariables; var++) {
          variables(var) = bl(var,N);
        }
      }
    }

    ZFSId p1 = getPointIdFromCell(start[0],j);
    ZFSId p2 = getPointIdFromPoint(p1,0,1);
    ZFSFloat yDelta = m_coordinates[1][p2]-m_coordinates[1][p1];
    theta += variables(PV->U)*variables(PV->RHO)/CV->rhoInfinity/PV->UInfinity*(F1 - variables(PV->U)/PV->UInfinity)*yDelta;

    fprintf(f_blasius, "%f", yCoord);
    fprintf(f_blasius, " %f", variables(0));
    fprintf(f_blasius, " %f", variables(1));
    fprintf(f_blasius, " %f", variables(2));
    fprintf(f_blasius, " %f", variables(3));
    fprintf(f_blasius, "\n");


    for(ZFSId i = start[0]; i<end[0]; i++) {
      const ZFSId cellId = cellIndex(i,j);
      for(ZFSId var=0; var<PV->noVariables; var++) {
        m_cells->pvariables[var][cellId] = variables(var);
      }
    }
  }
  fclose(f_blasius);
  cout << "Theta: " << theta << endl;
}

/** Rescaling inflow
 *
 *  put values from field to gc to avoid nan, only matters
 *  for first Runge-Kutta Step at first timeStep
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::initBc2510(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  for(ZFSId var=0; var < PV->noVariables; var++) {
    for(ZFSId i=start[0]; i<end[0]; i++) {
      for(ZFSId j=start[1]; j<end[1]; j++) {
        ZFSId cellId = cellIndex(i,j);
        ZFSId cellIdAdj = cellIndex(m_noGhostLayers, j);
        m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdAdj];
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
void ZFSStrctrdBndryCnd2D<isRans>::initBc2600(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  m_block->m_bc2600 = true;

  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
      {
        if(m_block->m_bc2600InitialStartup) {
          //First copy values from the field into the ghostcells
          for(ZFSInt i = start[0]; i<end[0]; i++) {
            for(ZFSInt j = start[1]; j<end[1]; j++) {
              const ZFSId cellId = cellIndex(i,j);
              const ZFSId cellIdAdj = cellIndex(m_noGhostLayers,j);
              for(ZFSId var=0; var<PV->noVariables; var++) {
                m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdAdj];
              }
            }
          }
        }

        //Then copy the values from the ghostcells into restart field
        for(ZFSInt i = start[0]; i<end[0]; i++) {
          for(ZFSInt j = m_noGhostLayers; j<end[1]-m_noGhostLayers; j++) {
            const ZFSId cellId = cellIndex(i,j);
            const ZFSId cellIdBc = i + (j - m_noGhostLayers)*m_noGhostLayers;
            for(ZFSId var=0; var<PV->noVariables; var++) {
              m_block->m_bc2600Variables[var][cellIdBc] = m_cells->pvariables[var][cellId];
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
void ZFSStrctrdBndryCnd2D<isRans>::initBc2021(ZFSId bcId)
{
  (void) bcId;
  m_bc2021Gradient = 1.0;
  m_bc2021Gradient = *(ZFSContext::getProperty("bc2021Gradient", m_blockId, __CALLING_FUNCTION__,&m_bc2021Gradient)->asFloat(0));
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::initBc3000(ZFSId bcId)
{
  (void) bcId;
}

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc1000(ZFSId bcId)
{
  TRACE();

  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
      {
        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;
        const ZFSInt cellShift = 2*m_noGhostLayers-1;

        for(ZFSId j=start[1]; j<end[1] ; j++)
          {
            for(ZFSId i=start[0]; i<end[0]; i++)
              {
                cellId=cellIndex(i,j);
                cellIdadj=cellIndex(cellShift-i,j);
                m_cells->pvariables[PV->RHO][cellId]  =m_cells->pvariables[PV->RHO][cellIdadj];
                m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj];
                m_cells->pvariables[PV->U][cellId]=-m_cells->pvariables[PV->U][cellIdadj];
                m_cells->pvariables[PV->V][cellId]=-m_cells->pvariables[PV->V][cellIdadj];
                if(isRans) {
                  m_cells->pvariables[PV->RANS_VAR[0]][cellId]=-m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
                }
              }
          }

        break;
      }
    case 1:
      {
        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;
        const ZFSInt cellShift = 2*(m_nCells[2]-1)-2*m_noGhostLayers+1;

        for(ZFSId j=start[1]; j<end[1] ; j++)
          {
            for(ZFSId i=start[0]; i<end[0]; i++)
              {
                cellId=cellIndex(i,j);
                cellIdadj=cellIndex(cellShift-i,j);
                m_cells->pvariables[PV->RHO][cellId]  =m_cells->pvariables[PV->RHO][cellIdadj];
                m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj];
                m_cells->pvariables[PV->U][cellId]=-m_cells->pvariables[PV->U][cellIdadj];
                m_cells->pvariables[PV->V][cellId]=-m_cells->pvariables[PV->V][cellIdadj];
                if(isRans) {
                  m_cells->pvariables[PV->RANS_VAR[0]][cellId]=-m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
                }
              }
          }

        break;
      }
    case 2:
      {

        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;
        const ZFSInt cellShift = 2*m_noGhostLayers-1;

        for(ZFSId j=start[1]; j<end[1] ; j++) {
          for(ZFSId i=start[0]; i<end[0]; i++) {
            cellId=cellIndex(i,j); //ghost
            cellIdadj=cellIndex(i,cellShift-j); // field
            m_cells->pvariables[PV->RHO][cellId] = m_cells->pvariables[PV->RHO][cellIdadj];
            m_cells->pvariables[PV->P][cellId]   = m_cells->pvariables[PV->P][cellIdadj];
            m_cells->pvariables[PV->U][cellId]   = -m_cells->pvariables[PV->U][cellIdadj];
            m_cells->pvariables[PV->V][cellId]   = -m_cells->pvariables[PV->V][cellIdadj];
            if(isRans) {
              m_cells->pvariables[PV->RANS_VAR[0]][cellId]=-m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
            }
          }
        }
        break;
      }

    case 3:
      {
        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;
        const ZFSInt cellShift = 2*(m_nCells[1]-1)-2*m_noGhostLayers+1;

        for(ZFSId j=start[1]; j<end[1] ; j++)
          {
            for(ZFSId i=start[0]; i<end[0]; i++)
              {
                cellId=cellIndex(i,j); //ghost
                cellIdadj=cellIndex(i,cellShift-j); //field
                m_cells->pvariables[PV->RHO][cellId]  =m_cells->pvariables[PV->RHO][cellIdadj];
                m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj];
                m_cells->pvariables[PV->U][cellId]=-m_cells->pvariables[PV->U][cellIdadj];
                m_cells->pvariables[PV->V][cellId]=-m_cells->pvariables[PV->V][cellIdadj];
                if(isRans) {
                  m_cells->pvariables[PV->RANS_VAR[0]][cellId]=-m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
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

// euler wall bc
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc1001(ZFSId bcId)
{
  TRACE();

  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face)
  {
  case 2:
  case 3:
  {
    ZFSInt cellShift;
    if (m_physicalBCMap[bcId]->face == 2) {
      cellShift = 2*m_noGhostLayers-1;
    } else if (m_physicalBCMap[bcId]->face == 3) {
      cellShift = 2*(m_nCells[0]-1)-2*m_noGhostLayers+1;
    }

    for(ZFSId j=start[1]; j<end[1] ; j++) {
      for(ZFSId i=start[0]; i<end[0]; i++) {
        ZFSId cellIdA1 = -1, pIJ =-1;
        if(m_physicalBCMap[bcId]->face==2) {
          pIJ = getPointIdFromCell(i,m_noGhostLayers);
          cellIdA1 = cellIndex(i,m_noGhostLayers);
        } else {
          pIJ = getPointIdFromCell(i,m_block->m_nPoints[0]-3);
          cellIdA1 = cellIndex(i,m_nCells[0]-3);
        }
        const ZFSFloat x1 = m_coordinates[0][pIJ];
        const ZFSFloat y1 = m_coordinates[1][pIJ];
        const ZFSFloat x2 = m_coordinates[0][pIJ + 1];
        const ZFSFloat y2 = m_coordinates[1][pIJ + 1];

        const ZFSFloat dydx = (y2 - y1) / (x2 - x1);
        const ZFSFloat alpha = -atan(dydx);

        const ZFSId cellId=cellIndex(i,j); //ghost
        const ZFSId cellIdadj=cellIndex(i,cellShift-j); //field

        const ZFSFloat rho = m_cells->pvariables[PV->RHO][cellIdA1]; //apply rho from first active cell to both ghost-cells
        const ZFSFloat u = m_cells->pvariables[PV->U][cellIdadj];
        const ZFSFloat v = m_cells->pvariables[PV->V][cellIdadj];
                
        const ZFSFloat uPrime =   u * cos(alpha) - v * sin(alpha);
        const ZFSFloat vPrime = -(u * sin(alpha) + v * cos(alpha));
                
        const ZFSFloat uGC =  uPrime * cos(-alpha) - vPrime * sin(-alpha);
        const ZFSFloat vGC =  uPrime * sin(-alpha) + vPrime * cos(-alpha);

        m_cells->pvariables[PV->RHO][cellId]  = rho;
        m_cells->pvariables[PV->U][cellId] = uGC;
        m_cells->pvariables[PV->V][cellId] = vGC;

        //apply pressure from first active cell to all cells
        m_cells->pvariables[PV->P][cellId]= m_cells->pvariables[PV->P][cellIdA1];
                
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellId]=-m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
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

//isothermal Wall
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc1003(ZFSId bcId){
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  const ZFSFloat temp = m_isothermalWallTemperature*PV->TInfinity;

  switch(m_physicalBCMap[bcId]->face){
  case 2:{
    const ZFSInt cellShift = 2*m_noGhostLayers-1;
    for(ZFSId j=start[1]; j<end[1] ; j++){
      for(ZFSId i=start[0]; i<end[0]; i++){
        const ZFSId cellId=cellIndex(i,j); //ghost
        const ZFSId cellIdAdj=cellIndex(i,cellShift-j); // field
        const ZFSFloat rhoActive = m_cells->pvariables[PV->RHO][cellIdAdj];
        const ZFSFloat u = (-1.0)*(m_cells->pvariables[PV->U][cellIdAdj]);
        const ZFSFloat v = (-1.0)*(m_cells->pvariables[PV->V][cellIdAdj]);
        const ZFSFloat pressure1=m_cells->pvariables[PV->P][cellIdAdj];
        const ZFSFloat rhoWall = pressure1*m_block->m_gamma/temp;
        const ZFSFloat rho = F2*rhoWall-rhoActive;

        m_cells->pvariables[PV->RHO][cellId]=rho;
        m_cells->pvariables[PV->U][cellId]=u;
        m_cells->pvariables[PV->V][cellId]=v;
        m_cells->pvariables[PV->P][cellId]=pressure1;

        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellId]=-m_cells->pvariables[PV->RANS_VAR[0]][cellIdAdj];
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

/** Subsonic Inflow <== tfs2001
 *
 *  rho=rho_inf, u=u_inf, v=v_inf, w=w_inf, dp/dn=0
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2001(ZFSId bcId)
{
  TRACE();

  //implemented for i-direction only for the moment
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
      {

        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;
            for(ZFSId j=start[1]; j<end[1] ; j++)
              {
                for(ZFSId i=start[0]; i<end[0]; i++)
                  {
                    cellId=cellIndex(m_noGhostLayers-1-i,j);
                    cellIdadj=cellIndex(m_noGhostLayers-i,j);

                    m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
                    m_cells->pvariables[PV->U][cellId]=PV->UInfinity;
                    m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
                    m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj];

                    if(isRans) {
                      m_cells->pvariables[PV->RANS_VAR[0]][cellId]=PV->ransInfinity[0];
                    }
                  }
              }
        break;
      }
    case 2:
      {
        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;
            for(ZFSId j=start[1]; j<end[1] ; j++)
              {
                for(ZFSId i=start[0]; i<end[0]; i++)
                  {
                    cellId=cellIndex(i,m_noGhostLayers-j-1);//ghost
                    cellIdadj=cellIndex(i,m_noGhostLayers-j);//field

                    m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
                    m_cells->pvariables[PV->U][cellId]=PV->UInfinity;
                    m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
                    m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj];

                    if(isRans) {
                      m_cells->pvariables[PV->RANS_VAR[0]][cellId]=PV->ransInfinity[0];
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


/** Subsonic Outflow, not really non-reflecting for face 0,1,3
 *  Simplified characteristic approach (Whitfield)
 *  Still reflects pressure waves at the outflow
 *  @author Pascal Meysonnat
 *  rho=rho, u=u, v=v, w=w, p=p_inf
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2004(ZFSId bcId)
{
  TRACE();
  /*
  const ZFSInt IJK[nDim] = {1,m_nCells[1]};
  const ZFSFloat gamma= m_block->m_gamma;
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  //Here we find out the normal direction of the
  //boundary and the tangential direction.
  //This way we can make a general formulation of
  //the boundary condition
  const ZFSInt face = m_physicalBCMap[bcId]->face;
  const ZFSInt normalDir = face/2;
  const ZFSInt firstTangentialDir = (normalDir+1)%nDim;
  const ZFSInt normalDirStart = start[normalDir];
  const ZFSInt firstTangentialStart = start[firstTangentialDir];
  const ZFSInt firstTangentialEnd = end[firstTangentialDir];
  const ZFSInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir]};
  //determine indices for direction help
  const ZFSInt n = (face%2)*2-1; //-1,+1
  const ZFSInt g1 = normalDirStart + (ZFSId)(0.5-(0.5*(ZFSFloat)n)); //+1,0
  const ZFSInt g2 = normalDirStart + (ZFSId)(0.5+(0.5*(ZFSFloat)n)); //0,+1
  const ZFSInt a1 = normalDirStart + (ZFSId)(0.5-(1.5*(ZFSFloat)n)); //+2,-1
  const ZFSInt a2 = normalDirStart + (ZFSId)(0.5-(2.5*(ZFSFloat)n)); //+3,-2
  
  for(ZFSId t1=firstTangentialStart; t1<firstTangentialEnd; t1++) {
    const ZFSId cellIdG1 = g1*inc[0] + t1*inc[1];
    const ZFSId cellIdG2 = g2*inc[0] + t1*inc[1];
    const ZFSId cellIdA1 = a1*inc[0] + t1*inc[1];
    const ZFSId cellIdA2 = a2*inc[0] + t1*inc[1];
    const ZFSFloat dxidx = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+0];
    const ZFSFloat dxidy = m_cells->surfaceMetrics[cellIdA1][normalDir*nDim+1];
    //multiply with n, so it will be -1 or +1 depending if we enter
    //or leave the domain of integration in positive direction
    const ZFSFloat gradxi = n*F1 / sqrt(dxidx*dxidx + dxidy*dxidy);
    const ZFSFloat dxHelp = dxidx*gradxi;
    const ZFSFloat dyHelp = dxidy*gradxi;
    //speed of sound
    const ZFSFloat cBC = sqrt(gamma*m_cells->pvariables[PV->P][cellIdG1]/m_cells->pvariables[PV->RHO][cellIdG1]);
    const ZFSFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];
    const ZFSFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
    const ZFSFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
    const ZFSFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
    const ZFSFloat pInner = m_cells->pvariables[PV->P][cellIdA1];
    const ZFSFloat maContravariant = (dxidx*uInner + dxidy*vInner - m_cells->dxt[normalDir][cellIdA1])*gradxi;
    if(maContravariant < F0) {
      //inflow
      const ZFSFloat p =F1B2*(pInner+PV->PInfinity+rhoBC*cBC*(dxHelp*(uInner - PV->UInfinity)+
                                                                dyHelp*(vInner - PV->VInfinity)));

      const ZFSFloat rho = CV->rhoInfinity + (p-PV->PInfinity)/POW2(cBC);
      const ZFSFloat help = (p-PV->PInfinity)/(rhoBC*cBC);

      m_cells->pvariables[PV->RHO][cellIdG1] = rho;
      m_cells->pvariables[PV->U][cellIdG1] = (PV->UInfinity + help*dxHelp);
      m_cells->pvariables[PV->V][cellIdG1] = (PV->VInfinity + help*dyHelp);
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
  */

  
  TRACE();
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face)
    {
    case 1:
      {
        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;
        ZFSId pIJ=0, pIJP=0;
        ZFSFloat xBC=F0;


        //fully new bc
        ZFSId i= start[0];
        ZFSFloat pBC=F0, rho=F0, u=F0, v=F0;
        ZFSFloat drho=F0, du=F0, dv=F0, dp=F0;
        pBC=PV->PInfinity;
        ZFSFloat pInner=F0, c02=F0, distance=F0;

        for(ZFSId j=start[1]; j<end[1] ; j++)
          {
            cellIdadj=cellIndex(i-1,j);

            //to determine the face coordinates!!!!!!
            pIJ=getPointIdFromCell(i,j);
            pIJP=getPointIdFromPoint(pIJ,0,1);

            xBC=F1B2*(m_coordinates[0][pIJ]+m_coordinates[0][pIJP]);
            pInner=m_cells->pvariables[PV->P][cellIdadj];
            c02=sqrt(m_block->m_gamma*pInner/m_cells->pvariables[PV->RHO][cellIdadj]);

            //first find out the values at the boundary
            //rho
            rho=m_cells->pvariables[PV->RHO][cellIdadj]+((pBC-pInner)/(c02*c02));
            //velocities u, v, w

            ZFSFloat dxidx = m_cells->surfaceMetrics[cellIdadj][0];
            ZFSFloat dxidy = m_cells->surfaceMetrics[cellIdadj][1];

            ZFSFloat gradxi = F1 / sqrt(dxidx*dxidx + dxidy*dxidy);

            ZFSFloat dxHelp = dxidx*gradxi;
            ZFSFloat dyHelp = dxidy*gradxi;

            u=m_cells->pvariables[PV->U][cellIdadj]+dxHelp*((pInner-pBC)/(m_cells->pvariables[PV->RHO][cellIdadj]*c02));
            v=m_cells->pvariables[PV->V][cellIdadj]+dyHelp*((pInner-pBC)/(m_cells->pvariables[PV->RHO][cellIdadj]*c02));

            //extrapolate the variables into the ghost cells
            //gradients
            distance=(xBC-m_cells->coordinates[0][cellIdadj]);

            drho=(rho-m_cells->pvariables[PV->RHO][cellIdadj])/distance;
            du=(u-m_cells->pvariables[PV->U][cellIdadj])/distance;
            dv=(v-m_cells->pvariables[PV->V][cellIdadj])/distance;
            dp=(pBC-m_cells->pvariables[PV->P][cellIdadj])/distance;

            //extrapolate:
            for(ZFSId ii=start[0]; ii<end[0]; ++ii)
              {
                cellId=cellIndex(ii,j);
                distance=(m_cells->coordinates[0][cellId]-m_cells->coordinates[0][cellIdadj]);
                m_cells->pvariables[PV->RHO][cellId]=m_cells->pvariables[PV->RHO][cellIdadj]+drho*distance;
                m_cells->pvariables[PV->U][cellId]=m_cells->pvariables[PV->U][cellIdadj]+du*distance;
                m_cells->pvariables[PV->V][cellId]=m_cells->pvariables[PV->V][cellIdadj]+dv*distance;
                m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj]+dp*distance;
                if(isRans) {
                  m_cells->pvariables[PV->RANS_VAR[0]][cellId]=m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
                }
              }
          }

        break;
      }

    case 3:
      {

        //fully new BC
        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;

        ZFSId j= start[1];
        ZFSFloat pBC=F0, rho=F0, u=F0, v=F0;
        ZFSFloat drho=F0, du=F0, dv=F0, dp=F0;
        pBC=PV->PInfinity;
        ZFSFloat pInner=F0, c02=F0, distance=F0;
        ZFSFloat yBC=F0;
        ZFSId pIJ=0, pIPJ=0;

        for(ZFSId i=start[0]; i<end[0] ; ++i)
          {
            cellIdadj=cellIndex(i,j-1);
            pIJ =getPointIdFromCell(i,j);
            pIPJ=getPointIdFromPoint(pIJ,1,0);

            yBC=F1B2*(m_coordinates[1][pIJ]+m_coordinates[1][pIPJ]);

            pInner=m_cells->pvariables[PV->P][cellIdadj];
            c02=sqrt(m_block->m_gamma*pInner/m_cells->pvariables[PV->RHO][cellIdadj]);

            //first find out the values at the boundary
            //rho
            rho=m_cells->pvariables[PV->RHO][cellIdadj]+((pBC-pInner)/(c02*c02));

            //velocities u, v, w
            ZFSFloat dxidx = m_cells->surfaceMetrics[cellIdadj][2];
            ZFSFloat dxidy = m_cells->surfaceMetrics[cellIdadj][3];

            ZFSFloat gradxi = F1 / sqrt(dxidx*dxidx + dxidy*dxidy);

            ZFSFloat dxHelp = dxidx*gradxi;
            ZFSFloat dyHelp = dxidy*gradxi;

            u=m_cells->pvariables[PV->U][cellIdadj]+dxHelp*((pInner-pBC)/(m_cells->pvariables[PV->RHO][cellIdadj]*c02));
            v=m_cells->pvariables[PV->V][cellIdadj]+dyHelp*((pInner-pBC)/(m_cells->pvariables[PV->RHO][cellIdadj]*c02));

            //extrapolate the variables into the ghost cells
            //gradients
            distance=(yBC-m_cells->coordinates[1][cellIdadj]);

            drho=(rho-m_cells->pvariables[PV->RHO][cellIdadj])/distance;
            du=(u-m_cells->pvariables[PV->U][cellIdadj])/distance;
            dv=(v-m_cells->pvariables[PV->V][cellIdadj])/distance;
            dp=(pBC-m_cells->pvariables[PV->P][cellIdadj])/distance;

            //extrapolate:
            for(ZFSId jj=start[1]; jj<end[1]; ++jj)
              {
                cellId=cellIndex(i,jj);
                distance=(m_cells->coordinates[1][cellId]-m_cells->coordinates[1][cellIdadj]);
                m_cells->pvariables[PV->RHO][cellId]=m_cells->pvariables[PV->RHO][cellIdadj]+drho*distance;
                m_cells->pvariables[PV->U][cellId]=m_cells->pvariables[PV->U][cellIdadj]+du*distance;
                m_cells->pvariables[PV->V][cellId]=m_cells->pvariables[PV->V][cellIdadj]+dv*distance;
                m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj]+dp*distance;
                if(isRans) {
                  m_cells->pvariables[PV->RANS_VAR[0]][cellId]=m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
                }
              }

          }
        break;
      }
    default:
      {
        zfsTerm(1, __CALLING_FUNCTION__, "Face direction not implemented");
      }
    }
  
}

/** Supersonic Inflow
 *
 *  rho=rho, u=u, v=v, w=w, p=p
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2002(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face){
  case 1:{
    for(ZFSId j=start[1]; j<end[1] ; j++){
      for(ZFSId i=start[0]; i<end[0]; i++){
        ZFSId cellId=cellIndex(i,j);
        m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
        m_cells->pvariables[PV->U][cellId]=PV->UInfinity;
        m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
        m_cells->pvariables[PV->P][cellId]=PV->PInfinity;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellId]=PV->ransInfinity[0];
        }
      }
    }
    break;
  }
  default:
    {
      //do nothing
    }
  }
}

/** Supersonic Outflow,
 *
 *  rho=rho, u=u, v=v, w=w, p=p
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2005(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face){
  case 1:{
     for(ZFSId j=start[1]; j<end[1] ; j++){
       for(ZFSId i=start[0]; i<end[0]; i++){
         ZFSId cellId=cellIndex(i,j);
         ZFSId cellIdadj=cellIndex(i-1,j);

         for(ZFSId var=0; var<PV->noVariables; var++){
           m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdadj];
         }
       }
     }
    break;
  }
  case 3:{
     for(ZFSId j=start[1]; j<end[1] ; j++){
       for(ZFSId i=start[0]; i<end[0]; i++){
         ZFSId cellId=cellIndex(i,j);
         ZFSId cellIdadj=cellIndex(i,j-1);
         for(ZFSId var=0; var<PV->noVariables; var++){
           m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdadj];
         }
       }
     }
     break;
  }
  default: {
    zfsTerm(1, __CALLING_FUNCTION__, "Face direction not implemented");
  }
  }
}


/** Subsonic Outflow
 * extrapolate all but pressure, prescribe p8
 *
 *  rho=rho, u=u, v=v, w=w, p=p_inf
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2007(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face)
  {
  case 1:
  {
    for(ZFSId j=start[1]; j<end[1] ; j++) {
      ZFSId cellId = cellIndex(start[0],j);
      ZFSId cellIdadj=cellIndex(start[0]-1,j);

      m_cells->pvariables[PV->RHO][cellId]  =m_cells->pvariables[PV->RHO][cellIdadj];
      m_cells->pvariables[PV->U][cellId]=m_cells->pvariables[PV->U][cellIdadj];
      m_cells->pvariables[PV->V][cellId]=m_cells->pvariables[PV->V][cellIdadj];
      m_cells->pvariables[PV->P][cellId]=PV->PInfinity;

      if(isRans) {
        m_cells->pvariables[PV->RANS_VAR[0]][cellId]=m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
      }

      for(ZFSId var=0; var<PV->noVariables; var++) {
        ZFSId cellIdP1 = cellIndex(start[0]+1,j);
        m_cells->pvariables[var][cellIdP1] = 2.0*m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
      }
    }
    break;
  }
  case 3:
  {
    for(ZFSId i=start[0]; i<end[0] ; i++) {
      ZFSId cellId = cellIndex(i,start[1]);
      ZFSId cellIdadj=cellIndex(i,start[1]-1);

      m_cells->pvariables[PV->RHO][cellId]  =m_cells->pvariables[PV->RHO][cellIdadj];
      m_cells->pvariables[PV->U][cellId]=m_cells->pvariables[PV->U][cellIdadj];
      m_cells->pvariables[PV->V][cellId]=m_cells->pvariables[PV->V][cellIdadj];
      m_cells->pvariables[PV->P][cellId]=PV->PInfinity;

      if(isRans) {
        m_cells->pvariables[PV->RANS_VAR[0]][cellId]=m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
      }

      for(ZFSId var=0; var<PV->noVariables; var++) {
        ZFSId cellIdP1 = cellIndex(i,start[1]+1);
        m_cells->pvariables[var][cellIdP1] = 2.0*m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
      }
    }

    break;
  }

  default:
  {
    zfsTerm(1, __CALLING_FUNCTION__, "Face direction not implemented");
  }
  }
}

/** BC: Symmetry plane
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc3000(ZFSId bcId)
{
  TRACE();
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
      {
        break;
      }
    case 1:
      {
        break;
      }
    case 2:
      {
        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;
        const ZFSInt cellShift = 2*m_noGhostLayers-1;

        for(ZFSId j=start[1]; j<end[1] ; j++) {
          for(ZFSId i=start[0]; i<end[0]; i++) {
            cellId=cellIndex(i,j); //ghost
            cellIdadj=cellIndex(i,cellShift-j); // field

            m_cells->pvariables[PV->RHO][cellId]=m_cells->pvariables[PV->RHO][cellIdadj];
            m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj];
            m_cells->pvariables[PV->U][cellId]=m_cells->pvariables[PV->U][cellIdadj];
            m_cells->pvariables[PV->V][cellId]=-m_cells->pvariables[PV->V][cellIdadj];
            if(isRans) {
              m_cells->pvariables[PV->RANS_VAR[0]][cellId]=m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
            }
          }
        }

        break;
      }
    case 3:
      {
        ZFSId cellId=-1;
        ZFSId cellIdadj=-1;
        const ZFSInt cellShift = 2*(m_nCells[1]-1)-2*m_noGhostLayers+1;

        for(ZFSId j=start[1]; j<end[1] ; j++) {
          for(ZFSId i=start[0]; i<end[0]; i++) {
            //mirroring
            cellId=cellIndex(i,j); //ghost
            cellIdadj=cellIndex(i,cellShift-j); // field

            m_cells->pvariables[PV->RHO][cellId]=m_cells->pvariables[PV->RHO][cellIdadj];
            m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj];
            m_cells->pvariables[PV->U][cellId]=m_cells->pvariables[PV->U][cellIdadj];
            m_cells->pvariables[PV->V][cellId]=-m_cells->pvariables[PV->V][cellIdadj];
            if(isRans) {
              m_cells->pvariables[PV->RANS_VAR[0]][cellId]=m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
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

// shear flow inflow
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2021(ZFSId bcId)
{
  TRACE();

  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face)
  {
  case 0:
  {
    for(ZFSId j=start[1]; j<end[1]; j++) {
      const ZFSId cellIdG1 = cellIndex(1,j);
      const ZFSId cellIdG2 = cellIndex(0,j);
      const ZFSId cellIdA1 = cellIndex(2,j);
      //const ZFSId cellIdA2 = cellIndex(3,j);

      const ZFSFloat dxidx = m_cells->surfaceMetrics[cellIdA1][0];
      const ZFSFloat dxidy = m_cells->surfaceMetrics[cellIdA1][1];
          
      //multiply with n, so it will be -1 or +1 depending if we enter
      //or leave the domain of integration in positive direction
      const ZFSFloat gradxi = -1*F1 / sqrt(dxidx*dxidx + dxidy*dxidy);
          
      const ZFSFloat dxHelp = dxidx*gradxi;
      const ZFSFloat dyHelp = dxidy*gradxi;
          
      const ZFSFloat cBC = sqrt(m_block->m_gamma*pressure(cellIdG1)/m_cells->pvariables[PV->RHO][cellIdG1]);
      const ZFSFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];

      const ZFSFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
      const ZFSFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
      const ZFSFloat pInner = pressure(cellIdA1);
          
      const ZFSFloat uInflow = PV->UInfinity * (m_cells->coordinates[1][cellIdG1] * m_bc2021Gradient);
      const ZFSFloat vInflow = 0.0;
          
      //inflow
      const ZFSFloat p =F1B2*(pInner+PV->PInfinity+rhoBC*cBC*(dxHelp*(uInner - uInflow)+
                                                              dyHelp*(vInner - vInflow)));

      const ZFSFloat rho = CV->rhoInfinity + (p-PV->PInfinity)/POW2(cBC);
      const ZFSFloat help = (p-PV->PInfinity)/(rhoBC*cBC);

      m_cells->pvariables[PV->RHO][cellIdG1] = rho;
      m_cells->pvariables[PV->U][cellIdG1] = uInflow + help*dxHelp;
      m_cells->pvariables[PV->V][cellIdG1] = vInflow + help*dyHelp;
      m_cells->pvariables[PV->P][cellIdG1] = p;

      if(isRans) {
        m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1]= PV->ransInfinity[0];
      }

      //extrapolate into second ghost cell
      for(ZFSId var=0; var<PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdG2] = F2*m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
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

template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2199(ZFSId bcId)
{
  TRACE();
  //!!!!!!!!!! check for the gamma stuff in tfs
  const ZFSFloat gamma=m_block->m_gamma;
  const ZFSFloat gammaMinusOne=m_block->m_gamma-F1;
  const ZFSFloat fgamma=F1/gamma;
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face)  {
  case 0:{
    ZFSFloat mflux=F0;
    ZFSFloat surface=F0;
    for(ZFSId j=start[1]; j<end[1]; j++){
      const ZFSId cellId= cellIndex(m_noGhostLayers,j);//first inner line
      const ZFSId pointId= getPointIdFromCell(m_noGhostLayers,j);
      const ZFSId pointIdP1=getPointIdFromPoint(pointId,0,1);
      //!!!!!!!!check the sign convention from tfs
      const ZFSFloat dx=m_coordinates[0][pointIdP1]-m_coordinates[0][pointId];
      const ZFSFloat dy=m_coordinates[1][pointIdP1]-m_coordinates[1][pointId];
      const ZFSFloat rhou=m_cells->pvariables[PV->RHO][cellId]*m_cells->pvariables[PV->U][cellId];
      const ZFSFloat rhov=m_cells->pvariables[PV->RHO][cellId]*m_cells->pvariables[PV->V][cellId];
      mflux+=dx*rhou+dy*rhov;
      surface+=sqrt(POW2(dx)+POW2(dy));
    }
    //get cellId at the domain top
    //!!!!!!! check for ii1
    ZFSId cellIdII1=cellIndex(m_noGhostLayers,end[1]-m_noGhostLayers);
    //get averaged flux
    const ZFSFloat rhouAVG=mflux/surface;
    //determine the pressure at the top at the iner line
    ZFSFloat pressure=min(F1, max(F0,m_cells->pvariables[PV->P][cellIdII1]*gamma));
    //fix poin iteration apparently converges fast
    //from schlichting and trockenbrot, page 155
    //ru=sqrt(2.*fgamm1)*pp**fgam*sqrt(f1-pp**(gamm1*fgam))
    for(ZFSId i=0; i<20; i++){
      pressure=pow((F1-pow((rhouAVG/(sqrt(F2*gammaMinusOne)*pow(pressure,fgamma))),F2)),(gammaMinusOne*fgamma));
    }
    pressure*=fgamma;
    
    //scale the velocity profile 

    break;
  }
  default:{ zfsTerm(1, __CALLING_FUNCTION__, "Face direction not implemented)");}
  }
}
   
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2300(ZFSId bcId)
{
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
      {
        for(ZFSId j=start[1]; j<end[1] ; j++) {
          ZFSId cellId = cellIndex(start[0]+1,j);
          ZFSId cellIdadj=cellIndex(start[0]+2,j);

          m_cells->pvariables[PV->P][cellId]=m_cells->pvariables[PV->P][cellIdadj];

          if(isRans) {
            m_cells->pvariables[PV->RANS_VAR[0]][cellId]=m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
          }

          for(ZFSId var=0; var<PV->noVariables; var++) {
            ZFSId cellIdP1 = cellIndex(start[0],j);
            m_cells->pvariables[var][cellIdP1] = 2.0*m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
          }
        }
        break;
      }
    default:
      {
        exit(1);
      }
    }
}

//Inlet station for RANS
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2510(ZFSId bcId){
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
  const ZFSFloat maxIntegrationHeight = 2.0*m_rescalingBLT;

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

  //compute the local moment thickness j=direction of integration
  for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; ++j) {
    const ZFSId cellId = cellIndex(i,j);
    const ZFSId pointIdM1 = getPointIdFromCell(i,j);
    const ZFSId pointIdP1 = getPointIdFromPoint(pointIdM1,0,1);

    if(m_coordinates[1][pointIdM1] > maxIntegrationHeight) {
      continue;
    }

    const ZFSFloat urat= m_cells->pvariables[PV->U][cellId]/PV->UInfinity;
    const ZFSFloat momThick=m_cells->pvariables[PV->U][cellId]*m_cells->pvariables[PV->RHO][cellId]*fabs(F1-urat)/(CV->rhoUInfinity);
    //integrate normal to the wall
    const ZFSFloat ydist= m_coordinates[1][pointIdP1] - m_coordinates[1][pointIdM1];
    thetaLocal(0) += momThick*ydist;
  }

  MPI_Allreduce(thetaLocal.begin(),thetaGlobal.begin(),2, MPI_DOUBLE, MPI_SUM,rescalingCommGrComm);

  thetaGlobal(0) = zfsMAX(thetaGlobal(0),0.0000001);
  thetaGlobal(1) = zfsMAX(thetaGlobal(1),0.0000001);

  if(globalTimeStep%50==0&&m_block->m_RKStep==0&&m_block->domainId() == rescalingCommGrRootGlobal) {
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
  ZFSFloatScratchSpace wallPropertiesLocal(noWallProperties, __CALLING_FUNCTION__, "wallPropertiesLocalInlet");
  ZFSFloatScratchSpace wallProperties(noWallProperties, __CALLING_FUNCTION__, "wallPropertiesInlet");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  MPI_Allreduce(wallPropertiesLocal.begin(),wallProperties.begin(),noWallProperties, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  //////////////////////////////////////////////////////////////
  ///////////////// GAMS, UTAUIN, INNER OUTER COORD ////////////
  //////////////////////////////////////////////////////////////

  ZFSFloat utauIn = F0;
  ZFSFloat gams = F0;

  ZFSFloat utauRe=wallProperties(0);
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

  gams =pow(thetaGlobal(1)/fabs(thetaGlobal(0)),facc);
  gams = zfsMIN(gams,1.8);
  utauIn= utauRe*gams;

  ZFSFloatScratchSpace coordInInner(m_nCells[0], __CALLING_FUNCTION__, "coordInInner");
  ZFSFloatScratchSpace coordInOuter(m_nCells[0], __CALLING_FUNCTION__, "coordInOuter");

  for(ZFSId j=0; j<m_nCells[0]; ++j){
    const ZFSId cellId=cellIndex(i,j);
    const ZFSFloat rho=m_cells->pvariables[PV->RHO][cellId];
    const ZFSFloat frho=F1/rho;
    const ZFSFloat p = m_cells->pvariables[PV->P][cellId];
    const ZFSFloat temp=p*gamma*frho;
    const ZFSFloat mu=zfsSUTHERLANDLAW(temp);

    coordInInner(j)=utauIn*rho*m_cells->coordinates[1][cellId]/(mu*sqrt(m_block->m_Re0));
    coordInOuter(j)=m_cells->coordinates[1][cellId]*rho/(m_rescalingBLT*CV->rhoInfinity);
  }

  const ZFSId wallLocalOffset=m_block->m_nOffsetCells[0]; //Offset in j-direction
  ZFSFloat tempWallInletLocal = F0;
  ZFSFloat tempWallInletGlobal = F0;
  //determine the wall stuff if wall is contained whithin the partition
  if(wallLocalOffset==0 && m_block->m_nActiveCells[0]>=m_noGhostLayers) {
    const ZFSId cellId=cellIndex(i, m_noGhostLayers);
    tempWallInletLocal= temperature(cellId);
  }

  MPI_Allreduce(&tempWallInletLocal,&tempWallInletGlobal,1, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  //////////////////////////////////////////////////////////////
  ///////////////// NOW EXCHANGE VAR SLICE /////////////////////
  //////////////////////////////////////////////////////////////
  const ZFSId noVariables = PV->noVariables+1;
  ZFSId totalCells= m_block->m_totalGridBlockDim[0][0];
  ZFSFloatScratchSpace varSliceLocal(noVariables, totalCells, __CALLING_FUNCTION__, "varSliceLocal");
  ZFSFloatScratchSpace varSlice(noVariables, totalCells, __CALLING_FUNCTION__, "varSlice");

  //we are at the inlet, only fill with zeros
  varSlice.fill(F0);
  varSliceLocal.fill(F0);

  MPI_Allreduce(varSliceLocal.begin(),varSlice.begin(),noVariables*totalCells, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  // if(globalTimeStep%5==0&&m_block->m_RKStep==0&&m_block->domainId() == rescalingCommGrRootGlobal) {
    // cout << m_block->domainId()<< " ThetaInflow " << thetaGlobal(0) <<" ThetaRecyclingStation " << thetaGlobal(1) << endl;

    // FILE* f_channel;
    // f_channel = fopen("./varslice.dat", "w");
    //          for(ZFSId jj=0; jj<totalCells-1; ++jj) {
    //            for(ZFSId var=0; var<6; var++) {
    //              fprintf(f_channel, " %f", varSlice(var, jj));
    //            }
    //            fprintf(f_channel, "\n");
    //          }
    // fclose(f_channel);
  // }

  ///////////////////////////////////////////////////////////////
  ///////////////// RESCALING ///////////////////////////////////
  ///////////////////////////////////////////////////////////////

  const ZFSFloat ctem1= (F1+ctema)*(F1-POW2(gams));
  const ZFSFloat ctem2= F2*ctema*gams*(F1-gams);
  const ZFSFloat ctem3= (F1-gams)*(F1+gams+F2*ctema*gams);

  ZFSId jStart = 0;
  ZFSId jEnd = m_nCells[0];

  if(m_block->m_nOffsetCells[0]==0) {
    jStart = m_noGhostLayers;
  }
  if(m_block->m_nOffsetCells[0]+m_block->m_nActiveCells[0]==m_block->m_totalGridBlockDim[0][0]-1) {
    jEnd = m_nCells[0]-m_noGhostLayers;
  }

  ZFSFloat blEdgeVValueLocal = F0;
  ZFSBool edgePointIsSet = false;
  ZFSId edgePointJ = 0;


  for(ZFSId j=jStart; j<jEnd; ++j) {
    const ZFSId cellId = cellIndex(i,j);

    if(j>0) {
      if( coordInOuter(j-1) < 1.05 && coordInOuter(j) >= 1.05){
        blEdgeVValueLocal = m_cells->pvariables[PV->V][cellIndex(i,j-1)];
      }
    }

    if( coordInOuter(j) < 1.05){
      ZFSFloat uInner=F0,vInner=F0,TInner=F0, mutInner=F0;
      ZFSFloat uOuter=F0,vOuter=F0,TOuter=F0, mutOuter=F0;
      const ZFSFloat count = alpha*(coordInOuter(j) -b);
      const ZFSFloat denom = (F1-F2*b)*coordInOuter(j) +b;
      const ZFSFloat ratio = count/denom;
      const ZFSFloat wfun  = F1B2*(F1+tanh(ratio)/tanh(alpha));

      for(ZFSId jj=0; jj<totalCells-1; ++jj) {
        const ZFSId localId   =  jj;
        const ZFSId localIdP1 =  jj+1;

        const ZFSFloat yInnerRe   =  varSlice(3, localId);
        const ZFSFloat yInnerReP1 =  varSlice(3, localIdP1);

        if((yInnerRe-coordInInner(j)) < rescalEPS && yInnerReP1>coordInInner(j)) {
          const ZFSFloat dy1 = coordInInner(j)- yInnerRe;
          const ZFSFloat dy2 = yInnerReP1 - coordInInner(j);
          const ZFSFloat dy  = yInnerReP1 - yInnerRe;

          const ZFSFloat u     = varSlice(0,localId);
          const ZFSFloat uP1   = varSlice(0,localIdP1);
          const ZFSFloat v     = varSlice(1,localId);
          const ZFSFloat vP1   = varSlice(1,localIdP1);
          const ZFSFloat t     = varSlice(2,localId);
          const ZFSFloat tP1   = varSlice(2,localIdP1);
          const ZFSFloat mut   = varSlice(5,localId);
          const ZFSFloat mutP1 = varSlice(5,localIdP1);
          uInner= (uP1*dy1+u*dy2)/dy;
          vInner= (vP1*dy1+v*dy2)/dy;
          TInner= (tP1*dy1+t*dy2)/dy;
          mutInner = (mutP1*dy1+mut*dy2)/dy;
        }
      }

      //outer region
      for(ZFSId jj=0; jj<totalCells-1; ++jj) {
        const ZFSId localId = jj;
        const ZFSId localIdP1 = jj+1;

        const ZFSFloat yOuterRe   =varSlice(4, localId);
        const ZFSFloat yOuterReP1 =varSlice(4, localIdP1);

        if((yOuterRe-coordInOuter(j))<rescalEPS && yOuterReP1>coordInOuter(j)) {
          const ZFSFloat dy1= coordInOuter(j) - yOuterRe;
          const ZFSFloat dy2= yOuterReP1 - coordInOuter(j);
          const ZFSFloat dy = yOuterReP1 - yOuterRe;

          const ZFSFloat u     = varSlice(0,localId);
          const ZFSFloat uP1   = varSlice(0,localIdP1);
          const ZFSFloat v     = varSlice(1,localId);
          const ZFSFloat vP1   = varSlice(1,localIdP1);
          const ZFSFloat t     = varSlice(2,localId);
          const ZFSFloat tP1   = varSlice(2,localIdP1);
          const ZFSFloat mut   = varSlice(5,localId);
          const ZFSFloat mutP1 = varSlice(5,localIdP1);
          uOuter= (uP1*dy1+u*dy2)/dy;
          vOuter= (vP1*dy1+v*dy2)/dy;
          TOuter= (tP1*dy1+t*dy2)/dy;
          mutOuter = (mutP1*dy1+mut*dy2)/dy;
        }
      }

      const ZFSFloat TInnerA= POW2(gams)*TInner+ctem1*PV->TInfinity;
      const ZFSFloat TOuterA= POW2(gams)*TOuter -(ctem2*(uOuter/PV->UInfinity)-ctem3)*PV->TInfinity;

      //van Driest transformation
      const ZFSFloat uvdInner= PV->UInfinity * asin(b_vd *uInner/PV->UInfinity)/b_vd;
      const ZFSFloat uvdOuter= PV->UInfinity * asin(b_vd *uOuter/PV->UInfinity)/b_vd;

      //scaling of transformed inner and outer velocities
      //use uvd8, the van Driest transformed u8 value
      uInner = gams*uvdInner;
      uOuter = gams*uvdOuter+(F1-gams)*uvd8;

      uInner= PV->UInfinity*sin(b_vd*uInner/PV->UInfinity)/b_vd;
      uOuter= PV->UInfinity*sin(b_vd*uOuter/PV->UInfinity)/b_vd;

      //turbulent viscosity
      const ZFSFloat tempWallInlet = tempWallInletGlobal;
      const ZFSFloat tempWallRecycling = wallProperties(2);

      const ZFSFloat viscWallInlet =  zfsSUTHERLANDLAW(tempWallInlet);
      const ZFSFloat viscWallRecycling = zfsSUTHERLANDLAW(tempWallRecycling);
      const ZFSFloat thetaInlet = thetaGlobal(0);
      const ZFSFloat thetaRecycling = thetaGlobal(1);
      mutInner = mutInner*(viscWallInlet/viscWallRecycling);
      mutOuter = mutOuter*gams*(thetaInlet/thetaRecycling);

      const ZFSFloat uMean= uInner*(F1-wfun)+uOuter*wfun;
      const ZFSFloat vMean= vInner*(F1-wfun)+vOuter*wfun;
      const ZFSFloat tMean= TInnerA*(F1-wfun)+TOuterA*wfun;
      ZFSFloat mutMean = mutInner*(F1-wfun)+mutOuter*wfun;

      const ZFSFloat clebf = 6.6;
      const ZFSFloat blt   = m_rescalingBLT;
      const ZFSFloat cleb=F1/(F1+pow((m_cells->coordinates[1][cellId]/(clebf*blt)), 6.0));

      const ZFSFloat pres = PV->PInfinity;
      const ZFSFloat rhoIn = gamma*pres/tMean;

      m_cells->pvariables[PV->RHO][cellId]=rhoIn*cleb;
      m_cells->pvariables[PV->U][cellId]=uMean;
      m_cells->pvariables[PV->V][cellId]=vMean;
      m_cells->pvariables[PV->P][cellId]= pres;
      m_cells->pvariables[PV->RANS_VAR[0]][cellId] = mutMean/rhoIn;

    }else{
      if(!edgePointIsSet) {
        edgePointJ = j;
        edgePointIsSet = true;
      }
      //const ZFSFloat pres = pressure(cellIndex(m_noGhostLayers,j,k)); //for supersonic PV->PInfinity
      const ZFSFloat pres = PV->PInfinity;
      const ZFSFloat rhoIn = gamma*pres/PV->TInfinity;

      const ZFSFloat uMean = PV->UInfinity;
      const ZFSFloat vMean = PV->VInfinity;

      m_cells->pvariables[PV->RHO][cellId]= rhoIn;
      m_cells->pvariables[PV->U][cellId]= uMean;
      m_cells->pvariables[PV->V][cellId]= vMean;//m_cells->pvariables[PV->V][cellIndex(i,j-1)]*rhoIn;
      m_cells->pvariables[PV->P][cellId]= pres;
      m_cells->pvariables[PV->RANS_VAR[0]][cellId]= PV->ransInfinity[0];
    }
  }

  ZFSFloat blEdgeVValueGlobal = F0;
  MPI_Allreduce(&blEdgeVValueLocal,&blEdgeVValueGlobal,1, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  if(edgePointIsSet) {
    for(ZFSId j=edgePointJ; j<jEnd; ++j) {
      const ZFSId cellId = cellIndex(i,j);
      m_cells->pvariables[PV->V][cellId] = blEdgeVValueGlobal;
      const ZFSFloat pres = PV->PInfinity;
      m_cells->pvariables[PV->P][cellId]= pres;
    }
  }

  for(ZFSId j=0; j<m_nCells[0]; ++j) {
    //extrapolation for second GC
    const ZFSId cellId = cellIndex(1,j);
    const ZFSId cellIdM1 = cellIndex(0,j);
    const ZFSId cellIdadj = cellIndex(2,j);

    for(ZFSId var=0; var<PV->noVariables; var++) {
      m_cells->pvariables[var][cellIdM1] = 2.0*m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
    }
  }
}




//Recycling station for RANS
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2511(ZFSId bcId){
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  cout.precision(8);

  MPI_Comm rescalingCommGrComm = *m_block->m_rescalingCommGrComm;

  //scaling between delta0 and delta2
  const ZFSFloat F727= 72.0/7.0; //8.5, 8.0
  const ZFSId i= start[0];
  const ZFSFloat yWall=F0; //this has been fixed else method does not works
  const ZFSFloat maxIntegrationHeight = 2.0*m_rescalingBLT;

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE THETA AND EXCHANGE /////////////////
  //////////////////////////////////////////////////////////////

  //thetaLocal.fill(F0); //initialize scratch space to zero // only for parallel use
  ZFSFloatScratchSpace thetaLocal(2, __CALLING_FUNCTION__, "thetaLocalRe");
  ZFSFloatScratchSpace thetaGlobal(2, __CALLING_FUNCTION__, "thetaGlobalRe");
  thetaLocal.fill(F0); //initialize scratch space
  thetaGlobal.fill(F0);

  //the offest position in k-direction is the offset
  //compute the local moment thickness j=direction of integration
  for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; ++j) {
    const ZFSId cellId = cellIndex(i,j);
    const ZFSId pointIdM1 = getPointIdFromCell(i,j);
    const ZFSId pointIdP1 = getPointIdFromPoint(pointIdM1,0,1);

    if(m_coordinates[1][pointIdM1] > maxIntegrationHeight) {
      continue;
    }

    const ZFSFloat urat= m_cells->pvariables[PV->U][cellId]/PV->UInfinity;
    const ZFSFloat momThick=m_cells->pvariables[PV->U][cellId]*m_cells->pvariables[PV->RHO][cellId]*fabs(F1-urat)/(CV->rhoUInfinity);

    //integrate normal to the wall
    const ZFSFloat ydist= m_coordinates[1][pointIdP1] - m_coordinates[1][pointIdM1];
    thetaLocal(1) += momThick*ydist;
  }


 //communicate the Thickness across the plane
  MPI_Allreduce(thetaLocal.begin(), thetaGlobal.begin(), 2 ,MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE WALL PROPERTIES ///////
  //////////////////////////////////////////////////////////////

  const ZFSFloat delta = F727*thetaGlobal(1);
  const ZFSId noVar =3; //for more variables if wanted
  const ZFSId wallLocalOffset=m_block->m_nOffsetCells[0]; //Offset in j-direction

  ZFSFloatScratchSpace wallPropertiesLocal(noVar, __CALLING_FUNCTION__, "wallPropertiesLocalRe");
  ZFSFloatScratchSpace wallProperties(noVar, __CALLING_FUNCTION__,"wallPropertiesRe");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  //determine the wall stuff if wall is contained whithin the partition
  if(wallLocalOffset==0 && m_block->m_nActiveCells[1]>=m_noGhostLayers) {
    const ZFSId cellId=cellIndex(i, m_noGhostLayers);
    const ZFSFloat rho = m_cells->pvariables[PV->RHO][cellId];
    const ZFSFloat p = m_cells->pvariables[PV->P][cellId];
    const ZFSFloat t=p*m_block->m_gamma/rho;
    const ZFSFloat mu=zfsSUTHERLANDLAW(t);
    const ZFSFloat uWall = fabs(m_cells->pvariables[PV->U][cellId]);
    const ZFSFloat ydist=m_cells->coordinates[1][cellId]-yWall;
    const ZFSFloat uTau= sqrt(uWall*mu/(ydist*rho));

    wallPropertiesLocal(0)= uTau;
    wallPropertiesLocal(1)= rho;
    wallPropertiesLocal(2)= t;
  }

  MPI_Allreduce(wallPropertiesLocal.begin(),wallProperties.begin(),noVar, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  ZFSFloat tempWallInletLocal = F0;
  ZFSFloat tempWallInletGlobal = F0;
  MPI_Allreduce(&tempWallInletLocal,&tempWallInletGlobal,1, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE VAR SLICE /////////////
  //////////////////////////////////////////////////////////////

  ZFSId totalCells=m_block->m_totalGridBlockDim[0][0];//-1;

  const ZFSId noVariables = PV->noVariables+1;
  ZFSFloatScratchSpace varSliceLocal(noVariables, totalCells, __CALLING_FUNCTION__, "varSliceLocal");
  ZFSFloatScratchSpace varSlice(noVariables, totalCells, __CALLING_FUNCTION__, "varSlice");

  for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; ++j) {
    const ZFSId cellId=cellIndex(i,j);
    const ZFSFloat rho=m_cells->pvariables[PV->RHO][cellId];
    const ZFSFloat frho=F1/rho;
    const ZFSFloat p = pressure(cellId);
    const ZFSFloat temp=p*m_block->m_gamma*frho;
    const ZFSFloat mu=zfsSUTHERLANDLAW(temp);
    const ZFSFloat uTauRe=wallProperties(0);
    const ZFSFloat yIn=(m_cells->coordinates[1][cellId]-yWall)*uTauRe*rho/(mu*sqrt(m_block->m_Re0));
    const ZFSFloat yOut=(m_cells->coordinates[1][cellId]-yWall)*rho/(delta*CV->rhoInfinity);
    const ZFSFloat u=m_cells->pvariables[PV->U][cellId];
    const ZFSFloat v=m_cells->pvariables[PV->V][cellId];

    //>RANS
    const ZFSFloat mut = m_cells->pvariables[PV->RANS_VAR[0]][cellId]*rho;
    //<RANS

    const ZFSId localId = m_block->m_nOffsetCells[0] + j-m_noGhostLayers + 1;

    //save the variables u,v,w,t,yI,yO
    varSliceLocal(0,localId) = u ;
    varSliceLocal(1,localId) = v ;
    varSliceLocal(2,localId) = temp ;
    varSliceLocal(3,localId) = yIn ;
    varSliceLocal(4,localId) = yOut ;
    varSliceLocal(5,localId) = mut ;
  }

  //set first value at the wall manually
  if(m_block->m_nOffsetCells[0]==0) {
    varSliceLocal(0,0) = 0.0 ;
    varSliceLocal(1,0) = 0.0 ;
    varSliceLocal(2,0) = varSliceLocal(2,1);
    varSliceLocal(3,0) = 0.0 ;
    varSliceLocal(4,0) = 0.0 ;
    varSliceLocal(5,0) = varSliceLocal(5,1) ;
  }

  //communicate the slice
  MPI_Allreduce(varSliceLocal.begin(),varSlice.begin(),noVariables*totalCells, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);

  ZFSFloat blEdgeVValueLocal = F0;
  ZFSFloat blEdgeVValueGlobal = F0;
  MPI_Allreduce(&blEdgeVValueLocal,&blEdgeVValueGlobal,1, MPI_DOUBLE, MPI_SUM, rescalingCommGrComm);
}

/** Prescribe given profile BC
 *
 *  Precribes a profile from the restart file
 *  extrapolate pressure from computational domain
 */
template <ZFSBool isRans>
void ZFSStrctrdBndryCnd2D<isRans>::bc2600(ZFSId bcId){
  ZFSInt* start = m_physicalBCMap[bcId]->start1;
  ZFSInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face)
    {
    case 0:
    {
      for(ZFSInt i = start[0]; i<end[0]; i++) {
        for(ZFSInt j = start[1]; j<end[1]; j++) {
          ZFSId cellId=cellIndex(m_noGhostLayers-1-i,j);
          ZFSId cellIdadj=cellIndex(m_noGhostLayers-i,j);

          //extrapolate pressure to ghost cells
          m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
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
inline ZFSFloat ZFSStrctrdBndryCnd2D<isRans>::pressure(ZFSId cellId){
  return m_cells->pvariables[PV->P][cellId]; 
}

template <ZFSBool isRans>
inline ZFSFloat ZFSStrctrdBndryCnd2D<isRans>::temperature(ZFSId cellId){
  const ZFSFloat gamma = m_block->m_gamma;
  ZFSFloat t = gamma*m_cells->pvariables[PV->P][cellId]/m_cells->pvariables[PV->RHO][cellId];
  return t;
}


/**
 *     New function to compute the skin friction and pressure coefficient
 *     and the part for the lift and drag coefficient
 *     works in all directions and offers line averaging in one direction
 *     @author: Pascal Meysonnat, Marian Albers
 *     @date: 01.01.1010  
 */

/*commented out but needs to be done 
template <ZFSBool isRans>
template <ZFSBool computePower>
void ZFSStrctrdBndryCnd2D<isRans>::computeFrictionPressureCoef(){
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
  for(ZFSId i =0; i< (ZFSInt)nDim*noWalls; ++i){
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
    //const ZFSInt mapOffsetPower = m_cells->cfOffsets[map];
    ZFSInt* start = m_auxDataMap[map]->start1;
    ZFSInt* end = m_auxDataMap[map]->end1;
    
    ZFSInt n=0;
    ZFSId normal = 0;
    //area
    ZFSFloat area = F0;    
    ZFSFloat tempArea = F0;
    //skin-friction
    ZFSFloat cf[2]={F0,F0};
    ZFSFloat tempCf[2]={F0,F0};
    //pressure coefficient
    ZFSFloat cp[2]={F0,F0};
    ZFSFloat tempCp[2]={F0,F0};
    //power computation
    ZFSFloat powerp[2]={F0,F0};
    ZFSFloat powerv[2]={F0,F0};
    ZFSFloat tempPowerV[2]={F0,F0};
    ZFSFloat tempPowerP[2]={F0,F0};

    //indices
    ZFSId n10[2]={0,0};
    ZFSId n01[2]={0,0};
    ZFSInt n1m1[2]={0,0};
    ZFSId pCoordDir[9]={0,0,0,0,0,0,0,0,0}; //be changed to 2d
    ZFSId firstTangential = -1, secondTangential = -1;
    ZFSId i=0, j=0;
    ZFSInt jj = 0, kk = 0; //be changed to 2d
    ZFSId *reali=nullptr, *realj=nullptr;
    ZFSInt *realjj=nullptr, *realkk=nullptr;
    ZFSInt sizeJ = 0, sizeK = 0; //be changed to 2d

    ZFSFloat dxidx=-F1, dxidy=-F1;
    ZFSFloat detadx=-F1, detady=-F1;
    ZFSFloat dzetadx=-F1, dzetady=-F1;

    ZFSFloat supportVec[2] = {F0,F0};
    ZFSFloat firstVec[2] = {F0,F0};
    ZFSFloat secondVec[2] = {F0,F0};
    ZFSFloat normalVec[2] = {F0,F0};
    ZFSFloat cellVec[2] = {F0,F0};

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
        const ZFSFloat dx=sqrt(POW2(m_cells->coordinates[normal][cellIdP1]-m_cells->coordinates[normal][cellId])+
                               POW2(m_cells->coordinates[firstTangential][cellIdP1]-m_cells->coordinates[firstTangential][cellId])+
                               POW2(m_cells->coordinates[secondTangential][cellIdP1]-m_cells->coordinates[secondTangential][cellId]));
        const ZFSFloat dx1=sqrt(POW2(m_cells->coordinates[normal][cellId]-xRef)+
                                POW2(m_cells->coordinates[firstTangential][cellId]-yRef)+
                                POW2(m_cells->coordinates[secondTangential][cellId]-zRef));
        //pressures
        const ZFSFloat p1 = m_cells->pvariables[PV->P][cellId];
        const ZFSFloat p2 = m_cells->pvariables[PV->P][cellIdP1];
        //extrapolation to the face
        const ZFSFloat pW=((p1-p2)/dx)*dx1+p1;
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
          m_cells->powerPres[mapOffsetCf+0*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_px/stagnationEnergy;
          m_cells->powerPres[mapOffsetCf+1*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_py/stagnationEnergy;
          m_cells->powerPres[mapOffsetCf+2*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_pz/stagnationEnergy;
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
        const ZFSFloat tBc=((T1-T2)/dx)*dx1+T1;
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
        const ZFSFloat det = 1.0 / m_cells->cellJac[cellId];

        //compute contravariant velocities parallel to the wall
        const ZFSFloat du = (dxidx*u1 + dxidy*v1 + dxidz*w1) - m_cells->dxt[0][cellId];
        const ZFSFloat dv = (detadx*u1 + detady*v1 + detadz*w1) - m_cells->dxt[1][cellId];
        const ZFSFloat dw = (dzetadx*u1 + dzetady*v1 + dzetadz*w1) - m_cells->dxt[2][cellId];
       

        //compute gradient
        const ZFSFloat dudx=du/orthDist;
        const ZFSFloat dvdx=dv/orthDist;
        const ZFSFloat dwdx=dw/orthDist;

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
          //Power due to viscous forces in the computational space ( tau_w*n*u)
          const ZFSFloat P_cfx = (P_cfuu*dxdxi + P_cfvv*dxdeta + P_cfww*dxdzeta) * det*gridVel[0];//PV->UInfinity;
          const ZFSFloat P_cfy = (P_cfuu*dydxi + P_cfvv*dydeta + P_cfww*dydzeta) * det*gridVel[1];//PV->VInfinity;
          const ZFSFloat P_cfz = (P_cfuu*dzdxi + P_cfvv*dzdeta + P_cfww*dzdzeta) * det*gridVel[2];//PV->WInfinity;
          //save power to map
          m_cells->powerVisc[mapOffsetCf+0*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_cfx/stagnationEnergy;
          m_cells->powerVisc[mapOffsetCf+1*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_cfy/stagnationEnergy;
          m_cells->powerVisc[mapOffsetCf+2*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]=P_cfz/stagnationEnergy;
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
        tempPowerP[0]+=m_cells->powerPres[mapOffsetCf+0*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*considerValue;
        tempPowerP[1]+=m_cells->powerPres[mapOffsetCf+1*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*considerValue;
        tempPowerP[2]+=m_cells->powerPres[mapOffsetCf+2*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*considerValue;
        //viscous power
        //this should be the valid method as in TFS (this is the t_w*n*u*da)
        tempPowerV[0]+=  m_cells->powerVisc[mapOffsetCf+0*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*(sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz))*considerValue;
        tempPowerV[1]+=  m_cells->powerVisc[mapOffsetCf+1*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*(sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz))*considerValue;
        tempPowerV[2]+= m_cells->powerVisc[mapOffsetCf+2*sizeJ*sizeK+(*realjj)+(*realkk)*sizeJ]*(sqrt(dxidx*dxidx + dxidy*dxidy + dxidz*dxidz))*considerValue;

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
      powerv[0]+=tempPowerV[0];
      powerv[1]+=tempPowerV[1];
      powerv[2]+=tempPowerV[2];
      powerp[0]+=tempPowerP[0];
      powerp[1]+=tempPowerP[1];
      powerp[2]+=tempPowerP[2];


      area+=tempArea;

      for(int m = 0; m<nDim; ++m) {
        tempCp[m] = F0;
        tempCf[m] = F0;
        tempPowerV[m] =F0;
        tempPowerP[m] =F0;
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
template void ZFSStrctrdBndryCnd2D<true>::computeFrictionPressureCoef<true>();
template void ZFSStrctrdBndryCnd2D<false>::computeFrictionPressureCoef<true>();
template void ZFSStrctrdBndryCnd2D<true>::computeFrictionPressureCoef<false>();
template void ZFSStrctrdBndryCnd2D<false>::computeFrictionPressureCoef<false>();

*/


template class ZFSStrctrdBndryCnd2D<true>;
template class ZFSStrctrdBndryCnd2D<false>;
