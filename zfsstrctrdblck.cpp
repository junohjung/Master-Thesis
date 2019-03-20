#include "zfsstrctrdblck.h"
#include "zfsglobals.h"
#include "zfsiolib.h"
#include "zfsstrctrdblckpartition.h"
#include "zfsstrctrdblckwindowinfo.h"
#include <cstdlib>
#include <cmath>
#if not defined(ZFS_MS_COMPILER)
#include <unistd.h>
#endif
//temporaray
#include <vector>
//#include <algorithm>

/** \brief Constructor of the structured block
 * \author Pascal Meysonnat
 * \date 01.01.1010
 */
template <ZFSInt nDim>
ZFSStrctrdBlck<nDim>::ZFSStrctrdBlck ( ZFSId blockId, ZFSBool* propertiesGroups, const MPI_Comm comm)
  : ZFSBlock(blockId, comm)
  , ZFSStrctrdPostprocessing < nDim, ZFSStrctrdBlck<nDim> > ()
  , m_eps(std::numeric_limits<ZFSFloat>::epsilon())
{
  (void) propertiesGroups;
  const ZFSLong oldAllocatedBytes = allocatedBytes();

  if(domainId() == 0) {cout << "Initializing Structured Block..." << endl;}

  // for the zonal approach we will have to split the communicator
  m_zfsStrctrdComm = comm;

  // intialize the most important properties
  initializeStrctrdBlck(propertiesGroups);

  // numericals method:
  setNumericalProperties();

  // set testcase parameters
  setTestcaseProperties();

  // set moving grid parameters
  setMovingGridProperties();

  // initialize all timers
  initializeTimers();

  // initialize cell
  m_cells = new ZFSStrctrdCell;

  ///////////////////////////////////////////////////
  ///////////////// OPEN GRID FILE //////////////////
  ///////////////////////////////////////////////////

  /*! \page propertyPage1
    \section gridInputFileName
    <code>ZFSString ZFSStrctrdBlck::m_gridInputFileName </code>\n
    default = <code>""</code>\n \n
    Name of the grid file.\n
    Keywords: <i>GRID, STRCTRD</i>
  */
  m_gridInputFileName=*(ZFSContext::getProperty("gridInputFileName", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString(0));
  ZFSId file_id = -1;
  //unique identifier needed to associate grid and solution in case of restart
  const char* aUID= new char[18];
  file_id = io_openfile("hdf5", m_gridInputFileName.c_str(), "collective", m_zfsStrctrdComm);
  ZFSFloat dummyBlocks=F0;
  io_read_dattribute1(file_id, "", "noBlocks", &dummyBlocks);
  m_noInputBlocks=dummyBlocks;
  io_read_sattribute1(file_id, "", "UID", aUID );
  m_uID=aUID;
  delete[] aUID;

  ///////////////////////////////////////////////////
  /////////////// GRID DECOMPOSITIONING /////////////
  ///////////////////////////////////////////////////
  if(domainId() == 0) {cout << "Doing  block decomposition..." << endl;}
  RECORD_TIMER_START(m_tdecomposition);
  m_partition = new ZFSStrctrdDecomposition<nDim>(m_noInputBlocks, m_gridInputFileName, noDomains());

  /*! \page propertyPage1
    \section readPartitionFromFile
    <code>ZFSInt ZFSStrctrdBlck::m_readDecompositionFromFile </code>\n
    default = <code>0</code>\n \n
    Trigger the use to read the MPI partitioning from a file (faster).\n
    possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>PARALLEL, PARTITIONING, STRCTRD</i>
  */
  m_readDecompositionFromFile=0;
  m_readDecompositionFromFile=*(ZFSContext::getProperty("readPartitionFromFile", m_blockId, __CALLING_FUNCTION__, &m_readDecompositionFromFile)->asInt(0));
  if(m_readDecompositionFromFile){
    zfs_log << "reading repartition from File ...." << endl;
    int success= m_partition->readFromFile();
    if(!success){
      zfs_log << "..... Reading from partition file FAILED --> new decomposition activated" << endl;
      m_partition->decompose();
    }
    zfs_log << ".... reading in successful " << endl;
  }else{
    m_partition->decompose();
  }
  m_totalGridCells=0;
  zfsAlloc(m_totalGridBlockCells, m_noInputBlocks, "m_totalGridBlockCells", (ZFSLong)-1, __CALLING_FUNCTION__);
  zfsAlloc(m_totalGridBlockDim, m_noInputBlocks, 3, "m_totalGridBlockDim", -1, __CALLING_FUNCTION__);
  for(ZFSId i=0; i< m_noInputBlocks; i++) {
    ZFSLong temp =1;
    for(ZFSId dim=0; dim<nDim; dim++) {
      temp*=(m_partition->inputBoxInfo[i]->DirLast[dim]);
      m_totalGridBlockDim[i][dim]=(m_partition->inputBoxInfo[i]->DirLast[dim])+1;//number of points in the grid File
    }
    m_totalGridBlockCells[i]=temp;
    if(temp!=1) m_totalGridCells+=temp;
  }
  RECORD_TIMER_STOP(m_tdecomposition);
  if(domainId() == 0) {cout << "Doing  block decomposition... SUCCESSFUL!" << endl;}

  zfs_log << "Setting zonal properties..." << endl;
  setZonalProperties();
  zfs_log << "Setting zonal properties... SUCCESSFUL!" << endl;

  zfs_log <<  "Allocating large fields..." << endl;
  allocateVariables();
  zfs_log << "Allocating large fields... SUCCESSFUL!" << endl;

  zfs_log << "Reading input output properties..." << endl;
  setInputOutputProperties();
  zfs_log << "Reading input output properties... SUCCESSFUL!" << endl;

  //////////////////////////////////////////////////////////////
  ///////////////// ASSIGN WINDOW INFORMATION //////////////////
  //////////////////////////////////////////////////////////////
  //4) assign window winformation (3D-2D possible because of dimension)
  //  but be careful, only 3D has been tested and implemented
  m_inputBlockId = m_partition->outputBoxInfo[getBoxId(domainId())]->inputBoxID;
  m_windowInfo = new ZFSStrctrdBlckWindowInfo<nDim>(m_partition, file_id,
                                                    m_noInputBlocks,
                                                    m_inputBlockId,
                                                    m_noGhostLayers,
                                                    m_zfsStrctrdComm,
                                                    noDomains(),
                                                    domainId());
                                                  //4.1) read in the original window information

  zfsAlloc(m_periodicDisplacements, nDim*nDim, "m_periodicDisplacement", F0, __CALLING_FUNCTION__);
  m_windowInfo->readNewWindowInfo(m_gridInputFileName);
  m_windowInfo->initGlobals();
  m_windowInfo->readWindowCoordinates(domainId(),
                                      m_periodicDisplacements,
                                      m_gridInputFileName);

  createMPIGroups();

  readAndSetSpongeLayerProperties();
  if(m_zonal || m_rans) {m_windowInfo->setWallInformation();}  //junoh

  m_singularity =new SingularInformation[30];

  //Create the mapping for the boundary and domain windows
  //the windows are needed for boundary condition application
  //and data exchange between the domains
  m_windowInfo->createWindowMapping(m_commChannelIn, m_commChannelOut,
                                    m_commChannelWorld, m_channelRoots,
                                    m_commStg, m_commStgRoot,
                                    m_commStgRootGlobal,
				    m_commBC2600, m_commBC2600Root,   //junoh
				    m_commBC2600RootGlobal,
                                    m_rescalingCommGrComm,
                                    m_rescalingCommGrRoot,
                                    m_rescalingCommGrRootGlobal,
                                    m_commPerRotOne, m_commPerRotTwo,
                                    m_commPerRotWorld, m_commPerRotRoots,
                                    m_commPerRotGroup,
                                    m_singularity, &m_hasSingularity);

  //Creating the communication flags needed
  //for all inter-domain exchange operations
  m_cmnctnFlag = new ZFSStrctrdCommunicationHandle(nDim);
  m_windowInfo->createCommunicationExchangeFlags(m_cmnctnFlag, CV->noVariables);

  //set the send and recieve buffer sizes;
  m_cmnctnFlag->setBufferSizes();

  if(m_nonBlockingComm) prepareNonBlockingCommunication();


  m_periodicConnection=0;
  //set flag for the periodic exchange
  for(ZFSId i=0; i<(ZFSInt)m_windowInfo->rcvMap.size(); i++) {
    if(m_windowInfo->rcvMap[i]->BC>=4000 && m_windowInfo->rcvMap[i]->BC<=4999) {
      m_periodicConnection=1;
      break;
    }
  }

  ///////////////////////////////////////////////////
  ///////////////// READ GRID ///////////////////////
  ///////////////////////////////////////////////////
  //5.1) calculate the no of grid points of the partition and the ijk_max

  zfsAlloc(m_nPoints, nDim, "m_nPoints", -1 , __CALLING_FUNCTION__);
  zfsAlloc(m_nActivePoints, nDim, "m_nActivePoints", -1 , __CALLING_FUNCTION__);
  zfsAlloc(m_nCells, nDim, "m_nCells", -1 , __CALLING_FUNCTION__);
  zfsAlloc(m_nActiveCells, nDim, "m_nActiveCells", -1 , __CALLING_FUNCTION__);
  zfsAlloc(m_nOffsetCells, nDim, "m_nOffsetCells", -1 , __CALLING_FUNCTION__);
  zfsAlloc(m_nInputBlockCells, nDim, "m_nInputBlockCells", -1, __CALLING_FUNCTION__);

  m_noGridPoints=1;
  m_noStrctrdActiveCells=1;
  m_noStrctrdCells=1;

  for (ZFSInt i=0;i<nDim;i++) {
    m_nOffsetCells[i]=m_partition->outputBoxInfo[getBoxId(domainId())]->offset[i];
    m_nActivePoints[i]=m_partition->outputBoxInfo[getBoxId(domainId())]->DirLast[i];
    m_nPoints[i]=m_nActivePoints[i]+2*m_noGhostLayers;
    //change dimension to that of cell centered!!
    m_nActiveCells[i]= m_nActivePoints[i]-1;
    m_nCells[i]=m_nPoints[i]-1;
    m_noGridPoints *= m_nPoints[i];
    m_nInputBlockCells[i] = m_partition->inputBoxInfo[m_inputBlockId]->DirLast[i];
    m_noStrctrdActiveCells *= m_nActiveCells[i];
    m_noStrctrdCells *= m_nCells[i];
  }


  allocateAndInitBlockMemory();

  ZFSInt averageCellsPerDomain = (ZFSInt) (m_totalGridCells / noDomains());
  ZFSFloat localDeviation = ((ZFSFloat) averageCellsPerDomain - (ZFSFloat)m_noStrctrdActiveCells ) / ((ZFSFloat) averageCellsPerDomain);
  ZFSFloat globalMaxDeviation = F0;
  ZFSFloat globalMinDeviation = F0;
  ZFSFloat globalAvgDeviation = F0;
  ZFSInt globalMaxNoCells = 0;
  ZFSInt globalMinNoCells = 0;
  MPI_Allreduce(&localDeviation, &globalMaxDeviation, 1, MPI_DOUBLE, MPI_MAX, m_zfsStrctrdComm);
  MPI_Allreduce(&localDeviation, &globalMinDeviation, 1, MPI_DOUBLE, MPI_MIN, m_zfsStrctrdComm);
  MPI_Allreduce(&localDeviation, &globalAvgDeviation, 1, MPI_DOUBLE, MPI_SUM, m_zfsStrctrdComm);
  MPI_Allreduce(&m_noStrctrdActiveCells, &globalMaxNoCells, 1, MPI_INT, MPI_MAX, m_zfsStrctrdComm);
  MPI_Allreduce(&m_noStrctrdActiveCells, &globalMinNoCells, 1, MPI_INT, MPI_MIN, m_zfsStrctrdComm);
  globalAvgDeviation /= noDomains();

  if(domainId() == 0) {
    cout << "///////////////////////////////////////////////////////////////////" << endl
         << "Total no. of grid cells: " << m_totalGridCells << endl
         << "Average cells per domain: " << averageCellsPerDomain <<  endl
         << "Max no of cells per domain: " << globalMaxNoCells << endl
         << "Min no of cells per domain: " << globalMinNoCells << endl
         << "Average deviation from average: " <<  globalAvgDeviation * 100.0 << " percent" << endl
         << "Maximum deviation from average: +" <<  globalMaxDeviation * 100.0 << " / "
         << globalMinDeviation * 100.0 << " percent" << endl
         << "///////////////////////////////////////////////////////////////////" << endl;
  }

  //read the grid from partition and
  //and move coordinates to right position
  if(domainId() == 0) { cout << "Reading Grid..." << endl; }
  zfs_log << "->reading the grid file" << endl;
  RECORD_TIMER_START(m_treadGrid);
  readGrid(file_id);
  RECORD_TIMER_STOP(m_treadGrid);
  zfs_log << "------------- Grid read successfully! -------------- " << endl;
  if(domainId() == 0) { cout << "Reading Grid SUCCESSFUL!" << endl; }
  io_closefile(file_id);
  moveCellPoints();


  //get the zonal BC information
  m_windowInfo->setZonalBCInformation();    //junoh
  // set properties for synthetic turbulence generation method
  setSTGProperties();

  // set the properties for bc2600 (if existing)
  setProfileBCProperties();
  
  //initialize the postprocessing class
  initStrctrdPostprocessing();

  //print allocated scratch memory
  if(domainId() == 0) {;
    ZFSFloat scratchMemory = (ZFSScratch::getTotalMemory()/1024.0)*noDomains();
    ZFSString memoryUnit = " KB";
    if(scratchMemory > 1024.0) {
      scratchMemory /= 1024.0;
      memoryUnit = " MB";
    }
    if(scratchMemory > 1024.0) {
      scratchMemory /= 1024.0;
      memoryUnit = " GB";
    }
    cout << "=== Total global scratch space memory: " << setprecision(2) << fixed
         << scratchMemory << memoryUnit << " ===" << endl;
  }
  printAllocatedMemory( oldAllocatedBytes, "ZFSStrctrdBlck", m_zfsStrctrdComm );

}


template <ZFSInt nDim>
ZFSStrctrdBlck<nDim>::~ZFSStrctrdBlck()
{
  delete m_commChannelIn;
  delete m_commChannelOut;
  delete m_commChannelWorld;
  delete [] m_channelRoots;

  delete m_commStg;
  delete m_commStgRoot;
  delete m_commStgRootGlobal;
  
  
  delete m_commBC2600; //junoh
  delete m_commBC2600Root;
  delete m_commBC2600RootGlobal;


  delete m_commPerRotOne;
  delete m_commPerRotTwo;
  delete m_commPerRotWorld;
  delete [] m_commPerRotRoots;

  delete m_rescalingCommGrComm;
  delete m_rescalingCommGrRoot;
  delete m_rescalingCommGrRootGlobal;

  delete m_cmnctnFlag;
  delete CV;
  delete PV;
  delete FQ;
  delete m_partition;
  delete m_windowInfo;
  delete m_cells;
}

/**
 * Counts the number of necessary FQ fields, allocates them
 * and corrects the indexes of the FQ variable pointers
 * \author Marian Albers
 * \changes Pascal Meysonnat
 * \date 17.09.2015
 */

ZFSStrctrdZonalBC::~ZFSStrctrdZonalBC() //junoh
{
  for(ZFSInt i=0; i<m_noSndNghbrDomains; i++) {
    delete[] m_bufferSndZonal[i];
    delete[] m_bufferSndMapCellId[i];
  }
  
  for(ZFSInt i=0; i<m_noRcvNghbrDomains; i++) {
    delete[] m_bufferRcvZonal[i];
    delete[] m_interpolatedVars[i];
    delete[]m_bufferRcvMapCellId[i];
  }
  
  delete [] m_bufferSndZonal;
  delete [] m_bufferRcvZonal;
  delete [] m_bufferSndMapCellId;
  delete [] m_bufferRcvMapCellId;

  delete []m_localBufferSndSize;
  delete []m_localBufferRcvSize;
  delete []m_localRcvId;
  delete []m_localSndId;
  delete []m_globalRcvZonalId;
  delete []m_globalSndZonalId;
  delete []m_localBufferMapCellId;
  delete []m_localBufferIndexCellId;

  delete []m_localCommReceiverIds;

  delete []m_interpolatedVars;
  delete []m_interpolatedVarsAV;
  delete []m_globalReceiverIds;
  // delete []m_globalLocalCellIds;
  delete []m_globalLocalMapCellIds;
  // delete []m_hasPartnerLocalBC;
  // delete []m_hasPartnerGlobalBC;
  // delete []m_localHasPartnerGlobalBC;
  delete []m_localMapCellsId;
  
  delete []m_zonalBCCells;

  delete []mpi_sndRequest;
  delete []mpi_rcvRequest;
  delete []mpi_sndStatus;
  delete []mpi_rcvStatus;

}

/**
 * Delete for allocating Zonal grid creation member variables
 * \author Junoh Jung
 * \date 07.2018
 */




template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::initializeFQField()
{
  TRACE();
  //count the number of needed FQ fields and allocate
  ZFSId noFQFieldsNeeded = 0;
  for(ZFSId i=0; i<FQ->maxNoFQVariables; i++) {
    noFQFieldsNeeded += FQ->neededFQVariables[i];
  }

  FQ->noFQVariables = noFQFieldsNeeded;

  zfs_log << "Allocating " << noFQFieldsNeeded << " FQ fields..." << endl;
  zfsAlloc(m_cells->fq, noFQFieldsNeeded, m_noStrctrdCells, "m_cells->fq", F0, __CALLING_FUNCTION__);
  zfsAlloc(FQ->loadedFromRestartFile, noFQFieldsNeeded, "FQ->loadedFromRestartFile", false, __CALLING_FUNCTION__);
  zfs_log << "Allocating " << noFQFieldsNeeded << " FQ fields... SUCCESSFUL" << endl;

  ZFSId currentPos = 0;
  for(ZFSId i=0; i<FQ->maxNoFQVariables; i++) {
    if(FQ->neededFQVariables[i] ==0)
      continue;

    FQ->activateFQField(i,currentPos,FQ->outputFQVariables[i], FQ->planeOutputFQVariables[i], FQ->boxOutputFQVariables[i]);

    FQ->noFQPlaneOutput += (ZFSId) FQ->planeOutputFQVariables[i];
    FQ->noFQBoxOutput   += (ZFSId) FQ->boxOutputFQVariables[i];
    currentPos++;
  }

  
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::allocateAndInitBlockMemory() {
  zfsAlloc(pointProperties, 2, m_noGridPoints, "pointPropertiesDummy", __CALLING_FUNCTION__);
  zfsAlloc(m_coordinates, nDim, m_noGridPoints, "m_coordinates", -1.01010101 , __CALLING_FUNCTION__ );

  if(m_movingGrid) {
    zfsAlloc(m_mgOldCoordinates, nDim, m_noGridPoints, "m_mgOldCoordinates", -1.01010101 , __CALLING_FUNCTION__ );
    zfsAlloc(m_mgInitCoordinates, nDim, m_noGridPoints, "m_mgInitCoordinates", -1.01010101 , __CALLING_FUNCTION__ );

    if (m_travelingWave) {
      zfsAlloc(m_tempWaveSample, (PV->noVariables+(2*nDim-3)), m_noStrctrdCells, "m_tempWaveSample", F0, FUN_);
    }
  }

  if(m_localTimeStep) {
    zfsAlloc(m_cells->localTimeStep, m_noStrctrdCells, "m_cells->localTimeStep", -1.01010101 , __CALLING_FUNCTION__ );
  }

  zfsAlloc(m_cells->coordinates, nDim, m_noStrctrdCells, "m_cells->coordinates", __CALLING_FUNCTION__);

  zfsAlloc(m_cells->variables, m_maxNoVariables, m_noStrctrdCells, "m_cells->variables", __CALLING_FUNCTION__);//junoh
  zfsAlloc(m_cells->pvariables, m_maxNoVariables, m_noStrctrdCells, "m_cells->pvariables", __CALLING_FUNCTION__);//junoh
  zfsAlloc(m_cells->temperature, m_noStrctrdCells, "m_cells->temperature", __CALLING_FUNCTION__);
  zfsAlloc(m_cells->lamvisc, m_noStrctrdCells, "m_cells->lamvisc", __CALLING_FUNCTION__);
  zfsAlloc(m_cells->oldVariables, m_maxNoVariables, m_noStrctrdCells, "m_cells->oldVariables", __CALLING_FUNCTION__); //junoh
  zfsAlloc(m_cells->ds, m_noStrctrdCells, "m_cells->ds", F0,  __CALLING_FUNCTION__);
  zfsAlloc(m_cells->rightHandSide, CV->noVariables, m_noStrctrdCells, "m_cells->rhs", F0,  __CALLING_FUNCTION__);
  zfsAlloc(m_cells->flux, CV->noVariables*m_noStrctrdCells, "m_cells->flux", 10000.0, __CALLING_FUNCTION__);
  zfsAlloc(m_cells->cellVolume, m_noStrctrdCells, "m_cells->cellVolume", F0, __CALLING_FUNCTION__);

#ifdef ZFS_EXTRA_DEBUG
  zfsAlloc(viscFluxOut, nDim, m_noStrctrdCells*CV->noVariables, "viscousFlux Debug", __CALLING_FUNCTION__);
  zfsAlloc(convFluxOut, nDim, m_noStrctrdCells*CV->noVariables, "convectiveFlux Debug", __CALLING_FUNCTION__);
#endif

  // always allocate moving grid volume fluxes
  // but leave them zero if not moving grid
  zfsAlloc(m_cells->dxt, nDim, m_noStrctrdCells, "m_cells->dxt", F0, __CALLING_FUNCTION__);

  zfsAlloc(m_cells->area, nDim, m_noStrctrdCells, "m_cells->area", 100000.0, __CALLING_FUNCTION__);

  // allocate viscous flux computation variables
  ZFSInt noVars = (CV->noVariables-1)*m_noStrctrdCells;
  ZFSInt noVarsFluxes = (CV->noVariables-1 + m_rans)*m_noStrctrdCells;

  zfs_log << "Allocating fFlux with " << (CV->noVariables-1 + m_rans) << " variables for " << m_noStrctrdCells << " cells"  << endl;

  zfsAlloc(m_cells->eFlux, noVarsFluxes, "m_cells->eFlux", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_cells->fFlux, noVarsFluxes, "m_cells->fFlux", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_cells->gFlux, noVarsFluxes, "m_cells->gFlux", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_cells->uVWT, noVars, "m_cells->uVWT", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_QLeft, m_maxNoVariables, "m_QLeft", F0, __CALLING_FUNCTION__); //junoh
  zfsAlloc(m_QRight, m_maxNoVariables, "m_QRight", F0, __CALLING_FUNCTION__); //junoh

  zfsAlloc(m_cells->viscousFlux, nDim*m_noStrctrdCells, "m_cells->viscousFlux", 1000.0, __CALLING_FUNCTION__);
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::prepareNonBlockingCommunication(){
  for(ZFSId i=0; i<m_cmnctnFlag->noNghbrDomainsNormal; i++){
    if(m_cmnctnFlag->mpi_rcvRequestCells[i] != MPI_REQUEST_NULL) {
      MPI_Request_free(&m_cmnctnFlag->mpi_rcvRequestCells[i]);
    }
    MPI_Recv_init(m_cmnctnFlag->m_bufferCellsRcv[i], m_cmnctnFlag->m_noNghbrDomainCellBufferSizeRcv[i], MPI_DOUBLE,m_cmnctnFlag->m_rcvNghbrId[i], 11111, m_zfsStrctrdComm, &m_cmnctnFlag->mpi_rcvRequestCells[i]);
    if(m_cmnctnFlag->mpi_sndRequestCells[i] != MPI_REQUEST_NULL) {
      MPI_Request_free(&m_cmnctnFlag->mpi_sndRequestCells[i]);
    }
    MPI_Send_init(m_cmnctnFlag->m_bufferCellsSnd[i],m_cmnctnFlag->m_noNghbrDomainCellBufferSizeSnd[i], MPI_DOUBLE,  m_cmnctnFlag->m_sndNghbrId[i] ,11111, m_zfsStrctrdComm,&m_cmnctnFlag->mpi_sndRequestCells[i]);
  }

  MPI_Startall(m_cmnctnFlag->noNghbrDomainsNormal, m_cmnctnFlag->mpi_rcvRequestCells);
  MPI_Startall(m_cmnctnFlag->noNghbrDomainsNormal, m_cmnctnFlag->mpi_sndRequestCells);
  //see comment of chris in zfsfvblock.cpp if in trouble
  MPI_Waitall(m_cmnctnFlag->noNghbrDomainsNormal, m_cmnctnFlag->mpi_rcvRequestCells, MPI_STATUSES_IGNORE);
  MPI_Startall(m_cmnctnFlag->noNghbrDomainsNormal, m_cmnctnFlag->mpi_rcvRequestCells );
}


/*! \fn void ZFSStrctrdBlck<nDim>::setInputOutputProperties()
 * \brief Reads properties and initializes variables associated with input/output.
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::setInputOutputProperties() {
  TRACE();
  ZFSInt tmpFalse = 0; 

  m_outputIterationNumber = 0;
  m_outputFormat = ".hdf5";
  m_lastOutputTimeStep = -1;

  /*! \page propertyPage1
    \section solutionOutput
    <code>ZFSString ZFSStrctrdBlck::m_solutionOutput </code>\n
    default = <code>"./out"</code>\n \n
    Solution output folder.\n
    Keywords: <i>SOLUTION, IO, STRCTRD</i>
  */
  m_solutionOutput = "./out";
  m_solutionOutput = *(ZFSContext::getProperty("solutionOutput", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString(0));

  /*! \page propertyPage1
    \section dragOutputInterval
    <code>ZFSInt ZFSStrctrdBlck::m_dragOutputInterval </code>\n
    default = <code>"./out"</code>\n \n
    Interval in which auxDataFiles (containing drag and cf/cp etc.).\n
    should be written.\n
    Possible values are:\n
    <ul>
      <li>Integer >= 0</li>
    </ul>
    Keywords: <i>DRAG, IO, STRCTRD</i>
  */
  m_dragOutputInterval = 0;
  m_dragOutputInterval = *(ZFSContext::getProperty("dragOutputInterval", m_blockId, __CALLING_FUNCTION__, &m_dragOutputInterval)->asInt(0));

  if(m_dragOutputInterval > 0) {

    /*! \page propertyPage1
      \section auxOutputDir
      <code>ZFSString ZFSStrctrdBlck::m_auxOutputDir </code>\n
      default = <code> solutionOutput</code>\n \n
      Folder for auxData files.\n
      Keywords: <i>DRAG, IO, STRCTRD</i>
    */
    m_auxOutputDir = m_solutionOutput;
    if(ZFSContext::propertyExists("auxOutputDir", m_blockId )){
      m_auxOutputDir = *(ZFSContext::getProperty("auxOutputDir", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString(0));

      ZFSString comparator = "/";
      if(strcmp((m_auxOutputDir.substr(m_auxOutputDir.length()-1,1)).c_str(),comparator.c_str()) != 0) {
        m_auxOutputDir = m_auxOutputDir + "/";
      }
    }
  }

  /*! \page propertyPage1
    \section dragAsciiOutputInterval
    <code>ZFSInt ZFSStrctrdBlck::m_dragAsciiOutputInterval </code>\n
    default = <code>"./out"</code>\n \n
    Interval in which the integrated drag .\n
    should be written to an ASCII file.\n
    Possible values are:\n
    <ul>
      <li>Integer >= 0</li>
    </ul>
    Keywords: <i>DRAG, IO, STRCTRD</i>
  */
  m_dragAsciiOutputInterval = 0;
  m_dragAsciiOutputInterval = *(ZFSContext::getProperty("dragAsciiOutputInterval", m_blockId, __CALLING_FUNCTION__, &m_dragAsciiOutputInterval)->asInt(0));

  /*! \page propertyPage1
    \section outputOffset
    <code>ZFSInt ZFSStrctrdBlck::m_outputOffset </code>\n
    default = <code> 0 </code>\n \n
    Time step before which no output should be written.\n
    Possible values are:\n
    <ul>
      <li>Integer >= 0</li>
    </ul>
    Keywords: <i>IO, STRCTRD</i>
  */
  m_outputOffset= 0;
  m_outputOffset= *(ZFSContext::getProperty("outputOffset", m_blockId, __CALLING_FUNCTION__, &m_outputOffset)->asInt(0));

  /*! \page propertyPage1
    \section ignoreUID
    <code>ZFSInt ZFSStrctrdBlck::m_ignoreUID </code>\n
    default = <code> 0</code>\n \n
    Switch to override the UID check for
    restart files
    Possible values are:\n
    <ul>
      <li>0 = UID is checked</li>
      <li>1 = UID is not checked</li>
    </ul>
    Keywords: <i>RESTART, IO, STRCTRD</i>
  */
  m_ignoreUID = 0;
  if(ZFSContext::propertyExists("ignoreUID", m_blockId )){
    m_ignoreUID = *(ZFSContext::getProperty("ignoreUID", m_blockId, __CALLING_FUNCTION__, (ZFSInt*) &m_ignoreUID )->asInt(0));
    zfs_log << "WARNING!!!!!!!!!!!!!!: UID was not checked. Solution and grid might not fit together" << endl;
  }

  /*! \page propertyPage1
    \section restartFile
    <code>ZFSInt ZFSStrctrdBlck::m_restartFile </code>\n
    default = <code> 0</code>\n \n
    Possible values are:\n
    <ul>
      <li>0 = initial start</li>
      <li>1 = start from restart file</li>
    </ul>
    Keywords: <i>RESTART, IO, STRCTRD</i>
  */
  m_restart = 0;
  m_restart = *(ZFSContext::getProperty("restartFile", m_blockId, __CALLING_FUNCTION__, (ZFSInt*) &tmpFalse )->asInt(0));

  /*! \page propertyPage1
    \section useNonSpecifiedRestartFile
    <code>ZFSString ZFSStrctrdBlck::m_useNonSpecifiedRestartFile </code>\n
    default = <code> 0</code>\n \n
    Keywords: <i>RESTART, IO, STRCTRD</i>
  */
  m_useNonSpecifiedRestartFile = 0;
  m_useNonSpecifiedRestartFile = *(ZFSContext::getProperty("useNonSpecifiedRestartFile", m_blockId, __CALLING_FUNCTION__, &tmpFalse )->asInt(0));

  /*! \page propertyPage1
    \section changeMa
    <code>ZFSInt ZFSStrctrdBlck::m_changeMa </code>\n
    default = <code> 0</code>\n \n
    Specify whether the variables should be transformed to
    a changed Ma number
    Possible values are:\n
    <ul>
      <li>0 = no conversion</li>
      <li>1 = convert all variables to new Ma number</li>
    </ul>
    Keywords: <i>RESTART, IO, STRCTRD</i>
  */
  m_changeMa = false;
  if(m_restart){
    m_changeMa = *(ZFSContext::getProperty("changeMa", m_blockId, __CALLING_FUNCTION__,  &m_changeMa )->asInt(0));
  }

  /*! \page propertyPage1
    \section debugOutput
    <code>ZFSInt ZFSStrctrdBlck::m_debugOutput </code>\n
    default = <code> 0</code>\n \n
    Write out debug information as blockId and cellId to solution file.
    Possible values are:\n
    <ul>
      <li>0 = no debug output</li>
      <li>1 = write blockId, cellId</li>
    </ul>
    Keywords: <i>RESTART, IO, STRCTRD</i>
  */
  m_debugOutput = false;
  m_debugOutput = *(ZFSContext::getProperty("debugOutput", m_blockId, __CALLING_FUNCTION__, (ZFSInt*) &tmpFalse )->asInt(0));

  if(m_debugOutput){
    FQ->neededFQVariables[FQ->BLOCKID] = 1;
    FQ->neededFQVariables[FQ->CELLID] = 1;
  }

  /*! \page propertyPage1
    \section writeSlopes
    <code>ZFSInt ZFSStrctrdBlck::m_writeSlopes </code>\n
    default = <code> 0</code>\n \n
    Write out the slopes at each rungeKutta step for debug
    Possible values are:\n
    <ul>
      <li>0 = nothing written</li>
      <li>1 = writing all slopes</li>
    </ul>
    Keywords: <i>RESTART, IO, STRCTRD</i>
  */
  m_writeSlopes = false;
  m_writeSlopes = *(ZFSContext::getProperty("writeSlopes", m_blockId, __CALLING_FUNCTION__, (ZFSInt*) &tmpFalse )->asInt(0));

  if(m_writeSlopes) {
    FQ->neededFQVariables[FQ->SLOPEX] = 1;
    FQ->neededFQVariables[FQ->SLOPEY] = 1;
    FQ->neededFQVariables[FQ->SLOPEZ] = 1;
  }

  if(!m_restart){
    m_restartTimeStep = 0;
  } else {
    if(!m_useNonSpecifiedRestartFile){
      m_restartTimeStep = -1;
      if(domainId()==0){//faster version if many processors are used especially on juqueen
        stringstream restartFileName;

        /*! \page propertyPage1
          \section restartVariablesFileName
          <code>ZFSString ZFSStrctrdBlck::m_restartFile </code>\n
          default = <code>""</code>\n \n
          Name of the restart file.\n
          Keywords: <i>GRID, STRCTRD</i>
        */
        ZFSString restartFile = *(ZFSContext::getProperty("restartVariablesFileName", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL )->asString(0));
        restartFileName << outputDir()<<  restartFile;
        //open the file
        ZFSId restartId = io_openfile("hdf5" ,(restartFileName.str()).c_str(),"collective", MPI_COMM_SELF);
        io_read_iattribute1(restartId,"", "globalTimeStep", &m_restartTimeStep);
        io_closefile(restartId);
	
	//send data to other processors
        MPI_Bcast((void*)&m_restartTimeStep, 1, MPI_INT, 0, m_zfsStrctrdComm);
      }else{
        MPI_Bcast((void*)&m_restartTimeStep, 1, MPI_INT, 0, m_zfsStrctrdComm);
      }
      if(m_restartTimeStep == -1) {
        zfsTerm(1, __CALLING_FUNCTION__, "Could not read from restart file! Is a restart desired and does the file given in the property file exist?");
      }
    }
  }

  /*! \page propertyPage1
    \section savePartitionOutput
    <code>ZFSInt ZFSStrctrdBlck::m_savePartitionOutput </code>\n
    default = <code> 0</code>\n \n
    Save also partitioned solution file with current number of domains.\n
    Possible values are:\n
    <ul>
      <li>0 = no output</li>
      <li>1 = write partitioned file</li>
    </ul>
    Keywords: <i>RESTART, IO, STRCTRD</i>
  */
  m_savePartitionOutput = false;
  m_savePartitionOutput = *(ZFSContext::getProperty("savePartitionOutput", m_blockId, __CALLING_FUNCTION__, &tmpFalse )->asInt(0));

  /*! \page propertyPage1
    \section computeCf
    <code>ZFSInt ZFSStrctrdBlck::m_bCfCpCoeff </code>\n
    default = <code> 0</code>\n \n
    Compute and write skin friction and pressure coefficient.\n
    Possible values are:\n
    <ul>
      <li>0 = no output</li>
      <li>1 = write cf values</li>
    </ul>
    Keywords: <i>DRAG, IO, STRCTRD</i>
  */
  m_bCfCpCoeff=0;
  m_bCfCpCoeff = *(ZFSContext::getProperty("computeCfCp", m_blockId, __CALLING_FUNCTION__, &m_bCfCpCoeff)->asInt(0));
  if(m_bCfCpCoeff) zfs_log << "<<<<< Skin-friction and Pressure Coefficient Computation: ENABLED" << endl;
  
  /*! \page propertyPage1
    \section computePower
    <code>ZFSInt ZFSStrctrdBlck::m_bPower </code>\n
    default = <code> 0</code>\n \n
    Compute and write out the power spent for actuation.\n
    Possible values are:\n
    <ul>
      <li>0 = no output</li>
      <li>1 = write power values</li>
    </ul>
    Keywords: <i>DRAG, IO, POWER, STRCTRD</i>
  */
  m_bPower=0;
  m_bPower = *(ZFSContext::getProperty("computePower", m_blockId, __CALLING_FUNCTION__, &m_bPower)->asInt(0));
  if(m_bPower) zfs_log << "<<<<< Power Computation: ENABLED" << endl;


  /*! \page propertyPage1
    \section computeCl
    <code>ZFSInt ZFSStrctrdBlck::m_bCl </code>\n
    default = <code> 0</code>\n \n
    Compute and write lift coefficient.\n
    Possible values are:\n
    <ul>
      <li>0 = no output</li>
      <li>1 = write cl values</li>
    </ul>
    Keywords: <i>DRAG, IO, STRCTRD</i>
  */
  m_bCl=0;
  m_bCl = *(ZFSContext::getProperty("computeCl", m_blockId, __CALLING_FUNCTION__, &m_bCl)->asInt(0));
  if(m_bCl) zfs_log << "<<<<< Lift Computation: ENABLED" << endl;
  /*! \page propertyPage1
    \section computeCd
    <code>ZFSInt ZFSStrctrdBlck::m_bCd </code>\n
    default = <code> 0</code>\n \n
    Compute and write drag coefficient.\n
    Possible values are:\n
    <ul>
      <li>0 = no output</li>
      <li>1 = write cd values</li>
    </ul>
    Keywords: <i>DRAG, IO, STRCTRD</i>
  */
  m_bCd=0;
  m_bCd = *(ZFSContext::getProperty("computeCd", m_blockId, __CALLING_FUNCTION__, &m_bCd)->asInt(0));
  if(m_bCd) zfs_log << "<<<<< Drag Computation: ENABLED" << endl;
  /*! \page propertyPage1
    \section detailAuxData
    <code>ZFSInt ZFSStrctrdBlck::m_detailAuxData </code>\n
    default = <code> 0</code>\n \n
    Write additional information into auxData files (areas, coordiantes).\n
    Possible values are:\n
    <ul>
      <li>0 = no output</li>
      <li>1 = write additional info</li>
    </ul>
    Keywords: <i>DRAG, IO, STRCTRD</i>
  */
  m_detailAuxData=0;
  m_detailAuxData = *(ZFSContext::getProperty("detailAuxData", m_blockId, __CALLING_FUNCTION__, &m_detailAuxData)->asInt(0));

  /*! \page propertyPage1
    \section computeCpLineAverage
    <code>ZFSInt ZFSStrctrdBlck::m_bCpLineAveraging </code>\n
    default = <code> 0</code>\n \n
    Trigger the line averaging of the cp/cf value.\n
    Possible values are:\n
    <ul>
      <li>0 = no averaging</li>
      <li>1 = compute line average</li>
    </ul>
    Keywords: <i>DRAG, IO, STRCTRD</i>
  */
  m_bCpLineAveraging=0;
  m_bCpLineAveraging = *(ZFSContext::getProperty("computeCpLineAverage", m_blockId, __CALLING_FUNCTION__, &m_bCpLineAveraging)->asInt(0));

  /*! \page propertyPage1
    \section cpAveragingDir
    <code>ZFSInt ZFSStrctrdBlck::m_cpAveragingDir </code>\n
    default = <code> 0</code>\n \n
    Direction in which to compute the cp/cf line average\n
    Possible values are:\n
    <ul>
      <li>0 = i-direction</li>
      <li>1 = j-direction</li>
      <li>2 = k-direction</li>
    </ul>
    Keywords: <i>DRAG, IO, STRCTRD</i>
  */
  m_cpAveragingDir=0;
  m_cpAveragingDir = *(ZFSContext::getProperty("cpAveragingDir", m_blockId, __CALLING_FUNCTION__, &m_cpAveragingDir)->asInt(0));

  /*! \page propertyPage1
    \section auxDataCoordinateLimits
    <code>ZFSInt ZFSStrctrdBlck::m_auxDataCoordinateLimits </code>\n
    default = <code> 0</code>\n \n
    Trigger the limitation of the cp/cf computation region.\n
    Possible values are:\n
    <ul>
      <li>0 = no limits</li>
      <li>1 = read in limits</li>
    </ul>
    Keywords: <i>DRAG, IO, STRCTRD</i>
  */
  m_auxDataCoordinateLimits=0;
  m_auxDataCoordinateLimits = *(ZFSContext::getProperty("auxDataCoordinateLimits", m_blockId, __CALLING_FUNCTION__, &m_auxDataCoordinateLimits)->asInt(0));

  //verify the conditions for auxilary data computation, i.e., when dragoutput is higher than 1 then cp and cf should be enabled
  if(m_dragOutputInterval || m_dragAsciiOutputInterval){
    if(!m_bCfCpCoeff){
      m_dragOutputInterval=0;
      m_dragAsciiOutputInterval=0;
      zfs_log <<"\t   _            _" << endl;
      zfs_log <<"\t _/ \\ _______ / \\_ "<<endl;
      zfs_log <<"\t /  __/__   __\\__              \\ " << endl;
      zfs_log <<"\t\\_/____ \\ / ____\\_/" << endl; 
      zfs_log <<"\t //    \\___/    \\             \\ " << endl;
      zfs_log <<"\t | |   o _ _ o   | | " << endl;
      zfs_log <<"\t\\_\\___/ V \\___/_/ " << endl;
      zfs_log <<"\t  < `    | |    ' > " << endl;
      zfs_log <<"\t  \\__.  \\^/  .__/ " << endl;
      zfs_log <<"\t     >         < " << endl;
      zfs_log <<"\t     [_I_I_I_I_] " << endl;
      zfs_log <<"\t    /  /   \\  \\ " << endl;
      zfs_log <<"\t  _/  /     \\  \\_" << endl;
      zfs_log <<"\t /   <       >   \\ " << endl;
      zfs_log <<"\t \\_.  |     |  ._/ " << endl;
      zfs_log <<"\t  \\_/       \\_/ " << endl;
      zfs_log << endl;
      zfs_log << " ########## WARNING::::::::: m_dragOutputInterval or m_dragAsciiOutputInterval was set to zero as no computation of cp cf was enabled" << endl;
      
      if(domainId()==0) cout << "\033[0;31m########## WARNING:::::::::\033[0m m_dragOutputInterval or m_dragAsciiOutputInterval was set to zero as no computation of cp cf was enabled" << endl;
    }
  }


  if(m_auxDataCoordinateLimits) {
    m_auxDataLimits = NULL;
    ZFSInt noAuxDataLimits = 4;
    zfsAlloc(m_auxDataLimits, noAuxDataLimits, "m_auxDataLimits", F0, __CALLING_FUNCTION__);

    for(ZFSInt i=0; i<noAuxDataLimits; ++i){
      /*! \page propertyPage1
        \section auxDataLimits
        <code>ZFSFloat ZFSStrctrdBlck::m_auxDataLimits </code>\n
        default = <code> 0</code>\n \n
        Limiting coordinates for 2D rectangle\n
        for c_d computation. \n
        Possible values are:\n
        <ul>
        <li>0.0,2.0,0.0,2.0 = limits a rectangle of 2.0 x 2.0 </li>
        </ul>
        Keywords: <i>DRAG, IO, STRCTRD</i>
      */
      m_auxDataLimits[i]=-99999.9;
      m_auxDataLimits[i] =*(ZFSContext::getProperty("auxDataLimits", m_blockId, __CALLING_FUNCTION__,(ZFSFloat*)  &m_auxDataLimits[i] , 4 )->asFloat(i));
    }

    zfs_log  << "AuxData limited area"
             << ", lower x-limit: " << m_auxDataLimits[0]
             << ", upper x-limit: " << m_auxDataLimits[1]
             << ", lower z-limit: " << m_auxDataLimits[2]
             << ", upper z-limit: " << m_auxDataLimits[3] << endl;
  }

  /*! \page propertyPage1
    \section computeLambda2
    <code>ZFSInt ZFSStrctrdBlock::m_computeLamda2 </code>\n
    default = <code>0</code>\n \n
    this property triggers the lambda2 calculation. It will then be written out with
    the same frequency as used for saveOutput \n
    possible values are:
    <ul>
    <li>Non-negative integervalues: <0> := off; <1>:=on</li>
    </ul>
    Keywords: <i>TURBULENCE, VORTICES, STRCTRD</i>
  */
  m_computeLambda2 = 0 ;
  m_computeLambda2 = *(ZFSContext::getProperty("computeLambda2", m_blockId, __CALLING_FUNCTION__, &m_computeLambda2)->asInt(0));
  if(m_computeLambda2) {
    FQ->neededFQVariables[FQ->LAMBDA2] = 1;
    FQ->planeOutputFQVariables[FQ->LAMBDA2] = 1;
    FQ->boxOutputFQVariables[FQ->LAMBDA2] = 1;
  }

  /*! \page propertyPage1
    \section primitive output
    <code>ZFSInt ZFSStrctrdBlock::m_primitiveOutput </code>\n
    default = <code>0</code>\n \n
    Trigger primitive variable output\n
    possible values are:
    <ul>
    <li>Non-negative integervalues: <0> := off; <1>:=on</li>
    </ul>
    Keywords: <i>IO, STRCTRD</i>
  */
  m_primitiveOutput = 1;
  m_primitiveOutput = *(ZFSContext::getProperty("primitiveOutput", m_blockId, __CALLING_FUNCTION__,&m_primitiveOutput)->asInt(0));
  if(m_primitiveOutput){ zfs_log << "Primitive Output is enabled" << endl;}
  else{ zfs_log << "Conservative Output is enabled" << endl;}

  /*! \page propertyPage1
    \section vorticityOutput
    <code>ZFSInt ZFSStrctrdBlock::m_vorticityOutput </code>\n
    default = <code>0</code>\n \n
    Trigger vorticity computation and output\n
    possible values are:
    <ul>
    <li>Non-negative integervalues: <0> := off; <1>:=on</li>
    </ul>
    Keywords: <i>TURBULENCE, VORTICITY, STRCTRD</i>
  */
  m_vorticityOutput = 0;
  m_vorticityOutput = *(ZFSContext::getProperty("vorticityOutput", m_blockId, __CALLING_FUNCTION__,&tmpFalse)->asInt(0));

  /*! \page propertyPage1
    \section averageVorticity
    <code>ZFSInt ZFSStrctrdBlock::m_averageVorticity </code>\n
    default = <code>0</code>\n \n
    Average the vorticity in the postprocessing\n
    possible values are:
    <ul>
    <li>Non-negative integervalues: <0> := off; <1>:=on</li>
    </ul>
    Keywords: <i>POSTPROCESSING, VORTICITY, STRCTRD</i>
  */
  m_averageVorticity = 0;
  m_averageVorticity = *(ZFSContext::getProperty("pp_averageVorticity", m_blockId, __CALLING_FUNCTION__,&tmpFalse)->asInt(0));

  if(m_vorticityOutput || m_averageVorticity) {
    for(ZFSId v=0; v<nDim; v++) {FQ->neededFQVariables[FQ->VORTX+v] = 1;}
  }
  //if no vorticity output is activated but averaged vorticity is on then output
  //is deactivated for solution output 
  if(!m_vorticityOutput){
    for(ZFSId v=0; v<nDim; v++) {FQ->outputFQVariables[FQ->VORTX+v]=false;}
  }
      



  ///////////////////////////////////////////////////
  ///////////////// Variable Names //////////////////
  ///////////////////////////////////////////////////
  zfsAlloc(m_variableNames,CV->noVariables, "m_variableNames", __CALLING_FUNCTION__);
  switch(nDim){
  case 1:{
    m_variableNames[CV->RHO_U] = "rhoU";
    m_variableNames[CV->RHO_E] = "rhoE";
    m_variableNames[CV->RHO] =   "rho";
    break;
  }
  case 2:{
    m_variableNames[CV->RHO_U] = "rhoU";
    m_variableNames[CV->RHO_V] = "rhoV";
    m_variableNames[CV->RHO_E] = "rhoE";
    m_variableNames[CV->RHO] =   "rho";
    break;
  }
  case 3:{
    m_variableNames[CV->RHO_U] = "rhoU";
    m_variableNames[CV->RHO_V] = "rhoV";
    m_variableNames[CV->RHO_W] = "rhoW";
    m_variableNames[CV->RHO_E] = "rhoE";
    m_variableNames[CV->RHO] =   "rho";
    break;
  }
  default:{
      zfsTerm(1, __CALLING_FUNCTION__, "spatial dimension not implemented for m_variableNames");
    break;
  }
  }

  //fill up the rest of the variales for the species
    if(nDim+2<CV->noVariables){
      m_variableNames[nDim+2]="rhoZ";
      for(ZFSId i=nDim+2; i<CV->noVariables; ++i){
      stringstream number;
      number<< i-(nDim+2);
      ZFSString varName="rho" + number.str();
      m_variableNames[i]=varName;
    }
  }

    zfsAlloc(m_pvariableNames,m_maxNoVariables, "m_pvariableNames", __CALLING_FUNCTION__); //junoh
    switch(nDim){
  case 1:{
    m_pvariableNames[PV->U] = "u";
    m_pvariableNames[PV->P] = "p";
    m_pvariableNames[PV->RHO] = "rho";
    break;
  }
  case 2:{
    m_pvariableNames[PV->U] = "u";
    m_pvariableNames[PV->V] = "v";
    m_pvariableNames[PV->P] = "p";
    m_pvariableNames[PV->RHO] = "rho";
    break;
  }
  case 3:{
    m_pvariableNames[PV->U] = "u";
    m_pvariableNames[PV->V] = "v";
    m_pvariableNames[PV->W] = "w";
    m_pvariableNames[PV->P] = "p";
    m_pvariableNames[PV->RHO] = "rho";
    break;
  }
  default:{
      zfsTerm(1, __CALLING_FUNCTION__, "spatial dimension not implemented for m_variableNames");
    break;
  }
  }

  //fill up the rest of the variales for the species
    if(nDim+2<m_maxNoVariables){  //junoh
    for(ZFSId i=nDim+2; i<m_maxNoVariables; ++i){
      stringstream number;
      number<< i-(nDim+2);
      ZFSString varName="rans" + number.str();
      m_pvariableNames[i]=varName;
    }
  }

  /*! \page propertyPage1
    \section residualOutputInterval
    <code>ZFSInt ZFSStrctrdBlck::m_residualOutputInterval </code>\n
    default = <code> 0</code>\n \n
    Interval for the Residual computation
    Possible values are:\n
    <ul>
      <li>Integer >= 1</li>
    </ul>
    Keywords: <i>RESIDUAL, IO, STRCTRD</i>
  */
  m_residualOutputInterval = 1;
  m_residualOutputInterval = *(ZFSContext::getProperty("residualOutputInterval", m_blockId, __CALLING_FUNCTION__, &m_residualOutputInterval)->asInt(0));

  m_residualFileExist=false;


  /////////////////////////////////////////////////
  ///////////////// Plane Output //////////////////
  /////////////////////////////////////////////////

  m_planeBlock= NULL;
  m_planeOffset = NULL;
  m_planeNormal = NULL;
  m_planeWriteCoordinates = false;

  //planeOutput
  m_noPlaneOutput=0;

  /*! \page propertyPage1
    \section planeOutputInterval
    <code>ZFSInt ZFSStrctrdBlck::m_planeOutputInterval </code>\n
    default = <code> 0</code>\n \n
    Interval for output of plane solution files.\n
    Possible values are:\n
    <ul>
      <li>Integer >= 0</li>
    </ul>
    Keywords: <i>PLANES, IO, STRCTRD</i>
  */
  m_planeOutputInterval = 0;
  m_planeOutputInterval= *(ZFSContext::getProperty("planeOutputInterval", m_blockId, __CALLING_FUNCTION__, &m_planeOutputInterval)->asInt(0));
  if(m_planeOutputInterval>0){
    //find out the number of planes

    /*! \page propertyPage1
      \section planeBlock
      <code>ZFSInt ZFSStrctrdBlck::planeBlocks </code>\n
      default = <code> 0</code>\n \n
      Block in which the plane is containted.\n
      Possible values are:\n
      <ul>
      <li>Integer >= 0</li>
      </ul>
      Keywords: <i>PLANES, IO, STRCTRD</i>
    */
    ZFSProperty * planeBlocks  = ZFSContext::getProperty("planeBlock", m_blockId, __CALLING_FUNCTION__, (ZFSId*) NULL, 1);
    m_noPlaneOutput = (ZFSInt)planeBlocks->count(); //number of Planes to be written out

    /*! \page propertyPage1
      \section planeOffset
      <code>ZFSInt ZFSStrctrdBlck::planeOffsets </code>\n
      default = <code> 0</code>\n \n
      Offset of the plane in normal direction in.\n
      computational coordinates.\n
      Possible values are:\n
      <ul>
      <li>Integer >= 0</li>
      </ul>
      Keywords: <i>PLANES, IO, STRCTRD</i>
    */
    ZFSProperty * planeOffsets = ZFSContext::getProperty("planeOffset", m_blockId, __CALLING_FUNCTION__, (ZFSId*) NULL, 1);

    /*! \page propertyPage1
      \section planeNormal
      <code>ZFSInt ZFSStrctrdBlck::planeNormals </code>\n
      default = <code> 0</code>\n \n
      Normal of the plane in computationl direction.\n
      Possible values are:\n
      <ul>
      <li>0: Normal in K-direction</li>
      <li>1: Normal in J-direction</li>
      <li>2: Normal in I-direction</li>
      </ul>
      Keywords: <i>PLANES, IO, STRCTRD</i>
    */
    ZFSProperty * planeNormals = ZFSContext::getProperty("planeNormal", m_blockId, __CALLING_FUNCTION__, (ZFSId*) NULL, 1);

    /*! \page propertyPage1
      \section planeOutputDir
      <code>ZFSInt ZFSStrctrdBlck::planeOutputDir </code>\n
      default = <code> m_solutionOutput</code>\n \n
      Output folder for plane output files.\n
      Possible values are:\n
      <ul>
      <li>String containing path</li>
      </ul>
      Keywords: <i>PLANES, IO, STRCTRD</i>
    */
    m_planeOutputDir = m_solutionOutput;
    if(ZFSContext::propertyExists("planeOutputDir", m_blockId )){
      m_planeOutputDir = *(ZFSContext::getProperty("planeOutputDir", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString(0));

      ZFSString comparator = "/";
      if(strcmp((m_planeOutputDir.substr(m_planeOutputDir.length()-1,1)).c_str(),comparator.c_str()) != 0) {
        m_planeOutputDir = m_planeOutputDir + "/";
      }
    }

   /*! \page propertyPage1
      \section planeWriteCoordinates
      <code>ZFSInt ZFSStrctrdBlck::m_planeWriteCoordinates </code>\n
      default = <code>0</code>\n \n
      Write the cell-centered coordinates\n
      also into the plane files.\n
      Possible values are:\n
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>PLANES, IO, STRCTRD</i>
    */
    m_planeWriteCoordinates = *(ZFSContext::getProperty("planeWriteCoordinates", m_blockId, __CALLING_FUNCTION__, &tmpFalse )->asInt(0));

    if(m_noPlaneOutput!=planeOffsets->count()|| m_noPlaneOutput!=planeNormals->count()){
      zfsTerm(1, __CALLING_FUNCTION__, "The number of Entries for 'planeNormal' and 'planeOffset' and 'planeBlock' do not coincide!! Please check");
    }
    zfsAlloc(m_planeBlock, m_noPlaneOutput, "m_planeBlock", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_planeOffset, m_noPlaneOutput, "m_planeOffset", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_planeNormal, m_noPlaneOutput, "m_planeNormal", 0, __CALLING_FUNCTION__);

    //read in the values for the Blocks, offset and, normals

    for(ZFSInt i=0; i<m_noPlaneOutput; ++i){
      m_planeBlock[i]=-1;
      m_planeOffset[i]=-1;
      m_planeNormal[i]=-1;
      m_planeBlock[i] =*(ZFSContext::getProperty("planeBlock", m_blockId, __CALLING_FUNCTION__, &m_planeBlock[i] ,m_noPlaneOutput )->asInt(i));
      m_planeOffset[i] =*(ZFSContext::getProperty("planeOffset", m_blockId, __CALLING_FUNCTION__, &m_planeOffset[i] ,m_noPlaneOutput )->asInt(i));

      m_planeNormal[i] =*(ZFSContext::getProperty("planeNormal", m_blockId, __CALLING_FUNCTION__, &m_planeNormal[i] ,m_noPlaneOutput )->asInt(i));
    }
    //finished reading in the plane output properties
  }

  /////////////////////////////////////////////////
  ///////////////// Line Output //////////////////
  /////////////////////////////////////////////////

  m_lineStart= NULL;
  m_lineDelta= NULL;
  m_lineDelta2d= NULL;
  m_lineNoPoints= NULL;
  m_lineNoPoints2d= NULL;
  m_pointCoordinates= NULL;;
  m_hasPartnerGlobal= NULL;
  m_hasPartnerLocal= NULL;
  m_interpolatedVarsGlobal= NULL;
  m_interpolatedVarsLocal= NULL;
  m_noFieldPointsTotal= 0;
  //lineOutput
  m_noLineOutput=0;
  m_no2dLines=0;

  /*! \page propertyPage1
    \section lineOutputInterval
    <code>ZFSInt ZFSStrctrdBlck::m_lineOutputInterval </code>\n
    default = <code> 0</code>\n \n
    Interval of the line interpolation output.\n
    Possible values are:\n
    <ul>
    <li>Integer >= 0</li>
    </ul>
    Keywords: <i>LINES, INTERPOLATION, IO, STRCTRD</i>
  */
  m_lineOutputInterval = 0;
  m_lineOutputInterval = *(ZFSContext::getProperty("lineOutputInterval", m_blockId, __CALLING_FUNCTION__, &m_lineOutputInterval)->asInt(0));

  if(m_lineOutputInterval>0){

    /*! \page propertyPage1
      \section lineOutputDir
      <code>ZFSInt ZFSStrctrdBlck::m_lineOutputDir </code>\n
      default = <code> m_solutionOutput</code>\n \n
      Folder to write the line output files.\n
      Possible values are:\n
      <ul>
      <li>String with path</li>
      </ul>
      Keywords: <i>LINES, INTERPOLATION, IO, STRCTRD</i>
    */
    m_lineOutputDir = m_solutionOutput;
    if(ZFSContext::propertyExists("lineOutputDir", m_blockId )){
      m_lineOutputDir = *(ZFSContext::getProperty("lineOutputDir", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString(0));

      ZFSString comparator = "/";
      if(strcmp((m_lineOutputDir.substr(m_lineOutputDir.length()-1,1)).c_str(),comparator.c_str()) != 0) {
        m_lineOutputDir = m_lineOutputDir + "/";
      }
    }

    /*! \page propertyPage1
      \section lineStartX
      <code>ZFSInt ZFSStrctrdBlck::lineStartX </code>\n
      default = <code> 0</code>\n \n
      Point in x-dir to start line distribution.\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>LINES, INTERPOLATION, IO, STRCTRD</i>
    */
    ZFSProperty * lineStartX = ZFSContext::getProperty("lineStartX", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);

    /*! \page propertyPage1
      \section lineStartY
      <code>ZFSInt ZFSStrctrdBlck::lineStartY </code>\n
      default = <code> 0</code>\n \n
      Point in y-dir to start line distribution.\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>LINES, INTERPOLATION, IO, STRCTRD</i>
    */
    ZFSProperty * lineStartY = ZFSContext::getProperty("lineStartY", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);

    /*! \page propertyPage1
      \section lineStartZ
      <code>ZFSInt ZFSStrctrdBlck::lineStartZ </code>\n
      default = <code> 0</code>\n \n
      Point in z-dir to start line distribution.\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>LINES, INTERPOLATION, IO, STRCTRD</i>
    */
    ZFSProperty * lineStartZ = ZFSContext::getProperty("lineStartZ", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);

    /*! \page propertyPage1
      \section lineDeltaX
      <code>ZFSInt ZFSStrctrdBlck::lineDeltaX </code>\n
      default = <code> 0</code>\n \n
      The delta in x-dir between the line points.\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>LINES, INTERPOLATION, IO, STRCTRD</i>
    */
    ZFSProperty * lineDeltaX = ZFSContext::getProperty("lineDeltaX", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);

    /*! \page propertyPage1
      \section lineDeltaY
      <code>ZFSInt ZFSStrctrdBlck::lineDeltaY </code>\n
      default = <code> 0</code>\n \n
      The delta in y-dir between the line points.\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>LINES, INTERPOLATION, IO, STRCTRD</i>
    */
    ZFSProperty * lineDeltaY = ZFSContext::getProperty("lineDeltaY", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);

    /*! \page propertyPage1
      \section lineDeltaZ
      <code>ZFSInt ZFSStrctrdBlck::lineDeltaZ </code>\n
      default = <code> 0</code>\n \n
      The delta in z-dir between the line points.\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>LINES, INTERPOLATION, IO, STRCTRD</i>
    */
    ZFSProperty * lineDeltaZ = ZFSContext::getProperty("lineDeltaZ", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);

    /*! \page propertyPage1
      \section lineNoPoints
      <code>ZFSInt ZFSStrctrdBlck::lineNoPoints </code>\n
      default = <code> 0</code>\n \n
      Number of the points to distribute.\n
      Possible values are:\n
      <ul>
      <li>Integer > 0</li>
      </ul>
      Keywords: <i>LINES, INTERPOLATION, IO, STRCTRD</i>
    */
    ZFSProperty * lineNoPoints = ZFSContext::getProperty("lineNoPoints", m_blockId, __CALLING_FUNCTION__, (ZFSId*) NULL, 1);

    m_noLineOutput = lineStartX->count();

    if(m_noLineOutput!=lineStartX->count()|| m_noLineOutput!=lineStartY->count() ||
       m_noLineOutput!=lineStartZ->count()|| m_noLineOutput!=lineDeltaX->count() ||
       m_noLineOutput!=lineDeltaY->count()|| m_noLineOutput!=lineDeltaZ->count() ||
       m_noLineOutput!=lineNoPoints->count()){
      zfsTerm(1, __CALLING_FUNCTION__, "The number of Entries for 'lineStartX', 'lineStartY', 'lineStartZ', 'lineDeltaX', 'lineDeltaY' and 'lineDeltaZ' do not coincide!! Please check");
    }

    //set number of lines in second direction
    if(ZFSContext::propertyExists("lineDeltaX_2d", m_blockId )&&ZFSContext::propertyExists("lineDeltaY_2d", m_blockId )&&ZFSContext::propertyExists("lineDeltaZ_2d", m_blockId )&&ZFSContext::propertyExists("lineNoPoints_2d", m_blockId )){

      /*! \page propertyPage1
        \section lineDeltaX_2d
        <code>ZFSInt ZFSStrctrdBlck::lineDeltaX2d </code>\n
        default = <code> 0</code>\n \n
        The delta in x-dir for the second\n
        dimension if 2d field interpolation is desired.\n
        Possible values are:\n
        <ul>
        <li>Float</li>
        </ul>
        Keywords: <i>FIELDS, INTERPOLATION, IO, STRCTRD</i>
      */
      ZFSProperty * lineDeltaX2d = ZFSContext::getProperty("lineDeltaX_2d", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);

      /*! \page propertyPage1
        \section lineDeltaY_2d
        <code>ZFSInt ZFSStrctrdBlck::lineDeltaY2d </code>\n
        default = <code> 0</code>\n \n
        The delta in y-dir for the second\n
        dimension if 2d field interpolation is desired.\n
        Possible values are:\n
        <ul>
        <li>Float</li>
        </ul>
        Keywords: <i>FIELDS, INTERPOLATION, IO, STRCTRD</i>
      */
      ZFSProperty * lineDeltaY2d = ZFSContext::getProperty("lineDeltaY_2d", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);

      /*! \page propertyPage1
        \section lineDeltaZ_2d
        <code>ZFSInt ZFSStrctrdBlck::lineDeltaZ2d </code>\n
        default = <code> 0</code>\n \n
        The delta in z-dir for the second\n
        dimension if 2d field interpolation is desired.\n
        Possible values are:\n
        <ul>
        <li>Float</li>
        </ul>
        Keywords: <i>FIELDS, INTERPOLATION, IO, STRCTRD</i>
      */
      ZFSProperty * lineDeltaZ2d = ZFSContext::getProperty("lineDeltaZ_2d", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);

      /*! \page propertyPage1
        \section lineNoPoints_2d
        <code>ZFSInt ZFSStrctrdBlck::lineNoPoints2d </code>\n
        default = <code> 0</code>\n \n
        The number of points for the second\n
        dimension if 2d field interpolation is desired.\n
        Possible values are:\n
        <ul>
        <li>Integer > 0</li>
        </ul>
        Keywords: <i>FIELDS, INTERPOLATION, IO, STRCTRD</i>
      */
      ZFSProperty * lineNoPoints2d = ZFSContext::getProperty("lineNoPoints_2d", m_blockId, __CALLING_FUNCTION__, (ZFSId*) NULL, 1);

      m_no2dLines= lineNoPoints2d->count();

      if(m_no2dLines!=lineDeltaX2d->count()|| m_no2dLines!=lineDeltaY2d->count()|| m_no2dLines!=lineDeltaZ2d->count() || m_no2dLines!=lineNoPoints2d->count()){
        cout << "no2dLines: " << m_no2dLines << " lineDeltaY2d: " << lineDeltaY2d->count() << "  lineDeltaX2d: " << lineDeltaX2d->count() << "  lineDeltaZ2d: " << lineDeltaZ2d->count() << endl;
        zfsTerm(1, __CALLING_FUNCTION__, "The number of Entries for 'lineDeltaX_2d', 'lineDeltaY_2d' and 'lineDeltaZ_2d' do not coincide!! Please check");
      }

    } else{
      m_no2dLines=m_noLineOutput; //just to solve allocating problems
    }

    zfsAlloc(m_lineStart, nDim, m_noLineOutput, "m_lineStart", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_lineDelta, nDim, m_noLineOutput, "m_lineDelta", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_lineNoPoints, m_noLineOutput, "m_lineNoPoints", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_lineNoPoints2d, m_no2dLines, "m_lineNoPoints2d", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_lineDelta2d, nDim, m_no2dLines, "m_lineDelta2d", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_fieldOffset, m_noLineOutput, "m_fieldOffset", 0, __CALLING_FUNCTION__);

    // for every "first direction" line there is only one field although there are more possibilities to combine the "first and second direction" lines
    // therefore the fieldOffset is equally sized to the number of "first direction" lines (m_noLineOutput)

    //read in the values for the startpoints and directionVectors

    //first direction
    for(ZFSInt i=0; i<m_noLineOutput; ++i){
      //initialize with unrealistic values
      for(ZFSId dim = 0; dim < nDim; dim++) {
        m_lineStart[dim][i]=-99999.9;
        m_lineDelta[dim][i]=-99999.9;
      }
      m_lineNoPoints[i]=-1;

      m_lineStart[0][i] =*(ZFSContext::getProperty("lineStartX", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) &m_lineStart[0][i] ,m_noLineOutput )->asFloat(i));
      m_lineStart[1][i] =*(ZFSContext::getProperty("lineStartY", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) &m_lineStart[1][i] ,m_noLineOutput )->asFloat(i));
      m_lineStart[2][i] =*(ZFSContext::getProperty("lineStartZ", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) &m_lineStart[2][i] ,m_noLineOutput )->asFloat(i));
      m_lineDelta[0][i] =*(ZFSContext::getProperty("lineDeltaX", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) &m_lineDelta[0][i] ,m_noLineOutput )->asFloat(i));
      m_lineDelta[1][i] =*(ZFSContext::getProperty("lineDeltaY", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) &m_lineDelta[1][i] ,m_noLineOutput )->asFloat(i));
      m_lineDelta[2][i] =*(ZFSContext::getProperty("lineDeltaZ", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) &m_lineDelta[2][i] ,m_noLineOutput )->asFloat(i));
      m_lineNoPoints[i] =*(ZFSContext::getProperty("lineNoPoints", m_blockId, __CALLING_FUNCTION__, (ZFSInt*) &m_lineNoPoints[i] ,m_noLineOutput )->asInt(i));
    }

    m_fieldInterpolation = false;
    //second direction
    if(ZFSContext::propertyExists("lineDeltaX_2d", m_blockId )&&ZFSContext::propertyExists("lineDeltaY_2d", m_blockId )&&ZFSContext::propertyExists("lineDeltaZ_2d", m_blockId )&&ZFSContext::propertyExists("lineNoPoints_2d", m_blockId )){

      m_fieldInterpolation = true;

      for(ZFSInt i=0; i<m_no2dLines; ++i){
        //initialize with unrealistic values
        for(ZFSId dim = 0; dim < nDim; dim++) {
          m_lineDelta2d[dim][i]=-99999.9;
        }
        m_lineNoPoints2d[i]=-1;

        m_lineDelta2d[0][i] =*(ZFSContext::getProperty("lineDeltaX_2d", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) &m_lineDelta2d[0][i] ,m_no2dLines )->asFloat(i));
        m_lineDelta2d[1][i] =*(ZFSContext::getProperty("lineDeltaY_2d", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) &m_lineDelta2d[1][i] ,m_no2dLines )->asFloat(i));
        m_lineDelta2d[2][i] =*(ZFSContext::getProperty("lineDeltaZ_2d", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) &m_lineDelta2d[2][i] ,m_no2dLines )->asFloat(i));
        m_lineNoPoints2d[i] =*(ZFSContext::getProperty("lineNoPoints_2d", m_blockId, __CALLING_FUNCTION__, (ZFSInt*) &m_lineNoPoints2d[i] ,m_no2dLines )->asInt(i));
      }

    } else{
      for(ZFSInt i=0; i<m_no2dLines; i++){
        m_lineNoPoints2d[i]=1;
      }
      for(ZFSInt fieldId=0; fieldId<m_noLineOutput; fieldId++){
        for(ZFSInt dim=0; dim< nDim; dim++){
          m_lineDelta2d[dim][fieldId]=0;
        }
      }
    }

    for(ZFSInt i=0; i<m_noLineOutput; i++){
      m_noFieldPointsTotal += m_lineNoPoints[i]*m_lineNoPoints2d[i];
    }

    zfsAlloc(m_pointCoordinates, nDim, m_noFieldPointsTotal, "m_pointCoordinates",F0,__CALLING_FUNCTION__);
    zfsAlloc(m_hasPartnerGlobal, m_noFieldPointsTotal, "m_hasPartnerGlobal",0,__CALLING_FUNCTION__);
    zfsAlloc(m_hasPartnerLocal, m_noFieldPointsTotal, "m_hasPartnerLocal",0,__CALLING_FUNCTION__);
    zfsAlloc(m_interpolatedVarsLocal, CV->noVariables, m_noFieldPointsTotal, "m_interpolatedVarsLocal",F0, __CALLING_FUNCTION__);
    zfsAlloc(m_interpolatedVarsGlobal, CV->noVariables, m_noFieldPointsTotal, "m_interpolatedVarsGlobal",F0, __CALLING_FUNCTION__);

    //all fields are in ONE array (m_pointCoordinates)
    ZFSId offset = 0;

    //putting startpoints at the right place and generating the field
    for(ZFSInt fieldId=0; fieldId<m_noLineOutput; fieldId++){
      m_fieldOffset[fieldId] = offset;

      for(ZFSInt pointId2d=0; pointId2d<m_lineNoPoints2d[fieldId]; pointId2d++){
        for(ZFSInt pointId=0; pointId<m_lineNoPoints[fieldId]; pointId++){
          for(ZFSInt dim=0; dim < nDim; dim++) {
            m_pointCoordinates[dim][offset]=m_lineStart[dim][fieldId]+pointId*m_lineDelta[dim][fieldId]+pointId2d*m_lineDelta2d[dim][fieldId];
          } offset++;
        }
      }
    }

    //Finished reading in the line output properties and generating the lines
  }


  /////////////////////////////////////////////////
  ///////////////// Box Output ////////////////////
  /////////////////////////////////////////////////

  /*! \page propertyPage1
    \section boxOutputInterval
    <code>ZFSInt ZFSStrctrdBlck::boxOutputInterval </code>\n
    default = <code> 0</code>\n \n
    Interval to write out the box output files.\n
    Possible values are:\n
    <ul>
    <li>Integer >= 0</li>
    </ul>
    Keywords: <i>BOXES, INTERPOLATION, IO, STRCTRD</i>
  */
  m_boxOutputInterval = 0;
  m_boxOutputInterval= *(ZFSContext::getProperty("boxOutputInterval", m_blockId, __CALLING_FUNCTION__, &m_boxOutputInterval)->asInt(0));

  if(m_boxOutputInterval>0) {

    /*! \page propertyPage1
      \section boxBlock
      <code>ZFSInt ZFSStrctrdBlck::boxBlocks </code>\n
      default = <code> 0</code>\n \n
      Blocks in which the box is contained.\n
      Possible values are:\n
      <ul>
      <li>Integer >= 0</li>
      </ul>
      Keywords: <i>BOXES, INTERPOLATION, IO, STRCTRD</i>
    */
    ZFSProperty * boxBlocks  = ZFSContext::getProperty("boxBlock", m_blockId, __CALLING_FUNCTION__, (ZFSId*) NULL, 1);
    m_noBoxOutput = (ZFSInt)boxBlocks->count(); //number of Boxes to be written out

    /*! \page propertyPage1
      \section boxWriteCoordinates
      <code>ZFSInt ZFSStrctrdBlck::m_boxWriteCoordinates </code>\n
      default = <code> 0</code>\n \n
      Write cell-center coordinates into the boxes.\n
      Possible values are:\n
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>BOXES, INTERPOLATION, IO, STRCTRD</i>
    */
    m_boxWriteCoordinates = 0;
    m_boxWriteCoordinates = *(ZFSContext::getProperty("boxWriteCoordinates", m_blockId, __CALLING_FUNCTION__, &tmpFalse )->asInt(0));

    /*! \page propertyPage1
      \section boxOutputDir
      <code>ZFSInt ZFSStrctrdBlck::m_boxOutputDir </code>\n
      default = <code> m_solutionOutput</code>\n \n
      Output folder to write the box output files.\n
      Possible values are:\n
      <ul>
      <li>String containing path</li>
      </ul>
      Keywords: <i>BOXES, INTERPOLATION, IO, STRCTRD</i>
    */
    m_boxOutputDir = m_solutionOutput;
    if(ZFSContext::propertyExists("boxOutputDir", m_blockId )){
      m_boxOutputDir = *(ZFSContext::getProperty("boxOutputDir", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString(0));

      ZFSString comparator = "/";
      if(strcmp((m_boxOutputDir.substr(m_boxOutputDir.length()-1,1)).c_str(),comparator.c_str()) != 0) {
        m_boxOutputDir = m_boxOutputDir + "/";
      }
    }

    zfsAlloc(m_boxBlock, m_noBoxOutput, "m_boxBlock", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_boxOffset, m_noBoxOutput, nDim, "m_boxOffset", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_boxSize, m_noBoxOutput, nDim, "m_boxSize", 0, __CALLING_FUNCTION__);

    for(ZFSInt i=0; i<m_noBoxOutput; ++i){
      m_boxBlock[i]=-1;

      for(ZFSId dim=0; dim<nDim; dim++) {
        m_boxOffset[i][dim]=-1;
        m_boxSize[i][dim]=-1;
      }

      m_boxBlock[i] =*(ZFSContext::getProperty("boxBlock", m_blockId, __CALLING_FUNCTION__, &m_boxBlock[i] ,m_noBoxOutput )->asInt(i));

      /*! \page propertyPage1
        \section boxOffsetK
        <code>ZFSInt ZFSStrctrdBlck::m_boxOffsetK </code>\n
        default = <code> 0</code>\n \n
        Offset of the box in K-dir.\n
        Possible values are:\n
        <ul>
        <li>Integer >= 0</li>
        </ul>
        Keywords: <i>BOXES, INTERPOLATION, IO, STRCTRD</i>
      */
      m_boxOffset[i][0] =*(ZFSContext::getProperty("boxOffsetK", m_blockId, __CALLING_FUNCTION__, &m_boxOffset[i][0] ,m_noBoxOutput )->asInt(i));

      /*! \page propertyPage1
        \section boxOffsetJ
        <code>ZFSInt ZFSStrctrdBlck::m_boxOffsetJ </code>\n
        default = <code> 0</code>\n \n
        Offset of the box in J-dir.\n
        Possible values are:\n
        <ul>
        <li>Integer >= 0</li>
        </ul>
        Keywords: <i>BOXES, INTERPOLATION, IO, STRCTRD</i>
      */
      m_boxOffset[i][1] =*(ZFSContext::getProperty("boxOffsetJ", m_blockId, __CALLING_FUNCTION__, &m_boxOffset[i][1] ,m_noBoxOutput )->asInt(i));

      /*! \page propertyPage1
        \section boxOffsetI
        <code>ZFSInt ZFSStrctrdBlck::m_boxOffsetI </code>\n
        default = <code> 0</code>\n \n
        Offset of the box in I-dir.\n
        Possible values are:\n
        <ul>
        <li>Integer >= 0</li>
        </ul>
        Keywords: <i>BOXES, INTERPOLATION, IO, STRCTRD</i>
      */
      m_boxOffset[i][2] =*(ZFSContext::getProperty("boxOffsetI", m_blockId, __CALLING_FUNCTION__, &m_boxOffset[i][2] ,m_noBoxOutput )->asInt(i));

      /*! \page propertyPage1
        \section boxSizeK
        <code>ZFSInt ZFSStrctrdBlck::m_boxSizeK </code>\n
        default = <code> 0</code>\n \n
        Size of the box in K-dir.\n
        Possible values are:\n
        <ul>
        <li>Integer > 0</li>
        </ul>
        Keywords: <i>BOXES, INTERPOLATION, IO, STRCTRD</i>
      */
      m_boxSize[i][0] =*(ZFSContext::getProperty("boxSizeK", m_blockId, __CALLING_FUNCTION__, &m_boxSize[i][0] ,m_noBoxOutput )->asInt(i));

      /*! \page propertyPage1
        \section boxSizeJ
        <code>ZFSInt ZFSStrctrdBlck::m_boxSizeJ </code>\n
        default = <code> 0</code>\n \n
        Size of the box in J-dir.\n
        Possible values are:\n
        <ul>
        <li>Integer > 0</li>
        </ul>
        Keywords: <i>BOXES, INTERPOLATION, IO, STRCTRD</i>
      */
      m_boxSize[i][1] =*(ZFSContext::getProperty("boxSizeJ", m_blockId, __CALLING_FUNCTION__, &m_boxSize[i][1] ,m_noBoxOutput )->asInt(i));

      /*! \page propertyPage1
        \section boxSizeI
        <code>ZFSInt ZFSStrctrdBlck::m_boxSizeI </code>\n
        default = <code> 0</code>\n \n
        Size of the box in I-dir.\n
        Possible values are:\n
        <ul>
        <li>Integer > 0</li>
        </ul>
        Keywords: <i>BOXES, INTERPOLATION, IO, STRCTRD</i>
      */
      m_boxSize[i][2] =*(ZFSContext::getProperty("boxSizeI", m_blockId, __CALLING_FUNCTION__, &m_boxSize[i][2] ,m_noBoxOutput )->asInt(i));
    }
  }


  /////////////////////////////////////////////////
  ////////// Convective Unit Output ///////////////
  /////////////////////////////////////////////////

  /*! \page propertyPage1
    \section useConvectiveUnitWrite
    <code>ZFSInt ZFSStrctrdBlck::m_useConvectiveUnitWrite </code>\n
    default = <code> 0</code>\n \n
    Solution interval controlled by fraction or multiple of convective unit\n
    instead of globalTimeStep.\n
    Possible values are:\n
    <ul>
      <li>0: off</li>
      <li>1: on</li>
    </ul>
    Keywords: <i>IO, STRCTRD</i>
  */
  m_useConvectiveUnitWrite = false;
  m_useConvectiveUnitWrite = *(ZFSContext::getProperty("useConvectiveUnitWrite", m_blockId, __CALLING_FUNCTION__, &tmpFalse )->asInt(0));

  if(m_useConvectiveUnitWrite) {
    /*! \page propertyPage1
      \section convectiveUnitInterval
      <code>ZFSInt ZFSStrctrdBlck::convectiveUnitInterval </code>\n
      default = <code> 0</code>\n \n
      Solution interval controlled by fraction or multiple of convective unit\n
      instead of globalTimeStep.\n
      Possible values are:\n
      <ul>
      <li>floating point number > 0.0 </li>
      </ul>
      Keywords: <i>IO, STRCTRD</i>
    */
    m_convectiveUnitInterval = 1.0;
    m_convectiveUnitInterval = *(ZFSContext::getProperty("convectiveUnitInterval", m_blockId, __CALLING_FUNCTION__,&m_convectiveUnitInterval )->asFloat(0));

    m_noConvectiveOutputs = 0;

    /*! \page propertyPage1
      \section sampleSolutionFiles
      <code>ZFSInt ZFSStrctrdBlck::sampleSolutionFiles </code>\n
      default = <code> 0</code>\n \n
      Trigger the output of solution files controlled\n
      by the convectiveUnitInterval.\n
      Possible values are:\n
      <ul>
      <li>0: off</li>
      <li>1: on</li>
      </ul>
      Keywords: <i>IO, STRCTRD</i>
    */
    m_sampleSolutionFiles = false;
    m_sampleSolutionFiles = *(ZFSContext::getProperty("sampleSolutionFiles", m_blockId, __CALLING_FUNCTION__, &tmpFalse )->asInt(0));
  }

  /////////////////////////////////////////////////
  ////////// Interpolation restart  ///////////////
  /////////////////////////////////////////////////
  /*! \page propertyPage1
    \section restartInterpolation
    <code>ZFSInt ZFSStrctrdBlck::restartInterpolation </code>\n
    default = <code> 0</code>\n \n
    Restart the computation with an interpolated field\n
    from a given donorVars/Grid.\n
    Possible values are:\n
    <ul>
    <li>0: off</li>
    <li>1: on</li>
    </ul>
    Keywords: <i>INTERPOLATION, IO, STRCTRD</i>
  */
  m_restartInterpolation = *(ZFSContext::getProperty("restartInterpolation", m_blockId, __CALLING_FUNCTION__, (ZFSInt*) &tmpFalse )->asInt(0));

  //set the correct function pointers for the proper output 
  if(m_primitiveOutput){
    saveSolution=[&](ZFSBool a) {this->saveOutputSolution<true>(a);};
    savePartitions=[&]() {this->saveOutputPartitions<true>();};
    if(nDim==3){
      savePlanes = [&] { this->saveOutputPlanes<true>(); };
      saveBoxes=[&]{ this->saveOutputBoxes<true>();};
    }else{
      savePlanes = nullptr;
      saveBoxes = nullptr;
    }
  }else{
    saveSolution=[&](ZFSBool a) {this->saveOutputSolution<false>(a);};
    savePartitions=[&]() {this->saveOutputPartitions<false>();};
    if(nDim==3){
      savePlanes = [&] { this->saveOutputPlanes<false>(); };
      saveBoxes=[&]{ this->saveOutputBoxes<false>();};
    }else{
      savePlanes = nullptr;
      saveBoxes = nullptr;
    }
  }

  //correct output if not wanted
  if(!m_boxOutputInterval){saveBoxes=nullptr;}
  if(!m_planeOutputInterval){savePlanes=nullptr;}

}

/* \brief Initializes the properties of the FvBlock.
// \author Pascal Meysonnat
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::initializeStrctrdBlck(bool* propertiesGroups){
  TRACE();

  (void) propertiesGroups;

  /*! \page propertyPage1
    \section noSpecies
    <code>ZFSInt ZFSStrctrdBlck::m_noSpecies </code>\n
    default = <code> 0</code>\n 
    Number of species for future species computation.\n
    Possible values are:\n
    <ul>
    <li>0: off</li>
    <li>1: on</li>
    </ul>
    Keywords: <i>SPECIES, STRCTRD</i>
  */
  m_noSpecies=0;
  m_noSpecies=*(ZFSContext::getProperty("noSpecies", m_blockId, __CALLING_FUNCTION__, &m_noSpecies)->asInt(0));

  FQ = new ZFSStrctrdFQVariables();

  // allocate the array for counting the needed fq fields
  zfsAlloc(FQ->neededFQVariables, FQ->maxNoFQVariables, "FQ->neededFQVariables", 0, __CALLING_FUNCTION__);
  zfsAlloc(FQ->outputFQVariables, FQ->maxNoFQVariables, "FQ->outputFQVariables", true, __CALLING_FUNCTION__);
  zfsAlloc(FQ->planeOutputFQVariables, FQ->maxNoFQVariables, "FQ->planeOutputFQVariables", false, __CALLING_FUNCTION__);
  zfsAlloc(FQ->boxOutputFQVariables, FQ->maxNoFQVariables, "FQ->boxOutputFQVariables", false, __CALLING_FUNCTION__);

  m_timeStep = 0;
  m_periodicConnection= 0;
}

/// \fn void ZFSStrctrdBlck<nDim>:: setTestcaseProperties()
/// \brief Reads and initializes properties associated with the Testcase.
///
/// \author Pascal Meysonnat
///
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>:: setTestcaseProperties(){
  ZFSInt tmpFalse = 0;

  /*! \page propertyPage1
    \section referenceLength
    <code>ZFSFloat ZFSStrctrdBlck::m_Pr </code>\n
    default = <code>1.0</code>\n \n
    WARNING: Do NOT use any value different than 1.0 - The correct implementation of this is not checked,
    so it probably will not do what you think it does/should do. Don't use it unless you REALLY know what you are doing.
    Reference Length L - The length = 1.0 of the grid is scaled with L.
    Possible values are:
    <ul>
    <li>1.0 +- eps</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_referenceLength = F1;
  m_referenceLength = *(ZFSContext::getProperty("referenceLength", m_blockId, __CALLING_FUNCTION__,&m_referenceLength)->asFloat(0));

  if( fabs(m_referenceLength - F1) > m_eps){
    zfs_log << "WARNING: referenceLength != 1.0. The correct implementation of this is not checked. Don't use it unless you REALLY know what you are doing." << endl;
  }

  /*! \page propertyPage1
    \section physicalReferenceLength
    <code>ZFSFloat ZFSStrctrdBlck::m_Pr </code>\n
    default = <code>1.0</code>\n \n
    Physical Reference Length L - TODO !
    Possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_physicalReferenceLength = F1;
  m_physicalReferenceLength = *(ZFSContext::getProperty("physicalReferenceLength", m_blockId, __CALLING_FUNCTION__,&m_physicalReferenceLength )->asFloat(0));

  /*! \page propertyPage1
    \section Re
    <code>ZFSFloat ZFSStrctrdBlck::m_Re </code>\n
    default = <code>no default value</code>\n \n
    Reynolds number is defined with your infinity variables. \n
    In the code the Reynolds number is nondimensionalized to a Reynolds number based on the stagnation variables a_0, mu_0, rho_0 \n
    \f$ Re_{0} = Re_{\infty}  \frac{\mu_{\infty}}{\rho_{\infty} Ma \sqrt{T_{\infty}} } = \frac{\rho_0 a_0 l}{\mu_{0}}\f$:
    <ul>
    <li> \f$ mu_{\infty} \f$,  \f$ mu_{0} \f$ - viscosity  by the infinity, stagnation temperature </li>
    <li> \f$ Ma \f$ is the mach number </li>
    <li> \f$ T_{\infty} \f$ is the infinity temperature (free stream temperature)</li>
    </ul>
    possible values are:
    <ul>
    <li>Non-negative floating point values of the order of 0.1</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_Re = *(ZFSContext::getProperty("Re", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0))/m_referenceLength;

  /*! \page propertyPage1
    \section Pr
    <code>ZFSFloat ZFSStrctrdBlck::m_Pr </code>\n
    default = <code>0.72</code>\n \n
    Prandtl number  - non-dimensionalized with stagnant flow conditions
    \f$ \mu_{0},  \lambda_{0}, c_{p} \f$:
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_Pr = 0.72;
  m_Pr = *(ZFSContext::getProperty("Pr", m_blockId, __CALLING_FUNCTION__, &m_Pr)->asFloat(0));
  m_rPr = 1./m_Pr;

  /*! \page propertyPage1
    \section ReTau
    <code>ZFSFloat ZFSStrctrdBlck::m_ReTau </code>\n
    default = <code>no default value</code>\n \n
    ReTau value for certain flows as channel or pipe flow\n
    to control the pressure gradient.\n
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_ReTau = *(ZFSContext::getProperty("ReTau", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0))/m_referenceLength;

  /*! \page propertyPage1
    \section Ma
    <code>ZFSFloat ZFSStrctrdBlck::m_Ma </code>\n
    default = <code>no default value</code>\n \n
    Mach's number - \f$ M_{\infty} = \frac{u_\infty}{a_\infty} \f$:
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_Ma = *(ZFSContext::getProperty("Ma", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));

  /*! \page propertyPage1
    \section angle
    <code>ZFSFloat* ZFSStrctrdBlck::m_angle </code>\n
    default = <code>no default value</code>\n \n
    m_angle[nDim] - Angles of rotation around the z, and y axes.
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>FINITE_VOLUME,BOUNDARY CONDITION</i>
  */
  zfsAlloc( m_angle, nDim, "m_angle", F0, __CALLING_FUNCTION__ );
  if(ZFSContext::propertyExists("angle", m_blockId )){
    for ( int i = 0; i < (nDim - 1); i++ ) {
      m_angle[i] = *(ZFSContext::getProperty("angle", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL, (nDim-1) )->asFloat(i));
      m_angle[i] *= PI / 180.0;
    }
  }

  /*! \page propertyPage1
    \section gamma
    <code>ZFSFloat ZFSStrctrdBlck::m_gamma </code>\n
    default = <code>1.4</code>\n \n
    Ratio of specific heats - \f$ \gamma = c_p / c_v \f$
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_gamma = 1.4;
  m_gamma = *(ZFSContext::getProperty("gamma", m_blockId, __CALLING_FUNCTION__, &m_gamma )->asFloat(0));

  m_gammaMinusOne = m_gamma-F1;
  m_fgammaMinusOne = F1/(m_gammaMinusOne);

  /*! \page propertyPage1
    \section initialCondition
    <code>ZFSId ZFSStrctrdBlck::m_initialCondition </code>\n
    default = <code>no default value</code>\n \n
    Selects the initial condition.
    possible values are:
    <ul>
    <li>See the initialCondition() function of the corresponding block</li>
    </ul>
    Keywords: <i>INITIAL_CONDITION, FINITE_VOLUME</i>
  */
  m_initialCondition = *(ZFSContext::getProperty("initialCondition", m_blockId, __CALLING_FUNCTION__,(ZFSInt*) NULL)->asInt(0));


  m_channelHeight=-F1;
  m_channelWidth =-F1;
  m_channelLength=-F1;
  m_channelInflowPlaneCoordinate=-111111.1111111;
  m_channelC1=5.0;
  m_channelC2=-3.05;
  m_channelC3= 2.5;
  m_channelC4= 5.5;
  switch(m_initialCondition){
  case 1233:
  case 1234:{//we are dealing with channel flows

    /*! \page propertyPage1
      \section channelHeight
      <code>ZFSInt ZFSStrctrdBlck::m_channelHeigth </code>\n
      default = <code> 1.0 </code>\n \n
      Height of the half of the channel, necessary\n
      to compute correct pressure gradient\n
      Possible values are:\n
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelHeight = *(ZFSContext::getProperty("channelHeight", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));

    /*! \page propertyPage1
      \section channelWidth
      <code>ZFSInt ZFSStrctrdBlck::m_channelWidth </code>\n
      default = <code> 1.0 </code>\n \n
      Width of the channel, necessary\n
      to compute correct pressure gradient\n
      Possible values are:\n
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelWidth = *(ZFSContext::getProperty("channelWidth", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));

    /*! \page propertyPage1
      \section channelLength
      <code>ZFSInt ZFSStrctrdBlck::m_channelLength </code>\n
      default = <code> 1.0 </code>\n \n
      Length of the channel, necessary\n
      to compute correct pressure gradient\n
      Possible values are:\n
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelLength = *(ZFSContext::getProperty("channelLength", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));

    /*! \page propertyPage1
      \section channelInflowCoordinate
      <code>ZFSInt ZFSStrctrdBlck::m_channelInflowCoordinate </code>\n
      default = <code> 1.0 </code>\n \n
      Coordinate of the channel inflow plane.\n
      Possible values are:\n
      <ul>
      <li>Floating point</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelInflowPlaneCoordinate = *(ZFSContext::getProperty("channelInflowCoordinate", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));

    /*! \page propertyPage1
      \section loglawC1
      <code>ZFSInt ZFSStrctrdBlck::m_channelC1 </code>\n
      default = <code> 5.0 </code>\n \n
      First parameter for log-law for channel initialization\n
      Possible values are:\n
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelC1 = *(ZFSContext::getProperty("loglawC1", m_blockId, __CALLING_FUNCTION__, &m_channelC1)->asFloat(0));

    /*! \page propertyPage1
      \section loglawC2
      <code>ZFSInt ZFSStrctrdBlck::m_channelC2 </code>\n
      default = <code> -3.05 </code>\n \n
      Second parameter for log-law for channel initialization\n
      Possible values are:\n
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelC2 = *(ZFSContext::getProperty("loglawC2", m_blockId, __CALLING_FUNCTION__, &m_channelC1)->asFloat(0));

    /*! \page propertyPage1
      \section loglawC3
      <code>ZFSInt ZFSStrctrdBlck::m_channelC3 </code>\n
      default = <code> 2.5 </code>\n \n
      Third parameter for log-law for channel initialization\n
      Possible values are:\n
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelC3 = *(ZFSContext::getProperty("loglawC3", m_blockId, __CALLING_FUNCTION__, &m_channelC1)->asFloat(0));

    /*! \page propertyPage1
      \section loglawC4
      <code>ZFSInt ZFSStrctrdBlck::m_channelC4 </code>\n
      default = <code> 5.5 </code>\n \n
      Fourth parameter for log-law for channel initialization\n
      Possible values are:\n
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelC4 = *(ZFSContext::getProperty("loglawC4", m_blockId, __CALLING_FUNCTION__, &m_channelC1)->asFloat(0));
    zfs_log << "============= Channel flow activated =============" << endl;
    zfs_log << "-> channelHeight: " << m_channelHeight << endl;
    zfs_log << "-> channelWidth:  " << m_channelWidth  << endl;
    zfs_log << "-> channelLength: " << m_channelLength << endl;
    zfs_log << "-> channelInflowPlaneCoordinate: " << m_channelInflowPlaneCoordinate << endl;
    zfs_log << "-> Log law properties: " << endl;
    zfs_log << "--> C1: " << m_channelC1 << endl;
    zfs_log << "--> C2: " << m_channelC2 << endl;
    zfs_log << "--> C3: " << m_channelC3 << endl;
    zfs_log << "--> C4: " << m_channelC4 << endl;
    zfs_log << "============= Channel flow summary finished =============" << endl;
    break;
  }
  case 1236:{//we are dealing with pipe flows
    //comment:
    //since channel and pipe have the same boundary conditions, the variable names  are "recycled".
    //  They stay distinguishable in the input property file

    /*! \page propertyPage1
      \section pipeDiameter
      <code>ZFSInt ZFSStrctrdBlck::m_pipeDiameter </code>\n
      default = <code> 1.0 </code>\n \n
      Diameter of the pipe.\n
      Possible values are:\n
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelHeight = *(ZFSContext::getProperty("pipeDiameter", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));
    m_channelWidth = m_channelHeight;

    /*! \page propertyPage1
      \section pipeLength
      <code>ZFSInt ZFSStrctrdBlck::m_pipeLength </code>\n
      default = <code> 1.0 </code>\n \n
      Length of the pipe.\n
      Possible values are:\n
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelLength = *(ZFSContext::getProperty("pipeLength", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));

    /*! \page propertyPage1
      \section pipeInflowCoordinate
      <code>ZFSInt ZFSStrctrdBlck::m_channelInflowCoordinate </code>\n
      default = <code> 1.0 </code>\n \n
      Coordinate of the channel inflow plane.\n
      Possible values are:\n
      <ul>
      <li>Floating point</li>
      </ul>
      Keywords: <i>CHANNEL, IO, STRCTRD</i>
    */
    m_channelInflowPlaneCoordinate = *(ZFSContext::getProperty("pipeInflowCoordinate", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));
    m_channelC1 = *(ZFSContext::getProperty("loglawC1", m_blockId, __CALLING_FUNCTION__, &m_channelC1)->asFloat(0));
    m_channelC2 = *(ZFSContext::getProperty("loglawC2", m_blockId, __CALLING_FUNCTION__, &m_channelC1)->asFloat(0));
    m_channelC3 = *(ZFSContext::getProperty("loglawC3", m_blockId, __CALLING_FUNCTION__, &m_channelC1)->asFloat(0));
    m_channelC4 = *(ZFSContext::getProperty("loglawC4", m_blockId, __CALLING_FUNCTION__, &m_channelC1)->asFloat(0));
    zfs_log << "============= Pipe flow activated =============" << endl;
    zfs_log << "-> pipeRadius: " << m_channelHeight << endl;
    zfs_log << "-> pipeLength: " << m_channelLength << endl;
    zfs_log << "-> pipeInflowPlaneCoordinate: " << m_channelInflowPlaneCoordinate << endl;
    zfs_log << "-> Log law properties: " << endl;
    zfs_log << "--> C1: " << m_channelC1 << endl;
    zfs_log << "--> C2: " << m_channelC2 << endl;
    zfs_log << "--> C3: " << m_channelC3 << endl;
    zfs_log << "--> C4: " << m_channelC4 << endl;
    zfs_log << "============= Pipe flow summary finished =============" << endl;
    break;
  }
  default:
    break;
  }


  /*! \page propertyPage1
    \section referenceTemperature
    <code>ZFSFloat ZFSStrctrdBlck::m_referenceTemperature </code>\n
    default = <code>273.15</code>\n \n
    Reference temperature \f$ T_{\mathrm{ref}}\f$
    Used to scale the Sutherland's constant as follows: \f$ S/T_{\mathrm{ref}} \f$
    Also used for the computation of the reference sound speed and combustion (TF) related quantities
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_referenceTemperature = 273.15;
  m_referenceTemperature = *(ZFSContext::getProperty("referenceTemperature", m_blockId, __CALLING_FUNCTION__, &m_referenceTemperature)->asFloat(0));

  /*! \page propertyPage1
    \section sutherlandConstant
    <code>ZFSFloat ZFSStrctrdBlck::m_sutherlandConstant </code>\n
    default = <code>110.4 K</code>\n \n
    Sutherland's constant. Used by Sutherland's law.
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_sutherlandConstant = 110.4;
  m_sutherlandConstant = *(ZFSContext::getProperty("sutherlandConstant", m_blockId, __CALLING_FUNCTION__, &m_sutherlandConstant)->asFloat(0));
  m_sutherlandConstant /= m_referenceTemperature;
  m_sutherlandPlusOne = m_sutherlandConstant + F1;

  /*! \page propertyPage1
    \section tripSandpaper
    <code>ZFSInt ZFSStrctrdBlck::m_useSandpaperTrip </code>\n
    default = <code> false </code>\n \n
    Activate sandpaper trip forcing.\n
    Possible values are:\n
    <ul>
    <li>Bool: True/False</li>
    </ul>
    Keywords: <i>TRIP, BOUNDARYLAYER, STRCTRD</i>
  */
  m_useSandpaperTrip = false;
  if(ZFSContext::propertyExists("tripSandpaper", m_blockId )){
    m_useSandpaperTrip = *(ZFSContext::getProperty("tripSandpaper", m_blockId, __CALLING_FUNCTION__, &m_useSandpaperTrip )->asInt(0));
  }

  /*! \page propertyPage1
    \section considerVolumeForces
    <code>ZFSInt ZFSStrctrdBlck::m_considerVolumeForces </code>\n
    default = <code> 0 </code>\n \n
    Trigger to use volume forces\n
    Possible values are:\n
    <ul>
    <li>0: off</li>
    <li>1: on</li>
    </ul>
    Keywords: <i>FORCING, STRCTRD</i>
  */
  zfsAlloc(m_volumeForce, nDim, "m_volumeAcceleration", F0, __CALLING_FUNCTION__);
  m_considerVolumeForces = *(ZFSContext::getProperty("considerVolumeForces", m_blockId, __CALLING_FUNCTION__, &tmpFalse )->asInt(0));

  if(m_considerVolumeForces) {
    for( ZFSId i=0; i<nDim; i++ ) {
      /*! \page propertyPage1
        \section volumeForce
        <code>ZFSInt ZFSStrctrdBlck::volumeForce </code>\n
        default = <code> 0 </code>\n \n
        Numerical value of the volume force\n
        in each space direction.\n
        Possible values are:\n
        <ul>
        <li>floating point number</li>
        </ul>
        Keywords: <i>FORCING, STRCTRD</i>
      */
      m_volumeForce[ i ] = *(ZFSContext::getProperty("volumeForce", m_blockId, __CALLING_FUNCTION__, m_volumeForce, nDim )->asFloat(i));
    }
  }

  ZFSString govEqs = "NAVIER_STOKES";
  govEqs = *(ZFSContext::getProperty("govEqs", m_blockId, AT_, &govEqs)->asString(0));

  m_euler = false;
  if(string2enum(govEqs) == EULER) {
    m_euler = true;
  }
}

/*! \brief Reads and initializes properties associated with the Moving Grid Methods
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::setMovingGridProperties(){
  TRACE();  

  m_movingGridTimeOffset = F0;
  m_movingGridStepOffset = 0;
  m_movingGridInitialStart = true;
  m_gridMovingMethod = 0;
  m_movingGrid = false;
  m_wallVel = F0;
  m_movingGridTimeOffset = F0;

  m_travelingWave = false;
  m_waveTimeStepComputed = false;
  m_waveSpeed = 0.0;
  m_waveBeginTransition = 0.0;
  m_waveEndTransition = 0.0;
  m_waveRestartFadeIn = false;
  m_waveLength = 0.0;
  m_waveAmplitude = 0.0;
  m_waveTime = 0.0;
  m_waveCellsPerWaveLength = -1;
  m_waveNoStepsPerCell = -1;
  m_synchronizedMGOutput=0;


  /*! \page propertyPage1
    \section movingGrid
    <code>ZFSInt ZFSStrctrdBlck::m_movingGrid </code>\n
    default = <code> 0 </code>\n \n
    Trigger to use moving grid methods\n
    Possible values are:\n
    <ul>
    <li>0: off</li>
    <li>1: on</li>
    </ul>
    Keywords: <i>MOVING, STRCTRD</i>
  */
  if(ZFSContext::propertyExists("movingGrid", m_blockId )){
    m_movingGrid = *(ZFSContext::getProperty("movingGrid", m_blockId,
                                             __CALLING_FUNCTION__,&m_movingGrid)->asInt(0));
  }

  /*! \page propertyPage1
    \section gridMovingMethod
    <code>ZFSInt ZFSStrctrdBlck::m_gridMovingMethod </code>\n
    default = <code> 0 </code>\n \n
    Number of the moving grid function\n
    specified in the switch in zfsstrctrdblck3d\n
    Possible values are:\n
    <ul>
    <li>Integer >= 0</li>
    </ul>
    Keywords: <i>MOVING, STRCTRD</i>
  */
  if (m_movingGrid) {
    m_gridMovingMethod = *(ZFSContext::getProperty("gridMovingMethod", m_blockId,
                                                   __CALLING_FUNCTION__,(ZFSInt*) NULL)->asInt(0));

    /*! \page propertyPage1
      \section wallVel
      <code>ZFSInt ZFSStrctrdBlck::m_wallVel </code>\n
      default = <code> 0 </code>\n \n
      Value of the wall velocity for\n
      certain moving grid methods.\n
      Possible values are:\n
      <ul>
      <li>Floating point</li>
      </ul>
      Keywords: <i>MOVING, STRCTRD</i>
    */
    if(ZFSContext::propertyExists("wallVel", m_blockId )){
      m_wallVel = *(ZFSContext::getProperty("wallVel", m_blockId,
                                            __CALLING_FUNCTION__,(ZFSFloat*) NULL)->asFloat(0));
    }

    /*! \page propertyPage1
      \section movingGridTimeOffset
      <code>ZFSInt ZFSStrctrdBlck::m_movingGridTimeOffset </code>\n
      default = <code> 0 </code>\n \n
      Numerical value of the time offset\n
      for moving grid functions.\n
      Possible values are:\n
      <ul>
      <li>Floating point</li>
      </ul>
      Keywords: <i>MOVING, STRCTRD</i>
    */
    if(ZFSContext::propertyExists("movingGridTimeOffset", m_blockId )){
      m_movingGridTimeOffset = *(ZFSContext::getProperty("movingGridTimeOffset", m_blockId,
                                                         __CALLING_FUNCTION__, &m_movingGridTimeOffset)->asFloat(0));
    }

    /*! \page propertyPage1
      \section movingGridSaveGrid
      <code>ZFSInt ZFSStrctrdBlck::m_movingGridSaveGrid </code>\n
      default = <code> 0 </code>\n \n
      Trigger to write out the moving grid\n
      at each solution time step.\n
      Possible values are:\n
      <ul>
      <li>0: off</li>
      <li>1: on</li>
      </ul>
      Keywords: <i>MOVING, STRCTRD</i>
    */
    m_movingGridSaveGrid = 0;
    if(ZFSContext::propertyExists("movingGridSaveGrid", m_blockId )){
      m_movingGridSaveGrid = *(ZFSContext::getProperty("movingGridSaveGrid", m_blockId,
                                                       __CALLING_FUNCTION__, &m_movingGridSaveGrid)->asInt(0));
    }


    /*! \page propertyPage1
      \section synchronizedMGOutput
      <code>ZFSInt ZFSStrctrdBlck::m_synchronizedMGOutput </code>\n
      default = <code> 0 </code>\n \n
      Trigger to write out the moving grid synchronized with the moving grid, i.e.,\n
      the solution is written out in an integer of timesteps to move one cell further.\n
      Possible values are:\n
      <ul>
      <li>0: off</li>
      <li>1: on</li>
      </ul>
      Keywords: <i>MOVING, STRCTRD, IO</i>
    */
    m_synchronizedMGOutput = 0;
    m_synchronizedMGOutput = *(ZFSContext::getProperty("synchronizedMGOutput", m_blockId,
                                                       __CALLING_FUNCTION__, &m_synchronizedMGOutput)->asInt(0));
    
    zfs_log<< "synchronizedMGOutput is activated? " << m_synchronizedMGOutput << endl;
    

    if(m_gridMovingMethod == 9 || 
       m_gridMovingMethod == 10 ||
       m_gridMovingMethod == 11) {
      m_travelingWave = true;
    }
  }
}

/*! \fn void ZFSStrctrdBlck<nDim>::setNumericalProperties()
 * \brief Reads and initializes properties associated with the numerical method.
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::setNumericalProperties() {
  TRACE();   

  /*! \page propertyPage1
    \section nonBlockingComm
    <code>ZFSInt ZFSStrctrdBlck::m_nonBlockingComm </code>\n
    default = <code> 0 </code>\n \n
    Switch on non-blocking communication (not working yet)\n
    Possible values are:\n
    <ul>
    <li>0: off</li>
    <li>1: on</li>
    </ul>
    Keywords: <i>COMMUNICATION, STRCTRD</i>
  */
  m_nonBlockingComm =0; //-->false
  m_nonBlockingComm = *(ZFSContext::getProperty("nonBlockingComm", m_blockId, __CALLING_FUNCTION__, &m_nonBlockingComm)->asInt(0));

  /*! \page propertyPage1
    \section constantTimeStep
    <code>ZFSInt ZFSStrctrdBlck::m_constantTimeStep </code>\n
    default = <code> 1 </code>\n \n
    Trigger the use of a constant time step\n
    (only computed once at startup)\n
    Possible values are:\n
    <ul>
    <li>0: off</li>
    <li>1: on</li>
    </ul>
    Keywords: <i>TIMESTEP, STRCTRD</i>
  */
  m_constantTimeStep= 1; //default is set to true
  m_constantTimeStep = *(ZFSContext::getProperty("constantTimeStep", m_blockId, __CALLING_FUNCTION__, &m_constantTimeStep)->asInt(0));

  /*! \page propertyPage1
    \section localTimeStep
    <code>ZFSInt ZFSStrctrdBlck::m_localTimeStep </code>\n
    default = <code> 0 </code>\n \n
    Trigger the use of local time-stepping.
    Possible values are:\n
    <ul>
    <li>0: off</li>
    <li>1: on</li>
    </ul>
    Keywords: <i>TIMESTEP, STRCTRD</i>
  */
  m_localTimeStep= 0;
  m_localTimeStep= *(ZFSContext::getProperty("localTimeStep", m_blockId, __CALLING_FUNCTION__, &m_localTimeStep)->asInt(0));

  /*! \page propertyPage1
    \section timeStepComputationInterval
    <code>ZFSInt ZFSStrctrdBlck::m_timeStepComputationInterval </code>\n
    default = <code> 1 </code>\n \n
    Set the interval for the recomputation of\n
    the time step.\n
    Possible values are:\n
    <ul>
    <li>Integer >= 1</li>
    </ul>
    Keywords: <i>TIMESTEP, STRCTRD</i>
  */
  m_timeStepComputationInterval=1; //time step computation interval
  if(!m_constantTimeStep) {
    m_timeStepComputationInterval = *(ZFSContext::getProperty("timeStepComputationInterval", m_blockId, __CALLING_FUNCTION__, &m_timeStepComputationInterval)->asInt(0));
  }

  /*! \page propertyPage1
    \section noGhostLayers
    <code>ZFSInt ZFSStrctrdBlck::m_noGhostLayers </code>\n
    default = <code> 2 </code>\n \n
    Number of ghost-layers around the active mesh.\n
    the time step.\n
    Possible values are:\n
    <ul>
    <li>Integer >= 1</li>
    </ul>
    Keywords: <i>GRID, STRCTRD</i>
  */
  m_noGhostLayers = *(ZFSContext::getProperty("noGhostLayers", m_blockId, __CALLING_FUNCTION__,(ZFSInt*) NULL)->asInt(0));


  /*! \page propertyPage1
    \section cfl
    <code>ZFSId ZFSStrctrdBlck::m_cfl </code>\n
    default = <code>no default</code>\n \n
    Courant number C - Factor of the CFL condition \n \n
    possible values are:
    <ul>
    <li> positive floating point values < stability limit of the time-stepping method </li>
    <li> For default RK5 scheme the theoretical stability limit is C<4. Due to cut and small cells C<1.5 or even C <= 1.0 is recommended. </li>
    </ul>
    Keywords: <i> FINITE_VOLUME, STABILITY, TIME_INTEGRATION </i>
  */
  m_cfl = *(ZFSContext::getProperty("cfl", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));

  /*! \page propertyPage1
    \section limiter
    <code>ZFSId ZFSStrctrdBlck::m_limiter </code>\n
    default = <code>0</code>\n \n
    Trigger the use of the limiter.\n
    possible values are:
    <ul>
    <li> 0 </li>
    </ul>
    Keywords: <i> FINITE_VOLUME, STABILITY, LIMITER, STRCTRD </i>
  */
  m_limiter = 0;
  m_limiter = *(ZFSContext::getProperty("limiter", m_blockId, __CALLING_FUNCTION__, &m_limiter )->asInt(0));

  if(m_limiter) {
    /*! \page propertyPage1
      \section limiterMethod
      <code>ZFSId ZFSStrctrdBlck::m_limiterMethod </code>\n
      default = <code>ALBADA</code>\n \n
      Name of the limiter to use.\n
      possible values are:
      <ul>
      <li> ALBADA, VENKATAKRISHNAN, MINMOD</li>
      </ul>
      Keywords: <i> FINITE_VOLUME, STABILITY, LIMITER, STRCTRD </i>
    */
    m_limiterMethod = "ALBADA";
    m_limiterMethod = *(ZFSContext::getProperty("limiterMethod", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString());
  }
  
  /*! \page propertyPage1
      \section musclScheme
      <code>ZFSId ZFSStrctrdBlck::m_musclScheme </code>\n
      default = <code>Standard</code>\n \n
      Sets the MUSCL scheme for the structured block.
      possible values are:
      <ul>
      <li> Other valid MUSCL schemes </li>
      </ul>
      Keywords: <i> STRCTRD, FV, MUSCL</i>
  */
  /*! \page propertyPage1
    \section musclScheme
    <code>ZFSId ZFSStrctrdBlck::m_musclScheme </code>\n
    default = <code>Standard</code>\n \n
    Name of the MUSCL scheme to be used\n
    possible values are:
    <ul>
    <li> Standard, Stretched</li>
    </ul>
    Keywords: <i> FINITE_VOLUME, STABILITY, LIMITER, STRCTRD </i>
  */
  m_musclScheme="Standard";
  m_musclScheme = *(ZFSContext::getProperty("musclScheme", m_blockId, AT_, &m_musclScheme)->asString(0));
  
  /*! \page propertyPage1
    \section ausmScheme
    <code>ZFSId ZFSStrctrdBlck::m_ausmScheme </code>\n
    default = <code>Standard</code>\n \n
    Sets the AUSM scheme for the structured block.
    possible values are:
    <ul>
    <li> Other valid AUSM schemes </li>
    </ul>
    Keywords: <i> STRCTRD, FV, MUSCL</i>
*/
  /*! \page propertyPage1
    \section ausmScheme
    <code>ZFSId ZFSStrctrdBlck::m_ausmScheme </code>\n
    default = <code>Standard</code>\n \n
    Name of the AUSM scheme to be used.\n
    possible values are:
    <ul>
    <li>Standard, PTHRC, AUSMDV</li>
    </ul>
    Keywords: <i> FINITE_VOLUME, STABILITY, LIMITER, STRCTRD </i>
  */
  m_ausmScheme="Standard";
  m_ausmScheme = *(ZFSContext::getProperty("ausmScheme", m_blockId, AT_, &m_ausmScheme)->asString(0));

  /*! \page propertyPage1
    \section convergenceCriterion
    <code>ZFSId ZFSStrctrdBlck::convergenceCriterion </code>\n
    default = <code>Standard</code>\n \n
    Set the convergence criterion to stop computation.\n
    possible values are:
    <ul>
    <li>Float < 1.0</li>
    </ul>
    Keywords: <i> FINITE_VOLUME, CONVERGENCE, STRCTRD </i>
  */
  m_convergenceCriterion = 1e-12;
  m_convergenceCriterion = *(ZFSContext::getProperty("convergenceCriterion", m_blockId, __CALLING_FUNCTION__, &m_convergenceCriterion )->asFloat(0));

  /*! \page propertyPage1
    \section upwindCoefficient
    <code>ZFSId ZFSStrctrdBlck::m_chi</code>\n
    default = <code>0.0</code>\n \n
    Chi for AUSM pressure splitting.\n
    possible values are:
    <ul>
    <li>Float >= 0.0</li>
    </ul>
    Keywords: <i> FINITE_VOLUME, CONVERGENCE, STRCTRD </i>
  */
  m_chi = 0.0;
  m_chi = *(ZFSContext::getProperty("upwindCoefficient", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0));
}


template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::initializeTimers(){
  TRACE();
  NEW_TIMER_GROUP_NOCREATE(m_tgStr, "STRUCTURED_FV");
  NEW_TIMER_NOCREATE(m_tcomm, "Communication", m_tgStr);
  NEW_TIMER_NOCREATE(m_tdecomposition, "Grid Decomposition", m_tgStr);
  NEW_TIMER_NOCREATE(m_treadGrid, "Grid Reading", m_tgStr);
  NEW_TIMER_NOCREATE(m_tloadRestart, "Load Restart File", m_tgStr);
  NEW_TIMER_NOCREATE(m_tcomputeMetrics, "Compute metrics", m_tgStr);
  NEW_TIMER_NOCREATE(m_tcomputeJacobian, "Compute jacobian", m_tgStr);
  NEW_TIMER_NOCREATE(m_tConvectiveFlux, "Convective Flux", m_tgStr);
  NEW_TIMER_NOCREATE(m_tViscousFlux, "Viscous Flux", m_tgStr);
  NEW_TIMER_NOCREATE(m_tVolumeFlux, "Volume Flux", m_tgStr);

  NEW_TIMER_NOCREATE(m_tbuildSponge, "Build up sponge", m_tgStr);
  NEW_TIMER_NOCREATE(m_tspanwiseReorder, "Wave Spanwise Reordering", m_tgStr);

  NEW_SUB_TIMER_NOCREATE(m_tloadRestartVars, "Load conservative varialbes", m_tloadRestart);
  NEW_SUB_TIMER_NOCREATE(m_tloadRestartSponge, "Load sponge data", m_tloadRestart);
  NEW_SUB_TIMER_NOCREATE(m_tloadRestartStg, "Load STG restart data", m_tloadRestart);
  NEW_SUB_TIMER_NOCREATE(m_tloadRestartStgEddies, "Read STG Eddies", m_tloadRestartStg);
  NEW_SUB_TIMER_NOCREATE(m_tloadRestartStgRead, "Read STG Values", m_tloadRestartStg);
  NEW_SUB_TIMER_NOCREATE(m_tloadRestartStgBcast, "Broadcast STG Values", m_tloadRestartStg);
  NEW_SUB_TIMER_NOCREATE(m_tloadRestartShift, "Shift restart variables", m_tloadRestart);

  NEW_SUB_TIMER_NOCREATE(m_texchange, "exchange",m_tcomm);
  NEW_SUB_TIMER_NOCREATE(m_tgatherAndSend, "gatherAndSend",m_texchange);
  NEW_SUB_TIMER_NOCREATE(m_tgatherAndSendWait, "gatherAndSendWait",m_tgatherAndSend);

  NEW_SUB_TIMER_NOCREATE(m_tscatterAndReceive, "scatterAndReceive",m_texchange);
  NEW_SUB_TIMER_NOCREATE(m_tscatterWaitSome, "gatherAndSendWait",m_tscatterAndReceive);
  NEW_SUB_TIMER_NOCREATE(m_tgather, "gather",m_texchange);
  NEW_SUB_TIMER_NOCREATE(m_tsend, "send",m_texchange);
  NEW_SUB_TIMER_NOCREATE(m_tgatherAndSendWait, "wait",m_texchange);
  NEW_SUB_TIMER_NOCREATE(m_treceive, "receive",m_texchange);
  NEW_SUB_TIMER_NOCREATE(m_tscatter, "scatter",m_texchange);
  NEW_SUB_TIMER_NOCREATE(m_tsendWait, "sendWait",m_texchange);
  NEW_SUB_TIMER_NOCREATE(m_treceiveWait, "receiveWait",m_texchange);

  NEW_SUB_TIMER_NOCREATE(m_texchangeDt, "exchange time-step",m_tcomm);
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::createMPIGroups() {
  //Channel Communication
  m_commChannelIn= new MPI_Comm;//MPI_COMM_NULL;
  m_commChannelOut= new MPI_Comm;//MPI_COMM_NULL;
  *m_commChannelOut = MPI_COMM_NULL;
  *m_commChannelIn = MPI_COMM_NULL;
  m_commChannelWorld = new MPI_Comm;
  *m_commChannelWorld = MPI_COMM_NULL;
  m_channelRoots=new ZFSInt[4];

  //STG Communication
  m_commStg = new MPI_Comm;
  *m_commStg = MPI_COMM_NULL;
  m_commStgRoot = new ZFSInt;
  *m_commStgRoot = -1;
  m_commStgRootGlobal = new ZFSInt;
  *m_commStgRootGlobal = -1;

  // //BC2600 for zonal //junoh
  m_commBC2600 = new MPI_Comm;
  *m_commBC2600 = MPI_COMM_NULL;
  m_commBC2600Root = new ZFSInt;
  *m_commBC2600Root = -1;
  m_commBC2600RootGlobal = new ZFSInt;
  *m_commBC2600RootGlobal = -1;
  //zonal Communication
  
  // m_commZonal = new MPI_Comm;
  // *m_commZonal = MPI_COMM_NULL;
  // m_commZonalRoot = new ZFSInt;
  // *m_commZonalRoot = -1;
  // m_commZonalRootGlobal = new ZFSInt;
  // *m_commZonalRootGlobal = -1;
  
  
  //Periodic rotation boundary condtions communicators
  m_commPerRotOne = new MPI_Comm;
  m_commPerRotTwo = new MPI_Comm;
  m_commPerRotWorld = new MPI_Comm;
  m_commPerRotRoots= new ZFSInt[4];
  m_commPerRotGroup = -1;

  //RESCALING CommunicationGroup
  m_rescalingCommGrComm = new MPI_Comm;
  *m_rescalingCommGrComm = MPI_COMM_NULL;
  m_rescalingCommGrRoot = new ZFSInt;
  *m_rescalingCommGrRoot = -1;
  m_rescalingCommGrRootGlobal = new ZFSInt;
  *m_rescalingCommGrRootGlobal = -1;
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////MPI Zonal///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

  m_commZonal = new MPI_Comm[m_noInputBlocks];
  for(ZFSId i=0; i<m_noInputBlocks; i++) {
    m_commZonal[i] = MPI_COMM_NULL;
  }
  zfsAlloc(m_commZonalRoot, m_noInputBlocks, "m_commZonalRoot", -1, __CALLING_FUNCTION__);
  zfsAlloc(m_commZonalRootGlobal, m_noInputBlocks, "m_commZonalRoot", -1, __CALLING_FUNCTION__);
  m_commZonalMyRank = -1;
  m_zonalRootRank = false;

  
  //first collect the ranks of each input block in one array and count the number
  for(ZFSId i=0; i<m_noInputBlocks; i++) {
    MPI_Group groupZonal, newGroupZonal;
    // ZFSInt tmpSize = -1;//no of domains in inputBlock[i]
    // ZFSIntScratchSpace tmpPartitions(tmpSize, __CALLING_FUNCTION__, "tmpPartitions"); //this one needs to be filled with all partitions belonging to inputBlock[i]

    vector<ZFSInt> tmpPartitions;
    ZFSInt blockDomainId =0;
    ZFSInt hasBlockDomain = 0;
    ZFSInt hasBlockDomainGlobal=0;
    ZFSIntScratchSpace nblockDomainArray(noDomains(), __CALLING_FUNCTION__, "nblockDomainArray");
    ZFSIntScratchSpace nblockDomainOffset(noDomains(), __CALLING_FUNCTION__, "nblockDomainOffset");
      if(m_inputBlockId==i) {
	blockDomainId = domainId();
	hasBlockDomain=1;
	// cout<<"m_partition->outputBoxInfo["<<j<<"]->inputBoxID:"<<m_partition->outputBoxInfo[j]->inputBoxID<<"domainId:"<<domainId()<<endl;
      }
      MPI_Allreduce(&hasBlockDomain,&hasBlockDomainGlobal,1,MPI_INT,MPI_SUM,m_zfsStrctrdComm);
      MPI_Allgather(&hasBlockDomain,1,MPI_INT,&nblockDomainArray[0], 1, MPI_INT,m_zfsStrctrdComm);
      // cout<<"hasBlockDomainGlobal:"<<hasBlockDomainGlobal<<" m_inputBlockId:"<<m_inputBlockId<<endl;
      nblockDomainOffset[0]=0;
      for (ZFSId j= 1; j< noDomains(); j++ ){ 
	nblockDomainOffset[j]= nblockDomainOffset[j-1] + nblockDomainArray[j-1];
      } 
      ZFSIntScratchSpace zonalRanks(hasBlockDomainGlobal, __CALLING_FUNCTION__, "zonalRanks");
      // cout<<"blockDomainId:"<<blockDomainId<<endl;
      // for(ZFSId j=0; j<noDomains(); j++){
      // 	cout<<"nblockDomainArray["<<j<<"]:"<<nblockDomainArray[j]<<endl;
      // 	cout<<"nblockDomainOffset["<<j<<"]:"<<nblockDomainOffset[j]<<endl;
      // }
      MPI_Allgatherv(&blockDomainId, nblockDomainArray[domainId()], MPI_INT, &zonalRanks[0], &nblockDomainArray[0], &nblockDomainOffset[0], MPI_INT, m_zfsStrctrdComm);


    MPI_Barrier(m_zfsStrctrdComm);
      int zonalcommsize = hasBlockDomainGlobal;   //no of domains in inputBlock[i]

    MPI_Comm_group(m_zfsStrctrdComm, &groupZonal);
    MPI_Group_incl(groupZonal, zonalcommsize, &zonalRanks[0], &newGroupZonal);
    MPI_Comm_create(m_zfsStrctrdComm, newGroupZonal, &m_commZonal[i]);

    if(domainId() == zonalRanks[0]) {
      MPI_Comm_rank(m_commZonal[i], &m_commZonalRoot[i]);
      MPI_Comm_rank(m_zfsStrctrdComm, &m_commZonalRootGlobal[i]);
      m_zonalRootRank = true;
    }
    

    MPI_Barrier(m_zfsStrctrdComm);
    MPI_Bcast(&m_commZonalRoot[0], m_noInputBlocks, MPI_INT, zonalRanks[0],m_zfsStrctrdComm);
    MPI_Bcast(&m_commZonalRootGlobal[0], m_noInputBlocks, MPI_INT, zonalRanks[0], m_zfsStrctrdComm);

    // for(ZFSId j=0; j<m_noInputBlocks; j++){
    //   cout<<"m_commZonalRoot["<<j<<"]:"<<m_commZonalRoot[j]<<endl;
    //   cout<<"m_commZonalRootGlobal["<<j<<"]:"<<m_commZonalRootGlobal[j]<<endl;
    // }
  }
}
  //   MPI_Barrier(m_zfsStrctrdComm);
  //   cout<<"check1"<<endl;
  //   ZFSInt* zonalranks;

  //   if(tmpPartitions.size() != 0) {
  //     int counterzonal = 0;
  //     zonalranks = new ZFSInt[tmpPartitions.size()];
  //     for(ZFSId j = 0; j < (ZFSInt)tmpPartitions.size(); j++) {
  // 	for(ZFSId k = 0; k < noDomains(); k++) {
  // 	  if(tmpPartitions[j] == k) {
  // 	    zonalranks[counterzonal] = k;
  // 	    counterzonal++;
  // 	    cout<<"zonalranks["<<counterzonal<<"]:"<<zonalranks[counterzonal]<<"domainId:"<<domainId()<<endl;
  //         }
  // 	}
  //     }
  //   }
  //   cout<<"check2"<<endl;
  //     int zonalcommsize = (int)(tmpPartitions.size());   //no of domains in inputBlock[i]

  //   MPI_Comm_group(m_zfsStrctrdComm, &groupZonal);
  //   MPI_Group_incl(groupZonal, zonalcommsize, zonalranks, &newGroupZonal);
  //   MPI_Comm_create(m_zfsStrctrdComm, newGroupZonal, &m_commZonal[i]);
  //   cout<<"check3"<<endl;
  //   if(domainId() == tmpPartitions[0]) {
  //     MPI_Comm_rank(m_commZonal[i], &m_commZonalRoot[i]);
  //     MPI_Comm_rank(m_zfsStrctrdComm, &m_commZonalRootGlobal[i]);
  //     m_zonalRootRank = true;
  //     cout<<"check domainId:"<<domainId()<<endl;
  //   }
  //   cout<<"check4"<<endl;
  //   MPI_Barrier(m_zfsStrctrdComm);
  //   // MPI_Bcast(&m_commZonalRoot[0], m_noInputBlocks, MPI_INT, tmpPartitions[0], m_commZonal[i]);
  //   // MPI_Bcast(&m_commZonalRootGlobal[0], m_noInputBlocks, MPI_INT, tmpPartitions[0], m_zfsStrctrdComm);
  // }
  // /////////////////////////////////////////////////////////////////////////////////////////////////

  
// }


template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::readAndSetSpongeLayerProperties(){
   TRACE();

  //initialize the values for the Sponge!
  m_spongeLayerThickness=NULL;
  m_betaSponge=NULL;
  m_sigmaSponge=NULL;
  m_noSpongeDomainInfos=0; //number of bc/windows for bc
  m_spongeLayerType=1;
  m_targetDensityFactor=F0;
  m_computeSpongeFactor=true;
  ZFSInt readSpongeFromBc=1;

  /*! \page propertyPage1
    \section useSponge
    <code>ZFSInt ZFSStrctrdBlck::m_useSponge </code>\n
    default = <code> 0</code>\n \n
    Trigger to use the sponge.\n
    Possible values are:\n
    <ul>
    <li>true/false</li>
    </ul>
    Keywords: <i>SPONGE, IO, STRCTRD</i>
  */
  m_useSponge =0;
  m_useSponge = *(ZFSContext::getProperty("useSponge", m_blockId, __CALLING_FUNCTION__, &m_useSponge )->asInt(0));

  if(m_useSponge){//yes
    FQ->neededFQVariables[FQ->SPONGE_FACTOR] = 1;

    /*! \page propertyPage1
      \section readSpongeFromBC
      <code>ZFSInt ZFSStrctrdBlck::m_readSpongeFromBC </code>\n
      default = <code> 0</code>\n \n
      Use the given BC numbers to apply sponge\n
      to each window connected with this BC number,\n
      otherwise use given windows IDs.\n
      Possible values are:\n
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>SPONGE, IO, STRCTRD</i>
    */
    readSpongeFromBc = *(ZFSContext::getProperty("readSpongeFromBC", m_blockId, __CALLING_FUNCTION__, &readSpongeFromBc )->asInt(0));


    ZFSProperty * spongeIds;
    if(readSpongeFromBc){
      /*! \page propertyPage1
        \section spongeBndryCndIds
        <code>ZFSInt ZFSStrctrdBlck::spongeIds </code>\n
        default = <code> 0</code>\n \n
        Use the given BC numbers to apply sponge\n
        to each window connected with this BC number.\n
        Possible values are:\n
        <ul>
        <li>BC Id</li>
        </ul>
        Keywords: <i>SPONGE, IO, STRCTRD</i>
      */
      spongeIds = ZFSContext::getProperty("spongeBndryCndIds", m_blockId, __CALLING_FUNCTION__, (ZFSId*) NULL, 1);
    } else {

      /*! \page propertyPage1
        \section spongeWindowIds
        <code>ZFSInt ZFSStrctrdBlck::m_spongeIds </code>\n
        default = <code> 0</code>\n \n
        Use the given window IDs to apply\n
        a sponge to them.\n
        <ul>
        <li>Window ID</li>
        </ul>
        Keywords: <i>SPONGE, IO, STRCTRD</i>
      */
      spongeIds = ZFSContext::getProperty("spongeWindowIds", m_blockId, __CALLING_FUNCTION__, (ZFSId*) NULL, 1);
    }
    m_noSpongeDomainInfos = spongeIds->elements;
    m_noSpongeDomainInfos = (ZFSInt)spongeIds->count();

    /*! \page propertyPage1
      \section spongeLayerType
      <code>ZFSInt ZFSStrctrdBlck::m_spongeLayerType </code>\n
      default = <code> 0</code>\n \n
      Type of the sponge layer, i.e., \n
      sponge to pressure, density, both,\n
      infinity values or predefined field etc.\n
      <ul>
      <li>Integer of Sponge type</li>
      </ul>
      Keywords: <i>SPONGE, IO, STRCTRD</i>
    */
    m_spongeLayerType = *(ZFSContext::getProperty("spongeLayerType", m_blockId, __CALLING_FUNCTION__, &m_spongeLayerType )->asInt(0));


    ZFSProperty* spongeBeta= ZFSContext::getProperty("betaSponge", m_blockId , __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);
    ZFSProperty* spongeSigma= ZFSContext::getProperty("sigmaSponge", m_blockId , __CALLING_FUNCTION__, (ZFSFloat*) NULL, 1);

    //The programm aborts if number of sponge properties does not fit together
    if(m_noSpongeDomainInfos!= spongeBeta->count() || m_noSpongeDomainInfos != spongeSigma->count()){
      zfsTerm(1,__CALLING_FUNCTION__,  "The number of sponge properties does not match");
    }

    //else we can allocate the memory necessary for the sponge etc'
    zfsAlloc(m_spongeLayerThickness, m_noSpongeDomainInfos, "m_spongeLayerThicknesses", F0, __CALLING_FUNCTION__);
    zfsAlloc(m_sigmaSponge, m_noSpongeDomainInfos, "m_sigmaSponge",  __CALLING_FUNCTION__);
    zfsAlloc(m_betaSponge, m_noSpongeDomainInfos, "m_betaSponge",  __CALLING_FUNCTION__);
    zfsAlloc(m_spongeBcWindowInfo, m_noSpongeDomainInfos, "m_spongeBcWindowInfo", __CALLING_FUNCTION__);

    //read all the parameters
    for(ZFSInt i=0; i<m_noSpongeDomainInfos; ++i){

      /*! \page propertyPage1
        \section spongeLayerThickness
        <code>ZFSFloat ZFSStrctrdBlck::m_spongeLayerThickness </code>\n
        default = <code>-1.0</code>\n \n
        The property controls the thickness of the sponge layer in which the sponge layer forcing is applied. The sponge forcing term added to the rhs of a cell inside the sponge layer is given by \n
        \f$ \Delta L(\phi) = V \sigma \frac{\Delta x_{sp}^2}{L_s^2} \Delta \phi  \f$,  \n
        where \f$ V \f$ is the cell volume, \f$  \sigma  \f$ is the forcing amplitude, \f$  x_{sp} \f$ is the inner sponge layer boundary, \f$  L_{sp} \f$ is the sponge layer thickness and \f$  \Delta \phi = \phi - \phi_{target}  \f$ is the difference between the local and the freesteam values of \f$  \phi  \f$.\n
        possible values are:
        <ul>
        <li> Any non-negative floating point value. </li>
        </ul>
        Keywords: <i> STRCTRD, SPONGE </i>
      */
      m_spongeLayerThickness[i]=-1.0;
      m_spongeLayerThickness[i] =*(ZFSContext::getProperty("spongeLayerThickness", m_blockId, __CALLING_FUNCTION__, &m_spongeLayerThickness[i] ,m_noSpongeDomainInfos )->asFloat(i));

      /*! \page propertyPage1
        \section betaSponge
        <code>ZFSFloat ZFSStrctrdBlock::m_betaSponge  </code>\n
        default = <code>0</code>\n \n
        The property controls the sponge function. Linear sponge spongeBeta = 1, quadratic spongeBeta = 2, etc. \n
        possible values are:
        <ul>
        <li>floating point values.</li>
        </ul>
        Keywords: <i>STRCTRD, SPONGE, BETA, PROFIL</i>
      */
      m_betaSponge[i]=F0;
      m_betaSponge[i] =*(ZFSContext::getProperty("betaSponge", m_blockId, __CALLING_FUNCTION__, &m_spongeLayerThickness[i] ,m_noSpongeDomainInfos )->asFloat(i));

      /*! \page propertyPage1
        \section sigmaSponge
        <code>ZFSFloat ZFSStrctrdBlock::m_sigmaSponge  </code>\n
        default = <code>0</code>\n \n
        Controls the sigma of the sponge, i.e., the strength of the sponge.\n
        possible values are:
        <ul>
        <li>floating point values.</li>
        </ul>
        Keywords: <i>STRCTRD, SPONGE, BETA, PROFIL</i>
      */
      m_sigmaSponge[i]=F0;
      m_sigmaSponge[i] =*(ZFSContext::getProperty("sigmaSponge", m_blockId, __CALLING_FUNCTION__, &m_spongeLayerThickness[i] ,m_noSpongeDomainInfos )->asFloat(i));
    }


    if(readSpongeFromBc){
      for(ZFSInt i=0; i<m_noSpongeDomainInfos; ++i){
        m_spongeBcWindowInfo[i]=-1;
        m_spongeBcWindowInfo[i] =*(ZFSContext::getProperty("spongeBndryCndIds", m_blockId, __CALLING_FUNCTION__, &m_spongeBcWindowInfo[i] ,m_noSpongeDomainInfos )->asInt(i));
      }
    }else {
      for(ZFSInt i=0; i<m_noSpongeDomainInfos; ++i){
        m_spongeBcWindowInfo[i]=-1;
        m_spongeBcWindowInfo[i] =*(ZFSContext::getProperty("spongeWindowIds", m_blockId, __CALLING_FUNCTION__, &m_spongeBcWindowInfo[i] ,m_noSpongeDomainInfos )->asInt(i));
      }
    }

    /*! \page propertyPage1
      \section targetDensityFactor
      <code>ZFSFloat ZFSStrctrdBlck::m_targetDensityFactor </code>\n
      default = <code>1.0</code>\n \n
      The property controls the Intensity of the sponge layer correction regarding the density forcing term for some values of the spongeLayerType. The sponge forcing term added to the rhs of a cell inside the sponge layer is given by \n
      \f$ \Delta L(\phi) = V \sigma \frac{\Delta x_{sp}^2}{L_s^2} \Delta \phi  \f$,  \n
      where \f$ V \f$ is the cell volume, \f$  \sigma  \f$ is the forcing amplitude, \f$  x_{sp} \f$ is the inner sponge layer boundary, \f$  L_{sp} \f$ is the sponge layer thickness and \f$  \Delta \phi = \phi - \phi_{target}  \f$ is the difference between the local and the freesteam values of \f$  \phi  \f$.\n
      The density target value is in these cases given as:    \n
      <code> deltaRho =a_pvariable( cellId ,  PV->RHO ) - m_rhoInfinity * m_targetDensityFactor;</code>\n
      See also spongeLayerType. Only meaningful and required with certain values for \ref spongeLayerType and if both \ref spongeLayerThickness and \ref sigmaSponge are specified and nonzero! ! \n \n
      possible values are:
      <ul>
      <li>Non-negative floating point values.</li>
      </ul>
      Keywords: <i>STRCTRD, SPONGE</i>
    */
    m_targetDensityFactor=F1;
    m_targetDensityFactor = *(ZFSContext::getProperty("targetDensityFactor", m_blockId, __CALLING_FUNCTION__, &m_targetDensityFactor )->asFloat(0));

    m_windowInfo->setSpongeInformation(m_noSpongeDomainInfos, m_betaSponge, m_sigmaSponge, m_spongeLayerThickness, m_spongeBcWindowInfo, readSpongeFromBc );

    /*! \page propertyPage1
      \section computeSpongeFactor
      <code>ZFSFloat ZFSStrctrdBlock::m_computeSpongeFactor  </code>\n
      default = <code>0</code>\n \n
      Trigger the sponge computation, if set to false.\n
      the sponge values will be read from the restart file\n
      which may be much faster due to the slow sponge computation.\n
      possible values are:
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>STRCTRD, SPONGE</i>
    */
    m_computeSpongeFactor = true;
    if(ZFSContext::propertyExists("computeSpongeFactor", m_blockId )){
      m_computeSpongeFactor = *(ZFSContext::getProperty("computeSpongeFactor", m_blockId, __CALLING_FUNCTION__, &m_computeSpongeFactor )->asInt(0));
    }

    //we now have the spongeProperties save in m_windowInfo->m_spongeInfoMap
    //allocating and other stuff can be handled in the boundaryCondition constructor

    if(m_spongeLayerType == 2) {
      //active the two sponge layer FQ fields for rho and rhoE
      FQ->neededFQVariables[FQ->SPONGE_RHO] = 1;
      FQ->neededFQVariables[FQ->SPONGE_RHO_E] = 1;
    }

    if(m_spongeLayerType == 4) {
      //active the two sponge layer FQ fields for rho and rhoE
      FQ->neededFQVariables[FQ->SPONGE_RHO] = 1;
    }
  }
}


///  \fn void ZFSStrctrdBlck<nDim>::setRungeKuttaProperties()
///  \brief This function reads the properties required for Runge Kutta time stepping.
///  \author Pascal S. Meysonnat (see FV-> Gonzalo G-B)
///  \date July, 2013
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::setRungeKuttaProperties() {
  TRACE();  

  //! Allocate and initialize Runge-Kutta coefficients:
  /*! \page propertyPage1
    \section noRKSteps
    <code>ZFSId ZFSStrctrdBlck::m_noRKSteps </code>\n
    default = <code>5</code>\n \n
    Number of steps in the Runge-Kutta time-stepping method.
    possible values are:
    <ul>
    <li>Positive integers.</li>
    </ul>
    Keywords: <i>FV, RUNGE KUTTA, TIME STEPPING</i>
  */
  m_noRKSteps = 5;
  m_noRKSteps = *(ZFSContext::getProperty("noRKSteps", m_blockId, __CALLING_FUNCTION__, (ZFSId*) &m_noRKSteps)->asInt(0));

  zfsAlloc( m_RKalpha, m_noRKSteps, "m_RKalpha", -F1, __CALLING_FUNCTION__ );

  ZFSFloat RK5DefaultCoeffs[5] = {0.25, 0.16666666666, 0.375, 0.5, 1.};
  /*! \page propertyPage1
    \section rkalpha-step
    <code>ZFSFloat ZFSStrctrdBlck::m_RKalpha[m_noRKSteps] </code>\n
    default = <code>0.25, 0.16666666666, 0.375, 0.5, 1</code> IF noRKSteps is 5.\n \n
    Coeffients of the Runge-Kutta time-stepping method.
    possible values are:
    <ul>
    <li>Floating point numbers (as many as Runge-Kutta steps).</li>
    </ul>
    Keywords: <i>FV, RUNGE KUTTA, TIME STEPPING</i>
  */
  if( m_noRKSteps == 5 ){ // default only valid m_noRKSteps == 5
    for(ZFSInt i = 0; i < m_noRKSteps; i++) {
      m_RKalpha[i] = *(ZFSContext::getProperty("rkalpha-step", m_blockId, __CALLING_FUNCTION__,(ZFSFloat*) &RK5DefaultCoeffs[i],m_noRKSteps)->asFloat(i));
    }
  }else{ // otherwise currently no default is set  -> can be extended later
    for(ZFSInt i = 0; i < m_noRKSteps; i++) {
      m_RKalpha[i] = *(ZFSContext::getProperty("rkalpha-step", m_blockId, __CALLING_FUNCTION__,(ZFSFloat*) NULL,m_noRKSteps)->asFloat(i));
    }
  }

  /*! \page propertyPage1
    \section rungeKuttaOrder
    <code>ZFSFloat ZFSStrctrdBlck::m_rungeKuttaOrder </code>\n
    default = <code>2</code>\n \n
    Defines the runge kutta method (order).
    possible values are:
    <ul>
    <li>2 - second order</li>
    <li>3 - third order</li>
    </ul>
    Keywords: <i>TIME_INTEGRATION, RUNGE_KUTTA</i>
  */
  m_rungeKuttaOrder = 2;
  m_rungeKuttaOrder = *(ZFSContext::getProperty("rungeKuttaOrder", m_blockId, __CALLING_FUNCTION__, &m_rungeKuttaOrder)->asInt(0));
}



template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::initializeRungeKutta()
{
  TRACE();
  setRungeKuttaProperties();

  if(!m_restart) {
    m_time = F0;
    m_physicalTime = F0;
    globalTimeStep=0;
  }

  m_RKStep = 0;

  // this is used to initialize the counter the very first time this is called
  // and to initialize the very first residual
  if( approx(m_time, F0, m_eps)) {
    m_workload = 0;
    m_workloadIncrement = 1;
    m_firstMaxResidual = F0;
    m_firstAvrgResidual =F0;
  }
}


template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::setSTGProperties()
{
  TRACE();
  m_stgIsActive = false;
  m_stgInitialStartup = false;
  m_stgFace = 0;

  if(nDim<3){
    return;
  }

  /*! \page propertyPage1
    \section useSTG
    <code>ZFSFloat ZFSStrctrdBlck::m_stgIsActive </code>\n
    default = <code>0</code>\n \n
    Trigger the use of the STG BC.\n
    possible values are:
    <ul>
    <li>true/false</li>
    </ul>
    Keywords: <i>STG, STRCTRD</i>
  */
  m_stgIsActive = false;
  if(ZFSContext::propertyExists("useSTG", m_blockId )){
    m_stgIsActive = *(ZFSContext::getProperty("useSTG", m_blockId, __CALLING_FUNCTION__, &m_stgIsActive )->asInt(0));
  }
    //switch this on for zonal computation   //junoh
  if(m_zonal) {
    m_stgIsActive = true;
  }
  if(domainId()==0){
  cout << "zonal: " << m_zonal << " stgIsActive:  " << m_stgIsActive << endl;
  }
  if(m_stgIsActive) {
    m_stgLocal = false;
    m_stgRootRank = false;

    for(ZFSInt i=0; i < abs((ZFSInt)m_windowInfo->physicalBCMap.size()); ++i) {
      if((m_windowInfo->physicalBCMap[i]->BC == 7909)||(m_windowInfo->physicalBCMap[i]->BC == 2221)) { //junoh      if(m_windowInfo->physicalBCMap[i]->BC == 7909) {
        m_stgLocal = true;
        m_stgFace = m_windowInfo->physicalBCMap[i]->face;
        break;
      }
    }

    if(m_stgLocal) {
      MPI_Comm_rank(*m_commStg, &m_commStgMyRank);
    } else {
      m_commStgMyRank = -1;
    }

    /*! \page propertyPage1
      \section stgSubSup
      <code>ZFSFloat ZFSStrctrdBlck::m_stgSubSup </code>\n
      default = <code>0</code>\n \n
      Use mixed subsonics/subsonic formulation\n
      of the STG boundary.\n
      possible values are:
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgSubSup = false;
    if(ZFSContext::propertyExists("stgSubSup", m_blockId )){
      m_stgSubSup = *(ZFSContext::getProperty("stgSubSup", m_blockId, __CALLING_FUNCTION__, &m_stgSubSup )->asInt(0));
    }

    /*! \page propertyPage1
      \section stgSupersonic
      <code>ZFSFloat ZFSStrctrdBlck::m_stgSupersonic </code>\n
      default = <code>0</code>\n \n
      Use supersonic STG boundary formulation.\n
      possible values are:
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgSupersonic = false;
    if(ZFSContext::propertyExists("stgSupersonic", m_blockId )){
      m_stgSupersonic = *(ZFSContext::getProperty("stgSupersonic", m_blockId, __CALLING_FUNCTION__, &m_stgSupersonic )->asInt(0));

      if(m_stgSupersonic && m_stgSubSup) {
        m_stgSubSup = false;
        if(domainId()==0) {
          cout << "WARNING: You activated conflicting properties stgSubSup "
               << "and stgSupersonic, thus only the pure supersonic formulation will be used. "
               << "Switch off stgSupersonic to get the mixed formulation" << endl;
        }
      }
    }

    /*! \page propertyPage1
      \section BLT1
      <code>ZFSFloat ZFSStrctrdBlck::m_stgBLT1 </code>\n
      default = <code>1.0</code>\n \n
      Defines the size of the STG virtual box\n
      in the x-direction as a fraction of the\n
      delta0 specified.\n
      possible values are:
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgBLT1 = 1.0;
    m_stgBLT1 = *(ZFSContext::getProperty("BLT1", m_blockId, __CALLING_FUNCTION__,&m_stgBLT1)->asFloat(0));

    /*! \page propertyPage1
      \section BLT2
      <code>ZFSFloat ZFSStrctrdBlck::m_stgBLT2 </code>\n
      default = <code>2</code>\n \n
      Defines the size of the STG virtual box\n
      in the y-direction as a fraction of the\n
      delta0 specified.\n
      possible values are:
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgBLT2 = 1.1;
    m_stgBLT2 = *(ZFSContext::getProperty("BLT2", m_blockId, __CALLING_FUNCTION__,&m_stgBLT2)->asFloat(0));

    /*! \page propertyPage1
      \section BLT3
      <code>ZFSFloat ZFSStrctrdBlck::m_stgBLT3 </code>\n
      default = <code>1.1</code>\n \n
      Defines the size of the STG virtual box\n
      in the z-direction as a fraction of the\n
      delta0 specified.\n
      possible values are:
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgBLT3 = 1.0;
    m_stgBLT3 = *(ZFSContext::getProperty("BLT3", m_blockId, __CALLING_FUNCTION__,&m_stgBLT3)->asFloat(0));

    /*! \page propertyPage1
      \section deltaIn
      <code>ZFSFloat ZFSStrctrdBlck::m_stgDelta99Inflow</code>\n
      default = <code>-1.0</code>\n \n
      Defines the delta0 thickness at the inflow for the STG.\n
      possible values are:
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgDelta99Inflow = -1.0;
    m_stgDelta99Inflow = *(ZFSContext::getProperty("deltaIn", m_blockId, __CALLING_FUNCTION__,&m_stgDelta99Inflow)->asFloat(0));

    m_stgBLT1 = m_stgBLT1 * m_stgDelta99Inflow;
    m_stgBLT2 = m_stgBLT2 * m_stgDelta99Inflow;

    zfsAlloc(m_stgLengthFactors, 3, "m_block->m_stgLengthFactors", F0, __CALLING_FUNCTION__);
    m_stgLengthFactors[0] = 1.0;
    m_stgLengthFactors[1] = 0.6;
    m_stgLengthFactors[2] = 1.5;

   /*! \page propertyPage1
      \section stgLengthFactors
      <code>ZFSFloat ZFSStrctrdBlck::m_stgLengthFactors </code>\n
      default = <code>1.0, 0.6, 1.5</code>\n \n
      The factor to scale the length scales\n
      in each coordinate direction with. For higher\n
      Reynolds number the values [1.0, 0.5, 1.4] \n
      produce better results.\n
      possible values are:
      <ul>
      <li>Float<3> > 0.0</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    if(ZFSContext::propertyExists("stgLengthFactors", m_blockId )){
      for(ZFSInt i = 0; i < 3; i++) {
        m_stgLengthFactors[i] = *(ZFSContext::getProperty("stgLengthFactors", m_blockId, __CALLING_FUNCTION__,
                                                          &m_stgLengthFactors[i], 3)->asFloat(i));
      }
    }

   /*! \page propertyPage1
      \section stgMaxNoEddies
      <code>ZFSFloat ZFSStrctrdBlck::m_stgMaxNoEddies </code>\n
      default = <code>200</code>\n \n
      Number of Eddies in the STG virtual box.\n
      possible values are:
      <ul>
      <li>Integer > 0</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgMaxNoEddies = 200;
    m_stgMaxNoEddies = *(ZFSContext::getProperty("stgMaxNoEddies", m_blockId, __CALLING_FUNCTION__,&m_stgMaxNoEddies)->asInt(0));

   /*! \page propertyPage1
      \section stgExple
      <code>ZFSFloat ZFSStrctrdBlck::m_stgExple </code>\n
      default = <code>0.5</code>\n \n
      Exponent of the STG LengthScale law.\n
      possible values are:
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgExple = 0.5;
    m_stgExple = *(ZFSContext::getProperty("stgExple", m_blockId, __CALLING_FUNCTION__,&m_stgExple)->asFloat(0));

    /*! \page propertyPage1
      \section stgEddieDistribution
      <code>ZFSFloat ZFSStrctrdBlck::m_stgEddieDistribution </code>\n
      default = <code>1.0</code>\n \n
      Shift die eddie distribution more to the wall\n
      or boundary layer edge.\n
      possible values are:
      <ul>
      <li>Floating point > 0.0</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgEddieDistribution = 1.0;
    if(ZFSContext::propertyExists("stgEddieDistribution", m_blockId )){
      m_stgEddieDistribution = *(ZFSContext::getProperty("stgEddieDistribution", m_blockId, __CALLING_FUNCTION__,&m_stgEddieDistribution)->asFloat(0));
    }

    /*! \page propertyPage1
      \section stgCreateNewEddies
      <code>ZFSFloat ZFSStrctrdBlck::m_stgCreateNewEddies </code>\n
      default = <code>0</code>\n \n
      Enforces the creation of all new eddies in STG virtual box\n
      or boundary layer edge.\n
      possible values are:
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgCreateNewEddies = false;
    if(ZFSContext::propertyExists("stgCreateNewEddies", m_blockId )){
      m_stgCreateNewEddies = *(ZFSContext::getProperty("stgCreateNewEddies", m_blockId, __CALLING_FUNCTION__, &m_stgCreateNewEddies )->asInt(0));
    }

    /*! \page propertyPage1
      \section stgInitialStartup
      <code>ZFSFloat ZFSStrctrdBlck::m_stgInitialStartup </code>\n
      default = <code>0</code>\n \n
      Initialize STG Method at Startup\n
      possible values are:
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgInitialStartup = false;
    m_stgInitialStartup = *(ZFSContext::getProperty("stgInitialStartup", m_blockId, __CALLING_FUNCTION__, &m_stgInitialStartup )->asInt(0));

    /*! \page propertyPage1
      \section stgEddieLengthScales
      <code>ZFSFloat ZFSStrctrdBlck::m_stgEddieLengthScales </code>\n
      default = <code>0</code>\n \n
      Connect length scales to eddies, not cells.\n
      possible values are:
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgEddieLengthScales = false;
    if(ZFSContext::propertyExists("stgEddieLengthScales", m_blockId )){
      m_stgEddieLengthScales = *(ZFSContext::getProperty("stgEddieLengthScales", m_blockId, __CALLING_FUNCTION__, &m_stgEddieLengthScales )->asInt(0));
    }

    /*! \page propertyPage1
      \section stgShapeFunction
      <code>ZFSFloat ZFSStrctrdBlck::m_stgFunction </code>\n
      default = <code>4</code>\n \n
      Shape function to be used in STG method.\n
      possible values are:
      <ul>
      <li>Integer >= 0</li>
      </ul>
      Keywords: <i>STG, STRCTRD</i>
    */
    m_stgShapeFunction = 4;
    if(ZFSContext::propertyExists("stgShapeFunction", m_blockId )){
      m_stgShapeFunction = *(ZFSContext::getProperty("stgShapeFunction", m_blockId, __CALLING_FUNCTION__, &m_stgShapeFunction )->asInt(0));
    }

    if(m_stgInitialStartup) {
      //activate the nut FQ field
      FQ->neededFQVariables[FQ->NU_T] = 1;
    }

    m_stgNoEddieProperties = 6;
    if(m_stgEddieLengthScales) {m_stgNoEddieProperties = 9;}
    zfsAlloc(m_stgEddies, m_stgMaxNoEddies, m_stgNoEddieProperties, "m_block->m_stgEddies", -F1, __CALLING_FUNCTION__);

    m_stgNoVariables = 20;
    ZFSInt noSTGCells = m_nCells[0]*m_nCells[1]*3;
    zfsAlloc(m_cells->stg_fq, m_stgNoVariables, noSTGCells, "m_cells->stg_fq", F0, __CALLING_FUNCTION__);

    m_stgBoxSize[0] = 0;
    m_stgBoxSize[1] = 0;
    m_stgBoxSize[2] = 0;

    zfs_log << "===========================================================" << endl
            << "                    STG PROPERTIES " << endl
            << "===========================================================" << endl
            << "Initial Start: " << m_stgInitialStartup << endl
            << "SubSup (Mixed subsonic/supersonic bc): " << m_stgSubSup << endl
            << "Supersonic BC: " << m_stgSupersonic << endl
            << "BLT 1,2,3: " << m_stgBLT1 << ", " << m_stgBLT2 << ", " << m_stgBLT3  << endl
            << "Delta0 inflow: " << m_stgDelta99Inflow << endl
            << "Length factors: " << m_stgLengthFactors[0] << ", "
            << m_stgLengthFactors[1] << ", "
            << m_stgLengthFactors[2] << endl
            << "Number of eddies: " << m_stgMaxNoEddies << endl
            << "Length scale exponent: " << m_stgExple << endl
            << "Eddie distribution: " << m_stgEddieDistribution << endl
            << "Create new eddies: " << m_stgCreateNewEddies << endl
            << "Eddie lengthscales: " << m_stgEddieLengthScales << endl
            << "Shape function: " << m_stgShapeFunction << endl
            << "Number of eddie properties: " << m_stgNoEddieProperties << endl
            << "Number of stg variables: " << m_stgNoVariables << endl
            << "===========================================================" << endl;

    switch(m_stgFace) {
    case 0:
    case 1:
      m_stgBoxSize[0] = m_nCells[0];
      m_stgBoxSize[1] = m_nCells[1];
      m_stgBoxSize[2] = 3;
      break;
    default:
      zfsTerm(1, __CALLING_FUNCTION__, "STG Method is not prepared for faces different than 0 or 1!");
    }
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::setProfileBCProperties()
{
  TRACE();

  ////////////////////////////////////////////////////////
  //////////////// BC 2600 ///////////////////////////////
  ////////////////////////////////////////////////////////

  m_bc2600IsActive = false;
  m_bc2600 = false;
  m_bc2600RootRank = -1;
  for(ZFSInt i=0; i < abs((ZFSInt)m_windowInfo->globalStrctrdBndryCndMaps.size()); ++i) {
    if(m_windowInfo->globalStrctrdBndryCndMaps[i]->BC == 2600) {
      m_bc2600IsActive = true;
    }
  }

  if(m_bc2600IsActive) {
    //now look if this domain contains parts of this bc
    ZFSInt localRank = 9999999;
    for(ZFSInt i=0; i < abs((ZFSInt)m_windowInfo->physicalBCMap.size()); ++i) {
      if(m_windowInfo->physicalBCMap[i]->BC == 2600) {
        m_bc2600 = true;
        localRank = domainId();
        break;
      }
    }

    MPI_Allreduce(&localRank, &m_bc2600RootRank, 1, MPI_INT, MPI_MIN, m_zfsStrctrdComm);
    
    //junoh
    if(m_bc2600) {
      MPI_Comm_rank(*m_commBC2600, &m_commBC2600MyRank);
    } else {
      m_commBC2600MyRank = -1;
    }

    /*! \page propertyPage1
      \section initialStartup2600
      <code>ZFSFloat ZFSStrctrdBlck::m_bc2600InitialStartup </code>\n
      default = <code>0</code>\n \n
      Trigger to indicate the initial start of the BC 2600\n
      load the value from the field into the BC field.\n
      possible values are:
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>BC2600, STRCTRD</i>
    */
    m_bc2600InitialStartup = false;
    m_bc2600InitialStartup = *(ZFSContext::getProperty("initialStartup2600", m_blockId, __CALLING_FUNCTION__,&m_bc2600InitialStartup)->asInt(0));

    zfsAlloc(m_bc2600noCells, nDim, "m_bc2600noCells", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_bc2600noActiveCells, nDim, "m_bc2600noCells", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_bc2600noOffsetCells, nDim, "m_bc2600noCells", 0, __CALLING_FUNCTION__);
    if(m_bc2600) {
      for(ZFSId dim=0; dim<nDim; dim++) {
        m_bc2600noOffsetCells[dim] = m_nOffsetCells[dim];
        m_bc2600noCells[dim] = m_nCells[dim];
        m_bc2600noActiveCells[dim] = m_bc2600noCells[dim]-2*m_noGhostLayers;
      }
      m_bc2600noOffsetCells[nDim-1] = 0;
      m_bc2600noCells[nDim-1] = m_noGhostLayers;
      m_bc2600noActiveCells[nDim-1] = m_noGhostLayers;
      ZFSInt noCellsBC = 1;
      for(ZFSId dim=0; dim<nDim; dim++) {
        noCellsBC *= m_bc2600noCells[dim];
      }
      zfsAlloc(m_bc2600Variables, m_maxNoVariables, noCellsBC, "m_bc2600Variables", -123.456, __CALLING_FUNCTION__);
    }
  }

  ////////////////////////////////////////////////////////
  //////////////// BC 2601 ///////////////////////////////
  ////////////////////////////////////////////////////////
  m_bc2601IsActive = false;
  m_bc2601 = false;
  for(ZFSInt i=0; i < abs((ZFSInt)m_windowInfo->globalStrctrdBndryCndMaps.size()); ++i) {
    if(m_windowInfo->globalStrctrdBndryCndMaps[i]->BC == 2601) {
      m_bc2601IsActive = true;
    }
  }

  if(m_bc2601IsActive) {

    /*! \page propertyPage1
      \section initialStartup2601
      <code>ZFSFloat ZFSStrctrdBlck::m_bc2601InitialStartup </code>\n
      default = <code>0</code>\n \n
      Trigger to indicate the initial start of the BC 2601\n
      load the value from the field into the BC field.\n
      possible values are:
      <ul>
      <li>true/false</li>
      </ul>
      Keywords: <i>BC2601, STRCTRD</i>
    */
    m_bc2601InitialStartup = false;
    m_bc2601InitialStartup = *(ZFSContext::getProperty("initialStartup2601", m_blockId, __CALLING_FUNCTION__,&m_bc2601InitialStartup)->asInt(0));

    /*! \page propertyPage1
      \section gammaEpsilon2601
      <code>ZFSFloat ZFSStrctrdBlck::m_bc2601GammaEpsilon </code>\n
      default = <code>0</code>\n \n
      Gamma Epsilon value for the BC2601\n
      Keywords: <i>BC2601, STRCTRD</i>
    */
    m_bc2601GammaEpsilon = 0.12;
    m_bc2601GammaEpsilon = *(ZFSContext::getProperty("gammaEpsilon2601", m_blockId, __CALLING_FUNCTION__,&m_bc2601GammaEpsilon)->asFloat());

    zfsAlloc(m_bc2601noCells, nDim, "m_bc2601noCells", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_bc2601noActiveCells, nDim, "m_bc2601noCells", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_bc2601noOffsetCells, nDim, "m_bc2601noCells", 0, __CALLING_FUNCTION__);
    if(m_bc2601) {
      for(ZFSId dim=0; dim<nDim; dim++) {
        m_bc2601noOffsetCells[dim] = m_nOffsetCells[dim];
        m_bc2601noCells[dim] = m_nCells[dim];
        m_bc2601noActiveCells[dim] = m_bc2601noCells[dim]-2*m_noGhostLayers;
      }
      m_bc2601noOffsetCells[nDim-2] = 0;
      m_bc2601noCells[nDim-2] = m_noGhostLayers;
      m_bc2601noActiveCells[nDim-2] = m_noGhostLayers;
      ZFSInt noCellsBC = 1;
      for(ZFSId dim=0; dim<nDim; dim++) {
        noCellsBC *= m_bc2601noCells[dim];
      }
      zfsAlloc(m_bc2601Variables, CV->noVariables, noCellsBC, "m_bc2601Variables", -123.456, __CALLING_FUNCTION__);
    }
  }
}


//! Set the zonal properties for RANS/LES calculations
/** \brief Set which zones are RANS and which are LES or if full LES or full RANS
 *  \author Marian Albers/Bhalaji
 *  \date Jan, 2015
 **/
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::setZonalProperties()
{
  TRACE();  

  ZFSBool fullRANS = false;
  m_zoneType = "LES";
  m_rans = false;

  /*! \page propertyPage1
    \section zonal
    <code>ZFSFloat ZFSStrctrdBlck::m_zonal </code>\n
    default = <code>0</code>\n \n
    Trigger a zonal computation.\n
    possible values are:
    <ul>
    <li>true/false</li>
    </ul>
    Keywords: <i>ZONAL, STRCTRD</i>
  */
  m_zonal = false;
  if(ZFSContext::propertyExists("zonal", m_blockId )){
    m_zonal = (ZFSBool)*(ZFSContext::getProperty("zonal", m_blockId, __CALLING_FUNCTION__, (ZFSInt*)NULL )->asInt(0));
  }

  /*! \page propertyPage1
    \section fullRANS
    <code>ZFSFloat ZFSStrctrdBlck::m_fullRANS </code>\n
    default = <code>0</code>\n \n
    Trigger a zonal computation.\n
    possible values are:
    <ul>
    <li>true/false</li>
    </ul>
    Keywords: <i>RANS, ZONAL, STRCTRD</i>
  */
  fullRANS = 0;
  if(ZFSContext::propertyExists("fullRANS", m_blockId )){
    fullRANS = (ZFSBool)*(ZFSContext::getProperty("fullRANS", 0, __CALLING_FUNCTION__, &fullRANS)->asInt(0));
  }

  if(m_zonal) {   //junoh
    //activate the nut FQ field
    FQ->neededFQVariables[FQ->NU_T] = 1;
    FQ->neededFQVariables[FQ->MU_T] = 1;
    FQ->neededFQVariables[FQ->WALLDISTANCE] = 1;
    FQ->neededFQVariables[FQ->AVG_RHO] = 1;
    FQ->neededFQVariables[FQ->AVG_U] = 1;
    FQ->neededFQVariables[FQ->AVG_V] = 1;
    FQ->neededFQVariables[FQ->AVG_W] = 1;
    FQ->neededFQVariables[FQ->AVG_P] = 1;


    // FQ->neededFQVariables[FQ->FLUC_U] = 1;
    // FQ->neededFQVariables[FQ->FLUC_V] = 1;
    // FQ->neededFQVariables[FQ->FLUC_W] = 1;
    // FQ->neededFQVariables[FQ->FLUC_UU] = 1;
    // FQ->neededFQVariables[FQ->FLUC_VV] = 1;
    // FQ->neededFQVariables[FQ->FLUC_WW] = 1;
    // FQ->neededFQVariables[FQ->FLUC_UV] = 1;
    // FQ->neededFQVariables[FQ->FLUC_UW] = 1;
    // FQ->neededFQVariables[FQ->FLUC_VW] = 1;
    
    // FQ->neededFQVariables[FQ->RECONST_NUT] = 1;
    // FQ->neededFQVariables[FQ->RECONST_NUTILDE] = 1;
    // FQ->neededFQVariables[FQ->NUTILDE] = 1;
    
      
    
    m_zonalExchangeInterval = 5;
    m_zonalExchangeInterval = *(ZFSContext::getProperty("zonalExchangeInterval", m_blockId, __CALLING_FUNCTION__, (ZFSId*) &m_zonalExchangeInterval)->asInt(0));
    
    m_zonalAveragingFactor = 128.0;
    m_zonalAveragingFactor = *(ZFSContext::getProperty("zonalAvgFactor", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) &m_zonalAveragingFactor)->asFloat(0));


    /*! \page propertyPage1
      \section noRansZones
      <code>ZFSFloat ZFSStrctrdBlck::m_noRansZones </code>\n
      default = <code>0</code>\n \n
      Number of zones (blocks) with RANS.\n
      possible values are:
      <ul>
      <li>Integer >= 0</li>
      </ul>
      Keywords: <i>RANS, ZONAL, STRCTRD</i>
    */
    ZFSInt noRansZones = 0;
    noRansZones = *(ZFSContext::getProperty("noRansZones", 0, __CALLING_FUNCTION__, &noRansZones)->asInt(0));
    ZFSIntScratchSpace ransZones(noRansZones, __CALLING_FUNCTION__, "ransZones");
    ransZones.fill(0);

    /*! \page propertyPage1
      \section ransZone
      <code>ZFSFloat ZFSStrctrdBlck::ransZones </code>\n
      default = <code>0</code>\n \n
      IDs of the RANS blocks.\n
      possible values are:
      <ul>
      <li>Integer >= 0</li>
      </ul>
      Keywords: <i>RANS, ZONAL, STRCTRD</i>
    */
    for(ZFSInt RANS =0;RANS<noRansZones;RANS++){
      ransZones[RANS] =*(ZFSContext::getProperty("ransZones", m_blockId, __CALLING_FUNCTION__, &ransZones[RANS] , noRansZones)->asInt(RANS));
    }

    //Find out if own partition is RANS
    for(ZFSInt RANS =0;RANS<noRansZones;RANS++){
      ZFSId inputBoxID = m_partition->outputBoxInfo[getBoxId(domainId())]->inputBoxID;
      if(inputBoxID == ransZones[RANS]) {
        m_zoneType = "RANS";
        m_rans = 1;
      }
    }

    MPI_Barrier(m_zfsStrctrdComm);

    //Count the zones
    ZFSInt NOZONES[2] = {0,0};
    ZFSInt NOZONESH[2] = {0,0};

    if(m_zoneType == "RANS") {
      NOZONES[0] = 1;
      NOZONES[1] = 0;
    } else {
      NOZONES[0] = 0;
      NOZONES[1] = 1;
    }

    MPI_Allreduce(NOZONES, NOZONESH, 2, MPI_INT, MPI_SUM, m_zfsStrctrdComm);

    //Give overview over RANS/LES zones
    if(domainId() == 0) {
      cout << "////////////////////////////////////////////////////////////////////////" << endl;
      cout << "No of RANS partitions: " << NOZONESH[0] << " , No of LES partitions: " << NOZONESH[1] << endl;
      cout << "////////////////////////////////////////////////////////////////////////" << endl;
    }
    zfs_log << "No of RANS partitions: " << NOZONESH[0] << " , No of LES partitions: " << NOZONESH[1] << endl;
  }

  if(fullRANS) {
    zfs_log << "Starting a full RANS computation" << endl;
    m_zoneType = "RANS";
    m_rans = 1;
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::allocateVariables(){
  //We will decide whether we will use RANS Variables or not depending on the block

  if(m_zoneType == "RANS"){
    //number of RANS equations
    //when we use another (two-equation) model
    //we need to alter this
    m_noRansEquations = 1;

    CV = new ZFSConservativeVariables<nDim>(m_noSpecies, m_noRansEquations);
    PV = new ZFSPrimitiveVariables<nDim>(m_noSpecies, m_noRansEquations );

    // FQ->neededFQVariables[FQ->WALLDISTANCE] = 1;
    // FQ->neededFQVariables[FQ->NU_T] = 1;
    // FQ->neededFQVariables[FQ->MU_T] = 1;
  } else {
    //we only have LES blocks
    m_noRansEquations = 0;
    CV = new ZFSConservativeVariables<nDim>(m_noSpecies );
    PV = new ZFSPrimitiveVariables<nDim>(m_noSpecies );
  }
  //junoh
  m_maxNoVariables = -1;
  MPI_Allreduce(&PV->noVariables, &m_maxNoVariables, 1, MPI_INT, MPI_MAX, m_zfsStrctrdComm);
  zfs_log << "Max number of variables: " << m_maxNoVariables << endl;
}


/** \brief Reset the right hand side to zero (floating point)
 * \author Pascal Meysonnat
 * \date 01.01.1010
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::resetRHS()
{
  //as only 1D array is used the ghostcells will also need to have a right hand side variables
  //if storage need to be saved go over the inner loops and remove storge of ghost rhs
  const ZFSId noVars = CV->noVariables;
  for( ZFSId cellId = 0; cellId < m_noStrctrdCells; cellId++ ) {
    for( ZFSId varId = 0;  varId < noVars;  varId++ ) {
      m_cells->rightHandSide[varId][ cellId ] = F0;
    }
  }
}

/** \brief Read in the grid file in serial/parallel
 * \author Pascal Meysonnat
 * \date 01.01.1010
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::readGrid(ZFSInt file_id)
{
  //this function reads in the grid !!!
  if(file_id<0) zfsTerm(1, __CALLING_FUNCTION__, "ERROR: grid File could not be opened :-(");

  ZFSInt start[3]; //max is three
  for (ZFSInt i=0; i<nDim;i++)
  {
    //get offset in file
    start[i]=m_partition->outputBoxInfo[getBoxId(domainId())]->offset[i];
  }
  //create the string to contain the datasetname in the file
  ZFSString sBlockName = "/block";
  stringstream dummy1;
  dummy1 << m_inputBlockId <<"/";
  sBlockName += dummy1.str();

  ZFSString varNames[]= {"x","y","z"};

  for(ZFSId dim=0; dim<nDim; dim++){
    io_read_ddataset_part1d1(file_id, sBlockName.c_str(), varNames[dim].c_str(), nDim, start, m_nActivePoints, m_coordinates[dim]);
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::computePV()
{
  //don't do this here, because it
  //will be executed before exchange() and
  //then wrong primitve values are in the GC
  //computePrimitiveVariables();
}

/** \brief Saves coordinates for whole grid without ghost points
 * \author Pascal Meysonnat
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::writeGrid() {
  ZFSId writeOutGhostLayers=1;
  //This function writes the grid only with or without ghost cells
  if(writeOutGhostLayers==1) {
    stringstream fileName;
    fileName << m_solutionOutput << "Grid" << globalTimeStep << m_outputFormat;

    ZFSInt file = io_openfile("hdf5", (fileName.str()).c_str(), "collective", m_zfsStrctrdComm);
    ZFSInt test[3] = {0,0,0};
    ZFSInt dummy= m_noInputBlocks;
    io_create_iattribute(file, "", "noBlocks", 1);
    io_write_iattribute1(file, "", "noBlocks", 1, &dummy);
    io_create_sattribute(file, "", "filetype", 6);
    io_write_sattribute1(file, "", "filetype", 6, "grid");
    io_create_sattribute(file, "", "gridType", 12);
    io_write_sattribute1(file, "", "gridType", 12, "structured");
    io_create_sattribute(file, "", "UID", m_uID.length());
    io_write_sattribute1(file, "", "UID", m_uID.length(), m_uID.c_str());
    //for moving grids we need the time in the header
    writePropertiesAsAttributes(file, "");
    ZFSInt allPoints[3];//not nDim because of compiler warning0
    for (ZFSId i=0; i<m_noInputBlocks; ++i){
      for(ZFSInt j =0; j<nDim; ++j){
        allPoints[j]=m_partition->inputBoxInfo[i]->DirLast[j]+1;//+2*m_noGhostLayers;
      }
      stringstream path;
      path << i;//m_inputBlockId;
      ZFSString blockPathStr = "block";
      blockPathStr += path.str();
      const char* blockPath = blockPathStr.c_str();

      io_create_ddataset(file, blockPath , "x", nDim,  allPoints);
      io_create_ddataset(file, blockPath , "y", nDim,  allPoints);
      if(nDim==3) {
        io_create_ddataset(file, blockPath , "z", nDim,  allPoints);
      }
    }

    stringstream path;
    path << m_inputBlockId;
    ZFSString blockPathStr = "block";
    blockPathStr += path.str();
    int start[3]={0,0,0};
    for (ZFSInt i=0; i<nDim;++i) {
      //get offset in file
      start[i]=m_partition->outputBoxInfo[getBoxId(domainId())]->offset[i];
    }

    for(ZFSId j=0; j<nDim; j++) {
      test[j]=m_partition->outputBoxInfo[getBoxId(domainId())]->DirLast[j];
    }

    ZFSInt ghostArray[3]={m_noGhostLayers, m_noGhostLayers, m_noGhostLayers};

    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "x", nDim,  test, start, ghostArray ,&m_coordinates[0][0]);
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "y", nDim,  test, start, ghostArray ,&m_coordinates[1][0]);

    if(nDim==3) {
      io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "z", nDim,  test, start, ghostArray ,&m_coordinates[2][0]);
    }

    io_closefile(file);
  }
}


/** \brief Saves coordinates for partitioned grid with ghost points
 * \author Pascal Meysonnat
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::writeGridPointsWithGhostPoints()
{
  TRACE();
  //first every process needs to create the datasets and the structure where the
  //data will be stored!
  cout << "writing the partitionedGrid.hdf5 File" << endl;
  const char* fileName = "partitionedGrid.hdf5";
  ZFSInt file = io_openfile("hdf5",fileName, "collective", m_zfsStrctrdComm);
  ZFSInt dummy = noDomains();
  io_create_iattribute(file, "", "noBlocks", 1 );
  io_write_iattribute1(file, "", "noBlocks", 1, &dummy);
  ZFSInt test[3] = {0,0,0};
  ZFSString gridVarNames[3]={"x","y","z"};
  for(ZFSInt i=0; i<noDomains(); i++) {
    //create datasets for the io library
    for(ZFSId j=0; j<nDim; j++) {
      test[j]=m_partition->outputBoxInfo[getBoxId(i)]->DirLast[j]+2*m_noGhostLayers;
    }
    stringstream path;
    path << getBoxId(i);
    ZFSString blockPathStr = "block";
    blockPathStr += path.str();
    const char* blockPath = blockPathStr.c_str();
    for(ZFSId dim=0; dim<nDim;++dim){
      io_create_ddataset(file, blockPath , gridVarNames[dim].c_str(),nDim,  test);
    }

    ////debugging only
    io_create_ddataset(file, blockPath , "pointId", nDim, test);
    io_create_ddataset(file, blockPath , "cpuId"  , nDim, test);
  }

  //write the values into the array so that we can visualize it
  ZFSInt offset[3]={0,0,0};
  for(ZFSId j=0; j<nDim; j++)
  {
    test[j]=m_partition->outputBoxInfo[getBoxId(domainId())]->DirLast[j]+2*m_noGhostLayers;
  }
  stringstream path;
  path << getBoxId(domainId());
  ZFSString blockPathStr = "block";
  blockPathStr += path.str();
  for(ZFSId dim=0; dim<nDim;++dim){
    io_write_ddataset_part(file, blockPathStr.c_str(), gridVarNames[dim].c_str(), nDim, test, offset, &m_coordinates[dim][0]);
  }

  io_write_ddataset_part(file, blockPathStr.c_str(), "pointId", nDim, test, offset, &pointProperties[0][0]);
  io_write_ddataset_part(file, blockPathStr.c_str(), "cpuId", nDim, test, offset, &pointProperties[1][0]);
  io_closefile(file);
}



template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::computeConservativeVariables()
{
  zfs_log << "we got into the general formulation but should be in the other one" << endl;
  const ZFSFloat FgammaMinusOne = F1 / (m_gamma - 1.0);
  ZFSFloat** const RESTRICT cvars = m_cells->variables;
  ZFSFloat** const RESTRICT pvars = m_cells->pvariables;
  ZFSFloat rhoVelocity[3] = {0.0, 0.0, 0.0};

  for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
    ZFSFloat velPOW2 = F0;
    //copy the density
    const ZFSFloat rho=pvars[PV->RHO][cellId];
    //compute the rho * velocity
    for(ZFSId i=0; i<nDim; i++) {
      rhoVelocity[i] = pvars[PV->VV[i]][cellId]*rho;
      velPOW2 += POW2(pvars[PV->VV[i]][cellId]);
    }

    //regular conservative variables
    cvars[CV->RHO][cellId] = rho;
    for(ZFSId i=0; i<nDim; i++) {
      cvars[CV->RHO_VV[i]][cellId] = rhoVelocity[i];
    }

    cvars[ CV->RHO_E ][cellId] =  pvars[ PV->P ][cellId] * FgammaMinusOne + F1B2*rho*velPOW2;

    //rans
    for(ZFSId ransVar=0; ransVar < m_noRansEquations; ransVar++) {
      // pvars[PV->RANS_VAR[ransVar]][cellId] = zfsMAX(pvars[PV->RANS_VAR[ransVar]][cellId], F0);
      cvars[CV->RANS_VAR[ransVar]][cellId] = zfsMAX(pvars[PV->RANS_VAR[ransVar]][cellId], F0)*rho;
    }

    //species
    for( ZFSId s=0; s<m_noSpecies; s++ ) {
      cvars[ CV->RHO_Y[ s ] ][cellId] = pvars[ PV->Y[ s ] ][cellId] * rho;
    }
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::saveVarToPrimitive(ZFSId cellId, ZFSId varId, ZFSFloat var) {
  m_cells->pvariables[varId][cellId] = var;
}


/**
 * @author Marian Albers, Pascal Meysonnat
 *    @date 01.01.1010
 *
 *    Saves variables in a plane 
 *    Plane start, end indices and normals
 *    can be given in the property file
 */
template <ZFSInt nDim>
template <ZFSBool primitiveOutput>
void ZFSStrctrdBlck<nDim>::saveOutputPlanes(){

  //create a file
  //a) all planes for each time in one file
  //b) all planes in one file and one file per output
  stringstream filename;
  filename << m_planeOutputDir << "planeOutput" << m_outputIterationNumber << m_outputFormat;  

  zfs_log << "Writing files to " << filename.str() << endl;
  ZFSInt file =io_openfile("hdf5", (filename.str()).c_str(), "collective", m_zfsStrctrdComm);

  ZFSString gridFile="../"+m_gridInputFileName;
  writeHeaderAttributes(file, gridFile, "planes");

  io_create_iattribute(file,"", "noPlanes", 1);
  io_write_iattribute1(file,"", "noPlanes", 1, &m_noPlaneOutput);

  //create datasets
  ZFSInt dataSize[2];
  ZFSInt dataOffset[2];

  stringstream timeStepPath;
  timeStepPath << "t" << m_outputIterationNumber; 
 						   

  writePropertiesAsAttributes(file,timeStepPath.str() );

  for(ZFSId i=0; i<m_noInputBlocks; ++i){
    for(ZFSId p=0; p<m_noPlaneOutput; ++p){
      if(m_planeBlock[p]==i){
        //create a dataset for the block
        int temp=0;
        for(int dim=0; dim<nDim; ++dim){
          if(m_planeNormal[p]!=dim){
            dataSize[temp]=m_totalGridBlockDim[i][dim]-1;
            dataOffset[temp]=m_nOffsetCells[dim];
            ++temp;
          }
        }

        if(m_planeOffset[p] > m_totalGridBlockDim[i][m_planeNormal[p]]){
          cout << m_planeOffset[p] << " " << m_totalGridBlockDim[i][m_planeNormal[p]] << endl;
          zfsTerm(1, __CALLING_FUNCTION__, "planeOffset is bigger than the dimension of the block");
        }

        //specific attributes only for the plane
        stringstream pathName;
        pathName << "t" << m_outputIterationNumber << "/plane" << p;
        io_create_iattribute(file,(pathName.str()).c_str(), "normal",1);
        io_write_iattribute1(file,(pathName.str()).c_str(), "normal",1, &m_planeNormal[p]);
        io_create_iattribute(file,(pathName.str()).c_str(), "blockId",1);
        io_write_iattribute1(file,(pathName.str()).c_str(), "blockId",1, &i);
        io_create_iattribute(file,(pathName.str()).c_str(), "index",1);
        io_write_iattribute1(file,(pathName.str()).c_str(), "index",1, &m_planeOffset[p]);
        int hasCoordinates=0;
        io_create_iattribute(file, (pathName.str()).c_str(), "hasCoordinates",1);

        /////////////////////////////
        /////create the datasets/////
        /////////////////////////////
        //-->primitive/conservative variables
        if(primitiveOutput){
          for(ZFSId v=0; v< m_maxNoVariables; ++v){ //junoh
            io_create_ddataset(file, (pathName.str()).c_str(), m_pvariableNames[v].c_str(), 2, dataSize);
          }
        }else{
          for(ZFSId v=0; v< CV->noVariables; ++v){
            io_create_ddataset(file, (pathName.str()).c_str(), m_variableNames[v].c_str(), 2, dataSize);
          }
        }

        //-->fq field
        for(ZFSId v=0; v<FQ->noFQVariables; ++v){
          if(FQ->fqWriteOutputPlanes[v]) {
            io_create_ddataset(file, (pathName.str()).c_str() , FQ->fqNames[v].c_str(), 2, dataSize);
          }
        }
        //-->plane coordinates
        if(m_planeWriteCoordinates) {
          hasCoordinates=1;
          io_create_ddataset(file, (pathName.str()).c_str(), "x", 2, dataSize);
          io_create_ddataset(file, (pathName.str()).c_str(), "y", 2, dataSize);
          io_create_ddataset(file, (pathName.str()).c_str(), "z", 2, dataSize);
        }
        //write output to check if coordinates are contained within the variable list
        io_write_iattribute1(file,(pathName.str()).c_str(), "hasCoordinates",1, &hasCoordinates);
      }
    }
  }

  for(ZFSId p=0; p<m_noPlaneOutput; ++p){
    const ZFSInt normal = m_planeNormal[p];
    //check if the plane is contained your inputblockId
    if(m_planeBlock[p]==m_inputBlockId && ( m_nOffsetCells[normal]  <= m_planeOffset[p] && m_planeOffset[p]<m_nOffsetCells[normal] + m_nActiveCells[normal] ) && m_planeOffset[p] != -1){ //the plane is contained
      //get the size of the plane!!!
      ZFSId temp=0;
      for(ZFSId dim=0; dim<nDim; ++dim){
        if(normal!=dim){
          dataSize[temp]=m_nActiveCells[dim];
          dataOffset[temp]=m_nOffsetCells[dim];
          ++temp;
        }
      }
      
      stringstream pathName;
      pathName << "t" << m_outputIterationNumber << "/plane" << p;
      ZFSInt totalLocalSize =1;
      for(int dim=0; dim<nDim-1; ++dim){totalLocalSize*=dataSize[dim];}
 
      ZFSInt noFields = CV->noVariables;
      if(FQ->noFQPlaneOutput>0) noFields+=FQ->noFQPlaneOutput;
      if(m_planeWriteCoordinates) noFields += nDim;
      ZFSFloatScratchSpace localPlaneVar(noFields*totalLocalSize, __CALLING_FUNCTION__, "local Plane Variables");
      //determine the local offset
      ZFSInt localOffset=m_planeOffset[p]-m_nOffsetCells[normal];
      //create one small section to write out everything in every direction independent of the normal
      ZFSId cellId =0, localId=0, offset=0;
      ZFSId *reali=nullptr, *realj=nullptr, *realk=nullptr;
      ZFSInt ii=0, jj=0;
      ZFSInt endCellIndexJ[3]={m_nCells[1], m_nCells[0], m_nCells[0]};
      ZFSInt endCellIndexI[3]={m_nCells[2], m_nCells[2], m_nCells[1]};
      ZFSInt activeCells[3]={m_nActiveCells[2], m_nActiveCells[2], m_nActiveCells[1]};
      ZFSInt kkOffset=localOffset + m_noGhostLayers;
      //assign to pointer to the right index
      switch(normal){
      case 0 :{
        reali=&ii;
        realj=&jj;
        realk=&kkOffset;
        break;
      }
      case 1:{
        reali=&ii;
        realj=&kkOffset;
        realk=&jj;
        break;
      }
      case 2:{
        reali=&kkOffset;
        realj=&ii;
        realk=&jj;
        break;
      }
      default: zfsTerm(1, __CALLING_FUNCTION__, "no such plane normal");
      }
      //careful the index i is only one direction, j the other, independend of the real direction
      //======================== GENERAL APPROACH ===============================
      for(ZFSId var=0; var<m_maxNoVariables; ++var){ //junoh
        for(jj=m_noGhostLayers; jj<endCellIndexJ[normal]-m_noGhostLayers; ++jj){
          for(ii=m_noGhostLayers;ii<endCellIndexI[normal]-m_noGhostLayers; ++ii){
            cellId=*reali+(*realj+*realk*m_nCells[1])*m_nCells[2];
            localId=var*totalLocalSize+(ii-m_noGhostLayers+(jj-m_noGhostLayers)*activeCells[normal]);
            if(primitiveOutput){
              localPlaneVar[localId]= m_cells->pvariables[var][cellId];
            }else{
              localPlaneVar[localId]= m_cells->variables[var][cellId];
            }                      
          }
        }
      }
      offset=CV->noVariables;
      for(ZFSId v=0; v<FQ->noFQVariables; ++v){
        if(FQ->fqWriteOutputPlanes[v]) {
          for(jj=m_noGhostLayers; jj<endCellIndexJ[normal]-m_noGhostLayers; ++jj){
            for(ii=m_noGhostLayers;ii<endCellIndexI[normal]-m_noGhostLayers; ++ii){
              cellId=*reali+(*realj+*realk*m_nCells[1])*m_nCells[2];
              localId= (offset)*totalLocalSize+(ii-m_noGhostLayers+(jj-m_noGhostLayers)*activeCells[normal]);
              localPlaneVar[localId]= m_cells->fq[v][cellId];
            }
          } 
          ++offset;
        }
      }   
      if(m_planeWriteCoordinates){
        for(jj=m_noGhostLayers; jj<endCellIndexJ[normal]-m_noGhostLayers; ++jj){
          for(ii=m_noGhostLayers;ii<endCellIndexI[normal]-m_noGhostLayers; ++ii){
            cellId=*reali+(*realj+*realk*m_nCells[1])*m_nCells[2];
            for(ZFSId dim=0; dim<nDim; dim++) {
              localId= (offset+dim)*totalLocalSize+(ii-m_noGhostLayers+(jj-m_noGhostLayers)*activeCells[normal]);
              localPlaneVar[localId]= m_cells->coordinates[dim][cellId];
            }
          }
        }
        offset += nDim;
      }
      //======================== GENERAL APPROACH END ===============================
      offset = 0;
      /////////////////////////////
      /////write out the data!/////
      /////////////////////////////
      //-->/primitive/conservative variables
      if(primitiveOutput){
        for(ZFSId v=0; v<m_maxNoVariables; ++v){ //junoh
          io_write_ddataset_part(file, (pathName.str()).c_str(), m_pvariableNames[v].c_str(), nDim-1, dataSize, dataOffset, &localPlaneVar[v*totalLocalSize]);
        }
      }else{
        for(ZFSId v=0; v<CV->noVariables; ++v){
          io_write_ddataset_part(file, (pathName.str()).c_str(), m_variableNames[v].c_str(), nDim-1, dataSize, dataOffset, &localPlaneVar[v*totalLocalSize]);
        }
      }
      offset=CV->noVariables;
      //-->fq field
      for(ZFSId v=0; v<FQ->noFQVariables; ++v){
        if(FQ->fqWriteOutputPlanes[v]) {
          io_write_ddataset_part(file, (pathName.str()).c_str(), FQ->fqNames[v].c_str(), nDim-1, dataSize, dataOffset, &localPlaneVar[offset*totalLocalSize]);
          ++offset;
        }
      }
      //-->plane coordinates
      if(m_planeWriteCoordinates) {
        io_write_ddataset_part(file, (pathName.str()).c_str(), "x", nDim-1, dataSize, dataOffset, &localPlaneVar[offset*totalLocalSize]);
        io_write_ddataset_part(file, (pathName.str()).c_str(), "y", nDim-1, dataSize, dataOffset, &localPlaneVar[(offset+1)*totalLocalSize]);
       io_write_ddataset_part(file, (pathName.str()).c_str(), "z", nDim-1, dataSize, dataOffset, &localPlaneVar[(offset+2)*totalLocalSize]);
       offset += nDim;
      }
    }else{//write out nothing as plane is not contained
      stringstream pathName;
      pathName << "t" << m_outputIterationNumber << "/plane" << p;
      for(int dim=0; dim<nDim-1; ++dim){
        dataSize[dim]=0;
        dataOffset[dim]=0;
      }
      if(primitiveOutput){
        for(ZFSId v=0; v< m_maxNoVariables; ++v){ //junoh
          io_write_ddataset_part(file, (pathName.str()).c_str(),m_pvariableNames[v].c_str(),2, dataSize, dataOffset, NULL);
        }
      }else{
        for(ZFSId v=0; v< CV->noVariables; ++v){
          io_write_ddataset_part(file, (pathName.str()).c_str(),m_variableNames[v].c_str(),2, dataSize, dataOffset, NULL);
        }
      }
      for(ZFSId v=0; v<FQ->noFQVariables; ++v){
        if(FQ->fqWriteOutputPlanes[v]) {
          io_write_ddataset_part(file, (pathName.str()).c_str(), FQ->fqNames[v].c_str(),2, dataSize, dataOffset, NULL);
        }
      }

      if(m_planeWriteCoordinates){
        io_write_ddataset_part(file, (pathName.str()).c_str(), "x",2, dataSize, dataOffset, NULL);
        io_write_ddataset_part(file, (pathName.str()).c_str(), "y",2, dataSize, dataOffset, NULL);
        io_write_ddataset_part(file, (pathName.str()).c_str(), "z",2, dataSize, dataOffset, NULL);
      }
    }
  }

  io_closefile(file);
}

template void ZFSStrctrdBlck<2>::saveOutputPlanes<true>();
template void ZFSStrctrdBlck<2>::saveOutputPlanes<false>();
template void ZFSStrctrdBlck<3>::saveOutputPlanes<true>();
template void ZFSStrctrdBlck<3>::saveOutputPlanes<false>();




/**
 *
 * @author Marian Albers, Pascal Meysonnat (modifications)
 *    @date:  Nov 1, 2015,
 *    Saves variables of a given box instead
 *    of whole domain. Box start and end indices
 *    can be given in the property file
 */
template <ZFSInt nDim>
template <ZFSBool primitiveOutput>
void ZFSStrctrdBlck<nDim>::saveOutputBoxes(){
  //create a file
  //a) all boxes for each time in one file
  //b) all boxes in one file and one file per output
  stringstream filename;
  filename << m_boxOutputDir << "boxOutput" << m_outputIterationNumber << m_outputFormat;

  ZFSInt file =io_openfile("hdf5", (filename.str()).c_str(), "collective", m_zfsStrctrdComm);

  ZFSString gridFileName= "../"+m_gridInputFileName;
  writeHeaderAttributes(file, gridFileName, "boxes");

  io_create_iattribute(file,"", "noBoxes", 1);
  io_write_iattribute1(file,"", "noBoxes", 1, &m_noBoxOutput);

  //create datasets
  ZFSInt localBoxSize[3];
  ZFSInt localBoxOffset[3];
  ZFSInt globalBoxSize[3];
  ZFSInt localDomainBoxOffset[3];

  stringstream timeStepPath;
  timeStepPath << "t" << m_outputIterationNumber;

  writePropertiesAsAttributes(file, timeStepPath.str());

  for(ZFSId i=0; i<m_noInputBlocks; ++i){
    for(ZFSId b=0; b<m_noBoxOutput; ++b){
      if(m_boxBlock[b]==i){
        //create a dataset for the block

        for(int dim=0; dim<nDim; ++dim){
          globalBoxSize[dim]=m_boxSize[b][dim];
        }

        stringstream pathName;
        pathName << "t" << m_outputIterationNumber <<  "/box" << b;

        io_create_iattribute(file,(pathName.str()).c_str(), "offseti",1);
        io_write_iattribute1(file,(pathName.str()).c_str(), "offseti",1, &m_boxOffset[b][2]);
        io_create_iattribute(file,(pathName.str()).c_str(), "offsetj",1);
        io_write_iattribute1(file,(pathName.str()).c_str(), "offsetj",1, &m_boxOffset[b][1]);
        io_create_iattribute(file,(pathName.str()).c_str(), "offsetk",1);
        io_write_iattribute1(file,(pathName.str()).c_str(), "offsetk",1, &m_boxOffset[b][0]);

        io_create_iattribute(file,(pathName.str()).c_str(), "sizei",1);
        io_write_iattribute1(file,(pathName.str()).c_str(), "sizei",1, &m_boxSize[b][2]);
        io_create_iattribute(file,(pathName.str()).c_str(), "sizej",1);
        io_write_iattribute1(file,(pathName.str()).c_str(), "sizej",1, &m_boxSize[b][1]);
        io_create_iattribute(file,(pathName.str()).c_str(), "sizek",1);
        io_write_iattribute1(file,(pathName.str()).c_str(), "sizek",1, &m_boxSize[b][0]);

        io_create_iattribute(file,(pathName.str()).c_str(), "blockId",1);
        io_write_iattribute1(file,(pathName.str()).c_str(), "blockId",1, &i);

        int hasCoordinates=0;
        io_create_iattribute(file, (pathName.str()).c_str(), "hasCoordinates",1);

        //create primitive/conservative variable datasets
        if(primitiveOutput){
          for(ZFSId v=0; v<m_maxNoVariables; v++){io_create_ddataset(file, (pathName.str()).c_str(), m_pvariableNames[v].c_str(), 3, globalBoxSize);}  //junoh
        }else{
          for(ZFSId v=0; v<CV->noVariables; v++){io_create_ddataset(file, (pathName.str()).c_str(), m_variableNames[v].c_str(), 3, globalBoxSize);}
        }
        //create datasets for fq-field
        for(ZFSId v=0; v<FQ->noFQVariables; v++){
          if(FQ->fqWriteOutputBoxes[v]) {
            io_create_ddataset(file, (pathName.str()).c_str(), FQ->fqNames[v].c_str(), 3, globalBoxSize);
          }
        }
        //create datasets for the variables
        if(m_boxWriteCoordinates) {
          hasCoordinates=1;
          io_create_ddataset(file, (pathName.str()).c_str(), "x", 3, globalBoxSize);
          io_create_ddataset(file, (pathName.str()).c_str(), "y", 3, globalBoxSize);
          io_create_ddataset(file, (pathName.str()).c_str(), "z", 3, globalBoxSize);
        }
        //write output to check if coordinates are contained within the variable list
        io_write_iattribute1(file,(pathName.str()).c_str(), "hasCoordinates",1, &hasCoordinates);
      }
    }
  }

  for(ZFSId b=0; b<m_noBoxOutput; ++b){
    //check if the box is contained your inputblockId
    if(m_boxBlock[b]==m_inputBlockId
       && (( m_nOffsetCells[2]  <= m_boxOffset[b][2] && m_boxOffset[b][2] < m_nOffsetCells[2] + m_nActiveCells[2] ) || ( m_boxOffset[b][2] <= m_nOffsetCells[2] && m_nOffsetCells[2] < m_boxOffset[b][2] + m_boxSize[b][2] ))
       && (( m_nOffsetCells[1]  <= m_boxOffset[b][1] && m_boxOffset[b][1] < m_nOffsetCells[1] + m_nActiveCells[1] ) || ( m_boxOffset[b][1] <= m_nOffsetCells[1] && m_nOffsetCells[1] < m_boxOffset[b][1] + m_boxSize[b][1] ))
       && (( m_nOffsetCells[0]  <= m_boxOffset[b][0] && m_boxOffset[b][0] < m_nOffsetCells[0] + m_nActiveCells[0] ) || ( m_boxOffset[b][0] <= m_nOffsetCells[0] && m_nOffsetCells[0] < m_boxOffset[b][0] + m_boxSize[b][0] ))){ //the box is contained
      //get the size of the box!!!

      for(ZFSId dim=0; dim<nDim; ++dim){
        if( m_nOffsetCells[dim]  <= m_boxOffset[b][dim] && m_boxOffset[b][dim] + m_boxSize[b][dim] < m_nOffsetCells[dim] + m_nActiveCells[dim]) {
          localBoxSize[dim]=m_boxSize[b][dim];
          localBoxOffset[dim]=0;
          localDomainBoxOffset[dim]=m_boxOffset[b][dim]-m_nOffsetCells[dim];
        } else if ( m_nOffsetCells[dim]  <= m_boxOffset[b][dim]) {
          localBoxSize[dim]= (m_nOffsetCells[dim]+m_nActiveCells[dim]) - m_boxOffset[b][dim];
          localBoxOffset[dim]=0;
          localDomainBoxOffset[dim]=m_boxOffset[b][dim]-m_nOffsetCells[dim];
        } else if (  m_boxOffset[b][dim] <= m_nOffsetCells[dim]  && m_nOffsetCells[dim]+m_nActiveCells[dim] < m_boxOffset[b][dim] + m_boxSize[b][dim] ) {
          localBoxSize[dim]= m_nActiveCells[dim];
          localBoxOffset[dim]=m_nOffsetCells[dim]-m_boxOffset[b][dim];
          localDomainBoxOffset[dim]=0;
        } else {
          localBoxSize[dim]= (m_boxOffset[b][dim]+m_boxSize[b][dim]) - m_nOffsetCells[dim] ;
          localBoxOffset[dim]=m_nOffsetCells[dim]-m_boxOffset[b][dim];
          localDomainBoxOffset[dim]=0;
        }
      }

      stringstream pathName;
      pathName << "t" << m_outputIterationNumber << "/box" << b;
      ZFSInt totalLocalSize =1;
      for(int dim=0; dim<nDim; ++dim){totalLocalSize*=localBoxSize[dim];}
      ZFSInt noFields = m_maxNoVariables; //junoh
      if(FQ->noFQBoxOutput>0) noFields += FQ->noFQBoxOutput;
      if(m_boxWriteCoordinates) noFields += nDim;


      ZFSFloatScratchSpace localBoxVar(noFields*totalLocalSize, __CALLING_FUNCTION__, "local Box Variables");

      ZFSId cellId =0;
      ZFSId localId = 0;
      ZFSId offset = 0;

      ZFSFloat noVars=CV->noVariables;
      if(primitiveOutput){noVars=m_maxNoVariables;} //junoh
      for(ZFSId var=0; var<noVars; ++var){
        for(ZFSId k=m_noGhostLayers+localDomainBoxOffset[0]; k<m_noGhostLayers+localDomainBoxOffset[0]+localBoxSize[0]; ++k){
          for(ZFSId j=m_noGhostLayers+localDomainBoxOffset[1]; j<m_noGhostLayers+localDomainBoxOffset[1]+localBoxSize[1]; ++j){
            for(ZFSId i=m_noGhostLayers+localDomainBoxOffset[2]; i<m_noGhostLayers+localDomainBoxOffset[2]+localBoxSize[2]; ++i){
              cellId=i+(j+k*m_nCells[1])*m_nCells[2];
              ZFSId boxI = i-m_noGhostLayers-localDomainBoxOffset[2];
              ZFSId boxJ = j-m_noGhostLayers-localDomainBoxOffset[1];
              ZFSId boxK = k-m_noGhostLayers-localDomainBoxOffset[0];
              localId = var*totalLocalSize+(boxI + (boxJ + boxK*localBoxSize[1])*localBoxSize[2]);
              if(primitiveOutput){localBoxVar[localId]= m_cells->pvariables[var][cellId];
              }else{
                localBoxVar[localId]= m_cells->variables[var][cellId];
              }
            }
          }
        }
      }
      if(primitiveOutput){offset += m_maxNoVariables; //junoh
      }else{offset +=CV->noVariables;}
      
      if(FQ->noFQBoxOutput>0) {
        for(ZFSId v=0; v<FQ->noFQVariables; ++v){
          if(FQ->fqWriteOutputBoxes[v]) {
            for(ZFSId k=m_noGhostLayers+localDomainBoxOffset[0]; k<m_noGhostLayers+localDomainBoxOffset[0]+localBoxSize[0]; ++k){
              for(ZFSId j=m_noGhostLayers+localDomainBoxOffset[1]; j<m_noGhostLayers+localDomainBoxOffset[1]+localBoxSize[1]; ++j){
                for(ZFSId i=m_noGhostLayers+localDomainBoxOffset[2]; i<m_noGhostLayers+localDomainBoxOffset[2]+localBoxSize[2]; ++i){
                  cellId=i+(j+k*m_nCells[1])*m_nCells[2];
                  ZFSId boxI = i-m_noGhostLayers-localDomainBoxOffset[2];
                  ZFSId boxJ = j-m_noGhostLayers-localDomainBoxOffset[1];
                  ZFSId boxK = k-m_noGhostLayers-localDomainBoxOffset[0];
                  localId = (offset)*totalLocalSize+(boxI + (boxJ + boxK*localBoxSize[1])*localBoxSize[2]);
                  localBoxVar[localId]= m_cells->fq[v][cellId];
                }
              }
            }
            ++offset;
          }
        }
      }

      if(m_boxWriteCoordinates) {
        for(ZFSId dim=0; dim<nDim; dim++){
          for(ZFSId k=m_noGhostLayers+localDomainBoxOffset[0]; k<m_noGhostLayers+localDomainBoxOffset[0]+localBoxSize[0]; ++k){
            for(ZFSId j=m_noGhostLayers+localDomainBoxOffset[1]; j<m_noGhostLayers+localDomainBoxOffset[1]+localBoxSize[1]; ++j){
              for(ZFSId i=m_noGhostLayers+localDomainBoxOffset[2]; i<m_noGhostLayers+localDomainBoxOffset[2]+localBoxSize[2]; ++i){
                cellId=i+(j+k*m_nCells[1])*m_nCells[2];
                ZFSId boxI = i-m_noGhostLayers-localDomainBoxOffset[2];
                ZFSId boxJ = j-m_noGhostLayers-localDomainBoxOffset[1];
                ZFSId boxK = k-m_noGhostLayers-localDomainBoxOffset[0];
                localId = (offset+dim)*totalLocalSize+(boxI + (boxJ + boxK*localBoxSize[1])*localBoxSize[2]);
                localBoxVar[localId]= m_cells->coordinates[dim][cellId];
              }
            }
          }
        }
        offset += nDim;
      }
      ///////////////////////////////////
      ///////write out the data!/////////
      ///////////////////////////////////
      //--> /primitive/conservative variables
      if(primitiveOutput){
        for(ZFSId v=0; v<m_maxNoVariables; v++){ //junoh
          io_write_ddataset_part(file, (pathName.str()).c_str(), m_pvariableNames[v].c_str(), nDim, localBoxSize, localBoxOffset, &localBoxVar[v*totalLocalSize]);
        }
      }else{
        for(ZFSId v=0; v<CV->noVariables; v++){
          io_write_ddataset_part(file, (pathName.str()).c_str(), m_variableNames[v].c_str(), nDim, localBoxSize, localBoxOffset, &localBoxVar[v*totalLocalSize]);
        }
      }
      offset=m_maxNoVariables; //junoh
      //--> fq field
      for(ZFSId v=0; v<FQ->noFQVariables; ++v){
        if(FQ->fqWriteOutputBoxes[v]) {
          io_write_ddataset_part(file, (pathName.str()).c_str(), FQ->fqNames[v].c_str(), nDim, localBoxSize, localBoxOffset, &localBoxVar[(offset)*totalLocalSize]);
          offset++;
        }
      }

      if(m_boxWriteCoordinates) {
        io_write_ddataset_part(file, (pathName.str()).c_str(), "x", nDim, localBoxSize, localBoxOffset, &localBoxVar[(offset+0)*totalLocalSize]);
        io_write_ddataset_part(file, (pathName.str()).c_str(), "y", nDim, localBoxSize, localBoxOffset, &localBoxVar[(offset+1)*totalLocalSize]);
        io_write_ddataset_part(file, (pathName.str()).c_str(), "z", nDim, localBoxSize, localBoxOffset, &localBoxVar[(offset+2)*totalLocalSize]);
        offset += nDim;
      }

    }else{//write out nothing as box is not contained
      stringstream pathName;
      pathName << "t" << m_outputIterationNumber << "/box" << b;
      for(int dim=0; dim<nDim; ++dim){
        localBoxSize[dim]=0;
        localBoxOffset[dim]=0;
      }
      if(primitiveOutput){
        for(ZFSId v=0; v<m_maxNoVariables; ++v){io_write_ddataset_part(file, (pathName.str()).c_str(), m_pvariableNames[v].c_str(), nDim, localBoxSize, localBoxOffset, NULL);} //junoh
      }else{
        for(ZFSId v=0; v<CV->noVariables; ++v){io_write_ddataset_part(file, (pathName.str()).c_str(), m_variableNames[v].c_str(), nDim, localBoxSize, localBoxOffset, NULL);
        }
      }
      for(ZFSId v=0; v<FQ->noFQVariables; ++v){
        if(FQ->fqWriteOutputBoxes[v]) {
          io_write_ddataset_part(file, (pathName.str()).c_str(), FQ->fqNames[v].c_str(),nDim, localBoxSize, localBoxOffset, NULL);
        }
      }

      if(m_boxWriteCoordinates){
        io_write_ddataset_part(file, (pathName.str()).c_str(), "x",nDim, localBoxSize, localBoxOffset, NULL);
        io_write_ddataset_part(file, (pathName.str()).c_str(), "y",nDim, localBoxSize, localBoxOffset, NULL);
        io_write_ddataset_part(file, (pathName.str()).c_str(), "z",nDim, localBoxSize, localBoxOffset, NULL);
      }
    }
  }
  io_closefile(file);
}

template void ZFSStrctrdBlck<2>::saveOutputBoxes<true>();
template void ZFSStrctrdBlck<2>::saveOutputBoxes<false>();
template void ZFSStrctrdBlck<3>::saveOutputBoxes<true>();
template void ZFSStrctrdBlck<3>::saveOutputBoxes<false>();



template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::writeHeaderAttributes(ZFSId fileId, ZFSString gridFile, ZFSString fileType){
  io_create_sattribute(fileId, "", "gridFile", 18);
  io_write_sattribute1(fileId, "", "gridFile", 18, gridFile.c_str());

  io_create_sattribute(fileId, "", "UID", m_uID.length());
  io_write_sattribute1(fileId, "", "UID", m_uID.length(), m_uID.c_str());

  io_create_sattribute(fileId,"","filetype",9);
  io_write_sattribute1(fileId,"","filetype",9, fileType.c_str());

  io_create_sattribute(fileId,"","blockType",18);
  io_write_sattribute1(fileId,"","blockType",18, "ZFS_STRUCTURED_FV");

  io_create_iattribute(fileId,"", "zonal", 1);
  io_write_iattribute1(fileId,"", "zonal", 1, &m_zonal);

  io_create_iattribute(fileId,"", "noBlocks", 1);
  io_write_iattribute1(fileId,"", "noBlocks", 1, &m_noInputBlocks);

  io_create_iattribute(fileId,"", "primitiveOutput",1);
  io_write_iattribute1(fileId,"", "primitiveOutput",1, &m_primitiveOutput);
}



template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::writePropertiesAsAttributes(ZFSId fileId, ZFSString path){
  io_create_dattribute(fileId, path.c_str(), "Ma",1);
  io_write_dattribute1(fileId, path.c_str(), "Ma",1, &m_Ma);

  io_create_dattribute(fileId, path.c_str(), "Re",1);
  io_write_dattribute1(fileId, path.c_str(), "Re",1, &m_Re);

  io_create_dattribute(fileId, path.c_str(), "Pr",1);
  io_write_dattribute1(fileId, path.c_str(), "Pr",1, &m_Pr);

  io_create_dattribute(fileId, path.c_str(), "timeStep",1);
  io_write_dattribute1(fileId, path.c_str(), "timeStep",1, &m_timeStep);

  io_create_dattribute(fileId, path.c_str(), "time",1);
  io_write_dattribute1(fileId, path.c_str(), "time",1, &m_time);

  io_create_dattribute(fileId, path.c_str(), "physicalTimeStep",1);
  io_write_dattribute1(fileId, path.c_str(), "physicalTimeStep",1, &m_physicalTimeStep);

  io_create_dattribute(fileId, path.c_str(), "physicalTime",1);
  io_write_dattribute1(fileId, path.c_str(), "physicalTime",1, &m_physicalTime);

  io_create_iattribute(fileId, path.c_str(), "globalTimeStep",1);
  io_write_iattribute1(fileId, path.c_str(), "globalTimeStep",1, &globalTimeStep);

  io_create_dattribute(fileId, path.c_str(), "firstMaxResidual",1);
  io_write_dattribute1(fileId, path.c_str(), "firstMaxResidual",1, &m_firstMaxResidual);

  io_create_dattribute(fileId, path.c_str(), "firstAvrgResidual",1);
  io_write_dattribute1(fileId, path.c_str(), "firstAvrgResidual",1, &m_firstAvrgResidual);


  //save the time(Step) at which the grid motion started
  //does only work safely for constant time step
  if(m_movingGrid) {
    io_create_iattribute(fileId, "", "movingGridStepOffset", 1);
    io_write_iattribute1(fileId, "", "movingGridStepOffset", 1, &m_movingGridStepOffset);

    io_create_dattribute(fileId, "", "movingGridTimeOffset", 1);
    io_write_dattribute1(fileId, "", "movingGridTimeOffset", 1, &m_movingGridTimeOffset);

    //save whether or not the wave time step has been computed
    if(m_travelingWave) {
      io_create_iattribute(fileId, "", "waveTimeStepComputed", 1);
      io_write_iattribute1(fileId, "", "waveTimeStepComputed", 1, &m_waveTimeStepComputed);

      io_create_iattribute(fileId, "", "waveNoStepsPerCell", 1);
      io_write_iattribute1(fileId, "", "waveNoStepsPerCell", 1, &m_waveNoStepsPerCell);
    }
  }

  if(m_stgIsActive){
    io_create_iattribute(fileId, "", "stgNRAN",1);
    io_write_iattribute1(fileId, "", "stgNRAN",1, &m_stgMaxNoEddies);
  }

}

/**
 *
 * @author Pascal S. Meysonnat, May 9, 2011
 * modified 03.03.2016
 */

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::saveOutput(ZFSBool forceOutput)
{
  //Function to write the solution to file with iolibrary
  ZFSBool writeSolution = false;
  ZFSBool writePlane = false;
  ZFSBool writeBox = false;
  ZFSBool writeLine = false;
  ZFSBool writeAux = false;
  ZFSBool writeLiftDrag = false;

  //first find out which writeOut-mode we use (iteration or convective unit intervals)
  //then check which functions should write out in this timestep
  if(m_useConvectiveUnitWrite) {
    //in this mode we check the convective unit intervals
    //and write out files of each type if activated
    //activation is done by setting the interval
    //to a value greater than 0 (planeOutputInterval = 1)
    if(m_physicalTime -  (ZFSFloat)(m_noConvectiveOutputs)*m_convectiveUnitInterval >= m_convectiveUnitInterval) {

      //restart file output is still triggered by iteration counter
      writeSolution = isInInterval(m_solutionInterval);
      forceOutput = writeSolution;

      //activate the desired outputs
      writeSolution = m_sampleSolutionFiles;
      writePlane = (m_planeOutputInterval > 0);
      writeBox = (m_boxOutputInterval > 0);
      writeLine = (m_lineOutputInterval > 0);
      writeAux = (m_dragOutputInterval > 0);
      writeLiftDrag = (m_dragAsciiOutputInterval > 0);

      m_noConvectiveOutputs++;
      m_outputIterationNumber = m_noConvectiveOutputs;
    }
  } else {
    //in this mode we check the iteration intervals
    //for each writeOut-type (solution, plane, box, line, aux)
    writeSolution = isInInterval(m_solutionInterval);
    writePlane = isInInterval(m_planeOutputInterval);
    writeBox = isInInterval(m_boxOutputInterval);
    writeLine = isInInterval(m_lineOutputInterval);
    writeAux = isInInterval(m_dragOutputInterval);
    writeLiftDrag = isInInterval(m_dragAsciiOutputInterval);

    // cout<<"domainId():"<<domainId()<<" and m_blockId:"<<m_blockId<<" writeSolution:"<<writeSolution<<"and writePlane:"<<writePlane<<"writeBox:"<<writeBox<<"writeLine"<<writeLine<<"writeAux:"<<writeAux<<"writeLiftDrag:"<<writeLiftDrag<<endl;
    



    m_outputIterationNumber = globalTimeStep;
  }

  //compute vorticity if necessary
  if(m_vorticityOutput && (writeSolution || writePlane || writeBox || forceOutput)) {
    computeVorticity();
  }

  //compute velocity if wanted
  if(m_computeLambda2 && (writeSolution || writePlane || writeBox || forceOutput)) {
    computeLambda2Criterion();
  }

  //boxes, planes, auxdata and lines
  //are only available for 3D checked by function pointer
  if(writePlane && savePlanes) { savePlanes(); }
  if(writeBox && saveBoxes) {saveBoxes();}
  
  if(nDim>2){ //needs also to be implemented for 2d !!!!!!!!!!!!!!
    if(writeAux) {saveAuxData();}
    if(writeLiftDrag){saveLiftDragToAsciiFile();}
    if(writeLine) {saveOutputLines();}
  }

  if( writeSolution || forceOutput ){
    //save out the partitions also if desired, i.e, for debugging purposses
    if(m_savePartitionOutput){savePartitions();}
    //save postprocessing variables if activated
    saveAverageRestart();
    //save solution/restart file
    saveSolution(forceOutput);
    m_lastOutputTimeStep = globalTimeStep;
  }
}

template <ZFSInt nDim>
ZFSBool ZFSStrctrdBlck<nDim>::isInInterval(ZFSInt interval) {
  if(interval > 0) {
    if((globalTimeStep-m_outputOffset) % interval == 0 && globalTimeStep-m_outputOffset >= 0) {
      return true;
    }
  }
  return false;
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::writeRestartFile(ZFSBool forceOutput) {
  if(forceOutput && m_lastOutputTimeStep != globalTimeStep) {
    saveSolution(forceOutput);      
    saveAverageRestart();
  }
}


/**
 *
 * @author Pascal S. Meysonnat
 *
 * saves the soution to hdf5 file
 * modified 17.08.2018
 */
template <ZFSInt nDim>
template <ZFSBool primitiveOutput>
void ZFSStrctrdBlck<nDim>::saveOutputSolution(ZFSBool forceOutput) {
  ZFSChar gridFile[25];
  stringstream fileName;
  stringstream gridFileName;
  ZFSString tempG;

  fileName << m_solutionOutput << m_outputIterationNumber << m_outputFormat;

  if(m_movingGrid){
    if(!forceOutput && m_movingGridSaveGrid) {
      writeGrid();
    }
    gridFileName << "Grid" << globalTimeStep << m_outputFormat;
    tempG = gridFileName.str();
  }else{
    tempG = "../"+m_gridInputFileName;
  }
  strcpy(gridFile,tempG.c_str());


  zfs_log << "writing Solution file " << fileName.str() << " ... forceOutput: " << forceOutput << endl;


  ZFSInt file =io_openfile("hdf5", (fileName.str()).c_str(), "collective", m_zfsStrctrdComm);

  writeHeaderAttributes(file, gridFile, "solution");
  writePropertiesAsAttributes(file, "");

  ZFSInt allCells[3] = {0,0,0};
  ZFSInt stgNoEddieFields = 1200;
  ZFSString stgGlobalPathStr = "stgGlobal";
  const char* stgGlobalPath = stgGlobalPathStr.c_str();


  for (ZFSId i=0; i<m_noInputBlocks; i++){
    for(ZFSInt j =0; j<nDim; j++){
      allCells[j]=m_partition->inputBoxInfo[i]->DirLast[j];
    }
    //create datasets for the io library
    stringstream path;
    path << i;
    ZFSString blockPathStr = "block";
    blockPathStr += path.str();
    const char* blockPath = blockPathStr.c_str();

    ////////////////////////////////////////////////
    ///////// Create Primitive/Conservative Variables ////////
    ////////////////////////////////////////////////
    if(primitiveOutput){
      zfs_log << "writing primitive Output" << endl;
      for(ZFSInt v=0; v<m_maxNoVariables; v++) {    //junoh
        io_create_ddataset(file,blockPath, (m_pvariableNames[v]).c_str(), nDim, allCells);
      }
    }else{
      zfs_log << "writing conservative output" << endl;
      for(ZFSInt v=0; v<m_maxNoVariables; v++) {    //junoh
        io_create_ddataset(file,blockPath, (m_variableNames[v]).c_str(), nDim, allCells);
      }
    }

    ////////////////////////////////////////////////
    ///////// Create FQ Information ////////////////
    ////////////////////////////////////////////////
    if(FQ->noFQVariables>0){
      for(ZFSInt v=0; v< FQ->noFQVariables; ++v){
        if(FQ->fqWriteOutput[v]) {
          io_create_ddataset(file, blockPath , FQ->fqNames[v].c_str(), nDim,  allCells);
        }
      }
    }
    ///  Create BC2600 Information ////////////
    
    if(m_bc2600IsActive){
      ZFSInt allCells2600[3] = {allCells[0], allCells[1], allCells[2]};
      allCells2600[nDim-1] = m_noGhostLayers;
      stringstream path2600Str;
      path2600Str << blockPath << "/bc2600" << endl;
      const char* path2600 = (path2600Str.str()).c_str();

      for(ZFSId var=0; var<m_maxNoVariables; var++) {  //junoh

        io_create_ddataset(file, path2600, m_pvariableNames[var].c_str() , nDim, allCells2600);
      }
     
    }
   
    ////////////////////////////////////////////////
    ///////// Create BC2601 Information ////////////
    ////////////////////////////////////////////////
    if(m_bc2601IsActive){
      ZFSInt allCells2601[3] = {allCells[0], allCells[1], allCells[2]};
      allCells2601[nDim-2] = m_noGhostLayers;
      stringstream path2601Str;
      path2601Str << blockPath << "/bc2601" << endl;
      const char* path2601 = (path2601Str.str()).c_str();

      for(ZFSId var=0; var<m_maxNoVariables; var++) {
        io_create_ddataset(file, path2601, m_pvariableNames[var].c_str() , nDim, allCells2601);
      }
    }

    ////////////////////////////////////////////////
    ///////// Create STG Information ///////////////
    ////////////////////////////////////////////////
    if(m_stgIsActive){
      ZFSInt allCells7909[3];
      allCells7909[0] = allCells[0];
      allCells7909[1] = allCells[1];
      allCells7909[2] = 3;
      for(ZFSId var = 0; var < m_stgNoVariables; var++){
        stringstream stgPath;
        stgPath << blockPathStr << "/stg";
        stringstream fieldName;
        fieldName << "stgFQ" <<  var;
        io_create_ddataset(file, (stgPath.str()).c_str(), (fieldName.str()).c_str(), nDim, allCells7909);
      }
    }
  }
  ////////////////////////////////////////////////
  //////// Create Sandpaper Tripping Info ////////
  ////////////////////////////////////////////////
  if(m_useSandpaperTrip) {
    stringstream tripPath;
    tripPath << "/trip";
    ZFSInt dataSize = 2*m_tripNoModes;

    io_create_ddataset(file, (tripPath.str()).c_str(), "tripModesG", 1, &dataSize);
    io_create_ddataset(file, (tripPath.str()).c_str(), "tripModesH1", 1, &dataSize);
    io_create_ddataset(file, (tripPath.str()).c_str(), "tripModesH2", 1, &dataSize);
  }

  if(m_stgIsActive){
    stgNoEddieFields = int(m_stgNoEddieProperties*m_stgMaxNoEddies);
    io_create_ddataset(file, stgGlobalPath, "FQeddies", 1, &stgNoEddieFields);
  }

 //  MPI_Barrier(m_zfsStrctrdComm);
 // cout<<"check 16 save"<<endl;
  ////////////////////////////////////////////////
  ///////// Write Primitive/Conservative Variables /////////
  ////////////////////////////////////////////////
  stringstream path;
  path << m_inputBlockId;
  ZFSString blockPathStr = "block";
  blockPathStr += path.str();
  //write out the data
  ZFSInt ghostArray[3]={m_noGhostLayers,m_noGhostLayers,m_noGhostLayers};
  if(primitiveOutput){
    for(ZFSId v=0; v<m_maxNoVariables; ++v){  //junoh
      io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), (m_pvariableNames[v]).c_str(), nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&m_cells->pvariables[v][0]);
    }
  }else{
    for(ZFSId v=0; v<m_maxNoVariables; ++v){  //junoh
      io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), (m_variableNames[v]).c_str(), nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&m_cells->variables[v][0]);
    }
  }
  ////////////////////////////////////////////////
  ///////// Write FQ Data ////////////////////////
  ////////////////////////////////////////////////
  if(FQ->noFQVariables>0){
    for(ZFSInt v=0; v< FQ->noFQVariables; ++v){
      if(FQ->fqWriteOutput[v]) {
        io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), FQ->fqNames[v].c_str(), nDim,  m_nActiveCells, m_nOffsetCells,ghostArray ,&m_cells->fq[v][0]);
      }
    }
  }
 //    MPI_Barrier(m_zfsStrctrdComm);
  ////////////////////////////////////////////////
  ///////// Write BC2600 Data ////////////////////
  ////////////////////////////////////////////////
  if(m_bc2600IsActive) {
    stringstream path2600Str;
    path2600Str << blockPathStr << "/bc2600" << endl;
    const char* path2600 = (path2600Str.str()).c_str();

    ZFSInt noGhostLayers2600[3] = {m_noGhostLayers,m_noGhostLayers,m_noGhostLayers};
    noGhostLayers2600[nDim-1] = 0;
    if(m_bc2600) {
      for(ZFSId var=0; var<m_maxNoVariables; var++) {   //junoh
        io_write_ddataset_part_ghost_array(file, path2600, m_pvariableNames[var].c_str(),  nDim,  m_bc2600noActiveCells,  m_bc2600noOffsetCells, noGhostLayers2600, &m_bc2600Variables[var][0]); //m_bc2600Variables
      }
    }else{
      for(ZFSId var=0; var<m_maxNoVariables; var++) {   //junoh
        io_write_ddataset_part_ghost_array(file, path2600, m_pvariableNames[var].c_str(),  nDim,  m_bc2600noActiveCells,  m_bc2600noOffsetCells, noGhostLayers2600, NULL);
      }
    }
  }
  ////////////////////////////////////////////////
  ///////// Write BC2601 Data ////////////////////
  ////////////////////////////////////////////////
  if(m_bc2601IsActive) {
    stringstream path2601Str;
    path2601Str << blockPathStr << "/bc2601" << endl;
    const char* path2601 = (path2601Str.str()).c_str();

    ZFSInt noGhostLayers2601[3] = {m_noGhostLayers,m_noGhostLayers,m_noGhostLayers};
    noGhostLayers2601[nDim-2] = 0;
    if(m_bc2601) {
      for(ZFSId var=0; var<m_maxNoVariables; var++) {  //junoh
        io_write_ddataset_part_ghost_array(file, path2601, m_pvariableNames[var].c_str(),  nDim,  m_bc2601noActiveCells,  m_bc2601noOffsetCells, noGhostLayers2601, &m_bc2601Variables[var][0]);
      }
    }else{
      for(ZFSId var=0; var<m_maxNoVariables; var++) { //junoh
        io_write_ddataset_part_ghost_array(file, path2601, m_pvariableNames[var].c_str(),  nDim,  m_bc2601noActiveCells,  m_bc2601noOffsetCells, noGhostLayers2601, NULL);
      }
    }
  }

  ////////////////////////////////////////////////
  ///////// Write STG Information ////////////////
  ////////////////////////////////////////////////
  if(m_stgIsActive){
    ZFSInt VBOffset = 0;
    if(m_stgRootRank){
      io_write_ddataset_part(file, stgGlobalPath, "FQeddies", 1,  &stgNoEddieFields,  &VBOffset, &m_stgEddies[0][0]);
    }else{
      stgNoEddieFields = 0;
      io_write_ddataset_part(file, stgGlobalPath, "FQeddies", 1,  &stgNoEddieFields,  &VBOffset, NULL);
    }

    // MPI_Barrier(m_zfsStrctrdComm);
    // cout<<"check 22 save"<<endl;

    //if this domain has part of the STG bc write value, otherwise only write NULL
    if(m_stgLocal){
      ZFSId noActiveStgCells = (m_stgBoxSize[0]-2*m_noGhostLayers)*(m_stgBoxSize[1]-2*m_noGhostLayers)*3;
      ZFSFloatScratchSpace stgFqDummy(m_stgNoVariables, noActiveStgCells, __CALLING_FUNCTION__, "stgFqDummy");

      for(ZFSId var=0; var < m_stgNoVariables; var++) {
        for(ZFSId k=m_noGhostLayers; k<m_stgBoxSize[0]-m_noGhostLayers; k++){
          for(ZFSId j=m_noGhostLayers; j<m_stgBoxSize[1]-m_noGhostLayers; j++){
            for(ZFSId i=0; i<m_stgBoxSize[2]; i++){

              ZFSId cellIdBC = i+(j + k*m_stgBoxSize[1])*3;
              ZFSId cellIdDummy = i+((j-m_noGhostLayers) + (k-m_noGhostLayers)*(m_stgBoxSize[1]-2*m_noGhostLayers))*3;
              stgFqDummy(var, cellIdDummy) = m_cells->stg_fq[var][cellIdBC];
            }
          }
        }
      }

      ZFSInt bcOffset[3] = {m_nOffsetCells[0], m_nOffsetCells[1], 0};
      ZFSInt bcCells[3] = {m_stgBoxSize[0]-2*m_noGhostLayers, m_stgBoxSize[1]-2*m_noGhostLayers, m_stgBoxSize[2]};

      for(ZFSId var = 0; var < m_stgNoVariables; var++){
        stringstream fieldName;
        stringstream stgPath;
        stgPath << blockPathStr << "/stg";
        fieldName << "stgFQ" <<  var;
        io_write_ddataset_part(file, (stgPath.str()).c_str(), (fieldName.str()).c_str(),  nDim,  bcCells,  bcOffset, &stgFqDummy(var,0));
      }
    }else{
      ZFSInt bcOffset[3] = {0, 0, 0};
      ZFSInt bcCells[3] = {0, 0, 0};

      for(ZFSId var = 0; var < m_stgNoVariables; var++){
        stringstream fieldName;
        stringstream stgPath;
        stgPath << blockPathStr << "/stg";
        fieldName << "stgFQ" <<  var;
        io_write_ddataset_part(file, (stgPath.str()).c_str(), (fieldName.str()).c_str(),  nDim,  bcCells,  bcOffset, NULL);
      }
    }

  }

  ////////////////////////////////////////////////
  ///////// Sandpaper Tripping ///////////////////
  ////////////////////////////////////////////////
  
  if(m_useSandpaperTrip) {
    stringstream tripPath;
    tripPath << "/trip";

    if(domainId()==0) {
      ZFSInt offset = 0;
      ZFSInt dataSize = m_tripNoModes*2;

      io_write_ddataset_part(file, (tripPath.str()).c_str(), "tripModesG",  1,  &dataSize,  &offset, m_tripModesG);
      io_write_ddataset_part(file, (tripPath.str()).c_str(), "tripModesH1",  1,  &dataSize,  &offset, m_tripModesH1);
      io_write_ddataset_part(file, (tripPath.str()).c_str(), "tripModesH2",  1,  &dataSize,  &offset, m_tripModesH2);
    } else {
      ZFSInt offset = 0;
      ZFSInt dataSize = 0;
      io_write_ddataset_part(file, (tripPath.str()).c_str(), "tripModesG",  1,  &dataSize,  &offset, NULL);
      io_write_ddataset_part(file, (tripPath.str()).c_str(), "tripModesH1",  1,  &dataSize,  &offset, NULL);
      io_write_ddataset_part(file, (tripPath.str()).c_str(), "tripModesH2",  1,  &dataSize,  &offset, NULL);
    }
  }

  io_closefile(file);
  zfs_log << "...-> OK " << endl;
}
template void ZFSStrctrdBlck<2>::saveOutputSolution<true>(ZFSBool);
template void ZFSStrctrdBlck<2>::saveOutputSolution<false>(ZFSBool);
template void ZFSStrctrdBlck<3>::saveOutputSolution<true>(ZFSBool);
template void ZFSStrctrdBlck<3>::saveOutputSolution<false>(ZFSBool);


/**
 *
 * @author: Pascal Meysonnat
 * @date: 01.01.1010 
 * 
 * saves the partition into file, good for debugging
 *
 */
template <ZFSInt nDim>
template <ZFSBool primitiveOutput>
void ZFSStrctrdBlck<nDim>::saveOutputPartitions()
{
  //Function to write the solution to file with iolibrary
  stringstream fileName;
  ZFSChar gridFile[25];
  stringstream gridFileName;
  ZFSString tempG;

  ZFSInt noCells[3] = {0,0,0};
	// cout<<"FQ->noFQVariables:"<<FQ->noFQVariables<<" domainId:"<<domainId()<<endl;
	
#ifndef ZFS_EXTRA_DEBUG
    fileName <<  m_solutionOutput << "partitioned" << globalTimeStep << m_outputFormat;
#endif

#ifdef ZFS_EXTRA_DEBUG
  fileName << m_solutionOutput << globalTimeStep<< "par_RK" << m_RKStep << m_outputFormat ;
#endif
  writeGridPointsWithGhostPoints();

  if(m_movingGrid)
  {
    gridFileName << "Grid" << globalTimeStep << ".hdf5";
    tempG = gridFileName.str();
  }
  else
  {
    tempG = "../" + m_gridInputFileName;
  }

  strcpy(gridFile,tempG.c_str());

  ZFSInt file =io_openfile("hdf5", (fileName.str()).c_str(), "collective", m_zfsStrctrdComm);
  writeHeaderAttributes(file, gridFile, "solution");
  writePropertiesAsAttributes(file, "");

  //save with ghostcells ==> multiple blocks;
  for(ZFSInt i=0; i<noDomains(); i++) {
    //create datasets for the io library
    for(ZFSId j=0; j<nDim; j++) {
      noCells[j]=m_partition->outputBoxInfo[getBoxId(i)]->DirLast[j]-1+2*m_noGhostLayers;
    }
    stringstream path;
    path << getBoxId(i);
    ZFSString blockPathStr = "block";
    blockPathStr += path.str();
    const char* blockPath = blockPathStr.c_str();
    //create dataset for primitive/conservative variables
    if(primitiveOutput){
      for(ZFSId v=0; v<m_maxNoVariables; v++){ //junoh
        io_create_ddataset(file, blockPath , m_pvariableNames[v].c_str(), nDim,  noCells);
      }
    } else {
      for(ZFSId v=0; v<CV->noVariables; v++){
        io_create_ddataset(file, blockPath , m_variableNames[v].c_str(), nDim,  noCells);
      }
    }
    //create dataset for fq field variables
    for(ZFSId v=0; v< FQ->noFQVariables; v++){
      io_create_ddataset(file, blockPath , FQ->fqNames[v].c_str(), nDim, noCells);
    }
#ifdef ZFS_EXTRA_DEBUG
    io_create_ddataset(file, blockPath , "FluxRhoI_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoJ_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoK_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoUI_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoUJ_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoUK_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoVI_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoVJ_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoVK_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoWI_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoWJ_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoWK_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoEI_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoEJ_conv", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoEK_conv", nDim,  noCells);

    io_create_ddataset(file, blockPath , "FluxRhoI_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoJ_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoK_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoUI_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoUJ_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoUK_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoVI_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoVJ_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoVK_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoWI_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoWJ_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoWK_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoEI_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoEJ_visc", nDim,  noCells);
    io_create_ddataset(file, blockPath , "FluxRhoEK_visc", nDim,  noCells);
#endif

    if(m_stgIsActive){
      ZFSInt allCells7909[3];
      allCells7909[0] = m_partition->outputBoxInfo[getBoxId(i)]->DirLast[0]-1+2*m_noGhostLayers;
      allCells7909[1] = m_partition->outputBoxInfo[getBoxId(i)]->DirLast[1]-1+2*m_noGhostLayers;
      allCells7909[2] = 3;

      for(ZFSId var = 0; var < m_stgNoVariables; var++){
        stringstream fieldName;
        stringstream stgPath;
        stgPath << blockPathStr << "/stg";
        fieldName << "stgFQ" <<  var;
        io_create_ddataset(file, (stgPath.str()).c_str(), (fieldName.str()).c_str(), nDim, allCells7909);
      }
    }
  }
  //write the values into the array so that we can visualize it
  ZFSInt offset[3]={0,0,0};
  for(ZFSId j=0; j<nDim; j++)
  {
    noCells[j]=m_partition->outputBoxInfo[getBoxId(domainId())]->DirLast[j]-1+2*m_noGhostLayers;
  }
  stringstream path;
  path << getBoxId(domainId());
  ZFSString blockPathStr = "block";
  blockPathStr += path.str();

  //write primitive/conservative variables
  if(primitiveOutput){
    for(ZFSId v=0; v<m_maxNoVariables; ++v){ //junoh
      io_write_ddataset_part(file,blockPathStr.c_str(), m_pvariableNames[v].c_str() , nDim, noCells, offset, &m_cells->pvariables[v][0]);
    }
  }else{
    for(ZFSId v=0; v<CV->noVariables; ++v){
      io_write_ddataset_part(file,blockPathStr.c_str(), m_variableNames[v].c_str() , nDim, noCells, offset, &m_cells->variables[v][0]);
    }
  }

  for(ZFSId v=0; v<FQ->noFQVariables; v++){
    io_write_ddataset_part(file,blockPathStr.c_str(), FQ->fqNames[v].c_str(), nDim, noCells, offset, &m_cells->fq[v][0]);
  }

#ifdef ZFS_EXTRA_DEBUG
  //fluxes:
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoI_conv" , nDim, noCells, offset, &convFluxOut[0][CV->RHO*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoJ_conv" , nDim, noCells, offset, &convFluxOut[1][CV->RHO*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoK_conv" , nDim, noCells, offset, &convFluxOut[2][CV->RHO*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoUI_conv" , nDim, noCells, offset, &convFluxOut[0][CV->RHO_U*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoUJ_conv" , nDim, noCells, offset, &convFluxOut[1][CV->RHO_U*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoUK_conv" , nDim, noCells, offset, &convFluxOut[2][CV->RHO_U*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoVI_conv" , nDim, noCells, offset, &convFluxOut[0][CV->RHO_V*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoVJ_conv" , nDim, noCells, offset, &convFluxOut[1][CV->RHO_V*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoVK_conv" , nDim, noCells, offset, &convFluxOut[2][CV->RHO_V*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoWI_conv" , nDim, noCells, offset, &convFluxOut[0][CV->RHO_W*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoWJ_conv" , nDim, noCells, offset, &convFluxOut[1][CV->RHO_W*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoWK_conv" , nDim, noCells, offset, &convFluxOut[2][CV->RHO_W*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoEI_conv" , nDim, noCells, offset, &convFluxOut[0][CV->RHO_E*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoEJ_conv" , nDim, noCells, offset, &convFluxOut[1][CV->RHO_E*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoEK_conv" , nDim, noCells, offset, &convFluxOut[2][CV->RHO_E*m_noStrctrdCells]);

  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoI_visc" , nDim, noCells, offset, &viscFluxOut[0][CV->RHO*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoJ_visc" , nDim, noCells, offset, &viscFluxOut[1][CV->RHO*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoK_visc" , nDim, noCells, offset, &viscFluxOut[2][CV->RHO*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoUI_visc" , nDim, noCells, offset, &viscFluxOut[0][CV->RHO_U*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoUJ_visc" , nDim, noCells, offset, &viscFluxOut[1][CV->RHO_U*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoUK_visc" , nDim, noCells, offset, &viscFluxOut[2][CV->RHO_U*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoVI_visc" , nDim, noCells, offset, &viscFluxOut[0][CV->RHO_V*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoVJ_visc" , nDim, noCells, offset, &viscFluxOut[1][CV->RHO_V*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoVK_visc" , nDim, noCells, offset, &viscFluxOut[2][CV->RHO_V*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoWI_visc" , nDim, noCells, offset, &viscFluxOut[0][CV->RHO_W*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoWJ_visc" , nDim, noCells, offset, &viscFluxOut[1][CV->RHO_W*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoWK_visc" , nDim, noCells, offset, &viscFluxOut[2][CV->RHO_W*m_noStrctrdCells]);
  io_write_ddataset_ppppart(file,blockPathStr.c_str(), "FluxRhoEI_visc" , nDim, noCells, offset, &viscFluxOut[0][CV->RHO_E*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoEJ_visc" , nDim, noCells, offset, &viscFluxOut[1][CV->RHO_E*m_noStrctrdCells]);
  io_write_ddataset_part(file,blockPathStr.c_str(), "FluxRhoEK_visc" , nDim, noCells, offset, &viscFluxOut[2][CV->RHO_E*m_noStrctrdCells]);
#endif

  ////////////////////////////////////////////////
  ///////// Write STG Information ////////////////
  ////////////////////////////////////////////////
  if(m_stgIsActive){
    ZFSInt bcOffset[3] ={0,0,0};
    ZFSInt bcCells[3] = {m_nCells[0], m_nCells[1], 3};

    for(ZFSId var = 0; var < m_stgNoVariables; var++){
      stringstream fieldName;
      stringstream stgPath;
      stgPath << blockPathStr << "/stg";
      fieldName << "stgFQ" <<  var;
      io_write_ddataset_part(file, (stgPath.str()).c_str(), (fieldName.str()).c_str(),  nDim,  bcCells,  bcOffset, &m_cells->stg_fq[var][0]);
    }
  }

  io_closefile(file);
}
template void ZFSStrctrdBlck<2>::saveOutputPartitions<true>();
template void ZFSStrctrdBlck<2>::saveOutputPartitions<false>();
template void ZFSStrctrdBlck<3>::saveOutputPartitions<true>();
template void ZFSStrctrdBlck<3>::saveOutputPartitions<false>();


/** \brief Load Restart File (primitive and conservative output)
 *         general formulation 
 *  
 * \author Pascal Meysonnat
 * \date: 01.01.1010
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::loadRestartFile(){
  TRACE();
  RECORD_TIMER_START(m_tloadRestart);
  zfs_log << "loading Restart file ... " << endl;
  stringstream restartFileName;
  ZFSInt restartFileId =-1;
  
  if(!m_useNonSpecifiedRestartFile) {
    ZFSString restartFile = "restart.hdf5";
    //! \page propertyPage1
    //\section restartVariableFileName
    //<code>ZFSInt ZFSStrctrdBlck::loadRestartFile </code>\n
    //default = <code> 0 </code>\n \n
    //Name of the specific restart file.\n
    //Possible values are:\n
    //<ul>
    //<li>string</li>
    //</ul>
    //Keywords: <i>RESTART, STRCTRD</i>
    //
    restartFile = *(ZFSContext::getProperty("restartVariablesFileName", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL )->asString(0)); //this should be removed 
    restartFileName << outputDir() << restartFile;
  }else{
    ZFSString restartFile = "restart";
    //! \page propertyPage1
    //\section restartTimeStep
    //<code>ZFSInt ZFSStrctrdBlck::loadRestartFile </code>\n
    //default = <code> 0 </code>\n \n
    //Start Iteration of the specific restart file.\n
    //Possible values are:\n
    //<ul>
    //<li>integer</li>
    //</ul>
    //Keywords: <i>RESTART, STRCTRD</i>
    //
    ZFSInt restartTimeStep = *(ZFSContext::getProperty("restartTimeStep", m_blockId, __CALLING_FUNCTION__, (ZFSInt*) NULL )->asInt(0));
    restartFileName << outputDir() << restartFile << restartTimeStep << ".hdf5";
  }

  restartFileId = io_openfile("hdf5" ,(restartFileName.str()).c_str(),"collective", m_zfsStrctrdComm); //not collective on juqueen for the moement
  
  //check if restart and grid do fit together through UID
  if(!m_ignoreUID) {
    const char* aUID=new char[18];
    io_read_sattribute1(restartFileId,"", "UID", aUID);
    if(strcmp(aUID, m_uID.c_str())!=0) {
      zfsTerm(1, __CALLING_FUNCTION__, "FATAL: the files do not match each other according to the attribute UID");
    }
    delete[] aUID;
  }
  //check general attributes
  io_read_iattribute1(restartFileId,"", "globalTimeStep", &globalTimeStep);
  io_read_dattribute1(restartFileId,"", "time", &m_time);
  io_read_dattribute1(restartFileId,"", "physicalTime", &m_physicalTime);
  io_read_dattribute1(restartFileId,"", "physicalTimeStep", &m_physicalTimeStep);
  io_read_dattribute1(restartFileId,"", "firstMaxResidual",&m_firstMaxResidual);
  io_read_dattribute1(restartFileId,"", "firstAvrgResidual",&m_firstAvrgResidual);
  io_read_dattribute1(restartFileId,"", "timeStep",&m_timeStep);
  //check if primitive or conservative output is given
  ZFSInt isPrimitiveOutput=0;
  if(io_checkObj(restartFileId, "", "primitiveOutput")){
    io_read_iattribute1(restartFileId,"","primitiveOutput",&isPrimitiveOutput);
  }
  //if moving Grid is actived read the moving grid time offset (if it exists)
  //otherwise assume this is an initial start and set the current restart time
  //as the moving grid time offset
  if(m_movingGrid) {
    if(io_checkObj(restartFileId, "", "movingGridStepOffset")) {
      io_read_dattribute1(restartFileId, "", "movingGridTimeOffset", &m_movingGridTimeOffset);
      io_read_iattribute1(restartFileId, "", "movingGridStepOffset", &m_movingGridStepOffset);
      m_movingGridInitialStart = false;
    } else {
      m_movingGridTimeOffset = m_time;
      m_movingGridStepOffset = m_restartTimeStep;
      m_movingGridInitialStart = true;
    }
    //check if the wave time step has already been computed
    if(io_checkObj(restartFileId, "", "waveTimeStepComputed")) {
      io_read_iattribute1(restartFileId, "", "waveTimeStepComputed", &m_waveTimeStepComputed);
    }else{ m_waveTimeStepComputed = false;}

    if(io_checkObj(restartFileId, "", "waveNoStepsPerCell")) {
      io_read_iattribute1(restartFileId, "", "waveNoStepsPerCell", &m_waveNoStepsPerCell);
    }else{m_waveNoStepsPerCell = 1;}
  }
  //check for convective unit output
  if(m_useConvectiveUnitWrite) {
    m_noConvectiveOutputs = (ZFSInt)(m_physicalTime / m_convectiveUnitInterval);
    zfs_log << "Convective unit output iteration counter: " << m_noConvectiveOutputs << endl;
  }
  //check for moving grid initial start
  if(m_movingGridInitialStart) {
    m_movingGridTimeOffset = m_time;
  }
  if(domainId() == 0) { cout << "Restarting at GlobalTimeStep " << globalTimeStep << endl;}
  m_restartTimeStep=globalTimeStep;
  //now read in the data!
  zfs_log << "-> reading in the data ... " << endl;
  ZFSInt start[nDim];
  ZFSInt Id = m_partition->outputBoxInfo[domainId()]->cpu;
  //the offset in the file
  for (ZFSInt i=0; i<nDim;i++) {
    //get offset in file
    start[i]=m_partition->outputBoxInfo[Id]->offset[i];
  }
  zfs_log << "Loading restart variables..." << endl;
  if(domainId()==0) {cout << "Loading restart variables..." << endl;}
  stringstream blockNumber;
  blockNumber << m_inputBlockId;
  ZFSString blockPathStr = "/block";
  blockPathStr += blockNumber.str();
  const char* blockPath = blockPathStr.c_str();
  //record tim to load restart file members
  RECORD_TIMER_START(m_tloadRestartVars);
  //check for primitive input or conservative input
  if(isPrimitiveOutput){  //junoh
    for(ZFSId var=0; var<m_maxNoVariables; var++) { io_read_ddataset_part1d1(restartFileId, blockPath, m_pvariableNames[var].c_str(), nDim, start, m_nActiveCells, m_cells->pvariables[var] ); zfs_log << "Reading " << m_pvariableNames[var] << endl;}
    //compute the conservative variables so that we can proceed as usual
    //computeConservativeVariables();
  }else{
    for(ZFSId var=0; var<CV->noVariables; var++) {io_read_ddataset_part1d1(restartFileId, blockPath, m_variableNames[var].c_str(), nDim, start, m_nActiveCells, m_cells->variables[var] );}
    //computePrimitiveVariables();
  }
  if(m_zonal){   //junoh 
    io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->AVG_U].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->AVG_U]);
    io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->AVG_V].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->AVG_V]);
    io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->AVG_W].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->AVG_W]);
    io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->AVG_RHO].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->AVG_RHO]);
    io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->AVG_P].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->AVG_P]);
    io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->NU_T].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->NU_T]);
    
    
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->FLUC_U].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->FLUC_U]);
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->FLUC_V].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->FLUC_V]);
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->FLUC_W].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->FLUC_W]);
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->FLUC_W].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->FLUC_W]);
    
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->FLUC_UU].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->FLUC_UU]);
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->FLUC_VV].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->FLUC_VV]);
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->FLUC_WW].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->FLUC_WW]);
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->FLUC_UV].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->FLUC_UV]);
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->FLUC_UW].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->FLUC_UW]);
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->FLUC_VW].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->FLUC_VW]);
    
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->RECONST_NUT].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->RECONST_NUT]);
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->RECONST_NUTILDE].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->RECONST_NUTILDE]);
    // io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->NUTILDE].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->NUTILDE]);


    
     FQ->loadedFromRestartFile[FQ->AVG_U] = true;
     FQ->loadedFromRestartFile[FQ->AVG_V] = true;
     FQ->loadedFromRestartFile[FQ->AVG_W] = true;
     FQ->loadedFromRestartFile[FQ->AVG_RHO] = true;
     FQ->loadedFromRestartFile[FQ->AVG_P] = true;
     FQ->loadedFromRestartFile[FQ->NU_T] = true;
     
     // FQ->loadedFromRestartFile[FQ->FLUC_U] = true;
     // FQ->loadedFromRestartFile[FQ->FLUC_V] = true;
     // FQ->loadedFromRestartFile[FQ->FLUC_W] = true;
     // FQ->loadedFromRestartFile[FQ->FLUC_UU] = true;
     // FQ->loadedFromRestartFile[FQ->FLUC_VV] = true;
     // FQ->loadedFromRestartFile[FQ->FLUC_WW] = true;
     // FQ->loadedFromRestartFile[FQ->FLUC_UV] = true;
     // FQ->loadedFromRestartFile[FQ->FLUC_UW] = true;
     // FQ->loadedFromRestartFile[FQ->FLUC_VW] = true;
    
     // FQ->loadedFromRestartFile[FQ->RECONST_NUT] = true;
     // FQ->loadedFromRestartFile[FQ->RECONST_NUTILDE] = true;
     // FQ->loadedFromRestartFile[FQ->NUTILDE] = true;

  }

  RECORD_TIMER_STOP(m_tloadRestartVars);
  zfs_log << "Loading restart variables... SUCCESSFUL!" << endl;
  if(domainId()==0) {cout << "Loading restart variables... SUCCESSFUL!" << endl;}
  zfs_log << "-> reading in auxilliary data for restart ..." << endl;
  zfs_log << "--> ... sponge ..." << endl;
  ////////////////////////////////////////////////////////////
  ///////////////// SPONGE ///////////////////////////////////
  ////////////////////////////////////////////////////////////

  RECORD_TIMER_START(m_tloadRestartSponge);
  if(m_useSponge) {
    zfs_log << "--> ... sponge ..." << endl;
    if(domainId()==0) {cout << "Loading sponge data..." << endl;}
    if(m_computeSpongeFactor==false) {
      io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->SPONGE_FACTOR].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->SPONGE_FACTOR] );
      FQ->loadedFromRestartFile[FQ->SPONGE_FACTOR] = true;
    }
    if(m_spongeLayerType == 2) {
      io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->SPONGE_RHO].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->SPONGE_RHO] );
      FQ->loadedFromRestartFile[FQ->SPONGE_RHO] = true;
      io_read_ddataset_part1d1(restartFileId,blockPath, FQ->fqNames[FQ->SPONGE_RHO_E].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->SPONGE_RHO_E]);
      FQ->loadedFromRestartFile[FQ->SPONGE_RHO_E] = true;
    }
    if(m_spongeLayerType == 4) {
      io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->SPONGE_RHO].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->SPONGE_RHO] );
      FQ->loadedFromRestartFile[FQ->SPONGE_RHO] = true;
    }
    if(domainId()==0) {cout << "Loading sponge data... SUCCESSFUL!" << endl;}
    zfs_log << "--> ... sponge ... SUCCESSFUL" << endl;
  }
  RECORD_TIMER_STOP(m_tloadRestartSponge);
   
  ////////////////////////////////////////////////////////////
  ///////////////// CELL SHIFTING ////////////////////////////
  ////////////////////////////////////////////////////////////
  RECORD_TIMER_START(m_tloadRestartShift);
  if(isPrimitiveOutput){
    this->shiftCellValuesRestart<true>();
    computeConservativeVariables();
  }else{
    this->shiftCellValuesRestart<false>();
    computePrimitiveVariables();
  }
  if(m_zonal){
    averagedFillGhostCells();  //junoh
  }

  RECORD_TIMER_STOP(m_tloadRestartShift);

  //load special variables
  loadRestartBC2600();
  if(nDim==3){//only implemented for 3d or does not work in 2d.
    loadRestartBC2601();
    if(m_stgIsActive) loadRestartSTG(isPrimitiveOutput);
  }

  ////////////////////////////////////////////////////////////
  ///////////////// CHANGE MACH NUMBER ///////////////////////
  ////////////////////////////////////////////////////////////
  if(m_changeMa){
    ZFSFloat oldMa=F0;
    io_read_dattribute1(restartFileId,"", "Ma", &oldMa);
    convertRestartVariables(oldMa);
    //if we use the STG also convert the vars in the STG fields
    if(m_stgIsActive) {
      convertRestartVariablesSTG(oldMa);
    }
  }
 
  zfs_log << "-> reading in auxilliary data for restart ...SUCCESSFUL" << endl;
  io_closefile(restartFileId);
  zfs_log << "loading Restart file ... SUCCESSFUL " << endl;
  //writeGridPointsWithGhostPoints();
  //saveOutput(true);
  //saveOutputPartitions<nDim>();
  //zfsTerm(-1, __CALLING_FUNCTION__, "Pascal killed it");
  RECORD_TIMER_STOP(m_tloadRestart);
}      

/**
 *     function to shift the values in the cell after restart
 *     to correct position (reading in does not take into account
 *     the ghost cells)
 *     /author Pascal Meysonnat
 *     /date   01.01.1010
 *
 */
template <ZFSInt nDim>
template <ZFSBool isPrimitive>
void ZFSStrctrdBlck<nDim>::shiftCellValuesRestart()
{
  TRACE(); 
  if(nDim==3){    
    //accounting for the ghost layers and shift the values to the right place
    for(ZFSId k=(m_nActiveCells[0]-1); k>=0; k-- ){
      for(ZFSId j=(m_nActiveCells[1]-1); j>=0; j-- ){
        for(ZFSId i=(m_nActiveCells[2]-1); i>=0; i--){
          const ZFSId cellId_org=i+(j+k*m_nActiveCells[1])*m_nActiveCells[2];
          const ZFSId i_new=i+m_noGhostLayers;
          const ZFSId j_new=j+m_noGhostLayers;
          const ZFSId k_new=k+m_noGhostLayers;
          const ZFSId cellId=i_new+(j_new+k_new*m_nCells[1])*m_nCells[2];
          if(!isPrimitive){
            for(ZFSId var=0; var<CV->noVariables; var++) {
              m_cells->variables[var][cellId]=m_cells->variables[var][cellId_org];
              m_cells->variables[var][cellId_org]=F0;
            }
          }else{
            for(ZFSId var=0; var<m_maxNoVariables; var++) { //junoh
              m_cells->pvariables[var][cellId]=m_cells->pvariables[var][cellId_org];
              m_cells->pvariables[var][cellId_org]=F0;
            }
          }
          //also shift values in the FQ field
	  // cout<<"FQ->noFQVariables:"<<FQ->noFQVariables<<endl;
          if(FQ->noFQVariables > 0) {
            for(ZFSId var = 0; var < FQ->noFQVariables; var++) {
              if(FQ->loadedFromRestartFile[var]) {
                m_cells->fq[var][cellId]=m_cells->fq[var][cellId_org];
                m_cells->fq[var][cellId_org]=F0;
              }
            }
          }
        }
      }
    }
  }  
  if(nDim==2){
    for(ZFSId j=(m_nActiveCells[0]-1); j>=0; j-- ){
      for(ZFSId i=(m_nActiveCells[1]-1); i>=0; i--){
        const ZFSId cellId_org=i+j*m_nActiveCells[1];
        const ZFSId i_new=i+m_noGhostLayers;
        const ZFSId j_new=j+m_noGhostLayers;
        const ZFSId cellId=i_new+j_new*m_nCells[1];
        if(!isPrimitive){
          for(ZFSId var=0; var<CV->noVariables; var++) {
            m_cells->variables[var][cellId]=m_cells->variables[var][cellId_org];
            m_cells->variables[var][cellId_org]=F0;
          }
        }else{
          for(ZFSId var=0; var<m_maxNoVariables; var++) { //junoh
            m_cells->pvariables[var][cellId]=m_cells->pvariables[var][cellId_org];
            m_cells->pvariables[var][cellId_org]=F0;
          }
        }
        //also shift values in the FQ field
        if(FQ->noFQVariables > 0) {
          for(ZFSId var = 0; var < FQ->noFQVariables; var++) {
            if(FQ->loadedFromRestartFile[var]) {
              m_cells->fq[var][cellId]=m_cells->fq[var][cellId_org];
              m_cells->fq[var][cellId_org]=F0;
            }
          }
        }
      }
    }
  }
}

template void ZFSStrctrdBlck<2>::shiftCellValuesRestart<true>();
template void ZFSStrctrdBlck<2>::shiftCellValuesRestart<false>();
template void ZFSStrctrdBlck<3>::shiftCellValuesRestart<true>();
template void ZFSStrctrdBlck<3>::shiftCellValuesRestart<false>();



template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::saveAuxData(){
  computeAuxData();
  stringstream fileName;
  ZFSChar gridFile[25];
  stringstream gridFileName;
  ZFSString tempG;

  fileName << m_auxOutputDir <<  "auxData" << m_outputIterationNumber << m_outputFormat;

  if(m_movingGrid) {
    gridFileName << "Grid" << globalTimeStep << ".hdf5";
    tempG = gridFileName.str();
  } else {
    tempG = "../" + m_gridInputFileName;
  }

  strcpy(gridFile,tempG.c_str());
  ZFSInt file =io_openfile("hdf5", (fileName.str()).c_str(), "collective", m_zfsStrctrdComm);
  writeHeaderAttributes(file, gridFile, "auxdata");
  writePropertiesAsAttributes(file,"" );
  const ZFSString powerNamesVisc[3]={"Pxv","Pyv","Pzv"};
  const ZFSString powerNamesPres[3]={"Pxp","Pyp","Pzp"};
  
  const ZFSString dataNames[9] = {"cfx","cfy","cfz","ax","ay","az","x","y","z"};
  ZFSId noFields = 3;
  if(m_detailAuxData) {
    noFields = 9;
  }

  for(ZFSUint i=0; i<m_windowInfo->globalStrctrdBndryCndMaps.size(); ++i) {
    const ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_windowInfo->globalStrctrdBndryCndMaps[i]->BC)/1000.0);
    if(firstDigit==1){
      ZFSInt datasetSize[2];
      ZFSId dim1=0;
      for(int dim=0; dim<nDim; ++dim){
        if(m_windowInfo->globalStrctrdBndryCndMaps[i]->end2[dim]==m_windowInfo->globalStrctrdBndryCndMaps[i]->start2[dim]) {
          continue;
        }
        datasetSize[dim1]=m_windowInfo->globalStrctrdBndryCndMaps[i]->end2[dim]-m_windowInfo->globalStrctrdBndryCndMaps[i]->start2[dim];
        dim1++;
      }
      int datasetSize1[2]={datasetSize[1], datasetSize[0]};

      if(m_bCfCpCoeff){
        stringstream datasetId;
        datasetId<< m_windowInfo->globalStrctrdBndryCndMaps[i]->Id2;
        ZFSString pathName="window";
        pathName+=datasetId.str();
        ZFSString dataNamesCp="cp";
        for(ZFSInt j=0; j<noFields; j++) {io_create_ddataset(file, pathName.c_str() , dataNames[j].c_str() , nDim-1, datasetSize1 );}
        io_create_ddataset(file, pathName.c_str() , dataNamesCp.c_str() , nDim-1, datasetSize1);
      }

      if(m_bPower){
        stringstream datasetId;
        datasetId<< m_windowInfo->globalStrctrdBndryCndMaps[i]->Id2;
        ZFSString pathName="window";
        pathName+=datasetId.str();
        for(ZFSInt j=0; j<nDim; j++) {
          io_create_ddataset(file, pathName.c_str() , powerNamesVisc[j].c_str() , nDim-1, datasetSize1 );
          io_create_ddataset(file, pathName.c_str() , powerNamesPres[j].c_str() , nDim-1, datasetSize1 );
        }
      }
    }
  }

  for(ZFSUint i=0; i<m_windowInfo->globalStrctrdBndryCndMaps.size(); ++i) {
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_windowInfo->globalStrctrdBndryCndMaps[i]->BC)/1000.0);
    if(firstDigit==1) {
      stringstream datasetname;
      datasetname<< m_windowInfo->globalStrctrdBndryCndMaps[i]->Id2;

      ZFSBool isLocalAuxMap = false;
      ZFSId localAuxMapId = 0;
      const ZFSUint noAuxDataMaps= m_windowInfo->physicalAuxDataMap.size();

      //loop through all local auxDataMaps
      //to find if this partition shares a part of it
      //and then save the local id
      for(ZFSUint j=0; j<noAuxDataMaps; ++j) {
        stringstream datasetId;
        datasetId<< m_windowInfo->physicalAuxDataMap[j]->Id2;
        if(datasetId.str() == datasetname.str()) {
          isLocalAuxMap = true;
          localAuxMapId = j;
          break;
        }
      }

      if(isLocalAuxMap) {
        stringstream datasetId;
        datasetId<< m_windowInfo->physicalAuxDataMap[localAuxMapId]->Id2;

        ZFSId offset[2];
        ZFSId dataSize[2];
        ZFSId dataSetSize=1;
        ZFSId dim1=0;
        for(ZFSId j=0; j<nDim; ++j){
          if( m_windowInfo->physicalAuxDataMap[localAuxMapId]->start2[j]==m_windowInfo->physicalAuxDataMap[localAuxMapId]->end2[j]) continue;
          dataSize[dim1]=m_windowInfo->physicalAuxDataMap[localAuxMapId]->end2[j]- m_windowInfo->physicalAuxDataMap[localAuxMapId]->start2[j];
          offset[dim1]=m_windowInfo->physicalAuxDataMap[localAuxMapId]->start2[j];
          dataSetSize*=dataSize[dim1];
          dim1++;
        }

        ZFSInt dataSet1[2]={dataSize[1], dataSize[0]};
        ZFSInt offset1[2]={offset[1], offset[0]};
        if(m_bCfCpCoeff) {
          const ZFSId mapOffset = m_cells->cfOffsets[localAuxMapId];
          ZFSString pathName= "window" + datasetId.str() ;
          for(ZFSInt j=0; j<noFields; j++) {io_write_ddataset_part(file, pathName.c_str() ,dataNames[j].c_str() , nDim-1, dataSet1, offset1, &m_cells->cf[mapOffset + dataSetSize*j]);}
          const ZFSId mapOffsetCp = m_cells->cpOffsets[localAuxMapId];
          ZFSString dataname = "cp";
          io_write_ddataset_part(file, pathName.c_str() ,dataname.c_str() , nDim-1, dataSet1, offset1, &m_cells->cp[mapOffsetCp]);
        }

        if(m_bPower) {
          const ZFSId mapOffset = m_cells->powerOffsets[localAuxMapId];
          ZFSString pathName= "window" + datasetId.str() ;
          for(ZFSInt j=0; j<nDim; j++) {
            io_write_ddataset_part(file, pathName.c_str() ,powerNamesVisc[j].c_str() , nDim-1, dataSet1, offset1, &m_cells->powerVisc[mapOffset + dataSetSize*j]);
            io_write_ddataset_part(file, pathName.c_str() ,powerNamesPres[j].c_str() , nDim-1, dataSet1, offset1, &m_cells->powerPres[mapOffset + dataSetSize*j]);
          }
        }
      } else {
        ZFSId offset[2];
        ZFSId dataSize[2];
        ZFSInt ghostArray[3]={m_noGhostLayers,m_noGhostLayers,m_noGhostLayers};
        for(ZFSId j=0; j<nDim-1; ++j){
          offset[j]=0;
          dataSize[j]=0;
        }
        //skin-friction and pressure coefficient
        if(m_bCfCpCoeff){
          ZFSString pathName= "window" + datasetname.str() ;
          for(ZFSInt j=0; j<noFields; j++) {io_write_ddataset_part_ghost_array(file, pathName.c_str(), dataNames[j].c_str() , nDim-1, dataSize, offset, ghostArray, NULL);}
          ZFSString dataname = "cp";
          io_write_ddataset_part_ghost_array(file, pathName.c_str(), dataname.c_str() , nDim-1, dataSize, offset, ghostArray, NULL);
        }
        //power 
        if(m_bPower){
          ZFSString pathName= "window" + datasetname.str() ;
          for(ZFSInt j=0; j<nDim; j++) {
            io_write_ddataset_part_ghost_array(file, pathName.c_str(), powerNamesVisc[j].c_str() , nDim-1, dataSize, offset, ghostArray, NULL);
            io_write_ddataset_part_ghost_array(file, pathName.c_str(), powerNamesPres[j].c_str() , nDim-1, dataSize, offset, ghostArray, NULL);
          }
        }


      }
    }
  }

  if(m_bCl){
    saveLiftCoefficient(file);
  }

  if(m_bCd){
    saveDragCoefficient(file);
  }

  if(m_bPower){
    savePowerCoefficient(file);
  }

  io_closefile(file);
}

/**
 *
 * @author Marian Albers, Nov 15, 2015
 * modified 15.11.2015
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::saveAveragedVariables(ZFSString name, ZFSId noVars, ZFSFloat** summedVars)
{
  //Function to write the solution to file with iolibrary
  stringstream fileName;
  stringstream gridFileName;
  string tempG;

  zfs_log << "writing Averaged Variables to file " << name << m_outputFormat << " ... " << endl;
  fileName << name << m_outputFormat;

  ZFSInt file =io_openfile("hdf5", (fileName.str()).c_str(), "collective", m_zfsStrctrdComm);

  io_create_sattribute(file, "", "UID", m_uID.length());
  io_write_sattribute1(file, "", "UID", m_uID.length(), m_uID.c_str());

  io_create_sattribute(file,"","filetype",9);
  io_write_sattribute1(file,"","filetype",9, "solution");

  io_create_sattribute(file,"","blockType",18);
  io_write_sattribute1(file,"","blockType",18, "ZFS_STRUCTURED_FV");

  io_create_dattribute(file, "", "Ma",1);
  io_write_dattribute1(file, "", "Ma",1, &m_Ma);

  io_create_dattribute(file, "", "Re",1);
  io_write_dattribute1(file, "", "Re",1, &m_Re);

  io_create_dattribute(file, "", "Pr",1);
  io_write_dattribute1(file, "", "Pr",1, &m_Pr);

  io_create_dattribute(file, "", "timeStep",1);
  io_write_dattribute1(file, "", "timeStep",1, &m_timeStep);

  io_create_dattribute(file, "", "time",1);
  io_write_dattribute1(file, "", "time",1, &m_time);

  io_create_dattribute(file, "", "physicalTime",1);
  io_write_dattribute1(file, "", "physicalTime",1, &m_physicalTime);

  io_create_iattribute(file, "", "averageStartTimeStep",1);
  io_write_iattribute1(file, "", "averageStartTimeStep",1, &m_averageStartTimestep);

  io_create_iattribute(file, "", "averageSampleInterval",1);
  io_write_iattribute1(file, "", "averageSampleInterval",1, &m_averageInterval);

  io_create_iattribute(file, "", "noSamples",1);
  io_write_iattribute1(file, "", "noSamples",1, &m_noSamples);


  ZFSInt allCells[3];
  for (ZFSId i=0; i<m_noInputBlocks; i++){
    for(ZFSInt j =0; j<nDim; j++){
      allCells[j]=m_partition->inputBoxInfo[i]->DirLast[j];
    }
    //create datasets for the io library
    stringstream path;
    path << i;
    ZFSString blockPathStr = "block";
    blockPathStr += path.str();
    const char* blockPath = blockPathStr.c_str();

    //create dataset and write
    for(ZFSId var=0; var<noVars; var++) {
      io_create_ddataset(file, blockPath , (m_avgVariableNames[var]).c_str(), nDim, allCells);
    }
  }

  stringstream path;
  path << m_inputBlockId;
  ZFSString blockPathStr = "block";
  blockPathStr += path.str();
  ZFSInt ghostArray[3]={m_noGhostLayers,m_noGhostLayers,m_noGhostLayers};
  for(ZFSId var=0; var<noVars; var++) {
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), (m_avgVariableNames[var]).c_str(), nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&summedVars[var][0] );
  }

  io_closefile(file);
}

/**
 *
 * @author Frederik Temme, Dez 15, 2015
 * modified 12.01.2015
 * equal to saveAveragedVariables below
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::saveAverageRestart(ZFSString name, ZFSId noVars, ZFSFloat** summedVars, ZFSFloat** square, ZFSFloat** cube, ZFSFloat** fourth){

  (void) noVars;

  //Function to write the solution to file with iolibrary
  stringstream fileName;
  stringstream gridFileName;
  string tempG;

  zfs_log << "writing Averaged Restart Variables to file " << name << m_outputFormat  << " ... " << endl;
  fileName << name << m_outputFormat;

  ZFSInt file =io_openfile("hdf5", (fileName.str()).c_str(), "collective", m_zfsStrctrdComm);

  io_create_iattribute(file, "", "averageStartTimeStep",1);
  io_write_iattribute1(file, "", "averageStartTimeStep",1, &m_averageStartTimestep);

  io_create_iattribute(file, "", "averageSampleInterval",1);
  io_write_iattribute1(file, "", "averageSampleInterval",1, &m_averageInterval);

  io_create_iattribute(file, "", "noSamples",1);
  io_write_iattribute1(file, "", "noSamples",1, &m_noSamples);


  ZFSInt allCells[3];//not nDim because of compiler warning0
  for (ZFSId i=0; i<m_noInputBlocks; i++){
    for(ZFSInt j =0; j<nDim; j++){
      allCells[j]=m_partition->inputBoxInfo[i]->DirLast[j];
    }
    //create datasets for the io library
    stringstream path;
    path << i;//m_inputBlockId;
    ZFSString blockPathStr = "block";
    blockPathStr += path.str();
    const char* blockPath = blockPathStr.c_str();

    io_create_ddataset(file, blockPath , "u", 3, allCells);
    io_create_ddataset(file, blockPath , "v", 3, allCells);
    io_create_ddataset(file, blockPath , "w", 3, allCells);
    io_create_ddataset(file, blockPath , "rho", 3, allCells);
    io_create_ddataset(file, blockPath , "p", 3, allCells);
    // io_create_ddataset(file, blockPath , "rans0", 3, allCells); //junoh

    if(m_averageVorticity) {
      io_create_ddataset(file, blockPath , "vortx", 3, allCells);
      io_create_ddataset(file, blockPath , "vorty", 3, allCells);
      io_create_ddataset(file, blockPath , "vortz", 3, allCells);
    }

    io_create_ddataset(file, blockPath , "uu", 3, allCells);
    io_create_ddataset(file, blockPath , "vv", 3, allCells);
    io_create_ddataset(file, blockPath , "ww", 3, allCells);
    io_create_ddataset(file, blockPath , "uv", 3, allCells);
    io_create_ddataset(file, blockPath , "vw", 3, allCells);
    io_create_ddataset(file, blockPath , "uw", 3, allCells);
    io_create_ddataset(file, blockPath , "pp", 3, allCells);

    if(m_averageVorticity) {
      io_create_ddataset(file, blockPath , "vortxvortx", 3, allCells);
      io_create_ddataset(file, blockPath , "vortyvorty", 3, allCells);
      io_create_ddataset(file, blockPath , "vortzvortz", 3, allCells);
    }

    if(m_kurtosis || m_skewness){
      io_create_ddataset(file, blockPath , "uuu", 3, allCells);
      io_create_ddataset(file, blockPath , "vvv", 3, allCells);
      io_create_ddataset(file, blockPath , "www", 3, allCells);
    }

    if(m_kurtosis) {
      io_create_ddataset(file, blockPath , "uuuu", 3, allCells);
      io_create_ddataset(file, blockPath , "vvvv", 3, allCells);
      io_create_ddataset(file, blockPath , "wwww", 3, allCells);
    }
  }

  ////////////////////////////////////////////////
  ///////////// Write Variables //////////////////
  ////////////////////////////////////////////////
  stringstream path;
  path << m_inputBlockId;
  ZFSString blockPathStr = "block";
  blockPathStr += path.str();
  ZFSInt offset = 0;
  ZFSInt ghostArray[3]={m_noGhostLayers,m_noGhostLayers,m_noGhostLayers};
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "u", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&summedVars[0][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "v", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&summedVars[1][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "w", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&summedVars[2][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "rho", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&summedVars[3][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "p", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&summedVars[4][0] );
  // io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "rans0", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&summedVars[5][0] ); //junoh

  // offset = m_maxNoVariables; //junoh
  
  offset = noVariables();

  if(m_averageVorticity) {
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "vortx", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&summedVars[offset+0][0] );
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "vorty", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&summedVars[offset+1][0] );
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "vortz", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&summedVars[offset+2][0] );
  }

  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "uu", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&square[0][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "vv", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&square[1][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "ww", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&square[2][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "uv", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&square[3][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "vw", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&square[4][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "uw", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&square[5][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "pp", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&square[6][0] );

  if(m_averageVorticity) {
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "vortxvortx", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&square[7][0] );
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "vortyvorty", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&square[8][0] );
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "vortzvortz", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&square[9][0] );
  }

  if(m_kurtosis || m_skewness){
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "uuu", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&cube[0][0] );
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "vvv", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&cube[1][0] );
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "www", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&cube[2][0] );
  }

  if(m_kurtosis) {
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "uuuu", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&fourth[0][0] );
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "vvvv", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&fourth[1][0] );
    io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "wwww", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&fourth[2][0] );
  }
  io_closefile(file);
}


/**
 * Writes the production terms into a given file
 * @author Marian Albers, Jan 10, 2017 
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::saveProductionTerms(ZFSString name, ZFSFloat** production)
{
  //Function to write the solution to file with iolibrary
  stringstream fileName;
  stringstream gridFileName;
  string tempG;

  zfs_log << "writing production terms to file " << name << m_outputFormat << " ... " << endl;
  fileName << name << m_outputFormat;

  ZFSInt file =io_openfile("hdf5", (fileName.str()).c_str(), "collective", m_zfsStrctrdComm);

  ZFSInt allCells[3];
  for (ZFSId i=0; i<m_noInputBlocks; i++){
    for(ZFSInt j =0; j<nDim; j++){
      allCells[j]=m_partition->inputBoxInfo[i]->DirLast[j];
    }
    //create datasets for the io library
    stringstream path;
    path << i;
    ZFSString blockPathStr = "block";
    blockPathStr += path.str();
    const char* blockPath = blockPathStr.c_str();

    //create dataset and write
    io_create_ddataset(file, blockPath , "p1j", nDim, allCells);
    io_create_ddataset(file, blockPath , "p2j", nDim, allCells);
    io_create_ddataset(file, blockPath , "p3j", nDim, allCells);
  }

  stringstream path;
  path << m_inputBlockId;
  ZFSString blockPathStr = "block";
  blockPathStr += path.str();
  ZFSInt ghostArray[3]={m_noGhostLayers,m_noGhostLayers,m_noGhostLayers};
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "p1j", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&production[0][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "p2j", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&production[1][0] );
  io_write_ddataset_part_ghost_array(file, blockPathStr.c_str(), "p3j", nDim,  m_nActiveCells, m_nOffsetCells, ghostArray ,&production[2][0] );

  io_closefile(file);
}


/** The readGrid-routine only reads the points in subsequent order without
 *  the ghost layers. Thus the points need to be moved to the right
 *  position to account for the ghost (cells) points
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::moveCellPoints()
{
  //switch command for different treatment depending on the spatial dimensions
  ZFSId i_org, j_org, k_org;
  ZFSId i_new, j_new, k_new;
  ZFSId pointId_org, pointId_new;
  //go over the inner loops of the array without ghost points!!!
  switch(nDim)
  {
  case 1:
  {
    //one dimenionsal case
    zfs_log << "ZFSStrctrdBlck::moveCellsPoints(): 1D case not implemented with Ghost Points (/Cells) " << endl;
    exit(0);
    break;
  }
  case 2:
  {
    //two dimensional case
    for(ZFSId j=(m_nActivePoints[0]-1); j>=0; --j ) {
      for(ZFSId i=(m_nActivePoints[1]-1); i>=0; --i) {
        i_org=i;
        j_org=j;
        i_new=i+m_noGhostLayers;
        j_new=j+m_noGhostLayers;
        pointId_org = i_org+(j_org*m_nActivePoints[1]); //postition in Array
        pointId_new = i_new+(j_new*m_nPoints[1]); //new postion in Array
        for(ZFSId dim=0; dim<nDim; ++dim) {
          m_coordinates[dim][pointId_new] = m_coordinates[dim][pointId_org]; //copy value to the right place
          m_coordinates[dim][pointId_org]= -1000;//for test purposes only
        }
      }
    }
    break;
  }
  case 3:
  {
    for(ZFSId k=(m_nActivePoints[0]-1); k>=0; k-- ) {
      for(ZFSId j=(m_nActivePoints[1]-1); j>=0; j-- ) {
        for(ZFSId i=(m_nActivePoints[2]-1); i>=0; i--) {
          i_org=i;
          j_org=j;
          k_org=k;
          i_new=i_org+m_noGhostLayers;
          j_new=j_org+m_noGhostLayers;
          k_new=k_org+m_noGhostLayers;
          pointId_org = i_org+(j_org+k_org*m_nActivePoints[1])*m_nActivePoints[2]; //postition in Array

          pointId_new = i_new+(j_new+k_new*m_nPoints[1])*m_nPoints[2]; //new postion in Array

          for(ZFSId dim=0; dim<nDim; ++dim) {
            m_coordinates[dim][pointId_new] = m_coordinates[dim][pointId_org]; //copy value to the right place
            m_coordinates[dim][pointId_org] = F0;
          }
        }
      }
    }

    break;
  }
  default:
  {
    zfs_log<< "number of dimensions \"" << nDim << "\" has not been implemented" << endl;
    exit(0);
    break;
  }
  }
}


template <ZFSInt nDim> //junoh
void ZFSStrctrdBlck<nDim>::exchange(){
  zfsTerm(1,__CALLING_FUNCTION__, "Derived Function was not implemented ==> Have a look at Block type");
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::exchangeTimeStep()
{
  TRACE();  
  RECORD_TIMER_START(m_tcomm);
  RECORD_TIMER_START(m_texchangeDt);

  ZFSFloat globalTimeStepMin =F0;

  MPI_Allreduce(&m_timeStep, &globalTimeStepMin, 1, MPI_DOUBLE, MPI_MIN, m_zfsStrctrdComm);

  m_timeStep=globalTimeStepMin;


  RECORD_TIMER_STOP(m_texchangeDt);
  RECORD_TIMER_STOP(m_tcomm);
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::exchangePoints()
{
  gatherPoints();
  sendPoints();
  receivePoints();

  for(ZFSId i= 0; i< m_cmnctnFlag->noNghbrDomainsNormal-m_cmnctnFlag->noNghbrDomainsSingular; i++) {
    MPI_Wait(&(m_cmnctnFlag->mpi_sndRequest[i]), &(m_cmnctnFlag->mpi_sndStatus[i]));
  }

  for(ZFSId i=0; i< m_cmnctnFlag->noNghbrDomainsNormal-m_cmnctnFlag->noNghbrDomainsSingular; i++) {
    MPI_Wait(&(m_cmnctnFlag->mpi_rcvRequest[i]), &(m_cmnctnFlag->mpi_rcvStatus[i]));
  }

  scatterPoints();
}


/**
 *   Householder Reduction according to
 *     Numercial Recipies in C: The Art of Scientific Computing
 *     Authors:
 *
 *     function tred2:
 *     Hosholder reduction of a real, symmetric matrix A.
 *     On output A is replaced by the orthogonla matrix Q (omitted here)
 *     diag  returns the diagonla elements of the tridiagonal matrix
 *     offdiag the off-diagonal elements with offdiag[1]=0;
 */


//compute Householder reduction
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::tred2(ZFSFloatScratchSpace& A, ZFSInt dim, ZFSFloat* diag, ZFSFloat* offdiag )
{

  //formulation from book pp. 582
  ZFSInt n=dim;
  ZFSId l=0,k=0,j=0,i=0;
  ZFSFloat scale=F0, hh=F0, h=F0, g=F0, f=F0;

  for(i=dim-1; i>0; i--)
  {
    l=i-1;
    h=F0;
    scale=F0;
    if(l>0)
    {
      for(k=0; k<i; ++k) { scale += abs(A(i,k));}
      if(approx(scale,F0, m_eps)) { offdiag[i]=A(i,l);}
      else
      {
        for (k=0; k<i;++k)
        {
          A(i,k)/=scale;
          h+=A(i,k)*A(i,k);
        }
        f=A(i,l);
        g=(f>=F0 ? -sqrt(h) : sqrt(h));
        offdiag[i]=scale*g;
        h-=f*g;
        A(i,l)=f-g;
        f=F0;
        for(j=0;j<i;++j)
        {
          g=F0;
          for(k=0;k<j+1;++k) { g+= A(j,k)*A(i,k);}
          for(k=j+1;k<i;++k) { g+= A(k,j)*A(i,k);}
          offdiag[j]=g/h;
          f+= offdiag[j]*A(i,j);
        }
        hh=f/(h+h);
        for(j=0;j<i;j++)
        {
          f=A(i,j);
          g=offdiag[j]-hh*f;
          offdiag[j]=g;
          for(k=0;k<j+1;k++) {A(j,k)-=f*offdiag[k]+g*A(i,k);}
        }
      }
    }
    else
    {
      offdiag[i]=A(i,l);
    }
    diag[i]=h;
  }

  offdiag[0]=F0;
  for(i=0;i<n; ++i)
  {
    diag[i]=A(i,i);
  }

}


/**
 *   Compute Eigenvalues with implicit shift according to
 *     Numercial Recipies in C: The Art of Scientific Computing
 *     Authors:
 *
 *     function tqli2:
 *      QL algorithm with implicit shifts, to determine the eigenvalues of
 *      a real symmetric matrix previously reduced by tred2. diag is a vector
 *      of length np. On input, its first n elements are the diagonal
 *      elements of the tridiagonal matrix. On output, it returns the
 *      eigenvalues. The vector offdiag inputs the subdiagonal elements of the
 *      tridiagonal matrix, with offidag(0) arbitrary. On output offdiag
 *      is destroyed.
 *
 */

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::tqli2(ZFSFloat* diag, ZFSFloat* offdiag, ZFSInt dim)//, ZFSFloat** z)
{


  const ZFSFloat eps = numeric_limits<ZFSFloat>::epsilon();
  ZFSInt m,l,iter,i,k;
  ZFSFloat s,r,p,g,f,dd,c,b;
  ZFSInt n=dim;
  for (i=1;i<n;++i) offdiag[i-1]=offdiag[i];
  offdiag[n-1]=0.0;
  for (l=0;l<n;++l) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
        dd=fabs(diag[m])+fabs(diag[m+1]);
        if (fabs(offdiag[m])<=eps*dd) break;
      }
      if (m != l) {
        if (iter++ == 30)
        {
          for(k=0; k<dim; ++k)
          {
            diag[k]=F0;
          }
          return;
        }
        g=(diag[l+1]-diag[l])/(2.0*offdiag[l]);
        r=pythag(g,F1);
        r=abs(r);
        if(g<F0) r*=-F1;
        g=diag[m]-diag[l]+offdiag[l]/(g+r);
        s=c=F1;
        p=F0;
        for (i=m-1;i>=l;i--) {
          f=s*offdiag[i];
          b=c*offdiag[i];
          offdiag[i+1]=(r=pythag(f,g));
          if (approx(r, F0,eps)) {
            diag[i+1] -= p;
            offdiag[m]=F0;
            break;
          }
          s=f/r;
          c=g/r;
          g=diag[i+1]-p;
          r=(diag[i]-g)*s+2.0*c*b;
          diag[i+1]=g+(p=s*r);
          g=c*r-b;
        }
        if ( approx(r, F0, eps) && i >= l) continue;
        diag[l] -= p;
        offdiag[l]=g;
        offdiag[m]=F0;
      }
    } while (m != l);
  }
}

/**
 *  Sorting function to sort list in ascending order
 *
 *  function: insertSort
 *
 *  sorts an array list into ascending numerical order, by straight
 *  insertion. dim is input(size of list); list is replaced on
 *  output by its sorted rearrangement.
 *
 */

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::insertSort(ZFSId dim, ZFSFloat* list)
{
  ZFSFloat temp=F0;
  ZFSInt j=0;
  for(ZFSId i=1; i<dim; i++)
  {
    temp = list[i];
    j = i-1;
    while(j>=0 && temp<list[j])
    {
      list[j+1] = list[j];
      j = j-1;
    }
    list[j+1] = temp;
  }
}

template <ZFSInt nDim>
ZFSFloat ZFSStrctrdBlck<nDim>::pythag(ZFSFloat a, ZFSFloat b)
{
  ZFSFloat absa=F0,absb=F0;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb)
    return absa*sqrt(1.0 + POW2(absb/absa));
  else
    return (approx(absb, F0, m_eps) ? F0 : absb*sqrt(1.0 + POW2(absa/absb)));
}


template <ZFSInt nDim>
ZFSFloat ZFSStrctrdBlck<nDim>::bessJ0 (ZFSFloat x) {
  // This subroutine calculates the First Kind Bessel Function of
  // order 0, for any real number X. The polynomial approximation by
  // series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
  // REFERENCES:
  // M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
  // C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
  // VOL.5, 1962.

  const ZFSFloat
    p1=1.0, p2=-0.1098628627e-2, p3=0.2734510407e-4,
    p4=-0.2073370639e-5, p5= 0.2093887211e-6,
    q1=-0.1562499995e-1, q2= 0.1430488765e-3, q3=-0.6911147651e-5,
    q4= 0.7621095161e-6, q5=-0.9349451520e-7,
    r1= 57568490574.0, r2=-13362590354.0, r3=651619640.7,
    r4=-11214424.18, r5= 77392.33017, r6=-184.9052456,
    s1= 57568490411.0, s2=1029532985.0, s3=9494680.718,
    s4= 59272.64853, s5=267.8532712, s6=1.0;
  ZFSFloat ax,fr,fs,z,fp,fq1,xx,y, tmp;

  if (approx(x,0.0,m_eps)) { return 1.0; };

  ax = fabs(x);
  if (ax < 8.0) {
    y = x*x;
    fr = r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))));
    fs = s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))));
    tmp = fr/fs;
  } else {
    z = 8./ax;
    y = z*z;
    xx = ax-0.785398164;
    fp = p1+y*(p2+y*(p3+y*(p4+y*p5)));
    fq1 = q1+y*(q2+y*(q3+y*(q4+y*q5)));
    tmp = sqrt(0.636619772/ax)*(fp*cos(xx)-z*fq1*sin(xx));
  }
  return tmp;
}

template <ZFSInt nDim>
ZFSFloat ZFSStrctrdBlck<nDim>::sign(ZFSFloat x, ZFSFloat y) {
  if (y<0.0) return (-fabs(x));
  else return (fabs(x));
}

template <ZFSInt nDim>
ZFSFloat ZFSStrctrdBlck<nDim>::bessJ1 (ZFSFloat x) {
  // This subroutine calculates the First Kind Bessel Function of
  // order 1, for any real number X. The polynomial approximation by
  // series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
  // REFERENCES:
  // M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
  // C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
  // VOL.5, 1962.

  const ZFSFloat
    p1=1.0, p2=0.183105e-2, p3=-0.3516396496e-4, p4=0.2457520174e-5,
    p5=-0.240337019e-6,  p6=0.636619772,
    q1= 0.04687499995, q2=-0.2002690873e-3, q3=0.8449199096e-5,
    q4=-0.88228987e-6, q5= 0.105787412e-6,
    r1= 72362614232.0, r2=-7895059235.0, r3=242396853.1,
    r4=-2972611.439,   r5=15704.48260,  r6=-30.16036606,
    s1=144725228442.0, s2=2300535178.0, s3=18583304.74,
    s4=99447.43394,    s5=376.9991397,  s6=1.0;

  ZFSFloat ax,fr,fs,y,z,fp,fq1,xx, tmp;

  ax = fabs(x);
  if (ax < 8.0) {
    y = x*x;
    fr = r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))));
    fs = s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))));
    tmp = x*(fr/fs);
  } else {
    z = 8.0/ax;
    y = z*z;
    xx = ax-2.35619491;
    fp = p1+y*(p2+y*(p3+y*(p4+y*p5)));
    fq1 = q1+y*(q2+y*(q3+y*(q4+y*q5)));
    tmp = sqrt(p6/ax)*(cos(xx)*fp-z*sin(xx)*fq1)*sign(s6,x);
  }
  return tmp;
}

template <ZFSInt nDim>
ZFSFloat ZFSStrctrdBlck<nDim>::bessJ (ZFSInt n, ZFSFloat x) {
  // This subroutine calculates the first kind modified Bessel function
  // of integer order N, for any REAL X. We use here the classical
  // recursion formula, when X > N. For X < N, the Miller's algorithm
  // is used to avoid overflows.
  // REFERENCE:
  // C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
  // MATHEMATICAL TABLES, VOL.5, 1962.

  const ZFSInt iacc = 40;
  const ZFSFloat bigno = 1e10,  bigni = 1e-10;

  ZFSFloat tox,bjm,bj,bjp,sum,tmp;
  ZFSInt j, jsum, m;

  if (n == 0) {return bessJ0(x);};
  if (n == 1) {return bessJ1(x);};
  if (approx(x,0.0,m_eps)) {return 0.0;};

  tox = 2.0/x;
  if (x > 1.0*n) {
    bjm = bessJ0(x);
    bj  = bessJ1(x);
    for (j=1; j<n; j++) {
      bjp = j*tox*bj-bjm;
      bjm = bj;
      bj  = bjp;
    }
    return bj;
  }
  else {
    m = (ZFSInt) (2*((n+floor(sqrt(1.0*(iacc*n))))/2));
    tmp = 0.0;
    jsum = 0;
    sum = 0.0;
    bjp = 0.0;
    bj  = 1.0;
    for (j=m; j>0; j--) {
      bjm = j*tox*bj-bjp;
      bjp = bj;
      bj  = bjm;
      if (fabs(bj) > bigno) {
        bj  = bj*bigni;
        bjp = bjp*bigni;
        tmp = tmp*bigni;
        sum = sum*bigni;
      }
      if (jsum != 0)  sum += bj;
      jsum = 1-jsum;
      if (j == n)  tmp = bjp;
    }
    sum = 2.0*sum-bj;
    return (tmp/sum);
  }
}

template <ZFSInt nDim>
ZFSFloat ZFSStrctrdBlck<nDim>::bessI(ZFSInt n, ZFSFloat x) {

  // This subroutine calculates the first kind modified Bessel function
  // of integer order N, for any REAL X. We use here the classical
  // recursion formula, when X > N. For X < N, the Miller's algorithm
  // is used to avoid overflows.
  // REFERENCE:
  // C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
  //  MATHEMATICAL TABLES, VOL.5, 1962.


  ZFSInt iacc = 40;
  ZFSFloat bigno = 1e10, bigni = 1e-10;
  ZFSFloat tox, bim, bi, bip, bsi;
  ZFSInt j, m;

  if (n==0)  return (bessI0(x));
  if (n==1)  return (bessI1(x));
  if (approx(x,0.0,m_eps)) return 0.0;

  tox = 2.0/x;
  bip = 0.0;
  bi  = 1.0;
  bsi = 0.0;
  m = (ZFSInt) (2*((n+floor(sqrt(iacc*n)))));
  for (j = m; j>0; j--) {
    bim = bip+j*tox*bi;
    bip = bi;
    bi  = bim;
    if (fabs(bi) > bigno) {
      bi  = bi*bigni;
      bip = bip*bigni;
      bsi = bsi*bigni;
    }
    if (j==n)  bsi = bip;
  }
  return (bsi*bessI0(x)/bi);
}

//  Auxiliary Bessel functions for N=0
template <ZFSInt nDim>
ZFSFloat ZFSStrctrdBlck<nDim>::bessI0(ZFSFloat x) {
  ZFSFloat y,p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,ax,bx;
  p1=1.0; p2=3.5156229; p3=3.0899424; p4=1.2067492;
  p5=0.2659732; p6=0.360768e-1; p7=0.45813e-2;
  q1=0.39894228; q2=0.1328592e-1; q3=0.225319e-2;
  q4=-0.157565e-2; q5=0.916281e-2; q6=-0.2057706e-1;
  q7=0.2635537e-1; q8=-0.1647633e-1; q9=0.392377e-2;
  if (fabs(x) < 3.75) {
    y=(x/3.75)*(x/3.75);
    return (p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
  }
  else {
    ax=fabs(x);
    y=3.75/ax;
    bx=exp(ax)/sqrt(ax);
    ax=q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9)))))));
    return (ax*bx);
  }
}

//  Auxiliary Bessel functions for N=1
template <ZFSInt nDim>
ZFSFloat ZFSStrctrdBlck<nDim>::bessI1(ZFSFloat x) {
  ZFSFloat y,p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,ax,bx;
  p1=0.5; p2=0.87890594; p3=0.51498869; p4=0.15084934;
  p5=0.2658733e-1; p6=0.301532e-2; p7=0.32411e-3;
  q1=0.39894228; q2=-0.3988024e-1; q3=-0.362018e-2;
  q4=0.163801e-2; q5=-0.1031555e-1; q6=0.2282967e-1;
  q7=-0.2895312e-1; q8=0.1787654e-1; q9=-0.420059e-2;
  if (fabs(x) < 3.75) {
    y=(x/3.75)*(x/3.75);
    return(x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))));
  } else {
    ax=fabs(x);
    y=3.75/ax;
    bx=exp(ax)/sqrt(ax);
    ax=q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9)))))));
    return (ax*bx);
  }
}


template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::allocateMetrics()
{
  TRACE();
  zfs_log << "allocating metrics..." << endl;
  allocateSurfaceMetrics();
  allocateCornerMetrics();
  allocateCellMetrics();
  zfs_log << "allocating metrics... SUCCESSFUL" << endl;
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::allocateSurfaceMetrics()
{
  TRACE();
  zfs_log << "allocating surface metrics ..." << endl;

  ZFSId no_metrics=9;
  ZFSId no_cells=1; // else multiplying with zero in for-loop

  if(nDim <3){
    no_metrics= 8;
  }

  for(ZFSId i=0; i<nDim; i++){
    no_cells*=m_nCells[i];
  }

  // 8 or 9  metrics terms
  // initialize with unrealistic values
  zfsAlloc(this->m_cells->surfaceMetrics, no_cells, no_metrics, "surfaceMetrics", 123.123, __CALLING_FUNCTION__);

  zfs_log << "allocating surface metrics ... SUCCESSFUL" << endl;
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::allocateCornerMetrics()
{
  TRACE();
  zfs_log << "allocating corner metrics ..." << endl;

  ZFSId no_metrics=9;
  ZFSId no_cells=1;// else multiplying with zero in for-loop

  if(nDim <3){
    no_metrics= 4;
  }

  for(ZFSId i=0; i<nDim; i++){
    no_cells*=m_nCells[i];
  }

  // 4 or 9 metric terms per cell, all evaluated at the point ipjp(kp)
  // initialize with unrealistic values
  zfsAlloc(this->m_cells->cornerMetrics, no_cells*no_metrics, "cornerMetrics", 123.123, __CALLING_FUNCTION__);

  zfs_log << "allocating corner metrics ... SUCCESSFUL" << endl;
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::allocateCellMetrics()
{
  TRACE();
  zfs_log << "allocating cell metrics ..." << endl;

  ZFSId no_metrics=9;
  ZFSId no_cells=1;// else multiplying with zero in for-loop

  if(nDim <3){
    no_metrics= 4;
  }

  for(ZFSId i=0; i<nDim; i++){
    no_cells*=m_nCells[i];
  }

  // 4 or 9 metric terms per cell, all evaluated at the point ipjp(kp)
  // initialize with unrealistic values
  zfsAlloc(m_cells->cellMetrics, no_cells, no_metrics, "cellMetrics", 123.123, __CALLING_FUNCTION__);

  zfs_log << "allocating cell metrics ... SUCCESSFUL" << endl;
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::allocateJacobian(){

  TRACE();
  zfs_log << "allocating jacobian..." << endl;

  ZFSId no_cells=1;// else multiplying with zero in for-loop

  for(ZFSId dim=0; dim<nDim; dim++){
    no_cells*=m_nCells[dim];
  }


  zfsAlloc(m_cells->cornerJac, no_cells, "m_cells->cornerJac", __CALLING_FUNCTION__);
  zfsAlloc(m_cells->cellJac, no_cells, "m_cells->cellJac", __CALLING_FUNCTION__);
  zfsAlloc(m_cells->oldCellJac, no_cells , "m_cells->oldCellJac", __CALLING_FUNCTION__);

  zfs_log << "allocating jacobian... SUCCESSFUL" << endl;
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::computeMetrics(){

  TRACE();
  RECORD_TIMER_START(m_tcomputeMetrics);
  //computeSurfaceMetrics();
  computeModSurfaceMetrics();
  computeCornerMetrics();
  //computeModCornerMetrics();
  computeCellMetrics();
  RECORD_TIMER_STOP(m_tcomputeMetrics);
}


template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::computeJacobian(){

  TRACE();
  RECORD_TIMER_START(m_tcomputeJacobian);
  computeCornerJacobian();
  //computeCellJacobian();
  //computeModCornerJacobian();
  computeModCellJacobian();
  RECORD_TIMER_STOP(m_tcomputeJacobian);
}

/**
 * Checks whole domain for NaNs and adds the number of NaNs globally
 * \author Marian Albers
 * \date 19.08.2015
 */
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::checkNans()
{
  TRACE();
  ZFSInt noNansLocal = 0;
  for(ZFSId cellid=0; cellid<m_noStrctrdCells; cellid++) {
    //go through every cell
    for( ZFSId i = 0; i < CV->noVariables; i++ ) {
      if(std::isnan(m_cells->variables[i][cellid])) {
        noNansLocal++;

      }
    }
  }

  ZFSInt noNansGlobal = 0;
  MPI_Allreduce(&noNansLocal, &noNansGlobal, 1, MPI_INT, MPI_SUM, m_zfsStrctrdComm);

  if(domainId() == 0) {
    cout << "GlobalTimeStep: " << globalTimeStep << " , noNansGlobal: " << noNansGlobal << endl;
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::setTimeStep()
{
  TRACE();

  if((!m_constantTimeStep && globalTimeStep % m_timeStepComputationInterval == 0) || globalTimeStep == 0) {
    computeTimeStep();
    if(noDomains()>1) {
      exchangeTimeStep();
    }

    m_physicalTimeStep = m_timeStep*m_timeRef;
  }

  //if (m_movingGrid && m_travelingWave) {
  //  fixTimeStepTravelingWave();
  //}
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::fixTimeStepTravelingWave()
{
  const ZFSFloat waveTransSpeed = m_waveSpeed;

  ZFSFloat deltaZ = F0;
  if(m_gridMovingMethod == 9 || m_gridMovingMethod == 11) {
    deltaZ = abs(m_coordinates[2][0]-m_coordinates[2][m_nPoints[2]*m_nPoints[1]]);
  } else {
    zfs_log << "WARNING!!!!!!!!!!!!!!!!!!!!: dz was fixed and hardcoded !!!!!!!!!! " << endl;
    deltaZ = 0.33200866159432146;
  }

  if(!m_waveTimeStepComputed) {
    //first compute time step as usual and exchange it
    computeTimeStep();
    if(noDomains()>1) {
      exchangeTimeStep();
    }

    //number of timeSteps the wave needs to move one cell width further
    m_waveNoStepsPerCell = ceil((deltaZ/waveTransSpeed)/(m_timeStep)); 
    if (m_waveNoStepsPerCell < 2) {
      m_waveNoStepsPerCell = 2; // although highly improbable, see Nyquist-Shannon-Theorem
    }

    cout.precision(18);
    if(domainId()==0) {
      cout << "Old time step = " << m_timeStep << endl;
    }
    zfs_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
    zfs_log << "Old time step = " << m_timeStep << endl;

    m_timeStep = (deltaZ/waveTransSpeed)/(m_waveNoStepsPerCell);
    m_physicalTimeStep = m_timeStep*m_timeRef;

    ZFSBool syncSolutionInterval=false;
    ZFSBool syncDragInterval=false;
    ZFSBool syncPlaneInterval=false;
    ZFSBool syncBoxInterval=false;
    ZFSBool syncLineInterval=false;
    
    
    //fix output writing functions in case 
    if(m_synchronizedMGOutput){
      //change the output frequency of the different io routines to be in accordance with the new time step such that only in full number of iterations in which the moving grid has traveld one cell distance is ensured
      //if(m_outputOffset){
      //  zfs_log << "ERROR: m_outputOffset and synchronized output cannot be combined (not implemented yet)" << endl;
      //  zfsTerm(1, __CALLING_FUNCTION__, "ERROR: m_outputOffset and synchronized output cannot be combined (not implemented yet)");
      //}
      m_outputOffset=m_movingGridStepOffset; //for input output such that it works propberly
      if(m_useConvectiveUnitWrite){
        zfs_log << "ERROR: m_useConvectiveUnitWrite and synchronized output cannot be combined (not implemented yet)" << endl;
        zfsTerm(1, __CALLING_FUNCTION__, "ERROR: m_useConvetiveUnitWrite and synchronized output cannot be combined (not implemented yet)");
      }
      if(m_solutionInterval){
        m_solutionInterval = m_waveNoStepsPerCell*floor((ZFSFloat)m_solutionInterval/(ZFSFloat)m_waveNoStepsPerCell);
        syncSolutionInterval=true;
      }
      if(m_dragOutputInterval){
        m_dragOutputInterval = m_waveNoStepsPerCell*floor((ZFSFloat)m_dragOutputInterval/(ZFSFloat)m_waveNoStepsPerCell);
        syncDragInterval =true;
      }
      if(m_planeOutputInterval){
        m_planeOutputInterval = m_waveNoStepsPerCell*floor((ZFSFloat)m_planeOutputInterval/(ZFSFloat)m_waveNoStepsPerCell);
        syncPlaneInterval = true;
      }
      if(m_boxOutputInterval){
        m_boxOutputInterval = m_waveNoStepsPerCell*floor((ZFSFloat)m_boxOutputInterval/(ZFSFloat)m_waveNoStepsPerCell);
        syncBoxInterval = true;
      }
      if(m_lineOutputInterval){
        m_lineOutputInterval = m_waveNoStepsPerCell*floor((ZFSFloat)m_boxOutputInterval/(ZFSFloat)m_waveNoStepsPerCell);
        syncLineInterval=true;
      }
      
      
    }

    if(domainId()==0) {
      cout << "New time step: " << m_timeStep << " new physical time step: " << m_physicalTimeStep << endl;
      cout << "Number of steps to move one cell width: " << deltaZ/(m_waveSpeed*m_timeStep) << " time steps" << endl;
      cout << "Number of steps to move one physical time step: " << F1/m_physicalTimeStep << endl;
      cout << "Number of steps to move one wave length: "<< m_waveLength/(m_waveSpeed*m_timeStep) << endl;
      cout << "Cells per wavelength: " << m_waveLength/deltaZ << " deltaZ: " << deltaZ << endl;
      cout << "Time for wave to travel one wave length: " << m_waveLength/m_waveSpeed << endl;
      cout << "Physical time for wave to travel one wave length: " << m_waveLength/m_waveSpeed*m_timeRef << endl;
      cout << "solution output interval was reset: " << syncSolutionInterval << " and changed to " << m_solutionInterval << endl;
      cout << "plane output interval was reset: " << syncPlaneInterval << " and changed to " << m_planeOutputInterval << endl;
      cout << "box output interval was reset: " << syncBoxInterval << " and changed to " << m_boxOutputInterval << endl;
      cout << "drag output interval was reset: " << syncDragInterval << " and changed to " << m_dragOutputInterval << endl;
      cout << "line output interval was reset: " << syncLineInterval << " and changed to " << m_lineOutputInterval << endl;
    }

    zfs_log << "New time step: " << m_timeStep << " new physical time step: " << m_physicalTimeStep << endl;
    zfs_log << "Number of steps to move one cell width: " << deltaZ/(m_waveSpeed*m_timeStep) << " time steps" << endl;
    zfs_log << "Number of steps to move one physical time step: " << F1/m_physicalTimeStep << endl;
    zfs_log << "Number of steps to move one wave length: "<< m_waveLength/(m_waveSpeed*m_timeStep) << endl;
    zfs_log << "Cells per wavelength: " << m_waveLength/deltaZ << " deltaZ: " << deltaZ << endl;
    zfs_log << "Time for wave to travel one wave length: " << m_waveLength/m_waveSpeed << endl;
    zfs_log << "Physical time for wave to travel one wave length: " << m_waveLength/m_waveSpeed*m_timeRef << endl;
    zfs_log << "solution output interval was reset: " << syncSolutionInterval << " and changed to " << m_solutionInterval << endl;
    zfs_log << "plane output interval was reset: " << syncPlaneInterval << " and changed to " << m_planeOutputInterval << endl;
    zfs_log << "box output interval was reset: " << syncBoxInterval << " and changed to " << m_boxOutputInterval << endl;
    zfs_log << "drag output interval was reset: " << syncDragInterval << " and changed to " << m_dragOutputInterval << endl;
    zfs_log << "line output interval was reset: " << syncLineInterval << " and changed to " << m_lineOutputInterval << endl;
    zfs_log << "solution writing out will start at" << m_outputOffset << endl;
    zfs_log << "////////////////////////////////////////////////////////////////" << endl;

    m_waveTimeStepComputed = true;
  }

  if(globalTimeStep==0||globalTimeStep==m_restartTimeStep) {
    ZFSBool syncSolutionInterval=false;
    ZFSBool syncDragInterval=false;
    ZFSBool syncPlaneInterval=false;
    ZFSBool syncBoxInterval=false;
    ZFSBool syncLineInterval=false;
    
    //fix output writing functions in case 
    if(m_synchronizedMGOutput){
      //change the output frequency of the different io routines to be in accordance with the new time step such that only in full number of iterations in which the moving grid has traveld one cell distance is ensured
      //if(m_outputOffset!=m_movingGridStepOffset){
      //  zfs_log << "ERROR: m_outputOffset and synchronized output cannot be combined (not implemented yet)" << endl;
      //  zfsTerm(1, __CALLING_FUNCTION__, "ERROR: m_outputOffset and synchronized output cannot be combined (not implemented yet)");
      //}
      //shift the starting point of the outputs for the moving grid
      m_outputOffset=m_movingGridStepOffset; //for input output such that it works propberly
      if(m_useConvectiveUnitWrite){
        zfs_log << "ERROR: m_useConvectiveUnitWrite and synchronized output cannot be combined (not implemented yet)" << endl;
        zfsTerm(1, __CALLING_FUNCTION__, "ERROR: m_useConvetiveUnitWrite and synchronized output cannot be combined (not implemented yet)");
      }
      if(m_solutionInterval){
        m_solutionInterval = m_waveNoStepsPerCell*floor((ZFSFloat)m_solutionInterval/(ZFSFloat)m_waveNoStepsPerCell);
        syncSolutionInterval=true;
      }
      if(m_dragOutputInterval){
        m_dragOutputInterval = m_waveNoStepsPerCell*floor((ZFSFloat)m_dragOutputInterval/(ZFSFloat)m_waveNoStepsPerCell);
        syncDragInterval =true;
      }
      if(m_planeOutputInterval){
        m_planeOutputInterval = m_waveNoStepsPerCell*floor((ZFSFloat)m_planeOutputInterval/(ZFSFloat)m_waveNoStepsPerCell);
        syncPlaneInterval = true;
      }
      if(m_boxOutputInterval){
        m_boxOutputInterval = m_waveNoStepsPerCell*floor((ZFSFloat)m_boxOutputInterval/(ZFSFloat)m_waveNoStepsPerCell);
        syncBoxInterval = true;
      }
      if(m_lineOutputInterval){
        m_lineOutputInterval = m_waveNoStepsPerCell*floor((ZFSFloat)m_boxOutputInterval/(ZFSFloat)m_waveNoStepsPerCell);
        syncLineInterval=true;
      } 
    }

    zfs_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
    zfs_log << "New time step: " << m_timeStep << " new physical time step: " << m_physicalTimeStep << endl;
    zfs_log << "Number of steps to move one cell width: " << deltaZ/(m_waveSpeed*m_timeStep) << " time steps" << endl;
    zfs_log << "Number of steps to move one physical time step: " << F1/m_physicalTimeStep << endl;
    zfs_log << "Number of steps to move one wave length: "<< m_waveLength/(m_waveSpeed*m_timeStep) << endl;
    zfs_log << "Time for wave to travel one wave length: " << m_waveLength/m_waveSpeed << endl;
    zfs_log << "Cells per wavelength: " << m_waveLength/deltaZ << " deltaZ: " << deltaZ << endl;
    zfs_log << "solution output interval was reset: " << syncSolutionInterval << " and changed to " << m_solutionInterval << endl;
    zfs_log << "plane output interval was reset: " << syncPlaneInterval << " and changed to " << m_planeOutputInterval << endl;
    zfs_log << "box output interval was reset: " << syncBoxInterval << " and changed to " << m_boxOutputInterval << endl;
    zfs_log << "drag output interval was reset: " << syncDragInterval << " and changed to " << m_dragOutputInterval << endl;
    zfs_log << "line output interval was reset: " << syncLineInterval << " and changed to " << m_lineOutputInterval << endl;
    zfs_log << "solution writing out will start at" << m_outputOffset << endl;
    zfs_log << "////////////////////////////////////////////////////////////////" << endl;
    //correct Averaging Time Steps
    if(m_postprocessing) {
      const ZFSInt waveNoStepsPerCell = m_waveNoStepsPerCell;
      if(m_averageInterval%waveNoStepsPerCell!=0) {
        zfs_log << "Changed averageInterval from " << m_averageInterval;
        m_averageInterval = waveNoStepsPerCell*floor((ZFSFloat)m_averageInterval/(ZFSFloat)waveNoStepsPerCell);
        zfs_log << " to " << m_averageInterval << ", every " << m_averageInterval*m_physicalTimeStep << " convective units" <<  endl;
      }

      const ZFSInt waveStepOffset = m_movingGridStepOffset;
      if((m_averageStartTimestep-waveStepOffset)%m_averageInterval!=0) {
        zfs_log << "Changed averageStartTimeStep from " << m_averageStartTimestep;
        ZFSInt offsetCounter = ceil(((ZFSFloat)m_averageStartTimestep-(ZFSFloat)waveStepOffset)/(ZFSFloat)m_averageInterval);
        m_averageStartTimestep = waveStepOffset+offsetCounter*m_averageInterval;
        zfs_log << " to " << m_averageStartTimestep << ". Averaging every " << m_averageInterval << " time steps"<< endl;
      }

      if((m_averageStopTimestep-waveStepOffset)%m_averageInterval!=0) {
        zfs_log << "Changed averageStopTimeStep from " << m_averageStopTimestep;
        ZFSInt offsetCounter = ceil(((ZFSFloat)m_averageStopTimestep-(ZFSFloat)waveStepOffset)/(ZFSFloat)m_averageInterval);
        m_averageStopTimestep = waveStepOffset+offsetCounter*m_averageInterval;
        zfs_log << " to " << m_averageStopTimestep << ". Averaging every " << m_averageInterval << " time steps"<< endl;
      }
    }
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::convertRestartVariables(ZFSFloat oldMa){
  TRACE();
  const  ZFSFloat eps = pow( 10.0, -5.0 );
  if(ABS(oldMa-m_Ma)>eps){
    zfs_log << "converting restart variables from old Ma: " << oldMa
            << " to new Ma: " << m_Ma << " ..." << endl;
    const ZFSFloat gammaMinusOne = m_gamma - 1.0;
    //old references
    ZFSFloat T8old = 1.0 / ( 1.0 + 0.5 * gammaMinusOne * POW2(oldMa));
    ZFSFloat p8old = pow( T8old , (m_gamma / gammaMinusOne)) / m_gamma;
    ZFSFloat u8old = oldMa * sqrt( T8old );
    ZFSFloat rho8old = pow( T8old, ( 1.0 / gammaMinusOne ) );
    // ZFSFloat rhoU8old = rho8old*u8old;
     ZFSFloat rhoE8old = p8old / gammaMinusOne + rho8old*( F1B2*POW2(u8old) );
    //new references
    ZFSFloat T8new = 1.0 / ( 1.0 + 0.5 * gammaMinusOne * POW2(m_Ma));
    ZFSFloat p8new = pow( T8new , (m_gamma / gammaMinusOne)) / m_gamma;
    ZFSFloat u8new = m_Ma * sqrt( T8new );
    ZFSFloat rho8new = pow( T8new, ( 1.0 / gammaMinusOne ) );
    // ZFSFloat rhoU8new = rho8new*u8new;
     ZFSFloat rhoE8new = p8new / gammaMinusOne + rho8new*( F1B2*POW2(u8new) );
    //ratios
    ZFSFloat velRatio  = u8new/u8old;
    ZFSFloat rhoRatio  = rho8new/rho8old;
    ZFSFloat pRatio = p8new/p8old;
    ZFSFloat rhoERatio = rhoE8new/rhoE8old;

    //conversion
    for(ZFSId cellId=0; cellId < m_noStrctrdCells; ++cellId){
      //density
      m_cells->pvariables[PV->RHO][cellId]*=rhoRatio;
      //velocities
      for(ZFSId i=0; i<nDim; ++i){
        m_cells->pvariables[PV->VV[i]][cellId]*=velRatio;
      }
      //energy
      m_cells->pvariables[PV->P][cellId]*=pRatio;
    }

    if(m_useSponge) {
      if(m_spongeLayerType == 2) {
        for(ZFSId cellId=0; cellId < m_noStrctrdCells; ++cellId){
          m_cells->fq[FQ->SPONGE_RHO][cellId]*=rhoRatio;
          m_cells->fq[FQ->SPONGE_RHO_E][cellId]*=rhoERatio;
        }
      }

      if(m_spongeLayerType == 4) {
        for(ZFSId cellId=0; cellId < m_noStrctrdCells; ++cellId){
          m_cells->fq[FQ->SPONGE_RHO][cellId]*=rhoRatio;        
        }
      }
    }
    
    computeConservativeVariables();

    zfs_log << "converting restart variables from old Ma: " << oldMa << " to new Ma: " << m_Ma << " ... SUCCESSFUL!" << endl;
  }
  return;
}

/// Compute right-hand side.
/// moved from zfsmethods
/// \author Sven Berger
/// \tparam nDim
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::rhs(const ZFSId rhsTimer) {
  TRACE();

  resetRHS();

  NEW_SUB_TIMER_STATIC(t_muscl, "MUSCL", rhsTimer);
  NEW_SUB_TIMER_STATIC(t_flux, "FLUX", rhsTimer);
  NEW_SUB_TIMER_STATIC(t_misc, "MISC", rhsTimer);

  RECORD_TIMER_START(t_muscl);
  Muscl(t_muscl);
  RECORD_TIMER_STOP(t_muscl);

  RECORD_TIMER_START(t_flux);
  if(!m_euler) {
    viscousFlux();
  }
  RECORD_TIMER_STOP(t_flux);

  RECORD_TIMER_START(t_misc);
  if( m_considerVolumeForces ) {
    computeVolumeForces();
  }
  RECORD_TIMER_STOP(t_misc);
}

/// Apply boundary condition
/// moved from zfsmethods
/// \author Sven Berger
/// \tparam nDim
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::rhsBnd()
{
  updateSpongeLayer();
}

template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::lhsBnd(const ZFSId lhsBndTimerId){
  NEW_SUB_TIMER_STATIC(t_computePV, "computePV", lhsBndTimerId);
  NEW_SUB_TIMER_STATIC(t_bnd, "bnd", lhsBndTimerId);
  NEW_SUB_TIMER_STATIC(t_exchange, "Exchange", lhsBndTimerId);

  RECORD_TIMER_START(t_computePV);
  computePrimitiveVariables();
  RECORD_TIMER_STOP(t_computePV);

  RECORD_TIMER_START(t_exchange);
  exchange();
  RECORD_TIMER_STOP(t_exchange);
  
  //junoh
  if(m_zonal) {
    computeCumulativeAverage(false);
    if(globalTimeStep%m_zonalExchangeInterval==0 && m_RKStep==0) {    
      // averagedFillGhostCells();
      // reconstructTurbulentVariables();
      spanwiseAvgZonal();
      zonalExchange();
      
      // averagedFillGhostCells();
      // reconstructTurbulentVariables();
    }
  }
  // saveOutputSolution<true>(true);
  // saveOutputPartitions<true>();

  RECORD_TIMER_START(t_bnd);
  applyBoundaryCondition();
  RECORD_TIMER_STOP(t_bnd);
}

/// Perform solution step
/// moved from zfsmethods
/// \author Sven Berger
/// \tparam nDim
template <ZFSInt nDim>
ZFSBool ZFSStrctrdBlck<nDim>::solutionStep()
{
  bool step = rungeKuttaStep();
  if( step ) setTimeStep();
  return step;
}


/// Initialize solver
/// moved from zfsmethods
/// \author Sven Berger
/// \tparam nDim
template <ZFSInt nDim>
void ZFSStrctrdBlck<nDim>::initSolver() {
  TRACE();
  initializeRungeKutta();
  
  initSolutionStep();
  
  setTimeStep();
}

// Explicit instantiations for 2D and 3D
template class ZFSStrctrdBlck<2>;
template class ZFSStrctrdBlck<3>;
