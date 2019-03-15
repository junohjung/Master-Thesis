#ifndef ZFSSTRCTBLCK_H
#define ZFSSTRCTBLCK_H

#include <sys/stat.h>

#include "zfsblock.h"
#include "zfslist.h"
#include "zfsstrctrdcell.h"
#include "zfsstrctrdbndrycnd.h"
#include "zfsstrctrdblckpartition.h"
#include "zfsstrctrdblckcommunicationhandle.h"
#include "zfsstrctrdblckwindowinfo.h"
#include "zfsmpi.h"
#include "zfsvariables.h"
#include "zfsstrctrdinterpolation.h"
#include "zfsstrctrdfqvariables.h"
#include "zfsstrctrdpostprocessing.h"
#include <cstring>
#include <functional>

template <ZFSInt nDim>
class ZFSStrctrdBlck : public ZFSBlock,
                       public ZFSStrctrdPostprocessing < nDim, ZFSStrctrdBlck<nDim> >
{
  template <ZFSInt nDim_> friend class ZFSStrctrdBndryCnd;
  template <ZFSInt nDim_> friend class ZFSStrctrdInterpolation;
  friend class ZFSStrctrdBlck3DRans;
  friend class ZFSStrctrdBlck2DRans;

public:
  // Add fields used from template base class to avoid calling with 'this->'
  using ZFSStrctrdPostprocessing<nDim, ZFSStrctrdBlck<nDim>>::m_averageStartTimestep;
  using ZFSStrctrdPostprocessing<nDim, ZFSStrctrdBlck<nDim>>::m_averageStopTimestep;
  using ZFSStrctrdPostprocessing<nDim, ZFSStrctrdBlck<nDim>>::m_averageInterval;
  using ZFSStrctrdPostprocessing<nDim, ZFSStrctrdBlck<nDim>>::m_noSamples;
  using ZFSStrctrdPostprocessing<nDim, ZFSStrctrdBlck<nDim>>::m_avgVariableNames;
  using ZFSStrctrdPostprocessing<nDim, ZFSStrctrdBlck<nDim>>::m_averageVorticity;
  using ZFSStrctrdPostprocessing<nDim, ZFSStrctrdBlck<nDim>>::m_kurtosis;
  using ZFSStrctrdPostprocessing<nDim, ZFSStrctrdBlck<nDim>>::m_skewness;
  using ZFSStrctrdPostprocessing<nDim, ZFSStrctrdBlck<nDim>>::m_tempWaveSample;
  using ZFSStrctrdPostprocessing<nDim, ZFSStrctrdBlck<nDim>>::m_postprocessing;

  ZFSStrctrdBlck ( ZFSId blockId, bool* propertiesGroups, const MPI_Comm comm);
  void         readGrid(ZFSInt file_id=-1);
  void         writeGrid();
  void         moveCellPoints(); //move the Cells to the right storage position in the cell Array (because of ghost cells)
  void         writeGridPointsWithGhostPoints();
  void         initializeFQField();

  virtual void writeHeaderAttributes(ZFSId fileId, ZFSString gridFile, ZFSString fileType);
  virtual void writePropertiesAsAttributes(ZFSId fileId, ZFSString path);
  virtual void saveOutput(ZFSBool forceOutput);
  std::function<void(ZFSBool)> saveSolution;
  std::function<void()> savePartitions;
  std::function<void()> savePlanes;
  std::function<void()>saveBoxes;
  template<ZFSBool primitiveOutput> void saveOutputSolution(ZFSBool forceOutput); //junoh
  template<ZFSBool primitiveOutput> void saveOutputPartitions();
  template<ZFSBool primitiveOutput> void saveOutputPlanes();
  template<ZFSBool primitiveOutput>  void saveOutputBoxes();
  virtual void saveOutputLines() {};
  virtual void writeRestartFile(ZFSBool forceOutput);
  //needed for restart
  template <ZFSBool isPrimitive>
  void shiftCellValuesRestart();

  void loadRestartFile();
  virtual void loadRestartBC2600(){};
  virtual void loadRestartBC2601(){};
  virtual void loadRestartSTG(ZFSBool){zfsTerm(-1, __CALLING_FUNCTION__, "not implemented in basic block");};
  //end
  //sandPaper
  virtual void tripForceCoefficients(ZFSFloat*, ZFSFloat*,  ZFSFloat*, ZFSInt, ZFSInt, ZFSFloat){}
  virtual void saveLiftDragToAsciiFile(){};
  void saveAveragedVariables(ZFSString, ZFSId, ZFSFloat**);
  void saveProductionTerms(ZFSString, ZFSFloat**);
  void saveAverageRestart(ZFSString, ZFSId, ZFSFloat**, ZFSFloat**, ZFSFloat**, ZFSFloat**);




  void saveAuxData();
  virtual void saveLiftCoefficient(ZFSInt){};
  virtual void saveDragCoefficient(ZFSInt){};
  virtual void savePowerCoefficient(ZFSInt){};
  
  virtual void computeAuxData(){};
  virtual void computeAuxDataRoot(){};
  virtual      ~ZFSStrctrdBlck();

  //!Structured Block Constructor reads and allocate properties/variables
  void         initializeStrctrdBlck(ZFSBool* propertiesGroups);
  void         allocateAndInitBlockMemory();
  void         setRungeKuttaProperties();
  void         setSamplingProperties();
  void         setNumericalProperties();
  void         setInputOutputProperties();
  void         setZonalProperties();
  void         allocateVariables();
  void         setTestcaseProperties();
  void         setMovingGridProperties();
  void         readAndSetSpongeLayerProperties();
  void         setSTGProperties();
  void         setProfileBCProperties();
  void         createMPIGroups();
  void         prepareNonBlockingCommunication();

  //methods:
  void         computePV();
  void         partitionGrid();
  virtual void computePrimitiveVariables(){zfsTerm(-1,__CALLING_FUNCTION__,"not implemented in general block");};
  virtual void computeConservativeVariables();//{zfsTerm(-1,__CALLING_FUNCTION__,"not implemented in general block");};
  void         saveVarToPrimitive(ZFSId, ZFSId, ZFSFloat);
  virtual void ppFillGhostCells() {};
  /* virtual void ppAveragedFillGhostCells(); //junoh */
  void         computeSamplingInterval();
  void         checkNans();
  virtual void computeVolumeForces(){};
  void         allocateMetrics();
  void         allocateSurfaceMetrics();
  void         allocateCornerMetrics();
  void         allocateCellMetrics();
  void         allocateJacobian();
  void         computeMetrics();
  void         computeJacobian();
  virtual void computeSurfaceMetrics(){};
  virtual void computeModSurfaceMetrics(){};
  virtual void computeCornerMetrics(){};
  virtual void computeModCornerMetrics(){};
  virtual void computeCellMetrics(){};
  virtual void computeCornerJacobian(){};
  virtual void computeCellJacobian(){};
  virtual void computeModCornerJacobian(){};
  virtual void computeModCellJacobian(){};
  virtual void applyBoundaryCondtition(){};
  virtual void moveGrid(ZFSInt ){};
  virtual void initMovingGrid(){};
  virtual void moveGrid(ZFSBool, ZFSBool){};
  virtual void saveGrid(){};
  virtual void computeLambda2Criterion() {};
  virtual void computeVorticity() {};
  virtual void exchange();  
  /* virtual void exchange(); */
  /* virtual void zonalExchange(){}; //junoh  */
  virtual void zonalExchange(){}; //junoh 
  /* virtual void reconstructTurbulentVariables(){}; //junoh */
  virtual void spanwiseAvgZonal(){}; //junoh

  virtual void waveExchange(){};
  virtual void spanwiseWaveReorder(){};
  virtual void exchangePoints();
  virtual void scatterPoints(){};
  virtual void gatherPoints(){};
  virtual void sendPoints(){};
  virtual void receivePoints(){};
  void         setTimeStep();
  void         fixTimeStepTravelingWave();
  void         exchangeTimeStep();
  void         initializeRungeKutta();
  virtual void computeTimeStep(){};
  ZFSBool      isInInterval(ZFSInt);
  ZFSInt       getNoCells() {return m_noStrctrdCells;};
  ZFSId        noVariables() const override { return PV->noVariables; };
  ZFSInt       getNoActiveCells() {return m_noStrctrdActiveCells;};
  ZFSInt*      getActiveCells() {return m_nActiveCells;};
  ZFSInt       getNoGhostLayers() {return m_noGhostLayers;};
  ZFSInt       getWaveAvrgInterval() {return (m_waveNoStepsPerCell);};
  ZFSInt       getWaveStepOffset() {return (m_movingGridStepOffset);};
  ZFSInt*      getCellGrid() {return m_nCells;};
  ZFSBool      isMovingGrid() {return m_movingGrid;};
  ZFSInt       getGridMovingMethod() {return m_gridMovingMethod;};
  ZFSId        getBoxId(ZFSId id) {return  m_partition->outputBoxInfo[id]->cpu;};
  virtual void computeCumulativeAverage(ZFSBool){};   //junoh
  virtual void averagedFillGhostCells(){};//junoh
 


  virtual void loadSampleFile(ZFSString) {};
  virtual void getSampleVariables(ZFSId, ZFSFloat*) {};
  virtual ZFSFloat getSampleVorticity(ZFSId, ZFSId) {return 0;};
  virtual ZFSFloat dvardxyz(ZFSId, ZFSId, ZFSFloat*) {return 0;};
  virtual ZFSFloat dvardx(ZFSId, ZFSFloat*) {return 0;};
  virtual void loadAverageRestartFile(const char*,  ZFSFloat**, ZFSFloat**, ZFSFloat**,ZFSFloat**) {}
  virtual void loadAveragedVariables(const char*) {};
  void         convertRestartVariables(ZFSFloat oldMa);
  virtual  void convertRestartVariablesSTG(ZFSFloat oldMa){(void)oldMa; zfsTerm(-1, __CALLING_FUNCTION__, "not implemented for 0d,2d");};
  

  void         tred2(ZFSFloatScratchSpace& A, ZFSInt dim, ZFSFloat* diag, ZFSFloat* offdiag );
  void         tqli2(ZFSFloat* diag, ZFSFloat* offdiag, ZFSInt dim);
  void         insertSort(ZFSId dim, ZFSFloat* list);
  ZFSFloat     pythag(ZFSFloat a, ZFSFloat b);

  ZFSFloat sign(ZFSFloat x, ZFSFloat y);
  ZFSFloat bessJ (ZFSInt n, ZFSFloat x);
  ZFSFloat bessJ0 (ZFSFloat x);
  ZFSFloat bessJ1 (ZFSFloat x);
  ZFSFloat bessI(ZFSInt n, ZFSFloat x);
  ZFSFloat bessI0(ZFSFloat x);
  ZFSFloat bessI1(ZFSFloat x);

  void resetRHS();
  void rhs(const ZFSId rhsTimer) override;
  void rhsBnd() override;
  void lhsBnd(const ZFSId lhsBndTimerId) override;
  void initSolver() override;
  ZFSBool solutionStep() override;
  ZFSFloat time() const override {return m_time;}

  //if zonal then this communicator is the subworld of each block
  MPI_Comm     m_zfsStrctrdComm;

  ZFSInt       m_restartTimeStep;
  ZFSString    m_outputFormat;
  ZFSInt       m_lastOutputTimeStep;

protected:
  // epsilon
  ZFSFloat m_eps;

  // left/right States
  ZFSFloat* m_QLeft;
  ZFSFloat* m_QRight;

  //moving grids:
  ZFSId        m_movingGrid;
  ZFSInt       m_gridMovingMethod;
  ZFSInt       m_movingGridStepOffset;
  ZFSInt       m_synchronizedMGOutput;
  ZFSInt       m_waveNoStepsPerCell;
  ZFSFloat     m_wallVel;
  ZFSFloat     m_movingGridTimeOffset;
  ZFSBool      m_movingGridSaveGrid;
  ZFSFloat*    m_makosMaxAmplitude;
  ZFSFloat     m_makosActuatorRadius;
  ZFSFloat     m_makosPeakAmplitude;
  ZFSFloat     m_makosCenterX;
  ZFSFloat     m_makosCenterZ;
  ZFSBool      m_travelingWave;
  ZFSFloat     m_waveLengthPlus;
  ZFSFloat     m_waveAmplitudePlus;
  ZFSFloat     m_waveTimePlus;
  ZFSFloat     m_waveTime;
  ZFSFloat     m_waveLength;
  ZFSFloat     m_waveAmplitude;
  ZFSFloat     m_waveSpeed;
  ZFSFloat     m_waveSpeedPlus;
  ZFSFloat     m_waveBeginTransition;
  ZFSFloat     m_waveEndTransition;
  ZFSFloat     m_waveOutBeginTransition;
  ZFSFloat     m_waveOutEndTransition;
  ZFSBool      m_movingGridInitialStart;
  ZFSBool      m_waveRestartFadeIn;
  ZFSInt       m_waveTimeStepComputed;
  ZFSInt       m_waveCellsPerWaveLength;

  ZFSFloat*    m_rhs;

  //for IO
  ZFSInt       m_noInputBlocks;
  ZFSBool      m_useNonSpecifiedRestartFile;
  ZFSId        m_outputOffset;
  ZFSInt       m_vorticityOutput;
  ZFSBool      m_debugOutput;
  ZFSBool      m_writeSlopes;
  ZFSInt       m_outputIterationNumber;
  ZFSBool      m_sampleSolutionFiles;

  ZFSInt       m_nonBlockingComm;

  //structure for cell data
  ZFSStrctrdCell*    m_cells;
  ZFSString*         m_variableNames;
  ZFSString*         m_pvariableNames;

  ZFSFloat**   m_coordinates; //(nogridpoints, dimension)
  ZFSFloat**   m_mgInitCoordinates;//(contains the initial values of the grid if moving grid is used)
  ZFSFloat**   m_mgOrgCoordinates; //(contains the original values of the grid if moving grid is used)
  ZFSFloat**   m_mgOldCoordinates; //(contains the values at the old time step if moving grid is used)

  //for debugging only
  ZFSFloat**   pointProperties;
#ifdef ZFS_EXTRA_DEBUG
  ZFSFloat**   viscFluxOut;
  ZFSFloat**   convFluxOut;
#endif
  //end debugging

  ZFSString    m_uID; //unique identifier for the restart and gridfile
  ZFSInt       m_ignoreUID;

  ZFSInt       m_noGridPoints; //gridpoints in the partition
  ZFSStrctrdDecomposition<nDim>* m_partition; //contains info about partitioning
  ZFSStrctrdBlckWindowInfo<nDim>* m_windowInfo; //contains info about the window information
  ZFSInt       m_noGhostLayers; //number of GhostLayerst to be added
  ZFSInt*      m_nPoints; //stores the maximum dimension of the partition with ghost points
  ZFSInt*      m_nActivePoints; // stores the  maximum dimension of the partition without ghost points
  ZFSInt       m_noStrctrdCells; //stores number of structured cells
  ZFSInt       m_noStrctrdActiveCells;

  ZFSInt       m_noSurfaces;
  ZFSInt       m_noActiveSurfaces;
  ZFSFloat     m_referenceLength;
  ZFSFloat     m_physicalReferenceLength;
  ZFSFloat     m_Pr;
  ZFSFloat     m_rPr;
  ZFSFloat     m_cfl;
  ZFSId        m_orderOfReconstruction;
  ZFSFloat     m_inflowTemperatureRatio;
  ZFSBool      m_considerVolumeForces;
  ZFSBool      m_euler;
  ZFSFloat*    m_volumeForce;
  ZFSFloat     m_gamma;
  ZFSFloat     m_gammaMinusOne;
  ZFSFloat     m_fgammaMinusOne;
  ZFSFloat     m_Re0;
  ZFSFloat     m_ReTau;
  ZFSFloat*    m_angle;
  ZFSId        m_periodicConnection;
  ZFSFloat*    m_periodicDisplacements;
  ZFSId        m_channelFlow;
  ZFSFloat     m_sutherlandConstant;
  ZFSFloat     m_sutherlandPlusOne;
  ZFSFloat     m_TinfS;
  ZFSInt       m_computeLambda2;


  //
  ZFSFloat     m_hInfinity;
  ZFSFloat     m_referenceEnthalpy;

  //variables
  ZFSConservativeVariables<nDim>* CV;
  ZFSPrimitiveVariables<nDim>* PV;
  ZFSStrctrdFQVariables* FQ;
  ZFSInt       m_maxNoVariables; //junoh

  //ZFSInt      m_bCf;
  //ZFSInt      m_bCp;
  ZFSInt      m_bCfCpCoeff;
  ZFSInt      m_bPower;
  ZFSInt      m_bCl;
  ZFSInt      m_bCd;
  ZFSInt      m_detailAuxData;
  ZFSInt      m_bCpLineAveraging;
  ZFSInt      m_cpAveragingDir;
  ZFSFloat**  m_cl;
  ZFSFloat**  m_cd;
  ZFSInt      m_dragOutputInterval;
  ZFSInt      m_dragAsciiOutputInterval;
  ZFSFloat    m_globalDomainWidth;
  ZFSInt      m_auxDataCoordinateLimits;
  ZFSFloat*   m_auxDataLimits;

  ZFSInt      m_primitiveOutput;
  ZFSInt      m_noPlaneOutput;
  ZFSInt      m_planeOutputInterval;
  ZFSInt*     m_planeBlock;
  ZFSInt*     m_planeOffset;
  ZFSInt*     m_planeNormal;
  ZFSBool     m_planeWriteCoordinates;

  ZFSInt      m_lineOutputInterval;
  ZFSInt      m_noLineOutput;
  ZFSInt      m_no2dLines;
  ZFSInt      m_noFieldPointsTotal;
  ZFSInt*     m_lineNoPoints;
  ZFSInt*     m_lineNoPoints2d;
  ZFSInt*     m_fieldOffset;
  ZFSFloat**  m_lineStart;
  ZFSFloat**  m_lineDelta;
  ZFSFloat**  m_lineDelta2d;
  ZFSFloat**  m_pointCoordinates;
  ZFSInt*     m_hasPartnerGlobal;
  ZFSInt*     m_hasPartnerLocal;
  ZFSFloat**  m_interpolatedVarsGlobal;
  ZFSFloat**  m_interpolatedVarsLocal;
  ZFSBool     m_fieldInterpolation;

  ZFSInt      m_noBoxOutput;
  ZFSInt      m_boxOutputInterval;
  ZFSInt*     m_boxBlock;
  ZFSInt**    m_boxOffset;
  ZFSInt**    m_boxSize;
  ZFSBool     m_boxWriteCoordinates;

  //
  ZFSFloat     m_referenceTemperature;
  ZFSInt       m_noSpecies;
  ZFSFloat*    m_referenceComposition;
  ZFSFloat*    m_formationEnthalpy;
  ZFSId        m_noRKSteps;
  ZFSFloat*    m_RKalpha;
  ZFSFloat     m_time;
  ZFSId        m_RKStep;
  ZFSFloat     m_timeStep;
  ZFSBool      m_localTimeStep;
  ZFSId        m_rungeKuttaOrder;
  ZFSId        m_timeStepMethod;
  ZFSInt       m_timeStepComputationInterval;
  ZFSFloat     m_timeRef;
  ZFSId        m_dualTimeStepping;
  ZFSFloat     m_physicalTimeStep;
  ZFSFloat     m_physicalTime;


  //constant Timestepping
  ZFSInt       m_constantTimeStep;

  //for averaging
  ZFSFloat     m_startAvgTime;
  ZFSFloat     m_avgFactor;

  //rescaling
  ZFSInt*      m_rescalingCommGrRoot;
  ZFSInt*      m_rescalingCommGrRootGlobal;
  MPI_Comm*    m_rescalingCommGrComm;

  //synthetic turbulence generation
  //for synthetic turbulence generation
  MPI_Comm*    m_commStg;
  ZFSInt*      m_commStgRoot;
  ZFSInt*      m_commStgRootGlobal;
  ZFSInt       m_commStgMyRank;
  ZFSBool      m_stgIsActive;
  ZFSFloat     m_stgBLT1;
  ZFSFloat     m_stgBLT2;
  ZFSFloat     m_stgBLT3;
  ZFSFloat     m_stgDelta99Inflow;
  ZFSInt       m_stgInitialStartup;
  ZFSId        m_stgNoEddieProperties;
  ZFSFloat**   m_stgEddies;
  ZFSInt       m_stgNoEddies;
  ZFSInt       m_stgMaxNoEddies;
  ZFSFloat     m_stgExple;
  ZFSFloat     m_stgEddieDistribution;
//  ZFSInt       m_numcellsBC;
  ZFSInt       m_stgBoxSize[3];
//  ZFSBool      m_restarted;
  ZFSBool      m_stgLocal;
  ZFSBool      m_stgCreateNewEddies;
  ZFSBool      m_stgRootRank;
  ZFSInt       m_stgNoVariables;
  ZFSInt       m_stgShapeFunction;
  ZFSBool      m_stgEddieLengthScales;
  ZFSBool      m_stgSubSup;
  ZFSBool      m_stgSupersonic;
  ZFSInt       m_stgFace;
  ZFSFloat*    m_stgLengthFactors;
  ZFSInt       m_stgMyRank;

  //global Upwind Coefficient;
  ZFSFloat     m_chi;

  //Sponge Variables
  ZFSInt       m_noSpongeDomainInfos;
  ZFSId*       m_spongeBcWindowInfo;
  ZFSInt       m_useSponge;
  ZFSInt       m_spongeLayerType;
  ZFSFloat*    m_sigmaSponge;
  ZFSFloat*    m_betaSponge;
  ZFSFloat*    m_spongeLayerThickness;
  ZFSFloat     m_targetDensityFactor;
  ZFSBool      m_computeSpongeFactor;

  //repartitionFile
  ZFSInt       m_readDecompositionFromFile;

  //Workload
  ZFSFloat     m_workload;
  ZFSFloat     m_workloadIncrement;
  ZFSString    m_gridInputFileName;
  //to store the absolute path to the grid file
  ZFSString    m_gridFileLocation;

  //Residual
  //MPI_derived datatypes
  typedef struct {
    ZFSFloat   maxRes;
    ZFSFloat   avrgRes;
    ZFSInt*    maxCellIndex;
  } ZFSRes;

  ZFSRes*      m_residualSnd;
  ZFSRes       m_residualRcv;

  ZFSInt       m_residualOutputInterval;

  ZFSFloat     m_avrgResidual;
  ZFSFloat     m_firstMaxResidual;
  ZFSFloat     m_firstAvrgResidual;
  ZFSLong      m_totalGridCells;
  ZFSInt**     m_totalGridBlockDim;
  ZFSLong*     m_totalGridBlockCells;
  ZFSInt       m_inputBlockId;
  ZFSBool      m_residualFileExist;
  MPI_Op       m_resOp;
  MPI_Datatype m_mpiStruct;
  FILE*        m_resFile;

  //Convergence
  ZFSBool      m_convergence;
  ZFSFloat     m_convergenceCriterion;

  //Restart & initial Condition
  ZFSId        m_restart;
  ZFSId        m_initialCondition;
  ZFSFloat     m_deltaP;
  ZFSId        m_restartInterpolation;

  //Strings for Output specification
  ZFSString    m_solutionOutput;
  ZFSString    m_solutionFileName;
  ZFSString    m_planeOutputDir;
  ZFSString    m_boxOutputDir;
  ZFSString    m_lineOutputDir;
  ZFSString    m_auxOutputDir;
  ZFSBool      m_savePartitionOutput;

  ZFSInt*      m_nCells; //cell array dimension with ghost layer
  ZFSInt*      m_nActiveCells; // cell array dimension without ghost layers
  ZFSInt*      m_nOffsetCells;
  ZFSInt*      m_nInputBlockCells;

  //For Boundary Conditions
  ZFSInt*      m_noWindows;      //contains the number of windows of each face
  ZFSInt**     m_bndryCndInfo;  //contains start and length of each window
  ZFSInt**     m_bndryCnd;      //contains the boundary condition

  //limiter
  ZFSId        m_limiter; //for limiter choice
  ZFSString    m_limiterMethod; //reads the limiter from the properties
  ZFSFloat     m_venkFactor; //customizable factor for the modified Venkatakrishnan Limiter

  ZFSString    m_musclScheme;
  ZFSString    m_ausmScheme;

  //zonal
  ZFSInt       m_zonal;
  ZFSInt       m_rans;
  ZFSInt       m_noRansEquations;
  ZFSString    m_zoneType;
  ZFSBool      m_hasSTG;  //junoh
  ZFSInt       m_zonalExchangeInterval;   //junoh

  MPI_Comm*    m_commBC2600; //junoh
  ZFSInt*      m_commBC2600Root;
  ZFSInt*      m_commBC2600RootGlobal;
  ZFSInt       m_commBC2600MyRank;
  /* ZFSBool      m_bc2600RootRank; */
  
  MPI_Comm*    m_commZonal; //junoh
  ZFSInt*      m_commZonalRoot;
  ZFSInt*      m_commZonalRootGlobal;
  ZFSInt       m_commZonalMyRank;
  ZFSBool      m_zonalRootRank; 
  

  //zonal cumulative averaging
  ZFSBool      m_zonalExponentialAveraging;
  ZFSInt       m_zonalStartAvgTime;
  ZFSFloat     m_zonalAveragingFactor;


  ZFSFloat     m_rhoNuTildeInfinty;
  ZFSFloat     m_nutInfinity;
  ZFSFloat     m_mutInfinity;

  //parallel
  ZFSStrctrdCommunicationHandle* m_cmnctnFlag;
  ZFSStrctrdCommunicationHandle* m_ppCmnctnFlag;
  ZFSStrctrdCommunicationHandle* m_avCmnctnFlag; //junoh
  ZFSStrctrdWaveCommunicationHandle* m_waveCmnctnFlag;
  ZFSInt       m_noNghbrDomains;
  ZFSInt*      m_nghbrDomainId;
  ZFSInt*      m_noNghbrDomainBufferSize;
  ZFSFloat**   m_bufferCellsSndRcv; //comunicator cells==>cells to be communicated
  ZFSFloat**   m_bufferPointsSendRcv; //comunicator points ==> points to be communicated
  ZFSInt*      m_nghbrFaceId;
  ZFSInt**     m_nghbrFaceInfo;
  MPI_Request* mpi_sndRequest;
  MPI_Request* mpi_rcvRequest;
  MPI_Status*  mpi_sndRcvStatus;

  //channel boundary condition
  MPI_Comm*    m_commChannelIn;
  MPI_Comm*    m_commChannelOut;
  MPI_Comm*    m_commChannelWorld;
  ZFSInt*      m_channelRoots;
  ZFSFloat     m_channelPresInlet;
  ZFSFloat     m_channelPresOutlet;
  ZFSFloat     m_channelLength;
  ZFSFloat     m_channelHeight;
  ZFSFloat     m_channelWidth;
  ZFSFloat     m_channelInflowPlaneCoordinate;
  ZFSFloat     m_channelC1;
  ZFSFloat     m_channelC2;
  ZFSFloat     m_channelC3;
  ZFSFloat     m_channelC4;

  //for the periodic rotation Boundary
  MPI_Comm*    m_commPerRotOne;
  MPI_Comm*    m_commPerRotTwo;
  MPI_Comm*    m_commPerRotWorld;
  ZFSInt*      m_commPerRotRoots;
  ZFSInt       m_commPerRotGroup;

  //profile BC 2600
  ZFSBool      m_bc2600IsActive;
  ZFSInt       m_bc2600InitialStartup;
  ZFSFloat**   m_bc2600Variables;
  ZFSBool      m_bc2600;
  ZFSInt*      m_bc2600noCells;
  ZFSInt*      m_bc2600noActiveCells;
  ZFSInt*      m_bc2600noOffsetCells;
  ZFSInt       m_bc2600RootRank;

  //effective boundary condition 2601
  ZFSBool      m_bc2601IsActive;
  ZFSInt       m_bc2601InitialStartup;
  ZFSFloat**   m_bc2601Variables;
  ZFSFloat**   m_bc2601ZerothOrderSolution;
  ZFSBool      m_bc2601;
  ZFSInt*      m_bc2601noCells;
  ZFSInt*      m_bc2601noActiveCells;
  ZFSInt*      m_bc2601noOffsetCells;
  ZFSFloat     m_bc2601GammaEpsilon;

  //convective solution output
  ZFSBool      m_useConvectiveUnitWrite;
  ZFSFloat     m_convectiveUnitInterval;
  ZFSInt       m_noConvectiveOutputs;

  //timer
  void         initializeTimers();
  ZFSId        m_tgStr;
  ZFSId        m_tcomm;
  ZFSId        m_texchange;
  ZFSId        m_tgatherAndSend;
  ZFSId        m_tscatterAndReceive;
  ZFSId        m_tgatherAndSendWait;
  ZFSId        m_tscatterWaitSome;
  ZFSId        m_tsendWait;
  ZFSId        m_tgather;
  ZFSId        m_tsend;
  ZFSId        m_treceive;
  ZFSId        m_tscatter;
  ZFSId        m_treceiving;
  ZFSId        m_treceiveWait;
  ZFSId        m_texchangeDt;
  ZFSId        m_tdecomposition;
  ZFSId        m_treadGrid;
  ZFSId        m_tloadRestart;
  ZFSId        m_tloadRestartVars;
  ZFSId        m_tloadRestartSponge;
  ZFSId        m_tloadRestartStg;
  ZFSId        m_tloadRestartStgEddies;
  ZFSId        m_tloadRestartStgRead;
  ZFSId        m_tloadRestartStgBcast;
  ZFSId        m_tloadRestartShift;
  ZFSId        m_tbuildSponge;
  ZFSId        m_tspanwiseReorder;
  ZFSId        m_tcomputeMetrics;
  ZFSId        m_tcomputeJacobian;
  ZFSId        m_tConvectiveFlux;
  ZFSId        m_tViscousFlux;
  ZFSId        m_tVolumeFlux;


  //infinity Values
  ZFSFloat      m_UInfinity;
  ZFSFloat      m_VInfinity;
  ZFSFloat      m_WInfinity;
  ZFSFloat      m_PInfinity;
  ZFSFloat      m_TInfinity;
  ZFSFloat      m_DthInfinity;
  ZFSFloat      m_muInfinity;
  ZFSFloat      m_DInfinity;
  ZFSFloat      m_VVInfinity[3];//3 is the max value!!!!
  ZFSFloat      m_rhoUInfinity;
  ZFSFloat      m_rhoVInfinity;
  ZFSFloat      m_rhoWInfinity;
  ZFSFloat      m_rhoEInfinity;
  ZFSFloat      m_rhoInfinity;
  ZFSFloat      m_rhoVVInfinity[3];

  ZFSInt        m_changeMa;

  //for modes:
  ZFSId         m_restartBc2800;
  ZFSFloat      m_restartTimeBc2800;
  ZFSFloat*     m_adiabaticTemperature;
  ZFSInt        m_useAdiabaticRestartTemperature;

  //singularity
  SingularInformation *m_singularity;
  ZFSInt        m_hasSingularity;

  //sandpaper tripping
  ZFSBool   m_useSandpaperTrip;
  ZFSFloat  m_tripXOrigin;
  ZFSFloat  m_tripXLength;
  ZFSFloat  m_tripYOrigin;
  ZFSFloat  m_tripYHeight;
  ZFSFloat  m_tripMaxAmpSteady;
  ZFSFloat  m_tripMaxAmpFluc;
  ZFSInt    m_tripNoModes;
  ZFSFloat  m_tripDeltaTime;
  ZFSId     m_tripTimeStep;
  ZFSInt    m_tripSeed;
  ZFSFloat* m_tripG;
  ZFSFloat* m_tripH1;
  ZFSFloat* m_tripH2;
  ZFSFloat* m_tripModesG;
  ZFSFloat* m_tripModesH1;
  ZFSFloat* m_tripModesH2;
  ZFSFloat* m_tripCoords;
  ZFSInt    m_tripNoCells;
  ZFSFloat  m_tripDomainWidth;

private:

  /// ZFSPostProcessingBlock interface:
  virtual void initStrctrdPostprocessing() { ZFSStrctrdPostprocessing< nDim, ZFSStrctrdBlck<nDim> >::initStrctrdPostprocessing(); }
  virtual void saveAverageRestart()  { ZFSStrctrdPostprocessing < nDim, ZFSStrctrdBlck<nDim> >::saveAverageRestart();   }
};

class ZFSStrctrdZonalBC   //junoh
{
 public:
  ZFSStrctrdZonalBC() {m_noCellsGlobalBC=0; m_noCellsLocalBC=0; m_noGlobalSndDomains=0; m_noGlobalRcvDomains=0; m_noSndNghbrDomains=0; m_noRcvNghbrDomains=0; m_noZonalVariables=0;  m_interpolatedVars=NULL; m_interpolatedVarsAV=NULL; m_globalReceiverIds=NULL; /* m_globalLocalCellIds=NULL; */ m_globalLocalMapCellIds=NULL; /* m_hasPartnerLocalBC=NULL; m_hasPartnerGlobalBC=NULL; */ m_hasLocalBCMap=false; m_hasSTG=false; m_hasZonalwithoutSTG=false; m_bufferSndZonal=NULL; m_bufferRcvZonal=NULL; m_bufferSndMapCellId=NULL; m_bufferRcvMapCellId=NULL;  m_globalRcvZonalId=NULL; m_globalSndZonalId=NULL; /* m_globalBufferAverage=NULL; */ m_localBufferSndSize=NULL; m_localBufferRcvSize=NULL; /* m_localBufferSndCellId=NULL; */ m_localRcvId=NULL; m_localSndId=NULL; /* m_localHasPartnerGlobalBC=NULL; */ m_localMapCellsId=NULL; m_zonalBCCells=NULL; m_localBufferIndexCellId=NULL; m_localBufferMapCellId=NULL;};
 ~ZFSStrctrdZonalBC();


 ZFSInt    m_noCellsGlobalBC;
 ZFSInt    m_noCellsLocalBC;
 /* ZFSBool   m_hasZonalCells; */
 /* ZFSInt    m_noDonorDomainsGlobal; */
 /* ZFSInt    m_noRcvDomainsGlobal; */
 /* ZFSId *   m_donorDomainIdGlobal;    //DomainIds which will send the data */
 /* ZFSId *   m_rcvDomainIdGlobal; */      //DomainIds which will receive the data
 ZFSBool   m_hasLocalBCMap;
 ZFSBool   m_hasSTG;
 ZFSBool   m_hasZonalwithoutSTG;

 ZFSInt*   m_globalReceiverIds;
 /* ZFSInt*   m_globalLocalCellIds; */
 ZFSInt*   m_globalLocalMapCellIds;
 ZFSInt*   m_localMapCellsId; 

 ZFSInt*   m_hasPartnerLocalBC;
 ZFSInt*   m_hasPartnerGlobalBC;
 /* ZFSInt*  m_nOffsetCells; */
 ZFSFloat** m_interpolatedVars;
 ZFSFloat* m_interpolatedVarsAV;
 /* ZFSFloat* m_globalBufferAverage; */

 //communicator
 ZFSFloat**   m_bufferSndZonal; //comunicator cells==>cells to be communicated (send and received)
 ZFSFloat**   m_bufferRcvZonal;
 ZFSFloat**   m_bufferSndMapCellId;
 ZFSFloat**   m_bufferRcvMapCellId;
 ZFSInt*      m_localBufferSndSize;
 ZFSInt*      m_localBufferRcvSize;
 ZFSInt*      m_localRcvId;
 ZFSInt*      m_localSndId;
 ZFSInt    m_noGlobalRcvDomains;
 ZFSInt    m_noGlobalSndDomains;
 ZFSInt*   m_globalRcvZonalId;
 ZFSInt*   m_globalSndZonalId;
 ZFSInt    m_noSndNghbrDomains;
 ZFSInt    m_noRcvNghbrDomains;
 ZFSInt*   m_localBufferMapCellId;
 ZFSInt*   m_localBufferIndexCellId;
 ZFSInt    m_noBufferSndSize;
 ZFSInt*   m_localCommReceiverIds;
 /* ZFSInt*   m_localHasPartnerGlobalBC; */

 ZFSInt*   m_zonalBCCells;
 ZFSId     zonalCellIndex(ZFSId i, ZFSId j, ZFSId k) {return i + (j+k*m_zonalBCCells[1])*m_zonalBCCells[2];};
 ZFSInt m_noZonalVariables;

 
 MPI_Request* mpi_sndRequest;
 MPI_Request* mpi_rcvRequest;
 MPI_Status* mpi_sndStatus;
 MPI_Status* mpi_rcvStatus;

};



#endif


