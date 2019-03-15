#ifndef ZFSSTRCTRDBLCK3D_H
#define ZFSSTRCTRDBLCK3D_H

#include "zfsstrctrdblck.h"
#include "zfsstrctrdbndrycnd3d.h"
#include "zfsstrctrdblck3drans.h"
//fftw is needed for the channel flow
#include <fftw3.h>
#include <complex>

extern ZFSId globalTimeStep;

enum FLUXMETHOD{
  AUSM_LES,
  AUSM_LES_PTHRC,
  AUSMDV,
};

class ZFSStrctrdBlck3D : public ZFSStrctrdBlck<3>
{
  template<ZFSBool isRans> friend class ZFSStrctrdBndryCnd3D;
  friend class ZFSStrctrdBlck3DRans;
public:
  ZFSStrctrdBlck3D( ZFSId, bool *, const MPI_Comm comm);
  ~ZFSStrctrdBlck3D();

  void         ZFSResOperationFunction(ZFSRes *in, ZFSRes *out, int* len);
  void         addGhostPointCoordinateValues();
  void         extrapolateGhostPointCoordinates();
  void         extrapolateGhostPointCoordinatesBC();

  void         copyPointIntoOrder(ZFSInt nghbr);
  void         initSolutionStep();
  void         initialCondition() override;
  void         computeCellCentreCoordinates();
  void         assignBndryCells();
  void         nonReflectingBC();
  void         applyBoundaryCondition() override;
  void         calcSurfaceMetrics();
  void         initBndryCnds();
  void         allocateAuxDataMaps();
  void         allocateSingularities();
  void         computeVolumeForces() override;
  void         applySandpaperTrip();
  void         initSandpaperTrip();
  void         tripForceCoefficients(ZFSFloat*, ZFSFloat*, ZFSFloat*, ZFSInt, ZFSInt, ZFSFloat);
  void         tripFourierCoefficients(ZFSFloat*, ZFSInt);

  void         computeSurfaceMetrics() override;
  void         computeModSurfaceMetrics() override;
  void         computeCornerMetrics() override;
  void         computeModCornerMetrics() override;
  void         computeCellMetrics() override;
  void         computeCornerJacobian() override;
  void         computeCellJacobian() override;
  void         computeModCornerJacobian() override;
  void         computeModCellJacobian() override;

   //computeZonalConnections   //junoh
  void         computeZonalConnections();
  void         zonalExchange() override; 
  void         zonalAllreduce();
  void         zonalRealInterpolation();
  void         zonalGather();
  void         zonalSend();
  void         zonalReceive();
  void         zonalCheck();
  void         zonalBufferAverage();
  void         zonalScatter();
  /* void         zonalSTG(); */
  ZFSStrctrdZonalBC* m_zonalBC;
  ZFSInt m_noZonalBCMaps;
  /* void         reconstructTurbulentVariables() override; //junoh */
  void         spanwiseAvgZonal() override; //junoh
  
  

  /* ZFSStrctrdZonalCommunicationHandle* m_zonalCmnctnFlag; */
  
 //

  void         distributeFluxToCells();
  void         computeTimeStep() override;
  void         writeGridCellsWithGhostCells();
  void         computeVorticity() override;
  void         computeLambda2Criterion() override;
  ZFSFloat     computeTotalKineticEngergy();
  ZFSFloat     computeTotalPressure();
  void         initializeNeighbourCellIds();
  void         initFluxMethod();
  void         computePrimitiveVariables() override;
  ZFSFloat     pressure(ZFSId);
  inline ZFSFloat getPSI(ZFSId, ZFSId);

  void         Ausm();
  inline void  AusmLES(ZFSFloat* QLeft, ZFSFloat* QRight, const ZFSId dim, const ZFSId cellId);
  inline void  AusmLES_PTHRC(ZFSFloat* QLeft, ZFSFloat* QRight, ZFSId dim, ZFSId cellId);
  inline void  AusmDV(ZFSFloat* QLeft, ZFSFloat* QRight, ZFSId dim, ZFSId cellId);

  typedef void (ZFSStrctrdBlck3D::*fluxmethod)(ZFSFloat*, ZFSFloat*, ZFSId, ZFSId);
  template< ZFSId noVars> void Muscl_AusmLES();
  template< ZFSId noVars> void Muscl_AusmLES_PTHRC();
  template< ZFSId noVars> void Muscl_AusmDV();
  template< fluxmethod ausm, ZFSId noVars> void MusclStretched_();

  void         Muscl(ZFSId timerId = -1);
  void         MusclRANS();
  void         (ZFSStrctrdBlck3D::*reconstructSurfaceData)();
  //different Muscl schemes:
  void         MusclVenkatakrishnan3D();
  void         MusclAlbada();
  void         MusclMinModLimiter();
  void         computeCellLength();
  ZFSFloat**   m_cellLength;

  void         computeReconstructionConstantsSVD();
  ZFSFloat RBF(const ZFSFloat, const ZFSFloat);
  ZFSFloat     computeRecConstSVD(const ZFSId ijk, const ZFSId noNghbrIds,ZFSInt* nghbr,  ZFSInt ID, ZFSInt sID, ZFSFloatScratchSpace& tmpA, ZFSFloatScratchSpace& tmpC, ZFSFloatScratchSpace& weights, const ZFSId recDim);
  ZFSId        gsl_linalg_SV_decomp_jacobi (ZFSFloatScratchSpace& A, ZFSFloatScratchSpace& Q, ZFSFloatScratchSpace& S);
  ZFSFloat     svdPreSolve( ZFSFloatScratchSpace& A, ZFSFloatScratchSpace& weights, ZFSFloatScratchSpace& AInv, const ZFSId m, const ZFSId n );

  ZFSBool      QRdec( ZFSFloatScratchSpace& a, ZFSFloatScratchSpace& d, const ZFSId m, const ZFSId n );
  void         getQ( ZFSFloatScratchSpace& qr, ZFSFloatScratchSpace& q, const ZFSId m, const ZFSId n );
  void         getR( ZFSFloatScratchSpace& qr, ZFSFloatScratchSpace& d, ZFSFloatScratchSpace& r, const ZFSId m, const ZFSId n );


  //limiter functions for MusclVenkatakrishnan3D()
  void VENKATAKRISHNAN_MOD_FCT(ZFSFloat effNghbrDelta, ZFSFloat srfcDelta, ZFSFloat dxEpsSqr, ZFSId cellPos, ZFSId var, ZFSFloatScratchSpace& minPhi);
  void VENKATAKRISHNAN_FCT(ZFSFloat effNghbrDelta, ZFSFloat srfcDelta, ZFSFloat dxEpsSqr, ZFSId cellPos, ZFSId var, ZFSFloatScratchSpace& minPhi);
  void BARTH_JESPERSON_FCT(ZFSFloat effNghbrDelta, ZFSFloat srfcDelta, ZFSFloat dxEpsSqr, ZFSId cellPos, ZFSId var, ZFSFloatScratchSpace& minPhi);
  void         (ZFSStrctrdBlck3D::*Venkatakrishnan_function)(ZFSFloat, ZFSFloat, ZFSFloat, ZFSId, ZFSId, ZFSFloatScratchSpace&);


  void         computeDxt();
  bool         maxResidual();
  bool         rungeKuttaStep();
  void         updateSpongeLayer();
  void         addDisturbance();
  void         viscousFlux();
  void         viscousFluxRANS();
  void         viscousFluxLES();
  void         (ZFSStrctrdBlck3D::*viscFluxMethod)();
  void         moveGrid(ZFSBool isRestart, ZFSBool zeroPos) override;
  void         initMovingGrid() override;
  void         saveGrid() override;
  void         saveCellJacobian();
  //for the treatment of boundary condition
  void         applyViscousBoundaryCondition();
  void         applyInviscidBoundaryCondition();
  //save the skin-friction etc ...
  void         computeAuxData();
  void         saveLiftCoefficient(ZFSInt fileId);
  void         savePowerCoefficient(ZFSInt fileId);
  void         saveDragCoefficient(ZFSInt fileId);
  void         saveLiftDragToAsciiFile() override;

  //fftw is needed for the channel flow:
  void         initFFTW(fftw_complex* uPhysField, fftw_complex* vPhysField, fftw_complex* wPhysField, ZFSInt lx, ZFSInt ly, ZFSInt lz, ZFSInt noPeakModes);
  void         getFourierCoefficients(ZFSFloat* k, ZFSFloat k0, complex<ZFSFloat>* fourierCoefficient);
  ZFSFloat     randnormal(ZFSFloat mu, ZFSFloat sigma);
  void         loadRestartBC2600();
  void         loadRestartBC2601();
  void         loadRestartSTG(ZFSBool);

  void         gather();
  void         scatter();
  void         send();
  void         receive();
  void         exchange() override;
  void         nonBlockingExchange();

  void         gatherPoints() override;
  void         sendPoints() override;
  void         receivePoints() override;
  void         scatterPoints() override;
  void         exchangePointsSingularity();
  void         gatherPointsS();
  void         sendPointsS();
  void         receivePointsS();
  void         scatterPointsS();

  //for traveling wave averaging
  void         waveGather();
  void         waveScatter();
  void         waveSend();
  void         waveReceive();
  void         waveExchange() override;
  void         spanwiseWaveReorder() override;

  //for exchange of postprocessing field
  void         ppFillGhostCells();
  void         ppExtrapolateVariables();
  void         ppGather();
  void         ppScatter();
  void         ppSend();
  void         ppReceive();
  void         ppExchange();
  //for exchange of  averaged vars for zonal method
  void         averagedFillGhostCells() override; //junoh
  void         averagedExtrapolateVariables();
  void         averagedGather();
  void         averagedSend();
  void         averagedReceive();
  void         averagedScatter();
  void         averagedExchange();

  void         viscousFluxCorrection();

  //For time averaging of flow variables
  void         computeCumulativeAverage(ZFSBool forceReset) override; //junoh

  //Line and plane interpolation
  void         initLineInterpolation();
  void         saveOutputLines() override;

 protected:
 class ZFSStrctrdBndryCnd<3> *m_strctrdBndryCnd;

 //index variables
 static const ZFSId xsd = 0;
 static const ZFSId ysd = 1;
 static const ZFSId zsd = 2;

 ZFSFloat m_kineticEOld;

 //for debugging only: revert metric order for output
 void revertMetrics();

 inline void crossProduct( ZFSFloat*, ZFSFloat*, ZFSFloat*);
 inline ZFSId   cellIndex(ZFSInt i, ZFSInt j, ZFSInt k);
 inline ZFSId   pointIndex(ZFSInt i, ZFSInt j, ZFSInt k);
 inline ZFSFloat dist(ZFSFloat* a, ZFSFloat* b);
 //routine to determine point ID from given cell ( point (0,0,0) for unit cube )
 inline ZFSId getPointIdFromCell( ZFSInt i, ZFSInt j, ZFSInt k );
 //routine to determine point ID from given point
 inline ZFSId getPointIdfromPoint( ZFSId origin, ZFSInt incI, ZFSInt incJ, ZFSInt incK );
 inline ZFSId getCellIdfromCell( ZFSId origin, ZFSInt incI, ZFSInt incJ, ZFSInt incK );
 inline ZFSId surfId(ZFSId point, ZFSId isd, ZFSId dim);

 class ZFSStrctrdBlck3DRans* m_ransBlck;


 //interpolation
 ZFSStrctrdInterpolation<3>* m_strctrdInterpolation;
 ZFSStrctrdInterpolation<3>* m_pointInterpolation;
 vector<ZFSStrctrdInterpolation<3>*> m_zonalV; //junoh
 void manualInterpolationCorrection();
 void interpolateFromDonor();

 //postprocessing helper functions
 void loadSampleFile(ZFSString) override;
 void getSampleVariables(ZFSId, ZFSFloat*) override;
 ZFSFloat getSampleVorticity(ZFSId, ZFSId) override;
 void loadAverageRestartFile(const char*,  ZFSFloat**, ZFSFloat**, ZFSFloat**,ZFSFloat**) override;
 void loadAveragedVariables(const char*) override;
 void shiftAverageCellValuesRestart();
 void shiftAverageCellValues();

 ZFSFloat dvardxyz(ZFSId, ZFSId, ZFSFloat*) override;
 ZFSFloat dvardx(ZFSId, ZFSFloat*) override ;

 static constexpr const ZFSInt nDim = 3;
};





#endif
