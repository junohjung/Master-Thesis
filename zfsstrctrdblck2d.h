#ifndef ZFSSTRCTRDBLCK2D_H
#define ZFSSTRCTRDBLCK2D_H

#include "zfsstrctrdblck.h"
#include "zfsstrctrdbndrycnd2d.h"
#include "zfsstrctrdblck2drans.h"

extern ZFSId globalTimeStep;


class ZFSStrctrdBlck2D : public ZFSStrctrdBlck<2>
{

  template<ZFSBool isRans> friend class ZFSStrctrdBndryCnd2D;
  friend class ZFSStrctrdBlck2DRans;

 public:
  ZFSStrctrdBlck2D( ZFSId, bool *, const MPI_Comm comm );
  ~ZFSStrctrdBlck2D();


  virtual void initialCondition();
  virtual void initSolutionStep();

  void      computeCellCentreCoordinates();
  void      addGhostPointCoordinateValues();
  void      extrapolateGhostPointCoordinates();
  void      extrapolateGhostPointCoordinatesBC();


  void initMovingGrid();

  void         initFluxMethod();
  void         assignBndryCells();
  void         initBndryCnds();
  void         applyBoundaryCondition();
  void         loadRestartBC2600();
  void         computePrimitiveVariables() override;

  virtual void exchange();
  void         gather();
  void         scatter();
  void         send();
  void         receive();

  virtual void gatherPoints();
  virtual void sendPoints();
  virtual void receivePoints();
  virtual void scatterPoints();

  //computing
  virtual void computeSurfaceMetrics();
  virtual void computeModSurfaceMetrics();
  virtual void computeCornerMetrics();
  virtual void computeModCornerMetrics();
  virtual void computeCellMetrics();

  virtual void computeCornerJacobian();
  virtual void computeCellJacobian();
  virtual void computeModCornerJacobian();
  virtual void computeModCellJacobian();
  void computeCellLength();

  virtual void computeTimeStep();
  bool         maxResidual();
  bool         rungeKuttaStep();
  void         updateSpongeLayer();
  void         viscousFlux();
  void         viscousFluxRANS();
  void         viscousFluxLES();
  void         (ZFSStrctrdBlck2D::*viscFluxMethod)();

  inline ZFSId    getPointIdFromCell(ZFSInt i, ZFSInt j);
  inline ZFSId    cellIndex(ZFSInt i, ZFSInt j);
  inline ZFSId    pointIndex(ZFSInt i, ZFSInt j);
  inline ZFSId    getPointIdFromPoint( ZFSId origin, ZFSInt incI, ZFSInt incJ );
  inline ZFSId    getCellIdFromCell( ZFSId origin, ZFSInt incI, ZFSInt incJ );
  inline ZFSFloat crossProduct( ZFSFloat vec1[2], ZFSFloat vec2[2]);
         ZFSFloat pressure(ZFSId cellId);
  inline ZFSFloat getPSI(ZFSId, ZFSId);

  void         Ausm();
  inline void  AusmLES(ZFSFloat* QLeft, ZFSFloat* QRight, const ZFSId dim, const ZFSId cellId);
  inline void  AusmLES_PTHRC(ZFSFloat* QLeft, ZFSFloat* QRight, ZFSId dim, ZFSId cellId);
  inline void  AusmDV(ZFSFloat* QLeft, ZFSFloat* QRight, const ZFSId dim, const ZFSId cellId);
  inline void  AusmLES_MovingGrid(ZFSFloatScratchSpace& QLeft, ZFSFloatScratchSpace& QRight, ZFSId dim, ZFSId cellId);

  typedef void (ZFSStrctrdBlck2D::*fluxmethod)(ZFSFloat*, ZFSFloat*, ZFSId, ZFSId);
  template< fluxmethod ausm, ZFSId noVars> void Muscl_();
  template< fluxmethod ausm, ZFSId noVars> void MusclStretched_();

 //function pointer for the Muscl_scheme
  void         MusclRANS();
  void         (ZFSStrctrdBlck2D::*reconstructSurfaceData)();

  //different Muscl schemes
  void         Muscl(ZFSId timerId = -1);
  void         MusclAlbada();
  void         MusclNoLimiter();

 class ZFSStrctrdBlck2DRans* m_ransBlck;

 protected:

  class ZFSStrctrdBndryCnd<2> *m_strctrdBndryCnd;
  static constexpr const ZFSInt nDim = 2;

};

#endif
