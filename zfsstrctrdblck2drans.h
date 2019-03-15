#ifndef ZFSSTRCTRDBLCK_2DRANS
#define ZFSSTRCTRDBLCK_2DRANS

#include "zfstypes.h"
#include "zfsscratch.h"
#include "zfsstrctrdblck.h"
#include "zfsstrctrdblck2d.h"
#include "zfsransmodelconstants.h"

class ZFSStrctrdBlck2DRans
{
 public:
  ZFSStrctrdBlck2DRans(class ZFSStrctrdBlck2D* block);
  ~ZFSStrctrdBlck2DRans();

  //void setBndryCndObject(class ZFSStrctrdBndryCnd2D* object, ZFSId noSpecies);
  void initFluxMethod();

  typedef void (ZFSStrctrdBlck2DRans::*fluxmethod)(ZFSFloat*, ZFSFloat*, ZFSId, ZFSId);
  template< fluxmethod ausm, ZFSId noVars> void Muscl_();
  void         (ZFSStrctrdBlck2DRans::*reconstructSurfaceData)();
  void Muscl();

  void Ausm(ZFSFloat* QLeft, ZFSFloat* QRight, const ZFSId dim, const ZFSId cellId);
  void AusmDV(ZFSFloat* QLeft, ZFSFloat* QRight, const ZFSId dim, const ZFSId cellId);

  void viscousFluxRANS();
  void viscousFlux_SA(); //SPALART-ALLMARAS
  void viscousFlux_FS(); //FARES-SCHROEDER
  void viscousFlux_KOmega();     //K-OMEGA MODELL
  void (ZFSStrctrdBlck2DRans::*viscFluxMethod) ();
  void (ZFSStrctrdBlck2DRans::*convFluxMethod) (ZFSFloat*,ZFSFloat*, ZFSId, ZFSId);

 protected:
  class ZFSStrctrdBndryCnd2DRans *m_strctrdBndryCndRans;


  void computeTurbViscosity();
 private:
  ZFSString m_ransMethod;
  ZFSStrctrdBlck2D* m_block;

  MPI_Comm m_zfsStrctrdComm;
  ZFSId m_blockId;
  ZFSInt* m_nCells;
  ZFSInt* m_nPoints;
  ZFSInt m_noStrctrdCells;
  ZFSStrctrdCell* m_cells;
  ZFSFloat** m_coordinates;
  ZFSConservativeVariables<2> * CV;
  ZFSPrimitiveVariables<2> *    PV;
  ZFSStrctrdFQVariables *    FQ;
  ZFSInt                     m_noGhostLayers;
  ZFSFloat m_eps;
  ZFSFloat m_chi;
  static constexpr const ZFSInt nDim = 2;
  ZFSFloat m_sutherlandConstant;
  ZFSFloat m_sutherlandPlusOne;

  const ZFSId xsd = 0;
  const ZFSId ysd = 1;
  const ZFSId zsd = 2;

  SAModelConstants m_sa;

  inline ZFSId   cellIndex(ZFSInt i, ZFSInt j);
  inline ZFSId getCellIdfromCell( ZFSId origin, ZFSInt incI, ZFSInt incJ);
  inline ZFSFloat getPSI(ZFSId, ZFSId);
};

#endif
