#ifndef ZFSSTRCTRDBLCK_3DRANS
#define ZFSSTRCTRDBLCK_3DRANS

#include "zfstypes.h"
#include "zfsscratch.h"
#include "zfsstrctrdblck.h"
#include "zfsstrctrdblck3d.h"
#include "zfsransmodelconstants.h"

class ZFSStrctrdBlck3DRans
{
 public:
  ZFSStrctrdBlck3DRans(class ZFSStrctrdBlck3D* block);
  ~ZFSStrctrdBlck3DRans();

  void initFluxMethod();

  typedef void (ZFSStrctrdBlck3DRans::*fluxmethod)(ZFSFloat*, ZFSFloat*, ZFSId, ZFSId);
  template< fluxmethod ausm, ZFSId noVars> void Muscl_();
  void (ZFSStrctrdBlck3DRans::*reconstructSurfaceData)();
  void Muscl();

  void AusmDV(ZFSFloat* QLeft, ZFSFloat* QRight, const ZFSId dim, const ZFSId cellId);

  void viscousFluxRANS();
  void viscousFlux_SA(); //SPALART-ALLMARAS
  void viscousFlux_FS(); //FARES-SCHROEDER
  void viscousFlux_KOmega(); //K-OMEGA MODELL
  void (ZFSStrctrdBlck3DRans::*viscFluxMethod) ();
  void (ZFSStrctrdBlck3DRans::*convFluxMethod) (ZFSFloat*,ZFSFloat*, ZFSId, ZFSId);
  void computeTurbViscosity();

 protected:
  class ZFSStrctrdBndryCnd3DRans *m_strctrdBndryCndRans;

 private:
  ZFSString m_ransMethod;
  ZFSStrctrdBlck3D* m_block;

  MPI_Comm m_zfsStrctrdComm;
  ZFSId m_blockId;
  ZFSInt* m_nCells;
  ZFSInt* m_nPoints;
  ZFSInt m_noStrctrdCells;
  ZFSStrctrdCell* m_cells;
  ZFSFloat** m_coordinates;
  ZFSConservativeVariables<3> * CV;
  ZFSPrimitiveVariables<3> *    PV;
  ZFSStrctrdFQVariables *    FQ;
  ZFSInt                     m_noGhostLayers;
  ZFSFloat m_eps;
  ZFSFloat m_chi;
  static constexpr const ZFSInt nDim = 3;
  ZFSFloat m_sutherlandConstant;
  ZFSFloat m_sutherlandPlusOne;

  const ZFSId xsd = 0;
  const ZFSId ysd = 1;
  const ZFSId zsd = 2;

  SAModelConstants m_sa;

  inline ZFSId   cellIndex(ZFSInt i, ZFSInt j, ZFSInt k);
  inline ZFSId getCellIdfromCell( ZFSId origin, ZFSInt incI, ZFSInt incJ, ZFSInt incK );
  inline ZFSFloat getPSI(ZFSId, ZFSId);
};

#endif
