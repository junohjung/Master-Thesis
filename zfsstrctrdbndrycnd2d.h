#ifndef ZFSSTRCTRDBNDRYCND2D
#define ZFSSTRCTRDBNDRYCND2D

#include "zfsstrctrdblck.h"
#include "zfsstrctrdblck2d.h"
#include "zfsstrctrdbndrycnd.h"

template <ZFSBool isRans>
class ZFSStrctrdBndryCnd2D : public ZFSStrctrdBndryCnd<2>
{
 public:
  friend class ZFSStrctrdBlck2D;

  ZFSStrctrdBndryCnd2D(class ZFSStrctrdBlck<2>* block, ZFSId );
  ~ZFSStrctrdBndryCnd2D();

  void         correctBndryCndIndices();

  void         bc1000(ZFSId);//wall no slip
  void         bc1001(ZFSId);//euler wall
  void         bc1003(ZFSId);//isothermal wall
  void         bc2001(ZFSId);//subsonic inflow
  void         bc2004(ZFSId);//subsonic outflow
  void         bc2002(ZFSId);//supersonic inflow
  void         bc2005(ZFSId);//supersonic outflow
  void         bc2007(ZFSId);//subsonic outflow
  void         bc2021(ZFSId);//shear flow inflow
  void         bc2199(ZFSId);//subsonic inflow compressible bernoulli (tfs2099)
  void         bc2300(ZFSId);//blasius
  void         bc2510(ZFSId);//rescaling inlet
  void         bc2511(ZFSId);//rescaling recycling
  void         bc2600(ZFSId);//prescribing
  void         bc3000(ZFSId);//symmetrie

  void         initBc1000(ZFSId);
  void         initBc1001(ZFSId){};
  void         initBc1003(ZFSId);
  void         initBc2002(ZFSId);
  void         initBc2004(ZFSId);
  void         initBc2005(ZFSId);
  void         initBc2021(ZFSId);
  void         initBc2199(ZFSId){};
  void         initBc2300(ZFSId);//blasius solution
  void         initBc2500(ZFSId){};
  void         initBc2501(ZFSId){};
  void         initBc2510(ZFSId);
  void         initBc2511(ZFSId){};
  void         initBc2600(ZFSId);
  void         initBc3000(ZFSId);

  void           readAndDistributeSpongeCoordinates();
  void           updateSpongeLayer();
  void           computeWallDistances();
  inline ZFSId   cellIndex(ZFSInt i, ZFSInt j);
  inline ZFSId   getPointIdFromCell( ZFSInt i, ZFSInt j);
  inline ZFSId   getPointIdFromPoint( ZFSId origin, ZFSInt incI, ZFSInt incJ);
  inline ZFSFloat pressure(ZFSId);
  inline ZFSFloat temperature(ZFSId);

 protected:
  class ZFSStrctrdBlck2D*  m_block;
  ZFSId                    m_noSpecies;
  ZFSId                    m_noStrctrdCells;
  ZFSInt                   m_startCommPeriodic;
  ZFSInt                   m_endCommPeriodic;
  ZFSInt                   m_startCommChannel;
  ZFSInt                   m_endCommChannel;

  //2500 rescaling bc
  ZFSFloat m_rescalingBLT;
  ZFSFloat m_isothermalWallTemperature;
  ZFSFloat m_bc2021Gradient;

  static constexpr const ZFSInt nDim = 2;
};


#endif
