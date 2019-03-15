#ifndef ZFSSTRCTRDBNDRYCND3D
#define ZFSSTRCTRDBNDRYCND3D

#include "zfsstrctrdblck3d.h"
#include "zfsstrctrdbndrycnd.h"
#include "zfsstrctrdblck.h"

template <ZFSBool isRans>
class ZFSStrctrdBndryCnd3D : public ZFSStrctrdBndryCnd<3>
{
 public:
  friend class ZFSStrctrdBlck3D;

  ZFSStrctrdBndryCnd3D(class ZFSStrctrdBlck<3>* block, ZFSId );
  ~ZFSStrctrdBndryCnd3D();

  void         createSpongeAtBndryCnd();
  void         correctBndryCndIndices();

  ZFSInt       m_currentPeriodicDirection;
  ZFSBool      skipPeriodicDirection(ZFSId);
  virtual void periodicExchange();
  void         periodicNonBlockingExchange();
  void         (ZFSStrctrdBndryCnd3D::*gatherPeriodic)();
  void         gather();
  void         gatherPeriodicRotation();
  void         scatter();
  void         send();
  void         receive();

  void         exchangePointsPeriodic();
  void         gatherPoints();
  void         receivePoints();
  void         sendPoints();
  void         scatterPoints();
  void         periodicPointsChange(ZFSFloat* pt,ZFSId type);

  void         exchangePointsPeriodicS();
  void         gatherPointsS();
  void         receivePointsS();
  void         sendPointsS();
  void         scatterPointsS();

  //void         applyDirichletNeumannBC();
  void         readAndDistributeSpongeCoordinates();
  void         updateSpongeLayer();
  void         computeWallDistances();
  void         setRotationalBCProperties();



  inline ZFSFloat dist(ZFSFloat* a, ZFSFloat* b);

  void         bc1000(ZFSId);//wall no slip
  void         bc1003(ZFSId);//isothermal no slip wall
  void         bc1004(ZFSId);//moving adiabatic wall
  void         bc1006(ZFSId);//moving isothermal wall
  void         bc1007(ZFSId);//oscillating wall
  void         bc2001(ZFSId);//subsonic inflow
  void         bc2004(ZFSId);//subsonic outflow
  void         bc2003(ZFSId);//simple subsonic in/outflow
  void         bc2002(ZFSId);//supersonic inflow
  void         bc2005(ZFSId);//supersonic outflow
  void         bc2009(ZFSId);//supersonic outflow after shock
  void         bc2020(ZFSId);//poiseulle flow inflow
  void         bc2099(ZFSId);//subsonic inflow (u=(y/d)^(1/7))
  void         bc2012(ZFSId){};//characteristic inflow
  void         bc2013(ZFSId){};//characteristic outflow
  void         bc2014(ZFSId);//subsonic rotational inflow

  void         bc2222(ZFSId);//subsonic RANS outflow bc //junoh

  void         bc2300(ZFSId){};//blasius
  void         bc2402(ZFSId);//channel flow
  void         bc2500(ZFSId);//Rescaling: recycle station
  void         bc2501(ZFSId);//Rescaling: inlet station
  void         bc2510(ZFSId);//Rescaling: recycle station RANS
  void         bc2511(ZFSId);//Rescaling: inlet station RANS
  void         bc2600(ZFSId);//Prescribing profile
  void         bc2601(ZFSId);//Prescribing profile
  void         bc2700(ZFSId);//mode inflow
  void         bc2900(ZFSId);//Jet Inlet Freund
  void         bc3000(ZFSId);//symmetry
  void         bc4001(ZFSId){};//periodic rotation
  void         bc6000(ZFSId);//communication
  void         bc7909(ZFSId);//synthetic turbulence generation

  void         initBc1000(ZFSId);//wall no slip
  void         initBc1003(ZFSId);//isothermal no slip wall
  void         initBc1004(ZFSId);// moving adiabatic wall
  void         initBc1006(ZFSId);// moving isothermal wall
  void         initBc1007(ZFSId);// oscillating wall
  void         initBc2001(ZFSId){};//subsonic inflow
  void         initBc2002(ZFSId){};//supersonic inflow
  void         initBc2003(ZFSId){};//simple subsonic in/outflow
  void         initBc2004(ZFSId);//subsonic outflow
  void         initBc2005(ZFSId){};//supersonic outflow
  void         initBc2009(ZFSId);//supersonic outflow after shock
  void         initBc2012(ZFSId){};//characteristic inflow
  void         initBc2013(ZFSId){};//characteristic outflow
  void         initBc2014(ZFSId){};//subsonic rotational bc
  void         initBc2020(ZFSId){};//poiseulle flow
  void         initBc2099(ZFSId){};//subsonic inflow (u=(y/d)^(1/7))
  void         initBc2300(ZFSId);//blasius solution
  void         initBc2402(ZFSId);//channel flow
  void         initBc2500(ZFSId);//Rescaling: recycle station
  void         initBc2501(ZFSId){};//Rescaling: inlet station
  void         initBc2600(ZFSId);//Prescribing profile
  void         initBc2601(ZFSId);//Prescribing profile
  void         initBc2700(ZFSId);//mode inflow
  void         initBc2900(ZFSId){};//Jet Inlet Freund
  void         initBc3000(ZFSId){};//symmetry
  void         initBc4001(ZFSId);//periodic rotation
  void         initBc6000(ZFSId){};//communication
  void         initBc7909(ZFSId);//synthetic turbulence generation

  //empty BC if BC==-1
 void         bc9999(ZFSId){};
 void         initBc9999(ZFSId){};

 inline void  crossProduct(ZFSFloat*, ZFSFloat*, ZFSFloat*);
 ZFSId        cellIndex(ZFSInt i, ZFSInt j, ZFSInt k);
 ZFSId        cellIndexBC(ZFSInt i, ZFSInt j, ZFSInt k); //For STG, only the first three rows are used

  void         computeAuxData();
  void         computeAuxDataRoot();
  template<ZFSBool computePower>
  void         computeFrictionPressureCoef();
  void         computeLiftCoef();
  void         computeDragCoef();
  void         computeLiftCoefRoot();
  void         computeDragCoefRoot();
  void         computePowerCoefRoot();
  void         computeMomentCoef();

 protected:
  class ZFSStrctrdBlck3D* m_block;
  ZFSId                   m_noSpecies;
  ZFSId                   m_noStrctrdCells;
  ZFSInt                  m_startCommPeriodic;
  ZFSInt                  m_endCommPeriodic;
  ZFSInt                  m_periodicS;
  ZFSInt                  m_startCommChannel;
  ZFSInt                  m_endCommChannel;
  ZFSInt                  m_channelInflowRank;

  //periodic
  ZFSInt                  m_noPeriodicConnections;


  //7909 synthetic turbulence generation
//  ZFSFloat** m_stgEddies;
  ZFSFloat generate_rand();
  ZFSFloat generate_rand_weighted();
  ZFSFloat* m_stgVbStart;
  ZFSFloat* m_stgVbEnd;
  ZFSFloat* m_stgMaxVel;
  ZFSFloat** m_stgGlobalLengthScales;

  //2700 mode bc
  ZFSFloat   m_isothermalWallTemperature;
  ZFSFloat   m_Ma;
  ZFSFloat*  m_modeAmp;   //for mode inflow
  ZFSFloat   m_modeSr;
  ZFSInt*    m_modeType;
  ZFSFloat*  m_modePhi;
  ZFSFloat   m_modeAngle;
  ZFSInt*    m_nmbrOfModes;
  ZFSInt     m_modes;
  ZFSFloat*  m_modeOmega;
  ZFSFloat*  m_modeEtaMin;
  ZFSFloat** m_modeK;

  //2601
  ZFSFloat** m_2601effConst;
  ZFSFloat*  m_2601streamwisePos;
  ZFSBool    m_2601wave;
  ZFSInt     m_2601noSteps;
  ZFSInt     m_2601noCoeff;
  ZFSInt     m_2601noPos;

  //2500 rescaling bc
  ZFSFloat m_rescalingBLT;

  //routine to determine point ID from given cell
  inline ZFSId    getPointIdFromCell( ZFSInt i, ZFSInt j, ZFSInt k );
  inline ZFSId    pointIndex(ZFSInt i, ZFSInt j, ZFSInt k);
  //routine to determine point ID from given point
  inline ZFSId    getPointIdfromPoint( ZFSId origin, ZFSInt incI, ZFSInt incJ, ZFSInt incK );
  inline ZFSFloat pressure(ZFSId cellId);
  inline ZFSFloat pressure(ZFSId i, ZFSId j, ZFSId k);
  inline ZFSFloat temperature(ZFSId cellId);
  inline ZFSId getReverseCellId(ZFSId i, ZFSId j, ZFSId k, ZFSId face);
  inline ZFSId getExtrNghbrId(ZFSId cellId, ZFSId face);
  inline pair<ZFSId, ZFSId> getMirrorCellIdPair(ZFSId i, ZFSId j, ZFSId k, ZFSId face);

  static constexpr ZFSInt m_reverseCellIdDim[18] = {-1, 1, 1, 1,1,1,
                                                    1,-1, 1, 1,1,1,
                                                    1, 1,-1, 1,1,1};

  static constexpr ZFSInt m_reverseCellIdGC[18] = {1,0,0, 0,0,0,
                                                   0,1,0, 0,0,0,
                                                   0,0,1, 0,0,0};
  static constexpr const ZFSInt nDim = 3;
};

#endif
