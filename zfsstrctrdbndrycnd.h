#ifndef ZFSSTRCTRDBNDRYCND
#define ZFSSTRCTRDBNDRYCND

#include "zfsglobals.h"
#include "zfsstrctrdcell.h"
#include "zfsstrctrdblck.h"
#include "zfsstrctrdwindowmapping.h"
#include "zfsstrctrdblckwindowinfo.h"
#include "zfsstrctrdfqvariables.h"
//#include "zfslist.h"

struct rotationalBC {
  ZFSFloat*                  perNormals;//Normals of the in/outflow planes
  ZFSInt                     hasRotationPeriodicBC;
  ZFSFloat                   perAngle;//Angle between the in- and outflow planes
  ZFSInt                     initRun;//iterator that is incremented each run of initBc4001 (important for domanis with both sides of rotational bc)
  ZFSInt                     firstRunBcId;//bcId of the side that was taken in the first run of initBc4001
  ZFSFloat*                  prevNormal;//normal of the side from the first run
  ZFSFloat**                 rotationMatrix4001;//rotational matrix for bc 4001
  ZFSFloat**                 rotationMatrix4002;//rotational matrix for bc 4002 (the reverse one of 4001)
};

// Forward declarations
template <ZFSInt nDim> class ZFSStrctrdBlck;

template <ZFSInt nDim>
class ZFSStrctrdBndryCnd
{
  template <ZFSInt nDim_> friend class ZFSStrctrdBlck;
  friend class ZFSStrctrdBlck3D;
  template <ZFSInt nDim_> friend class ZFSStrctrdBlckWindowInfo;
 public:
  //functions
  ZFSStrctrdBndryCnd(ZFSStrctrdBlck<nDim>* block);
  virtual ~ZFSStrctrdBndryCnd();
  void applyNonReflectingBC();
  void assignBndryCnds();
  void mapCpy(ZFSStrctrdWindowMap* input, ZFSStrctrdWindowMap* output);
  void applyDirichletNeumannBC();
  virtual void correctBndryCndIndices();
  virtual void periodicExchange(){};

  //MPI Communicator
  MPI_Comm m_zfsStrctrdComm;

  //variables
  ZFSInt*                    m_nCells;
  ZFSInt*                    m_nPoints;
  ZFSStrctrdCell*            m_cells;
  ZFSFloat**                 m_coordinates;
  ZFSStrctrdBlck<nDim>*            m_block;
  ZFSConservativeVariables<nDim> * CV;
  ZFSPrimitiveVariables<nDim> *    PV;
  ZFSStrctrdFQVariables *    FQ;
  ZFSId                      m_blockId;
  //stores the number of boundary conditions;
  ZFSId                              m_noBndryCndIds;

  ZFSInt                     m_noGhostLayers;
  //for Boundary conditions

  vector<ZFSStrctrdWindowMap*> bcInfoMap;
  ZFSId*                       bcIdInfoMap;
  vector<ZFSStrctrdWindowMap*> m_physicalBCMap;
  vector<ZFSStrctrdWindowMap*> m_auxDataMap;
  vector<ZFSStrctrdWindowMap*> m_globalStrctrdBndryMaps;

  ZFSInt       m_noSpongeDomainInfos;
  ZFSId*       m_spongeBcWindowInfo;
  //ZFSInt       m_useSponge;
  ZFSInt      m_spongeLayerType;
  ZFSFloat*    m_spongeLayerThickness;
  ZFSFloat*    m_sigmaSponge;
  ZFSFloat*    m_betaSponge;
  ZFSFloat     m_targetDensityFactor;
  /*
  bool                       m_createSpongeBoundary;
  ZFSFloat*                  m_spongeDistance;
  ZFSFloat                   m_spongeLayerThickness;
  ZFSId                      m_noSpongeBndryCndIds;
  ZFSFloat*                  m_sigmaSpongeBndryId;
  ZFSId*                     m_spongeBndryCndIds;
  ZFSId                      m_spongeLayerLayout;
  ZFSId                      m_noCellsInsideSpongeLayer;
  ZFSId*                     m_cellsInsideSpongeLayer;
  */
  ZFSFloat                   m_sigma;

  void                       saveAuxData();

  virtual void               bc0(ZFSId){};//nothing to do bc
  virtual void               bc1000(ZFSId){};//wall no slip
  virtual void               bc1001(ZFSId){};//euler wall
  virtual void               bc1003(ZFSId){};//isothermal no slip wall
  virtual void               bc1004(ZFSId){};// moving adiabatic wall
  virtual void               bc1006(ZFSId){};// moving isothermal wall
  virtual void               bc1007(ZFSId){};// oscillating wall
  virtual void               bc2001(ZFSId){};//subsonic inflow
  virtual void               bc2002(ZFSId){};//supersonic inflow
  virtual void               bc2003(ZFSId){};//simple subsonic in/outflow
  virtual void               bc2004(ZFSId){};//subsonic outflow
  virtual void               bc2024(ZFSId){};//subsonic outflow
  virtual void               bc2005(ZFSId){};//supersonic outflow
  virtual void               bc2007(ZFSId){};//supersonic outflow
  virtual void               bc2009(ZFSId){};//supersonic outflow after shock
  virtual void               bc2020(ZFSId){};//poiseulle flow inflow
  virtual void               bc2021(ZFSId){};//shear flow inflow
  virtual void               bc2099(ZFSId){};//subsonic inflow (u=(y/d)^(1/7))

  virtual void               bc2221(ZFSId){};//zonal with STG //junoh
  virtual void               bc2222(ZFSId){};//zonal without STG   //junoh

  virtual void               bc2199(ZFSId){};//subsonic inflow compressible bernoulli (tfs2099)x
  virtual void               bc2300(ZFSId){};//blasius solution
  virtual void               bc2402(ZFSId){};//channel flow
  virtual void               bc3000(ZFSId){};//symmetry
  virtual void               bc4001(ZFSId){};//periodic rotation
  virtual void               bc6000(ZFSId){};//communication
  virtual void               bc2012(ZFSId){};//characteristic inflow
  virtual void               bc2014(ZFSId){};//subsonic rotational inflow
  virtual void               bc2015(ZFSId){};//non reflecting outflow poinsot lele
  virtual void               bc7909(ZFSId){};//synthetic turbulence generation
  virtual void               bc2013(ZFSId){};//characteristic outflow
  virtual void               bc2500(ZFSId){};//Rescaling: recycle station
  virtual void               bc2511(ZFSId){};//Rescaling: inlet station RANS
  virtual void               bc2510(ZFSId){};//Rescaling: recycle station RANS
  virtual void               bc2501(ZFSId){};//Rescaling: inlet station
  virtual void               bc2600(ZFSId){};//Prescribing profile
  virtual void               bc2601(ZFSId){};//Prescribing profile
  virtual void               bc2700(ZFSId){};//mode inflow
  virtual void               bc2900(ZFSId){};//Jet Freund inlet
  virtual void               initBc0(ZFSId){};//wall no slip
  virtual void               initBc1000(ZFSId){};//wall no slip
  virtual void               initBc1001(ZFSId){};//euler wall
  virtual void               initBc1003(ZFSId){};//isothermal wall no slip
  virtual void               initBc1004(ZFSId){};// moving adiabatic wall
  virtual void               initBc1006(ZFSId){};// moving isothermal wall
  virtual void               initBc1007(ZFSId){};// oscillating wall
  virtual void               initBc2001(ZFSId){};//subsonic inflow
  virtual void               initBc2003(ZFSId){};//simple in/outflow
  virtual void               initBc2004(ZFSId){};//subsonic outflow
  virtual void               initBc2024(ZFSId){};//subsonic outflow
  virtual void               initBc2002(ZFSId){};//supersonic inflow
  virtual void               initBc2005(ZFSId){};//supersonic outflow
  virtual void               initBc2007(ZFSId){};//supersonic outflow
  virtual void               initBc2009(ZFSId){};//supersonic outflow after shock
  virtual void               initBc2020(ZFSId){};//poiseulle flow inflow
  virtual void               initBc2021(ZFSId){};//shear flow inflow
  virtual void               initBc2099(ZFSId){};//subsonic inflow (u=(y/d)^(1/7))
  
  virtual void               initBc2221(ZFSId){}; //zonal with STG //junoh
  virtual void               initBc2222(ZFSId){}; //zonal without STG //junoh

  virtual void               initBc2199(ZFSId){};//subsonic inflow compressible bernoulli (tfs2099)
  virtual void               initBc2300(ZFSId){};//blasius solution
  virtual void               initBc2402(ZFSId){};//channel flow
  virtual void               initBc3000(ZFSId){};//symmetry
  virtual void               initBc4001(ZFSId){};//periodic rotation
  virtual void               initBc6000(ZFSId){};//communication
  virtual void               initBc2012(ZFSId){};//characteristic inflow
  virtual void               initBc7909(ZFSId){};//synthetic turbulence generation
  virtual void               initBc2013(ZFSId){};//characteristic outflow
  virtual void               initBc2014(ZFSId){};//subsonic rotational inflow
  virtual void               initBc2015(ZFSId){};//non reflecting outflow poinsot lele
  virtual void               initBc2500(ZFSId){};//Rescaling: recycle station
  virtual void               initBc2501(ZFSId){};//Rescaling: inlet station
  virtual void               initBc2510(ZFSId){};//Rescaling: inlet station
  virtual void               initBc2600(ZFSId){};//Prescribing profile
  virtual void               initBc2601(ZFSId){};//Prescribing profile
  virtual void               initBc2700(ZFSId){};//mode inflow
  virtual void               initBc2900(ZFSId){};//Jet Freund inlet
  virtual void               bc9999(ZFSId){};
  virtual void               initBc9999(ZFSId){};

  virtual void               exchangePointsPeriodic(){};
  virtual void               exchangePointsPeriodicS(){};


  virtual void computeAuxData(){};
  virtual void computeAuxDataRoot(){};
  virtual void computeWallDistances(){};
  virtual void updateSpongeLayer(){};


  //function pointer to the boundary method used;
  typedef void ( ZFSStrctrdBndryCnd::* BndryCndHandler ) (ZFSId);

  //Dirichlet Conditions
  BndryCndHandler* nonReflectingBoundaryCondition;
  BndryCndHandler* bndryCndHandler;
  BndryCndHandler* initBndryCndHandler;
  //compute cf
  BndryCndHandler* skinFrictionHandler;

 protected:

  //for the channel flow we need the surfaces:
  // sorted for in and out if one partition contains both parts of the channel


  ZFSId                      m_noStrctrdCells;
  ZFSFloat                   m_channelSurfaceIn;
  ZFSFloat                   m_channelSurfaceOut;
  ZFSInt*                    m_channelSurfacIndexMap;
  //periodic boundary conditions
  rotationalBC               m_rotBC;

  ZFSFloat                   m_sutherlandPlusOne;
  ZFSFloat                   m_sutherlandConstant;
  // ZFSInt                     m_bCf;
  //ZFSInt                     m_bCp;
  ZFSInt                     m_bCfCpCoeff;
  ZFSInt                     m_bPower;
  ZFSInt                     m_bCl;
  ZFSInt                     m_bCd;
  ZFSInt                     m_bCpLineAveraging;
  ZFSInt                     m_cpAveragingDir;
  ZFSFloat*                  m_cD;
  ZFSFloat*                  m_cL;
  ZFSFloat*                  m_cDp;
  ZFSFloat*                  m_cLp;
  ZFSFloat*                  m_cDv;
  ZFSFloat*                  m_cLv;
  ZFSFloat*                  m_Powerv;
  ZFSFloat*                  m_Powerp;
  ZFSFloat*                  m_cArea;

  //infintiy Values!!!!
  ZFSFloat m_UInfinity;
  ZFSFloat m_VInfinity;
  ZFSFloat m_WInfinity;
  ZFSFloat m_PInfinity;
  ZFSFloat m_TInfinity;
  ZFSFloat m_DthInfinity;
  ZFSFloat m_muInfinity;
  ZFSFloat m_DInfinity;
  ZFSFloat m_VVInfinity[3];
  ZFSFloat m_rhoUInfinity;
  ZFSFloat m_rhoVInfinity;
  ZFSFloat m_rhoWInfinity;
  ZFSFloat m_rhoEInfinity;
  ZFSFloat m_rhoInfinity;
  ZFSFloat m_rhoVVInfinity[3];
};

#endif
