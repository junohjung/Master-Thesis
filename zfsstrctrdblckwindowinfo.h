#ifndef ZFSFVBLOCKSTRUCTWINDOWSINFO
#define ZFSFVBLOCKSTRUCTWINDOWSINFO

#include <vector>
#include <set>
#include "zfstypes.h"
#include "zfsstrctrdblckpartition.h"
#include "zfsstrctrdwindowmapping.h"
#include "zfsstrctrdblckcommunicationhandle.h"

#include <cstring>
using namespace std;



typedef struct
{
  ZFSInt BC;
  ZFSInt blockId;
  ZFSInt pos[3];
  ZFSInt found;
} pointType;

struct SingularInformation
{
  ZFSInt start[3];
  ZFSInt end[3];
  ZFSInt Nstar;
  ZFSInt displacement[5][3];
  ZFSInt count;
  ZFSInt totalPoints;
  ZFSInt totalCells;
  ZFSFloat **coordinates;
  ZFSFloat **variables;
  ZFSFloat **ReconstructionConstants;
  ZFSInt Viscous[3];
  ZFSInt BC;
  ZFSInt SingularBlockId[4];
};

class connectionNode
{
public:
  ZFSId blockId1;
  ZFSInt pos1[3];
  ZFSId blockId2;
  ZFSInt pos2[3];
  ZFSInt BC;
  ZFSInt Nstar;
  connectionNode()
  {BC=-1; blockId1=-1;blockId2=-1;
    pos1[0]=-1;pos1[1]=-1;pos1[2]=-1;
    pos2[0]=-1;pos2[1]=-1;pos2[2]=-1;Nstar=-1;};
  // connectionNode(ZFSInt a,ZFSId b1, ZFSInt p1[3], ZFSId b2, ZFSInt p2[3])
  //  {BC=a; blockId1=b1;blockId2=b2;
  //    pos1[0]=p1[0];pos1[1]=p1[1];pos1[2]=p1[2];
  //    pos2[0]=p2[0];pos2[1]=p2[1];pos2[2]=p2[2];};

  connectionNode(ZFSInt a,ZFSId b1, ZFSInt* p1, ZFSId b2, ZFSInt* p2)
  {BC=a; blockId1=b1;blockId2=b2;
    pos1[0]=p1[0];pos1[1]=p1[1];pos1[2]=p1[2];
    pos2[0]=p2[0];pos2[1]=p2[1];pos2[2]=p2[2];Nstar=-1;};

  connectionNode(ZFSInt a,ZFSId b1, ZFSInt* p1, ZFSId b2, ZFSInt* p2,ZFSBool enableSwap)
  {
    ZFSBool swap =false;

    if (enableSwap) {
      //ensure a fixed order, i.e., lower blockId first or lower pos first if the blockIds match
      if (b2 < b1) swap = true;
      if (b1 == b2) {
        for(ZFSInt countDim = 0;countDim<3;++countDim)
        {
          if (p1[countDim] < p2[countDim]) break;
          if (p1[countDim] > p2[countDim])
          {
            swap = true;
            break;
          }
        }
      }
    }

    if (swap) {
      BC=a;
      blockId1=b2;blockId2=b1;
      pos1[0]=p2[0];pos1[1]=p2[1];pos1[2]=p2[2];
      pos2[0]=p1[0];pos2[1]=p1[1];pos2[2]=p1[2];
    }
    else {
      BC=a;
      blockId1=b1;blockId2=b2;
      pos1[0]=p1[0];pos1[1]=p1[1];pos1[2]=p1[2];
      pos2[0]=p2[0];pos2[1]=p2[1];pos2[2]=p2[2];
    }
    Nstar=-1;
  };

  void print()
  {
    cout<<"============================"<<endl;
    cout<<BC<<endl;
    cout<<blockId1<<" "<<pos1[0]<<" "<<pos1[1]<<" "<<pos1[2]<<endl;
    cout<<blockId2<<" "<<pos2[0]<<" "<<pos2[1]<<" "<<pos2[2]<<endl;
    cout<<"Nstar: "<<Nstar<<endl;
    cout<<"============================"<<endl;
  };
  ZFSBool operator< (const connectionNode& entry2) const {

    if (BC < entry2.BC) {return true;}
    if (BC > entry2.BC) {return false;}
    if (blockId1 < entry2.blockId1) {return true;}
    if (blockId1 > entry2.blockId1) {return false;}
    if (blockId2 < entry2.blockId2) {return true;}
    if (blockId2 > entry2.blockId2) {return false;}
    for(ZFSInt i=0;i<3;++i)
    {
      if (pos1[i] < entry2.pos1[i]) {return true;}
      if (pos1[i] > entry2.pos1[i]) {return false;}
      if (pos2[i] < entry2.pos2[i]) {return true;}
      if (pos2[i] > entry2.pos2[i]) {return false;}
    }
    return false;
  };
  ZFSBool operator== (const connectionNode& entry2) const {

    if (BC != entry2.BC) {return false;}
    if (blockId1 != entry2.blockId1) {return false;}
    if (blockId2 != entry2.blockId2) {return false;}
    for(ZFSInt i=0;i<3;++i)
    {
      if (pos1[i] != entry2.pos1[i]) {return false;}
      if (pos2[i] != entry2.pos2[i]) {return false;}
    }
    return true;
  };
};


struct connectivity {
  ZFSId window;
  ZFSInt* permutation;
  ZFSInt* order;
  ZFSInt* increments; //if HangingNodes will be used
};


class windowInformation
{
public:
  windowInformation(ZFSInt ndx){cpu =-1; BC=-1;startindex = new ZFSInt[ndx]; endindex= new ZFSInt[ndx];};
  ~windowInformation(){delete[] startindex; delete[] endindex;};
  //int index;
  connectivity* connec;
  ZFSId fixDir;
  ZFSId inDir1;
  ZFSId inDir2;
  ZFSId* startindex;
  ZFSId* endindex;
  //ZFSId* permutation;
  ZFSId cpu;
  ZFSId BC;
  ZFSId inputBlockId;
  ZFSId windowId;

  //private:
};

struct noBoxWindows {
  //for information about the order of the numbering see tfs manual
  vector<vector<windowInformation*> > sides;

};

template <ZFSInt nDim>
class  ZFSStrctrdBlckWindowInfo
{
  template <ZFSInt nDim_> friend class ZFSStrctrdBlck;
  friend class ZFSStrctrdBlck3D;
  template <ZFSInt nDim_> friend class ZFSStrctrdBndryCnd;
  template <ZFSInt nDim_> friend class ZFSStrctrdDecomposition;
public:
  ZFSStrctrdBlckWindowInfo(ZFSStrctrdDecomposition<nDim>* part, ZFSId file_ident, ZFSInt noInputBlocks_, ZFSInt inputBlockId_, ZFSInt GhostLayers, MPI_Comm strctrdCommunicator, const ZFSId noDomains_, const ZFSId domainId_);
  ~ZFSStrctrdBlckWindowInfo();
  void readNewWindowInfo();
  void readNewWindowInfo(ZFSString gridFileName);
  void readWindowInfo();
  void createWindowInfo();
  void mapCreate(ZFSId Id1, ZFSInt* start1, ZFSInt* end1, ZFSInt* step1, ZFSId Id2=-1, ZFSInt* start2=NULL, ZFSInt* end2=NULL, ZFSInt* step2=NULL, ZFSInt* order=NULL, ZFSInt BC=-1, ZFSStrctrdWindowMap* output=NULL);

  void setSpongeInformation(ZFSInt noSpongeInfo, ZFSFloat* beta, ZFSFloat* sigma, ZFSFloat* thickness, ZFSInt* bcInfo, ZFSInt informationType);
  void setWallInformation();
  void setZonalBCInformation(); //junoh

  ZFSBool mapCheck(ZFSStrctrdWindowMap* input);
  ZFSBool mapCheck0d(ZFSStrctrdWindowMap* input);
  ZFSBool mapCheck1d(ZFSStrctrdWindowMap* input);
  ZFSBool mapCheck3d(ZFSStrctrdWindowMap* input);
  ZFSBool mapCheck2d(ZFSStrctrdWindowMap* input);
  ZFSBool mapCheckWave(ZFSStrctrdWindowMap* input);
  void mapPrint(ZFSStrctrdWindowMap* input);
  void mapPrintSimple(ZFSStrctrdWindowMap* input);
  void mapZero(ZFSStrctrdWindowMap* output );
  int mapCompare(ZFSStrctrdWindowMap* map1, ZFSStrctrdWindowMap* map2);
  int mapCompare11(ZFSStrctrdWindowMap* map1, ZFSStrctrdWindowMap* map2);
  void mapInvert(ZFSStrctrdWindowMap* input,ZFSStrctrdWindowMap* output );
  void mapInvert1(ZFSStrctrdWindowMap* output );
  void mapNormalize1(ZFSStrctrdWindowMap* input, ZFSStrctrdWindowMap* output);
  void mapNormalize2(ZFSStrctrdWindowMap* input, ZFSStrctrdWindowMap* output);
  void mapNormalize3(ZFSStrctrdWindowMap* output);
  void mapCombine11(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2,ZFSStrctrdWindowMap* output);
  void mapCombine12(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2,ZFSStrctrdWindowMap* output);
  void mapCombine21(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2,ZFSStrctrdWindowMap* output);
  void mapCombine22(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2,ZFSStrctrdWindowMap* output);
  void mapCombineWave(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2,ZFSStrctrdWindowMap* output);
  void mapCombineCell11(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2,ZFSStrctrdWindowMap* output);
  void mapCombineCell12(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2,ZFSStrctrdWindowMap* output);
  void mapCombineCell21(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2,ZFSStrctrdWindowMap* output);
  void mapCombineCell22(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2,ZFSStrctrdWindowMap* output);

  void mapCpy(ZFSStrctrdWindowMap* input, ZFSStrctrdWindowMap* output );
  void readMapFromArray(ZFSStrctrdWindowMap* map, ZFSInt* array);
  void writeMapToArray(ZFSStrctrdWindowMap* map, ZFSInt* array);
  void initGlobals();
  void writeConnectionWindowInformation(ZFSFloat* periodicDisplacements);
  void readConnectionWindowInformation(ZFSFloat* periodicDisplacments, ZFSString fileName);
  void createWindowMapping(MPI_Comm* channelIn, MPI_Comm* channelOut, MPI_Comm* channelWorld, ZFSInt* channelRoots, MPI_Comm* commStg, ZFSInt* commStgRoot, ZFSInt* commStgRootGlobal, MPI_Comm* commBC2600, ZFSInt* commBC2600Root, ZFSInt* commBC2600RootGlobal,  MPI_Comm* rescalingCommGrComm, ZFSInt* rescalingCommGrRoot, ZFSInt* rescalingCommGrRootGlobal, MPI_Comm* commPerRotOne, MPI_Comm* commPerRotTwo, MPI_Comm* commPerRotWorld, ZFSInt* rotationRoots, ZFSInt& perRotGroup, SingularInformation *singularity,  ZFSInt* hasSingularity); //junoh
  void createWaveWindowMapping(ZFSInt);
  void createCommunicationExchangeFlags(ZFSStrctrdCommunicationHandle* comm, ZFSInt noVariables);
  void createWaveCommunicationExchangeFlags(ZFSStrctrdWaveCommunicationHandle* comm, ZFSInt noVariables);
  vector<ZFSStrctrdWindowMap*> localStrctrdBndryCndMaps;
  vector<ZFSStrctrdWindowMap*> channelSurfaceIndices;
  //for the sponge
  vector<ZFSStrctrdWindowMap*> m_spongeInfoMap;
  vector<ZFSStrctrdWindowMap*> m_wallDistInfoMap;
  //zonal
  vector<ZFSStrctrdWindowMap*> m_zonalBCMaps; //junoh
  ZFSBool checkZonalBCMaps(ZFSStrctrdWindowMap* map1, ZFSStrctrdWindowMap* map2); //junoh


  //=============================================================================================
  //new functions and variables
  void readWindowCoordinates(ZFSInt cpu, ZFSFloat* periodicDisplacements,ZFSString fileName);
  void periodicPointsChange(ZFSFloat* pt,ZFSId type,ZFSFloat* periodicDisplacements); //transform point coords

  //for mutilblock and periodic connections
  ZFSBool addConnection(ZFSInt connectiontype,ZFSInt b1,ZFSInt* p1,ZFSInt b2,ZFSInt* p2);
  ZFSBool findConnection(connectionNode a);
  void removeConnection(connectionNode a);

  //for special connecitons
  ZFSBool addConnection(ZFSInt connectiontype,ZFSInt b1,ZFSInt* p1,ZFSInt b2,ZFSInt* p2,ZFSInt Nstar);
  ZFSBool findConnection(connectionNode a,ZFSInt Nstar);
  void removeConnection(connectionNode a,ZFSInt Nstar);
  void multiBlockAssembling();
  void singularityAssembling();

  set<connectionNode> connectionset;
  set<connectionNode> singularconnectionset;
  vector<ZFSStrctrdWindowMap*> window2d;
  vector<ZFSStrctrdWindowMap*> window1d;
  vector<ZFSStrctrdWindowMap*> window0d;
  vector<ZFSStrctrdWindowMap*> singularwindow;
  vector<ZFSStrctrdWindowMap*> localSingularMap;
  //==============================================================================================


private:
  noBoxWindows* inputBoxWindowInformation;
  noBoxWindows* outputBoxWindowInformation;
  vector<ZFSStrctrdWindowMap*> globalStrctrdBndryCndMaps;
  //vector<ZFSStrctrdWindowMap*> localStrctrdBndryCndMaps;
  vector<ZFSStrctrdWindowMap*> localStrctrdDomainMaps;
  vector<ZFSStrctrdWindowMap*> localStrctrdDomainBndryMaps;

  //diagonal communication
  vector<ZFSStrctrdWindowMap*> m_partitionMapsWithGC;
  vector<ZFSStrctrdWindowMap*> m_partitionMapsWithoutGC;
  vector<ZFSStrctrdWindowMap*> rcvMap;
  vector<ZFSStrctrdWindowMap*> sndMap;
  vector<ZFSStrctrdWindowMap*> rcvMapPeriodic;
  vector<ZFSStrctrdWindowMap*> sndMapPeriodic;
  vector<ZFSStrctrdWindowMap*> physicalBCMap;
  vector<ZFSStrctrdWindowMap*> physicalAuxDataMap;
  vector<ZFSStrctrdWindowMap*> waveRcvMap;
  vector<ZFSStrctrdWindowMap*> waveSndMap;
  vector<windowInformation*> inputWindows;
  ZFSStrctrdDecomposition<nDim>* partition;
  ZFSInt noInputWindowInformation;
  ZFSInt noInputWindowConnections;
  ZFSInt noInputBndryCnds;
  ZFSId file_id;
  ZFSStrctrdWindowMap* m_myMapWithGC;
  ZFSStrctrdWindowMap* m_myMapWithoutGC;
  ZFSInt m_noPartitions;
  ZFSInt m_noGhostLayers;

  // Communication related
  ZFSId noDomains() const { return m_noDomains; }
  ZFSId domainId() const { return m_domainId; }

  MPI_Comm m_zfsStrctrdComm;
  const ZFSId m_noInputBlocks;
  const ZFSId m_inputBlockId;
  const ZFSId m_noDomains;
  const ZFSId m_domainId;

};

#endif
