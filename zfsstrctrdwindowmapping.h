#ifndef ZFSFVBLOCKSTRUCTWINDOWMAPPING
#define ZFSFVBLOCKSTRUCTWINDOWMAPPING

#include "zfstypes.h"

class ZFSStrctrdWindowMap
{
  //friend class ZFSStrctrdBlck;
 public:
  ZFSStrctrdWindowMap(ZFSInt dim);
  ~ZFSStrctrdWindowMap();
  //~ZFSStrctrdWindowMap();
  ZFSId   Id1;
  ZFSInt* start1;
  ZFSInt* end1;
  ZFSInt* step1;
  ZFSId   Id2;
  ZFSInt* start2;
  ZFSInt* end2;
  ZFSInt* step2;
  ZFSInt* order;
  ZFSInt  nDim;
  ZFSInt  BC;
  ZFSInt constIndex;
  ZFSInt face;
  ZFSInt dir;
  ZFSInt dc1;
  ZFSInt dc2;
  ZFSInt originShape;

  //sponge properties:
  ZFSBool hasSponge;
  ZFSFloat spongeThickness;
  ZFSFloat beta;
  ZFSFloat sigma;

  //singularity
  ZFSInt Nstar;
  ZFSInt SingularId;
  ZFSInt SingularBlockId[4];
};

#endif
