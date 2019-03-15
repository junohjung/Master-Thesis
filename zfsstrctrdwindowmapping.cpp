#include "zfsmpi.h"
#include "zfsstrctrdblckwindowinfo.h"
#include "zfsglobals.h"

ZFSStrctrdWindowMap::ZFSStrctrdWindowMap(ZFSInt dim)
{
  start1 = new ZFSInt[dim];
  step1 = new ZFSInt[dim];
  end1 = new ZFSInt[dim];
  start2 = new ZFSInt[dim];
  step2 = new ZFSInt[dim];
  end2 = new ZFSInt[dim];
  order = new ZFSInt[dim];
  nDim=dim;
  BC=-1;
  face=-1;
  hasSponge=false;
  spongeThickness=F0;
  beta=F0;
  sigma=F0;
}

ZFSStrctrdWindowMap::~ZFSStrctrdWindowMap()
{
  delete[] start1;
  start1 = NULL;
  delete[] step1;
  step1 = NULL;
  delete[] end1;
  end1 = NULL;
  delete[] start2;
  start2 = NULL;
  delete[] step2;
  step2 = NULL;
  delete[] end2;
  end2 = NULL;
  delete[] order;
  order = NULL;
}
