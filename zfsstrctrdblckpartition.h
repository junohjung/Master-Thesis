#ifndef ZFSFVBLOCKSTRUCTPARTITION_H
#define ZFSFVBLOCKSTRUCTPARTITION_H

#include "zfstypes.h"
#include "zfsscratch.h"

class nodeType
{
 public:
  nodeType(){strctrdLevel=0;NoStrctrdChild=1;NoLeaf=1;};
  ~nodeType();
  ZFSInt strctrdLevel;
  ZFSInt NoStrctrdChild;
  ZFSInt NoLeaf;
  nodeType** strctrdchild;
};

template <ZFSInt nDim>
class BoxInfoType
{
 public:
  BoxInfoType();
  ~BoxInfoType();
  ZFSInt inputBoxID;
  ZFSInt outputBoxID;
  ZFSInt cpu;
  ZFSInt* DirLast;
  ZFSInt* offset;
  ZFSInt totalPhysicalSize;
};

template <ZFSInt nDim>
class ZFSStrctrdDecomposition
{
  template <ZFSInt nDim_> friend class ZFSStrctrdBlck;
  friend class ZFSStrctrdBlck3D;
  friend class ZFSStrctrdBlck2D;
  template <ZFSInt nDim_> friend class ZFSStrctrdBlckWindowInfo;
 public:
  ZFSStrctrdDecomposition(ZFSInt noBox, ZFSString fileName, const ZFSId noDomains_);
  ZFSStrctrdDecomposition(ZFSInt noBox, ZFSId file_id, const ZFSId noDomains_);
  ~ZFSStrctrdDecomposition();
  ZFSInt rounding(ZFSFloat x);
  void initializeChilds(nodeType* &treePointer);
  void setOutputBoxInfoHelper(nodeType* &treePointer, ZFSInt inputBoxID, ZFSIntScratchSpace &level2dimension, ZFSIntScratchSpace &level2dimension1D, ZFSInt &countOutputBox, ZFSIntScratchSpace &beginBoxHelper, ZFSIntScratchSpace &endBoxHelper, ZFSIntScratchSpace &divisor);
  void setOutputBoxInfo(nodeType** &treeRoot, ZFSIntScratchSpace &level2dimension);
  ZFSInt quickSortPartition(ZFSInt left, ZFSInt right, ZFSIntScratchSpace &sizeSortedOutputBoxID);
  void quickSort(ZFSIntScratchSpace &sizeSortedOutputBoxID, ZFSInt low, ZFSInt high);
  ZFSBool loadBalancing();
  ZFSBool isBoxCountExceeded();
  void diveForInsertionPoint(nodeType* &treePointer, ZFSInt inputBoxID, ZFSIntScratchSpace &level2dimension1D, nodeType* &insertTreePointer);
  void destroyChilds(nodeType* &treePointer);
  ZFSInt sumLeaves(nodeType* &treePointer);
  void insertChild(nodeType* &insertTreePointer, ZFSInt inputBoxID, ZFSIntScratchSpace &level2dimension1D);
  void addLeaf(nodeType** &treeRoot, ZFSIntScratchSpace &level2dimension);
  ZFSInt decompose();
  ZFSInt readFromFile();
 private:
  ZFSId noDomains() const { return m_noDomains; }

  ZFSFloat m_eps;
  const ZFSId m_noDomains;
  ZFSInt m_oldNoOutputBoxes;
  ZFSInt m_noOutputBoxes;
  ZFSInt m_noInputBoxes;
  BoxInfoType<nDim>** inputBoxInfo;
  BoxInfoType<nDim>** outputBoxInfo;
};

#endif
