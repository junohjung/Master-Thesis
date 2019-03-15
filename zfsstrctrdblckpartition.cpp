#include <cstring>
#include "zfsglobals.h"
#include "zfsmpi.h"
#include "zfsiolib.h"
#include "zfsstrctrdblckpartition.h"
#include <cstring>

using namespace std;


nodeType::~nodeType()
{
}

template <ZFSInt nDim>
ZFSInt ZFSStrctrdDecomposition<nDim>::rounding(ZFSFloat x)
{
  ZFSFloat a;
  ((x - floor(x)) >= 0.5) ? a=ceil(x) : a=floor(x);
  return (ZFSInt)a;
}


template <ZFSInt nDim>
BoxInfoType<nDim>::BoxInfoType()
{
  DirLast = new ZFSInt[nDim];
  offset = new ZFSInt[nDim];
}

template <ZFSInt nDim>
BoxInfoType<nDim>::~BoxInfoType()
{
  delete [] DirLast;
  delete [] offset;
}

const ZFSInt m_maxRatioNoOutputBoxNoDomains = 1;
const ZFSFloat m_maxRelativeSizeDeviationFromAverage = 0.05;

template <ZFSInt nDim>
ZFSStrctrdDecomposition<nDim>::ZFSStrctrdDecomposition(ZFSInt noBox,
                                                       ZFSString fileName,
                                                       const ZFSId noDomains_):
  m_eps(std::numeric_limits<ZFSFloat>::epsilon()), m_noDomains(noDomains_)
{
  m_oldNoOutputBoxes = 0;
  m_noOutputBoxes=1;
  m_noInputBoxes=noBox;
  inputBoxInfo=NULL;
  outputBoxInfo=NULL;
  inputBoxInfo = new BoxInfoType<nDim> *[m_noInputBoxes];
  for(int j=0;j<m_noInputBoxes; j++) {
    inputBoxInfo[j]=new BoxInfoType<nDim>();
  }
  ZFSIntScratchSpace dataSetSize(m_noInputBoxes*nDim, __CALLING_FUNCTION__, "datasetsize");
  dataSetSize.fill(0);
  if(globalDomainId()==0){
    int file_id=io_openfile("hdf5", fileName.c_str(), "collective", MPI_COMM_SELF);
    for(ZFSId i=0; i<m_noInputBoxes; i++){
      const char* cBlockName;
      ZFSString sBlockName = "/block";
      stringstream dummy1;
      dummy1 << i;
      sBlockName += dummy1.str();
      cBlockName = sBlockName.c_str();
      io_getDatasetSize(file_id , cBlockName,"x", nDim , &(dataSetSize[i*nDim]));
    }
    io_closefile(file_id);
    MPI_Bcast(dataSetSize.getPointer(),m_noInputBoxes*nDim, MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    MPI_Bcast(dataSetSize.getPointer(),m_noInputBoxes*nDim, MPI_INT, 0, MPI_COMM_WORLD);
  }

  for(ZFSId i=0; i<m_noInputBoxes; i++){
    inputBoxInfo[i]->inputBoxID = i;
    for(ZFSId j=0; j<nDim; j++){
      inputBoxInfo[i]->DirLast[j]=dataSetSize[i*nDim+j]-1;
      inputBoxInfo[i]->offset[j]=0;
    }
  }
}

template <ZFSInt nDim>
ZFSStrctrdDecomposition<nDim>::~ZFSStrctrdDecomposition() {

  for(int j=0;j<m_noInputBoxes; j++) {
    delete inputBoxInfo[j];
  }

  delete [] inputBoxInfo;

  for(int j=0;j<m_noOutputBoxes; j++) {
    delete outputBoxInfo[j];
  }

  delete [] outputBoxInfo;
}

template <ZFSInt nDim>
ZFSStrctrdDecomposition<nDim>::ZFSStrctrdDecomposition(ZFSInt noBox,
                                                       ZFSId file_id,
                                                       const ZFSId noDomains_):
  m_eps(std::numeric_limits<ZFSFloat>::epsilon()), m_noDomains(noDomains_)
{
  m_oldNoOutputBoxes = 0;
  m_noOutputBoxes=1;
  m_noInputBoxes=noBox;
  inputBoxInfo=NULL;
  outputBoxInfo=NULL;
  inputBoxInfo = new BoxInfoType<nDim> *[m_noInputBoxes];
  for(int j=0;j<m_noInputBoxes; j++) {
    inputBoxInfo[j]=new BoxInfoType<nDim>();
  }
  ZFSInt* address_ijk_max;
  ZFSInt ijk_max[3];
  address_ijk_max =ijk_max;
  for(ZFSId i=0; i<m_noInputBoxes; i++) {
    inputBoxInfo[i]->inputBoxID = i;
    //create string containing the block number
    const char* cBlockName;

    ZFSString sBlockName = "/block";
    stringstream dummy1;
    //need to add one because of Fortran notation
    dummy1 << i;
    sBlockName += dummy1.str();
    cBlockName = sBlockName.c_str();
    //get ijk_max
    io_getDatasetSize(file_id , cBlockName,"x", nDim , address_ijk_max );
    //assign ijk_max to inputBoxInfo and also assign the offset
    for(ZFSId j=0; j<nDim; j++) {
      inputBoxInfo[i]->DirLast[j]=ijk_max[j]-1;
      inputBoxInfo[i]->offset[j]=0;
    }
  }
}


template <ZFSInt nDim>
void  ZFSStrctrdDecomposition<nDim>::initializeChilds(nodeType* &treePointer)
{
  nodeType* myTreePointer;
  ZFSInt countChild;
  if (treePointer->strctrdLevel < (nDim-1)) {
    treePointer->strctrdchild = new nodeType *[treePointer->NoStrctrdChild];
    for(ZFSInt i=0; i<treePointer->NoStrctrdChild; i++) {
      treePointer->strctrdchild[i]= new nodeType;
    }
    for(countChild=0;countChild<treePointer->NoStrctrdChild; countChild++) {
      treePointer->strctrdchild[countChild]->strctrdLevel = (treePointer->strctrdLevel)+1;
      myTreePointer = treePointer->strctrdchild[countChild];
      initializeChilds(myTreePointer);
    }
  }
}

template <ZFSInt nDim>
void  ZFSStrctrdDecomposition<nDim>::setOutputBoxInfoHelper(nodeType* &treePointer,
                                                            ZFSInt inputBoxID,
                                                            ZFSIntScratchSpace &level2dimension,
                                                            ZFSIntScratchSpace &level2dimension1D,
                                                            ZFSInt &countOutputBox,
                                                            ZFSIntScratchSpace &beginBoxHelper,
                                                            ZFSIntScratchSpace &endBoxHelper,
                                                            ZFSIntScratchSpace &divisor){

  nodeType* myTreePointer;
  ZFSIntScratchSpace beginBox(nDim, __CALLING_FUNCTION__, "beginBox");
  ZFSIntScratchSpace endBox(nDim, __CALLING_FUNCTION__, "endBox");

  beginBox.fill(0);
  endBox.fill(0);

  divisor[level2dimension1D[treePointer->strctrdLevel]] = treePointer->NoLeaf;

  if (treePointer->strctrdLevel<(nDim-1)) {
    beginBoxHelper[level2dimension1D[treePointer->strctrdLevel]] =0;
    for (ZFSInt i=0;i<treePointer->NoStrctrdChild;i++) {
      if(i>0) {
        beginBoxHelper[level2dimension1D[treePointer->strctrdLevel]]= beginBoxHelper[level2dimension1D[treePointer->strctrdLevel]]+ treePointer->strctrdchild[i-1]->NoLeaf;
      }
      endBoxHelper[level2dimension1D[treePointer->strctrdLevel]] = beginBoxHelper[level2dimension1D[treePointer->strctrdLevel]] + treePointer->strctrdchild[i]->NoLeaf;
      myTreePointer= treePointer->strctrdchild[i];
      setOutputBoxInfoHelper(myTreePointer, inputBoxID, level2dimension, level2dimension1D, countOutputBox, beginBoxHelper, endBoxHelper, divisor);
    }
  } else {
    for (ZFSInt i=0;i<treePointer->NoStrctrdChild;i++) {
      beginBoxHelper[level2dimension1D[treePointer->strctrdLevel]]=i;
      endBoxHelper[level2dimension1D[treePointer->strctrdLevel]]=i+1;
      for(ZFSInt j=0;j<nDim;j++) {
        beginBox[j] =rounding( ((ZFSFloat)inputBoxInfo[inputBoxID]->DirLast[j])*((ZFSFloat)beginBoxHelper[j])/((ZFSFloat)divisor[j]));
        endBox[j] =rounding( ((ZFSFloat)inputBoxInfo[inputBoxID]->DirLast[j])*((ZFSFloat)endBoxHelper[j])/((ZFSFloat)divisor[j]));
      }

      outputBoxInfo[countOutputBox]->totalPhysicalSize=1;
      for (ZFSInt j=0;j<nDim;j++) {
        outputBoxInfo[countOutputBox]->DirLast[j]=endBox[j]-beginBox[j]+1;
        outputBoxInfo[countOutputBox]->totalPhysicalSize *= (outputBoxInfo[countOutputBox]->DirLast[j]);
        outputBoxInfo[countOutputBox]->offset[j]=beginBox[j];
      }

      outputBoxInfo[countOutputBox]->inputBoxID = inputBoxID;
      outputBoxInfo[countOutputBox]->outputBoxID = countOutputBox;
      countOutputBox++;
    }
  }
}

template <ZFSInt nDim>
void  ZFSStrctrdDecomposition<nDim>::setOutputBoxInfo(nodeType** &treeRoot, ZFSIntScratchSpace &level2dimension){
  nodeType* myTreeRootPointer;
  ZFSIntScratchSpace beginBoxHelper(nDim, __CALLING_FUNCTION__, "beginBoxHelper");
  ZFSIntScratchSpace endBoxHelper(nDim, __CALLING_FUNCTION__, "endBoxHelper");
  ZFSIntScratchSpace divisor(nDim, __CALLING_FUNCTION__, "divisor");
  ZFSInt countOutputBox;

  beginBoxHelper.fill(0);
  endBoxHelper.fill(0);
  divisor.fill(0);

  for(int i=0;i<m_oldNoOutputBoxes;i++) {
    delete outputBoxInfo[i];
  }

  delete[] outputBoxInfo;
  outputBoxInfo = new BoxInfoType<nDim> *[m_noOutputBoxes];

  for(ZFSInt j=0;j<m_noOutputBoxes; j++){
    outputBoxInfo[j]=new BoxInfoType<nDim>();
  }

  m_oldNoOutputBoxes = m_noOutputBoxes;
  countOutputBox = 0;

  for(ZFSInt i=0;i<m_noInputBoxes;i++) {
    myTreeRootPointer=treeRoot[i];
    ZFSIntScratchSpace level2dimension1D(nDim, __CALLING_FUNCTION__, "level2dimension1D");

    for(ZFSInt j=0;j<nDim;j++) {
      level2dimension1D[j]=level2dimension(j,i);
    }

    setOutputBoxInfoHelper(myTreeRootPointer, i, level2dimension, level2dimension1D, countOutputBox, beginBoxHelper, endBoxHelper, divisor);
  }
}

template <ZFSInt nDim>
ZFSInt  ZFSStrctrdDecomposition<nDim>::quickSortPartition(ZFSInt left, ZFSInt right,
                                                          ZFSIntScratchSpace &sizeSortedOutputBoxID)
{
  ZFSInt myLeft, myRight, tmp, pivot;
  myLeft = left;
  myRight = right-1 ;
  pivot = outputBoxInfo[sizeSortedOutputBoxID[right]]->totalPhysicalSize;

  do {
    while(outputBoxInfo[sizeSortedOutputBoxID[myLeft]]->totalPhysicalSize >= pivot && myLeft < right) {
      myLeft++;
    }

    while(outputBoxInfo[sizeSortedOutputBoxID[myRight]]->totalPhysicalSize <=pivot && myRight > left) {
      myRight--;
    }

    if(myLeft < myRight) {
      tmp=sizeSortedOutputBoxID[myLeft];
      sizeSortedOutputBoxID[myLeft]=sizeSortedOutputBoxID[myRight];
      sizeSortedOutputBoxID[myRight]= tmp;
    } else {
      break;
    }

  } while(myLeft < myRight);

  if(outputBoxInfo[sizeSortedOutputBoxID[myLeft]]->totalPhysicalSize < pivot) {
    tmp = sizeSortedOutputBoxID[myLeft];
    sizeSortedOutputBoxID[myLeft] = sizeSortedOutputBoxID[right];
    sizeSortedOutputBoxID[right] = tmp;
  }

  return myLeft;
}

//new QuickSort

template <ZFSInt nDim>
void  ZFSStrctrdDecomposition<nDim>::quickSort(ZFSIntScratchSpace &sizeSortedOutputBoxID,
                                               ZFSInt low, ZFSInt high)
{
  ZFSInt index;

  if(low<high) {
    index = quickSortPartition(low, high, sizeSortedOutputBoxID);
    quickSort(sizeSortedOutputBoxID, low, index-1);
    quickSort(sizeSortedOutputBoxID, index+1, high);
  }
}

template <ZFSInt nDim>
ZFSBool  ZFSStrctrdDecomposition<nDim>::loadBalancing(){
  ZFSIntScratchSpace nodeCount(noDomains(), __CALLING_FUNCTION__, "nodeCount");
  ZFSInt countCPU, CPUWithLeastNodes;
  ZFSFloat averageNodeCount;

  if (m_noOutputBoxes<noDomains()) {
    return false;
  } else {
    ZFSIntScratchSpace sizeSortedOutputBoxID(m_noOutputBoxes,__CALLING_FUNCTION__, "sizeSortedOutputBoxID");
    for(ZFSInt i=0;i<m_noOutputBoxes;i++) {
      sizeSortedOutputBoxID[i]=i;
    }
    quickSort(sizeSortedOutputBoxID, 0, m_noOutputBoxes-1);
    for (countCPU=0; countCPU<noDomains(); countCPU++) {
      outputBoxInfo[sizeSortedOutputBoxID[countCPU]]->cpu = countCPU;
      nodeCount[countCPU]= outputBoxInfo[sizeSortedOutputBoxID[countCPU]]->totalPhysicalSize;
    }

    ZFSIntScratchSpace dummy(noDomains(), __CALLING_FUNCTION__, "dummy");
    ZFSIntScratchSpace minpos(noDomains(), __CALLING_FUNCTION__, "minpos");
    ZFSInt temp;

    for(ZFSInt j=noDomains();j<m_noOutputBoxes;j++) {
      //find the smalles index
      //copy for sorting algorithm
      for (ZFSInt i=0;i<noDomains();i++) {
        dummy[i]=nodeCount[i];
        minpos[i]=i;
      }

      //sort array and index
      for(ZFSInt i = 1; (i <= noDomains()); i++) {
        for (ZFSInt k=0; k < (noDomains() -1); k++) {
          // ascending order simply changes to <
          if (dummy[k+1] < dummy[k]) {
            // swap elements
            temp = dummy[k];
            dummy[k] = dummy[k+1];
            dummy[k+1] = temp;
            temp=minpos[k];
            minpos[k]=minpos[k+1];
            minpos[k+1]=temp;
          }
        }
      }

      CPUWithLeastNodes = minpos[0];
      outputBoxInfo[sizeSortedOutputBoxID[j]]->cpu =CPUWithLeastNodes;
      nodeCount[CPUWithLeastNodes] = nodeCount[CPUWithLeastNodes] + outputBoxInfo[sizeSortedOutputBoxID[j]]->totalPhysicalSize;
    }

    ZFSLong sum=0;
    for(ZFSInt i=0; i<noDomains(); i++){
      sum +=nodeCount[i];
    }

    averageNodeCount = ((ZFSFloat)sum)/((ZFSFloat)noDomains());
    ZFSInt min=nodeCount[0], max=nodeCount[0];

    for(ZFSInt i=0;i<noDomains();i++) {
      if(nodeCount[i]>max) {
        max=nodeCount[i];
      }
      if(nodeCount[i]<min) {
        min=nodeCount[i];
      }
    }

    if (abs(((ZFSFloat)max-averageNodeCount)/averageNodeCount) <= m_maxRelativeSizeDeviationFromAverage &&
        abs(((ZFSFloat)min-averageNodeCount)/averageNodeCount) <= m_maxRelativeSizeDeviationFromAverage) {
      return true;
    } else {
      return false;
    }
  }
}

template <ZFSInt nDim>
ZFSBool ZFSStrctrdDecomposition<nDim>:: isBoxCountExceeded(){
  if(((ZFSFloat)(m_noOutputBoxes+1))/((ZFSFloat)noDomains())>m_maxRatioNoOutputBoxNoDomains) {
    return true;
  } else {
    return false;
  }
}

template <ZFSInt nDim>
void  ZFSStrctrdDecomposition<nDim>::diveForInsertionPoint(nodeType* &treePointer, ZFSInt inputBoxID, ZFSIntScratchSpace &level2dimension1D, nodeType* &insertTreePointer)
{
  nodeType* myTreePointer;



  if ((((ZFSFloat)treePointer->NoStrctrdChild)/((ZFSFloat)inputBoxInfo[inputBoxID]->DirLast[level2dimension1D[treePointer->strctrdLevel]]+1.0)) < (((ZFSFloat)insertTreePointer->NoStrctrdChild)/((ZFSFloat)inputBoxInfo[inputBoxID]->DirLast[level2dimension1D[insertTreePointer->strctrdLevel]]+1.0)) || (approx(((ZFSFloat)(treePointer->NoStrctrdChild) / ((ZFSFloat)(inputBoxInfo[inputBoxID]->DirLast[level2dimension1D[treePointer->strctrdLevel]])+1.0)), ((ZFSFloat)(insertTreePointer->NoStrctrdChild) / ((ZFSFloat)(inputBoxInfo[inputBoxID]->DirLast[level2dimension1D[insertTreePointer->strctrdLevel]])+1.0)), m_eps) && (treePointer->strctrdLevel < insertTreePointer->strctrdLevel))) {
    insertTreePointer = treePointer;
  }

  if (treePointer->strctrdLevel<(nDim-1)) {
    myTreePointer =  treePointer->strctrdchild[0];
    for (ZFSInt i=1;i<treePointer->NoStrctrdChild;i++) {
      if (treePointer->strctrdchild[i]->NoLeaf < myTreePointer->NoLeaf) {
        myTreePointer= treePointer->strctrdchild[i];
      }
    }
    diveForInsertionPoint(myTreePointer, inputBoxID, level2dimension1D, insertTreePointer);
  }
}

template <ZFSInt nDim>
void  ZFSStrctrdDecomposition<nDim>::destroyChilds(nodeType* &treePointer)
{
  nodeType* myTreePointer;
  if (treePointer->strctrdLevel<(nDim-1)) {
    for (ZFSInt i=0; i<treePointer->NoStrctrdChild;i++) {
      myTreePointer = treePointer->strctrdchild[i];
      destroyChilds(myTreePointer);
    }
    delete[] treePointer->strctrdchild;
  }
}

template <ZFSInt nDim>
ZFSInt  ZFSStrctrdDecomposition<nDim>::sumLeaves(nodeType* &treePointer)
{
  nodeType* myTreePointer;
  ZFSInt sumLeaf;

  if (treePointer->strctrdLevel<(nDim-1)) {
    sumLeaf=0;
    for (ZFSInt i=0;i<treePointer->NoStrctrdChild;i++) {
      myTreePointer = treePointer->strctrdchild[i];
      sumLeaf += sumLeaves(myTreePointer);
    }
  } else {
    sumLeaf = treePointer->NoStrctrdChild;
  }
  treePointer->NoLeaf = sumLeaf;
  return sumLeaf;
}


template <ZFSInt nDim>
void  ZFSStrctrdDecomposition<nDim>::insertChild(nodeType* &insertTreePointer, ZFSInt inputBoxID, ZFSIntScratchSpace &level2dimension1D)
{
  nodeType* myInsertTreePointer;
  ZFSInt myNoLeaf;

  destroyChilds(insertTreePointer);
  insertTreePointer->NoStrctrdChild++;
  myNoLeaf = insertTreePointer->NoLeaf+1;
  initializeChilds(insertTreePointer);

  if (insertTreePointer->strctrdLevel <(nDim-1)) {
    for(ZFSInt i=insertTreePointer->NoStrctrdChild+1; i<=myNoLeaf; i++) {
      myInsertTreePointer=insertTreePointer;
      diveForInsertionPoint(insertTreePointer, inputBoxID, level2dimension1D, myInsertTreePointer);
      insertChild(myInsertTreePointer, inputBoxID, level2dimension1D);
      insertTreePointer->NoLeaf = sumLeaves(insertTreePointer);
    }
  }
}

template <ZFSInt nDim>
void  ZFSStrctrdDecomposition<nDim>::addLeaf(nodeType** &treeRoot, ZFSIntScratchSpace &level2dimension)
{
  nodeType* myTreeRootPointer;
  nodeType* myInsertTreePointer;
  ZFSInt inputBoxID = 0;
  ZFSInt maxTotalPhysicalSize=0;

  for(ZFSInt i=0;i<m_noOutputBoxes;i++) {
    if(outputBoxInfo[i]->totalPhysicalSize>maxTotalPhysicalSize) {
      inputBoxID = outputBoxInfo[i]->inputBoxID;
      maxTotalPhysicalSize = outputBoxInfo[i]->totalPhysicalSize;
    }
  }

  myTreeRootPointer = treeRoot[inputBoxID];
  myInsertTreePointer = treeRoot[inputBoxID];
  ZFSIntScratchSpace level2dimension1D(nDim, __CALLING_FUNCTION__, "level2dimension1D");

  for(ZFSInt i=0;i<nDim;i++) {
    level2dimension1D[i]=level2dimension(i,inputBoxID);
  }

  diveForInsertionPoint(myTreeRootPointer, inputBoxID, level2dimension1D, myInsertTreePointer);
  insertChild(myInsertTreePointer, inputBoxID, level2dimension1D);
  myTreeRootPointer->NoLeaf = sumLeaves(myTreeRootPointer);
  m_noOutputBoxes++;;
  setOutputBoxInfo(treeRoot, level2dimension);
}

template <ZFSInt nDim>
ZFSInt  ZFSStrctrdDecomposition<nDim>::decompose(){
  nodeType* myTreeRootPointer;
  nodeType** myTreeRoot;
  ZFSInt countInputBox;
  m_noOutputBoxes = m_noInputBoxes;
  myTreeRoot = new nodeType *[m_noInputBoxes];

  //need to call the constructor for each pointer
  for(ZFSInt i=0;i<m_noInputBoxes; i++) {
    myTreeRoot[i]=new nodeType;
  }
  for(countInputBox=0;countInputBox<m_noInputBoxes; countInputBox++) {
    myTreeRootPointer = myTreeRoot[countInputBox];
    initializeChilds(myTreeRootPointer);
  }

  ZFSIntScratchSpace level2dimensionA(nDim, m_noInputBoxes, __CALLING_FUNCTION__, "level2dimensionA");
  level2dimensionA.fill(0);
  ZFSIntScratchSpace dummy(nDim, __CALLING_FUNCTION__, "dummy");
  ZFSInt maxpos[MAX_SPACE_DIMENSIONS];

  //copy for sorting algorithm
  for(ZFSInt k=0;k<m_noInputBoxes;k++) {
    for (ZFSInt i=0;i<nDim;i++) {
      dummy[i]=inputBoxInfo[k]->DirLast[i];
      maxpos[i]=i;
    }

    ZFSInt temp;
    //sort array and index
    for(ZFSInt i = 1; (i <= nDim); i++) {
      for (ZFSInt j=0; j < (nDim -1); j++) {
        // ascending order simply changes to <
        if (dummy[j+1] > dummy[j]) {
          // swap elements
          temp = dummy[j];
          dummy[j] = dummy[j+1];
          dummy[j+1] = temp;
          temp=maxpos[j];
          maxpos[j]=maxpos[j+1];
          maxpos[j+1]=temp;
        }
      }
    }

    //copy back the sorted index
    for (ZFSInt i=0;i<nDim;i++) {
      level2dimensionA(i,k)=maxpos[i];
    }
  }

  setOutputBoxInfo(myTreeRoot, level2dimensionA);
  ZFSBool loadBalanced = loadBalancing();
  ZFSBool boxCountExceeded = isBoxCountExceeded();

  while(!loadBalanced&&!boxCountExceeded) {
    addLeaf(myTreeRoot, level2dimensionA);
    loadBalanced = loadBalancing();
    boxCountExceeded = isBoxCountExceeded();
  }

  for (ZFSInt i=0;i<m_noInputBoxes;i++) {
    myTreeRootPointer = myTreeRoot[i];
    destroyChilds(myTreeRootPointer);
  }


  for(int j=0;j<m_noInputBoxes; j++) {
    delete myTreeRoot[j];
  }

  delete[] myTreeRoot;

  return 0;
}

template <ZFSInt nDim>
ZFSInt ZFSStrctrdDecomposition<nDim>::readFromFile(){
  //temporary scratch to send/receive the partition Info
  ZFSIntScratchSpace partitionInfo((3+3+4)*noDomains(), "partitionInfo", __CALLING_FUNCTION__);
  partitionInfo.fill(0);
  ZFSId noCPUs=0;
  if(globalDomainId()==0){ //only root reads in the information
    //open the file and read information out of it
    ZFSFloat dumptime = MPI_Wtime();
    zfs_log << "opening partition file " << dumptime << endl;
    ZFSInt fileId=io_openfile("hdf5", "partitionFile.hdf5", "collective",MPI_COMM_SELF);
    dumptime = MPI_Wtime();
    zfs_log << "reading attribute " << dumptime << endl;
    io_read_iattribute1(fileId,"", "noCPUs", &noCPUs);
    dumptime = MPI_Wtime();
    zfs_log << "before bcase " << dumptime << endl;
    MPI_Bcast(&noCPUs, 1, MPI_INT, 0, MPI_COMM_WORLD); //broadcast the number of cpus
    if(noCPUs!=noDomains()){
      zfs_log << "WARNING no CPUs is not equal to ones given in partitionFile" << endl;
      return 0;
    }
    m_noOutputBoxes=noCPUs;
    //allocate the number of boxes
    outputBoxInfo = new BoxInfoType<nDim> *[noCPUs];
    //read the information out of the file
    ZFSInt asize[2]={noCPUs,10};
    dumptime = MPI_Wtime();
    zfs_log << "before reading data " << dumptime << endl;
    io_read_idataset1(fileId,"","partitionData", 2, asize, partitionInfo.getPointer());

    dumptime = MPI_Wtime();
    zfs_log << "before closing " << dumptime << endl;

    io_closefile(fileId);
    dumptime = MPI_Wtime();
    zfs_log << "after closing " << dumptime << endl;

    MPI_Bcast(partitionInfo.getPointer(), noCPUs*(10), MPI_INT, 0, MPI_COMM_WORLD );
  } else{ // all other processes
    //get the number of cpus
    MPI_Bcast(&noCPUs, 1, MPI_INT, 0, MPI_COMM_WORLD); //broadcast the number of cpus
    if(noCPUs!=noDomains()){
      zfs_log << "WARNING no CPUs is not equal to ones given in partitionFile" << endl;
      return 0;
    }
    m_noOutputBoxes=noCPUs;
    //allocate the number of boxes
    outputBoxInfo = new BoxInfoType<nDim> *[noCPUs];
    MPI_Bcast(partitionInfo.getPointer(), noCPUs*(10), MPI_INT, 0, MPI_COMM_WORLD );
  }
  ZFSFloat dumptime;
  MPI_Barrier(MPI_COMM_WORLD);
  dumptime = MPI_Wtime();
  zfs_log << "before distributing " << dumptime << endl;

  //put the information into the right place
  for(ZFSInt j=0;j<m_noOutputBoxes; j++){
    outputBoxInfo[j]=new BoxInfoType<nDim>();
    memcpy(outputBoxInfo[j]->DirLast, &(partitionInfo[j*10]),3*sizeof(ZFSInt));
    memcpy(outputBoxInfo[j]->offset, &(partitionInfo[j*10+3]),3*sizeof(ZFSInt));
    memcpy(&(outputBoxInfo[j]->inputBoxID), &(partitionInfo[j*10+6]), sizeof(ZFSInt));
    memcpy(&(outputBoxInfo[j]->outputBoxID), &(partitionInfo[j*10+7]), sizeof(ZFSInt));
    memcpy(&(outputBoxInfo[j]->cpu), &(partitionInfo[j*10+8]), sizeof(ZFSInt));
    memcpy(&(outputBoxInfo[j]->totalPhysicalSize), &(partitionInfo[j*10+9]), sizeof(ZFSInt));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  dumptime = MPI_Wtime();
  zfs_log << " fter distribution" << dumptime << endl;

  return 1;
}

// Explicit instantiations for 2D and 3D
template class ZFSStrctrdDecomposition<2>;
template class ZFSStrctrdDecomposition<3>;
