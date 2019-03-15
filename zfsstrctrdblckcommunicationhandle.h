#ifndef ZFSFVBLOCKSTRUCTCOMMUNICTATIONHANDLE
#define ZFSFVBLOCKSTRUCTCOMMUNICTATIONHANDLE

#include <vector>
#include "zfstypes.h"
#include "zfsmpi.h"

//This class contains the information needed for the block to communicate via MPI with other processes. It also contains the send and recieve buffers

class ZFSStrctrdCommunicationHandle
{
 public:
  ZFSStrctrdCommunicationHandle(ZFSInt Dimension){noNghbrDomains=0; m_noNghbrDomainCellBufferSizeSnd=NULL;m_noNghbrDomainPointBufferSizeSnd=NULL; m_noNghbrDomainCellBufferSizeRcv=NULL;m_noNghbrDomainPointBufferSizeRcv=NULL; m_bufferCellsSnd=NULL; m_bufferCellsRcv=NULL; m_bufferPointsSnd=NULL; m_bufferPointsRcv=NULL;nDim=Dimension;};
  ~ZFSStrctrdCommunicationHandle();
  void setBufferSizes();
  ZFSInt noNghbrDomains;//includes the periodic boundary condition (all maps)
  ZFSInt noNghbrDomainsNormal; // all excluding periodic boundary condition
  ZFSInt noNghbrDomainsPeriodic; //contains the number of periodic BCs
  ZFSInt noNghbrDomainsChannel; //contains the number of Channel nghbr
  ZFSInt noNghbrDomainsSingular; //singularity
  ZFSInt noNghbrDomainsMultiBlock;
  ZFSInt noNghbrDomainsPeriodicS;
  //for the wave exchange
  ZFSInt noNghbrDomainsRcv;
  ZFSInt noNghbrDomainsSnd;

  ZFSInt*      m_tagHelperRCV;
  ZFSInt*      m_tagHelperSND;
  ZFSInt*      m_sndNghbrId;
  ZFSInt*      m_rcvNghbrId;
  ZFSInt*      m_noNghbrDomainCellBufferSizeSnd;
  ZFSInt*      m_noNghbrDomainPointBufferSizeSnd;
  ZFSInt*      m_noNghbrDomainCellBufferSizeRcv;
  ZFSInt*      m_noNghbrDomainPointBufferSizeRcv;

  ZFSFloat**   m_bufferCellsSnd; //comunicator cells==>cells to be communicated (send and received)
  ZFSFloat**   m_bufferCellsRcv;
  ZFSFloat**   m_bufferPointsSnd;
  ZFSFloat**   m_bufferPointsRcv; //comunicator points ==> points to be communicated (send and received)
  ZFSInt*      m_nghbrFaceId;
  ZFSInt**     m_nghbrFaceInfo;
  //ZFSInt*      m_sndRcvBufferPos;
  MPI_Request* mpi_sndRequest;
  MPI_Request* mpi_rcvRequest;
  MPI_Request* mpi_sndRequestCells;
  MPI_Request* mpi_rcvRequestCells;
  MPI_Status*  mpi_sndStatus;
  MPI_Status*  mpi_rcvStatus;
  ZFSInt**     startInfoSNDcells;
  ZFSInt**     startInfoRCVcells;
  ZFSInt**     endInfoSNDcells;
  ZFSInt**     endInfoRCVcells;
  ZFSInt**     startInfoSNDpoints;
  ZFSInt**     startInfoRCVpoints;
  ZFSInt**     startInfo2;
  ZFSInt**     startInfo1;
  ZFSInt**     endInfo1;
  ZFSInt**     endInfoSNDpoints;
  ZFSInt**     endInfoRCVpoints;
  ZFSInt**     orderInfo;
  ZFSInt**     stepInfoSND;
  ZFSInt**     stepInfoRCV;
  ZFSInt*      constIndex;
  ZFSInt*      constIndexDir;
  ZFSInt       nDim;
  ZFSInt*      bcId;

  //singular communications
  ZFSInt*      singularID;
};

class ZFSStrctrdWaveCommunicationHandle
{
 public:
  ZFSStrctrdWaveCommunicationHandle(ZFSInt Dimension){noNghbrDomainsSnd=0; noNghbrDomainsRcv=0; m_noNghbrDomainCellBufferSizeSnd=NULL; m_noNghbrDomainCellBufferSizeRcv=NULL; m_bufferCellsSnd=NULL; m_bufferCellsRcv=NULL; nDim=Dimension;};
  ~ZFSStrctrdWaveCommunicationHandle();
  void setBufferSizes();

  ZFSInt nDim;
  ZFSInt noNghbrDomainsRcv;
  ZFSInt noNghbrDomainsSnd;

  ZFSInt*      m_tagHelperRCV;
  ZFSInt*      m_tagHelperSND;
  ZFSInt*      m_sndNghbrId;
  ZFSInt*      m_rcvNghbrId;
  ZFSInt*      m_noNghbrDomainCellBufferSizeSnd;
  ZFSInt*      m_noNghbrDomainCellBufferSizeRcv;

  ZFSFloat**   m_bufferCellsSnd; //comunicator cells==>cells to be communicated (send and received)
  ZFSFloat**   m_bufferCellsRcv;
  MPI_Request* mpi_sndRequest;
  MPI_Request* mpi_rcvRequest;
  MPI_Request* mpi_sndRequestCells;
  MPI_Request* mpi_rcvRequestCells;
  MPI_Status*  mpi_sndStatus;
  MPI_Status*  mpi_rcvStatus;
  ZFSInt**     startInfoSNDcells;
  ZFSInt**     startInfoRCVcells;
  ZFSInt**     endInfoSNDcells;
  ZFSInt**     endInfoRCVcells;
};



#endif








