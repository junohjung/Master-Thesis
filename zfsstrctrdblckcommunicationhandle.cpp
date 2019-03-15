#include "zfsglobals.h"
#include "zfstypes.h"
#include "zfsstrctrdblckcommunicationhandle.h"


void ZFSStrctrdCommunicationHandle::setBufferSizes()
{
  m_bufferCellsSnd= new ZFSFloat *[noNghbrDomains];
  m_bufferPointsSnd= new ZFSFloat *[noNghbrDomains];
  m_bufferCellsRcv= new ZFSFloat *[noNghbrDomains];
  m_bufferPointsRcv= new ZFSFloat *[noNghbrDomains];

  for(ZFSInt i=0; i<noNghbrDomains; i++) {
    m_bufferCellsSnd[i]  = new ZFSFloat[m_noNghbrDomainCellBufferSizeSnd[i]];
    m_bufferPointsSnd[i] =  new ZFSFloat[m_noNghbrDomainPointBufferSizeSnd[i]];
    m_bufferCellsRcv[i]  =  new ZFSFloat[m_noNghbrDomainCellBufferSizeRcv[i]];
    m_bufferPointsRcv[i] = new ZFSFloat[m_noNghbrDomainPointBufferSizeRcv[i]];
  }

  mpi_sndRequest= new  MPI_Request[noNghbrDomains];
  mpi_rcvRequest= new MPI_Request[noNghbrDomains];
  mpi_sndRequestCells= new  MPI_Request[noNghbrDomains];
  mpi_rcvRequestCells= new MPI_Request[noNghbrDomains];
  mpi_rcvStatus= new MPI_Status[noNghbrDomains];
  mpi_sndStatus= new MPI_Status[noNghbrDomains];

  for(ZFSId i=0; i<noNghbrDomains; i++){
    mpi_sndRequest[i]=MPI_REQUEST_NULL;
    mpi_rcvRequest[i]=MPI_REQUEST_NULL;
    mpi_sndRequestCells[i]=MPI_REQUEST_NULL;
    mpi_rcvRequestCells[i]=MPI_REQUEST_NULL;
  }
}


ZFSStrctrdCommunicationHandle::~ZFSStrctrdCommunicationHandle()
{
  //release all the objects
  for(ZFSInt i=0; i<noNghbrDomains; i++) {
    delete[]  m_bufferCellsRcv[i];
    delete[] m_bufferCellsSnd[i];
    delete[] m_bufferPointsSnd[i];
    delete[]  m_bufferPointsRcv[i];

    delete[] startInfoSNDcells[i];
    delete[] endInfoSNDcells[i];
    delete[] startInfoRCVcells[i];
    delete[] endInfoRCVcells[i];
    delete[] startInfoSNDpoints[i];
    delete[] endInfoSNDpoints[i];
    delete[] startInfoRCVpoints[i];
    delete[] endInfoRCVpoints[i];
    delete[] orderInfo[i];
    delete[] stepInfoRCV[i];
  }

  delete[] m_sndNghbrId;
  delete[] m_rcvNghbrId;
  delete[] m_noNghbrDomainCellBufferSizeSnd;
  delete[] m_noNghbrDomainPointBufferSizeSnd;
  delete[] m_noNghbrDomainCellBufferSizeRcv;
  delete[] m_noNghbrDomainPointBufferSizeRcv;
  delete[] m_bufferCellsRcv;
  delete[] m_bufferCellsSnd;
  delete[] m_bufferPointsSnd;
  delete[] m_bufferPointsRcv;
  delete[] mpi_sndRequest;
  delete[] mpi_rcvRequest;
  delete[] mpi_sndStatus;
  delete[] mpi_rcvStatus;
  delete[] startInfoSNDcells;
  delete[] endInfoSNDcells;
  delete[] startInfoRCVcells;
  delete[] endInfoRCVcells;
  delete[] startInfoSNDpoints;
  delete[] endInfoSNDpoints;
  delete[] startInfoRCVpoints;
  delete[] endInfoRCVpoints;
  delete[] orderInfo;
  delete[] stepInfoRCV;
  delete[] bcId;
  delete[] singularID;
  delete[] m_tagHelperSND;
  delete[] m_tagHelperRCV;
  delete[] mpi_sndRequestCells;
  delete[] mpi_rcvRequestCells;
}

void ZFSStrctrdWaveCommunicationHandle::setBufferSizes()
{
  mpi_sndRequest= new  MPI_Request[noNghbrDomainsSnd];
  mpi_rcvRequest= new MPI_Request[noNghbrDomainsRcv];
  mpi_sndRequestCells= new  MPI_Request[noNghbrDomainsSnd];
  mpi_rcvRequestCells= new MPI_Request[noNghbrDomainsRcv];
  mpi_sndStatus= new MPI_Status[noNghbrDomainsSnd];
  mpi_rcvStatus= new MPI_Status[noNghbrDomainsRcv];
  m_bufferCellsSnd= new ZFSFloat *[noNghbrDomainsSnd];
  m_bufferCellsRcv= new ZFSFloat *[noNghbrDomainsRcv];

  for(ZFSInt i=0; i<noNghbrDomainsSnd; i++) {
    mpi_sndRequest[i]=MPI_REQUEST_NULL;
    mpi_sndRequestCells[i]=MPI_REQUEST_NULL;
    m_bufferCellsSnd[i]  =  new ZFSFloat[m_noNghbrDomainCellBufferSizeSnd[i]];
  }

  for(ZFSInt i=0; i<noNghbrDomainsRcv; i++) {
    mpi_rcvRequest[i]=MPI_REQUEST_NULL;
    mpi_rcvRequestCells[i]=MPI_REQUEST_NULL;
    m_bufferCellsRcv[i]  =  new ZFSFloat[m_noNghbrDomainCellBufferSizeRcv[i]];
  }
}

ZFSStrctrdWaveCommunicationHandle::~ZFSStrctrdWaveCommunicationHandle()
{
  //release all the objects
  for(ZFSInt i=0; i<noNghbrDomainsSnd; i++) {
    delete [] m_bufferCellsSnd[i];
    delete [] startInfoSNDcells[i];
    delete [] endInfoSNDcells[i];
  }

  for(ZFSInt i=0; i<noNghbrDomainsRcv; i++) {
    delete [] m_bufferCellsRcv[i];
    delete [] startInfoRCVcells[i];
    delete [] endInfoRCVcells[i];
  }

  delete[] m_tagHelperRCV;
  delete[] m_tagHelperSND;
  delete[] m_sndNghbrId;
  delete[] m_rcvNghbrId;
  delete[] m_noNghbrDomainCellBufferSizeSnd;
  delete[] m_noNghbrDomainCellBufferSizeRcv;
  delete[] m_bufferCellsRcv;
  delete[] m_bufferCellsSnd;
  delete[] mpi_sndRequest;
  delete[] mpi_rcvRequest;
  delete[] mpi_sndStatus;
  delete[] mpi_rcvStatus;
  delete[] startInfoSNDcells;
  delete[] endInfoSNDcells;
  delete[] startInfoRCVcells;
  delete[] endInfoRCVcells;
  delete[] mpi_sndRequestCells;
  delete[] mpi_rcvRequestCells;
}

