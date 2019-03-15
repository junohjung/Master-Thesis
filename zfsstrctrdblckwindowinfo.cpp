#include "zfsmpi.h"
#include "zfsstrctrdblckwindowinfo.h"
#include "zfsglobals.h"
#include "zfsiolib.h"
#include "zfsstrctrdwindowmapping.h"
#include "vector"
#include "zfspointbox.h"
#include "zfskdtree.h"

template <ZFSInt nDim>
ZFSStrctrdBlckWindowInfo<nDim>::ZFSStrctrdBlckWindowInfo(ZFSStrctrdDecomposition<nDim>* part,
                                                         ZFSId file_ident,
                                                         ZFSInt noInputBlocks_,
                                                         ZFSInt inputBlockId_,
                                                         ZFSInt GhostLayers,
                                                         MPI_Comm strctrdCommunicator,
                                                         const ZFSId noDomains_,
                                                         const ZFSId domainId_)
: m_noInputBlocks(noInputBlocks_),
  m_inputBlockId(inputBlockId_),
  m_noDomains(noDomains_),
  m_domainId(domainId_)
{
  partition=part;
  file_id=file_ident;
  m_noPartitions=m_noInputBlocks;
  m_noGhostLayers=GhostLayers;

  if(noDomains()>m_noPartitions) {
    //this means we have more cpus than inputBlocks
    m_noPartitions=noDomains();
  }

  m_zfsStrctrdComm = strctrdCommunicator;
}

template <ZFSInt nDim>
ZFSStrctrdBlckWindowInfo<nDim>::~ZFSStrctrdBlckWindowInfo()
{
  delete m_myMapWithGC;
  delete m_myMapWithoutGC;

  for(ZFSUint i=0; i<inputWindows.size(); i++) {
    delete inputWindows[i];
  }

  for(ZFSUint i=0; i<m_partitionMapsWithoutGC.size(); i++) {
    delete m_partitionMapsWithoutGC[i];
  }

  for(ZFSUint i=0; i<m_partitionMapsWithGC.size(); i++) {
    delete m_partitionMapsWithGC[i];
  }

  for(ZFSUint i=0; i<sndMap.size(); i++) {
    delete sndMap[i];
  }

  for(ZFSUint i=0; i<rcvMap.size(); i++) {
    delete rcvMap[i];
  }

  for(ZFSUint i=0; i<sndMapPeriodic.size(); i++) {
    delete sndMapPeriodic[i];
  }

  for(ZFSUint i=0; i<rcvMapPeriodic.size(); i++) {
    delete rcvMapPeriodic[i];
  }

  for(ZFSUint i=0; i<physicalAuxDataMap.size(); i++) {
    delete physicalAuxDataMap[i];
  }

  for(ZFSUint i=0; i<waveRcvMap.size(); i++) {
    delete waveRcvMap[i];
  }

  for(ZFSUint i=0; i<waveSndMap.size(); i++) {
    delete waveSndMap[i];
  }

  for(ZFSUint i=0; i<m_wallDistInfoMap.size(); i++) {
    delete m_wallDistInfoMap[i];
  }

  for(ZFSUint i=0; i<m_spongeInfoMap.size(); i++) {
    delete m_spongeInfoMap[i];
  }

  for(ZFSUint i=0; i<globalStrctrdBndryCndMaps.size(); i++) {
    delete globalStrctrdBndryCndMaps[i];
  }

  for(ZFSUint i=0; i<localStrctrdDomainMaps.size(); i++) {
    delete localStrctrdDomainMaps[i];
  }

  for(ZFSUint i=0; i<localStrctrdDomainBndryMaps.size(); i++) {
    delete localStrctrdDomainBndryMaps[i];
  }

  for(ZFSUint i=0; i<window2d.size(); i++) {
    delete window2d[i];
  }

  for(ZFSUint i=0; i<window1d.size(); i++) {
    delete window1d[i];
  }

  for(ZFSUint i=0; i<window0d.size(); i++) {
    delete window0d[i];
  }

  for(ZFSUint i=0; i<singularwindow.size(); i++) {
    delete singularwindow[i];
  }

  for(ZFSUint i=0; i<localSingularMap.size(); i++) {
    delete localSingularMap[i];
  }

  for(ZFSUint i=0; i<localStrctrdBndryCndMaps.size(); i++) {
    delete localStrctrdBndryCndMaps[i];
  }

  for(ZFSUint i=0; i<channelSurfaceIndices.size(); i++) {
    delete channelSurfaceIndices[i];
  }


}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::readWindowCoordinates(ZFSInt cpu,
                                                           ZFSFloat* periodicDisplacements,
                                                           ZFSString fileName)
{
  //currently only works for 3D
  if(nDim==2) {
    return;
  }

  ZFSInt hasConnectionInfo = 0;
  if(globalDomainId()==0){
    ZFSId fid = io_openfile("hdf5", fileName.c_str(), "collective", MPI_COMM_SELF);
    ZFSBool attributeExists = io_checkObj(fid, "/Connectivity", "hasConnectionInfo");
    if(attributeExists) {
      zfs_log << "Grid file has connection info!" << endl;
      io_read_iattribute1(fid, "/Connectivity", "hasConnectionInfo", &hasConnectionInfo);
    }
    MPI_Bcast(&hasConnectionInfo, 1, MPI_INT, 0, m_zfsStrctrdComm);
    io_closefile(fid);
  } else {
    MPI_Bcast(&hasConnectionInfo, 1, MPI_INT, 0, m_zfsStrctrdComm);
  }

  if(hasConnectionInfo) {
    (void) fileName;
    zfs_log << "Connection info available, reading connection information from grid file!" << endl;
    readConnectionWindowInformation(periodicDisplacements, fileName);
  } else {
    if(domainId()==0) {
      cout<<"Starting grid connection info search..."<<endl;
    }
    ZFSInt offset[3],size[3],count=0;
    //3) read in the coordinates of the grid points
    //open file for reading the grid data
    pointType* nodeMap;
    ZFSInt totalGridCells=0;

    ZFSInt**  totalGridBlockDim;
    ZFSInt* totalGridBlockCells;

    zfsAlloc(totalGridBlockDim,  m_noInputBlocks, 3, "totalGridBlockDim", -1 , __CALLING_FUNCTION__ );
    zfsAlloc(totalGridBlockCells, m_noInputBlocks  , "totalGridBlockCells", -1 , __CALLING_FUNCTION__ );

    for(ZFSId i=0; i< m_noInputBlocks; ++i)
    {
      ZFSInt temp =1;
      for(ZFSId dim=0; dim<3; ++dim)
      {
        temp*=partition->inputBoxInfo[i]->DirLast[dim]+1;
        totalGridBlockDim[i][dim]=partition->inputBoxInfo[i]->DirLast[dim]+1;//number of points in the grid File
      }
      totalGridBlockCells[i]=temp;
      if(temp!=1) totalGridCells+=temp;
    }

    ZFSInt countNode=0;
    for (ZFSInt i=0;i< m_noInputBlocks; i++)
    {
      countNode = countNode + 2* totalGridBlockDim[i][1]*totalGridBlockDim[i][2]
        + 2*totalGridBlockDim[i][2]*(totalGridBlockDim[i][0]-2)
        + 2*(totalGridBlockDim[i][0]-2)*(totalGridBlockDim[i][1]-2);
    }

    ZFSFloatScratchSpace coordinates(3, countNode, __CALLING_FUNCTION__, "coordinates");
    coordinates.fill(-1.01010101);
    nodeMap=new pointType [countNode];

    ZFSInt memsize=0;
    //CPU 0 reads the coordinates for all 6 sides of all blocks
    if(domainId()==0) {
      cout<<"Reading in all coordinates of all block faces"<<endl;
    }
    for(ZFSId i=0; i< m_noInputBlocks; ++i) {
      ZFSString bName="/block";
      stringstream number;
      number <<partition->inputBoxInfo[i]->inputBoxID<<"/";
      bName += number.str();

      //in grid file 2:i  1:j   0:k
      //////////////////////////////
      ///////// FACE 1 /////////////
      //////////////////////////////
      //face +x  face 1
      if(cpu==0) {
        offset[2]=0; size[2]=1; //i
        offset[1]=0; size[1]=partition->inputBoxInfo[i]->DirLast[1]+1; //j
        offset[0]=0; size[0]=partition->inputBoxInfo[i]->DirLast[0]+1; //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, &coordinates(0,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, &coordinates(1,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, &coordinates(2,memsize));
      } else {
        offset[2]=0; size[2]=0; //i
        offset[1]=0; size[1]=0; //j
        offset[0]=0; size[0]=0; //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, NULL);
      }

      //write all computational coordinates
      //of the face to nodeMap object
      for(ZFSInt a=offset[2];a<offset[2]+size[2];++a) {
        for(ZFSInt c=offset[0];c<offset[0]+size[0];++c) {
          for(ZFSInt b=offset[1];b<offset[1]+size[1];++b) {
            nodeMap[count].blockId=partition->inputBoxInfo[i]->inputBoxID;
            nodeMap[count].pos[0]=a;
            nodeMap[count].pos[1]=b;
            nodeMap[count].pos[2]=c;
            count++;
          }
        }
      }
      memsize+=size[0]*size[1]*size[2];

      //////////////////////////////
      ///////// FACE 2 /////////////
      //////////////////////////////
      //face -x  face 2
      if(cpu==0){
        offset[2]=partition->inputBoxInfo[i]->DirLast[2]; size[2]=1; //i
        offset[1]=0; size[1]=partition->inputBoxInfo[i]->DirLast[1]+1; //j
        offset[0]=0; size[0]=partition->inputBoxInfo[i]->DirLast[0]+1; //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, &coordinates(0,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, &coordinates(1,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, &coordinates(2,memsize));
      } else {
        offset[2]=0; size[2]=0; //i
        offset[1]=0; size[1]=0; //j
        offset[0]=0; size[0]=0; //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, NULL);
      }

      for(ZFSInt a=offset[2];a<offset[2]+size[2];++a) {
        for(ZFSInt c=offset[0];c<offset[0]+size[0];++c) {
          for(ZFSInt b=offset[1];b<offset[1]+size[1];++b) {
            nodeMap[count].blockId=partition->inputBoxInfo[i]->inputBoxID;
            nodeMap[count].pos[0]=a;
            nodeMap[count].pos[1]=b;
            nodeMap[count].pos[2]=c;
            count++;
          }
        }
      }
      memsize+=size[0]*size[1]*size[2];

      //////////////////////////////
      ///////// FACE 3 /////////////
      //////////////////////////////
      //face +y  face 3
      if(cpu==0){
        offset[2]=1; size[2]=partition->inputBoxInfo[i]->DirLast[2]-1;  //i
        offset[1]=0; size[1]=1;    //j
        offset[0]=0; size[0]=partition->inputBoxInfo[i]->DirLast[0]+1;    //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, &coordinates(0,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, &coordinates(1,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, &coordinates(2,memsize));
      } else {
        offset[2]=0; size[2]=0; //i
        offset[1]=0; size[1]=0; //j
        offset[0]=0; size[0]=0; //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, NULL);
      }

      for(ZFSInt b=offset[1];b<offset[1]+size[1];++b) {
        for(ZFSInt c=offset[0];c<offset[0]+size[0];++c) {
          for(ZFSInt a=offset[2];a<offset[2]+size[2];++a) {
            nodeMap[count].blockId=partition->inputBoxInfo[i]->inputBoxID;
            nodeMap[count].pos[0]=a;
            nodeMap[count].pos[1]=b;
            nodeMap[count].pos[2]=c;
            count++;
          }
        }
      }
      memsize+=size[0]*size[1]*size[2];

      //////////////////////////////
      ///////// FACE 4 /////////////
      //////////////////////////////
      //face -y  face 4
      if(cpu==0){
        offset[2]=1; size[2]=partition->inputBoxInfo[i]->DirLast[2]-1; //i
        offset[1]=partition->inputBoxInfo[i]->DirLast[1]; size[1]=1; //j
        offset[0]=0; size[0]=partition->inputBoxInfo[i]->DirLast[0]+1; //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, &coordinates(0,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, &coordinates(1,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, &coordinates(2,memsize));
      } else {
        offset[2]=0; size[2]=0; //i
        offset[1]=0; size[1]=0; //j
        offset[0]=0; size[0]=0; //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, NULL);
      }

      for(ZFSInt b=offset[1];b<offset[1]+size[1];++b) {
        for(ZFSInt c=offset[0];c<offset[0]+size[0];++c) {
          for(ZFSInt a=offset[2];a<offset[2]+size[2];++a) {
            nodeMap[count].blockId=partition->inputBoxInfo[i]->inputBoxID;
            nodeMap[count].pos[0]=a;
            nodeMap[count].pos[1]=b;
            nodeMap[count].pos[2]=c;
            count++;
          }
        }
      }
      memsize+=size[0]*size[1]*size[2];

      //////////////////////////////
      ///////// FACE 5 /////////////
      //////////////////////////////
      //face +z  face 5
      if(cpu==0){
        offset[2]=1; size[2]=partition->inputBoxInfo[i]->DirLast[2]-1;  //i
        offset[1]=1; size[1]=partition->inputBoxInfo[i]->DirLast[1]-1;    //j
        offset[0]=0; size[0]=1;    //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, &coordinates(0,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, &coordinates(1,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, &coordinates(2,memsize));
      } else {
        offset[2]=0; size[2]=0;           //i
        offset[1]=0; size[1]=0; //j
        offset[0]=0; size[0]=0; //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, NULL);
      }

      for(ZFSInt c=offset[0];c<offset[0]+size[0];++c) {
        for(ZFSInt b=offset[1];b<offset[1]+size[1];++b) {
          for(ZFSInt a=offset[2];a<offset[2]+size[2];++a) {
            nodeMap[count].blockId=partition->inputBoxInfo[i]->inputBoxID;
            nodeMap[count].pos[0]=a;
            nodeMap[count].pos[1]=b;
            nodeMap[count].pos[2]=c;
            count++;
          }
        }
      }
      memsize+=size[0]*size[1]*size[2];

      //////////////////////////////
      ///////// FACE 6 /////////////
      //////////////////////////////
      //face -z  face 6
      if(cpu==0){
        offset[2]=1; size[2]=partition->inputBoxInfo[i]->DirLast[2]-1; //i
        offset[1]=1; size[1]=partition->inputBoxInfo[i]->DirLast[1]-1; //j
        offset[0]=partition->inputBoxInfo[i]->DirLast[0]; size[0]=1; //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, &coordinates(0,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, &coordinates(1,memsize));
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, &coordinates(2,memsize));
      } else {
        offset[2]=0; size[2]=0; //i
        offset[1]=0; size[1]=0; //j
        offset[0]=0; size[0]=0; //k
        io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, NULL);
        io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, NULL);
      }

      for(ZFSInt c=offset[0];c<offset[0]+size[0];++c) {
        for(ZFSInt b=offset[1];b<offset[1]+size[1];++b) {
          for(ZFSInt a=offset[2];a<offset[2]+size[2];++a){
            nodeMap[count].blockId=partition->inputBoxInfo[i]->inputBoxID;
            nodeMap[count].pos[0]=a;
            nodeMap[count].pos[1]=b;
            nodeMap[count].pos[2]=c;
            count++;
          }
        }
      }
      memsize+=size[0]*size[1]*size[2];
    }

    //now broadcast the surface coordinates to every processor
    MPI_Bcast(&count,  1, MPI_INT, 0, m_zfsStrctrdComm);
    MPI_Bcast(&memsize,1, MPI_INT, 0, m_zfsStrctrdComm);
    MPI_Bcast(&coordinates(0,0),memsize, MPI_DOUBLE, 0, m_zfsStrctrdComm);
    MPI_Bcast(&coordinates(1,0),memsize, MPI_DOUBLE, 0, m_zfsStrctrdComm);
    MPI_Bcast(&coordinates(2,0),memsize, MPI_DOUBLE, 0, m_zfsStrctrdComm);

    //also broadcast the nodeMap to every processor
    MPI_Datatype matrix;
    MPI_Type_contiguous(6, MPI_INT, &matrix);
    MPI_Type_commit(&matrix);
    MPI_Bcast(nodeMap,memsize, matrix, 0, m_zfsStrctrdComm);

    pointType* nodeMapP;
    ZFSFloat tmppoint[3];
    ZFSInt pcount=0,plocation=0,numofpoints[3];

    for(ZFSInt i=0; i<noInputWindowInformation; i++) {
      if(inputWindows[i]->BC>4000&&inputWindows[i]->BC<5000) {
        for(ZFSInt j=0;j<3;j++) {
          numofpoints[j]=inputWindows[i]->endindex[j]-inputWindows[i]->startindex[j]+1;
        }
        pcount+= numofpoints[0]*numofpoints[1]*numofpoints[2];
      }
    }

    ZFSFloatScratchSpace periodicCoordinates(3, pcount, __CALLING_FUNCTION__, "periodicCoordinates");
    periodicCoordinates.fill(-1.01010101);
    nodeMapP=new pointType [pcount];

    if(domainId()==0) {
      cout<<"Loading periodic face coordinates"<<endl;
    }
    ZFSInt pmemsize=0;
    for(ZFSInt i=0; i<noInputWindowInformation; i++) {
      if(inputWindows[i]->BC>=4000&&inputWindows[i]->BC<5000) {
        ZFSString bName="/block";
        stringstream number;
        number << inputWindows[i]->inputBlockId<<"/";
        bName += number.str();

        if(cpu==0) {
          for(ZFSInt j=0;j<3;j++) {
            offset[j]=inputWindows[i]->startindex[2-j];
            size[j]=inputWindows[i]->endindex[2-j]-inputWindows[i]->startindex[2-j]+1;
          }

          io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, &periodicCoordinates(0,plocation));
          io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, &periodicCoordinates(1,plocation));
          io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, &periodicCoordinates(2,plocation));
        } else {
          offset[2]=0; size[2]=0; //i
          offset[1]=0; size[1]=0; //j
          offset[0]=0; size[0]=0; //k
          io_read_ddataset_part1d1(file_id, bName.c_str(), "x", 3, offset, size, NULL);
          io_read_ddataset_part1d1(file_id, bName.c_str(), "y", 3, offset, size, NULL);
          io_read_ddataset_part1d1(file_id, bName.c_str(), "z", 3, offset, size, NULL);
        }

        for(ZFSInt c=offset[0];c<offset[0]+size[0];++c) {
          for(ZFSInt b=offset[1];b<offset[1]+size[1];++b) {
            for(ZFSInt a=offset[2];a<offset[2]+size[2];++a) {

              nodeMapP[plocation].BC=inputWindows[i]->BC;
              nodeMapP[plocation].blockId=inputWindows[i]->inputBlockId;
              nodeMapP[plocation].pos[0]=a;
              nodeMapP[plocation].pos[1]=b;
              nodeMapP[plocation].pos[2]=c;
              plocation++;
            }
          }
        }
        pmemsize+=size[0]*size[1]*size[2];
      }
    }

    MPI_Bcast(&plocation,  1, MPI_INT, 0, m_zfsStrctrdComm);
    MPI_Bcast(&pmemsize,   1, MPI_INT, 0, m_zfsStrctrdComm);

    MPI_Bcast(&periodicCoordinates(0,0),pmemsize, MPI_DOUBLE, 0, m_zfsStrctrdComm);
    MPI_Bcast(&periodicCoordinates(1,0),pmemsize, MPI_DOUBLE, 0, m_zfsStrctrdComm);
    MPI_Bcast(&periodicCoordinates(2,0),pmemsize, MPI_DOUBLE, 0, m_zfsStrctrdComm);

    MPI_Bcast(nodeMapP, pmemsize, matrix, 0, m_zfsStrctrdComm);


    /////////////////////////////////
    /////// CONNECTIONS /////////////
    /////////////////////////////////
    if(domainId()==0) {
      cout<<"Building up connections for multiblock connection search"<<endl;
    }
    vector < Point<3> > pts;
    for(ZFSId j=0; j<count; ++j) {
      Point<3> a(coordinates(0,j),coordinates(1,j),coordinates(2,j));
      pts.push_back(a);
      nodeMap[j].found=false;
    }

    ZFSFloat m_gridEps = 0.0000001;

    KDtree<3> tree(pts);
    ZFSBool add;
    ZFSInt numConnections = 0,numSingularConnections = 0,nfound,tempnum,tempcount;
    ZFSInt results[10];
    // results=new ZFSInt [10];
    for(ZFSInt i = 0;i<count;++i) {
      if (!nodeMap[i].found) {
        Point<3> a(coordinates(0,i),coordinates(1,i),coordinates(2,i));
        nfound=tree.locatenear(a,m_gridEps,results,10,false); //0.00001
        for(ZFSInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
          for(ZFSInt countNode3 = countNode2+1;countNode3 < nfound;++countNode3) {
            if (addConnection( 6000,
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               nodeMap[results[countNode3]].blockId,
                               nodeMap[results[countNode3]].pos)) {
              numConnections = numConnections + 1;
            }
          }
          nodeMap[results[countNode2]].found = true;
        }

        //if three points share the same
        //coordinate it must be a 3-star
        if(nfound==3) {
          add=false;
          tempcount=0;

          //check if it is a singularity 3-star or normal 3-point connection
          for(ZFSInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
            for(ZFSInt j=0;j< m_noInputBlocks; ++j) {
              if(nodeMap[results[countNode2]].blockId==partition->inputBoxInfo[j]->inputBoxID) {
                tempnum = j;
                break;
              }
            }
            for(ZFSInt j=0;j<3;++j) {
              //pos[]: ijk  DirLast[]: kji
              if(nodeMap[results[countNode2]].pos[j]==0||
                 nodeMap[results[countNode2]].pos[j]==partition->inputBoxInfo[tempnum]->DirLast[2-j]) {
                ++tempcount;
              }
            }
          }

          if(tempcount==9||tempcount==6) {
            add=true;
          }

          if(add) {
            for(ZFSInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
              if (addConnection( 6000,
                                 nodeMap[results[countNode2]].blockId,
                                 nodeMap[results[countNode2]].pos,
                                 nodeMap[results[countNode2]].blockId,
                                 nodeMap[results[countNode2]].pos,
                                 3)) {
                numSingularConnections = numSingularConnections + 1;
              }
            }
          }
        }

        //if three points share the same
        //coordinate it must be a 5-star
        if(nfound==5) {
          for(ZFSInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
            if (addConnection( 6000,
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               5)) {
              numSingularConnections = numSingularConnections + 1;
            }
          }
        }
        //coordinate it must be a 6-star
        if(nfound==6) {
          for(ZFSInt countNode2 = 0; countNode2 < nfound; ++countNode2) {
            if (addConnection( 6000,
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               nodeMap[results[countNode2]].blockId,
                               nodeMap[results[countNode2]].pos,
                               6)) {
              numSingularConnections = numSingularConnections + 1;
            }
          }
        }
      }
    }

    if(domainId()==0) {
      cout<<"Computing periodic window displacements"<<endl;
    }
    /////////////////////////////////////////////////
    //////// PERIODIC WINDOW DISPLACEMENTS //////////
    /////////////////////////////////////////////////
    ZFSFloatScratchSpace periodicWindowCenter(m_noInputBlocks, nDim, 2*nDim, __CALLING_FUNCTION__, "periodicWindowCOGS");
    ZFSIntScratchSpace periodicWindowNoNodes(m_noInputBlocks, 2*nDim, __CALLING_FUNCTION__, "periodicWindowNoNodes");
    cout.precision(10);
    periodicWindowCenter.fill(F0);
    periodicWindowNoNodes.fill(0);

    //compute displacement
    const ZFSInt periodicOffset = 4401;
    for(ZFSInt i = 0;i<pcount;++i) {
      if(nodeMapP[i].BC < 4401 || nodeMapP[i].BC > 4406) {
        continue;
      }

      //add up all coordinates
      for(ZFSId dim=0; dim<nDim; dim++) {
        periodicWindowCenter(nodeMapP[i].blockId, dim, nodeMapP[i].BC-periodicOffset) += periodicCoordinates(dim,i);
      }

      periodicWindowNoNodes(nodeMapP[i].blockId, nodeMapP[i].BC-periodicOffset) += 1;
    }

    //now compute the center
    //(not weighted, should be basically the same surface in another location)
    //and compute distance between the two centers
    for(ZFSId periodicWindowId = 0; periodicWindowId < 2*nDim; periodicWindowId++) {
      if(periodicWindowNoNodes(periodicWindowId) <= 0) {
        continue;
      }

      for(ZFSId block=0; block<m_noInputBlocks; block++) {
        for(ZFSId dim=0; dim<nDim; dim++) {
          periodicWindowCenter(block, dim, periodicWindowId) /= (ZFSFloat)periodicWindowNoNodes(block,periodicWindowId);
        }
      }

      if(periodicWindowId%2==1) {
        for(ZFSId dim=0; dim<nDim; dim++) {
          const ZFSId displacementId = (ZFSFloat)(periodicWindowId+1)/2.0 - 1;
          periodicDisplacements[dim*nDim+displacementId] = periodicWindowCenter(m_inputBlockId, dim, periodicWindowId) - periodicWindowCenter(m_inputBlockId, dim,periodicWindowId-1);
          zfs_log << "periodicWindowCenter for window: "
                  << periodicWindowId+periodicOffset << "  dim: " << dim
                  << " displacementId: " << displacementId << " displacement: "
                  << periodicDisplacements[dim*nDim+displacementId] << endl;
        }
      }
    }

    /////////////////////////////////////////////////
    //////// PERIODIC CONNECTIONS ///////////////////
    /////////////////////////////////////////////////

    for(ZFSInt i = 0;i<pcount;++i) {
      tmppoint[0]=periodicCoordinates(0,i);
      tmppoint[1]=periodicCoordinates(1,i);
      tmppoint[2]=periodicCoordinates(2,i);
      periodicPointsChange(tmppoint,nodeMapP[i].BC,periodicDisplacements);
      Point<3> aaa(tmppoint[0],tmppoint[1],tmppoint[2]);
      nfound=tree.locatenear(aaa,m_gridEps,results,10,false); //0.0001
      for(ZFSInt countNode2 = 0;countNode2 < nfound;++countNode2) {
        if (addConnection( nodeMapP[i].BC,
                           nodeMapP[i].blockId,
                           nodeMapP[i].pos,
                           nodeMap[results[countNode2]].blockId,
                           nodeMap[results[countNode2]].pos)) {
          numConnections = numConnections + 1;
        }
      }

      if(nfound==5) {
        if (addConnection( nodeMapP[i].BC,
                           nodeMapP[i].blockId,
                           nodeMapP[i].pos,
                           nodeMapP[i].blockId,
                           nodeMapP[i].pos,
                           5)) {
          numSingularConnections = numSingularConnections + 1;
        }
      }
      //6 star for periodic bc
      if(nfound==6) {
        if (addConnection( nodeMapP[i].BC,
                           nodeMapP[i].blockId,
                           nodeMapP[i].pos,
                           nodeMapP[i].blockId,
                           nodeMapP[i].pos,
                           6)) {
          numSingularConnections = numSingularConnections + 1;
        }
      }
    }

    if(domainId()==0) {
      cout<<"Assemble multiblock connections"<<endl;
    }
    ///////////////////////////////////////////////////
    /////// CREATE THE WINDOWS FROM CONNECTIONS ///////
    ///////////////////////////////////////////////////
    multiBlockAssembling();
    singularityAssembling();


    ZFSStrctrdWindowMap* windowMap;
    for(ZFSInt i=0;i<(ZFSInt)window2d.size();++i) {
      windowMap=new ZFSStrctrdWindowMap(nDim);
      mapCreate(window2d[i]->Id1, window2d[i]->start1, window2d[i]->end1, window2d[i]->step1,
                window2d[i]->Id2, window2d[i]->start2, window2d[i]->end2, window2d[i]->step2,
                window2d[i]->order, window2d[i]->BC, windowMap);
      windowMap->Nstar=-1;windowMap->SingularId=-1;
      windowMap->dc1=window2d[i]->dc1;
      windowMap->dc2=window2d[i]->dc2;

      //i-surface: change sign of I; same for J and K
      if(windowMap->dc1*windowMap->dc2>0) {
        ZFSInt tempdim=abs(windowMap->dc2)-1;
        windowMap->step2[tempdim]=-1;
      }

      //add the map to the list!
      globalStrctrdBndryCndMaps.push_back(windowMap);


      if( window2d[i]->BC==6000) {
        windowMap=new ZFSStrctrdWindowMap(nDim);
        mapCreate(window2d[i]->Id1, window2d[i]->start1, window2d[i]->end1, window2d[i]->step1,
                  window2d[i]->Id2, window2d[i]->start2, window2d[i]->end2 ,window2d[i]->step2,
                  window2d[i]->order, window2d[i]->BC, windowMap);
        windowMap->Nstar=-1;windowMap->SingularId=-1;
        windowMap->dc1=window2d[i]->dc1;
        windowMap->dc2=window2d[i]->dc2;
        //i-surface: change sign of I; same for J and K
        if(windowMap->dc1*windowMap->dc2>0) {
          ZFSInt tempdim=abs(windowMap->dc2)-1;
          windowMap->step2[tempdim]=-1;
        }

        mapInvert1(windowMap);
        mapNormalize3(windowMap);
        globalStrctrdBndryCndMaps.push_back(windowMap);
      }
    }


    for(ZFSInt i=0;i<(ZFSInt)window1d.size();++i) {
      ZFSInt test;

      windowMap=new ZFSStrctrdWindowMap(nDim);
      mapCreate(window1d[i]->Id1, window1d[i]->start1, window1d[i]->end1, window1d[i]->step1,
                window1d[i]->Id2,  window1d[i]->start2, window1d[i]->end2 ,window1d[i]->step2,
                window1d[i]->order, window1d[i]->BC, windowMap);
      windowMap->Nstar=-1;windowMap->SingularId=-1;
      windowMap->dc1=window1d[i]->dc1;
      windowMap->dc2=window1d[i]->dc2;

      //i-surface: change sign of I; same for J and K
      for(ZFSInt j=0;j<3;++j) {
        if(windowMap->start1[j]==windowMap->end1[j]) {
          if(windowMap->start1[j]==0&&windowMap->start2[windowMap->order[j]]==0) {
            windowMap->step2[windowMap->order[j]]=-1;
          }
          if(windowMap->start1[j]>0&&windowMap->start2[windowMap->order[j]]>0) {
            windowMap->step2[windowMap->order[j]]=-1;
          }
        }
      }

      //now for 5 star communicaitons (or ...... if needed)
      //3 star does not need additional information
      for (ZFSInt j=0; j<(ZFSInt)singularwindow.size();++j)
      {
        test=mapCompare11(windowMap,singularwindow[j]);
        if(test)
        {
          windowMap->Nstar=singularwindow[j]->Nstar;
          windowMap->SingularId=singularwindow[j]->SingularId;
        }
      }

      //add the map to the list!
      globalStrctrdBndryCndMaps.push_back(windowMap);

      if( window1d[i]->BC==6000) {
        windowMap=new ZFSStrctrdWindowMap(nDim);
        mapCreate(window1d[i]->Id1, window1d[i]->start1, window1d[i]->end1, window1d[i]->step1,
                  window1d[i]->Id2, window1d[i]->start2, window1d[i]->end2 ,window1d[i]->step2,
                  window1d[i]->order, window1d[i]->BC, windowMap);
        windowMap->Nstar=-1;windowMap->SingularId=-1;
        windowMap->dc1=window1d[i]->dc1;
        windowMap->dc2=window1d[i]->dc2;

        //i-surface: change sign of I; same for J and K
        for(ZFSInt j=0;j<3;++j) {
          if(windowMap->start1[j]==windowMap->end1[j]) {
            if(windowMap->start1[j]==0&&windowMap->start2[windowMap->order[j]]==0) {
              windowMap->step2[windowMap->order[j]]=-1;
            }
            if(windowMap->start1[j]>0&&windowMap->start2[windowMap->order[j]]>0) {
              windowMap->step2[windowMap->order[j]]=-1;
            }
          }
        }

        mapInvert1(windowMap);
        mapNormalize3(windowMap);

        //now for 5 star communicaitons (or ...... if needed)
        //3 star does not need additional information
        for (ZFSInt j=0; j<(ZFSInt)singularwindow.size();++j) {
          test=mapCompare11(windowMap,singularwindow[j]);
          if(test) {
            windowMap->Nstar=singularwindow[j]->Nstar;
            windowMap->SingularId=singularwindow[j]->SingularId;
          }
        }

        globalStrctrdBndryCndMaps.push_back(windowMap);
      }
    }


    for(ZFSInt i=0;i<(ZFSInt)window0d.size();++i) {
      ZFSInt test;

      //point communicaiton (rare)
      windowMap=new ZFSStrctrdWindowMap(nDim);
      mapCreate(window0d[i]->Id1, window0d[i]->start1, window0d[i]->end1, window0d[i]->step1,
                window0d[i]->Id2 ,  window0d[i]->start2, window0d[i]->end2 ,window0d[i]->step2,
                window0d[i]->order, window0d[i]->BC, windowMap);
      windowMap->Nstar=-1;windowMap->SingularId=-1;
      windowMap->dc1=window0d[i]->dc1;
      windowMap->dc2=window0d[i]->dc2;
      //i-surface: change sign of I; same for J and K
      for(ZFSInt j=0;j<3;++j) {
        if(windowMap->start1[j]==windowMap->end1[j]) {
          if(windowMap->start1[j]==0&&windowMap->start2[windowMap->order[j]]==0) {
            windowMap->step2[windowMap->order[j]]=-1;
          }
          if(windowMap->start1[j]>0&&windowMap->start2[windowMap->order[j]]>0) {
            windowMap->step2[windowMap->order[j]]=-1;
          }
        }
      }

      //now for 5 star communicaitons (or ...... if needed)
      //3 star does not need additional information
      for (ZFSInt j=0; j<(ZFSInt)singularwindow.size();++j) {
        test=mapCompare11(windowMap,singularwindow[j]);
        if(test) {
          windowMap->Nstar=singularwindow[j]->Nstar;
          windowMap->SingularId=singularwindow[j]->SingularId;
        }
      }

      //add the map to the list!
      globalStrctrdBndryCndMaps.push_back(windowMap);

      if( window0d[i]->BC==6000) {
        windowMap=new ZFSStrctrdWindowMap(nDim);
        mapCreate(window0d[i]->Id1, window0d[i]->start1, window0d[i]->end1, window0d[i]->step1,
                  window0d[i]->Id2, window0d[i]->start2, window0d[i]->end2 ,window0d[i]->step2,
                  window0d[i]->order, window0d[i]->BC, windowMap);
        windowMap->Nstar=-1;windowMap->SingularId=-1;
        windowMap->dc1=window0d[i]->dc1;
        windowMap->dc2=window0d[i]->dc2;

        //i-surface: change sign of I; same for J and K
        for(ZFSInt j=0;j<3;++j) {
          if(windowMap->start1[j]==windowMap->end1[j]) {
            if(windowMap->start1[j]==0&&windowMap->start2[windowMap->order[j]]==0) {
              windowMap->step2[windowMap->order[j]]=-1;
            }
            if(windowMap->start1[j]>0&&windowMap->start2[windowMap->order[j]]>0) {
              windowMap->step2[windowMap->order[j]]=-1;
            }
          }
        }

        mapInvert1(windowMap);
        mapNormalize3(windowMap);

        //now for 5 star communicaitons (or ...... if needed)
        //3 star does not need additional information
        for (ZFSInt j=0; j<(ZFSInt)singularwindow.size();++j)
        {
          test=mapCompare11(windowMap,singularwindow[j]);
          if(test)
          {
            windowMap->Nstar=singularwindow[j]->Nstar;
            windowMap->SingularId=singularwindow[j]->SingularId;
          }
        }
        globalStrctrdBndryCndMaps.push_back(windowMap);
      }
    }

    writeConnectionWindowInformation(periodicDisplacements);

    if(domainId()==0) {
      cout<<"Connection identification and window creation finished!"<<endl;
    }
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::multiBlockAssembling()
{
  ZFSInt countConnection=0;//numConnections;
  ZFSStrctrdWindowMap* newWindow;
  ZFSBool found,notexisted,labell;
  connectionNode temp1;
  ZFSInt numWindows[3]={0,0,0};

  set<connectionNode>::iterator it;
  //numConnections=connectionset.size();

  while(connectionset.size()!=0) {
    auto element=connectionset.begin();
    ZFSInt order[3]={-1,-1,-1};
    ZFSInt step1[3]={1,1,1};
    ZFSInt step2[3]={1,1,1};
    ZFSInt pos1[3],pos2[3],b1,b2;

    pos1[0]=element->pos1[0];
    pos1[1]=element->pos1[1];
    pos1[2]=element->pos1[2];
    pos2[0]=element->pos2[0];
    pos2[1]=element->pos2[1];
    pos2[2]=element->pos2[2];
    b1=element->blockId1;
    b2=element->blockId2;

    newWindow=new ZFSStrctrdWindowMap(nDim);
    mapCreate(b1, pos1, pos1, step1,
              b2, pos2, pos2, step2,
              order, element->BC, newWindow);
    newWindow->Nstar=-1;

    for (ZFSInt i=0;i<3;++i) {
      found=false;
      pos1[0]=newWindow->start1[0];
      pos1[1]=newWindow->start1[1];
      pos1[2]=newWindow->start1[2];

      if(newWindow->order[i]==-1) {
        //D1+

        pos1[i]=newWindow->start1[i] + 1;
        for (ZFSInt j=0;j<3;++j) {
          //D2
          pos2[0]=newWindow->start2[0];
          pos2[1]=newWindow->start2[1];
          pos2[2]=newWindow->start2[2];

          notexisted=true;
          for(ZFSInt k=0;k<3;++k) {
            if(newWindow->order[k]==j)notexisted=false;
          }

          if(notexisted) {
            pos2[j]=newWindow->start2[j]+1;

            temp1.BC=newWindow->BC;
            temp1.blockId1=newWindow->Id1; temp1.blockId2=newWindow->Id2;
            temp1.pos1[0]=pos1[0];temp1.pos1[1]=pos1[1];temp1.pos1[2]=pos1[2];
            temp1.pos2[0]=pos2[0];temp1.pos2[1]=pos2[1];temp1.pos2[2]=pos2[2];
            temp1.Nstar=-1;
            found=findConnection(temp1);

            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }

            //D2-
            pos2[j]=newWindow->start2[j]-1;
            temp1.BC=newWindow->BC;
            temp1.blockId1=newWindow->Id1;
            temp1.blockId2=newWindow->Id2;
            temp1.pos1[0]=pos1[0];temp1.pos1[1]=pos1[1];
            temp1.pos1[2]=pos1[2];temp1.pos2[0]=pos2[0];
            temp1.pos2[1]=pos2[1];temp1.pos2[2]=pos2[2];
            temp1.Nstar=-1;
            found=findConnection(temp1);
            if(found)
            {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }

        if(found) {
          continue;
        }

        //D1-
        pos1[i]=newWindow->start1[i] - 1;
        for (ZFSInt j=0;j<3;++j) {
          //D2+
          pos2[0]=newWindow->start2[0];
          pos2[1]=newWindow->start2[1];
          pos2[2]=newWindow->start2[2];
          notexisted=true;
          for(ZFSInt k=0;k<3;++k) {
            if(newWindow->order[k]==j) {
              notexisted=false;
            }
          }

          if(notexisted) {
            pos2[j]=newWindow->start2[j]+1;
            temp1.BC=newWindow->BC;
            temp1.blockId1=newWindow->Id1; temp1.blockId2=newWindow->Id2;
            temp1.pos1[0]=pos1[0];
            temp1.pos1[1]=pos1[1];
            temp1.pos1[2]=pos1[2];
            temp1.pos2[0]=pos2[0];
            temp1.pos2[1]=pos2[1];
            temp1.pos2[2]=pos2[2];
            temp1.Nstar=-1;
            found=findConnection(temp1);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }

            //D2-
            pos2[j]=newWindow->start2[j]-1;
            temp1.BC=newWindow->BC;
            temp1.blockId1=newWindow->Id1; temp1.blockId2=newWindow->Id2;
            temp1.pos1[0]=pos1[0];
            temp1.pos1[1]=pos1[1];
            temp1.pos1[2]=pos1[2];
            temp1.pos2[0]=pos2[0];
            temp1.pos2[1]=pos2[1];
            temp1.pos2[2]=pos2[2];
            temp1.Nstar=-1;
            found=findConnection(temp1);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }
        if(found) {
          continue;
        }
      }
    }

    ZFSInt ordercount=0;
    ZFSBool facewindow,directionInWindow[3];
    for (ZFSInt i=0;i<3;++i) {
      directionInWindow[i]=false;
      if(newWindow->order[i]!=-1) {
        ordercount++;directionInWindow[i]=true;
      }
    }

    if(ordercount>2) {
      cout<<"Invalid volume mapping found! Are your blocks overlapping or is the grid epsilon too large?"<<endl;

    }

    if(ordercount==2) {
      facewindow=true;
    } else {
      facewindow=false;
    }

    for (ZFSInt i=0;i<3;++i) {
      ZFSInt j=0;

      if(newWindow->order[i]==-1) {
        for(j=0;j<3;++j) {
          labell=true;
          for(ZFSInt k=0;k<3;++k) {
            if ( newWindow->order[k]==j ) {
              labell=false;
            }
          }

          if(labell==true) {
            newWindow->order[i]=j;
            break;
          }
        }

        if(facewindow) {
          if(newWindow->start1[i]==0) {
            newWindow->dc1=i+1;
          } else {
            newWindow->dc1=-i-1;
          }

          if(newWindow->start2[j]==0) {
            newWindow->dc2=j+1;
          } else {
            newWindow->dc2=-j-1;
          }
        } else {
          newWindow->dc1= 999;
          newWindow->dc2= 999;
        }
      }
    }

    ZFSInt start1[3],end1[3],start2[3],end2[3];
    ZFSBool goGo=true;
    ZFSInt countDim,countDim2;
    ZFSInt ii,jj,kk;
    while (goGo) {
      goGo=false;
      for ( countDim=0; countDim<3; ++countDim) {
        if(directionInWindow[countDim]) {
          countDim2=newWindow->order[countDim];

          for(ZFSInt i=0; i<3; ++i) {
            start1[i]=newWindow->start1[i];
            end1[i]=newWindow->end1[i];
            start2[i]=newWindow->start2[i];
            end2[i]=newWindow->end2[i];
          }
          end1[countDim]  = end1[countDim]  + newWindow->step1[countDim];
          end2[countDim2] = end2[countDim2] + newWindow->step2[countDim2];
          start1[countDim] = end1[countDim];
          start2[countDim2] = end2[countDim2];

          pos2[newWindow->order[2]] = start2[newWindow->order[2]];
          kk=start1[2];
          do {
            pos2[newWindow->order[1]] = start2[newWindow->order[1]];
            jj=start1[1];
            do {
              pos2[newWindow->order[0]] = start2[newWindow->order[0]];
              ii=start1[0];
              do {

                pos1[0]=ii;pos1[1]=jj;pos1[2]=kk;
                if (newWindow->BC==6000) {
                  labell=true;
                } else {
                  labell=false;
                }

                connectionNode temp2(newWindow->BC,newWindow->Id1,pos1,newWindow->Id2,pos2,labell);
                found=findConnection(temp2);

                if(!found&&domainId()==0&&newWindow->BC==4401&&newWindow->Id1==2&&newWindow->Id2!=2) {
                  // temp2.print();
                  connectionNode temp3(newWindow->BC,newWindow->Id2,pos2,newWindow->Id1,pos1,labell);
                  found=findConnection(temp3);

                  if(found) {
                    cout<<"wooooo!! somthing is  wrong and i do not know where and why!!!!!!"<<endl;
                  }
                }
                if (!found) {
                  break;
                }
                pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
                ii=ii+newWindow->step1[0];

              } while(ii>=start1[0]&&ii<=end1[0]);
              pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
              jj=jj+newWindow->step1[1];

            } while(jj>=start1[1]&&jj<=end1[1]);
            pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
            kk=kk+newWindow->step1[2];

          } while(kk>=start1[2]&&kk<=end1[2]);

          if (found) {
            //all connections have been found
            for(ZFSInt m=0; m<3; ++m) {
              newWindow->end1[m]=end1[m];
              newWindow->end2[m]=end2[m];
            }
            goGo = true;
          }

          for(ZFSInt i=0; i<3; ++i) {
            start1[i]=newWindow->start1[i];
            end1[i]=newWindow->end1[i];
            start2[i]=newWindow->start2[i];
            end2[i]=newWindow->end2[i];
          }
          start1[countDim]  = start1[countDim]  - newWindow->step1[countDim];
          start2[countDim2] = start2[countDim2] - newWindow->step2[countDim2];
          end1[countDim] = start1[countDim];
          end2[countDim2] = start2[countDim2];

          pos2[newWindow->order[2]] = start2[newWindow->order[2]];
          kk=start1[2];
          do {
            pos2[newWindow->order[1]] = start2[newWindow->order[1]];
            jj=start1[1];
            do {
              pos2[newWindow->order[0]] = start2[newWindow->order[0]];
              ii=start1[0];
              do {
                pos1[0]=ii;pos1[1]=jj;pos1[2]=kk;
                if (newWindow->BC==6000) {
                  labell=true;
                } else {
                  labell=false;
                }

                connectionNode temp2(newWindow->BC,newWindow->Id1,pos1,newWindow->Id2,pos2,labell);
                found=findConnection(temp2);

                if(!found&&domainId()==0&&newWindow->BC==4401&&newWindow->Id1==2&&newWindow->Id2!=2) {
                  connectionNode temp3(newWindow->BC,newWindow->Id2,pos2,newWindow->Id1,pos1,labell);
                  found=findConnection(temp3);
                  if(found) {
                    cout<<"wooooo!! somthing is  wrong and i do not know where and why!!!!!!"<<endl;
                  }
                }
                if (!found) {
                  break;
                }
                pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
                ii=ii+newWindow->step1[0];

              } while(ii>=start1[0]&&ii<=end1[0]);
              pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
              jj=jj+newWindow->step1[1];

            } while(jj>=start1[1]&&jj<=end1[1]);
            pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
            kk=kk+newWindow->step1[2];

          } while(kk>=start1[2]&&kk<=end1[2]);


          if (found) {
            //all connections have been found
            for(ZFSInt m=0; m<3; ++m) {
              newWindow->start1[m]=start1[m];
              newWindow->start2[m]=start2[m];
            }
            goGo = true;
          }
        }
      }
    }

    if (newWindow->BC == 6000 && newWindow->Id2 < newWindow->Id1) {
      mapInvert1(newWindow);
    }
    mapNormalize3(newWindow);

    //delete treated connections
    pos2[newWindow->order[2]] = newWindow->start2[newWindow->order[2]];
    for(ZFSInt k = newWindow->start1[2];k<=newWindow->end1[2]; k=k+newWindow->step1[2]) {
      pos2[newWindow->order[1]] = newWindow->start2[newWindow->order[1]];
      for(ZFSInt j = newWindow->start1[1]; j<=newWindow->end1[1]; j=j+newWindow->step1[1]) {
        pos2[newWindow->order[0]] = newWindow->start2[newWindow->order[0]];
        for(ZFSInt i = newWindow->start1[0]; i<=newWindow->end1[0]; i=i+newWindow->step1[0]) {
          pos1[0]=i;pos1[1]=j;pos1[2]=k;

          if (newWindow->BC==6000) {
            labell=true;
          } else {
            labell=false;
          }

          connectionNode temp2(newWindow->BC,newWindow->Id1,pos1,newWindow->Id2,pos2,labell);
          found=findConnection(temp2);

          if (!found) {
            cout<<"error!! can not delete the element!!!"<<endl;
            cout<<"BC: "<<newWindow->BC<<" id1:"<<newWindow->Id1<<" id1:"<<newWindow->Id2<<endl;
            cout<<"i: "<<i<<" iend:"<<newWindow->end1[0]<<" j: "<<j
                <<" jend:"<<newWindow->end1[1]<<" k: "<<k<<" kend:"<<newWindow->end1[2]<<endl;
          }

          if(found)  {
            removeConnection(temp2);
            countConnection = countConnection + 1;
          }

          pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
        }
        pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
      }
      pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
    }

    ZFSInt facecount;
    if (!mapCheck(newWindow)) {
      cout<<"invalid window!"<<endl;
    }

    facecount=0;
    for(ZFSInt i =0;i<3; ++i) {
      if(newWindow->end1[i]-newWindow->start1[i]!=0) {
        facecount++;
      }
    }

    switch(facecount)
    {
    case 2:
    {
      window2d.push_back(newWindow);
      numWindows[2]++;
    }
    break;

    case 1:
    {
      window1d.push_back(newWindow);
      numWindows[1]++;
    }
    break;

    case 0:
    {
      window0d.push_back(newWindow);
      numWindows[0]++;
    }
    break;

    default:
      cout<<"ERROR!!! face dim is wrong!!!"<<endl;
    }
  }

  zfs_log << "Connection Statistics: " << numWindows[2]<<" 2D-, "<< numWindows[1]
          << " 1D-, and "<< numWindows[0]<< " 0D-connections found"<<endl;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::singularityAssembling()
{

  // cout<<"number of connections: "<<numConnections<<endl;
  ZFSInt countConnection=0,numConnections=0,Nstar;
  ZFSStrctrdWindowMap* newWindow;
  ZFSBool found,notexisted,labell;
  connectionNode temp1;
  ZFSInt numWindows[3]={0,0,0};
  numConnections=singularconnectionset.size();

  set<connectionNode>::iterator it;

  while(singularconnectionset.size()!=0) {
    auto element=singularconnectionset.begin();
    ZFSInt order[3]={-1,-1,-1};
    ZFSInt step1[3]={1,1,1};
    ZFSInt step2[3]={1,1,1};
    ZFSInt pos1[3],pos2[3],b1,b2;

    pos1[0]=element->pos1[0];
    pos1[1]=element->pos1[1];
    pos1[2]=element->pos1[2];
    pos2[0]=element->pos2[0];
    pos2[1]=element->pos2[1];
    pos2[2]=element->pos2[2];
    b1=element->blockId1;
    b2=element->blockId2;
    Nstar=element->Nstar;

    newWindow=new ZFSStrctrdWindowMap(nDim);
    mapCreate(b1, pos1, pos1, step1,
              b2, pos2, pos2, step2,
              order, element->BC, newWindow);
    newWindow->Nstar=Nstar;

    for (ZFSInt i=0;i<3;++i) {
      found=false;
      pos1[0]=newWindow->start1[0];
      pos1[1]=newWindow->start1[1];
      pos1[2]=newWindow->start1[2];

      if(newWindow->order[i]==-1) {
        //D1+

        pos1[i]=newWindow->start1[i] + 1;
        for (ZFSInt j=0;j<3;++j) {
          //D2
          pos2[0]=newWindow->start2[0];
          pos2[1]=newWindow->start2[1];
          pos2[2]=newWindow->start2[2];

          notexisted=true;
          for(ZFSInt k=0;k<3;++k) {
            if(newWindow->order[k]==j) {
              notexisted=false;
            }
          }

          if(notexisted) {

            pos2[j]=newWindow->start2[j]+1;

            temp1.BC=newWindow->BC;
            temp1.blockId1=newWindow->Id1; temp1.blockId2=newWindow->Id2;
            temp1.pos1[0]=pos1[0];temp1.pos1[1]=pos1[1];temp1.pos1[2]=pos1[2];
            temp1.pos2[0]=pos2[0];temp1.pos2[1]=pos2[1];temp1.pos2[2]=pos2[2];
            temp1.Nstar=Nstar;
            found=findConnection(temp1,Nstar);

            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }

            //D2-
            pos2[j]=newWindow->start2[j]-1;
            // temp1(newWindow->BC,newWindow->Id1,pos1,newWindow->Id2,pos2);
            temp1.BC=newWindow->BC;
            temp1.blockId1=newWindow->Id1;
            temp1.blockId2=newWindow->Id2;
            temp1.pos1[0]=pos1[0];temp1.pos1[1]=pos1[1];
            temp1.pos1[2]=pos1[2];temp1.pos2[0]=pos2[0];
            temp1.pos2[1]=pos2[1];temp1.pos2[2]=pos2[2];
            temp1.Nstar=Nstar;
            found=findConnection(temp1,Nstar);

            if(found) {
              newWindow->step1[i] = 1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }
        if(found) {
          continue;
        }

        //D1-
        pos1[i]=newWindow->start1[i] - 1;
        for (ZFSInt j=0;j<3;++j) {
          //D2+
          pos2[0]=newWindow->start2[0];
          pos2[1]=newWindow->start2[1];
          pos2[2]=newWindow->start2[2];
          notexisted=true;

          for(ZFSInt k=0;k<3;++k) {
            if(newWindow->order[k]==j)notexisted=false;
          }

          if(notexisted) {
            pos2[j]=newWindow->start2[j]+1;
            temp1.BC=newWindow->BC;
            temp1.blockId1=newWindow->Id1; temp1.blockId2=newWindow->Id2;
            temp1.pos1[0]=pos1[0];
            temp1.pos1[1]=pos1[1];
            temp1.pos1[2]=pos1[2];
            temp1.pos2[0]=pos2[0];
            temp1.pos2[1]=pos2[1];
            temp1.pos2[2]=pos2[2];
            temp1.Nstar=Nstar;
            found=findConnection(temp1,Nstar);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = 1;
              newWindow->order[i] = j;
              break;
            }
            //D2-
            pos2[j]=newWindow->start2[j]-1;
            temp1.BC=newWindow->BC;
            temp1.blockId1=newWindow->Id1; temp1.blockId2=newWindow->Id2;
            temp1.pos1[0]=pos1[0];
            temp1.pos1[1]=pos1[1];
            temp1.pos1[2]=pos1[2];
            temp1.pos2[0]=pos2[0];
            temp1.pos2[1]=pos2[1];
            temp1.pos2[2]=pos2[2];
            temp1.Nstar=Nstar;
            found=findConnection(temp1,Nstar);

            if(found) {
              newWindow->step1[i] = -1;
              newWindow->step2[j] = -1;
              newWindow->order[i] = j;
              break;
            }
          }
        }
        if(found) {
          continue;
        }
      }
    }

    ZFSInt ordercount=0;
    ZFSBool facewindow,directionInWindow[3];

    for (ZFSInt i=0;i<3;++i) {
      directionInWindow[i]=false;
      if(newWindow->order[i]!=-1) {
        ordercount++;directionInWindow[i]=true;
      }
    }

    if(ordercount>2) {
      cout<< "Invalid volume mapping found! "
          << "Are your blocks overlapping or is the grid epsilon too large?"<<endl;
    }

    if(ordercount==2) {
      facewindow=true;
    } else {
      facewindow=false;
    }

    for (ZFSInt i=0;i<3;++i) {
      ZFSInt j=0;
      if(newWindow->order[i]==-1) {
        for(j=0;j<3;++j) {
          labell=true;
          for(ZFSInt k=0;k<3;++k) {
            if ( newWindow->order[k]==j ) {
              labell=false;
            }
          }
          if(labell==true) {
            newWindow->order[i]=j;break;
          }
        }

        if(facewindow) {
          if(newWindow->start1[i]==0) {
            newWindow->dc1=i+1;
          } else {
            newWindow->dc1=-i-1;
          }

          if(newWindow->start2[j]==0) {
            newWindow->dc2=j+1;
          } else {
            newWindow->dc2=-j-1;
          }
        } else {
          newWindow->dc1= 999;
          newWindow->dc2= 999;
        }
      }
    }

    ZFSInt start1[3],end1[3],start2[3],end2[3];
    ZFSBool goGo=true;
    ZFSInt countDim,countDim2;
    ZFSInt ii,jj,kk;

    while (goGo) {
      goGo=false;
      for ( countDim=0; countDim<3; ++countDim) {
        if(directionInWindow[countDim]) {
          countDim2=newWindow->order[countDim];
          for(ZFSInt i=0; i<3; ++i) {
            start1[i]=newWindow->start1[i];
            end1[i]=newWindow->end1[i];
            start2[i]=newWindow->start2[i];
            end2[i]=newWindow->end2[i];
          }
          end1[countDim]  = end1[countDim]  + newWindow->step1[countDim];
          end2[countDim2] = end2[countDim2] + newWindow->step2[countDim2];
          start1[countDim] = end1[countDim];
          start2[countDim2] = end2[countDim2];

          pos2[newWindow->order[2]] = start2[newWindow->order[2]];
          kk=start1[2];

          do {
            pos2[newWindow->order[1]] = start2[newWindow->order[1]];
            jj=start1[1];
            do {
              pos2[newWindow->order[0]] = start2[newWindow->order[0]];
              ii=start1[0];
              do {
                pos1[0]=ii;pos1[1]=jj;pos1[2]=kk;
                if (newWindow->BC==6000) {
                  labell=true;
                } else {
                  labell=false;
                }

                connectionNode temp2(newWindow->BC,newWindow->Id1,pos1,newWindow->Id2,pos2,labell);
                temp2.Nstar=Nstar;
                found=findConnection(temp2,Nstar);

                if (!found) {
                  break;
                }
                pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
                ii=ii+newWindow->step1[0];

              } while(ii>=start1[0]&&ii<=end1[0]);
              pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
              jj=jj+newWindow->step1[1];

            } while(jj>=start1[1]&&jj<=end1[1]);
            pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
            kk=kk+newWindow->step1[2];

          } while(kk>=start1[2]&&kk<=end1[2]);

          if (found) {
            //all connections have been found
            for(ZFSInt m=0; m<3; ++m) {
              newWindow->end1[m]=end1[m];
              newWindow->end2[m]=end2[m];
            }
            goGo = true;
          }

          for(ZFSInt i=0; i<3; ++i) {
            start1[i]=newWindow->start1[i];
            end1[i]=newWindow->end1[i];
            start2[i]=newWindow->start2[i];
            end2[i]=newWindow->end2[i];
          }

          start1[countDim]  = start1[countDim]  - newWindow->step1[countDim];
          start2[countDim2] = start2[countDim2] - newWindow->step2[countDim2];
          end1[countDim] = start1[countDim];
          end2[countDim2] = start2[countDim2];

          pos2[newWindow->order[2]] = start2[newWindow->order[2]];
          kk=start1[2];

          do {
            pos2[newWindow->order[1]] = start2[newWindow->order[1]];
            jj=start1[1];
            do {
              pos2[newWindow->order[0]] = start2[newWindow->order[0]];
              ii=start1[0];
              do {
                pos1[0]=ii;pos1[1]=jj;pos1[2]=kk;

                if (newWindow->BC==6000) {
                  labell=true;
                } else {
                  labell=false;
                }

                connectionNode temp2(newWindow->BC,newWindow->Id1,pos1,newWindow->Id2,pos2,labell);
                temp2.Nstar=Nstar;
                found=findConnection(temp2,Nstar);

                if (!found) {
                  break;
                }

                pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
                ii=ii+newWindow->step1[0];

              } while(ii>=start1[0]&&ii<=end1[0]);
              pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
              jj=jj+newWindow->step1[1];

            } while(jj>=start1[1]&&jj<=end1[1]);
            pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
            kk=kk+newWindow->step1[2];

          } while(kk>=start1[2]&&kk<=end1[2]);

          if (found) {
            //all connections have been found
            for(ZFSInt m=0; m<3; ++m) {
              newWindow->start1[m]=start1[m];
              newWindow->start2[m]=start2[m];
            }

            goGo = true;
          }
        }
      }
    }

    if (newWindow->BC == 6000 && newWindow->Id2 < newWindow->Id1) {
      mapInvert1(newWindow);
    }
    mapNormalize3(newWindow);

    //delete treated connections
    pos2[newWindow->order[2]] = newWindow->start2[newWindow->order[2]];
    for(ZFSInt k = newWindow->start1[2];k<=newWindow->end1[2]; k=k+newWindow->step1[2]) {
      pos2[newWindow->order[1]] = newWindow->start2[newWindow->order[1]];
      for(ZFSInt j = newWindow->start1[1]; j<=newWindow->end1[1]; j=j+newWindow->step1[1]) {
        pos2[newWindow->order[0]] = newWindow->start2[newWindow->order[0]];
        for(ZFSInt i = newWindow->start1[0]; i<=newWindow->end1[0]; i=i+newWindow->step1[0]) {
          pos1[0]=i;pos1[1]=j;pos1[2]=k;

          if (newWindow->BC==6000) {
            labell=true;
          } else {
            labell=false;
          }

          connectionNode temp2(newWindow->BC,newWindow->Id1,pos1,newWindow->Id2,pos2,labell);
          temp2.Nstar=Nstar;
          found=findConnection(temp2,Nstar);

          if (!found) {
            cout<<"singular  error!! can not delete the element!!!"<<endl;
          }

          removeConnection(temp2,Nstar);
          if (found) {
            countConnection = countConnection + 1;
          }
          pos2[newWindow->order[0]] = pos2[newWindow->order[0]] + newWindow->step2[newWindow->order[0]];
        }
        pos2[newWindow->order[1]] = pos2[newWindow->order[1]] + newWindow->step2[newWindow->order[1]];
      }
      pos2[newWindow->order[2]] = pos2[newWindow->order[2]] + newWindow->step2[newWindow->order[2]];
    }

    //set Id for singularmap
    newWindow->SingularId=numWindows[0];
    singularwindow.push_back(newWindow);
    numWindows[0]++;

    if(domainId()==0)  {
      cout<< numWindows[0]<< " singular connections found ("<< countConnection*100/numConnections<< "% done)"<<endl;
    }
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::periodicPointsChange(ZFSFloat* pt,
                                                          ZFSId type,
                                                          ZFSFloat* periodicDisplacements)
{
  ZFSFloat angle=0;
  /*
    case 4401 4402: first periodic direction
    case 4403 4404: second periodic direction
    case 4405 4405: third periodic direction
    case 4011 4012:  rotation X axis clockwise and anticlockwise
  */
  ZFSFloat rotationMatrix[3][3];
  ZFSFloat tmp[3];
  switch(type)
  {
  case 4401:
  case 4403:
  case 4405:
  {
    const ZFSId displacementId = (ZFSFloat)(type-4400+1)/2.0 - 1;
    for (ZFSInt dim=0;dim<3;++dim) {
      pt[dim]=pt[dim]+periodicDisplacements[dim*nDim+displacementId];
    }
    break;
  }

  case 4402:
  case 4404:
  case 4406:
  {
    const ZFSId displacementId = (ZFSFloat)(type-4400)/2.0 - 1;
    for (ZFSInt dim=0;dim<3;++dim) {
      pt[dim]=pt[dim]-periodicDisplacements[dim*nDim+displacementId];
    }
    break;
  }

  case 4011:
  {
    rotationMatrix[1][1]=cos(angle);
    rotationMatrix[1][2]=-sin(angle);
    rotationMatrix[2][1]=sin(angle);
    rotationMatrix[2][2]=cos(angle);
    tmp[0]=pt[0];
    tmp[1]=pt[1];
    tmp[2]=pt[2];
    for (ZFSInt i=0;i<3;++i)
    {
      pt[i]=rotationMatrix[i][0]*tmp[0]+rotationMatrix[i][1]*tmp[1]
        +rotationMatrix[i][2]*tmp[2];
    }
    break;
  }

  case 4012:
  {
    rotationMatrix[1][1]=cos(-angle);
    rotationMatrix[1][2]=-sin(-angle);
    rotationMatrix[2][1]=sin(-angle);
    rotationMatrix[2][2]=cos(-angle);
    tmp[0]=pt[0];
    tmp[1]=pt[1];
    tmp[2]=pt[2];
    for (ZFSInt i=0;i<3;++i)
    {
      pt[i]=rotationMatrix[i][0]*tmp[0]+rotationMatrix[i][1]*tmp[1]
        +rotationMatrix[i][2]*tmp[2];
    }
    break;
  }

  default:
  {
    cout<<"ERROR!!! periodic type is wrong!!! in windowinfoe"<<endl;
  }
  }
}

template <ZFSInt nDim>
ZFSBool ZFSStrctrdBlckWindowInfo<nDim>::addConnection(ZFSInt connectiontype,
                                                      ZFSInt b1, ZFSInt* p1,
                                                      ZFSInt b2, ZFSInt* p2)
{
  //add connections for 6000 and 6000+
  pair<set<connectionNode>::iterator,ZFSBool> ret;

  if(connectiontype==6000) {
    connectionNode a(connectiontype,b1,p1,b2,p2,true);
    ret=connectionset.insert(a);
  } else {
    connectionNode a(connectiontype,b1,p1,b2,p2,false);
    ret=connectionset.insert(a);
  }

  if(ret.second==false) {
    cout<<"error! same connection can not be added! check the connections!"<<endl;
    return false;
  }

  return true;
}

template <ZFSInt nDim>
ZFSBool ZFSStrctrdBlckWindowInfo<nDim>::addConnection(ZFSInt connectiontype,
                                                      ZFSInt b1, ZFSInt* p1,
                                                      ZFSInt b2, ZFSInt* p2,
                                                      ZFSInt Nstar)
{
  //add special connections (3 or 5-star)
  pair<set<connectionNode>::iterator,ZFSBool> ret;

  if(connectiontype==6000) {
    connectionNode a(connectiontype,b1,p1,b2,p2,true);
    a.Nstar=Nstar;
    ret=singularconnectionset.insert(a);
  } else if(connectiontype>=4000&&connectiontype<5000 ) {
    //periodic boundary condition only need to handle 5 star communications
    connectionNode a(connectiontype,b1,p1,b2,p2,false);
    a.Nstar=Nstar;
    ret=singularconnectionset.insert(a);
  } else {
    cout<<"singular point on BC "<<connectiontype<<" is not supported!! check it!!!"<<endl;
    exit(0);
  }

  if(ret.second==false) {
    cout<<"error! same connection can not be added! check the connections!"<<endl;
    return false;
  }

  return true;
}


template <ZFSInt nDim>
ZFSBool ZFSStrctrdBlckWindowInfo<nDim>::findConnection(connectionNode a)
{
  set<connectionNode>::iterator it;
  it=connectionset.find(a);

  if(it!=connectionset.end()) {
    return true;
  } else {
    return false;
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::removeConnection(connectionNode a)
{
  connectionset.erase(a);
}


template <ZFSInt nDim>
ZFSBool ZFSStrctrdBlckWindowInfo<nDim>::findConnection(connectionNode a, ZFSInt Nstar)
{
  set<connectionNode>::iterator it;
  a.Nstar=Nstar;
  it=singularconnectionset.find(a);

  if(it!=singularconnectionset.end()) {
    return true;
  } else {
    return false;
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::removeConnection(connectionNode a, ZFSInt Nstar)
{
  a.Nstar=Nstar;
  singularconnectionset.erase(a);
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::readNewWindowInfo(ZFSString gridFileName)
{
  zfs_log << "-->reading windowInformation from file ...(serial->parallel approach)" << endl;
  if(globalDomainId()==0) {
    ZFSInt fId = io_openfile("hdf5", gridFileName.c_str(), "collective", MPI_COMM_SELF);
    stringstream dummy1;
    //go to the windowInformation folder and get number of windows
    noInputWindowInformation=0;
    io_read_iattribute1(fId, "/WindowInformation", "noWindows", &noInputWindowInformation);
    //io_read_iattribute1(file_id, "/Connectivity", "noConnec", &noInputWindowConnections);
    noInputWindowConnections = 0;
    io_read_iattribute1(fId, "/BC", "noBndryCnds",&noInputBndryCnds);

    ZFSIntScratchSpace windowInfo(noInputWindowInformation*7, __CALLING_FUNCTION__, "widnow information");

    ZFSInt indices[MAX_SPACE_DIMENSIONS][2]; //stores the indices of the windowInformation
    //copy the data of the windowInformations into a Scratspace which is send to everyone
    for(ZFSInt i=0; i<noInputWindowInformation; i++){
      stringstream windowdummy;
      windowdummy << (i+1);
      ZFSString window = "window" + windowdummy.str();
      ZFSString window1 = "/WindowInformation/window" + windowdummy.str() ;
      io_read_iattribute1(fId, window1.c_str(), "block", &(windowInfo[i*7]));
      int size[2]={nDim ,2};
      io_read_idataset1(fId, "WindowInformation",window.c_str(),2,size, &(windowInfo[i*7+1]));
    }

    ZFSInt noBcWindows;
    ZFSInt** BCInformation = new ZFSInt *[noInputBndryCnds];
    ZFSInt* bcCndId = new ZFSInt[noInputBndryCnds];
    ZFSInt* noBcCndWindows = new ZFSInt[noInputBndryCnds];

    for(ZFSInt i =0; i<noInputBndryCnds; i++) {
      stringstream bcDummy;
      stringstream bcDummy1;
      bcDummy<<"/BC/BC" << (i+1);
      bcDummy1<< "BC" << (i+1);
      char bcIdWindows[80];
      char bcId[80];
      strcpy(bcIdWindows, &bcDummy.str()[0]);
      strcpy(bcId, &bcDummy.str()[0]);
      strcat(bcIdWindows, "/noWindows");
      io_read_iattribute1(file_id, (bcDummy.str()).c_str(), "noWindows",  &noBcWindows);
      noBcCndWindows[i]=noBcWindows;
      BCInformation[i] = new ZFSInt[noBcWindows];
      io_read_idataset1(fId, "BC", (bcDummy1.str()).c_str(),1,&noBcWindows, &BCInformation[i][0]);
      strcat(bcId, "/BC");
      io_read_iattribute1(fId, (bcDummy.str()).c_str(),"BC", &bcCndId[i]);
    }
    io_closefile(fId);

    ZFSInt noWindows=0;
    for(ZFSInt i =0; i<noInputBndryCnds; i++) { noWindows+=noBcCndWindows[i];}
    //put everything into a vector to send the data
    ZFSId count=0;
    ZFSIntScratchSpace bcInfo(noWindows+3*noInputBndryCnds, __CALLING_FUNCTION__, "send bc info");
    for(ZFSInt i =0; i<noInputBndryCnds; i++) {
      bcInfo[count]=bcCndId[i];
      count++;
      bcInfo[count]=noBcCndWindows[i];
      count++;
      memcpy(&(bcInfo[count]), &(BCInformation[i][0]), noBcCndWindows[i]*sizeof(ZFSInt));
      count+=noBcCndWindows[i];
      bcInfo[count]=bcCndId[i];
      count++;
    }

    ZFSInt countInfo[3]={noInputWindowInformation, noInputBndryCnds, noWindows};
    MPI_Bcast(countInfo, 3, MPI_INT, 0, m_zfsStrctrdComm);

    //broadcast the information
    MPI_Bcast(windowInfo.getPointer(),noInputWindowInformation*7, MPI_INT, 0, m_zfsStrctrdComm);
    MPI_Bcast(bcInfo.getPointer(), (noWindows+3*noInputBndryCnds), MPI_INT, 0, m_zfsStrctrdComm);

    //put the window information into the right place
    for(ZFSInt i=0; i<noInputWindowInformation; i++) {
      windowInformation* temp  = new windowInformation(nDim);
      //fill the windowInformation with life
      //first attribute the windowId to the object
      temp->windowId=i;
      temp->inputBlockId=windowInfo[i*7];
      memcpy(&indices[0][0], &(windowInfo[i*7+1]), MAX_SPACE_DIMENSIONS*2*sizeof(ZFSInt));
      for( int dim=0;dim<nDim;dim++) {
        temp->startindex[dim]=indices[dim][0]-1;
        temp->endindex[dim]=indices[dim][1]-1;
      }
      //set the conneciton to NULL
      temp->connec=NULL;
      //add the window to the vector
      inputWindows.push_back(temp);
    }
    //distribute the boundary conditions to the windows!
    for(ZFSInt i=0; i<noInputBndryCnds; ++i) {
      for(ZFSInt j=0; j<noBcCndWindows[i]; ++j) {
        ZFSInt windowId=BCInformation[i][j]-1; //correct the fortran notation
        inputWindows[windowId]->BC=bcCndId[i];
      }
    }

    delete[] bcCndId;
    delete[] noBcCndWindows;
    for(ZFSInt i =0; i<noInputBndryCnds; i++) {
      delete[] BCInformation[i];
    }
    delete[] BCInformation;
  }else{
    //receive the data first
    ZFSInt countInfo[3]={0, 0,0};
    MPI_Bcast(countInfo, 3, MPI_INT, 0, m_zfsStrctrdComm);
    noInputWindowInformation=countInfo[0];
    noInputBndryCnds=countInfo[1];
    ZFSInt noWindows=countInfo[2];
    ZFSIntScratchSpace bcInfo(noWindows+3*noInputBndryCnds, __CALLING_FUNCTION__, "send bc info");
    ZFSIntScratchSpace windowInfo(noInputWindowInformation*7, __CALLING_FUNCTION__, "widnow information");
    //receive the data from the other processes
    MPI_Bcast(windowInfo.getPointer(),noInputWindowInformation*7, MPI_INT, 0, m_zfsStrctrdComm);
    MPI_Bcast(bcInfo.getPointer(), (noWindows+3*noInputBndryCnds), MPI_INT, 0, m_zfsStrctrdComm);
    ZFSInt indices[MAX_SPACE_DIMENSIONS][2]; //stores the indices of the windowInformation
    //put the window information into the right place
    for(ZFSInt i=0; i<noInputWindowInformation; i++) {
      windowInformation* temp  = new windowInformation(nDim);
      //fill the windowInformation with life
      //first attribute the windowId to the object
      temp->windowId=i;
      temp->inputBlockId=windowInfo[i*7];
      memcpy(&indices[0][0], &(windowInfo[i*7+1]), MAX_SPACE_DIMENSIONS*2*sizeof(ZFSInt));
      for( int dim=0;dim<nDim;dim++) {
        temp->startindex[dim]=indices[dim][0]-1;
        temp->endindex[dim]=indices[dim][1]-1;
      }
      //set the conneciton to NULL
      temp->connec=NULL;
      //add the window to the vector
      inputWindows.push_back(temp);
    }

    //put the boundary condition info into the right place
    ZFSInt** BCInformation = new ZFSInt *[noInputBndryCnds];
    ZFSInt* bcCndId = new ZFSInt[noInputBndryCnds];
    ZFSInt* noBcCndWindows = new ZFSInt[noInputBndryCnds];
    ZFSInt count=0;
    for(ZFSInt i =0; i<noInputBndryCnds; i++) {
      bcCndId[i]=bcInfo[count];
      count++;
      noBcCndWindows[i]=bcInfo[count];
      count++;
      BCInformation[i] = new ZFSInt[noBcCndWindows[i]];
      memcpy(&(BCInformation[i][0]), &(bcInfo[count]),noBcCndWindows[i]*sizeof(ZFSInt));
      count+=noBcCndWindows[i];
      bcCndId[i]=bcInfo[count];
      count++;
    }

    //distribute the boundary conditions to the windows!
    for(ZFSInt i=0; i<noInputBndryCnds; ++i) {
      for(ZFSInt j=0; j<noBcCndWindows[i]; ++j) {
        ZFSInt windowId=BCInformation[i][j]-1; //correct the fortran notation
        inputWindows[windowId]->BC=bcCndId[i];
      }
    }

    delete[] bcCndId;
    delete[] noBcCndWindows;
    for(ZFSInt i =0; i<noInputBndryCnds; i++) {
      delete[] BCInformation[i];
    }
    delete[] BCInformation;

  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::readNewWindowInfo()
{
  zfs_log << "-->reading windowInformation from file ..." << endl;
  stringstream dummy1;

  //go to the windowInformation folder and get number of windows
  noInputWindowInformation=0;
  io_read_iattribute1(file_id, "/WindowInformation", "noWindows", &noInputWindowInformation);
  //io_read_iattribute1(file_id, "/Connectivity", "noConnec", &noInputWindowConnections);
  noInputWindowConnections = 0;
  io_read_iattribute1(file_id, "/BC", "noBndryCnds",&noInputBndryCnds);

  int indices[MAX_SPACE_DIMENSIONS][2]; //stores the indices of the windowInformation

  //first read in the windowInformations
  for(ZFSInt i=0; i<noInputWindowInformation; i++) {

    windowInformation* temp  = new windowInformation(nDim);
    //fill the windowInformation with life
    //first attribute the windowId to the object
    temp->windowId=i;

    //create path to the dataset of the windowInformation
    stringstream windowdummy;
    stringstream windowdummy1;

    windowdummy1 << "window" << (i+1);
    windowdummy << "/WindowInformation/window" << (i+1) ;
    //read the number of the inputBox
    char blockId[80];
    strcpy(blockId, &windowdummy.str()[0]);
    strcat(blockId, "/block");
    string windoMoep = windowdummy.str();
    io_read_iattribute1(file_id, windoMoep.c_str(), "block", &(temp->inputBlockId));

    //read the indices of the windowInformation
    char windowId[80];
    strcpy(windowId, &windowdummy.str()[0]);
    int size[2]={nDim ,2};
    io_read_idataset1(file_id, "WindowInformation",(windowdummy1.str()).c_str(),2,size, &indices[0][0]);

    for( int dim=0;dim<nDim;dim++) {
      temp->startindex[dim]=indices[dim][0]-1;
      temp->endindex[dim]=indices[dim][1]-1;
    }

    //set the conneciton to NULL
    temp->connec=NULL;
    //add the window to the vector
    inputWindows.push_back(temp);
  }

  //check the connectivities of the windows
  ZFSInt  connecIndices[5];
  for(ZFSId i=0; i<noInputWindowConnections; i++) {
    stringstream connecdummy;
    stringstream connecdummy1;
    connecdummy << "/Connectivity/connec" << (i+1) ;
    connecdummy1 << "connec" << (i+1);
    char connecId[80];
    strcpy(connecId, &connecdummy.str()[0]);
    int size=nDim+2;

    //check dimension in property file and grid
    ZFSInt dummySize;
    io_getDatasetSize(file_id,"Connectivity", (connecdummy1.str()).c_str(), 1, &dummySize);
    if((dummySize==4 && nDim==3) || (dummySize==5 && nDim==2)){
      zfsTerm(1, __CALLING_FUNCTION__, "nDim in properties does not fit the grid!! Please check");
    }


    io_read_idataset1(file_id, "Connectivity",(connecdummy1.str()).c_str(),1, &size, &connecIndices[0]);

    inputWindows[connecIndices[0]-1]->connec = new connectivity;
    inputWindows[connecIndices[0]-1]->connec->permutation = new ZFSInt[nDim];
    inputWindows[connecIndices[1]-1]->connec = new connectivity;
    inputWindows[connecIndices[1]-1]->connec->permutation = new ZFSInt[nDim];
    inputWindows[connecIndices[0]-1]->connec->order = new ZFSInt[nDim];
    inputWindows[connecIndices[1]-1]->connec->order = new ZFSInt[nDim];
    inputWindows[connecIndices[0]-1]->connec->window = connecIndices[1]-1;
    inputWindows[connecIndices[1]-1]->connec->window = connecIndices[0]-1;

    for(ZFSId j=2; j<5; j++) {
      if(connecIndices[j]>=0) {
        inputWindows[connecIndices[0]-1]->connec->order[j-2]=connecIndices[j]-1;
        inputWindows[connecIndices[1]-1]->connec->order[j-2]=connecIndices[j]-1;
        inputWindows[connecIndices[0]-1]->connec->permutation[j-2]=1;
        inputWindows[connecIndices[1]-1]->connec->permutation[j-2]=1;

        if(connecIndices[j]==0) {
          inputWindows[connecIndices[0]-1]->connec->order[j-2]=j-2;
          inputWindows[connecIndices[1]-1]->connec->order[j-2]=j-2;
        }
      } else {
        inputWindows[connecIndices[0]-1]->connec->order[j-2]=(-connecIndices[j])-1;
        inputWindows[connecIndices[1]-1]->connec->order[j-2]=(-connecIndices[j])-1;
        inputWindows[connecIndices[0]-1]->connec->permutation[j-2]=-1;
        inputWindows[connecIndices[1]-1]->connec->permutation[j-2]=-1;
      }
    }
  }

  //read in the boundary conditions and attribute them to the windows;
  ZFSInt noBcWindows;
  ZFSInt** BCInformation = new ZFSInt *[noInputBndryCnds];
  ZFSInt* bcCndId = new ZFSInt[noInputBndryCnds];
  ZFSInt* noBcCndWindows = new ZFSInt[noInputBndryCnds];

  for(ZFSInt i =0; i<noInputBndryCnds; i++) {
    stringstream bcDummy;
    stringstream bcDummy1;
    bcDummy<<"/BC/BC" << (i+1);
    bcDummy1<< "BC" << (i+1);
    char bcIdWindows[80];
    char bcId[80];
    strcpy(bcIdWindows, &bcDummy.str()[0]);
    strcpy(bcId, &bcDummy.str()[0]);
    strcat(bcIdWindows, "/noWindows");
    io_read_iattribute1(file_id, (bcDummy.str()).c_str(), "noWindows",  &noBcWindows);
    noBcCndWindows[i]=noBcWindows;
    BCInformation[i] = new ZFSInt[noBcWindows];
    io_read_idataset1(file_id, "BC", (bcDummy1.str()).c_str(),1,&noBcWindows, &BCInformation[i][0]);
    strcat(bcId, "/BC");
    io_read_iattribute1(file_id, (bcDummy.str()).c_str(),"BC", &bcCndId[i]);
  }

  //distribute the boundary conditions to the windows!
  for(ZFSInt i=0; i<noInputBndryCnds; ++i) {
    for(ZFSInt j=0; j<noBcCndWindows[i]; ++j) {
      ZFSInt windowId=BCInformation[i][j]-1; //correct the fortran notation
      inputWindows[windowId]->BC=bcCndId[i];
    }
  }

  //delete the pointers for the moment. if needed for later make a member variable out of it
  delete[] bcCndId;
  delete[] noBcCndWindows;
  for(ZFSInt i =0; i<noInputBndryCnds; i++) {
    delete[] BCInformation[i];
  }
  delete[] BCInformation;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::initGlobals()
{
  ZFSStrctrdWindowMap* windowMap;
  windowInformation* windowInfo;
  // ZFSInt noPeriodicWindows[3] = {0,0,0};
  //function creates a map of all the global boundary conditions
  //global means from Input file !!!

  //first look which periodic direction is not used
  //ZFSBool switchPeriodicIds = true;
  ZFSBool isPeriodicDirUsed[3] = {false,false,false};
  ZFSBool switchPeriodicBC[3] = {true, true, true};

  //check if periodic direction (4401,4403,4405) is already in use (find free slots for channel bc)
  //also check if distinct bc numbers (4401/4402, 4403/4404, 4405/4406) for both surfaces are used
  //or if both periodic surfaces have the same number, in that case swtich
  for(ZFSId windowId=0; windowId< noInputWindowInformation; windowId++) {
    switch(inputWindows[windowId]->BC) {
    case 4401:
      isPeriodicDirUsed[0] = true;
      break;
    case 4402:
      switchPeriodicBC[0] = false;
      break;
    case 4403:
      isPeriodicDirUsed[1] = true;
      break;
    case 4404:
      switchPeriodicBC[1] = false;
      break;
    case 4405:
      isPeriodicDirUsed[2] = true;
      break;
    case 4406:
      switchPeriodicBC[2] = false;
      break;
    default:
      //do nothing
      break;
    }
  }

  for(ZFSId noWin=0; noWin< noInputWindowInformation; noWin++) {
    windowInfo= inputWindows[noWin];
    windowMap=new ZFSStrctrdWindowMap(nDim);
    ZFSInt order[3]={0,1,2};
    ZFSInt step1[3]={1,1,1};

    ZFSInt start2[3]={0,0,0};
    ZFSInt end2[3]={0,0,0};
    ZFSInt step2[3]={1,1,1};
    if(windowInfo->connec!=NULL) {
      //if that is the case be carefull, if a channel BC or rescaling is
      //used (splitting in communication and physical BC): A
      //splitting of the maps has to be performed first for the correctness
      //of the map correction schemes.

      // connection to another window ==> check the permutation and order
      memcpy(order, windowInfo->connec->order, nDim*sizeof(ZFSInt));
      memcpy(step2, windowInfo->connec->permutation, nDim*sizeof(ZFSInt));
      ZFSInt adjwindow = windowInfo->connec->window;
      memcpy(start2, inputWindows[adjwindow]->startindex, nDim*sizeof(ZFSInt));
      memcpy(end2, inputWindows[adjwindow]->endindex, nDim*sizeof(ZFSInt));

      for(ZFSInt i=0; i<nDim; i++) {
        if(step2[i]<0) {
          ZFSInt tmp = end2[i];
          end2[i]= start2[i];
          start2[i]=tmp;
        }
      }

      mapCreate(windowInfo->inputBlockId, windowInfo->startindex, windowInfo->endindex, step1,
                inputWindows[adjwindow]->inputBlockId ,  start2, end2 , step2, order, windowInfo->BC, windowMap);
      globalStrctrdBndryCndMaps.push_back(windowMap);

    } else {
      //else use the standard mapping to itself
      if(windowInfo->BC!=6000&&(windowInfo->BC<4000||windowInfo->BC>=5000)&&windowInfo->BC!=2222&&windowInfo->BC!=2221) { //junoh
        mapCreate(windowInfo->inputBlockId, windowInfo->startindex, windowInfo->endindex, step1,
                  (windowInfo->windowId+1) ,  NULL, NULL , NULL, order,windowInfo->BC , windowMap );
        //add the map to the list!
        globalStrctrdBndryCndMaps.push_back(windowMap);
      }
            // special treatment for zonal BC 2222   junoh
     
      if(windowInfo->BC==2222) {
        mapCreate(windowInfo->inputBlockId, windowInfo->startindex, windowInfo->endindex, step1,
                  (windowInfo->windowId+1) ,  NULL, NULL , NULL, order,windowInfo->BC , windowMap );
        //add map for bc2222 to list
        globalStrctrdBndryCndMaps.push_back(windowMap);
      }
      //special treatment for zonal BC 2221
      if(windowInfo->BC==2221) {
        mapCreate(windowInfo->inputBlockId, windowInfo->startindex, windowInfo->endindex, step1,
                  (windowInfo->windowId+1) ,  NULL, NULL , NULL, order,windowInfo->BC , windowMap );
        //add map for bc2221 to list
        globalStrctrdBndryCndMaps.push_back(windowMap);
      
	//now add the same map for bc7909 to list
	windowMap=new ZFSStrctrdWindowMap(nDim); //junoh
        mapCreate(windowInfo->inputBlockId, windowInfo->startindex, windowInfo->endindex, step1,
                  (windowInfo->windowId+1) ,  NULL, NULL , NULL, order, 7909, windowMap );
        //add the map to the list!
        globalStrctrdBndryCndMaps.push_back(windowMap);
      }


      if(windowInfo->BC>4000&&windowInfo->BC<5000) {
        //add channel flow and rotating periodic BC to physical BC.
        //in order to build the comm groups for channel and rotation.
        //!!!important set numbers for channel and rotation.
        //4001 and 4002 for rotation
        //2401 and 2402 for channel

        if(windowInfo->BC==4011) {
          mapCreate(windowInfo->inputBlockId, windowInfo->startindex, windowInfo->endindex, step1,
                    (windowInfo->windowId+1) ,  NULL, NULL , NULL, order,4001 , windowMap );
          globalStrctrdBndryCndMaps.push_back(windowMap);
        }

        if(windowInfo->BC==4012){
          mapCreate(windowInfo->inputBlockId, windowInfo->startindex, windowInfo->endindex, step1,
                    (windowInfo->windowId+1) ,  NULL, NULL , NULL, order,4002 , windowMap );
          globalStrctrdBndryCndMaps.push_back(windowMap);
        }

        if(windowInfo->BC>4400&&windowInfo->BC<4407) {
          //if same BC number is used for both periodic surfaces
          //switch BC number of the surface at the end of the block
          ZFSInt periodicDirection = -1;
          for(ZFSId dim=0; dim<nDim; dim++) {
            if(windowInfo->startindex[dim] == windowInfo->endindex[dim]) {
              periodicDirection = dim;
              break;
            }
          }

          ZFSInt bigBlockId = windowInfo->inputBlockId;
          ZFSId bigBlockSize = partition->inputBoxInfo[bigBlockId]->DirLast[nDim-1-periodicDirection];
          ZFSInt periodicId = windowInfo->BC-4401;
          ZFSInt periodicIndex = periodicId/2;

          if(switchPeriodicBC[periodicIndex]&&
             bigBlockSize==windowInfo->startindex[periodicDirection]&&
             periodicId%2==0) {
            zfs_log << "Changing BC from " << windowInfo->BC << " to " << windowInfo->BC + 1 << endl;
            windowInfo->BC++;
          }
        }

        //Channel BCs
        //for BC 4005 and 4006 two periodic bcs 4401 and 4402
        //will be created and also physical bcs 2401 and 2402
        //to correct the pressure and density at both sides
        if(windowInfo->BC == 4005 || windowInfo->BC == 4006) {
          //use the first unused periodic direction
          ZFSInt freePeriodicBC = 4400;
          for(ZFSId dir=0; dir<nDim; dir++) {
            if(isPeriodicDirUsed[dir]==false) {
              freePeriodicBC += dir*2+1;
              break;
            }
          }

          zfs_log << "Using Periodic BC " << freePeriodicBC << " / "
                  << freePeriodicBC+1 << " for the channel/pipe flow" << endl;

          if(windowInfo->BC==4005) {
            zfs_log << "Identified channel inflow bc 4005, creating periodic map BC " << freePeriodicBC
                    << " and physicalMap BC " << 2401 << endl;
            windowInfo->BC=freePeriodicBC;
            mapCreate(windowInfo->inputBlockId, windowInfo->startindex, windowInfo->endindex, step1,
                      (windowInfo->windowId+1) ,  NULL, NULL , NULL, order,2401 , windowMap );
            globalStrctrdBndryCndMaps.push_back(windowMap);
          }

          if(windowInfo->BC==4006) {
            zfs_log << "Identified channel inflow bc 4006, creating periodic map BC " << freePeriodicBC+1
                    << " and physicalMap BC " << 2402 << endl;
            windowInfo->BC=freePeriodicBC+1;
            mapCreate(windowInfo->inputBlockId, windowInfo->startindex, windowInfo->endindex, step1,
                      (windowInfo->windowId+1) ,  NULL, NULL , NULL, order,2402 , windowMap );
            globalStrctrdBndryCndMaps.push_back(windowMap);
          }
        }
      }
    }
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::writeConnectionWindowInformation(ZFSFloat* periodicDisplacements)
{
  /////////////////////////////////
  /////// NORMAL CONNECTIONS //////
  /////////////////////////////////
  ZFSInt noConnectionMaps = 0;
  for(ZFSId mapId = 0; mapId<(ZFSId)globalStrctrdBndryCndMaps.size(); mapId++) {
    //only write out the connection maps
    if((globalStrctrdBndryCndMaps[mapId]->BC==6000) ||
       (globalStrctrdBndryCndMaps[mapId]->BC > 4400 && globalStrctrdBndryCndMaps[mapId]->BC < 4407)) {

      stringstream windowName;
      windowName << "map" << noConnectionMaps << endl;
      ZFSString windowNameStr = windowName.str();
      ZFSInt noWindowInformation = 7*nDim+16;
      io_create_idataset(file_id, "/Connectivity", windowNameStr.c_str(), 1, &noWindowInformation);

      ZFSIntScratchSpace windowInformation(noWindowInformation, __CALLING_FUNCTION__, "windowInformation");
      writeMapToArray(globalStrctrdBndryCndMaps[mapId],windowInformation.begin());

      ZFSInt offset = 0;
      if(domainId()==0) {
        io_write_idataset_part(file_id, "/Connectivity", windowNameStr.c_str(), 1,  &noWindowInformation,  &offset, windowInformation.begin());
      } else {
        ZFSInt zeroWindows = 0;
        io_write_idataset_part(file_id, "/Connectivity", windowNameStr.c_str(), 1,  &zeroWindows,  &offset, NULL);
      }
      noConnectionMaps++;
    }
  }

  /////////////////////////////////
  /////// SINGULARITIES ///////////
  /////////////////////////////////

  ZFSInt noSingularityMaps = 0;
  for(ZFSId mapId = 0; mapId<(ZFSId)singularwindow.size(); mapId++) {
    stringstream windowName;
    windowName << "singularMap" << noSingularityMaps << endl;
    ZFSString windowNameStr = windowName.str();
    ZFSInt noWindowInformation = 7*nDim+16;
    io_create_idataset(file_id, "/Connectivity", windowNameStr.c_str(), 1, &noWindowInformation);

    ZFSIntScratchSpace windowInformation(noWindowInformation, __CALLING_FUNCTION__, "windowInformation");
    writeMapToArray(singularwindow[mapId],windowInformation.begin());

    ZFSInt offset = 0;
    if(domainId()==0) {
      io_write_idataset_part(file_id, "/Connectivity", windowNameStr.c_str(), 1,  &noWindowInformation,  &offset, windowInformation.begin());
    } else {
      ZFSInt zeroWindows = 0;
      io_write_idataset_part(file_id, "/Connectivity", windowNameStr.c_str(), 1,  &zeroWindows,  &offset, NULL);
    }
    noSingularityMaps++;
  }

  /////////////////////////////////
  ///// PERIODIC DISPLACEMENTS ////
  /////////////////////////////////
  cout << "Starting write out " << endl;
  ZFSInt noPeriodicDisplacementInfo = nDim*nDim;
  ZFSFloatScratchSpace localPeriodicDisplacements(m_noInputBlocks*noPeriodicDisplacementInfo, __CALLING_FUNCTION__, "allPeriodicDisplacements");
  ZFSFloatScratchSpace globalPeriodicDisplacements(m_noInputBlocks*noPeriodicDisplacementInfo, __CALLING_FUNCTION__, "allPeriodicDisplacements");
  localPeriodicDisplacements.fill(-99999999.0);
  for(ZFSId periodicId=0; periodicId<noPeriodicDisplacementInfo; periodicId++) {
    localPeriodicDisplacements(m_inputBlockId*noPeriodicDisplacementInfo + periodicId) = periodicDisplacements[periodicId];
  }

  cout << "Allreduce noblock: " << m_noInputBlocks << " totalSize: " << m_noInputBlocks*noPeriodicDisplacementInfo << endl;
  MPI_Allreduce(&localPeriodicDisplacements(0), &globalPeriodicDisplacements(0), m_noInputBlocks*noPeriodicDisplacementInfo, MPI_DOUBLE, MPI_MAX, m_zfsStrctrdComm);

  cout << "Writing out periodic Window displacements" << endl;

  for(ZFSId blockId = 0; blockId<m_noInputBlocks; blockId++) {
    stringstream path;
    path << "periodicDisplacements" << blockId;
    ZFSString pathStr = path.str();
    if(!io_checkObj(file_id, "/Connectivity", pathStr.c_str())) {
      io_create_ddataset(file_id, "/Connectivity", pathStr.c_str(), 1, &noPeriodicDisplacementInfo);
    }
    ZFSInt offset = 0;
    if(domainId()==0) {
      io_write_ddataset_part(file_id, "/Connectivity", pathStr.c_str(), 1,  &noPeriodicDisplacementInfo,  &offset, &globalPeriodicDisplacements[blockId*noPeriodicDisplacementInfo]);
    } else {
      ZFSInt zeroWindows = 0;
      io_write_ddataset_part(file_id, "/Connectivity", pathStr.c_str(), 1,  &zeroWindows,  &offset, NULL);
    }
  }

  ZFSInt hasConnectionInfo = 1;
  zfs_log << "Connection info written to grid file, normalConnections: " << noConnectionMaps << " noSingularityMaps: " << noSingularityMaps << endl;

  if(!io_checkObj(file_id, "/Connectivity", "hasConnectionInfo")) { io_create_iattribute(file_id, "/Connectivity", "hasConnectionInfo",1);}
  io_write_iattribute1(file_id, "/Connectivity", "hasConnectionInfo",1, &hasConnectionInfo);

  if(!io_checkObj(file_id, "/Connectivity", "noConnections")) {io_create_iattribute(file_id, "/Connectivity", "noConnections",1);}
  io_write_iattribute1(file_id, "/Connectivity", "noConnections",1, &noConnectionMaps);

  if(!io_checkObj(file_id, "/Connectivity", "noSingularConnections")) { io_create_iattribute(file_id, "/Connectivity", "noSingularConnections",1);}
  io_write_iattribute1(file_id, "/Connectivity", "noSingularConnections",1, &noSingularityMaps);
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::readConnectionWindowInformation(ZFSFloat* periodicDisplacements, ZFSString fileName) {
  ZFSInt noConnectionMaps = 0;
  ZFSInt noSingularityMaps = 0;
  int fid=-1;
  ZFSInt noPeriodicDisplacementInfo = nDim*nDim;
  ZFSFloatScratchSpace globalPeriodicDisplacements(m_noInputBlocks*noPeriodicDisplacementInfo, __CALLING_FUNCTION__, "allPeriodicDisplacements");

  //first read the number of connection information
  //and broadcast it to all partitions
  if(globalDomainId()==0){
    fid=io_openfile("hdf5", fileName.c_str(), "collective", MPI_COMM_SELF);
    io_read_iattribute1(fid, "/Connectivity", "noConnections", &noConnectionMaps);
    io_read_iattribute1(fid, "/Connectivity", "noSingularConnections", &noSingularityMaps);
    io_closefile(fid);
    ZFSInt dummy[2]={noConnectionMaps,noSingularityMaps};
    MPI_Bcast(dummy, 2, MPI_INT, 0, m_zfsStrctrdComm);
  } else {
    ZFSInt dummy[2]={0,0};
    MPI_Bcast(dummy, 2, MPI_INT, 0, m_zfsStrctrdComm);
    noConnectionMaps=dummy[0];
    noSingularityMaps=dummy[1];
  }

  //if connections were looked for previously (hasConnectionInfo is set)
  //but there are no connections there is no need to read them
  if(noConnectionMaps==0 && noSingularityMaps==0) {
    if(globalDomainId()==0) {
      cout << "Connections were previously searched for, none found!" << endl
           << "---> hasConnectionInfo = 1" << endl
           << "---> noConnectionMaps = 0" << endl
           << "---> noSingularityMaps = 0" << endl;
    }
    return;
  }

  if(globalDomainId()==0){
    fid=io_openfile("hdf5", fileName.c_str(), "collective", MPI_COMM_SELF);
    ZFSInt noWindowInformation = 7*nDim+16;
    if(noConnectionMaps>0){
      ZFSIntScratchSpace windowInformation(noWindowInformation*noConnectionMaps, __CALLING_FUNCTION__, "windowInformation");
      for(ZFSId mapId = 0; mapId<noConnectionMaps; mapId++) {
        stringstream windowName;
        windowName << "map" << mapId << endl;
        ZFSString windowNameStr = windowName.str();
        io_read_idataset1(fid, "/Connectivity", windowNameStr.c_str(), 1,  &noWindowInformation, &(windowInformation[mapId*noWindowInformation]));
      }
      MPI_Bcast(windowInformation.getPointer(), noWindowInformation*noConnectionMaps, MPI_INT, 0, m_zfsStrctrdComm);

      for(ZFSId mapId = 0; mapId<noConnectionMaps; mapId++) {
        ZFSStrctrdWindowMap* temp = new ZFSStrctrdWindowMap(nDim);
        readMapFromArray(temp,&windowInformation[mapId*noWindowInformation]);
        globalStrctrdBndryCndMaps.push_back(temp);
      }
    }

    if(noSingularityMaps>0){
      ZFSIntScratchSpace singularityInformation(noWindowInformation*noSingularityMaps, __CALLING_FUNCTION__, "singularityInformation");
      for(ZFSId mapId = 0; mapId<noSingularityMaps; mapId++) {
        stringstream windowName;
        windowName << "singularMap" << mapId << endl;
        ZFSString windowNameStr = windowName.str();
        io_read_idataset1(fid, "/Connectivity", windowNameStr.c_str(), 1,  &noWindowInformation, &(singularityInformation[mapId*noWindowInformation]));
      }
      MPI_Bcast(singularityInformation.getPointer(), noWindowInformation*noSingularityMaps, MPI_INT, 0, m_zfsStrctrdComm);
      for(ZFSId mapId = 0; mapId<noSingularityMaps; mapId++) {
        ZFSStrctrdWindowMap* temp = new ZFSStrctrdWindowMap(nDim);
        readMapFromArray(temp,&(singularityInformation[mapId*noWindowInformation]));
        singularwindow.push_back(temp);
      }
    }

    for(ZFSId blockId = 0; blockId<m_noInputBlocks; blockId++) {
      stringstream path;
      path << "periodicDisplacements" << blockId;
      ZFSString pathStr = path.str();
      if(io_checkObj(fid, "/Connectivity", pathStr.c_str())) {
        //check if path for block individual displacements exists, if not read
        //same dataset for every block
        io_read_ddataset1(fid, "/Connectivity", pathStr.c_str(), 1,  &noPeriodicDisplacementInfo, &globalPeriodicDisplacements[blockId*noPeriodicDisplacementInfo]);
      } else {
        io_read_ddataset1(fid, "/Connectivity", "periodicDisplacements", 1,  &noPeriodicDisplacementInfo, &globalPeriodicDisplacements[blockId*noPeriodicDisplacementInfo]);
      }
    }

    MPI_Bcast(&globalPeriodicDisplacements(0), m_noInputBlocks*noPeriodicDisplacementInfo, MPI_DOUBLE, 0, m_zfsStrctrdComm);
    io_closefile(fid);
  }else{
    ZFSInt noWindowInformation = 7*nDim+16;
    if(noConnectionMaps>0){
      ZFSIntScratchSpace windowInformation(noWindowInformation*noConnectionMaps, __CALLING_FUNCTION__, "windowInformation");
      MPI_Bcast(windowInformation.getPointer(), noWindowInformation*noConnectionMaps, MPI_INT, 0, m_zfsStrctrdComm);
      for(ZFSId mapId = 0; mapId<noConnectionMaps; mapId++) {
        ZFSStrctrdWindowMap* temp = new ZFSStrctrdWindowMap(nDim);
        readMapFromArray(temp,&windowInformation[mapId*noWindowInformation]);
        globalStrctrdBndryCndMaps.push_back(temp);
      }
    }
    if(noSingularityMaps>0){
      ZFSIntScratchSpace singularityInformation(noWindowInformation*noSingularityMaps, __CALLING_FUNCTION__, "singularityInformation");
      MPI_Bcast(singularityInformation.getPointer(), noWindowInformation*noSingularityMaps, MPI_INT, 0, m_zfsStrctrdComm);
      for(ZFSId mapId = 0; mapId<noSingularityMaps; mapId++) {
        ZFSStrctrdWindowMap* temp = new ZFSStrctrdWindowMap(nDim);
        readMapFromArray(temp,&(singularityInformation[mapId*noWindowInformation]));
        singularwindow.push_back(temp);
      }
    }

    MPI_Bcast(&globalPeriodicDisplacements(0), m_noInputBlocks*noPeriodicDisplacementInfo, MPI_DOUBLE, 0, m_zfsStrctrdComm);
  }

  for(ZFSId periodicId=0; periodicId<noPeriodicDisplacementInfo; periodicId++) {
    periodicDisplacements[periodicId] = globalPeriodicDisplacements(m_inputBlockId*noPeriodicDisplacementInfo + periodicId);
  }
}


template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::setSpongeInformation(ZFSInt noSpongeInfo, ZFSFloat* beta, ZFSFloat* sigma, ZFSFloat* thickness, ZFSInt* bcInfo, ZFSInt informationType){
  //we could also use the window information for sponges
  //set through properties.
  ZFSStrctrdWindowMap* temp = new ZFSStrctrdWindowMap(nDim);

  //add the sponge Layter type to it too!!.
  if(informationType){
    for(ZFSUint bcId=0; bcId<globalStrctrdBndryCndMaps.size(); ++bcId ){
      for(ZFSId i=0; i<noSpongeInfo; ++i){
        if(globalStrctrdBndryCndMaps[bcId]->BC==bcInfo[i]){
          mapCpy(globalStrctrdBndryCndMaps[bcId], temp);
          temp->spongeThickness=thickness[i];
          temp->beta=beta[i];
          temp->sigma=sigma[i];
          m_spongeInfoMap.push_back(temp);
          temp=new ZFSStrctrdWindowMap(nDim);
        }
      }
    }
  } else {
    for(ZFSUint bcId=0; bcId<globalStrctrdBndryCndMaps.size(); ++bcId ){
      for(ZFSId i=0; i<noSpongeInfo; ++i){
        if(globalStrctrdBndryCndMaps[bcId]->Id2==bcInfo[i]){
          mapCpy(globalStrctrdBndryCndMaps[bcId], temp);
          temp->spongeThickness=thickness[i];
          temp->beta=beta[i];
          temp->sigma=sigma[i];
          m_spongeInfoMap.push_back(temp);
          temp=new ZFSStrctrdWindowMap(nDim);
        }
      }
    }
  }

  delete temp;
}


template <ZFSInt nDim>    //junoh
void ZFSStrctrdBlckWindowInfo<nDim>::setZonalBCInformation(){
  ZFSStrctrdWindowMap* temp = new ZFSStrctrdWindowMap(nDim);
  for(ZFSId i = 0; i < (ZFSInt)globalStrctrdBndryCndMaps.size(); i++) {
    if(globalStrctrdBndryCndMaps[i]->BC == 2221 || 
       globalStrctrdBndryCndMaps[i]->BC == 2222) {
      mapCpy(globalStrctrdBndryCndMaps[i], temp);
      m_zonalBCMaps.push_back(temp);
      temp=new ZFSStrctrdWindowMap(nDim);
    }
  }
  delete temp;
}

template <ZFSInt nDim>  //junoh
ZFSBool ZFSStrctrdBlckWindowInfo<nDim>::checkZonalBCMaps(ZFSStrctrdWindowMap* map1, ZFSStrctrdWindowMap* map2) {
  ZFSBool isValid = false;
  if(map1->BC==map2->BC) {
    if(map1->Id1==map2->Id1){
      isValid= true;
    }
    // ZFSStrctrdWindowMap* localMapDummy= new ZFSStrctrdWindowMap(nDim);
    // mapCombine11(map1, map2, localMapDummy);
    // mapPrint(localMapDummy);
    // mapPrint(map1);
    // mapPrint(map2);
    // isValid = mapCheck(localMapDummy);
    // delete localMapDummy;
  }
  return isValid;
}



//new function to create list of all wall bcs needed for SA RANS
template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::setWallInformation(){
  ZFSStrctrdWindowMap* temp = new ZFSStrctrdWindowMap(nDim);
  for(ZFSUint bcId=0; bcId<globalStrctrdBndryCndMaps.size(); ++bcId ){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) globalStrctrdBndryCndMaps[bcId]->BC)/1000.0);
    if(firstDigit==1 || globalStrctrdBndryCndMaps[bcId]->BC==2601){
      mapCpy(globalStrctrdBndryCndMaps[bcId], temp);
      m_wallDistInfoMap.push_back(temp);
      temp=new ZFSStrctrdWindowMap(nDim);
    }
  }

  delete temp;
}

template <ZFSInt nDim> //junoh
void ZFSStrctrdBlckWindowInfo<nDim>::createWindowMapping(MPI_Comm* channelIn, MPI_Comm* channelOut,
                                                         MPI_Comm* channelWorld, ZFSInt* channelRoots,
                                                         MPI_Comm* commStg, ZFSInt* commStgRoot, ZFSInt* commStgRootGlobal, MPI_Comm* commBC2600, ZFSInt* commBC2600Root,
                                                         ZFSInt* commBC2600RootGlobal, //junoh
                                                         MPI_Comm* rescalingCommGrComm,
                                                         ZFSInt* rescalingCommGrRoot,
                                                         ZFSInt* rescalingCommGrRootGlobal,
                                                         MPI_Comm* commPerRotOne, MPI_Comm* commPerRotTwo,
                                                         MPI_Comm* commPerRotWorld, ZFSInt* rotationRoots,
                                                         ZFSInt& perRotGroup, SingularInformation *singularity, ZFSInt* hasSingularity)
{
  //create a map of the partition and check whether it fits any of the globals
  //boundary Condition Maps
  //thus this functions relates the partition to the actual window and there-
  //fore also to the boundary condition.

  ZFSInt start1[3]={0,0,0};
  ZFSInt end1[3]={0,0,0};
  ZFSInt step1[3]={1,1,1};
  ZFSInt order[3]={0,1,2};
  ZFSInt start2[3]={0,0,0};
  ZFSInt end2[3]={0,0,0};

  ZFSInt Id = partition->outputBoxInfo[domainId()]->cpu;
  ZFSInt inputBoxId = partition->outputBoxInfo[Id]->inputBoxID;

  /////////////////////////////////////////////////////////////
  //////////////// CREATE OWN/PARTITION MAPS //////////////////
  /////////////////////////////////////////////////////////////
  //this is the map of the active cells (no ghostcells) in my own partition
  //shifted by the no of ghost-cells
  for(ZFSInt i=0; i<nDim; i++) {
    start1[i]=partition->outputBoxInfo[Id]->offset[nDim-1-i]+m_noGhostLayers;//shifted by the number of ghost layers
    start2[i]=m_noGhostLayers;//shifted by the number of ghost layers
    end1[i]=start1[i]+partition->outputBoxInfo[Id]->DirLast[nDim-1-i]-1;
    end2[i]=start2[i]+partition->outputBoxInfo[Id]->DirLast[nDim - 1 -i]-1;
  }

  m_myMapWithoutGC=new ZFSStrctrdWindowMap(nDim);
  mapCreate(inputBoxId, start1, end1, step1, Id , start2, end2, step1, order,-1, m_myMapWithoutGC);

  //this is the map of all the cells  (including ghost-cells) in my own partition
  for(ZFSInt i=0; i<nDim; i++) {
    if(partition->outputBoxInfo[Id]->offset[nDim-1-i]!=0) {
      start1[i]=partition->outputBoxInfo[Id]->offset[nDim-1-i];
      end1[i]=start1[i]+partition->outputBoxInfo[Id]->DirLast[nDim-1-i]-1+2*m_noGhostLayers;
    } else {
      start1[i]=0;
      end1[i]=start1[i]+partition->outputBoxInfo[Id]->DirLast[nDim-1-i]-1+2*m_noGhostLayers;
    }
    start2[i]=0;
    end2[i]=partition->outputBoxInfo[Id]->DirLast[nDim - 1 -i]-1+2*m_noGhostLayers;
  }
  m_myMapWithGC=new ZFSStrctrdWindowMap(nDim);
  mapCreate(inputBoxId, start1, end1, step1, Id , start2, end2, step1, order,-1, m_myMapWithGC);


  //maps of all partitions WITHOUT the ghostcells
  ZFSStrctrdWindowMap* localMapDummy = NULL;
  for(ZFSId j=0; j<noDomains(); j++) {
    localMapDummy= new ZFSStrctrdWindowMap(nDim);

    ZFSInt Id2 = partition->outputBoxInfo[j]->cpu;
    inputBoxId = partition->outputBoxInfo[Id2]->inputBoxID;
    for(ZFSInt i=0; i<nDim; i++) {
      start1[i]=partition->outputBoxInfo[Id2]->offset[nDim-1-i]+m_noGhostLayers;
      start2[i]=0;
      end1[i]=start1[i]+partition->outputBoxInfo[Id2]->DirLast[nDim-1-i]-1;
      end2[i]=partition->outputBoxInfo[Id2]->DirLast[nDim - 1 -i]-1;
    }
    mapCreate(inputBoxId, start1, end1, step1, Id2, start2, end2, step1, order, 6000, localMapDummy);
    m_partitionMapsWithoutGC.push_back(localMapDummy);
  }

  //maps of all partitions WITH the ghostcells
  for(ZFSId j=0; j<noDomains(); j++) {
    localMapDummy= new ZFSStrctrdWindowMap(nDim);

    ZFSInt Id2 = partition->outputBoxInfo[j]->cpu;
    inputBoxId = partition->outputBoxInfo[Id2]->inputBoxID;

    for(ZFSInt i=0; i<nDim; i++) {
      if(partition->outputBoxInfo[Id2]->offset[nDim-1-i]!=0) {
        start1[i]=partition->outputBoxInfo[Id2]->offset[nDim-1-i];
        end1[i]=start1[i]+partition->outputBoxInfo[Id2]->DirLast[nDim-1-i]-1+2*m_noGhostLayers;
      } else {
        start1[i]=0;
        end1[i]=start1[i]+partition->outputBoxInfo[Id2]->DirLast[nDim-1-i]-1+2*m_noGhostLayers;
      }

      start2[i]=0;
      end2[i]=partition->outputBoxInfo[Id2]->DirLast[nDim - 1 -i]-1+4;
    }
    mapCreate(inputBoxId, start1, end1, step1, Id2, start2, end2, step1, order, 6000, localMapDummy);
    m_partitionMapsWithGC.push_back(localMapDummy);
  }

  //////////////////////////////////////////////////////////////////////////////////
  ////////////////// USE MAPS TO CHECK FOR OVERLAPPING PARTS ///////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  //check overlapping of own partition with other partition and put into SND maps
  localMapDummy= new ZFSStrctrdWindowMap(nDim);
  for(ZFSId i=0; i<noDomains(); i++) {
    //skip if comparison to itself
    if(domainId()==i) {
      continue;
    }
    mapCombine11(m_myMapWithoutGC, m_partitionMapsWithGC[i], localMapDummy );
    ZFSBool test=false;
    test=mapCheck(localMapDummy);
    if(test==true) {
      Id = partition->outputBoxInfo[domainId()]->cpu;
      localMapDummy->BC=6000;
      sndMap.push_back(localMapDummy);
      localMapDummy=new ZFSStrctrdWindowMap(nDim);
    }
  }

  //check overlapping of own partition with other partition and put into RCV maps
  for(ZFSId i=0; i<noDomains(); i++) {
    //skip if comparison to itself
    if(domainId()==i) {
      continue;
    }
    mapCombine11(m_myMapWithGC, m_partitionMapsWithoutGC[i], localMapDummy );
    ZFSBool test=false;
    test=mapCheck(localMapDummy);
    if(test==true) {
      Id = partition->outputBoxInfo[domainId()]->cpu;
      localMapDummy->BC=6000;
      rcvMap.push_back(localMapDummy);
      localMapDummy=new ZFSStrctrdWindowMap(nDim);
    }
  }

  //shift globalBC maps by ghost-cells
  for(ZFSInt i=0; i<(ZFSInt)globalStrctrdBndryCndMaps.size(); i++) {
    for(ZFSId dim=0; dim <nDim; dim ++) {
      globalStrctrdBndryCndMaps[i]->start1[dim]+=m_noGhostLayers;
      globalStrctrdBndryCndMaps[i]->end1[dim]+=m_noGhostLayers;
      if(globalStrctrdBndryCndMaps[i]->BC==6000 || (globalStrctrdBndryCndMaps[i]->BC>=4000&&globalStrctrdBndryCndMaps[i]->BC<5000) ||  globalStrctrdBndryCndMaps[i]->BC==6001 ) {
        //for these special cases we also need to change the opposite side
        globalStrctrdBndryCndMaps[i]->start2[dim]+=m_noGhostLayers;
        globalStrctrdBndryCndMaps[i]->end2[dim]+=m_noGhostLayers;
      }
    }
  }

  //go over all the global boundary maps and check if there is any overlapping part with the local partition and store it in localStrctrdBndryCndMaps
  //==>find all the local boundaries of the partition (including the periodic part)

  //check for local boundary maps
  for(ZFSInt i=0; i<(ZFSInt)globalStrctrdBndryCndMaps.size(); i++) {
    mapCombine11(m_myMapWithoutGC, globalStrctrdBndryCndMaps[i], localMapDummy );
    ZFSBool test=false;
    test=mapCheck(localMapDummy);

    if(test==true) {
      localMapDummy->BC=globalStrctrdBndryCndMaps[i]->BC;
      localMapDummy->Nstar=globalStrctrdBndryCndMaps[i]->Nstar;
      localMapDummy->SingularId=globalStrctrdBndryCndMaps[i]->SingularId;
      localMapDummy->dc1=globalStrctrdBndryCndMaps[i]->dc1;
      localMapDummy->dc2=globalStrctrdBndryCndMaps[i]->dc2;
      localMapDummy->Id1=globalStrctrdBndryCndMaps[i]->Id1; //junoh
      localStrctrdBndryCndMaps.push_back(localMapDummy);
      localMapDummy=new ZFSStrctrdWindowMap(nDim);
    }
  }

  for(ZFSInt j=0; j<(ZFSInt)m_partitionMapsWithoutGC.size(); j++) {
    for(ZFSId i=0; i<nDim; i++) {
      m_partitionMapsWithoutGC[j]->start2[i]+=m_noGhostLayers;
      m_partitionMapsWithoutGC[j]->end2[i]+=m_noGhostLayers;
    }
  }

  //check overlapping with bc 6000 (multiblock)
  vector<ZFSStrctrdWindowMap*> SndMapBC6000;

  for(ZFSInt i=0; i<(ZFSInt)globalStrctrdBndryCndMaps.size(); i++){
    if(globalStrctrdBndryCndMaps[i]->BC==6000) {
      for(ZFSInt j=0; j<(ZFSInt)m_partitionMapsWithoutGC.size(); j++) {
        //skip own domain, but only if 6000 does not connect the same domain
        if(globalStrctrdBndryCndMaps[i]->Id1!=globalStrctrdBndryCndMaps[i]->Id2) {
          if(partition->outputBoxInfo[domainId()]->cpu==m_partitionMapsWithoutGC[j]->Id2) {
            continue;
          }
        }

        ZFSInt cpuid = partition->outputBoxInfo[domainId()]->cpu;
        ZFSInt ownbigblock = partition->outputBoxInfo[cpuid]->inputBoxID;
        if(globalStrctrdBndryCndMaps[i]->Id2==ownbigblock) {
          mapCombine11( globalStrctrdBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy );
          ZFSBool test=false;
          test=mapCheck(localMapDummy);

          if(test==true) {
            localMapDummy->BC=globalStrctrdBndryCndMaps[i]->BC;
            localMapDummy->Nstar=globalStrctrdBndryCndMaps[i]->Nstar;
            localMapDummy->SingularId=globalStrctrdBndryCndMaps[i]->SingularId;
            localMapDummy->dc1=globalStrctrdBndryCndMaps[i]->dc1;
            localMapDummy->dc2=globalStrctrdBndryCndMaps[i]->dc2;

            SndMapBC6000.push_back(localMapDummy);
            localMapDummy=new ZFSStrctrdWindowMap(nDim);
          }
        }
      }
    }
  }


  //now check overlapping with bc 4xxx (periodic bc)
  for(ZFSInt i=0; i<(ZFSInt)globalStrctrdBndryCndMaps.size(); i++) {
    if(globalStrctrdBndryCndMaps[i]->BC>=4000&&globalStrctrdBndryCndMaps[i]->BC<5000) {
      for(ZFSInt j=0; j<(ZFSInt)m_partitionMapsWithoutGC.size(); j++) {
        //self communication is allowed
        mapCombine11( globalStrctrdBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy );
        ZFSBool test=false;
        test=mapCheck(localMapDummy);

        if(test==true) {
          localMapDummy->BC=globalStrctrdBndryCndMaps[i]->BC;
          localMapDummy->Nstar=globalStrctrdBndryCndMaps[i]->Nstar;
          localMapDummy->SingularId=globalStrctrdBndryCndMaps[i]->SingularId;
          localMapDummy->dc1=globalStrctrdBndryCndMaps[i]->dc1;
          localMapDummy->dc2=globalStrctrdBndryCndMaps[i]->dc2;
          SndMapBC6000.push_back(localMapDummy);
          localMapDummy=new ZFSStrctrdWindowMap(nDim);
        }
      }
    }
  }

  //now 6000 and 4xxx communication bc are both in SndMpaBC6000

  ///////////////////////////////////////////////
  /////// SINGULARITY COMMUNICATION  ////////////
  ///////////////////////////////////////////////
  //special treatment of singularities
  ZFSInt singularPoint=(ZFSInt)singularwindow.size();
  for(ZFSInt i=0; i<singularPoint; ++i) {
    for(ZFSId dim=0; dim <nDim; ++dim) {
      singularwindow[i]->start1[dim]+=m_noGhostLayers;
      singularwindow[i]->end1[dim]+=m_noGhostLayers;
      //need also to change the opposite side
      singularwindow[i]->start2[dim]+=m_noGhostLayers;
      singularwindow[i]->end2[dim]+=m_noGhostLayers;
    }
  }

  for(ZFSInt i=0; i<singularPoint; ++i) {
    ZFSInt singularcount=0,limit=singularwindow[i]->Nstar-3;
    for(ZFSInt j=0; j<(ZFSInt)globalStrctrdBndryCndMaps.size(); j++) {
      if((globalStrctrdBndryCndMaps[j]->Nstar==5||globalStrctrdBndryCndMaps[j]->Nstar==6)&&globalStrctrdBndryCndMaps[j]->BC==6000&&globalStrctrdBndryCndMaps[j]->Nstar==singularwindow[i]->Nstar) {

        ZFSInt labell=0;
        labell = mapCompare11(singularwindow[i], globalStrctrdBndryCndMaps[j]);
        if(labell) {
          singularwindow[i]->SingularBlockId[singularcount]=globalStrctrdBndryCndMaps[j]->Id2;
          singularcount++;
        }
      }
      if(singularcount>limit) cout<<"Important!! "<<globalStrctrdBndryCndMaps[j]->Nstar<<" star check error!!!!!!!!!!!!!!!!!! check it !!!!!"<<endl;
    }
  }

  // if(domainId()==0)
  //   {
  //     for(ZFSInt i=0; i<singularPoint; ++i) {
  //     if(singularwindow[i]->BC==6000)
  //       cout<< "singularmaps  Id1: " << singularwindow[i]->Id1 << " Id2: " << singularwindow[i]->Id2 << " BC: " << singularwindow[i]->BC<<"  "<<  singularwindow[i]->start1[0] <<"-" <<  singularwindow[i]->end1[0]<< " "<<  singularwindow[i]->start1[1] <<"-" <<  singularwindow[i]->end1[1]<< " "<<  singularwindow[i]->start1[2] <<"-" <<  singularwindow[i]->end1[2]<< " ID: " <<singularwindow[i]->SingularBlockId[0]<<"  "<<singularwindow[i]->SingularBlockId[1]<< endl;

  //     }
  //   }

  for(ZFSInt i=0; i<singularPoint; ++i) {
    mapCombine11(m_myMapWithoutGC,singularwindow[i], localMapDummy );
    ZFSBool test=false;
    test=mapCheck(localMapDummy);
    localMapDummy->Nstar=singularwindow[i]->Nstar;
    localMapDummy->SingularId=singularwindow[i]->SingularId;
    localMapDummy->BC=singularwindow[i]->BC;
    localMapDummy->SingularBlockId[0]=singularwindow[i]->SingularBlockId[0];
    localMapDummy->SingularBlockId[1]=singularwindow[i]->SingularBlockId[1];
    localMapDummy->SingularBlockId[2]=singularwindow[i]->SingularBlockId[2];
    localMapDummy->SingularBlockId[3]=singularwindow[i]->SingularBlockId[3];

    if(test==true) {
      localSingularMap.push_back(localMapDummy);
      localMapDummy=new ZFSStrctrdWindowMap(nDim);
    }
  }

  *hasSingularity=(ZFSInt)localSingularMap.size();
  singularPoint=(ZFSInt)localSingularMap.size();

  for(ZFSInt i=0; i<(ZFSInt)localSingularMap.size(); ++i) {
    memcpy(singularity[i].start,localSingularMap[i]->start1,nDim*sizeof(ZFSInt));
    memcpy(singularity[i].end,  localSingularMap[i]->end1,  nDim*sizeof(ZFSInt));
    singularity[i].Nstar=localSingularMap[i]->Nstar;
    singularity[i].count=localSingularMap[i]->SingularId;
    singularity[i].BC=localSingularMap[i]->BC;

    singularity[i].SingularBlockId[0]=localSingularMap[i]->SingularBlockId[0];
    singularity[i].SingularBlockId[1]=localSingularMap[i]->SingularBlockId[1];
    singularity[i].SingularBlockId[2]=localSingularMap[i]->SingularBlockId[2];
    singularity[i].SingularBlockId[3]=localSingularMap[i]->SingularBlockId[3];

    //!!!!important!!!!
    //for local singular map
    //BC6000: line only (no point and face) BC4000: point only.

    if(localSingularMap[i]->BC>=4400&&localSingularMap[i]->BC<4410) {
      localSingularMap[i]->face=localSingularMap[i]->BC-4401;
    } else {
      for(ZFSInt k=0; k<3;++k) {
        if(localSingularMap[i]->start1[k]!=localSingularMap[i]->end1[k]) {
          localSingularMap[i]->face=k*2;
        }
      }
    }

    ZFSInt displace[3],dimcount=0;
    ZFSInt dimN,dimT1,dimT2;
    for(ZFSInt k=0; k<3;++k) {
      if(singularity[i].start[k]==singularity[i].end[k]) {
        if(singularity[i].start[k]==2) {
          displace[k]=-1;
        } else {
          displace[k]=1;
        }
        dimcount++;
      } else {
        displace[k]=0;
        dimN=k;
      }
    }

    if(dimcount==3) {
      if(localSingularMap[i]->face!=-1) {
        dimN=localSingularMap[i]->face/2;
        displace[dimN]=0;
        if(singularity[i].BC==6000) {
          singularity[i].BC=6003;
        }
      } else {
        cout<<"ERROR!!! point singular communication cannot decide the direction!! check the singularity!!!"<<endl;
        exit(0);
      }
    }

    if(dimN==0)     {dimT1=1;dimT2=2;}
    else if(dimN==1){dimT1=0;dimT2=2;}
    else            {dimT1=0;dimT2=1;}

    switch(singularity[i].Nstar)
    {
    case 3:
    {
      // bugs possibly exist
      singularity[i].displacement[0][dimN]=0;
      singularity[i].displacement[0][dimT1]=displace[dimT1];
      singularity[i].displacement[0][dimT2]=0;

      singularity[i].displacement[1][dimN]=0;
      singularity[i].displacement[1][dimT1]=0;
      singularity[i].displacement[1][dimT2]=displace[dimT2];

      break;
    }
    case 5:
    {
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // bugs  exist need to be fixed
      singularity[i].displacement[0][dimN]=0;
      singularity[i].displacement[0][dimT1]=displace[dimT1];
      singularity[i].displacement[0][dimT2]=0;

      singularity[i].displacement[1][dimN]=0;
      singularity[i].displacement[1][dimT1]=0;
      singularity[i].displacement[1][dimT2]=displace[dimT2];

      singularity[i].displacement[2][dimN]=0;
      singularity[i].displacement[2][dimT1]=displace[dimT1];
      singularity[i].displacement[2][dimT2]=displace[dimT2];

      singularity[i].displacement[3][dimN]=0;
      singularity[i].displacement[3][dimT1]=2*displace[dimT1];
      singularity[i].displacement[3][dimT2]=2*displace[dimT2];


      //communication map change start from 2 to 3
      singularity[i].count=2;

      break;
    }
    case 6:
    {
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // bugs  exist need to be fixed
      singularity[i].displacement[0][dimN]=0;
      singularity[i].displacement[0][dimT1]=displace[dimT1];
      singularity[i].displacement[0][dimT2]=0;

      singularity[i].displacement[1][dimN]=0;
      singularity[i].displacement[1][dimT1]=0;
      singularity[i].displacement[1][dimT2]=displace[dimT2];

      singularity[i].displacement[2][dimN]=0;
      singularity[i].displacement[2][dimT1]=displace[dimT1];
      singularity[i].displacement[2][dimT2]=displace[dimT2];

      singularity[i].displacement[3][dimN]=0;
      singularity[i].displacement[3][dimT1]=2*displace[dimT1];
      singularity[i].displacement[3][dimT2]=2*displace[dimT2];

      singularity[i].displacement[4][dimN]=0;
      singularity[i].displacement[4][dimT1]=displace[dimT1];
      singularity[i].displacement[4][dimT2]=2*displace[dimT2];


      //communication map change start from 2 to 3
      singularity[i].count=2;

      break;
    }
    default:
    {
      break;
    }
    }
  }

  //expand the singularmap
  for(ZFSInt i=0; i<(ZFSInt)localSingularMap.size(); ++i)
  {
    if(localSingularMap[i]->BC==6000)
    {
      for (ZFSInt k=0;k<3;++k)
      {
        if(localSingularMap[i]->start1[k]!=localSingularMap[i]->end1[k])
        {
          // if ( m_myMapWithoutGC->start2[k]==localSingularMap[i]->start1[k])
          {
            localSingularMap[i]->start1[k]-=2;
            localSingularMap[i]->start2[k]-=2;
          }
          // if ( m_myMapWithoutGC->end2[k]==localSingularMap[i]->end1[k])
          {
            localSingularMap[i]->end1[k]+=2;
            localSingularMap[i]->end2[k]+=2;
          }
        }//if !=
      }//for k
    }//if 6000
  }

  //first put multiblock communication maps into addCommBC,
  //afterwards also add the periodic communcation maps

  vector<ZFSStrctrdWindowMap*> periodicBC;
  vector<ZFSStrctrdWindowMap*> addCommBC;
  MPI_Barrier(m_zfsStrctrdComm);
  for(ZFSId i=0; i<(ZFSInt)localStrctrdBndryCndMaps.size(); i++)
  {
    switch(localStrctrdBndryCndMaps[i]->BC)
    {
    case 4011: //periodic rotational 1
    case 4012: //periodic rotational 2
    case 4401: //periodic surface 1
    case 4402: //periodic surface 2
    case 4403:
    case 4404:
    case 4405:
    case 4406:
    {
      break;
    }
    case 6000: //additional block communication (interblock communication)
    {
      //save the additional Communication
      addCommBC.push_back(localStrctrdBndryCndMaps[i]);
      break;
    }
    case 2401:
    case 2402: //this is for the channel flow
    {
      //we need a copy of the map to find the channel In/Outflow parti-
      //cipating processors.
      physicalBCMap.push_back(localStrctrdBndryCndMaps[i]);
      ZFSStrctrdWindowMap* copy = new ZFSStrctrdWindowMap(nDim);
      mapCpy(localStrctrdBndryCndMaps[i], copy);
      channelSurfaceIndices.push_back(copy); //needed for identifiaction of aereas etc.

      break;
    }
    default:
    {
      physicalBCMap.push_back(localStrctrdBndryCndMaps[i]);
      break;
    }
    }

  }

  for(ZFSId i=0; i<(ZFSInt)localStrctrdBndryCndMaps.size(); i++)
  {
    switch(localStrctrdBndryCndMaps[i]->BC)
    {
    case 4011: //periodic rotational 1
    case 4012: //periodic rotational 2
    case 4401: //periodic surface 1
    case 4402: //periodic surface 2
    case 4403:
    case 4404:
    case 4405:
    case 4406:
    {
      //save the periodic boundary condition to addcommBC, same position
      addCommBC.push_back(localStrctrdBndryCndMaps[i]);
      break;
    }
    default:
    {
      //do nothing
      break;
    }
    }
  }


  ///////////////////////////////////////////
  //////// CHANNEL BC ///////////////////////
  ///////////////////////////////////////////

  //creating the communicators
  vector<ZFSInt> channelface;
  vector<ZFSInt> channelpartitions;
  int anzahl =1;
  ZFSInt counterIn=0;
  ZFSInt counterOut=0;
  for(ZFSId i=0; i<(ZFSInt)globalStrctrdBndryCndMaps.size(); i++) {
    if(globalStrctrdBndryCndMaps[i]->BC==2401||globalStrctrdBndryCndMaps[i]->BC==2402) {
      //found a face containing the channel in/out flow
      //now check the partition it belongs to!!
      //cout << "found " << anzahl << endl;
      anzahl++;
      for(ZFSId j=0; j<(ZFSInt)m_partitionMapsWithoutGC.size(); j++) {
        mapCombine11( globalStrctrdBndryCndMaps[i] , m_partitionMapsWithoutGC[j], localMapDummy );
        ZFSBool test=false;
        test=mapCheck(localMapDummy);
        if(test==true) {
          channelpartitions.push_back(m_partitionMapsWithoutGC[j]->Id2);
          //partition contains part of the boundary
          //check the face

          //now check the face
          if(globalStrctrdBndryCndMaps[i]->BC == 2401 ) {
            channelface.push_back(globalStrctrdBndryCndMaps[i]->BC);
            counterIn++;
          } else {
            channelface.push_back(globalStrctrdBndryCndMaps[i]->BC);
            counterOut++;
          }
        }
      }
    }
  }

  //now assign the ranks
  ZFSInt* channelInflow=NULL;
  ZFSInt* channelOutflow=NULL;
  if(channelpartitions.size()!=0) {
    ZFSInt counterI=0, counterO=0;
    channelInflow = new ZFSInt[counterIn];
    channelOutflow = new ZFSInt[counterOut];
    for(ZFSInt i=0; i<(ZFSInt) channelpartitions.size() ; i++) {
      for(ZFSInt j=0; j<noDomains(); j++) {
        if(channelpartitions[i]==partition->outputBoxInfo[j]->cpu) {
          if(channelface[i]==2401) {
            channelInflow[counterI]=j;
            counterI++;
          }
          if(channelface[i]==2402) {
            channelOutflow[counterO]=j;
            counterO++;
          }
        }
      }
    }


    zfs_log << "ChannelIn: " << counterI << " ChannelOut: " << counterOut << endl;

    //create a communicator for the in- and the outflow
    //averaging process
    MPI_Group groupIn=MPI_GROUP_NULL;
    MPI_Group* newgroupIn= new MPI_Group;
    MPI_Group groupOut=MPI_GROUP_NULL;
    MPI_Group* newgroupOut=new MPI_Group;

    MPI_Comm_group(m_zfsStrctrdComm, &groupIn);
    MPI_Comm_group(m_zfsStrctrdComm, &groupOut);

    MPI_Group_incl(groupIn, (int)counterIn , channelInflow, newgroupIn);
    MPI_Comm_create(m_zfsStrctrdComm, newgroupIn[0], channelIn);

    MPI_Group_incl(groupOut, (int)counterOut , channelOutflow, newgroupOut);
    MPI_Comm_create(m_zfsStrctrdComm, newgroupOut[0], channelOut);

    //finally we need to create a overall channel surface communicator to distribute
    //the averaged values p & T

    //create group of all in- and outflow participants
    //but pay attention to duplicates (one domain shares in- and outflow)
    ZFSInt counterAll = 0;
    ZFSInt* channelAll= new ZFSInt[counterOut+counterIn];
    for(ZFSInt i=0; i<counterIn; i++) {
      channelAll[i]=channelInflow[i];
      counterAll++;
    }
    for(ZFSInt i=0; i<counterOut; i++) {
      ZFSBool check = false;
      for(ZFSInt j=0; j<counterIn; ++j) {
        if(channelAll[j] == channelOutflow[i]) {
          check = true;
          cout << "ATTENTION: Skipping duplicate generation in channelAll" << endl;
          break;
        }
      }
      //if this cpu is already in channelAll we must not insert it again
      if(!check) {
        channelAll[counterIn+i]=channelOutflow[i];
        counterAll++;
      }
    }
    MPI_Group groupAll = MPI_GROUP_NULL;
    MPI_Group* newgroupAll = new MPI_Group;
    MPI_Comm_group(m_zfsStrctrdComm, &groupAll);
    MPI_Group_incl(groupAll, (int)(counterAll), channelAll, newgroupAll);
    MPI_Comm_create(m_zfsStrctrdComm, newgroupAll[0], channelWorld);

    //for global exchange we need a root process
    if(domainId()==channelInflow[0])
    {
      MPI_Comm_rank(channelIn[0], &channelRoots[0]);
      MPI_Comm_rank(channelWorld[0], &channelRoots[2]);
    }
    MPI_Barrier(m_zfsStrctrdComm);
    MPI_Bcast(&channelRoots[0], 1, MPI_INT, channelInflow[0], m_zfsStrctrdComm);
    MPI_Bcast(&channelRoots[2], 1, MPI_INT, channelInflow[0], m_zfsStrctrdComm);

    if(domainId()==channelOutflow[0])
    {
      MPI_Comm_rank(channelOut[0], &channelRoots[1]);
      MPI_Comm_rank(channelWorld[0], &channelRoots[3]);
    }
    MPI_Barrier(m_zfsStrctrdComm);
    MPI_Bcast(&channelRoots[1], 1, MPI_INT, channelOutflow[0], m_zfsStrctrdComm);
    MPI_Bcast(&channelRoots[3], 1, MPI_INT, channelOutflow[0], m_zfsStrctrdComm);

    // if(domainId()==0) {
    //          cout<<"&&&&inlet cpu: ";
    //          for(ZFSInt i=0; i<counterIn; i++) {
    //            cout<<channelInflow[i]<<" ";
    //          }
    //          cout<<endl;
    //          cout<<"&&&&outlet cpu: ";
    //          for(ZFSInt i=0; i<counterOut; i++) {
    //            cout<<channelOutflow[i]<<" ";
    //          }
    //          cout<<endl;
    //          cout<<"&&&&channel cpu: ";
    //          for(ZFSInt i=0; i<counterAll; i++) {
    //            cout<<channelAll[i]<<" ";
    //          }
    //          cout<<endl;
    //          cout<<"roots: "<<channelRoots[0]<<" "<<channelRoots[1]<<" "<<channelRoots[2]<<" "<<channelRoots[3]<<endl;
    // }
  }

  //////////////////////////////////////////
  //////// PERIODIC ROTATION BC ////////////
  //////////////////////////////////////////

  vector<ZFSInt> rotatingface;
  vector<ZFSInt> rotatingpartitions;

  ZFSInt noRotPartitions =0;
  counterIn=0;
  counterOut=0;

  for(ZFSId i=0; i<(ZFSInt)globalStrctrdBndryCndMaps.size(); i++){
    if(globalStrctrdBndryCndMaps[i]->BC==4001 || globalStrctrdBndryCndMaps[i]->BC==4002){ //== containing rotation info
      //found a face containing the per. rotation BC ==> which partition?
      noRotPartitions++;
      for(ZFSId j=0; j<(ZFSInt)m_partitionMapsWithoutGC.size(); j++){
        mapCombine11( globalStrctrdBndryCndMaps[i] , m_partitionMapsWithoutGC[j], localMapDummy );
        ZFSBool test=false;
        test=mapCheck(localMapDummy);
        if(test==true){
          rotatingpartitions.push_back(m_partitionMapsWithoutGC[j]->Id2);
          //partition contains part of the boundary
          if(globalStrctrdBndryCndMaps[i]->BC == 4001){
            rotatingface.push_back(globalStrctrdBndryCndMaps[i]->BC);
            counterIn++;
          } else {
            rotatingface.push_back(globalStrctrdBndryCndMaps[i]->BC);
            counterOut++;
          }
        }
      }
    }
  }

  ZFSInt* rotationOne=NULL;
  ZFSInt* rotationTwo=NULL;
  if(rotatingpartitions.size()!=0){
    ZFSInt counterI=0, counterO=0;
    rotationOne = new ZFSInt[counterIn];
    rotationTwo = new ZFSInt[counterOut];
    for(ZFSInt i=0; i<(ZFSInt) rotatingpartitions.size() ; i++){
      for(ZFSInt j=0; j<noDomains(); j++){
        if(rotatingpartitions[i]==partition->outputBoxInfo[j]->cpu){
          if(rotatingface[i]==4001){
            rotationOne[counterI]=j;
            ++counterI;
          } else {
            rotationTwo[counterO]=j;
            counterO++;
          }
        }
      }
    }
    //check which group the processor is in
    perRotGroup=0;
    for(ZFSInt i=0;i<counterI;++i){
      if(domainId()==rotationOne[i]){ perRotGroup=1; break;}
    }
    for (ZFSInt i=0;i<counterO;++i){
      if(domainId()==rotationTwo[i]){
        //Maybe domain is in both groups, then set to 3
        if(perRotGroup == 1)
          perRotGroup=3;
        else
          perRotGroup=2;

        break;
      }
    }

    //create a communicator for the two sides of the rotation BC
    MPI_Group groupOne=MPI_GROUP_NULL;
    MPI_Group* newgroupOne= new MPI_Group;
    MPI_Group groupTwo=MPI_GROUP_NULL;
    MPI_Group* newgroupTwo=new MPI_Group;

    MPI_Comm_group(m_zfsStrctrdComm, &groupOne);
    MPI_Comm_group(m_zfsStrctrdComm, &groupTwo);

    MPI_Group_incl(groupOne, (int)counterIn , rotationOne, newgroupOne);
    MPI_Comm_create(m_zfsStrctrdComm, *newgroupOne, commPerRotOne);

    MPI_Group_incl(groupTwo, (int)counterOut , rotationTwo, newgroupTwo);
    MPI_Comm_create(m_zfsStrctrdComm, *newgroupTwo, commPerRotTwo);

    //we also need one cross communicator containing all processors involved
    //in the periodic rotation boundary condition

    ZFSInt counterAll = 0;
    ZFSInt* rotationAll= new ZFSInt[counterOut+counterIn];
    for(ZFSInt i=0; i<counterIn; ++i){
      rotationAll[i]=rotationOne[i];
      counterAll++;
    }

    for(ZFSInt i=0; i<counterOut; ++i){
      //check if this cpu is already in rotationAll
      ZFSBool check = false;
      for(ZFSInt j=0; j<counterIn; ++j){
        if(rotationAll[j] == rotationTwo[i]){
          check = true;
          cout << "ATTENTION: Skipping duplicate generation in rotationAll" << endl;
          break;
        }
      }

      //if this cpu is already in rotationAll we must not insert it again
      if(!check) {
        rotationAll[counterIn+i]=rotationTwo[i];
        counterAll++;
      }
    }

    MPI_Group groupAll = MPI_GROUP_NULL;
    MPI_Group* newgroupAll = new MPI_Group;

    MPI_Comm_group(m_zfsStrctrdComm, &groupAll);
    MPI_Group_incl(groupAll, (int)(counterAll), rotationAll, newgroupAll);
    MPI_Comm_create(m_zfsStrctrdComm, *newgroupAll, commPerRotWorld);
    //fix the root processes for the communication
    if(domainId()==rotationOne[0]){
      MPI_Comm_rank(commPerRotOne[0], &rotationRoots[0]);
      MPI_Comm_rank(commPerRotWorld[0], &rotationRoots[2]);
    }
    MPI_Barrier(m_zfsStrctrdComm);
    MPI_Bcast(&rotationRoots[0], 1, MPI_INT, rotationOne[0], m_zfsStrctrdComm);
    MPI_Bcast(&rotationRoots[2], 1, MPI_INT, rotationOne[0], m_zfsStrctrdComm);

    if(domainId()==rotationTwo[0]){
      MPI_Comm_rank(commPerRotTwo[0], &rotationRoots[1]);
      MPI_Comm_rank(commPerRotWorld[0], &rotationRoots[3]);
    }
    MPI_Barrier(m_zfsStrctrdComm);
    MPI_Bcast(&rotationRoots[1], 1, MPI_INT, rotationTwo[0], m_zfsStrctrdComm);
    MPI_Bcast(&rotationRoots[3], 1, MPI_INT, rotationTwo[0], m_zfsStrctrdComm);
  }



  /////////////////////////////////////////////////////////
  ///////////////// STG BC Communication //////////////////
  /////////////////////////////////////////////////////////

  /*
    Create the MPI Group for the STG methods that only includes those
    domains that are located at the STG inflow.
  */

  vector<ZFSInt> stgpartitions;

  for(ZFSId i = 0; i < (ZFSInt)globalStrctrdBndryCndMaps.size(); i++) {
    if(globalStrctrdBndryCndMaps[i]->BC == 7909) {
      // cout<<"domainId in stg comm:"<<domainId()<<endl;
      for(ZFSId j = 0; j < (ZFSInt)m_partitionMapsWithoutGC.size(); j++) {
        mapCombine11(globalStrctrdBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
        ZFSBool test = false;
        test = mapCheck(localMapDummy);
        if(test == true)
          stgpartitions.push_back(m_partitionMapsWithoutGC[j]->Id2);
      }
    }
  }
    // for(ZFSId j=0; j<noDomains(); j++) {
    //   cout<<"partition->outputBoxInfo["<<j<<"]->inputBoxID:"<<partition->outputBoxInfo[j]->inputBoxID<<endl;
    // }
  // cout<<"stgpartitions.size():"<<stgpartitions.size()<<endl;
  //   for(ZFSId j = 0; j < (ZFSInt)m_partitionMapsWithoutGC.size(); j++) {
  //     cout<<"m_partitionMapsWithoutGC["<<j<<"]->Id2:"<<m_partitionMapsWithoutGC[j]->Id2<<" m_partitionMapsWithoutGC["<<j<<"]->Id1"<<m_partitionMapsWithoutGC[j]->Id1<<" m_partitionMapsWithoutGC["<<j<<"]->inputBoxID:"<<m_partitionMapsWithoutGC[j]->inputBoxID<<endl;
  //   }

  MPI_Barrier(m_zfsStrctrdComm);

  ZFSInt* stgranks;
  if(stgpartitions.size() != 0) {
    zfs_log << "STG CPUs: ";
    int counterstg = 0;
    stgranks = new ZFSInt[stgpartitions.size()];
    for(ZFSId i = 0; i < (ZFSInt)stgpartitions.size(); i++) {
      for(ZFSId j = 0; j < noDomains(); j++) {
        if(stgpartitions[i] == partition->outputBoxInfo[j]->cpu) {
          stgranks[counterstg] = j;
          counterstg++;
          zfs_log << partition->outputBoxInfo[j]->cpu << ", ";
	    // cout<<"partition->outputBoxInfo[j]->cpu:"<<partition->outputBoxInfo[j]->cpu<<endl;

        }
      }
    }
    zfs_log << endl;
    zfs_log << "Total number of STG cpus: " << counterstg << endl;

    MPI_Barrier(m_zfsStrctrdComm);

    MPI_Group groupStg, newgroupStg;
    int stgcommsize = (int)(stgpartitions.size());

    MPI_Comm_group(m_zfsStrctrdComm, &groupStg);
    MPI_Group_incl(groupStg, stgcommsize, stgranks, &newgroupStg);
    MPI_Comm_create(m_zfsStrctrdComm, newgroupStg, commStg);

    if(domainId() == stgranks[0]) {
      MPI_Comm_rank(commStg[0], &commStgRoot[0]);
      MPI_Comm_rank(m_zfsStrctrdComm, &commStgRootGlobal[0]);
    }

    MPI_Barrier(m_zfsStrctrdComm);
    MPI_Bcast(commStgRoot, 1, MPI_INT, stgranks[0], m_zfsStrctrdComm);
    MPI_Bcast(commStgRootGlobal, 1, MPI_INT, stgranks[0], m_zfsStrctrdComm);
  }



  //////////////////////////////////////////////////////
  ////////////////BC 2600 Communication////////////////
  ///////////////////////////////////////////////////// //junoh
  

  vector<ZFSInt> bc2600partitions;

  for(ZFSId i = 0; i < (ZFSInt)globalStrctrdBndryCndMaps.size(); i++) {
    if(globalStrctrdBndryCndMaps[i]->BC == 2600) {
      // cout<<"domainId in 2600 comm:"<<domainId()<<endl;
      for(ZFSId j = 0; j < (ZFSInt)m_partitionMapsWithoutGC.size(); j++) {
        mapCombine11(globalStrctrdBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
        ZFSBool test = false;
        test = mapCheck(localMapDummy);
        if(test == true)
          bc2600partitions.push_back(m_partitionMapsWithoutGC[j]->Id2);
      }
    }
  }
  // cout<<"stgpartitions.size():"<<bc2600partitions.size()<<endl;
  //  for(ZFSId i = 0; i < (ZFSInt)globalStrctrdBndryCndMaps.size(); i++) {

  //    cout<<"globalStrctrdBndryCndMaps["<<i<<"]->Id1:"<<globalStrctrdBndryCndMaps[i]->Id1<<"globalStrctrdBndryCndMaps["<<i<<"]->Id2:"<<globalStrctrdBndryCndMaps[i]->Id2<<endl;
  //  }

  MPI_Barrier(m_zfsStrctrdComm);

  ZFSInt* bc2600ranks;
  if(bc2600partitions.size() != 0) {
    zfs_log << "bc2600 CPUs: ";
    int counterBC2600 = 0;
    bc2600ranks = new ZFSInt[bc2600partitions.size()];
    for(ZFSId i = 0; i < (ZFSInt)bc2600partitions.size(); i++) {
      for(ZFSId j = 0; j < noDomains(); j++) {
        if(bc2600partitions[i] == partition->outputBoxInfo[j]->cpu) {
          bc2600ranks[counterBC2600] = j;
          counterBC2600++;
          zfs_log << partition->outputBoxInfo[j]->cpu << ", ";
	    // cout<<"partition->outputBoxInfo[j]->cpu:"<<partition->outputBoxInfo[j]->cpu<<endl;

        }
      }
    }
    zfs_log << endl;
    zfs_log << "Total number of BC2600 cpus: " << counterBC2600 << endl;

    MPI_Barrier(m_zfsStrctrdComm);

    MPI_Group groupBC2600, newgroupBC2600;
    int bc2600commsize = (int)(bc2600partitions.size());

    MPI_Comm_group(m_zfsStrctrdComm, &groupBC2600);
    MPI_Group_incl(groupBC2600, bc2600commsize, bc2600ranks, &newgroupBC2600);
    MPI_Comm_create(m_zfsStrctrdComm, newgroupBC2600, commBC2600);

    if(domainId() == bc2600ranks[0]) {
      MPI_Comm_rank(commBC2600[0], &commBC2600Root[0]);
      MPI_Comm_rank(m_zfsStrctrdComm, &commBC2600RootGlobal[0]);
    }

    MPI_Barrier(m_zfsStrctrdComm);
    MPI_Bcast(commBC2600Root, 1, MPI_INT, bc2600ranks[0], m_zfsStrctrdComm);
    MPI_Bcast(commBC2600RootGlobal, 1, MPI_INT, bc2600ranks[0], m_zfsStrctrdComm);
  }







  /////////////////////////////////////////////////////////
  ///////////////// Rescaling BC Communication ////////////
  /////////////////////////////////////////////////////////


  vector<ZFSInt> rescalingCommGrPartitions;
  for(ZFSId i = 0; i < (ZFSInt)globalStrctrdBndryCndMaps.size(); i++) {
    if(globalStrctrdBndryCndMaps[i]->BC == 2500 || globalStrctrdBndryCndMaps[i]->BC==2501) {
      for(ZFSId j = 0; j < (ZFSInt)m_partitionMapsWithoutGC.size(); j++) {
        mapCombine11(globalStrctrdBndryCndMaps[i], m_partitionMapsWithoutGC[j], localMapDummy);
        ZFSBool test = false;
        test = mapCheck(localMapDummy);
        if(test == true) {
          rescalingCommGrPartitions.push_back(m_partitionMapsWithoutGC[j]->Id2);
        }
      }
    }
  }

  MPI_Barrier(m_zfsStrctrdComm);

  ZFSInt* rescalingCommGrRanks;

  if(rescalingCommGrPartitions.size() != 0) {
    int counterCommGrRescaling = 0;
    rescalingCommGrRanks = new ZFSInt[rescalingCommGrPartitions.size()];
    for(ZFSId i = 0; i < (ZFSInt)rescalingCommGrPartitions.size(); i++) {
      for(ZFSId j = 0; j < noDomains(); j++) {
        if(rescalingCommGrPartitions[i] == partition->outputBoxInfo[j]->cpu) {
          rescalingCommGrRanks[counterCommGrRescaling] = j;
          counterCommGrRescaling++;
        }
      }
    }

    MPI_Barrier(m_zfsStrctrdComm);
    MPI_Group groupRescalingCommGr, newgroupRescalingCommGr;
    ZFSInt rescalingCommGrCommsize = (int)(rescalingCommGrPartitions.size());

    MPI_Comm_group(m_zfsStrctrdComm, &groupRescalingCommGr);
    MPI_Group_incl(groupRescalingCommGr, rescalingCommGrCommsize , rescalingCommGrRanks, &newgroupRescalingCommGr);
    MPI_Comm_create(m_zfsStrctrdComm, newgroupRescalingCommGr,rescalingCommGrComm);

    if(domainId() == rescalingCommGrRanks[0]) {
      MPI_Comm_rank(rescalingCommGrComm[0], &rescalingCommGrRoot[0]);
      MPI_Comm_rank(m_zfsStrctrdComm, &rescalingCommGrRootGlobal[0]);
    }

    MPI_Barrier(m_zfsStrctrdComm);
    MPI_Bcast(rescalingCommGrRoot, 1, MPI_INT, rescalingCommGrRanks[0], m_zfsStrctrdComm);
    MPI_Bcast(rescalingCommGrRootGlobal, 1, MPI_INT, rescalingCommGrRanks[0], m_zfsStrctrdComm);
  }

  MPI_Barrier(m_zfsStrctrdComm);

  //////////////////////////////////////////////////////////////////
  ///////////////// MULTIBLOCK/PERIODIC COMMUNICATION  /////////////
  //////////////////////////////////////////////////////////////////

  vector<ZFSStrctrdWindowMap*> addComm6000Recv; //store the maps for the communication
  vector<ZFSStrctrdWindowMap*> addComm6000Snd;
  vector<ZFSInt> adjacentPartitionBC6000;
  //we do know from the addCommBC which side we possess but we need to search for the opposite side
  //and the partition containing the opposite side!!!

  //correct indices only for local multiblock bc comm maps
  for(ZFSId i=0; i<(ZFSInt)addCommBC.size(); i++) {
    ZFSStrctrdWindowMap* temp= new ZFSStrctrdWindowMap(nDim);
    ZFSInt dimcount=0;
    for(ZFSId dim=0; dim<nDim; dim++) {
      if(addCommBC[i]->start1[dim]!=addCommBC[i]->end1[dim]&&addCommBC[i]->BC==6000) {
        if(addCommBC[i]->step2[addCommBC[i]->order[dim]]>0) {
          addCommBC[i]->start1[dim]-=m_noGhostLayers;
          addCommBC[i]->end1[dim]+=m_noGhostLayers;
          addCommBC[i]->start2[addCommBC[i]->order[dim]]-=m_noGhostLayers;
          addCommBC[i]->end2[addCommBC[i]->order[dim]]+=m_noGhostLayers;
        } else {
          addCommBC[i]->start1[dim]-=m_noGhostLayers;
          addCommBC[i]->end1[dim]+=m_noGhostLayers;
          addCommBC[i]->start2[addCommBC[i]->order[dim]]+=m_noGhostLayers;
          addCommBC[i]->end2[addCommBC[i]->order[dim]]-=m_noGhostLayers;
        }
      }
      if(addCommBC[i]->start1[dim]==addCommBC[i]->end1[dim]&&addCommBC[i]->BC==6000)
        dimcount++;
    }
    if(dimcount==3) {cout<<"error!!!!!! 0d surface found!!  check it .... "<<endl;}

    //now compare the multiblock/periodic bc maps to all partitions
    for(ZFSInt j=0; j<(ZFSInt)m_partitionMapsWithoutGC.size(); j++) {
      mapCombine21(addCommBC[i], m_partitionMapsWithoutGC[j], temp);
      temp->BC=addCommBC[i]->BC;
      temp->Nstar=addCommBC[i]->Nstar;
      temp->SingularId=addCommBC[i]->SingularId;
      temp->dc1=addCommBC[i]->dc1;
      temp->dc2=addCommBC[i]->dc2;
      ZFSBool test=false;

      if(addCommBC[i]->BC==6000)
      {
        if(dimcount==1)
          test=mapCheck2d(temp);
        else if(dimcount==2)
          test=mapCheck1d(temp);
        else if(dimcount==3)
          test=mapCheck0d(temp);
      }
      else
      {
        test=mapCheck(temp); //map does exist
      }
      if(test==true) {
        addComm6000Recv.push_back(temp);
        adjacentPartitionBC6000.push_back(j);
        temp=new ZFSStrctrdWindowMap(nDim);
      }
    }
  }



  //first correct the SndMapBC6000 only for multiblock by GC layers
  for(ZFSId i=0; i<(ZFSInt)SndMapBC6000.size(); i++) {

    ZFSInt dimcount=0;
    for(ZFSId dim=0; dim<nDim; dim++) {
      if(SndMapBC6000[i]->start1[dim]!=SndMapBC6000[i]->end1[dim]&&SndMapBC6000[i]->BC==6000) {
        if(SndMapBC6000[i]->step2[SndMapBC6000[i]->order[dim]]>0) {
          SndMapBC6000[i]->start1[dim]-=m_noGhostLayers;
          SndMapBC6000[i]->end1[dim]+=m_noGhostLayers;
          SndMapBC6000[i]->start2[SndMapBC6000[i]->order[dim]]-=m_noGhostLayers;
          SndMapBC6000[i]->end2[SndMapBC6000[i]->order[dim]]+=m_noGhostLayers;
        } else {
          SndMapBC6000[i]->start1[dim]-=m_noGhostLayers;
          SndMapBC6000[i]->end1[dim]+=m_noGhostLayers;
          SndMapBC6000[i]->start2[SndMapBC6000[i]->order[dim]]+=m_noGhostLayers;
          SndMapBC6000[i]->end2[SndMapBC6000[i]->order[dim]]-=m_noGhostLayers;
        }
      }
      if(SndMapBC6000[i]->start1[dim]==SndMapBC6000[i]->end1[dim]&&SndMapBC6000[i]->BC==6000)
        dimcount++;
    }
    if(dimcount==3) {cout<<"error!!!!!! 0d surface found!!  check it .... "<<endl;}

    //now compare multiblock/periodic bc map (with ghost-cells) with own map (w/o ghost-cells)
    //for(ZFSId i=0; i<(ZFSInt)SndMapBC6000.size(); i++) {
    ZFSStrctrdWindowMap* temp= new ZFSStrctrdWindowMap(nDim);

    mapCombine11(m_myMapWithoutGC, SndMapBC6000[i], temp);
    temp->BC=SndMapBC6000[i]->BC;
    temp->Nstar=SndMapBC6000[i]->Nstar;
    temp->SingularId=SndMapBC6000[i]->SingularId;
    temp->dc1=SndMapBC6000[i]->dc1;
    temp->dc2=SndMapBC6000[i]->dc2;

    ZFSBool test=false;
    if(SndMapBC6000[i]->BC==6000)
    {
      if(dimcount==1)
        test=mapCheck2d(temp);
      else if(dimcount==2)
        test=mapCheck1d(temp);
      else if(dimcount==3)
        test=mapCheck0d(temp);
    }
    else
    {
      test=mapCheck(temp); //map does exist
    }
    if(test==true) {
      addComm6000Snd.push_back(temp);
      temp=new ZFSStrctrdWindowMap(nDim);
    }
  }

  //special treatment for certain cases where number of SND and RCV maps is unequal
  ZFSInt count1[4]={0,0,0,0},count2[4]={0,0,0,0};
  for (ZFSInt i=0;i<(ZFSInt)addComm6000Snd.size();++i) {
    if(addComm6000Snd[i]->BC==6000&&addComm6000Snd[i]->Nstar==-1)
      count1[0]++;
    else if(addComm6000Snd[i]->BC==6000&&addComm6000Snd[i]->Nstar!=-1)
      count1[1]++;
    else if(addComm6000Snd[i]->BC>=4000&&addComm6000Snd[i]->BC<5000&&addComm6000Snd[i]->Nstar==-1)
      count1[2]++;
    else if(addComm6000Snd[i]->BC>=4000&&addComm6000Snd[i]->BC<5000&&addComm6000Snd[i]->Nstar!=-1)
      count1[3]++;
  }
  for (ZFSInt i=0;i<(ZFSInt)addComm6000Recv.size();++i) {
    if(addComm6000Recv[i]->BC==6000&&addComm6000Recv[i]->Nstar==-1)
      count2[0]++;
    else if(addComm6000Recv[i]->BC==6000&&addComm6000Recv[i]->Nstar!=-1)
      count2[1]++;
    else if(addComm6000Recv[i]->BC>=4000&&addComm6000Recv[i]->BC<5000&&addComm6000Recv[i]->Nstar==-1)
      count2[2]++;
    else if(addComm6000Recv[i]->BC>=4000&&addComm6000Recv[i]->BC<5000&&addComm6000Recv[i]->Nstar!=-1)
      count2[3]++;
  }

  MPI_Barrier(m_zfsStrctrdComm);


  //////////////////////////////////////////////////////////////////
  ///////////////// CORRECTIONS FOR BC 6000 (MULTIBLOCK) ///////////
  //////////////////////////////////////////////////////////////////

  //determine faces of receiving maps
  //and make maps three-dimensional
  //RCV maps, no singularity
  for(ZFSInt i=0; i<(ZFSInt)addComm6000Recv.size(); ++i) {
    if(addComm6000Recv[i]->Nstar==-1&&addComm6000Recv[i]->BC==6000) {
      ZFSStrctrdWindowMap* tempRCV= new ZFSStrctrdWindowMap(nDim);
      mapCpy(addComm6000Recv[i], tempRCV);

      //IMPORTANT:
      //the faces are important for the unique send and rcv tags
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Recv[i]->step2[dim]<0) {
          ZFSInt haha= addComm6000Recv[i]->start2[dim];
          addComm6000Recv[i]->start2[dim]=addComm6000Recv[i]->end2[dim];
          addComm6000Recv[i]->end2[dim]=haha;
        }
      }

      //check the side for the receiving parts and save the face
      for(ZFSId part =0; part<(ZFSInt)m_partitionMapsWithoutGC.size(); ++part) {
        if(addComm6000Recv[i]->Id2==m_partitionMapsWithoutGC[part]->Id2) {
          for(ZFSId dim=0; dim<nDim; dim++) {
            if(addComm6000Recv[i]->start2[dim]==addComm6000Recv[i]->end2[dim]) {
              switch(dim){
              case 0:{//face 0 or 1
                      //check if face 0 or 1
                if(addComm6000Recv[i]->start2[dim]==m_partitionMapsWithoutGC[part]->start2[dim]){
                  //test of start2 instead of start1.
                  tempRCV->face=0;
                }else {
                  tempRCV->face=1;
                }
                break;
              }
              case 1:{//face 2 or 3
                      //check if face 2 or 3
                if(addComm6000Recv[i]->start2[dim]==m_partitionMapsWithoutGC[part]->start2[dim]){
                  tempRCV->face=2;
                }else{
                  tempRCV->face=3;
                }
                break;
              }
              case 2:{//face 4 or 5
                      //check if face 4 or 5
                if(addComm6000Recv[i]->start2[dim]==m_partitionMapsWithoutGC[part]->start2[dim]){
                  tempRCV->face=4;
                }else{
                  tempRCV->face=5;
                }
                break;
              }
              default:{
                cerr << "error no side could be attributed" << endl; exit(1);
              }
              }
            }
          }
        }
      }

      // //check how many dimensions the map has
      // ZFSInt mapDimension = 0;
      // for(ZFSId dim=0; dim<nDim; dim++) {
      //        if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
      //          mapDimension++;
      //        }
      // }
      // //this is a 2d face with one zero dimension
      // //in that case extend in both finite dimensions
      // if(mapDimension==1||mapDimension==2) {
      //        for(ZFSId dim=0; dim<nDim; dim++) {
      //          if(addComm6000Recv[i]->start1[dim]!=addComm6000Recv[i]->end1[dim]) {
      //              tempRCV->start1[dim]-=m_noGhostLayers;
      //              tempRCV->end1[dim]+=m_noGhostLayers;
      //          }
      //        }
      // }

      //make the RCV maps three-dimensional
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
          //values are the same:
          //change the send and receive values
          if(addComm6000Recv[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
            //we are on the start side of the domain
            tempRCV->start1[dim]-=m_noGhostLayers;
          } else {
            //we are at the end side of the domain
            tempRCV->end1[dim]+=m_noGhostLayers;
          }
        }
      }

      tempRCV->BC=addComm6000Recv[i]->BC;
      rcvMap.push_back(tempRCV);
    }
  }


  //determine faces of receiving maps
  //and make maps three-dimensional
  //RCV maps, with singularity
  for(ZFSInt i=0; i<(ZFSInt)addComm6000Recv.size(); ++i) {
    if(addComm6000Recv[i]->Nstar!=-1&&addComm6000Recv[i]->BC==6000) {
      //ZFSStrctrdWindowMap* tempSND= new ZFSStrctrdWindowMap(nDim);
      ZFSStrctrdWindowMap* tempRCV= new ZFSStrctrdWindowMap(nDim);
      mapCpy(addComm6000Recv[i], tempRCV);

      //IMPORTANT:
      //the faces are important for the unique send and rcv tags
      ZFSInt dimcount=0;
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
          dimcount++;
        }
      }

      if(dimcount==2) {
        for(ZFSId dim=0; dim<nDim; dim++) {
          if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
            //values are the same:
            //change the send and receive values
            if(addComm6000Recv[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //corner
              //rcv Maps:
              tempRCV->end1[dim]+=1;
            } else {
              //rcv Maps:
              tempRCV->start1[dim]-=1;
            }
          } else {
            //line
            //rcv Maps:
            // tempRCV->start1[dim]-=m_noGhostLayers;
            // tempRCV->end1[dim]  +=m_noGhostLayers;
          }
        }
      } else if(dimcount==3) {
        ZFSInt dimN;
        if(addComm6000Recv[i]->BC>=4400&&addComm6000Recv[i]->BC<4410) {
          addComm6000Recv[i]->face=addComm6000Recv[i]->BC-4401;
        }
        if(addComm6000Recv[i]->face!=-1) {
          dimN=addComm6000Recv[i]->face/2;
        } else {
          cout<<"ERROR!!! point singular communication cannot decide the direction!! check the singularity!!!"<<endl;
          exit(0);
        }

        for(ZFSId dim=0; dim<nDim; dim++) {
          if(dim!=dimN) {
            //values are the same:
            //change the send and receive values
            if(addComm6000Recv[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //corner
              //rcv Maps:
              tempRCV->end1[dim]+=1;

            } else {
              //rcv Maps:
              tempRCV->start1[dim]-=1;
            }

          } else {
            //rcv Maps:
            if(addComm6000Recv[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //rcv Maps:
              tempRCV->start1[dim]-=2;

            } else {
              //rcv Maps:
              tempRCV->end1[dim]+=2;
            }
          }
        }
      }

      for(ZFSInt j=0; j<singularPoint; ++j) {
        ZFSStrctrdWindowMap* temp= new ZFSStrctrdWindowMap(nDim);
        if(localSingularMap[j]->BC==addComm6000Recv[i]->BC)
          mapCombine11(localSingularMap[j], addComm6000Recv[i], temp);
        ZFSBool test=false;
        test=mapCheck(temp); //map does exist
        if(test==true) {

          //   inputBoxId = partition->outputBoxInfo[Id2]->inputBoxID;
          // ZFSInt RecvId = partition->outputBoxInfo[addComm6000Recv[i]->Id2]->cpu;
          ZFSInt Recvblock = partition->outputBoxInfo[addComm6000Recv[i]->Id2]->inputBoxID;
          ZFSInt singnumber =-1;
          for(ZFSInt k=0;k<addComm6000Recv[i]->Nstar-3;k++)
          {
            if (Recvblock== singularity[j].SingularBlockId[k] )
            {
              singnumber=k;
              break;
            }
          }

          // if (Recvblock== singularity[j].SingularBlockId[0] )
          //   singnumber=0;
          // else  if (Recvblock== singularity[j].SingularBlockId[1] )
          //   singnumber=1;
          if(singnumber==-1)
          {
            cout<<"ERROR!!!!!!!!can not find the correct the displacement!!!!!!"<<endl;
            cout<<"recvid2:"<<addComm6000Recv[i]->Id2<<" recblock:"<<Recvblock<<" id:"<<singularity[j].SingularBlockId[0]<<" "<<singularity[j].SingularBlockId[1]<<endl;
          }

          //ZFSInt singnumber=addComm6000Recv[i]->Nstar;
          for(ZFSId dim=0; dim<nDim; dim++) {
            tempRCV->start1[dim]+=singularity[j].displacement[singnumber+2][dim];
            tempRCV->end1[dim]  +=singularity[j].displacement[singnumber+2][dim];
          }
          //  singularity[j].count++;
          tempRCV->SingularId=j;
          break;
        }
      }

      tempRCV->face=100;

      if(addComm6000Recv[i]->BC>=4000&&addComm6000Recv[i]->BC<5000) {
        tempRCV->BC=addComm6000Recv[i]->BC;
      } else {
        tempRCV->BC=6333;
      }
      rcvMap.push_back(tempRCV);
    }
  }

  //determine faces of receiving maps
  //and make maps three-dimensional
  //SND maps, no singularity
  for(ZFSInt i=0; i<(ZFSInt)addComm6000Snd.size(); i++) {
    if(addComm6000Snd[i]->Nstar==-1&&addComm6000Snd[i]->BC==6000) {
      ZFSStrctrdWindowMap* tempSND= new ZFSStrctrdWindowMap(nDim);
      mapCpy(addComm6000Snd[i], tempSND);
      //IMPORTANT:
      //the faces are important for the unique send and rcv tags
      //first check the side for the sending parts
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Snd[i]->start1[dim]==addComm6000Snd[i]->end1[dim]){
          switch(dim){
          case 0:{//face 0 or 1
                  //check if face 0 or 1
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]){
              //test of start2 instead of start1.
              tempSND->face=0;
            }else {
              tempSND->face=1;
            }
            break;
          }
          case 1:{//face 2 or 3
                  //check if face 2 or 3
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]){
              tempSND->face=2;
            }else{
              tempSND->face=3;
            }
            break;
          }
          case 2:{//face 4 or 5
                  //check if face 4 or 5
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]){
              tempSND->face=4;
            }else{
              tempSND->face=5;
            }
            break;
          }
          default:{
            cerr << "error no side could be attributed" << endl; exit(1);
          }
          }
        }
      }

      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Snd[i]->start1[dim]==addComm6000Snd[i]->end1[dim]) {
          //values are the same:
          //change the send and receive values
          if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
            //we are on the start side of the domain
            tempSND->end1[dim]+=m_noGhostLayers;
          } else {
            //we are at the end side of the domain
            tempSND->start1[dim]-=m_noGhostLayers;
          }
        }
      }

      tempSND->BC=addComm6000Snd[i]->BC;
      sndMap.push_back(tempSND);
    }
  }

  //determine faces of receiving maps
  //and make maps three-dimensional
  //SND maps, with singularity
  for(ZFSInt i=0; i<(ZFSInt)addComm6000Snd.size(); i++) {
    if(addComm6000Snd[i]->Nstar!=-1&&addComm6000Snd[i]->BC==6000) {
      ZFSStrctrdWindowMap* tempSND= new ZFSStrctrdWindowMap(nDim);
      mapCpy(addComm6000Snd[i], tempSND);

      ZFSInt dimcount=0;
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
          dimcount++;
        }
      }

      if(dimcount==2) {
        for(ZFSId dim=0; dim<nDim; dim++) {
          if(addComm6000Snd[i]->start1[dim]==addComm6000Snd[i]->end1[dim]) {
            //values are the same:
            //change the send and receive values
            //corner
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //snd Maps:
              tempSND->end1[dim]+=1;
            } else {
              //snd Maps:
              tempSND->start1[dim]-=1;
            }
          } else {
            //line
            //snd Maps:
            // tempSND->start1[dim]-=m_noGhostLayers;
            // tempSND->end1[dim]+=m_noGhostLayers;
          }
        }
      } else if(dimcount==3) {
        ZFSInt dimN;
        if(addComm6000Snd[i]->BC>=4400&&addComm6000Snd[i]->BC<4410) {
          addComm6000Snd[i]->face=addComm6000Snd[i]->BC-4401;
        }
        if(addComm6000Snd[i]->face!=-1) {
          dimN=addComm6000Snd[i]->face/2;
        } else {
          cout<<"erroor!!! point singular communication cannot decide the direction!! check the singylairy!!!!!!"<<endl;
          exit(0);
        }

        for(ZFSId dim=0; dim<nDim; dim++) {
          if(dim!=dimN) {
            //values are the same:
            //change the send and receive values
            //corner
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //rcv Maps:
              tempSND->end1[dim]+=1;
            } else {
              //rcv Maps:
              tempSND->start1[dim]-=1;
            }
          } else{
            //line
            //rcv Maps:
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //rcv Maps:
              tempSND->end1[dim]+=2;

            } else {
              //rcv Maps:
              tempSND->start1[dim]-=2;
            }
          }
        }
      }

      tempSND->face=100;

      if(addComm6000Snd[i]->BC>=4000&&addComm6000Snd[i]->BC<5000) {
        tempSND->BC=addComm6000Snd[i]->BC;
      } else {
        tempSND->BC=6333;
      }

      sndMap.push_back(tempSND);
    }
  }

  //////////////////////////////////////////////////////////////////
  ///////////////// CORRECTIONS FOR BC 4xxx (PERIODIC) /////////////
  //////////////////////////////////////////////////////////////////

  //determine faces of receiving maps
  //and make maps three-dimensional
  //RCV maps, no singularity
  for(ZFSInt i=0; i<(ZFSInt)addComm6000Recv.size(); ++i) {
    if(addComm6000Recv[i]->Nstar==-1&&addComm6000Recv[i]->BC!=6000) {
      ZFSStrctrdWindowMap* tempRCV= new ZFSStrctrdWindowMap(nDim);
      mapCpy(addComm6000Recv[i], tempRCV);

      //IMPORTANT:
      //the faces are important for the unique send and rcv tags
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Recv[i]->step2[dim]<0) {
          ZFSInt tmp= addComm6000Recv[i]->start2[dim];
          addComm6000Recv[i]->start2[dim]=addComm6000Recv[i]->end2[dim];
          addComm6000Recv[i]->end2[dim]=tmp;
        }
      }
      //check the side for the receiving parts
      for(ZFSId part =0; part<(ZFSInt)m_partitionMapsWithoutGC.size(); ++part) {
        if(addComm6000Recv[i]->Id2==m_partitionMapsWithoutGC[part]->Id2) {
          for(ZFSId dim=0; dim<nDim; dim++) {
            if(addComm6000Recv[i]->start2[dim]==addComm6000Recv[i]->end2[dim]) {
              switch(dim){
              case 0:{//face 0 or 1
                      //check if face 0 or 1
                if(addComm6000Recv[i]->start2[dim]==m_partitionMapsWithoutGC[part]->start2[dim]){
                  //test of start2 instead of start1.
                  tempRCV->face=0;
                }else {
                  tempRCV->face=1;
                }
                break;
              }
              case 1:{//face 2 or 3
                      //check if face 2 or 3
                if(addComm6000Recv[i]->start2[dim]==m_partitionMapsWithoutGC[part]->start2[dim]){
                  tempRCV->face=2;
                }else{
                  tempRCV->face=3;
                }
                break;
              }
              case 2:{//face 4 or 5
                      //check if face 4 or 5
                if(addComm6000Recv[i]->start2[dim]==m_partitionMapsWithoutGC[part]->start2[dim]){
                  tempRCV->face=4;
                }else{
                  tempRCV->face=5;
                }
                break;
              }
              default:{
                cerr << "error no side could be attributed" << endl; exit(1);
              }
              }
            }
          }
        }
      }

      //check how many dimensions the map has
      ZFSInt mapDimension = 0;
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
          mapDimension++;
        }
      }

      //this is a 2d face with one zero dimension
      //in that case extend in both finite dimensions
      if(mapDimension==1) {
        for(ZFSId dim=0; dim<nDim; dim++) {
          if(addComm6000Recv[i]->start1[dim]!=addComm6000Recv[i]->end1[dim]) {
            tempRCV->start1[dim]-=m_noGhostLayers;
            tempRCV->end1[dim]+=m_noGhostLayers;
          }
        }
      }

      //now thicken the map:
      //point: thicken in all three dimensions
      //line: thicken in two zero dimensions
      //face: thicken in the only zero dimension
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
          //values are the same:
          //change the send and receive values
          if(addComm6000Recv[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
            //we are on the start side of the domain
            //rcv Maps:
            tempRCV->start1[dim]-=m_noGhostLayers;
          } else {
            //we are at the end side of the domain
            //rcv Maps:
            tempRCV->end1[dim]+=m_noGhostLayers;
          }
        }
      }

      tempRCV->BC=addComm6000Recv[i]->BC;
      rcvMap.push_back(tempRCV);
    }
  }

  //determine faces of receiving maps
  //and make maps three-dimensional
  //RCV maps, with singularity
  for(ZFSInt i=0; i<(ZFSInt)addComm6000Recv.size(); ++i) {
    if(addComm6000Recv[i]->Nstar!=-1&&addComm6000Recv[i]->BC!=6000) {
      ZFSStrctrdWindowMap* tempRCV= new ZFSStrctrdWindowMap(nDim);
      mapCpy(addComm6000Recv[i], tempRCV);

      //IMPORTANT:
      //the faces are important for the unique send and rcv tags
      ZFSInt dimcount=0;
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
          dimcount++;
        }
      }

      if(dimcount==2) {
        for(ZFSId dim=0; dim<nDim; dim++) {
          if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
            //values are the same:
            //change the send and receive values
            //corner
            if(addComm6000Recv[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //rcv Maps:
              tempRCV->end1[dim]+=1;
            } else {
              //rcv Maps:
              tempRCV->start1[dim]-=1;
            }
          } else {
            //line
            //rcv Maps:
            tempRCV->start1[dim]-=m_noGhostLayers;
            tempRCV->end1[dim]  +=m_noGhostLayers;
          }
        }
      } else if(dimcount==3) {
        ZFSInt dimN;
        if(addComm6000Recv[i]->BC>=4400&&addComm6000Recv[i]->BC<4410) {
          addComm6000Recv[i]->face=addComm6000Recv[i]->BC-4401;
        }
        if(addComm6000Recv[i]->face!=-1) {
          dimN=addComm6000Recv[i]->face/2;
        } else{
          cout<<"erroor!!! point singular communication cannot decide the direction!! check the singylairy!!!!!!"<<endl;
          exit(0);
        }

        for(ZFSId dim=0; dim<nDim; dim++) {
          if(dim!=dimN) {
            //values are the same:
            //change the send and receive values
            //corner
            if(addComm6000Recv[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //rcv Maps:
              tempRCV->end1[dim]+=1;
            } else {
              //rcv Maps:
              tempRCV->start1[dim]-=1;
            }
          } else {
            //line
            //rcv Maps:
            if(addComm6000Recv[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //corner
              //rcv Maps:
              tempRCV->start1[dim]-=2;
            } else {
              //rcv Maps:
              tempRCV->end1[dim]+=2;
            }
          }
        }
      }

      for(ZFSInt j=0; j<singularPoint; ++j) {
        ZFSStrctrdWindowMap* temp= new ZFSStrctrdWindowMap(nDim);
        if(localSingularMap[j]->BC==addComm6000Recv[i]->BC)
          mapCombine11(localSingularMap[j], addComm6000Recv[i], temp);
        ZFSBool test=false;
        test=mapCheck(temp); //map does exist
        if(test==true) {
          //communication map change start from 2 to 3
          ZFSInt singnumber=singularity[j].count;
          //ZFSInt singnumber=addComm6000Recv[i]->Nstar;
          for(ZFSId dim=0; dim<nDim; dim++) {
            tempRCV->start1[dim]+=singularity[j].displacement[singnumber][dim];
            tempRCV->end1[dim]  +=singularity[j].displacement[singnumber][dim];
          }
          singularity[j].count++;
          tempRCV->SingularId=j;
          break;
        }
      }

      tempRCV->face=100;

      if(addComm6000Recv[i]->BC>=4000&&addComm6000Recv[i]->BC<5000) {
        tempRCV->BC=addComm6000Recv[i]->BC;
      } else {
        tempRCV->BC=6333;
      }

      rcvMap.push_back(tempRCV);
    }
  }

  //determine faces of receiving maps
  //and make maps three-dimensional
  //SND maps, no singularity
  for(ZFSInt i=0; i<(ZFSInt)addComm6000Snd.size(); i++) {
    if(addComm6000Snd[i]->Nstar==-1&&addComm6000Snd[i]->BC!=6000) {
      ZFSStrctrdWindowMap* tempSND= new ZFSStrctrdWindowMap(nDim);
      mapCpy(addComm6000Snd[i], tempSND);
      //IMPORTANT:
      //the faces are important for the unique send and rcv tags

      //first check the side for the sending parts
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Snd[i]->start1[dim]==addComm6000Snd[i]->end1[dim]) {
          switch(dim){
          case 0:{//face 0 or 1
                  //check if face 0 or 1
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //test of start2 instead of start1.
              tempSND->face=0;
            }else {
              tempSND->face=1;
            }
            break;
          }
          case 1:{//face 2 or 3
                  //check if face 2 or 3
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              tempSND->face=2;
            }else{
              tempSND->face=3;
            }
            break;
          }
          case 2:{//face 4 or 5
                  //check if face 4 or 5
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              tempSND->face=4;
            }else{
              tempSND->face=5;
            }
            break;
          }
          default:{
            cerr << "error no side could be attributed" << endl; exit(1);
          }
          }
        }
      }

      //check how many dimensions the map has
      ZFSInt mapDimension = 0;
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Snd[i]->start1[dim]==addComm6000Snd[i]->end1[dim]) {
          mapDimension++;
        }
      }

      //this is a 2d face with one zero dimension
      //in that case extend in both finite dimensions
      if(mapDimension==1) {
        for(ZFSId dim=0; dim<nDim; dim++) {
          if(addComm6000Snd[i]->start1[dim]!=addComm6000Snd[i]->end1[dim]) {
            // if(tempSND->end1[dim]-tempSND->start1[dim]!=m_noGhostLayers) {
            tempSND->start1[dim]-=m_noGhostLayers;
            tempSND->end1[dim]+=m_noGhostLayers;
            // }
          }
        }
      }

      //now thicken the map:
      //point: thicken in all three dimensions
      //line: thicken in two zero dimensions
      //face: thicken in the only zero dimension
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Snd[i]->start1[dim]==addComm6000Snd[i]->end1[dim]) {
          //values are the same:
          //change the send and receive values
          if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
            //we are on the start side of the domain
            //snd Maps:
            tempSND->end1[dim]+=m_noGhostLayers;
          } else {
            //we are at the end side of the domain
            //snd Maps:
            tempSND->start1[dim]-=m_noGhostLayers;
          }
        }
      }

      tempSND->BC=addComm6000Snd[i]->BC;
      sndMap.push_back(tempSND);
    }
  }

  //determine faces of receiving maps
  //and make maps three-dimensional
  //SND maps, with singularity
  for(ZFSInt i=0; i<(ZFSInt)addComm6000Snd.size(); i++) {
    if(addComm6000Snd[i]->Nstar!=-1&&addComm6000Snd[i]->BC!=6000) {
      ZFSStrctrdWindowMap* tempSND= new ZFSStrctrdWindowMap(nDim);
      mapCpy(addComm6000Snd[i], tempSND);

      ZFSInt dimcount=0;
      for(ZFSId dim=0; dim<nDim; dim++) {
        if(addComm6000Recv[i]->start1[dim]==addComm6000Recv[i]->end1[dim]) {
          dimcount++;
        }
      }

      if(dimcount==2) {
        for(ZFSId dim=0; dim<nDim; dim++) {
          if(addComm6000Snd[i]->start1[dim]==addComm6000Snd[i]->end1[dim]) {
            //values are the same:
            //change the send and receive values
            //corner
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //snd Maps:
              tempSND->end1[dim]+=1;
            } else {
              //snd Maps:
              tempSND->start1[dim]-=1;
            }
          } else {
            //line
            //snd Maps:
            tempSND->start1[dim]-=m_noGhostLayers;
            tempSND->end1[dim]+=m_noGhostLayers;
          }
        }
      } else if(dimcount==3) {
        ZFSInt dimN;
        if(addComm6000Snd[i]->BC>=4400&&addComm6000Snd[i]->BC<4410) {
          addComm6000Snd[i]->face=addComm6000Snd[i]->BC-4401;
        }
        if(addComm6000Snd[i]->face!=-1) {
          dimN=addComm6000Snd[i]->face/2;
        } else {
          cout<<"erroor!!! point singular communication cannot decide the direction!! check the singylairy!!!!!!"<<endl;
          exit(0);
        }

        for(ZFSId dim=0; dim<nDim; dim++) {
          if(dim!=dimN) {
            //values are the same:
            //change the send and receive values
            //corner
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //rcv Maps:
              tempSND->end1[dim]+=1;
            } else {
              //rcv Maps:
              tempSND->start1[dim]-=1;
            }
          } else {
            //line
            //rcv Maps:
            if(addComm6000Snd[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
              //corner
              //rcv Maps:
              tempSND->end1[dim]+=2;
            } else {
              //rcv Maps:
              tempSND->start1[dim]-=2;
            }
          }
        }
      }

      tempSND->face=100;

      if(addComm6000Snd[i]->BC>=4000&&addComm6000Snd[i]->BC<5000) {
        tempSND->BC=addComm6000Snd[i]->BC;
      } else {
        tempSND->BC=6333;
      }

      sndMap.push_back(tempSND);
    }
  }

  ////////////////////////////////////////////////////
  /////////// SINGULARITY CORRECTION  ////////////////
  ////////////////////////////////////////////////////

  for(ZFSInt i=0; i<singularPoint; ++i) {
    for(ZFSId dim=0; dim<nDim; dim++) {
      singularity[i].Viscous[dim]=0;

      if(singularity[i].start[dim]==singularity[i].end[dim]) {
        //values are the same:
        if(singularity[i].start[dim]==m_myMapWithoutGC->start2[dim]) {
          //corner
          singularity[i].end[dim]+=1;
          singularity[i].Viscous[dim]=-1;
        } else {
          singularity[i].start[dim]-=1;
        }
      } else {
        singularity[i].start[dim]-=m_noGhostLayers;
        singularity[i].end[dim]  +=m_noGhostLayers;
      }
    }
  }

  //=============> now the snd and rcv maps contain the following order
  //1) all the communication because of partitioning
  //2) communication to other partition/block because of multiblock 6000 bc or periodic boundary condition 4xxx

  ///////////////////////////////////////////////////////////
  //////////////// AUX DATA MAP CORRECTION //////////////////
  ///////////////////////////////////////////////////////////

  //create the auxilary data maps
  //needed to write out cf,cp etc ...
  int counter=0;
  for(ZFSInt i=0; i<(ZFSInt)globalStrctrdBndryCndMaps.size(); i++) {
    int firstDigit=(int)(((double) globalStrctrdBndryCndMaps[i]->BC)/1000.0);
    if(firstDigit==1) {//all the wall should be covered
      mapCombine11(m_myMapWithoutGC, globalStrctrdBndryCndMaps[i], localMapDummy );
      ZFSBool test=false;
      test=mapCheck(localMapDummy);
      if(test==true){
        physicalAuxDataMap.push_back(localMapDummy);
        for(int j=0; j<nDim; ++j){
          //correct the global part for the output!!!
          physicalAuxDataMap[counter]->start2[j]-=globalStrctrdBndryCndMaps[i]->start2[j];
          physicalAuxDataMap[counter]->end2[j]-=globalStrctrdBndryCndMaps[i]->start2[j];
        }
        counter++;
        localMapDummy=new ZFSStrctrdWindowMap(nDim);
      }
    }
  }

  for(ZFSInt i=0; i<(ZFSInt)physicalAuxDataMap.size(); i++) {
    for(ZFSId dim=0; dim<nDim; dim++) {
      //check if the index is the same
      //MPI_Barrier(m_zfsStrctrdComm);
      if(physicalAuxDataMap[i]->start1[dim]==physicalAuxDataMap[i]->end1[dim]) {
        //check which face it is and save this information
        switch(dim)
        {
        case 0:{//face 0 or 1
          //check if face 0 or 1
          if(physicalAuxDataMap[i]->start1[dim]==m_myMapWithoutGC->start2[dim]){
            //test of start2 instead of start1.
            physicalAuxDataMap[i]->face=0;
          }else {
            physicalAuxDataMap[i]->face=1;
          }
          break;
        }
        case 1:{//face 2 or 3
          //check if face 2 or 3
          if(physicalAuxDataMap[i]->start1[dim]==m_myMapWithoutGC->start2[dim]){
            physicalAuxDataMap[i]->face=2;
          }else{
            physicalAuxDataMap[i]->face=3;
          }
          break;
        }
        case 2:{//face 4 or 5
          //check if face 4 or 5
          if(physicalAuxDataMap[i]->start1[dim]==m_myMapWithoutGC->start2[dim]){
            physicalAuxDataMap[i]->face=4;
          }else{
            physicalAuxDataMap[i]->face=5;
          }
          break;
        }
        default:{
          cerr << "error no side could be attributed" << endl; exit(1);
        }

        }
        continue;
      }
    }
  }



  ///////////////////////////////////////////////////////////
  //////////////// PHYSICAL BC MAP CORRECTION ///////////////
  ///////////////////////////////////////////////////////////
  //correct all the start2/end2 indices on the physical maps
  //to get correct results from mapCombine
  for(ZFSInt i=0; i<(ZFSInt)physicalBCMap.size(); i++) {
    for(ZFSId dim=0; dim<nDim; dim++) {
      physicalBCMap[i]->start2[dim] += m_noGhostLayers;
      physicalBCMap[i]->end2[dim] += m_noGhostLayers;
    }
  }

  // until now, all physical BCs are still 2D or 1D,
  // now thicken them to 3D and extend them if necessary
  for(ZFSInt i=0; i<(ZFSInt)physicalBCMap.size(); i++) {
    ZFSInt addDimensionCount = 0;
    for(ZFSId dim=0; dim<nDim; dim++) {
      if(physicalBCMap[i]->start1[dim]==physicalBCMap[i]->end1[dim]) {
        //check which face it is and save this information
        //also thicken in the face normal direction by two ghost-cells
        switch(dim)
        {
        case 0:
        {//face 0 or 1
          //check if face 0 or 1
          if(physicalBCMap[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
            physicalBCMap[i]->face=0;
            physicalBCMap[i]->start1[dim]-=m_noGhostLayers;
            physicalBCMap[i]->start2[dim]-=m_noGhostLayers;
          } else {
            physicalBCMap[i]->face=1;
            physicalBCMap[i]->end1[dim]+=m_noGhostLayers;
            physicalBCMap[i]->end2[dim]+=m_noGhostLayers;
          }
          break;
        }
        case 1:
        {//face 2 or 3
          if(physicalBCMap[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
            //check if face 2 or 3
            physicalBCMap[i]->face=2;
            physicalBCMap[i]->start1[dim]-=m_noGhostLayers;
            physicalBCMap[i]->start2[dim]-=m_noGhostLayers;
          } else {
            physicalBCMap[i]->face=3;
            physicalBCMap[i]->end1[dim]+=m_noGhostLayers;
            physicalBCMap[i]->end2[dim]+=m_noGhostLayers;
          }
          break;
        }
        case 2:
        {//face 4 or 5
          if(physicalBCMap[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
            //check if face 4 or 5
            physicalBCMap[i]->face=4;
            physicalBCMap[i]->start1[dim]-=m_noGhostLayers;
            physicalBCMap[i]->start2[dim]-=m_noGhostLayers;
          } else {
            physicalBCMap[i]->face=5;
            physicalBCMap[i]->end1[dim]+=m_noGhostLayers;
            physicalBCMap[i]->end2[dim]+=m_noGhostLayers;
          }
          break;
        }
        default:
        {
          cerr << "error no side could be attributed" << endl; exit(1);
        }
        }

        addDimensionCount++;
        continue;
      }

      //if the map starts at the begin of the domain, extend in negative direction
      if(physicalBCMap[i]->start1[dim]==m_myMapWithoutGC->start2[dim]) {
        physicalBCMap[i]->start1[dim]-=m_noGhostLayers;
        physicalBCMap[i]->start2[dim]-=m_noGhostLayers;
      }

      //if the map end at the end of the domain, extend in positive direction
      if(physicalBCMap[i]->end1[dim]==m_myMapWithoutGC->end2[dim]) {
        physicalBCMap[i]->end1[dim]+=m_noGhostLayers;
        physicalBCMap[i]->end2[dim]+=m_noGhostLayers;
      }
    }

    physicalBCMap[i]->originShape = addDimensionCount;
  }

  ///////////////////////////////////////////////////////////
  ///////// EXTENSION CORRECTION ////////////////////////////
  ///////////////////////////////////////////////////////////

  //remove the extension in the cases where a partition boundary
  //coincides with two different BCs
  //(if we don't remove it, there might be two overlapping BCs with
  // uncertainty which one is applied first and second, this
  // can cause oscillations and NANs eventually)

  //go over all physical maps
  for(ZFSInt i=0; i<(ZFSInt)physicalBCMap.size(); i++) {
    //first loop: only maps that were originally surfaces
    if(physicalBCMap[i]->originShape==2) {
      continue;
    }
    for(ZFSInt j=0; j<(ZFSInt)physicalBCMap.size(); j++) {
      //second loop: only maps that were originally lines
      if(physicalBCMap[j]->originShape!=2) {
        continue;
      }
      //skip if BC number is equal (we can leave the extension then)
      if(physicalBCMap[i]->BC==physicalBCMap[j]->BC) {
        continue;
      }

      mapCombine11(physicalBCMap[i], physicalBCMap[j],localMapDummy);
      ZFSBool is3dMap = mapCheck3d(localMapDummy);

      if(is3dMap) {
        ZFSInt faceDim = physicalBCMap[i]->face/2;

        for(ZFSId dim=0; dim<nDim; dim++) {
          //skip if is the face normal direction
          if(dim==faceDim) {
            continue;
          }

          //start of domain
          if(physicalBCMap[j]->start1[dim]==physicalBCMap[i]->start1[dim]&&
             physicalBCMap[j]->end1[dim]==m_myMapWithoutGC->start2[dim]) {
            physicalBCMap[i]->start1[dim] += m_noGhostLayers;
            cout << "DomainID: " << domainId() << " cutting extension at start of partition because of partition/bc conflict with BCs " << physicalBCMap[i]->BC << " and " << physicalBCMap[j]->BC << endl;
          }

          //end of domain
          if(physicalBCMap[j]->start1[dim]==m_myMapWithoutGC->end2[dim]&&
             physicalBCMap[j]->end1[dim]==physicalBCMap[i]->end1[dim]) {
            physicalBCMap[i]->end1[dim] -= m_noGhostLayers;
            cout << "DomainID: " << domainId() << " cutting extension at end of partition because of partition/bc conflict with BCs " << physicalBCMap[i]->BC << " and " << physicalBCMap[j]->BC << endl;
          }
        }
      }
    }
  }


  //correct all the start2/end2 indices on the rcv maps
  Id = partition->outputBoxInfo[domainId()]->cpu;
  for(ZFSInt j=0; j<(ZFSInt)rcvMap.size(); j++) {
    for(ZFSId dim=0; dim<nDim; dim++) {
      rcvMap[j]->start2[dim] = rcvMap[j]->start1[dim] + partition->outputBoxInfo[Id]->offset[nDim-1-dim];
      rcvMap[j]->end2[dim] = rcvMap[j]->end1[dim] + partition->outputBoxInfo[Id]->offset[nDim-1-dim];
    }
  }

  for(ZFSInt i=0; i<(ZFSInt)physicalBCMap.size(); i++) {
    if(physicalBCMap[i]->BC==2501) {
      continue;
    }
    for(ZFSInt j=0; j<(ZFSInt)rcvMap.size(); j++) {
      mapCombine11(rcvMap[j],physicalBCMap[i],localMapDummy);
      ZFSBool is3dMap = mapCheck3d(localMapDummy);
      if(is3dMap&&rcvMap[j]->BC==6000) {
        // cout << "############ ATTENTION: BC MAP IS OVERLAPPING WITH 6000 RCV MAP!!! #############" << endl;
        // cout << "receiveMap: " << endl;
        // mapPrint(rcvMap[j]);
        // cout << "physicalBCMap: " << endl;
        // mapPrint(physicalBCMap[i]);
        // cout << "combined21: " << endl;
        // mapPrint(localMapDummy);
      }
    }
  }

  delete localMapDummy;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::createWaveWindowMapping(ZFSInt waveZeroPos)
{
  for(ZFSId j=0; j<(ZFSId)waveRcvMap.size(); j++) {
    delete waveRcvMap[j];
    waveRcvMap[j] = NULL;
  }

  for(ZFSId j=0; j<(ZFSId)waveSndMap.size(); j++) {
    delete waveSndMap[j];
    waveSndMap[j] = NULL;
  }

  waveRcvMap.clear();
  waveSndMap.clear();

  ZFSInt start1[3]={0,0,0};
  ZFSInt end1[3]={0,0,0};
  ZFSInt step1[3]={1,1,1};
  ZFSInt order[3]={0,1,2};
  ZFSInt start2[3]={0,0,0};
  ZFSInt end2[3]={0,0,0};
  ZFSInt offsetCells[3]={0,0,0};
  ZFSInt activeCells[3]={0,0,0};

  ZFSInt Id = partition->outputBoxInfo[domainId()]->cpu;
  ZFSInt inputBoxId = partition->outputBoxInfo[Id]->inputBoxID;

  /////////////////////////////////////////////////////////////
  //////////////// CREATE OWN/PARTITION MAPS //////////////////
  /////////////////////////////////////////////////////////////
  //this is the map of the active cells (no ghostcells) in my own partition
  //shifted by the no of ghost-cells
  for(ZFSInt i=0; i<nDim; i++) {
    start1[i]=partition->outputBoxInfo[Id]->offset[nDim-1-i];//shifted by the number of ghost layers
    start2[i]=0;//shifted by the number of ghost layers
    end1[i]=start1[i]+partition->outputBoxInfo[Id]->DirLast[nDim-1-i]-1;
    end2[i]=start2[i]+partition->outputBoxInfo[Id]->DirLast[nDim - 1 -i]-1;
  }

  ZFSStrctrdWindowMap* waveMyMap = new ZFSStrctrdWindowMap(nDim);
  mapCreate(inputBoxId, start1, end1, step1, Id , start2, end2, step1, order,-1, waveMyMap);

  ZFSInt allCells[3] = {0,0,0};
  for(ZFSInt dim=0; dim<nDim; dim++) {
    allCells[dim]    = partition->inputBoxInfo[0]->DirLast[nDim-1-dim];
  }

  vector<ZFSStrctrdWindowMap*> waveSndPartitionMaps;
  ZFSStrctrdWindowMap* localMapDummy = NULL;
  for(ZFSId j=0; j<noDomains(); j++) {
    localMapDummy= new ZFSStrctrdWindowMap(nDim);
    const ZFSInt Id2 = partition->outputBoxInfo[j]->cpu;
    inputBoxId = partition->outputBoxInfo[Id2]->inputBoxID;

    for(ZFSInt dim=0; dim<nDim; dim++) {
      offsetCells[dim] = partition->outputBoxInfo[Id2]->offset[nDim-1-dim];
      activeCells[dim] = partition->outputBoxInfo[Id2]->DirLast[nDim-1-dim]-1;
    }

    for(ZFSInt dim=0; dim<nDim; dim++) {
      start1[dim]=offsetCells[dim];
      end1[dim]=start1[dim]+activeCells[dim];
    }

    mapCreate(inputBoxId, start1, end1, step1, Id2, start2, end2, step1, order, 0, localMapDummy);
    waveSndPartitionMaps.push_back(localMapDummy);
  }

  vector<ZFSStrctrdWindowMap*> wavePartitionMaps;
  //maps of all the partitions moved to the relative system
  for(ZFSId j=0; j<noDomains(); j++) {
    localMapDummy= new ZFSStrctrdWindowMap(nDim);

    //preparation
    ZFSInt shiftedOffset[3] = {0,0,0};
    ZFSInt croppedActiveCells[3]={0,0,0};
    ZFSInt overhangCells[3] = {0,0,0};
    ZFSInt overhangOffset[3] = {0,0,0};
    const ZFSInt Id2 = partition->outputBoxInfo[j]->cpu;
    inputBoxId = partition->outputBoxInfo[Id2]->inputBoxID;
    for(ZFSInt dim=0; dim<nDim; dim++) {
      offsetCells[dim] = partition->outputBoxInfo[Id2]->offset[nDim-1-dim];
      activeCells[dim] = partition->outputBoxInfo[Id2]->DirLast[nDim-1-dim]-1;
    }

    //first fill with the unmoved values
    //important for 0- and 1-direction
    for(ZFSInt dim=0; dim<nDim; dim++) {
      shiftedOffset[dim] = offsetCells[dim];
      croppedActiveCells[dim] = activeCells[dim];
      overhangCells[dim] = activeCells[dim];
      overhangOffset[dim] = offsetCells[dim];
    }

    if(offsetCells[2]+activeCells[2]<=waveZeroPos) {
      //domain is before wave zero
      //only shift domain
      shiftedOffset[2] = allCells[2]-waveZeroPos+offsetCells[2];
      for(ZFSInt dim=0; dim<nDim; dim++) {
        start1[dim]=shiftedOffset[dim];
        end1[dim]=start1[dim]+activeCells[dim];
      }
      mapCreate(inputBoxId, start1, end1, step1, Id2, start2, end2, step1, order, 0, localMapDummy);
      wavePartitionMaps.push_back(localMapDummy);
    } else if(offsetCells[2]<waveZeroPos && offsetCells[2]+activeCells[2]>waveZeroPos) {
      //create first map
      croppedActiveCells[2] = waveZeroPos-offsetCells[2];
      shiftedOffset[2] = allCells[2]-croppedActiveCells[2];
      for(ZFSInt dim=0; dim<nDim; dim++) {
        start1[dim]=shiftedOffset[dim];
        end1[dim]=start1[dim]+croppedActiveCells[dim];
      }
      mapCreate(inputBoxId, start1, end1, step1, Id2, start2, end2, step1, order, 0, localMapDummy);
      wavePartitionMaps.push_back(localMapDummy);
      localMapDummy= new ZFSStrctrdWindowMap(nDim);
      //create second map
      overhangCells[2] = activeCells[2]-croppedActiveCells[2];
      overhangOffset[2] = 0;
      for(ZFSInt dim=0; dim<nDim; dim++) {
        start1[dim]=overhangOffset[dim];
        end1[dim]=start1[dim]+overhangCells[dim];
      }
      mapCreate(inputBoxId, start1, end1, step1, Id2, start2, end2, step1, order, 1, localMapDummy);
      wavePartitionMaps.push_back(localMapDummy);
    } else {

      shiftedOffset[2] = offsetCells[2]-waveZeroPos;
      for(ZFSInt dim=0; dim<nDim; dim++) {
        start1[dim]=shiftedOffset[dim];
        end1[dim]=start1[dim]+activeCells[dim];
      }
      mapCreate(inputBoxId, start1, end1, step1, Id2, start2, end2, step1, order, 1, localMapDummy);
      wavePartitionMaps.push_back(localMapDummy);
    }
  }

  //////////////////////////////////////////////////////////////////////////////////
  ////////////////// USE MAPS TO CHECK FOR OVERLAPPING PARTS ///////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  //check overlapping of own non-gc partition with other gc partitions and put into SND maps
  localMapDummy= new ZFSStrctrdWindowMap(nDim);
  for(ZFSId i=0; i<(ZFSId)wavePartitionMaps.size(); i++) {
    //only do this for own maps
    if(wavePartitionMaps[i]->Id2 != Id) {
      continue;
    }

    for(ZFSId j=0; j<(ZFSId)waveSndPartitionMaps.size(); j++) {
      mapCombineWave(wavePartitionMaps[i], waveSndPartitionMaps[j], localMapDummy );
      ZFSBool test=false;
      test=mapCheckWave(localMapDummy);
      if(test==true) {
        localMapDummy->BC = wavePartitionMaps[i]->BC;
        if(localMapDummy->start1[2] < allCells[2]-waveZeroPos) {
          const ZFSInt translationK = allCells[2]-waveZeroPos-localMapDummy->start1[2];
          const ZFSInt sizeK = localMapDummy->end1[2] - localMapDummy->start1[2];
          localMapDummy->start1[2] = allCells[2]-translationK;
          localMapDummy->end1[2] = localMapDummy->start1[2]+sizeK;
        } else {
          const ZFSInt translationK = allCells[2]-waveZeroPos;
          localMapDummy->start1[2] -= translationK;
          localMapDummy->end1[2] -= translationK;
        }
        waveSndMap.push_back(localMapDummy);
        localMapDummy=new ZFSStrctrdWindowMap(nDim);
      }
    }
  }

  //check overlapping of own gc partition with other non-gc partitions and put into RCV maps
  for(ZFSId i=0; i<(ZFSId)wavePartitionMaps.size(); i++) {
    mapCombineWave(waveMyMap, wavePartitionMaps[i], localMapDummy );
    ZFSBool test=false;
    test=mapCheckWave(localMapDummy);
    if(test==true) {
      waveRcvMap.push_back(localMapDummy);
      localMapDummy=new ZFSStrctrdWindowMap(nDim);
    }
  }

  for(ZFSId i=0; i<(ZFSInt)waveRcvMap.size(); i++) {
    for(ZFSInt dim=0; dim<nDim; dim++) {
      const ZFSId offset = partition->outputBoxInfo[Id]->offset[nDim-1-dim];
      waveRcvMap[i]->start1[dim] = waveRcvMap[i]->start1[dim] - offset + m_noGhostLayers;
      waveRcvMap[i]->end1[dim] = waveRcvMap[i]->end1[dim] - offset + m_noGhostLayers;
    }
  }

  for(ZFSId i=0; i<(ZFSInt)waveSndMap.size(); i++) {
    for(ZFSInt dim=0; dim<nDim; dim++) {
      const ZFSId offset = partition->outputBoxInfo[Id]->offset[nDim-1-dim];
      waveSndMap[i]->start1[dim] = waveSndMap[i]->start1[dim] - offset + m_noGhostLayers;
      waveSndMap[i]->end1[dim] = waveSndMap[i]->end1[dim] - offset + m_noGhostLayers;
    }
  }

  for(ZFSId j=0; j<(ZFSId)waveSndPartitionMaps.size(); j++) {
    delete waveSndPartitionMaps[j];
    waveSndPartitionMaps[j] = NULL;
  }

  for(ZFSId j=0; j<(ZFSId)wavePartitionMaps.size(); j++) {
    delete wavePartitionMaps[j];
    wavePartitionMaps[j] = NULL;
  }

  delete localMapDummy;
  delete waveMyMap;
  waveSndPartitionMaps.clear();
  wavePartitionMaps.clear();
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::createWaveCommunicationExchangeFlags(ZFSStrctrdWaveCommunicationHandle* comm,
                                                                          ZFSInt noVariables )
{
  const ZFSInt noNghbrsRcv = waveRcvMap.size();
  const ZFSInt noNghbrsSnd = waveSndMap.size();

  // cout << "noNghbrsRcv: " << noNghbrsRcv << " noNghbrsSnd: " << noNghbrsSnd << endl;
  // ZFSInt totalRcv = 0;
  // ZFSInt totalSnd = 0;
  // MPI_Allreduce(&noNghbrsRcv, &totalRcv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  // MPI_Allreduce(&noNghbrsSnd, &totalSnd, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  // cout << " Total rcv: " << totalRcv << " total snd: " << totalSnd << endl;

  comm->noNghbrDomainsRcv=noNghbrsRcv;
  comm->noNghbrDomainsSnd=noNghbrsSnd;

  comm->m_sndNghbrId=new ZFSInt[comm->noNghbrDomainsSnd]; //partitionNumber of neighbour Domain
  comm->m_rcvNghbrId= new ZFSInt[comm->noNghbrDomainsRcv]; //rank of neighbour Domain
  comm->m_noNghbrDomainCellBufferSizeSnd = new ZFSInt[comm->noNghbrDomainsSnd]; //cell size Snd
  comm->m_noNghbrDomainCellBufferSizeRcv = new ZFSInt[comm->noNghbrDomainsRcv]; //cell size Rcv
  comm->startInfoSNDcells = new ZFSInt *[comm->noNghbrDomainsSnd]; //start cells for snd
  comm->endInfoSNDcells = new ZFSInt *[comm->noNghbrDomainsSnd]; //end cells for snd
  comm->startInfoRCVcells = new ZFSInt *[comm->noNghbrDomainsRcv]; //start cells for rcv
  comm->endInfoRCVcells = new ZFSInt *[comm->noNghbrDomainsRcv]; //end cells for rcv

  for(ZFSInt i=0; i<comm->noNghbrDomainsSnd; i++) {
    comm->startInfoSNDcells[i]= new ZFSInt [nDim];
    comm->endInfoSNDcells[i]=new ZFSInt[nDim];
  }

  for(ZFSInt i=0; i<comm->noNghbrDomainsRcv; i++) {
    comm->startInfoRCVcells[i]= new ZFSInt [nDim];
    comm->endInfoRCVcells[i]=new ZFSInt[nDim];
  }

  comm->m_tagHelperSND=NULL;
  comm->m_tagHelperRCV=NULL;

  ZFSInt tempId = 0;
  comm->m_tagHelperSND=new ZFSInt[comm->noNghbrDomainsSnd];
  comm->m_tagHelperRCV=new ZFSInt[comm->noNghbrDomainsRcv];

  ZFSInt accBufferSnd = 0;
  ZFSInt accBufferRcv = 0;

  //SND Maps
  for(ZFSInt i=0; i<comm->noNghbrDomainsSnd; i++) {
    tempId=waveSndMap[i]->Id2;
    //rank != cpuId ==> find the corresponding rank
    for(ZFSInt j=0; j<m_noPartitions; j++) {
      if(tempId==partition->outputBoxInfo[j]->cpu) {
        comm->m_sndNghbrId[i]=j;
        break;
      }
    }

    memcpy(comm->startInfoSNDcells[i], waveSndMap[i]->start1, nDim*sizeof(ZFSInt));
    memcpy(comm->endInfoSNDcells[i], waveSndMap[i]->end1, nDim*sizeof(ZFSInt));

    //compute the buffersizes for sending and receiving
    ZFSInt cellproductSnd=1;
    for(ZFSInt j=0; j<nDim; j++) {
      ZFSInt cellsizesSnd=(comm->endInfoSNDcells[i][j] - comm->startInfoSNDcells[i][j]);
      if(cellsizesSnd!=0) {
        cellproductSnd*=(cellsizesSnd);
      }
    }

    comm->m_noNghbrDomainCellBufferSizeSnd[i]=cellproductSnd*noVariables;
    accBufferSnd += comm->m_noNghbrDomainCellBufferSizeSnd[i];
    comm->m_tagHelperSND[i]=waveSndMap[i]->BC;
  }

  //RCV Maps
  for(ZFSInt i=0; i<comm->noNghbrDomainsRcv; i++) {
    tempId=waveRcvMap[i]->Id2;
    //rank != cpuId ==> find the corresponding rank
    for(ZFSInt j=0; j<m_noPartitions; j++) {
      if(tempId==partition->outputBoxInfo[j]->cpu) {
        comm->m_rcvNghbrId[i]=j;
        break;
      }
    }

    memcpy(comm->startInfoRCVcells[i], waveRcvMap[i]->start1, nDim*sizeof(ZFSInt));
    memcpy(comm->endInfoRCVcells[i], waveRcvMap[i]->end1, nDim*sizeof(ZFSInt));

    //compute the buffersizes for sending and receiving
    ZFSInt cellproductRcv=1;
    for(ZFSInt j=0; j<nDim; j++) {
      ZFSInt cellsizesRcv=(comm->endInfoRCVcells[i][j] - comm->startInfoRCVcells[i][j]);
      if(cellsizesRcv!=0) {
        cellproductRcv*=(cellsizesRcv);
      }
    }

    comm->m_noNghbrDomainCellBufferSizeRcv[i]=cellproductRcv*noVariables;
    accBufferRcv += comm->m_noNghbrDomainCellBufferSizeRcv[i];
    comm->m_tagHelperRCV[i]=waveRcvMap[i]->BC;
  }

  // ZFSInt totalRcv = 0;
  // ZFSInt totalSnd = 0;
  // MPI_Allreduce(&accBufferRcv, &totalRcv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  // MPI_Allreduce(&accBufferSnd, &totalSnd, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  // cout << " Total rcv: " << totalRcv << " total snd: " << totalSnd << endl;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::createCommunicationExchangeFlags(ZFSStrctrdCommunicationHandle* comm,
                                                                  ZFSInt noVariables )
{

  ZFSInt label=0;
  if(rcvMap.size()!=sndMap.size()) {
    ZFSInt count1[4]={0,0,0,0},count2[4]={0,0,0,0};
    for (ZFSInt i=0;i<(ZFSInt)sndMap.size();++i) {
      if(sndMap[i]->BC==6000&&sndMap[i]->Nstar==-1)
        count1[0]++;
      else if(sndMap[i]->BC==6000&&sndMap[i]->Nstar!=-1)
        count1[1]++;
      else if(sndMap[i]->BC>=4000&&sndMap[i]->BC<5000&&sndMap[i]->Nstar==-1)
        count1[2]++;
      else if(sndMap[i]->BC>=4000&&sndMap[i]->BC<5000&&sndMap[i]->Nstar!=-1)
        count1[3]++;
    }
    for (ZFSInt i=0;i<(ZFSInt)rcvMap.size();++i) {
      if(rcvMap[i]->BC==6000&&rcvMap[i]->Nstar==-1)
        count2[0]++;
      else if(rcvMap[i]->BC==6000&&rcvMap[i]->Nstar!=-1)
        count2[1]++;
      else if(rcvMap[i]->BC>=4000&&rcvMap[i]->BC<5000&&rcvMap[i]->Nstar==-1)
        count2[2]++;
      else if(rcvMap[i]->BC>=4000&&rcvMap[i]->BC<5000&&rcvMap[i]->Nstar!=-1)
        count2[3]++;
    }

    //sending and receiving does not correspond ==> error
    cout<<"***********************************WARNING************************************** "<<endl;
    cout<<"this error may be caused by the partition when you have a Step in the grid, try to change number of CPUs and run again!"<<endl;
    cout<<"******************************************************************************** "<<endl;
    //zfsTerm(1,__CALLING_FUNCTION__, "number of sending and receiving processes does not match");
    label=1;
  }

  if(label==1) {
    cout<<"CPU: "<<domainId()<<"  2@@@@@@@@@@!!"<<endl;
  }

  //go through the snd and rcv domain and exclude the bc 4444/4441/4442 (no point exchange for the periodic domain boundary
  ZFSInt noNghbrDomains=0;
  ZFSInt noNghbrDomainsPeriodic=0;
  ZFSInt noNghbrDomainsPeriodicS=0;
  ZFSInt noNghbrDomainsChannel=0;
  ZFSInt noNghbrDomainsSingular=0;
  ZFSInt noNghbrDomainsMultiBlock=0;

  for(ZFSInt i=0; i<(ZFSInt)rcvMap.size(); i++) {
    switch(rcvMap[i]->BC)
    {
    case 6000: // normal partition communication
    case 6333: // singular communications
    case 6666: // additional partition communication
    {
      noNghbrDomains++;
      if(rcvMap[i]->BC==6333) noNghbrDomainsSingular++;
      if(rcvMap[i]->BC==6000) noNghbrDomainsMultiBlock++;

      break;
    }
    case 4011:
    case 4012:
    case 4401:
    case 4402:
    case 4403:
    case 4404:
    case 4405:
    case 4406:
    {
      noNghbrDomainsPeriodic++;
      if(rcvMap[i]->Nstar!=-1) noNghbrDomainsPeriodicS++;
      break;
    }
    case 6001:
    {
      noNghbrDomainsChannel++;
      break;
    }
    default:
    {
      //do nothing;
    }
    }
  }
  comm->noNghbrDomainsNormal=noNghbrDomains;
  comm->noNghbrDomainsChannel=noNghbrDomainsChannel;
  comm->noNghbrDomainsPeriodic=noNghbrDomainsPeriodic;
  comm->noNghbrDomainsPeriodicS=noNghbrDomainsPeriodicS;
  comm->noNghbrDomainsSingular=noNghbrDomainsSingular;
  comm->noNghbrDomainsMultiBlock=noNghbrDomainsMultiBlock;

  // cout << "rank " << domainId()
  //      << " noDomains6000: " << noNghbrDomains
  //      << " noDomainsnormal: " << noNghbrDomainsMultiBlock
  //      << " noDomainsingular: " <<  noNghbrDomainsSingular
  //      << " noNghbrDomainsPeriodicTotal: "<< noNghbrDomainsPeriodic
  //      << " noNghbrDomainsPeriodicS: "<< noNghbrDomainsPeriodicS
  //      << " all noNghbrDomains: " << sndMap.size() << endl;

  //from the snd- and rcvMap it should be now known where communication takes place and especially the number of cells which are needed
  //==> create a communication handle containing all the infos. ==> copy info from map to communication handle

  comm->noNghbrDomains=sndMap.size(); //number of neighbour Domains

  comm->m_sndNghbrId=new ZFSInt[comm->noNghbrDomains]; //rank of snd neighbour domain
  comm->m_rcvNghbrId= new ZFSInt[comm->noNghbrDomains]; //rank of rcv neighbour domain
  comm->m_noNghbrDomainCellBufferSizeSnd = new ZFSInt[comm->noNghbrDomains]; //cell size Snd
  comm->m_noNghbrDomainPointBufferSizeSnd = new ZFSInt[comm->noNghbrDomains]; // point size Snd
  comm->m_noNghbrDomainCellBufferSizeRcv = new ZFSInt[comm->noNghbrDomains]; //cell size Rcv
  comm->m_noNghbrDomainPointBufferSizeRcv = new ZFSInt[comm->noNghbrDomains]; // point size Rcv
  // comm->mpi_sndRequest = new MPI_Request[comm->noNghbrDomains]; //sending request
  // comm->mpi_rcvRequest = new MPI_Request[comm->noNghbrDomains]; //receiving request
  // comm->mpi_rcvStatus = new MPI_Status[comm->noNghbrDomains]; //receiving status
  // comm->mpi_sndStatus = new MPI_Status[comm->noNghbrDomains]; //sending status
  comm->startInfoSNDcells = new ZFSInt *[comm->noNghbrDomains]; //start cells for snd
  comm->endInfoSNDcells = new ZFSInt *[comm->noNghbrDomains]; //end cells for snd
  comm->startInfoRCVcells = new ZFSInt *[comm->noNghbrDomains]; //start cells for rcv
  comm->endInfoRCVcells = new ZFSInt *[comm->noNghbrDomains]; //end cells for rcv
  comm->startInfoSNDpoints = new ZFSInt *[comm->noNghbrDomains]; //start points for snd
  comm->endInfoSNDpoints = new ZFSInt *[comm->noNghbrDomains]; //end points for snd
  comm->startInfoRCVpoints = new ZFSInt *[comm->noNghbrDomains]; //start points for rcv
  comm->endInfoRCVpoints = new ZFSInt *[comm->noNghbrDomains]; //end points for rcv
  comm->orderInfo= new ZFSInt *[comm->noNghbrDomains]; // information of the order i,j,k
  comm->stepInfoRCV= new ZFSInt *[comm->noNghbrDomains]; // information of the steps
  comm->bcId= new ZFSInt [comm->noNghbrDomains]; // information of the comm bc

  //singularity
  comm->singularID= new ZFSInt[comm->noNghbrDomainsSingular];

  for(ZFSInt i=0; i<comm->noNghbrDomains; i++) {
    comm->startInfoSNDcells[i]= new ZFSInt [nDim];
    comm->endInfoSNDcells[i]=new ZFSInt[nDim];
    comm->startInfoRCVcells[i]= new ZFSInt [nDim];
    comm->endInfoRCVcells[i]=new ZFSInt[nDim];
    comm->startInfoSNDpoints[i]= new ZFSInt [nDim];
    comm->endInfoSNDpoints[i]=new ZFSInt[nDim];
    comm->startInfoRCVpoints[i]= new ZFSInt [nDim];
    comm->endInfoRCVpoints[i]=new ZFSInt[nDim];
    comm->orderInfo[i]= new ZFSInt[nDim];
    comm->stepInfoRCV[i]= new ZFSInt[nDim];
  }

  comm->m_tagHelperSND=NULL;
  comm->m_tagHelperRCV=NULL;

  if(comm->noNghbrDomains!=0) {
    ZFSInt tempId;
    comm->m_tagHelperSND=new ZFSInt[comm->noNghbrDomains];
    comm->m_tagHelperRCV=new ZFSInt[comm->noNghbrDomains];
    //go through the communication domains
    for(ZFSInt i=0; i<comm->noNghbrDomains; i++) {
      tempId=rcvMap[i]->Id2;
      //rank != cpuId ==> find the corresponding rank
      for(ZFSInt j=0; j<m_noPartitions; j++) {
        if(tempId==partition->outputBoxInfo[j]->cpu) {
          comm->m_rcvNghbrId[i]=j;
          break;
        }
      }

      tempId=sndMap[i]->Id2;
      //rank != cpuId ==> find the corresponding rank
      for(ZFSInt j=0; j<m_noPartitions; j++) {
        if(tempId==partition->outputBoxInfo[j]->cpu) {
          comm->m_sndNghbrId[i]=j;
          break;
        }
      }

      memcpy(comm->startInfoSNDcells[i], sndMap[i]->start1, nDim*sizeof(ZFSInt));
      memcpy(comm->endInfoSNDcells[i], sndMap[i]->end1, nDim*sizeof(ZFSInt));
      memcpy(comm->startInfoRCVcells[i], rcvMap[i]->start1, nDim*sizeof(ZFSInt));
      memcpy(comm->endInfoRCVcells[i], rcvMap[i]->end1, nDim*sizeof(ZFSInt));
      memcpy(comm->startInfoSNDpoints[i], sndMap[i]->start1, nDim*sizeof(ZFSInt));
      memcpy(comm->endInfoSNDpoints[i], sndMap[i]->end1, nDim*sizeof(ZFSInt));
      memcpy(comm->startInfoRCVpoints[i], rcvMap[i]->start1, nDim*sizeof(ZFSInt));
      memcpy(comm->endInfoRCVpoints[i], rcvMap[i]->end1, nDim*sizeof(ZFSInt));
      memcpy(comm->orderInfo[i], rcvMap[i]->order, nDim*sizeof(ZFSInt));
      memcpy(comm->stepInfoRCV[i], rcvMap[i]->step2, nDim*sizeof(ZFSInt));

      //compute the buffersizes for sending and receiving
      ZFSInt cellproductSnd=1, cellproductRcv=1;
      ZFSInt pointproductSnd=1, pointproductRcv=1;
      for(ZFSInt j=0; j<nDim; j++) {
        ZFSInt cellsizes=(comm->endInfoSNDcells[i][j] - comm->startInfoSNDcells[i][j]);
        ZFSInt cellsizesRcv=(comm->endInfoRCVcells[i][j] - comm->startInfoRCVcells[i][j]);
        if(cellsizes!=0&&cellsizesRcv!=0) {
          //cout << "length is given as " << cellsizes << endl;
          cellproductSnd*=(cellsizes);
          pointproductSnd*=(cellsizes+1);
          cellproductRcv*=(cellsizesRcv);
          pointproductRcv*=(cellsizesRcv+1);
        }
      }

      // if(cellproductSnd!=1) comm->m_noNghbrDomainCellBufferSizeSnd[i]=cellproductSnd*noVariables;
      // if(cellproductRcv!=1)  comm->m_noNghbrDomainCellBufferSizeRcv[i]=cellproductRcv*noVariables;
      // if(pointproductSnd!=1) comm->m_noNghbrDomainPointBufferSizeSnd[i]=pointproductSnd*nDim;
      // if(pointproductRcv!=1) comm->m_noNghbrDomainPointBufferSizeRcv[i]=pointproductRcv*nDim;

      comm->m_noNghbrDomainCellBufferSizeSnd[i]=cellproductSnd*noVariables;
      comm->m_noNghbrDomainCellBufferSizeRcv[i]=cellproductRcv*noVariables;
      comm->m_noNghbrDomainPointBufferSizeSnd[i]=pointproductSnd*nDim;
      comm->m_noNghbrDomainPointBufferSizeRcv[i]=pointproductRcv*nDim;

      comm->m_tagHelperRCV[i]=rcvMap[i]->face;
      if(comm->m_tagHelperRCV[i]<0) {
        comm->m_tagHelperRCV[i]=0;
      }
      comm->m_tagHelperSND[i]=sndMap[i]->face;
      if(comm->m_tagHelperSND[i]<0) {
        comm->m_tagHelperSND[i]=0;
      }


      comm->bcId[i] = sndMap[i]->BC;
    }
  }

  // for(ZFSId i=0; i<comm->noNghbrDomains; i++) {
  //   if(comm->bcId[i]==4403&&sndMap[i]->Nstar==-1)//if(comm->bcId[i]==6000)
  //     cout << "rank " << domainId() << " send to rank " << comm->m_nghbrDomainCPU[i]<<" i:"<< comm->startInfoSNDcells[i][0] <<"-"<<comm->endInfoSNDcells[i][0] <<" j:"<< comm->startInfoSNDcells[i][1] <<"-"<<comm->endInfoSNDcells[i][1] <<" k:"<< comm->startInfoSNDcells[i][2] <<"-"<<comm->endInfoSNDcells[i][2] <<" tag:"<<comm->m_tagHelperSND[i]<< " BC:"<<comm->bcId[i]<< endl;

  // }
  // for(ZFSId i=0; i<comm->noNghbrDomains; i++) {
  //   if(comm->bcId[i]==4403&&rcvMap[i]->Nstar==-1)
  //     cout << "rank " << domainId() << " recv from rank " << comm->m_nghbrDomainRankId[i]<< " i:"<< comm->startInfoRCVcells[i][0] <<"-"<<comm->endInfoRCVcells[i][0] <<" j:"<< comm->startInfoRCVcells[i][1] <<"-"<<comm->endInfoRCVcells[i][1] <<" k:"<< comm->startInfoRCVcells[i][2] <<"-"<<comm->endInfoRCVcells[i][2]<< " order: "<< comm->orderInfo[i][0] <<" "<<comm->orderInfo[i][1] <<" "<< comm->orderInfo[i][2] <<" step:"<<comm->stepInfoRCV[i][0] <<" "<< comm->stepInfoRCV[i][1] <<" "<<comm->stepInfoRCV[i][2] <<" tag:"<<comm->m_tagHelperRCV[i]<< " BC:"<<comm->bcId[i]<< endl;
  // }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::createWindowInfo()
{
  //noBoxWindows*
  outputBoxWindowInformation = new noBoxWindows;
  int par_size = noDomains();
  int rank= domainId();

  /*this function creates the window_information for each partition (actually only the one of mpi_rank), making use of the information conatined in the inputBoxInfo and the information about the partition properties stored in outputBoxInfo

    coordinates ==> indices
    //    face 1/2: const i: j, j, k, k ==> 0, 1, 1, 2, 2
    //    face 3/4: const j: k, k, i, i ==> 1, 2, 2, 0, 0
    //    face 5/6: const k: i, i, j, j ==> 2, 0, 0, 1, 1

    -----------------------------
    |\                          |\
    | \                         | \
    |  \                        |  \
    |   \             (4)       |   \
    |    \                      |    \
    |     \      (6)            |     \
    |      \ __________________________\
    |       |                   |      |
    |  (1)  |                   |  (2) |
    |       |                   |      |
    --------|--------------------      |
    \       |          (5)      \      |
     \      |                    \     |
    j \     |                     \    |
    |  \    |         (3)          \   |
    |__i\   |                       \  |
    \    \  |                        \ |
     \    \ |                         \|
      k     ----------------------------

  */


  int noFace=2*nDim;
  for(int face=0;face<noFace;face++)
  {
    int inputID = partition->outputBoxInfo[rank]->inputBoxID;
    vector<windowInformation*> side;
    outputBoxWindowInformation->sides.push_back(side);
    //this is for indeces starting with one
    //int fixed=(face+2)/2;
    //int dir1=(fixed%nDim)+1;
    //int dir2=(dir1%nDim)+1;
    //==> need to subtract one unit (C++ starts with 0)!!
    // fixed--;
    //dir1--;
    //dir2--;
    int fixed=-1, dir1=-1, dir2=-1;
    switch(face)
    {
    case 0:
    {
      fixed=2;
      dir1=1;
      dir2=0;
      break;
    }
    case 1:
    {
      fixed=2;
      dir1=1;
      dir2=0;
      break;
    }
    case 2:
    {
      fixed=1;
      dir1=2;
      dir2=0;
      break;
    }
    case 3:
    {
      fixed=1;
      dir1=2;
      dir2=0;
      break;
    }
    case 4:
    {
      fixed=0;
      dir1=2;
      dir2=1;
      break;
    }
    case 5:
    {
      fixed=0;
      dir1=2;
      dir2=1;
      break;
    }
    default:
    {
      zfsTerm(1,__CALLING_FUNCTION__, "No matching case for face id");
    }
    }

    /*
      first check if current side is part of the inputBox->WindowBoundary!
      If yes attribute the Boundary condition as given for the file else search for corresponding pair.

      if modulo 2 is odd then the side starts with index zero (see TFS side convention)
    */
    int inputBoxFixDir, outputBoxFixDir;
    if(((face+1)%2)==1)//on odd sides
    {
      inputBoxFixDir=0;
      outputBoxFixDir=0;
    }
    else
    {
      inputBoxFixDir=partition->inputBoxInfo[inputID]->DirLast[fixed];
      outputBoxFixDir=partition->outputBoxInfo[rank]->DirLast[fixed];

    }


    int BoxStartDim1=partition->outputBoxInfo[rank]->offset[dir1]+1;
    int BoxStartDim2=partition->outputBoxInfo[rank]->offset[dir2]+1;
    int BoxEndDim1=partition->outputBoxInfo[rank]->offset[dir1]+partition->outputBoxInfo[rank]->DirLast[dir1];
    int BoxEndDim2=partition->outputBoxInfo[rank]->offset[dir2]+partition->outputBoxInfo[rank]->DirLast[dir2];
    int BoxL1=(BoxEndDim1-BoxStartDim1)+1;
    //int BoxL2=(BoxEndDim2-BoxStartDim2)+1;
    //int Box1off=BoxStartDim1;
    //int Box2off=BoxStartDim2;

    if ((outputBoxFixDir + partition->outputBoxInfo[rank]->offset[fixed]) == inputBoxFixDir)
    {
      //At the edge of the inputBox!
      //RESTRICTION:: for the moment the partitioned box can only be part of one original block!!
      //This needs to be changed later!!!!!


      int noInputFaceSegments= inputBoxWindowInformation[inputID].sides[face].size();
      printf("no of Segments %d on side %d\n", noInputFaceSegments, (face+1));

      //search through all input segments of the side if one does correspod to the boxside

      for(int i=0; i<noInputFaceSegments; i++)
      {
        //is segment part of the box side???
        int SegStartDim1=inputBoxWindowInformation[inputID].sides[face][i]->startindex[0];
        int SegStartDim2=inputBoxWindowInformation[inputID].sides[face][i]->startindex[1];
        int SegEndDim1=inputBoxWindowInformation[inputID].sides[face][i]->endindex[0];
        int SegEndDim2=inputBoxWindowInformation[inputID].sides[face][i]->endindex[1];
        int BoundaryCond = inputBoxWindowInformation[inputID].sides[face][i]->BC;
        int SegL1=(SegEndDim1-SegStartDim1)+1;
        //int SegL2=(SegEndDim2-SegStartDim2)+1;
        //int Seg1off=SegStartDim1;
        //int Seg2off=SegStartDim2;
        printf("\n---------\n");
        printf("side %d\n",face);
        printf("BoxStartDim1 %d\n", BoxStartDim1);
        printf("BoxStartDim2 %d\n", BoxStartDim2);
        printf("BoxEndDim1 %d\n", BoxEndDim1);
        printf("BoxEndDim2 %d\n\n", BoxEndDim2);
        printf("startDim1 %d\n", SegStartDim1);
        printf("startDim2 %d\n", SegStartDim2);
        printf("endDim1 %d\n", SegEndDim1);
        printf("endDim2 %d\n", SegEndDim2);
        //do the search if a part of the segment is contained in the box with the modulo function
        //check all cornerpoints: is one containded?
        if((BoxStartDim1>=SegStartDim1)&&(BoxStartDim1<=SegEndDim1)&&(BoxStartDim2>=SegStartDim2)&&(BoxStartDim2<=SegEndDim2))
        {
          //bottom left corner
          windowInformation* segment =new windowInformation(nDim);
          segment->inDir1=dir1;
          segment->inDir2=dir2;
          segment->fixDir=fixed;
          segment->startindex[0]=BoxStartDim1;
          segment->startindex[1]=BoxStartDim2;
          //search in both direction how far one can go until boarder is reached
          //start with direction one
          if(BoxEndDim1>=SegEndDim1)
          {
            segment->endindex[0]=SegEndDim1;
          }
          else
          {
            segment->endindex[0]=BoxEndDim1;
          }
          //search in direction two
          if(BoxEndDim2>=SegEndDim2)
          {
            segment->endindex[1]=SegEndDim2;
          }
          else
          {
            segment->endindex[1]=BoxEndDim2;
          }
          outputBoxWindowInformation->sides[face].push_back(segment);
          segment->BC=BoundaryCond;
        }
        else if((BoxEndDim1>=SegStartDim1)&&(BoxEndDim1<=SegEndDim1)&&(BoxStartDim2>=SegStartDim2)&&(BoxStartDim2<=SegEndDim2))
        {
          //bottom right corner
          windowInformation* segment =new windowInformation(nDim);
          segment->inDir1=dir1;
          segment->inDir2=dir2;
          segment->fixDir=fixed;
          segment->endindex[0]=BoxEndDim1;
          segment->startindex[1]=BoxStartDim2;
          //search in both direction how far one can go until boarder is reached
          //start with direction one
          if(BoxStartDim1>=SegStartDim1)
          {
            segment->startindex[0]=BoxStartDim1;
          }
          else
          {
            segment->startindex[0]=SegStartDim1;
          }
          //search in direction two
          if(BoxEndDim2>=SegEndDim2)
          {
            segment->endindex[1]=SegEndDim2;
          }
          else
          {
            segment->endindex[1]=BoxEndDim2;
          }
          segment->BC=BoundaryCond;
          outputBoxWindowInformation->sides[face].push_back(segment);
        }
        else if((BoxStartDim1>=SegStartDim1)&&(BoxStartDim1<=SegEndDim1)&&(BoxEndDim2>=SegStartDim2)&&(BoxEndDim2<=SegEndDim2))
        {
          //top left corner
          windowInformation* segment =new windowInformation(nDim);
          segment->inDir1=dir1;
          segment->inDir2=dir2;
          segment->fixDir=fixed;
          segment->startindex[0]=BoxStartDim1;
          segment->endindex[1]=BoxEndDim2;
          //search in both direction how far one can go until boarder is reached
          //start with direction one
          if(BoxEndDim1>=SegEndDim1)
          {
            segment->endindex[0]=SegEndDim1;
          }
          else
          {
            segment->endindex[0]=BoxEndDim1;
          }
          //search in direction two
          if(BoxStartDim2>=SegStartDim2)
          {
            segment->startindex[1]=BoxStartDim2;
          }
          else
          {
            segment->endindex[1]=SegStartDim2;
          }
          segment->BC=BoundaryCond;
          outputBoxWindowInformation->sides[face].push_back(segment);
        }
        else if((BoxEndDim1>=SegStartDim1)&&(BoxEndDim1<=SegEndDim1)&&(BoxEndDim2>=SegStartDim2)&&(BoxEndDim2<=SegEndDim2))
        {
          //top right corner
          windowInformation* segment =new windowInformation(nDim);
          segment->inDir1=dir1;
          segment->inDir2=dir2;
          segment->fixDir=fixed;
          segment->endindex[0]=BoxEndDim1;
          segment->endindex[1]=BoxEndDim2;
          //search in both direction how far one can go until boarder is reached
          //start with direction one
          if(BoxStartDim1>=SegStartDim1)
          {
            segment->startindex[0]=SegEndDim1;
          }
          else
          {
            segment->startindex[0]=BoxStartDim1;
          }
          //search in direction two
          if(BoxStartDim2>=SegStartDim2)
          {
            segment->startindex[1]=BoxStartDim2;
          }
          else
          {
            segment->endindex[1]=SegStartDim2;
          }
          segment->BC=BoundaryCond;
          outputBoxWindowInformation->sides[face].push_back(segment);
        }
        else if((BoxStartDim1<SegStartDim1)&&(BoxEndDim1>SegEndDim1)&&(SegL1<BoxL1))
        {
          //segment is completely contained in output box
          windowInformation* segment =new windowInformation(nDim);
          segment->inDir1=dir1;
          segment->inDir2=dir2;
          segment->fixDir=fixed;
          segment->startindex[0]=SegStartDim1;
          segment->endindex[0]=SegEndDim1;
          //check second direction only
          if(BoxStartDim2>=SegStartDim2)
          {
            segment->startindex[1]=BoxStartDim2;
          }
          else
          {
            segment->startindex[1]=SegStartDim2;
          }
          if(BoxEndDim2>=SegEndDim2)
          {
            segment->endindex[1]=SegEndDim2;
          }
          else
          {
            segment->endindex[1]=BoxEndDim2;
          }
          segment->BC=BoundaryCond;
          outputBoxWindowInformation->sides[face].push_back(segment);
        }
      }
    }
    else
    {
      //inside the input Box

      int dir;

      if(((face+1)%2)==1)//on odd sides
      {
        //on odd sides
        dir= partition->outputBoxInfo[rank]->offset[fixed];
      }
      else
      {
        dir= partition->outputBoxInfo[rank]->offset[fixed]+partition->outputBoxInfo[rank]->DirLast[fixed];
      }

      //check with even sides of all other blocks if side is in common
      for(int j=0; j<par_size; j++)
      {
        //go through all domains an search if the indices are common
        if(((partition->outputBoxInfo[j]->offset[fixed]+partition->outputBoxInfo[j]->DirLast[fixed])==dir||partition->outputBoxInfo[j]->offset[fixed]==dir)&&j!=rank)
        {
          //found a pair==>same procedure as above! find the common points
          //start with direction one and end with the second one.
          //search for every corner point and the possibility of being
          //completely contained!!!
          int Box1StartDim1=partition->outputBoxInfo[j]->offset[dir1]+1;
          int Box1StartDim2=partition->outputBoxInfo[j]->offset[dir2]+1;
          int Box1EndDim1=partition->outputBoxInfo[j]->offset[dir1]+partition->outputBoxInfo[j]->DirLast[dir1];
          int Box1EndDim2=partition->outputBoxInfo[j]->offset[dir2]+partition->outputBoxInfo[j]->DirLast[dir2];
          int Box1L1=(Box1EndDim1-Box1StartDim1)+1;
          //same as before
          if((Box1StartDim1>=BoxStartDim1)&&(Box1StartDim1<=BoxEndDim1)&&(Box1StartDim2>=BoxStartDim2)&&(Box1StartDim2<=BoxEndDim2))
          {
            //bottom left corner
            windowInformation* segment =new windowInformation(nDim);
            segment->inDir1=dir1;
            segment->inDir2=dir2;
            segment->fixDir=fixed;
            segment->startindex[0]=Box1StartDim1;
            segment->startindex[1]=Box1StartDim2;
            //search in both direction how far one can go until boarder is reached
            //start with direction one
            if(Box1EndDim1>=BoxEndDim1)
            {
              segment->endindex[0]=BoxEndDim1;
            }
            else
            {
              segment->endindex[0]=Box1EndDim1;
            }
            //search in direction two
            if(Box1EndDim2>=BoxEndDim2)
            {
              segment->endindex[1]=BoxEndDim2;
            }
            else
            {
              segment->endindex[1]=Box1EndDim2;
            }
            outputBoxWindowInformation->sides[face].push_back(segment);
            segment->BC=0;
            segment->cpu=j;
          }
          else if((Box1EndDim1>=BoxStartDim1)&&(Box1EndDim1<=BoxEndDim1)&&(Box1StartDim2>=BoxStartDim2)&&(Box1StartDim2<=BoxEndDim2))
          {
            //bottom right corner
            windowInformation* segment =new windowInformation(nDim);
            segment->inDir1=dir1;
            segment->inDir2=dir2;
            segment->fixDir=fixed;
            segment->endindex[0]=Box1EndDim1;
            segment->startindex[1]=Box1StartDim2;
            //search in both direction how far one can go until boarder is reached
            //start with direction one
            if(Box1StartDim1>=BoxStartDim1)
            {
              segment->startindex[0]=Box1StartDim1;
            }
            else
            {
              segment->startindex[0]=BoxStartDim1;
            }
            //search in direction two
            if(Box1EndDim2>=BoxEndDim2)
            {
              segment->endindex[1]=BoxEndDim2;
            }
            else
            {
              segment->endindex[1]=Box1EndDim2;
            }
            segment->BC=0;
            segment->cpu=j;
            outputBoxWindowInformation->sides[face].push_back(segment);
          }
          else if((Box1StartDim1>=BoxStartDim1)&&(Box1StartDim1<=BoxEndDim1)&&(Box1EndDim2>=BoxStartDim2)&&(Box1EndDim2<=BoxEndDim2))
          {
            //top left corner
            windowInformation* segment =new windowInformation(nDim);
            segment->inDir1=dir1;
            segment->inDir2=dir2;
            segment->fixDir=fixed;
            segment->startindex[0]=Box1StartDim1;
            segment->endindex[1]=Box1EndDim2;
            //search in both direction how far one can go until boarder is reached
            //start with direction one
            if(Box1EndDim1>=BoxEndDim1)
            {
              segment->endindex[0]=BoxEndDim1;
            }
            else
            {
              segment->endindex[0]=Box1EndDim1;
            }
            //search in direction two
            if(Box1StartDim2>=BoxStartDim2)
            {
              segment->startindex[1]=Box1StartDim2;
            }
            else
            {
              segment->endindex[1]=BoxStartDim2;
            }
            segment->BC=0;
            segment->cpu=j;
            outputBoxWindowInformation->sides[face].push_back(segment);
          }
          else if((Box1EndDim1>=BoxStartDim1)&&(Box1EndDim1<=BoxEndDim1)&&(Box1EndDim2>=BoxStartDim2)&&(Box1EndDim2<=BoxEndDim2))
          {
            //top right corner
            windowInformation* segment =new windowInformation(nDim);
            segment->inDir1=dir1;
            segment->inDir2=dir2;
            segment->fixDir=fixed;
            segment->endindex[0]=Box1EndDim1;
            segment->endindex[1]=Box1EndDim2;
            //search in both direction how far one can go until boarder is reached
            //start with direction one
            if(Box1StartDim1>=BoxStartDim1)
            {
              segment->startindex[0]=BoxEndDim1;
            }
            else
            {
              segment->startindex[0]=Box1StartDim1;
            }
            //search in direction two
            if(Box1StartDim2>=BoxStartDim2)
            {
              segment->startindex[1]=Box1StartDim2;
            }
            else
            {
              segment->endindex[1]=BoxStartDim2;
            }
            segment->BC=0;
            segment->cpu=j;
            outputBoxWindowInformation->sides[face].push_back(segment);
          }
          else if((Box1StartDim1<BoxStartDim1)&&(Box1EndDim1>BoxEndDim1)&&(Box1L1<BoxL1))
          {
            //segment is completely contained in output box
            windowInformation* segment =new windowInformation(nDim);
            segment->inDir1=dir1;
            segment->inDir2=dir2;
            segment->fixDir=fixed;
            segment->startindex[0]=BoxStartDim1;
            segment->endindex[0]=BoxEndDim1;
            //check second direction only
            if(Box1StartDim2>=BoxStartDim2)
            {
              segment->startindex[1]=Box1StartDim2;
            }
            else
            {
              segment->startindex[1]=BoxStartDim2;
            }
            if(Box1EndDim2>=BoxEndDim2)
            {
              segment->endindex[1]=BoxEndDim2;
            }
            else
            {
              segment->endindex[1]=Box1EndDim2;
            }
            segment->BC=0;
            segment->cpu=j;
            outputBoxWindowInformation->sides[face].push_back(segment);
          }
        }
      }
    }
  }

  //just a test: write out all the window information data to see if correct
  if(domainId()==1)
  {
    for(int face=0;face<noFace;face++)
    {
      printf("outputboxside %i contains %i segments, <==with rank %i\n", (face+1),(ZFSInt)outputBoxWindowInformation->sides[face].size(), rank);
      for(int i=0; i<(int)outputBoxWindowInformation->sides[face].size(); i++)
      {
        printf("startDim1 %d \n", outputBoxWindowInformation->sides[face][i]->startindex[0]);
        printf("startDim2 %d \n", outputBoxWindowInformation->sides[face][i]->startindex[1]);
        printf("endDim1 %d \n", outputBoxWindowInformation->sides[face][i]->endindex[0]);
        printf("endDim2 %d \n", outputBoxWindowInformation->sides[face][i]->endindex[1]);
        printf("Boundary Condition ==> %d\n", outputBoxWindowInformation->sides[face][i]->BC);
        printf("corresponding cpu==> %d\n", outputBoxWindowInformation->sides[face][i]->cpu);

      }
    }
  }

  //=======NOTICE==========

  //The start and end values are given in absolute coordinates and not in local ones
  //therefore need to subtract the local box infos to get everything in good order

  for(int face=0;face<noFace;face++)
  {
    printf("outputboxside %i contains %i segments, <==with rank %i\n", (face+1), (ZFSInt)outputBoxWindowInformation->sides[face].size(), rank);
    for(int i=0; i<(int) outputBoxWindowInformation->sides[face].size(); i++)
    {
      switch(face)
      {
      case 0:
      {
        outputBoxWindowInformation->sides[face][i]->startindex[0]=outputBoxWindowInformation->sides[face][i]->startindex[0]-partition->outputBoxInfo[rank]->offset[1];
        outputBoxWindowInformation->sides[face][i]->startindex[1]=outputBoxWindowInformation->sides[face][i]->startindex[1]-partition->outputBoxInfo[rank]->offset[0];
        outputBoxWindowInformation->sides[face][i]->endindex[0]=outputBoxWindowInformation->sides[face][i]->endindex[0]-partition->outputBoxInfo[rank]->offset[1];
        outputBoxWindowInformation->sides[face][i]->endindex[1]=outputBoxWindowInformation->sides[face][i]->endindex[1]-partition->outputBoxInfo[rank]->offset[0];
        break;
      }
      case 1:
      {
        outputBoxWindowInformation->sides[face][i]->startindex[0]=outputBoxWindowInformation->sides[face][i]->startindex[0]-partition->outputBoxInfo[rank]->offset[1];
        outputBoxWindowInformation->sides[face][i]->startindex[1]=outputBoxWindowInformation->sides[face][i]->startindex[1]-partition->outputBoxInfo[rank]->offset[0];
        outputBoxWindowInformation->sides[face][i]->endindex[0]=outputBoxWindowInformation->sides[face][i]->endindex[0]-partition->outputBoxInfo[rank]->offset[1];
        outputBoxWindowInformation->sides[face][i]->endindex[1]=outputBoxWindowInformation->sides[face][i]->endindex[1]-partition->outputBoxInfo[rank]->offset[0];
        break;
      }
      case 2:
      {
        outputBoxWindowInformation->sides[face][i]->startindex[0]=outputBoxWindowInformation->sides[face][i]->startindex[0]-partition->outputBoxInfo[rank]->offset[2];
        outputBoxWindowInformation->sides[face][i]->startindex[1]=outputBoxWindowInformation->sides[face][i]->startindex[1]-partition->outputBoxInfo[rank]->offset[0];
        outputBoxWindowInformation->sides[face][i]->endindex[0]=outputBoxWindowInformation->sides[face][i]->endindex[0]-partition->outputBoxInfo[rank]->offset[2];
        outputBoxWindowInformation->sides[face][i]->endindex[1]=outputBoxWindowInformation->sides[face][i]->endindex[1]-partition->outputBoxInfo[rank]->offset[0];
        break;
      }
      case 3:
      {
        outputBoxWindowInformation->sides[face][i]->startindex[0]=outputBoxWindowInformation->sides[face][i]->startindex[0]-partition->outputBoxInfo[rank]->offset[2];
        outputBoxWindowInformation->sides[face][i]->startindex[1]=outputBoxWindowInformation->sides[face][i]->startindex[1]-partition->outputBoxInfo[rank]->offset[0];
        outputBoxWindowInformation->sides[face][i]->endindex[0]=outputBoxWindowInformation->sides[face][i]->endindex[0]-partition->outputBoxInfo[rank]->offset[2];
        outputBoxWindowInformation->sides[face][i]->endindex[1]=outputBoxWindowInformation->sides[face][i]->endindex[1]-partition->outputBoxInfo[rank]->offset[0];
        break;
      }
      case 4:
      {
        outputBoxWindowInformation->sides[face][i]->startindex[0]=outputBoxWindowInformation->sides[face][i]->startindex[0]-partition->outputBoxInfo[rank]->offset[2];
        outputBoxWindowInformation->sides[face][i]->startindex[1]=outputBoxWindowInformation->sides[face][i]->startindex[1]-partition->outputBoxInfo[rank]->offset[1];
        outputBoxWindowInformation->sides[face][i]->endindex[0]=outputBoxWindowInformation->sides[face][i]->endindex[0]-partition->outputBoxInfo[rank]->offset[2];
        outputBoxWindowInformation->sides[face][i]->endindex[1]=outputBoxWindowInformation->sides[face][i]->endindex[1]-partition->outputBoxInfo[rank]->offset[1];
        break;
      }
      case 5:
      {
        outputBoxWindowInformation->sides[face][i]->startindex[0]=outputBoxWindowInformation->sides[face][i]->startindex[0]-partition->outputBoxInfo[rank]->offset[2];
        outputBoxWindowInformation->sides[face][i]->startindex[1]=outputBoxWindowInformation->sides[face][i]->startindex[1]-partition->outputBoxInfo[rank]->offset[1];
        outputBoxWindowInformation->sides[face][i]->endindex[0]=outputBoxWindowInformation->sides[face][i]->endindex[0]-partition->outputBoxInfo[rank]->offset[2];
        outputBoxWindowInformation->sides[face][i]->endindex[1]=outputBoxWindowInformation->sides[face][i]->endindex[1]-partition->outputBoxInfo[rank]->offset[1];
        break;
      }
      default:
      {
        zfsTerm(1,__CALLING_FUNCTION__, "No matching case for face id");
      }
      }
    }
  }
}

template <ZFSInt nDim>
int ZFSStrctrdBlckWindowInfo<nDim>::mapCompare(ZFSStrctrdWindowMap* map1, ZFSStrctrdWindowMap* map2)
{
  if(map1->Id1!=map2->Id1) return 0;
  if(map1->Id2!=map2->Id2) return 0;
  for(ZFSId i=0; i<nDim; i++)
    {
      if(map1->start1[i]!=map2->start1[i]) return 0;
      if(map1->end1[i]!=map2->end1[i]) return 0;
      //if(map1->start2[i]!=map2->start2[i]) return 0;
      //if(map1->end2[i]!=map2->end2[i]) return 0;
      if(map1->step1[i]!=map2->step1[i]) return 0;
      //if(map1->step2[i]!=map2->step2[i]) return 0;
      if(map1->order[i]!=map2->order[i]) return 0;
    }
  return 1;
}

template <ZFSInt nDim>
int ZFSStrctrdBlckWindowInfo<nDim>::mapCompare11(ZFSStrctrdWindowMap* map1, ZFSStrctrdWindowMap* map2)
{
  //only compare the Id1, start1 and end1
  if(map1->Id1!=map2->Id1) return 0;
  for(ZFSId i=0; i<nDim; i++)
    {
      if(map1->start1[i]!=map2->start1[i]) return 0;
      if(map1->end1[i]!=map2->end1[i]) return 0;
    }
  return 1;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCreate(ZFSId Id1, ZFSInt* start1, ZFSInt* end1, ZFSInt* step1, ZFSId Id2, ZFSInt* start2, ZFSInt* end2, ZFSInt* step2, ZFSInt* order, ZFSInt BC, ZFSStrctrdWindowMap* output)
{
  output->Id1=Id1;
  output->BC=BC;
  output->Nstar=-1;
  output->dc1=888;
  output->dc2=888;
  output->SingularId=-1;

  if(Id2!=-1)
    {
      output->Id2=Id2;
    }
  else
    {
      output->Id2=Id1;
    }

  memcpy(output->start1, start1, nDim*sizeof(ZFSInt));
  if(start2!=NULL)
    {
      memcpy(output->start2, start2, nDim*sizeof(ZFSInt));
    }
  else
    {
      memcpy(output->start2, start1, nDim*sizeof(ZFSInt));
    }

  memcpy(output->end1, end1, nDim*sizeof(ZFSInt));
  if(end2!=NULL)
    {
      memcpy(output->end2, end2, nDim*sizeof(ZFSInt));
    }
  else
    {
      memcpy(output->end2, end1, nDim*sizeof(ZFSInt));
    }

  if(step1!=NULL)
    {
      memcpy(output->step1, step1, nDim*sizeof(ZFSInt));
    }
  else
    {
      for(ZFSInt i=0; i<nDim; i++) output->step1[i]=1;
    }

  if(step2!=NULL)
    {
      memcpy(output->step2, step2, nDim*sizeof(ZFSInt));
    }
  else
    {
      for(ZFSInt i=0; i<nDim; i++) output->step2[i]=1;
    }

  if(order!=NULL)
    {
      memcpy(output->order, order, nDim*sizeof(ZFSInt));
    }
  else
    {
      for(ZFSInt i=1; i<=nDim; i++) output->order[i]=i;
    }

  output->SingularId=-1;
  output->Nstar=-1;
  output->dc1=-999;
  output->dc2=-999;

  /*
  bool check = mapCheck(output);
  if(!check)
    {
      mapPrint(output);
      zfsTerm(1, __CALLING_FUNCTION__, "Invalid Map");
    }
  */
}

template <ZFSInt nDim>
bool ZFSStrctrdBlckWindowInfo<nDim>::mapCheck(ZFSStrctrdWindowMap* input)
{
  // ZFSInt pointorline = nDim-1;
  //ZFSInt dummy=0;
  bool test;
  for(ZFSInt i =0; i<nDim; i++)
    {
      test = false;
      if(input->step1[i]==0) break;
      if(input->step2[i]==0) break;
      if(input->step1[i]>0)
        {
          if(input->end1[i]< input->start1[i]) break;
        }
      else
        {
          if(input->start1[i] < input->end1[i]) break;
        }
      if(input->step2[i]>0)
        {
          if(input->end2[i]< input->start2[i]) break;
        }
      else
        {
          if(input->start2[i] < input->end2[i]) break;
        }

      if(((input->end1[i]-input->start1[i])%(input->step1[i]))!=0) break;
      if(((input->end2[i]-input->start2[i])%(input->step2[i]))!=0) break;

      if(input->order[i] <0 || input->order[i]>(nDim-1)) break;

      if(((input->end1[i]-input->start1[i])/input->step1[i])!=((input->end2[input->order[i]]-input->start2[input->order[i]])/input->step2[input->order[i]])) break;

      //check if the map is a line (3D) or a point (2D-problem)
      //if it is a line then two of the three indicies are the same else for a point
      //only one is
      //if(input->start1[i]==input->end1[i]) dummy++;
       //if(dummy==pointorline) break;

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // line or point is acceptable!!!!

      test=true;

    }
  //check on Boundary Condition not possible
  return test;
}

template <ZFSInt nDim>
bool ZFSStrctrdBlckWindowInfo<nDim>::mapCheck0d(ZFSStrctrdWindowMap* input)
{
  // ZFSInt pointorline = 2;
  ZFSInt dummy=0;
  bool test;
  for(ZFSInt i =0; i<nDim; i++)
    {
      test = false;
      if(input->step1[i]==0) break;
      if(input->step2[i]==0) break;
      if(input->step1[i]>0)
        {
          if(input->end1[i]< input->start1[i]) break;
        }
      else
        {
          if(input->start1[i] < input->end1[i]) break;
        }
      if(input->step2[i]>0)
        {
          if(input->end2[i]< input->start2[i]) break;
        }
      else
        {
          if(input->start2[i] < input->end2[i]) break;
        }

      if(((input->end1[i]-input->start1[i])%(input->step1[i]))!=0) break;
      if(((input->end2[i]-input->start2[i])%(input->step2[i]))!=0) break;

      if(input->order[i] <0 || input->order[i]>(nDim-1)) break;

      if(((input->end1[i]-input->start1[i])/input->step1[i])!=((input->end2[input->order[i]]-input->start2[input->order[i]])/input->step2[input->order[i]])) break;

      //check if the map is a line (3D) or a point (2D-problem)
      //if it is a line then two of the three indicies are the same else for a point
      //only one is
      if(input->start1[i]==input->end1[i]) dummy++;


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // line or point is acceptable!!!!

      test=true;

    }
  if(dummy!=3) test=false;
  //check on Boundary Condition not possible
  return test;
}

template <ZFSInt nDim>
bool ZFSStrctrdBlckWindowInfo<nDim>::mapCheck1d(ZFSStrctrdWindowMap* input)
{
  // ZFSInt pointorline = 2;
  ZFSInt dummy=0;
  bool test;
  for(ZFSInt i =0; i<nDim; i++)
    {
      test = false;
      if(input->step1[i]==0) break;
      if(input->step2[i]==0) break;
      if(input->step1[i]>0)
        {
          if(input->end1[i]< input->start1[i]) break;
        }
      else
        {
          if(input->start1[i] < input->end1[i]) break;
        }
      if(input->step2[i]>0)
        {
          if(input->end2[i]< input->start2[i]) break;
        }
      else
        {
          if(input->start2[i] < input->end2[i]) break;
        }

      if(((input->end1[i]-input->start1[i])%(input->step1[i]))!=0) break;
      if(((input->end2[i]-input->start2[i])%(input->step2[i]))!=0) break;

      if(input->order[i] <0 || input->order[i]>(nDim-1)) break;

      if(((input->end1[i]-input->start1[i])/input->step1[i])!=((input->end2[input->order[i]]-input->start2[input->order[i]])/input->step2[input->order[i]])) break;

      //check if the map is a line (3D) or a point (2D-problem)
      //if it is a line then two of the three indicies are the same else for a point
      //only one is
      if(input->start1[i]==input->end1[i]) dummy++;


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // line or point is acceptable!!!!

      test=true;

    }
  if(dummy!=2) test=false;
  //check on Boundary Condition not possible
  return test;
}

template <ZFSInt nDim>
bool ZFSStrctrdBlckWindowInfo<nDim>::mapCheck2d(ZFSStrctrdWindowMap* input)
{
  //ZFSInt pointorline = nDim-1;
  ZFSInt dummy=0;
  ZFSBool result=false;
  for(ZFSInt i =0; i<nDim; i++)
    {
      result = false;
      if(input->step1[i]==0) break;
      if(input->step2[i]==0) break;
      if(input->step1[i]>0)
        {
          if(input->end1[i]< input->start1[i]) break;
        }
      else
        {
          if(input->start1[i] < input->end1[i]) break;
        }
      if(input->step2[i]>0)
        {
          if(input->end2[i]< input->start2[i]) break;
        }
      else
        {
          if(input->start2[i] < input->end2[i]) break;
        }

      if(((input->end1[i]-input->start1[i])%(input->step1[i]))!=0) break;
      if(((input->end2[i]-input->start2[i])%(input->step2[i]))!=0) break;

      if(input->order[i] <0 || input->order[i]>(nDim-1)) break;

      if(((input->end1[i]-input->start1[i])/input->step1[i])!=((input->end2[input->order[i]]-input->start2[input->order[i]])/input->step2[input->order[i]])) break;


      //check if the map is a line (3D) or a point (2D-problem)
      //if it is a line then two of the three indicies are the same else for a point
      //only one is
      if(input->start1[i]==input->end1[i]) dummy++;

      result = true;
    }

  if(dummy!=1||result==false) {
    return false;
  } else {
    return true;
  }
}

template <ZFSInt nDim>
bool ZFSStrctrdBlckWindowInfo<nDim>::mapCheck3d(ZFSStrctrdWindowMap* input)
{
  //ZFSInt pointorline = nDim-1;
  ZFSInt dummy=0;
  ZFSBool result = false;
  for(ZFSInt i =0; i<nDim; i++)
    {
      result = false;
      if(input->step1[i]==0) break;
      if(input->step2[i]==0) break;
      if(input->step1[i]>0)
        {
          if(input->end1[i]< input->start1[i]) break;
        }
      else
        {
          if(input->start1[i] < input->end1[i]) break;
        }
      if(input->step2[i]>0)
        {
          if(input->end2[i]< input->start2[i]) break;
        }
      else
        {
          if(input->start2[i] < input->end2[i]) break;
        }

      if(((input->end1[i]-input->start1[i])%(input->step1[i]))!=0) break;
      if(((input->end2[i]-input->start2[i])%(input->step2[i]))!=0) break;

      if(input->order[i] <0 || input->order[i]>(nDim-1)) break;

      if(((input->end1[i]-input->start1[i])/input->step1[i])!=((input->end2[input->order[i]]-input->start2[input->order[i]])/input->step2[input->order[i]])) break;

      //check if the map is a line (3D) or a point (2D-problem)
      //if it is a line then two of the three indicies are the same else for a point
      //only one is
      if(input->start1[i]==input->end1[i]) dummy++;

      result=true;
    }

  if(dummy!=0||result==false) {
    return false;
  } else {
    return true;
  }
}

template <ZFSInt nDim>
ZFSBool ZFSStrctrdBlckWindowInfo<nDim>::mapCheckWave(ZFSStrctrdWindowMap* input)
{
  //ZFSInt pointorline = nDim-1;
  ZFSInt dummy=0;
  ZFSBool result = false;
  for(ZFSInt i =0; i<nDim; i++)
    {
      result = false;
      if(input->step1[i]==0) break;
      if(input->step1[i]>0) {
        if(input->end1[i]< input->start1[i]) break;
      } else {
        if(input->start1[i] < input->end1[i]) break;
      }

      if(((input->end1[i]-input->start1[i])%(input->step1[i]))!=0) break;

      //check if the map is a line (3D) or a point (2D-problem)
      //if it is a line then two of the three indicies are the same else for a point
      //only one is
      if(input->start1[i]==input->end1[i]) dummy++;
      result=true;
    }

  if(dummy!=0||result==false) {
    return false;
  } else {
    return true;
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapPrint(ZFSStrctrdWindowMap* input)
{
  cout << "======== MAP INFO =======" << endl;
  cout << "Id1: " << input->Id1 << endl;
  stringstream start1;
  start1 << "start1: ";
  stringstream start2;
  start2 << "start2: ";
  stringstream end1;
  end1 << "end1: ";
  stringstream end2;
  end2 << "end2: ";
  stringstream step1;
  step1 << "step1: ";
  stringstream step2;
  step2 << "step2: ";
  stringstream order;
  order << "order: ";
  for(ZFSInt i=0; i<nDim; i++)
    {
      start1 << input->start1[i] <<" ";
      start2 << input->start2[i] <<" ";
      end1 << input->end1[i] << " ";
      end2 << input->end2[i] << " ";
      step1 << input->step1[i] << " ";
      step2 << input->step2[i] << " ";
      order << input->order[i] << " ";
    }
  cout << start1.str() << endl;
  cout << end1.str() << endl;
  cout << step1.str() << endl;
  cout << order.str() << endl;
  cout << "Id2: " << input->Id2 << endl;
  cout << start2.str() << endl;
  cout << end2.str() << endl;
  cout << step2.str() << endl;
  cout << "BC: " << input->BC << endl;
  cout << "Face: " << input->face << endl;

  cout << "spongInfos:" << endl;
  cout << "hasSponge: " << input->hasSponge << endl;
  cout << "spongeThickness: " << input->spongeThickness << endl;
  cout << "beta: "<< input->beta << endl;
  cout << "sigma: " << input->sigma << endl;

  cout<<"Nstar: "<<input->Nstar<<endl;
  cout<<"DC1: "<<input->dc1<<"  DC2: "<<input->dc2<<endl;
  cout << "======== MAP INFO END=======" << endl;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapPrintSimple(ZFSStrctrdWindowMap* input)
{
  stringstream start1;
  start1 << "start1: ";
  stringstream start2;
  start2 << "start2: ";

  for(ZFSInt i=0; i<nDim; i++)
    {
      start1 << input->start1[i] <<"-" << input->end1[i]<< " ";
      start2 << input->start2[i] <<"-" << input->end2[i]<< " ";
    }
  cout << "Id1: " << input->Id1 << " Id2: " << input->Id2 << " BC: " << input->BC<<"  " << start1.str()<<"  " << start2.str()<<" step:"<<input->step2[0]<<" "<< input->step2[1]<<" "<< input->step2[2]<<" order:"<<input->order[0]<<" "<< input->order[1]<<" "<< input->order[2]<<" BC:"<<input->BC<<" Nstar:"<<input->Nstar<< endl;
  // cout << start1.str() << endl;
  // cout << start2.str() << endl;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapZero(ZFSStrctrdWindowMap* output)
{
  //output=new ZFSStrctrdWindowMap(nDim);
  output->Id1=0;
  output->Id2=0;
  for(ZFSInt i=0; i<nDim; i++)
    {
      output->start1[i]=0;
      output->end1[i]=-1;
      output->step1[i]=1;
      output->start2[i]=0;
      output->end2[i]=-1;
      output->step2[i]=1;
      output->order[i]=0;
      output->BC=-1;
    }
  output->hasSponge=false;
  output->spongeThickness=F0;
  output->beta=F0;
  output->sigma=F0;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapInvert(ZFSStrctrdWindowMap* input,ZFSStrctrdWindowMap* output)
{
  //output= new ZFSStrctrdWindowMap(nDim);
  output->Id1=input->Id2;
  output->Id2=input->Id1;

  memcpy(output->start2, input->start1, nDim*sizeof(ZFSInt));

  memcpy(output->end2, input->end1, nDim*sizeof(ZFSInt));
  memcpy(output->step2, input->step1, nDim*sizeof(ZFSInt));
  memcpy(output->start1, input->start2, nDim*sizeof(ZFSInt));
  memcpy(output->end1, input->end2, nDim*sizeof(ZFSInt));
  memcpy(output->step1, input->step2, nDim*sizeof(ZFSInt));

  for(ZFSInt i=0; i<nDim; i++)
    {

      output->order[input->order[i]]=i;
    }

}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapInvert1(ZFSStrctrdWindowMap* output)
{
  //output= new ZFSStrctrdWindowMap(nDim);
  ZFSStrctrdWindowMap* input;
  input = new ZFSStrctrdWindowMap(nDim);
  mapCpy(output, input);
  output->Id1=input->Id2;
  output->Id2=input->Id1;

  output->dc1=input->dc2;
  output->dc2=input->dc1;

  memcpy(output->start2, input->start1, nDim*sizeof(ZFSInt));

  memcpy(output->end2, input->end1, nDim*sizeof(ZFSInt));
  memcpy(output->step2, input->step1, nDim*sizeof(ZFSInt));
  memcpy(output->start1, input->start2, nDim*sizeof(ZFSInt));
  memcpy(output->end1, input->end2, nDim*sizeof(ZFSInt));
  memcpy(output->step1, input->step2, nDim*sizeof(ZFSInt));

  for(ZFSInt i=0; i<nDim; i++)
    {
      output->order[input->order[i]]=i;
    }

  delete input;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCpy(ZFSStrctrdWindowMap* input, ZFSStrctrdWindowMap* output)
{
  output->Id1=input->Id1;
  output->Id2=input->Id2;
  memcpy(output->start1, input->start1, nDim*sizeof(ZFSInt));
  memcpy(output->start2, input->start2, nDim*sizeof(ZFSInt));
  memcpy(output->end1, input->end1, nDim*sizeof(ZFSInt));
  memcpy(output->end2, input->end2, nDim*sizeof(ZFSInt));
  memcpy(output->step1, input->step1, nDim*sizeof(ZFSInt));
  memcpy(output->step2, input->step2, nDim*sizeof(ZFSInt));
  memcpy(output->order, input->order, nDim*sizeof(ZFSInt));
  output->BC=input->BC;
  output->Nstar=input->Nstar;
  output->SingularId=input->SingularId;
  output->dc1=input->dc1;
  output->dc2=input->dc2;
  memcpy(output->SingularBlockId, input->SingularBlockId, 4*sizeof(ZFSInt));
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapNormalize3(ZFSStrctrdWindowMap* output)
{
  //output = new ZFSStrctrdWindowMap(nDim);
  ZFSStrctrdWindowMap* temp;
  temp = new ZFSStrctrdWindowMap(nDim);
  mapCpy(output, temp);
  for(ZFSInt i=0; i<nDim; i++)
    {
      if(temp->step1[i] < 0)
        {
          output->step1[i]=-temp->step1[i];
          output->start1[i]=temp->end1[i];
          output->end1[i]=temp->start1[i];
          output->step2[output->order[i]]=-temp->step2[temp->order[i]];
          output->start2[output->order[i]]=temp->end2[temp->order[i]];
          output->end2[output->order[i]]=temp->start2[temp->order[i]];
        }
    }
  delete temp;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapNormalize1(ZFSStrctrdWindowMap* input, ZFSStrctrdWindowMap* output)
{
  //output = new ZFSStrctrdWindowMap(nDim);
  mapCpy(input, output);
  for(ZFSInt i=0; i<nDim; i++)
    {
      if(input->step1[i] < 0)
        {
          output->step1[i]=-input->step1[i];
          output->start1[i]=input->end1[i];
          output->end1[i]=input->start1[i];
          output->step2[output->order[i]]=-input->step2[input->order[i]];
          output->start2[output->order[i]]=input->end2[input->order[i]];
          output->end2[output->order[i]]=input->start2[input->order[i]];
        }
    }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapNormalize2(ZFSStrctrdWindowMap* input, ZFSStrctrdWindowMap* output)
{
  ZFSStrctrdWindowMap out(nDim);
  ZFSStrctrdWindowMap out1(nDim);
  mapInvert(input, &out);
  mapNormalize1(&out, &out1);
  mapInvert(&out1, output);
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCombine11(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2, ZFSStrctrdWindowMap* output)
{
  ZFSStrctrdWindowMap tmp(nDim);
  ZFSInt shift, shift2;
  if(input1->Id1!=input2->Id1)
    {
      mapZero(output);
    }
  else
    {
      mapNormalize1(input1, output);
      mapNormalize1(input2, &tmp);
      for(ZFSInt i=0; i<nDim; i++)
        {
          shift=0;
          while(output->start1[i]+shift*output->step1[i]<=output->end1[i])
            {
              if((output->start1[i]+shift*output->step1[i]>=tmp.start1[i])&&(((output->start1[i]+shift*output->step1[i])-(tmp.start1[i]))%(tmp.step1[i]))==0) break;
              shift=shift+1;
            }
          output->start1[i]=output->start1[i]+shift*output->step1[i];
          output->start2[output->order[i]]=output->start2[output->order[i]]+shift* output->step2[output->order[i]];

          shift=0;

          while((output->end1[i]-shift*output->step1[i])>=output->start1[i])
            {
              if(((output->end1[i]-shift*output->step1[i])<= tmp.end1[i])&& ((output->end1[i]-shift*output->step1[i]-tmp.end1[i])%tmp.step1[i])==0) break;
              ++shift;
            }

          output->end1[i]=output->end1[i]-shift*output->step1[i];
          output->end2[output->order[i]]=output->end2[output->order[i]]-shift*output->step2[output->order[i]];

          shift=0;

          while(tmp.start1[i]+shift*tmp.step1[i]<=tmp.end1[i])
            {
              if((tmp.start1[i]+shift*tmp.step1[i]>=output->start1[i])&& ((tmp.start1[i]+(shift*tmp.step1[i])-output->start1[i])%output->step1[i])==0) break;
              ++shift;
            }

          tmp.start1[i]=tmp.start1[i]+shift*tmp.step1[i];
          tmp.start2[tmp.order[i]]= tmp.start2[tmp.order[i]]+shift*tmp.step2[tmp.order[i]];

          shift=0;
          while(tmp.end1[i]-shift*tmp.step1[i]>=tmp.start1[i])
            {
             if(((tmp.end1[i]-shift*tmp.step1[i])<= output->end1[i])&& ((tmp.end1[i]-shift*tmp.step1[i]-output->end1[i])%output->step1[i])==0) break;
             ++shift;
            }
          tmp.end1[i]=tmp.end1[i]-shift*tmp.step1[i];
          tmp.end2[tmp.order[i]]=tmp.end2[tmp.order[i]]-shift*tmp.step2[tmp.order[i]];


          shift=1;
          shift2=1;

          while(output->step1[i]*shift!=tmp.step1[i]*shift2)
            {
              if(output->step1[i]*shift< tmp.step1[i]*shift2)
                {
                  ++shift;
                }
              else
                {
                  ++shift2;
                }
            }
          output->step1[i]=shift*output->step1[i];
          output->step2[output->order[i]]=shift*output->step2[output->order[i]];
          tmp.step1[i]=shift2*tmp.step1[i];
          tmp.step2[tmp.order[i]]=shift2*tmp.step2[tmp.order[i]];

        }

      ZFSStrctrdWindowMap tmp1(nDim);
      mapCpy(output, &tmp1);
      mapInvert(&tmp1, output);
      for(ZFSInt i=0; i<nDim; i++)
        {
          output->order[i]=tmp.order[output->order[i]];
        }
      output->Id2=tmp.Id2;
      memcpy(output->start2, tmp.start2, nDim*sizeof(ZFSInt));
      memcpy(output->end2, tmp.end2, nDim*sizeof(ZFSInt));
      memcpy(output->step2, tmp.step2, nDim*sizeof(ZFSInt));

      mapCpy(output, &tmp1);
      mapNormalize1(&tmp1, output);
      if(input1->BC==-1 && input2->BC>0)
        {
          output->BC=input2->BC;
        }
      else
        {
          output->BC=input1->BC;
        }
    }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCombine12(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2, ZFSStrctrdWindowMap* output)
{
  ZFSStrctrdWindowMap tmp(nDim);
  mapInvert(input2, &tmp);
  mapCombine11(input1, &tmp, output);
  //return tmp;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCombine21(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2, ZFSStrctrdWindowMap* output)
{
  ZFSStrctrdWindowMap tmp(nDim);
  mapInvert(input1, &tmp);
  mapCombine11(&tmp, input2, output);
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCombine22(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2, ZFSStrctrdWindowMap* output)
{
  ZFSStrctrdWindowMap tmp(nDim);
  ZFSStrctrdWindowMap tmp1(nDim);
  mapInvert(input1, &tmp);
  mapInvert(input2, &tmp1);
  mapCombine11(&tmp, &tmp1, output);
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCombineCell11(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2, ZFSStrctrdWindowMap* output)
{
  ZFSStrctrdWindowMap tmp(nDim);
  ZFSInt shift, shift2;
  if(input1->Id1!=input2->Id1)
    {
      mapZero(output);
    }
  else
    {
      mapNormalize1(input1, output);
      mapNormalize1(input2, &tmp);
      for(ZFSInt i=0; i<nDim; i++)
        {
          if(output->start1[i]==output->end1[i]||tmp.start1[i]==tmp.end1[i])
            {
              shift=0;
              while(output->start1[i]+shift*output->step1[i]<=output->end1[i])
                {
                  if((output->start1[i]+shift*output->step1[i]>=tmp.start1[i])&&(((output->start1[i]+shift*output->step1[i])-(tmp.start1[i]))%(tmp.step1[i]))==0) break;
                  shift=shift+1;
                }
              output->start1[i]=output->start1[i]+shift*output->step1[i];
              output->start2[output->order[i]]=output->start2[output->order[i]]+shift* output->step2[output->order[i]];

              shift=0;

              while((output->end1[i]-shift*output->step1[i])>=output->start1[i])
                {
                  if(((output->end1[i]-shift*output->step1[i])<= tmp.end1[i])&& ((output->end1[i]-shift*output->step1[i]-tmp.end1[i])%tmp.step1[i])==0) break;
                  ++shift;
                }

              output->end1[i]=output->end1[i]-shift*output->step1[i];
              output->end2[output->order[i]]=output->end2[output->order[i]]-shift*output->step2[output->order[i]];

              shift=0;

              while(tmp.start1[i]+shift*tmp.step1[i]<=tmp.end1[i])
                {
                  if((tmp.start1[i]+shift*tmp.step1[i]>=output->start1[i])&& ((tmp.start1[i]+(shift*tmp.step1[i])-output->start1[i])%output->step1[i])==0) break;
                  ++shift;
                }

              tmp.start1[i]=tmp.start1[i]+shift*tmp.step1[i];
              tmp.start2[tmp.order[i]]= tmp.start2[tmp.order[i]]+shift*tmp.step2[tmp.order[i]];

              shift=0;
              while(tmp.end1[i]-shift*tmp.step1[i]>=tmp.start1[i])
                {
                  if(((tmp.end1[i]-shift*tmp.step1[i])<= output->end1[i])&& ((tmp.end1[i]-shift*tmp.step1[i]-output->end1[i])%output->step1[i])==0) break;
                  ++shift;
                }
              tmp.end1[i]=tmp.end1[i]-shift*tmp.step1[i];
              tmp.end2[tmp.order[i]]=tmp.end2[tmp.order[i]]-shift*tmp.step2[tmp.order[i]];
            }



          else
            {

              shift=1;
              output->end1[i]=output->end1[i]-shift*output->step1[i];
              output->end2[output->order[i]]=output->end2[output->order[i]]-shift*output->step2[output->order[i]];
              tmp.end1[i]=tmp.end1[i]-shift*tmp.step1[i];
              tmp.end2[tmp.order[i]]=tmp.end2[tmp.order[i]]-shift*tmp.step2[tmp.order[i]];



              shift=0;
              while(output->start1[i]+shift*output->step1[i]<=output->end1[i])
                {
                  if((output->start1[i]+shift*output->step1[i]>=tmp.start1[i])&&(((output->start1[i]+shift*output->step1[i])-(tmp.start1[i]))%(tmp.step1[i]))==0) break;
                  shift=shift+1;
                }
              output->start1[i]=output->start1[i]+shift*output->step1[i];
              output->start2[output->order[i]]=output->start2[output->order[i]]+shift* output->step2[output->order[i]];

              shift=0;

              while((output->end1[i]-shift*output->step1[i])>=output->start1[i])
                {
                  if(((output->end1[i]-shift*output->step1[i])<= tmp.end1[i])&& ((output->end1[i]-shift*output->step1[i]-tmp.end1[i])%tmp.step1[i])==0) break;
                  ++shift;
                }

              output->end1[i]=output->end1[i]-shift*output->step1[i];
              output->end2[output->order[i]]=output->end2[output->order[i]]-shift*output->step2[output->order[i]];


              shift=0;

              while(tmp.start1[i]+shift*tmp.step1[i]<=tmp.end1[i])
                {
                  if((tmp.start1[i]+shift*tmp.step1[i]>=output->start1[i])&& ((tmp.start1[i]+(shift*tmp.step1[i])-output->start1[i])%output->step1[i])==0) break;
                  ++shift;
                }

              tmp.start1[i]=tmp.start1[i]+shift*tmp.step1[i];
              tmp.start2[tmp.order[i]]= tmp.start2[tmp.order[i]]+shift*tmp.step2[tmp.order[i]];

              shift=0;
              while(tmp.end1[i]-shift*tmp.step1[i]>=tmp.start1[i])
                {
                  if(((tmp.end1[i]-shift*tmp.step1[i])<= output->end1[i])&& ((tmp.end1[i]-shift*tmp.step1[i]-output->end1[i])%output->step1[i])==0) break;
                  ++shift;
                }
              tmp.end1[i]=tmp.end1[i]-shift*tmp.step1[i];
              tmp.end2[tmp.order[i]]=tmp.end2[tmp.order[i]]-shift*tmp.step2[tmp.order[i]];


              if(output->start1[i]<=output->end1[i]&&tmp.start1[i]<=tmp.end1[i])
                {
                  shift=-1;
                  output->end1[i]=output->end1[i]-shift*output->step1[i];
                  output->end2[output->order[i]]=output->end2[output->order[i]]-shift*output->step2[output->order[i]];
                  tmp.end1[i]=tmp.end1[i]-shift*tmp.step1[i];
                  tmp.end2[tmp.order[i]]=tmp.end2[tmp.order[i]]-shift*tmp.step2[tmp.order[i]];
                }

            }

          shift=1;
          shift2=1;

          while(output->step1[i]*shift!=tmp.step1[i]*shift2)
            {
              if(output->step1[i]*shift< tmp.step1[i]*shift2)
                {
                  ++shift;
                }
              else
                {
                  ++shift2;
                }
            }
          output->step1[i]=shift*output->step1[i];
          output->step2[output->order[i]]=shift*output->step2[output->order[i]];
          tmp.step1[i]=shift2*tmp.step1[i];
          tmp.step2[tmp.order[i]]=shift2*tmp.step2[tmp.order[i]];

        }

      ZFSStrctrdWindowMap tmp1(nDim);
      mapCpy(output, &tmp1);
      mapInvert(&tmp1, output);
      for(ZFSInt i=0; i<nDim; i++)
        {
          output->order[i]=tmp.order[output->order[i]];
        }
      output->Id2=tmp.Id2;
      memcpy(output->start2, tmp.start2, nDim*sizeof(ZFSInt));
      memcpy(output->end2, tmp.end2, nDim*sizeof(ZFSInt));
      memcpy(output->step2, tmp.step2, nDim*sizeof(ZFSInt));

      mapCpy(output, &tmp1);
      mapNormalize1(&tmp1, output);
      if(input1->BC==-1 && input2->BC>0)
        {
          output->BC=input2->BC;
        }
      else
        {
          output->BC=input1->BC;
        }
    }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCombineWave(ZFSStrctrdWindowMap*input1,ZFSStrctrdWindowMap* input2, ZFSStrctrdWindowMap* output)
{
  ZFSStrctrdWindowMap tmp(nDim);
  ZFSInt shift;
  if(input1->Id1!=input2->Id1)
    {
      mapZero(output);
    }
  else
    {
      mapNormalize1(input1, output);
      mapNormalize1(input2, &tmp);
      for(ZFSInt i=0; i<nDim; i++)
    {
      shift=0;
      //go from map1.start1 to map1.end1 and find starting point of
      //overlapping part
      while(output->start1[i]+shift*output->step1[i]<=output->end1[i])
        {
          //stop incrementing if map1.start1 >= map2.start2

if((output->start1[i]+shift*output->step1[i]>=tmp.start1[i])&&(((output->start1[i]+shift*output->step1[i])-(tmp.start1[i]))%(tmp.step1[i]))==0)
break;
          shift=shift+1;
        }

      output->start1[i]=output->start1[i]+shift*output->step1[i];

output->start2[output->order[i]]=output->start2[output->order[i]]+shift*
output->step2[output->order[i]];

      shift=0;

      //go from map1.end1 to map1.start1 and find ending point of
      //overlapping part
      while((output->end1[i]-shift*output->step1[i])>=output->start1[i])
        {
          //stop incrementing if map1.end1 <= map2.end2
          if(((output->end1[i]-shift*output->step1[i])<= tmp.end1[i])&&
((output->end1[i]-shift*output->step1[i]-tmp.end1[i])%tmp.step1[i])==0)
break;
          ++shift;
        }

      output->end1[i]=output->end1[i]-shift*output->step1[i];

output->end2[output->order[i]]=output->end2[output->order[i]]-shift*output->step2[output->order[i]];

      shift=0;

    }

      output->Id1=output->Id2;
      output->Id2=tmp.Id2;
      output->BC=tmp.BC;
    }
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCombineCell12(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2, ZFSStrctrdWindowMap* output)
{
  ZFSStrctrdWindowMap tmp(nDim);
  mapInvert(input2, &tmp);
  mapCombineCell11(input1, &tmp, output);
  //return tmp;
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCombineCell21(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2, ZFSStrctrdWindowMap* output)
{
  ZFSStrctrdWindowMap tmp(nDim);
  mapInvert(input1, &tmp);
  mapCombineCell11(&tmp, input2, output);
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::mapCombineCell22(ZFSStrctrdWindowMap* input1,ZFSStrctrdWindowMap* input2, ZFSStrctrdWindowMap* output)
{
  ZFSStrctrdWindowMap tmp(nDim);
  ZFSStrctrdWindowMap tmp1(nDim);
  mapInvert(input1, &tmp);
  mapInvert(input2, &tmp1);
  mapCombineCell11(&tmp, &tmp1, output);
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::readMapFromArray(ZFSStrctrdWindowMap* map, ZFSInt* dataField) {
    for(ZFSId dim=0; dim<nDim; dim++) {
      map->start1[dim] = dataField[0*nDim+dim];
      map->end1[dim] = dataField[1*nDim+dim];
      map->step1[dim] = dataField[2*nDim+dim];
      map->start2[dim] = dataField[3*nDim+dim];
      map->end2[dim] = dataField[4*nDim+dim];
      map->step2[dim] = dataField[5*nDim+dim];
      map->order[dim] = dataField[6*nDim+dim];
    }

    map->Id1 = dataField[7*nDim+0];
    map->Id2 = dataField[7*nDim+1];
    map->nDim = dataField[7*nDim+2];
    map->BC = dataField[7*nDim+3];
    map->constIndex = dataField[7*nDim+4] ;
    map->face = dataField[7*nDim+5];
    map->dir = dataField[7*nDim+6];
    map->dc1 = dataField[7*nDim+7];
    map->dc2 = dataField[7*nDim+8];
    map->originShape = dataField[7*nDim+9];
    map->hasSponge = dataField[7*nDim+10];
    map->spongeThickness = dataField[7*nDim+11];
    map->beta = dataField[7*nDim+12];
    map->sigma = dataField[7*nDim+13];
    map->Nstar = dataField[7*nDim+14];
    map->SingularId = dataField[7*nDim+15];
}

template <ZFSInt nDim>
void ZFSStrctrdBlckWindowInfo<nDim>::writeMapToArray(ZFSStrctrdWindowMap* map, ZFSInt* dataField) {
      for(ZFSId dim=0; dim<nDim; dim++) {
        dataField[0*nDim+dim] = map->start1[dim];
        dataField[1*nDim+dim] = map->end1[dim];
        dataField[2*nDim+dim] = map->step1[dim];
        dataField[3*nDim+dim] = map->start2[dim];
        dataField[4*nDim+dim] = map->end2[dim];
        dataField[5*nDim+dim] = map->step2[dim];
        dataField[6*nDim+dim] = map->order[dim];
      }

      dataField[7*nDim+0] = map->Id1;
      dataField[7*nDim+1] = map->Id2;
      dataField[7*nDim+2] = map->nDim;
      dataField[7*nDim+3] = map->BC;
      dataField[7*nDim+4] = map->constIndex;
      dataField[7*nDim+5] = map->face;
      dataField[7*nDim+6] = map->dir;
      dataField[7*nDim+7] = map->dc1;
      dataField[7*nDim+8] = map->dc2;
      dataField[7*nDim+9] = map->originShape;
      dataField[7*nDim+10] = map->hasSponge;
      dataField[7*nDim+11] = map->spongeThickness;
      dataField[7*nDim+12] = map->beta;
      dataField[7*nDim+13] = map->sigma;
      dataField[7*nDim+14] = map->Nstar;
      dataField[7*nDim+15] = map->SingularId;
}

// Explicit instantiations for 2D and 3D
template class ZFSStrctrdBlckWindowInfo<2>;
template class ZFSStrctrdBlckWindowInfo<3>;
