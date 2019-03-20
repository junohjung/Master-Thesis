#include "zfsstrctrdblck3d.h"
#include "zfsglobals.h"
#include "zfsconstants.h"
#include "zfsiolib.h"
#include <cmath>


ZFSStrctrdBlck3D::ZFSStrctrdBlck3D( ZFSId blockId, bool* propertiesGroups, const MPI_Comm comm): ZFSStrctrdBlck<3>( blockId, propertiesGroups, comm )
{
  TRACE();
  const ZFSLong oldAllocatedBytes = allocatedBytes();

  // count the no of necessary FQ fields and allocate
  initializeFQField();

  // compute the cell center coordinates from point coordinates
  computeCellCentreCoordinates();

  if(m_rans) {
    m_strctrdBndryCnd = new ZFSStrctrdBndryCnd3D<true>( this, m_noSpecies );
  } else {
    m_strctrdBndryCnd = new ZFSStrctrdBndryCnd3D<false>( this, m_noSpecies );
  }

  allocateSingularities();

  // allocate memory for aux data maps (cf,cp)
  allocateAuxDataMaps();

  // assign coordinates to all ghost points
  addGhostPointCoordinateValues();

  //if we are Rans we should allocate a new RANS block
  if(m_rans==true){m_ransBlck = new ZFSStrctrdBlck3DRans(this);}

  initFluxMethod();

  //allocate and compute metrics and jacobian
  allocateMetrics();
  allocateJacobian();
  computeCellCentreCoordinates();
  computeMetrics();
  computeJacobian();
   
  m_convergence=false;

  computeZonalConnections();     //junoh
  
  if(m_zonal || m_rans) {m_strctrdBndryCnd->computeWallDistances();}
  // cout<<"globalTimeStep:"<<globalTimeStep<<" m_restartTimeStep:"<<m_restartTimeStep<<endl;
  //Assign handlers to the correct boundary conditions
  assignBndryCells();

  if(m_lineOutputInterval>0){initLineInterpolation();}

  printAllocatedMemory( oldAllocatedBytes, "ZFSStrctrdBlck3D",
                        m_zfsStrctrdComm );
}

ZFSStrctrdBlck3D::~ZFSStrctrdBlck3D()
{
  TRACE();

  delete m_strctrdBndryCnd;

  if(m_rans) {
    delete m_ransBlck;
  }

  if(m_hasSingularity) {
    delete [] m_singularity;
  }
}
///////////////////////////////////////////////////////////////
///////////////////////junoh's work for Master Thesis//////////
///////////////////////////////////////////////////////////////
void ZFSStrctrdBlck3D::computeZonalConnections() {   
  if(domainId()==0){
    cout<<" ZonalConnections starting... "<<endl;
    cout << "////////////////////////////////////////////////" << endl;
    cout << "//////////creating  zonal connections //////////" << endl;
    cout << "////////////////////////////////////////////////" << endl;

  }
  zfs_log<<" ZonalConnections starting... "<<endl;
  //////////////////////////////////////
  ////Find the number of ZonalBC////////
  //////////////////////////////////////
  const ZFSInt noZonalBCMaps = m_windowInfo->m_zonalBCMaps.size();

  ZFSIntScratchSpace receiverInputBlockId(noZonalBCMaps, __CALLING_FUNCTION__, "receiverInputBlockId");

  for(ZFSInt id=0; id<noZonalBCMaps; id++) {
    receiverInputBlockId[id] = m_windowInfo->m_zonalBCMaps[id]->Id1;
  }

 
  const ZFSInt nDim = 3; 
  m_zonalBC = new ZFSStrctrdZonalBC[noZonalBCMaps];    //create the object for saving all the data 

  ///////////////////////////////////
  ///// Start ZonalBC Loop //////////
  ///////////////////////////////////
  for(ZFSInt id=0; id<noZonalBCMaps; id++){
      if(domainId()==0){
	cout<<" ZonalBC "<<id<<" Loop  starting... noBlock: "<<receiverInputBlockId[id]<<endl;
  }
    //set up the data of the zonal domain and noCellsBC
    ZFSIntScratchSpace noCellsBCArray(noDomains(), __CALLING_FUNCTION__, "noCellsBCArray"); 
    ZFSInt noCellsGlobalBC=0;
    ZFSInt noCellsLocalBC=0;
    ZFSInt localBCMapId = -1;
    ZFSBool hasLocalBCMap = false;
    ZFSInt stgIP=0;
    ZFSBool hasSTG = false;
    ZFSBool hasZonalwithoutSTG = false;
    ZFSInt hasRcvDomain = 0;
    // if(m_zonalBC[id].m_hasSTG){
    m_zonalBC[id].m_noZonalVariables = 6;
    //}
    m_zonalBC[id].m_zonalBCCells = new ZFSInt[nDim];
    m_zonalBC[id].m_zonalBCCells[0]=m_nCells[0];
    m_zonalBC[id].m_zonalBCCells[1]=m_nCells[1];
    
    for(ZFSInt bcId=0; bcId < abs((ZFSInt)m_windowInfo->physicalBCMap.size()); ++bcId) {

      hasLocalBCMap = m_windowInfo->checkZonalBCMaps(m_windowInfo->m_zonalBCMaps[id], m_windowInfo->physicalBCMap[bcId]);

      if(hasLocalBCMap){
	if(m_windowInfo->physicalBCMap[bcId]->BC == 2221){
	  m_zonalBC[id].m_zonalBCCells[2]=3;
	  stgIP=1;
	  hasSTG = true;
	}
	else if(m_windowInfo->physicalBCMap[bcId]->BC == 2222){
	  hasZonalwithoutSTG= true;
	  m_zonalBC[id].m_zonalBCCells[2]=m_nCells[2];
	  // stgIP=0;      if you want to collect more than 2rows(ghostscell), set up stgIP to the number you want   
	    }
	ZFSInt* start = m_strctrdBndryCnd->m_physicalBCMap[bcId]->start1;
	ZFSInt* end = m_strctrdBndryCnd->m_physicalBCMap[bcId]->end1;  
	noCellsLocalBC += (end[0]+stgIP-start[0])*(end[1]-start[1])*(end[2]-start[2]);
	localBCMapId = bcId;
	hasRcvDomain++;
	break;
      }
    }

      if(hasLocalBCMap == false) {
	noCellsLocalBC = 1;
      }
      ZFSFloatScratchSpace localCoordinatesBCX(noCellsLocalBC, __CALLING_FUNCTION__, "localCoordinatesBCX");
      ZFSFloatScratchSpace localCoordinatesBCY(noCellsLocalBC, __CALLING_FUNCTION__, "localCoordinatesBCY");
      ZFSFloatScratchSpace localCoordinatesBCZ(noCellsLocalBC, __CALLING_FUNCTION__, "localCoordinatesBCZ");
      ZFSIdScratchSpace localReceiverIds(noCellsLocalBC, __CALLING_FUNCTION__, "localReceiverIds"); 
      ZFSIdScratchSpace localMapCellsId(noCellsLocalBC,__CALLING_FUNCTION__,"localMapCellsId");
      
      if(hasLocalBCMap == false) {
	noCellsLocalBC=0;
      }

      if(hasLocalBCMap) {
	// cout<<"stgIP:"<<stgIP<<"id:"<<id<<"domainId():"<<domainId()<<endl;
	ZFSInt* start = m_strctrdBndryCnd->m_physicalBCMap[localBCMapId]->start1;
	ZFSInt* end = m_strctrdBndryCnd->m_physicalBCMap[localBCMapId]->end1;
	ZFSId cellIdBC = 0;
	for(ZFSId k=start[2]; k<end[2] ; k++) {
	  for(ZFSId j=start[1]; j<end[1] ; j++) {
	    for(ZFSId i=start[0]; i<end[0]+stgIP; i++) {
	      ZFSId cellId = cellIndex(i,j,k);
	      ZFSId bcCellId = m_zonalBC[id].zonalCellIndex(i,j,k);
	      localCoordinatesBCX(cellIdBC) = m_cells->coordinates[0][cellId];
	      localCoordinatesBCY(cellIdBC) = m_cells->coordinates[1][cellId];
	      localCoordinatesBCZ(cellIdBC) = m_cells->coordinates[2][cellId];
	      localMapCellsId(cellIdBC) = bcCellId;   //in order to use zonalScatter
	      localReceiverIds(cellIdBC)=domainId();
	      cellIdBC++;
	    }
	  }
	}     
	// cout << "this BC starts at " << start[0] << " , " << start[1] << " , " << start[2] << endl;
	// cout << "this BC ends at " << end[0]+stgIP << " , " << end[1] << " , " << end[2] << endl;
	// cout << "I'm partition " << domainId() << " and I have part of BC 2222 or 2221 " <<"id: "<<id<<endl; 
}

      MPI_Allreduce(&noCellsLocalBC,&noCellsGlobalBC,1,MPI_INT,MPI_SUM,m_zfsStrctrdComm);
      MPI_Allgather(&noCellsLocalBC,1,MPI_INT,&noCellsBCArray[0], 1, MPI_INT,m_zfsStrctrdComm);


 if(domainId()==0){
	cout<<" cellCoordinates gathering..."<<endl;
  }

      //////////////////////////////////////////////////////////////////////////
      ////// Gathering the cellCoordinatesGlobalBC from CoordinatesLocalBC//////
      //////////////////////////////////////////////////////////////////////////
      ZFSIdScratchSpace noCellsBCOffsets(noDomains(), __CALLING_FUNCTION__, "noCellsBCOffsets");
      ZFSFloatScratchSpace globalCoordinatesBCX(noCellsGlobalBC,__CALLING_FUNCTION__, "globalCoordinatesBCX"); 
      ZFSFloatScratchSpace globalCoordinatesBCY(noCellsGlobalBC,__CALLING_FUNCTION__, "globalCoordinatesBCY");
      ZFSFloatScratchSpace globalCoordinatesBCZ(noCellsGlobalBC,__CALLING_FUNCTION__, "globalCoordinatesBCZ");
      ZFSIdScratchSpace globalReceiverIds(noCellsGlobalBC, __CALLING_FUNCTION__, "globalReceiverIds");
      ZFSIdScratchSpace globalLocalMapCellIds(noCellsGlobalBC,__CALLING_FUNCTION__,"globalLocalMapCellIds");

      noCellsBCOffsets[0]=0;
      for (ZFSId i= 1; i< noDomains(); i++ ){ 
	noCellsBCOffsets[i]= noCellsBCOffsets[i-1] + noCellsBCArray[i-1];
      }

      MPI_Allgatherv(&localCoordinatesBCX[0], noCellsBCArray[domainId()], MPI_DOUBLE, &globalCoordinatesBCX[0], &noCellsBCArray[0], &noCellsBCOffsets[0], MPI_DOUBLE, m_zfsStrctrdComm);
      MPI_Allgatherv(&localCoordinatesBCY[0], noCellsBCArray[domainId()], MPI_DOUBLE, &globalCoordinatesBCY[0], &noCellsBCArray[0], &noCellsBCOffsets[0], MPI_DOUBLE, m_zfsStrctrdComm);
      MPI_Allgatherv(&localCoordinatesBCZ[0], noCellsBCArray[domainId()], MPI_DOUBLE, &globalCoordinatesBCZ[0], &noCellsBCArray[0], &noCellsBCOffsets[0], MPI_DOUBLE, m_zfsStrctrdComm);
      MPI_Allgatherv(&localReceiverIds[0], noCellsBCArray[domainId()], MPI_INT, &globalReceiverIds[0], &noCellsBCArray[0], &noCellsBCOffsets[0], MPI_INT, m_zfsStrctrdComm);   
      MPI_Allgatherv(&localMapCellsId[0],noCellsBCArray[domainId()], MPI_INT, &globalLocalMapCellIds[0], &noCellsBCArray[0], &noCellsBCOffsets[0], MPI_INT, m_zfsStrctrdComm);


 if(domainId()==0){
	cout<<" Interpolation Partner finding..."<<endl;
  }

      /////////////////////////////////////////
      ////Set Up PrepareInterpolation  ////////
      /////////////////////////////////////////

      ZFSFloat **  m_cellCoordinatesBC;
      ZFSInt *     m_hasPartnerLocalBC;
      ZFSBool      m_hasInterpolationPartnerDomain=true;
      zfsAlloc(m_cellCoordinatesBC, 3, noCellsGlobalBC,"m_cellCoordinatesBC",F0,__CALLING_FUNCTION__);
      zfsAlloc(m_hasPartnerLocalBC, noCellsGlobalBC, "m_hasPartnerLocalBC",0,__CALLING_FUNCTION__);
     
      for(ZFSId i=0; i<noCellsGlobalBC; i++){
	m_cellCoordinatesBC[0][i]=globalCoordinatesBCX(i);
	m_cellCoordinatesBC[1][i]=globalCoordinatesBCY(i);
	m_cellCoordinatesBC[2][i]=globalCoordinatesBCZ(i);
      } 
      
	if(m_inputBlockId==receiverInputBlockId[id]){ 
	  m_hasInterpolationPartnerDomain=false;
	}
      //PrepareInterpolation
      ZFSStrctrdInterpolation<3>* m_zonalInterpolation = new ZFSStrctrdInterpolation<3>(m_nCells, m_cells->coordinates, m_zfsStrctrdComm); 
      m_zonalInterpolation->prepareZonalInterpolation(noCellsGlobalBC, m_cellCoordinatesBC, m_hasPartnerLocalBC,m_hasInterpolationPartnerDomain);
      m_zonalV.push_back(m_zonalInterpolation);   // m_zonalV is the vector in order to put the Interpolation data in each zonal BC

      /////////////////////////////////////////////////////////////////////////////////////////	
       
	// ZFSInt *     m_hasPartnerGlobalmodified;      
	// zfsAlloc(m_hasPartnerGlobalmodified, noCellsGlobalBC, "m_hasPartnerGlobalmodified",0,__CALLING_FUNCTION__);
	// MPI_Allreduce(m_hasPartnerLocalBC, m_hasPartnerGlobalmodified, noCellsGlobalBC, MPI_INT, MPI_SUM, m_zfsStrctrdComm);

 if(domainId()==0){
	cout<<" Snd and Rcv Domains finding..."<<endl;
  }

	/////////////////////////////////////////////////////////////////////////////////////////
	//////////////////// global sndDomain,rcvDomain Info in each zonalBC///////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////
	// to find out which domains will be rcvDomains
	ZFSInt noGlobalRcvDomains =0;
	MPI_Allreduce(&hasRcvDomain, &noGlobalRcvDomains,1,MPI_INT,MPI_SUM,m_zfsStrctrdComm);
	ZFSIdScratchSpace globalRcvZonalId(noGlobalRcvDomains,__CALLING_FUNCTION__,"globalRcvZonalId");
   
	ZFSInt pos=0;
	for(ZFSId j=0; j<noDomains(); j++){
	  if(noCellsBCArray[j]>0){
	    globalRcvZonalId[pos++]=j; 
	  }
	}
     
	//to find out which domains will be sndDomains
	ZFSInt hasSndDomain=false;
	for(ZFSId cellIdBC=0; cellIdBC<noCellsGlobalBC; cellIdBC++){
	  if(m_hasPartnerLocalBC[cellIdBC]){
	    hasSndDomain=true;
	    break;
	  }
	}
     
      ZFSInt noGlobalSndDomains=0;
      MPI_Allreduce(&hasSndDomain,&noGlobalSndDomains,1,MPI_INT,MPI_SUM,m_zfsStrctrdComm);
    
      ZFSIdScratchSpace hasSndDomainInfo(noDomains(),__CALLING_FUNCTION__,"hasSndDomainInfo");
      MPI_Allgather(&hasSndDomain,1,MPI_INT,&hasSndDomainInfo[0],1,MPI_INT,m_zfsStrctrdComm);

      ZFSIdScratchSpace globalSndZonalId(noGlobalSndDomains,__CALLING_FUNCTION__,"globalSndZonalId");
      pos=0; 
	for(ZFSId j=0; j<noDomains(); j++){
	  if(hasSndDomainInfo[j]>0){
	    globalSndZonalId[pos++]=j;
	  }
	}
	/////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////local send and receive domains Infomation//////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	//get SndSize for zonalexchange        
	ZFSInt noBufferSndSize=0;
	for(ZFSId cellIdBC=0; cellIdBC<noCellsGlobalBC; cellIdBC++){
	  if(m_hasPartnerLocalBC[cellIdBC]){
	    noBufferSndSize++;
	  }
	}
    

	if(hasSndDomain==false){
	  noBufferSndSize=1;
	}
	ZFSIdScratchSpace localBufferMapCellId(noBufferSndSize,__CALLING_FUNCTION__,"localBufferMapCellId");
	ZFSIdScratchSpace localBufferIndexCellId(noBufferSndSize,__CALLING_FUNCTION__,"localBufferIndexCellId");
	ZFSIdScratchSpace localCommReceiverIds(noBufferSndSize,__CALLING_FUNCTION__,"localCommReceiverIds");


       	if(hasSndDomain==false){
	  noBufferSndSize=0;
	}

	//get noSndNghbrDomains
	//first we get localCommReceiverIds
	pos=0;
	for(ZFSId cellIdBC=0; cellIdBC<noCellsGlobalBC; cellIdBC++){
	  if(m_hasPartnerLocalBC[cellIdBC]){
	    localCommReceiverIds[pos]=globalReceiverIds[cellIdBC];
	    pos++;
	  }
	}

	ZFSInt noSndNghbrDomains=0;
	for(ZFSId i=0; i<noGlobalSndDomains; i++){
	  if(globalSndZonalId[i]==domainId()){
	    for(ZFSId Id=0; Id<noDomains(); Id++){
	      for(ZFSId j=0; j<noBufferSndSize; j++){
		if(localCommReceiverIds[j]==Id){
		  noSndNghbrDomains++;
		  break;
		}
	      }
	    }
	  }
	}

       	if(hasSndDomain==false){
	  noSndNghbrDomains=1;
	}

	ZFSIdScratchSpace localRcvId(noSndNghbrDomains,__CALLING_FUNCTION__,"localRcvId");
	ZFSIdScratchSpace localBufferSndSize(noSndNghbrDomains,__CALLING_FUNCTION__,"localBufferSndSize");
  	if(hasSndDomain==false){
	  noSndNghbrDomains=0;
	}


	////////get rcvId in local domains///////////////// 

	pos=0;
	for(ZFSId i=0; i<noGlobalSndDomains; i++){
	  if(globalSndZonalId[i]==domainId()){
	    for(ZFSId Id=0; Id<noDomains(); Id++){
	      for(ZFSId j=0; j<noBufferSndSize; j++){
		if(localCommReceiverIds[j]==Id){
		  localRcvId[pos]=Id;
		  pos++;
		  break;
		}
	      }
	    }
	  }
	}


	//get localBufferMapCellId, localBufferIndexCellId
	pos=0;
	for(ZFSId cellIdBC=0; cellIdBC<noCellsGlobalBC; cellIdBC++){
	  if(m_hasPartnerLocalBC[cellIdBC]){
	    ZFSId MapCellId= globalLocalMapCellIds[cellIdBC];
	    localBufferMapCellId[pos]=MapCellId;
	    localBufferIndexCellId[pos]=cellIdBC;            // this is needed to get InterpolatedVars in zonalGather
	    pos++;
	  }
	}
	
	//////////get localbufferSndSize//////////////////
	
	for(ZFSId i=0; i<noSndNghbrDomains; i++){
	  for(ZFSId j=0; j<noBufferSndSize; j++){
	    if(localRcvId[i]==localCommReceiverIds[j]){
	      localBufferSndSize[i]++;   
	    }
	  }
	}
	//get noRcvNghbrDomains and SndId
	ZFSIdScratchSpace sndBufferRcvSize(noDomains(),__CALLING_FUNCTION__,"sndBufferRcvSize");

	for(ZFSId noRcv=0; noRcv<noGlobalRcvDomains; noRcv++){
	  ZFSInt bufferRcvSize=0;
	  for(ZFSId cellIdBC=0; cellIdBC<noCellsGlobalBC; cellIdBC++){
	    if(globalRcvZonalId[noRcv]==globalReceiverIds[cellIdBC]){
	      if(m_hasPartnerLocalBC[cellIdBC]){
		bufferRcvSize++;
	      }
	    }
	  }
	  //get the corresponding rcvDomainId
	  ZFSId rcvDomainId=globalRcvZonalId[noRcv];
	  //get the bufferSndSize
	  sndBufferRcvSize[rcvDomainId]=bufferRcvSize;
	}
	ZFSIdScratchSpace rcvBufferRcvSize(noDomains(),__CALLING_FUNCTION__,"rcvBufferRcvSize");
	MPI_Alltoall(&sndBufferRcvSize[0],1,MPI_INT,&rcvBufferRcvSize[0],1,MPI_INT,m_zfsStrctrdComm);

//get noRcvNghbrDomains, localSndId, localbufferRcvSize. These are from which domains each rcvDoamin will receive and  how much size it is.

    	ZFSInt noRcvNghbrDomains=0;
    	if(hasLocalBCMap){
    	  for(ZFSId i=0; i<noDomains(); i++){
    	    if(rcvBufferRcvSize[i]>0){                 //get noRcvNghbrDoamins
    	      noRcvNghbrDomains++;
    	    }
    	  }
    	}

    	if(hasLocalBCMap==false){
    	  noRcvNghbrDomains=1;
    	}

    	ZFSIdScratchSpace localSndId(noRcvNghbrDomains,__CALLING_FUNCTION__,"localSndId");
    	ZFSIdScratchSpace localBufferRcvSize(noRcvNghbrDomains,__CALLING_FUNCTION__,"localBufferRcvSize");

    	if(hasLocalBCMap==false){
    	  noRcvNghbrDomains=0;
    	}
	//get localSndId from which domains each rcvDoamin will receive data  
    	pos=0;
    	if(hasLocalBCMap){
        for(ZFSId i=0; i<noDomains(); i++){
    	  if(rcvBufferRcvSize[i]>0){                 //sndId is the info the receiver Domains have from which domains I will receive
    	    localSndId[pos]=i;
    	    pos++;
    	  }
    	}  
    	}
    	//get localbufferRcvSize
	pos=0;
    	if(hasLocalBCMap){
    	  for(ZFSId i=0; i<noDomains(); i++){
    	    if(rcvBufferRcvSize[i]>0){
    	      localBufferRcvSize[pos]=rcvBufferRcvSize[i];
    	      pos++;
    	    }
    	  }
    	}
	



      //////////////////////////////////////////////////////
      ///////////// SET UP COMMUNCATION ////////////////////
      //////////////////////////////////////////////////////

      // // // member Vars dynamic memory
      m_zonalBC[id].m_localMapCellsId = new ZFSInt [noCellsLocalBC];
      m_zonalBC[id].m_globalReceiverIds = new ZFSInt [noCellsGlobalBC];
      m_zonalBC[id].m_globalLocalMapCellIds= new ZFSInt [noCellsGlobalBC];
      m_zonalBC[id].m_globalRcvZonalId = new ZFSInt [noGlobalRcvDomains];
      m_zonalBC[id].m_globalSndZonalId = new ZFSInt [noGlobalSndDomains];
      m_zonalBC[id].mpi_sndRequest= new  MPI_Request[noSndNghbrDomains];
      m_zonalBC[id].mpi_rcvRequest= new MPI_Request[noRcvNghbrDomains];
      m_zonalBC[id].mpi_rcvStatus= new MPI_Status[noRcvNghbrDomains];
      m_zonalBC[id].mpi_sndStatus= new MPI_Status[noSndNghbrDomains];     
      
      m_zonalBC[id].m_localBufferSndSize= new ZFSInt [noSndNghbrDomains];
      m_zonalBC[id].m_localBufferRcvSize= new ZFSInt [noRcvNghbrDomains];
      m_zonalBC[id].m_bufferSndZonal= new ZFSFloat *[noSndNghbrDomains];
      m_zonalBC[id].m_bufferRcvZonal= new ZFSFloat *[noRcvNghbrDomains];
      m_zonalBC[id].m_localRcvId= new ZFSInt [noSndNghbrDomains];
      m_zonalBC[id].m_localSndId= new ZFSInt [noRcvNghbrDomains];
      m_zonalBC[id].m_localBufferMapCellId = new ZFSInt [noBufferSndSize];
      m_zonalBC[id].m_localBufferIndexCellId = new ZFSInt [noBufferSndSize];
      m_zonalBC[id].m_localCommReceiverIds= new ZFSInt [noBufferSndSize];

      m_zonalBC[id].m_bufferSndMapCellId= new ZFSFloat*[noSndNghbrDomains];
      m_zonalBC[id].m_bufferRcvMapCellId= new ZFSFloat*[noRcvNghbrDomains];
      


      for(ZFSInt i=0; i<noSndNghbrDomains; i++){
      	m_zonalBC[id].m_localRcvId[i]= localRcvId[i];
      	m_zonalBC[id].m_localBufferSndSize[i]= localBufferSndSize[i];
	m_zonalBC[id].m_bufferSndZonal[i]= new ZFSFloat [localBufferSndSize[i]*6];
	m_zonalBC[id].m_bufferSndMapCellId[i]= new ZFSFloat [localBufferSndSize[i]];
      }
      
      for(ZFSInt i=0; i<noRcvNghbrDomains; i++){
      	m_zonalBC[id].m_localSndId[i]= localSndId[i];
      	m_zonalBC[id].m_localBufferRcvSize[i]=localBufferRcvSize[i];
	m_zonalBC[id].m_bufferRcvZonal[i]= new ZFSFloat [localBufferRcvSize[i]*6];
	m_zonalBC[id].m_bufferRcvMapCellId[i]= new ZFSFloat [localBufferRcvSize[i]];
      }
      

      m_zonalBC[id].m_noSndNghbrDomains=noSndNghbrDomains;
      m_zonalBC[id].m_noRcvNghbrDomains=noRcvNghbrDomains;
      m_zonalBC[id].m_noBufferSndSize=noBufferSndSize;
      m_noZonalBCMaps= noZonalBCMaps;
      
      for(ZFSInt i=0; i<noBufferSndSize; i++){
      	m_zonalBC[id].m_localBufferIndexCellId[i]=localBufferIndexCellId[i];
	m_zonalBC[id].m_localBufferMapCellId[i]=localBufferMapCellId[i];
	m_zonalBC[id].m_localCommReceiverIds[i]=localCommReceiverIds[i];
      }
      
      // Store the Globaldata
      m_zonalBC[id].m_noCellsGlobalBC = noCellsGlobalBC;
      m_zonalBC[id].m_noCellsLocalBC = noCellsLocalBC;
      m_zonalBC[id].m_hasLocalBCMap = hasLocalBCMap;
      m_zonalBC[id].m_noGlobalRcvDomains = noGlobalRcvDomains;
      m_zonalBC[id].m_noGlobalSndDomains = noGlobalSndDomains;
      m_zonalBC[id].m_hasSTG = hasSTG;
      m_zonalBC[id].m_hasZonalwithoutSTG = hasZonalwithoutSTG;

      for(ZFSInt cellId=0; cellId<noCellsGlobalBC; cellId++){
      	m_zonalBC[id].m_globalReceiverIds[cellId] = globalReceiverIds(cellId);
      	m_zonalBC[id].m_globalLocalMapCellIds[cellId] = globalLocalMapCellIds(cellId);
      }
 

      
      for(ZFSId i=0; i<noCellsLocalBC; i++){
	m_zonalBC[id].m_localMapCellsId[i]=localMapCellsId[i];
      }
      
      for(ZFSInt i=0; i<noGlobalSndDomains; i++){
      	m_zonalBC[id].m_globalSndZonalId[i]= globalSndZonalId(i);
      }
      for(ZFSInt i=0; i<noGlobalRcvDomains; i++){
      	m_zonalBC[id].m_globalRcvZonalId[i]= globalRcvZonalId(i);
      }

      for(ZFSId i=0; i<noSndNghbrDomains; i++){
      	m_zonalBC[id].mpi_sndRequest[i]=MPI_REQUEST_NULL;
      }

       for(ZFSId i=0; i<noRcvNghbrDomains; i++){
       	m_zonalBC[id].mpi_rcvRequest[i]=MPI_REQUEST_NULL;
       }


 if(domainId()==0){
	cout<<" MapCellId exchanging ..."<<endl;
  }
 

        ///////////////////////////////////////////////////////////////
      /////////////////////MapCellId exchange////////////////////////
      ///////////////////////////////////////////////////////////////
     // cout<<"check mapcellId first starting point"<<endl; 
    for(ZFSInt noSnd=0; noSnd<m_zonalBC[id].m_noSndNghbrDomains; noSnd++){
      ZFSFloat* bufferSnd=m_zonalBC[id].m_bufferSndMapCellId[noSnd];
      pos=0;
      for(ZFSId cellId=0; cellId<m_zonalBC[id].m_noBufferSndSize; cellId++){
	if(m_zonalBC[id].m_localCommReceiverIds[cellId]==m_zonalBC[id].m_localRcvId[noSnd]){
	  bufferSnd[pos]=m_zonalBC[id].m_localBufferMapCellId[cellId];
	  pos++;
	}
      }
    }  
    for(ZFSId i=0; i<m_zonalBC[id].m_noGlobalSndDomains; i++){
      if(m_zonalBC[id].m_globalSndZonalId[i]==domainId()){
	for(ZFSId noSnd=0; noSnd<m_zonalBC[id].m_noSndNghbrDomains; noSnd++){ 
	  ZFSInt err = MPI_Isend((void*)&m_zonalBC[id].m_bufferSndMapCellId[noSnd][0],m_zonalBC[id].m_localBufferSndSize[noSnd],MPI_DOUBLE, m_zonalBC[id].m_localRcvId[noSnd],0,MPI_COMM_WORLD,&m_zonalBC[id].mpi_sndRequest[noSnd]);
	  if(err) cout << "rank " << domainId() << " zonal sending throws error " << endl;
	}
      }
    }  
    for(ZFSId i=0; i<m_zonalBC[id].m_noGlobalRcvDomains; i++){
      if(m_zonalBC[id].m_globalRcvZonalId[i]==domainId()){
	for(ZFSId noRcv=0; noRcv<m_zonalBC[id].m_noRcvNghbrDomains; noRcv++){
	  ZFSInt err = MPI_Irecv((void*)&m_zonalBC[id].m_bufferRcvMapCellId[noRcv][0], m_zonalBC[id].m_localBufferRcvSize[noRcv], MPI_DOUBLE, m_zonalBC[id].m_localSndId[noRcv],0,MPI_COMM_WORLD,&m_zonalBC[id].mpi_rcvRequest[noRcv]);
	  if(err) cout << "rank " << domainId() << " zonal receiving throws error " << endl;
	}
      }
    }
    MPI_Waitall(m_zonalBC[id].m_noSndNghbrDomains,m_zonalBC[id].mpi_sndRequest,m_zonalBC[id].mpi_sndStatus);
    MPI_Waitall(m_zonalBC[id].m_noRcvNghbrDomains,m_zonalBC[id].mpi_rcvRequest,m_zonalBC[id].mpi_rcvStatus);


      if(domainId()==0){
	cout<<" ZonalBC "<<id<<" Loop  SUCCESSFUL! "<<endl;
  }
 
  }

      if(domainId()==0){
	cout<<" ZonalConnections SUCCESSFUL! "<<endl;
  }
      zfs_log<<" ZonalConnections SUCCESSFUL!... "<<endl;
    // saveOutputPartitions(true);
    // zfsTerm(1, __CALLING_FUNCTION__, "Leaving program at end of compute zonal connections!");
} 

/**
 * Zonal Grid Connection
 * \author Junoh Jung, Master Thesis
 * \ date 07.2018 
 */

/////Zonal Exchange Data////////         
void ZFSStrctrdBlck3D::zonalExchange()
{
  // RECORD_TIMER_START(m_tcomm);
  // RECORD_TIMER_START(m_tzonalExchange);
  if(m_zonal){
    // RECORD_TIMER_START(m_tzonalGather);
    zonalGather();
    // RECORD_TIMER_STOP(m_tzonalGather);

    // RECORD_TIMER_START(m_tzonalSend);
    zonalSend();
    // RECORD_TIMER_STOP(m_tzonalSend);

    // RECORD_TIMER_START(m_tzonalReceive);
    zonalReceive();
    // RECORD_TIMER_STOP(m_tzonalReceive);

    // RECORD_TIMER_START(m_tzonalScatter);
    zonalScatter();
    // RECORD_TIMER_STOP(m_tzonalScatter);    
    
  }

  // if(m_periodicConnection) m_strctrdBndryCnd->periodicExchange();

   // RECORD_TIMER_STOP(m_tzonalExchange);
   // RECORD_TIMER_STOP(m_tcomm);
}

/**
 * Zonal Exchange Procedure
 * \author Junoh Jung, Master Thesis
 * \ date 07.2018 
 */


void ZFSStrctrdBlck3D::zonalGather()
{ 
  //gather sending variables in buffer
  for(ZFSId id=0; id<(ZFSInt)m_windowInfo->m_zonalBCMaps.size(); id++){
    for(ZFSInt noSnd=0; noSnd<m_zonalBC[id].m_noSndNghbrDomains; noSnd++){
      ZFSFloat* bufferSnd=m_zonalBC[id].m_bufferSndZonal[noSnd];
      ZFSInt pos=0;
      for(ZFSId cellId=0; cellId<m_zonalBC[id].m_noBufferSndSize; cellId++){
	if(m_zonalBC[id].m_localCommReceiverIds[cellId]==m_zonalBC[id].m_localRcvId[noSnd]){
	  ZFSInt cellIdBC=m_zonalBC[id].m_localBufferIndexCellId[cellId];      
	  bufferSnd[pos+PV->U*m_zonalBC[id].m_localBufferSndSize[noSnd]]=m_zonalV.at(id)->interpolateVariableZonal(m_cells->pvariables[PV->U], cellIdBC);
	  bufferSnd[pos+PV->V*m_zonalBC[id].m_localBufferSndSize[noSnd]]=m_zonalV.at(id)->interpolateVariableZonal(m_cells->pvariables[PV->V], cellIdBC);
	  bufferSnd[pos+PV->W*m_zonalBC[id].m_localBufferSndSize[noSnd]]=m_zonalV.at(id)->interpolateVariableZonal(m_cells->pvariables[PV->W], cellIdBC);
	  bufferSnd[pos+PV->RHO*m_zonalBC[id].m_localBufferSndSize[noSnd]]=m_zonalV.at(id)->interpolateVariableZonal(m_cells->pvariables[PV->RHO], cellIdBC);
	  bufferSnd[pos+PV->P*m_zonalBC[id].m_localBufferSndSize[noSnd]]=m_zonalV.at(id)->interpolateVariableZonal(m_cells->pvariables[PV->P], cellIdBC);
	  bufferSnd[pos+FQ->NU_T*m_zonalBC[id].m_localBufferSndSize[noSnd]]=m_zonalV.at(id)->interpolateVariableZonal(m_cells->fq[FQ->NU_T], cellIdBC);
	  pos++;
	}    

	
      }
    }
  }
}

/**
 * Zonal Gather in zonal exchange
 * \author Junoh Jung, Master Thesis
 * \ date 07.2018 
 */

void ZFSStrctrdBlck3D::zonalSend()     //junoh
{
  for(ZFSId id=0; id<(ZFSInt)m_windowInfo->m_zonalBCMaps.size(); id++){
    //access only sndDomains
    for(ZFSId i=0; i<m_zonalBC[id].m_noGlobalSndDomains; i++){
      if(m_zonalBC[id].m_globalSndZonalId[i]==domainId()){
	for(ZFSId noSnd=0; noSnd<m_zonalBC[id].m_noSndNghbrDomains; noSnd++){ 
	  ZFSInt err = MPI_Isend((void*)&m_zonalBC[id].m_bufferSndZonal[noSnd][0],m_zonalBC[id].m_localBufferSndSize[noSnd]*6,MPI_DOUBLE, m_zonalBC[id].m_localRcvId[noSnd],0,MPI_COMM_WORLD,&m_zonalBC[id].mpi_sndRequest[noSnd]);
	  if(err) cout << "rank " << domainId() << " zonal sending throws error " << endl;
	}
      }
    }  
  }
}
/**
 * Zonal Send in zonal exchange
 * \author Junoh Jung, Master Thesis
 * \ date 07.2018 
 */

void ZFSStrctrdBlck3D::zonalReceive() //junoh
{
  for(ZFSId id=0; id<(ZFSInt)m_windowInfo->m_zonalBCMaps.size(); id++){
    for(ZFSId i=0; i<m_zonalBC[id].m_noGlobalRcvDomains; i++){
      //access only rcvDomains
      if(m_zonalBC[id].m_globalRcvZonalId[i]==domainId()){
	for(ZFSId noRcv=0; noRcv<m_zonalBC[id].m_noRcvNghbrDomains; noRcv++){
	  ZFSInt err = MPI_Irecv((void*)&m_zonalBC[id].m_bufferRcvZonal[noRcv][0], m_zonalBC[id].m_localBufferRcvSize[noRcv]*6, MPI_DOUBLE, m_zonalBC[id].m_localSndId[noRcv],0,MPI_COMM_WORLD,&m_zonalBC[id].mpi_rcvRequest[noRcv]);
	  if(err) cout << "rank " << domainId() << " zonal receiving throws error " << endl;
	}
      }
    }
    MPI_Waitall(m_zonalBC[id].m_noSndNghbrDomains,m_zonalBC[id].mpi_sndRequest,m_zonalBC[id].mpi_sndStatus);
    MPI_Waitall(m_zonalBC[id].m_noRcvNghbrDomains,m_zonalBC[id].mpi_rcvRequest,m_zonalBC[id].mpi_rcvStatus);
  }
}



void ZFSStrctrdBlck3D::zonalScatter()
{
  for(ZFSId id=0; id<(ZFSInt)m_windowInfo->m_zonalBCMaps.size(); id++){
    for(ZFSId noRcv=0; noRcv<m_zonalBC[id].m_noRcvNghbrDomains; noRcv++){
      for(ZFSId i=0; i<m_zonalBC[id].m_localBufferRcvSize[noRcv]; i++){ 
	//get cellId and localCellId
	ZFSId localMapCellId=m_zonalBC[id].m_bufferRcvMapCellId[noRcv][i];
	if(m_zonalBC[id].m_hasZonalwithoutSTG){
	  m_cells->fq[FQ->AVG_U][localMapCellId]=m_zonalBC[id].m_bufferRcvZonal[noRcv][PV->U*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	  m_cells->fq[FQ->AVG_V][localMapCellId]=m_zonalBC[id].m_bufferRcvZonal[noRcv][PV->V*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	  m_cells->fq[FQ->AVG_W][localMapCellId]=m_zonalBC[id].m_bufferRcvZonal[noRcv][PV->W*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	  m_cells->fq[FQ->AVG_RHO][localMapCellId]=m_zonalBC[id].m_bufferRcvZonal[noRcv][PV->RHO*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	  m_cells->fq[FQ->AVG_P][localMapCellId]=m_zonalBC[id].m_bufferRcvZonal[noRcv][PV->P*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	  m_cells->fq[FQ->NU_T][localMapCellId]=m_zonalBC[id].m_bufferRcvZonal[noRcv][FQ->NU_T*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	  } 
	  else if(m_zonalBC[id].m_hasSTG){
	    m_cells->stg_fq[PV->U][localMapCellId]= m_zonalBC[id].m_bufferRcvZonal[noRcv][PV->U*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	    m_cells->stg_fq[PV->V][localMapCellId]= m_zonalBC[id].m_bufferRcvZonal[noRcv][PV->V*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	    m_cells->stg_fq[PV->W][localMapCellId]= m_zonalBC[id].m_bufferRcvZonal[noRcv][PV->W*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	    m_cells->stg_fq[PV->RHO][localMapCellId]= m_zonalBC[id].m_bufferRcvZonal[noRcv][PV->RHO*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	    m_cells->stg_fq[PV->P][localMapCellId]= m_zonalBC[id].m_bufferRcvZonal[noRcv][PV->P*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];
	    m_cells->stg_fq[FQ->NU_T][localMapCellId]= m_zonalBC[id].m_bufferRcvZonal[noRcv][FQ->NU_T*m_zonalBC[id].m_localBufferRcvSize[noRcv]+i];

	}
      }
    }
  }
}
/**
 * Zonal Scatter in zonal exchange
 * \author Junoh Jung, Master Thesis
 * \ date 07.2018 
 */


void ZFSStrctrdBlck3D::initFluxMethod() {
  //choose the convective flux method
  //for performance issues, it is templated
  //(we do not need to write the muscl stuff too often
  //we need to put the limiter also into the template
  //in order to reduce this function!!
  if(m_rans==true){
    reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclRANS;
    viscFluxMethod=&ZFSStrctrdBlck3D::viscousFluxRANS;
  }else{
    viscFluxMethod=&ZFSStrctrdBlck3D::viscousFluxLES;
    //check if limiter
    if(m_limiter)
    {
      switch(string2enum(m_limiterMethod))
      {
      case VENKATAKRISHNAN_MOD:
      {

        /*! \page propertyPage1
          \section venkFactor
          <code>ZFSInt ZFSStrctrdBlck::m_venkFactor </code>\n
          default = <code> 0 </code>\n \n
          Factor for the Venkatakrishnan Limiter.\n
          Possible values are:\n
          <ul>
          <li>Float > 0.0</li>
          </ul>
          Keywords: <i>LIMITER, STRCTRD</i>
        */
        m_venkFactor = *(ZFSContext::getProperty("venkFactor", m_blockId, __CALLING_FUNCTION__, (ZFSFloat*) NULL )->asFloat(0)); //reads the customizable parameter from properties

        zfs_log << "Using VENKATAKRISHNAN MOD limiter with VENK factor of " << m_venkFactor << " !" << endl;
        reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclVenkatakrishnan3D;
        Venkatakrishnan_function = &ZFSStrctrdBlck3D::VENKATAKRISHNAN_MOD_FCT;
        break;
      }

      case VENKATAKRISHNAN:
      {
        zfs_log << "Using VENKATAKRISHNAN limiter!" << endl;
        reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclVenkatakrishnan3D;
        Venkatakrishnan_function = &ZFSStrctrdBlck3D::VENKATAKRISHNAN_FCT;
        break;
      }

      case BARTH_JESPERSON :
      {
        zfs_log << "Using BARTH JESPERSON limiter!" << endl;
        reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclVenkatakrishnan3D;
        Venkatakrishnan_function = &ZFSStrctrdBlck3D::BARTH_JESPERSON_FCT;
        break;
      }
      case MINMOD :
      {
        zfs_log << "Using MINMOD limiter!" << endl;
        reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclMinModLimiter;
        break;
      }
      case ALBADA:
      {
        zfs_log << "Using VAN ALBADA limiter!" << endl;
        reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclAlbada;
        break;
      }
      default:
      {
        stringstream errorMessage;
        errorMessage << "Limiter function " << m_limiterMethod << " not implemented!" << endl;
        zfsTerm(1, __CALLING_FUNCTION__, errorMessage.str());
      }
      }
    }else
      if(m_musclScheme=="Standard"){
        zfs_log << "Using unlimited MUSCL! (standard Formulation)" << endl;
        if(m_ausmScheme=="Standard"){
          zfs_log << "Using standard AUSM central"<<endl;
          switch(CV->noVariables){
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::Muscl_AusmLES<5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::Muscl_AusmLES<6>;break;}
          case 7:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::Muscl_AusmLES<7>;break;}
          default:{
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            zfsTerm(1, __CALLING_FUNCTION__);
          }
          }
        }else if(m_ausmScheme=="PTHRC"){
          zfs_log << "Using AUSM PTHRC"<<endl;
          switch(CV->noVariables){
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::Muscl_AusmLES_PTHRC<5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::Muscl_AusmLES_PTHRC<6>;break;}
          case 7:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::Muscl_AusmLES_PTHRC<7>;break;}
          default:{
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            zfsTerm(1, __CALLING_FUNCTION__);
          }
          }
        }else if(m_ausmScheme=="AUSMDV"){
          zfs_log << "Using AUSMDV"<<endl;
          switch(CV->noVariables){
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::Muscl_AusmDV<5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::Muscl_AusmDV<6>;break;}
          case 7:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::Muscl_AusmDV<7>;break;}
          default:{
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            zfsTerm(1, __CALLING_FUNCTION__);
          }
          }
        }
      } else if(m_musclScheme=="Stretched"){
        zfs_log << "Using unlimited MUSCL (streched Grids)";
        zfsAlloc(m_cells->cellLength, nDim, m_noStrctrdCells, "m_cells->cellLength", -F1, __CALLING_FUNCTION__);
        computeCellLength();
        zfs_log << "Using standard AUSM central"<<endl;
        if(m_ausmScheme=="Standard"){
          switch(CV->noVariables){
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES,5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES,6>;break;}
          case 7:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES,7>;break;}
          default:{
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            zfsTerm(1, __CALLING_FUNCTION__);
          }
          }
        }else if(m_ausmScheme=="PTHRC"){
          zfs_log << "Using AUSM PTHRC"<<endl;
          switch(CV->noVariables){
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES_PTHRC,5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES_PTHRC,6>;break;}
          case 7:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES_PTHRC,7>;break;}
          default:{
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            zfsTerm(1, __CALLING_FUNCTION__);
          }
          }
        }else if(m_ausmScheme=="AUSMDV"){
          zfs_log << "Using AUSMDV"<<endl;
          switch(CV->noVariables){
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmDV,5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmDV,6>;break;}
          case 7:{ reconstructSurfaceData=&ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmDV,7>;break;}
          default:{
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            zfsTerm(1, __CALLING_FUNCTION__);
          }
          }
        }
      }
  }
}


void ZFSStrctrdBlck3D::computeCellLength(){
  //this function can be moved into the MusclSchemeStreched later but for testing it is easier
  // REMEMBER: FOR MOVINg GRIDS THIS NEEDS TO BE CALLED EACH TIME

  ZFSId P1=-1, P2=-1, P3=-1, P4=-1, P5=-1, P6=-1, P7=-1, P8=-1;
  ZFSFloat f1x=F0,f1y=F0,f1z=F0, f2x=F0, f2y=F0,f2z=F0;
  for(ZFSId k=0; k<m_nCells[0]; k++){
    for(ZFSId j=0; j<m_nCells[1]; j++){
      for(ZFSId i=0; i<m_nCells[2];i++){
        const ZFSId cellId=cellIndex(i,j,k);
        P1 = getPointIdFromCell(i,j,k);
        P2 = getPointIdfromPoint( P1, 1, 0, 0 );
        P3 = getPointIdfromPoint( P1, 1, 1, 0 );
        P4 = getPointIdfromPoint( P1, 1, 0, 1 );
        P5 = getPointIdfromPoint( P1, 1, 1, 1 );
        P6 = getPointIdfromPoint( P1, 0, 1, 0 );
        P7 = getPointIdfromPoint( P1, 0, 1, 1 );
        P8 = getPointIdfromPoint( P1, 0, 0, 1 );
        //----------Idirection
        //face 1
        f1x=F1B4*(m_coordinates[0][P1]+m_coordinates[0][P6]+m_coordinates[0][P7]+m_coordinates[0][P8]);
        f1y=F1B4*(m_coordinates[1][P1]+m_coordinates[1][P6]+m_coordinates[1][P7]+m_coordinates[1][P8]);
        f1z=F1B4*(m_coordinates[2][P1]+m_coordinates[2][P6]+m_coordinates[2][P7]+m_coordinates[2][P8]);
        //face 2
        f2x=F1B4*(m_coordinates[0][P2]+m_coordinates[0][P3]+m_coordinates[0][P4]+m_coordinates[0][P5]);
        f2y=F1B4*(m_coordinates[1][P2]+m_coordinates[1][P3]+m_coordinates[1][P4]+m_coordinates[1][P5]);
        f2z=F1B4*(m_coordinates[2][P2]+m_coordinates[2][P3]+m_coordinates[2][P4]+m_coordinates[2][P5]);
        m_cells->cellLength[0][cellId]=sqrt(POW2(f2x-f1x)+POW2(f2y-f1y)+POW2(f2z-f1z));
        //----------Jdirection
        //face 1
        f1x=F1B4*(m_coordinates[0][P1]+m_coordinates[0][P2]+m_coordinates[0][P4]+m_coordinates[0][P8]);
        f1y=F1B4*(m_coordinates[1][P1]+m_coordinates[1][P2]+m_coordinates[1][P4]+m_coordinates[1][P8]);
        f1z=F1B4*(m_coordinates[2][P1]+m_coordinates[2][P2]+m_coordinates[2][P4]+m_coordinates[2][P8]);
        //face 2
        f2x=F1B4*(m_coordinates[0][P3]+m_coordinates[0][P5]+m_coordinates[0][P6]+m_coordinates[0][P7]);
        f2y=F1B4*(m_coordinates[1][P3]+m_coordinates[1][P5]+m_coordinates[1][P6]+m_coordinates[1][P7]);
        f2z=F1B4*(m_coordinates[2][P3]+m_coordinates[2][P5]+m_coordinates[2][P6]+m_coordinates[2][P7]);
        m_cells->cellLength[1][cellId]=sqrt(POW2(f2x-f1x)+POW2(f2y-f1y)+POW2(f2z-f1z));
        //----------Kdirection
        //face 1
        f1x=F1B4*(m_coordinates[0][P1]+m_coordinates[0][P2]+m_coordinates[0][P3]+m_coordinates[0][P6]);
        f1y=F1B4*(m_coordinates[1][P1]+m_coordinates[1][P2]+m_coordinates[1][P3]+m_coordinates[1][P6]);
        f1z=F1B4*(m_coordinates[2][P1]+m_coordinates[2][P2]+m_coordinates[2][P3]+m_coordinates[2][P6]);
        //face 2
        f2x=F1B4*(m_coordinates[0][P4]+m_coordinates[0][P5]+m_coordinates[0][P7]+m_coordinates[0][P8]);
        f2y=F1B4*(m_coordinates[1][P4]+m_coordinates[1][P5]+m_coordinates[1][P7]+m_coordinates[1][P8]);
        f2z=F1B4*(m_coordinates[2][P5]+m_coordinates[2][P5]+m_coordinates[2][P7]+m_coordinates[2][P8]);
        m_cells->cellLength[2][cellId]=sqrt(POW2(f2x-f1x)+POW2(f2y-f1y)+POW2(f2z-f1z));
      }
    }
  }
}


void ZFSStrctrdBlck3D::initSolutionStep()
{
  TRACE();
  //Compute infinity values from property file
  //and (if no restart) fill cells according
  //to the initialCondition property
  initialCondition();

  // if(m_restartInterpolation) {     //junoh
  //   interpolateFromDonor();
  // } else if ( m_restart ) {
  //   loadRestartFile();
  // }     

  if ( m_restart ) {
    loadRestartFile();   
  }

  if(m_restartInterpolation) {
    interpolateFromDonor();
  }

  //timestep will be computed and
  //set if globalTimeStep == 0 or
  //constantTimeStep = 0
  //function needs to be after initMovingGrid
  //else no parameters are known
  setTimeStep();

  //initialize moving grid
  //functions and move grid
  //to correct position
  if (m_movingGrid) {
    initMovingGrid();
  }

  computeCumulativeAverage(false);

  if(m_useSandpaperTrip) {
    initSandpaperTrip();
  }

  //Get the correct values
  //in the exchange ghostcells
  exchange();
  if(m_rans) {
    m_ransBlck->computeTurbViscosity();
  }

  if(m_zonal) {
    // averagedFillGhostCells();
    // reconstructTurbulentVariables();
    spanwiseAvgZonal();
    zonalExchange(); //junoh
   
   // reconstructTurbulentVariables();
  }

  //Call the init function of each BC
  initBndryCnds();

  //Apply boundary conditions
  //and fill the non-exchange ghostcells
  applyBoundaryCondition();
  if(m_rans) {
    m_ransBlck->computeTurbViscosity();
  }


  // globalTimeStep++;
  // saveOutputSolution<true>(false);

  //Check for NaNs
  checkNans();

  computeConservativeVariables();
  // globalTimeStep++;
  // saveOutputPartitions<true>();
  // saveOutputSolution<true>(false);


}

/**
 * Instead of a restart from a restart file
 * this method will interpolate from a donor solution
 * onto the flow field
 * \author Marian Albers 
 */
void ZFSStrctrdBlck3D::interpolateFromDonor() {
  m_strctrdInterpolation = new ZFSStrctrdInterpolation<3>(m_zfsStrctrdComm);
  m_strctrdInterpolation->prepareInterpolationField(m_nCells, m_cells->coordinates);

  if(m_zonal){  //junoh
    m_strctrdInterpolation->interpolateField(m_pvariableNames[PV->U], m_cells->pvariables[PV->U]);
    m_strctrdInterpolation->interpolateField(m_pvariableNames[PV->V], m_cells->pvariables[PV->V]);
    m_strctrdInterpolation->interpolateField(m_pvariableNames[PV->W], m_cells->pvariables[PV->W]);
    m_strctrdInterpolation->interpolateField(m_pvariableNames[PV->P], m_cells->pvariables[PV->P]);
    m_strctrdInterpolation->interpolateField(m_pvariableNames[PV->RHO], m_cells->pvariables[PV->RHO]);
    m_strctrdInterpolation->interpolateField(m_pvariableNames[PV->RANS_FIRST], m_cells->pvariables[PV->RANS_FIRST]);
    m_strctrdInterpolation->interpolateField(FQ->fqNames[FQ->NU_T], m_cells->fq[FQ->NU_T]);
  } else{
    for(ZFSId var=0; var<PV->noVariables; var++) { 
      m_strctrdInterpolation->interpolateField(m_pvariableNames[var], m_cells->pvariables[var]);
    }
  }

  if(m_useSponge && m_spongeLayerType == 2) {
    m_strctrdInterpolation->interpolateField("rho", m_cells->fq[FQ->SPONGE_RHO]);
    m_strctrdInterpolation->interpolateField("rhoE", m_cells->fq[FQ->SPONGE_RHO_E]);
  }

  if(m_useSponge && m_spongeLayerType == 4) {
    m_strctrdInterpolation->interpolateField("rho", m_cells->fq[FQ->SPONGE_RHO]);
  }

  manualInterpolationCorrection();

  delete m_strctrdInterpolation;
  m_strctrdInterpolation = 0;
}



/**
 * Computation of infinity values for the conservative and primitive variables
 * Initialization ot the entire flow field
 * \author Pascal S. Meysonnat
 * \date 16.02.2011
 * modified last 25.02.2011
 */
void ZFSStrctrdBlck3D::initialCondition()
{
  TRACE();
  const ZFSFloat gammaMinusOne = m_gamma - 1.0;
  ZFSFloat UT;
  ZFSFloat pressureCH=F0;

  PV->TInfinity = 1.0 / ( 1.0 + F1B2 * gammaMinusOne * POW2(m_Ma));
  UT = m_Ma * sqrt( PV->TInfinity );
  PV->UInfinity = UT * cos( m_angle[ 0 ] ) * cos( m_angle[ 1 ] );
  PV->VInfinity = UT * sin( m_angle[ 0 ] ) * cos( m_angle[ 1 ] );
  PV->WInfinity = UT * sin( m_angle[ 1 ] );
  PV->VVInfinity[ 0 ] = PV->UInfinity;
  PV->VVInfinity[ 1 ] = PV->VInfinity;
  PV->VVInfinity[ 2 ] = PV->WInfinity;
  PV->PInfinity = pow( PV->TInfinity, (m_gamma / gammaMinusOne)) / m_gamma;

  // compute conservative variables
  CV->rhoInfinity = pow( PV->TInfinity, ( 1.0 / gammaMinusOne ) );
  CV->rhoUInfinity = CV->rhoInfinity * PV->UInfinity;
  CV->rhoVInfinity = CV->rhoInfinity * PV->VInfinity;
  CV->rhoWInfinity = CV->rhoInfinity * PV->WInfinity;
  CV->rhoVVInfinity[ 0 ] = CV->rhoUInfinity;
  CV->rhoVVInfinity[ 1 ] = CV->rhoVInfinity;
  CV->rhoVVInfinity[ 2 ] = CV->rhoWInfinity;
  CV->rhoEInfinity = PV->PInfinity / gammaMinusOne + CV->rhoInfinity
    * ( F1B2 * POW2(UT));

  // internal Reynolds number Re0 = Re / ( rho8*M*sqrt(T8)/T8^F072)
  m_Re0 = m_Re * zfsSUTHERLANDLAW(PV->TInfinity) / (CV->rhoInfinity*m_Ma*sqrt(PV->TInfinity));

  // reference enthalpies (needed for combustion computations)
  m_hInfinity = PV->PInfinity / CV->rhoInfinity * m_gamma / gammaMinusOne;

  // reference time (convection time)
  m_timeRef = UT / m_referenceLength;

  m_deltaP = F0;
  // pressure loss per unit length dp = rho_00 u_tau^2 L / D ) here: D=1.0, L=1;
  // channel: dp = rho_00 u_tau^2 L / D )
  // m_deltaP = POW2( m_Ma * m_ReTau * sqrt(PV->TInfinity) / m_Re  ) * CV->rhoInfinity / m_referenceLength;x
  // result is obtained by making deltap dimensionless with a_0^2 and rho_0
  // pipe: dp = lambda * L/D * rho/2 * u^2, lambda = 0.3164 Re^(-1/4) (Blasius)

  if(m_rans) {
    const ZFSFloat lamVisc = zfsSUTHERLANDLAW(PV->TInfinity);
    const ZFSFloat chi = 0.1 ;
    CV->ransInfinity[0] = chi*(lamVisc);
    PV->ransInfinity[0] = chi*(lamVisc/CV->rhoInfinity);
  }

  zfs_log << "=================================================" << endl;
  zfs_log << "           INITIAL CONDITION SUMMARY" << endl;
  zfs_log << "=================================================" << endl;
  zfs_log << "Re = " << m_Re << endl;
  zfs_log << "Re0 = " << m_Re0 << endl;
  zfs_log << "Ma = " << m_Ma << endl;
  zfs_log << "TInfinity = " << PV->TInfinity << endl;
  zfs_log << "UInfinity = " << PV->UInfinity << endl;
  zfs_log << "VInfinity = " << PV->VInfinity << endl;
  zfs_log << "WInfinity = " << PV->WInfinity << endl;
  zfs_log << "PInfinity = " << PV->PInfinity << endl;
  zfs_log << "rhoInfinity = " << CV->rhoInfinity << endl;
  zfs_log << "rhoEInfinity = " << CV->rhoEInfinity << endl;
  zfs_log << "referenceTime = " << m_timeRef << endl;

  zfs_log << "=================================================" << endl;
  zfs_log << "           ZONAL CONNECTION SUMMARY" << endl;
  zfs_log << "=================================================" << endl;
  zfs_log << " noZonalBC: " <<m_noZonalBCMaps<< endl;
  
  for(ZFSInt id=0; id<m_noZonalBCMaps; id++){
    zfs_log<<" zonalBC: "<<id<< " zonal exchange cells: " << m_zonalBC[id].m_noCellsGlobalBC<<endl;
    // zfs_log<<" has STG: "<<m_zonalBC[id].m_hasSTG << endl;
    zfs_log<< "Number of Donor domains :"<< m_zonalBC[id].m_noGlobalSndDomains<<endl;
    zfs_log<< "Number of Receiver domains :"<<m_zonalBC[id].m_noGlobalRcvDomains<<endl;	
  }

  if(domainId() == 0 )
  {
    cout << "////////////////////////////////////////////////" << endl;
    cout << "////////// Initial Condition summary ///////////" << endl;
    cout << "////////////////////////////////////////////////" << endl;
    cout << "Re = " << m_Re << endl;
    cout << "Re0 = " << m_Re0 << endl;
    cout << "Ma = " << m_Ma << endl;
    cout << "TInfinity = " << PV->TInfinity << endl;
    cout << "UInfinity = " << PV->UInfinity << endl;
    cout << "VInfinity = " << PV->VInfinity << endl;
    cout << "WInfinity = " << PV->WInfinity << endl;
    cout << "Angle = "<<m_angle[ 0 ]<< "  "<<m_angle[ 1 ]<<endl;
    cout << "PInfinity = " << PV->PInfinity << endl;
    cout << "rhoInfinity = " << CV->rhoInfinity << endl;
    cout << "rhoEInfinity = " << CV->rhoEInfinity << endl;
    cout << "referenceTime = " << m_timeRef << endl;
    cout << " zonal = "<<m_zonal<< endl;

    if(m_zonal){
      cout << "////////////////////////////////////////////////" << endl;
      cout << "//////////Zonal Connections Summary/////////////" << endl;
      cout << "////////////////////////////////////////////////" << endl;
      cout << " noZonalBC: " <<m_noZonalBCMaps<< endl;
      cout << "////////////////////////////////////////////////" << endl;
  
      for(ZFSInt id=0; id<m_noZonalBCMaps; id++){
	cout<<" zonalBC: "<<id<< " zonal exchange cells: " << m_zonalBC[id].m_noCellsGlobalBC<<endl;
	// cout<<" has STG: "<<m_zonalBC[id].m_hasSTG << endl;
	cout<< "No. of Donor domains :"<< m_zonalBC[id].m_noGlobalSndDomains<<endl;
	cout<< "No. of Receiver domains :"<<m_zonalBC[id].m_noGlobalRcvDomains<<endl;
	cout << "////////////////////////////////////////////////" << endl;
      }
    }
  }

  if( !m_restart ) {
    // initialize random number generator
    srand(0);
    // inflow condition
    // ----------------
    switch (m_initialCondition)
    {
    case 0 :
    {
      //parallel inflow field
      for(ZFSId cellid=0; cellid<m_noStrctrdCells; cellid++) {
        //go through every cell
        m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
        for( ZFSId i = 0; i < nDim; i++ ) {
          m_cells->pvariables[PV->VV[i]][cellid] = PV->VVInfinity[i];
        }

        m_cells->pvariables[PV->P][cellid]= PV->PInfinity;

        if(m_rans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellid] = PV->ransInfinity[0];
        }
      }
      break;
    }
    case 43 :
    {
      //parallel inflow field with pressure peak in the middle of the domain
      for(ZFSId cellid=0; cellid<m_noStrctrdCells; cellid++) {
        //go through every cell
        m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
        for( ZFSId i = 0; i < nDim; i++ ) {
          m_cells->pvariables[PV->VV[i]][cellid] = F0;
        }

        m_cells->pvariables[PV->P][cellid]= PV->PInfinity;

        ZFSFloat radius = sqrt(POW2(m_cells->coordinates[0][cellid] - 0.5) + 
                               POW2(m_cells->coordinates[1][cellid] - 0.5));

        //impose pressure peak in the middle of the domain
        if(radius <= 0.05) {
          ZFSFloat pAmp = 0.005;
          ZFSFloat pressureSignal = sin(radius/0.05*PI)*pAmp + PV->PInfinity;
          m_cells->pvariables[PV->P][cellid]= pressureSignal;
        }
      }
      break;
    }
    case 333 :
    {
      //parallel inflow field
      for(ZFSId cellid=0; cellid<m_noStrctrdCells; cellid++) {
        //go through every cell
        m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
        for( ZFSId i = 0; i < nDim; i++ ) {
          m_cells->pvariables[PV->VV[i]][cellid] = PV->VVInfinity[i];
        }

        m_cells->pvariables[PV->P][cellid]= PV->PInfinity;

        //impose pressure peak in the middle of the domain
        if(m_cells->coordinates[0][cellid] > 0.4 && m_cells->coordinates[0][cellid] < 0.5) {
          ZFSFloat pAmp = 0.005;
          ZFSFloat xCoordinate = m_cells->coordinates[0][cellid] - 0.4;
          ZFSFloat pressureSignal = sin(xCoordinate/0.1*PI)*pAmp + PV->PInfinity;
          m_cells->pvariables[PV->P][cellid]= pressureSignal;
        }
      }
      break;
    }
    case 314 : //stagnating flow field
    {
      for(ZFSId cellid=0; cellid<m_noStrctrdCells; cellid++) {
        m_cells->pvariables[PV->RHO][cellid] = 1.0;
        for( ZFSId i = 0; i < nDim; i++ ) {
          m_cells->pvariables[PV->VV[i]][cellid] = F0;
        }

        m_cells->pvariables[PV->P][cellid]= PV->PInfinity;
      }
      cout << "I.C. stagnating flow field was applied! "<< endl;
      break;
    }
    case 315 : // Poiseuille flow
    {
      ZFSFloat x=F0, y=F0;// p=F0, T=F0, T0=F0;
      ZFSFloat y_max = F1;  // channel height

      for(ZFSId cellid=0; cellid<m_noStrctrdCells; cellid++) {
        x=m_cells->coordinates[0][cellid];
        y=m_cells->coordinates[1][cellid];

        for( ZFSId i = 0; i < nDim; i++ ) {
          m_cells->pvariables[PV->VV[i]][cellid] = F0;
        }
        // all velocity components are 0 except u, Poiseuille distribution:
        m_cells->pvariables[PV->VV[0]][cellid] = - (F3/F2)*PV->UInfinity*(POW2(y-y_max/F2)-POW2(y_max/F2))/POW2(y_max/F2);

        // the pressure is defined through the axial pressure gradient:
        m_cells->pvariables[ PV->P ][cellid] = PV->PInfinity - F3*(x+15.0)*zfsSUTHERLANDLAW(PV->TInfinity)*PV->UInfinity*POW2(F2/y_max)/m_Re0;

        // compressible Poiseuille Temperature distribution
        //T = T0 - (m_gamma-F1)*m_Pr * POW2(PV->UInfinity*F3/F2) * (F1B2*(F1+ pow( (y-y_max/F2)/(y_max/F2), F4 ) )-POW2((y-y_max/F2)/(y_max/F2)));
        m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;//m_gamma * m_cells->variables[ PV->P ][cellid] / T;//
      }

      break;
    }
    case 101 : //TAYLOR_GREEN_VORTEX
    {
      //domain boundaries are all 2*pi
      //rho=1.0;
      //u=A*SIN(x)*COS(y)*COS(z)
      //v=- A*COS(x)*SIN(y)*COS(z)
      //w=0.0
      //p=A*A*rho*( 1./(Ms*Ms*kappa) + 1./16.*(COS(2*x)*COS(2.*z)+ 2.*COS(2.*y) +2.*COS(2.*x) +COS(2*y)*COS(2.*z)))
      //Ms =0.1 maximum Mach number
      //A= speed magnitude set to 1.0
      ZFSInt cellId=0;
      ZFSFloat A=PV->UInfinity;
      ZFSFloat x=F0;
      ZFSFloat y=F0;
      ZFSFloat z=F0;
      for(ZFSId k=0; k<m_nCells[0]; k++){
        for(ZFSId j=0; j<m_nCells[1]; j++){
          for(ZFSId i=0; i<m_nCells[2]; i++){
            cellId=cellIndex(i,j,k);
            x=m_cells->coordinates[0][cellId];
            y=m_cells->coordinates[1][cellId];
            z=m_cells->coordinates[2][cellId];
            m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
            m_cells->pvariables[PV->VV[0]][cellId]=A*sin(x)*cos(y)*cos(z);
            m_cells->pvariables[PV->VV[1]][cellId]=-A*cos(x)*sin(y)*cos(z);
            m_cells->pvariables[PV->VV[2]][cellId]=0.0;
            m_cells->pvariables[PV->P][cellId]=PV->PInfinity+F1B16*(POW2(A)*CV->rhoInfinity)*(cos(2.0*x)+cos(2.0*y))*(2.0+cos(2.0*z));
          }
        }
      }
      break;
    }
    case 1234:
    {
      //laminar channel flow
      ZFSId cellId =0;
      m_channelHeight = 2.0;
      m_channelLength = 6.2831;
      // m_deltaP=32.0*zfsSUTHERLANDLAW(PV->TInfinity)*m_Ma*sqrt(PV->TInfinity)*m_channelLength/(m_Re0*m_channelHeight);
      m_deltaP = -12.0*PV->UInfinity*zfsSUTHERLANDLAW(PV->TInfinity)*m_channelLength/(POW2(m_channelHeight)*m_Re0);

      m_channelPresInlet=PV->PInfinity;
      m_channelPresOutlet=PV->PInfinity+m_deltaP;
      ZFSFloat u=PV->UInfinity; //POW2(m_channelHeight)*m_deltaP*m_Re0*m_referenceLength*m_referenceLength/m_channelLength;
      for(ZFSInt k =m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
        for(ZFSInt j =m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
          for(ZFSInt i =m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
            cellId=cellIndex(i,j,k);

            //channel height is in j direction
            //so prescribe a mean profile in u(y)
            ZFSFloat x = m_cells->coordinates[0][cellId];
            pressureCH=m_channelPresInlet+(m_channelPresOutlet-m_channelPresInlet)*x/m_channelLength;
            //density:
            ZFSFloat rho = pressureCH*m_gamma/PV->TInfinity;
            m_cells->pvariables[PV->RHO][cellId]=rho;
            m_cells->pvariables[PV->U][cellId]=u;
            m_cells->pvariables[PV->V][cellId]=F0;
            m_cells->pvariables[PV->W][cellId]=F0;
            m_cells->pvariables[PV->P][cellId]=pressureCH;
          }
        }
      }

      break;
    }
    case 1233:
    {
      //turbulent channel with perturbations
      //calculate the Pressure loss;

      //for the law of the wall
      const ZFSFloat C1      = m_channelC1;
      const ZFSFloat C2      = m_channelC2;
      const ZFSFloat C3      = m_channelC3;
      const ZFSFloat C4      = m_channelC4;
      ZFSFloat xINIT=m_channelInflowPlaneCoordinate;
      ZFSId cellId=0;
      ZFSFloat yplus=F0;
      ZFSFloat uTau= m_ReTau*m_Ma*sqrt(PV->TInfinity)/m_Re;
      ZFSFloat prefactor=m_ReTau;//uTau*m_Re0/zfsSUTHERLANDLAW(PV->TInfinity);

      m_deltaP = -CV->rhoInfinity*POW2(uTau)*F2*(m_channelLength)/m_channelHeight;

      zfs_log << "uTau: " << uTau << " channelLength: " << m_channelLength << endl;
      //mean velocity profile
      m_channelPresInlet = PV->PInfinity;
      m_channelPresOutlet = PV->PInfinity - CV->rhoInfinity*POW2(uTau)*F2*(xINIT+m_channelLength)/m_channelHeight;
      ZFSFloat deltaP = m_channelPresInlet - m_channelPresOutlet;
      ZFSFloat cfTheo = 2*POW2(uTau/PV->UInfinity);
      ZFSFloat cdTheo = cfTheo*2.0*m_channelWidth*m_channelLength;
      zfs_log << "deltaP: " << deltaP << " cfTheoretisch: " << cfTheo << " cdTheo: " << cdTheo << endl;


      for(ZFSInt k =0; k<m_nCells[0]; k++) {
        for(ZFSInt j =0; j<m_nCells[1]; j++) {
          for(ZFSInt i =0; i<m_nCells[2]; i++) {
            //channel height is in j direction
            //so prescribe a mean profile in u(y)
            cellId=cellIndex(i,j,k);
            //channel starts at x=0 or we need to prescribe an offset
            ZFSFloat y = m_cells->coordinates[1][cellId];
            ZFSFloat x = m_cells->coordinates[0][cellId];
            pressureCH=m_channelPresInlet+((m_channelPresOutlet-m_channelPresInlet)*(x-xINIT)/m_channelLength);

            ZFSFloat rho = pressureCH*m_gamma/PV->TInfinity;
            ZFSFloat velFactor = 2.0;
            if(y<m_channelHeight/2) {
              yplus=prefactor*y;
            } else {
              yplus=prefactor*(m_channelHeight-y);
            }

            //C1 etc are defined in zfsconstants.h
            if(yplus<=5.0) {
              m_cells->pvariables[PV->U][cellId]=0.5*uTau*yplus*velFactor;
            } else if (yplus<=30 && yplus>5.0) {
              m_cells->pvariables[PV->U][cellId]=0.5*uTau*(C1*log(yplus)+C2)*velFactor;
            } else if(yplus>30) {
              m_cells->pvariables[PV->U][cellId]=0.5*uTau*(C3*log(yplus)+C4)*velFactor;
            }

            m_cells->pvariables[PV->RHO][cellId]=rho;
            m_cells->pvariables[PV->V][cellId]=F0;
            m_cells->pvariables[PV->W][cellId]=F0;
            m_cells->pvariables[PV->P][cellId]=pressureCH;
          }
        }
      }
      //create the fluctuations:


      //this way of creating the
      //fluctuations is only possible for
      //quite small mounts of cells
      //if too large this approach will not work anymore
      //125000000 is a random integer number which has to be tested
      //else the approach with white noise will be employed
      if(m_totalGridBlockCells[m_inputBlockId]<=0) {
        fftw_complex *uPhysField, *vPhysField, *wPhysField;

        //we need the total number of points
        ZFSInt lx=m_totalGridBlockDim[m_inputBlockId][2],ly=m_totalGridBlockDim[m_inputBlockId][1],lz=m_totalGridBlockDim[m_inputBlockId][0];

        // field of velocities from positve frequencies
        uPhysField = (fftw_complex*) fftw_malloc(lx*ly*lz * sizeof(fftw_complex));
        vPhysField = (fftw_complex*) fftw_malloc(lx*ly*lz * sizeof(fftw_complex));
        wPhysField = (fftw_complex*) fftw_malloc(lx*ly*lz * sizeof(fftw_complex));

        //no real parallel computation is done
        //should be implemented later !!!!

        if(noDomains()>1) {
          ZFSFloatScratchSpace sendRecvBufferU(lx*ly*lz, __CALLING_FUNCTION__, "sendRecvBufferU");
          ZFSFloatScratchSpace sendRecvBufferV(lx*ly*lz, __CALLING_FUNCTION__, "sendRecvBufferV");
          ZFSFloatScratchSpace sendRecvBufferW(lx*ly*lz, __CALLING_FUNCTION__, "sendRecvBufferW");
          if(domainId()==0) {
            ZFSInt m_noPeakModes=100;
            initFFTW(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes);
            //copy values into the sendRCVbuffer

            for(ZFSInt id=0; id<lz*ly*lz; id++) {
              sendRecvBufferU[id]=uPhysField[id][0];
              sendRecvBufferV[id]=vPhysField[id][0];
              sendRecvBufferW[id]=wPhysField[id][0];
            }
          }
          MPI_Bcast(&sendRecvBufferU[0], lx*ly*lz, MPI_DOUBLE, 0, m_zfsStrctrdComm);
          MPI_Bcast(&sendRecvBufferV[0], lx*ly*lz, MPI_DOUBLE, 0, m_zfsStrctrdComm);
          MPI_Bcast(&sendRecvBufferW[0], lx*ly*lz, MPI_DOUBLE, 0, m_zfsStrctrdComm);

          if(domainId()!=0) {
            for(ZFSInt id=0; id<lz*ly*lz; id++) {
              uPhysField[id][0]=sendRecvBufferU[id];
              vPhysField[id][0]=sendRecvBufferV[id];
              wPhysField[id][0]=sendRecvBufferW[id];
            }
          }
        } else {
          // create velocity field
          ZFSInt m_noPeakModes=100;
          initFFTW(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes);
        }
        //now we need to distribute the
      } else {
        ZFSFloat amp=0.15;
        for(ZFSInt k =m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
          for(ZFSInt j =m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
            for(ZFSInt i =m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
              cellId=cellIndex(i,j,k);
              m_cells->pvariables[PV->V][cellId]+=amp * 2.0 * (0.5 - (ZFSFloat) (1.0*rand()/(RAND_MAX+1.0)))* m_cells->pvariables[PV->U][cellId];
              m_cells->pvariables[PV->W][cellId]+=amp * 2.0 * (0.5 - (ZFSFloat) (1.0*rand()/(RAND_MAX+1.0)))* m_cells->pvariables[PV->U][cellId];
              m_cells->pvariables[PV->U][cellId]+=amp * 2.0 * (0.5 - (ZFSFloat) (1.0*rand()/(RAND_MAX+1.0)))* m_cells->pvariables[PV->U][cellId];
            }
          }
        }
      }

      break;
    }
    case 1236:
    {
      // pipe with perturbations
      //calculate the Pressure loss;
      //for the law of the wall
      //const ZFSFloat C1      = m_channelC1;
      //const ZFSFloat C2      = m_channelC2;
      const ZFSFloat C3      = m_channelC3;
      const ZFSFloat C4      = m_channelC4;
      ZFSId cellId=0;
      ZFSFloat uTau= m_ReTau*m_Ma*sqrt(PV->TInfinity)/m_Re;
      //ZFSFloat prefactor=m_ReTau;//uTau*m_Re0/zfsSUTHERLANDLAW(PV->TInfinity);

      m_deltaP = -4.0*CV->rhoInfinity*POW2(uTau)*(m_channelLength)/m_channelHeight;
      m_channelPresInlet=PV->PInfinity;
      m_channelPresOutlet=PV->PInfinity+m_deltaP;
      const ZFSFloat bulkVel = uTau *(C3*log(m_ReTau/2)+C4);

      zfs_log << "=========== Turb. Pipe Flow Inital Condition Summary =========== " <<endl;
      zfs_log << "-->Turbulent pipe flow deltaP: " << m_deltaP << endl;
      zfs_log << "-->pipe friciton velocity: " << uTau << endl;
      zfs_log << "-->pipe pressure inflow: " << m_channelPresInlet << endl;
      zfs_log << "-->pipe pressure outflow: " << m_channelPresOutlet << endl;
      zfs_log << "--> bulk velocity (u_max)" << bulkVel << endl;
      zfs_log << "=========== Turb. Pipe Flow Initial Condition Summary Finished =========== " <<endl;

      for(ZFSInt k =0; k<m_nCells[0]; k++) {
        for(ZFSInt j =0; j<m_nCells[1]; j++) {
          for(ZFSInt i =0; i<m_nCells[2]; i++) {
            //IMPORTANT PIPE LENGTH HAS TO GO INTO X-DIRECTION
            //so prescribe a mean profile in u(r)
            //centerline is assumed at y=0.0,z=0.0

            cellId=cellIndex(i,j,k);
            //determine the radius
            ZFSFloat r = sqrt(POW2(m_cells->coordinates[1][cellId])+POW2(m_cells->coordinates[2][cellId]));

            ZFSFloat x = m_cells->coordinates[0][cellId];
            //determine the pressure drop
            pressureCH = m_deltaP/m_channelLength*(x - m_channelInflowPlaneCoordinate) + PV->PInfinity;
            ZFSFloat rho = pressureCH*m_gamma/PV->TInfinity;
            //important: viscous sublayer is not build explicitly. The flow has to build it by itself
            ZFSFloat vel= bulkVel + C3*log(F1-min(((F2*r)/m_channelHeight),0.999999999))*uTau;
            m_cells->pvariables[PV->RHO][cellId]=rho;
            m_cells->pvariables[PV->U][cellId]=vel;
            m_cells->pvariables[PV->V][cellId]=F0;
            m_cells->pvariables[PV->W][cellId]=F0;
            m_cells->pvariables[PV->P][cellId]=pressureCH;
          }
        }
      }

      //create the fluctuations:
      //this way of creating the
      //fluctuations is only possible for
      //quite small mounts of cells
      //if too large this approach will not work anymore
      //125000000 is a random integer number which has to be tested
      //else the approach with white noise will be employed
      if(m_totalGridBlockCells[m_inputBlockId]<=0) {
        fftw_complex *uPhysField, *vPhysField, *wPhysField;

        //we need the total number of points
        ZFSInt lx=m_totalGridBlockDim[m_inputBlockId][2],ly=m_totalGridBlockDim[m_inputBlockId][1],lz=m_totalGridBlockDim[m_inputBlockId][0];

        // field of velocities from positve frequencies
        uPhysField = (fftw_complex*) fftw_malloc(lx*ly*lz * sizeof(fftw_complex));
        vPhysField = (fftw_complex*) fftw_malloc(lx*ly*lz * sizeof(fftw_complex));
        wPhysField = (fftw_complex*) fftw_malloc(lx*ly*lz * sizeof(fftw_complex));

        //no real parallel computation is done
        //should be implemented later !!!!

        if(noDomains()>1) {
          ZFSFloatScratchSpace sendRecvBufferU(lx*ly*lz, __CALLING_FUNCTION__, "sendRecvBufferU");
          ZFSFloatScratchSpace sendRecvBufferV(lx*ly*lz, __CALLING_FUNCTION__, "sendRecvBufferV");
          ZFSFloatScratchSpace sendRecvBufferW(lx*ly*lz, __CALLING_FUNCTION__, "sendRecvBufferW");
          if(domainId()==0) {
            ZFSInt m_noPeakModes=100;
            initFFTW(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes);
            //copy values into the sendRCVbuffer

            for(ZFSInt id=0; id<lz*ly*lz; id++) {
              sendRecvBufferU[id]=uPhysField[id][0];
              sendRecvBufferV[id]=vPhysField[id][0];
              sendRecvBufferW[id]=wPhysField[id][0];
            }
          }
          MPI_Bcast(&sendRecvBufferU[0], lx*ly*lz, MPI_DOUBLE, 0, m_zfsStrctrdComm);
          MPI_Bcast(&sendRecvBufferV[0], lx*ly*lz, MPI_DOUBLE, 0, m_zfsStrctrdComm);
          MPI_Bcast(&sendRecvBufferW[0], lx*ly*lz, MPI_DOUBLE, 0, m_zfsStrctrdComm);

          if(domainId()!=0) {
            for(ZFSInt id=0; id<lz*ly*lz; id++) {
              uPhysField[id][0]=sendRecvBufferU[id];
              vPhysField[id][0]=sendRecvBufferV[id];
              wPhysField[id][0]=sendRecvBufferW[id];
            }
          }
        } else {
          // create velocity field
          ZFSInt m_noPeakModes=100;
          initFFTW(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes);
        }
        //now we need to distribute the
      } else {
        ZFSFloat amp=0.15;
        for(ZFSInt k =m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
          for(ZFSInt j =m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
            for(ZFSInt i =m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
              cellId=cellIndex(i,j,k);
              m_cells->pvariables[PV->V][cellId]+=amp * 2.0 * (0.5 - (ZFSFloat) (1.0*rand()/(RAND_MAX+1.0)))* m_cells->pvariables[PV->U][cellId];
              m_cells->pvariables[PV->W][cellId]+=amp * 2.0 * (0.5 - (ZFSFloat) (1.0*rand()/(RAND_MAX+1.0)))* m_cells->pvariables[PV->U][cellId];
              m_cells->pvariables[PV->U][cellId]+=amp * 2.0 * (0.5 - (ZFSFloat) (1.0*rand()/(RAND_MAX+1.0)))* m_cells->pvariables[PV->U][cellId];
            }
          }
        }
      }

      break;
    }
    case 79092:{
      //approximate mean turbulent boundary layer profile
      const ZFSFloat epss = 1e-10;
      const ZFSFloat reTheta=1000.000;
      const ZFSFloat theta = 1.0;
      const ZFSFloat delta0 = 72.0/7.0 * theta;
      const ZFSFloat kappa=0.4;
      const ZFSFloat C1=3.573244189003983;
      const ZFSFloat PI1 = 0.55;
      const ZFSFloat cf= 0.024/pow(reTheta, 0.25);
      const ZFSFloat uTau=sqrt(cf/2.0)*PV->UInfinity;
      const ZFSFloat nu = zfsSUTHERLANDLAW(PV->TInfinity);
      for(ZFSId k=0; k<m_nCells[0]; k++){
        for(ZFSId j=0; j<m_nCells[1]; j++){
          for(ZFSId i=0; i<m_nCells[2]; i++){
            const ZFSId cellId=cellIndex(i,j,k);
            m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
            const ZFSFloat y=m_cells->coordinates[1][cellId];
            ZFSFloat yPlus= zfsMAX(uTau*y*m_Re0/nu, F0);

            if(y>delta0){
              m_cells->pvariables[PV->U][cellId]= PV->UInfinity;
            } else if(yPlus <= 10.0){
              m_cells->pvariables[PV->U][cellId]= yPlus*uTau;
            } else if ( yPlus <=30 && yPlus>10.0){
              m_cells->pvariables[PV->U][cellId]=uTau*((F1/kappa)*log(max(yPlus,epss))+C1);
            }else if( yPlus >30.0 && y<=delta0){
              m_cells->pvariables[PV->U][cellId]= uTau * ((F1/kappa)*log(max(yPlus,epss)) + C1 + 2*(PI1/kappa)*(3*y*y - 2*y*y*y));
            }

            m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
            m_cells->pvariables[PV->W][cellId]=PV->WInfinity;
            m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
            m_cells->pvariables[PV->P][cellId]=PV->PInfinity;

            if(m_rans) {
              m_cells->pvariables[PV->RANS_VAR[0]][cellId] = PV->ransInfinity[0];
            }
          }
        }
      }

      break;
    }

    case 79091: //Turbulent plate
    {
      const ZFSFloat epss = 1e-10;
      const ZFSFloat reTheta = 1000.0;
      const ZFSFloat theta = 1.0;
      const ZFSFloat delta0 = 72.0/7.0 * theta;
      const ZFSFloat K = 0.4;
      const ZFSFloat C1 = 3.573244189003983;  //With coles
      const ZFSFloat PI1 = 0.55;
      const ZFSFloat cf = 0.024 / pow(reTheta,0.25);

      for(ZFSId k=0; k<m_nCells[0]; k++) {
        for(ZFSId j=0; j<m_nCells[1]; j++) {
          for(ZFSId i=0; i<m_nCells[2]; i++) {
            const ZFSId cellId=cellIndex(i,j,k);
            const ZFSFloat mu = zfsSUTHERLANDLAW(PV->TInfinity);
            const ZFSFloat utau = sqrt(cf/2.0) * m_Ma * sqrt(PV->TInfinity);
            const ZFSFloat yplus = m_cells->coordinates[1][cellId] * sqrt(cf/2.) * CV->rhoUInfinity/mu*m_Re0;
            const ZFSFloat eta = m_cells->coordinates[1][cellId] / delta0; //y/delta

            //1-7th profile
            //log-law + wake
            if(m_cells->coordinates[1][cellId] > delta0) {
              m_cells->pvariables[PV->U][cellId] = PV->UInfinity; //Outside BL
            } else if(yplus < 10) {
              m_cells->pvariables[PV->U][cellId] = utau * yplus;
            } else {
              m_cells->pvariables[PV->U][cellId] = zfsMIN(utau * ((1./K)*log(max(yplus,epss)) + C1
                                                                     + 2*PI1/K*(3*eta*eta - 2*eta*eta*eta)) , PV->UInfinity);
            }

            m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
            m_cells->pvariables[PV->W][cellId]=PV->WInfinity;
            m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
            m_cells->pvariables[PV->P][cellId] = PV->PInfinity;

            if(m_rans) {
              m_cells->pvariables[PV->RANS_VAR[0]][cellId] = PV->ransInfinity[0];
            }
          }
        }
      }

      break;
    }
    case 11 : //point source in the middle
    {
      //contain an initial perturbation in the middle
      for(ZFSId cellid=0; cellid<m_noStrctrdCells; cellid++) {
        //go through every cell
        m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
        for( ZFSId i = 0; i < nDim; i++ ) {
          m_cells->pvariables[PV->VV[i]][cellid] = PV->VVInfinity[i];
        }

        ZFSFloat amp= 0.0001;
        ZFSFloat x = m_cells->coordinates[0][cellid];
        ZFSFloat y = m_cells->coordinates[1][cellid];
        ZFSFloat z = m_cells->coordinates[2][cellid];
        ZFSFloat r = 0.025;
        ZFSFloat a1= POW2(x-0.5);
        ZFSFloat a2= POW2(y-0.5);
        ZFSFloat a3= POW2(z-0.5);
        ZFSFloat disturb = amp*exp(-(a1+a2+a3)/POW2(r)/2.0);
        m_cells->pvariables[PV->P][cellid]= PV->PInfinity+disturb;
      }
      break;
    }
    /** TESTCASE from C. Bogey, C. Bailly, Three-dimensional non-reflective boundary conditions
     *  for acoustic simulations: far field formulation and validation test cases
     *  Acta acustica united with acustica Volume 88 (2002), 463-471
     */
    case 111 :
    {
      ZFSFloat amp=0.01;
      ZFSFloat alpha=log(2)/9;
      for(ZFSId cellid=0; cellid<m_noStrctrdCells; cellid++) {
        //go through every cell
        m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
        ZFSFloat x = m_cells->coordinates[0][cellid];
        ZFSFloat y = m_cells->coordinates[1][cellid];
        ZFSFloat z = m_cells->coordinates[2][cellid];
        ZFSFloat fluc = amp*exp(-alpha*(POW2(x)+POW2(y)+POW2(z)));
        m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity+fluc;
        m_cells->pvariables[PV->P][cellid]=(PV->PInfinity+fluc);
        for( ZFSId i = 0; i < nDim; i++ ) {
          m_cells->pvariables[PV->VV[i]][cellid] = PV->VVInfinity[i];
        }
      }
      break;
    }
    /** TESTCASE from C. Bogey, C. Bailly, Three-dimensional non-reflective boundary conditions
     *  for acoustic simulations: far field formulation and validation test cases
     *  Acta acustica united with acustica Volume 88 (2002), 463-471 the vortex
     */
    case 112:
    {
      ZFSFloat amp=0.003;
      ZFSFloat b=0.025;
      ZFSFloat alpha=log(2)/POW2(b);
      ZFSFloat r0=0.1; //changed!!!
      ZFSFloat r=F0;
      ZFSFloat phi=F0;
      ZFSFloat v=F0;
      for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
        //go through every cell
        m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
        m_cells->pvariables[PV->P][cellId]= PV->PInfinity;
        m_cells->pvariables[PV->U][cellId]=PV->UInfinity;
        m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
        m_cells->pvariables[PV->W][cellId]=PV->WInfinity;
        ZFSFloat x = m_cells->coordinates[0][cellId]-0.5;
        ZFSFloat y = m_cells->coordinates[1][cellId]-0.5;
        ZFSFloat z = m_cells->coordinates[2][cellId]-0.25;
        r=sqrt(POW2(y)+POW2(z));
        phi=atan(z/y);
        m_cells->pvariables[PV->U][cellId]+=amp*(r0/r)*(r-r0)*exp(-1.0*alpha*(POW2(x)+POW2(r-r0)));
        v= -1.0*amp*(r0/r)*x*exp(-1.0*alpha*(POW2(x)+POW2(r-r0)));
        m_cells->pvariables[PV->V][cellId]+=v*r*sin(phi);
        m_cells->pvariables[PV->W][cellId]+=v*r*cos(phi);
      }
      saveOutput(true);
      break;
    }
    case 113:
    {
      ZFSFloat amp=0.003;
      ZFSFloat b=0.25/4;
      ZFSFloat alpha=log(2)/POW2(b);
      ZFSFloat r0=0.25; //changed!!!
      ZFSFloat r=F0;
      ZFSFloat phi=F0;
      ZFSFloat v=F0;
      for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
        //go through every cell
        m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
        m_cells->pvariables[PV->P][cellId]= PV->PInfinity;
        m_cells->pvariables[PV->U][cellId]=PV->UInfinity;
        m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
        m_cells->pvariables[PV->W][cellId]=PV->WInfinity;
        ZFSFloat x = m_cells->coordinates[0][cellId]+0.5;
        ZFSFloat y = m_cells->coordinates[1][cellId];
        ZFSFloat z = m_cells->coordinates[2][cellId]-0.25;
        r=sqrt(POW2(y)+POW2(z));
        phi=atan(z/y);
        m_cells->pvariables[PV->U][cellId]+=amp*(r0/r)*(r-r0)*exp(-1.0*alpha*(POW2(x)+POW2(r-r0)));
        v= -1.0*amp*(r0/r)*x*exp(-1.0*alpha*(POW2(x)+POW2(r-r0)));
        m_cells->pvariables[PV->V][cellId]+=v*r*sin(phi);
        m_cells->pvariables[PV->W][cellId]+=v*r*cos(phi);
      }
      break;
    }
    case 2: //shear flow is prescribed
    {
      ZFSInt cellId=0;
      ZFSFloat x=F0;
      for(ZFSId k=0; k<m_nCells[0]; k++) {
        for(ZFSId j=0; j<m_nCells[1]; j++) {
          for(ZFSId i=0; i<m_nCells[2]; i++) {
            cellId=cellIndex(i,j,k);
            x=m_cells->coordinates[0][cellId];
            m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
            if(x<0.5) {
              m_cells->pvariables[PV->U][cellId] = F0;
              m_cells->pvariables[PV->V][cellId] = 0.5*PV->UInfinity;
              m_cells->pvariables[PV->W][cellId] = F0;
            } else {
              m_cells->pvariables[PV->U][cellId] = F0;
              m_cells->pvariables[PV->V][cellId] = PV->UInfinity;
              m_cells->pvariables[PV->W][cellId] = F0;
            }
            m_cells->pvariables[PV->P][cellId]=  PV->PInfinity;
          }
        }
      }
      break;
    }
    case 900:{ //jet Freund
      for(ZFSId k=0; k<m_nCells[0]; k++) {
        for(ZFSId j=0; j<m_nCells[1]; j++) {
          for(ZFSId i=0; i<m_nCells[2]; i++) {
            ZFSId cellId = cellIndex(i,j,k);
            ZFSFloat r=sqrt(POW2(m_cells->coordinates[0][cellId])+POW2(m_cells->coordinates[1][cellId])+POW2(m_cells->coordinates[2][cellId]));
            ZFSFloat u= F1B2*(F1-tanh(12.5*(fabs(r/0.5)-fabs(0.5/r))))*PV->VVInfinity[0];
            m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
            m_cells->pvariables[PV->U][cellId]=u;
            m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
            m_cells->pvariables[PV->W][cellId]=PV->WInfinity;
            m_cells->pvariables[PV->P][cellId]=PV->PInfinity;
          }
        }
      }
      break;
    }
    case 4001:{//test periodic rotation boundary conditions
      for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; ++k){
        for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; ++j){
          for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; ++i){
            ZFSId cellId = cellIndex(i,j,k);
            //ZFSFloat x= m_cells->coordinates[0][cellId];
            ZFSFloat y= m_cells->coordinates[1][cellId];
            ZFSFloat z= m_cells->coordinates[2][cellId];
            ZFSFloat phi=atan2(y,z);
            ZFSFloat r=sqrt(POW2(y)+POW2(z));
            ZFSFloat rmax = 10.0;
            m_cells->pvariables[PV->RHO][cellId]=CV->rhoInfinity;
            m_cells->pvariables[PV->U][cellId]=PV->UInfinity;
            m_cells->pvariables[PV->V][cellId]=-(r/rmax)*cos(phi)* 0.1* PV->UInfinity;
            m_cells->pvariables[PV->W][cellId]=(r/rmax)*sin(phi) * 0.1* PV->UInfinity;
            m_cells->pvariables[PV->P][cellId]=PV->PInfinity;
          }
        }
      }
      break;
    }
    case 42: {
      for(ZFSId cellId = 0; cellId < m_noStrctrdCells; cellId++) {
        m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
        m_cells->pvariables[PV->V][cellId] = F0;
        m_cells->pvariables[PV->W][cellId] = F0;

        if(m_cells->coordinates[1][cellId] < 0.5) {
          m_cells->pvariables[PV->U][cellId] = 1.001*PV->UInfinity;
        } else {
          m_cells->pvariables[PV->U][cellId] = 0.999*PV->UInfinity;
        }

        m_cells->pvariables[PV->P][cellId]=  PV->PInfinity;
      }
      break;
    }
    case 44: {
      ZFSFloat amp=0.003;
      ZFSFloat b=0.25/4;
      ZFSFloat alpha=log(2)/POW2(b);
      ZFSFloat r0=0.25; //changed!!!
      ZFSFloat r=F0;
      ZFSFloat phi=F0;
      ZFSFloat v=F0;
      for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
        //go through every cell
        m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
        m_cells->pvariables[PV->P][cellId]= PV->PInfinity;
        m_cells->pvariables[PV->U][cellId]=PV->UInfinity;
        m_cells->pvariables[PV->V][cellId]=PV->VInfinity;
        m_cells->pvariables[PV->W][cellId]=PV->WInfinity;
        ZFSFloat x = m_cells->coordinates[0][cellId]-0.61;
        ZFSFloat y = m_cells->coordinates[1][cellId]-0.55;
        ZFSFloat z = m_cells->coordinates[2][cellId]-0.43;
        r=sqrt(POW2(y)+POW2(z));
        phi=atan(z/y);
        m_cells->pvariables[PV->U][cellId]+=amp*(r0/zfsMAX(r,0.00001))*(r-r0)*exp(-1.0*alpha*(POW2(x)+POW2(r-r0)));
        v=-1.0*amp*(r0/zfsMAX(r,0.00001))*x*exp(-1.0*alpha*(POW2(x)+POW2(r-r0)));
        m_cells->pvariables[PV->V][cellId]+=v*r*sin(phi);
        m_cells->pvariables[PV->W][cellId]+=v*r*cos(phi);
      }
      break;
    }
    default:
    {
      //put the parallel flow field input in here
      //force output that no specific initial condition was chosen
      zfsTerm(1, __CALLING_FUNCTION__, "No (correct) initial Condition is given!");
      break;
    }
    }
  }
}


/* Compute all the interpolation coefficients
   necessary for the line/field output
*/
void ZFSStrctrdBlck3D::initLineInterpolation(){
  TRACE();
  if(!m_movingGrid) {
    m_pointInterpolation = new ZFSStrctrdInterpolation<3>(m_nCells, m_cells->coordinates, m_cells->pvariables, m_zfsStrctrdComm);
    //allocate domains points of the lines
    m_pointInterpolation->prepareInterpolation(m_noFieldPointsTotal, m_pointCoordinates, m_hasPartnerLocal);
    MPI_Allreduce(m_hasPartnerLocal, m_hasPartnerGlobal, m_noFieldPointsTotal, MPI_INT, MPI_SUM, m_zfsStrctrdComm);
  }
}

/* Manually correct errors made by the restart interpolation
   In case some cells did not get correct values in
   the interpolation process (due to no matching donor cells)
   Put your routines that will correct these cells here.
*/
void ZFSStrctrdBlck3D::manualInterpolationCorrection(){

  /*! \page propertyPage1
    \section interpolationCorrection
    <code>ZFSInt ZFSStrctrdBlck::m_interpolationCorrection </code>\n
    default = <code> 0 </code>\n \n
    Trigger a manual correction of the interpolation\n
    from a donor grid.\n
    Possible values are:\n
    <ul>
    <li>true/false</li>
    </ul>
    Keywords: <i>LIMITER, STRCTRD</i>
  */
  ZFSBool interpolationCorrection = false;
  if(ZFSContext::propertyExists("interpolationCorrection", m_blockId )){
    interpolationCorrection = *(ZFSContext::getProperty("interpolationCorrection", m_blockId, __CALLING_FUNCTION__,&interpolationCorrection)->asInt(0));
  }
  if(m_nOffsetCells[1] == 0 && m_rans && interpolationCorrection) {
    cout << "Correcting interpolation" << endl;
    for(ZFSId i=0; i<m_nCells[2]; i++) {
      for(ZFSId k=0; k<m_nCells[0]; k++) {
        ZFSId cellIdA1 = cellIndex(i,2,k);
        ZFSId cellIdA2 = cellIndex(i,3,k);
        ZFSId cellIdA3 = cellIndex(i,4,k);

        for(ZFSId var=0; var<PV->noVariables; var++) {
          m_cells->pvariables[var][cellIdA1] = m_cells->pvariables[var][cellIdA2] +
            (m_cells->coordinates[1][cellIdA1]-m_cells->coordinates[1][cellIdA2])/(m_cells->coordinates[1][cellIdA2]-m_cells->coordinates[1][cellIdA3])*
            (m_cells->variables[var][cellIdA2]-m_cells->pvariables[var][cellIdA3]);
        }
      }
    }
  }
}


/* brief returns a normal distributed random-
 * number with mu=mean and sigma=standard deviation
 *
 */
ZFSFloat ZFSStrctrdBlck3D::randnormal(ZFSFloat mu, ZFSFloat sigma) {
  TRACE();
  static bool deviateAvailable=false;        //        flag
  static float storedDeviate;                        //        deviate from previous calculation
  double polar, rsquared, var1, var2;

  // If no deviate has been stored, the polar Box-Muller transformation is
  // performed, producing two independent normally-distributed random
  // deviates.  One is stored for the next round, and one is returned.
  if (!deviateAvailable) {
    // choose pairs of uniformly distributed deviates, discarding those
    // that don't fall within the unit circle
    do {
      var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
      var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
      rsquared=var1*var1+var2*var2;
    }
    while ( rsquared>=1.0 || approx(rsquared, F0, m_eps));
      
    // calculate polar tranformation for each deviate
    polar=sqrt(-2.0*log(rsquared)/rsquared);
      
    // store first deviate and set flag
    storedDeviate=var1*polar;
    deviateAvailable=true;
      
    // return second deviate
    return var2*polar*sigma + mu;

    // If a deviate is available from a previous call to this function, it is
    // eturned, and the flag is set to false.
  } else {
    deviateAvailable=false;
    return storedDeviate*sigma + mu;
  }
}

void ZFSStrctrdBlck3D::nonReflectingBC()
{
  TRACE();
  m_strctrdBndryCnd->applyNonReflectingBC();
}

void ZFSStrctrdBlck3D::applyBoundaryCondition()
{
  TRACE();
  // treat Dirichlet and Neumann BC in one go!!!
  m_strctrdBndryCnd->applyDirichletNeumannBC();
}


void ZFSStrctrdBlck3D::MusclRANS(){
  m_ransBlck->Muscl();
}


void ZFSStrctrdBlck3D::MusclMinModLimiter()
{
  TRACE();
  //stencil identifier
  const ZFSId IJK[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};
  // switch for order
  const ZFSInt sword = 2;

  //reduce to onedimensional arrays
  ZFSFloat* __restrict x=&m_cells->coordinates[0][0];
  ZFSFloat* __restrict y=&m_cells->coordinates[1][0];
  ZFSFloat* __restrict z=&m_cells->coordinates[2][0];
  ZFSFloat* RESTRICT flux = ALIGNED_F(m_cells->flux);

  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSFloat *const RESTRICT cellVariables= ALIGNED_F(m_cells->pvariables[0]);

  for(ZFSId dim=0; dim<nDim; ++dim) {
    for(ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers; ++k) {
      for(ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers; ++j) {
        for(ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers; ++i) {
          //cell ids
          const ZFSId I=cellIndex(i,j,k);
          const ZFSId IP1=I+IJK[dim];
          const ZFSId IM1=I-IJK[dim];
          const ZFSId IP2=I+2*IJK[dim];

          //distances q_i+1 - q_i
          const ZFSFloat DS=sqrt(POW2(x[IP1]-x[I]) + POW2(y[IP1]-y[I]) + POW2(z[IP1]-z[I]));
          //distances q_i - q_i-1
          const ZFSFloat DSM1=sqrt(POW2(x[I]-x[IM1]) + POW2(y[I]-y[IM1]) + POW2(z[I]-z[IM1]));
          const ZFSFloat DSP1=sqrt(POW2(x[IP2]-x[IP1]) + POW2(y[IP2]-y[IP1]) + POW2(z[IP2]-z[IP1]));
          const ZFSFloat DSP=DS/POW2(DSP1+DS);
          const ZFSFloat DSM=DS/POW2(DSM1+DS);

          if (sword == 2) {
            for(ZFSId var=0; var<PV->noVariables; ++var) {
              const ZFSUint offset = var*noCells;
              const ZFSFloat *const RESTRICT vars = ALIGNED_F(cellVariables+offset);

              const ZFSFloat DQ   = vars[IP1]-vars[I];
              const ZFSFloat DQP1 = vars[IP2]-vars[IP1];
              const ZFSFloat DQM1 = vars[I]-vars[IM1];

              const ZFSFloat ri = DQM1/DQ;
              const ZFSFloat rip = DQ/DQP1;

              const ZFSFloat phii = zfsMAX(F0,zfsMIN(F1,ri));
              const ZFSFloat phiip = zfsMAX(F0,zfsMIN(F1,rip));

              m_QLeft[var]  = vars[I]   + ( DQ*DSM1 + DQM1*DS )*DSM*phii;
              m_QRight[var] = vars[IP1] - ( DQP1*DS + DQ*DSP1 )*DSP*phiip;
            }
          }

          AusmLES(m_QLeft, m_QRight, dim, I);
        }
      }
    }

    //FLUX BALANCE
    for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; ++k) {
      for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; ++j) {
        for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; ++i) {
          const ZFSId I=cellIndex(i,j,k);
          const ZFSId IM1=I-IJK[dim];

          for(ZFSId v=0; v<CV->noVariables; ++v) {
            m_cells->rightHandSide[v][I]+=flux[IM1+noCells*v]-flux[I+noCells*v];
          }
        }
      }
    }
  }
}

//Muscl reconstruction with Albada limiter
void ZFSStrctrdBlck3D::MusclAlbada()
{
  TRACE();
  //stencil identifier
  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSId IJK[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};

  //reduce to onedimensional arrays
  ZFSFloat* RESTRICT x=&m_cells->coordinates[0][0];
  ZFSFloat* RESTRICT y=&m_cells->coordinates[1][0];
  ZFSFloat* RESTRICT z=&m_cells->coordinates[2][0];
  ZFSFloat* RESTRICT flux = ALIGNED_F(m_cells->flux);
  ZFSFloat** RESTRICT pvars = m_cells->pvariables;

  /////////IMPORTANT PARAMETER
  //ZFSFloat epsi=F1;
  //ZFSFloat kappa=F1B3;
  /////////END IMPORTANT PARAMETER
  for(ZFSId dim=0; dim<nDim; dim++) {
    for(ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers; k++) {
      for(ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers; j++) {
        for(ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers; i++) {
          //cell ids
          const ZFSId I=cellIndex(i,j,k);
          const ZFSId IP1=I+IJK[dim];
          const ZFSId IM1=I-IJK[dim];
          const ZFSId IP2=I+2*IJK[dim];

          //distances q_i+1 - q_i
          const ZFSFloat DS=sqrt(POW2(x[IP1]-x[I]) + POW2(y[IP1]-y[I]) + POW2(z[IP1]-z[I]));
          //distances q_i - q_i-1
          const ZFSFloat DSM1=sqrt(POW2(x[I]-x[IM1]) + POW2(y[I]-y[IM1]) + POW2(z[I]-z[IM1]));
          const ZFSFloat DSP1=sqrt(POW2(x[IP2]-x[IP1]) + POW2(y[IP2]-y[IP1]) + POW2(z[IP2]-z[IP1]));
          const ZFSFloat DSP=DS/POW2(DSP1+DS);
          const ZFSFloat DSM=DS/POW2(DSM1+DS);

          const ZFSFloat pIM2 = pvars[PV->P][IM1];
          const ZFSFloat pIM1 = pvars[PV->P][I];
          const ZFSFloat pIP2 = pvars[PV->P][IP2];
          const ZFSFloat pIP1 = pvars[PV->P][IP1];

          const ZFSFloat smps=DS*DSP1;
          const ZFSFloat dummy=fabs(pIM2-F2*pIM1+pIP1)/(pIM2+F2*pIM1+pIP1);
          const ZFSFloat dummy1=fabs(pIM1-F2*pIP1+pIP2)/(pIM1+F2*pIP1+pIP2);
          const ZFSFloat psi=zfsMIN(F1,F6*zfsMAX(dummy,dummy1));
          const ZFSFloat epsLim=zfsMAX(m_eps, pow(F1B2*smps, F5));

          for(ZFSId var=0; var< PV->noVariables; ++var){
            const ZFSFloat DQ   = pvars[var][IP1]-pvars[var][I];
            const ZFSFloat DQP1 = pvars[var][IP2]-pvars[var][IP1];
            const ZFSFloat DQM1 = pvars[var][I]-pvars[var][IM1];
            const ZFSFloat phi=F1B2-(F1B2-zfsMAX(F0,(DQP1*DQM1*smps+F1B2*epsLim)/(POW2(DQP1*DS)+POW2(DQM1*DSP1)+epsLim)))*psi;

            m_QLeft[var]  = pvars[var][I]+DSM*(DSM1*DQ+DS*DQM1)*phi;
            m_QRight[var] = pvars[var][IP1]-DSP*(DS*DQP1+DSP1*DQ)*phi;
          }

          AusmLES(m_QLeft, m_QRight, dim, I); //Flux balance in AUSM
        }
      }
    }

    //FLUX BALANCE
    for(ZFSId v=0; v<CV->noVariables; v++) {
      for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
        for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
          for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
            const ZFSId I=cellIndex(i,j,k);
            const ZFSId IM1=I-IJK[dim];
            m_cells->rightHandSide[v][I]+=flux[IM1+noCells*v]-flux[I+noCells*v];
          }
        }
      }
    }
  }
}


/** \brief MUSCL with Venkatakrishan limiter
 * Here, MUSCL and AUSM are run through separately
 * The values of QLeft and QRight (of every cell) from the MUSCL are stored in the scratchspace
 * before passing it to the AUSM scheme.
 * Pros: the limiter can include the values of the neighbouring cells, better results for waves
 * Cons: higher computational effort
 * /author Leo Hoening, september 2015
 */
void ZFSStrctrdBlck3D::MusclVenkatakrishnan3D()
{
  TRACE();
  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSId IJK[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};
  //reduce to onedimensional arrays
  const ZFSFloat *const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const ZFSFloat *const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const ZFSFloat *const RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);
  const ZFSFloat *const RESTRICT cellVariables= ALIGNED_F(m_cells->pvariables[0]);
  ZFSFloat *const RESTRICT flux = ALIGNED_F(m_cells->flux);

  ZFSFloatScratchSpace QLeft(m_noStrctrdCells,PV->noVariables,3, __CALLING_FUNCTION__, "QLeft" );
  ZFSFloatScratchSpace QRight(m_noStrctrdCells,PV->noVariables,3, __CALLING_FUNCTION__, "QRight" );
  ZFSFloatScratchSpace minPhi(PV->noVariables, 2, __CALLING_FUNCTION__, "minPhi");

  QLeft.fill(F0);
  QRight.fill(F0);

  for(ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers; k++){
    for(ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers; j++){
      for(ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers; i++){
        for(ZFSId var=0; var<nDim+2; var++) {
          const ZFSUint offset = var*noCells;
          const ZFSFloat *const RESTRICT pvars = ALIGNED_F(cellVariables+offset);

          ZFSFloat minNghbrDelta = F0;
          ZFSFloat maxNghbrDelta = F0;
          ZFSFloat effNghbrDelta = F0;

          //1. get the abs min value of the max and min value from the reconstruction neighbours
          for(ZFSId dim1=0; dim1<nDim; dim1++) {
            for(ZFSId rcnstructnNghbr = 0; rcnstructnNghbr < 2 ; rcnstructnNghbr++) {
              const ZFSId cellId=cellIndex(i,j,k);
              const ZFSId IP1=cellId+IJK[dim1];
              const ZFSId IM1=cellId-IJK[dim1];

              ZFSFloat rcnstrctnNghbrValue = F0;
              if (rcnstructnNghbr==0){
                rcnstrctnNghbrValue = pvars[IP1]; //(i/j/k)+1
              }else{
                rcnstrctnNghbrValue = pvars[IM1]; //(i/j/k)-1
              }

              const ZFSFloat tmpDelta = rcnstrctnNghbrValue - pvars[cellId]; //i
              maxNghbrDelta = zfsMAX(maxNghbrDelta, tmpDelta);
              minNghbrDelta = zfsMIN(minNghbrDelta, tmpDelta);
            }
          }

          effNghbrDelta = zfsMIN(maxNghbrDelta, abs(minNghbrDelta));

          ZFSFloat srfcDelta = F0;
          ZFSFloat dxEpsSqr = F1;
          for(ZFSId dim1=0; dim1<nDim; dim1++) {
            const ZFSId cellId=cellIndex(i,j,k);
            const ZFSId IP1=cellId+IJK[dim1];
            const ZFSId IM1=cellId-IJK[dim1];

            const ZFSFloat DS=sqrt(POW2(x[IP1]-x[cellId]) + POW2(y[IP1]-y[cellId]) + POW2(z[IP1]-z[cellId]));
            //distances q_i - q_i-1
            const ZFSFloat DSM1=sqrt(POW2(x[cellId]-x[IM1]) + POW2(y[cellId]-y[IM1]) + POW2(z[cellId]-z[IM1]));
            const ZFSFloat DSM=DS/POW2(DSM1+DS);

            //2. get srfcDelta and compute the minimum phi
            const ZFSFloat dx1 = DSM * (DSM1*sqrt(POW2(x[IP1]-x[cellId]) + POW2(y[IP1]-y[cellId]) + POW2(z[IP1]-z[cellId])) + DS*sqrt(POW2(x[cellId]-x[IM1]) + POW2(y[cellId]-y[IM1]) + POW2(z[cellId]-z[IM1])));
            const ZFSFloat DQ=(pvars[IP1]-pvars[IM1])/sqrt(POW2(x[IP1]-x[IM1]) + POW2(y[IP1]-y[IM1]) + POW2(z[IP1]-z[IM1]));
            srfcDelta += abs(DQ*dx1);
            dxEpsSqr *= dx1;
          }

          ZFSId cellPos = 0;

          //calling limiter function
          (this->*Venkatakrishnan_function)(effNghbrDelta, srfcDelta, dxEpsSqr, cellPos, var, minPhi);

          minNghbrDelta = F0;
          maxNghbrDelta = F0;
          effNghbrDelta = F0;

          //1. get the abs min value of the max and min value from the reconstruction neighbours
          for(ZFSId dim1=0; dim1<nDim; dim1++) {
            for(ZFSId rcnstructnNghbr = 0; rcnstructnNghbr < 2 ; rcnstructnNghbr++) {
              const ZFSId cellId=cellIndex(i,j,k);
              const ZFSId IP1=cellId+IJK[dim1];
              const ZFSId IP2=cellId+2*IJK[dim1];

              ZFSFloat rcnstrctnNghbrValue = F0;
              if (rcnstructnNghbr==0){
                rcnstrctnNghbrValue = pvars[IP2]; //(i/j/k)+2
              } else {
                rcnstrctnNghbrValue = pvars[cellId]; //(i/j/k)
              }

              const ZFSFloat tmpDelta = rcnstrctnNghbrValue - pvars[IP1]; //(i/j/k)+1

              maxNghbrDelta = zfsMAX(maxNghbrDelta, tmpDelta);
              minNghbrDelta = zfsMIN(minNghbrDelta, tmpDelta);
            }
          }

          effNghbrDelta = zfsMIN(maxNghbrDelta, abs(minNghbrDelta));

          srfcDelta = F0;
          dxEpsSqr = F1;
          for(ZFSId dim1=0; dim1<nDim; dim1++) {
            const ZFSId cellId=cellIndex(i,j,k);
            const ZFSId IP1=cellId+IJK[dim1];
            const ZFSId IP2=cellId+2*IJK[dim1];

            const ZFSFloat DS=sqrt(POW2(x[IP1]-x[cellId]) + POW2(y[IP1]-y[cellId]) + POW2(z[IP1]-z[cellId]));
            //distances q_i - q_i-1
            const ZFSFloat DSP1=sqrt(POW2(x[IP2]-x[IP1]) + POW2(y[IP2]-y[IP1]) + POW2(z[IP2]-z[IP1]));
            const ZFSFloat DSP=DS/POW2(DSP1+DS);

            //2. get srfcDelta and compute the minimum phi

            const ZFSFloat dx2 = DSP * (DS*sqrt(POW2(x[IP2]-x[IP1]) + POW2(y[IP2]-y[IP1]) + POW2(z[IP2]-z[IP1])) + DSP1*sqrt(POW2(x[IP1]-x[cellId]) + POW2(y[IP1]-y[cellId]) + POW2(z[IP1]-z[cellId])));
            const ZFSFloat DQ=(pvars[IP2]-pvars[cellId])/sqrt(POW2(x[IP2]-x[cellId]) + POW2(y[IP2]-y[cellId]) + POW2(z[IP2]-z[cellId]));

            srfcDelta += abs(DQ*dx2);
            dxEpsSqr *= dx2;
          }

          cellPos = 1;
          (this->*Venkatakrishnan_function)(effNghbrDelta, srfcDelta, dxEpsSqr, cellPos, var, minPhi); //calling Venk

          for(ZFSId dim1=0; dim1<nDim; dim1++) {
            const ZFSId cellId=cellIndex(i,j,k);
            const ZFSId IP1=cellId+IJK[dim1];
            const ZFSId IM1=cellId-IJK[dim1];
            const ZFSId IP2=cellId+2*IJK[dim1];

            const ZFSFloat DS=sqrt(POW2(x[IP1]-x[cellId]) + POW2(y[IP1]-y[cellId]) + POW2(z[IP1]-z[cellId]));
            //distances q_i - q_i-1
            const ZFSFloat DSM1=sqrt(POW2(x[cellId]-x[IM1]) + POW2(y[cellId]-y[IM1]) + POW2(z[cellId]-z[IM1]));
            const ZFSFloat DSP1=sqrt(POW2(x[IP2]-x[IP1]) + POW2(y[IP2]-y[IP1]) + POW2(z[IP2]-z[IP1]));
            const ZFSFloat DSP=DS/POW2(DSP1+DS);
            const ZFSFloat DSM=DS/POW2(DSM1+DS);

            QLeft(cellId,var,dim1)  = pvars[cellId] + DSM*(DSM1*(pvars[IP1]-pvars[cellId])+DS*(pvars[cellId]-pvars[IM1]))*minPhi(var,0);
            QRight(cellId,var,dim1) = pvars[IP1]    - DSP*(DS*(pvars[IP2]-pvars[IP1])     +DSP1*(pvars[IP1]-pvars[cellId]))*minPhi(var,1);

          }
        }
      }
    }
  }

  ////// AUSM /////////
  for(ZFSId dim=0; dim<nDim; dim++) {
    for(ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers; k++) {
      for(ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers; j++) {
        for(ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers; i++) {
          //cell ids
          const ZFSId cellId=cellIndex(i,j,k);

          for(ZFSId v=0; v<PV->noVariables; v++) {
            m_QLeft[v] = QLeft(cellId,v,dim);
            m_QRight[v] = QRight(cellId,v,dim);
          }

          AusmLES(m_QLeft, m_QRight, dim, cellId);
        }
      }
    }


    //FLUX BALANCE
    for(ZFSId v=0; v<CV->noVariables; v++) {
      for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
        for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
          for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
            const ZFSId I=cellIndex(i,j,k);
            const ZFSId IM1=I-IJK[dim];
            m_cells->rightHandSide[v][I]+=flux[IM1+noCells*v]-flux[I+noCells*v];
          }
        }
      }
    }
  }
}


/** Venkatakrishnan limiter, modified for better results
 *  \author Leo Hoening, September 2015
 */
void ZFSStrctrdBlck3D::VENKATAKRISHNAN_MOD_FCT(ZFSFloat effNghbrDelta, ZFSFloat srfcDelta, ZFSFloat dxEpsSqr, ZFSId cellPos, ZFSId var, ZFSFloatScratchSpace& minPhi)
{
  ZFSFloat epsSqr = pow(m_venkFactor,F3)*dxEpsSqr;
  minPhi(var,cellPos) = ( pow(effNghbrDelta,F2) + epsSqr + F2*effNghbrDelta*srfcDelta ) / \
    ( pow(effNghbrDelta,F2) + F2*pow(srfcDelta,F2) + effNghbrDelta*srfcDelta + epsSqr );
}

/** Standard Venkatakrishnan limiter
 *  \author Leo Hoening, September 2015
 */
void ZFSStrctrdBlck3D::VENKATAKRISHNAN_FCT(ZFSFloat effNghbrDelta, ZFSFloat srfcDelta, ZFSFloat dxEpsSqr, ZFSId cellPos, ZFSId var, ZFSFloatScratchSpace& minPhi)
{
  (void) dxEpsSqr;
  const ZFSFloat eps = 1e-12;
  ZFSFloat yps1 = effNghbrDelta/(srfcDelta+eps);
  minPhi(var,cellPos) = zfsMIN((yps1*yps1 + F2*yps1)/(yps1*yps1 + yps1 + F2),F1);
}

/** Barth-Jesperson Limiter
 *  \author Leo Hoening, September 2015
 */
void ZFSStrctrdBlck3D::BARTH_JESPERSON_FCT(ZFSFloat effNghbrDelta, ZFSFloat srfcDelta, ZFSFloat dxEpsSqr, ZFSId cellPos, ZFSId var, ZFSFloatScratchSpace& minPhi)
{
  (void) dxEpsSqr;
  const ZFSFloat eps = 1e-12;
  ZFSFloat phi_max = effNghbrDelta/(srfcDelta+eps);
  minPhi(var,cellPos) = zfsMIN(phi_max, F1);
}

/** \brief AUSM CENTRAL as in TFS
 *  can be used for moving grids, dxt term is included
 */
//inline void ZFSStrctrdBlck3D::AusmNew(ZFSFloat* QLeft, ZFSFloat* QRight, const ZFSId dim, const ZFSId cellId)
inline void ZFSStrctrdBlck3D::AusmLES(ZFSFloat* RESTRICT QLeft, ZFSFloat* RESTRICT QRight, const ZFSId dim, const ZFSId I)
{
  // ZFSFloat pFactor[3]={F0,F0,F0};
  const ZFSFloat gamma = m_gamma;
  const ZFSFloat gammaMinusOne = gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_cells->surfaceMetrics[ I ]);

  const ZFSFloat dxdtau = m_cells->dxt[dim][I];

  //calculate pressure
  const ZFSFloat PL= QLeft[PV->P];
  const ZFSFloat UL   = QLeft[ PV->U ];
  const ZFSFloat VL   = QLeft[ PV->V ];
  const ZFSFloat WL   = QLeft[ PV->W ];
  const ZFSFloat RHOL = QLeft[ PV->RHO ];

  const ZFSFloat PR= QRight[PV->P];
  const ZFSFloat UR   = QRight[ PV->U ];
  const ZFSFloat VR   = QRight[ PV->V ];
  const ZFSFloat WR   = QRight[ PV->W ];
  const ZFSFloat RHOR = QRight[ PV->RHO ];

  // compute lenght of metric vector for normalization
  const ZFSFloat metricLength = sqrt(POW2( surf[dim*3+0] ) + POW2(surf[dim*3+1]) + POW2(surf[dim*3+2]));
  const ZFSFloat fMetricLength = F1 / metricLength;

  //scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const ZFSFloat UUL = ((UL * surf[ dim*3+0 ] +
                         VL * surf[ dim*3+1 ] +
                         WL * surf[ dim*3+2 ]) - dxdtau) * fMetricLength;


  const ZFSFloat UUR = ((UR * surf[ dim*3+0 ] +
                         VR * surf[ dim*3+1 ] +
                         WR * surf[ dim*3+2 ]) - dxdtau) * fMetricLength;


  //speed of sound
  const ZFSFloat AL = sqrt(gamma * zfsMAX(m_eps, PL / zfsMAX(m_eps, RHOL)));
  const ZFSFloat AR = sqrt(gamma * zfsMAX(m_eps, PR / zfsMAX(m_eps, RHOR)));

  const ZFSFloat MAL = UUL / AL;
  const ZFSFloat MAR = UUR / AR;

  const ZFSFloat MALR = F1B2*(MAL+MAR);
  const ZFSFloat PLR = PL*(F1B2+ m_chi*MAL) + PR*(F1B2- m_chi*MAR) ;

  const ZFSFloat RHO_AL = RHOL*AL;
  const ZFSFloat RHO_AR = RHOR*AR;

  const ZFSFloat PLfRHOL = PL/RHOL;
  const ZFSFloat PRfRHOR = PR/RHOR;

  const ZFSFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
  const ZFSFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;

  const ZFSFloat RHOU = F1B2 * ( MALR * (RHO_AL+RHO_AR) + fabs(MALR) * (RHO_AL-RHO_AR) ) * metricLength;
  const ZFSFloat RHOU2 = F1B2*RHOU;
  // multiply by metric length to take surface area into account
  const ZFSFloat AbsRHO_U2 = fabs(RHOU2);

  const ZFSFloat pFactor0 = surf[ 3 * dim + 0 ];
  const ZFSFloat pFactor1 = surf[ 3 * dim + 1 ];
  const ZFSFloat pFactor2 = surf[ 3 * dim + 2 ];
  const ZFSInt noCells = m_noStrctrdCells;

  //==>fluxes:
  ZFSFloat* RESTRICT flux = ALIGNED_F(m_cells->flux);
  flux[I+noCells*CV->RHO_U] = RHOU2 * ( UL + UR ) + AbsRHO_U2 * ( UL - UR ) + PLR * pFactor0;
  flux[I+noCells*CV->RHO_V] = RHOU2 * ( VL + VR ) + AbsRHO_U2 * ( VL - VR ) + PLR * pFactor1;
  flux[I+noCells*CV->RHO_W] = RHOU2 * ( WL + WR ) + AbsRHO_U2 * ( WL - WR ) + PLR * pFactor2;
  flux[I+noCells*CV->RHO_E] = RHOU2 * ( e0 + e1)  + AbsRHO_U2 * ( e0 - e1 ) + PLR * dxdtau;
  flux[I+noCells*CV->RHO]   = RHOU;
}


/**
 *  Same AUSM scheme as AusmNew with additional damping controlled
 *  by the 4th order pressure derivative. Pressure needs to computed
 *  beforehand.
 *
 */
inline void ZFSStrctrdBlck3D::AusmLES_PTHRC(ZFSFloat* QLeft, ZFSFloat* QRight, ZFSId dim, ZFSId I)
{
  ZFSFloat pFactor[3]={F0,F0,F0};
  const ZFSFloat gamma = m_gamma;
  const ZFSFloat FgammaMinusOne = m_fgammaMinusOne;

  const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_cells->surfaceMetrics[ I ]);
  const ZFSFloat *const RESTRICT p = ALIGNED_F(m_cells->pvariables[ PV->P ]);

  const ZFSFloat dxdtau = m_cells->dxt[dim][I];

  //calculate pressure
  const ZFSFloat PL= QLeft[PV->P];
  const ZFSFloat UL   = QLeft[ PV->U ];
  const ZFSFloat VL   = QLeft[ PV->V ];
  const ZFSFloat WL   = QLeft[ PV->W ];
  const ZFSFloat RHOL = QLeft[ PV->RHO ];

  const ZFSFloat PR= QRight[PV->P];
  const ZFSFloat UR   = QRight[ PV->U ];
  const ZFSFloat VR   = QRight[ PV->V ];
  const ZFSFloat WR   = QRight[ PV->W ];
  const ZFSFloat RHOR = QRight[ PV->RHO ];

  // compute lenght of metric vector for normalization
  const ZFSFloat metricLength = sqrt(POW2( surf[dim*3+0] ) + POW2(surf[dim*3+1]) + POW2(surf[dim*3+2]));
  const ZFSFloat fMetricLength = F1 / metricLength;

  //scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const ZFSFloat UUL = ((UL * surf[ dim*3+0 ] +
                         VL * surf[ dim*3+1 ] +
                         WL * surf[ dim*3+2 ]) - dxdtau) * fMetricLength;


  const ZFSFloat UUR = ((UR * surf[ dim*3+0 ] +
                         VR * surf[ dim*3+1 ] +
                         WR * surf[ dim*3+2 ]) - dxdtau) * fMetricLength;


  //speed of sound
  const ZFSFloat AL = sqrt(gamma * zfsMAX(m_eps, PL / zfsMAX(m_eps, RHOL)));
  const ZFSFloat AR = sqrt(gamma * zfsMAX(m_eps, PR / zfsMAX(m_eps, RHOR)));

  const ZFSFloat MAL = UUL / AL;
  const ZFSFloat MAR = UUR / AR;

  const ZFSFloat MALR = F1B2*(MAL+MAR);

  //4th order pressure damping
  const ZFSId IPJK  = getCellIdfromCell(I,  1,0,0);
  const ZFSId IMJK  = getCellIdfromCell(I, -1,0,0);
  const ZFSId IP2JK = getCellIdfromCell(I,  2,0,0);
  const ZFSId IM2JK = getCellIdfromCell(I, -2,0,0);

  const ZFSId IJPK = getCellIdfromCell(I, 0,1,0);
  const ZFSId IJMK = getCellIdfromCell(I, 0,-1,0);
  const ZFSId IJP2K = getCellIdfromCell(I, 0,2,0);
  const ZFSId IJM2K = getCellIdfromCell(I, 0,-2,0);

  const ZFSId IJKP = getCellIdfromCell(I, 0,0,1);
  const ZFSId IJKM = getCellIdfromCell(I, 0,0,-1);
  const ZFSId IJKP2 = getCellIdfromCell(I, 0,0,2);
  const ZFSId IJKM2 = getCellIdfromCell(I, 0,0,-2);

  const ZFSFloat p4I4 = F4 * (p[IPJK] + p[IMJK]) -F6 * (p[I]) - p[IP2JK] - p[IM2JK];
  const ZFSFloat p4J4 = F4 * (p[IJPK] + p[IJMK]) -F6 * (p[I]) - p[IJP2K] - p[IJM2K];
  const ZFSFloat p4K4 = F4 * (p[IJKP] + p[IJKM]) -F6 * (p[I]) - p[IJKP2] - p[IJKM2];

  const ZFSFloat cfac = 1.0/1.3;
  const ZFSFloat pfac = fabs(p4I4) + fabs(p4J4) + fabs(p4K4);
  ZFSFloat fac = cfac*pfac;
  fac = min(1/64.0,fac*5.0);

  const ZFSFloat PLR= PL*(F1B2+ fac * MAL) + PR*(F1B2- fac * MAR);

  const ZFSFloat RHO_AL = RHOL*AL;
  const ZFSFloat RHO_AR = RHOR*AR;

  const ZFSFloat PLfRHOL = PL/RHOL;
  const ZFSFloat PRfRHOR = PR/RHOR;

  const ZFSFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
  const ZFSFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;

  const ZFSFloat RHOU = F1B2 * ( MALR * (RHO_AL+RHO_AR) + fabs(MALR) * (RHO_AL-RHO_AR) ) * metricLength;
  const ZFSFloat RHOU2 = F1B2*RHOU;
  // multiply by metric length to take surface area into account
  const ZFSFloat AbsRHO_U2 = fabs(RHOU2);

  // setup pressure factors, include metric terms (not normalized due to needed surface area)
  for( ZFSInt isd = xsd; isd < nDim; isd ++ ) {
    pFactor[ isd ] = surf[ 3 * dim + isd ];
  }

  const ZFSInt noCells = m_noStrctrdCells;
  ZFSFloat* RESTRICT flux = ALIGNED_F(m_cells->flux);
  flux[I+noCells*CV->RHO_U] = RHOU2 * ( UL + UR ) + AbsRHO_U2 * ( UL - UR ) + PLR * pFactor[ 0 ];
  flux[I+noCells*CV->RHO_V] = RHOU2 * ( VL + VR ) + AbsRHO_U2 * ( VL - VR ) + PLR * pFactor[ 1 ];
  flux[I+noCells*CV->RHO_W] = RHOU2 * ( WL + WR ) + AbsRHO_U2 * ( WL - WR ) + PLR * pFactor[ 2 ];
  flux[I+noCells*CV->RHO_E] = RHOU2 * ( e0 + e1)  + AbsRHO_U2 * ( e0 - e1 ) + PLR * dxdtau;
  flux[I+noCells*CV->RHO]   = RHOU;
}

void ZFSStrctrdBlck3D::AusmDV(ZFSFloat* QLeft, ZFSFloat*  QRight, const ZFSId dim, const ZFSId I){
  ZFSFloat pFactor[3]={F0,F0,F0};
  const ZFSFloat gamma = m_gamma;
  const ZFSFloat gammaMinusOne = gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_cells->surfaceMetrics[I]);
  const ZFSFloat dxdtau = m_cells->dxt[dim][I];

  //left side
  const ZFSFloat RHOL = QLeft[PV->RHO];
  const ZFSFloat FRHOL = F1/RHOL;
  ZFSFloat UL = QLeft[PV->U];
  ZFSFloat VL = QLeft[PV->V];
  ZFSFloat WL = QLeft[PV->W];
  const ZFSFloat PL = QLeft[PV->P];

  //right side
  const ZFSFloat RHOR = QRight[PV->RHO];
  const ZFSFloat FRHOR = F1/RHOR;
  ZFSFloat UR = QRight[PV->U];
  ZFSFloat VR = QRight[PV->V];
  ZFSFloat WR = QRight[PV->W];
  const ZFSFloat PR = QRight[PV->P];

  const ZFSFloat PLfRHOL = PL/RHOL;
  const ZFSFloat PRfRHOR = PR/RHOR;
  const ZFSFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
  const ZFSFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;


  // compute lenght of metric vector for normalization
  const ZFSFloat DGRAD = sqrt(POW2( surf[dim*3+0] ) + POW2(surf[dim*3+1]) + POW2(surf[dim*3+2]));
  const ZFSFloat FDGRAD = F1 / DGRAD;

  //scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const ZFSFloat UUL = ((UL * surf[ dim*3+0 ] +
                         VL * surf[ dim*3+1 ] +
                         WL * surf[ dim*3+2 ]) - dxdtau) * FDGRAD;


  const ZFSFloat UUR = ((UR * surf[ dim*3+0 ] +
                         VR * surf[ dim*3+1 ] +
                         WR * surf[ dim*3+2 ]) - dxdtau) * FDGRAD;

  ZFSFloat AL = FRHOL*PL;
  ZFSFloat AR = FRHOR*PR;

  const ZFSFloat FALR = 2.0/(AL + AR);
  const ZFSFloat ALPHAL = AL*FALR;
  const ZFSFloat ALPHAR = AR*FALR;

  AL = sqrt(gamma*AL);
  AR = sqrt(gamma*AR);
  AL = zfsMAX(AL,AR);
  AR = AL;

  const ZFSFloat XMAL = UUL/AL;
  const ZFSFloat XMAR = UUR/AR;

  AL = AL*DGRAD;
  AR = AR*DGRAD;

  const ZFSFloat RHOAL = AL*RHOL;
  const ZFSFloat RHOAR = AR*RHOR;

  const ZFSId IJK[2]={1,m_nCells[1]};
  const ZFSId IP1=I+IJK[dim];

  const ZFSFloat FDV = 0.3;
  const ZFSFloat DXDXEZ= m_cells->coordinates[0][IP1] - m_cells->coordinates[0][I];
  const ZFSFloat DYDXEZ= m_cells->coordinates[1][IP1] - m_cells->coordinates[1][I];
  const ZFSFloat DZDXEZ= m_cells->coordinates[2][IP1] - m_cells->coordinates[2][I];
  ZFSFloat SV = 2.0*DGRAD/(m_cells->cellJac[I]+m_cells->cellJac[IP1])*
    (FDV+(F1-FDV)*getPSI(I,dim));
  const ZFSFloat SV1 = F0*SV*DXDXEZ;
  const ZFSFloat SV2 = F0*SV*DYDXEZ;
  const ZFSFloat SV3 = F0*SV*DZDXEZ;

  const ZFSFloat XMAL1 = zfsMIN(F1,zfsMAX(-F1,XMAL));
  const ZFSFloat XMAR1 = zfsMIN(F1,zfsMAX(-F1,XMAR));

  ZFSFloat FXMA = F1B2*(XMAL1+fabs(XMAL1));
  const ZFSFloat XMALP = ALPHAL*(F1B4*POW2(XMAL1+F1)-FXMA)+FXMA+(zfsMAX(F1,XMAL)-F1);
  FXMA = F1B2*(XMAR1-fabs(XMAR1));
  const ZFSFloat XMARM = ALPHAR*(-F1B4*POW2(XMAR1-F1)-FXMA)+FXMA+(zfsMIN(-F1,XMAR)+F1);

  const ZFSFloat FLP = PL*((F2-XMAL1)*POW2(F1+XMAL1));
  const ZFSFloat FRP = PR*((F2+XMAR1)*POW2(F1-XMAR1));
  const ZFSFloat PLR = F1B4*(FLP+FRP);

  const ZFSFloat RHOUL = XMALP*RHOAL;
  const ZFSFloat RHOUR = XMARM*RHOAR;
  const ZFSFloat RHOU = RHOUL+RHOUR;
  const ZFSFloat RHOU2 = F1B2*RHOU;
  const ZFSFloat ARHOU2 = fabs(RHOU2);

  const ZFSFloat UUL2 = SV1*UUL;
  const ZFSFloat UUR2 = SV1*UUR;
  UL  = UL-UUL2;
  UR  = UR-UUR2;
  const ZFSFloat UUL3 = SV2*UUL;
  const ZFSFloat UUR3 = SV2*UUR;
  VL  = VL-UUL3;
  VR  = VR-UUR3;
  const ZFSFloat UUL4 = SV3*UUL;
  const ZFSFloat UUR4 = SV3*UUR;
  WL  = WL-UUL4;
  WR  = WR-UUR4;

  // setup pressure factors, include metric terms (not normalized due to needed surface area)
  for( ZFSInt isd = xsd; isd < nDim; isd ++ ) {
    pFactor[ isd ] = surf[ 3 * dim + isd ];
  }

  const ZFSInt noCells = m_noStrctrdCells;
  ZFSFloat *const RESTRICT flux = m_cells->flux;

  flux[I+noCells*CV->RHO_U] = RHOU2 * ( UL + UR) + ARHOU2 * ( UL - UR) + PLR*pFactor[0] + RHOUL*UUL2+RHOUR*UUR2;
  flux[I+noCells*CV->RHO_V] = RHOU2 * ( VL + VR) + ARHOU2 * ( VL - VR) + PLR*pFactor[1] + RHOUL*UUL3+RHOUR*UUR3;
  flux[I+noCells*CV->RHO_W] = RHOU2 * ( WL + WR) + ARHOU2 * ( WL - WR) + PLR*pFactor[2] + RHOUL*UUL4+RHOUR*UUR4;
  flux[I+noCells*CV->RHO_E] = RHOU2 * ( e0 + e1) + ARHOU2 * ( e0 - e1) + PLR*dxdtau;
  flux[I+noCells*CV->RHO]   = RHOU;
}

template<ZFSId noVars>
void ZFSStrctrdBlck3D::Muscl_AusmLES(){
  TRACE();

  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSId IJK[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};

  const ZFSFloat *const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const ZFSFloat *const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const ZFSFloat *const RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);
  const ZFSFloat *const *const RESTRICT vars= ALIGNED_F(m_cells->pvariables);
  ZFSFloat *const RESTRICT ds = ALIGNED_F(m_cells->ds);  
  ZFSFloat *const RESTRICT cellRhs= ALIGNED_MF(m_cells->rightHandSide[0]);
  ZFSFloat *const RESTRICT flux = ALIGNED_F(m_cells->flux);

  const ZFSFloat gamma = m_gamma;
  const ZFSFloat gammaMinusOne = gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  const ZFSInt noCellsI = m_nCells[2]-2;
  const ZFSInt noCellsJ = m_nCells[1]-2;
  const ZFSInt noCellsK = m_nCells[0]-2;

  const ZFSInt noCellsIP1 = m_nCells[2]-1;
  const ZFSInt noCellsJP1 = m_nCells[1]-1;
  const ZFSInt noCellsKP1 = m_nCells[0]-1;

  for(ZFSId dim=0; dim<nDim; dim++){
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(ZFSId k=0; k<noCellsKP1; k++) {
      for(ZFSId j=0; j<noCellsJP1; j++) {
        for(ZFSId i=0; i<noCellsIP1; i++) {
          const ZFSId I=cellIndex(i,j,k);
          const ZFSId IP1=I+IJK[dim];
          ds[I] = sqrt(POW2(x[IP1]-x[I]) + POW2(y[IP1]-y[I]) + POW2(z[IP1]-z[I]));
        }
      }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(ZFSId k=1; k<noCellsK; k++) {
      for(ZFSId j=1; j<noCellsJ; j++) {
#if defined(ZFS_INTEL_COMPILER)
#pragma ivdep 
#pragma vector always
#endif
        for(ZFSId i=1; i<noCellsI; i++) {
          const ZFSId I=cellIndex(i,j,k);
          const ZFSId IP1=I+IJK[dim];
          const ZFSId IM1=I-IJK[dim];
          const ZFSId IP2=I+2*IJK[dim];

          const ZFSFloat DS = ds[I];
          const ZFSFloat DSM1 = ds[IM1];
          const ZFSFloat DSP1 = ds[IP1];

          const ZFSFloat DSP=DS/POW2(DSP1+DS);
          const ZFSFloat DSM=DS/POW2(DSM1+DS);

          //unrolled the loop so the compiler
          //can optimize better
          const ZFSFloat DQU   = vars[PV->U][IP1]-vars[PV->U][I];
          const ZFSFloat DQPU = vars[PV->U][IP2]-vars[PV->U][IP1];
          const ZFSFloat DQMU = vars[PV->U][I]-vars[PV->U][IM1];
          const ZFSFloat UL  = vars[PV->U][I]+DSM*(DSM1*DQU+DS*DQMU);
          const ZFSFloat UR= vars[PV->U][IP1]-DSP*(DS*DQPU+DSP1*DQU);

          const ZFSFloat DQV   = vars[PV->V][IP1]-vars[PV->V][I];
          const ZFSFloat DQPV = vars[PV->V][IP2]-vars[PV->V][IP1];
          const ZFSFloat DQMV = vars[PV->V][I]-vars[PV->V][IM1];
          const ZFSFloat VL = vars[PV->V][I]+DSM*(DSM1*DQV+DS*DQMV);
          const ZFSFloat VR = vars[PV->V][IP1]-DSP*(DS*DQPV+DSP1*DQV);

          const ZFSFloat DQW   = vars[PV->W][IP1]-vars[PV->W][I];
          const ZFSFloat DQPW = vars[PV->W][IP2]-vars[PV->W][IP1];
          const ZFSFloat DQMW = vars[PV->W][I]-vars[PV->W][IM1];
          const ZFSFloat WL = vars[PV->W][I]+DSM*(DSM1*DQW+DS*DQMW);
          const ZFSFloat WR = vars[PV->W][IP1]-DSP*(DS*DQPW+DSP1*DQW);

          const ZFSFloat DQP   = vars[PV->P][IP1]-vars[PV->P][I];
          const ZFSFloat DQPP = vars[PV->P][IP2]-vars[PV->P][IP1];
          const ZFSFloat DQMP = vars[PV->P][I]-vars[PV->P][IM1];
          const ZFSFloat PL = vars[PV->P][I]+DSM*(DSM1*DQP+DS*DQMP);
          const ZFSFloat PR = vars[PV->P][IP1]-DSP*(DS*DQPP+DSP1*DQP);

          const ZFSFloat DQRHO   = vars[PV->RHO][IP1]-vars[PV->RHO][I];
          const ZFSFloat DQPRHO = vars[PV->RHO][IP2]-vars[PV->RHO][IP1];
          const ZFSFloat DQMRHO = vars[PV->RHO][I]-vars[PV->RHO][IM1];
          const ZFSFloat RHOL = vars[PV->RHO][I]+DSM*(DSM1*DQRHO+DS*DQMRHO);
          const ZFSFloat RHOR = vars[PV->RHO][IP1]-DSP*(DS*DQPRHO+DSP1*DQRHO);

          const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_cells->surfaceMetrics[ I ]);
          const ZFSFloat dxdtau = m_cells->dxt[dim][I];

          // compute length of metric vector for normalization
          const ZFSFloat metricLength = sqrt(POW2( surf[dim*3+0] ) + POW2(surf[dim*3+1]) + POW2(surf[dim*3+2]));
          const ZFSFloat fMetricLength = F1 / metricLength;

          //scale by metric length to get velocity in the new basis (get normalized basis vectors)
          const ZFSFloat UUL = ((UL * surf[ dim*3+0 ] +
                                 VL * surf[ dim*3+1 ] +
                                 WL * surf[ dim*3+2 ]) - dxdtau) * fMetricLength;


          const ZFSFloat UUR = ((UR * surf[ dim*3+0 ] +
                                 VR * surf[ dim*3+1 ] +
                                 WR * surf[ dim*3+2 ]) - dxdtau) * fMetricLength;


          //speed of sound
          const ZFSFloat AL = sqrt(gamma * max(m_eps,(PL/max(m_eps,RHOL))));
          const ZFSFloat AR = sqrt(gamma * max(m_eps,(PR/max(m_eps,RHOR))));

          const ZFSFloat MAL = UUL / AL;
          const ZFSFloat MAR = UUR / AR;

          const ZFSFloat MALR = F1B2*(MAL+MAR);
          const ZFSFloat PLR = PL*(F1B2+ m_chi*MAL) + PR*(F1B2- m_chi*MAR) ;

          const ZFSFloat RHO_AL = RHOL*AL;
          const ZFSFloat RHO_AR = RHOR*AR;

          const ZFSFloat PLfRHOL = PL/RHOL;
          const ZFSFloat PRfRHOR = PR/RHOR;

          const ZFSFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
          const ZFSFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;

          const ZFSFloat RHOU = F1B2 * ( MALR * (RHO_AL+RHO_AR) + fabs(MALR) * (RHO_AL-RHO_AR) ) * metricLength;
          const ZFSFloat RHOU2 = F1B2*RHOU;
          // multiply by metric length to take surface area into account
          const ZFSFloat AbsRHO_U2 = fabs(RHOU2);

          // setup pressure factors, include metric terms (not normalized due to needed surface area)
          const ZFSFloat pFactor0 = surf[ 3 * dim + 0 ];
          const ZFSFloat pFactor1 = surf[ 3 * dim + 1 ];
          const ZFSFloat pFactor2 = surf[ 3 * dim + 2 ];
  
          flux[I+noCells*CV->RHO_U] = RHOU2 * ( UL + UR ) + AbsRHO_U2 * ( UL - UR ) + PLR * pFactor0;
          flux[I+noCells*CV->RHO_V] = RHOU2 * ( VL + VR ) + AbsRHO_U2 * ( VL - VR ) + PLR * pFactor1;
          flux[I+noCells*CV->RHO_W] = RHOU2 * ( WL + WR ) + AbsRHO_U2 * ( WL - WR ) + PLR * pFactor2;
          flux[I+noCells*CV->RHO_E] = RHOU2 * ( e0 + e1)  + AbsRHO_U2 * ( e0 - e1 ) + PLR * dxdtau;
          flux[I+noCells*CV->RHO]   = RHOU;
        }
      }
    }

    //FLUX BALANCE
    for(ZFSUint v=0; v<noVars; v++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
        for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
#if defined(ZFS_INTEL_COMPILER)
#pragma ivdep 
#pragma vector always
#endif
          for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
            const ZFSId I = cellIndex(i,j,k);
            const ZFSId IM1=I-IJK[dim];
            const ZFSUint offset = v*noCells;
            ZFSFloat *const RESTRICT rhs = ALIGNED_F(cellRhs+offset);
            rhs[I]+=flux[IM1+noCells*v]-flux[I+noCells*v];
          }
        }
      }
    }
  }
}

template void ZFSStrctrdBlck3D::Muscl_AusmLES<5>();
template void ZFSStrctrdBlck3D::Muscl_AusmLES<6>();
template void ZFSStrctrdBlck3D::Muscl_AusmLES<7>();

template<ZFSId noVars>
void ZFSStrctrdBlck3D::Muscl_AusmLES_PTHRC(){
  TRACE();

  //stencil identifier
  const ZFSId IJK[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};

  const ZFSFloat *const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const ZFSFloat *const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const ZFSFloat *const RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);
  const ZFSFloat *const *const  RESTRICT vars= ALIGNED_F(m_cells->pvariables);
  const ZFSFloat *const RESTRICT p = ALIGNED_F(m_cells->pvariables[ PV->P ]);
  ZFSFloat *const RESTRICT ds = ALIGNED_F(m_cells->ds);  
  ZFSFloat *const RESTRICT cellRhs= ALIGNED_MF(m_cells->rightHandSide[0]);
  ZFSFloat *const RESTRICT flux = ALIGNED_F(m_cells->flux);

  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSFloat gamma = m_gamma;
  const ZFSFloat gammaMinusOne = gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  const ZFSInt noCellsI = m_nCells[2]-2;
  const ZFSInt noCellsJ = m_nCells[1]-2;
  const ZFSInt noCellsK = m_nCells[0]-2;

  const ZFSInt noCellsIP1 = m_nCells[2]-1;
  const ZFSInt noCellsJP1 = m_nCells[1]-1;
  const ZFSInt noCellsKP1 = m_nCells[0]-1;

  for(ZFSId dim=0; dim<nDim; dim++){
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(ZFSId k=0; k<noCellsKP1; k++) {
      for(ZFSId j=0; j<noCellsJP1; j++) {
        for(ZFSId i=0; i<noCellsIP1; i++) {
          const ZFSId I=cellIndex(i,j,k);
          const ZFSId IP1=I+IJK[dim];
          ds[I] = sqrt(POW2(x[IP1]-x[I]) + POW2(y[IP1]-y[I]) + POW2(z[IP1]-z[I]));
        }
      }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(ZFSId k=1; k<noCellsK; k++) {
      for(ZFSId j=1; j<noCellsJ; j++) {
#if defined(ZFS_INTEL_COMPILER)
#pragma ivdep 
#pragma vector always
#endif
        for(ZFSId i=1; i<noCellsI; i++) {
          const ZFSId I=cellIndex(i,j,k);
          const ZFSId IP1=I+IJK[dim];
          const ZFSId IM1=I-IJK[dim];
          const ZFSId IP2=I+2*IJK[dim];

          const ZFSFloat DS = ds[I];
          const ZFSFloat DSM1 = ds[IM1];
          const ZFSFloat DSP1 = ds[IP1];

          const ZFSFloat DSP=DS/POW2(DSP1+DS);
          const ZFSFloat DSM=DS/POW2(DSM1+DS);

          //unrolled the loop so the compiler
          //can optimize better
          const ZFSFloat DQU   = vars[PV->U][IP1]-vars[PV->U][I];
          const ZFSFloat DQPU = vars[PV->U][IP2]-vars[PV->U][IP1];
          const ZFSFloat DQMU = vars[PV->U][I]-vars[PV->U][IM1];
          const ZFSFloat UL  = vars[PV->U][I]+DSM*(DSM1*DQU+DS*DQMU);
          const ZFSFloat UR= vars[PV->U][IP1]-DSP*(DS*DQPU+DSP1*DQU);

          const ZFSFloat DQV   = vars[PV->V][IP1]-vars[PV->V][I];
          const ZFSFloat DQPV = vars[PV->V][IP2]-vars[PV->V][IP1];
          const ZFSFloat DQMV = vars[PV->V][I]-vars[PV->V][IM1];
          const ZFSFloat VL = vars[PV->V][I]+DSM*(DSM1*DQV+DS*DQMV);
          const ZFSFloat VR = vars[PV->V][IP1]-DSP*(DS*DQPV+DSP1*DQV);

          const ZFSFloat DQW   = vars[PV->W][IP1]-vars[PV->W][I];
          const ZFSFloat DQPW = vars[PV->W][IP2]-vars[PV->W][IP1];
          const ZFSFloat DQMW = vars[PV->W][I]-vars[PV->W][IM1];
          const ZFSFloat WL = vars[PV->W][I]+DSM*(DSM1*DQW+DS*DQMW);
          const ZFSFloat WR = vars[PV->W][IP1]-DSP*(DS*DQPW+DSP1*DQW);

          const ZFSFloat DQP   = vars[PV->P][IP1]-vars[PV->P][I];
          const ZFSFloat DQPP = vars[PV->P][IP2]-vars[PV->P][IP1];
          const ZFSFloat DQMP = vars[PV->P][I]-vars[PV->P][IM1];
          const ZFSFloat PL = vars[PV->P][I]+DSM*(DSM1*DQP+DS*DQMP);
          const ZFSFloat PR = vars[PV->P][IP1]-DSP*(DS*DQPP+DSP1*DQP);

          const ZFSFloat DQRHO   = vars[PV->RHO][IP1]-vars[PV->RHO][I];
          const ZFSFloat DQPRHO = vars[PV->RHO][IP2]-vars[PV->RHO][IP1];
          const ZFSFloat DQMRHO = vars[PV->RHO][I]-vars[PV->RHO][IM1];
          const ZFSFloat RHOL = vars[PV->RHO][I]+DSM*(DSM1*DQRHO+DS*DQMRHO);
          const ZFSFloat RHOR = vars[PV->RHO][IP1]-DSP*(DS*DQPRHO+DSP1*DQRHO);

          const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_cells->surfaceMetrics[ I ]);
          const ZFSFloat dxdtau = m_cells->dxt[dim][I];

          // compute lenght of metric vector for normalization
          const ZFSFloat metricLength = sqrt(POW2( surf[dim*3+0] ) + POW2(surf[dim*3+1]) + POW2(surf[dim*3+2]));
          const ZFSFloat fMetricLength = F1 / metricLength;

          //scale by metric length to get velocity in the new basis (get normalized basis vectors)
          const ZFSFloat UUL = ((UL * surf[ dim*3+0 ] +
                                 VL * surf[ dim*3+1 ] +
                                 WL * surf[ dim*3+2 ]) - dxdtau) * fMetricLength;


          const ZFSFloat UUR = ((UR * surf[ dim*3+0 ] +
                                 VR * surf[ dim*3+1 ] +
                                 WR * surf[ dim*3+2 ]) - dxdtau) * fMetricLength;


          //speed of sound
          const ZFSFloat AL = sqrt(gamma * max(m_eps,(PL/max(m_eps,RHOL))));
          const ZFSFloat AR = sqrt(gamma * max(m_eps,(PR/max(m_eps,RHOR))));

          const ZFSFloat MAL = UUL / AL;
          const ZFSFloat MAR = UUR / AR;

          const ZFSFloat MALR = F1B2*(MAL+MAR);

          //4th order pressure damping
          const ZFSId IPJK  = getCellIdfromCell(I,  1,0,0);
          const ZFSId IMJK  = getCellIdfromCell(I, -1,0,0);
          const ZFSId IP2JK = getCellIdfromCell(I,  2,0,0);
          const ZFSId IM2JK = getCellIdfromCell(I, -2,0,0);

          const ZFSId IJPK = getCellIdfromCell(I, 0,1,0);
          const ZFSId IJMK = getCellIdfromCell(I, 0,-1,0);
          const ZFSId IJP2K = getCellIdfromCell(I, 0,2,0);
          const ZFSId IJM2K = getCellIdfromCell(I, 0,-2,0);

          const ZFSId IJKP = getCellIdfromCell(I, 0,0,1);
          const ZFSId IJKM = getCellIdfromCell(I, 0,0,-1);
          const ZFSId IJKP2 = getCellIdfromCell(I, 0,0,2);
          const ZFSId IJKM2 = getCellIdfromCell(I, 0,0,-2);

          const ZFSFloat p4I4 = F4 * (p[IPJK] + p[IMJK]) -F6 * (p[I]) - p[IP2JK] - p[IM2JK];
          const ZFSFloat p4J4 = F4 * (p[IJPK] + p[IJMK]) -F6 * (p[I]) - p[IJP2K] - p[IJM2K];
          const ZFSFloat p4K4 = F4 * (p[IJKP] + p[IJKM]) -F6 * (p[I]) - p[IJKP2] - p[IJKM2];

          const ZFSFloat cfac = 1.0/1.3;
          const ZFSFloat pfac = fabs(p4I4) + fabs(p4J4) + fabs(p4K4);
          ZFSFloat fac = cfac*pfac;
          fac = min(1/64.0,fac*5.0);

          const ZFSFloat PLR= PL*(F1B2+ fac * MAL) + PR*(F1B2- fac * MAR);

          const ZFSFloat RHO_AL = RHOL*AL;
          const ZFSFloat RHO_AR = RHOR*AR;

          const ZFSFloat PLfRHOL = PL/RHOL;
          const ZFSFloat PRfRHOR = PR/RHOR;

          const ZFSFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
          const ZFSFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;

          const ZFSFloat RHOU = F1B2 * ( MALR * (RHO_AL+RHO_AR) + fabs(MALR) * (RHO_AL-RHO_AR) ) * metricLength;
          const ZFSFloat RHOU2 = F1B2*RHOU;
          // multiply by metric length to take surface area into account
          const ZFSFloat AbsRHO_U2 = fabs(RHOU2);

          // setup pressure factors, include metric terms (not normalized due to needed surface area)
          const ZFSFloat pFactor0 = surf[ 3 * dim + 0 ];
          const ZFSFloat pFactor1 = surf[ 3 * dim + 1 ];
          const ZFSFloat pFactor2 = surf[ 3 * dim + 2 ];
  
          flux[I+noCells*CV->RHO_U] = RHOU2 * ( UL + UR ) + AbsRHO_U2 * ( UL - UR ) + PLR * pFactor0;
          flux[I+noCells*CV->RHO_V] = RHOU2 * ( VL + VR ) + AbsRHO_U2 * ( VL - VR ) + PLR * pFactor1;
          flux[I+noCells*CV->RHO_W] = RHOU2 * ( WL + WR ) + AbsRHO_U2 * ( WL - WR ) + PLR * pFactor2;
          flux[I+noCells*CV->RHO_E] = RHOU2 * ( e0 + e1)  + AbsRHO_U2 * ( e0 - e1 ) + PLR * dxdtau;
          flux[I+noCells*CV->RHO]   = RHOU;
        }
      }
    }

    //FLUX BALANCE
    for(ZFSUint v=0; v<noVars; v++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
        for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
#if defined(ZFS_INTEL_COMPILER)
#pragma ivdep 
#pragma vector always
#endif
          for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
            const ZFSId I = cellIndex(i,j,k);
            const ZFSId IM1=I-IJK[dim];
            const ZFSUint offset = v*noCells;
            ZFSFloat *const RESTRICT rhs = ALIGNED_F(cellRhs+offset);
            rhs[I]+=flux[IM1+noCells*v]-flux[I+noCells*v];
          }
        }
      }
    }
  }
}

template void ZFSStrctrdBlck3D::Muscl_AusmLES_PTHRC<5>();
template void ZFSStrctrdBlck3D::Muscl_AusmLES_PTHRC<6>();
template void ZFSStrctrdBlck3D::Muscl_AusmLES_PTHRC<7>();


template<ZFSId noVars>
void ZFSStrctrdBlck3D::Muscl_AusmDV(){
  TRACE();

  //stencil identifier
  const ZFSId IJK[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};

  const ZFSFloat *const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const ZFSFloat *const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const ZFSFloat *const RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);
  const ZFSFloat *const *const RESTRICT vars= ALIGNED_F(m_cells->pvariables);
  ZFSFloat *const RESTRICT ds = ALIGNED_F(m_cells->ds);
  ZFSFloat *const RESTRICT flux = ALIGNED_F(m_cells->flux);
  ZFSFloat *const RESTRICT cellRhs= ALIGNED_MF(m_cells->rightHandSide[0]);

  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSFloat gamma = m_gamma;
  const ZFSFloat gammaMinusOne = gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  const ZFSInt noCellsI = m_nCells[2]-2;
  const ZFSInt noCellsJ = m_nCells[1]-2;
  const ZFSInt noCellsK = m_nCells[0]-2;

  const ZFSInt noCellsIP1 = m_nCells[2]-1;
  const ZFSInt noCellsJP1 = m_nCells[1]-1;
  const ZFSInt noCellsKP1 = m_nCells[0]-1;

  for(ZFSId dim=0; dim<nDim; dim++){
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(ZFSId k=0; k<noCellsKP1; k++) {
      for(ZFSId j=0; j<noCellsJP1; j++) {
        for(ZFSId i=0; i<noCellsIP1; i++) {
          const ZFSId I=cellIndex(i,j,k);
          const ZFSId IP1=I+IJK[dim];
          ds[I] = sqrt(POW2(x[IP1]-x[I]) + POW2(y[IP1]-y[I]) + POW2(z[IP1]-z[I]));
        }
      }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(ZFSId k=1; k<noCellsK; k++) {
      for(ZFSId j=1; j<noCellsJ; j++) {
#if defined(ZFS_INTEL_COMPILER)
#pragma ivdep 
#pragma vector always
#endif
        for(ZFSId i=1; i<noCellsI; i++) {
          const ZFSId I=cellIndex(i,j,k);
          const ZFSId IP1=I+IJK[dim];
          const ZFSId IM1=I-IJK[dim];
          const ZFSId IP2=I+2*IJK[dim];

          const ZFSFloat DS = ds[I];
          const ZFSFloat DSM1 = ds[IM1];
          const ZFSFloat DSP1 = ds[IP1];

          const ZFSFloat DSP=DS/POW2(DSP1+DS);
          const ZFSFloat DSM=DS/POW2(DSM1+DS);

          //unrolled the loop so the compiler
          //can optimize better
          const ZFSFloat DQU   = vars[PV->U][IP1]-vars[PV->U][I];
          const ZFSFloat DQPU = vars[PV->U][IP2]-vars[PV->U][IP1];
          const ZFSFloat DQMU = vars[PV->U][I]-vars[PV->U][IM1];
          ZFSFloat UL  = vars[PV->U][I]+DSM*(DSM1*DQU+DS*DQMU);
          ZFSFloat UR= vars[PV->U][IP1]-DSP*(DS*DQPU+DSP1*DQU);

          const ZFSFloat DQV   = vars[PV->V][IP1]-vars[PV->V][I];
          const ZFSFloat DQPV = vars[PV->V][IP2]-vars[PV->V][IP1];
          const ZFSFloat DQMV = vars[PV->V][I]-vars[PV->V][IM1];
          ZFSFloat VL = vars[PV->V][I]+DSM*(DSM1*DQV+DS*DQMV);
          ZFSFloat VR = vars[PV->V][IP1]-DSP*(DS*DQPV+DSP1*DQV);

          const ZFSFloat DQW   = vars[PV->W][IP1]-vars[PV->W][I];
          const ZFSFloat DQPW = vars[PV->W][IP2]-vars[PV->W][IP1];
          const ZFSFloat DQMW = vars[PV->W][I]-vars[PV->W][IM1];
          ZFSFloat WL = vars[PV->W][I]+DSM*(DSM1*DQW+DS*DQMW);
          ZFSFloat WR = vars[PV->W][IP1]-DSP*(DS*DQPW+DSP1*DQW);

          const ZFSFloat DQP   = vars[PV->P][IP1]-vars[PV->P][I];
          const ZFSFloat DQPP = vars[PV->P][IP2]-vars[PV->P][IP1];
          const ZFSFloat DQMP = vars[PV->P][I]-vars[PV->P][IM1];
          const ZFSFloat PL = vars[PV->P][I]+DSM*(DSM1*DQP+DS*DQMP);
          const ZFSFloat PR = vars[PV->P][IP1]-DSP*(DS*DQPP+DSP1*DQP);

          const ZFSFloat DQRHO   = vars[PV->RHO][IP1]-vars[PV->RHO][I];
          const ZFSFloat DQPRHO = vars[PV->RHO][IP2]-vars[PV->RHO][IP1];
          const ZFSFloat DQMRHO = vars[PV->RHO][I]-vars[PV->RHO][IM1];
          const ZFSFloat RHOL = vars[PV->RHO][I]+DSM*(DSM1*DQRHO+DS*DQMRHO);
          const ZFSFloat RHOR = vars[PV->RHO][IP1]-DSP*(DS*DQPRHO+DSP1*DQRHO);

          const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_cells->surfaceMetrics[ I ]);
          const ZFSFloat dxdtau = m_cells->dxt[dim][I];

          const ZFSFloat FRHOL = F1/RHOL;
          const ZFSFloat FRHOR = F1/RHOR;

          const ZFSFloat PLfRHOL = PL/RHOL;
          const ZFSFloat PRfRHOR = PR/RHOR;
          const ZFSFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
          const ZFSFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;


          // compute lenght of metric vector for normalization
          const ZFSFloat DGRAD = sqrt(POW2( surf[dim*3+0] ) + POW2(surf[dim*3+1]) + POW2(surf[dim*3+2]));
          const ZFSFloat FDGRAD = F1 / DGRAD;

          //scale by metric length to get velocity in the new basis (get normalized basis vectors)
          const ZFSFloat UUL = ((UL * surf[ dim*3+0 ] +
                                 VL * surf[ dim*3+1 ] +
                                 WL * surf[ dim*3+2 ]) - dxdtau) * FDGRAD;


          const ZFSFloat UUR = ((UR * surf[ dim*3+0 ] +
                                 VR * surf[ dim*3+1 ] +
                                 WR * surf[ dim*3+2 ]) - dxdtau) * FDGRAD;

          ZFSFloat AL = FRHOL*PL;
          ZFSFloat AR = FRHOR*PR;

          const ZFSFloat FALR = 2.0/(AL + AR);
          const ZFSFloat ALPHAL = AL*FALR;
          const ZFSFloat ALPHAR = AR*FALR;

          AL = sqrt(gamma*AL);
          AR = sqrt(gamma*AR);
          AL = zfsMAX(AL,AR);
          AR = AL;

          const ZFSFloat XMAL = UUL/AL;
          const ZFSFloat XMAR = UUR/AR;

          AL = AL*DGRAD;
          AR = AR*DGRAD;

          const ZFSFloat RHOAL = AL*RHOL;
          const ZFSFloat RHOAR = AR*RHOR;

          const ZFSFloat FDV = 0.3;
          const ZFSFloat DXDXEZ= m_cells->coordinates[0][IP1] - m_cells->coordinates[0][I];
          const ZFSFloat DYDXEZ= m_cells->coordinates[1][IP1] - m_cells->coordinates[1][I];
          const ZFSFloat DZDXEZ= m_cells->coordinates[2][IP1] - m_cells->coordinates[2][I];
          ZFSFloat SV = 2.0*DGRAD/(m_cells->cellJac[I]+m_cells->cellJac[IP1])*
            (FDV+(F1-FDV)*getPSI(I,dim));
          const ZFSFloat SV1 = F0*SV*DXDXEZ;
          const ZFSFloat SV2 = F0*SV*DYDXEZ;
          const ZFSFloat SV3 = F0*SV*DZDXEZ;

          const ZFSFloat XMAL1 = zfsMIN(F1,zfsMAX(-F1,XMAL));
          const ZFSFloat XMAR1 = zfsMIN(F1,zfsMAX(-F1,XMAR));

          ZFSFloat FXMA = F1B2*(XMAL1+fabs(XMAL1));
          const ZFSFloat XMALP = ALPHAL*(F1B4*POW2(XMAL1+F1)-FXMA)+FXMA+(zfsMAX(F1,XMAL)-F1);
          FXMA = F1B2*(XMAR1-fabs(XMAR1));
          const ZFSFloat XMARM = ALPHAR*(-F1B4*POW2(XMAR1-F1)-FXMA)+FXMA+(zfsMIN(-F1,XMAR)+F1);

          const ZFSFloat FLP = PL*((F2-XMAL1)*POW2(F1+XMAL1));
          const ZFSFloat FRP = PR*((F2+XMAR1)*POW2(F1-XMAR1));
          const ZFSFloat PLR = F1B4*(FLP+FRP);

          const ZFSFloat RHOUL = XMALP*RHOAL;
          const ZFSFloat RHOUR = XMARM*RHOAR;
          const ZFSFloat RHOU = RHOUL+RHOUR;
          const ZFSFloat RHOU2 = F1B2*RHOU;
          const ZFSFloat ARHOU2 = fabs(RHOU2);

          const ZFSFloat UUL2 = SV1*UUL;
          const ZFSFloat UUR2 = SV1*UUR;
          UL  = UL-UUL2;
          UR  = UR-UUR2;
          const ZFSFloat UUL3 = SV2*UUL;
          const ZFSFloat UUR3 = SV2*UUR;
          VL  = VL-UUL3;
          VR  = VR-UUR3;
          const ZFSFloat UUL4 = SV3*UUL;
          const ZFSFloat UUR4 = SV3*UUR;
          WL  = WL-UUL4;
          WR  = WR-UUR4;

          // setup pressure factors, include metric terms (not normalized due to needed surface area)
          const ZFSFloat pFactor0 = surf[ 3 * dim + 0 ];
          const ZFSFloat pFactor1 = surf[ 3 * dim + 1 ];
          const ZFSFloat pFactor2 = surf[ 3 * dim + 2 ];
  
          flux[I+noCells*CV->RHO_U] = RHOU2 * ( UL + UR) + ARHOU2 * ( UL - UR) + PLR*pFactor0 + RHOUL*UUL2+RHOUR*UUR2;
          flux[I+noCells*CV->RHO_V] = RHOU2 * ( VL + VR) + ARHOU2 * ( VL - VR) + PLR*pFactor1 + RHOUL*UUL3+RHOUR*UUR3;
          flux[I+noCells*CV->RHO_W] = RHOU2 * ( WL + WR) + ARHOU2 * ( WL - WR) + PLR*pFactor2 + RHOUL*UUL4+RHOUR*UUR4;
          flux[I+noCells*CV->RHO_E] = RHOU2 * ( e0 + e1) + ARHOU2 * ( e0 - e1) + PLR*dxdtau;
          flux[I+noCells*CV->RHO]   = RHOU;
        }
      }
    }

    //FLUX BALANCE
    for(ZFSUint v=0; v<noVars; v++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
        for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
#if defined(ZFS_INTEL_COMPILER)
#pragma ivdep 
#pragma vector always
#endif
          for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
            const ZFSId I = cellIndex(i,j,k);
            const ZFSId IM1=I-IJK[dim];
            const ZFSUint offset = v*noCells;
            ZFSFloat *const RESTRICT rhs = ALIGNED_F(cellRhs+offset);
            rhs[I]+=flux[IM1+noCells*v]-flux[I+noCells*v];
          }
        }
      }
    }
  }
}

template void ZFSStrctrdBlck3D::Muscl_AusmDV<5>();
template void ZFSStrctrdBlck3D::Muscl_AusmDV<6>();
template void ZFSStrctrdBlck3D::Muscl_AusmDV<7>();


template<ZFSStrctrdBlck3D::fluxmethod ausm,ZFSId noVars>
void ZFSStrctrdBlck3D::MusclStretched_(){
TRACE();

  //stencil identifier
  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSId IJK[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};
  const ZFSFloat *const RESTRICT cellVariables= ALIGNED_F(m_cells->pvariables[0]);
  const ZFSFloat *const RESTRICT cellLength= ALIGNED_F(m_cells->cellLength[0]);
  ZFSFloat *const RESTRICT cellRhs= ALIGNED_MF(m_cells->rightHandSide[0]);
  ZFSFloat *const RESTRICT qleft= ALIGNED_MF(m_QLeft);
  ZFSFloat *const RESTRICT qright= ALIGNED_MF(m_QRight);
  ZFSFloat *const RESTRICT flux = ALIGNED_F(m_cells->flux);
  /////////IMPORTANT PARAMETER
  //ZFSFloat epsi=F1;
  const ZFSFloat phi =F1;
  const ZFSFloat kappa=F0;//F1B3;
  /////////END IMPORTANT PARAMETER
  for(ZFSId dim=0; dim<nDim; dim++) {
    const ZFSUint dimOffset = dim*m_noStrctrdCells;
    const ZFSFloat *const RESTRICT length = ALIGNED_F(cellLength+dimOffset);

    for(ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers; k++) {
      for(ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers; j++) {
        for(ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers; i++) {
          const ZFSId I=cellIndex(i,j,k);
          const ZFSId IP1=I+IJK[dim];
          const ZFSId IM1=I-IJK[dim];
          const ZFSId IP2=I+2*IJK[dim];

          const ZFSFloat rp=(length[I]+length[IP1])/(F2*length[I]);
          const ZFSFloat rm=(length[I]+length[IM1])/(F2*length[I]);
          const ZFSFloat f=phi/(F2*(rp+rm));
          const ZFSFloat f1=(rm+kappa*phi)/rp;
          const ZFSFloat f2=(rp-kappa*phi)/rm;

          const ZFSFloat rp1=(length[IP1]+length[IP2])/(F2*length[IP1]);
          const ZFSFloat rm1=(length[IP1]+length[I])/(F2*length[IP1]);
          const ZFSFloat fa=phi/(F2*(rp1+rm1));
          const ZFSFloat fb=(rm1-kappa*phi)/rp1;
          const ZFSFloat fc=(rp1+kappa*phi)/rm1;

          for(ZFSUint v=0; v<noVars; v++) {
            const ZFSUint offset = v*m_noStrctrdCells;
            const ZFSFloat *const RESTRICT vars = ALIGNED_F(cellVariables+offset);
            //left variables
            const ZFSFloat DQ=(vars[IP1]-vars[I]);
            const ZFSFloat DQM1=(vars[I]-vars[IM1]);
            qleft[v] = vars[I]+f*(f1*DQ+f2*DQM1);

            //right variables
            const ZFSFloat DQP1=(vars[IP2]-vars[IP1]);
            const ZFSFloat DQ1=(vars[IP1]-vars[I]);
            qright[v]= vars[IP1]-fa*(fb*DQP1+fc*DQ1);
          }

          (this->*ausm)(m_QLeft, m_QRight, dim, I);
        }
      }
    }

    //FLUX BALANCE
    for(ZFSUint v=0; v<noVars; v++) {
      for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
        for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
          for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
            const ZFSId I=cellIndex(i,j,k);
            const ZFSId IM1=I-IJK[dim];
            const ZFSUint offset = v*noCells;
            ZFSFloat *const RESTRICT rhs = ALIGNED_F(cellRhs+offset);
            rhs[I]+=flux[IM1+noCells*v]-flux[I+noCells*v];
          }
        }
      }
    }
  }
}
//standard Ausm
template void ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES,5>();
template void ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES,6>();
template void ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES,7>();
//pthrc
template void ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES_PTHRC,5>();
template void ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES_PTHRC,6>();
template void ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmLES_PTHRC,7>();
//AusmDV
template void ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmDV,5>();
template void ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmDV,6>();
template void ZFSStrctrdBlck3D::MusclStretched_<&ZFSStrctrdBlck3D::AusmDV,7>();





void ZFSStrctrdBlck3D::Muscl(ZFSId zfsNotUsed(timerId))
{
  TRACE();
  if (m_movingGrid) {
    if(m_RKStep == 0) {
      saveGrid();
      saveCellJacobian();
    }

    moveGrid(false,false);

    //compute the volume fluxes
    RECORD_TIMER_START(m_tVolumeFlux);
    computeDxt();
    RECORD_TIMER_STOP(m_tVolumeFlux);
  }

  RECORD_TIMER_START(m_tConvectiveFlux);
  (this->*reconstructSurfaceData)();
  RECORD_TIMER_STOP(m_tConvectiveFlux);

  if(m_useSandpaperTrip) {
    applySandpaperTrip();
  }
}

void ZFSStrctrdBlck3D::Ausm()
{
  //Ausm routines have been moved and are called from inside Muscl (better performance)
}

void ZFSStrctrdBlck3D::computeVolumeForces()
{
  TRACE();
  
  for(ZFSId dim=0; dim<nDim; dim++) {
    for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
      m_cells->rightHandSide[CV->RHO_VV[dim]][cellId] += m_cells->variables[CV->RHO][cellId] * m_volumeForce[ dim ] * m_cells->cellJac[ cellId ];
      m_cells->rightHandSide[CV->RHO_E][cellId] += m_cells->variables[CV->RHO_VV[dim]][cellId] * m_volumeForce[ dim ] * m_cells->cellJac[ cellId ];
    }
  }
}

void ZFSStrctrdBlck3D::initSandpaperTrip() {
  TRACE();

  /*! \page propertyPage1
    \section tripXOrigin
    <code>ZFSInt ZFSStrctrdBlck::m_tripXOrigin </code>\n
    default = <code> 30.0 </code>\n \n
    Streamwise center position of the trip forcing.\n
    Possible values are:\n
    <ul>
    <li>Float <> 0.0</li>
    </ul>
    Keywords: <i>TRIP, BOUNDARYLAYER, STRCTRD</i>
  */
  m_tripXOrigin = 30.0;
  m_tripXOrigin = *(ZFSContext::getProperty("tripXOrigin", m_blockId, __CALLING_FUNCTION__, &m_tripXOrigin)->asFloat(0));

  /*! \page propertyPage1
    \section tripXLength
    <code>ZFSInt ZFSStrctrdBlck::m_tripXLength </code>\n
    default = <code> 1.0 </code>\n \n
    Streamwise extent of the trip forcing.\n
    Possible values are:\n
    <ul>
    <li>Float > 0.0</li>
    </ul>
    Keywords: <i>TRIP, BOUNDARYLAYER, STRCTRD</i>
  */
  m_tripXLength = 1.0;
  m_tripXLength = *(ZFSContext::getProperty("tripXLength", m_blockId, __CALLING_FUNCTION__, &m_tripXLength)->asFloat(0));

  /*! \page propertyPage1
    \section tripYOrigin
    <code>ZFSInt ZFSStrctrdBlck::m_tripYOrigin </code>\n
    default = <code> 0.5 </code>\n \n
    Wall-normal center position of the trip forcing.\n
    Possible values are:\n
    <ul>
    <li>Float > 0.0</li>
    </ul>
    Keywords: <i>TRIP, BOUNDARYLAYER, STRCTRD</i>
  */
  m_tripYOrigin = 0.5;
  m_tripYOrigin = *(ZFSContext::getProperty("tripYOrigin", m_blockId, __CALLING_FUNCTION__, &m_tripYOrigin)->asFloat(0));

  /*! \page propertyPage1
    \section tripYHeight
    <code>ZFSInt ZFSStrctrdBlck::m_tripYHeight </code>\n
    default = <code> 0.5 </code>\n \n
    Wall-normal extent of the trip forcing.\n
    Possible values are:\n
    <ul>
    <li>Float > 0.0</li>
    </ul>
    Keywords: <i>TRIP, BOUNDARYLAYER, STRCTRD</i>
  */
  m_tripYHeight = 0.5;
  m_tripYHeight = *(ZFSContext::getProperty("tripYHeight", m_blockId, __CALLING_FUNCTION__, &m_tripYHeight)->asFloat(0));

  /*! \page propertyPage1
    \section tripMaxAmpSteady
    <code>ZFSInt ZFSStrctrdBlck::m_tripMaxAmpSteady </code>\n
    default = <code> 0.0 </code>\n \n
    Strength of the steady forcing amplitude.\n
    Possible values are:\n
    <ul>
    <li>Float > 0.0</li>
    </ul>
    Keywords: <i>TRIP, BOUNDARYLAYER, STRCTRD</i>
  */
  m_tripMaxAmpSteady = 0.0;
  m_tripMaxAmpSteady = *(ZFSContext::getProperty("tripMaxAmpSteady", m_blockId, __CALLING_FUNCTION__, &m_tripMaxAmpSteady)->asFloat(0));

  /*! \page propertyPage1
    \section tripMaxAmpFluc
    <code>ZFSInt ZFSStrctrdBlck::m_tripMaxAmpFluc </code>\n
    default = <code> 0.005 </code>\n \n
    Strength of the fluctuating forcing amplitude.\n
    Possible values are:\n
    <ul>
    <li>Float > 0.0</li>
    </ul>
    Keywords: <i>TRIP, BOUNDARYLAYER, STRCTRD</i>
  */
  m_tripMaxAmpFluc = 0.005;
  m_tripMaxAmpFluc = *(ZFSContext::getProperty("tripMaxAmpFluc", m_blockId, __CALLING_FUNCTION__, &m_tripMaxAmpFluc)->asFloat(0));

  /*! \page propertyPage1
    \section tripNoModes
    <code>ZFSInt ZFSStrctrdBlck::m_tripNoModes </code>\n
    default = <code> 30 </code>\n \n
    Number of Fourier modes to use for\n
    the trip forcing.\n
    Possible values are:\n
    <ul>
    <li>Int > 0</li>
    </ul>
    Keywords: <i>TRIP, BOUNDARYLAYER, STRCTRD</i>
  */
  m_tripNoModes = 30;
  m_tripNoModes = *(ZFSContext::getProperty("tripNoModes", m_blockId, __CALLING_FUNCTION__, &m_tripNoModes)->asInt(0));

  /*! \page propertyPage1
    \section tripDeltaTime
    <code>ZFSInt ZFSStrctrdBlck::m_tripDeltaTime </code>\n
    default = <code> 2.0 </code>\n \n
    Delta t of the tripping, i.e. the\n
    time step to change the Fourier coeffients.\n
    Possible values are:\n
    <ul>
    <li>Float > 0.0</li>
    </ul>
    Keywords: <i>TRIP, BOUNDARYLAYER, STRCTRD</i>
  */
  m_tripDeltaTime = 2.0;
  m_tripDeltaTime = *(ZFSContext::getProperty("tripDeltaTime", m_blockId, __CALLING_FUNCTION__, &m_tripDeltaTime)->asFloat(0));


  m_tripTimeStep = (ZFSInt)(m_time/m_tripDeltaTime);
  m_tripSeed = 70;
  srand(m_tripSeed);
  m_tripDomainWidth = 0.0;
  m_tripNoCells = m_nActiveCells[0];

  zfsAlloc(m_tripCoords, m_tripNoCells, "m_tripCoords", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_tripG, m_tripNoCells, "m_tripG", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_tripH1, m_tripNoCells, "m_tripH1", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_tripH2, m_tripNoCells, "m_tripH2", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_tripModesG, 2*m_tripNoModes, "m_tripModesG", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_tripModesH1, 2*m_tripNoModes, "m_tripModesH1", F0, __CALLING_FUNCTION__);
  zfsAlloc(m_tripModesH2, 2*m_tripNoModes, "m_tripModesH2", F0, __CALLING_FUNCTION__);
  
  ZFSFloat localDomainWidth = m_cells->coordinates[2][cellIndex(0,0,m_nCells[0]-2)];
  MPI_Allreduce(&localDomainWidth, &m_tripDomainWidth, 1, MPI_DOUBLE, MPI_MAX, m_zfsStrctrdComm);

  zfs_log << "=================================================" << endl
          << "           SANDPAPER TRIP PROPERTIES             " << endl
          << "=================================================" << endl
          << "tripXOrigin: " << m_tripXOrigin << endl
          << "tripXLength: " << m_tripXLength << endl
          << "tripYOrigin: " << m_tripYOrigin << endl
          << "tripYHeight: " << m_tripYHeight << endl
          << "tripMaxAmpSteady: " << m_tripMaxAmpSteady << endl
          << "tripMaxAmpFluc: " << m_tripMaxAmpFluc << endl
          << "tripNoModes: " << m_tripNoModes << endl
          << "tripDeltaTime: " << m_tripDeltaTime << endl
          << "tripTimeStep: " << m_tripTimeStep << endl
          << "tripDomainWidth: " << m_tripDomainWidth << endl
          << "=================================================" <<endl;

  for(ZFSId k=0; k<m_tripNoCells; k++) {
    m_tripCoords[k] = m_cells->coordinates[2][cellIndex(0,0,m_noGhostLayers+k)];
  }

  if(m_restart) {
      stringstream tripPath;
      tripPath << "/trip";
      ZFSInt dummyOffset = 0;
      
      if(domainId()==0) {
        stringstream restartFileName;
        ZFSString restartFile = *(ZFSContext::getProperty("restartVariablesFileName", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL )->asString(0));
        restartFileName << outputDir() << restartFile;
        ZFSInt fid=io_openfile("hdf5" ,(restartFileName.str()).c_str(),"collective", MPI_COMM_SELF);

        ZFSInt dataSize = m_tripNoModes*2;
        io_read_ddataset_part1d1(fid, (tripPath.str()).c_str(), "tripModesG", 1, &dummyOffset, &dataSize, m_tripModesG);
        io_read_ddataset_part1d1(fid, (tripPath.str()).c_str(), "tripModesH1", 1, &dummyOffset, &dataSize, m_tripModesH1);
        io_read_ddataset_part1d1(fid, (tripPath.str()).c_str(), "tripModesH2", 1, &dummyOffset, &dataSize, m_tripModesH2);
        io_closefile(fid);
      }

      MPI_Bcast(m_tripModesG, 2*m_tripNoModes, MPI_DOUBLE, 0, m_zfsStrctrdComm);
      MPI_Bcast(m_tripModesH1, 2*m_tripNoModes, MPI_DOUBLE, 0, m_zfsStrctrdComm);
      MPI_Bcast(m_tripModesH2, 2*m_tripNoModes, MPI_DOUBLE, 0, m_zfsStrctrdComm);
  } else {
    tripFourierCoefficients(m_tripModesG, m_tripNoModes);
    tripFourierCoefficients(m_tripModesH1, m_tripNoModes);
    tripFourierCoefficients(m_tripModesH2, m_tripNoModes);
  }

  tripForceCoefficients(m_tripModesG, m_tripG, m_tripCoords, m_tripNoCells, m_tripNoModes, m_tripDomainWidth);
  tripForceCoefficients(m_tripModesH1, m_tripH1, m_tripCoords, m_tripNoCells, m_tripNoModes, m_tripDomainWidth);
  tripForceCoefficients(m_tripModesH2, m_tripH2, m_tripCoords, m_tripNoCells, m_tripNoModes, m_tripDomainWidth);
}

void ZFSStrctrdBlck3D::applySandpaperTrip() {
  TRACE();
  const ZFSFloat t = m_time + m_timeStep*m_RKalpha[m_RKStep];
  const ZFSInt tripTime = (ZFSInt)(t/m_tripDeltaTime);
  const ZFSFloat p = t/m_tripDeltaTime - tripTime;
  const ZFSFloat b = 3*pow(p,2) - 2*pow(p,3);
  
  if(tripTime > m_tripTimeStep) {
    m_tripTimeStep = tripTime;

    //copy old values from H2 to H1
    for(ZFSId k=0; k<m_tripNoCells; k++) {
      m_tripH1[k] = m_tripH2[k];
    }

    //also copy the old mode coefficients
    for(ZFSId n=0; n<2*m_tripNoModes; n++) {
      m_tripModesH1[n] = m_tripModesH2[n];
    }

    //compute new fourier coefficients
    tripFourierCoefficients(m_tripModesH2, m_tripNoModes);
    tripForceCoefficients(m_tripModesH2, m_tripH2, m_tripCoords, m_tripNoCells, m_tripNoModes, m_tripDomainWidth);
  }
    
  for(ZFSId k=0; k<m_nActiveCells[0]; k++) {
    const ZFSFloat forceStrength = (m_tripMaxAmpSteady*m_tripG[k] + 
                                    m_tripMaxAmpFluc*((1.0-b)*m_tripH1[k] + b*m_tripH2[k]));

    for(ZFSId j=0; j<m_nActiveCells[1]; j++) {
      for(ZFSId i=0; i<m_nActiveCells[2];i++) {
        const ZFSId cellId = cellIndex(i+m_noGhostLayers,j+m_noGhostLayers,k+m_noGhostLayers);
        const ZFSFloat x = m_cells->coordinates[0][cellId];
        const ZFSFloat y = m_cells->coordinates[1][cellId];

        ZFSFloat force = 0.0; 
        if(x > m_tripXOrigin-m_tripXLength &&
           x < m_tripXOrigin+m_tripXLength &&
           y > m_tripYOrigin-m_tripYHeight &&
           y < m_tripYOrigin+m_tripYHeight) {
          force = exp(POW2((x-m_tripXOrigin)/m_tripXLength) + POW2((y-m_tripYOrigin)/m_tripYHeight))*forceStrength;
        }

        m_cells->rightHandSide[CV->RHO_V][cellId] += m_cells->variables[CV->RHO][cellId] * force * m_cells->cellJac[ cellId ];
        m_cells->rightHandSide[CV->RHO_E][cellId] += m_cells->variables[CV->RHO_V][cellId] * force * m_cells->cellJac[ cellId ];
      }
    }
  }
}

void ZFSStrctrdBlck3D::tripForceCoefficients(ZFSFloat* modes, ZFSFloat* forceCoef,
                                             ZFSFloat* coords, ZFSInt noCells, ZFSInt noModes,
                                             ZFSFloat maxWaveLength) {
  ZFSFloat maxLocalValue = 0.0;
  ZFSFloat* ak = &modes[0];
  ZFSFloat* phik = &modes[noModes];
  
  for(ZFSId k=0; k<noCells; k++) {
    const ZFSFloat z = coords[k];
    forceCoef[k] = ak[0]*F1B2;
    for(ZFSId n=1; n<noModes; n++) {
      forceCoef[k] += ak[n]*cos(n*z/maxWaveLength*2*PI - phik[n]);
    }

    maxLocalValue = zfsMAX(maxLocalValue, fabs(forceCoef[k]));
  }

  ZFSFloat maxGlobalValue = 0.0;
  MPI_Allreduce(&maxLocalValue, &maxGlobalValue, 1, MPI_DOUBLE, MPI_MAX, m_zfsStrctrdComm);

  //normalize the series
  for(ZFSId k=0; k<noCells; k++) {
    forceCoef[k] = forceCoef[k]/maxGlobalValue;
  }
}

void ZFSStrctrdBlck3D::tripFourierCoefficients(ZFSFloat* modes, ZFSInt noModes) {
  ZFSFloat* ak = &modes[0];
  ZFSFloat* phik = &modes[noModes];
  if(domainId()==0) {
    for(ZFSId n=0; n<noModes; n++) {
      ak[n] = rand()/double(RAND_MAX);
      phik[n] = rand()/double(RAND_MAX)*2*PI;
    }
  }
  
  MPI_Bcast(&ak[0], noModes, MPI_DOUBLE, 0, m_zfsStrctrdComm);
  MPI_Bcast(&phik[0], noModes, MPI_DOUBLE, 0, m_zfsStrctrdComm);
}



void ZFSStrctrdBlck3D::computeCellCentreCoordinates()
{
  TRACE();
  //function to compute the coordinates at cell centre
  //do it over I, J, K loop but change to one array
  for(ZFSId k=0; k<m_nCells[0]; k++)
  {
    for(ZFSId j=0; j<m_nCells[1]; j++)
    {
      for(ZFSId i=0; i<m_nCells[2];i++)
      {
        ZFSId pointId = i+(j+k*m_nPoints[1])*m_nPoints[2];
        ZFSId IJK = pointId;
        ZFSId IP1JK= pointId+1;
        ZFSId IJP1K= pointId+m_nPoints[2];
        ZFSId IP1JP1K= IJP1K+1;
        ZFSId IJKP1= pointId+m_nPoints[2]*m_nPoints[1];
        ZFSId IP1JKP1= IJKP1+1;
        ZFSId IJP1KP1= pointId+m_nPoints[2]+m_nPoints[2]*m_nPoints[1];
        ZFSId IP1JP1KP1= IJP1KP1+1;
        ZFSId cellId = i+(j+k*m_nCells[1])*m_nCells[2];
              for(ZFSId dim=0; dim<nDim; dim++)
        {
          //average the coordinates for cell centre data
          m_cells->coordinates[dim][cellId]=F1B8*(m_coordinates[dim][IJK]+m_coordinates[dim][IP1JK]+m_coordinates[dim][IJP1K]+m_coordinates[dim][IP1JP1K]+m_coordinates[dim][IJKP1]+m_coordinates[dim][IP1JKP1]+m_coordinates[dim][IJP1KP1]+m_coordinates[dim][IP1JP1KP1]);
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::distributeFluxToCells()
{
  TRACE();
}

void ZFSStrctrdBlck3D::computeCumulativeAverage(ZFSBool forceReset)
{
  TRACE();
  if(m_RKStep==0 && m_zoneType == "LES") {
    if(globalTimeStep == 0 || 
       forceReset) {
      // cout << "///////// Starting averaging, resetting AVG variables ////////// domainId(): " << domainId() << " m_zoneType: " << m_zoneType << endl;
      for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
	m_cells->fq[FQ->AVG_RHO][cellId] = m_cells->pvariables[PV->RHO][cellId];
	m_cells->fq[FQ->AVG_U][cellId] = m_cells->pvariables[PV->U][cellId];
	m_cells->fq[FQ->AVG_V][cellId] = m_cells->pvariables[PV->V][cellId];
	m_cells->fq[FQ->AVG_W][cellId] = m_cells->pvariables[PV->W][cellId];
	m_cells->fq[FQ->AVG_P][cellId] = m_cells->pvariables[PV->P][cellId];
      }
    } else {
      const ZFSFloat timeFacRho = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
      const ZFSFloat timeFacU = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
      const ZFSFloat timeFacV = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
      const ZFSFloat timeFacW = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
      const ZFSFloat timeFacE = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
      // cout<<"timeFacRho:"<<timeFacRho<<" and timeFacU:"<<timeFacU<<"and timeFacV:"<<timeFacV<<"and m_timeStep:"<<m_timeStep<<"sqrt(PV->TInfinity):"<<sqrt(PV->TInfinity)<<endl;

      //! do the time average of the flow variables for LES and store them.
      for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
	//Exponential averaging
	m_cells->fq[FQ->AVG_RHO][cellId] = timeFacRho* m_cells->pvariables[PV->RHO][cellId] + (1.0 - timeFacRho)*m_cells->fq[FQ->AVG_RHO][cellId];
	m_cells->fq[FQ->AVG_U][cellId] = timeFacU* m_cells->pvariables[PV->U][cellId] + (1.0 - timeFacU)*m_cells->fq[FQ->AVG_U][cellId];
	m_cells->fq[FQ->AVG_V][cellId] = timeFacV* m_cells->pvariables[PV->V][cellId] + (1.0 - timeFacV)*m_cells->fq[FQ->AVG_V][cellId];
	m_cells->fq[FQ->AVG_W][cellId] = timeFacW* m_cells->pvariables[PV->W][cellId] + (1.0 - timeFacW)*m_cells->fq[FQ->AVG_W][cellId];
	m_cells->fq[FQ->AVG_P][cellId] = timeFacE* m_cells->pvariables[PV->P][cellId] + (1.0 - timeFacE)*m_cells->fq[FQ->AVG_P][cellId];
      }
      
      	// //! COMPUTE THE FLUCTUATIONS AND STORE THEM
	// // ZFSFloat fRHO = 0.0;
	// // ZFSFloat fRHO_avg = 0.0;
	// const ZFSFloat timeFacUU = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
	// const ZFSFloat timeFacVV = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
	// const ZFSFloat timeFacWW = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
	// const ZFSFloat timeFacUV = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
	// const ZFSFloat timeFacUW = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
	// const ZFSFloat timeFacVW = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
	

	// for (ZFSId k=0; k<m_nActiveCells[0]; k++)
	//   {
	//     for (ZFSId j=0; j<m_nActiveCells[1]; j++)
	//       {
	// 	for (ZFSId i=0; i<m_nActiveCells[2]; i++)
	// 	  {			   
	// 	    ZFSId localCellId3D= i+m_noGhostLayers+((j+m_noGhostLayers)+(k+m_noGhostLayers)*m_nCells[1])*m_nCells[2];
	// 	    // COMPUTING the fluctuating part of each velocity components 
	// 	    m_cells->fq[FQ->FLUC_U][localCellId3D]=m_cells->pvariables[PV->U][localCellId3D] - m_cells->fq[FQ->AVG_U][localCellId3D];
	// 	    m_cells->fq[FQ->FLUC_V][localCellId3D]=m_cells->pvariables[PV->V][localCellId3D] - m_cells->fq[FQ->AVG_V][localCellId3D];
	// 	    m_cells->fq[FQ->FLUC_W][localCellId3D]=m_cells->pvariables[PV->W][localCellId3D] - m_cells->fq[FQ->AVG_W][localCellId3D];

	// 	    // computing the momemts with time averaging 1 
	// 	    m_cells->fq[FQ->FLUC_UU][localCellId3D]=timeFacUU*(m_cells->fq[FQ->FLUC_U][localCellId3D]*m_cells->fq[FQ->FLUC_U][localCellId3D])+(F1-timeFacUU)*m_cells->fq[FQ->FLUC_UU][localCellId3D];
	// 	    m_cells->fq[FQ->FLUC_VV][localCellId3D]=timeFacVV*(m_cells->fq[FQ->FLUC_V][localCellId3D]*m_cells->fq[FQ->FLUC_V][localCellId3D])+(F1-timeFacVV)*m_cells->fq[FQ->FLUC_VV][localCellId3D];
	// 	    m_cells->fq[FQ->FLUC_WW][localCellId3D]=timeFacWW*(m_cells->fq[FQ->FLUC_W][localCellId3D]*m_cells->fq[FQ->FLUC_W][localCellId3D])+(F1-timeFacWW)*m_cells->fq[FQ->FLUC_WW][localCellId3D];
	// 	    m_cells->fq[FQ->FLUC_UV][localCellId3D]=timeFacUV*(m_cells->fq[FQ->FLUC_U][localCellId3D]*m_cells->fq[FQ->FLUC_V][localCellId3D])+(F1-timeFacUV)*m_cells->fq[FQ->FLUC_UV][localCellId3D];
	// 	    m_cells->fq[FQ->FLUC_UW][localCellId3D]=timeFacUW*(m_cells->fq[FQ->FLUC_U][localCellId3D]*m_cells->fq[FQ->FLUC_W][localCellId3D])+(F1-timeFacUW)*m_cells->fq[FQ->FLUC_UW][localCellId3D];
	// 	    m_cells->fq[FQ->FLUC_VW][localCellId3D]=timeFacVW*(m_cells->fq[FQ->FLUC_V][localCellId3D]*m_cells->fq[FQ->FLUC_W][localCellId3D])+(F1-timeFacVW)*m_cells->fq[FQ->FLUC_VW][localCellId3D];

	// 	  }
	//       }
	//   }

      
    }
  }
}

void ZFSStrctrdBlck3D::spanwiseAvgZonal() {
  if(m_zonal == true)
    {
      if(m_zoneType=="LES"){ 	

	const ZFSInt LESBlock = m_inputBlockId; 

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////Spanwise averaging/////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	ZFSFloat totalNoCellsIJ=(m_totalGridBlockDim[LESBlock][2]-1)*(m_totalGridBlockDim[LESBlock][1]-1);

	ZFSFloatScratchSpace localSpannwiseVars(totalNoCellsIJ,__CALLING_FUNCTION__,"localSpannwiseVars");
	ZFSFloatScratchSpace globalSpannwiseVars(totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseVars");

	ZFSFloatScratchSpace globalSpannwiseVarsSumIJ(totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseVarsSumIJ");
	ZFSFloatScratchSpace globalSpannwiseVarsSumK(totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseVarsSumK");
	ZFSFloatScratchSpace globalSpannwiseAveragedVars(totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseAveragedVars");

	
     
	for(ZFSId k=0; k<m_nActiveCells[0]; k++){
	  for(ZFSId j=0; j<m_nActiveCells[1]; j++){
	    for(ZFSId i=0; i<m_nActiveCells[2]; i++){
	      ZFSId localCellId3D= i+m_noGhostLayers+((j+m_noGhostLayers)+(k+m_noGhostLayers)*m_nCells[1])*m_nCells[2];
	      ZFSId globalCellId2D =(i+m_nOffsetCells[2])+(j+m_nOffsetCells[1])*(m_totalGridBlockDim[LESBlock][2]-1);
	      
		localSpannwiseVars(globalCellId2D) +=m_cells->fq[FQ->AVG_P][localCellId3D]; // adding the vars in spanwise direction
	      
	    }
	  }
	}

	MPI_Allreduce(&localSpannwiseVars(0),&globalSpannwiseVars(0),totalNoCellsIJ,MPI_DOUBLE,MPI_SUM,m_commZonal[LESBlock]);
	for(ZFSId i=0; i<totalNoCellsIJ; i++){
	    globalSpannwiseVars(i)/=(m_totalGridBlockDim[LESBlock][0]-1);
	}

	
	for (ZFSId k=0; k<m_nActiveCells[0]; k++)
	  {
	    for (ZFSId j=0; j<m_nActiveCells[1]; j++)
	      {
		for (ZFSId i=0; i<m_nActiveCells[2]; i++)
		  {
		    ZFSId localCellId3D= i+m_noGhostLayers+((j+m_noGhostLayers)+(k+m_noGhostLayers)*m_nCells[1])*m_nCells[2];
		    ZFSId globalCellId2D =(i+m_nOffsetCells[2])+(j+m_nOffsetCells[1])*(m_totalGridBlockDim[LESBlock][2]-1);

		    m_cells->fq[FQ->AVG_P][localCellId3D]   = globalSpannwiseVars(globalCellId2D);
		   
		  }
	      }
	  }
      }
    }
}




// void ZFSStrctrdBlck3D::reconstructTurbulentVariables() {
//   if(m_zonal == true)
//     {
//       if(m_zoneType=="LES"){ 
// 	const ZFSInt LESBlock = m_inputBlockId; 
// 	// cout<<"m_totalGridBlockDim[blocknumber][0]:"<<m_totalGridBlockDim[LESBlock][0]<<" m_totalGridBlockDim[blocknumber][1]"<<m_totalGridBlockDim[LESBlock][1]<<" m_totalGridBlockDim[blocknumber][2]"<<m_totalGridBlockDim[LESBlock][2]<<endl;



// 	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	/////////////////////////////////////Spanwise averaging/////////////////////////////////////////////////////////////
// 	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 	ZFSFloat totalNoCellsIJ=(m_totalGridBlockDim[LESBlock][2]-1)*(m_totalGridBlockDim[LESBlock][1]-1);

// 	ZFSFloatScratchSpace localSpannwiseVars(PV->noVariables,totalNoCellsIJ,__CALLING_FUNCTION__,"localSpannwiseVars");
// 	ZFSFloatScratchSpace globalSpannwiseVars(PV->noVariables,totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseVars");

// 	ZFSFloatScratchSpace globalSpannwiseVarsSumIJ(PV->noVariables,totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseVarsSumIJ");
// 	ZFSFloatScratchSpace globalSpannwiseVarsSumK(PV->noVariables,totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseVarsSumK");
// 	ZFSFloatScratchSpace globalSpannwiseAveragedVars(PV->noVariables,totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseAveragedVars");

	
     
// 	// cout<<"PV->noVariables*totalNoCellsIJ:"<<PV->noVariables*totalNoCellsIJ<<endl;
// 	// cout<<"totalNoCellsIJ:"<<totalNoCellsIJ<<endl;
// 	// cout<<"m_nCells[0]:"<<m_nCells[0]<<"m_nCells[1]:"<<m_nCells[1]<<" m_nCells[2]:"<<m_nCells[2]<<" m_nActiveCells[0]:"<<m_nActiveCells[0]<<" m_nActiveCells[1]:"<<m_nActiveCells[1]<<" m_nActiveCells[2]:"<<m_nActiveCells[2]<<" m_nOffsetCells[0]:"<<m_nOffsetCells[0]<<" m_nOffsetCells[1]:"<<m_nOffsetCells[1]<<" m_nOffsetCells[2]"<<m_nOffsetCells[2]<<"domainId:"<<domainId()<<endl;
// 	// for(ZFSId k=m_nOffsetCells[0]; k<m_nOffsetCells[0]+m_nActiveCells[0]; k++){
// 	//   for(ZFSId j=m_nOffsetCells[0]; j<m_nOffsetCells[1]+m_nActiveCells[1]; j++){
// 	//     for(ZFSId i=m_nOffsetCells[0]; i<m_nOffsetCells[2]+m_nActiveCells[2]; i++){
// 	//       ZFSIad cellId3D= i +(j+k*m_nActiveCells[1])*m_nActiveCells[2];
// 	//       ZFSId cellId2D= i +(j*m_nActiveCells[2]);

// 	for(ZFSId k=0; k<m_nActiveCells[0]; k++){
// 	  for(ZFSId j=0; j<m_nActiveCells[1]; j++){
// 	    for(ZFSId i=0; i<m_nActiveCells[2]; i++){
// 	      ZFSId localCellId3D= i+m_noGhostLayers+((j+m_noGhostLayers)+(k+m_noGhostLayers)*m_nCells[1])*m_nCells[2];
// 	      ZFSId globalCellId2D =(i+m_nOffsetCells[2])+(j+m_nOffsetCells[1])*(m_totalGridBlockDim[LESBlock][2]-1);
// 	      for(ZFSId vars=0; vars<PV->noVariables; vars++){
// 		localSpannwiseVars(vars,globalCellId2D) +=m_cells->pvariables[vars][localCellId3D]; // adding the vars in spanwise direction
// 	      }
// 	    }
// 	  }
// 	}

// 	MPI_Allreduce(&localSpannwiseVars(0,0),&globalSpannwiseVars(0,0),PV->noVariables*totalNoCellsIJ,MPI_DOUBLE,MPI_SUM,m_commZonal[LESBlock]);
// 	for(ZFSId i=0; i<totalNoCellsIJ; i++){
// 	  for(ZFSId vars=0; vars<PV->noVariables; vars++){
// 	    globalSpannwiseVars(vars,i)/=(m_totalGridBlockDim[LESBlock][0]-1);
// 	  }
// 	}

// 	// if(domainId()==m_commZonalRootGlobal[LESBlock]){
// 	//   for(ZFSId i=0; i<totalNoCellsIJ; i++){
// 	//     for(ZFSId vars=0; vars<PV->noVariables; vars++){
// 	// 	// cout<<"localSpannwiseVars["<<vars<<"]["<<i<<"]:"<<localSpannwiseVars[vars][i]<<endl;
// 	// 	cout<<"globalSpannwiseVars["<<vars<<"]["<<i<<"]:"<<globalSpannwiseVars[vars][i]<<" domainId:"<<domainId()<<endl; 
// 	//     }
// 	//   }
// 	// }
	
// 	for (ZFSId k=0; k<m_nActiveCells[0]; k++)
// 	  {
// 	    for (ZFSId j=0; j<m_nActiveCells[1]; j++)
// 	      {
// 		for (ZFSId i=0; i<m_nActiveCells[2]; i++)
// 		  {
// 		    ZFSId localCellId3D= i+m_noGhostLayers+((j+m_noGhostLayers)+(k+m_noGhostLayers)*m_nCells[1])*m_nCells[2];
// 		    ZFSId globalCellId2D =(i+m_nOffsetCells[2])+(j+m_nOffsetCells[1])*(m_totalGridBlockDim[LESBlock][2]-1);

// 		    m_cells->fq[FQ->AVG_U][localCellId3D]   = globalSpannwiseVars(0,globalCellId2D);
// 		    m_cells->fq[FQ->AVG_V][localCellId3D] = globalSpannwiseVars(1,globalCellId2D);
// 		    m_cells->fq[FQ->AVG_W][localCellId3D] = globalSpannwiseVars(2,globalCellId2D);
// 		    m_cells->fq[FQ->AVG_RHO][localCellId3D] = globalSpannwiseVars(3,globalCellId2D);
// 		    m_cells->fq[FQ->AVG_P][localCellId3D] = globalSpannwiseVars(4,globalCellId2D);

// 		  }
// 	      }
// 	  }


// 	// //! COMPUTE THE FLUCTUATIONS AND STORE THEM
// 	// // ZFSFloat fRHO = 0.0;
// 	// // ZFSFloat fRHO_avg = 0.0;
// 	// const ZFSFloat timeFacUU = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
// 	// const ZFSFloat timeFacVV = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
// 	// const ZFSFloat timeFacWW = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
// 	// const ZFSFloat timeFacUV = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
// 	// const ZFSFloat timeFacUW = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
// 	// const ZFSFloat timeFacVW = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
	

// 	// for (ZFSId k=0; k<m_nActiveCells[0]; k++)
// 	//   {
// 	//     for (ZFSId j=0; j<m_nActiveCells[1]; j++)
// 	//       {
// 	// 	for (ZFSId i=0; i<m_nActiveCells[2]; i++)
// 	// 	  {			   
// 	// 	    ZFSId localCellId3D= i+m_noGhostLayers+((j+m_noGhostLayers)+(k+m_noGhostLayers)*m_nCells[1])*m_nCells[2];
// 	// 	    // COMPUTING the fluctuating part of each velocity components 
// 	// 	    m_cells->fq[FQ->FLUC_U][localCellId3D]=m_cells->pvariables[PV->U][localCellId3D] - m_cells->fq[FQ->AVG_U][localCellId3D];
// 	// 	    m_cells->fq[FQ->FLUC_V][localCellId3D]=m_cells->pvariables[PV->V][localCellId3D] - m_cells->fq[FQ->AVG_V][localCellId3D];
// 	// 	    m_cells->fq[FQ->FLUC_W][localCellId3D]=m_cells->pvariables[PV->W][localCellId3D] - m_cells->fq[FQ->AVG_W][localCellId3D];

// 	// 	    // computing the momemts with time averaging 1 
// 	// 	    m_cells->fq[FQ->FLUC_UU][localCellId3D]=timeFacUU*(m_cells->fq[FQ->FLUC_U][localCellId3D]*m_cells->fq[FQ->FLUC_U][localCellId3D])+(F1-timeFacUU)*m_cells->fq[FQ->FLUC_UU][localCellId3D];
// 	// 	    m_cells->fq[FQ->FLUC_VV][localCellId3D]=timeFacVV*(m_cells->fq[FQ->FLUC_V][localCellId3D]*m_cells->fq[FQ->FLUC_V][localCellId3D])+(F1-timeFacVV)*m_cells->fq[FQ->FLUC_VV][localCellId3D];
// 	// 	    m_cells->fq[FQ->FLUC_WW][localCellId3D]=timeFacWW*(m_cells->fq[FQ->FLUC_W][localCellId3D]*m_cells->fq[FQ->FLUC_W][localCellId3D])+(F1-timeFacWW)*m_cells->fq[FQ->FLUC_WW][localCellId3D];
// 	// 	    m_cells->fq[FQ->FLUC_UV][localCellId3D]=timeFacUV*(m_cells->fq[FQ->FLUC_U][localCellId3D]*m_cells->fq[FQ->FLUC_V][localCellId3D])+(F1-timeFacUV)*m_cells->fq[FQ->FLUC_UV][localCellId3D];
// 	// 	    m_cells->fq[FQ->FLUC_UW][localCellId3D]=timeFacUW*(m_cells->fq[FQ->FLUC_U][localCellId3D]*m_cells->fq[FQ->FLUC_W][localCellId3D])+(F1-timeFacUW)*m_cells->fq[FQ->FLUC_UW][localCellId3D];
// 	// 	    m_cells->fq[FQ->FLUC_VW][localCellId3D]=timeFacVW*(m_cells->fq[FQ->FLUC_V][localCellId3D]*m_cells->fq[FQ->FLUC_W][localCellId3D])+(F1-timeFacVW)*m_cells->fq[FQ->FLUC_VW][localCellId3D];

// 	// 	  }
// 	//       }
// 	//   }





// 	ZFSInt noAveragedVars = 6;
// 	ZFSFloatScratchSpace localSpannwiseVarsFLUC(noAveragedVars,totalNoCellsIJ,__CALLING_FUNCTION__,"localSpannwiseVars");
// 	ZFSFloatScratchSpace globalSpannwiseVarsFLUC(noAveragedVars,totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseVars");

// 	ZFSFloatScratchSpace globalSpannwiseVarsSumIJFLUC(noAveragedVars,totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseVarsSumIJ");
// 	ZFSFloatScratchSpace globalSpannwiseVarsSumKFLUC(noAveragedVars,totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseVarsSumK");
// 	ZFSFloatScratchSpace globalSpannwiseAveragedVarsFLUC(noAveragedVars,totalNoCellsIJ,__CALLING_FUNCTION__,"globalSpannwiseAveragedVars");

	
     
// 	// cout<<"PV->noVariables*totalNoCellsIJ:"<<PV->noVariables*totalNoCellsIJ<<endl;
// 	// cout<<"totalNoCellsIJ:"<<totalNoCellsIJ<<endl;
// 	// cout<<"m_nCells[0]:"<<m_nCells[0]<<"m_nCells[1]:"<<m_nCells[1]<<" m_nCells[2]:"<<m_nCells[2]<<" m_nActiveCells[0]:"<<m_nActiveCells[0]<<" m_nActiveCells[1]:"<<m_nActiveCells[1]<<" m_nActiveCells[2]:"<<m_nActiveCells[2]<<" m_nOffsetCells[0]:"<<m_nOffsetCells[0]<<" m_nOffsetCells[1]:"<<m_nOffsetCells[1]<<" m_nOffsetCells[2]"<<m_nOffsetCells[2]<<"domainId:"<<domainId()<<endl;
// 	// for(ZFSId k=m_nOffsetCells[0]; k<m_nOffsetCells[0]+m_nActiveCells[0]; k++){
// 	//   for(ZFSId j=m_nOffsetCells[0]; j<m_nOffsetCells[1]+m_nActiveCells[1]; j++){
// 	//     for(ZFSId i=m_nOffsetCells[0]; i<m_nOffsetCells[2]+m_nActiveCells[2]; i++){
// 	//       ZFSIad cellId3D= i +(j+k*m_nActiveCells[1])*m_nActiveCells[2];
// 	//       ZFSId cellId2D= i +(j*m_nActiveCells[2]);

// 	for(ZFSId k=0; k<m_nActiveCells[0]; k++){
// 	  for(ZFSId j=0; j<m_nActiveCells[1]; j++){
// 	    for(ZFSId i=0; i<m_nActiveCells[2]; i++){
// 	      ZFSId localCellId3D= i+m_noGhostLayers+((j+m_noGhostLayers)+(k+m_noGhostLayers)*m_nCells[1])*m_nCells[2];
// 	      ZFSId globalCellId2D =(i+m_nOffsetCells[2])+(j+m_nOffsetCells[1])*(m_totalGridBlockDim[LESBlock][2]-1);
	      
// 		localSpannwiseVarsFLUC(0,globalCellId2D) +=m_cells->fq[FQ->FLUC_UU][localCellId3D]; // adding the vars in spanwise direction
// 		localSpannwiseVarsFLUC(1,globalCellId2D) +=m_cells->fq[FQ->FLUC_VV][localCellId3D]; // adding the vars in spanwise direction
// 	    	localSpannwiseVarsFLUC(2,globalCellId2D) +=m_cells->fq[FQ->FLUC_WW][localCellId3D]; // adding the vars in spanwise direction
// 	    	localSpannwiseVarsFLUC(3,globalCellId2D) +=m_cells->fq[FQ->FLUC_UV][localCellId3D]; // adding the vars in spanwise direction
// 	    	localSpannwiseVarsFLUC(4,globalCellId2D) +=m_cells->fq[FQ->FLUC_UW][localCellId3D]; // adding the vars in spanwise direction
// 	    	localSpannwiseVarsFLUC(5,globalCellId2D) +=m_cells->fq[FQ->FLUC_VW][localCellId3D]; // adding the vars in spanwise direction
	    
	      
// 	    }
// 	  }
// 	}

// 	MPI_Allreduce(&localSpannwiseVarsFLUC(0,0),&globalSpannwiseVarsFLUC(0,0),noAveragedVars*totalNoCellsIJ,MPI_DOUBLE,MPI_SUM,m_commZonal[LESBlock]);
// 	for(ZFSId i=0; i<totalNoCellsIJ; i++){
// 	  for(ZFSId vars=0; vars<noAveragedVars; vars++){
// 	    globalSpannwiseVarsFLUC(vars,i)/=(m_totalGridBlockDim[LESBlock][0]-1);
// 	  }
// 	}

// 	for (ZFSId k=0; k<m_nActiveCells[0]; k++)
// 	  {
// 	    for (ZFSId j=0; j<m_nActiveCells[1]; j++)
// 	      {
// 		for (ZFSId i=0; i<m_nActiveCells[2]; i++)
// 		  {
// 		    ZFSId localCellId3D= i+m_noGhostLayers+((j+m_noGhostLayers)+(k+m_noGhostLayers)*m_nCells[1])*m_nCells[2];
// 		    ZFSId globalCellId2D =(i+m_nOffsetCells[2])+(j+m_nOffsetCells[1])*(m_totalGridBlockDim[LESBlock][2]-1);

// 		    m_cells->fq[FQ->FLUC_UU][localCellId3D]   = globalSpannwiseVarsFLUC(0,globalCellId2D);
// 		    m_cells->fq[FQ->FLUC_VV][localCellId3D] = globalSpannwiseVarsFLUC(1,globalCellId2D);
// 		    m_cells->fq[FQ->FLUC_WW][localCellId3D] = globalSpannwiseVarsFLUC(2,globalCellId2D);
// 		    m_cells->fq[FQ->FLUC_UV][localCellId3D] = globalSpannwiseVarsFLUC(3,globalCellId2D);
// 		    m_cells->fq[FQ->FLUC_UW][localCellId3D] = globalSpannwiseVarsFLUC(4,globalCellId2D);
// 		    m_cells->fq[FQ->FLUC_VW][localCellId3D] = globalSpannwiseVarsFLUC(5,globalCellId2D);

// 		  }
// 	      }
// 	  }





// 	//! COMPUTIMG NUT FOR RANS
// 	ZFSFloat dudxi=0.0, dudeta=0.0, dudzeta=0.0, dvdxi=0.0, dvdeta=0.0, dvdzeta=0.0, dwdxi=0.0, dwdeta=0.0, dwdzeta=0.0;
// 	ZFSFloat dudx=0.0,dudy=0.0,dudz=0.0,dvdx=0.0,dvdy=0.0,dvdz=0.0,dwdx=0.0,dwdy=0.0,dwdz=0.0;
// 	ZFSFloat dxidx, dxidy, dxidz, detadx, detady, detadz, dzetadx, dzetady, dzetadz;
// 	ZFSFloat s11=0.0,s12=0.0,s13=0.0,s21=0.0,s22=0.0,s23=0.0,s31=0.0,s32=0.0,s33=0.0,SijSij=0.0,SijSijLim=0.0;
// 	ZFSFloat tke = 0.0,tuarg=0.0,tu=0.0,tufac=0.0,tulim=0.0;
// 	ZFSFloat omega=0.0,eps=0.0000000000001; //Clebf=1.1, blt=4.5,Cleb=0.0
// 	ZFSFloat maxNut = 0.0;
// 	ZFSId IJK =0, IPJK=0, IJPK=0, IJKP=0, IMJK=0, IJMK=0, IJKM=0;
// 	//ZFSFloat Clebf = 0.0, blt = 0.0, Cleb = 0.0;
// 	const ZFSFloat rRe = F1 / m_Re0;
// 	ZFSFloat   epss = 1e-34;
      
// 	// first compute the velocity derivatives
// 	for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++)
// 	  {
// 	    for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++)
// 	      {		      
// 		for (ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++)
// 		  {

// 		     IJK   = cellIndex(i,j,k);	  
// 		     IPJK  = cellIndex((i+1),j, k);
// 		     IJPK  = cellIndex(i,(j+1), k);
// 		     IJKP  = cellIndex(i,j, (k+1));
// 		     IMJK  = cellIndex((i-1),j, k);
// 		     IJMK  = cellIndex(i,(j-1), k);
// 		     IJKM  = cellIndex(i,j,(k-1));
	      
// 		    dudxi = F1B2*(m_cells->fq[FQ->AVG_U][IPJK]-m_cells->fq[FQ->AVG_U][IMJK]);
// 		    dudeta = F1B2*(m_cells->fq[FQ->AVG_U][IJPK]-m_cells->fq[FQ->AVG_U][IJMK]);     
// 		    dudzeta = F1B2*(m_cells->fq[FQ->AVG_U][IJKP]-m_cells->fq[FQ->AVG_U][IJKM]);     
	 
// 		    dvdxi = F1B2*(m_cells->fq[FQ->AVG_V][IPJK]-m_cells->fq[FQ->AVG_V][IMJK]);
// 		    dvdeta = F1B2*(m_cells->fq[FQ->AVG_V][IJPK]-m_cells->fq[FQ->AVG_V][IJMK]);     
// 		    dvdzeta = F1B2*(m_cells->fq[FQ->AVG_V][IJKP]-m_cells->fq[FQ->AVG_V][IJKM]);
	 
// 		    dwdxi = F1B2*(m_cells->fq[FQ->AVG_W][IPJK]-m_cells->fq[FQ->AVG_W][IMJK]);
// 		    dwdeta = F1B2*(m_cells->fq[FQ->AVG_W][IJPK]-m_cells->fq[FQ->AVG_W][IJMK]);     
// 		    dwdzeta = F1B2*(m_cells->fq[FQ->AVG_W][IJKP]-m_cells->fq[FQ->AVG_W][IJKM]);

// 		    dxidx = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][0];
// 		    dxidy = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][1];
// 		    dxidz = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][2];
	
// 		    detadx = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][3 + 0];
// 		    detady = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][3 + 1];
// 		    detadz = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][3 + 2];
	
// 		    dzetadx = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][6 + 0];
// 		    dzetady = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][6 + 1];
// 		    dzetadz = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][6 + 2];

	
// 		    dudx = dudxi*dxidx + dudeta*detadx + dudzeta*dzetadx;
// 		    dudy = dudxi*dxidy + dudeta*detady + dudzeta*dzetady;
// 		    dudz = dudxi*dxidz + dudeta*detadz + dudzeta*dzetadz;
	
// 		    dvdx = dvdxi*dxidx + dvdeta*detadx + dvdzeta*dzetadx;
// 		    dvdy = dvdxi*dxidy + dvdeta*detady + dvdzeta*dzetady;
// 		    dvdz = dvdxi*dxidz + dvdeta*detadz + dvdzeta*dzetadz;
	
// 		    dwdx = dwdxi*dxidx + dwdeta*detadx + dwdzeta*dzetadx;
// 		    dwdy = dwdxi*dxidy + dwdeta*detady + dwdzeta*dzetady;
// 		    dwdz = dwdxi*dxidz + dwdeta*detadz + dwdzeta*dzetadz;
		
// 		    // compute the strain components
// 		    s11 = 2.0*dudx;
// 		    s12 = dvdx+dudy;
// 		    s13 = dwdx+dudz;
	
// 		    s21 = dudy+dvdx;
// 		    s22 = 2.0*dvdy;
// 		    s23 = dwdy+dvdz;
	
// 		    s31 = dudz+dwdx;
// 		    s32 = dvdz+dwdy;
// 		    s33 = 2.0*dwdz;
	
// 		    // norm of the tstrain tensors
// 		    SijSij = F1B4*( s11*s11 + s12*s12 + s13*s13 +
// 				    s21*s21 + s22*s22 + s23*s23 +
// 				    s31*s31 + s32*s32 + s33*s33 );

// 		    // limiter for strain 		 
// 		    SijSijLim = (4.0*m_Ma*m_Ma)*0.0001;
// 		    SijSij = max(SijSij,SijSijLim);
	
// 		    // compute turbulent kinetic energy
// 		    tke = F1B2*(m_cells->fq[FQ->FLUC_UU][IJK]+
// 				m_cells->fq[FQ->FLUC_VV][IJK]+
// 				m_cells->fq[FQ->FLUC_WW][IJK]);
// 		    tuarg =  max((m_cells->fq[FQ->FLUC_UU][IJK]+
// 				  m_cells->fq[FQ->FLUC_VV][IJK]+
// 				  m_cells->fq[FQ->FLUC_WW][IJK]),0.0000001);

// 		    tu = sqrt(tuarg/3.0)/m_Ma;       
// 		    tufac = tanh((tu-0.01)/0.01);    
// 		    tulim=max(tu,0.01);              
// 		    tulim=tu/tulim;                  

// 		    tke = tke*tufac*tulim;
	 
// 		    // compute dissipation
// 		    omega = sqrt(2.0*SijSij)/sqrt(0.09);
// 		    // 0.09 is c_mu
// 		    omega = max(omega,eps);
// 		    //Clebf = 1.1;
// 		    //blt = 10.0;
// 		    //Cleb = F1 / (F1 + pow( (m_cells->coordinates[1][IJK] / (Clebf * blt)), 6.0));
// 		    // Finally computing the turb viscosity reconstructed

		   
// 		    m_cells->fq[FQ->RECONST_NUT][IJK] = max(tke/(omega*rRe),0.1); ///max(Cleb, 0.5);  
// 		    maxNut = max(maxNut, m_cells->fq[FQ->RECONST_NUT][IJK]);
// 		  }
// 	      }
// 	  }

// 	cout << "MAX NUT: " << maxNut << endl;

// 	//Now compute the SA Transport Variable Rho_Nyutilde from Nu_t
// 	//START_TIMER(m_block->m_treconstruction);
// 	// compute / reconstruct SA transport variable
// 	// fRHO = 0.0;
// 	// fRHO_avg = 0.0;
// 	ZFSFloat c9to3 = 9.1*9.1*9.1,chi=0.0; //,xacc=0.01
// 	ZFSFloat dummy11 = 0.0, dummy12 = 0.0, dummy21 =0.0, dummy22 = 0.0, dummy23 =0.0;
// 	ZFSFloat x1 =0.0,x2 = 500.0,dxx=0.0;
// 	ZFSId NewtonIter = 50;
// 	ZFSFloat rtnewt = 0.0,visc_lam=0.0,visc_tur=0.0,pressure1=0.0,temperature = 0.0;
// 	// ZFSFloat gammaMinusOne = m_gamma - F1;

// 	for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++)
// 	  {
// 	    for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++)
// 	      { 
// 		for (ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++)
// 		  {
	      
// 		    IJK   = cellIndex(i,j,k);
// 		    IPJK = cellIndex(i+1,j,k);
// 		    IMJK = cellIndex(i-1,j,k);
// 		    ZFSId IM2JK = cellIndex(i-2,j,k);
// 		    ZFSId IP2JK = cellIndex(i+2,j,k);


// 		    for (ZFSId n=0; n<NewtonIter; n++)
// 		      {
			    
// 			dummy11 = F1B2*(x1+x2);
// 			dummy12 =m_cells->fq[FQ->RECONST_NUT][IJK]/max(m_cells->fq[FQ->RECONST_NUT][IJK],0.0000000000001);
// 			dummy21 =0.0;
// 			dummy22 =0.0;
// 			dummy22 =0.0; 	

// 			rtnewt = dummy11;
// 			dummy21 = rtnewt;
	  
// 			// COMPUTE LAMINAR VISCOSITY 
// 			// pressure1 = gammaMinusOne*(m_cells->fq[FQ->AVG_P][IJK]-F1B2*m_cells->fq[FQ->AVG_RHO][IJK]*(
// 			// 											       POW2(m_cells->fq[FQ->AVG_U][IJK])+
// 			// 											       POW2(m_cells->fq[FQ->AVG_V][IJK])+
// 			// 											       POW2(m_cells->fq[FQ->AVG_W][IJK])	            
// 			// 											       ));
			
// 			// COMPUTE LAMINAR VISCOSITY
		       
// 			pressure1=m_cells->fq[FQ->AVG_P][IJK];
// 			temperature = m_gamma*pressure1/(m_cells->fq[FQ->AVG_RHO][IJK]);
// 			visc_lam = zfsSUTHERLANDLAW(temperature);
// 			visc_lam = visc_lam/m_cells->fq[FQ->AVG_RHO][IJK];
// 			visc_tur =  min(m_cells->fq[FQ->RECONST_NUT][IJK],1500.0);  
// 			chi = rtnewt/visc_lam;
// 			ZFSFloat chi3 = chi*chi*chi;                               
// 			// F
// 			dummy22=((rtnewt*chi3/(chi3+c9to3))-visc_tur)*dummy12;          
// 			// DF

// 			dummy23=chi3/(chi3+9.1)+                
// 			  3.0*chi3*(1.0/(chi3+c9to3)-
// 				    chi3/(POW2(chi3+c9to3)));
// 			dxx =  dummy22/ dummy23;
	      
// 			dummy11 =  dummy21-dxx;
// 			m_cells->fq[FQ->RECONST_NUTILDE][IJK] = max(dummy11*dummy12,0.1);   
// 			// m_cells->fq[FQ->NU_T][IJK] = max(dummy11*dummy12,0.1);
			
// 			// m_cells->fq[FQ->RECONST_NUTILDE][IJK] = max(m_cells->fq[FQ->AVG_RHO][IJK]*dummy11*dummy12,0.1);
// 			// m_cells->fq[FQ->NU_T][IJK] = max(m_cells->fq[FQ->AVG_RHO][IJK]*dummy11*dummy12,0.1);
// 		      }

// 		    m_cells->fq[FQ->RECONST_NUTILDE][IMJK] =  m_cells->fq[FQ->RECONST_NUTILDE][IJK];
// 		    m_cells->fq[FQ->RECONST_NUTILDE][IPJK] =  m_cells->fq[FQ->RECONST_NUTILDE][IJK];	
// 		    m_cells->fq[FQ->RECONST_NUTILDE][IM2JK] =  m_cells->fq[FQ->RECONST_NUTILDE][IJK];
// 		    m_cells->fq[FQ->RECONST_NUTILDE][IP2JK] =  m_cells->fq[FQ->RECONST_NUTILDE][IJK];		
// 		  }
// 	      }
// 	  }
//       }
//       // saveOutputPartitions<true>();
//       // saveOutputSolution<true>(false);
//       // zfsTerm(1, __CALLING_FUNCTION__, "Leaving program at end of reconstructed!");

//      }
//  }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




	


// 	  cout<<"check1 allreduce"<<endl;
// 	  MPI_Allreduce(&localSpannwiseVars[0][0],&globalSpannwiseVars[0][0],PV->noVariables*totalNoCellsIJ,MPI_DOUBLE,MPI_SUM,m_commZonal[LESBlock]);
// 	  cout<<"check2 allreduce"<<endl;
// 	  for(ZFSId j=0; j<m_totalGridBlockDim[LESBlock][1]-1; j++){
// 	    for(ZFSId i=0; i<m_totalGridBlockDim[LESBlock][2]-1; i++){
// 	      ZFSId cellIndex = i +(j*(m_totalGridBlockDim[LESBlock][2]-1));
// 	      for(ZFSId vars=0; vars<PV->noVariables; vars++){
// 		globalSpannwiseVarsSumIJ[vars][cellIndex]+=globalSpannwiseVars[vars][cellIndex];
// 	      }
// 	    }
// 	  }
// 	}
// 	cout<<"check3 allreduce"<<endl;
//       	MPI_Allreduce(&globalSpannwiseVarsSumIJ[0][0],&globalSpannwiseVarsSumK[0][0],PV->noVariables*totalNoCellsIJ,MPI_DOUBLE,MPI_SUM,m_commZonal[LESBlock]);
// 	cout<<"check4 allreduce"<<endl;
// 	for(ZFSId j=0; j<m_totalGridBlockDim[LESBlock][1]-1; j++){
// 	  for(ZFSId i=0; i<m_totalGridBlockDim[LESBlock][2]-1; i++){
// 	    ZFSId cellId2DGlobal =i +(j*(m_totalGridBlockDim[LESBlock][2]-1));
// 	    for(ZFSId vars=0; vars<PV->noVariables; vars++){
// 	      globalSpannwiseAveragedVars[vars][cellId2DGlobal]=globalSpannwiseVarsSumK[vars][cellId2DGlobal]/m_totalGridBlockDim[LESBlock][0];
// 	     //  if(domainId()==m_commZonalRootGlobal[LESBlock]){
// 	    // 	cout<<"globalSpannwiseAveragedVars["<<vars<<"]["<<cellId2DGlobal<<"]:"<<globalSpannwiseAveragedVars[vars][cellId2DGlobal]<<endl;}
// 	    // }
// 	  }
// 	}
//       }
    
//     }
// }





//       for(ZFSId avgPos = 0; avgPos < m_noSpanwiseAvgPositions; avgPos++)
// 	{
// 	  ZFSId interpolationWindowId = m_spanwiseAvgPositions[avgPos];

// 	  //Do the reconstruction only for the interpolation area
// 	  if(m_nOffsetCells[2] <= interpolationWindowId && m_nOffsetCells[2] + m_nActiveCells[2] > interpolationWindowId && m_zoneType == "LES" && m_useSpanwiseAveraging == true)
// 	    {
// 	      ZFSId localInterpolationPos = interpolationWindowId - m_nOffsetCells[2] + m_noGhostLayers;

// 	      //First do the spanwise averaging for all participating domains
// 	      ZFSId spanwiseAvgLocalDim[3] = {0,0,0};
// 	      ZFSId spanwiseAvgOffset[3] = {0,0,0};
// 	      ZFSId spanwiseAvgGlobalDim[3] = {0,0,0};


// 	      for(ZFSId dim = 0; dim < m_spaceDimensions; dim++)
// 		{
// 		  spanwiseAvgLocalDim[dim] = m_nActiveCells[m_spaceDimensions - 1 - dim];
// 		  spanwiseAvgOffset[dim] = m_nOffsetCells[m_spaceDimensions - 1 - dim];
// 		  spanwiseAvgGlobalDim[dim] = m_totalGridBlockDim[m_inputBlockId][m_spaceDimensions - 1 - dim]-1;
// 		}

// 	      ZFSInt myrank, numprocs;
// 	      MPI_Comm_size(m_spanwiseAvgComm[avgPos], &numprocs);
// 	      MPI_Comm_rank(m_spanwiseAvgComm[avgPos], &myrank);

// 	      ZFSIntScratchSpace spanwiseLocalSizesRcv(numprocs, m_spaceDimensions, __CALLING_FUNCTION__, "spanwiseLocalSize");
// 	      ZFSIntScratchSpace spanwiseOffsetsRcv(numprocs, m_spaceDimensions, __CALLING_FUNCTION__, "spanwiseOffsetsRcv");

// 	      MPI_Gather(spanwiseAvgLocalDim, m_spaceDimensions, MPI_INT, spanwiseLocalSizesRcv.begin(), m_spaceDimensions, MPI_INT, m_spanwiseAvgRoot[avgPos], m_spanwiseAvgComm[avgPos]);
// 	      MPI_Gather(spanwiseAvgOffset, m_spaceDimensions, MPI_INT, spanwiseOffsetsRcv.begin(), m_spaceDimensions, MPI_INT, m_spanwiseAvgRoot[avgPos], m_spanwiseAvgComm[avgPos]);

// 	      {
// 		ZFSId numberElementsSend = spanwiseAvgLocalDim[1];
// 		ZFSId noVariablesAvg = 5;
// 		ZFSId varOffset = 8;

// 		ZFSFloatScratchSpace spanwiseLocalSndBuf(noVariablesAvg*numberElementsSend, __CALLING_FUNCTION__, "spanwiseLocalSndBuf");
	
// 		ZFSFloat lineAvg = F0;
// 		for(ZFSId var = 0; var < noVariablesAvg; var++)
// 		  {
// 		    for(ZFSId j=m_noGhostLayers; j< spanwiseAvgLocalDim[1] + m_noGhostLayers; j++)
// 		      {
// 			for(ZFSId k=m_noGhostLayers; k< spanwiseAvgLocalDim[2] + m_noGhostLayers; k++)
// 			  {
// 			    ZFSId cellId = cellIndex(localInterpolationPos,j,k);
// 			    lineAvg += m_cells->variables[var + varOffset][cellId];
// 			  }

// 			ZFSId localId = var*numberElementsSend  + (j-2);
// 			spanwiseLocalSndBuf[localId] = lineAvg;
// 			lineAvg = F0;
// 		      }
// 		  }

// 		ZFSFloatScratchSpace spanwiseGlobalVariablesAveraged(noVariablesAvg, spanwiseAvgGlobalDim[1], __CALLING_FUNCTION__, "spanwiseGlobalVariablesAveraged");

// 		if(myrank != m_spanwiseAvgRoot[avgPos])
// 		  {
// 		    MPI_Send(spanwiseLocalSndBuf.begin(), numberElementsSend*noVariablesAvg, MPI_DOUBLE, m_spanwiseAvgRoot[avgPos], 1, m_spanwiseAvgComm[avgPos]); 
// 		  }
// 		else
// 		  {
// 		    spanwiseGlobalVariablesAveraged.fill(F0);

// 		    for(ZFSId proc = 0; proc<numprocs; proc++)
// 		      {
// 			ZFSId numberElementsRcv = spanwiseLocalSizesRcv(proc,1);
// 			ZFSFloatScratchSpace spanwiseRcvBuf(noVariablesAvg*numberElementsRcv, __CALLING_FUNCTION__, "spanwiseLocalRcvBuf");

// 			if(proc == 0)
// 			  {
// 			    copy(spanwiseLocalSndBuf.begin(), spanwiseLocalSndBuf.end(), spanwiseRcvBuf.begin());
// 			  }
// 			else
// 			  {
// 			    MPI_Status    status;
// 			    MPI_Recv(spanwiseRcvBuf.begin(), noVariablesAvg*numberElementsRcv, MPI_DOUBLE, proc, 1, m_spanwiseAvgComm[avgPos], &status);
// 			  }

// 			//Now put them in the global array
// 			for(ZFSId var = 0; var < noVariablesAvg; var++)
// 			  {
// 			    for(ZFSId j=0; j<spanwiseLocalSizesRcv(proc,1); j++)
// 			      {
// 				ZFSId localId  = var*numberElementsRcv + j;
// 				ZFSId globalId = j+spanwiseOffsetsRcv(proc,1);
				
// 				//Add at the global pos of this slice and divide by spanwise dim
// 				spanwiseGlobalVariablesAveraged(var,globalId) += spanwiseRcvBuf[localId] / spanwiseAvgGlobalDim[2];
// 			      }    
// 			  }
// 		      }
// 		  }

// 		//Broadcast the spanwise average to everyone
// 		MPI_Bcast(spanwiseGlobalVariablesAveraged.begin(), noVariablesAvg*spanwiseAvgGlobalDim[1], MPI_DOUBLE, m_spanwiseAvgRoot[avgPos], m_spanwiseAvgComm[avgPos]);	  

// 		ZFSFloatScratchSpace spanwiseLocalAveraged(noVariablesAvg,  spanwiseAvgLocalDim[1] + 2*m_noGhostLayers, __CALLING_FUNCTION__, "spanwiseLocalAveraged");
	
// 		//Now put the global spanwise average back in a local array 
// 		//and fill the GCs accordingly
// 		for(ZFSId var = 0; var < noVariablesAvg; var++)
// 		  {
// 		    ZFSId lowerBound = 0;
// 		    ZFSId upperBound = 0;

// 		    if(spanwiseAvgOffset[1] == 0)
// 		      lowerBound = m_noGhostLayers;
// 		    if(spanwiseAvgOffset[1] + spanwiseAvgLocalDim[1] == spanwiseAvgGlobalDim[1])
// 		      upperBound = -m_noGhostLayers;

// 		    for(ZFSId j=lowerBound; j<spanwiseAvgLocalDim[1]+ 2*m_noGhostLayers + upperBound; j++)
// 		      {
// 			ZFSId globalId = spanwiseAvgOffset[1] - m_noGhostLayers + j;
// 			spanwiseLocalAveraged(var, j) = spanwiseGlobalVariablesAveraged(var,globalId);
// 		      }

// 		    //In case offset is 0 add the wall BC
// 		    if(spanwiseAvgOffset[1] == 0)
// 		      {
// 			if(var == 0 || var == 1 || var == 2 || var == 5)
// 			  {
// 			    spanwiseLocalAveraged(var,1) = -spanwiseLocalAveraged(var,2);
// 			    spanwiseLocalAveraged(var,0) = -spanwiseLocalAveraged(var,3);
// 			  }
// 			else
// 			  {
// 			    spanwiseLocalAveraged(var,1) = spanwiseLocalAveraged(var,2);
// 			    spanwiseLocalAveraged(var,0) = spanwiseLocalAveraged(var,3);
// 			  }
// 		      }

// 		    //In case we're at the end of the domain add top GC
// 		    if(spanwiseAvgOffset[1] + spanwiseAvgLocalDim[1] == spanwiseAvgGlobalDim[1])
// 		      {
// 			spanwiseLocalAveraged(var,spanwiseAvgLocalDim[1]+m_noGhostLayers) = spanwiseLocalAveraged(var,spanwiseAvgLocalDim[1]+m_noGhostLayers - 1);
// 			spanwiseLocalAveraged(var,spanwiseAvgLocalDim[1]+m_noGhostLayers+1) = spanwiseLocalAveraged(var,spanwiseAvgLocalDim[1]+m_noGhostLayers);
// 		      }
// 		  }
  
   
// 		// assign the spanwise averaged values to all the cells	
// 		for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++)
// 		  {
// 		    for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++)
// 		      {
// 			for (ZFSId i=localInterpolationPos-2; i<localInterpolationPos+3; i++)
// 			  {
// 			    IJK   = cellIndex(i,j,k);
// 			    m_cells->variables[FQ->SPAVG_RHO][IJK]   = spanwiseLocalAveraged(0,j);
// 			    m_cells->variables[CV->AVG_RHO_U][IJK] = spanwiseLocalAveraged(1,j);
// 			    m_cells->variables[CV->AVG_RHO_V][IJK] = spanwiseLocalAveraged(2,j);
// 			    m_cells->variables[CV->AVG_RHO_W][IJK] = spanwiseLocalAveraged(3,j);
// 			    m_cells->variables[CV->AVG_RHO_E][IJK] = spanwiseLocalAveraged(4,j);
// 			  }
// 		      }
// 		  }
// 	      }
	
// 	      //! COMPUTE THE FLUCTUATIONS AND STORE THEM
// 	      ZFSFloat fRHO = 0.0;
// 	      ZFSFloat fRHO_avg = 0.0;

// 	      for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++)
// 		{
// 		  for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++)
// 		    {
// 		      IJK   = cellIndex(localInterpolationPos,j,k);
// 		      fRHO = F1/m_cells->variables[CV->RHO][IJK];
// 		      fRHO_avg = F1/m_cells->variables[FQ->SPAVG_RHO][IJK];

// 		      // COMPUTING the fluctuating part of each velocity components 
// 		      m_cells->variables[CV->FLUC_U][IJK]=m_cells->variables[CV->RHO_U][IJK]*fRHO - m_cells->variables[CV->AVG_RHO_U][IJK]*fRHO_avg;
// 		      m_cells->variables[CV->FLUC_V][IJK]=m_cells->variables[CV->RHO_V][IJK]*fRHO - m_cells->variables[CV->AVG_RHO_V][IJK]*fRHO_avg;
// 		      m_cells->variables[CV->FLUC_W][IJK]=m_cells->variables[CV->RHO_W][IJK]*fRHO - m_cells->variables[CV->AVG_RHO_W][IJK]*fRHO_avg;

// 		      // computing the momemts with time averaging 
// 		      m_cells->variables[CV->FLUC_UU][IJK]=timeFacUU*(m_cells->variables[CV->FLUC_U][IJK]*m_cells->variables[CV->FLUC_U][IJK])+(F1-timeFacUU)*m_cells->variables[CV->FLUC_UU][IJK];
// 		      m_cells->variables[CV->FLUC_VV][IJK]=timeFacVV*(m_cells->variables[CV->FLUC_V][IJK]*m_cells->variables[CV->FLUC_V][IJK])+(F1-timeFacVV)*m_cells->variables[CV->FLUC_VV][IJK];
// 		      m_cells->variables[CV->FLUC_WW][IJK]=timeFacWW*(m_cells->variables[CV->FLUC_W][IJK]*m_cells->variables[CV->FLUC_W][IJK])+(F1-timeFacWW)*m_cells->variables[CV->FLUC_WW][IJK];
// 		      m_cells->variables[CV->FLUC_UV][IJK]=timeFacUV*(m_cells->variables[CV->FLUC_U][IJK]*m_cells->variables[CV->FLUC_V][IJK])+(F1-timeFacUV)*m_cells->variables[CV->FLUC_UV][IJK];
// 		      m_cells->variables[CV->FLUC_UW][IJK]=timeFacUW*(m_cells->variables[CV->FLUC_U][IJK]*m_cells->variables[CV->FLUC_W][IJK])+(F1-timeFacUW)*m_cells->variables[CV->FLUC_UW][IJK];
// 		      m_cells->variables[CV->FLUC_VW][IJK]=timeFacVW*(m_cells->variables[CV->FLUC_V][IJK]*m_cells->variables[CV->FLUC_W][IJK])+(F1-timeFacVW)*m_cells->variables[CV->FLUC_VW][IJK];	
// 		    }
// 		}

// 	      ZFSId numberElementsSend = spanwiseAvgLocalDim[1]*spanwiseAvgLocalDim[2];
// 	      ZFSId noVariablesAvg = 9;
// 	      ZFSId varOffset = 15;

// 	      ZFSFloatScratchSpace spanwiseLocalSndBuf(noVariablesAvg*numberElementsSend, __CALLING_FUNCTION__, "spanwiseLocalSndBuf");

// 	      ZFSFloat lineAvg = F0;
// 	      for(ZFSId var = 0; var < noVariablesAvg; var++)
// 		{
// 		  for(ZFSId j=m_noGhostLayers; j< spanwiseAvgLocalDim[1] + m_noGhostLayers; j++)
// 		    {
// 		      for(ZFSId k=m_noGhostLayers; k< spanwiseAvgLocalDim[2] + m_noGhostLayers; k++)
// 			{
// 			  ZFSId cellId = cellIndex(localInterpolationPos,j,k);
// 			  lineAvg += m_cells->variables[var + varOffset][cellId];
// 			}

// 		      ZFSId localId = var*numberElementsSend  + (j-2);
// 		      spanwiseLocalSndBuf[localId] = lineAvg;
// 		      lineAvg = F0;
// 		    }
// 		}

// 	      ZFSFloatScratchSpace spanwiseGlobalVariablesAveraged(noVariablesAvg, spanwiseAvgGlobalDim[1], __CALLING_FUNCTION__, "spanwiseGlobalVariablesAveraged");

// 	      if(myrank != m_spanwiseAvgRoot[avgPos])
// 		{
// 		  MPI_Send(spanwiseLocalSndBuf.begin(), numberElementsSend*noVariablesAvg, MPI_DOUBLE, m_spanwiseAvgRoot[avgPos], 1, m_spanwiseAvgComm[avgPos]); 
// 		}
// 	      else
// 		{
// 		  spanwiseGlobalVariablesAveraged.fill(F0);

// 		  for(ZFSId proc = 0; proc<numprocs; proc++)
// 		    {
// 		      ZFSId numberElementsRcv = spanwiseLocalSizesRcv(proc,1) * spanwiseLocalSizesRcv(proc,2);
// 		      ZFSFloatScratchSpace spanwiseRcvBuf(noVariablesAvg*numberElementsRcv, __CALLING_FUNCTION__, "spanwiseLocalRcvBuf");

// 		      if(proc == 0)
// 			{
// 			  copy(spanwiseLocalSndBuf.begin(), spanwiseLocalSndBuf.end(), spanwiseRcvBuf.begin());
// 			}
// 		      else
// 			{
// 			  MPI_Status    status;
// 			  MPI_Recv(spanwiseRcvBuf.begin(), noVariablesAvg*numberElementsRcv, MPI_DOUBLE, proc, 1, m_spanwiseAvgComm[avgPos], &status);
// 			}

// 		      //Now put them in the global array

// 		      for(ZFSId var = 0; var < noVariablesAvg; var++)
// 			{
// 			  for(ZFSId j=0; j<spanwiseLocalSizesRcv(proc,1); j++)
// 			    {
// 			      ZFSId localId  = var*numberElementsRcv + j;
// 			      ZFSId globalId = j+spanwiseOffsetsRcv(proc,1);
				
// 			      spanwiseGlobalVariablesAveraged(var,globalId) += spanwiseRcvBuf[localId]  / spanwiseAvgGlobalDim[2];
// 			    }    
// 			}
// 		    }
// 		}

// 	      //Broadcast the spanwise average to everyone
// 	      MPI_Bcast(spanwiseGlobalVariablesAveraged.begin(), noVariablesAvg*spanwiseAvgGlobalDim[1], MPI_DOUBLE, m_spanwiseAvgRoot[avgPos], m_spanwiseAvgComm[avgPos]);	  

// 	      ZFSFloatScratchSpace spanwiseLocalAveraged(noVariablesAvg,  spanwiseAvgLocalDim[1] + 2*m_noGhostLayers, __CALLING_FUNCTION__, "spanwiseLocalAveraged");
	    
// 	      for(ZFSId var = 0; var < noVariablesAvg; var++)
// 		{
// 		  ZFSId lowerBound = 0;
// 		  ZFSId upperBound = 0;

// 		  if(spanwiseAvgOffset[1] == 0)
// 		    lowerBound = m_noGhostLayers;
// 		  if(spanwiseAvgOffset[1] + spanwiseAvgLocalDim[1] == spanwiseAvgGlobalDim[1])
// 		    upperBound = -m_noGhostLayers;

// 		  for(ZFSId j=lowerBound; j<spanwiseAvgLocalDim[1]+ 2*m_noGhostLayers + upperBound; j++)
// 		    {
// 		      ZFSId globalId = spanwiseAvgOffset[1] - m_noGhostLayers + j;
// 		      spanwiseLocalAveraged(var, j) = spanwiseGlobalVariablesAveraged(var,globalId);
// 		    }

// 		  //In case offset is 0 add the wall BC
// 		  if(spanwiseAvgOffset[1] == 0)
// 		    {
// 		      spanwiseLocalAveraged(var,1) = -spanwiseLocalAveraged(var,2);
// 		      spanwiseLocalAveraged(var,0) = -spanwiseLocalAveraged(var,3);
// 		    }

// 		  //In case we're at the end of the domain add top GC
// 		  if(spanwiseAvgOffset[1] + spanwiseAvgLocalDim[1] == spanwiseAvgGlobalDim[1])
// 		    {
// 		      spanwiseLocalAveraged(var,spanwiseAvgLocalDim[1]+m_noGhostLayers) = spanwiseLocalAveraged(var,spanwiseAvgLocalDim[1]+m_noGhostLayers - 1);
// 		      spanwiseLocalAveraged(var,spanwiseAvgLocalDim[1]+m_noGhostLayers+1) = spanwiseLocalAveraged(var,spanwiseAvgLocalDim[1]+m_noGhostLayers);
// 		    }
// 		}

   
// 	      // assign the spanwise averaged values to all the cells	
// 	      for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++)
// 		{
// 		  for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++)
// 		    {  
// 		      for (ZFSId i=localInterpolationPos-2; i<localInterpolationPos+3; i++)
// 			{
// 			  IJK   = cellIndex(i,j,k);
// 			  m_cells->variables[CV->FLUC_UU][IJK] = spanwiseLocalAveraged(0,j);
// 			  m_cells->variables[CV->FLUC_VV][IJK] = spanwiseLocalAveraged(1,j);
// 			  m_cells->variables[CV->FLUC_WW][IJK] = spanwiseLocalAveraged(2,j);
// 			  m_cells->variables[CV->FLUC_UV][IJK] = spanwiseLocalAveraged(3,j);
// 			  m_cells->variables[CV->FLUC_UW][IJK] = spanwiseLocalAveraged(4,j);
// 			  m_cells->variables[CV->FLUC_VW][IJK] = spanwiseLocalAveraged(5,j);	      
// 			}
// 		    }   
// 		}

// 	      //! COMPUTIMG NUT FOR RANS
// 	      ZFSFloat dudxi=0.0, dudeta=0.0, dudzeta=0.0, dvdxi=0.0, dvdeta=0.0, dvdzeta=0.0, dwdxi=0.0, dwdeta=0.0, dwdzeta=0.0;
// 	      ZFSFloat dudx=0.0,dudy=0.0,dudz=0.0,dvdx=0.0,dvdy=0.0,dvdz=0.0,dwdx=0.0,dwdy=0.0,dwdz=0.0;
// 	      ZFSFloat dxidx, dxidy, dxidz, detadx, detady, detadz, dzetadx, dzetady, dzetadz;
// 	      ZFSFloat s11=0.0,s12=0.0,s13=0.0,s21=0.0,s22=0.0,s23=0.0,s31=0.0,s32=0.0,s33=0.0,SijSij=0.0,SijSijLim=0.0;
// 	      ZFSFloat tke = 0.0,tuarg=0.0,tu=0.0,tufac=0.0,tulim=0.0;
// 	      ZFSFloat omega=0.0,eps=0.0000000000001; //Clebf=1.1, blt=4.5,Cleb=0.0
// 	      ZFSFloat maxNut = 0.0;
// 	      //ZFSFloat Clebf = 0.0, blt = 0.0, Cleb = 0.0;
// 	      const ZFSFloat rRe = F1 / m_Re0;
// 	      ZFSFloat   epss = 1e-34;
      
// 	      // first compute the velocity derivatives
// 	      for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++)
// 		{
// 		  for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++)
// 		    {
// 		      ZFSId i =  localInterpolationPos;

// 		      IJK   = cellIndex(i,j,k);	  
// 		      IPJK  = cellIndex((i+1),j, k);
// 		      IJPK  = cellIndex(i,(j+1), k);
// 		      IJKP  = cellIndex(i,j, (k+1));
// 		      IMJK  = cellIndex((i-1),j, k);
// 		      IJMK  = cellIndex(i,(j-1), k);
// 		      IJKM  = cellIndex(i,j,(k-1));
	      
// 		      dudxi = F1B2*((m_cells->variables[CV->AVG_RHO_U][IPJK]/m_cells->variables[CV->AVG_RHO][IPJK])-(m_cells->variables[CV->AVG_RHO_U][IMJK]/m_cells->variables[CV->AVG_RHO][IMJK]));
// 		      dudeta = F1B2*((m_cells->variables[CV->AVG_RHO_U][IJPK]/m_cells->variables[CV->AVG_RHO][IJPK])-(m_cells->variables[CV->AVG_RHO_U][IJMK]/m_cells->variables[CV->AVG_RHO][IJMK]));     
// 		      dudzeta = F1B2*((m_cells->variables[CV->AVG_RHO_U][IJKP]/m_cells->variables[CV->AVG_RHO][IJKP])-(m_cells->variables[CV->AVG_RHO_U][IJKM]/m_cells->variables[CV->AVG_RHO][IJKM]));     
	 
// 		      dvdxi = F1B2*((m_cells->variables[CV->AVG_RHO_V][IPJK]/m_cells->variables[CV->AVG_RHO][IPJK])-(m_cells->variables[CV->AVG_RHO_V][IMJK]/m_cells->variables[CV->AVG_RHO][IMJK]));
// 		      dvdeta = F1B2*((m_cells->variables[CV->AVG_RHO_V][IJPK]/m_cells->variables[CV->AVG_RHO][IJPK])-(m_cells->variables[CV->AVG_RHO_V][IJMK]/m_cells->variables[CV->AVG_RHO][IJMK]));     
// 		      dvdzeta = F1B2*((m_cells->variables[CV->AVG_RHO_V][IJKP]/m_cells->variables[CV->AVG_RHO][IJKP])-(m_cells->variables[CV->AVG_RHO_V][IJKM]/m_cells->variables[CV->AVG_RHO][IJKM]));
	 
// 		      dwdxi = F1B2*((m_cells->variables[CV->AVG_RHO_W][IPJK]/m_cells->variables[CV->AVG_RHO][IPJK])-(m_cells->variables[CV->AVG_RHO_W][IMJK]/m_cells->variables[CV->AVG_RHO][IMJK]));
// 		      dwdeta = F1B2*((m_cells->variables[CV->AVG_RHO_W][IJPK]/m_cells->variables[CV->AVG_RHO][IJPK])-(m_cells->variables[CV->AVG_RHO_W][IJMK]/m_cells->variables[CV->AVG_RHO][IJMK]));     
// 		      dwdzeta = F1B2*((m_cells->variables[CV->AVG_RHO_W][IJKP]/m_cells->variables[CV->AVG_RHO][IJKP])-(m_cells->variables[CV->AVG_RHO_W][IJKM]/m_cells->variables[CV->AVG_RHO][IJKM]));

// 		      dxidx = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][0];
// 		      dxidy = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][1];
// 		      dxidz = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][2];
	
// 		      detadx = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][3 + 0];
// 		      detady = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][3 + 1];
// 		      detadz = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][3 + 2];
	
// 		      dzetadx = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][6 + 0];
// 		      dzetady = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][6 + 1];
// 		      dzetadz = (1./max(m_cells->cellJac[IJK], epss))*m_cells->cellMetrics[IJK][6 + 2];

	
// 		      dudx = dudxi*dxidx + dudeta*detadx + dudzeta*dzetadx;
// 		      dudy = dudxi*dxidy + dudeta*detady + dudzeta*dzetady;
// 		      dudz = dudxi*dxidz + dudeta*detadz + dudzeta*dzetadz;
	
// 		      dvdx = dvdxi*dxidx + dvdeta*detadx + dvdzeta*dzetadx;
// 		      dvdy = dvdxi*dxidy + dvdeta*detady + dvdzeta*dzetady;
// 		      dvdz = dvdxi*dxidz + dvdeta*detadz + dvdzeta*dzetadz;
	
// 		      dwdx = dwdxi*dxidx + dwdeta*detadx + dwdzeta*dzetadx;
// 		      dwdy = dwdxi*dxidy + dwdeta*detady + dwdzeta*dzetady;
// 		      dwdz = dwdxi*dxidz + dwdeta*detadz + dwdzeta*dzetadz;
		
// 		      // compute the strain components
// 		      s11 = 2.0*dudx;
// 		      s12 = dvdx+dudy;
// 		      s13 = dwdx+dudz;
	
// 		      s21 = dudy+dvdx;
// 		      s22 = 2.0*dvdy;
// 		      s23 = dwdy+dvdz;
	
// 		      s31 = dudz+dwdx;
// 		      s32 = dvdz+dwdy;
// 		      s33 = 2.0*dwdz;
	
// 		      // norm of the tstrain tensors
// 		      SijSij = F1B4*( s11*s11 + s12*s12 + s13*s13 +
// 				      s21*s21 + s22*s22 + s23*s23 +
// 				      s31*s31 + s32*s32 + s33*s33 );

// 		      // limiter for strain 		 
// 		      SijSijLim = (4.0*m_Ma*m_Ma)*0.0001;
// 		      SijSij = max(SijSij,SijSijLim);
	
// 		      // compute turbulent kinetic energy
// 		      tke = F1B2*(m_cells->variables[CV->FLUC_UU][IJK]+
// 				  m_cells->variables[CV->FLUC_VV][IJK]+
// 				  m_cells->variables[CV->FLUC_WW][IJK]);
// 		      tuarg =  max((m_cells->variables[CV->FLUC_UU][IJK]+
// 				    m_cells->variables[CV->FLUC_VV][IJK]+
// 				    m_cells->variables[CV->FLUC_WW][IJK]),0.0000001);

// 		      tu = sqrt(tuarg/3.0)/m_Ma;
// 		      tufac = tanh((tu-0.01)/0.01);
// 		      tulim=max(tu,0.01);
// 		      tulim=tu/tulim;

// 		      tke = tke*tufac*tulim;
	 
// 		      // compute dissipation
// 		      omega = sqrt(2.0*SijSij)/sqrt(0.09);
// 		      // 0.09 is c_mu
// 		      omega = max(omega,eps);
// 		      //Clebf = 1.1;
// 		      //blt = 10.0;
// 		      //Cleb = F1 / (F1 + pow( (m_cells->coordinates[1][IJK] / (Clebf * blt)), 6.0));
// 		      // Finally computing the turb viscosity reconstructed
// 		      m_cells->variables[CV->RECONST_NYUT][IJK] = max(tke/(omega*rRe),0.1); ///max(Cleb, 0.5);
// 		      maxNut = max(maxNut, m_cells->variables[CV->RECONST_NYUT][IJK]);
// 		    }
// 		}

// 	      cout << "MAX NUT: " << maxNut << endl;

// 	      //Now compute the SA Transport Variable Rho_Nyutilde from Nu_t
// 	      //START_TIMER(m_block->m_treconstruction);
// 	      // compute / reconstruct SA transport variable
// 	      fRHO = 0.0;
// 	      fRHO_avg = 0.0;
// 	      ZFSFloat c9to3 = 9.1*9.1*9.1,chi=0.0; //,xacc=0.01
// 	      ZFSFloat dummy11 = 0.0, dummy12 = 0.0, dummy21 =0.0, dummy22 = 0.0, dummy23 =0.0;
// 	      ZFSFloat x1 =0.0,x2 = 500.0,dxx=0.0;
// 	      ZFSId NewtonIter = 50;
// 	      ZFSFloat rtnewt = 0.0,visc_lam=0.0,visc_tur=0.0,pressure1=0.0,temperature = 0.0;
// 	      ZFSFloat gammaMinusOne = m_gamma - F1;

// 	      for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++)
// 		{
// 		  for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++)
// 		    {
// 		      ZFSId i =  localInterpolationPos;
	      
// 		      IJK   = cellIndex(i,j,k);
// 		      IPJK = cellIndex(i+1,j,k);
// 		      IMJK = cellIndex(i-1,j,k);
// 		      ZFSId IM2JK = cellIndex(i-2,j,k);
// 		      ZFSId IP2JK = cellIndex(i+2,j,k);


// 		      for (ZFSId n=0; n<NewtonIter; n++)
// 			{
			    
// 			  dummy11 = F1B2*(x1+x2);
// 			  dummy12 =m_cells->variables[CV->RECONST_NYUT][IJK]/max(m_cells->variables[CV->RECONST_NYUT][IJK],0.0000000000001);
// 			  dummy21 =0.0;
// 			  dummy22 =0.0;
// 			  dummy22 =0.0; 	

// 			  rtnewt = dummy11;
// 			  dummy21 = rtnewt;
	  
// 			  // COMPUTE LAMINAR VISCOSITY 
// 			  pressure1 = gammaMinusOne*(m_cells->variables[CV->AVG_RHO_E][IJK]-F1B2*m_cells->variables[CV->AVG_RHO][IJK]*(
// 																       POW2(m_cells->variables[CV->AVG_RHO_U][IJK])+
// 																       POW2(m_cells->variables[CV->AVG_RHO_V][IJK])+
// 																       POW2(m_cells->variables[CV->AVG_RHO_W][IJK])	            
// 																       ));

// 			  temperature = m_gamma*pressure1/m_cells->variables[CV->AVG_RHO][IJK];
// 			  visc_lam = zfsSUTHERLANDLAW(temperature);
// 			  visc_lam = visc_lam/m_cells->variables[CV->AVG_RHO][IJK];
// 			  visc_tur =  min(m_cells->variables[CV->RECONST_NYUT][IJK],1500.0);
// 			  chi = rtnewt/visc_lam;
// 			  ZFSFloat chi3 = pow(3.0, chi);
// 			  // F
// 			  dummy22=((rtnewt*chi3/(chi3+c9to3))-visc_tur)*dummy12;
// 			  // DF

// 			  dummy23=chi3/(chi3+9.1)+
// 			    3.0*chi3*(1.0/(chi3+c9to3)-
// 				      chi3/(POW2(chi3+c9to3)));
// 			  dxx =  dummy22/ dummy23;
	      
// 			  dummy11 =  dummy21-dxx;
// 			  m_cells->variables[CV->RECONST_NYUTILDE][IJK] = max(m_cells->variables[CV->AVG_RHO][IJK]*dummy11*dummy12,0.1);
// 			  m_cells->variables[CV->RHO_NYUTILDE][IJK] = max(m_cells->variables[CV->AVG_RHO][IJK]*dummy11*dummy12,0.1);
// 			}

// 		      m_cells->variables[CV->RHO_NYUTILDE][IMJK] =  m_cells->variables[CV->RHO_NYUTILDE][IJK];
// 		      m_cells->variables[CV->RHO_NYUTILDE][IPJK] =  m_cells->variables[CV->RHO_NYUTILDE][IJK];	
// 		      m_cells->variables[CV->RHO_NYUTILDE][IM2JK] =  m_cells->variables[CV->RHO_NYUTILDE][IJK];
// 		      m_cells->variables[CV->RHO_NYUTILDE][IP2JK] =  m_cells->variables[CV->RHO_NYUTILDE][IJK];		
// 		    }
// 		}
// 	    }
// 	}
//     }
// }



// void ZFSStrctrdBlck3D::computeCumulativeAverage(ZFSBool forceReset)
// {
//   TRACE();
//   if(m_RKStep==0 && m_zoneType == "LES") {
//     if(globalTimeStep == m_zonalStartAvgTime || forceReset) {
//       cout << "///////// Starting averaging, resetting AVG variables //////////" << endl;
//       for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
// 	m_cells->fq[FQ->AVG_RHO][cellId] = m_cells->pvariables[PV->RHO][cellId];
// 	m_cells->fq[FQ->AVG_U][cellId] = m_cells->pvariables[PV->U][cellId];
// 	m_cells->fq[FQ->AVG_V][cellId] = m_cells->pvariables[PV->V][cellId];
// 	m_cells->fq[FQ->AVG_W][cellId] = m_cells->pvariables[PV->W][cellId];
// 	m_cells->fq[FQ->AVG_P][cellId] = m_cells->pvariables[PV->P][cellId];
//       }
//     } else {
//       if(m_zonalExponentialAveraging) {
// 	const ZFSFloat timeFacRho = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
// 	const ZFSFloat timeFacU = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
// 	const ZFSFloat timeFacV = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
// 	const ZFSFloat timeFacW = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;
// 	const ZFSFloat timeFacE = m_Ma *(1.0/m_zonalAveragingFactor)* sqrt(PV->TInfinity)*m_timeStep;

// 	//! do the time average of the flow variables for LES and store them.
// 	for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
// 	  //Exponential averaging
// 	  m_cells->fq[FQ->AVG_RHO][cellId] = timeFacRho* m_cells->pvariables[PV->RHO][cellId] + (1.0 - timeFacRho)*m_cells->fq[FQ->AVG_RHO][cellId];
// 	  m_cells->fq[FQ->AVG_U][cellId] = timeFacU* m_cells->pvariables[PV->U][cellId] + (1.0 - timeFacU)*m_cells->fq[FQ->AVG_U][cellId];
// 	  m_cells->fq[FQ->AVG_V][cellId] = timeFacV* m_cells->pvariables[PV->V][cellId] + (1.0 - timeFacV)*m_cells->fq[FQ->AVG_V][cellId];
// 	  m_cells->fq[FQ->AVG_W][cellId] = timeFacW* m_cells->pvariables[PV->W][cellId] + (1.0 - timeFacW)*m_cells->fq[FQ->AVG_W][cellId];
// 	  m_cells->fq[FQ->AVG_P][cellId] = timeFacE* m_cells->pvariables[PV->P][cellId] + (1.0 - timeFacE)*m_cells->fq[FQ->AVG_P][cellId];
// 	}
//       } else {
// 	ZFSInt averagingStep = 0;

// 	if(globalTimeStep > m_startAvgTime) {
// 	  averagingStep = globalTimeStep - m_zonalStartAvgTime;
// 	} else {
// 	  averagingStep = globalTimeStep;
// 	}

// 	for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
// 	  //Normal cumulative averaging
// 	  m_cells->fq[FQ->AVG_RHO][cellId] = (m_cells->pvariables[PV->RHO][cellId]   + m_cells->fq[FQ->AVG_RHO][cellId]   * averagingStep) / (averagingStep + 1) ;
// 	  m_cells->fq[FQ->AVG_U][cellId]   = (m_cells->pvariables[PV->U][cellId] + m_cells->fq[FQ->AVG_U][cellId] * averagingStep) / (averagingStep + 1) ;
// 	  m_cells->fq[FQ->AVG_V][cellId]   = (m_cells->pvariables[PV->V][cellId] + m_cells->fq[FQ->AVG_V][cellId] * averagingStep) / (averagingStep + 1) ;
// 	  m_cells->fq[FQ->AVG_W][cellId]   = (m_cells->pvariables[PV->W][cellId] + m_cells->fq[FQ->AVG_W][cellId] * averagingStep) / (averagingStep + 1) ;
// 	  m_cells->fq[FQ->AVG_P][cellId]   = (m_cells->pvariables[PV->P][cellId] + m_cells->fq[FQ->AVG_P][cellId] * averagingStep) / (averagingStep + 1) ;
// 	}
//       }
//     }
//   }
// }


// void ZFSStrctrdBlck3D::flucFillGhostCells() {
//   flucExtrapolateVariables();
//   m_flucCmnctnFlag = new ZFSStrctrdCommunicationHandle(nDim);
//   m_windowInfo->createCommunicationExchangeFlags(m_flucCmnctnFlag, 17);
//   m_flucCmnctnFlag->setBufferSizes();
//   flucExchange();
//   delete m_flucCmnctnFlag;
// }

// void ZFSStrctrdBlck3D::flucExtrapolateVariables() //junoh
// {
//   TRACE();
//   //i-direction
//   ZFSId cellId, cellIdAdj1, cellIdAdj2;
//   for(ZFSId k=0; k<m_nCells[0]; k++) {
//     for(ZFSId j=0; j<m_nCells[1]; j++ ) {
//       for(ZFSId i=0; i<m_noGhostLayers; i++) {
//         cellId = cellIndex(m_noGhostLayers-1-i,j,k); //pointId in Array
//         cellIdAdj1 = cellIndex(m_noGhostLayers-i,j,k);
//         cellIdAdj2 = cellIndex(m_noGhostLayers+1-i,j,k);

// 	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);	  
// 	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);
// 	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);

// 	m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

// 	m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
// 	m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	


	  
//         cellId = cellIndex(m_nCells[2]-m_noGhostLayers+i,j,k); //pointId in Array
//         cellIdAdj1 = cellIndex(m_nCells[2]-m_noGhostLayers-1+i,j,k);
//         cellIdAdj2 = cellIndex(m_nCells[2]-m_noGhostLayers-2+i,j,k);

// 	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);
// 	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);
	
// 	m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

// 	m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
// 	m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	
//       }
//     }
//   }

//   //j-direction
//   for(ZFSId k=0; k<m_nCells[0]; k++) {
//     for(ZFSId i=0; i<m_nCells[2]; i++ ) {
//       for(ZFSId j=0; j<m_noGhostLayers; j++) {
//         cellId = cellIndex(i,m_noGhostLayers-1-j,k); //pointId in Array
//         cellIdAdj1 = cellIndex(i,m_noGhostLayers-j,k);
//         cellIdAdj2 = cellIndex(i,m_noGhostLayers+1-j,k);
	
// 	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);
// 	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);

//       	m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

// 	m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
// 	m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	
//         cellId = cellIndex(i,m_nCells[1]-m_noGhostLayers+j,k); //pointId in Array
//         cellIdAdj1 = cellIndex(i,m_nCells[1]-m_noGhostLayers-1+j,k);
//         cellIdAdj2 = cellIndex(i,m_nCells[1]-m_noGhostLayers-2+j,k);
        
// 	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);	
// 	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);
	
// 	m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

// 	m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
// 	m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	
//       }
//     }
//   }

//   //k-direction
//   for(ZFSId j=0; j<m_nCells[1]; j++) {
//     for(ZFSId i=0; i<m_nCells[2]; i++ ) {
//       for(ZFSId k=0; k<m_noGhostLayers; k++) {
//         cellId = cellIndex(i,j,m_noGhostLayers-1-k); //pointId in Array
//         cellIdAdj1 = cellIndex(i,j,m_noGhostLayers-k);
//         cellIdAdj2 = cellIndex(i,j,m_noGhostLayers+1-k);

// 	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);
// 	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);

// 	m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

// 	m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
// 	m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	
	
//         cellId = cellIndex(i,j,m_nCells[0]-m_noGhostLayers+k); //pointId in Array
//         cellIdAdj1 = cellIndex(i,j,m_nCells[0]-m_noGhostLayers-1+k);
//         cellIdAdj2 = cellIndex(i,j,m_nCells[0]-m_noGhostLayers-2+k);
	
// 	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);
// 	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);
// 	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);
	
// 	m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
// 	m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

// 	m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
// 	m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	
//       }
//     }
//   }
// }

// void ZFSStrctrdBlck3D::flucExchange()
// {
//   if(noDomains()>1){
//     flucGather();
//     flucSend();
//     flucReceive();
//     MPI_Waitall(m_flucCmnctnFlag->noNghbrDomainsNormal,m_flucCmnctnFlag->mpi_sndRequest,m_flucCmnctnFlag->mpi_sndStatus);
//     MPI_Waitall(m_flucCmnctnFlag->noNghbrDomainsNormal,m_flucCmnctnFlag->mpi_rcvRequest, m_flucCmnctnFlag->mpi_rcvStatus);
//     flucScatter();
//   }
// }

// void ZFSStrctrdBlck3D::flucGather()
// {
//   ZFSId cellId;
//   for(ZFSInt nghbr=0; nghbr<m_flucCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
//     ZFSInt* startInfo=m_flucCmnctnFlag->startInfoSNDcells[nghbr];
//     ZFSInt* endInfo= m_flucCmnctnFlag->endInfoSNDcells[nghbr];
//     ZFSFloat* bufferSnd = m_flucCmnctnFlag->m_bufferCellsSnd[nghbr];
//     ZFSInt pos=0;
//     ZFSInt noCells = (endInfo[0]-startInfo[0])*(endInfo[1]-startInfo[1])*(endInfo[2]-startInfo[2]);
//       for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
//         for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
//           for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
//             cellId = i +(j+k*m_nCells[1])*m_nCells[2];
// 	    // bufferSnd[pos]=m_cells->fq[var][cellId];
//             bufferSnd[pos+FQ->AVG_U*noCells]=m_cells->fq[FQ->AVG_U][cellId];
// 	    bufferSnd[pos+FQ->AVG_V*noCells]=m_cells->fq[FQ->AVG_V][cellId];
// 	    bufferSnd[pos+FQ->AVG_W*noCells]=m_cells->fq[FQ->AVG_W][cellId];
// 	    bufferSnd[pos+FQ->AVG_RHO*noCells]=m_cells->fq[FQ->AVG_RHO][cellId];
// 	    bufferSnd[pos+FQ->AVG_P*noCells]=m_cells->fq[FQ->AVG_P][cellId];
// 	    bufferSnd[pos+FQ->NU_T*noCells]=m_cells->fq[FQ->NU_T][cellId];
	    
// 	    bufferSnd[pos+6*noCells]=m_cells->fq[FQ->FLUC_U][cellId];
// 	    bufferSnd[pos+7*noCells]=m_cells->fq[FQ->FLUC_V][cellId];
// 	    bufferSnd[pos+8*noCells]=m_cells->fq[FQ->FLUC_W][cellId];
// 	    bufferSnd[pos+9*noCells]=m_cells->fq[FQ->FLUC_UU][cellId];
// 	    bufferSnd[pos+10*noCells]=m_cells->fq[FQ->FLUC_VV][cellId];
// 	    bufferSnd[pos+11*noCells]=m_cells->fq[FQ->FLUC_WW][cellId];
// 	    bufferSnd[pos+12*noCells]=m_cells->fq[FQ->FLUC_UV][cellId];
// 	    bufferSnd[pos+13*noCells]=m_cells->fq[FQ->FLUC_UW][cellId];
// 	    bufferSnd[pos+14*noCells]=m_cells->fq[FQ->FLUC_VW][cellId];
// 	    bufferSnd[pos+15*noCells]=m_cells->fq[FQ->RECONST_NUT][cellId];
// 	    bufferSnd[pos+16*noCells]=m_cells->fq[FQ->RECONST_NUTILDE][cellId];
// 	    bufferSnd[pos+17*noCells]=m_cells->fq[FQ->NUTILDE][cellId];
//             pos++;
          
//         }
//       }
//     }
//   }
// }

// void ZFSStrctrdBlck3D::flucSend()
// {
//   for(ZFSId nghbr=0; nghbr<m_avCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
//     ZFSInt tag= domainId()+(m_avCmnctnFlag->m_tagHelperSND[nghbr])*noDomains();
//     ZFSInt err = MPI_Isend((void*)&m_avCmnctnFlag->m_bufferCellsSnd[nghbr][0], m_avCmnctnFlag->m_noNghbrDomainCellBufferSizeSnd[nghbr], MPI_DOUBLE, m_avCmnctnFlag->m_sndNghbrId[nghbr], tag, MPI_COMM_WORLD, &m_avCmnctnFlag->mpi_sndRequest[nghbr]);
//     if(err) cout << "rank " << domainId() << " sending throws error " << endl;
//   }
// }

// void ZFSStrctrdBlck3D::averagedReceive()
// {
//   for(ZFSId nghbr=0; nghbr<m_avCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
//     ZFSInt tag=m_avCmnctnFlag->m_rcvNghbrId[nghbr]+(m_avCmnctnFlag->m_tagHelperRCV[nghbr])*noDomains();
//     ZFSInt err = MPI_Irecv((void*)&m_avCmnctnFlag->m_bufferCellsRcv[nghbr][0],m_avCmnctnFlag->m_noNghbrDomainCellBufferSizeRcv[nghbr], MPI_DOUBLE,m_avCmnctnFlag->m_rcvNghbrId[nghbr] ,tag, MPI_COMM_WORLD, &m_avCmnctnFlag->mpi_rcvRequest[nghbr]);
//     if(err) cout << "rank " << domainId() << " sending throws error " << endl;
//   }
// }

// void ZFSStrctrdBlck3D::averagedScatter()
// {
//   ZFSId cellId;
//   //the ordering of the grid points can be different from
//   //sending instance ==> reorder it and copy it to the
//   //right place

//   //for(ZFSInt nghbr=0; nghbr< m_ppCmnctnFlag->noNghbrDomainsNormal-m_ppCmnctnFlag->noNghbrDomainsSingular; nghbr++)
//   for(ZFSInt nghbr=0; nghbr< m_avCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
//     ZFSInt k2, j2, i2, id2_U, id2_V, id2_W, id2_RHO, id2_P, id2_NU_T, id2_FLUC_U, id2_FLUC_V, id2_FLUC_W, id2_FLUC_UU, id2_FLUC_VV, id2_FLUC_WW, id2_FLUC_UV, id2_FLUC_UW, id2_FLUC_VW, id2_RECONST_NUT, id2_RECONST_NUTILDE, id2_NUTILDE; //id2 junoh
//     ZFSInt* step1 = m_avCmnctnFlag->stepInfoRCV[nghbr];
//     ZFSInt  step2[3];
//     ZFSInt* order = m_avCmnctnFlag->orderInfo[nghbr];
//     ZFSInt start1[3];
//     ZFSInt start2[3];
//     ZFSInt end2[3];
//     ZFSInt len2[3];
//     ZFSInt totalCells=1;
//     ZFSInt len1[3];

//     for(ZFSInt j=0; j<nDim; j++) {
//       len1[j]=m_avCmnctnFlag->endInfoRCVcells[nghbr][j] - m_avCmnctnFlag->startInfoRCVcells[nghbr][j];
//       if(len1[j]!=0) totalCells*=len1[j];
//       //added    check the step for RCV part !!!!!!!!important
//       step2[order[j]]=step1[j];
//       //*******************************************************
//     }

//     for(ZFSInt j=0; j<nDim; j++) {
//       start2[j]=0;
//       end2[j]=len1[j]-1;
//       len2[order[j]]=len1[j];
//       if(step2[j]<0) {
//         ZFSInt dummy=start2[j];
//         start2[j]=end2[j];
//         end2[j]=dummy;
//       }
//     }

//     ZFSInt* startInfo=m_avCmnctnFlag->startInfoRCVcells[nghbr];
//     ZFSInt* endInfo= m_avCmnctnFlag->endInfoRCVcells[nghbr];
//     ZFSFloat* bufferRcv = m_avCmnctnFlag->m_bufferCellsRcv[nghbr];
//     ZFSInt pos=0;
//     // ZFSInt noCells = (endInfo[0]-startInfo[0])*(endInfo[1]-startInfo[1])*(endInfo[2]-startInfo[2]);
//      // for(ZFSId var=0; var<5; var++){
//       k2=start2[2];
//       for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
//         j2=start2[1];
//         for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
//           i2=start2[0];
//           for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
//             start1[order[0]]=i2;
//             start1[order[1]]=j2;
//             start1[order[2]]=k2;
	    
// 	     // id2=var*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];

//             id2_U=FQ->AVG_U*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_V=FQ->AVG_V*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_W=FQ->AVG_W*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_RHO=FQ->AVG_RHO*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_P=FQ->AVG_P*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];	    
// 	    id2_NU_T=FQ->NU_T*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_FLUC_U=6*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_FLUC_V=7*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_FLUC_W=8*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_FLUC_UU=9*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_FLUC_VV=10*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_FLUC_WW=11*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_FLUC_UV=12*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_FLUC_UW=13*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_FLUC_VW=14*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_RECONST_NUT=15*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_RECONST_NUTILDE=16*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
// 	    id2_NUTILDE=17*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    




//             cellId = i +(j+k*m_nCells[1])*m_nCells[2];

// 	    // m_cells->fq[var][cellId]= bufferRcv[id2];

//             m_cells->fq[FQ->AVG_U][cellId]= bufferRcv[id2_U];
// 	    m_cells->fq[FQ->AVG_V][cellId]= bufferRcv[id2_V];
// 	    m_cells->fq[FQ->AVG_W][cellId]= bufferRcv[id2_W];
// 	    m_cells->fq[FQ->AVG_RHO][cellId]= bufferRcv[id2_RHO];
// 	    m_cells->fq[FQ->AVG_P][cellId]= bufferRcv[id2_P];
// 	    m_cells->fq[FQ->NU_T][cellId]= bufferRcv[id2_NU_T];

// 	    m_cells->fq[FQ->FLUC_U][cellId]= bufferRcv[id2_FLUC_U];
// 	    m_cells->fq[FQ->FLUC_V][cellId]= bufferRcv[id2_FLUC_V];
// 	    m_cells->fq[FQ->FLUC_W][cellId]= bufferRcv[id2_FLUC_W];
// 	    m_cells->fq[FQ->FLUC_UU][cellId]= bufferRcv[id2_FLUC_UU];
// 	    m_cells->fq[FQ->FLUC_VV][cellId]= bufferRcv[id2_FLUC_VV];
// 	    m_cells->fq[FQ->FLUC_WW][cellId]= bufferRcv[id2_FLUC_WW];
// 	    m_cells->fq[FQ->FLUC_UV][cellId]= bufferRcv[id2_FLUC_UV];
// 	    m_cells->fq[FQ->FLUC_UW][cellId]= bufferRcv[id2_FLUC_UW];
// 	    m_cells->fq[FQ->FLUC_VW][cellId]= bufferRcv[id2_FLUC_VW];
// 	    m_cells->fq[FQ->RECONST_NUT][cellId]= bufferRcv[id2_RECONST_NUT];
// 	    m_cells->fq[FQ->RECONST_NUTILDE][cellId]= bufferRcv[id2_RECONST_NUTILDE];
// 	    m_cells->fq[FQ->NUTILDE][cellId]= bufferRcv[id2_NUTILDE];

//             i2+=step2[0];
//             pos++;
//           }
//           j2+=step2[1];
//         }
//         k2+=step2[2];
//       }
//       // }
//   }
// }






void ZFSStrctrdBlck3D::computeTimeStep()
{
  TRACE();

  m_timeStep=1000.0;
  const ZFSFloat *const RESTRICT dxtx = ALIGNED_F(m_cells->dxt[0]);
  const ZFSFloat *const RESTRICT dxty = ALIGNED_F(m_cells->dxt[1]);
  const ZFSFloat *const RESTRICT dxtz = ALIGNED_F(m_cells->dxt[2]);

  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
      for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers;i++) {
        const ZFSId cellId = cellIndex(i,j,k);
        const ZFSFloat *const RESTRICT metric = ALIGNED_F(m_cells->cellMetrics[ cellId ]);
        const ZFSFloat Frho = F1 / m_cells->pvariables[ PV->RHO ][cellId];

        // compute the speed of sound
        const ZFSFloat speedOfSound = sqrt (m_gamma * m_cells->pvariables[PV->P][cellId] * Frho);

        // no need for simplified metrics, since information is already contained
        // in cell metrics
        const ZFSFloat lenXi = sqrt( POW2( metric[0] ) +
                                     POW2( metric[1] ) +
                                     POW2( metric[2] ) );

        const ZFSFloat lenEt = sqrt( POW2( metric[3] ) +
                                     POW2( metric[4] ) +
                                     POW2( metric[5] ) );

        const ZFSFloat lenZe = sqrt( POW2( metric[6] ) +
                                     POW2( metric[7] ) +
                                     POW2( metric[8] ) );

        // contravariant velocities
        ZFSFloat U_c = F0;
        ZFSFloat V_c = F0;
        ZFSFloat W_c = F0;

        for( ZFSInt isd = xsd; isd < nDim; isd ++ ) {
          U_c += m_cells->pvariables[ PV->VV[isd] ][ cellId ] * metric[ xsd * nDim + isd ];
          V_c += m_cells->pvariables[ PV->VV[isd] ][ cellId ] * metric[ ysd * nDim + isd ];
          W_c += m_cells->pvariables[ PV->VV[isd] ][ cellId ] * metric[ zsd * nDim + isd ];
        }

        // subtract grid velocity
        U_c -= dxtx[cellId];
        V_c -= dxty[cellId];
        W_c -= dxtz[cellId];

        U_c = fabs(U_c);
        V_c = fabs(V_c);
        W_c = fabs(W_c);

        // has area information in it due to metric terms
        const ZFSFloat eigenvalue = U_c + V_c + W_c + speedOfSound * ( lenXi + lenEt + lenZe );

        // divide volume information (jacobian) through area to get characteristic length for CFL
        const ZFSFloat deltaT = m_cfl * m_cells->cellJac[cellId] / eigenvalue;

        if(m_localTimeStep) {
          m_cells->localTimeStep[cellId] = deltaT;
          m_timeStep = F1;
          m_timeRef = F1;
        } else {
          m_timeStep = zfsMIN(m_timeStep, deltaT);
        }
      }
    }
  }
}


void ZFSStrctrdBlck3D::updateSpongeLayer()
{
  TRACE();
  if(m_useSponge) m_strctrdBndryCnd->updateSpongeLayer();
}


bool ZFSStrctrdBlck3D::rungeKuttaStep()
{
  TRACE();
  const ZFSId noVars = CV->noVariables;
  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSFloat rkAlpha = m_RKalpha[ m_RKStep ];
  const ZFSFloat rkFactor = rkAlpha*m_timeStep;

  ZFSFloat *const RESTRICT oldVars = ALIGNED_F(m_cells->oldVariables[0]);
  ZFSFloat *const RESTRICT vars = ALIGNED_F(m_cells->variables[0]);
  ZFSFloat *const RESTRICT oldCellJac= ALIGNED_MF(m_cells->oldCellJac);
  const ZFSFloat *const RESTRICT cellJac= ALIGNED_MF(m_cells->cellJac);
  const ZFSFloat *const RESTRICT rhs= ALIGNED_MF(m_cells->rightHandSide[0]);

//#ifdef ZFS_EXTRA_DEBUG
//  savePartitions(); //testing only
//#endif

  // set old variables
  if( m_RKStep == 0 ) {
    for(ZFSId v=0; v< noVars; v++) {
      const ZFSUint offset = v*noCells;
      ZFSFloat *const RESTRICT oldCellVars = ALIGNED_F(oldVars + offset);
      const ZFSFloat *const RESTRICT cellVars = ALIGNED_F(vars + offset);
      for(ZFSId cellId =0; cellId<m_noStrctrdCells; cellId++) {
        oldCellVars[cellId] = cellVars[cellId];
      }
    }
  }

  switch( m_rungeKuttaOrder ) {
  case 2: {
    //for moving grids we take the old Jacobian into account
    if (m_localTimeStep) {
      for(ZFSId v=0; v<noVars; v++) {
        const ZFSUint cellOffset = v*noCells;
        ZFSFloat *const RESTRICT cellVars = vars + cellOffset;
        const ZFSFloat *const RESTRICT oldCellVars = oldVars + cellOffset;
        const ZFSFloat *const RESTRICT cellRhs = rhs + cellOffset;

        for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
          for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
            for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers;i++) {
              const ZFSId cellId = cellIndex(i,j,k);
              const ZFSId localRkFactor = rkAlpha*m_cells->localTimeStep[cellId];
              const ZFSFloat factor= localRkFactor / m_cells->cellJac[cellId];
              cellVars[cellId]=oldCellVars[cellId]+factor*cellRhs[cellId];
            }
          }
        }
      }
    } else if (m_movingGrid) {
      for(ZFSId v=0; v<noVars; v++) {
        const ZFSUint cellOffset = v*noCells;
        ZFSFloat *const RESTRICT cellVars = vars + cellOffset;
        const ZFSFloat *const RESTRICT oldCellVars = oldVars + cellOffset;
        const ZFSFloat *const RESTRICT cellRhs = rhs + cellOffset;

        for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
          for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
            for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers;i++) {
              const ZFSId cellId = cellIndex(i,j,k);
              cellVars[cellId]=(oldCellVars[cellId]*oldCellJac[cellId]+rkFactor*cellRhs[cellId])/cellJac[cellId];
            }
          }
        }
      }
    } else {
      for(ZFSId v=0; v<noVars; v++) {
        const ZFSUint cellOffset = v*noCells;
        ZFSFloat *const RESTRICT cellVars = vars + cellOffset;
        const ZFSFloat *const RESTRICT oldCellVars = oldVars + cellOffset;
        const ZFSFloat *const RESTRICT cellRhs = rhs + cellOffset;

        for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
          for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
            for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers;i++) {
              const ZFSId cellId = cellIndex(i,j,k);
              const ZFSFloat factor= rkFactor / m_cells->cellJac[cellId];
              cellVars[cellId]=oldCellVars[cellId]+factor*cellRhs[cellId];
            }
          }
        }
      }
    }
    break;
  }
  case 3: {
    for(ZFSId v=0; v<noVars; v++) {
      const ZFSUint cellOffset = v*noCells;
      ZFSFloat *const RESTRICT cellVars = vars + cellOffset;
      const ZFSFloat *const RESTRICT oldCellVars = oldVars + cellOffset;
      const ZFSFloat *const RESTRICT cellRhs = rhs + cellOffset;

      for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
        for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
          for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers;i++){
            const ZFSId cellId = cellIndex(i,j,k);
            const ZFSFloat factor= rkFactor / m_cells->cellJac[cellId];
            cellVars[cellId] = rkAlpha*cellVars[cellId]+(F1-rkAlpha)*oldCellVars[cellId]-factor*cellRhs[cellId];
          }
        }
      }
    }
    break;
  }
  default:
  {
    stringstream errorMessage;
    errorMessage << "Given RungeKutta Order " << m_rungeKuttaOrder << " not implemented! " << endl;
    zfsTerm(1, __CALLING_FUNCTION__, errorMessage.str());
  }
  }

  ++m_RKStep;



  if( m_RKStep == m_noRKSteps ) {
    globalTimeStep++;
    m_physicalTime += m_timeStep * m_timeRef;
    m_time += m_timeStep;

    m_RKStep = 0;

    return true;
  } else {
    return false;
  }
}

void ZFSStrctrdBlck3D::addDisturbance()
{
  TRACE();
  cout << "enterin addDisturbance " << endl;
  cout << " m_time = " << m_time << endl;
  cout << "m_timeStep = " << m_timeStep << endl;
  cout << " global t= " << globalTimeStep << endl;
  if(m_time>=0)
  {

    ZFSFloat pi2= 8.0*atan(1);;
    ZFSFloat period =0.1;
    ZFSFloat amp1= 0.00001*sin((pi2/period)*m_time);

    //ZFSFloat T=(m_cells->variables[CV->RHO_E][364])*(m_gamma-F1)/287.1500;
    //ZFSFloat fluc_p = rhsdist*287.15000*T;
    //ZFSFloat fluc_rhoE= (F1/(m_gamma-1))*fluc_p;
    //rhsdist=fluc_rhoE;
    ZFSId cellId=0;
    cout << "sin " << amp1 << endl;
    ZFSFloat centerloc = 3.1415;
    //go through all cells and adopt disturbance smoothly!!
    //m_cells->rightHandSide[CV->RHO_E][2191941]+=rhsdist;
    for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++)
    {
      for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++)
      {
        for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers;i++)
        {
          cellId = i+(j+k*m_nCells[1])*m_nCells[2];
          ZFSFloat x = m_cells->coordinates[0][cellId];
          ZFSFloat y = m_cells->coordinates[1][cellId];
          ZFSFloat z = m_cells->coordinates[2][cellId];
          ZFSFloat r = 0.025;
          ZFSFloat a1= POW2(x-centerloc);
          ZFSFloat a2= POW2(y-centerloc);
          ZFSFloat a3= POW2(z-centerloc);
          ZFSFloat disturb = amp1*exp(-(a1+a2+a3)/POW2(r)/2.0);
          m_cells->rightHandSide[CV->RHO_E][cellId]+=disturb;
          /*    ZFSFloat x=m_cells->coordinates[0][cellId]-centreloc;
                ZFSFloat y=m_cells->coordinates[1][cellId]-centreloc;
                ZFSFloat z=m_cells->coordinates[2][cellId]-centreloc;
                m_cells->rightHandSide[CV->RHO_E][cellId]+=(F1+tanh(x/0.009))*(F1-tanh(x/0.009))*(F1+tanh(y/0.009))*(F1-tanh(y/0.009))*(F1+tanh(z/0.009))*(F1-tanh(z/0.009))*rhsdist;*/
        }
      }
    }
  }

  cout << "leaving addDisturbance " << endl;
}

void ZFSStrctrdBlck3D::assignBndryCells()
{
  TRACE();
  m_strctrdBndryCnd->assignBndryCnds();
}

void ZFSStrctrdBlck3D::initBndryCnds()
{
  TRACE();
  m_strctrdBndryCnd->correctBndryCndIndices();
}

/** \brief Extrapolates and exchanges ghost point coordinates
 * \author Pascal Meysonnat
 * \date 01.01.1010
 */
void ZFSStrctrdBlck3D::addGhostPointCoordinateValues()
{
  TRACE();
  // > for debugging only: save cellId/pointId and domainId for all cells and points
  if(m_debugOutput){
    for(ZFSId k=0; k<(m_nCells[0]); k++ ){
      for(ZFSId j=0; j<(m_nCells[1]); j++ ){
        for(ZFSId i=0; i<(m_nCells[2]); i++){
          ZFSId cellId = cellIndex(i,j,k);
          m_cells->fq[FQ->CELLID][cellId]=cellId;
          m_cells->fq[FQ->BLOCKID][cellId]=domainId();
        }
      }
    }


    //for the points also
    for(ZFSId k=0; k<(m_nPoints[0]); k++ ) {
      for(ZFSId j=0; j<(m_nPoints[1]); j++ ) {
        for(ZFSId i=0; i<(m_nPoints[2]); i++) {
          ZFSId pointId = i+(j+k*m_nPoints[1])*m_nPoints[2];
          pointProperties[0][pointId]=pointId;
          pointProperties[1][pointId]=getBoxId(domainId());
        }
      }
    }
  }
  //< end debugging only

  extrapolateGhostPointCoordinates();

  if(noDomains()>1) exchangePoints();
  if(m_periodicConnection) m_strctrdBndryCnd->exchangePointsPeriodic();

  extrapolateGhostPointCoordinatesBC();

  computeCellCentreCoordinates();

  //MUST be done after cell center computation!!!
  if(noDomains()>1&&m_cmnctnFlag->noNghbrDomainsSingular>0) exchangePointsSingularity();
  if(m_periodicConnection&&m_cmnctnFlag->noNghbrDomainsPeriodicS>0) m_strctrdBndryCnd->exchangePointsPeriodicS();

  if(m_hasSingularity>0) {
    computeReconstructionConstantsSVD();
  }

  if(m_savePartitionOutput)
  {
    //3) write the totalGridFile with GhostPoints
    writeGridPointsWithGhostPoints();
  }
}

void ZFSStrctrdBlck3D::extrapolateGhostPointCoordinatesBC()
{
  for(ZFSInt bcId=0;bcId<(ZFSInt)m_strctrdBndryCnd->m_physicalBCMap.size(); ++bcId) {
    // all the periodic BCs are NOT included.
    // also skip the channel bc
    if( m_strctrdBndryCnd->m_physicalBCMap[bcId]->BC==2401||m_strctrdBndryCnd->m_physicalBCMap[bcId]->BC==2402) {
      continue;
    }

    ZFSInt* start = m_strctrdBndryCnd->m_physicalBCMap[bcId]->start1;
    ZFSInt* end = m_strctrdBndryCnd->m_physicalBCMap[bcId]->end1;
    ZFSInt index= m_strctrdBndryCnd->m_physicalBCMap[bcId]->face/2;
    ZFSInt step=  m_strctrdBndryCnd->m_physicalBCMap[bcId]->face%2;
    ZFSInt pos[3],fix[3],mirror[3],ijk[3],extendijk[3];
    ZFSInt pointId,FixPointId,MirrorPointId;

    extendijk[0]=1;extendijk[1]=1;extendijk[2]=1;
    extendijk[index]=0;

    for(ijk[2]=start[2]; ijk[2]<end[2]+extendijk[2]; ++ijk[2]) {
      for(ijk[1]=start[1]; ijk[1]<end[1]+extendijk[1]; ++ijk[1]) {
        for(ijk[0]=start[0]; ijk[0]<end[0]+extendijk[0]; ++ijk[0]) {
          for(ZFSInt m=0;m<3;++m) {
            if(index==m) {
              if(step==1) {
                pos[m]=ijk[m]+1;
                fix[m]=start[m];
                mirror[m]=2*fix[m]-pos[m];
              } else {
                pos[m]=ijk[m];
                fix[m]=end[m];
                mirror[m]=2*fix[m]-pos[m];
              }
            } else {
              pos[m]=ijk[m];
              fix[m]=ijk[m];
              mirror[m]=ijk[m];
            }
          }//m

          pointId       =pointIndex(pos[0],pos[1],pos[2]);
          FixPointId    =pointIndex(fix[0],fix[1],fix[2]);
          MirrorPointId =pointIndex(mirror[0],mirror[1],mirror[2]);

          for(ZFSId dim =0; dim < nDim; dim++) {
            m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
          }
        }//ijk
      }
    }
  }//bcid
}

/** \brief Extrapolates ghost point coordinates from active point coordinates
 * \author Pascal Meysonnat
 */
void ZFSStrctrdBlck3D::extrapolateGhostPointCoordinates()
{
  TRACE();
  //This function mirrors the grid points on the faces
  //Can this be written in a more efficient way????

  //i-direction
  ZFSId pointId, FixPointId, MirrorPointId;

  for(ZFSId k=m_noGhostLayers; k<(m_nPoints[0]-m_noGhostLayers); k++)
  {
    for(ZFSId j=m_noGhostLayers; j<(m_nPoints[1]-m_noGhostLayers); j++ )
    {
      for(ZFSId i=0; i<m_noGhostLayers; i++)
      {
        pointId = (m_noGhostLayers -1 -i)+(j+k*m_nPoints[1])*m_nPoints[2]; //pointId in Array
        FixPointId = (m_noGhostLayers-i)+(j+k*m_nPoints[1])*m_nPoints[2]; //point about which everything is mirrored
        MirrorPointId = (m_noGhostLayers+1 - i)+(j+k*m_nPoints[1])*m_nPoints[2];
        for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }
        //coordinates at the other end!!

        pointId = (m_nPoints[2]-i-1)+(j+k*m_nPoints[1])*m_nPoints[2];
        FixPointId = (m_nPoints[2]-m_noGhostLayers-1)+(j+k*m_nPoints[1])*m_nPoints[2];
        MirrorPointId = (m_nPoints[2]-1-(2*m_noGhostLayers - i))+(j+k*m_nPoints[1])*m_nPoints[2];
        for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }

  //j-direction

  for(ZFSId k=m_noGhostLayers; k<(m_nPoints[0]-m_noGhostLayers); k++)
  {
    for(ZFSId j=0; j<m_noGhostLayers; j++ )
    {
      for(ZFSId i=m_noGhostLayers; i<(m_nPoints[2]-m_noGhostLayers); i++)
      {
        pointId = i+((m_noGhostLayers-j-1)+k*m_nPoints[1])*m_nPoints[2]; //pointId in Array
        FixPointId = i+((m_noGhostLayers-j)+k*m_nPoints[1])*m_nPoints[2]; //point about which everything is mirrored
        MirrorPointId = i+((m_noGhostLayers+1-j)+k*m_nPoints[1])*m_nPoints[2];
              for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }
        //coordinates at the other end!!
        pointId = i+((m_nPoints[1]-j-1)+k*m_nPoints[1])*m_nPoints[2];
        FixPointId =i+ ((m_nPoints[1]-m_noGhostLayers-1)+k*m_nPoints[1])*m_nPoints[2];
        MirrorPointId =i+ ((m_nPoints[1]-1-(2*m_noGhostLayers - j))+k*m_nPoints[1])*m_nPoints[2];
              for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }

  //k-direction
  for(ZFSId k=0; k<m_noGhostLayers; k++)
  {
    for(ZFSId j=0; j<(m_nPoints[1]-m_noGhostLayers); j++ )
    {
      for(ZFSId i=m_noGhostLayers; i<(m_nPoints[2]-m_noGhostLayers); i++)
      {
        pointId = i+(j+(m_noGhostLayers-1-k)*m_nPoints[1])*m_nPoints[2]; //pointId in Array
        FixPointId = i+(j+(m_noGhostLayers-k)*m_nPoints[1])*m_nPoints[2]; //point about which everything is mirrored
        MirrorPointId = i+(j+(m_noGhostLayers+1-k)*m_nPoints[1])*m_nPoints[2]; //m_noGhostLayers+(m_noGhostLayers-i)
              for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }


        //coordinates at the other end!!
        pointId = i+(j+(m_nPoints[0]-k-1)*m_nPoints[1])*m_nPoints[2];
        FixPointId =i+ (j+(m_nPoints[0]-m_noGhostLayers-1)*m_nPoints[1])*m_nPoints[2];
        MirrorPointId =i+ (j+(m_nPoints[0]-1-(2*m_noGhostLayers - k))*m_nPoints[1])*m_nPoints[2];
              for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }

  //corner points missing yet!! They only need to be calculated if a Visualisation tool is used to
  //show the ghost points and the grid


  //in i-direction

  for(ZFSId k=0; k<(m_nPoints[0]); k++)
  {
    for(ZFSId j=0; j<(m_nPoints[1]); j++ )
    {
      for(ZFSId i=0; i<m_noGhostLayers; i++)
      {

        pointId = (m_noGhostLayers -1 -i)+(j+k*m_nPoints[1])*m_nPoints[2]; //pointId in Array
        FixPointId = (m_noGhostLayers-i)+(j+k*m_nPoints[1])*m_nPoints[2]; //point about which everything is mirrored
        MirrorPointId = (m_noGhostLayers+1 - i)+(j+k*m_nPoints[1])*m_nPoints[2]; //m_noGhostLayers+(m_noGhostLayers-i)
              for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }
        //coordinates at the other end!!

        pointId = (m_nPoints[2]-i-1)+(j+k*m_nPoints[1])*m_nPoints[2];
        FixPointId = (m_nPoints[2]-m_noGhostLayers-1)+(j+k*m_nPoints[1])*m_nPoints[2];
        MirrorPointId = (m_nPoints[2]-1-(2*m_noGhostLayers - i))+(j+k*m_nPoints[1])*m_nPoints[2];
              for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }

      }
    }
  }

  for(ZFSId k=0; k<(m_nPoints[0]); k++)
  {
    for(ZFSId j=0; j<m_noGhostLayers; j++ )
    {
      for(ZFSId i=0; i<(m_nPoints[2]); i++)
      {

        pointId = i+((m_noGhostLayers-j-1)+k*m_nPoints[1])*m_nPoints[2]; //pointId in Array
        FixPointId = i+((m_noGhostLayers-j)+k*m_nPoints[1])*m_nPoints[2]; //point about which everything is mirrored
        MirrorPointId = i+((m_noGhostLayers+1-j)+k*m_nPoints[1])*m_nPoints[2];
              for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }
        //coordinates at the other end!!
        pointId = i+((m_nPoints[1]-j-1)+k*m_nPoints[1])*m_nPoints[2];
        FixPointId =i+ ((m_nPoints[1]-m_noGhostLayers-1)+k*m_nPoints[1])*m_nPoints[2];
        MirrorPointId =i+ ((m_nPoints[1]-1-(2*m_noGhostLayers - j))+k*m_nPoints[1])*m_nPoints[2];
              for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }

  for(ZFSId k=0; k<m_noGhostLayers; k++)
  {
    for(ZFSId j=0; j<(m_nPoints[1]); j++ )
    {
      for(ZFSId i=0; i<(m_nPoints[2]); i++)
      {

        pointId = i+(j+(m_noGhostLayers-1-k)*m_nPoints[1])*m_nPoints[2]; //pointId in Array
        FixPointId = i+(j+(m_noGhostLayers-k)*m_nPoints[1])*m_nPoints[2]; //point about which everything is mirrored
        MirrorPointId = i+(j+(m_noGhostLayers+1-k)*m_nPoints[1])*m_nPoints[2];
              for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }
        //coordinates at the other end!!
        pointId = i+(j+(m_nPoints[0]-k-1)*m_nPoints[1])*m_nPoints[2];
        FixPointId =i+ (j+(m_nPoints[0]-m_noGhostLayers-1)*m_nPoints[1])*m_nPoints[2];
        MirrorPointId =i+ (j+(m_nPoints[0]-1-(2*m_noGhostLayers - k))*m_nPoints[1])*m_nPoints[2];
              for(ZFSId dim =0; dim < nDim; dim++)
        {
          m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }
}

/** \brief Saves the partitioned cell coordinates with ghost cells
 * \author Pascal Meysonnat
 */
void ZFSStrctrdBlck3D::writeGridCellsWithGhostCells()
{
  TRACE();
  //first every process needs to create the datasets and the structure where the
  //data will be stored!
  const char* fileName = "totalGridCells";
  ZFSInt file = io_openfile("hdf5", fileName, "collective", m_zfsStrctrdComm);
  ZFSInt noCells[3] = {0,0,0};
  for(ZFSInt i=0; i<noDomains(); i++) {
    //create datasets for the io library
    for(ZFSId j=0; j<nDim; j++) {
      noCells[j]=m_partition->outputBoxInfo[getBoxId(i)]->DirLast[j]-1+2*m_noGhostLayers;
    }
    stringstream path;
    path << getBoxId(i);
    ZFSString solutionpath = "cpu";
    solutionpath += path.str();
    const char* dsetname = solutionpath.c_str();
    io_create_ddataset(file, dsetname , "x", 3,  noCells);
    io_create_ddataset(file, dsetname , "y", 3,  noCells);
    io_create_ddataset(file, dsetname , "z", 3,  noCells);
  }
  //write the values into the array so that we can visualize it
  ZFSInt offset[3]={0,0,0};
  for(ZFSId j=0; j<nDim; j++) {
    noCells[j]=m_partition->outputBoxInfo[getBoxId(domainId())]->DirLast[j]-1+2*m_noGhostLayers;
  }
  stringstream path;
  path << getBoxId(domainId());
  ZFSString solutionpath = "cpu";
  solutionpath += path.str();
  ZFSString x= solutionpath+ "/x";
  ZFSString y= solutionpath+ "/y";
  ZFSString z= solutionpath+ "/z";

  io_write_ddataset_part(file, solutionpath.c_str(), "x", nDim, noCells, offset, &m_cells->coordinates[0][0]);
  io_write_ddataset_part(file, solutionpath.c_str(), "y", nDim, noCells, offset, &m_cells->coordinates[1][0]);
  io_write_ddataset_part(file, solutionpath.c_str(), "z", nDim, noCells, offset, &m_cells->coordinates[2][0]);
  io_closefile(file);
}

//============================================================================================================
//====================================COMMUNICATIONS==========================================================





void ZFSStrctrdBlck3D::gatherPoints()
{
  ZFSId pointId;
  for(ZFSInt nghbr=0; nghbr< m_cmnctnFlag->noNghbrDomainsNormal-m_cmnctnFlag->noNghbrDomainsSingular; nghbr++)
  {
    ZFSInt* startInfo=m_cmnctnFlag->startInfoSNDpoints[nghbr];
    ZFSInt* endInfo= m_cmnctnFlag->endInfoSNDpoints[nghbr];
    ZFSFloat* bufferSnd = m_cmnctnFlag->m_bufferPointsSnd[nghbr];
    ZFSInt pos=0;
    for(ZFSId dim=0; dim<nDim; dim++)
    {
      for(ZFSInt k=startInfo[2]; k<endInfo[2]+1; k++)
      {
        for(ZFSInt j=startInfo[1]; j<endInfo[1]+1; j++)
        {
          for(ZFSInt i=startInfo[0]; i<endInfo[0]+1; i++)
          {
            pointId = i +(j+k*m_nPoints[1])*m_nPoints[2];
            bufferSnd[pos]=m_coordinates[dim][pointId];
            pos++;
          }
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::sendPoints()
{
  // MPI_Waitall( m_cmnctnFlag->noNghbrDomainsNormal, m_cmnctnFlag->mpi_sndRequest, m_cmnctnFlag->mpi_sndStatus);
  ZFSInt tag=-1;
  for(ZFSId nghbr=0; nghbr<m_cmnctnFlag->noNghbrDomainsNormal-m_cmnctnFlag->noNghbrDomainsSingular; nghbr++)
  {
    tag= domainId()+(m_cmnctnFlag->m_tagHelperSND[nghbr])*noDomains();
    ZFSInt err= MPI_Isend((void*)&m_cmnctnFlag->m_bufferPointsSnd[nghbr][0], m_cmnctnFlag->m_noNghbrDomainPointBufferSizeSnd[nghbr], MPI_DOUBLE, m_cmnctnFlag->m_sndNghbrId[nghbr], tag, MPI_COMM_WORLD, &m_cmnctnFlag->mpi_sndRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}


void ZFSStrctrdBlck3D::receivePoints()
{
  ZFSInt tag=-1;
  for(ZFSId nghbr=0; nghbr<m_cmnctnFlag->noNghbrDomainsNormal-m_cmnctnFlag->noNghbrDomainsSingular; nghbr++)
  {
    tag=m_cmnctnFlag->m_rcvNghbrId[nghbr]+(m_cmnctnFlag->m_tagHelperRCV[nghbr])*noDomains();
    MPI_Irecv((void*)&m_cmnctnFlag->m_bufferPointsRcv[nghbr][0],m_cmnctnFlag->m_noNghbrDomainPointBufferSizeRcv[nghbr], MPI_DOUBLE,m_cmnctnFlag->m_rcvNghbrId[nghbr] ,tag, MPI_COMM_WORLD, &m_cmnctnFlag->mpi_rcvRequest[nghbr]);
  }
}


void ZFSStrctrdBlck3D::scatterPoints()
{
  ZFSId pointId;
  //  ZFSInt cellId;
  //the ordering of the grid points can be different from
  //sending instance ==> reorder it and copy it to the
  //right place
  for(ZFSInt nghbr=0; nghbr< m_cmnctnFlag->noNghbrDomainsNormal-m_cmnctnFlag->noNghbrDomainsSingular; nghbr++)
  {
    ZFSInt k2, j2, i2, id2;

    ZFSInt* step1 = m_cmnctnFlag->stepInfoRCV[nghbr];
    ZFSInt  step2[3];

    ZFSInt* order = m_cmnctnFlag->orderInfo[nghbr];
    ZFSInt start1[3];
    //    ZFSInt end1[3];
    ZFSInt start2[3];
    ZFSInt end2[3];
    ZFSInt len2[3];
    ZFSInt totalPoints=1;
    ZFSInt len1[3];
    for(ZFSInt j=0; j<nDim; j++)
    {
      len1[j]=m_cmnctnFlag->endInfoRCVpoints[nghbr][j] - m_cmnctnFlag->startInfoRCVpoints[nghbr][j]+1;
      totalPoints*=len1[j];
      //added    check the step for RCV part !!!!!!!!important
      step2[order[j]]=step1[j];
    }
    for(ZFSInt j=0; j<nDim; j++)
    {
      start2[j]=0;
      end2[j]=len1[j]-1;
      len2[order[j]]=len1[j];
      if(step2[j]<0)
      {
        ZFSInt dummy=start2[j];
        start2[j]=end2[j];
        end2[j]=dummy;
      }
    }
    ZFSInt* startInfo=m_cmnctnFlag->startInfoRCVpoints[nghbr];
    ZFSInt* endInfo= m_cmnctnFlag->endInfoRCVpoints[nghbr];

    ZFSFloat* bufferRcv = m_cmnctnFlag->m_bufferPointsRcv[nghbr];
    for(ZFSId dim=0; dim<nDim; dim++)
    {
      k2=start2[2];
      for(ZFSInt k=startInfo[2]; k<endInfo[2]+1; k++)
      {
        j2=start2[1];
        for(ZFSInt j=startInfo[1]; j<endInfo[1]+1; j++)
        {
          i2=start2[0];
          for(ZFSInt i=startInfo[0]; i<endInfo[0]+1; i++)
          {
            start1[order[0]]=i2;
            start1[order[1]]=j2;
            start1[order[2]]=k2;
            id2=dim*totalPoints+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
            pointId = i +(j+k*m_nPoints[1])*m_nPoints[2];
            m_coordinates[dim][pointId]= bufferRcv[id2];
            i2+=step2[0];
          }
          j2+=step2[1];
        }
        k2+=step2[2];
      }
    }
  }
}

void ZFSStrctrdBlck3D::exchangePointsSingularity()
{
  gatherPointsS();

  sendPointsS();

  receivePointsS();

  for(ZFSId i=m_cmnctnFlag->noNghbrDomainsNormal- m_cmnctnFlag->noNghbrDomainsSingular; i< m_cmnctnFlag->noNghbrDomainsNormal; i++) {
    MPI_Wait(&(m_cmnctnFlag->mpi_sndRequest[i]), &(m_cmnctnFlag->mpi_sndStatus[i]));
  }

  for(ZFSId i=m_cmnctnFlag->noNghbrDomainsNormal- m_cmnctnFlag->noNghbrDomainsSingular; i< m_cmnctnFlag->noNghbrDomainsNormal; i++) {
    MPI_Wait(&(m_cmnctnFlag->mpi_rcvRequest[i]), &(m_cmnctnFlag->mpi_rcvStatus[i]));
  }

  MPI_Waitall(m_cmnctnFlag->noNghbrDomainsSingular,m_cmnctnFlag->mpi_rcvRequest, m_cmnctnFlag->mpi_rcvStatus);

  scatterPointsS();
}


void ZFSStrctrdBlck3D::gatherPointsS()
{
  //  ZFSId pointId;
  ZFSId cellId;

  for(ZFSInt nghbr= m_cmnctnFlag->noNghbrDomainsNormal- m_cmnctnFlag->noNghbrDomainsSingular; nghbr< m_cmnctnFlag->noNghbrDomainsNormal; nghbr++)
  {
    ZFSInt* startInfo=m_cmnctnFlag->startInfoSNDpoints[nghbr];
    ZFSInt* endInfo= m_cmnctnFlag->endInfoSNDpoints[nghbr];
    ZFSFloat* bufferSnd = m_cmnctnFlag->m_bufferPointsSnd[nghbr];

    ZFSInt pos=0;
    for(ZFSId dim=0; dim<nDim; dim++)
    {
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++)
      {
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++)
        {
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++)
          {
            //pointId = i +(j+k*m_nPoints[1])*m_nPoints[2];
            // bufferSnd[pos]=m_coordinates[dim][pointId];
            cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            bufferSnd[pos]= m_cells->coordinates[dim][cellId];
            pos++;
          }
        }
      }
    }
  }

}

void ZFSStrctrdBlck3D::sendPointsS()
{
  for(ZFSInt nghbr= m_cmnctnFlag->noNghbrDomainsNormal- m_cmnctnFlag->noNghbrDomainsSingular; nghbr< m_cmnctnFlag->noNghbrDomainsNormal; nghbr++){
    ZFSInt tag=domainId()+(m_cmnctnFlag->m_tagHelperSND[nghbr])*noDomains();
    ZFSInt err= MPI_Isend((void*)&m_cmnctnFlag->m_bufferPointsSnd[nghbr][0], m_cmnctnFlag->m_noNghbrDomainPointBufferSizeSnd[nghbr], MPI_DOUBLE, m_cmnctnFlag->m_sndNghbrId[nghbr], tag, MPI_COMM_WORLD, &m_cmnctnFlag->mpi_sndRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck3D::receivePointsS()
{
  for(ZFSInt nghbr= m_cmnctnFlag->noNghbrDomainsNormal- m_cmnctnFlag->noNghbrDomainsSingular; nghbr< m_cmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt tag=m_cmnctnFlag->m_rcvNghbrId[nghbr]+(m_cmnctnFlag->m_tagHelperRCV[nghbr])*noDomains();
    MPI_Irecv((void*)&m_cmnctnFlag->m_bufferPointsRcv[nghbr][0],m_cmnctnFlag->m_noNghbrDomainPointBufferSizeRcv[nghbr], MPI_DOUBLE,m_cmnctnFlag->m_rcvNghbrId[nghbr] ,tag, MPI_COMM_WORLD, &m_cmnctnFlag->mpi_rcvRequest[nghbr]);
  }
}

void ZFSStrctrdBlck3D::scatterPointsS()
{
  //  ZFSId pointId;
  ZFSInt cellId;
  //the ordering of the grid points can be different from
  //sending instance ==> reorder it and copy it to the
  //right place
  for(ZFSInt nghbr= m_cmnctnFlag->noNghbrDomainsNormal- m_cmnctnFlag->noNghbrDomainsSingular; nghbr< m_cmnctnFlag->noNghbrDomainsNormal; nghbr++)
  {
    ZFSInt k2, j2, i2, id2;
    ZFSInt* step1 = m_cmnctnFlag->stepInfoRCV[nghbr];
    ZFSInt  step2[3];

    ZFSInt* order = m_cmnctnFlag->orderInfo[nghbr];
    ZFSInt start1[3];
    //      ZFSInt end1[3];
    ZFSInt start2[3];
    ZFSInt end2[3];
    ZFSInt len2[3];
    ZFSInt totalCells=1;
    ZFSInt len1[3];
    for(ZFSInt j=0; j<nDim; j++)
    {
      len1[j]=m_cmnctnFlag->endInfoRCVpoints[nghbr][j] - m_cmnctnFlag->startInfoRCVpoints[nghbr][j];
      if(len1[j]!=0)  totalCells*=len1[j];

      //added    check the step for RCV part !!!!!!!!important
      step2[order[j]]=step1[j];
      //*******************************************************
    }
    for(ZFSInt j=0; j<nDim; j++)
    {
      start2[j]=0;
      end2[j]=len1[j]-1;
      len2[order[j]]=len1[j];
      if(step2[j]<0)
      {
        ZFSInt dummy=start2[j];
        start2[j]=end2[j];
        end2[j]=dummy;
      }
    }

    ZFSInt* startInfo=m_cmnctnFlag->startInfoRCVpoints[nghbr];
    ZFSInt* endInfo= m_cmnctnFlag->endInfoRCVpoints[nghbr];

    ZFSFloat* bufferRcv = m_cmnctnFlag->m_bufferPointsRcv[nghbr];
    for(ZFSId dim=0; dim<nDim; dim++)
    {
      k2=start2[2];
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++)
      {
        j2=start2[1];
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++)
        {
          i2=start2[0];
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++)
          {
            start1[order[0]]=i2;
            start1[order[1]]=j2;
            start1[order[2]]=k2;
            id2=dim*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
            cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            m_cells->coordinates[dim][cellId]=bufferRcv[id2];

            i2+=step2[0];
          }
          j2+=step2[1];
        }
        k2+=step2[2];
      }
    }
  }


}

//=======================================================================================
//=======================================================================================



/** brief parallel: coordinates the communication (exchange)
 *
 *
 */
void ZFSStrctrdBlck3D::exchange()
{
  RECORD_TIMER_START(m_tcomm);
  RECORD_TIMER_START(m_texchange);
  if(noDomains()>1){
    if(!m_nonBlockingComm){
      RECORD_TIMER_START(m_tgather);
      gather();
      RECORD_TIMER_STOP(m_tgather);

      RECORD_TIMER_START(m_tsend);
      send();
      RECORD_TIMER_STOP(m_tsend);

      RECORD_TIMER_START(m_treceive);
      receive();
      RECORD_TIMER_STOP(m_treceive);

      RECORD_TIMER_START(m_tsendWait);
      MPI_Waitall(m_cmnctnFlag->noNghbrDomainsNormal,m_cmnctnFlag->mpi_sndRequest,m_cmnctnFlag->mpi_sndStatus);
      RECORD_TIMER_STOP(m_tsendWait);

      RECORD_TIMER_START(m_treceiveWait);
      MPI_Waitall(m_cmnctnFlag->noNghbrDomainsNormal,m_cmnctnFlag->mpi_rcvRequest, m_cmnctnFlag->mpi_rcvStatus);
      RECORD_TIMER_STOP(m_treceiveWait);

      RECORD_TIMER_START(m_tscatter);
      scatter();
      RECORD_TIMER_STOP(m_tscatter);
    }else{
      nonBlockingExchange();
    }
  }

  if(m_periodicConnection) m_strctrdBndryCnd->periodicExchange();

  RECORD_TIMER_STOP(m_texchange);
  RECORD_TIMER_STOP(m_tcomm);
}


void ZFSStrctrdBlck3D::nonBlockingExchange(){
  //gather and send
  RECORD_TIMER_START(m_tgatherAndSend);
  RECORD_TIMER_START(m_tgatherAndSendWait);
  MPI_Waitall(m_cmnctnFlag->noNghbrDomainsNormal, m_cmnctnFlag->mpi_sndRequestCells, MPI_STATUSES_IGNORE);
  RECORD_TIMER_STOP(m_tgatherAndSendWait);
  for(ZFSId nghbr = 0; nghbr < m_cmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt* startInfo=m_cmnctnFlag->startInfoSNDcells[nghbr];
    ZFSInt* endInfo= m_cmnctnFlag->endInfoSNDcells[nghbr];
    ZFSFloat* bufferSnd = m_cmnctnFlag->m_bufferCellsSnd[nghbr];
    ZFSInt pos=0;
    for(ZFSId var=0; var<PV->noVariables; var++){
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++){
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++){
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++){
            const ZFSId cellId = cellIndex(i,j,k);
            bufferSnd[pos]=m_cells->pvariables[var][cellId];
            pos++;
          }
        }
      }
    }
    MPI_Start(&m_cmnctnFlag->mpi_sndRequestCells[nghbr]);
  }
  RECORD_TIMER_STOP(m_tgatherAndSend);
  //scatter
  RECORD_TIMER_START(m_tscatterAndReceive);
  ZFSIntScratchSpace array_of_indices(m_cmnctnFlag->noNghbrDomainsNormal, __CALLING_FUNCTION__, "array_of_indices");
  ZFSInt outcount;
  while(true){ //exit only in case of completed_index == MPI_UNDEFINED
    RECORD_TIMER_START(m_tscatterWaitSome);
    MPI_Waitsome(m_cmnctnFlag->noNghbrDomainsNormal, m_cmnctnFlag->mpi_rcvRequestCells, &outcount,array_of_indices.getPointer(), MPI_STATUSES_IGNORE);
    RECORD_TIMER_STOP(m_tscatterWaitSome);
    if(outcount == MPI_UNDEFINED) { break; }
    for(ZFSInt a = 0; a < outcount; a++) {
      ZFSInt nghbr = array_of_indices[a];
      ZFSInt k2, j2, i2, id2;
      ZFSInt* step1 = m_cmnctnFlag->stepInfoRCV[nghbr];
      ZFSInt  step2[3];
      ZFSInt* order = m_cmnctnFlag->orderInfo[nghbr];
      ZFSInt start1[3];
      ZFSInt start2[3];
      ZFSInt end2[3];
      ZFSInt len2[3];
      ZFSInt totalCells=1;
      ZFSInt len1[3];
      for(ZFSInt j=0; j<nDim; j++){
        len1[j]=m_cmnctnFlag->endInfoRCVcells[nghbr][j] - m_cmnctnFlag->startInfoRCVcells[nghbr][j];
        if(len1[j]!=0) totalCells*=len1[j];
        //added    check the step for RCV part !!!!!!!!important
        step2[order[j]]=step1[j];
        //*******************************************************
      }
      for(ZFSInt j=0; j<nDim; j++){
        start2[j]=0;
        end2[j]=len1[j]-1;
        len2[order[j]]=len1[j];
        if(step2[j]<0){
          ZFSInt dummy=start2[j];
          start2[j]=end2[j];
          end2[j]=dummy;
        }
      }
      ZFSInt* startInfo=m_cmnctnFlag->startInfoRCVcells[nghbr];
      ZFSInt* endInfo= m_cmnctnFlag->endInfoRCVcells[nghbr];
      ZFSFloat* bufferRcv = m_cmnctnFlag->m_bufferCellsRcv[nghbr];
      ZFSInt pos=0;
      for(ZFSId var=0; var<PV->noVariables; var++){
        k2=start2[2];
        for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++){
          j2=start2[1];
          for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++){
            i2=start2[0];
            for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++){
              start1[order[0]]=i2;
              start1[order[1]]=j2;
              start1[order[2]]=k2;
              id2=var*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
              ZFSId cellId = i +(j+k*m_nCells[1])*m_nCells[2];
              m_cells->variables[var][cellId]= bufferRcv[id2];

              i2+=step2[0];
              pos++;
            }
            j2+=step2[1];
          }
          k2+=step2[2];
        }
      }
    }
  }
  RECORD_TIMER_STOP(m_tscatterAndReceive);
  //open new receive
  RECORD_TIMER_START(m_treceive);
  MPI_Startall(m_cmnctnFlag->noNghbrDomainsNormal, m_cmnctnFlag->mpi_rcvRequestCells);
  RECORD_TIMER_STOP(m_treceive);
}



void ZFSStrctrdBlck3D::gather()
{
  ZFSId cellId;
  for(ZFSInt nghbr=0; nghbr<m_cmnctnFlag->noNghbrDomainsNormal; nghbr++)
  {
    ZFSInt* startInfo=m_cmnctnFlag->startInfoSNDcells[nghbr];
    ZFSInt* endInfo= m_cmnctnFlag->endInfoSNDcells[nghbr];
    ZFSFloat* bufferSnd = m_cmnctnFlag->m_bufferCellsSnd[nghbr];
    ZFSInt pos=0;

    for(ZFSId var=0; var<PV->noVariables; var++) {
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            bufferSnd[pos]=m_cells->pvariables[var][cellId];
            pos++;
          }
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::send()
{

  for(ZFSId nghbr=0; nghbr<m_cmnctnFlag->noNghbrDomainsNormal; nghbr++)
  {
    ZFSInt tag;
    tag= domainId()+(m_cmnctnFlag->m_tagHelperSND[nghbr])*noDomains();
    ZFSInt err = MPI_Isend((void*)&m_cmnctnFlag->m_bufferCellsSnd[nghbr][0], m_cmnctnFlag->m_noNghbrDomainCellBufferSizeSnd[nghbr], MPI_DOUBLE, m_cmnctnFlag->m_sndNghbrId[nghbr], tag, MPI_COMM_WORLD, &m_cmnctnFlag->mpi_sndRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }

}
void ZFSStrctrdBlck3D::receive()
{
  for(ZFSId nghbr=0; nghbr<m_cmnctnFlag->noNghbrDomainsNormal; nghbr++)
  {
    ZFSInt tag;
    tag=m_cmnctnFlag->m_rcvNghbrId[nghbr]+(m_cmnctnFlag->m_tagHelperRCV[nghbr])*noDomains();
    ZFSInt err = MPI_Irecv((void*)&m_cmnctnFlag->m_bufferCellsRcv[nghbr][0],m_cmnctnFlag->m_noNghbrDomainCellBufferSizeRcv[nghbr], MPI_DOUBLE,m_cmnctnFlag->m_rcvNghbrId[nghbr] ,tag, MPI_COMM_WORLD, &m_cmnctnFlag->mpi_rcvRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck3D::scatter()
{
  ZFSId cellId;
  //the ordering of the grid points can be different from
  //sending instance ==> reorder it and copy it to the
  //right place

  //for(ZFSInt nghbr=0; nghbr< m_cmnctnFlag->noNghbrDomainsNormal-m_cmnctnFlag->noNghbrDomainsSingular; nghbr++)
  for(ZFSInt nghbr=0; nghbr< m_cmnctnFlag->noNghbrDomainsNormal; nghbr++)
  {
    ZFSInt k2, j2, i2, id2;
    ZFSInt* step1 = m_cmnctnFlag->stepInfoRCV[nghbr];
    ZFSInt  step2[3];
    ZFSInt* order = m_cmnctnFlag->orderInfo[nghbr];
    ZFSInt start1[3];
    ZFSInt start2[3];
    ZFSInt end2[3];
    ZFSInt len2[3];
    ZFSInt totalCells=1;
    ZFSInt len1[3];

    for(ZFSInt j=0; j<nDim; j++) {
      len1[j]=m_cmnctnFlag->endInfoRCVcells[nghbr][j] - m_cmnctnFlag->startInfoRCVcells[nghbr][j];
      if(len1[j]!=0) totalCells*=len1[j];
      //added    check the step for RCV part !!!!!!!!important
      step2[order[j]]=step1[j];
      //*******************************************************
    }

    for(ZFSInt j=0; j<nDim; j++) {
      start2[j]=0;
      end2[j]=len1[j]-1;
      len2[order[j]]=len1[j];
      if(step2[j]<0) {
        ZFSInt dummy=start2[j];
        start2[j]=end2[j];
        end2[j]=dummy;
      }
    }

    ZFSInt* startInfo=m_cmnctnFlag->startInfoRCVcells[nghbr];
    ZFSInt* endInfo= m_cmnctnFlag->endInfoRCVcells[nghbr];

    ZFSFloat* bufferRcv = m_cmnctnFlag->m_bufferCellsRcv[nghbr];
    ZFSInt pos=0;
    for(ZFSId var=0; var<PV->noVariables; var++) {
      k2=start2[2];
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        j2=start2[1];
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          i2=start2[0];
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            start1[order[0]]=i2;
            start1[order[1]]=j2;
            start1[order[2]]=k2;

            id2=var*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
            cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            m_cells->pvariables[var][cellId]= bufferRcv[id2];

            i2+=step2[0];
            pos++;
          }
          j2+=step2[1];
        }
        k2+=step2[2];
      }
    }
  }
}

void ZFSStrctrdBlck3D::waveExchange()
{
  waveGather();
  waveSend();
  waveReceive();
  MPI_Waitall(m_waveCmnctnFlag->noNghbrDomainsSnd, m_waveCmnctnFlag->mpi_sndRequest, m_waveCmnctnFlag->mpi_sndStatus);
  MPI_Waitall(m_waveCmnctnFlag->noNghbrDomainsRcv, m_waveCmnctnFlag->mpi_rcvRequest, m_waveCmnctnFlag->mpi_rcvStatus);
  waveScatter();
}


void ZFSStrctrdBlck3D::waveGather()
{
  for(ZFSInt nghbr=0; nghbr<m_waveCmnctnFlag->noNghbrDomainsSnd; nghbr++) {
    ZFSInt* startInfo=m_waveCmnctnFlag->startInfoSNDcells[nghbr];
    ZFSInt* endInfo= m_waveCmnctnFlag->endInfoSNDcells[nghbr];
    ZFSFloat* bufferSnd = m_waveCmnctnFlag->m_bufferCellsSnd[nghbr];
    ZFSInt pos=0;

    for(ZFSId var=0; var<PV->noVariables; var++) {
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            const ZFSId cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            bufferSnd[pos]=m_cells->pvariables[var][cellId];
            pos++;
          }
        }
      }
    }

    if(m_averageVorticity) {
      for(ZFSId var=0; var<nDim; var++) {
        for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
          for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
            for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
              const ZFSId cellId = i +(j+k*m_nCells[1])*m_nCells[2];
              bufferSnd[pos]=m_cells->fq[FQ->VORTICITY[var]][cellId];
              pos++;
            }
          }
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::waveSend()
{
  for(ZFSId nghbr=0; nghbr<m_waveCmnctnFlag->noNghbrDomainsSnd; nghbr++) {
    const ZFSId tag= domainId()+(m_waveCmnctnFlag->m_tagHelperSND[nghbr])*noDomains();
    ZFSInt err = MPI_Isend((void*)&m_waveCmnctnFlag->m_bufferCellsSnd[nghbr][0],
                           m_waveCmnctnFlag->m_noNghbrDomainCellBufferSizeSnd[nghbr], MPI_DOUBLE,
                           m_waveCmnctnFlag->m_sndNghbrId[nghbr], tag, MPI_COMM_WORLD,
                           &m_waveCmnctnFlag->mpi_sndRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck3D::waveReceive()
{
  for(ZFSId nghbr=0; nghbr<m_waveCmnctnFlag->noNghbrDomainsRcv; nghbr++) {
    const ZFSId tag = m_waveCmnctnFlag->m_rcvNghbrId[nghbr]+(m_waveCmnctnFlag->m_tagHelperRCV[nghbr])*noDomains();
    ZFSInt err = MPI_Irecv((void*)&m_waveCmnctnFlag->m_bufferCellsRcv[nghbr][0],
                           m_waveCmnctnFlag->m_noNghbrDomainCellBufferSizeRcv[nghbr], MPI_DOUBLE,
                           m_waveCmnctnFlag->m_rcvNghbrId[nghbr] ,tag, MPI_COMM_WORLD,
                           &m_waveCmnctnFlag->mpi_rcvRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck3D::waveScatter()
{
  for(ZFSInt nghbr=0; nghbr< m_waveCmnctnFlag->noNghbrDomainsRcv; nghbr++) {
    ZFSInt* startInfo=m_waveCmnctnFlag->startInfoRCVcells[nghbr];
    ZFSInt* endInfo= m_waveCmnctnFlag->endInfoRCVcells[nghbr];

    ZFSFloat* bufferRcv = m_waveCmnctnFlag->m_bufferCellsRcv[nghbr];
    ZFSInt pos=0;

    for(ZFSId var=0; var<PV->noVariables; var++) {
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            const ZFSId cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            m_tempWaveSample[var][cellId]= bufferRcv[pos];
            pos++;
          }
        }
      }
    }

    if(m_averageVorticity) {
      for(ZFSId var=0; var<(2*nDim-3); var++) {
        for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
          for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
            for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
              const ZFSId cellId = i +(j+k*m_nCells[1])*m_nCells[2];
              m_tempWaveSample[PV->noVariables+var][cellId]= bufferRcv[pos];
              pos++;
            }
          }
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::spanwiseWaveReorder(){
  RECORD_TIMER_START(m_tspanwiseReorder);

  ZFSInt allCellsK=m_partition->inputBoxInfo[0]->DirLast[0];
  const ZFSInt waveZeroPos = ((globalTimeStep-m_movingGridStepOffset)/m_waveNoStepsPerCell)%allCellsK;
  m_windowInfo->createWaveWindowMapping(waveZeroPos);
  m_waveCmnctnFlag = new ZFSStrctrdWaveCommunicationHandle(nDim);

  ZFSInt noVars = PV->noVariables;
  if(m_averageVorticity) {
    noVars += (2*nDim-3);
    computeVorticity();
  }

  m_windowInfo->createWaveCommunicationExchangeFlags(m_waveCmnctnFlag, noVars);
  m_waveCmnctnFlag->setBufferSizes();
  waveExchange();
  delete m_waveCmnctnFlag;

  RECORD_TIMER_STOP(m_tspanwiseReorder);
}




void ZFSStrctrdBlck3D::ppFillGhostCells() {
  ppExtrapolateVariables();
  m_ppCmnctnFlag = new ZFSStrctrdCommunicationHandle(nDim);
  m_windowInfo->createCommunicationExchangeFlags(m_ppCmnctnFlag, getNoPPVars());
  m_ppCmnctnFlag->setBufferSizes();
  ppExchange();
  delete m_ppCmnctnFlag;
}


void ZFSStrctrdBlck3D::ppExtrapolateVariables()
{
  TRACE();

  //i-direction
  ZFSId cellId, cellIdAdj1, cellIdAdj2;
  for(ZFSId k=0; k<m_nCells[0]; k++) {
    for(ZFSId j=0; j<m_nCells[1]; j++ ) {
      for(ZFSId i=0; i<m_noGhostLayers; i++) {
        cellId = cellIndex(m_noGhostLayers-1-i,j,k); //pointId in Array
        cellIdAdj1 = cellIndex(m_noGhostLayers-i,j,k);
        cellIdAdj2 = cellIndex(m_noGhostLayers+1-i,j,k);
        for(ZFSId var =0; var < getNoPPVars(); var++) {
          m_summedVars[var][cellId]=(F2*m_summedVars[var][cellIdAdj1]-m_summedVars[var][cellIdAdj2]);
        }

        cellId = cellIndex(m_nCells[2]-m_noGhostLayers+i,j,k); //pointId in Array
        cellIdAdj1 = cellIndex(m_nCells[2]-m_noGhostLayers-1+i,j,k);
        cellIdAdj2 = cellIndex(m_nCells[2]-m_noGhostLayers-2+i,j,k);
        for(ZFSId var =0; var < getNoPPVars(); var++) {
          m_summedVars[var][cellId]=(F2*m_summedVars[var][cellIdAdj1]-m_summedVars[var][cellIdAdj2]);
        }
      }
    }
  }

  //j-direction
  for(ZFSId k=0; k<m_nCells[0]; k++) {
    for(ZFSId i=0; i<m_nCells[2]; i++ ) {
      for(ZFSId j=0; j<m_noGhostLayers; j++) {
        cellId = cellIndex(i,m_noGhostLayers-1-j,k); //pointId in Array
        cellIdAdj1 = cellIndex(i,m_noGhostLayers-j,k);
        cellIdAdj2 = cellIndex(i,m_noGhostLayers+1-j,k);
        for(ZFSId var =0; var < getNoPPVars(); var++) {
          m_summedVars[var][cellId]=(F2*m_summedVars[var][cellIdAdj1]-m_summedVars[var][cellIdAdj2]);
        }

        cellId = cellIndex(i,m_nCells[1]-m_noGhostLayers+j,k); //pointId in Array
        cellIdAdj1 = cellIndex(i,m_nCells[1]-m_noGhostLayers-1+j,k);
        cellIdAdj2 = cellIndex(i,m_nCells[1]-m_noGhostLayers-2+j,k);
        for(ZFSId var =0; var < getNoPPVars(); var++) {
          m_summedVars[var][cellId]=(F2*m_summedVars[var][cellIdAdj1]-m_summedVars[var][cellIdAdj2]);
        }
      }
    }
  }

  //k-direction
  for(ZFSId j=0; j<m_nCells[1]; j++) {
    for(ZFSId i=0; i<m_nCells[2]; i++ ) {
      for(ZFSId k=0; k<m_noGhostLayers; k++) {
        cellId = cellIndex(i,j,m_noGhostLayers-1-k); //pointId in Array
        cellIdAdj1 = cellIndex(i,j,m_noGhostLayers-k);
        cellIdAdj2 = cellIndex(i,j,m_noGhostLayers+1-k);
        for(ZFSId var =0; var < getNoPPVars(); var++) {
          m_summedVars[var][cellId]=(F2*m_summedVars[var][cellIdAdj1]-m_summedVars[var][cellIdAdj2]);
        }

        cellId = cellIndex(i,j,m_nCells[0]-m_noGhostLayers+k); //pointId in Array
        cellIdAdj1 = cellIndex(i,j,m_nCells[0]-m_noGhostLayers-1+k);
        cellIdAdj2 = cellIndex(i,j,m_nCells[0]-m_noGhostLayers-2+k);
        for(ZFSId var =0; var < getNoPPVars(); var++) {
          m_summedVars[var][cellId]=(F2*m_summedVars[var][cellIdAdj1]-m_summedVars[var][cellIdAdj2]);
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::ppExchange()
{
  if(noDomains()>1){
    ppGather();
    ppSend();
    ppReceive();
    MPI_Waitall(m_ppCmnctnFlag->noNghbrDomainsNormal,m_ppCmnctnFlag->mpi_sndRequest,m_ppCmnctnFlag->mpi_sndStatus);
    MPI_Waitall(m_ppCmnctnFlag->noNghbrDomainsNormal,m_ppCmnctnFlag->mpi_rcvRequest, m_ppCmnctnFlag->mpi_rcvStatus);
    ppScatter();
  }
}

void ZFSStrctrdBlck3D::ppGather()
{
  ZFSId cellId;
  for(ZFSInt nghbr=0; nghbr<m_ppCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt* startInfo=m_ppCmnctnFlag->startInfoSNDcells[nghbr];
    ZFSInt* endInfo= m_ppCmnctnFlag->endInfoSNDcells[nghbr];
    ZFSFloat* bufferSnd = m_ppCmnctnFlag->m_bufferCellsSnd[nghbr];
    ZFSInt pos=0;

    for(ZFSId var=0; var<getNoPPVars(); var++) {
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            bufferSnd[pos]=m_summedVars[var][cellId];
            pos++;
          }
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::ppSend()
{
  for(ZFSId nghbr=0; nghbr<m_ppCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt tag= domainId()+(m_ppCmnctnFlag->m_tagHelperSND[nghbr])*noDomains();
    ZFSInt err = MPI_Isend((void*)&m_ppCmnctnFlag->m_bufferCellsSnd[nghbr][0], m_ppCmnctnFlag->m_noNghbrDomainCellBufferSizeSnd[nghbr], MPI_DOUBLE, m_ppCmnctnFlag->m_sndNghbrId[nghbr], tag, MPI_COMM_WORLD, &m_ppCmnctnFlag->mpi_sndRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck3D::ppReceive()
{
  for(ZFSId nghbr=0; nghbr<m_ppCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt tag=m_ppCmnctnFlag->m_rcvNghbrId[nghbr]+(m_ppCmnctnFlag->m_tagHelperRCV[nghbr])*noDomains();
    ZFSInt err = MPI_Irecv((void*)&m_ppCmnctnFlag->m_bufferCellsRcv[nghbr][0],m_ppCmnctnFlag->m_noNghbrDomainCellBufferSizeRcv[nghbr], MPI_DOUBLE,m_ppCmnctnFlag->m_rcvNghbrId[nghbr] ,tag, MPI_COMM_WORLD, &m_ppCmnctnFlag->mpi_rcvRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck3D::ppScatter()
{
  ZFSId cellId;
  //the ordering of the grid points can be different from
  //sending instance ==> reorder it and copy it to the
  //right place

  //for(ZFSInt nghbr=0; nghbr< m_ppCmnctnFlag->noNghbrDomainsNormal-m_ppCmnctnFlag->noNghbrDomainsSingular; nghbr++)
  for(ZFSInt nghbr=0; nghbr< m_ppCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt k2, j2, i2, id2;
    ZFSInt* step1 = m_ppCmnctnFlag->stepInfoRCV[nghbr];
    ZFSInt  step2[3];
    ZFSInt* order = m_ppCmnctnFlag->orderInfo[nghbr];
    ZFSInt start1[3];
    ZFSInt start2[3];
    ZFSInt end2[3];
    ZFSInt len2[3];
    ZFSInt totalCells=1;
    ZFSInt len1[3];

    for(ZFSInt j=0; j<nDim; j++) {
      len1[j]=m_ppCmnctnFlag->endInfoRCVcells[nghbr][j] - m_ppCmnctnFlag->startInfoRCVcells[nghbr][j];
      if(len1[j]!=0) totalCells*=len1[j];
      //added    check the step for RCV part !!!!!!!!important
      step2[order[j]]=step1[j];
      //*******************************************************
    }

    for(ZFSInt j=0; j<nDim; j++) {
      start2[j]=0;
      end2[j]=len1[j]-1;
      len2[order[j]]=len1[j];
      if(step2[j]<0) {
        ZFSInt dummy=start2[j];
        start2[j]=end2[j];
        end2[j]=dummy;
      }
    }

    ZFSInt* startInfo=m_ppCmnctnFlag->startInfoRCVcells[nghbr];
    ZFSInt* endInfo= m_ppCmnctnFlag->endInfoRCVcells[nghbr];
    ZFSFloat* bufferRcv = m_ppCmnctnFlag->m_bufferCellsRcv[nghbr];
    ZFSInt pos=0;

    for(ZFSId var=0; var<getNoPPVars(); var++) {
      k2=start2[2];
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        j2=start2[1];
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          i2=start2[0];
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            start1[order[0]]=i2;
            start1[order[1]]=j2;
            start1[order[2]]=k2;

            id2=var*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
            cellId = i +(j+k*m_nCells[1])*m_nCells[2];
            m_summedVars[var][cellId]= bufferRcv[id2];

            i2+=step2[0];
            pos++;
          }
          j2+=step2[1];
        }
        k2+=step2[2];
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////zonal ghost cell averaged variables exchange//////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
void ZFSStrctrdBlck3D::averagedFillGhostCells() {
  averagedExtrapolateVariables();
  m_avCmnctnFlag = new ZFSStrctrdCommunicationHandle(nDim);
  m_windowInfo->createCommunicationExchangeFlags(m_avCmnctnFlag, 17);
  m_avCmnctnFlag->setBufferSizes();
  averagedExchange();
  delete m_avCmnctnFlag;
}
/**
 * Output from zonal methods
 * \author Junoh Jung, Master Thesis
 * \ date 07.2018 
 */


void ZFSStrctrdBlck3D::averagedExtrapolateVariables() //junoh
{
  TRACE();
  //i-direction
  ZFSId cellId, cellIdAdj1, cellIdAdj2;
  for(ZFSId k=0; k<m_nCells[0]; k++) {
    for(ZFSId j=0; j<m_nCells[1]; j++ ) {
      for(ZFSId i=0; i<m_noGhostLayers; i++) {
        cellId = cellIndex(m_noGhostLayers-1-i,j,k); //pointId in Array
        cellIdAdj1 = cellIndex(m_noGhostLayers-i,j,k);
        cellIdAdj2 = cellIndex(m_noGhostLayers+1-i,j,k);
	for(ZFSId var=0; var<m_maxNoVariables; var++){
        m_cells->pvariables[var][cellId]=(F2*m_cells->pvariables[var][cellIdAdj1]-m_cells->pvariables[var][cellIdAdj2]);
	}
	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);	  
	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);
	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);

	// m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

	// m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
	// m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	


	  
        cellId = cellIndex(m_nCells[2]-m_noGhostLayers+i,j,k); //pointId in Array
        cellIdAdj1 = cellIndex(m_nCells[2]-m_noGhostLayers-1+i,j,k);
        cellIdAdj2 = cellIndex(m_nCells[2]-m_noGhostLayers-2+i,j,k);

	for(ZFSId var=0; var<m_maxNoVariables; var++){
        m_cells->pvariables[var][cellId]=(F2*m_cells->pvariables[var][cellIdAdj1]-m_cells->pvariables[var][cellIdAdj2]);
	}
	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);
	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);
	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);
	
	// m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

	// m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
	// m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	// m_cells->fq[FQ->NUTILDE][cellId]=(F2*m_cells->fq[FQ->NUTILDE][cellIdAdj1]-m_cells->fq[FQ->NUTILDE][cellIdAdj2]);
	
      }
    }
  }

  //j-direction
  for(ZFSId k=0; k<m_nCells[0]; k++) {
    for(ZFSId i=0; i<m_nCells[2]; i++ ) {
      for(ZFSId j=0; j<m_noGhostLayers; j++) {
        cellId = cellIndex(i,m_noGhostLayers-1-j,k); //pointId in Array
        cellIdAdj1 = cellIndex(i,m_noGhostLayers-j,k);
        cellIdAdj2 = cellIndex(i,m_noGhostLayers+1-j,k);
	
	for(ZFSId var=0; var<m_maxNoVariables; var++){
        m_cells->pvariables[var][cellId]=(F2*m_cells->pvariables[var][cellIdAdj1]-m_cells->pvariables[var][cellIdAdj2]);
	}
	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);
	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);
	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);

	  
	// m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

	// m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
	// m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	


        cellId = cellIndex(i,m_nCells[1]-m_noGhostLayers+j,k); //pointId in Array
        cellIdAdj1 = cellIndex(i,m_nCells[1]-m_noGhostLayers-1+j,k);
        cellIdAdj2 = cellIndex(i,m_nCells[1]-m_noGhostLayers-2+j,k);
        
	for(ZFSId var=0; var<m_maxNoVariables; var++){
        m_cells->pvariables[var][cellId]=(F2*m_cells->pvariables[var][cellIdAdj1]-m_cells->pvariables[var][cellIdAdj2]);
	}
	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);
	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);	
	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);
	
	// m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

	// m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
	// m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	
      }
    }
  }

  //k-direction
  for(ZFSId j=0; j<m_nCells[1]; j++) {
    for(ZFSId i=0; i<m_nCells[2]; i++ ) {
      for(ZFSId k=0; k<m_noGhostLayers; k++) {
        cellId = cellIndex(i,j,m_noGhostLayers-1-k); //pointId in Array
        cellIdAdj1 = cellIndex(i,j,m_noGhostLayers-k);
        cellIdAdj2 = cellIndex(i,j,m_noGhostLayers+1-k);

	for(ZFSId var=0; var<m_maxNoVariables; var++){
        m_cells->pvariables[var][cellId]=(F2*m_cells->pvariables[var][cellIdAdj1]-m_cells->pvariables[var][cellIdAdj2]);
	}
	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);
	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);
	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);

	// m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

	// m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
	// m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	
	
        cellId = cellIndex(i,j,m_nCells[0]-m_noGhostLayers+k); //pointId in Array
        cellIdAdj1 = cellIndex(i,j,m_nCells[0]-m_noGhostLayers-1+k);
        cellIdAdj2 = cellIndex(i,j,m_nCells[0]-m_noGhostLayers-2+k);
	
	for(ZFSId var=0; var<m_maxNoVariables; var++){
        m_cells->pvariables[var][cellId]=(F2*m_cells->pvariables[var][cellIdAdj1]-m_cells->pvariables[var][cellIdAdj2]);
	}
	m_cells->fq[FQ->AVG_U][cellId]=(F2*m_cells->fq[FQ->AVG_U][cellIdAdj1]-m_cells->fq[FQ->AVG_U][cellIdAdj2]);
	m_cells->fq[FQ->AVG_V][cellId]=(F2*m_cells->fq[FQ->AVG_V][cellIdAdj1]-m_cells->fq[FQ->AVG_V][cellIdAdj2]);
	m_cells->fq[FQ->AVG_W][cellId]=(F2*m_cells->fq[FQ->AVG_W][cellIdAdj1]-m_cells->fq[FQ->AVG_W][cellIdAdj2]);
	m_cells->fq[FQ->AVG_RHO][cellId]=(F2*m_cells->fq[FQ->AVG_RHO][cellIdAdj1]-m_cells->fq[FQ->AVG_RHO][cellIdAdj2]);
	m_cells->fq[FQ->AVG_P][cellId]=(F2*m_cells->fq[FQ->AVG_P][cellIdAdj1]-m_cells->fq[FQ->AVG_P][cellIdAdj2]);
	m_cells->fq[FQ->NU_T][cellId]=(F2*m_cells->fq[FQ->NU_T][cellIdAdj1]-m_cells->fq[FQ->NU_T][cellIdAdj2]);
	
	// m_cells->fq[FQ->FLUC_U][cellId]=(F2*m_cells->fq[FQ->FLUC_U][cellIdAdj1]-m_cells->fq[FQ->FLUC_U][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_V][cellId]=(F2*m_cells->fq[FQ->FLUC_V][cellIdAdj1]-m_cells->fq[FQ->FLUC_V][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_W][cellId]=(F2*m_cells->fq[FQ->FLUC_W][cellIdAdj1]-m_cells->fq[FQ->FLUC_W][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UU][cellId]=(F2*m_cells->fq[FQ->FLUC_UU][cellIdAdj1]-m_cells->fq[FQ->FLUC_UU][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VV][cellId]=(F2*m_cells->fq[FQ->FLUC_VV][cellIdAdj1]-m_cells->fq[FQ->FLUC_VV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_WW][cellId]=(F2*m_cells->fq[FQ->FLUC_WW][cellIdAdj1]-m_cells->fq[FQ->FLUC_WW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UV][cellId]=(F2*m_cells->fq[FQ->FLUC_UV][cellIdAdj1]-m_cells->fq[FQ->FLUC_UV][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_UW][cellId]=(F2*m_cells->fq[FQ->FLUC_UW][cellIdAdj1]-m_cells->fq[FQ->FLUC_UW][cellIdAdj2]);
	// m_cells->fq[FQ->FLUC_VW][cellId]=(F2*m_cells->fq[FQ->FLUC_VW][cellIdAdj1]-m_cells->fq[FQ->FLUC_VW][cellIdAdj2]);

	// m_cells->fq[FQ->RECONST_NUT][cellId]=(F2*m_cells->fq[FQ->RECONST_NUT][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUT][cellIdAdj2]);
	// m_cells->fq[FQ->RECONST_NUTILDE][cellId]=(F2*m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj1]-m_cells->fq[FQ->RECONST_NUTILDE][cellIdAdj2]);
	
      }
    }
  }
}

/**
 * Output from zonal methods
 * \author Junoh Jung, Master Thesis
 * \ date 07.2018 
 */

void ZFSStrctrdBlck3D::averagedExchange()
{
  if(noDomains()>1){
    averagedGather();
    averagedSend();
    averagedReceive();
    MPI_Waitall(m_avCmnctnFlag->noNghbrDomainsNormal,m_avCmnctnFlag->mpi_sndRequest,m_avCmnctnFlag->mpi_sndStatus);
    MPI_Waitall(m_avCmnctnFlag->noNghbrDomainsNormal,m_avCmnctnFlag->mpi_rcvRequest, m_avCmnctnFlag->mpi_rcvStatus);
    averagedScatter();
  }
}

void ZFSStrctrdBlck3D::averagedGather()
{
  ZFSId cellId;
  for(ZFSInt nghbr=0; nghbr<m_avCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt* startInfo=m_avCmnctnFlag->startInfoSNDcells[nghbr];
    ZFSInt* endInfo= m_avCmnctnFlag->endInfoSNDcells[nghbr];
    ZFSFloat* bufferSnd = m_avCmnctnFlag->m_bufferCellsSnd[nghbr];
    ZFSInt pos=0;
    ZFSInt noCells = (endInfo[0]-startInfo[0])*(endInfo[1]-startInfo[1])*(endInfo[2]-startInfo[2]);
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            cellId = i +(j+k*m_nCells[1])*m_nCells[2];
	    // bufferSnd[pos]=m_cells->fq[var][cellId];
            bufferSnd[pos+FQ->AVG_U*noCells]=m_cells->fq[FQ->AVG_U][cellId];
	    bufferSnd[pos+FQ->AVG_V*noCells]=m_cells->fq[FQ->AVG_V][cellId];
	    bufferSnd[pos+FQ->AVG_W*noCells]=m_cells->fq[FQ->AVG_W][cellId];
	    bufferSnd[pos+FQ->AVG_RHO*noCells]=m_cells->fq[FQ->AVG_RHO][cellId];
	    bufferSnd[pos+FQ->AVG_P*noCells]=m_cells->fq[FQ->AVG_P][cellId];
	    bufferSnd[pos+FQ->NU_T*noCells]=m_cells->fq[FQ->NU_T][cellId];
	    
	    // bufferSnd[pos+6*noCells]=m_cells->fq[FQ->FLUC_U][cellId];
	    // bufferSnd[pos+7*noCells]=m_cells->fq[FQ->FLUC_V][cellId];
	    // bufferSnd[pos+8*noCells]=m_cells->fq[FQ->FLUC_W][cellId];
	    // bufferSnd[pos+9*noCells]=m_cells->fq[FQ->FLUC_UU][cellId];
	    // bufferSnd[pos+10*noCells]=m_cells->fq[FQ->FLUC_VV][cellId];
	    // bufferSnd[pos+11*noCells]=m_cells->fq[FQ->FLUC_WW][cellId];
	    // bufferSnd[pos+12*noCells]=m_cells->fq[FQ->FLUC_UV][cellId];
	    // bufferSnd[pos+13*noCells]=m_cells->fq[FQ->FLUC_UW][cellId];
	    // bufferSnd[pos+14*noCells]=m_cells->fq[FQ->FLUC_VW][cellId];
	    // bufferSnd[pos+15*noCells]=m_cells->fq[FQ->RECONST_NUT][cellId];
	    // bufferSnd[pos+16*noCells]=m_cells->fq[FQ->RECONST_NUTILDE][cellId];
	    pos++;
          
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::averagedSend()
{
  for(ZFSId nghbr=0; nghbr<m_avCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt tag= domainId()+(m_avCmnctnFlag->m_tagHelperSND[nghbr])*noDomains();
    ZFSInt err = MPI_Isend((void*)&m_avCmnctnFlag->m_bufferCellsSnd[nghbr][0], m_avCmnctnFlag->m_noNghbrDomainCellBufferSizeSnd[nghbr], MPI_DOUBLE, m_avCmnctnFlag->m_sndNghbrId[nghbr], tag, MPI_COMM_WORLD, &m_avCmnctnFlag->mpi_sndRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck3D::averagedReceive()
{
  for(ZFSId nghbr=0; nghbr<m_avCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt tag=m_avCmnctnFlag->m_rcvNghbrId[nghbr]+(m_avCmnctnFlag->m_tagHelperRCV[nghbr])*noDomains();
    ZFSInt err = MPI_Irecv((void*)&m_avCmnctnFlag->m_bufferCellsRcv[nghbr][0],m_avCmnctnFlag->m_noNghbrDomainCellBufferSizeRcv[nghbr], MPI_DOUBLE,m_avCmnctnFlag->m_rcvNghbrId[nghbr] ,tag, MPI_COMM_WORLD, &m_avCmnctnFlag->mpi_rcvRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck3D::averagedScatter()
{
  ZFSId cellId;
  //the ordering of the grid points can be different from
  //sending instance ==> reorder it and copy it to the
  //right place

  //for(ZFSInt nghbr=0; nghbr< m_ppCmnctnFlag->noNghbrDomainsNormal-m_ppCmnctnFlag->noNghbrDomainsSingular; nghbr++)
  for(ZFSInt nghbr=0; nghbr< m_avCmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    // ZFSInt k2, j2, i2, id2_U, id2_V, id2_W, id2_RHO, id2_P, id2_NU_T, id2_FLUC_U, id2_FLUC_V, id2_FLUC_W, id2_FLUC_UU, id2_FLUC_VV, id2_FLUC_WW, id2_FLUC_UV, id2_FLUC_UW, id2_FLUC_VW, id2_RECONST_NUT, id2_RECONST_NUTILDE ; //id2 junoh
    
    ZFSInt k2, j2, i2, id2_U, id2_V, id2_W, id2_RHO, id2_P, id2_NU_T;  //id2 junoh
    ZFSInt* step1 = m_avCmnctnFlag->stepInfoRCV[nghbr];
    ZFSInt  step2[3];
    ZFSInt* order = m_avCmnctnFlag->orderInfo[nghbr];
    ZFSInt start1[3];
    ZFSInt start2[3];
    ZFSInt end2[3];
    ZFSInt len2[3];
    ZFSInt totalCells=1;
    ZFSInt len1[3];

    for(ZFSInt j=0; j<nDim; j++) {
      len1[j]=m_avCmnctnFlag->endInfoRCVcells[nghbr][j] - m_avCmnctnFlag->startInfoRCVcells[nghbr][j];
      if(len1[j]!=0) totalCells*=len1[j];
      //added    check the step for RCV part !!!!!!!!important
      step2[order[j]]=step1[j];
      //*******************************************************
    }

    for(ZFSInt j=0; j<nDim; j++) {
      start2[j]=0;
      end2[j]=len1[j]-1;
      len2[order[j]]=len1[j];
      if(step2[j]<0) {
        ZFSInt dummy=start2[j];
        start2[j]=end2[j];
        end2[j]=dummy;
      }
    }

    ZFSInt* startInfo=m_avCmnctnFlag->startInfoRCVcells[nghbr];
    ZFSInt* endInfo= m_avCmnctnFlag->endInfoRCVcells[nghbr];
    ZFSFloat* bufferRcv = m_avCmnctnFlag->m_bufferCellsRcv[nghbr];
    ZFSInt pos=0;
    // ZFSInt noCells = (endInfo[0]-startInfo[0])*(endInfo[1]-startInfo[1])*(endInfo[2]-startInfo[2]);
     // for(ZFSId var=0; var<5; var++){
      k2=start2[2];
      for(ZFSInt k=startInfo[2]; k<endInfo[2]; k++) {
        j2=start2[1];
        for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
          i2=start2[0];
          for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
            start1[order[0]]=i2;
            start1[order[1]]=j2;
            start1[order[2]]=k2;
	    
	     // id2=var*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];

            id2_U=FQ->AVG_U*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    id2_V=FQ->AVG_V*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    id2_W=FQ->AVG_W*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    id2_RHO=FQ->AVG_RHO*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    id2_P=FQ->AVG_P*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];	    
	    id2_NU_T=FQ->NU_T*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_FLUC_U=6*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_FLUC_V=7*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_FLUC_W=8*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_FLUC_UU=9*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_FLUC_VV=10*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_FLUC_WW=11*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_FLUC_UV=12*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_FLUC_UW=13*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_FLUC_VW=14*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_RECONST_NUT=15*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    // id2_RECONST_NUTILDE=16*totalCells+start1[0]+(start1[1]+start1[2]*len2[1])*len2[0];
	    




            cellId = i +(j+k*m_nCells[1])*m_nCells[2];

	    // m_cells->fq[var][cellId]= bufferRcv[id2];

            m_cells->fq[FQ->AVG_U][cellId]= bufferRcv[id2_U];
	    m_cells->fq[FQ->AVG_V][cellId]= bufferRcv[id2_V];
	    m_cells->fq[FQ->AVG_W][cellId]= bufferRcv[id2_W];
	    m_cells->fq[FQ->AVG_RHO][cellId]= bufferRcv[id2_RHO];
	    m_cells->fq[FQ->AVG_P][cellId]= bufferRcv[id2_P];
	    m_cells->fq[FQ->NU_T][cellId]= bufferRcv[id2_NU_T];

	    // m_cells->fq[FQ->FLUC_U][cellId]= bufferRcv[id2_FLUC_U];
	    // m_cells->fq[FQ->FLUC_V][cellId]= bufferRcv[id2_FLUC_V];
	    // m_cells->fq[FQ->FLUC_W][cellId]= bufferRcv[id2_FLUC_W];
	    // m_cells->fq[FQ->FLUC_UU][cellId]= bufferRcv[id2_FLUC_UU];
	    // m_cells->fq[FQ->FLUC_VV][cellId]= bufferRcv[id2_FLUC_VV];
	    // m_cells->fq[FQ->FLUC_WW][cellId]= bufferRcv[id2_FLUC_WW];
	    // m_cells->fq[FQ->FLUC_UV][cellId]= bufferRcv[id2_FLUC_UV];
	    // m_cells->fq[FQ->FLUC_UW][cellId]= bufferRcv[id2_FLUC_UW];
	    // m_cells->fq[FQ->FLUC_VW][cellId]= bufferRcv[id2_FLUC_VW];
	    // m_cells->fq[FQ->RECONST_NUT][cellId]= bufferRcv[id2_RECONST_NUT];
	    // m_cells->fq[FQ->RECONST_NUTILDE][cellId]= bufferRcv[id2_RECONST_NUTILDE];
	   
            i2+=step2[0];
            pos++;
          }
          j2+=step2[1];
        }
        k2+=step2[2];
      }
      // }
  }
}


//====================================COMMUNICATIONS==========================================================
//============================================================================================================
//============================================================================================================




inline ZFSFloat ZFSStrctrdBlck3D::dist(ZFSFloat* a, ZFSFloat* b)
{
  ZFSFloat dist1=F0;
  for(ZFSId dim=0; dim<3; dim++)
  {
    dist1+=POW2(a[dim*m_noStrctrdCells]-b[dim*m_noStrctrdCells]);
  }
  return sqrt(dist1);
}

inline ZFSId ZFSStrctrdBlck3D::cellIndex(ZFSInt i, ZFSInt j, ZFSInt k)
{
  return i+(j+k*m_nCells[1])*m_nCells[2];
}

inline ZFSId ZFSStrctrdBlck3D::getCellIdfromCell( ZFSId origin, ZFSInt incI, ZFSInt incJ, ZFSInt incK )
{
  return origin + incI + incJ * m_nCells[2] + incK * m_nCells[2] * m_nCells[1];
}

inline ZFSId ZFSStrctrdBlck3D::pointIndex(ZFSInt i, ZFSInt j, ZFSInt k)
{
  return i+(j+k*m_nPoints[1])*m_nPoints[2];
}

void ZFSStrctrdBlck3D::viscousFlux(){
  RECORD_TIMER_START(m_tViscousFlux);
  (this->*viscFluxMethod)();
  RECORD_TIMER_STOP(m_tViscousFlux);
}

void ZFSStrctrdBlck3D::viscousFluxRANS(){
  m_ransBlck->viscousFluxRANS();
}

/*new viscous Flux calculation (exactly the same as in tfs) sss
 */
void ZFSStrctrdBlck3D::viscousFluxLES()
{
  TRACE();
  const ZFSFloat rPr = F1/m_Pr;
  const ZFSFloat rRe = F1/m_Re0;
  const ZFSFloat gammaMinusOne = m_gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;
  const ZFSInt noCells = m_noStrctrdCells;

  ZFSFloat* const RESTRICT u = ALIGNED_F(&m_cells->pvariables[PV->U][0]);
  ZFSFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  ZFSFloat* const RESTRICT w = &m_cells->pvariables[PV->W][0];
  ZFSFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  ZFSFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  ZFSFloat* const RESTRICT T = &m_cells->temperature[0];
  ZFSFloat* const RESTRICT lamvisc = &m_cells->lamvisc[0];

  ZFSFloat *const RESTRICT eflux= ALIGNED_MF(m_cells->eFlux);
  ZFSFloat *const RESTRICT fflux= ALIGNED_MF(m_cells->fFlux);
  ZFSFloat *const RESTRICT gflux= ALIGNED_MF(m_cells->gFlux);
  ZFSFloat *const RESTRICT vflux= ALIGNED_MF(m_cells->viscousFlux);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers+1; k++) {
    for(ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers+1; j++) {
      for(ZFSId ii=m_noGhostLayers-1; ii<m_nCells[2]-m_noGhostLayers+1; ii++) {
        const ZFSId I=cellIndex(ii,j,k);
        T[I] = m_gamma*p[I]/rho[I];
        lamvisc[I] = zfsSUTHERLANDLAW(T[I]);
      }
    }
  }

  ZFSFloat tau1, tau2, tau3, tau4, tau5, tau6;
  ZFSFloat dTdx, dTdy, dTdz;
  //ZFSFloat mueOverRe;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers+1; k++) {
    for(ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers+1; j++) {
      for(ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers+1; i++) {
        //get the adjacent cells;
        const ZFSId IJK   = cellIndex(i,j,k);
        const ZFSId IPJK  = cellIndex((i+1),j, k);
        const ZFSId IPJPK = cellIndex((i+1),(j+1), k);
        const ZFSId IJPK  = cellIndex(i,(j+1), k);
        const ZFSId IJKP  = cellIndex(i,j,(k+1));
        const ZFSId IPJKP = cellIndex((i+1),j,(k+1));
        const ZFSId IPJPKP= cellIndex((i+1),(j+1),(k+1));
        const ZFSId IJPKP = cellIndex(i,(j+1),(k+1));

	const ZFSFloat cornerMetrics[9] = {m_cells->cornerMetrics[0*noCells+IJK],
					   m_cells->cornerMetrics[1*noCells+IJK],
					   m_cells->cornerMetrics[2*noCells+IJK],
					   m_cells->cornerMetrics[3*noCells+IJK],
					   m_cells->cornerMetrics[4*noCells+IJK],
					   m_cells->cornerMetrics[5*noCells+IJK],
					   m_cells->cornerMetrics[6*noCells+IJK],
					   m_cells->cornerMetrics[7*noCells+IJK],
					   m_cells->cornerMetrics[8*noCells+IJK]};


        const ZFSFloat dudxi=F1B4*(u[IPJPKP]+u[IPJPK]+u[IPJKP]+u[IPJK]-u[IJPKP]-u[IJPK]-u[IJKP]-u[IJK]);
        const ZFSFloat dudet=F1B4*(u[IPJPKP]+u[IJPKP]+u[IPJPK]+u[IJPK]-u[IPJKP]-u[IJKP]-u[IPJK]-u[IJK]);
        const ZFSFloat dudze=F1B4*(u[IPJPKP]+u[IJPKP]+u[IPJKP]+u[IJKP]-u[IPJPK]-u[IJPK]-u[IPJK]-u[IJK]);

        const ZFSFloat dvdxi=F1B4*(v[IPJPKP]+v[IPJPK]+v[IPJKP]+v[IPJK]-v[IJPKP]-v[IJPK]-v[IJKP]-v[IJK]);
        const ZFSFloat dvdet=F1B4*(v[IPJPKP]+v[IJPKP]+v[IPJPK]+v[IJPK]-v[IPJKP]-v[IJKP]-v[IPJK]-v[IJK]);
        const ZFSFloat dvdze=F1B4*(v[IPJPKP]+v[IJPKP]+v[IPJKP]+v[IJKP]-v[IPJPK]-v[IJPK]-v[IPJK]-v[IJK]);

        const ZFSFloat dwdxi=F1B4*(w[IPJPKP]+w[IPJPK]+w[IPJKP]+w[IPJK]-w[IJPKP]-w[IJPK]-w[IJKP]-w[IJK]);
        const ZFSFloat dwdet=F1B4*(w[IPJPKP]+w[IJPKP]+w[IPJPK]+w[IJPK]-w[IPJKP]-w[IJKP]-w[IPJK]-w[IJK]);
        const ZFSFloat dwdze=F1B4*(w[IPJPKP]+w[IJPKP]+w[IPJKP]+w[IJKP]-w[IPJPK]-w[IJPK]-w[IPJK]-w[IJK]);

        const ZFSFloat dTdxi=F1B4*(T[IPJPKP]+T[IPJPK]+T[IPJKP]+T[IPJK]-T[IJPKP]-T[IJPK]-T[IJKP]-T[IJK]);
        const ZFSFloat dTdet=F1B4*(T[IPJPKP]+T[IJPKP]+T[IPJPK]+T[IJPK]-T[IPJKP]-T[IJKP]-T[IPJK]-T[IJK]);
        const ZFSFloat dTdze=F1B4*(T[IPJPKP]+T[IJPKP]+T[IPJKP]+T[IJKP]-T[IPJPK]-T[IJPK]-T[IPJK]-T[IJK]);

        const ZFSFloat uAvg=F1B8*(u[IPJPKP]+u[IJPKP]+u[IJPK]+u[IPJPK]+u[IPJKP]+u[IJKP]+u[IJK]+u[IPJK]);
        const ZFSFloat vAvg=F1B8*(v[IPJPKP]+v[IJPKP]+v[IJPK]+v[IPJPK]+v[IPJKP]+v[IJKP]+v[IJK]+v[IPJK]);
        const ZFSFloat wAvg=F1B8*(w[IPJPKP]+w[IJPKP]+w[IJPK]+w[IPJPK]+w[IPJKP]+w[IJKP]+w[IJK]+w[IPJK]);

	const ZFSFloat mue = F1B8*(lamvisc[IPJPKP]+lamvisc[IJPKP]+lamvisc[IJPK]+lamvisc[IPJPK]+
                                   lamvisc[IPJKP]+lamvisc[IJKP]+lamvisc[IJK]+lamvisc[IPJK]);

        // compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

        // tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
        //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
        //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
        tau1 = F4B3 * ( dudxi * cornerMetrics[ xsd * 3 + xsd ] +
                        dudet * cornerMetrics[ ysd * 3 + xsd ] +
                        dudze * cornerMetrics[ zsd * 3 + xsd ] ) -

          F2B3 * ( dvdxi * cornerMetrics[ xsd * 3 + ysd ] +
                   dvdet * cornerMetrics[ ysd * 3 + ysd ] +
                   dvdze * cornerMetrics[ zsd * 3 + ysd ] ) -

          F2B3 * ( dwdxi * cornerMetrics[ xsd * 3 + zsd ] +
                   dwdet * cornerMetrics[ ysd * 3 + zsd ] +
                   dwdze * cornerMetrics[ zsd * 3 + zsd ] );

        // compute tau2 = du/dy + dv/dx

        // tau_xy = du/dxi * dxi/dy + du/deta * deta/dy + du/dzeta * dzeta/dy
        //        + dv/dxi * dxi/dx + dv/deta * deta/dx + dv/dzeta * dzeta/dx
        tau2 = dudxi * cornerMetrics[ xsd * 3 + ysd ] +
          dudet * cornerMetrics[ ysd * 3 + ysd ] +
          dudze * cornerMetrics[ zsd * 3 + ysd ] +

          dvdxi * cornerMetrics[ xsd * 3 + xsd ] +
          dvdet * cornerMetrics[ ysd * 3 + xsd ] +
          dvdze * cornerMetrics[ zsd * 3 + xsd ];

        // compute tau3 = du/dz + dw/dx

        // tau_xz = du/dxi * dxi/dz + du/deta * deta/dz + du/dzeta * dzeta/dz
        //        + dw/dxi * dxi/dx + dw/deta * deta/dx + dw/dzeta * dzeta/dx
        tau3 = dudxi * cornerMetrics[ xsd * 3 + zsd ] +
          dudet * cornerMetrics[ ysd * 3 + zsd ] +
          dudze * cornerMetrics[ zsd * 3 + zsd ] +

          dwdxi * cornerMetrics[ xsd * 3 + xsd ] +
          dwdet * cornerMetrics[ ysd * 3 + xsd ] +
          dwdze * cornerMetrics[ zsd * 3 + xsd ];

        // compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

        // tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
        //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
        //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
        tau4 = F4B3 * ( dvdxi * cornerMetrics[ xsd * 3 + ysd ] +
                        dvdet * cornerMetrics[ ysd * 3 + ysd ] +
                        dvdze * cornerMetrics[ zsd * 3 + ysd ] ) -

          F2B3 * ( dudxi * cornerMetrics[ xsd * 3 + xsd ] +
                   dudet * cornerMetrics[ ysd * 3 + xsd ] +
                   dudze * cornerMetrics[ zsd * 3 + xsd ] ) -

          F2B3 * ( dwdxi * cornerMetrics[ xsd * 3 + zsd ] +
                   dwdet * cornerMetrics[ ysd * 3 + zsd ] +
                   dwdze * cornerMetrics[ zsd * 3 + zsd ] );

        // compute tau5 = dv/dz + dw/dy

        // tau_yz = dv/dxi * dxi/dz + dv/deta * deta/dz + dv/dzeta * dzeta/dz
        //        + dw/dxi * dxi/dy + dw/deta * deta/dy + dw/dzeta * dzeta/dy
        tau5 = dvdxi * cornerMetrics[ xsd * 3 + zsd ] +
          dvdet * cornerMetrics[ ysd * 3 + zsd ] +
          dvdze * cornerMetrics[ zsd * 3 + zsd ] +

          dwdxi * cornerMetrics[ xsd * 3 + ysd ] +
          dwdet * cornerMetrics[ ysd * 3 + ysd ] +
          dwdze * cornerMetrics[ zsd * 3 + ysd ];

        // compute tau6 = 2 dw/dz - 2/3 ( du/dx + dv/dy + dw/dz )

        // tau_zz = 4/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
        //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
        //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
        tau6 = F4B3 * ( dwdxi * cornerMetrics[ xsd * 3 + zsd ] +
                        dwdet * cornerMetrics[ ysd * 3 + zsd ] +
                        dwdze * cornerMetrics[ zsd * 3 + zsd ] ) -

          F2B3 * ( dudxi * cornerMetrics[ xsd * 3 + xsd ] +
                   dudet * cornerMetrics[ ysd * 3 + xsd ] +
                   dudze * cornerMetrics[ zsd * 3 + xsd ] ) -

          F2B3 * ( dvdxi * cornerMetrics[ xsd * 3 + ysd ] +
                   dvdet * cornerMetrics[ ysd * 3 + ysd ] +
                   dvdze * cornerMetrics[ zsd * 3 + ysd ] );


        dTdx = dTdxi * cornerMetrics[ xsd * 3 + xsd ] +
          dTdet * cornerMetrics[ ysd * 3 + xsd ] +
          dTdze * cornerMetrics[ zsd * 3 + xsd ];

        dTdy = dTdxi * cornerMetrics[ xsd * 3 + ysd ] +
          dTdet * cornerMetrics[ ysd * 3 + ysd ] +
          dTdze * cornerMetrics[ zsd * 3 + ysd ];

        dTdz = dTdxi * cornerMetrics[ xsd * 3 + zsd ] +
          dTdet * cornerMetrics[ ysd * 3 + zsd ] +
          dTdze * cornerMetrics[ zsd * 3 + zsd ];

        const ZFSFloat mueOverRe = mue * rRe / m_cells->cornerJac[IJK]; // divide by Jacobian
        tau1*=mueOverRe;
        tau2*=mueOverRe;
        tau3*=mueOverRe;
        tau4*=mueOverRe;
        tau5*=mueOverRe;
        tau6*=mueOverRe;
        const ZFSFloat mueH=FgammaMinusOne*mueOverRe*rPr;

        const ZFSFloat qx=mueH*dTdx+uAvg*tau1+vAvg*tau2+wAvg*tau3;
        const ZFSFloat qy=mueH*dTdy+uAvg*tau2+vAvg*tau4+wAvg*tau5;
        const ZFSFloat qz=mueH*dTdz+uAvg*tau3+vAvg*tau5+wAvg*tau6;


        //efluxes
        eflux[ 0*noCells+IJK ]    = tau1 * cornerMetrics[ xsd * 3 + xsd ] +
          tau2 * cornerMetrics[ xsd * 3 + ysd ] +
          tau3 * cornerMetrics[ xsd * 3 + zsd ];

        eflux[ 1*noCells+IJK ] = tau2 * cornerMetrics[ xsd * 3 + xsd ] +
          tau4 * cornerMetrics[ xsd * 3 + ysd ] +
          tau5 * cornerMetrics[ xsd * 3 + zsd ];

        eflux[ 2*noCells+IJK ] = tau3 * cornerMetrics[ xsd * 3 + xsd ] +
          tau5 * cornerMetrics[ xsd * 3 + ysd ] +
          tau6 * cornerMetrics[ xsd * 3 + zsd ];

        eflux[ 3*noCells+IJK ] = qx * cornerMetrics[ xsd * 3 + xsd ] +
          qy * cornerMetrics[ xsd * 3 + ysd ] +
          qz * cornerMetrics[ xsd * 3 + zsd ];

        //ffluxes
        fflux[ 0*noCells+IJK ]    = tau1 * cornerMetrics[ ysd * 3 + xsd ] +
          tau2 * cornerMetrics[ ysd * 3 + ysd ] +
          tau3 * cornerMetrics[ ysd * 3 + zsd ];

        fflux[ 1*noCells+IJK ] = tau2 * cornerMetrics[ ysd * 3 + xsd ] +
          tau4 * cornerMetrics[ ysd * 3 + ysd ] +
          tau5 * cornerMetrics[ ysd * 3 + zsd ];

        fflux[ 2*noCells+IJK ] = tau3 * cornerMetrics[ ysd * 3 + xsd ] +
          tau5 * cornerMetrics[ ysd * 3 + ysd ] +
          tau6 * cornerMetrics[ ysd * 3 + zsd ];

        fflux[ 3*noCells+IJK ] = qx * cornerMetrics[ ysd * 3 + xsd ] +
          qy * cornerMetrics[ ysd * 3 + ysd ] +
          qz * cornerMetrics[ ysd * 3 + zsd ];

        //gfluxes
        gflux[ 0*noCells+IJK ]    = tau1 * cornerMetrics[ zsd * 3 + xsd ] +
          tau2 * cornerMetrics[ zsd * 3 + ysd ] +
          tau3 * cornerMetrics[ zsd * 3 + zsd ];

        gflux[ 1*noCells+IJK ] = tau2 * cornerMetrics[ zsd * 3 + xsd ] +
          tau4 * cornerMetrics[ zsd * 3 + ysd ] +
          tau5 * cornerMetrics[ zsd * 3 + zsd ];

        gflux[ 2*noCells+IJK ] = tau3 * cornerMetrics[ zsd * 3 + xsd ] +
          tau5 * cornerMetrics[ zsd * 3 + ysd ] +
          tau6 * cornerMetrics[ zsd * 3 + zsd ];

        gflux[ 3*noCells+IJK ] = qx * cornerMetrics[ zsd * 3 + xsd ] +
          qy * cornerMetrics[ zsd * 3 + ysd ] +
          qz * cornerMetrics[ zsd * 3 + zsd ];

      }
    }
  }


  //viscous flux correction for the singular points
  //m_hasSingularity=0 means no singular points in this block, otherwise do flux correction
  if(m_hasSingularity>0){
    viscousFluxCorrection();
  }

  for(ZFSId var=0; var<(CV->noVariables)-1; var++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
      for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
        for(ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers; i++) {
          const ZFSId IJK    = cellIndex(i,j,k);
          const ZFSId IJMK   = cellIndex(i,(j-1),k);
          const ZFSId IJKM   = cellIndex(i,j,(k-1));
          const ZFSId IJMKM  = cellIndex(i,(j-1),(k-1));

          vflux[0*noCells+IJK]=F1B4*(eflux[var*noCells+IJK]+
                                     eflux[var*noCells+IJKM]+
                                     eflux[var*noCells+IJMK]+
                                     eflux[var*noCells+IJMKM]);

#ifdef ZFS_EXTRA_DEBUG
          viscFluxOut[0][var*noCells+IJK]= vflux[3*IJK];
#endif
        }
      }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
      for(ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers; j++) {
        for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
          const ZFSId IJK   = cellIndex(i,j,k);
          const ZFSId IMJK  = cellIndex((i-1),j,k);
          const ZFSId IJKM  = cellIndex(i,j,(k-1));
          const ZFSId IMJKM = cellIndex((i-1),j,(k-1));

          vflux[1*noCells+IJK]=F1B4*(fflux[var*noCells+IJK]+
                                     fflux[var*noCells+IJKM]+
                                     fflux[var*noCells+IMJK]+
                                     fflux[var*noCells+IMJKM]);

#ifdef ZFS_EXTRA_DEBUG
          viscFluxOut[1][var*noCells+IJK]= vflux[3*IJK+1];
#endif
        }
      }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers; k++) {
      for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
        for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
          const ZFSId IJK   = cellIndex(i,j,k);
          const ZFSId IMJK  = cellIndex((i-1),j,k);
          const ZFSId IJMK  = cellIndex(i,(j-1),k);
          const ZFSId IMJMK = cellIndex((i-1),(j-1),k);

          vflux[2*noCells+IJK]=F1B4*(gflux[var*noCells+IJK]+
                                     gflux[var*noCells+IMJK]+
                                     gflux[var*noCells+IJMK]+
                                     gflux[var*noCells+IMJMK]);

#ifdef ZFS_EXTRA_DEBUG
          viscFluxOut[2][var*m_noStrctrdCells+IJK]= vflux[3*IJK+2];
#endif
        }
      }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
      for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
        for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
          const ZFSId IJK=cellIndex(i,j,k);
          const ZFSId IMJK=cellIndex(i-1,j,k);
          const ZFSId IJMK=cellIndex(i,j-1,k);
          const ZFSId IJKM=cellIndex(i,j,k-1);
          m_cells->rightHandSide[var][IJK]+=  vflux[0*noCells+IJK]-vflux[0*noCells+IMJK]
                                             +vflux[1*noCells+IJK]-vflux[1*noCells+IJMK]
                                             +vflux[2*noCells+IJK]-vflux[2*noCells+IJKM];
        }
      }
    }
  }
}


void ZFSStrctrdBlck3D::viscousFluxCorrection()
{
  const ZFSInt noCells = m_noStrctrdCells;
  const ZFSFloat rPr = F1/m_Pr;
  const ZFSFloat rRe = F1 / m_Re0;
  const ZFSFloat gammaMinusOne = m_gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  ZFSFloat* __restrict rhou=&m_cells->variables[CV->RHO_U][0];
  ZFSFloat* __restrict rhov=&m_cells->variables[CV->RHO_V][0];
  ZFSFloat* __restrict rhow=&m_cells->variables[CV->RHO_W][0];
  ZFSFloat* __restrict rhoE=&m_cells->variables[CV->RHO_E][0];
  ZFSFloat* __restrict rho=&m_cells->variables[CV->RHO][0];

  ZFSFloat* __restrict eflux=&m_cells->eFlux[0];
  ZFSFloat* __restrict fflux=&m_cells->fFlux[0];
  ZFSFloat* __restrict gflux=&m_cells->gFlux[0];

  ZFSInt dim;
  ZFSInt start[3],end[3],nghbr[20];
  ZFSInt len1[3];
  ZFSInt totalCells;

  for(ZFSInt i=0; i<m_hasSingularity; ++i) {
    //only correct for bc 6000 not for bc 4000-5000
    if(m_singularity[i].BC==6000) {
      totalCells=1;
      for(ZFSInt j=0; j<nDim; j++) {
        len1[j]=m_singularity[i].end[j]-m_singularity[i].start[j];
        if(len1[j]!=0)  totalCells*=len1[j];
      }

      for( ZFSInt n = 0; n < 3 ; ++n ) {
        if(m_singularity[i].end[n]-m_singularity[i].start[n]>1) {
          dim=n;
          // start[n]=m_singularity[i].start[n]+1;
          start[n]=m_singularity[i].start[n]+1;
          end[n]=m_singularity[i].end[n]-1;
        } else {
          start[n]=m_singularity[i].start[n];
          end[n]=m_singularity[i].end[n];
        }
      }

      ZFSFloat u[20],v[20],w[20],T[20];
      ZFSFloat U,V,W,t,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dTdx,dTdy,dTdz;
      ZFSFloat qx, qy, qz;
      ZFSFloat tau1, tau2, tau3, tau4, tau5, tau6;
      ZFSFloat mueOverRe,mue,mueH;
      for( ZFSInt kk = start[2]; kk <end[2]; ++kk ) {
        for( ZFSInt jj = start[1]; jj <end[1]; ++jj ) {
          for( ZFSInt ii = start[0]; ii <end[0]; ++ii ) {
            ZFSInt count=0;
            ZFSInt temp[3]={0,0,0};
            ZFSInt IJK = cellIndex( ii+ m_singularity[i].Viscous[0], jj+ m_singularity[i].Viscous[1], kk+ m_singularity[i].Viscous[2] );

            const ZFSFloat cornerMetrics[9] = {m_cells->cornerMetrics[0*noCells+IJK],
                                               m_cells->cornerMetrics[1*noCells+IJK],
                                               m_cells->cornerMetrics[2*noCells+IJK],
                                               m_cells->cornerMetrics[3*noCells+IJK],
                                               m_cells->cornerMetrics[4*noCells+IJK],
                                               m_cells->cornerMetrics[5*noCells+IJK],
                                               m_cells->cornerMetrics[6*noCells+IJK],
                                               m_cells->cornerMetrics[7*noCells+IJK],
                                               m_cells->cornerMetrics[8*noCells+IJK]};
	    
            temp[dim]=1;
            nghbr[count++]=cellIndex(ii,jj,kk);
            nghbr[count++]=cellIndex(ii+temp[0],jj+temp[1],kk+temp[2]);

            for(ZFSInt m=0; m< m_singularity[i].Nstar-1; ++m) {
              ZFSInt *change=  m_singularity[i].displacement[m];
              nghbr[count++]=cellIndex(ii+change[0],jj+change[1],kk+change[2]);
              nghbr[count++]=cellIndex(ii+temp[0]+change[0],jj+temp[1]+change[1],kk+temp[2]+change[2]);
            }

            if(count!=m_singularity[i].Nstar*2) {
              cout<<"what the hell! it is wrong!!!"<<endl;
            }

            for(ZFSInt m=0;m<m_singularity[i].Nstar*2;++m) {
              u[m]=rhou[nghbr[m]]/rho[nghbr[m]];
              v[m]=rhov[nghbr[m]]/rho[nghbr[m]];
              w[m]=rhow[nghbr[m]]/rho[nghbr[m]];
              T[m]=(m_gamma*gammaMinusOne*(rhoE[nghbr[m]]-F1B2*rho[nghbr[m]]*(POW2(u[m])+POW2(v[m])+POW2(w[m]))))/rho[nghbr[m]];
            }

            U=F0;V=F0;W=F0;t=F0;
            dudx=F0;dudy=F0;dudz=F0;
            dvdx=F0;dvdy=F0;dvdz=F0;
            dwdx=F0;dwdy=F0;dwdz=F0;
            dTdx=F0;dTdy=F0;dTdz=F0;

            ZFSInt id2=ii-start[0]+((jj-start[1])+(kk-start[2])*len1[1])*len1[0];

            for(ZFSInt n=0; n<count;n++) {
              ZFSInt ID=id2*count+n;
              U   += m_singularity[i].ReconstructionConstants[ 0 ][ ID ]*u[n];
              dudx+= m_singularity[i].ReconstructionConstants[ 1 ][ ID ]*u[n];
              dudy+= m_singularity[i].ReconstructionConstants[ 2 ][ ID ]*u[n];
              dudz+= m_singularity[i].ReconstructionConstants[ 3 ][ ID ]*u[n];

              V   += m_singularity[i].ReconstructionConstants[ 0 ][ ID ]*v[n];
              dvdx+= m_singularity[i].ReconstructionConstants[ 1 ][ ID ]*v[n];
              dvdy+= m_singularity[i].ReconstructionConstants[ 2 ][ ID ]*v[n];
              dvdz+= m_singularity[i].ReconstructionConstants[ 3 ][ ID ]*v[n];

              W   += m_singularity[i].ReconstructionConstants[ 0 ][ ID ]*w[n];
              dwdx+= m_singularity[i].ReconstructionConstants[ 1 ][ ID ]*w[n];
              dwdy+= m_singularity[i].ReconstructionConstants[ 2 ][ ID ]*w[n];
              dwdz+= m_singularity[i].ReconstructionConstants[ 3 ][ ID ]*w[n];

              t   += m_singularity[i].ReconstructionConstants[ 0 ][ ID ]*T[n];
              dTdx+= m_singularity[i].ReconstructionConstants[ 1 ][ ID ]*T[n];
              dTdy+= m_singularity[i].ReconstructionConstants[ 2 ][ ID ]*T[n];
              dTdz+= m_singularity[i].ReconstructionConstants[ 3 ][ ID ]*T[n];
            }

            tau1 = 2*dudx - 2/3*( dudx + dvdy + dwdz );
            tau2 = dudy + dvdx;
            tau3 = dudz + dwdx;
            tau4 = 2*dvdy - 2/3*( dudx + dvdy + dwdz );
            tau5 = dvdz + dwdy;
            tau6 = 2*dwdz - 2/3*( dudx + dvdy + dwdz );

            mue=zfsSUTHERLANDLAW(t);
            mueOverRe = mue * rRe;
            tau1*=mueOverRe;
            tau2*=mueOverRe;
            tau3*=mueOverRe;
            tau4*=mueOverRe;
            tau5*=mueOverRe;
            tau6*=mueOverRe;
            mueH=FgammaMinusOne*mueOverRe*rPr;

            qx=mueH*dTdx+U*tau1+V*tau2+W*tau3;
            qy=mueH*dTdy+U*tau2+V*tau4+W*tau5;
            qz=mueH*dTdz+U*tau3+V*tau5+W*tau6;


            //efluxes
            eflux[ 0*noCells+IJK ]    = tau1 * cornerMetrics[ xsd * 3 + xsd ] +
              tau2 * cornerMetrics[ xsd * 3 + ysd ] +
              tau3 * cornerMetrics[ xsd * 3 + zsd ];

            eflux[ 1*noCells+IJK ] = tau2 * cornerMetrics[ xsd * 3 + xsd ] +
              tau4 * cornerMetrics[ xsd * 3 + ysd ] +
              tau5 * cornerMetrics[ xsd * 3 + zsd ];

            eflux[ 2*noCells+IJK ] = tau3 * cornerMetrics[ xsd * 3 + xsd ] +
              tau5 * cornerMetrics[ xsd * 3 + ysd ] +
              tau6 * cornerMetrics[ xsd * 3 + zsd ];

            eflux[ 3*noCells+IJK ] = qx * cornerMetrics[ xsd * 3 + xsd ] +
              qy * cornerMetrics[ xsd * 3 + ysd ] +
              qz * cornerMetrics[ xsd * 3 + zsd ];

            //ffluxes
            fflux[ 0*noCells+IJK ]    = tau1 * cornerMetrics[ ysd * 3 + xsd ] +
              tau2 * cornerMetrics[ ysd * 3 + ysd ] +
              tau3 * cornerMetrics[ ysd * 3 + zsd ];

            fflux[ 1*noCells+IJK ] = tau2 * cornerMetrics[ ysd * 3 + xsd ] +
              tau4 * cornerMetrics[ ysd * 3 + ysd ] +
              tau5 * cornerMetrics[ ysd * 3 + zsd ];

            fflux[ 2*noCells+IJK ] = tau3 * cornerMetrics[ ysd * 3 + xsd ] +
              tau5 * cornerMetrics[ ysd * 3 + ysd ] +
              tau6 * cornerMetrics[ ysd * 3 + zsd ];

            fflux[ 3*noCells+IJK ] = qx * cornerMetrics[ ysd * 3 + xsd ] +
              qy * cornerMetrics[ ysd * 3 + ysd ] +
              qz * cornerMetrics[ ysd * 3 + zsd ];

            //gfluxes
            gflux[ 0*noCells+IJK ]    = tau1 * cornerMetrics[ zsd * 3 + xsd ] +
              tau2 * cornerMetrics[ zsd * 3 + ysd ] +
              tau3 * cornerMetrics[ zsd * 3 + zsd ];

            gflux[ 1*noCells+IJK ] = tau2 * cornerMetrics[ zsd * 3 + xsd ] +
              tau4 * cornerMetrics[ zsd * 3 + ysd ] +
              tau5 * cornerMetrics[ zsd * 3 + zsd ];

            gflux[ 2*noCells+IJK ] = tau3 * cornerMetrics[ zsd * 3 + xsd ] +
              tau5 * cornerMetrics[ zsd * 3 + ysd ] +
              tau6 * cornerMetrics[ zsd * 3 + zsd ];

            gflux[ 3*noCells+IJK ] = qx * cornerMetrics[ zsd * 3 + xsd ] +
              qy * cornerMetrics[ zsd * 3 + ysd ] +
              qz * cornerMetrics[ zsd * 3 + zsd ];
          }
        }
      }
    }
  }
}


void ZFSStrctrdBlck3D::computeSurfaceMetrics()
{
    TRACE();
    zfs_log << "computing surface metrics ... " << endl;

  for( ZFSInt k = 0; k < this->m_nCells[0]; k++ )
  {
    for( ZFSInt j = 0; j < this->m_nCells[1]; j++ )
    {
      for( ZFSInt i = 0; i < this->m_nCells[2]; i++ )
      {
        // determine global cell ID
        ZFSId cellId = this->cellIndex(i,j,k); //i + ( k * m_nCells[1] + j ) * m_nCells[2];
        // determine global point ID for local cell IDs
        ZFSId ijk    = getPointIdFromCell( i, j, k );
        ZFSId ipjk   = getPointIdfromPoint( ijk, 1, 0, 0 );
        ZFSId ipjpk  = getPointIdfromPoint( ijk, 1, 1, 0 );
        ZFSId ipjkp  = getPointIdfromPoint( ijk, 1, 0, 1 );
        ZFSId ipjpkp = getPointIdfromPoint( ijk, 1, 1, 1 );
        ZFSId ijpk   = getPointIdfromPoint( ijk, 0, 1, 0 );
        ZFSId ijpkp  = getPointIdfromPoint( ijk, 0, 1, 1 );
        ZFSId ijkp   = getPointIdfromPoint( ijk, 0, 0, 1 );

        //auxilliary variables
        ZFSFloat DcoordDxi[ 3 ];
        ZFSFloat DcoordDeta[ 3 ];
        ZFSFloat DcoordDzeta[ 3 ];

        ZFSFloat metricTmp[ 3 ];

        //////////////////////////////////////////////////////
        ////////////////////// FACE I ////////////////////////
        //////////////////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          DcoordDeta[isd]  = ( ( m_coordinates[isd][ipjpkp] + m_coordinates[isd][ipjpk] ) -
                               ( m_coordinates[isd][ipjkp] + m_coordinates[isd][ipjk] ) ) *
            F1B2;
          DcoordDzeta[isd] = ( ( m_coordinates[isd][ipjpkp] + m_coordinates[isd][ipjkp] ) -
                               ( m_coordinates[isd][ipjpk] + m_coordinates[isd][ipjk] ) ) *
            F1B2;
        }

        // compute Dxi and store
        this->crossProduct( metricTmp, DcoordDeta, DcoordDzeta );

        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          m_cells->surfaceMetrics[cellId][ xsd * nDim + isd ] = metricTmp[isd];
        }

        ///////////////////////////////////////////////////////
        ////////////////////// FACE J /////////////////////////
        ///////////////////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          DcoordDxi[isd]   = ( ( m_coordinates[isd][ipjpkp] + m_coordinates[isd][ipjpk] ) -
                               ( m_coordinates[isd][ijpkp] + m_coordinates[isd][ijpk] ) ) *
            F1B2;

          DcoordDzeta[isd] = ( ( m_coordinates[isd][ijpkp] + m_coordinates[isd][ipjpkp] ) -
                               ( m_coordinates[isd][ipjpk] + m_coordinates[isd][ijpk] ) ) *
            F1B2;
        }

        // compute Deta and store
        this->crossProduct( metricTmp, DcoordDzeta, DcoordDxi );

        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          m_cells->surfaceMetrics[cellId][ ysd * nDim + isd ] = metricTmp[isd];
        }

        ///////////////////////////////////////////////////////
        ////////////////////// FACE K /////////////////////////
        ///////////////////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          DcoordDxi[isd]   = ( ( m_coordinates[isd][ipjpkp] + m_coordinates[isd][ipjkp] ) -
                               ( m_coordinates[isd][ijpkp] + m_coordinates[isd][ijkp] ) ) *
            F1B2;

          DcoordDeta[isd]  = ( ( m_coordinates[isd][ipjpkp] + m_coordinates[isd][ijpkp] ) -
                               ( m_coordinates[isd][ipjkp] + m_coordinates[isd][ijkp] ) ) *
            F1B2;
        }

        // compute Dzeta and store
        this->crossProduct( metricTmp, DcoordDxi, DcoordDeta );

        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          m_cells->surfaceMetrics[cellId][ zsd * nDim + isd ] = metricTmp[isd];
        }
      }
    }
  }
    zfs_log << "computing surface metrics ... SUCCESSFUL" << endl;
}

/** brief Computes surface metrics
 *  similar routine as DXNC in TFS adapted
 *  for a cell-centered FV scheme
 *  more exact calculation of the metrics
 *  than with standard method by splitting each
 *  surface into 4 sub-surfaces.
 *
 */
void ZFSStrctrdBlck3D::computeModSurfaceMetrics()
{
  TRACE();
  ZFSFloat diagonal1[3];
  ZFSFloat diagonal2[3];

  ZFSFloat subMetric1[3];
  ZFSFloat subMetric2[3];
  ZFSFloat subMetric3[3];
  ZFSFloat subMetric4[3];

  ZFSFloat p1[3];
  ZFSFloat p2[3];
  ZFSFloat p3[3];
  ZFSFloat p4[3];
  ZFSFloat center[3];

  for( ZFSInt k = 0; k < this->m_nCells[0]; k++ )
  {
    for( ZFSInt j = 0; j < this->m_nCells[1]; j++ )
    {
      for( ZFSInt i = 0; i < this->m_nCells[2]; i++ )
      {
        // determine global cell ID
        const ZFSId cellId = this->cellIndex(i,j,k); //i + ( k * m_nCells[1] + j ) * m_nCells[2];
        // determine global point ID for local cell IDs
        const ZFSId ijk    = getPointIdFromCell( i, j, k );
        const ZFSId ipjk   = getPointIdfromPoint( ijk, 1, 0, 0 );
        const ZFSId ipjpk  = getPointIdfromPoint( ijk, 1, 1, 0 );
        const ZFSId ipjkp  = getPointIdfromPoint( ijk, 1, 0, 1 );
        const ZFSId ipjpkp = getPointIdfromPoint( ijk, 1, 1, 1 );
        const ZFSId ijpk   = getPointIdfromPoint( ijk, 0, 1, 0 );
        const ZFSId ijpkp  = getPointIdfromPoint( ijk, 0, 1, 1 );
        const ZFSId ijkp   = getPointIdfromPoint( ijk, 0, 0, 1 );

        ///////////////////////////////////////////////////////
        ////////////////////// FACE I /////////////////////////
        ///////////////////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          //edge centers
          p1[isd] = (m_coordinates[isd][ipjpk] + m_coordinates[isd][ipjk])/F2;
          p2[isd] = (m_coordinates[isd][ipjpkp] + m_coordinates[isd][ipjkp])/F2;
          p3[isd] = (m_coordinates[isd][ipjpkp] + m_coordinates[isd][ipjpk])/F2;
          p4[isd] = (m_coordinates[isd][ipjkp] + m_coordinates[isd][ipjk])/F2;

          //surface centroid
          center[isd] = (m_coordinates[isd][ipjk] + m_coordinates[isd][ipjpk] + m_coordinates[isd][ipjkp] + m_coordinates[isd][ipjpkp])/F4;
        }

        ////////////////
        // Submetric 1//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = p1[isd] - p4[isd];
          diagonal2[isd] = center[isd] - m_coordinates[isd][ipjk];
        }

        // compute first subMetric
        this->crossProduct( subMetric1, diagonal1, diagonal2 );

        ////////////////
        // Submetric 2//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = m_coordinates[isd][ipjpk] - center[isd];
          diagonal2[isd] = p3[isd] - p1[isd];
        }

        // compute second subMetric
        this->crossProduct( subMetric2, diagonal1, diagonal2 );

        ////////////////
        // Submetric 3//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = p3[isd] - p2[isd];
          diagonal2[isd] = m_coordinates[isd][ipjpkp] - center[isd];
        }

        // compute third subMetric
        this->crossProduct( subMetric3, diagonal1, diagonal2 );

        ////////////////
        // Submetric 4//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = center[isd] - m_coordinates[isd][ipjkp];
          diagonal2[isd] = p2[isd] - p4[isd];
        }

        // compute fourth subMetric
        this->crossProduct( subMetric4, diagonal1, diagonal2 );

        // assemble subMetrics
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          // ZFSFloat oldSurfMetric = m_cells->surfaceMetrics[cellId][isd];
          // ZFSFloat newSurfMetric = (subMetric1[isd] + subMetric2[isd] + subMetric3[isd] + subMetric4[isd])/F2;
          // if(fabs(oldSurfMetric - newSurfMetric) > m_eps) {
          //   cout.precision(18);
          //   cout << "Diff between metrics, easy metric: " << oldSurfMetric << " mod metric: " << newSurfMetric << " diff: " << fabs(oldSurfMetric - newSurfMetric) << endl;
          // }
          m_cells->surfaceMetrics[cellId][isd] = (subMetric1[isd] + subMetric2[isd] + subMetric3[isd] + subMetric4[isd])/F2;
        }


        ///////////////////////////////////////////////////////
        ////////////////////// FACE J /////////////////////////
        //////////////////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          //edge centers
          p1[isd] = (m_coordinates[isd][ipjpk] + m_coordinates[isd][ijpk])/F2;
          p2[isd] = (m_coordinates[isd][ipjpkp] + m_coordinates[isd][ijpkp])/F2;
          p3[isd] = (m_coordinates[isd][ipjpkp] + m_coordinates[isd][ipjpk])/F2;
          p4[isd] = (m_coordinates[isd][ijpkp] + m_coordinates[isd][ijpk])/F2;

          //surface centroid
          center[isd] = (m_coordinates[isd][ijpk] + m_coordinates[isd][ipjpk] + m_coordinates[isd][ipjpkp] + m_coordinates[isd][ijpkp])/F4;
        }

        ////////////////
        // Submetric 1//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = center[isd] - m_coordinates[isd][ijpk];
          diagonal2[isd] = p1[isd] - p4[isd];
        }

        // compute first subMetric
        this->crossProduct( subMetric1, diagonal1, diagonal2 );

        ////////////////
        // Submetric 2//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = p3[isd] - p1[isd];
          diagonal2[isd] = m_coordinates[isd][ipjpk] - center[isd];
        }

        // compute second subMetric
        this->crossProduct( subMetric2, diagonal1, diagonal2 );

        ////////////////
        // Submetric 3//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = m_coordinates[isd][ipjpkp] - center[isd];
          diagonal2[isd] = p3[isd] - p2[isd];
        }

        // compute third subMetric
        this->crossProduct( subMetric3, diagonal1, diagonal2 );

        ////////////////
        // Submetric 4//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = p2[isd] - p4[isd];
          diagonal2[isd] = center[isd] - m_coordinates[isd][ijpkp];
        }

        // compute fourth subMetric
        this->crossProduct( subMetric4, diagonal1, diagonal2 );

        // assemble subMetrics
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          // ZFSFloat oldSurfMetric = m_cells->surfaceMetrics[cellId][3+isd];
          // ZFSFloat newSurfMetric = (subMetric1[isd] + subMetric2[isd] + subMetric3[isd] + subMetric4[isd])/F2;
          // if(fabs(oldSurfMetric - newSurfMetric) > m_eps) {
          //   cout.precision(18);
          //   cout << "Diff between metrics, easy metric: " << oldSurfMetric << " mod metric: " << newSurfMetric << " diff: " << fabs(oldSurfMetric - newSurfMetric) << endl;
          // }
          m_cells->surfaceMetrics[cellId][3 + isd] = (subMetric1[isd] + subMetric2[isd] + subMetric3[isd] + subMetric4[isd])/F2;
        }


        ///////////////////////////////////////////////////////
        ////////////////////// FACE K /////////////////////////
        //////////////////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          //edge centers
          p1[isd] = (m_coordinates[isd][ipjpkp] + m_coordinates[isd][ipjkp])/F2;
          p2[isd] = (m_coordinates[isd][ijpkp] + m_coordinates[isd][ijkp])/F2;
          p3[isd] = (m_coordinates[isd][ipjpkp] + m_coordinates[isd][ijpkp])/F2;
          p4[isd] = (m_coordinates[isd][ipjkp] + m_coordinates[isd][ijkp])/F2;

          //surface centroid
          center[isd] = (m_coordinates[isd][ijkp] + m_coordinates[isd][ipjpkp] + m_coordinates[isd][ipjkp] + m_coordinates[isd][ijpkp])/F4;
        }

        ////////////////
        // Submetric 1//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = p1[isd] - p4[isd];
          diagonal2[isd] = center[isd] - m_coordinates[isd][ipjkp];
        }

        // compute first subMetric
        this->crossProduct( subMetric1, diagonal1, diagonal2 );

        ////////////////
        // Submetric 2//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = m_coordinates[isd][ipjpkp] - center[isd];
          diagonal2[isd] = p3[isd] - p1[isd];
        }

        // compute second subMetric
        this->crossProduct( subMetric2, diagonal1, diagonal2 );

        ////////////////
        // Submetric 3//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = p3[isd] - p2[isd];
          diagonal2[isd] = m_coordinates[isd][ijpkp] - center[isd];
        }

        // compute third subMetric
        this->crossProduct( subMetric3, diagonal1, diagonal2 );

        ////////////////
        // Submetric 4//
        ////////////////
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          diagonal1[isd] = center[isd] - m_coordinates[isd][ijkp];
          diagonal2[isd] = p2[isd] - p4[isd];
        }

        // compute fourth subMetric
        this->crossProduct( subMetric4, diagonal1, diagonal2 );

        // assemble subMetricsp
        for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          // ZFSFloat oldSurfMetric = m_cells->surfaceMetrics[cellId][6+isd];
          // ZFSFloat newSurfMetric = (subMetric1[isd] + subMetric2[isd] + subMetric3[isd] + subMetric4[isd])/F2;
          // if(fabs(oldSurfMetric - newSurfMetric) > m_eps) {
          //   cout.precision(18);
          //   cout << "Diff between metrics, easy metric: " << oldSurfMetric << " mod metric: " << newSurfMetric << " diff: " << fabs(oldSurfMetric - newSurfMetric) << endl;
          // }
          m_cells->surfaceMetrics[cellId][6 + isd] = (subMetric1[isd] + subMetric2[isd] + subMetric3[isd] + subMetric4[isd])/F2;
        }
      }
    }
  }
}



void ZFSStrctrdBlck3D::computeCellMetrics()
{
  TRACE();

  for( ZFSInt k = m_noGhostLayers - 1; k < this->m_nCells[0] - 1; k++ ) {
    for( ZFSInt j = m_noGhostLayers - 1; j < this->m_nCells[1] - 1; j++ ) {
      for( ZFSInt i = m_noGhostLayers -1; i < this->m_nCells[2] - 1; i++ ) {
        ZFSId cellId = this->cellIndex( i, j, k );
        // auxilliary variables
        ZFSFloat   DcoordDxi[ 3 ];
        ZFSFloat  DcoordDeta[ 3 ];
        ZFSFloat DcoordDzeta[ 3 ];

        ZFSFloat metricTmp[ 3 ];

        for( ZFSInt isd = xsd; isd < nDim; isd++) {
          DcoordDxi[isd]   = ( m_cells->coordinates[ isd ][ cellIndex( i + 1, j, k ) ] -
                               m_cells->coordinates[ isd ][ cellIndex( i - 1, j, k ) ] ) * F1B2;

          DcoordDeta[isd]  = ( m_cells->coordinates[ isd ][ cellIndex( i, j + 1, k ) ] -
                               m_cells->coordinates[ isd ][ cellIndex( i, j - 1, k ) ] ) * F1B2;

          DcoordDzeta[isd] = ( m_cells->coordinates[ isd ][ cellIndex( i, j, k + 1 ) ] -
                               m_cells->coordinates[ isd ][ cellIndex( i, j, k - 1 ) ] ) * F1B2;
        }

        //compute metric terms and store them

        //dxi
        this->crossProduct( metricTmp, DcoordDeta, DcoordDzeta );
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          m_cells->cellMetrics[cellId][ xsd * nDim + isd ] = metricTmp[isd];
        }
        //deta
        this->crossProduct( metricTmp, DcoordDzeta, DcoordDxi );
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          m_cells->cellMetrics[cellId][ ysd * nDim + isd ] = metricTmp[isd];
        }
        //dzeta
        this->crossProduct( metricTmp, DcoordDxi, DcoordDeta );
        for( ZFSInt isd = xsd; isd < nDim; isd ++ ) {
          m_cells->cellMetrics[cellId][ zsd * nDim + isd ] = metricTmp[isd];
        }
      }
    }
  }
}


void ZFSStrctrdBlck3D::computeCornerMetrics()
{
  TRACE();
  const ZFSInt noCells = m_noStrctrdCells;

  for( ZFSInt k = m_noGhostLayers - 1; k < this->m_nCells[0] - m_noGhostLayers; k++ ) {
    for( ZFSInt j = m_noGhostLayers - 1; j < this->m_nCells[1] - m_noGhostLayers; j++ ) {
      for( ZFSInt i = m_noGhostLayers -1; i < this->m_nCells[2] - m_noGhostLayers; i++ ) {
        // determine global cell ID
        ZFSId cellId = this->cellIndex(i,j,k); //i + ( k * m_nCells[1] + j ) * m_nCells[2];

        // auxilliary variables
        ZFSFloat   DcoordDxi[3 ];
        ZFSFloat  DcoordDeta[3];
        ZFSFloat DcoordDzeta[3];

        ZFSFloat metricTmp[3];

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          // looks complicated, but what happens is that we always catch the point Id of ipjpkp
          // from the neighboring cell and build the centered difference

          //compute d(x,y,z)/dxi
          DcoordDxi[isd] = F1B2 *
            ( m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k), 1, 1, 1 ) ] -
              m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i-1, j, k), 1, 1, 1 ) ] );

          //compute d(x,y,z)/deta
          DcoordDeta[isd] = F1B2 *
            ( m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k), 1, 1, 1 ) ] -
              m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j-1, k), 1, 1, 1 ) ] );

          //compute d(x,y,z)/dzeta
          DcoordDzeta[isd] = F1B2 *
            ( m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k+1), 1, 1, 1 ) ] -
              m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k-1), 1, 1, 1 ) ] );
        }

        //compute metric terms and store them

        //dxi
        this->crossProduct( metricTmp, DcoordDeta, DcoordDzeta );
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          m_cells->cornerMetrics[(xsd * nDim+isd)*noCells + cellId] = metricTmp[isd];
        }
        //deta
        this->crossProduct( metricTmp, DcoordDzeta, DcoordDxi );
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          m_cells->cornerMetrics[(ysd * nDim+isd)*noCells + cellId] = metricTmp[isd];
        }
        //dzeta
        this->crossProduct( metricTmp, DcoordDxi, DcoordDeta );
        for( ZFSInt isd = xsd; isd < nDim; isd ++ ) {
          m_cells->cornerMetrics[(zsd * nDim+isd)*noCells + cellId] = metricTmp[isd];
        }
      }
    }
  }
}


void ZFSStrctrdBlck3D::computeModCornerMetrics()
{
  TRACE();
  zfs_log << "computing corner metrics ..." << endl;
  const ZFSInt noCells = m_noStrctrdCells;

  for( ZFSInt k = m_noGhostLayers - 1; k < this->m_nCells[0] - m_noGhostLayers; k++ ) {
    for( ZFSInt j = m_noGhostLayers - 1; j < this->m_nCells[1] - m_noGhostLayers; j++ ) {
      for( ZFSInt i = m_noGhostLayers -1; i < this->m_nCells[2] - m_noGhostLayers; i++ ) {
        // determine global cell ID
        ZFSId cellId = this->cellIndex(i,j,k); //i + ( k * m_nCells[1] + j ) * m_nCells[2];
        ZFSFloat metricTmp[3];

        ZFSFloat p1[3];
        ZFSFloat p2[3];
        ZFSFloat p3[3];
        ZFSFloat p4[3];

        ZFSFloat diag1[3];
        ZFSFloat diag2[3];


        ////////////////////////////////
        ////////// DXI /////////////////
        ////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          // looks complicated, but what happens is that we always catch the point Id of ipjpkp
          // from the neighboring cell and build the centered difference

          p1[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k+1), 1, 0, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k+1), 1, 0, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k+1), 1, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k+1), 1, 1, 1 ) ]);


          p2[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 1, 0, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 1, 0, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 1, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 1, 1, 1 ) ]);


          p3[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k+1), 1, 0, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k+1), 1, 0, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k+1), 1, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k+1), 1, 1, 1 ) ]);

          p4[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k), 1, 0, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k), 1, 0, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k), 1, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k), 1, 1, 1 ) ]);

          diag1[isd] = p1[isd] - p2[isd];
          diag2[isd] = p3[isd] - p4[isd];
        }

        this->crossProduct( metricTmp, diag1, diag2 );
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          m_cells->cornerMetrics[(xsd * nDim+isd)*noCells + cellId] = F1B2*metricTmp[isd];
        }


        ////////////////////////////////
        ////////// DETA ////////////////
        ////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          // looks complicated, but what happens is that we always catch the point Id of ipjpkp
          // from the neighboring cell and build the centered difference

          p1[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k+1), 0, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k+1), 1, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k+1), 0, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k+1), 1, 1, 1 ) ]);


          p2[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 0, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 1, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 0, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 1, 1, 1 ) ]);

          p3[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k), 0, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k), 1, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k), 0, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k), 1, 1, 1 ) ]);

          p4[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k+1), 0, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k+1), 1, 1, 0 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k+1), 0, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k+1), 1, 1, 1 ) ]);

          diag1[isd] = p1[isd] - p2[isd];
          diag2[isd] = p3[isd] - p4[isd];
        }

        this->crossProduct( metricTmp, diag1, diag2 );
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          m_cells->cornerMetrics[(ysd * nDim+isd)*noCells + cellId] = F1B2*metricTmp[isd];
        }

        ////////////////////////////////
        ////////// DZETA ///////////////
        ////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {

          p1[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j+1, k), 1, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j+1, k), 0, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j+1, k), 1, 0, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j+1, k), 0, 0, 1 ) ]);


          p2[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 1, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 0, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 1, 0, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j, k), 0, 0, 1 ) ]);

          p3[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k), 1, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k), 0, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k), 1, 0, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i, j+1, k), 0, 0, 1 ) ]);

          p4[isd] = F1B4*(
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k), 1, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k), 0, 1, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k), 1, 0, 1 ) ] +
            m_coordinates[isd][ getPointIdfromPoint( getPointIdFromCell(i+1, j, k), 0, 0, 1 ) ]);

          diag1[isd] = p1[isd] - p2[isd];
          diag2[isd] = p3[isd] - p4[isd];
        }

        this->crossProduct( metricTmp, diag1, diag2 );
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          m_cells->cornerMetrics[(zsd * nDim+isd)*noCells + cellId] = F1B2*metricTmp[isd];
        }
      }
    }
  }

  zfs_log << "computing mod corner metrics ... SUCCESSFUL " << endl;
}


void ZFSStrctrdBlck3D::computeCellJacobian(){
  TRACE();

  for( ZFSInt k = m_noGhostLayers - 1; k < this->m_nCells[0] - m_noGhostLayers; k++ ) {
    for( ZFSInt j = m_noGhostLayers - 1; j < this->m_nCells[1] - m_noGhostLayers; j++ ) {
      for( ZFSInt i = m_noGhostLayers -1; i < this->m_nCells[2] - m_noGhostLayers; i++ ) {
        ZFSId cellId = cellIndex( i, j, k );
        ZFSFloat invJac = m_cells->cellMetrics[ cellId ][ xsd * 3 + xsd ] *
          (
            m_cells->cellMetrics[ cellId ][ ysd * 3 + ysd] *
            m_cells->cellMetrics[ cellId ][ zsd * 3 + zsd] -
            m_cells->cellMetrics[ cellId ][ ysd * 3 + zsd] *
            m_cells->cellMetrics[ cellId ][ zsd * 3 + ysd]
            ) -
          m_cells->cellMetrics[ cellId ][ ysd * 3 + xsd ] *
          (
            m_cells->cellMetrics[ cellId ][ xsd * 3 + ysd] *
            m_cells->cellMetrics[ cellId ][ zsd * 3 + zsd] -
            m_cells->cellMetrics[ cellId ][ xsd * 3 + zsd] *
            m_cells->cellMetrics[ cellId ][ zsd * 3 + ysd]
            ) +
          m_cells->cellMetrics[ cellId ][ zsd * 3 + xsd ] *
          (
            m_cells->cellMetrics[ cellId ][ xsd * 3 + ysd] *
            m_cells->cellMetrics[ cellId ][ ysd * 3 + zsd] -
            m_cells->cellMetrics[ cellId ][ xsd * 3 + zsd] *
            m_cells->cellMetrics[ cellId ][ ysd * 3 + ysd]
            );

        // since metric terms are with omitted jacobian
        // there is factor of J^3; multiplied with J^-1 (invJac) we get J^2
        // --> take square root to get J
        this->m_cells->cellJac[cellId] = sqrt( invJac );
      }
    }
  }

}


/** brief computes cell Jacobian
 *  similar routine as JACGP in TFS
 *  more exact calculation of the Jacobian
 *  than with metrics by using 8 subJacobians
 *  for each cell.
 *  accuracy
 *
 */
void ZFSStrctrdBlck3D::computeModCellJacobian()
{
  ZFSFloat subJ[8] = {F0,F0,F0,F0,F0,F0,F0,F0};

  ZFSFloat** __restrict coords = m_coordinates;


  // auxilliary variables for surface values
  ZFSFloat S1[3];
  ZFSFloat S2[3];
  ZFSFloat S3[3];
  ZFSFloat S1P[3];
  ZFSFloat S2P[3];
  ZFSFloat S3P[3];
  
  // tmp storage for vectors
  ZFSFloat CP1[3];
  ZFSFloat CP2[3];
  ZFSFloat CP3[3];
  ZFSFloat CP4[3];
  ZFSFloat CP5[3];
  ZFSFloat CP6[3];
  ZFSFloat CP7[3];
  ZFSFloat CP8[3];
  ZFSFloat CP9[3];
  ZFSFloat CP10[3];
  ZFSFloat CP11[3];
  ZFSFloat CP12[3];
  // vectors for metrics
  ZFSFloat tmpX1[3];
  ZFSFloat tmpX2[3];
  ZFSFloat tmpX3[3];
  ZFSFloat tmpX4[3];
  ZFSFloat tmpX5[3];
  ZFSFloat tmpX6[3];
  //tmp metric storage dxi
  ZFSFloat DX1[3];
  //tmp metric storage deta
  ZFSFloat DX2[3];
  //tmp metric storage dzeta
  ZFSFloat DX3[3];
  
  for( ZFSInt k = m_noGhostLayers - 1; k < this->m_nCells[0] - m_noGhostLayers; k++ ) {
    for( ZFSInt j = m_noGhostLayers - 1; j < this->m_nCells[1] - m_noGhostLayers; j++ ) {
      for( ZFSInt i = m_noGhostLayers -1; i < this->m_nCells[2] - m_noGhostLayers; i++ ) {
        const ZFSId cellId = cellIndex( i, j, k );

        const ZFSId ijk    = getPointIdFromCell( i, j, k );
        const ZFSId ijpk   = getPointIdfromPoint( ijk, 0, 1, 0 );
        const ZFSId ijkp   = getPointIdfromPoint( ijk, 0, 0, 1 );
        const ZFSId ijpkp  = getPointIdfromPoint( ijk, 0, 1, 1 );
        const ZFSId ipjk   = getPointIdfromPoint( ijk, 1, 0 ,0 );
        const ZFSId ipjpk  = getPointIdfromPoint( ijk, 1, 1, 0 );
        const ZFSId ipjkp  = getPointIdfromPoint( ijk, 1, 0, 1 );
        const ZFSId ipjpkp = getPointIdfromPoint( ijk, 1, 1, 1 );


        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          //averaging of the grid points of surface 1 (j+1/2,k+1/2) around corner point
          S1[isd] = F1B4 * ( coords[isd][ijpkp] + coords[isd][ijk]
                             + coords[isd][ijkp] + coords[isd][ijpk] );
          //averaging of the grid points of surface 2 (i+1/2,k+1/2) around corner point
          S2[isd] = F1B4 * ( coords[isd][ipjk] + coords[isd][ijk]
                             + coords[isd][ipjkp] + coords[isd][ijkp] );
          //averaging of the grid points of surface 3 (i+1/2,j+1/2) around corner point
          S3[isd] = F1B4 * ( coords[isd][ipjk] + coords[isd][ijk]
                             + coords[isd][ipjpk] + coords[isd][ijpk] );
          //averaging of the grid points of surface 1p (j+1/2,k+1/2) around corner point
          S1P[isd] = F1B4 * ( coords[isd][ipjpkp] + coords[isd][ipjk]
                              + coords[isd][ipjkp] + coords[isd][ipjpk] );
          //averaging of the grid points of surface 2p (i+1/2,k+1/2)
          S2P[isd] = F1B4 * ( coords[isd][ipjpk] + coords[isd][ijpk]
                              + coords[isd][ipjpkp] + coords[isd][ijpkp] );
          //averaging of the grid oints of surface 3p (i+1/2,j+1/2)
          S3P[isd] = F1B4 * ( coords[isd][ipjkp] + coords[isd][ijkp]
                              + coords[isd][ipjpkp] + coords[isd][ijpkp] );
        }



        ///////////////////////////////////////////
        ////////// subjacobian 1 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          //averaging corner points
          CP1[isd] = F1B2 * ( coords[isd][ipjk] + coords[isd][ijk] );
          CP2[isd] = F1B2 * ( coords[isd][ijpk] + coords[isd][ijk] );
          CP3[isd] = F1B2 * ( coords[isd][ijkp] + coords[isd][ijk] );

          //setting up vectors for new metric terms
          tmpX1[isd] = ( CP2[isd] - CP3[isd] ) ;
          tmpX2[isd] = ( S1[isd] - coords[isd][ijk] ) ;

          tmpX3[isd] = ( CP3[isd] - CP1[isd] ) ;
          tmpX4[isd] = ( S2[isd] - coords[isd][ijk] ) ;

          tmpX5[isd] = ( S3[isd] - coords[isd][ijk] ) ;
          tmpX6[isd] = ( CP2[isd] - CP1[isd] ) ;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ[0] = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ[0] += ( m_cells->coordinates[isd][cellId] - coords[isd][ijk] )
            * F1B2*( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 2 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP4[isd] = F1B2 * ( coords[isd][ipjk] + coords[isd][ipjpk] );
          CP5[isd] = F1B2 * ( coords[isd][ipjk] + coords[isd][ipjkp] );

          tmpX1[isd] = ( S3[isd] - S2[isd] ) ;
          tmpX2[isd] = ( m_cells->coordinates[isd][cellId] - CP1[isd] ) ;

          tmpX3[isd] = ( S2[isd] - coords[isd][ipjk] ) ;
          tmpX4[isd] = ( CP5[isd] - CP1[isd] ) ;

          tmpX5[isd] = ( CP4[isd] - CP1[isd] ) ;
          tmpX6[isd] = ( S3[isd] - coords[isd][ipjk] ) ;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ[1] = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ[1] += ( S1P[isd] - CP1[isd] ) * F1B2*( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 3 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP6[isd] = F1B2 * ( coords[isd][ipjpk] + coords[isd][ijpk] );
          CP7[isd] = F1B2 * ( coords[isd][ijpkp] + coords[isd][ijpk] );

          tmpX1[isd] = ( coords[isd][ijpk] - S1[isd] ) ;
          tmpX2[isd] = ( CP7[isd] - CP2[isd] ) ;

          tmpX3[isd] = ( S1[isd] - S3[isd] ) ;
          tmpX4[isd] = ( m_cells->coordinates[isd][cellId] - CP2[isd] ) ;

          tmpX5[isd] = ( CP6[isd] - CP2[isd] ) ;
          tmpX6[isd] = ( coords[isd][ijpk] - S3[isd] ) ;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ[2] = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ[2] += ( S2P[isd] - CP2[isd] ) * F1B2*( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 4 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP8[isd] = F1B2 * ( coords[isd][ipjpkp] + coords[isd][ipjpk] );

          tmpX1[isd] = ( CP6[isd] - m_cells->coordinates[isd][cellId] ) ;
          tmpX2[isd] = ( S2P[isd] - S3[isd] ) ;

          tmpX3[isd] = ( m_cells->coordinates[isd][cellId] - CP4[isd] ) ;
          tmpX4[isd] = ( S1P[isd] - S3[isd] ) ;

          tmpX5[isd] = ( coords[isd][ipjpk] - S3[isd] ) ;
          tmpX6[isd] = ( CP6[isd] - CP4[isd] ) ;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ[3] = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ[3] += ( CP8[isd] - S3[isd] ) * F1B2*( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 5 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP9[isd]  = F1B2 * ( coords[isd][ipjkp] + coords[isd][ijkp] );
          CP10[isd] = F1B2 * ( coords[isd][ijpkp] + coords[isd][ijkp] );

          tmpX1[isd] = ( S1[isd] - coords[isd][ijkp] ) ;
          tmpX2[isd] = ( CP10[isd] - CP3[isd] ) ;

          tmpX3[isd] = ( coords[isd][ijkp] - S2[isd] ) ;
          tmpX4[isd] = ( CP9[isd] - CP3[isd] ) ;

          tmpX5[isd] = ( m_cells->coordinates[isd][cellId] - CP3[isd] ) ;
          tmpX6[isd] = ( S1[isd] - S2[isd] ) ;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ[4] = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ[4] += ( S3P[isd] - CP3[isd] ) * F1B2*( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 6 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP11[isd] = F1B2 * ( coords[isd][ipjkp] + coords[isd][ipjpkp] );

          tmpX1[isd] = ( m_cells->coordinates[isd][cellId] - CP9[isd] ) ;
          tmpX2[isd] = ( S3P[isd] - S2[isd] ) ;

          tmpX3[isd] = ( CP9[isd] - CP5[isd] ) ;
          tmpX4[isd] = ( coords[isd][ipjkp] - S2[isd] ) ;

          tmpX5[isd] = ( S1P[isd] - S2[isd] ) ;
          tmpX6[isd] = ( m_cells->coordinates[isd][cellId] - CP5[isd] ) ;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ[5] = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ[5] += ( CP11[isd] - S2[isd] ) * F1B2*( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 7 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP12[isd] = F1B2 * ( coords[isd][ipjpkp] + coords[isd][ijpkp] );

          tmpX1[isd] = ( CP7[isd] - CP10[isd] ) ;
          tmpX2[isd] = ( coords[isd][ijpkp] - S1[isd] ) ;

          tmpX3[isd] = ( CP10[isd] - m_cells->coordinates[isd][cellId] ) ;
          tmpX4[isd] = ( S3P[isd] - S1[isd] ) ;

          tmpX5[isd] = ( S2P[isd] - S1[isd] ) ;
          tmpX6[isd] = ( CP7[isd] - m_cells->coordinates[isd][cellId] ) ;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ[6] = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ[6] += ( CP12[isd] - S1[isd] ) * F1B2*( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 8 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          tmpX1[isd] = ( S2P[isd] - S3P[isd] ) ;
          tmpX2[isd] = ( CP12[isd] - m_cells->coordinates[isd][cellId] ) ;

          tmpX3[isd] = ( S3P[isd] - S1P[isd] ) ;
          tmpX4[isd] = ( CP11[isd] - m_cells->coordinates[isd][cellId] ) ;

          tmpX5[isd] = ( CP8[isd] - m_cells->coordinates[isd][cellId] ) ;
          tmpX6[isd] = ( S2P[isd] - S1P[isd] ) ;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ[7] = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ[7] += ( coords[isd][ipjpkp] - m_cells->coordinates[isd][cellId] )
            * F1B2*( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        //////////////////////////////////////////////
        ///// assemble subjacobians //////////////////
        //////////////////////////////////////////////
        m_cells->cellJac[cellId] = F0;

        for( ZFSInt jacId = 0; jacId < 8; jacId ++) {
          m_cells->cellJac[cellId] += subJ[jacId];
        }

        m_cells->cellJac[cellId] = F1B3 * fabs(m_cells->cellJac[cellId]);
      }
    }
  }
}



// jacobian for the viscous fluxes
void ZFSStrctrdBlck3D::computeCornerJacobian(){
  TRACE();
  //jacobian in the physical space is the inverse of the jacobian in computational space
  const ZFSInt noCells = m_noStrctrdCells;

  for( ZFSInt k = m_noGhostLayers - 1; k < this->m_nCells[0] - m_noGhostLayers; k++ ) {
    for( ZFSInt j = m_noGhostLayers - 1; j < this->m_nCells[1] - m_noGhostLayers; j++ ) {
      for( ZFSInt i = m_noGhostLayers -1; i < this->m_nCells[2] - m_noGhostLayers; i++ ) {
        const ZFSId cellId = cellIndex( i, j, k );

        const ZFSFloat invJac = m_cells->cornerMetrics[ (xsd * nDim+xsd)*noCells + cellId ] *
          (
            m_cells->cornerMetrics[ (ysd * nDim+ysd)*noCells + cellId] *
            m_cells->cornerMetrics[ (zsd * nDim+zsd)*noCells + cellId] -
            m_cells->cornerMetrics[ (ysd * nDim+zsd)*noCells + cellId] *
            m_cells->cornerMetrics[ (zsd * nDim+ysd)*noCells + cellId]
            ) -
          m_cells->cornerMetrics[ (ysd * nDim+xsd)*noCells + cellId ] *
          (
            m_cells->cornerMetrics[ (xsd * nDim+ysd)*noCells + cellId] *
            m_cells->cornerMetrics[ (zsd * nDim+zsd)*noCells + cellId] -
            m_cells->cornerMetrics[ (xsd * nDim+zsd)*noCells + cellId] *
            m_cells->cornerMetrics[ (zsd * nDim+ysd)*noCells + cellId]
            ) +
          m_cells->cornerMetrics[ (zsd * nDim+xsd)*noCells + cellId ] *
          (
            m_cells->cornerMetrics[ (xsd * nDim+ysd)*noCells + cellId] *
            m_cells->cornerMetrics[ (ysd * nDim+zsd)*noCells + cellId] -
            m_cells->cornerMetrics[ (xsd * nDim+zsd)*noCells + cellId] *
            m_cells->cornerMetrics[ (ysd * nDim+ysd)*noCells + cellId]
            );

        // since metric terms are with omitted jacobian
        // there is factor of J^3; multiplied with J^-1 (invJac) we get J^2
        // --> take square root to get J
        this->m_cells->cornerJac[cellId] = sqrt( invJac );
      }
    }
  }
}


// same as TFS routine JACGP, more exact calculation of control volume
// for corner point than with corner metrics
void ZFSStrctrdBlck3D::computeModCornerJacobian()
{  
  TRACE();
  ZFSFloatScratchSpace subJ(m_noStrctrdCells, 8, __CALLING_FUNCTION__, "subJ");
  ZFSFloatScratchSpace subJtmp(m_noStrctrdCells, 8, __CALLING_FUNCTION__, "subJtmp");

  ZFSFloat** __restrict coords = m_coordinates;

  // tmp storage for vectors
  ZFSFloat CP1[3];
  ZFSFloat CP2[3];
  ZFSFloat CP3[3];
  ZFSFloat CP4[3];
  ZFSFloat CP5[3];
  ZFSFloat CP6[3];
  ZFSFloat CP7[3];
  ZFSFloat CP8[3];
  ZFSFloat CP9[3];
  ZFSFloat CP10[3];
  ZFSFloat CP11[3];
  ZFSFloat CP12[3];
  // vectors for metrics
  ZFSFloat tmpX1[3];
  ZFSFloat tmpX2[3];
  ZFSFloat tmpX3[3];
  ZFSFloat tmpX4[3];
  ZFSFloat tmpX5[3];
  ZFSFloat tmpX6[3];
  //tmp metric storage dxi
  ZFSFloat DX1[3];
  //tmp metric storage deta
  ZFSFloat DX2[3];
  //tmp metric storage dzeta
  ZFSFloat DX3[3];

  // auxilliary variables for surface values
  ZFSFloat S1[3];
  ZFSFloat S2[3];
  ZFSFloat S3[3];
  ZFSFloat S1P[3];
  ZFSFloat S2P[3];
  ZFSFloat S3P[3];

  subJ.fill(1234.56);
  subJtmp.fill(5678.9);

  const ZFSFloat sqrttwo = sqrt(2);
  const ZFSFloat fsqrttwo = F1/sqrttwo;

  for( ZFSInt k = m_noGhostLayers - 2; k < this->m_nCells[0] - m_noGhostLayers; k++ ) {
    for( ZFSInt j = m_noGhostLayers - 2; j < this->m_nCells[1] - m_noGhostLayers; j++ ) {
      for( ZFSInt i = m_noGhostLayers -2; i < this->m_nCells[2] - m_noGhostLayers; i++ ) {
        const ZFSId cellId = cellIndex( i, j, k );
        const ZFSId centCellId = cellIndex( i+1, j+1, k+1 );
        const ZFSId tmpId = getPointIdFromCell( i, j, k );

        const ZFSId ijk    = getPointIdfromPoint( tmpId, 1, 1, 1 );
        const ZFSId ijpk   = getPointIdfromPoint( ijk, 0, 1, 0 );
        const ZFSId ijkp   = getPointIdfromPoint( ijk, 0, 0, 1 );
        const ZFSId ijpkp  = getPointIdfromPoint( ijk, 0, 1, 1 );
        const ZFSId ipjk   = getPointIdfromPoint( ijk, 1, 0 ,0 );
        const ZFSId ipjpk  = getPointIdfromPoint( ijk, 1, 1, 0 );
        const ZFSId ipjkp  = getPointIdfromPoint( ijk, 1, 0, 1 );
        const ZFSId ipjpkp = getPointIdfromPoint( ijk, 1, 1, 1 );

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          //averaging of the grid points of surface 1 (j+1/2,k+1/2) around corner point
          S1[isd] = F1B4 * ( coords[isd][ijpkp] + coords[isd][ijk]
                             + coords[isd][ijkp] + coords[isd][ijpk] );
          //averaging of the grid points of surface 2 (i+1/2,k+1/2) around corner point
          S2[isd] = F1B4 * ( coords[isd][ipjk] + coords[isd][ijk]
                             + coords[isd][ipjkp] + coords[isd][ijkp] );
          //averaging of the grid points of surface 3 (i+1/2,j+1/2) around corner point
          S3[isd] = F1B4 * ( coords[isd][ipjk] + coords[isd][ijk]
                             + coords[isd][ipjpk] + coords[isd][ijpk] );
          //averaging of the grid points of surface 1p (j+1/2,k+1/2) around corner point
          S1P[isd] = F1B4 * ( coords[isd][ipjpkp] + coords[isd][ipjk]
                              + coords[isd][ipjkp] + coords[isd][ipjpk] );
          //averaging of the grid points of surface 2p (i+1/2,k+1/2)
          S2P[isd] = F1B4 * ( coords[isd][ipjpk] + coords[isd][ijpk]
                              + coords[isd][ipjpkp] + coords[isd][ijpkp] );
          //averaging of the grid oints of surface 3p (i+1/2,j+1/2)
          S3P[isd] = F1B4 * ( coords[isd][ipjkp] + coords[isd][ijkp]
                              + coords[isd][ipjpkp] + coords[isd][ijpkp] );
        }


        ///////////////////////////////////////////
        ////////// subjacobian 1 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          //averaging corner points
          CP1[isd] = F1B2 * ( coords[isd][ipjk] + coords[isd][ijk] );
          CP2[isd] = F1B2 * ( coords[isd][ijpk] + coords[isd][ijk] );
          CP3[isd] = F1B2 * ( coords[isd][ijkp] + coords[isd][ijk] );

          //setting up vectors for new metric terms
          tmpX1[isd] = ( CP2[isd] - CP3[isd] ) *fsqrttwo;
          tmpX2[isd] = ( S1[isd] - coords[isd][ijk] ) *fsqrttwo;

          tmpX3[isd] = ( CP3[isd] - CP1[isd] ) *fsqrttwo;
          tmpX4[isd] = ( S2[isd] - coords[isd][ijk] ) *fsqrttwo;

          tmpX5[isd] = ( S3[isd] - coords[isd][ijk] ) *fsqrttwo;
          tmpX6[isd] = ( CP2[isd] - CP1[isd] ) *fsqrttwo;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ(cellId,0) = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ(cellId,0) += ( m_cells->coordinates[isd][centCellId] - coords[isd][ijk] )
            * ( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 2 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP4[isd] = F1B2 * ( coords[isd][ipjk] + coords[isd][ipjpk] );
          CP5[isd] = F1B2 * ( coords[isd][ipjk] + coords[isd][ipjkp] );

          tmpX1[isd] = ( S3[isd] - S2[isd] ) *fsqrttwo;
          tmpX2[isd] = ( m_cells->coordinates[isd][centCellId] - CP1[isd] ) *fsqrttwo;

          tmpX3[isd] = ( S2[isd] - coords[isd][ipjk] ) *fsqrttwo;
          tmpX4[isd] = ( CP5[isd] - CP1[isd] ) *fsqrttwo;

          tmpX5[isd] = ( CP4[isd] - CP1[isd] ) *fsqrttwo;
          tmpX6[isd] = ( S3[isd] - coords[isd][ipjk] ) *fsqrttwo;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ(cellId,1) = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ(cellId,1) += ( S1P[isd] - CP1[isd] ) * ( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 3 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP6[isd] = F1B2 * ( coords[isd][ipjpk] + coords[isd][ijpk] );
          CP7[isd] = F1B2 * ( coords[isd][ijpkp] + coords[isd][ijpk] );

          tmpX1[isd] = ( coords[isd][ijpk] - S1[isd] ) *fsqrttwo;
          tmpX2[isd] = ( CP7[isd] - CP2[isd] ) *fsqrttwo;

          tmpX3[isd] = ( S1[isd] - S3[isd] ) *fsqrttwo;
          tmpX4[isd] = ( m_cells->coordinates[isd][centCellId] - CP2[isd] ) *fsqrttwo;

          tmpX5[isd] = ( CP6[isd] - CP2[isd] ) *fsqrttwo;
          tmpX6[isd] = ( coords[isd][ijpk] - S3[isd] ) *fsqrttwo;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ(cellId,2) = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ(cellId,2) += ( S2P[isd] - CP2[isd] ) * ( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 4 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP8[isd] = F1B2 * ( coords[isd][ipjpkp] + coords[isd][ipjpk] );

          tmpX1[isd] = ( CP6[isd] - m_cells->coordinates[isd][centCellId] ) *fsqrttwo;
          tmpX2[isd] = ( S2P[isd] - S3[isd] ) *fsqrttwo;

          tmpX3[isd] = ( m_cells->coordinates[isd][centCellId] - CP4[isd] ) *fsqrttwo;
          tmpX4[isd] = ( S1P[isd] - S3[isd] ) *fsqrttwo;

          tmpX5[isd] = ( coords[isd][ipjpk] - S3[isd] ) *fsqrttwo;
          tmpX6[isd] = ( CP6[isd] - CP4[isd] ) *fsqrttwo;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ(cellId,3) = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ(cellId,3) += ( CP8[isd] - S3[isd] ) * ( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 5 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP9[isd]  = F1B2 * ( coords[isd][ipjkp] + coords[isd][ijkp] );
          CP10[isd] = F1B2 * ( coords[isd][ijpkp] + coords[isd][ijkp] );

          tmpX1[isd] = ( S1[isd] - coords[isd][ijkp] ) *fsqrttwo;
          tmpX2[isd] = ( CP10[isd] - CP3[isd] ) *fsqrttwo;

          tmpX3[isd] = ( coords[isd][ijkp] - S2[isd] ) *fsqrttwo;
          tmpX4[isd] = ( CP9[isd] - CP3[isd] ) *fsqrttwo;

          tmpX5[isd] = ( m_cells->coordinates[isd][centCellId] - CP3[isd] ) *fsqrttwo;
          tmpX6[isd] = ( S1[isd] - S2[isd] ) *fsqrttwo;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ(cellId,4) = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ(cellId,4) += ( S3P[isd] - CP3[isd] ) * ( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 6 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP11[isd] = F1B2 * ( coords[isd][ipjkp] + coords[isd][ipjpkp] );

          tmpX1[isd] = ( m_cells->coordinates[isd][centCellId] - CP9[isd] ) *fsqrttwo;
          tmpX2[isd] = ( S3P[isd] - S2[isd] ) *fsqrttwo;

          tmpX3[isd] = ( CP9[isd] - CP5[isd] ) *fsqrttwo;
          tmpX4[isd] = ( coords[isd][ipjkp] - S2[isd] ) *fsqrttwo;

          tmpX5[isd] = ( S1P[isd] - S2[isd] ) *fsqrttwo;
          tmpX6[isd] = ( m_cells->coordinates[isd][centCellId] - CP5[isd] ) *fsqrttwo;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ(cellId,5) = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ(cellId,5) += ( CP11[isd] - S2[isd] ) * ( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 7 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          CP12[isd] = F1B2 * ( coords[isd][ipjpkp] + coords[isd][ijpkp] );

          tmpX1[isd] = ( CP7[isd] - CP10[isd] ) *fsqrttwo;
          tmpX2[isd] = ( coords[isd][ijpkp] - S1[isd] ) *fsqrttwo;

          tmpX3[isd] = ( CP10[isd] - m_cells->coordinates[isd][centCellId] ) *fsqrttwo;
          tmpX4[isd] = ( S3P[isd] - S1[isd] ) *fsqrttwo;

          tmpX5[isd] = ( S2P[isd] - S1[isd] ) *fsqrttwo;
          tmpX6[isd] = ( CP7[isd] - m_cells->coordinates[isd][centCellId] ) *fsqrttwo;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ(cellId,6) = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ(cellId,6) += ( CP12[isd] - S1[isd] ) * ( DX1[isd] + DX2[isd] + DX3[isd] );
        }

        ///////////////////////////////////////////
        ////////// subjacobian 8 //////////////////
        ///////////////////////////////////////////

        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          tmpX1[isd] = ( S2P[isd] - S3P[isd] ) *fsqrttwo;
          tmpX2[isd] = ( CP12[isd] - m_cells->coordinates[isd][centCellId] ) *fsqrttwo;

          tmpX3[isd] = ( S3P[isd] - S1P[isd] ) *fsqrttwo;
          tmpX4[isd] = ( CP11[isd] - m_cells->coordinates[isd][centCellId] ) *fsqrttwo;

          tmpX5[isd] = ( CP8[isd] - m_cells->coordinates[isd][centCellId] ) *fsqrttwo;
          tmpX6[isd] = ( S2P[isd] - S1P[isd] ) *fsqrttwo;
        }

        this->crossProduct( DX1, tmpX1, tmpX2 );
        this->crossProduct( DX2, tmpX3, tmpX4 );
        this->crossProduct( DX3, tmpX5, tmpX6 );

        subJ(cellId,7) = F0;
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          subJ(cellId,7) += ( coords[isd][ipjpkp] - m_cells->coordinates[isd][centCellId] )
            * ( DX1[isd] + DX2[isd] + DX3[isd] );
        }
      }
    }
  }

//////////////////////////////////////////////
///// assemble subjacobians //////////////////
//////////////////////////////////////////////

//copy into dummy array
  for( ZFSInt i = 0; i < m_noStrctrdCells; i++ ) {
    for( ZFSInt j = 0; j < 8; j++ ) {
      subJtmp(i,j) = subJ(i,j);
    }
  }

//shift subjacobians
  for( ZFSInt k = m_noGhostLayers-1; k < this->m_nCells[0] - m_noGhostLayers; k++ ) {
    for( ZFSInt j = m_noGhostLayers-1; j < this->m_nCells[1] - m_noGhostLayers; j++ ) {
      for( ZFSInt i = m_noGhostLayers-1; i < this->m_nCells[2] - m_noGhostLayers; i++ ) {
        const ZFSId cellId = cellIndex( i, j, k );

        subJ(cellId,0) = subJ( cellIndex(i-1, j-1, k-1) , 7 );
        subJ(cellId,1) = subJ( cellIndex(i  , j-1, k-1) , 6 );
        subJ(cellId,2) = subJ( cellIndex(i-1, j  , k-1) , 5 );
        subJ(cellId,3) = subJ( cellIndex(i  , j  , k-1) , 4 );
      }
    }
  }

  for( ZFSInt k = m_noGhostLayers-1; k < this->m_nCells[0] - m_noGhostLayers; k++ ) {
    for( ZFSInt j = m_noGhostLayers-1; j < this->m_nCells[1] - m_noGhostLayers; j++ ) {
      for( ZFSInt i = m_noGhostLayers-1; i < this->m_nCells[2] - m_noGhostLayers; i++ ) {
        const ZFSId cellId = cellIndex( i, j, k );
        subJ(cellId,4) = subJtmp( cellIndex(i-1, j-1, k) , 3 );
        subJ(cellId,5) = subJtmp( cellIndex(i  , j-1, k) , 2 );
        subJ(cellId,6) = subJtmp( cellIndex(i-1, j  , k) , 1 );
        subJ(cellId,7) = subJtmp( cellId , 0 );
      }
    }
  }

//finally jacobian at corner point!
  for( ZFSInt k = m_noGhostLayers-1; k < this->m_nCells[0] - m_noGhostLayers; k++ ) {
    for( ZFSInt j = m_noGhostLayers-1; j < this->m_nCells[1] - m_noGhostLayers; j++ ) {
      for( ZFSInt i = m_noGhostLayers-1; i < this->m_nCells[2] - m_noGhostLayers; i++ ) {
        const ZFSId cellId = cellIndex( i, j, k );

        //ZFSFloat oldJacobian = m_cells->cornerJac[cellId];
        m_cells->cornerJac[cellId] = F0;
        for( ZFSInt jacId = 0; jacId < 8; jacId ++) {
          m_cells->cornerJac[cellId] += subJ(cellId,jacId);
        }

        m_cells->cornerJac[cellId] = F1B3 * m_cells->cornerJac[cellId];
      }
    }
  }
}

/** \brief help function for computeDxt
 */
inline ZFSId ZFSStrctrdBlck3D::surfId(ZFSId point, ZFSId isd, ZFSId dim)
{
  return point + (isd + 3*dim)*9;
}

/** \brief compute volume fluxes
 *  same routine as dxtnc in TFS
 *  splits the xi, eta and zeta cell surface into 4 subparts
 *  by computing 9 points on the surface and then computing
 *  the 4 subjacobians on this surface. In general, the procedure
 *  is exactly the same as for one big surface, but with higher
 *  accuracy
 * \author Marian Albers, 2015
 */
void ZFSStrctrdBlck3D::computeDxt()
{
  TRACE();
  const ZFSFloat frk = F1/(m_timeStep*m_RKalpha[m_RKStep]);
  ZFSFloat surfCoordNew[9*3*3];
  ZFSFloat surfCoordOld[9*3*3];
  const ZFSInt IJK[nDim] = {m_nCells[2], m_nCells[1], m_nCells[0]};

  const ZFSFloat *const *const RESTRICT coords = m_coordinates;
  const ZFSFloat *const *const RESTRICT oldCoords = m_mgOldCoordinates;


  for( ZFSInt k = 0; k < IJK[2]; ++k ) {
    for( ZFSInt j = 0; j < IJK[1]; ++j ) {
      for( ZFSInt i = 0; i < IJK[0]; ++i ) {
        // determine global cell ID
        const ZFSId cellId = cellIndex(i,j,k);

        // determine global point ID for local cell IDs
        const ZFSId ijk    = getPointIdFromCell( i, j, k );
        const ZFSId ipjk   = getPointIdfromPoint( ijk, 1, 0, 0 );
        const ZFSId ipjpk  = getPointIdfromPoint( ijk, 1, 1, 0 );
        const ZFSId ipjkp  = getPointIdfromPoint( ijk, 1, 0, 1 );
        const ZFSId ipjpkp = getPointIdfromPoint( ijk, 1, 1, 1 );
        const ZFSId ijpk   = getPointIdfromPoint( ijk, 0, 1, 0 );
        const ZFSId ijpkp  = getPointIdfromPoint( ijk, 0, 1, 1 );
        const ZFSId ijkp   = getPointIdfromPoint( ijk, 0, 0, 1 );

        //compute values for the 9 points of the surface
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          //coordinates at corner point (i+1/2,j-1/2,k-1/2)
          surfCoordNew[surfId(0,isd,0)] = coords[isd][ipjk];
          surfCoordOld[surfId(0,isd,0)] = oldCoords[isd][ipjk];

          //coordinates at surface edge center (i+1/2,j,k-1/2)
          surfCoordNew[surfId(1,isd,0)] = F1B2 * (coords[isd][ipjpk] + coords[isd][ipjk]);
          surfCoordOld[surfId(1,isd,0)] = F1B2 * (oldCoords[isd][ipjpk] + oldCoords[isd][ipjk]);

          //coordinates at corner point (i+1/2,j+1/2,k-1/2)
          surfCoordNew[surfId(2,isd,0)] = coords[isd][ipjpk];
          surfCoordOld[surfId(2,isd,0)] = oldCoords[isd][ipjpk];

          //coordinates at surface edge center (i+1/2,j,k-1/2)
          surfCoordNew[surfId(3,isd,0)] = F1B2 * (coords[isd][ipjkp] + coords[isd][ipjk]);
          surfCoordOld[surfId(3,isd,0)] = F1B2 * (oldCoords[isd][ipjkp] + oldCoords[isd][ipjk]);

          //coordinates at surface centroid (i+1/2,j,k)
          surfCoordNew[surfId(4,isd,0)] = F1B4 * (coords[isd][ipjk]        + coords[isd][ipjpk]        + coords[isd][ipjkp]        + coords[isd][ipjpkp]);
          surfCoordOld[surfId(4,isd,0)] = F1B4 * (oldCoords[isd][ipjk]   + oldCoords[isd][ipjpk]   + oldCoords[isd][ipjkp]   + oldCoords[isd][ipjpkp]);

          //coordinates at surface edge center (i+1/2,j+1/2,k)
          surfCoordNew[surfId(5,isd,0)] = F1B2 * (coords[isd][ipjpkp] + coords[isd][ipjpk]);
          surfCoordOld[surfId(5,isd,0)] = F1B2 * (oldCoords[isd][ipjpkp] + oldCoords[isd][ipjpk]);

          //coordinates at corner point (i+1/2,j-1/2,k+1/2)
          surfCoordNew[surfId(6,isd,0)] = coords[isd][ipjkp];
          surfCoordOld[surfId(6,isd,0)] = oldCoords[isd][ipjkp];

          //coordinates at surface edge center (i+1/2,j,k+1/2)
          surfCoordNew[surfId(7,isd,0)] = F1B2 * (coords[isd][ipjpkp] + coords[isd][ipjkp]);
          surfCoordOld[surfId(7,isd,0)] = F1B2 * (oldCoords[isd][ipjpkp] + oldCoords[isd][ipjkp]);

          //coordinates at corner point (i+1/2,j+1/2,k+1/2)
          surfCoordNew[surfId(8,isd,0)] = coords[isd][ipjpkp];
          surfCoordOld[surfId(8,isd,0)] = oldCoords[isd][ipjpkp];
        }

        //compute values for the 9 points of the surface
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          //coordinates at corner point (i-1/2,j+1/2,k-1/2)
          surfCoordNew[surfId(0,isd,1)] = coords[isd][ijpk];
          surfCoordOld[surfId(0,isd,1)] = oldCoords[isd][ijpk];

          //coordinates at surface edge center (i-1/2,j+1/2,k)
          surfCoordNew[surfId(1,isd,1)] = F1B2 * (coords[isd][ijpkp] + coords[isd][ijpk]);
          surfCoordOld[surfId(1,isd,1)] = F1B2 * (oldCoords[isd][ijpkp] + oldCoords[isd][ijpk]);

          //coordinates at corner point (i-1/2,j+1/2,k+1/2)
          surfCoordNew[surfId(2,isd,1)] = coords[isd][ijpkp];
          surfCoordOld[surfId(2,isd,1)] = oldCoords[isd][ijpkp];

          //coordinates at surface edge center (i,j+1/2,k-1/2)
          surfCoordNew[surfId(3,isd,1)] = F1B2 * (coords[isd][ipjpk] + coords[isd][ijpk]);
          surfCoordOld[surfId(3,isd,1)] = F1B2 * (oldCoords[isd][ipjpk] + oldCoords[isd][ijpk]);

          //coordinates at surface centroid (i,j+1/2,k)
          surfCoordNew[surfId(4,isd,1)] = F1B4 * (coords[isd][ijpk] + coords[isd][ijpkp] + coords[isd][ipjpk] + coords[isd][ipjpkp]);
          surfCoordOld[surfId(4,isd,1)] = F1B4 * (oldCoords[isd][ijpk]   + oldCoords[isd][ijpkp]   + oldCoords[isd][ipjpk]   + oldCoords[isd][ipjpkp]);

          //coordinates at surface edge center (i,j+1/2,k+1/2)
          surfCoordNew[surfId(5,isd,1)] = F1B2 * (coords[isd][ipjpkp] + coords[isd][ijpkp]);
          surfCoordOld[surfId(5,isd,1)] = F1B2 * (oldCoords[isd][ipjpkp] + oldCoords[isd][ijpkp]);

          //coordinates at corner point (i+1/2,j+1/2,k-1/2)
          surfCoordNew[surfId(6,isd,1)] = coords[isd][ipjpk];
          surfCoordOld[surfId(6,isd,1)] = oldCoords[isd][ipjpk];

          //coordinates at surface edge center (i+1/2,j+1/2,k)
          surfCoordNew[surfId(7,isd,1)] = F1B2 * (coords[isd][ipjpkp] + coords[isd][ipjpk]);
          surfCoordOld[surfId(7,isd,1)] = F1B2 * (oldCoords[isd][ipjpkp] + oldCoords[isd][ipjpk]);

          //coordinates at corner point (i+1/2,j+1/2,k)
          surfCoordNew[surfId(8,isd,1)] = coords[isd][ipjpkp];
          surfCoordOld[surfId(8,isd,1)] = oldCoords[isd][ipjpkp];
        }

        //compute values for the 9 points of the surface
        for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
          //coordinates at corner point (i-1/2,j-1/2,k+1/2)
          surfCoordNew[surfId(0,isd,2)] = coords[isd][ijkp];
          surfCoordOld[surfId(0,isd,2)] = oldCoords[isd][ijkp];

          //coordinates at surface edge center (i,j-1/2,k+1/2)
          surfCoordNew[surfId(1,isd,2)] = F1B2 * (coords[isd][ipjkp] + coords[isd][ijkp]);
          surfCoordOld[surfId(1,isd,2)] = F1B2 * (oldCoords[isd][ipjkp] + oldCoords[isd][ijkp]);

          //coordinates at corner point (i+1/2,j-1/2,k+1/2)
          surfCoordNew[surfId(2,isd,2)] = coords[isd][ipjkp];
          surfCoordOld[surfId(2,isd,2)] = oldCoords[isd][ipjkp];

          //coordinates at surface edge center (i-1/2,j,k+1/2)
          surfCoordNew[surfId(3,isd,2)] = F1B2 * (coords[isd][ijpkp] + coords[isd][ijkp]);
          surfCoordOld[surfId(3,isd,2)] = F1B2 * (oldCoords[isd][ijpkp] + oldCoords[isd][ijkp]);

          //coordinates at surface centroid (i,j,k+1/2)
          surfCoordNew[surfId(4,isd,2)] = F1B4 * (coords[isd][ijkp] + coords[isd][ijpkp] + coords[isd][ipjkp] + coords[isd][ipjpkp]);
          surfCoordOld[surfId(4,isd,2)] = F1B4 * (oldCoords[isd][ijkp]   + oldCoords[isd][ijpkp]   + oldCoords[isd][ipjkp]   + oldCoords[isd][ipjpkp]);

          //coordinates at surface edge center (i+1/2,j,k+1/2)
          surfCoordNew[surfId(5,isd,2)] = F1B2 * (coords[isd][ipjpkp] + coords[isd][ipjkp]);
          surfCoordOld[surfId(5,isd,2)] = F1B2 * (oldCoords[isd][ipjpkp] + oldCoords[isd][ipjkp]);

          //coordinates at corner point (i-1/2,j+1/2,k+1/2)
          surfCoordNew[surfId(6,isd,2)] = coords[isd][ijpkp];
          surfCoordOld[surfId(6,isd,2)] = oldCoords[isd][ijpkp];

          //coordinates at surface edge center (i,j+1/2,k+1/2)
          surfCoordNew[surfId(7,isd,2)] = F1B2 * (coords[isd][ipjpkp] + coords[isd][ijpkp]);
          surfCoordOld[surfId(7,isd,2)] = F1B2 * (oldCoords[isd][ipjpkp] + oldCoords[isd][ijpkp]);

          //coordinates at corner point (i+1/2,j+1/2,k+1/2)
          surfCoordNew[surfId(8,isd,2)] = coords[isd][ipjpkp];
          surfCoordOld[surfId(8,isd,2)] = oldCoords[isd][ipjpkp];
        }

        for( ZFSInt dim = 0; dim < nDim; ++dim) {
          //reset the volume flux
          m_cells->dxt[dim][cellId] = F0;

          ZFSFloat subJacobian1, subJacobian2, subJacobian3, subJacobian4;
          ZFSFloat diag1[3];
          ZFSFloat diag2[3];
          ZFSFloat oldNormal[3];
          ZFSFloat newNormal1[3];
          ZFSFloat newNormal2[3];

          //
          //+1/2    X-----x-----X
          //        |     |/////|
          //        |     |/////|  <----  1 sub-Jacobian
          //0       x     o-----x
          //        |           |
          //        |           |
          //-1/2    X-----x-----X
          //
          //      -1/2    0   +1/2
          //
          //
          // The big X are points on the grid, no calculation needed
          // The small x are points on the edge, halfway between two grid points
          // The o is at the surface centroid, take the average from all four surrounding grid points
          // The picture also shows one sub-Jacobian


          //compute 4 subjacobians for 4 subsurfaces of the surface

          //////////////////////////////
          /////// SUB JACOBIAN 1 ///////
          //////////////////////////////

          //(i+1/2,j-1/4,k-1/4)

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordOld[surfId(4,isd,dim)] - surfCoordOld[surfId(0,isd,dim)]);
            diag2[isd] = (surfCoordOld[surfId(3,isd,dim)] - surfCoordOld[surfId(1,isd,dim)]);
          }

          this->crossProduct(oldNormal, diag1, diag2);

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordNew[surfId(3,isd,dim)] - surfCoordOld[surfId(0,isd,dim)]);
            diag2[isd] = (surfCoordNew[surfId(0,isd,dim)] - surfCoordOld[surfId(3,isd,dim)]);
          }


          this->crossProduct(newNormal1, diag1, diag2);

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordNew[surfId(0,isd,dim)] - surfCoordOld[surfId(1,isd,dim)]);
            diag2[isd] = (surfCoordNew[surfId(1,isd,dim)] - surfCoordOld[surfId(0,isd,dim)]);
          }

          this->crossProduct(newNormal2, diag1, diag2);

          subJacobian1 = F1B3 * (  (surfCoordNew[surfId(4,0,dim)] - surfCoordOld[surfId(0,0,dim)])*(oldNormal[0] + newNormal1[0] + newNormal2[0])*F1B2
                                   + (surfCoordNew[surfId(4,1,dim)] - surfCoordOld[surfId(0,1,dim)])*(oldNormal[1] + newNormal1[1] + newNormal2[1])*F1B2
                                   + (surfCoordNew[surfId(4,2,dim)] - surfCoordOld[surfId(0,2,dim)])*(oldNormal[2] + newNormal1[2] + newNormal2[2])*F1B2);



          //////////////////////////////
          /////// SUB JACOBIAN 2 ///////
          //////////////////////////////

          //(i+1/2,j+1/4,k-1/4)

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordOld[surfId(5,isd,dim)] - surfCoordOld[surfId(1,isd,dim)]);
            diag2[isd] = (surfCoordOld[surfId(4,isd,dim)] - surfCoordOld[surfId(2,isd,dim)]);
          }

          this->crossProduct(oldNormal, diag1, diag2);

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordNew[surfId(4,isd,dim)] - surfCoordOld[surfId(1,isd,dim)]);
            diag2[isd] = (surfCoordNew[surfId(1,isd,dim)] - surfCoordOld[surfId(4,isd,dim)]);
          }

          this->crossProduct(newNormal1, diag1, diag2);

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordNew[surfId(1,isd,dim)] - surfCoordOld[surfId(2,isd,dim)]);
            diag2[isd] = (surfCoordNew[surfId(2,isd,dim)] - surfCoordOld[surfId(1,isd,dim)]);
          }

          this->crossProduct(newNormal2, diag1, diag2);

          subJacobian2 = F1B3 * (  (surfCoordNew[surfId(5,0,dim)] - surfCoordOld[surfId(1,0,dim)])*(oldNormal[0] + newNormal1[0] + newNormal2[0])*F1B2
                                   + (surfCoordNew[surfId(5,1,dim)] - surfCoordOld[surfId(1,1,dim)])*(oldNormal[1] + newNormal1[1] + newNormal2[1])*F1B2
                                   + (surfCoordNew[surfId(5,2,dim)] - surfCoordOld[surfId(1,2,dim)])*(oldNormal[2] + newNormal1[2] + newNormal2[2])*F1B2);


          //////////////////////////////
          /////// SUB JACOBIAN 3 ///////
          //////////////////////////////

          //(i+1/2,j-1/4,k+1/4)

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordOld[surfId(7,isd,dim)] - surfCoordOld[surfId(3,isd,dim)]);
            diag2[isd] = (surfCoordOld[surfId(6,isd,dim)] - surfCoordOld[surfId(4,isd,dim)]);
          }

          this->crossProduct(oldNormal, diag1, diag2);

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordNew[surfId(6,isd,dim)] - surfCoordOld[surfId(3,isd,dim)]);
            diag2[isd] = (surfCoordNew[surfId(3,isd,dim)] - surfCoordOld[surfId(6,isd,dim)]);
          }

          this->crossProduct(newNormal1, diag1, diag2);

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordNew[surfId(3,isd,dim)] - surfCoordOld[surfId(4,isd,dim)]);
            diag2[isd] = (surfCoordNew[surfId(4,isd,dim)] - surfCoordOld[surfId(3,isd,dim)]);
          }

          this->crossProduct(newNormal2, diag1, diag2);

          subJacobian3 = F1B3 * (  (surfCoordNew[surfId(7,0,dim)] - surfCoordOld[surfId(3,0,dim)])*(oldNormal[0] + newNormal1[0] + newNormal2[0])*F1B2
                                   + (surfCoordNew[surfId(7,1,dim)] - surfCoordOld[surfId(3,1,dim)])*(oldNormal[1] + newNormal1[1] + newNormal2[1])*F1B2
                                   + (surfCoordNew[surfId(7,2,dim)] - surfCoordOld[surfId(3,2,dim)])*(oldNormal[2] + newNormal1[2] + newNormal2[2])*F1B2);



          //////////////////////////////
          /////// SUB JACOBIAN 4 ///////
          //////////////////////////////

          //(i+1/2,j+1/4,k+1/4)

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordOld[surfId(8,isd,dim)] - surfCoordOld[surfId(4,isd,dim)]);
            diag2[isd] = (surfCoordOld[surfId(7,isd,dim)] - surfCoordOld[surfId(5,isd,dim)]);
          }

          this->crossProduct(oldNormal, diag1, diag2);

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordNew[surfId(7,isd,dim)] - surfCoordOld[surfId(4,isd,dim)]);
            diag2[isd] = (surfCoordNew[surfId(4,isd,dim)] - surfCoordOld[surfId(7,isd,dim)]);
          }

          this->crossProduct(newNormal1, diag1, diag2);

          for( ZFSInt isd = xsd; isd < nDim; isd++ ) {
            diag1[isd] = (surfCoordNew[surfId(4,isd,dim)] - surfCoordOld[surfId(5,isd,dim)]);
            diag2[isd] = (surfCoordNew[surfId(5,isd,dim)] - surfCoordOld[surfId(4,isd,dim)]);
          }

          this->crossProduct(newNormal2, diag1, diag2);

          subJacobian4 = F1B3 * (  (surfCoordNew[surfId(8,0,dim)] - surfCoordOld[surfId(4,0,dim)])*(oldNormal[0] + newNormal1[0] + newNormal2[0])*F1B2
                                   + (surfCoordNew[surfId(8,1,dim)] - surfCoordOld[surfId(4,1,dim)])*(oldNormal[1] + newNormal1[1] + newNormal2[1])*F1B2
                                   + (surfCoordNew[surfId(8,2,dim)] - surfCoordOld[surfId(4,2,dim)])*(oldNormal[2] + newNormal1[2] + newNormal2[2])*F1B2);

          //Compute volume-flux from 4 sub-Jacobians
          m_cells->dxt[dim][cellId] = (subJacobian1 + subJacobian2 + subJacobian3 + subJacobian4)*frk;

        }
      }
    }
  }

}


void ZFSStrctrdBlck3D::allocateAuxDataMaps()
{
  TRACE();
  if(m_bCl || m_bCd){
    if(!m_bCfCpCoeff){
      m_bCfCpCoeff=1;
      zfs_log << "WARNING:: Cf & Cp compuatation were enabeled but not activated in property file (required for lift and drag computation)" << endl;
    }
  }

  if(m_bCfCpCoeff){
    ZFSId noFields = 3;
    if(m_detailAuxData) {
      noFields = 9;
    }

    ZFSId noAuxDataMaps = m_windowInfo->physicalAuxDataMap.size();
    zfsAlloc(m_cells->cfOffsets, noAuxDataMaps, "m_cells->cfOffsets", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_cells->cpOffsets, noAuxDataMaps, "m_cells->cpOffsets", 0, __CALLING_FUNCTION__);
    zfsAlloc(m_cells->powerOffsets, noAuxDataMaps, "m_cells->powerOffsets", 0, __CALLING_FUNCTION__);
    ZFSInt totalSizeCf = 0, totalSizeCp=0, totalSizePower=0;
    for(ZFSId i=0; i<noAuxDataMaps; ++i){
      ZFSInt dataSize=1;
      for(ZFSId j=0; j<nDim; ++j){
        if(m_windowInfo->physicalAuxDataMap[i]->end1[j]==m_windowInfo->physicalAuxDataMap[i]->start1[j]) {
          continue;
        }
        dataSize*=m_windowInfo->physicalAuxDataMap[i]->end1[j] - m_windowInfo->physicalAuxDataMap[i]->start1[j];
      }
      m_cells->cfOffsets[i] = totalSizeCf;
      m_cells->cpOffsets[i] = totalSizeCp;
      m_cells->powerOffsets[i] = totalSizePower;
      totalSizeCf += dataSize*noFields;
      totalSizeCp += dataSize;
      totalSizePower +=dataSize*nDim;
    }

    if(m_bCfCpCoeff) {
      zfsAlloc(m_cells->cf, totalSizeCf, "m_cells->cf", -1.23456123456, __CALLING_FUNCTION__);
      zfsAlloc(m_cells->cp, totalSizeCp, "m_cells->cp", -1.23456123456, __CALLING_FUNCTION__);
    }
    
    if(m_bPower){
     zfsAlloc(m_cells->powerVisc, totalSizePower, "m_cells->power", -1.23456123456, __CALLING_FUNCTION__);
     zfsAlloc(m_cells->powerPres, totalSizePower, "m_cells->power", -1.23456123456, __CALLING_FUNCTION__);
    } 
  }

  //compute Domain Width for Averaging
  if(m_bCpLineAveraging) {
    //if we only compute cd,cl for a part of the domain (m_auxDataCoordinateLimits == true)
    //only compute the average with the width of this section
    if(m_auxDataCoordinateLimits) {
      m_globalDomainWidth = fabs(m_auxDataLimits[3]-m_auxDataLimits[2]);
    } else {
      ZFSFloat minCoordinate = F0, maxCoordinate = F0, minCoordinateGlobal = F0, maxCoordinateGlobal = F0;
      ZFSId lowPoint = -1, highPoint = -1;

      if(m_cpAveragingDir == 0) {
        highPoint = getPointIdFromCell(m_nCells[2]-m_noGhostLayers,0,0);
        lowPoint = getPointIdFromCell(m_noGhostLayers,0,0);
      } else if(m_cpAveragingDir == 1) {
        highPoint = getPointIdFromCell(0,m_nCells[1]-m_noGhostLayers,0);
        lowPoint = getPointIdFromCell(0,m_noGhostLayers,0);
      } else {
        highPoint = getPointIdFromCell(m_noGhostLayers,m_noGhostLayers,m_nCells[0]-m_noGhostLayers);
        lowPoint = getPointIdFromCell(m_noGhostLayers,m_noGhostLayers,m_noGhostLayers);
      }

      minCoordinate = m_coordinates[m_cpAveragingDir][lowPoint];
      maxCoordinate = m_coordinates[m_cpAveragingDir][highPoint];
      MPI_Allreduce(&minCoordinate, &minCoordinateGlobal, 1, MPI_DOUBLE, MPI_MIN, m_zfsStrctrdComm);
      MPI_Allreduce(&maxCoordinate, &maxCoordinateGlobal, 1, MPI_DOUBLE, MPI_MAX, m_zfsStrctrdComm);
      m_globalDomainWidth = fabs(maxCoordinateGlobal - minCoordinateGlobal);
      zfs_log << "Global domain width: " << m_globalDomainWidth << endl;
    }
  }
}

//! Computes the maxResiduum for all cells
/** This function computes the maxResidual using
 * Res = deltaT/(CFL*VolOfCell) * |RHS|
 *
 * with deltaT depending on local or global time stepping
 * is used.
 * checks if the computed max density residual
 * is below the convergence criterion and returns
 * boolean variable
 *
 * method in analogy to fvblock3d but sligthly optimized
 * for the structured block usage
 *
 * Last change: Pascal Meysonnat, April 16, 2012
 *
 */

bool ZFSStrctrdBlck3D::maxResidual()
{
  TRACE();
   
  if( globalTimeStep % m_residualInterval != 0) return true;
  ZFSFloat epsilon = pow( 10.0, -10.0 );
  m_avrgResidual=F0;
  ZFSId cellId=F0;
  ZFSFloat tmpResidual = F0;
  ZFSFloat maxResidual1 =F0;
  ZFSId maxResIndex[3];
  //ZFSId localCounter=F0;
  ZFSFloat maxResidualOrg=F0;
  ZFSFloat localMaxResidual=F0;
  ZFSFloat localAvrgResidual=F0;
  ZFSFloat accumAvrgResidual=F0;
  ZFSFloat globalMaxResidual=F0;
  m_workload = F0;
  //ZFSInt accumCounter=0;
  for(ZFSId dim=0; dim<nDim; dim++) {
    maxResIndex[dim]=F0;
  }

  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
      for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
        cellId=cellIndex(i,j,k);
        //cerr << cellId << endl;
        tmpResidual = m_timeStep / (m_cfl *m_cells->cellJac[cellId] )* fabs( m_cells->rightHandSide[CV->RHO][cellId]);
        m_avrgResidual += tmpResidual;
        
        if( tmpResidual > maxResidual1 ) {
          maxResIndex[0]=i-m_noGhostLayers;
          maxResIndex[1]=j-m_noGhostLayers;
          maxResIndex[2]=k-m_noGhostLayers;
          maxResidual1=tmpResidual;
        }
      }
    }
  }

  //localCounter = counter;
  localMaxResidual=maxResidual1;
  localAvrgResidual=m_avrgResidual;
  //reset average Residual
  m_avrgResidual=F0;

  ZFSFloat localTotalEnergy =F0;
  ZFSFloat globalTotalEnergy =F0;
  ZFSFloat globalPressure =F0;
  ZFSFloat localPressure = F0;
  if(m_initialCondition==101){
    localTotalEnergy=computeTotalKineticEngergy();
    MPI_Allreduce(&localTotalEnergy, &globalTotalEnergy, 1, MPI_DOUBLE, MPI_SUM, m_zfsStrctrdComm);
    globalTotalEnergy/=((16.0*pow(4*atan(1),3.0)));//divided by the overall volume
    localPressure=computeTotalPressure();
    MPI_Allreduce(&localPressure, &globalPressure, 1, MPI_DOUBLE, MPI_SUM, m_zfsStrctrdComm);
    globalPressure/=((16.0*pow(4*atan(1),3.0)));//divided by the overall volume
  }

  //if( noDomains()>1 )
  //  {
  //MPI_Allreduce(m_residualSnd, &m_residualRcv, 1, m_mpiStruct, m_resOp, m_zfsStrctrdComm);
  MPI_Allreduce(&localAvrgResidual, &accumAvrgResidual, 1, MPI_DOUBLE, MPI_SUM, m_zfsStrctrdComm);
  MPI_Allreduce(&localMaxResidual, &globalMaxResidual, 1, MPI_DOUBLE, MPI_MAX, m_zfsStrctrdComm);

  m_avrgResidual = accumAvrgResidual;//m_residualRcv.avrgRes;
  maxResidualOrg=globalMaxResidual;
  //globalMaxResidual=globalMaxResidual;//m_residualRcv.maxRes;
  //for(ZFSId i=0; i<3; i++)
  //{
  //  maxResIndex[i]=m_residualRcv.maxCellIndex[i];
  //}

  // }

  //cout << "m_avrgResidual = " << m_avrgResidual<< " | totalCells " << m_totalGridCells <<endl;
  m_avrgResidual=m_avrgResidual/m_totalGridCells;
  // write first residuals;
  if( ABS(m_firstMaxResidual) < epsilon )
  {
    m_firstMaxResidual = zfsMAX(epsilon,globalMaxResidual);
    m_firstAvrgResidual =zfsMAX(epsilon,m_avrgResidual);
    if(m_initialCondition != 101){ //we need an extra treatment of TGV because of symmetie
      if(approx(localMaxResidual, maxResidualOrg, m_eps)){//so only cpu with the max writes out ==> no need to communicate the max index[i]

        //write out values into residual file
        FILE* f_residual;
        f_residual = fopen("./Residual", "a+");
        fprintf(f_residual, "#MaxRes_1: %1.10e \n", m_firstMaxResidual);
        fprintf(f_residual, "#MaxAvgRes_1: %1.10e \n",m_firstAvrgResidual );
        fprintf(f_residual, "#iter, physTime, time, dT, wLoad, avrgRes, maxRes, blockId, i, j, k ");
        fclose(f_residual);
      }
    }
    else{
      if(domainId()==0){
        //write out values into residual file
        FILE* f_residual;
        f_residual = fopen("./Residual", "a+");
        fprintf(f_residual, "#MaxRes_1: %1.10e \n", m_firstMaxResidual);
        fprintf(f_residual, "#MaxAvgRes_1: %1.10e \n",m_firstAvrgResidual );
        fprintf(f_residual, "#iter, physTime, time, dT, wLoad, avrgRes, maxRes, blockId, i, j, k, k_mean, p_mean, pProbe ");
      }
    }
  }

  // normalize residuals
  globalMaxResidual = globalMaxResidual / m_firstMaxResidual;
  m_avrgResidual = (m_avrgResidual / m_firstAvrgResidual);

  //question if "( m_avrgResidual >= F0 || m_avrgResidual < F0 ) {} else {"
  //is better to capture the also inf???

  if(std::isnan(m_avrgResidual))
  {
    cerr << "Solution diverged, average residual is nan " << endl;
    zfs_log << "Solution diverged, average residual is nan " << endl;
    saveOutput(true);
    zfsTerm(1,__CALLING_FUNCTION__,"Solution diverged, average residual is nan ");
  }

  //convergence Check

  m_convergence=false;
  if( maxResidual1 < m_convergenceCriterion )
  {
    m_convergence = true;
  }

  //need again special treatment for TGV due to symmetrie many processors would write out (!!!!IMPORTANT NO CORRECT INDEX WILL BE WRITTEN OUT )
  if(m_initialCondition!=101){
    //processor with the highest Residual writes out!!! saves communication;
    if(approx(localMaxResidual, maxResidualOrg, m_eps)){//so only cpu with the max writes out ==> no need to communicate the max index[i]
      //write out values into residual file
      FILE* f_residual;
      f_residual = fopen("./Residual", "a+");
      fprintf(f_residual, "%d", globalTimeStep);
      fprintf(f_residual, " %f", m_physicalTime);
      fprintf(f_residual, " %f", m_time);
      fprintf(f_residual, " %f", m_timeStep);
      fprintf(f_residual, " %f", m_workload);
      fprintf(f_residual, " %1.10e", m_avrgResidual);
      fprintf(f_residual, " %1.10e", globalMaxResidual);
      fprintf(f_residual, " %d", m_inputBlockId);
      fprintf(f_residual, " %d", m_nOffsetCells[2]+maxResIndex[0]);//i
      fprintf(f_residual, " %d", m_nOffsetCells[1]+maxResIndex[1]);//j
      fprintf(f_residual, " %d", m_nOffsetCells[0]+maxResIndex[2]);//k
      fprintf(f_residual, "\n");
      fclose(f_residual);
    }
  }else{
    if(domainId()==0){
      ZFSFloat dissip =F0;
      if(globalTimeStep==1){m_kineticEOld=1.2474901617e-03;}
      else{
        dissip=(globalTotalEnergy-m_kineticEOld)/m_timeStep;
        m_kineticEOld=globalTotalEnergy;
      }
      //write out values into residual file
      //compute the dissipation rate

      FILE* f_residual;
      f_residual = fopen("./Residual", "a+");
      fprintf(f_residual, "%d", globalTimeStep);
      fprintf(f_residual, " %f", m_physicalTime);
      fprintf(f_residual, " %f", m_time);
      fprintf(f_residual, " %f", m_timeStep);
      fprintf(f_residual, " %f", m_workload);
      fprintf(f_residual, " %1.10e", m_avrgResidual);
      fprintf(f_residual, " %1.10e", globalMaxResidual);
      fprintf(f_residual, " %d", m_inputBlockId);
      fprintf(f_residual, " %d", m_nOffsetCells[2]+maxResIndex[0]);//i Will be wrong
      fprintf(f_residual, " %d", m_nOffsetCells[1]+maxResIndex[1]);//j Will be wrong
      fprintf(f_residual, " %d", m_nOffsetCells[0]+maxResIndex[2]);//k Will be wrong
      fprintf(f_residual, " %1.10e", globalTotalEnergy); //kinetic Energy
      fprintf(f_residual, " %1.10e", globalPressure); //averaged pressure
      fprintf(f_residual, " %1.10e", dissip); //dissipation rate
      fprintf(f_residual, "\n");
      fclose(f_residual);
    }
  }

  if( maxResidual1 < m_convergenceCriterion ) {
    return true;
  } else {
    return false;
  }
}

ZFSFloat ZFSStrctrdBlck3D::computeTotalKineticEngergy(){
  TRACE();
  ZFSFloat localEnergy=F0;
   for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++){
     for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++){
       for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++){
	 const ZFSId cellId=cellIndex(i,j,k);
	 localEnergy+= (POW2(m_cells->pvariables[PV->U][cellId])+
                        POW2(m_cells->pvariables[PV->V][cellId])+
                        POW2(m_cells->pvariables[PV->W][cellId]))*m_cells->cellJac[cellId];
       }
     }
   }
   return localEnergy;
}

ZFSFloat ZFSStrctrdBlck3D::computeTotalPressure(){
  TRACE();
  ZFSFloat localPressure=F0;
  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++){
    for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++){
      for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++){
        const ZFSId cellId=cellIndex(i,j,k);
        localPressure+= m_cells->pvariables[PV->P][cellId]*m_cells->cellJac[cellId];
      }
    }
  }
  return localPressure;
}

inline ZFSFloat ZFSStrctrdBlck3D::pressure(ZFSId cellId){
  return m_cells->pvariables[PV->P][cellId];
}


void ZFSStrctrdBlck3D::ZFSResOperationFunction(ZFSRes* in, ZFSRes* out, int* len)
{
  TRACE();
  ZFSRes res;
  res.maxRes=F0;
  res.avrgRes=F0;
  res.maxCellIndex = new ZFSInt[3];
  //cout << "nDim " << nDim;
  for(ZFSInt i=0; i<3; i++)
  {
    res.maxCellIndex[i]=0;
  }

  for(ZFSInt j=0; j<*len; j++)
  {

    if(res.maxRes< (*in).maxRes)
    {
      res.maxRes=(*in).maxRes;
      for(ZFSInt i=0; i<3; i++)
      {
        res.maxCellIndex[i]=in->maxCellIndex[i];
      }
    }

    res.avrgRes=in->avrgRes+out->avrgRes;
    out->avrgRes=res.avrgRes;
    out->maxRes=res.maxRes;
    for(ZFSInt i=0; i<3; i++)
    {
      out->maxCellIndex[i]=res.maxCellIndex[i];
    }
    in++;
    out++;
  }
   delete[] res.maxCellIndex;
}


inline void ZFSStrctrdBlck3D::crossProduct( ZFSFloat* result, ZFSFloat* vec1,
                                            ZFSFloat* vec2)
{
  result[xsd] = vec1[ysd] * vec2[zsd] - vec1[zsd] * vec2[ysd];
  result[ysd] = vec1[zsd] * vec2[xsd] - vec1[xsd] * vec2[zsd];
  result[zsd] = vec1[xsd] * vec2[ysd] - vec1[ysd] * vec2[xsd];
}

inline ZFSId ZFSStrctrdBlck3D::getPointIdFromCell( ZFSInt i, ZFSInt j, ZFSInt k ){
  return i + ( k * ( m_nCells[1] + 1 ) + j ) * ( m_nCells[2] + 1 );
}

inline ZFSId ZFSStrctrdBlck3D::getPointIdfromPoint( ZFSId origin, ZFSInt incI,
                                                    ZFSInt incJ, ZFSInt incK )
{
  return origin + incI + incJ * m_nPoints[2] + incK * m_nPoints[2] * m_nPoints[1];
}

void ZFSStrctrdBlck3D:: revertMetrics()
{
  TRACE();
  ZFSFloat **tmpMetrics;

  zfsAlloc(tmpMetrics, 9, m_nCells[0] * m_nCells[1] * m_nCells[2], "revertSurfaceMetrics", 234.234, __CALLING_FUNCTION__);

  for( ZFSInt i = 0; i < m_noStrctrdCells; i++ )
  {
    for( ZFSInt j = 0; j < 9; j++ )
    {
      tmpMetrics[j][i] = this->m_cells->surfaceMetrics[i][j];
    }
  }

  zfsDeallocate(m_cells->surfaceMetrics );
  m_cells->surfaceMetrics = tmpMetrics;
  tmpMetrics = NULL;
}

void ZFSStrctrdBlck3D::moveGrid(ZFSBool isRestart, ZFSBool zeroPos)
{
  TRACE();
  ZFSFloat pi=4.0*atan(1);
  ZFSFloat t = m_time + m_timeStep*m_RKalpha[m_RKStep];

  if(isRestart) {
    t = m_time;
  }

  switch(m_gridMovingMethod)
  {
  case 0: //travelling wave   ?
  {
    //we need some relaxation function in wall-normal direction (here take j)
    //choose for example tanh(y)
    cout << "time = " << t << endl;
    //ZFSFloat x=F0;
    ZFSFloat y=F0;
    ZFSFloat z=F0;
    //ZFSFloat p=F0;
    ZFSId pointId=0;
    ZFSFloat amp=0;
    ZFSFloat l=F0;
    cout << "moving the grid" << endl;
    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId j=0; j<m_nPoints[1]; j++)
      {
        for(ZFSId i=0; i<m_nPoints[2]; i++)
        {
          pointId=pointIndex(i,j,k);
          //x=m_coordinates[0][cellId];
          y=m_mgOrgCoordinates[0][j];
          if(y>2.0) break;
          l=(m_mgOrgCoordinates[0][j+1]-y);
          z=m_mgOrgCoordinates[1][k];
          amp=(1.0-((tanh(16*y-3.75)+1.0)/2.0));

          m_coordinates[1][pointId]=y+sin(4*pi*z-0.5*t)*amp*0.025;
          pointId=pointIndex(i,j+1,k);
          ZFSFloat m=1.0/(4*pi*cos(4*pi*z-0.5*t)*amp*0.025);
          if(m>F0)
          {
            m_coordinates[2][pointId]=z-sqrt(POW2(l)/(1+POW2(m)));
          }
          else
          {
            m_coordinates[2][pointId]=z+sqrt(POW2(l)/(1+POW2(m)));
          }

        }
      }
    }

    break;
  }
  case 1:  // oscillating grid in x-direction
  {
    //we need some relaxation function in wall-normal direction
    // here we take a linear function: ( y_max - y ) / y_max
    ZFSId pointId = 0;
    ZFSFloat x = F0;
    ZFSFloat y_max = 10.0;  // y of upper domain boundary
    ZFSFloat frequency = F1 / ( F2*m_timeStep*1000.0 );
    ZFSFloat amp = PV->UInfinity / ( F2*pi*frequency );

    //cout<<"frequency = "<<frequency<<" amp= "<<amp<<" PV->UInfinity = "<<PV->UInfinity<<endl;

    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId j=0; j<(m_nPoints[1]); j++)
      {
        for(ZFSId i=0; i<m_nPoints[2]; i++)
        {// its important to move also the ghostcell's points in this case since these are needed for periodic bc (don't mirror in this case!)
          pointId = pointIndex(i,j,k);
          //x=m_coordinates[0][cellId];
          x = m_mgOrgCoordinates[0][i];
          //m_coordinates[1][pointId]=m_coordinates[1][pointId]+(0.5-y)*l;
          m_coordinates[0][pointId] = x + ( 1 - m_coordinates[1][pointId] / y_max )*amp*sin( 2*pi*t*frequency );
        }
      }
    }

    break;
  }
  case 2:  // channel with moving indentation
  {
    //break;//comment this and compVolFlux in / out for mg /not mg
    //}
    //case 222:
    //{
    ZFSId pointId = 0;
    ZFSFloat y = F0;
    ZFSFloat h = F0, g = F0, beta = 4.14l;
    ZFSFloat x2, x3, x4, x5;
    ZFSFloat StrNum = m_wallVel; // Strouhal Number, set by wallVel
    ZFSFloat ver = 0.0l; // ver can be used to translate x coordinate
    ZFSFloat y_max = 1.0l;  //can be used to scale (old mesh had y_max = m_mgOrgCoordinates[0][32] = 1/30)

    // values from Ralph and Pedley:
    x2 = y_max*( ver - 11.75l );
    x3 = y_max*( ver - 9.25l );
    x4 = y_max*( ver - 1.25 );
    x5 = y_max*( ver + 1.25 );
    // 1000.0l*
    h = 0.38l*( 1.0l - cos( PV->UInfinity*StrNum*2.0l*pi*t/y_max) ) / 2.0l; // eps = 0.38 like in Kwak

    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId j=0; j<m_nPoints[1]; j++)
      {
        y=m_mgOrgCoordinates[0][j];        // original coordinates y
        for(ZFSId i=0; i<m_nPoints[2]; i++)
        {
          pointId = pointIndex(i,j,k);

          g = F0;
          if (m_coordinates[0][pointId] > x2 && m_coordinates[0][pointId] <x5)
          {
            g = ( (m_coordinates[0][pointId] < x3) ? (1.0l+tanhl(beta * (m_coordinates[0][pointId]-(x2+x3)/2.0l)/y_max) )/2.0l : ((m_coordinates[0][pointId] < x4) ? 1.0l : (1.0l-tanhl(beta * (m_coordinates[0][pointId]-(x4+x5)/2.0l)/y_max) )/2.0l ) );
          }

          m_coordinates[1][pointId]=y*(1.0l-h*g);
        }
      }
    }

    break;
  }
  case 3:  // piston moving in x-direction
  {
    ZFSId pointId = 0;
    ZFSFloat x = F0;

    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId j=0; j<m_nPoints[1]; j++)
      {
        for(ZFSId i=0; i<m_nPoints[2]; i++)
        {
          pointId = pointIndex(i,j,k);
          x = m_mgInitCoordinates[0][pointId];

          m_coordinates[0][pointId] = x*( 1 + t*m_wallVel );//x+t*m_wallVel*(i-2)/(m_nPoints[2]-5);
        }
      }
    }

    break;
  }
  case 4:  // inner grid movement
  {/*
     break;//comment this and compVolFlux in / out for mg /not mg
     }
     case 222:
     {*/
    ZFSId pointId = 0;
    ZFSFloat x = F0, y = F0;
    ZFSFloat h = F0, g = F0, beta = 16.0l;
    ZFSFloat x2, x3, x4, x5;
    ZFSFloat StrNum = m_wallVel; // Strouhal Number, set by wallVel
    ZFSFloat ver = 0.0l; // ver can be used to translate x coordinate
    ZFSFloat y_max = 1.0l; // can be used to scale

    //for Square:
    x2 = y_max*( ver+0.1l );
    x3 = y_max*( ver+0.5l );
    x4 = y_max*( ver+0.5l );
    x5 = y_max*( ver+0.9l );

    h = 0.35l*( 1.0l - cos( m_Ma * sqrt( PV->TInfinity )*StrNum*2.0l*pi*t/y_max ) ) / 2.0l;

    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId j=0; j<m_nPoints[1]; j++)
      {
        for(ZFSId i=0; i<m_nPoints[2]; i++)
        {
          pointId=pointIndex(i,j,k);
          // original coordinates
          //x=m_mgOrgCoordinates[0][i];
          //y=m_mgOrgCoordinates[1][j];

          x = m_mgInitCoordinates[0][pointId];
          y = m_mgInitCoordinates[1][pointId];

          g=F0;
          if (y > x2 && y < x5 && x > x2 && x < x5)
          {
            g = ( (y < x3) ? (1.0l+tanhl(beta * (y-(x2+x3)/2.0l)/y_max) )/2.0l : ((y < x4) ? 1.0l : (1.0l-tanhl(beta * (y-(x4+x5)/2.0l)/y_max) )/2.0l ) );
          }
          m_coordinates[0][pointId]=x*(1.0l-h*g*(1.0l- x));

          g=F0;
          if (x > x2 && x < x5 && y > x2 && y < x5)
          {
            g = ( (x < x3) ? (1.0l+tanhl(beta * (x-(x2+x3)/2.0l)/y_max) )/2.0l : ((x < x4) ? 1.0l : (1.0l-tanhl(beta * (x-(x4+x5)/2.0l)/y_max) )/2.0l ) );
          }
          m_coordinates[1][pointId]=y*(1.0l-h*g*(1.0l- y));
        }
      }
    }

    break;
  }
  case 5:  // makos point-actuation - gauss distribution
  {
    // for ReTheta = 1000
    // experiment: ReTheta = 10.000, Ma = 0.04, amplitude = 0.07mm, frequency = 2000, A+ = 24.3, T+ = 44.17;
    // numerical: Ma = 0.2, amplitude = 0.5260, frequency = 0.04831 for same A+ and T+ (based on non-convection time  t = m_time / m_timeRef)
    // gauss distribution function with sigma = 4 to be reached at a radius of 2cm
    ZFSFloat t_offset = t - m_movingGridTimeOffset;
    ZFSFloat frequency = 0.04831; // for T+ = 44.17

    ZFSFloat pos = (F1-cos(F2*pi*t_offset*frequency))/F2;

    ZFSFloat maxGridMoveHeight = 10.0; // max height of grid points which are affected by grid movement
    ZFSFloat maxGridMoveRadius = 10.0; // max radius of grid points which are affected by grid movement
    ZFSFloat centerX = 0.0;
    ZFSFloat centerZ = 10.71576;
    ZFSFloat amplitude = 0.5260;
    ZFSFloat sigma = m_wallVel;

    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId j=0; j<m_nPoints[1]; j++)
      {
        for(ZFSId i=0; i<m_nPoints[2]; i++)
        {
          ZFSId pointId=pointIndex(i,j,k);

          ZFSFloat distX = centerX - m_mgInitCoordinates[0][pointId];
          ZFSFloat distZ = centerZ - m_mgInitCoordinates[2][pointId];
          ZFSFloat radius = sqrt(POW2(distX) + POW2(distZ));
          ZFSFloat height = amplitude*exp((-0.5)*(POW2(radius/sigma))) * pos;

          if(m_mgInitCoordinates[1][pointId] < maxGridMoveHeight && radius < maxGridMoveRadius) {
            //move grid point in y-direction with scaling factor
            //only grid points that are roughly inside the BL are moved
            m_coordinates[1][pointId] = m_mgInitCoordinates[1][pointId] + (1.0 - (m_mgInitCoordinates[1][pointId]/maxGridMoveHeight))*height;
          }
        }
      }
    }

    break;
  }
  case 6: // makos line actuation - gauss distribution
  {
    // for ReTheta = 1000
    // experiment: ReTheta = 10.000, Ma = 0.04, amplitude = 0.07mm, frequency = 2000, A+ = 24.3, T+ = 44.17;
    // numerical: Ma = 0.2, amplitude = 0.5260, frequency = 0.04831 for same A+ and T+ (based on non-convection time  t = m_time / m_timeRef)
    // gauss distribution function with sigma = 4 to be reached at a radius of 2cm
    ZFSFloat t_offset = t - m_movingGridTimeOffset;
    ZFSFloat frequency = 0.04831; // for T+ = 44.17

    ZFSFloat pos = (F1-cos(F2*pi*t_offset*frequency))/F2;

    ZFSFloat maxGridMoveHeight = 10.0; // max height of grid points which are affected by grid movement
    ZFSFloat maxGridMoveRadius = 10.0; // max radius of grid points which are affected by grid movement
    ZFSFloat centerX = 0.0;
    ZFSFloat amplitude = 0.5260; //A+ = 24.3
    ZFSFloat sigma = m_wallVel;

    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId j=0; j<m_nPoints[1]; j++)
      {
        for(ZFSId i=0; i<m_nPoints[2]; i++)
        {
          ZFSId pointId=pointIndex(i,j,k);

          ZFSFloat distX = centerX - m_mgInitCoordinates[0][pointId];
          ZFSFloat radius = sqrt(POW2(distX));
          ZFSFloat height = amplitude*exp((-0.5)*(POW2(radius/sigma))) * pos;

          if(m_mgInitCoordinates[1][pointId] < maxGridMoveHeight && radius < maxGridMoveRadius) {
            //move grid point in y-direction with scaling factor
            //only grid points that are roughly inside the BL are moved
            m_coordinates[1][pointId] = m_mgInitCoordinates[1][pointId] + (1.0 - (m_mgInitCoordinates[1][pointId]/maxGridMoveHeight))*height;
          }
        }
      }
    }
    break;
  }
  case 7:  // makos point-actuation - realistic distribution with Bessel functions
  {
    // for ReTheta = 1000
    // experiment: ReTheta = 10.000, Ma = 0.04, amplitude = 0.07mm, frequency = 2000, A+ = 24.3, T+ = 44.17;
    // numerical: Ma = 0.2, amplitude = 0.5260, frequency = 0.04831 for same A+ and T+ (based on non-convection time  t = m_time / m_timeRef)
    // gauss distribution function with sigma = 4 to be reached at a radius of 2cm
    ZFSFloat t_offset = t - m_movingGridTimeOffset;
    ZFSFloat frequency = 0.04831; // for T+ = 44.17

    //ZFSFloat pos = (F1-cos(F2*pi*t_offset*frequency))/F2;
    ZFSFloat pos = cos(pi/F2 + F2*pi*t_offset*frequency);

    const ZFSFloat maxGridMoveHeight = 10.0; // max height of grid points which are affected by grid movement

    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId j=0; j<m_nPoints[1]; j++)
      {
        for(ZFSId i=0; i<m_nPoints[2]; i++)
        {
          ZFSId pointId=pointIndex(i,j,k);

          ZFSFloat distX = m_makosCenterX - m_mgInitCoordinates[0][pointId];
          ZFSFloat distZ = m_makosCenterZ - m_mgInitCoordinates[2][pointId];
          ZFSFloat radius = sqrt(POW2(distX) + POW2(distZ));

          if(radius <= m_makosActuatorRadius && m_mgInitCoordinates[1][pointId] < maxGridMoveHeight) {
            //move grid point in y-direction with scaling factor
            //only grid points that are roughly inside the BL are moved
            ZFSFloat height = m_makosMaxAmplitude[i+k*m_nPoints[2]]*pos;

            m_coordinates[1][pointId] = m_mgInitCoordinates[1][pointId] + (1.0 - (m_mgInitCoordinates[1][pointId]/maxGridMoveHeight))*height;
          }
        }
      }
    }
    break;
  }
  case 8:   // makos line-actuation - realistic distribution with Bessel functions
  {
    // for ReTheta = 1000
    // experiment: ReTheta = 10.000, Ma = 0.04, amplitude = 0.07mm, frequency = 2000, A+ = 24.3, T+ = 44.17;
    // numerical: Ma = 0.2, amplitude = 0.5260, frequency = 0.04831 for same A+ and T+ (based on non-convection time  t = m_time / m_timeRef)
    // gauss distribution function with sigma = 4 to be reached at a radius of 2cm
    ZFSFloat t_offset = t - m_movingGridTimeOffset;
    ZFSFloat frequency = 0.04831; // for T+ = 44.17

    //ZFSFloat pos = (F1-cos(F2*pi*t_offset*frequency))/F2;
    ZFSFloat pos = cos(pi/F2 + F2*pi*t_offset*frequency);

    const ZFSFloat maxGridMoveHeight = 10.0; // max height of grid points which are affected by grid movement

    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId j=0; j<m_nPoints[1]; j++)
      {
        for(ZFSId i=0; i<m_nPoints[2]; i++)
        {
          ZFSId pointId=pointIndex(i,j,k);

          ZFSFloat distX = m_makosCenterX - m_mgInitCoordinates[0][pointId];
          ZFSFloat radius = sqrt(POW2(distX));

          //if(radius <= m_makosActuatorRadius && m_mgInitCoordinates[1][pointId] < maxGridMoveHeight) {
          if(radius <= m_makosActuatorRadius && approx(m_mgInitCoordinates[1][pointId],F0,m_eps)) {
            //move grid point in y-direction with scaling factor
            //only grid points that are roughly inside the BL are moved
            ZFSFloat height = m_makosMaxAmplitude[i+k*m_nPoints[2]]*pos;

            m_coordinates[1][pointId] = m_mgInitCoordinates[1][pointId] + (1.0 - (m_mgInitCoordinates[1][pointId]/maxGridMoveHeight))*height;
          }
        }
      }
    }
    break;
  }
  case 9: {
    //traveling wave case
    ZFSFloat t_offset = t-m_movingGridTimeOffset;
    if(zeroPos) {t_offset = F0;};
    const ZFSFloat transitionLength = m_waveEndTransition - m_waveBeginTransition;
    const ZFSFloat transitionOutLength = m_waveOutEndTransition - m_waveOutBeginTransition;

    ZFSFloat fadeInFactor = 0;
    const ZFSFloat timeRelaxation = 80.0;

    if (t_offset < timeRelaxation) {
      fadeInFactor = (1.0-cos(t_offset/timeRelaxation*pi))*F1B2;
    } else {
      fadeInFactor = 1.0;
    }

    if(zeroPos) {fadeInFactor = F1;};

    for(ZFSId k=0; k<m_nPoints[0]; k++) {
      for(ZFSId j=0; j<m_nPoints[1]; j++) {
        for(ZFSId i=0; i<m_nPoints[2]; i++) {
          const ZFSId pointId=pointIndex(i,j,k);
          const ZFSFloat xInit = m_mgInitCoordinates[0][pointId];
          ZFSFloat transitionFactor = F0;
          if(xInit <= m_waveBeginTransition) {
            transitionFactor = F0;
          } else if (xInit > m_waveBeginTransition && xInit< m_waveEndTransition) {
            transitionFactor = (1-cos((xInit-m_waveBeginTransition)/transitionLength*pi))*F1B2;
          } else if (m_waveEndTransition <= xInit && xInit <= m_waveOutBeginTransition) {
            transitionFactor = F1;
          } else if (xInit > m_waveOutBeginTransition && xInit< m_waveOutEndTransition) {
            transitionFactor = (1+cos((xInit-m_waveOutBeginTransition)/transitionOutLength*pi))*F1B2;
          } else {
            transitionFactor = F0;
          }

          const ZFSFloat zInit = m_mgInitCoordinates[2][pointId];
          const ZFSFloat yInit = m_mgInitCoordinates[1][pointId];
          m_coordinates[1][pointId] = fadeInFactor*(m_waveAmplitude*transitionFactor*cos((F2*pi)/m_waveLength*(zInit - m_waveSpeed*t_offset))) + yInit;
        }
      }
    }

    break;
  }
  case 10: {
    ZFSFloat t_offset = t-m_movingGridTimeOffset;
    if(zeroPos) {t_offset = F0;};
    const ZFSFloat transitionLength = m_waveEndTransition - m_waveBeginTransition;

    ZFSFloat fadeInFactor = 0;
    const ZFSFloat timeRelaxation = 80.0;

    if (t_offset < timeRelaxation) {
      fadeInFactor = (1.0-cos(t_offset/timeRelaxation*pi))*F1B2;
    } else {
      fadeInFactor = 1.0;
    }

    if(zeroPos) {fadeInFactor = F1;};

    for(ZFSId k=0; k<m_nPoints[0]; k++) {
      for(ZFSId j=0; j<m_nPoints[1]; j++) {
        for(ZFSId i=0; i<m_nPoints[2]; i++) {
          const ZFSId pointId=pointIndex(i,j,k);
          const ZFSFloat xInit = m_mgInitCoordinates[0][pointId];
          ZFSFloat transitionFactor = F0;
          if(xInit <= m_waveBeginTransition) {
            transitionFactor = F0;
          } else if (xInit > m_waveBeginTransition && xInit< m_waveEndTransition) {
            transitionFactor = (1-cos((xInit-m_waveBeginTransition)/transitionLength*pi))*F1B2;
          } else {
            transitionFactor = F1;
          }

          if(xInit > 63.0) {
            transitionFactor = F0;
          } else if (xInit >= 60.0 && xInit <= 63.0) {
            transitionFactor = F1-(F1-cos((xInit-60.0)/3.0*pi))*F1B2;
          }

          const ZFSFloat zInit = m_mgInitCoordinates[2][pointId];
          const ZFSFloat yInit = m_mgInitCoordinates[1][pointId];
          m_coordinates[1][pointId] = fadeInFactor*(m_waveAmplitude*transitionFactor*(cos((F2*pi)/m_waveLength*(zInit - m_waveSpeed*t_offset)) + 1.0)) + yInit;
        }
      }
    }
    break;
  }
  case 11: {
    //traveling wave channel (Tomiyama & Fukagata 2013)
    ZFSFloat t_offset = t-m_movingGridTimeOffset;
    if(zeroPos) {t_offset = F0;};

    ZFSFloat fadeInFactor = 0;
    const ZFSFloat timeRelaxation = 50.0;

    if (t_offset < timeRelaxation) {
      fadeInFactor = (1.0-cos(t_offset/timeRelaxation*pi))*F1B2;
    } else {
      fadeInFactor = 1.0;
    }

    if(zeroPos) {fadeInFactor = F1;};

    for(ZFSId k=0; k<m_nPoints[0]; k++) {
      for(ZFSId j=0; j<m_nPoints[1]; j++) {
        for(ZFSId i=0; i<m_nPoints[2]; i++) {
          const ZFSId pointId=pointIndex(i,j,k);
          const ZFSFloat yInit = m_mgInitCoordinates[1][pointId];
          const ZFSFloat zInit = m_mgInitCoordinates[2][pointId];
          ZFSFloat yRelaxation = F0;

          if(yInit <= F0) {
            yRelaxation = F1;
          } else if(yInit > F0 && yInit < F1) {
            yRelaxation = F1-yInit;
          } else if(yInit > F1 && yInit < F2) {
            yRelaxation = yInit - F1;
          } else {
            yRelaxation = F1;
          }

          if(yInit <= F1) {
            m_coordinates[1][pointId] = yInit + fadeInFactor*(m_waveAmplitude*yRelaxation*cos((F2*pi)/m_waveLength*(zInit - m_waveSpeed*t_offset)));
          } else {
            m_coordinates[1][pointId] = yInit - fadeInFactor*(m_waveAmplitude*yRelaxation*cos((F2*pi)/m_waveLength*(zInit - m_waveSpeed*t_offset)));
          }
        }
      }
    }

    break;
  }
  default:
  {
    zfsTerm(1, __CALLING_FUNCTION__, "Grid Moving Method not implemented!");
  }
  }

  if (m_gridMovingMethod != 1) {extrapolateGhostPointCoordinates();};
  if(noDomains()>1) {exchangePoints();};
  computeCellCentreCoordinates();
  computeMetrics();
  computeJacobian();

}

void ZFSStrctrdBlck3D::initMovingGrid()
{
  TRACE();
  zfs_log  << "in initMovingGrid " << endl;


  //First approach: save whole mesh in m_mgInitCoordinates (for analytical channel with indentation)
  for(ZFSId k=0; k<m_nPoints[0]; ++k) {
    for(ZFSId j=0; j<m_nPoints[1]; ++j) {
      for(ZFSId i=0; i<m_nPoints[2]; ++i) {
        const ZFSId pointId=pointIndex(i,j,k);
        for( ZFSInt isd = xsd; isd < nDim; ++isd) {
          m_mgInitCoordinates[isd][pointId] = m_coordinates[isd][pointId];
        }
      }
    }
  }

  // Second approach: save only parts of the mesh depending on moving grid case
  switch(m_gridMovingMethod)
  {
  case 1: 
  case 3 : // grid moves only in x-direction
  {
    //we need to save the original grid distribution in some direction
    //here only the x-direction is stored, as the points will be moved only in x-direction
    //store one original line of x (it is a rectangular domain so we do not need to store all points)
    ZFSFloat* dummy = new ZFSFloat[m_nPoints[2]+m_nPoints[1]];
    m_mgOrgCoordinates = new ZFSFloat *[2];
    m_mgOrgCoordinates[0]=&dummy[0];
    m_mgOrgCoordinates[1]=&dummy[m_nPoints[1]];

    ZFSInt i=m_noGhostLayers;
    ZFSInt j=0;
    ZFSInt k=m_noGhostLayers;

    //store one i-line
    for(i=0;i<m_nPoints[2];i++) {
      const ZFSId pointId=pointIndex(i,j,k);
      m_mgOrgCoordinates[0][i]=m_coordinates[0][pointId];
    }
    break;
  }
  case 2 : // grid moves only in y-direction
  {
    //we need to save the original grid distribution in some direction
    //here only the y-direction is stored, as the points will be moved only in y-direction
    //store one original line of y (it is a rectangular domain so we do not need to store all points)
    ZFSFloat* dummy = new ZFSFloat[m_nPoints[2]+m_nPoints[1]];
    m_mgOrgCoordinates = new ZFSFloat *[2];
    m_mgOrgCoordinates[0]=&dummy[0];
    m_mgOrgCoordinates[1]=&dummy[m_nPoints[1]];

    ZFSInt i=m_noGhostLayers;
    ZFSInt j=0;
    ZFSInt k=m_noGhostLayers;
    cout<<"Channel with moving indentation, St = "<<m_wallVel<<endl;
    //store one j-line
    for(j=0;j<m_nPoints[1];j++)
    {
      const ZFSId pointId=pointIndex(i,j,k);
      m_mgOrgCoordinates[0][j]=m_coordinates[1][pointId];
    }

    break;
  }
  case 4 : // grid moves in x and y-direction
  {
    //store one original line of i and j (it is a rectangular domain so we do not need to store all points
    ZFSFloat* dummy = new ZFSFloat[m_nPoints[2]+m_nPoints[1]];
    m_mgOrgCoordinates = new ZFSFloat *[2];
    m_mgOrgCoordinates[0]=&dummy[0];
    m_mgOrgCoordinates[1]=&dummy[m_nPoints[2]];


    ZFSInt i=m_noGhostLayers;
    ZFSInt j=0;
    ZFSInt k=m_noGhostLayers;

    //store one i-line
    for(i=0;i<m_nPoints[2];i++)
    {
      const ZFSId pointId=pointIndex(i,j,k);
      m_mgOrgCoordinates[0][i]=m_coordinates[0][pointId];
    }

    //store one j-line
    for(j=0;j<m_nPoints[1];j++)
    {
      const ZFSId pointId=pointIndex(i,j,k);
      m_mgOrgCoordinates[1][j+1]=m_coordinates[1][pointId];
    }

    break;
  }
  case 5:
  case 6:
  {
    break;
  }
  case 7:
  {
    //compute max amplitude height for all cells - point actuation
    m_makosActuatorRadius = 11.275;
    m_makosPeakAmplitude = 0.5260/2.0;
    const ZFSFloat lambda = 3.19622;

    m_makosCenterX = 0.0;
    m_makosCenterZ = 15.4996668;

    zfsAlloc(m_makosMaxAmplitude, m_nPoints[2]*m_nPoints[0], "m_makosMaxAmplitude", F0, __CALLING_FUNCTION__);

    ZFSFloat maxFunctionValue = bessJ(0,0) - (bessJ(0,lambda)/bessI(0,lambda)) * bessI(0,0);

    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId i=0; i<m_nPoints[2]; i++)
      {
        ZFSId j = m_noGhostLayers;
        const ZFSId pointId=pointIndex(i,j,k);

        ZFSFloat distX = m_makosCenterX - m_mgInitCoordinates[0][pointId];
        ZFSFloat distZ = m_makosCenterZ - m_mgInitCoordinates[2][pointId];
        ZFSFloat radius = sqrt(POW2(distX) + POW2(distZ));

        if(radius <= m_makosActuatorRadius) {
          m_makosMaxAmplitude[i+k*m_nPoints[2]] = (bessJ(0,(lambda/m_makosActuatorRadius)*radius) - (bessJ(0,lambda)/bessI(0,lambda)) * bessI(0,(lambda/m_makosActuatorRadius)*radius))/maxFunctionValue * m_makosPeakAmplitude;
        }
      }
    }

    break;
  }
  case 8:
  {
    //compute max amplitude height for all cells - line actuation
    m_makosActuatorRadius = 11.72;
    m_makosPeakAmplitude = 0.5260/2.0;
    const ZFSFloat lambda = 3.19622;

    m_makosCenterX = 0.0;
    m_makosCenterZ = 10.71576;

    zfsAlloc(m_makosMaxAmplitude, m_nPoints[2]*m_nPoints[0], "m_makosMaxAmplitude", F0, __CALLING_FUNCTION__);

    ZFSFloat maxFunctionValue = bessJ(0,0) - (bessJ(0,lambda)/bessI(0,lambda)) * bessI(0,0);

    for(ZFSId k=0; k<m_nPoints[0]; k++)
    {
      for(ZFSId i=0; i<m_nPoints[2]; i++)
      {
        ZFSId j = m_noGhostLayers;
        const ZFSId pointId=pointIndex(i,j,k);

        ZFSFloat distX = m_makosCenterX - m_mgInitCoordinates[0][pointId];
        ZFSFloat radius = sqrt(POW2(distX));

        if(radius <= m_makosActuatorRadius) {
          m_makosMaxAmplitude[i+k*m_nPoints[2]] = (bessJ(0,(lambda/m_makosActuatorRadius)*radius) - (bessJ(0,lambda)/bessI(0,lambda)) * bessI(0,(lambda/m_makosActuatorRadius)*radius))/maxFunctionValue * m_makosPeakAmplitude;
        }
      }
    }

    break;
  }
  case 9:
    //traveling wave defined by viscous units
    //used Smits formula to compute friction velocity
  {
    m_travelingWave = true;
    if(!m_restart){ m_waveTimeStepComputed = false;}
    m_waveSpeed = 0.0;
    m_waveLength = 0.0;
    m_waveAmplitude = 0.0;
    m_waveCellsPerWaveLength = 1;
    if(!m_restart){m_waveNoStepsPerCell = 1;}

    //time needs to be constant for traveling wave
    m_constantTimeStep = true;

    /*! \page propertyPage1
      \section waveLengthPlus
      <code>ZFSInt ZFSStrctrdBlck::m_waveLengthPlus </code>\n
      default = <code> 1.0 </code>\n \n
      Wavelength of the traveling wave in inner units.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>WAVE, MOVING, STRCTRD</i>
    */
    m_waveLengthPlus = 1.0;
    m_waveLengthPlus = *(ZFSContext::getProperty("waveLengthPlus", m_blockId, __CALLING_FUNCTION__, &m_waveLengthPlus)->asFloat(0));

    /*! \page propertyPage1
      \section waveAmplitudePlus
      <code>ZFSInt ZFSStrctrdBlck::m_waveAmplitudePlus </code>\n
      default = <code> 1.0 </code>\n \n
      Amplitude of the traveling wave in inner units.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>WAVE, MOVING, STRCTRD</i>
    */
    m_waveAmplitudePlus = 0.0;
    m_waveAmplitudePlus = *(ZFSContext::getProperty("waveAmplitudePlus", m_blockId, __CALLING_FUNCTION__, &m_waveAmplitudePlus)->asFloat(0));

    /*! \page propertyPage1
      \section waveTimePlus
      <code>ZFSInt ZFSStrctrdBlck::m_waveTimePlus </code>\n
      default = <code> 1.0 </code>\n \n
      Period time of the traveling wave in inner units.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>WAVE, MOVING, STRCTRD</i>
    */
    m_waveTimePlus = 0.0;
    m_waveTimePlus = *(ZFSContext::getProperty("waveTimePlus", m_blockId, __CALLING_FUNCTION__, &m_waveTimePlus)->asFloat(0));

    /*! \page propertyPage1
      \section waveBeginTransition
      <code>ZFSInt ZFSStrctrdBlck::m_waveBeginTransition </code>\n
      default = <code> 1.0 </code>\n \n
      Start of the transition from flat to wave in x-dir.\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>WAVE, MOVING, STRCTRD</i>
    */
    m_waveBeginTransition = 0.0;
    m_waveBeginTransition = *(ZFSContext::getProperty("waveBeginTransition", m_blockId, __CALLING_FUNCTION__, &m_waveBeginTransition)->asFloat(0));

    /*! \page propertyPage1
      \section waveEndTransition
      <code>ZFSInt ZFSStrctrdBlck::m_waveEndTransition </code>\n
      default = <code> 1.0 </code>\n \n
      End of the transition from flat to wave in x-dir.\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>WAVE, MOVING, STRCTRD</i>
    */
    m_waveEndTransition = 0.0;
    m_waveEndTransition = *(ZFSContext::getProperty("waveEndTransition", m_blockId, __CALLING_FUNCTION__, &m_waveEndTransition)->asFloat(0));

    /*! \page propertyPage1
      \section waveOutBeginTransition
      <code>ZFSInt ZFSStrctrdBlck::m_waveOutBeginTransition </code>\n
      default = <code> 1.0 </code>\n \n
      Start of the transition from wave to flat in x-dir.\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>WAVE, MOVING, STRCTRD</i>
    */
    m_waveOutBeginTransition = 1000000.0;
    m_waveOutBeginTransition = *(ZFSContext::getProperty("waveOutBeginTransition", m_blockId, __CALLING_FUNCTION__, &m_waveOutBeginTransition)->asFloat(0));

    /*! \page propertyPage1
      \section waveOutEndTransition
      <code>ZFSInt ZFSStrctrdBlck::m_waveOutEndTransition </code>\n
      default = <code> 1.0 </code>\n \n
      End of the transition from wave to flat in x-dir.\n
      Possible values are:\n
      <ul>
      <li>Float</li>
      </ul>
      Keywords: <i>WAVE, MOVING, STRCTRD</i>
    */
    m_waveOutEndTransition = 2000000.0;
    m_waveOutEndTransition = *(ZFSContext::getProperty("waveOutEndTransition", m_blockId, __CALLING_FUNCTION__, &m_waveOutEndTransition)->asFloat(0));

    //compute Wave parameters
    ZFSFloat deltaS = -1.0;
    ZFSFloat cf = -1.0;
    if(m_Re<5000) {
      // use Smits formula for smaller Reynolds numbers
      deltaS = pow(m_Re, -7.0 / 8.0)*sqrt(2.0 / 0.024);
      cf = 0.024*pow(m_Re, -F1B4);
    } else {
      // Coles-Fernholz is better for larger Re
      deltaS = ((log(m_Re))/0.384 + 4.127)/m_Re;
      cf = 2.0*pow((log(m_Re)/0.384 + 4.127), -2.0);
    }
    const ZFSFloat uTau = sqrt(POW2(PV->UInfinity)*CV->rhoInfinity*cf*F1B2);

    m_waveLength = m_waveLengthPlus*deltaS;
    m_waveAmplitude = m_waveAmplitudePlus*deltaS;
    m_waveSpeedPlus = m_waveLengthPlus/m_waveTimePlus;
    m_waveSpeed = m_waveSpeedPlus*uTau;
    
    //assume equidistant grid in z-direction
    const ZFSFloat deltaZ = abs(m_coordinates[2][0]-m_coordinates[2][m_nPoints[2]*m_nPoints[1]]);
    m_waveCellsPerWaveLength = round(m_waveLength/deltaZ);

    zfs_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
    zfs_log << "Re: " << m_Re << " c_f: " << cf << " u_tau: " << uTau << " deltaZ: " << deltaZ << " viscousUnit delta_v: " << deltaS << endl;
    zfs_log << "wavelengthPlus: " << m_waveLengthPlus << " AmplitudePlus: " << m_waveAmplitudePlus << " SpeedPlus: " << m_waveSpeedPlus << endl;
    zfs_log << "Wavelength: " << m_waveLength << " Amplitude: " << m_waveAmplitude << " Speed: " << m_waveSpeed << endl;
    zfs_log << "Max up/down speed: " << m_waveSpeed*m_waveAmplitude << endl;
    if(m_waveSpeed*m_waveAmplitude > 0.3) {
      zfs_log << "WARNING: max up/down speed is above 0.3, beware of compressibility effects!" << endl;
      if(domainId()==0) {
        cout << "///////////////////////////////////////////////////////////////////////////" << endl;
        cout << "WARNING: max up/down speed is above 0.3, beware of compressibility effects!" << endl;
        cout << "///////////////////////////////////////////////////////////////////////////" << endl;
      }
    }
    zfs_log << "////////////////////////////////////////////////////////////////" << endl;
    fixTimeStepTravelingWave();
    break;
  }
  case 10:
    //travelling wave defined by non-plus units
  {
    m_travelingWave = true;
    m_waveTimeStepComputed = false;
    m_waveSpeed = 0.0;
    m_waveBeginTransition = 0.0;
    m_waveEndTransition = 0.0;
    m_waveCellsPerWaveLength = 1;
    m_waveNoStepsPerCell = 1;

    //time needs to be constant for traveling wave
    m_constantTimeStep = true;

    /*! \page propertyPage1
      \section waveLength
      <code>ZFSInt ZFSStrctrdBlck::m_waveLength </code>\n
      default = <code> 1.0 </code>\n \n
      Wavelength of the traveling wave.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>WAVE, MOVING, STRCTRD</i>
    */
    m_waveLength = 0.0;
    m_waveLength = *(ZFSContext::getProperty("waveLength", m_blockId, __CALLING_FUNCTION__, &m_waveLengthPlus)->asFloat(0));

    /*! \page propertyPage1
      \section waveAmplitude
      <code>ZFSInt ZFSStrctrdBlck::m_waveAmplitude </code>\n
      default = <code> 1.0 </code>\n \n
      Amplitude of the traveling wave.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>WAVE, MOVING, STRCTRD</i>
    */
    m_waveAmplitude = 0.0;
    m_waveAmplitude = *(ZFSContext::getProperty("waveAmplitude", m_blockId, __CALLING_FUNCTION__, &m_waveAmplitudePlus)->asFloat(0));

    /*! \page propertyPage1
      \section waveTime
      <code>ZFSInt ZFSStrctrdBlck::m_waveTime </code>\n
      default = <code> 1.0 </code>\n \n
      Period time of the traveling wave.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>WAVE, MOVING, STRCTRD</i>
    */
    m_waveTime = 0.0;
    m_waveTime = *(ZFSContext::getProperty("waveTime", m_blockId, __CALLING_FUNCTION__, &m_waveTimePlus)->asFloat(0));
    m_waveBeginTransition = *(ZFSContext::getProperty("waveBeginTransition", m_blockId, __CALLING_FUNCTION__, &m_waveBeginTransition)->asFloat(0));
    m_waveEndTransition = *(ZFSContext::getProperty("waveEndTransition", m_blockId, __CALLING_FUNCTION__, &m_waveEndTransition)->asFloat(0));

    m_waveSpeed = m_waveLength/m_waveTime;
    const ZFSFloat deltaZ = abs(m_coordinates[2][0]-m_coordinates[2][m_nPoints[2]*m_nPoints[1]]);
    m_waveCellsPerWaveLength = round(m_waveLength/deltaZ);

    zfs_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
    zfs_log << "Wavelength: " << m_waveLength << " Amplitude: " << m_waveAmplitude << " Speed: " << m_waveSpeed << endl;
    zfs_log << "Max up/down speed: " << m_waveSpeed*m_waveAmplitude << endl;
    if(m_waveSpeed*m_waveAmplitude > 0.3) {
      zfs_log << "WARNING: max up/down speed is above 0.3, beware of compressibility effects!" << endl;
      if(domainId()==0) {
        cout << "///////////////////////////////////////////////////////////////////////////" << endl;
        cout << "WARNING: max up/down speed is above 0.3, beware of compressibility effects!" << endl;
        cout << "///////////////////////////////////////////////////////////////////////////" << endl;
      }
    }
    zfs_log << "////////////////////////////////////////////////////////////////" << endl;
    fixTimeStepTravelingWave();
    break;
  }
  case 11:
    //traveling wave channel (Tomiyama & Fukagata 2013)
  {
    m_travelingWave = true;
    m_waveTimeStepComputed = false;
    m_waveSpeed = 0.0;
    m_waveLength = 0.0;
    m_waveAmplitude = 0.0;
    m_waveLengthPlus = 0.0;
    m_waveAmplitudePlus = 0.0;
    m_waveTimePlus = 0.0;
    m_waveCellsPerWaveLength = 1;
    m_waveNoStepsPerCell = 1;

    //time needs to be constant for traveling wave
    m_constantTimeStep = true;
    m_waveLengthPlus = *(ZFSContext::getProperty("waveLengthPlus", m_blockId, __CALLING_FUNCTION__, &m_waveLengthPlus)->asFloat(0));
    m_waveAmplitudePlus = *(ZFSContext::getProperty("waveAmplitudePlus", m_blockId, __CALLING_FUNCTION__, &m_waveAmplitudePlus)->asFloat(0));
    m_waveTimePlus = *(ZFSContext::getProperty("waveTimePlus", m_blockId, __CALLING_FUNCTION__, &m_waveTimePlus)->asFloat(0));

    //compute Wave parameters
    // const ZFSFloat cf = 0.008185; //0.024*pow(m_Re, -F1B4);
    // const ZFSFloat uTau = sqrt(POW2(PV->UInfinity)*CV->rhoInfinity*cf*F1B2);
    const ZFSFloat uTau= m_ReTau*m_Ma*sqrt(PV->TInfinity)/m_Re;
    const ZFSFloat cf = 2*POW2(uTau/PV->UInfinity);

    const ZFSFloat mu8 = zfsSUTHERLANDLAW(PV->TInfinity);
    m_waveLength = m_waveLengthPlus/(sqrt(cf/2.0)*PV->UInfinity*m_Re0*CV->rhoInfinity/mu8);
    m_waveAmplitude = m_waveAmplitudePlus/(sqrt(cf/2.0)*PV->UInfinity*m_Re0*CV->rhoInfinity/mu8);
    m_waveSpeed = (m_waveLengthPlus/m_waveTimePlus)*uTau;///PV->UInfinity;
    const ZFSFloat deltaZ = abs(m_coordinates[2][0]-m_coordinates[2][m_nPoints[2]*m_nPoints[1]]);
    m_waveCellsPerWaveLength = round(m_waveLength/deltaZ);

    zfs_log << "WaveLength: " << m_waveLength << " waveAmplitude: " << m_waveAmplitude << " waveSpeed: " << m_waveSpeed << endl;
    fixTimeStepTravelingWave();
    break;
  }
  default:
  {
    zfsTerm(1, __CALLING_FUNCTION__, "Grid Moving Method not implemented!");
  }
  }

  saveGrid();

  //now move the grid to the correct position
  if(m_restart) {
    if (m_movingGrid && !m_stgInitialStartup) {
      if(m_movingGridInitialStart) {
        //if this is an initial start of the
        //grid movement, just move to initial pos
        moveGrid(true,false);
        saveGrid();
      } else {
        //move to last pos before restart,
        //save and move to current pos again
        //this way the grid velocity is computed
        //correctly in the BC
        m_time -= m_timeStep;
        moveGrid(true,false);
        saveGrid();
        m_time += m_timeStep;
        moveGrid(true,false);
      }
    }
  }
}

void ZFSStrctrdBlck3D::saveGrid()  // saves the whole grid
{
  TRACE();
  ZFSId pointId=0;
  for(ZFSId k=0; k<m_nPoints[0]; ++k) {
    for(ZFSId j=0; j<m_nPoints[1]; ++j) {
      for(ZFSId i=0; i<m_nPoints[2]; ++i) {
        pointId=pointIndex(i,j,k);
        for( ZFSInt isd = xsd; isd < nDim; ++isd) {
          m_mgOldCoordinates[isd][pointId] = m_coordinates[isd][pointId];
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::saveCellJacobian()  //saves the old cell Jacobian
{
  TRACE();
  ZFSFloat *const RESTRICT oldCellJac= ALIGNED_MF(m_cells->oldCellJac);
  const ZFSFloat *const RESTRICT cellJac= ALIGNED_MF(m_cells->cellJac);
  for(ZFSId cellId =0; cellId<m_noStrctrdCells; cellId++) {
    oldCellJac[cellId] = cellJac[cellId];
  }
}


void ZFSStrctrdBlck3D::loadRestartBC2600() {
  if(m_bc2600IsActive && !m_bc2600InitialStartup) {
   
    if(domainId()==0) {
      cout << "Loading BC2600 values..." << endl;
    }
    if(m_bc2600){  //junoh   
    // ZFSInt nInputBlockCellsRANS[3]={0,0,0};   //junoh   
    // MPI_Allreduce(&m_nInputBlockCells[0], &nInputBlockCellsRANS[0], 1, MPI_INT, MPI_MIN, m_zfsStrctrdComm);
    // MPI_Allreduce(&m_nInputBlockCells[1], &nInputBlockCellsRANS[1], 1, MPI_INT, MPI_MIN, m_zfsStrctrdComm);
    // ZFSInt bcCellsRANS[3] = {nInputBlockCellsRANS[0], nInputBlockCellsRANS[1] ,m_noGhostLayers};
    ZFSInt bcCells[3] = {m_nInputBlockCells[0],m_nInputBlockCells[1],m_noGhostLayers};    
    ZFSInt noCellsBC = bcCells[0]*bcCells[1]*bcCells[2];
    // ZFSInt noCellsBCRANS = bcCellsRANS[0]*bcCellsRANS[1]*bcCellsRANS[2];   
    ZFSInt bcOffset[3] = {0,0,0};
    ZFSFloatScratchSpace tmpRestartVars(noCellsBC*m_maxNoVariables, __CALLING_FUNCTION__, "m_tmpRestartVars2600"); //junoh
    // cout<<"noCellsBC_RANS*m_maxNoVariables:"<<noCellsBCRANS*m_maxNoVariables<<"domainId:"<<domainId()<<endl;
    // cout<<"bcCells_RANS[0]:"<<bcCellsRANS[0]<<"bcCells_RANS[1]:"<<bcCellsRANS[1]<<"bcCells_RANS[2]"<<bcCellsRANS[2]<<endl;
    if(m_commBC2600MyRank==0) {// m_zonalRootRank
      stringstream restartFileName;
      ZFSInt restartFileId =-1;
      ZFSString restartFile = *(ZFSContext::getProperty("restartVariablesFileName", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL )->asString(0));
      restartFileName << outputDir() << restartFile;
      restartFileId = io_openfile("hdf5" ,(restartFileName.str()).c_str(),"collective", MPI_COMM_SELF);
      stringstream pathStr;
      pathStr << "/block" << m_inputBlockId << "/bc2600" << endl;
      const char* path = (pathStr.str()).c_str();

      for(ZFSId var=0; var<m_maxNoVariables; var++) { //junoh
        io_read_ddataset_part1d1(restartFileId, path, m_pvariableNames[var].c_str(), nDim, bcOffset, bcCells, &tmpRestartVars[var*noCellsBC]); //junoh
      }

      io_closefile(restartFileId);
        }
  
    MPI_Bcast(&tmpRestartVars[0], noCellsBC*m_maxNoVariables, MPI_DOUBLE, 0, *m_commBC2600); //m_commZonal[m_inputBlockId] m_zfsStrctrdComm
      

    if(m_commBC2600MyRank==0) {
      cout << "Loading BC2600 values... SUCCESSFUL!" << endl;
    }

    // if(m_bc2600) {   junoh
      ZFSId startGC[3] = {0,0,0};
      ZFSId endGC[3] = {0,0,0};

      if(m_bc2600noOffsetCells[1] == 0) { startGC[1] = m_noGhostLayers; }
      if(m_bc2600noOffsetCells[0] == 0) { startGC[0] = m_noGhostLayers; }
      if(m_bc2600noOffsetCells[1]+m_bc2600noActiveCells[1] == bcCells[1]) { endGC[1] = m_noGhostLayers; }//junoh
      if(m_bc2600noOffsetCells[0]+m_bc2600noActiveCells[0] == bcCells[0]) { endGC[0] = m_noGhostLayers; } //junoh

      for(ZFSInt i = 0; i< m_noGhostLayers; i++) {
        for(ZFSInt j = startGC[1]; j<m_bc2600noCells[1]-endGC[1]; j++) {
          for(ZFSInt k = startGC[0]; k<m_bc2600noCells[0]-endGC[0]; k++) {
            ZFSId cellId = cellIndex(i,j,k);
            ZFSId globalI = i;
            ZFSId globalJ = m_bc2600noOffsetCells[1]-m_noGhostLayers+j;
            ZFSId globalK = m_bc2600noOffsetCells[0]-m_noGhostLayers+k;
            ZFSId cellIdBC = globalI+(globalJ + globalK*bcCells[1])*bcCells[2];    //junoh

            //load values from restart field
            for(ZFSId var=0; var<m_maxNoVariables; var++) {
	      // m_bc2600Variables[var][cellId] = tmpRestartVars[var*noCellsBCRANS+cellIdBC]; //junoh
	      m_cells->pvariables[var][cellId] = tmpRestartVars[var*noCellsBC+cellIdBC]; //junoh

	    }
	    // m_cells->pvariables[PV->RANS_FIRST][cellId] = tmpRestartVars[PV->RANS_FIRST*noCellsBC+cellIdBC]; //junoh
          }
        }
      }


      //Fix diagonal cells at end of domain
      if(m_bc2600noOffsetCells[1] + m_bc2600noActiveCells[1] == m_nInputBlockCells[1]) { //junoh m_nInputBlockCells[1]
        for(ZFSInt i = 0; i<m_noGhostLayers; i++) {
          for(ZFSInt k = 0; k<m_bc2600noCells[0]; k++) {
            const ZFSId cellIdA2 = cellIndex(i,m_noGhostLayers+m_bc2600noActiveCells[1]-2,k);
            const ZFSId cellIdA1 = cellIndex(i,m_noGhostLayers+m_bc2600noActiveCells[1]-1,k);
            const ZFSId cellIdG1 = cellIndex(i,m_noGhostLayers+m_bc2600noActiveCells[1],  k);
            for(ZFSId var=0; var<m_maxNoVariables; var++) { //junoh
              const ZFSFloat distA1A2 = sqrt(POW2(m_cells->coordinates[0][cellIdA1]-m_cells->coordinates[0][cellIdA2])+
                                             POW2(m_cells->coordinates[1][cellIdA1]-m_cells->coordinates[1][cellIdA2])+
                                             POW2(m_cells->coordinates[2][cellIdA1]-m_cells->coordinates[2][cellIdA2]));
	    
	      // const ZFSFloat slope = (m_bc2600Variables[var][cellIdA1]- m_bc2600Variables[var][cellIdA2])/distA1A2; 
              const ZFSFloat slope = (m_cells->pvariables[var][cellIdA1]-m_cells->pvariables[var][cellIdA2])/distA1A2;


              const ZFSFloat distG1A1 = sqrt(POW2(m_cells->coordinates[0][cellIdG1]-m_cells->coordinates[0][cellIdA1])+
                                             POW2(m_cells->coordinates[1][cellIdG1]-m_cells->coordinates[1][cellIdA1])+
                                             POW2(m_cells->coordinates[2][cellIdG1]-m_cells->coordinates[2][cellIdA1]));
	      // m_bc2600Variables[var][cellIdG1] = m_bc2600Variables[var][cellIdA1] + distG1A1*slope;
	   
              m_cells->pvariables[var][cellIdG1] = m_cells->pvariables[var][cellIdA1] + distG1A1*slope;
           
	    }
          }
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::loadRestartBC2601() {
  //if we have the prescribing boundary,
  //put the values from the restart file into the ghostcells
  if(m_bc2601IsActive && !m_bc2601InitialStartup) {
    if(domainId()==0) {cout << "Loading restart values 2601" << endl;}
    ZFSInt bcCells[3] = {0,0,0};
    ZFSInt bcOffset[3] = {0,0,0};

    bcCells[0] = m_partition->inputBoxInfo[m_inputBlockId]->DirLast[0];
    bcCells[1] = m_noGhostLayers;
    bcCells[2] = m_partition->inputBoxInfo[m_inputBlockId]->DirLast[2];

    ZFSId noCellsBC = bcCells[0]*bcCells[1]*bcCells[2];
    ZFSFloatScratchSpace tmpRestartVars(noCellsBC*PV->noVariables, __CALLING_FUNCTION__, "m_tmpRestartVars2600");

    if(domainId()==0) {
      stringstream restartFileName;
      ZFSInt restartFileId =-1;
      ZFSString restartFile = *(ZFSContext::getProperty("restartVariablesFileName", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL )->asString(0));
      restartFileName << outputDir() << restartFile;

      restartFileId = io_openfile("hdf5" ,(restartFileName.str()).c_str(),"collective", MPI_COMM_SELF);
      stringstream pathStr;
      pathStr << "/block" << m_inputBlockId << "/bc2601" << endl;
      const char* path = (pathStr.str()).c_str();

      for(ZFSId var=0; var<PV->noVariables; var++) {
        io_read_ddataset_part1d1(restartFileId, path, m_pvariableNames[var].c_str(), nDim, bcOffset, bcCells, &tmpRestartVars[var*noCellsBC]);
      }
      io_closefile(restartFileId);
    }

    MPI_Bcast(&tmpRestartVars[0], noCellsBC*PV->noVariables, MPI_DOUBLE, 0, m_zfsStrctrdComm);

    if(m_bc2601) {
      ZFSId startGC[3] = {0,0,0};
      ZFSId endGC[3] = {0,0,0};

      if(m_nOffsetCells[2] == 0) { startGC[2] = m_noGhostLayers; }
      if(m_nOffsetCells[0] == 0) { startGC[0] = m_noGhostLayers; }
      if(m_nOffsetCells[2]+m_nActiveCells[2] == m_partition->inputBoxInfo[m_inputBlockId]->DirLast[2]) { endGC[2] = m_noGhostLayers; }
      if(m_nOffsetCells[0]+m_nActiveCells[0] == m_partition->inputBoxInfo[m_inputBlockId]->DirLast[0]) { endGC[0] = m_noGhostLayers; }

      for(ZFSInt i = startGC[2]; i<m_nCells[2]-endGC[2]; i++) {
        for(ZFSInt j = 0; j< m_noGhostLayers; j++) {
          for(ZFSInt k = startGC[0]; k<m_nCells[0]-endGC[0]; k++) {
            ZFSId cellId = cellIndex(i,j,k);
            ZFSId globalI = m_nOffsetCells[2]-m_noGhostLayers+i;
            ZFSId globalJ = j;
            ZFSId globalK = m_nOffsetCells[0]-m_noGhostLayers+k;
            ZFSId cellIdBC = globalI+(globalJ + globalK*m_noGhostLayers)*m_partition->inputBoxInfo[m_inputBlockId]->DirLast[2];

            //load values from restart field
            for(ZFSId var=0; var<PV->noVariables; var++) {
              m_cells->pvariables[var][cellId] = tmpRestartVars[var*noCellsBC+cellIdBC];
            }
          }
        }
      }

      //Fix diagonal cells at start of domain
      if(m_nOffsetCells[2] == 0) {
        for(ZFSInt j = 0; j<m_noGhostLayers; j++) {
          for(ZFSInt k = 0; k<m_nCells[0]; k++) {
            const ZFSId cellIdA2 = cellIndex(3,j,k);
            const ZFSId cellIdA1 = cellIndex(2,j,k);
            const ZFSId cellIdG1 = cellIndex(1,j,k);
            for(ZFSId var=0; var<PV->noVariables; var++) {
              const ZFSFloat slope = (m_cells->pvariables[var][cellIdA2]-m_cells->pvariables[var][cellIdA1])/
                (m_cells->coordinates[0][cellIdA2]-m_cells->coordinates[0][cellIdA1]);
              m_cells->pvariables[var][cellIdG1] = m_cells->pvariables[var][cellIdA1] + (m_cells->coordinates[0][cellIdG1]-m_cells->coordinates[0][cellIdA1])*slope;
            }
          }
        }
      }

      //Fix diagonal cells at end of domain
      if(m_nOffsetCells[2] + m_nActiveCells[2] == m_partition->inputBoxInfo[m_inputBlockId]->DirLast[2]) {
        for(ZFSInt j = 0; j<m_noGhostLayers; j++) {
          for(ZFSInt k = 0; k<m_nCells[0]; k++) {
            const ZFSId cellIdA2 = cellIndex(m_noGhostLayers+m_nActiveCells[2]-2,j,k);
            const ZFSId cellIdA1 = cellIndex(m_noGhostLayers+m_nActiveCells[2]-1,j,k);
            const ZFSId cellIdG1 = cellIndex(m_noGhostLayers+m_nActiveCells[2],  j,k);
            for(ZFSId var=0; var<PV->noVariables; var++) {
              const ZFSFloat slope = (m_cells->pvariables[var][cellIdA1]-m_cells->pvariables[var][cellIdA2])/
                (m_cells->coordinates[0][cellIdA1]-m_cells->coordinates[0][cellIdA2]);
              m_cells->pvariables[var][cellIdG1] = m_cells->pvariables[var][cellIdA1] + (m_cells->coordinates[0][cellIdG1]-m_cells->coordinates[0][cellIdA1])*slope;
            }
          }
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::loadRestartSTG(ZFSBool isPrimitiveOutput) {
  RECORD_TIMER_START(m_tloadRestartStg);
  if(m_stgIsActive) {
    if(!isPrimitiveOutput && domainId()==0) {
      cout << "Restart file has conservative variables, converting STG variables to primitive!" << endl;
    }

    // if(m_zonal){  //junoh
    //   m_stgNoVariables=6;
      
    // } 
    stringstream restartFileName;
    ZFSString restartFile = *(ZFSContext::getProperty("restartVariablesFileName", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL )->asString(0));
    restartFileName << outputDir() << restartFile;
    ZFSId restartFileId = io_openfile("hdf5" ,(restartFileName.str()).c_str(),"collective", m_zfsStrctrdComm);

    stringstream blockNumber;
    blockNumber << m_inputBlockId;
    ZFSString blockPathStr = "/block";
    blockPathStr += blockNumber.str();
    const char* blockPath = blockPathStr.c_str();

    ZFSInt start[ 3 ];
    ZFSInt Id = m_partition->outputBoxInfo[domainId()]->cpu;
    //the offset in the file
    for (ZFSInt i=0; i<nDim;i++) {
      //get offset in file
      start[i]=m_partition->outputBoxInfo[Id]->offset[i];
    }

    ZFSInt restartNoEddies = 0;
    io_read_iattribute1(restartFileId,"", "stgNRAN",&restartNoEddies);

    if(restartNoEddies != m_stgMaxNoEddies) {
      zfs_log << "STG: NRAN in restart file (" << restartNoEddies
              << ") not the same as given in property file ("
              << m_stgMaxNoEddies << "), creating new random distribution of eddies!" << endl;
      m_stgCreateNewEddies = true;
    } else {
      zfs_log << "STG: Reading in " << restartNoEddies << " eddies from restart file" << endl;
    }

    //if this is an initialStartup the new
    //eddies need to be created in any case
    if(m_stgInitialStartup) {
      m_stgCreateNewEddies = true;
      //also load nu_t into fq field
      io_read_ddataset_part1d1(restartFileId, blockPath, FQ->fqNames[FQ->NU_T].c_str(), nDim, start, m_nActiveCells, m_cells->fq[FQ->NU_T]);
      FQ->loadedFromRestartFile[FQ->NU_T] = true;
    } else {
      // has to be set manually for restart from RANS profile
      ZFSInt ninmax = int(m_stgNoEddieProperties*m_stgMaxNoEddies);
      ZFSInt *VBStart = new ZFSInt;
      *VBStart = 0;


      //////////////////////////////////////////////////
      ////////////// LOAD EDDIES ///////////////////////
      //////////////////////////////////////////////////
      if(m_stgCreateNewEddies) {
        if(domainId()==0) {
          cout << "NRAN in property differs from NRAN in restart file"
               << " NOT READING EDDIES FROM RESTART!" << endl;
        }
      } else {
        ZFSString stgGlobalPathStr = "stgGlobal";
	
        if(io_checkObj(restartFileId, stgGlobalPathStr.c_str(), "FQeddies")) {
          if(domainId()==0) {cout << "FQeddies field is at new position /stgGlobal/FQeddies" << endl;}
        } else {
          if(domainId()==0) {cout << "FQeddies field is NOT at new position! Using old path within block..." << endl;}
          stgGlobalPathStr = blockPathStr;
          if(domainId()==0) {cout << "Loading FQeddies from path " << stgGlobalPathStr << endl;}
        }

	RECORD_TIMER_START(m_tloadRestartStgEddies);
        const char* stgGlobalPath = stgGlobalPathStr.c_str();
        //do this only serially
	if(globalDomainId()==0){ cout << "Loading STG Eddies..." << endl;}
        if(globalDomainId()==0){
	  ZFSInt fid=io_openfile("hdf5" ,(restartFileName.str()).c_str(),"collective", MPI_COMM_SELF);
          io_read_ddataset_part1d1(fid, stgGlobalPath, "FQeddies", 1, VBStart, &ninmax, m_stgEddies[0]);
	  io_closefile(fid);
        }

        MPI_Bcast(m_stgEddies[0], ninmax, MPI_DOUBLE, 0, m_zfsStrctrdComm);
	RECORD_TIMER_STOP(m_tloadRestartStgEddies);
	if(globalDomainId()==0){ cout << "Loading STG Eddies... SUCCESSFUL!" << endl;}
      }


      //////////////////////////////////////////////////
      ////////////// LOAD STG VARIABLES ////////////////
      //////////////////////////////////////////////////

      if(m_stgLocal) {
        ZFSInt bcCells[3] = {0,0,0};
        ZFSInt bcOffset[3] = {0,0,0};
       
        bcCells[0] = m_partition->inputBoxInfo[m_inputBlockId]->DirLast[0];
        bcCells[1] = m_partition->inputBoxInfo[m_inputBlockId]->DirLast[1];
        bcCells[2] = 3;
        const ZFSId noCellsBC = bcCells[0]*bcCells[1]*bcCells[2];
        ZFSFloatScratchSpace tmpRestartVars(noCellsBC*m_stgNoVariables, __CALLING_FUNCTION__, "tmpRestartVars");
        RECORD_TIMER_START(m_tloadRestartStgRead);
        if(m_commStgMyRank==0){ cout << "Loading STG Datasets..." << endl;}
        if(m_commStgMyRank==0) {
          ZFSInt fid=io_openfile("hdf5" ,(restartFileName.str()).c_str(),"collective", MPI_COMM_SELF);
          for(ZFSId var = 0; var < m_stgNoVariables; var++) {
            stringstream fieldName;
            stringstream stgPath;
            stgPath << blockPathStr << "/stg";
            fieldName << "stgFQ" <<  var;
            io_read_ddataset_part1d1(fid, (stgPath.str()).c_str(), (fieldName.str()).c_str(),  nDim, bcOffset, bcCells, &tmpRestartVars[var*noCellsBC]);
          }
          io_closefile(fid);
	}
        RECORD_TIMER_STOP(m_tloadRestartStgRead);
       
        RECORD_TIMER_START(m_tloadRestartStgBcast);
        MPI_Bcast(&tmpRestartVars[0], noCellsBC*m_stgNoVariables, MPI_DOUBLE, 0, *m_commStg);
        RECORD_TIMER_STOP(m_tloadRestartStgBcast);
        if(m_commStgMyRank==0){ cout << "Loading STG Datasets... SUCCESSFUL!" << endl;}

        //////////////////////////////////////////////////
        ////////// DISTRIBUTE STG VARIABLES //////////////
        //////////////////////////////////////////////////

        ZFSId startGC[3] = {0,0,0};
        ZFSId endGC[3] = {0,0,0};
        if(m_nOffsetCells[1] == 0) { startGC[1] = m_noGhostLayers; }
        if(m_nOffsetCells[0] == 0) { startGC[0] = m_noGhostLayers; }
        if(m_nOffsetCells[1]+m_nActiveCells[1] == m_partition->inputBoxInfo[m_inputBlockId]->DirLast[1]) { endGC[1] = m_noGhostLayers; }
        if(m_nOffsetCells[0]+m_nActiveCells[0] == m_partition->inputBoxInfo[m_inputBlockId]->DirLast[0]) { endGC[0] = m_noGhostLayers; }


        for(ZFSInt k = startGC[0]; k<m_nCells[0]-endGC[0]; k++) {
          for(ZFSInt j = startGC[1]; j<m_nCells[1]-endGC[1]; j++) {
            for(ZFSInt i = 0; i< 3; i++) {
              ZFSId cellId = cellIndex(i,j,k);
              ZFSId globalI = i;
              ZFSId globalJ = m_nOffsetCells[1]-m_noGhostLayers+j;
              ZFSId globalK = m_nOffsetCells[0]-m_noGhostLayers+k;
              ZFSId cellIdBCGlobal = globalI+(globalJ + globalK*m_partition->inputBoxInfo[m_inputBlockId]->DirLast[1])*3;
              ZFSId cellIdBC = i + (j+k*m_nCells[1])*3;

              //load values from restart field
              for(ZFSId var = 0; var < m_stgNoVariables; var++) {
                m_cells->stg_fq[var][cellIdBC] = tmpRestartVars[var*noCellsBC+cellIdBCGlobal];
              }

              if(!isPrimitiveOutput) {
                const ZFSFloat rho = m_cells->stg_fq[0][cellIdBC];
                const ZFSFloat rhoU = m_cells->stg_fq[1][cellIdBC];
                const ZFSFloat rhoV = m_cells->stg_fq[2][cellIdBC];
                const ZFSFloat rhoW = m_cells->stg_fq[3][cellIdBC];
                const ZFSFloat rhoE = m_cells->stg_fq[4][cellIdBC];
                
                const ZFSFloat gammaMinusOne = m_gamma-1.0;
                const ZFSFloat u = rhoU/rho;
                const ZFSFloat v = rhoV/rho;
                const ZFSFloat w = rhoW/rho;
                const ZFSFloat p = gammaMinusOne*(rhoE-F1B2*rho*(POW2(u)+POW2(v)+POW2(w)));

                m_cells->stg_fq[PV->RHO][cellIdBC] = rho;
                m_cells->stg_fq[PV->U][cellIdBC] = u;
                m_cells->stg_fq[PV->V][cellIdBC] = v;
                m_cells->stg_fq[PV->W][cellIdBC] = w;
                m_cells->stg_fq[PV->P][cellIdBC] = p;
              }

              if(i<2) {
                m_cells->pvariables[PV->RHO][cellId] = m_cells->stg_fq[PV->RHO][cellIdBC];
                m_cells->pvariables[PV->U][cellId] = m_cells->stg_fq[PV->U][cellIdBC];
                m_cells->pvariables[PV->V][cellId] = m_cells->stg_fq[PV->V][cellIdBC];
                m_cells->pvariables[PV->W][cellId] = m_cells->stg_fq[PV->W][cellIdBC];
                m_cells->pvariables[PV->P][cellId] = m_cells->stg_fq[PV->P][cellIdBC];
              }
            }
          }
        }
      }
    }

    io_closefile(restartFileId);
  }
  RECORD_TIMER_STOP(m_tloadRestartStg);
}


void ZFSStrctrdBlck3D::computeVorticity(){
  TRACE();
  ZFSFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  ZFSFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  ZFSFloat* const RESTRICT w = &m_cells->pvariables[PV->W][0];
  ZFSFloat* const RESTRICT vortx=&m_cells->fq[FQ->VORTX][0];
  ZFSFloat* const RESTRICT vorty=&m_cells->fq[FQ->VORTY][0];
  ZFSFloat* const RESTRICT vortz=&m_cells->fq[FQ->VORTZ][0];
  ZFSFloat* const RESTRICT jac=&m_cells->cellJac[0];


  for(ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++){
    for(ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++){
      for(ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++){

        const ZFSId IJK   = i+(j+k*m_nCells[1])*m_nCells[2];
        const ZFSId IPJK  = (i+1)+(j+k*m_nCells[1])*m_nCells[2];
        const ZFSId IMJK  = (i-1)+(j+k*m_nCells[1])*m_nCells[2];
        const ZFSId IJPK  = i+((j+1)+k*m_nCells[1])*m_nCells[2];
        const ZFSId IJMK  = i+((j-1)+k*m_nCells[1])*m_nCells[2];
        const ZFSId IJKP  = i+(j+(k+1)*m_nCells[1])*m_nCells[2];
        const ZFSId IJKM  = i+(j+(k-1)*m_nCells[1])*m_nCells[2];

        const ZFSFloat dudxi=u[IPJK]-u[IMJK];
        const ZFSFloat dudet=u[IJPK]-u[IJMK];
        const ZFSFloat dudze=u[IJKP]-u[IJKM];

        const ZFSFloat dvdxi=v[IPJK]-v[IMJK];
        const ZFSFloat dvdet=v[IJPK]-v[IJMK];
        const ZFSFloat dvdze=v[IJKP]-v[IJKM];

        const ZFSFloat dwdxi=w[IPJK]-w[IMJK];
        const ZFSFloat dwdet=w[IJPK]-w[IJMK];
        const ZFSFloat dwdze=w[IJKP]-w[IJKM];

        const ZFSFloat dvdz=  dvdxi * m_cells->cellMetrics[IJK][ xsd * 3 + zsd ] +
          dvdet * m_cells->cellMetrics[IJK][ ysd * 3 + zsd ] +
          dvdze * m_cells->cellMetrics[IJK][ zsd * 3 + zsd ];
        const ZFSFloat dwdy=  dwdxi * m_cells->cellMetrics[IJK][ xsd * 3 + ysd ] +
          dwdet * m_cells->cellMetrics[IJK][ ysd * 3 + ysd ] +
          dwdze * m_cells->cellMetrics[IJK][ zsd * 3 + ysd ];

        const ZFSFloat dudz=  dudxi * m_cells->cellMetrics[IJK][ xsd * 3 + zsd ] +
          dudet * m_cells->cellMetrics[IJK][ ysd * 3 + zsd ] +
          dudze * m_cells->cellMetrics[IJK][ zsd * 3 + zsd ];
        const ZFSFloat dwdx=  dwdxi * m_cells->cellMetrics[IJK][ xsd * 3 + xsd ] +
          dwdet * m_cells->cellMetrics[IJK][ ysd * 3 + xsd ] +
          dwdze * m_cells->cellMetrics[IJK][ zsd * 3 + xsd ];

        const ZFSFloat dvdx=  dvdxi * m_cells->cellMetrics[IJK][ xsd * 3 + xsd ] +
          dvdet * m_cells->cellMetrics[IJK][ ysd * 3 + xsd ] +
          dvdze * m_cells->cellMetrics[IJK][ zsd * 3 + xsd ];
        const ZFSFloat dudy=  dudxi * m_cells->cellMetrics[IJK][ xsd * 3 + ysd ] +
          dudet * m_cells->cellMetrics[IJK][ ysd * 3 + ysd ] +
          dudze * m_cells->cellMetrics[IJK][ zsd * 3 + ysd ];

        vortx[IJK] = F1B2*(dwdy-dvdz)/ jac[IJK];
        vorty[IJK] = F1B2*(dudz-dwdx)/ jac[IJK];
        vortz[IJK] = F1B2*(dvdx-dudy)/ jac[IJK];
      }
    }
  }
}


/**
 *     function to compute the lambda_2 criterion 
 *     /author Pascal Meysonnat
 *     /date   01.01.1010
 *
 */
void ZFSStrctrdBlck3D::computeLambda2Criterion()
{
  TRACE();
  ZFSFloatScratchSpace J(nDim,nDim,__CALLING_FUNCTION__,"J");
  ZFSFloat d[3] = {F0,F0,F0};
  ZFSFloat e[3] = {F0,F0,F0};
  J.fill(F0);

  //ZFSId IMJK, IJMK, IMJMK, IJKM, IJMKM, IMJKM;
  ZFSFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  ZFSFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  ZFSFloat* const RESTRICT w = &m_cells->pvariables[PV->W][0];
  ZFSFloat* const RESTRICT lambda2=&m_cells->fq[FQ->LAMBDA2][0];

  //compute the lambda2 criterion.
  for(ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers+1; k++) {
    for(ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers+1; j++) {
      for(ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers+1; i++) {
        const ZFSId cellId = cellIndex(i,j,k);
        const ZFSId IMJK=cellIndex(i-1,j,k);
        const ZFSId IPJK=cellIndex(i+1,j,k);
        const ZFSId IJMK=cellIndex(i,j-1,k);
        const ZFSId IJPK=cellIndex(i,j+1,k);
        const ZFSId IJKM=cellIndex(i,j,k-1);
        const ZFSId IJKP=cellIndex(i,j,k+1);
        const ZFSFloat FcellJac=F1/m_cells->cellJac[cellId];

        const ZFSFloat dudx=FcellJac*(m_cells->cellMetrics[cellId][0]*(u[IPJK]-u[IMJK])+
                                      m_cells->cellMetrics[cellId][3]*(u[IJPK]-u[IJMK])+
                                      m_cells->cellMetrics[cellId][6]*(u[IJKP]-u[IJKM]));
        const ZFSFloat dudy=FcellJac*(m_cells->cellMetrics[cellId][1]*(u[IPJK]-u[IMJK])+
                                      m_cells->cellMetrics[cellId][4]*(u[IJPK]-u[IJMK])+
                                      m_cells->cellMetrics[cellId][7]*(u[IJKP]-u[IJKM]));
        const ZFSFloat dudz=FcellJac*(m_cells->cellMetrics[cellId][2]*(u[IPJK]-u[IMJK])+
                                      m_cells->cellMetrics[cellId][5]*(u[IJPK]-u[IJMK])+
                                      m_cells->cellMetrics[cellId][8]*(u[IJKP]-u[IJKM]));

        const ZFSFloat dvdx=FcellJac*(m_cells->cellMetrics[cellId][0]*(v[IPJK]-v[IMJK])+
                                      m_cells->cellMetrics[cellId][3]*(v[IJPK]-v[IJMK])+
                                      m_cells->cellMetrics[cellId][6]*(v[IJKP]-v[IJKM]));
        const ZFSFloat dvdy=FcellJac*(m_cells->cellMetrics[cellId][1]*(v[IPJK]-v[IMJK])+
                                      m_cells->cellMetrics[cellId][4]*(v[IJPK]-v[IJMK])+
                                      m_cells->cellMetrics[cellId][7]*(v[IJKP]-v[IJKM]));
        const ZFSFloat dvdz=FcellJac*(m_cells->cellMetrics[cellId][2]*(v[IPJK]-v[IMJK])+
                                      m_cells->cellMetrics[cellId][5]*(v[IJPK]-v[IJMK])+
                                      m_cells->cellMetrics[cellId][8]*(v[IJKP]-v[IJKM]));

        const ZFSFloat dwdx=FcellJac*(m_cells->cellMetrics[cellId][0]*(w[IPJK]-w[IMJK])+
                                      m_cells->cellMetrics[cellId][3]*(w[IJPK]-w[IJMK])+
                                      m_cells->cellMetrics[cellId][6]*(w[IJKP]-w[IJKM]));
        const ZFSFloat dwdy=FcellJac*(m_cells->cellMetrics[cellId][1]*(w[IPJK]-w[IMJK])+
                                      m_cells->cellMetrics[cellId][4]*(w[IJPK]-w[IJMK])+
                                      m_cells->cellMetrics[cellId][7]*(w[IJKP]-w[IJKM]));
        const ZFSFloat dwdz=FcellJac*(m_cells->cellMetrics[cellId][2]*(w[IPJK]-w[IMJK])+
                                      m_cells->cellMetrics[cellId][5]*(w[IJPK]-w[IJMK])+
                                      m_cells->cellMetrics[cellId][8]*(w[IJKP]-w[IJKM]));

        //Compute the matrix (S^2+Omega^2)/2
        J(0,0)=POW2(dudx)+dvdx*dudy+dwdx*dudz;
        J(0,1)=F1B2*(dudx*dudy+dudx*dvdx+dvdy*dudy+dvdy*dvdx+dwdy*dudz+dvdz*dwdx);
        J(0,2)=F1B2*(dudx*dudz+dudx*dwdx+dudy*dvdz+dvdx*dwdy+dwdz*dudz+dwdz*dwdx);
        J(1,0)=J(0,1);
        J(1,1)=POW2(dvdy)+dvdx*dudy+dvdz*dwdy;
        J(1,2)=F1B2*(dudz*dvdx+dwdx*dudy+dvdy*dvdz+dvdy*dwdy+dwdz*dvdz+dwdz*dwdy);
        J(2,0)=J(0,2);
        J(2,1)=J(1,2);
        J(2,2)=POW2(dwdz)+dwdx*dudz+dvdz*dwdy;

        //perform householder tridiagonalization of
        //symmetric real matrix
        tred2(J, nDim, d, e);
        //compute eigenvalues
        tqli2(d,e,nDim);
        //sort eigenvalues
        insertSort(nDim,d);

        lambda2[cellId]=d[1];
      }
    }
  }
}


/**
 *     call function to compute the auxillary data (lift, cf, cp ...)
 *     /author Pascal Meysonnat
 *     /date   01.01.1010
 *
 */
void ZFSStrctrdBlck3D::computeAuxData(){
  TRACE();
  m_strctrdBndryCnd->computeAuxData();
}


/**
 *     function to save the power of actuation to file (auxData)
 *     /author Pascal Meysonnat
 *     /date   01.01.1010
 *
 */
void ZFSStrctrdBlck3D::savePowerCoefficient(ZFSInt fileId){
  TRACE();
  ZFSInt count =0;
  for(ZFSUint i=0; i<m_windowInfo->globalStrctrdBndryCndMaps.size(); ++i){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_windowInfo->globalStrctrdBndryCndMaps[i]->BC)/1000.0);
    if(firstDigit==1){
      count ++;
    }
  }
  ZFSInt noWalls= count;
  ZFSFloatScratchSpace cPower(noWalls, 3, __CALLING_FUNCTION__, "Ptot");
  ZFSFloatScratchSpace cPowerVisc(noWalls, 3, __CALLING_FUNCTION__, "Pvisc");
  ZFSFloatScratchSpace cPowerPres(noWalls, 3, __CALLING_FUNCTION__, "Ppres");

  for(ZFSId i=0; i<noWalls; ++i){
    for(ZFSId j=0; j<3; j++){
      cPower(i,j)=m_strctrdBndryCnd->m_Powerp[i*3+j] + m_strctrdBndryCnd->m_Powerv[i*3+j];
      cPowerVisc(i,j)=m_strctrdBndryCnd->m_Powerv[i*3+j];
      cPowerPres(i,j)=m_strctrdBndryCnd->m_Powerp[i*3+j];
    }
  }

  count=0;
  for(ZFSUint i=0; i<m_windowInfo->globalStrctrdBndryCndMaps.size(); ++i){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_windowInfo->globalStrctrdBndryCndMaps[i]->BC)/1000.0);
    if(firstDigit==1){
      stringstream datasetname;
      datasetname<< m_windowInfo->globalStrctrdBndryCndMaps[i]->Id2;
      ZFSString pathName= "window" + datasetname.str() ;

      io_create_dattribute(fileId,pathName.c_str(),"Ptot",1);
      io_create_dattribute(fileId,pathName.c_str(),"Pvisc",1);
      io_create_dattribute(fileId,pathName.c_str(),"Ppres",1);

      io_write_dattribute1(fileId,pathName.c_str(),"Ptot",1,  &cPower(count,1));
      io_write_dattribute1(fileId,pathName.c_str(),"Pvisc",1, &cPowerVisc(count,1));
      io_write_dattribute1(fileId,pathName.c_str(),"Ppres",1, &cPowerPres(count,1));

      count ++;
    }
  }
}


/**
 *     function to save the lift coefficient to file (auxData)
 *     /author Marian Albers/Pascal Meysonnat
 *     /date   01.01.1010
 *
 */

void ZFSStrctrdBlck3D::saveLiftCoefficient(ZFSInt fileId){
  TRACE();
  ZFSInt count =0;
  for(ZFSUint i=0; i<m_windowInfo->globalStrctrdBndryCndMaps.size(); ++i){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_windowInfo->globalStrctrdBndryCndMaps[i]->BC)/1000.0);
    if(firstDigit==1){
      count ++;
    }
  }
  ZFSInt noWalls= count;
  ZFSFloatScratchSpace cL(noWalls, 3, __CALLING_FUNCTION__, "cL");
  ZFSFloatScratchSpace cLv(noWalls, 3, __CALLING_FUNCTION__, "cLv");
  ZFSFloatScratchSpace cLp(noWalls, 3, __CALLING_FUNCTION__, "cLp");

  for(ZFSId i=0; i<noWalls; ++i){
    for(ZFSId j=0; j<3; j++){
      cL(i,j)=m_strctrdBndryCnd->m_cLv[i*3+j] + m_strctrdBndryCnd->m_cLp[i*3+j];
      cLv(i,j)=m_strctrdBndryCnd->m_cLv[i*3+j];
      cLp(i,j)=m_strctrdBndryCnd->m_cLp[i*3+j];
    }
  }

  count=0;
  for(ZFSUint i=0; i<m_windowInfo->globalStrctrdBndryCndMaps.size(); ++i){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_windowInfo->globalStrctrdBndryCndMaps[i]->BC)/1000.0);
    if(firstDigit==1){
      stringstream datasetname;
      datasetname<< m_windowInfo->globalStrctrdBndryCndMaps[i]->Id2;
      ZFSString pathName= "window" + datasetname.str() ;

      io_create_dattribute(fileId,pathName.c_str(),"Cly",1);
      io_create_dattribute(fileId,pathName.c_str(),"Clyv",1);
      io_create_dattribute(fileId,pathName.c_str(),"Clyp",1);

      io_write_dattribute1(fileId,pathName.c_str(),"Cly",1,  &cL(count,1));
      io_write_dattribute1(fileId,pathName.c_str(),"Clyv",1, &cLv(count,1));
      io_write_dattribute1(fileId,pathName.c_str(),"Clyp",1, &cLp(count,1));

      count ++;
    }
  }
}


/**
 *     function to save the drag coefficient to file (auxData)
 *     /author Marian Albers Pascal Meysonnat
 *     /date   01.01.1010
 *
 */
void ZFSStrctrdBlck3D::saveDragCoefficient(ZFSInt fileId){
  TRACE();
  ZFSInt count =0;
  for(ZFSUint i=0; i<m_windowInfo->globalStrctrdBndryCndMaps.size(); ++i){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_windowInfo->globalStrctrdBndryCndMaps[i]->BC)/1000.0);
    if(firstDigit==1){
      count ++;
    }
  }
  ZFSInt noWalls= count;
  ZFSFloatScratchSpace cD(noWalls, 3, __CALLING_FUNCTION__, "cD");
  ZFSFloatScratchSpace cDv(noWalls, 3, __CALLING_FUNCTION__, "cDv");
  ZFSFloatScratchSpace cDp(noWalls, 3, __CALLING_FUNCTION__, "cDp");

  for(ZFSId i=0; i<noWalls; ++i){
    for(ZFSId j=0; j<3; j++){
      cD(i,j)=m_strctrdBndryCnd->m_cDv[i*3+j] + m_strctrdBndryCnd->m_cDp[i*3+j];
      cDv(i,j)=m_strctrdBndryCnd->m_cDv[i*3+j];
      cDp(i,j)=m_strctrdBndryCnd->m_cDp[i*3+j];
    }
  }

  count=0;
  for(ZFSUint i=0; i<m_windowInfo->globalStrctrdBndryCndMaps.size(); ++i){
    ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_windowInfo->globalStrctrdBndryCndMaps[i]->BC)/1000.0);
    if(firstDigit==1){
      stringstream datasetname;
      datasetname<< m_windowInfo->globalStrctrdBndryCndMaps[i]->Id2;
      ZFSString pathName= "window" + datasetname.str() ;

      io_create_dattribute(fileId,pathName.c_str(),"Cdx",1);
      io_create_dattribute(fileId,pathName.c_str(),"Cdxv",1);
      io_create_dattribute(fileId,pathName.c_str(),"Cdxp",1);

      io_write_dattribute1(fileId,pathName.c_str(),"Cdx",1, &cD(count,0));
      io_write_dattribute1(fileId,pathName.c_str(),"Cdxv",1, &cDv(count,0));
      io_write_dattribute1(fileId,pathName.c_str(),"Cdxp",1, &cDp(count,0));

      count ++;
    }
  }
}


/**
 *     function to save the lift/drag coefficient and power to file (ascii)
 *     /author Marian Albers
 *     /date   01.01.1010
 *
 */
void ZFSStrctrdBlck3D::saveLiftDragToAsciiFile(){
  TRACE();
  //compute the data
  if(m_bCl && m_bCd){
    m_strctrdBndryCnd->computeAuxDataRoot();

    //write all the data to files if cpuRank==0
    if(domainId()==0){
      //we need to decide how many files to open
      ZFSInt count =0;
      for(ZFSUint i=0; i<m_windowInfo->globalStrctrdBndryCndMaps.size(); ++i){
	ZFSInt firstDigit=(ZFSInt)(((ZFSFloat) m_windowInfo->globalStrctrdBndryCndMaps[i]->BC)/1000.0);
	if(firstDigit==1){
	  count ++;
	}
      }//--> number of walls
      ZFSInt noWalls= count;
      
      //collect the data
      ZFSFloatScratchSpace cD(noWalls, 3, __CALLING_FUNCTION__, "cD");
      ZFSFloatScratchSpace cDv(noWalls, 3, __CALLING_FUNCTION__, "cDv");
      ZFSFloatScratchSpace cDp(noWalls, 3, __CALLING_FUNCTION__, "cDp");
      //drag
      for(ZFSId i=0; i<noWalls; ++i){
	for(ZFSId j=0; j<3; j++){
	  cD(i,j)=m_strctrdBndryCnd->m_cDv[i*3+j] + m_strctrdBndryCnd->m_cDp[i*3+j];
	  cDv(i,j)=m_strctrdBndryCnd->m_cDv[i*3+j];
	  cDp(i,j)=m_strctrdBndryCnd->m_cDp[i*3+j];
	}
      }
      //lift
      ZFSFloatScratchSpace cL(noWalls, 3, __CALLING_FUNCTION__, "cL");
      ZFSFloatScratchSpace cLv(noWalls, 3, __CALLING_FUNCTION__, "cLv");
      ZFSFloatScratchSpace cLp(noWalls, 3, __CALLING_FUNCTION__, "cLp");
      
      for(ZFSId i=0; i<noWalls; ++i){
	for(ZFSId j=0; j<3; j++){
	  cL(i,j)=m_strctrdBndryCnd->m_cLv[i*3+j] + m_strctrdBndryCnd->m_cLp[i*3+j];
	  cLv(i,j)=m_strctrdBndryCnd->m_cLv[i*3+j];
	  cLp(i,j)=m_strctrdBndryCnd->m_cLp[i*3+j];
	}
      }

      //area
      ZFSFloatScratchSpace cArea(noWalls, __CALLING_FUNCTION__, "cArea");
      for(ZFSId i=0; i<noWalls; ++i){
	cArea(i) = m_strctrdBndryCnd->m_cArea[i];
      }
      
      //power consumption
      ZFSFloatScratchSpace cPower(noWalls, 3, __CALLING_FUNCTION__, "Ptot");
      ZFSFloatScratchSpace cPowerVisc(noWalls, 3, __CALLING_FUNCTION__, "Pvisc");
      ZFSFloatScratchSpace cPowerPres(noWalls, 3, __CALLING_FUNCTION__, "Ppres");
      
      if(m_bPower){
        for(ZFSId i=0; i<noWalls; ++i){
          for(ZFSId j=0; j<3; j++){
            cPower(i,j)=m_strctrdBndryCnd->m_Powerp[i*3+j] + m_strctrdBndryCnd->m_Powerv[i*3+j];
            cPowerVisc(i,j)=m_strctrdBndryCnd->m_Powerv[i*3+j];
            cPowerPres(i,j)=m_strctrdBndryCnd->m_Powerp[i*3+j];
          }
        }
      }

      for(ZFSId i=0; i<noWalls; ++i){
        stringstream iWall;
        iWall << i;
        ZFSString filename = "./forces." + iWall.str() + ".dat";
        FILE* f_forces;
        f_forces = fopen(filename.c_str(), "a+");
        fprintf(f_forces, "%d", globalTimeStep);//#1
        fprintf(f_forces, " %f", m_time);//#2
        fprintf(f_forces, " %f", m_physicalTime);//#3
        fprintf(f_forces, " %f", cL(i,1));//#4
        fprintf(f_forces, " %f", cLv(i,1));//#5
        fprintf(f_forces, " %f", cLp(i,1));//#6
        fprintf(f_forces, " %f", cD(i,0));//#7
        fprintf(f_forces, " %f", cDv(i,0));//#8
        fprintf(f_forces, " %f", cDp(i,0));//#9
        fprintf(f_forces, " %f", cArea(i));//#10
        if(m_bPower){
          fprintf(f_forces, " %f", cPower(i,0));//#11
          fprintf(f_forces, " %f", cPowerVisc(i,0));//#12
          fprintf(f_forces, " %f", cPowerPres(i,0));//#13
          fprintf(f_forces, " %f", cPower(i,1));//#14
          fprintf(f_forces, " %f", cPowerVisc(i,1));//#15
          fprintf(f_forces, " %f", cPowerPres(i,1));//#16
          fprintf(f_forces, " %f", cPower(i,2));//#17
          fprintf(f_forces, " %f", cPowerVisc(i,2));//#18
          fprintf(f_forces, " %f", cPowerPres(i,2));//#19
        }
        fprintf(f_forces, "\n");
        fclose(f_forces);
      }
    }
  }
}

void ZFSStrctrdBlck3D::saveOutputLines() {
  TRACE();
  ////////////////////////////////
  ////////// INTERPOLATION ///////
  ////////////////////////////////

  //if it a moving grid the interpolation
  //coefficients need to be computed every time
  if(m_movingGrid) {
    for(ZFSInt pointId=0; pointId < m_noFieldPointsTotal; pointId++){
      m_hasPartnerLocal[pointId] = 0;
      m_hasPartnerGlobal[pointId] = 0;
      for(ZFSInt var=0; var<PV->noVariables; var++){
        m_interpolatedVarsLocal[var][pointId] = F0;
        m_interpolatedVarsGlobal[var][pointId] = F0;
      }
    }

    m_pointInterpolation = new ZFSStrctrdInterpolation<3>(m_nCells, m_cells->coordinates, m_cells->pvariables, m_zfsStrctrdComm);
    //allocate domains points of the lines
    m_pointInterpolation->prepareInterpolation(m_noFieldPointsTotal, m_pointCoordinates, m_hasPartnerLocal);
    MPI_Allreduce(m_hasPartnerLocal, m_hasPartnerGlobal, m_noFieldPointsTotal, MPI_INT, MPI_SUM, m_zfsStrctrdComm);
  }

  // calculation of interpolated variables only for the points in domain
  for(ZFSInt pointId=0;pointId < m_noFieldPointsTotal; pointId++){
    if(m_hasPartnerLocal[pointId]){
      for(ZFSInt var=0; var< PV->noVariables; var++){
        m_interpolatedVarsLocal[var][pointId]=m_pointInterpolation->getInterpolatedVariable(pointId,var);
      }
    }
  }

  MPI_Allreduce(&m_interpolatedVarsLocal[0][0], &m_interpolatedVarsGlobal[0][0], m_noFieldPointsTotal*PV->noVariables, MPI_DOUBLE, MPI_SUM, m_zfsStrctrdComm);

  // calculation of right value, if a point is assigned to more than one domain
  for(ZFSInt pointId=0; pointId < m_noFieldPointsTotal; pointId++){
    if(m_hasPartnerGlobal[pointId]>1 ){
      for(ZFSInt var=0; var<PV->noVariables; var++){
        m_interpolatedVarsGlobal[var][pointId]=m_interpolatedVarsGlobal[var][pointId]/ (ZFSFloat)m_hasPartnerGlobal[pointId];
      }
    }
  }

  if(m_movingGrid) {
    delete m_pointInterpolation;
    m_pointInterpolation = NULL;
  }

  ////////////////////////////////
  ////////// WRITE HDF5 //////////
  ////////////////////////////////

  ZFSChar gridFile[25];
  stringstream fileName;
  stringstream GridFileName;
  ZFSString tempG;
  fileName << m_lineOutputDir << "fieldOutput" << globalTimeStep << m_outputFormat;

  if(m_movingGrid){
    GridFileName << "Grid" << globalTimeStep << ".hdf5";
    tempG = GridFileName.str();
  }else{
    tempG = "../"+m_gridInputFileName;
  }
  strcpy(gridFile,tempG.c_str());

  ZFSInt file =io_openfile("hdf5", (fileName.str()).c_str(), "collective", m_zfsStrctrdComm);

  writeHeaderAttributes(file, gridFile, "field");
  writePropertiesAsAttributes(file, "");

  io_create_iattribute(file,"", "noFields", 1);
  io_write_iattribute1(file,"", "noFields", 1, &m_noLineOutput);

  for(ZFSId lineId=0; lineId<m_noLineOutput; lineId++) {
    ZFSInt dataOffset[2] = {0,0};
    ZFSInt dataSize[2] = {m_lineNoPoints2d[lineId],m_lineNoPoints[lineId]};
    ZFSInt fieldOffset = m_fieldOffset[lineId];

    ZFSId noDims = 1;
    if(m_fieldInterpolation) {
      noDims = 2;
    }

    stringstream path;
    path << lineId;
    ZFSString solutionpath = "field";
    solutionpath += path.str();
    const char* dsetname = solutionpath.c_str();

    for(ZFSId v=0; v<PV->noVariables; v++) {
      io_create_ddataset(file,dsetname, (m_pvariableNames[v]).c_str(), noDims, dataSize);
    }

    io_create_ddataset(file,dsetname, "x", noDims, dataSize);
    io_create_ddataset(file,dsetname, "y", noDims, dataSize);
    io_create_ddataset(file,dsetname, "z", noDims, dataSize);

    if(domainId() == 0) {
      for(ZFSId v=0; v<PV->noVariables; v++) {
        io_write_ddataset_part(file, dsetname, m_pvariableNames[v].c_str(), noDims, dataSize, dataOffset, &m_interpolatedVarsGlobal[v][fieldOffset]);
      }

      io_write_ddataset_part(file, dsetname, "x",noDims, dataSize, dataOffset, &m_pointCoordinates[0][fieldOffset]);
      io_write_ddataset_part(file, dsetname, "y",noDims, dataSize, dataOffset, &m_pointCoordinates[1][fieldOffset]);
      io_write_ddataset_part(file, dsetname, "z",noDims, dataSize, dataOffset, &m_pointCoordinates[2][fieldOffset]);
    } else {
      dataSize[0] = 0; dataSize[1] = 0;
      for(ZFSId v=0; v<PV->noVariables; v++) {
        io_write_ddataset_part(file, dsetname, m_pvariableNames[v].c_str(),1, dataSize, dataOffset, NULL);
      }

      io_write_ddataset_part(file, dsetname, "x",noDims, dataSize, dataOffset, NULL);
      io_write_ddataset_part(file, dsetname, "y",noDims, dataSize, dataOffset, NULL);
      io_write_ddataset_part(file, dsetname, "z",noDims, dataSize, dataOffset, NULL);
    }
  }

  io_closefile(file);
}

void ZFSStrctrdBlck3D::applyInviscidBoundaryCondition()
{
  TRACE();
}

void ZFSStrctrdBlck3D::applyViscousBoundaryCondition()
{
  TRACE();
}

void ZFSStrctrdBlck3D::getSampleVariables(ZFSId cellId, ZFSFloat* cellVars) {
  cellVars[PV->U]  = m_cells->pvariables[PV->U][cellId];
  cellVars[PV->V]  = m_cells->pvariables[PV->V][cellId];
  cellVars[PV->W]  = m_cells->pvariables[PV->W][cellId];
  cellVars[PV->RHO]= m_cells->pvariables[PV->RHO][cellId];
  cellVars[PV->P]  = m_cells->pvariables[PV->P][cellId];
}

ZFSFloat ZFSStrctrdBlck3D::getSampleVorticity(ZFSId cellId, ZFSId dim) {
  return m_cells->fq[FQ->VORTICITY[dim]][cellId];
}


/**
 *
 * @author Frederik Temme, Jan 14, 2015
 * modified 14.1.2016
 *
 */
void ZFSStrctrdBlck3D::loadSampleFile(ZFSString fileName){ // loading files for averaging (pre- and postsolve)
  TRACE();
  ZFSInt FileId =-1;

  //open the file
  FileId = io_openfile("hdf5" ,fileName.c_str(),"collective", m_zfsStrctrdComm);

  //check whether the file exists
  if(FileId==-1){
    zfsTerm(1, __CALLING_FUNCTION__, "Sample file not found");
  }

  //now read in the data!
  ZFSInt start[ 3 ];
  ZFSInt Id = m_partition->outputBoxInfo[domainId()]->cpu;
  //the offset in the file
  for (ZFSInt i=0; i<nDim;i++){
    //get offset in file
    start[i]=m_partition->outputBoxInfo[Id]->offset[i];
  }

  stringstream blockNumber;
  blockNumber << m_inputBlockId;
  ZFSString blockPathStr = "/block";
  blockPathStr += blockNumber.str();
  const char* blockPath = blockPathStr.c_str();

  for(ZFSId var=0; var<PV->noVariables; var++) {
    io_read_ddataset_part1d1(FileId,blockPath, m_pvariableNames[var].c_str(), nDim, start, m_nActiveCells, m_cells->pvariables[var]);
  }

  io_closefile(FileId);
  shiftCellValuesRestart<true>();
}


ZFSFloat ZFSStrctrdBlck3D::dvardxyz(ZFSId IJK,ZFSId dir,ZFSFloat* var) {
  const ZFSId INC[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};

  const ZFSId IPJK  = IJK+INC[0];
  const ZFSId IMJK  = IJK-INC[0];
  const ZFSId IJPK  = IJK+INC[1];
  const ZFSId IJMK  = IJK-INC[1];
  const ZFSId IJKP  = IJK+INC[2];
  const ZFSId IJKM  = IJK-INC[2];

  const ZFSFloat dvardxi=var[IPJK]-var[IMJK];
  const ZFSFloat dvardet=var[IJPK]-var[IJMK];
  const ZFSFloat dvardze=var[IJKP]-var[IJKM];

  const ZFSFloat ddxyz=  dvardxi * m_cells->cellMetrics[IJK][ xsd * 3 + dir ] +
    dvardet * m_cells->cellMetrics[IJK][ ysd * 3 + dir ] +
    dvardze * m_cells->cellMetrics[IJK][ zsd * 3 + dir ];

  return ddxyz/m_cells->cellJac[IJK];
}

ZFSFloat ZFSStrctrdBlck3D::dvardx(ZFSId IJK, ZFSFloat* var) {
  const ZFSId INC[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};

  const ZFSId IPJK  = IJK+INC[0];
  const ZFSId IMJK  = IJK-INC[0];
  const ZFSId IJPK  = IJK+INC[1];
  const ZFSId IJMK  = IJK-INC[1];
  const ZFSId IJKP  = IJK+INC[2];
  const ZFSId IJKM  = IJK-INC[2];

  const ZFSFloat dvardxi=var[IPJK]-var[IMJK];
  const ZFSFloat dvardet=var[IJPK]-var[IJMK];
  const ZFSFloat dvardze=var[IJKP]-var[IJKM];

  const ZFSFloat ddx =  dvardxi * m_cells->cellMetrics[IJK][ xsd * 3 + xsd ] +
    dvardet * m_cells->cellMetrics[IJK][ ysd * 3 + xsd ] +
    dvardze * m_cells->cellMetrics[IJK][ zsd * 3 + xsd ];

  return ddx/m_cells->cellJac[IJK];
}

/**
 * Loads the postprocessing restart file to continue the postprocessing
 *
 * @author Frederik Temme, Dez 17, 2015
 * modified 06.01.2016
 */

void ZFSStrctrdBlck3D::loadAverageRestartFile(const char* fileName, ZFSFloat** sum, ZFSFloat** square, ZFSFloat** cube, ZFSFloat** fourth){
  TRACE();
  zfs_log << "loading average restart file ... " << endl;
  ZFSInt restartFileId =-1;
  if(domainId()==0) {cout << "Opening file: " << fileName << endl;}

  //open the file
  restartFileId = io_openfile("hdf5", fileName,"collective", m_zfsStrctrdComm);

  //check whether the file exists
  if(restartFileId==-1){
    zfsTerm(1, __CALLING_FUNCTION__, "AverageRestartFile cannot be found");
  }

  //now read in the data!
  zfs_log << "-> reading in the data ... " << endl;
  io_read_iattribute1(restartFileId,"", "noSamples", &m_noSamples);
  zfs_log << "Current number of postprocessing samples: " << m_noSamples << endl;

  ZFSInt start[ 3 ];
  ZFSInt Id = m_partition->outputBoxInfo[domainId()]->cpu;
  //the offset in the file
  for (ZFSInt i=0; i<nDim;i++) {
    //get offset in file
    start[i]=m_partition->outputBoxInfo[Id]->offset[i];
  }

  stringstream blockNumber;
  blockNumber << m_inputBlockId;
  ZFSString blockPathStr = "/block";
  blockPathStr += blockNumber.str();
  const char* blockPath = blockPathStr.c_str();
  ZFSId offset = 0;

  io_read_ddataset_part1d1(restartFileId, blockPath, "u", nDim, start, m_nActiveCells, sum[0]);
  io_read_ddataset_part1d1(restartFileId, blockPath, "v", nDim, start, m_nActiveCells, sum[1]);
  io_read_ddataset_part1d1(restartFileId, blockPath, "w", nDim, start, m_nActiveCells, sum[2]);
  io_read_ddataset_part1d1(restartFileId, blockPath, "rho", nDim, start, m_nActiveCells, sum[3]);
  io_read_ddataset_part1d1(restartFileId, blockPath, "p", nDim, start, m_nActiveCells, sum[4]);
  offset = noVariables();

  if(m_averageVorticity) {
    io_read_ddataset_part1d1(restartFileId,blockPath, "vortx", nDim, start, m_nActiveCells, sum[offset+0]);
    io_read_ddataset_part1d1(restartFileId,blockPath, "vorty", nDim, start, m_nActiveCells, sum[offset+1]);
    io_read_ddataset_part1d1(restartFileId,blockPath, "vortz", nDim, start, m_nActiveCells, sum[offset+2]);
  }

  io_read_ddataset_part1d1(restartFileId,blockPath, "uu", nDim, start, m_nActiveCells, square[0]);
  io_read_ddataset_part1d1(restartFileId,blockPath, "vv", nDim, start, m_nActiveCells, square[1]);
  io_read_ddataset_part1d1(restartFileId,blockPath, "ww", nDim, start, m_nActiveCells, square[2]);
  io_read_ddataset_part1d1(restartFileId,blockPath, "uv", nDim, start, m_nActiveCells, square[3]);
  io_read_ddataset_part1d1(restartFileId,blockPath, "vw", nDim, start, m_nActiveCells, square[4]);
  io_read_ddataset_part1d1(restartFileId,blockPath, "uw", nDim, start, m_nActiveCells, square[5]);

  io_read_ddataset_part1d1(restartFileId,blockPath, "pp", nDim, start, m_nActiveCells, square[6]);

  if(m_averageVorticity) {
    io_read_ddataset_part1d1(restartFileId,blockPath, "vortxvortx", nDim, start, m_nActiveCells, square[7]);
    io_read_ddataset_part1d1(restartFileId,blockPath, "vortyvorty", nDim, start, m_nActiveCells, square[8]);
    io_read_ddataset_part1d1(restartFileId,blockPath, "vortzvortz", nDim, start, m_nActiveCells, square[9]);
  }

  if(m_kurtosis || m_skewness){
    io_read_ddataset_part1d1(restartFileId,blockPath, "uuu", nDim, start, m_nActiveCells, cube[0]);
    io_read_ddataset_part1d1(restartFileId,blockPath, "vvv", nDim, start, m_nActiveCells, cube[1]);
    io_read_ddataset_part1d1(restartFileId,blockPath, "www", nDim, start, m_nActiveCells, cube[2]);
  }

  if(m_kurtosis){
    io_read_ddataset_part1d1(restartFileId,blockPath, "uuuu", nDim, start, m_nActiveCells, fourth[0]);
    io_read_ddataset_part1d1(restartFileId,blockPath, "vvvv", nDim, start, m_nActiveCells, fourth[1]);
    io_read_ddataset_part1d1(restartFileId,blockPath, "wwww", nDim, start, m_nActiveCells, fourth[2]);
  }

  io_closefile(restartFileId);
  zfs_log << "loading Restart file ... SUCCESSFUL " << endl;

  shiftAverageCellValuesRestart();
}

/*
 * Loads the averaged variables again to do further postprocessing
 *
 * @author Marian Albers, Apr 20, 2016
 */
void ZFSStrctrdBlck3D::loadAveragedVariables(const char* fileName){
  TRACE();
  zfs_log << "loading averaged variables file ... " << endl;
  ZFSInt restartFileId =-1;

  //open the file
  restartFileId = io_openfile("hdf5" ,fileName,"collective", m_zfsStrctrdComm);

  //check whether the file exists
  if(restartFileId==-1){
    zfsTerm(1, __CALLING_FUNCTION__, "AverageRestartFile cannot be found");
  }

  ZFSInt start[ 3 ];
  ZFSInt Id = m_partition->outputBoxInfo[domainId()]->cpu;
  //the offset in the file
  for (ZFSInt i=0; i<nDim;i++) {
    //get offset in file
    start[i]=m_partition->outputBoxInfo[Id]->offset[i];
  }

  stringstream blockNumber;
  blockNumber << m_inputBlockId;
  ZFSString blockPathStr = "/block";
  blockPathStr += blockNumber.str();
  const char* blockPath = blockPathStr.c_str();

  for(ZFSId var=0; var<getNoPPVars(); var++) {
    io_read_ddataset_part1d1(restartFileId, blockPath, (m_avgVariableNames[var]).c_str(), nDim, start, m_nActiveCells, m_summedVars[var]);
  }

  io_closefile(restartFileId);
  zfs_log << "loading Restart file ... SUCCESSFUL " << endl;

  shiftAverageCellValues();
}

/**
 *
 * @author Frederik Temme, Dez 17, 2015
 * modified 12.01.2015
 */
void ZFSStrctrdBlck3D::shiftAverageCellValuesRestart(){
  TRACE();
  ZFSId cellId_org=0;
  ZFSId cellId=0;
  ZFSInt i_new, j_new, k_new;

  //accounting for the ghost layers and shift the values to the right place
  for(ZFSId k=(m_nActiveCells[0]-1); k>=0; k-- ){
    for(ZFSId j=(m_nActiveCells[1]-1); j>=0; j-- ){
      for(ZFSId i=(m_nActiveCells[2]-1); i>=0; i--){
        cellId_org=i+(j+k*m_nActiveCells[1])*m_nActiveCells[2];
        i_new=i+m_noGhostLayers;
        j_new=j+m_noGhostLayers;
        k_new=k+m_noGhostLayers;
        cellId=i_new+(j_new+k_new*m_nCells[1])*m_nCells[2];

        for(ZFSId var = 0; var < getNoPPVars(); var++) {
          m_summedVars[var][cellId]=m_summedVars[var][cellId_org];
          m_summedVars[var][cellId_org]=F0;
        }

        for(ZFSId var = 0; var < getNoPPSquareVars(); var++) {
          m_square[var][cellId]=m_square[var][cellId_org];
          m_square[var][cellId_org]=F0;
        }

        if(m_kurtosis || m_skewness){
          for(ZFSId var = 0; var < nDim; var++) {
            m_cube[var][cellId]=m_cube[var][cellId_org];
            m_cube[var][cellId_org]=F0;
          }
        }

        if(m_kurtosis){
          for(ZFSId var = 0; var < nDim; var++) {
            m_fourth[var][cellId]=m_fourth[var][cellId_org];
            m_fourth[var][cellId_org]=F0;
          }
        }
      }
    }
  }
}

/**
 * Shifts the averaged variables
 * @author Marian Albers, Apr 20, 2016
 */
void ZFSStrctrdBlck3D::shiftAverageCellValues(){
  TRACE();
  ZFSId cellId_org=0;
  ZFSId cellId=0;
  ZFSInt i_new, j_new, k_new;

  //accounting for the ghost layers and shift the values to the right place
  for(ZFSId k=(m_nActiveCells[0]-1); k>=0; k-- ){
    for(ZFSId j=(m_nActiveCells[1]-1); j>=0; j-- ){
      for(ZFSId i=(m_nActiveCells[2]-1); i>=0; i--){
        cellId_org=i+(j+k*m_nActiveCells[1])*m_nActiveCells[2];
        i_new=i+m_noGhostLayers;
        j_new=j+m_noGhostLayers;
        k_new=k+m_noGhostLayers;
        cellId=i_new+(j_new+k_new*m_nCells[1])*m_nCells[2];

        for(ZFSId var = 0; var < getNoPPVars(); var++) {
          m_summedVars[var][cellId]=m_summedVars[var][cellId_org];
          m_summedVars[var][cellId_org]=F0;
        }
      }
    }
  }
}

void ZFSStrctrdBlck3D::computePrimitiveVariables() {
  const ZFSFloat gammaMinusOne = m_gamma - 1.0;

  ZFSFloat** const RESTRICT cvars = m_cells->variables;
  ZFSFloat** const RESTRICT pvars = m_cells->pvariables;

  for(ZFSId k=m_noGhostLayers; k < m_nCells[0]-m_noGhostLayers; ++k) {
    for(ZFSId j=m_noGhostLayers; j < m_nCells[1]-m_noGhostLayers; ++j) {
      for(ZFSId i=m_noGhostLayers; i < m_nCells[2]-m_noGhostLayers; ++i) {
        const ZFSId cellId = cellIndex(i,j,k);
        const ZFSFloat fRho = F1 / cvars[CV->RHO][cellId];
        ZFSFloat velPOW2 = F0;
        for(ZFSId vel = 0;  vel < nDim; ++vel) { // compute velocity
          pvars[vel][cellId] = cvars[vel][cellId] * fRho;
          velPOW2 += POW2(pvars[vel][cellId]);
        }

        // density and pressure:
        pvars[PV->RHO][cellId] = cvars[CV->RHO][cellId]; // density
        pvars[PV->P][cellId] = gammaMinusOne * (cvars[CV->RHO_E][cellId]
                                                - F1B2 * pvars[PV->RHO][cellId] * velPOW2);

        for(ZFSId ransVar=0; ransVar < m_noRansEquations; ransVar++) {
          cvars[CV->RANS_VAR[ransVar]][cellId] = zfsMAX(cvars[CV->RANS_VAR[ransVar]][cellId], F0);
          pvars[PV->RANS_VAR[ransVar]][cellId] = cvars[CV->RANS_VAR[ransVar]][cellId]*fRho;
        }
      }
    }
  }
}


void ZFSStrctrdBlck3D::allocateSingularities() {
  for(ZFSInt i=0; i<m_hasSingularity; ++i) {
    ZFSInt len[3];
    m_singularity[i].totalPoints=1;
    m_singularity[i].totalCells=1;

    for(ZFSInt j=0; j<nDim; j++) {
      len[j]=m_singularity[i].end[j]-m_singularity[i].start[j];
      m_singularity[i].totalPoints*=(len[j]+1);
      m_singularity[i].totalCells *=len[j];
    }

    m_singularity[i].ReconstructionConstants=  new ZFSFloat* [4];
    //4 unknowns and 2*Nstar cells
    for(ZFSInt j=0; j< 4; j++) {
      m_singularity[i].ReconstructionConstants[j]=new ZFSFloat [m_singularity[i].totalCells*m_singularity[i].Nstar*2];
    }
  }
}

//see also the function in zfslbmblckdxqy.cpp for this function
//implemented originally by Georg Eitel Amor for LBM

void ZFSStrctrdBlck3D::initFFTW(fftw_complex* uPhysField, fftw_complex* vPhysField, fftw_complex* wPhysField, ZFSInt lx, ZFSInt ly, ZFSInt lz, ZFSInt noPeakModes)
{
  TRACE();
  //first check for odd numbers of cells in each direction
  if (lx%2 != 0 || ly%2 != 0 || lz%2 != 0){
    stringstream errorMessage;
    errorMessage << " FFTInit: Domain size must NOT be an odd number!: (lx)x(ly)x(lz)-> "<< lx << "x"<< ly <<"x"<<lz<< endl;
    zfsTerm(1, __CALLING_FUNCTION__, errorMessage.str());
  }

  zfs_log <<" --- initializing FFTW --- " <<endl;
  zfs_log <<" domain size = "<<lx<<"x"<<ly<<"x"<<lz<<endl;

  ZFSFloat waveVector[3], k0;

  complex<double>* fourierCoefficient = new complex<double>[3];

  fftw_complex *uHatField, *vHatField, *wHatField;

  fftw_plan planU, planV, planW;

  //1) Allocation of Fourier coefficients
  uHatField = (fftw_complex*) fftw_malloc(lx*ly*lz * sizeof(fftw_complex));
  vHatField = (fftw_complex*) fftw_malloc(lx*ly*lz * sizeof(fftw_complex));
  wHatField = (fftw_complex*) fftw_malloc(lx*ly*lz * sizeof(fftw_complex));

  //2) Creation of the plans for the FFTW
  planU = fftw_plan_dft_3d(lx, ly, lz, uHatField, uPhysField, FFTW_BACKWARD, FFTW_MEASURE);
  planV = fftw_plan_dft_3d(lx, ly, lz, vHatField, vPhysField, FFTW_BACKWARD, FFTW_MEASURE);
  planW = fftw_plan_dft_3d(lx, ly, lz, wHatField, wPhysField, FFTW_BACKWARD, FFTW_MEASURE);

  for (ZFSId p=0; p<lx; p++)
  {
    for (ZFSId q=0; q<ly; q++)
    {
      for (ZFSId r=0; r<lz; r++)
      {
        ZFSId cellId=r+lz*(q+ly*p);
        //u-component
        uHatField[cellId][0]  = 0.0;
        uHatField[cellId][1]  = 0.0;
        uPhysField[cellId][0]  = 0.0;
        uPhysField[cellId][1] = 0.0;

        //v-component
        vHatField[cellId][0]  = 0.0;
        vHatField[cellId][1]  = 0.0;
        vPhysField[cellId][0] = 0.0;
        vPhysField[cellId][1] = 0.0;

        //w-component
        wHatField[cellId][0]  = 0.0;
        wHatField[cellId][1]  = 0.0;
        wPhysField[cellId][0] = 0.0;
        wPhysField[cellId][1] = 0.0;

      }
    }
  }

  //IMPORTANT NOTICE ON USE OF FFTW:
  //comment from Georg Eitel Amor:
  // FFTW stores the coefficients for positive wavenumbers in the first half of the array,
  // and those for negative wavenumbers in reverse order in the second half.
  // [0, 1, ... , N/2-1, N/2, ... , N-1]
  //  - the entry for zero-wavenumber is at position 0
  //  - the k-th entry and the (N-k)th entry correspond to wavenumbers with opposite sign
  //  - the entry at position N/2 corresponds to the Nyquist wavenumber and appears only once

  // peak wave number of energy spectrum
  k0 = 2.0 * PI / (lx/noPeakModes);

  for (ZFSId p=0; p<=lx/2; p++)
  {
    for (ZFSId q=0; q<=ly/2; q++)
    {
      for (ZFSId r=0; r<=lz/2; r++)
      {

        // wave-vector: k(p,q,r) = (2 \pi p / lx, 2 \pi q / ly, 2 \pi r / lz)
        waveVector[0] = (p) * 2.0 * PI / lx;
        waveVector[1] = (q) * 2.0 * PI / ly;
        waveVector[2] = (r) * 2.0 * PI / lz;

        getFourierCoefficients(waveVector, k0, fourierCoefficient);

        // 1. Positive frequencies:
        uHatField[r+lz*(q+ly*p)][0] = real(fourierCoefficient[0]);
        uHatField[r+lz*(q+ly*p)][1] = imag(fourierCoefficient[0]);

        vHatField[r+lz*(q+ly*p)][0] = real(fourierCoefficient[1]);
        vHatField[r+lz*(q+ly*p)][1] = imag(fourierCoefficient[1]);

        wHatField[r+lz*(q+ly*p)][0] = real(fourierCoefficient[2]);
        wHatField[r+lz*(q+ly*p)][1] = imag(fourierCoefficient[2]);

        // 2. Negative frequencies:
        if(p>1 && q>1 && r>1)
        {
          if(p<lx/2 && q<ly/2 && r<lz/2)
          {
            // since the physical velocity field is real, the coefficients for negative frequencies
            // are the complex conjugate of those for positive frequencies
            uHatField[(lz-r)+lz*((ly-q)+ly*(lx-p))][0] = uHatField[r+lz*(q+ly*p)][0];
            uHatField[(lz-r)+lz*((ly-q)+ly*(lx-p))][1] = -uHatField[r+lz*(q+ly*p)][1];

            vHatField[(lz-r)+lz*((ly-q)+ly*(lx-p))][0] = vHatField[r+lz*(q+ly*p)][0];
            vHatField[(lz-r)+lz*((ly-q)+ly*(lx-p))][1] = -vHatField[r+lz*(q+ly*p)][1];

            wHatField[(lz-r)+lz*((ly-q)+ly*(lx-p))][0] = wHatField[r+lz*(q+ly*p)][0];
            wHatField[(lz-r)+lz*((ly-q)+ly*(lx-p))][1] = -wHatField[r+lz*(q+ly*p)][1];
          }
        }

      }
    }
  }

  // Do Fourier transform (backward, see plan definition)
  // Definition in one dimension:
  // u(x) = \sum_{j=0}^{lx-1} \hat{u}_j exp(i 2 \pi j x / lx)

  fftw_execute(planU);
  fftw_execute(planV);
  fftw_execute(planW);


  // normalize (this preserves the norm of the basis functions)
  for (ZFSId p=0; p<lx; p++)
  {
    for (ZFSId q=0; q<ly; q++)
    {
      for (ZFSId r=0; r<lz; r++)
      {

        uPhysField[r+lz*(q+ly*p)][0] /= sqrt(double(lx*ly*lz));
        vPhysField[r+lz*(q+ly*p)][0] /= sqrt(double(lx*ly*lz));
        wPhysField[r+lz*(q+ly*p)][0] /= sqrt(double(lx*ly*lz));

        uPhysField[r+lz*(q+ly*p)][1] /= sqrt(double(lx*ly*lz));
        vPhysField[r+lz*(q+ly*p)][1] /= sqrt(double(lx*ly*lz));
        wPhysField[r+lz*(q+ly*p)][1] /= sqrt(double(lx*ly*lz));

      }
    }
  }

  fftw_destroy_plan(planU);
  fftw_destroy_plan(planV);
  fftw_destroy_plan(planW);
  fftw_free(uHatField);
  fftw_free(vHatField);
  fftw_free(wHatField);
}

/** brief Generates a single complex coefficient of Fourier series
 *  Original Implementation by Georg Eitel Amor
 *  for a given wavenumber k, and a certain energy spectrum
 *  (see Appendix of Orszag, 1969)
 *
 */
void ZFSStrctrdBlck3D::getFourierCoefficients(ZFSFloat* k, ZFSFloat k0, complex<ZFSFloat>* fourierCoefficient) {
  TRACE();
  ZFSFloat r[6], s[6], kAbs, energy;
  complex<ZFSFloat> uHat, vHat, wHat;
  //complex<double>* fourierCoefficient;

    kAbs = sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
    
    // the zero-frequency component is always set to zero, so there is no offset
    if(approx(kAbs, F0, m_eps))
      {
	
	//fourierCoefficient = new complex<double>[3];
	
	fourierCoefficient[0] = complex<double>(0,0);
	fourierCoefficient[1] = complex<double>(0,0);
	fourierCoefficient[2] = complex<double>(0,0);
      }
    else
      {
	
	//energy = (kAbs/k0)*(kAbs/k0)*(kAbs/k0)*(kAbs/k0) * exp(-2.0*(kAbs/k0)*(kAbs/k0));
        //energy = pow(kAbs/k0,8.0) * exp(-4.0*(kAbs/k0)*(kAbs/k0)); // set spectral distribution
        energy = pow(kAbs/k0,4.0) * exp(-2.0*(kAbs/k0)*(kAbs/k0)); // set spectral distribution
        energy *= exp(2.0) * 0.499 * ( m_Ma/F1BCS) * ( m_Ma/F1BCS ); // set maximal fluctuation amplitude to 20% of the freestream velocity (for 128^3: 0.88)
	
	// determine Fourier coefficients:
	// r and s are Independant random vector fields with independant 
	// components (zero mean and rms according to energy spectrum).
	// Each vector has three components for k and another three for -k.

    for (ZFSId i=0; i<6; i++){
      r[i] = randnormal( 0.0, PI*sqrt(energy)/(SQRT2*kAbs) );
      //r[i] = randNumGen.randNorm(0.0, PI*sqrt(energy)/(SQRT2*kAbs));
      s[i] = randnormal( 0.0, PI*sqrt(energy)/(SQRT2*kAbs) );
      //s[i] = randNumGen.randNorm(0.0, PI*sqrt(energy)/(SQRT2*kAbs));
    }

    uHat = (1.0 - k[0]*k[0]/(kAbs*kAbs)) * complex<double>(r[0] + r[3], s[0] - s[3])
      - k[0]*k[1]/(kAbs*kAbs) * complex<double>(r[1] + r[4], s[1] - s[4])
      - k[0]*k[2]/(kAbs*kAbs) * complex<double>(r[2] + r[5], s[2] - s[5]);


    vHat = -k[1]*k[0]/(kAbs*kAbs) * complex<double>(r[0] + r[3], s[0] - s[3])
      + (1.0 - k[1]*k[1]/(kAbs*kAbs)) * complex<double>(r[1] + r[4], s[1] - s[4])
      - k[1]*k[2]/(kAbs*kAbs) * complex<double>(r[2] + r[5], s[2] - s[5]);


    wHat = -k[2]*k[0]/(kAbs*kAbs) * complex<double>(r[0] + r[3], s[0] - s[3])
      - k[2]*k[1]/(kAbs*kAbs) * complex<double>(r[1] + r[4], s[1] - s[4])
      + (1.0 - k[2]*k[2]/(kAbs*kAbs)) * complex<double>(r[2] + r[5], s[2] - s[5]);

    //fourierCoefficient = new complex<double>[3];

    // fourierCoefficient[0] = complex<double>(sqrt(2*energy)/SQRT2,sqrt(2*energy)/SQRT2);// uHat;
    // fourierCoefficient[1] = complex<double>(sqrt(2*energy)/SQRT2,sqrt(2*energy)/SQRT2);//vHat;
    // fourierCoefficient[2] = complex<double>(sqrt(2*energy)/SQRT2,sqrt(2*energy)/SQRT2);//wHat;

    fourierCoefficient[0] = uHat;
    fourierCoefficient[1] = vHat;
    fourierCoefficient[2] = wHat;
    //return fourierCoefficient;
    }
}

void ZFSStrctrdBlck3D::computeReconstructionConstantsSVD()
{
  ZFSInt nghbr[30],dim;
  ZFSInt start[3],end[3];
  m_orderOfReconstruction = 1;
  const ZFSId recDim = ( m_orderOfReconstruction == 2 ) ? (IPOW2[nDim]+1) : nDim+1;
  ZFSInt maxNoSingularityRecNghbrIds = 14;
  ZFSFloatScratchSpace tmpA(maxNoSingularityRecNghbrIds,recDim,__CALLING_FUNCTION__,"tmpA");
  ZFSFloatScratchSpace tmpC(recDim,maxNoSingularityRecNghbrIds,__CALLING_FUNCTION__,"tmpC");
  ZFSFloatScratchSpace weights(maxNoSingularityRecNghbrIds,__CALLING_FUNCTION__,"weights");
  ZFSFloat counter = F0;
  ZFSFloat avg = F0;
  ZFSFloat maxc = F0;

  for(ZFSInt i=0; i<m_hasSingularity; ++i) {
    if(m_singularity[i].BC==6000) {
      ZFSInt totalCells=1;
      ZFSInt len1[3];

      //(p)reset the reconstruction constants
      for( ZFSId n = 0; n < 4 ; ++n ) {
        for(ZFSId m = 0; m < m_singularity[i].totalCells*m_singularity[i].Nstar*2 ; ++m ) {
          m_singularity[i].ReconstructionConstants[ n ][ m ] = -999;
        }
      }

      for(ZFSInt j=0; j<nDim; j++) {
        len1[j]=m_singularity[i].end[j]-m_singularity[i].start[j];
        if(len1[j]!=0)  totalCells*=len1[j];
      }

      for( ZFSInt n = 0; n < 3 ; ++n ) {
        if(m_singularity[i].end[n]-m_singularity[i].start[n]>1) {
          dim=n;
          start[n]=m_singularity[i].start[n]+1;
          end[n]=m_singularity[i].end[n]-1;
        } else {
          start[n]=m_singularity[i].start[n];
          end[n]=m_singularity[i].end[n];
        }
      }

      for( ZFSInt kk = start[2]; kk <end[2]; ++kk ) {
        for( ZFSInt jj = start[1]; jj <end[1]; ++jj ) {
          for( ZFSInt ii = start[0]; ii <end[0]; ++ii ) {
            ZFSInt count=0;
            ZFSInt temp[3]={0,0,0};
            temp[dim]=1;

            nghbr[count++]=cellIndex(ii,jj,kk);
            nghbr[count++]=cellIndex(ii+temp[0],jj+temp[1],kk+temp[2]);

            //the coordinates of the corner where the viscousflux should be corrected.
            ZFSInt ijk    = getPointIdFromCell( ii+ m_singularity[i].Viscous[0], jj+ m_singularity[i].Viscous[1], kk+ m_singularity[i].Viscous[2] );
            ijk= getPointIdfromPoint(ijk,1,1,1);

            for(ZFSInt m=0; m< m_singularity[i].Nstar-1; ++m) {
              ZFSInt *change=  m_singularity[i].displacement[m];
              nghbr[count++]=cellIndex(ii+change[0],jj+change[1],kk+change[2]);
              nghbr[count++]=cellIndex(ii+temp[0]+change[0],jj+temp[1]+change[1],kk+temp[2]+change[2]);
            }

            if(count!=m_singularity[i].Nstar*2) {
              cerr << "Something wrong with the singularities in the LS coeffiecient computation" << endl;
            }

            //weighted Least square
            weights.fill(F0);

            //Compute weights with RBF (take mean distance as R0)
            for( ZFSId n = 0; n < count; n++ ) {
              ZFSId nghbrId = nghbr[n];
              ZFSFloat dxdx = F0;
              for ( ZFSId m = 0; m < nDim; ++m ) {
                dxdx += POW2( m_cells->coordinates[m][nghbrId] -  m_coordinates[m][ijk] );
              }

              weights[n] =1/dxdx;   // RBF( dxdx, POW2( dist) );
            }

            ZFSInt id2=ii-start[0]+((jj-start[1])+(kk-start[2])*len1[1])*len1[0];
            ZFSInt ID=id2*m_singularity[i].Nstar*2;

            ZFSFloat condNum = computeRecConstSVD(ijk, count, nghbr, ID, i, tmpA, tmpC, weights, recDim);
            avg += condNum;
            maxc = zfsMAX(maxc,condNum);
            counter += F1;
            if ( condNum < F0 || condNum > 1e7 || std::isnan(condNum) ) {
              cerr << domainId() << " SVD decomposition for pointId " << ijk
                   << " with large condition number: "
                   << condNum << " num of neighbor" << count << "x" << recDim << " "
                   << " coords " << m_coordinates[0][ijk] << ", " <<  m_coordinates[1][ijk]
                   << ", " <<  m_coordinates[2][ijk] << endl;
            }
          }
        }
      }
    }
  }
}

//this one is for the singularities
ZFSFloat ZFSStrctrdBlck3D::computeRecConstSVD(const ZFSId ijk,
                                              const ZFSId noNghbrIds,
                                              ZFSInt* nghbr,  ZFSInt ID,
                                              ZFSInt sID,
                                              ZFSFloatScratchSpace& tmpA,
                                              ZFSFloatScratchSpace& tmpC,
                                              ZFSFloatScratchSpace& weights,
                                              const ZFSId recDim)
{
  if ( noNghbrIds == 0 ) return F0;

  const ZFSFloat normalizationFactor = 1/0.01; // reduces the condition number of the eq system

  for( ZFSId n = 0; n < noNghbrIds; n++ ) {
    ZFSId nghbrId = nghbr[n];
    ZFSFloat dx[3];
    for ( ZFSId i = 0; i < nDim; i++ ) {
      dx[i] = ( m_cells->coordinates[i][nghbrId] - m_coordinates[i][ijk] ) * normalizationFactor;
    }

    tmpA(n,0)=F1 * normalizationFactor;
    for ( ZFSId i = 0; i < nDim; i++ ) {
      tmpA(n,i+1) = dx[i];
    }
  }
  ZFSFloat condNum = svdPreSolve( tmpA, weights, tmpC, noNghbrIds, recDim );

  if ( condNum < F1 ) {
    cerr << domainId() << ": SVD failed for point " << ijk  << " " << endl;
  }

  for( ZFSId n = 0;  n < noNghbrIds;  n++ ) {
    for ( ZFSId i = 0; i < nDim+1; i++ ) {
      m_singularity[sID].ReconstructionConstants[i][ID+n]  = tmpC(i,n) * normalizationFactor;

    }
  }

  return condNum;
}


/**
 * \brief Compute a weighted pseudo inverse of the m*n matrix 'A' by singular value decomposition, also for underdetermined systems
 * \author Lennart Schneiders
 */
ZFSFloat ZFSStrctrdBlck3D::svdPreSolve( ZFSFloatScratchSpace& A,
                                        ZFSFloatScratchSpace& weights,
                                        ZFSFloatScratchSpace& AInv,
                                        const ZFSId m, const ZFSId n )
{
  const ZFSId p = min( m, n );

  ZFSFloatScratchSpace U( m, p, __CALLING_FUNCTION__, "U" );
  ZFSFloatScratchSpace V( n, p, __CALLING_FUNCTION__, "V" );
  ZFSFloatScratchSpace S( p, __CALLING_FUNCTION__, "S" );
  ZFSId mrank = -1;
  if( m >= n ) {
    for ( ZFSId i = 0; i < m; i++ ) {
      for ( ZFSId j = 0; j < n; j++ ) {
        U(i,j) = weights[i] * A(i,j);
      }
    }
#define SVD_QR_PRECONDITIONING
#ifdef SVD_QR_PRECONDITIONING
    ZFSBool precon = (m > (n+n/2+1)); //only improves performance for significantly overdetermined matrices
    if ( precon ) {
      if ( QRdec(U, S, m, n ) ) { //A=U=Q*R, returns true if not rank deficient
        ZFSFloatScratchSpace Q( m, p, __CALLING_FUNCTION__, "Q" );
        ZFSFloatScratchSpace R( p, p, __CALLING_FUNCTION__, "R" );
        getQ(U, Q, m, n); //Q
        getR(U, S, R, m, n); //R
        mrank = gsl_linalg_SV_decomp_jacobi (R, V, S); //R=R*S*V^t
        for ( ZFSId i = 0; i < m; i++ ) {    //U^t=R^t*Q^t ==> U=QR
          for ( ZFSId j = 0; j < p; j++ ) {
            U(i,j) = F0;
            for ( ZFSId k = 0; k < p; k++ ) {
              U(i,j) += Q(i,k)*R(k,j);
            }
          }
        }
      }
      else { //if A rank deficient don't use QR preconditioning
        for ( ZFSId i = 0; i < m; i++ ) {
          for ( ZFSId j = 0; j < n; j++ ) {
            U(i,j) = weights[i] * A(i,j);
          }
        }
        precon = false;
      }
    }
    if ( !precon ) {
      mrank = gsl_linalg_SV_decomp_jacobi (U, V, S); // A = U*S*V^t
    }
#else
    mrank = gsl_linalg_SV_decomp_jacobi (U, V, S); // A = U*S*V^t
#endif
  }
  else {
    for ( ZFSId i = 0; i < m; i++ ) {
      for ( ZFSId j = 0; j < n; j++ ) {
        V(j,i) = weights[i] * A(i,j);
      }
    }
    mrank = gsl_linalg_SV_decomp_jacobi (V, U, S); // A^t = V*S*U^t
  }

  if ( mrank < 0 ) {
    cerr << domainId() << " Jacobi SVD algorithm for matrix '" << A.m_variable_name << "' (" << m << "x" << n << ") did not reach desired tolerance" << endl;
    stringstream wgt;
    wgt << "weights: ";
    for ( ZFSId i = 0; i < m; i++ ) wgt << weights[i] << " ";
    cerr << wgt.str() << endl;
  }

  for ( ZFSId i = 0; i < n; i++ ) {
    for ( ZFSId j = 0; j < m; j++ ) {
      AInv(i,j) = F0;
      for ( ZFSId k = 0; k < p; k++ ) {
        if ( fabs(S(k)) < m_eps ) continue;
        AInv(i,j) += V(i,k) * U(j,k) / S(k);
      }
      AInv(i,j) *= weights[j];
    }
  }

  ZFSFloat condNum = S(0) / S(S.size()-1);
  if ( condNum > 100000.0 || std::isnan( condNum ) ) {
/*    cerr << "S: "; for ( ZFSId i = 0; i < p; i++ ) { cerr << S[i] << " "; } cerr << endl << endl;
      cerr << "Matrix '" << A.m_variable_name << "' (" << m << "x" << n << ") condition number: " << condNum
      << ", rank: " << mrank << ", regularized condition number: " << S(0) / S(mrank-1) << "." << endl;

      cerr << "weights: ";
      for ( ZFSId i = 0; i < m; i++ ) { cerr << weights[i] << " "; } cerr << endl << endl;
      cerr << "singular values: ";
      for ( ZFSId i = 0; i < p; i++ ) { cerr << S[i] << " "; } cerr << endl << endl;
      cerr << "U: " << m << "x" << p << endl;
      for ( ZFSId i = 0; i < m; i++ ) { for ( ZFSId j = 0; j < p; j++ ) { cerr << U(i,j) << " "; } cerr << endl; } cerr << endl << endl;
      cerr << "V: " << n << "x" << p << endl;
      for ( ZFSId i = 0; i < n; i++ ) { for ( ZFSId j = 0; j < p; j++ ) { cerr << V(i,j) << " "; } cerr << endl; } cerr << endl << endl;
      cerr << "A: " << m << "x" << n << endl;
      for ( ZFSId i = 0; i < m; i++ ) { for ( ZFSId j = 0; j < n; j++ ) { cerr << A(i,j) << " "; } cerr << endl; } cerr << endl << endl;
      cerr << "AInv: " << n << "x" << m << endl;
      for ( ZFSId i = 0; i < n; i++ ) { for ( ZFSId j = 0; j < m; j++ ) { cerr << AInv(i,j) << " "; } cerr << endl; } cerr << endl << endl;
*/
  }
  if ( mrank <= 0 ) {
    condNum = -F1;
  }
  else if ( mrank < p ) {
    condNum = S(0) / S(mrank-1); //spectral norm
    for ( ZFSId i = mrank; i < p; i++ ) {
      S[i] = F0;
    }
  }
  if ( mrank < p ) {
    cerr << "Matrix rank reduced to " << mrank << " in SVD for matrix " << A.m_variable_name << endl;
    for ( ZFSId i = 0; i < m; i++ ) {
      for ( ZFSId j = 0; j < n; j++ ) { cerr << A(i,j) << " "; } cerr << endl;
    }
    condNum = -F2;
  }

  return condNum;
}




/* linalg/svd.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007, 2010 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 */

/* This is a the jacobi version */
/* Author:  G. Jungman */

/*
 * Algorithm due to J.C. Nash, Compact Numerical Methods for
 * Computers (New York: Wiley and Sons, 1979), chapter 3.
 * See also Algorithm 4.1 in
 * James Demmel, Kresimir Veselic, "Jacobi's Method is more
 * accurate than QR", Lapack Working Note 15 (LAWN15), October 1989.
 * Available from netlib.
 *
 * Based on code by Arthur Kosowsky, Rutgers University
 *                  kosowsky@physics.rutgers.edu
 *
 * Another relevant paper is, P.P.M. De Rijk, "A One-Sided Jacobi
 * Algorithm for computing the singular value decomposition on a
 * vector computer", SIAM Journal of Scientific and Statistical
 * Computing, Vol 10, No 2, pp 359-371, March 1989.
 *
 */
ZFSId ZFSStrctrdBlck3D::gsl_linalg_SV_decomp_jacobi (ZFSFloatScratchSpace& A,
                                                     ZFSFloatScratchSpace& Q,
                                                     ZFSFloatScratchSpace& S)
{
  const ZFSFloat EPS = 1e-15;
  const ZFSFloat EPS_REG = 1e-12;
  const ZFSId M = A.size0();
  const ZFSId N = A.size1();
  ZFSId rank = N;
  if ( M < N ) {
    zfsTerm(1,__CALLING_FUNCTION__, "svd of MxN matrix, M<N, is not implemented. Just take the transpose of A before calling this function.");
  }
  else if (Q.size0() != A.size1()) {
    zfsTerm(1,__CALLING_FUNCTION__, "square matrix Q must match second dimension of matrix A");
  }
  else if (Q.size0() != Q.size1()) {
    zfsTerm(1,__CALLING_FUNCTION__, "matrix Q must be square");
  }
  else if (S.size0() != A.size1()) {
    zfsTerm(1,__CALLING_FUNCTION__, "length of vector S must match second dimension of matrix A");
  }
  else
  {
    ZFSId i, j, k;

    /* Initialize the rotation counter and the sweep counter. */
    ZFSId count = 1;
    ZFSId sweep = 0;
    //ZFSId sweepmax = 20*N;//5*N;
    ZFSId sweepmax = 200;
    //ZFSFloat tolerance = 10 * (ZFSFloat)M * EPS;
    ZFSFloat tolerance = 1e-12;

    /* Always do at least 12 sweeps. */
    sweepmax = zfsMAX(sweepmax, 12);

    /* Set Q to the identity matrix. */
    Q.fill(F0);
    for ( ZFSId t = 0; t < N; t++ ) Q(t,t) = F1;

    /* Store the column error estimates in S, for use during the
       orthogonalization */

    for (j = 0; j < N; j++)
    {
      // gsl_vector_view cj = gsl_matrix_column (A, j);
      //double sj = gsl_blas_dnrm2 (&cj.vector);
      //gsl_vector_set(S, j, GSL_DBL_EPSILON * sj);
      ZFSFloat sj = F0;
      for ( ZFSId t = 0; t < M; t++ ) sj += POW2( A(t,j) );
      S(j) = EPS * sqrt(sj);
    }

    /* Orthogonalize A by plane rotations. */

    while (count > 0 && sweep <= sweepmax)
    {
      /* Initialize rotation counter. */
      count = N * (N - 1) / 2;

      for (j = 0; j < N - 1; j++)
      {
        for (k = j + 1; k < N; k++)
        {
          ZFSFloat a = 0.0;
          ZFSFloat b = 0.0;
          ZFSFloat p = 0.0;
          ZFSFloat q = 0.0;
          ZFSFloat cosine, sine;
          ZFSFloat v;
          ZFSFloat abserr_a, abserr_b;
          ZFSId sorted, orthog, noisya, noisyb;

          //gsl_vector_view cj = gsl_matrix_column (A, j);
          //gsl_vector_view ck = gsl_matrix_column (A, k);
          //gsl_blas_ddot (&cj.vector, &ck.vector, &p);
          for ( ZFSId t = 0; t < M; t++ ) p += A(t,j)*A(t,k);

          p *= 2.0 ;  /* equation 9a:  p = 2 x.y */

          //a = gsl_blas_dnrm2 (&cj.vector);
          //b = gsl_blas_dnrm2 (&ck.vector);
          for ( ZFSId t = 0; t < M; t++ ) { a += POW2( A(t,j) ); } a = sqrt(a);
          for ( ZFSId t = 0; t < M; t++ ) { b += POW2( A(t,k) ); } b = sqrt(b);

          q = a * a - b * b;
          v = hypot(p, q);

          /* test for columns j,k orthogonal, or dominant errors */

          abserr_a = S(j);//gsl_vector_get(S,j);
          abserr_b = S(k);//gsl_vector_get(S,k);

          //sorted = !(a<b);
          sorted = a >= b;//(GSL_COERCE_DBL(a) >= GSL_COERCE_DBL(b));
          //orthog = fabs(p) < tolerance*a*b;
          orthog = fabs(p) <= tolerance*a*b; //(fabs (p) <= tolerance * GSL_COERCE_DBL(a * b));
          noisya = (a < abserr_a);
          noisyb = (b < abserr_b);

          if (sorted && (orthog || noisya || noisyb))
          {
            count--;
            continue;
          }

          /* calculate rotation angles */
          if (fabs(v) < 1e-15 || !sorted)
            //if (v == 0 || !sorted)
          {
            cosine = 0.0;
            sine = 1.0;
          }
          else
          {
            cosine = sqrt((v + q) / (2.0 * v));
            sine = p / (2.0 * v * cosine);
          }

          /* apply rotation to A */
          for (i = 0; i < M; i++)
          {
            const ZFSFloat Aik = A(i,k);//gsl_matrix_get (A, i, k);
            const ZFSFloat Aij = A(i,j);//gsl_matrix_get (A, i, j);
            A(i,j) = Aij * cosine + Aik * sine; // gsl_matrix_set (A, i, j, Aij * cosine + Aik * sine);
            A(i,k) = -Aij * sine + Aik * cosine; //gsl_matrix_set (A, i, k, -Aij * sine + Aik * cosine);
          }

          S(j) = fabs(cosine) * abserr_a + fabs(sine) * abserr_b; //gsl_vector_set(S, j, fabs(cosine) * abserr_a + fabs(sine) * abserr_b);
          S(k) = fabs(sine) * abserr_a + fabs(cosine) * abserr_b; //gsl_vector_set(S, k, fabs(sine) * abserr_a + fabs(cosine) * abserr_b);

          /* apply rotation to Q */
          for (i = 0; i < N; i++)
          {
            const ZFSFloat Qij = Q(i,j);//gsl_matrix_get (Q, i, j);
            const ZFSFloat Qik = Q(i,k);//gsl_matrix_get (Q, i, k);
            Q(i,j) = Qij * cosine + Qik * sine;//gsl_matrix_set (Q, i, j, Qij * cosine + Qik * sine);
            Q(i,k) =  -Qij * sine + Qik * cosine;//gsl_matrix_set (Q, i, k, -Qij * sine + Qik * cosine);
          }
        }
      }

      /* Sweep completed. */
      sweep++;
    }

    /*
     * Orthogonalization complete. Compute singular values.
     */

    {
      ZFSFloat prev_norm = -1.0;

      for (j = 0; j < N; j++)
      {
        //gsl_vector_view column = gsl_matrix_column (A, j);
        //ZFSFloat norm = gsl_blas_dnrm2 (&column.vector);
        ZFSFloat norm = F0;
        for ( ZFSId t = 0; t < M; t++ ) { norm += POW2( A(t,j) ); } norm = sqrt(norm);

        /* Determine if singular value is zero, according to the
           criteria used in the main loop above (i.e. comparison
           with norm of previous column). */
        if ( fabs(norm) < EPS_REG || fabs(prev_norm) < EPS_REG )
          //if (norm == 0.0 || prev_norm == 0.0 || (j > 0 && norm <= tolerance * prev_norm))
        {
          S(j) = F0;//gsl_vector_set (S, j, 0.0);     /* singular */
          for ( ZFSId t = 0; t < M; t++ ) A(t,j) = F0; //gsl_vector_set_zero (&column.vector);   /* annihilate column */

          rank = zfsMIN(rank, j);
          prev_norm = 0.0;
        }
        else
        {
          S(j) = norm;//gsl_vector_set (S, j, norm);    /* non-singular */
          for ( ZFSId t = 0; t < M; t++ ) A(t,j) /= norm; //gsl_vector_scale (&column.vector, 1.0 / norm);  /* normalize column */

          prev_norm = norm;
        }
      }
    }

    if (count > 0)
    {
      /* reached sweep limit */
      //zfsTerm(1,__CALLING_FUNCTION__, "Jacobi iterations did not reach desired tolerance");
      return -1;
    }

    return rank;
  }
  return -1;
}



/**
 * Create a QR factorization for a.
 * Copyright (c) 2008-2011 Zhang Ming (M. Zhang), zmjerry@163.com
 */
ZFSBool ZFSStrctrdBlck3D::QRdec( ZFSFloatScratchSpace& a,
                                 ZFSFloatScratchSpace& d,
                                 const ZFSId m, const ZFSId n )
{
  ASSERT( m >= n, "" );
  const ZFSId p = min(m,n);
  ZFSId rank = p;
  for( ZFSId k=0; k<p; ++k ) {
    ZFSFloat nrm = F0;
    for( ZFSId i=k; i<m; ++i ) {
      nrm = hypot( nrm, a(i,k) );
    }
    if ( fabs(nrm) > m_eps ) {
      if( a(k,k) < F0 ) nrm = -nrm;
      for( ZFSId i=k; i<m; ++i ) {
        a(i,k) /= nrm;
      }
      a(k,k) += F1;
      for( ZFSId j=k+1; j<n; ++j ) {
        ZFSFloat s = F0;
        for( ZFSId i=k; i<m; ++i ) {
          s += a(i,k)*a(i,j);
        }
        s = -s/a(k,k);
        for( ZFSId i=k; i<m; ++i ) {
          a(i,j) += s*a(i,k);
        }
      }
      d(k) = -nrm;
    }
    else {
      d(k) = F0;
      rank = zfsMIN( rank, k );
    }
  }
  //cerr << rank << " " << p << endl; cin.get();
  if ( rank < p ) {
    return false;
  }
  return true;
}


//-----------------------------------------------------------------------------


/**
 * \brief retrieves the 'Q' from the given QR factorization
 * \author Lennart Schneiders
 */
void ZFSStrctrdBlck3D::getQ( ZFSFloatScratchSpace& qr, ZFSFloatScratchSpace& q, const ZFSId m, const ZFSId n )
{
  ASSERT( m >= n, "" );
  const ZFSId p = min(m,n);
  for( ZFSId k=p-1; k>=0; --k ) {
    for( ZFSId i=0; i<m; ++i ) q(i,k) = F0;
    q(k,k) = F1;
    for( ZFSId j=k; j<p; ++j ) {
      ZFSFloat s = F0;
      for( ZFSId i=k; i<m; ++i ) s += qr(i,k) * q(i,j);
      s = -s / qr(k,k);
      for( ZFSId i=k; i<m; ++i ) q(i,j) += s*qr(i,k);
    }
  }
}


//----------------------------------------------------------------------------


/**
 * \brief retrieves the 'R' from the given QR factorization
 * \author Lennart Schneiders
 */
void ZFSStrctrdBlck3D::getR( ZFSFloatScratchSpace& qr, ZFSFloatScratchSpace& d, ZFSFloatScratchSpace& r, const ZFSId m, const ZFSId n )
{
  ASSERT( m >= n, "" );
  const ZFSId p = min(m,n);
  for( ZFSId i=0; i<p; ++i ) {
    for( ZFSId j=0; j<n; ++j ) {
      if( i < j ) r(i,j) = qr(i,j);
      else if( i == j ) r(i,j) = d(i);
      else r(i,j) = F0;
    }
  }
}


inline ZFSFloat ZFSStrctrdBlck3D::getPSI(ZFSId I, ZFSId dim) {
  const ZFSFloat FK = 18.0;
  const ZFSId IJK[3]={1,m_nCells[2],m_nCells[1]*m_nCells[2]};
  const ZFSId IP1=I+IJK[dim];
  const ZFSId IM1=I-IJK[dim];
  const ZFSId IP2=I+2*IJK[dim];

  const ZFSFloat PIM2 = m_cells->pvariables[PV->P][IM1];
  const ZFSFloat PIM1 = m_cells->pvariables[PV->P][I];
  const ZFSFloat PIP2 = m_cells->pvariables[PV->P][IP2];
  const ZFSFloat PIP1 = m_cells->pvariables[PV->P][IP1];

  const ZFSFloat PSI = zfsMIN(F1,FK*zfsMAX(zfsMAX(
                                           fabs((PIM2-PIM1)/zfsMIN(PIM2,PIM1)),
                                           fabs((PIM1-PIP1)/zfsMIN(PIM1,PIP1))),
                                           fabs((PIP1-PIP2)/zfsMIN(PIP1,PIP2))));
  return PSI;
}
