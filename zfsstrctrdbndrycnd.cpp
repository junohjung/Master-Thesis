#include "zfsstrctrdbndrycnd.h"

using namespace std;

template <ZFSInt nDim>
ZFSStrctrdBndryCnd<nDim>::ZFSStrctrdBndryCnd(ZFSStrctrdBlck<nDim>* block)
  :m_zfsStrctrdComm(block->m_zfsStrctrdComm),
   m_nCells(block->m_nCells),
   m_nPoints(block->m_nPoints),
   m_cells(block->m_cells),
   CV(block->CV),
   PV(block->PV),
   FQ(block->FQ),
   m_blockId(block->m_blockId),
   m_noGhostLayers(block->m_noGhostLayers),
   m_physicalBCMap(block->m_windowInfo->physicalBCMap),
   m_auxDataMap(block->m_windowInfo->physicalAuxDataMap),
   m_globalStrctrdBndryMaps(block->m_windowInfo->globalStrctrdBndryCndMaps),
   m_noSpongeDomainInfos(block->m_noSpongeDomainInfos),
   m_spongeBcWindowInfo(block->m_spongeBcWindowInfo),
   m_spongeLayerType(block->m_spongeLayerType),
   m_spongeLayerThickness(block->m_spongeLayerThickness),
   m_sigmaSponge(block->m_sigmaSponge),
   m_betaSponge(block->m_betaSponge),
   m_targetDensityFactor(block->m_targetDensityFactor),
   m_noStrctrdCells(block->m_noStrctrdCells),
   m_sutherlandPlusOne(block->m_sutherlandPlusOne),
   m_sutherlandConstant(block->m_sutherlandConstant),
   m_bCfCpCoeff(block->m_bCfCpCoeff),
   m_bPower(block->m_bPower),
   m_bCl(block->m_bCl),
   m_bCd(block->m_bCd)
{
  TRACE();


  m_coordinates=new ZFSFloat*[nDim];
  for(ZFSId i=0; i<nDim; i++) {
    m_coordinates[i]=&(block->m_coordinates[i][0]);
  }

  m_block=block;

  //the surfaces for the channel flow calculation
  m_channelSurfaceIn=F0;
  m_channelSurfaceOut=F0;
}

template <ZFSInt nDim>
void ZFSStrctrdBndryCnd<nDim>::applyDirichletNeumannBC()
{
  for(ZFSId bcId=0; bcId < m_noBndryCndIds; bcId++) {
    (this->*bndryCndHandler[bcId])(bcId);
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBndryCnd<nDim>::applyNonReflectingBC()
{
  TRACE();  
  //not implemented yet
  //for( ZFSId bcId = 0;  bcId < m_noBndryCndIds;  bcId++ )
  //  (this->*nonReflectingBoundaryCondition[bcId]) (bcId);
}

template <ZFSInt nDim>
void ZFSStrctrdBndryCnd<nDim>::mapCpy(ZFSStrctrdWindowMap* input, ZFSStrctrdWindowMap* output)
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
  output->face=input->face;
}

template <ZFSInt nDim>
void ZFSStrctrdBndryCnd<nDim>::assignBndryCnds()
{
  TRACE();

  m_noBndryCndIds=m_physicalBCMap.size();

  bndryCndHandler=NULL;
  bndryCndHandler=new BndryCndHandler [m_physicalBCMap.size()];
  initBndryCndHandler= new BndryCndHandler [m_physicalBCMap.size()];

  //relation between surface index map and the bcId
  zfsAlloc(m_channelSurfacIndexMap, m_physicalBCMap.size(), "m_channelSurfacIndexMap", -1, __CALLING_FUNCTION__);
  ZFSInt counter=0;
  //assign the function pointers
  for(ZFSUint bcId=0; bcId<m_physicalBCMap.size(); bcId++) {
    switch(m_physicalBCMap[bcId]->BC){
    case 0:{
      // euler wall
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc0;//empty boundary condition
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc0; //empty boundary condition
      break;
    }
    case 1000:
    case 1004: {
      // adiabatic wall
      if(m_block->m_movingGrid) {
        bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc1004;
        initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc1004;
      } else {
        bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc1000;
        initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc1000;
      }
      break;
    }
    case 1003:
    case 1006: {
      // isothermal wall
      if(m_block->m_movingGrid) {
        bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc1006;
        initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc1006;
      } else {
        bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc1003;
        initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc1003;
      }
      break;
    }
    case 1001:{
      // euler wall
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc1001;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc1001;
      break;
    }
    case 1007:{
      // oscillating wall
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc1007;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc1007;
      break;
    }
    case 2001:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2001;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2001;
      break;
    }
    case 2002:
    case 2010:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2002;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2002;
      break;
    }
    case 2003:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2003;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2003;
      break;
    }
    case 2004:
    case 2024:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2004;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2004;
      break;
    }
    case 2005:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2005;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2005;
      break;
    }
    case 2007:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2007;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2007;
      break;
    }
    case 2009:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2009;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2009;
      break;
    }
    case 2012:{ //characteristic inflow
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2012;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2012;
      break;
    }
    case 2013:{ //characteristic outflow
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2013;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2013;
      break;
    }
    case 2014:{ //characteristic outflow
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2014;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2014;
      break;
    }
    case 2020:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2020;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2020;
      break;
    }
    case 2021:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2021;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2021;
      break;
    }
    case 2099:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2099;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2099;
      break;
    }
    case 2221:{        //junoh //zonal with STG
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2221;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2221;
      break;
    }
    case 2222:{        //junoh //zonal without STG
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2222;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2222;
      break;
    }
    case 2199:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2199;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2199;
      break;
    }
    case 2300:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2300;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2300;
      counter++;
      break;
    }
    case 2401:
    case 2402:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2402;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2402;
      m_channelSurfacIndexMap[bcId]=counter;
      counter++;
      break;
    }
    case 2500:{//Rescaling
      if(m_block->m_rans) {
        bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2510;
        initBndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::initBc2510;
      } else {
        bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2500;
        initBndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::initBc2500;
      }
      break;
    }
    case 2501:{//Rescaling
      if(m_block->m_rans) {
        bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2511;
        initBndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::initBc2501;
      } else {
        bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2501;
        initBndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::initBc2501;
      }
      break;
    }
    case 2600:{//Prescribing profile
      m_block->m_bc2600 = true;
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2600;
      initBndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::initBc2600;
      m_block->m_bc2600 = true;
      break;
    }
    case 2601:{//Prescribing profile
      m_block->m_bc2601 = true;
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2601;
      initBndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::initBc2601;
      m_block->m_bc2601 = true;
      break;
    }
    case 2700:{//mode inflow
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2700;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc2700;
      break;
    }
    case 2900:{//Jet Freund Inlet
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc2900;
      initBndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::initBc2900;
      break;
    }
    case 3000:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc3000;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc3000;
      break;
    }
    case 4001:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc4001;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc4001;
      break;
    }
    case 4002:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc4001;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc4001;
      break;
    }
    case 6000:{
      bndryCndHandler[bcId]=&ZFSStrctrdBndryCnd::bc6000;
      initBndryCndHandler[bcId]= &ZFSStrctrdBndryCnd::initBc6000;
      break;
    }
    case 7909:
      {
        bndryCndHandler[bcId] = &ZFSStrctrdBndryCnd::bc7909;
        initBndryCndHandler[bcId] = &ZFSStrctrdBndryCnd::initBc7909;
        break;
      }
      //empty BC!!!!!!!!!!!!!!!!
    case 4401:
    case 4402:
    case 4403:
    case 4404:
    case 4405:
    case 4406:
      {
        bndryCndHandler[bcId] = &ZFSStrctrdBndryCnd::bc9999;
        initBndryCndHandler[bcId] = &ZFSStrctrdBndryCnd::initBc9999;
        break;
      }

    default:{
      cout << "boundary condtition is missing" <<m_physicalBCMap[bcId]->BC << endl;
      zfsTerm(1, __CALLING_FUNCTION__, "Boundary Condition is not implemented" );
      break;
    }
    }
  }
}

template <ZFSInt nDim>
void ZFSStrctrdBndryCnd<nDim>::correctBndryCndIndices()
{
  //in correcting cell Information
  for(ZFSId bcId=0; bcId < m_noBndryCndIds; bcId++)
    {
      (this->*initBndryCndHandler[bcId])(bcId);
    }
}

template <ZFSInt nDim>
ZFSStrctrdBndryCnd<nDim>::~ZFSStrctrdBndryCnd()
{
  delete[] bndryCndHandler;
  delete[] initBndryCndHandler;
  delete[] m_coordinates;
}

template <ZFSInt nDim>
void ZFSStrctrdBndryCnd<nDim>::saveAuxData(){

}

// Explicit instantiations for 2D and 3D
template class ZFSStrctrdBndryCnd<2>;
template class ZFSStrctrdBndryCnd<3>;
