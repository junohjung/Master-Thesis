#include "zfsstrctrdblck2drans.h"
#include "zfsglobals.h"
#include "zfsiolib.h"

#include "zfsstrctrdblckpartition.h"
#include "zfsstrctrdblckwindowinfo.h"
#include <cstdlib>
#if not defined(ZFS_MS_COMPILER)
#include <unistd.h>
#endif
//temporaray
#include <vector>


ZFSStrctrdBlck2DRans::ZFSStrctrdBlck2DRans(class ZFSStrctrdBlck2D* block ):
  m_zfsStrctrdComm(block->m_zfsStrctrdComm),
  m_blockId(block->m_blockId),
  m_nCells(block->m_nCells),
  m_nPoints(block->m_nPoints),
  m_noStrctrdCells(block->m_noStrctrdCells),
  m_cells(block->m_cells),
  m_coordinates(block->m_coordinates),
  CV(block->CV),
  PV(block->PV),
  FQ(block->FQ),
  m_noGhostLayers(block->m_noGhostLayers),
  m_eps(block->m_eps),
  m_chi(block->m_chi),
  m_sutherlandConstant(block->m_sutherlandConstant),
  m_sutherlandPlusOne(block->m_sutherlandPlusOne)
{
  m_block=block;
  //set all the methods for the 2d RANS code here

  initFluxMethod();
}

ZFSStrctrdBlck2DRans::~ZFSStrctrdBlck2DRans(){
}

void ZFSStrctrdBlck2DRans::initFluxMethod() {
  //set pointer to the right AUSM for the RANS equations;

  /*! \page propertyPage1
    \section ransMethod
    <code>ZFSInt ZFSStrctrdBlck3DRans::m_ransMethod </code>\n
    default = <code> 1.0 </code>\n \n
    Name of the RANS method to be used.\n
    Possible values are:\n
    <ul>
    <li>RANS_SA_DV</li>
    </ul>
    Keywords: <i>RANS, STRCTRD</i>
  */
  m_ransMethod = *(ZFSContext::getProperty("ransMethod", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL)->asString(0));
  switch (string2enum(m_ransMethod))
    {
    case RANS_SA :
      {
        zfsAlloc(m_cells->saFlux1, nDim*m_noStrctrdCells, "m_cells->saFlux1", -999999.9, __CALLING_FUNCTION__);
        zfsAlloc(m_cells->saFlux2, nDim*m_noStrctrdCells, "m_cells->saFlux2", -999999.9, __CALLING_FUNCTION__);
        zfsAlloc(m_cells->prodDest, m_noStrctrdCells, "m_cells->prodDest", -999999.9, __CALLING_FUNCTION__);
        viscFluxMethod = &ZFSStrctrdBlck2DRans::viscousFlux_SA;
        switch(m_block->CV->noVariables){
        case 5: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm, 5>; break;}
        case 6: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm, 6>; break;}
        case 7: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm, 7>; break;}
        default:{
          stringstream errorMessage;
          errorMessage << "Number of Variables " << m_block->CV->noVariables << " not implemented! in temlate Rans AUSM " << endl;
          zfsTerm(1, __CALLING_FUNCTION__, errorMessage.str());
        }
        }
        break;
      }
    case RANS_SA_DV :
      {
        zfsAlloc(m_cells->saFlux1, nDim*m_noStrctrdCells, "m_cells->saFlux1", -999999.9, __CALLING_FUNCTION__);
        zfsAlloc(m_cells->saFlux2, nDim*m_noStrctrdCells, "m_cells->saFlux2", -999999.9, __CALLING_FUNCTION__);
        zfsAlloc(m_cells->prodDest, m_noStrctrdCells, "m_cells->prodDest", -999999.9, __CALLING_FUNCTION__);
        viscFluxMethod = &ZFSStrctrdBlck2DRans::viscousFlux_SA;
        switch(m_block->CV->noVariables){
        case 5: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::AusmDV, 5>; break;}
        case 6: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::AusmDV, 6>; break;}
        case 7: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::AusmDV, 7>; break;}
        default:{
          stringstream errorMessage;
          errorMessage << "Number of Variables " << m_block->CV->noVariables << " not implemented! in temlate Rans AUSM " << endl;
          zfsTerm(1, __CALLING_FUNCTION__, errorMessage.str());
        }
        }
        break;
      }
    case RANS_FS :
      {
        viscFluxMethod = &ZFSStrctrdBlck2DRans::viscousFlux_FS;
        switch(m_block->CV->noVariables){
        case 5: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm, 5>; break;}
        case 6: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm, 6>; break;}
        case 7: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm, 7>; break;}
        default:{
          stringstream errorMessage;
          errorMessage << "Number of Variables " << m_block->CV->noVariables << " not implemented! in temlate Rans AUSM " << endl;
          zfsTerm(1, __CALLING_FUNCTION__, errorMessage.str());
        }
        }
        break;
      }
    case RANS_KOMEGA :
      {
        viscFluxMethod = &ZFSStrctrdBlck2DRans::viscousFlux_KOmega;
        switch(m_block->CV->noVariables){
        case 5: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm, 5>; break;}
        case 6: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm, 6>; break;}
        case 7: { reconstructSurfaceData=&ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm, 7>; break;}
        default:{
          stringstream errorMessage;
          errorMessage << "Number of Variables " << m_block->CV->noVariables << " not implemented! in temlate Rans AUSM " << endl;
          zfsTerm(1, __CALLING_FUNCTION__, errorMessage.str());
        }
        }
        break;
      }
    default:
      { zfsTerm(1, __CALLING_FUNCTION__, "RANS METHOD wsa not specified in properties"); break; }
    }
}

void ZFSStrctrdBlck2DRans::Muscl(){
  (this->*reconstructSurfaceData)();
}

// void ZFSStrctrdBlck2DRans::AusmRANS(ZFSFloat* QLeft, ZFSFloat* QRight, const ZFSId dim, const ZFSId cellId){
//   (this->*convFluxMethod)(QLeft,QRight, dim, cellId);
// }

void ZFSStrctrdBlck2DRans::Ausm(ZFSFloat* QLeft, ZFSFloat*  QRight, const ZFSId dim, const ZFSId I){
  ZFSFloat pFactor[2]={F0,F0};
  const ZFSFloat gamma = m_block->m_gamma;
  const ZFSFloat gammaMinusOne = gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_cells->surfaceMetrics[ I ]);

  const ZFSFloat dxdtau = m_cells->dxt[dim][I];

  //calculate pressure
  const ZFSFloat PL= QLeft[PV->P];
  const ZFSFloat UL   = QLeft[ PV->U ];
  const ZFSFloat VL   = QLeft[ PV->V ];
  const ZFSFloat RHOL = QLeft[ PV->RHO ];

  const ZFSFloat PR= QRight[PV->P];
  const ZFSFloat UR   = QRight[ PV->U ];
  const ZFSFloat VR   = QRight[ PV->V ];
  const ZFSFloat RHOR = QRight[ PV->RHO ];

  // compute lenght of metric vector for normalization
  const ZFSFloat metricLength = sqrt(POW2( surf[dim*2+0] ) + POW2(surf[dim*2+1]));
  const ZFSFloat fMetricLength = F1 / metricLength;

  //scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const ZFSFloat UUL = ((UL * surf[ dim*2+0 ] +
                         VL * surf[ dim*2+1 ]) - dxdtau) * fMetricLength;


  const ZFSFloat UUR = ((UR * surf[ dim*2+0 ] +
                         VR * surf[ dim*2+1 ]) - dxdtau) * fMetricLength;


  //speed of sound
  const ZFSFloat AL = sqrt(gamma * zfsMAX(m_eps, PL / zfsMAX(m_eps, RHOL)));
  const ZFSFloat AR = sqrt(gamma * zfsMAX(m_eps, PR / zfsMAX(m_eps, RHOR)));

  const ZFSFloat MAL = UUL / AL;
  const ZFSFloat MAR = UUR / AR;

  const ZFSFloat MALR = F1B2*(MAL+MAR);
  const ZFSFloat PLR = F1B2*(PL+PR);

  const ZFSFloat RHO_AL = RHOL*AL;
  const ZFSFloat RHO_AR = RHOR*AR;

  const ZFSFloat PLfRHOL = PL/RHOL;
  const ZFSFloat PRfRHOR = PR/RHOR;

  //sa-transport variable
  const ZFSFloat nutildeR = QRight[PV->RANS_FIRST];
  const ZFSFloat nutildeL = QLeft[PV->RANS_FIRST];

  const ZFSFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + PLfRHOL;
  const ZFSFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + PRfRHOR;

  const ZFSFloat RHOU = F1B2 * ( MALR * (RHO_AL+RHO_AR) + fabs(MALR) * (RHO_AL-RHO_AR) ) * metricLength;
  const ZFSFloat RHOU2 = F1B2*RHOU;
  // multiply by metric length to take surface area into account
  const ZFSFloat AbsRHO_U2 = fabs(RHOU2);

  // setup pressure factors, include metric terms (not normalized due to needed surface area)
  for( ZFSInt isd = xsd; isd < nDim; isd ++ ) {
    pFactor[ isd ] = surf[ 2 * dim + isd ];
  }

  const ZFSInt noCells = m_noStrctrdCells;
  ZFSFloat *const RESTRICT flux = m_cells->flux;

  flux[I+noCells*CV->RHO_U] = RHOU2 * ( UL + UR ) + AbsRHO_U2 * ( UL - UR ) + PLR * pFactor[ 0 ];
  flux[I+noCells*CV->RHO_V] = RHOU2 * ( VL + VR ) + AbsRHO_U2 * ( VL - VR ) + PLR * pFactor[ 1 ];
  flux[I+noCells*CV->RHO_E] = RHOU2 * ( e0 + e1)  + AbsRHO_U2 * ( e0 - e1 ) + PLR * dxdtau;
  flux[I+noCells*CV->RHO]   = RHOU;
  flux[I+noCells*CV->RANS_FIRST] = RHOU2 * (nutildeL+ nutildeR ) + AbsRHO_U2 * ( nutildeL- nutildeR);
}

void ZFSStrctrdBlck2DRans::AusmDV(ZFSFloat* QLeft, ZFSFloat*  QRight, const ZFSId dim, const ZFSId I)
{
  ZFSFloat pFactor[2]={F0,F0};
  const ZFSFloat gamma = m_block->m_gamma;
  const ZFSFloat gammaMinusOne = gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_block->m_cells->surfaceMetrics[I]);

  //left side
  const ZFSFloat RHOL = QLeft[PV->RHO];
  const ZFSFloat FRHOL = F1/RHOL;
  ZFSFloat UL = QLeft[PV->U];
  ZFSFloat VL = QLeft[PV->V];
  const ZFSFloat PL = QLeft[PV->P];
  const ZFSFloat NUTILDEL = QLeft[PV->RANS_VAR[0]];

  //right side
  const ZFSFloat RHOR = QRight[PV->RHO];
  const ZFSFloat FRHOR = F1/RHOR;
  ZFSFloat UR = QRight[PV->U];
  ZFSFloat VR = QRight[PV->V];
  const ZFSFloat PR = QRight[PV->P];
  const ZFSFloat NUTILDER = QRight[PV->RANS_VAR[0]];

  const ZFSFloat PLfRHOL = PL/RHOL;
  const ZFSFloat PRfRHOR = PR/RHOR;
  const ZFSFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + PLfRHOL;
  const ZFSFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + PRfRHOR;


  // compute lenght of metric vector for normalization
  const ZFSFloat DGRAD = sqrt(POW2( surf[dim*2+0] ) + POW2(surf[dim*2+1]));
  const ZFSFloat FDGRAD = F1 / DGRAD;

  //scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const ZFSFloat UUL = ((UL * surf[ dim*2+0 ] +
                         VL * surf[ dim*2+1 ])) * FDGRAD;


  const ZFSFloat UUR = ((UR * surf[ dim*2+0 ] +
                         VR * surf[ dim*2+1 ])) * FDGRAD;

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
  ZFSFloat SV = 2.0*DGRAD/(m_cells->cellJac[I]+m_cells->cellJac[IP1])*
    (FDV+(F1-FDV)*getPSI(I,dim));
  const ZFSFloat SV1 = SV*DXDXEZ;
  const ZFSFloat SV2 = SV*DYDXEZ;

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

  // setup pressure factors, include metric terms (not normalized due to needed surface area)
  for( ZFSInt isd = xsd; isd < nDim; isd ++ ) {
    pFactor[ isd ] = surf[ 2 * dim + isd ];
  }

  const ZFSInt noCells = m_noStrctrdCells;
  ZFSFloat *const RESTRICT flux = m_cells->flux;

  flux[I+noCells*CV->RHO_U] = RHOU2 * ( UL + UR) + ARHOU2 * ( UL - UR) + PLR*pFactor[0] + RHOUL*UUL2+RHOUR*UUR2;
  flux[I+noCells*CV->RHO_V] = RHOU2 * ( VL + VR) + ARHOU2 * ( VL - VR) + PLR*pFactor[1] + RHOUL*UUL3+RHOUR*UUR3;
  flux[I+noCells*CV->RHO_E] = RHOU2 * ( e0 + e1) + ARHOU2 * ( e0 - e1);
  flux[I+noCells*CV->RHO]   = RHOU;
  flux[I+noCells*CV->RANS_VAR[0]] = RHOU2*(NUTILDEL+NUTILDER) + ARHOU2*(NUTILDEL-NUTILDER);
}

void ZFSStrctrdBlck2DRans::viscousFluxRANS(){
  (this->*viscFluxMethod)();
}


void ZFSStrctrdBlck2DRans::viscousFlux_SA(){
  computeTurbViscosity();

  //OTHER variables required to calculate the laminar viscous fluxes
  const ZFSInt noCells = m_noStrctrdCells;
  const ZFSFloat rPrL = F1/m_block->m_Pr;
  const ZFSFloat rPrT = F1/0.9;
  const ZFSFloat rRe = F1 / m_block->m_Re0;
  const ZFSFloat gammaMinusOne = m_block->m_gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;


  ZFSFloat* __restrict eflux=&m_cells->eFlux[0];
  ZFSFloat* __restrict fflux=&m_cells->fFlux[0];
  ZFSFloat* __restrict vflux=&m_cells->viscousFlux[0];
  ZFSFloat* __restrict sa_1flux=&m_cells->saFlux1[0];
  ZFSFloat* __restrict sa_2flux=&m_cells->saFlux2[0];

  ZFSFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  ZFSFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  ZFSFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  ZFSFloat* const RESTRICT nuTilde=&m_cells->pvariables[PV->RANS_VAR[0]][0];
  ZFSFloat* const RESTRICT T = &m_cells->temperature[0];
  ZFSFloat* __restrict mut = &m_cells->fq[FQ->MU_T][0];

  for (ZFSId j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers+1; j++) {
    for (ZFSId i=m_noGhostLayers-1; i<m_nCells[1]-m_noGhostLayers+1; i++) {
      //get the adjacent cells;

      const ZFSId IJ   = cellIndex(i,j);
      const ZFSId IPJ  = cellIndex((i+1),j);
      const ZFSId IPJP = cellIndex((i+1),(j+1));
      const ZFSId IJP  = cellIndex(i,(j+1));
      const ZFSId IMJ  = cellIndex((i-1),j);
      const ZFSId IJM  = cellIndex(i,(j-1));

      const ZFSFloat cornerMetrics[9] = {m_cells->cornerMetrics[0*noCells+IJ],
                                         m_cells->cornerMetrics[1*noCells+IJ],
                                         m_cells->cornerMetrics[2*noCells+IJ],
                                         m_cells->cornerMetrics[3*noCells+IJ]};

      const ZFSFloat dudxi=F1B2*(u[IPJP]+u[IPJ]-u[IJP]-u[IJ]);
      const ZFSFloat dudet=F1B2*(u[IPJP]+u[IJP]-u[IPJ]-u[IJ]);

      const ZFSFloat dvdxi=F1B2*(v[IPJP]+v[IPJ]-v[IJP]-v[IJ]);
      const ZFSFloat dvdet=F1B2*(v[IPJP]+v[IJP]-v[IPJ]-v[IJ]);

      const ZFSFloat dTdxi=F1B2*(T[IPJP]+T[IPJ]-T[IJP]-T[IJ]);
      const ZFSFloat dTdet=F1B2*(T[IPJP]+T[IJP]-T[IPJ]-T[IJ]);

      const ZFSFloat dnutdxi=F1B2*(nuTilde[IPJP]+nuTilde[IPJ]-nuTilde[IJP]-nuTilde[IJ]);
      const ZFSFloat dnutdet=F1B2*(nuTilde[IPJP]+nuTilde[IJP]-nuTilde[IPJ]-nuTilde[IJ]);

      const ZFSFloat uAvg=F1B4*(u[IJP]+u[IPJP]+u[IJ]+u[IPJ]);
      const ZFSFloat vAvg=F1B4*(v[IJP]+v[IPJP]+v[IJ]+v[IPJ]);
      const ZFSFloat nutldAvg=F1B4*(nuTilde[IJP]+nuTilde[IPJP]+nuTilde[IJ]+nuTilde[IPJ]);

      const ZFSFloat nuAvg=F1B4*((zfsSUTHERLANDLAW(abs(T[IJ]))/rho[IJ])
                                 +(zfsSUTHERLANDLAW(abs(T[IPJ]))/rho[IPJ])
                                 +(zfsSUTHERLANDLAW(abs(T[IJP]))/rho[IJP])
                                 +(zfsSUTHERLANDLAW(abs(T[IPJP]))/rho[IPJP]));

      const ZFSFloat muAvg=F1B4*(zfsSUTHERLANDLAW(T[IJP])+zfsSUTHERLANDLAW(T[IPJP])+zfsSUTHERLANDLAW(T[IJ])+zfsSUTHERLANDLAW(T[IPJ]));

      const ZFSFloat mutAvg = F1B4*(mut[IPJP]+mut[IPJ]+mut[IJP]+mut[IJ]);


      const ZFSFloat dnutldx = dnutdxi * m_cells->cornerMetrics[(xsd * 2 + xsd)*noCells+IJ ] +
                               dnutdet * m_cells->cornerMetrics[(ysd * 2 + xsd)*noCells+IJ ];

      const ZFSFloat dnutldy = dnutdxi * m_cells->cornerMetrics[(xsd * 2 + ysd)*noCells+IJ ] +
                               dnutdet * m_cells->cornerMetrics[(ysd * 2 + ysd)*noCells+IJ ];

      const ZFSFloat Frj = rRe/m_cells->cornerJac[IJ];

      const ZFSFloat sax1 = Frj*(nuAvg+(1.0+m_sa.cb2)*nutldAvg)*
        (dnutldx * cornerMetrics[ xsd * 2 + xsd ]+
         dnutldy * cornerMetrics[ xsd * 2 + ysd ]);
      const ZFSFloat sax2 = -Frj*m_sa.cb2*
        (dnutldx * cornerMetrics[ xsd * 2 + xsd ]+
         dnutldy * cornerMetrics[ xsd * 2 + ysd ]);
      const ZFSFloat say1 = Frj*(nuAvg+(1.0+m_sa.cb2)*nutldAvg)*
        (dnutldx * cornerMetrics[ ysd * 2 + xsd ]+
         dnutldy * cornerMetrics[ ysd * 2 + ysd ]);
      const ZFSFloat say2 = -Frj*m_sa.cb2*
        (dnutldx * cornerMetrics[ ysd * 2 + xsd ]+
         dnutldy * cornerMetrics[ ysd * 2 + ysd ]);



      //compute vorticity
      const ZFSFloat du1 = u[IPJ]- u[IMJ];
      const ZFSFloat du2 = u[IJP]- u[IJM];
      const ZFSFloat dv1 = v[IPJ]- v[IMJ];
      const ZFSFloat dv2 = v[IJP]- v[IJM];
      const ZFSFloat vortk = (m_cells->cellMetrics[IJ][ xsd * 2 + xsd ]*dv1)+
                             (m_cells->cellMetrics[IJ][ ysd * 2 + xsd ]*dv2)-
                             (m_cells->cellMetrics[IJ][ xsd * 2 + ysd ]*du1)-
                             (m_cells->cellMetrics[IJ][ ysd * 2 + ysd ]*du2);
      const ZFSFloat s = F1B2*fabs(vortk)/m_cells->cellJac[IJ];


      const ZFSFloat distance = m_cells->fq[FQ->WALLDISTANCE][IJ];
      const ZFSFloat Fdist2 = 1.0/(distance*distance);
      const ZFSFloat chi = nuTilde[IJ]*rho[IJ]/(zfsSUTHERLANDLAW(T[IJ]));
      const ZFSFloat chip3 = chi*chi*chi;
      const ZFSFloat Fv1 =  chip3/(chip3+m_sa.cv1to3);
      const ZFSFloat Fv2 = F1 - (chi/(F1+chi*Fv1));

      const ZFSFloat term = nuTilde[IJ]*Fdist2*m_sa.Fkap2;
      const ZFSFloat stilde = s+term*Fv2*rRe;
      const ZFSFloat r = min(10.0, rRe*term/stilde);

      const ZFSFloat g = r+m_sa.cw2*(pow(r,6)-r);
      const ZFSFloat Fwterm = (1+m_sa.cw3to6)/(pow(g,6)+m_sa.cw3to6);
      const ZFSFloat Fw = g*pow(Fwterm,(1.0/6.0));
      const ZFSFloat prodValue = rho[IJ]*m_sa.cb1*(F1-m_sa.Ft2)*stilde*nuTilde[IJ];
      const ZFSFloat destValue = rRe* rho[IJ]*(m_sa.cw1*Fw-m_sa.cb1*m_sa.Fkap2*m_sa.Ft2)*pow( nuTilde[IJ],2.0)*Fdist2;


      m_cells->prodDest[IJ] = (prodValue-destValue)*m_cells->cellJac[IJ];

      // compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

      // tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
      //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
      //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
      ZFSFloat tau1 = F4B3 * ( dudxi * cornerMetrics[ xsd * 2 + xsd ]   +
                               dudet * cornerMetrics[ ysd * 2 + xsd ] ) -
        F2B3 * ( dvdxi * cornerMetrics[ xsd * 2 + ysd ]   +
                 dvdet * cornerMetrics[ ysd * 2 + ysd ] );

      // compute tau2 = du/dy + dv/dx

      // tau_xy = du/dxi * dxi/dy + du/deta * deta/dy + du/dzeta * dzeta/dy
      //        + dv/dxi * dxi/dx + dv/deta * deta/dx + dv/dzeta * dzeta/dx
      ZFSFloat tau2 = dudxi * cornerMetrics[ xsd * 2 + ysd ] +
        dudet * cornerMetrics[ ysd * 2 + ysd ] +
        dvdxi * cornerMetrics[ xsd * 2 + xsd ] +
        dvdet * cornerMetrics[ ysd * 2 + xsd ];

      // compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

      // tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
      //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
      //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
      ZFSFloat tau4 = F4B3 * ( dvdxi * cornerMetrics[ xsd * 2 + ysd ] +
                               dvdet * cornerMetrics[ ysd * 2 + ysd ] ) -
        F2B3 * ( dudxi * cornerMetrics[ xsd * 2 + xsd ] +
                 dudet * cornerMetrics[ ysd * 2 + xsd ]);


      const ZFSFloat dTdx = dTdxi * cornerMetrics[ xsd * 2 + xsd ] +
        dTdet * cornerMetrics[ ysd * 2 + xsd ];

      const ZFSFloat dTdy = dTdxi * cornerMetrics[ xsd * 2 + ysd ] +
        dTdet * cornerMetrics[ ysd * 2 + ysd ];



      const ZFSFloat muOverRe = (muAvg+mutAvg) * rRe / m_cells->cornerJac[IJ]; // divide by Jacobian
      tau1*=muOverRe;
      tau2*=muOverRe;
      tau4*=muOverRe;

      const ZFSFloat muH = rRe*FgammaMinusOne/ m_cells->cornerJac[IJ]*((muAvg*rPrL)+(mutAvg*rPrT));

      const ZFSFloat qx=muH*dTdx+uAvg*tau1+vAvg*tau2;
      const ZFSFloat qy=muH*dTdy+uAvg*tau2+vAvg*tau4;


      //efluxes
      eflux[ 5 * IJ ]     = tau1 * cornerMetrics[ xsd * 2 + xsd ] +
        tau2 * cornerMetrics[ xsd * 2 + ysd ];
      eflux[ 5 * IJ + 1 ] = tau2 * cornerMetrics[ xsd * 2 + xsd ] +
        tau4 * cornerMetrics[ xsd * 2 + ysd ];
      eflux[ 5 * IJ + 2 ] = qx * cornerMetrics[ xsd * 2 + xsd ] +
        qy * cornerMetrics[ xsd * 2 + ysd ];
      eflux[ 5* IJ + 3 ]=  sax1;// diffusion of nutilde for every cell
      eflux[ 5* IJ + 4 ]=  sax2;// diffusion of nutilde for every cell


      //ffluxes
      fflux[ 5 * IJ ]     = tau1 * cornerMetrics[ ysd * 2 + xsd ] +
        tau2 * cornerMetrics[ ysd * 2 + ysd ];
      fflux[ 5 * IJ + 1 ] = tau2 * cornerMetrics[ ysd * 2 + xsd ] +
        tau4 * cornerMetrics[ ysd * 2 + ysd ];
      fflux[ 5 * IJ + 2 ] = qx * cornerMetrics[ ysd * 2 + xsd ] +
        qy * cornerMetrics[ ysd * 2 + ysd ];
      fflux[ 5* IJ + 3 ]= say1;// diffusion of nutilde for every cell
      fflux[ 5* IJ + 4 ]= say2;// diffusion of nutilde for every cell
    }
  }


  for(ZFSId var = 0; var < ( CV->noVariables ) - 2; ++var) {
    ////////////////////////////////////////
    /////////// UVT ////////////////////////
    ////////////////////////////////////////
    for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; ++j){
      for(ZFSId i=m_noGhostLayers-1; i<m_nCells[1]-m_noGhostLayers; ++i) {
        const ZFSId IJ    = cellIndex(i,j);
        const ZFSId IJM   = cellIndex(i,(j-1));
        vflux[ 2 * IJ ] = F1B2 * ( eflux[ 5 * IJ + var ] + eflux[ 5 * IJM + var ]);
      }
    }

    for(ZFSId j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers; ++j) {
      for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; ++i) {
        const ZFSId IJ    = cellIndex(i,j);
        const ZFSId IMJ   = cellIndex((i-1),j);
        vflux[ 2 * IJ + 1 ] = F1B2 * ( fflux[ 5 * IJ + var ] + fflux[ 5 * IMJ + var ]);
      }
    }

    for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; ++j) {
      for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; ++i) {
        const ZFSId IJ    = cellIndex(i,j);
        const ZFSId IMJ   = cellIndex((i-1),j);
        const ZFSId IJM   = cellIndex(i,(j-1));
        m_cells->rightHandSide[var][IJ]+= vflux[2*IJ] - vflux[2*IMJ] + vflux[2*IJ+1] - vflux[2*IJM+1];
      }
    }
  }



  ///////////////////////////////////////////
  //////////// SA1/SA2 //////////////////////
  ///////////////////////////////////////////
  for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; ++j){
    for(ZFSId i=m_noGhostLayers-1; i<m_nCells[1]-m_noGhostLayers; ++i) {
      const ZFSId IJ    = cellIndex(i,j);
      const ZFSId IJM   = cellIndex(i,(j-1));

      sa_1flux[ 2 * IJ ] = F1B2 * ( eflux[ 5 * IJ + 3 ] + eflux[ 5 * IJM + 3 ]);
      sa_2flux[ 2 * IJ ] = F1B2 * ( eflux[ 5 * IJ + 4 ] + eflux[ 5 * IJM + 4 ]);
    }
  }

  for(ZFSId j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers; ++j) {
    for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; ++i) {
      const ZFSId IJ    = cellIndex(i,j);
      const ZFSId IMJ   = cellIndex((i-1),j);

      sa_1flux[ 2 * IJ + 1 ] = F1B2 * ( fflux[ 5 * IJ + 3 ] + fflux[ 5 * IMJ + 3 ]);
      sa_2flux[ 2 * IJ + 1 ] = F1B2 * ( fflux[ 5 * IJ + 4 ] + fflux[ 5 * IMJ + 4 ]);
    }
  }


  // separate loop for adding the prodn nad destrn terms for tur kin viscosity transport variable
  for (ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++) {
    for (ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; i++) {
      const ZFSId IJ    = cellIndex(i,j);
      const ZFSId IMJ   = cellIndex((i-1),j);
      const ZFSId IJM   = cellIndex(i,(j-1));

      const ZFSFloat dissipation_term = ((sa_1flux[2*IJ]  -sa_1flux[2*IMJ])  +((sa_2flux[2*IJ]  -sa_2flux[2*IMJ])  *nuTilde[IJ]))*rho[IJ]*m_sa.Fsigma +
        ((sa_1flux[2*IJ+1]-sa_1flux[2*IJM+1])+((sa_2flux[2*IJ+1]-sa_2flux[2*IJM+1])*nuTilde[IJ]))*rho[IJ]*m_sa.Fsigma;

      m_cells->rightHandSide[CV->RANS_VAR[0]][IJ]+=dissipation_term;
      m_cells->rightHandSide[CV->RANS_VAR[0]][IJ]+=m_cells->prodDest[IJ];
    }
  }
}

void ZFSStrctrdBlck2DRans::viscousFlux_FS(){
  cout << "back";
}
void ZFSStrctrdBlck2DRans::viscousFlux_KOmega(){
  cout << "back";
}

void ZFSStrctrdBlck2DRans::computeTurbViscosity()
{
  // OTHER variables required to calculate the laminar viscous fluxes
  ZFSFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  ZFSFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  ZFSFloat* const RESTRICT nuTilde=&m_cells->pvariables[PV->RANS_VAR[0]][0];
  ZFSFloat* const RESTRICT T = &m_cells->temperature[0];

  for (ZFSId i=0; i< m_noStrctrdCells; i++) {
    T[i] = m_block->m_gamma*p[i]/rho[i];
    // decode the kinematic turbulent viscosity from the turb dynamic visc arrays
    const ZFSFloat nuLaminar = zfsSUTHERLANDLAW(T[i])/rho[i];
    const ZFSFloat chi = nuTilde[i]/(nuLaminar);
    const ZFSFloat fv1  =  pow(chi,3)/(pow(chi,3)+m_sa.cv1to3);
    m_cells->fq[FQ->NU_T][i] = fv1*nuTilde[i];
    m_cells->fq[FQ->MU_T][i] = rho[i]*fv1*nuTilde[i];
  }
}


template<ZFSStrctrdBlck2DRans::fluxmethod ausm,ZFSId noVars>
void ZFSStrctrdBlck2DRans::Muscl_(){
  TRACE();

  //stencil identifier
  const ZFSId IJK[2]={1,m_nCells[1]};
  const ZFSFloat *const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const ZFSFloat *const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);

  // const ZFSUint noVars = CV->noVariables;
  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSFloat *const RESTRICT cellVariables= ALIGNED_F(m_cells->pvariables[0]);
  ZFSFloat *const RESTRICT cellRhs= ALIGNED_MF(m_cells->rightHandSide[0]);
  ZFSFloat* RESTRICT flux = ALIGNED_F(m_cells->flux);
  ZFSFloat *const RESTRICT qleft= ALIGNED_MF(m_block->m_QLeft);
  ZFSFloat *const RESTRICT qright= ALIGNED_MF(m_block->m_QRight);

  /////////IMPORTANT PARAMETER
  //ZFSFloat epsi=F1;
  //ZFSFloat kappa=F1B3;
  /////////END IMPORTANT PARAMETER

  for(ZFSId dim=0; dim<nDim; dim++){
    for(ZFSId j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers; j++) {
      for(ZFSId i=m_noGhostLayers-1; i<m_nCells[1]-m_noGhostLayers; i++) {
        const ZFSId I=cellIndex(i,j);
        const ZFSId IP1=I+IJK[dim];
        const ZFSId IM1=I-IJK[dim];
        const ZFSId IP2=I+2*IJK[dim];

        //distances q_i+1 - q_i
        const ZFSFloat DS=sqrt(POW2(x[IP1]-x[I]) + POW2(y[IP1]-y[I]));
        //distances q_i - q_i-1
        const ZFSFloat DSM1=sqrt(POW2(x[I]-x[IM1]) + POW2(y[I]-y[IM1]));
        const ZFSFloat DSP1=sqrt(POW2(x[IP2]-x[IP1]) + POW2(y[IP2]-y[IP1]));
        const ZFSFloat DSP=DS/POW2(DSP1+DS);
        const ZFSFloat DSM=DS/POW2(DSM1+DS);

        for(ZFSUint v=0; v<noVars; ++v){
          const ZFSUint offset = v*noCells;
          const ZFSFloat *const RESTRICT vars = ALIGNED_F(cellVariables+offset);
          const ZFSFloat DQ   = vars[IP1]-vars[I];
          const ZFSFloat DQP1 = vars[IP2]-vars[IP1];
          const ZFSFloat DQM1 = vars[I]-vars[IM1];
          qleft[v] = vars[I]+DSM*(DSM1*DQ+DS*DQM1);
          qright[v]= vars[IP1]-DSP*(DS*DQP1+DSP1*DQ);
        }

        (this->*ausm)(qleft, qright, dim, I); //Flux balance in AUSM
      }
    }

    //FLUX BALANCE
    for(ZFSUint v=0; v<noVars; v++) {
      for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++) {
        for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; i++) {
          const ZFSId I=cellIndex(i,j);
          const ZFSId IM1=I-IJK[dim];
          const ZFSUint offset = v*noCells;
          ZFSFloat *const RESTRICT rhs = ALIGNED_F(cellRhs+offset);
          rhs[I]+= flux[IM1+noCells*v]-flux[I+noCells*v];
        }
      }
    }
  }
}

template void ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm,5>();
template void ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm,6>();
template void ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::Ausm,7>();
template void ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::AusmDV,5>();
template void ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::AusmDV,6>();
template void ZFSStrctrdBlck2DRans::Muscl_<&ZFSStrctrdBlck2DRans::AusmDV,7>();


inline ZFSId ZFSStrctrdBlck2DRans::cellIndex(ZFSInt i, ZFSInt j)
{
  return i+j*m_nCells[1];
}

inline ZFSId ZFSStrctrdBlck2DRans::getCellIdfromCell( ZFSId origin, ZFSInt incI, ZFSInt incJ)
{
  return origin + incI + incJ * m_nCells[1];
}

inline ZFSFloat ZFSStrctrdBlck2DRans::getPSI(ZFSId I, ZFSId dim) {
  const ZFSFloat FK = 18.0;
  const ZFSId IJK[2]={1,m_nCells[1]};
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
