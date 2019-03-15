#include "zfsstrctrdblck3drans.h"
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


ZFSStrctrdBlck3DRans::ZFSStrctrdBlck3DRans(class ZFSStrctrdBlck3D* block ):
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

  initFluxMethod();
}

ZFSStrctrdBlck3DRans::~ZFSStrctrdBlck3DRans(){
}

void ZFSStrctrdBlck3DRans::initFluxMethod() {
  //m_strctrdBndryCndRans = new ZFSStrctrdBndryCnd3DRans(m_block, m_block->m_noSpecies);
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
    case RANS_SA_DV :
      {
        zfsAlloc(m_cells->saFlux1, nDim*m_noStrctrdCells, "m_cells->saFlux1", -999999.9, __CALLING_FUNCTION__);
        zfsAlloc(m_cells->saFlux2, nDim*m_noStrctrdCells, "m_cells->saFlux2", -999999.9, __CALLING_FUNCTION__);
        zfsAlloc(m_cells->prodDest, m_noStrctrdCells, "m_cells->prodDest", -999999.9, __CALLING_FUNCTION__);
        viscFluxMethod = &ZFSStrctrdBlck3DRans::viscousFlux_SA;
        switch(m_block->CV->noVariables){
        case 5: { reconstructSurfaceData=&ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV, 5>; break;}
        case 6: { reconstructSurfaceData=&ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV, 6>; break;}
        case 7: { reconstructSurfaceData=&ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV, 7>; break;}
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
        viscFluxMethod = &ZFSStrctrdBlck3DRans::viscousFlux_FS;
        switch(m_block->CV->noVariables){
        case 5: { reconstructSurfaceData=&ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV, 5>; break;}
        case 6: { reconstructSurfaceData=&ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV, 6>; break;}
        case 7: { reconstructSurfaceData=&ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV, 7>; break;}
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
        viscFluxMethod = &ZFSStrctrdBlck3DRans::viscousFlux_KOmega;
        switch(m_block->CV->noVariables){
        case 5: { reconstructSurfaceData=&ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV, 5>; break;}
        case 6: { reconstructSurfaceData=&ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV, 6>; break;}
        case 7: { reconstructSurfaceData=&ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV, 7>; break;}
        default:{
          stringstream errorMessage;
          errorMessage << "Number of Variables " << m_block->CV->noVariables << " not implemented! in template Rans AUSM " << endl;
          zfsTerm(1, __CALLING_FUNCTION__, errorMessage.str());
        }
        }
        break;
      }
    default:
      { zfsTerm(1, __CALLING_FUNCTION__, "RANS METHOD wsa not specified in properties"); break; }
    }
}

void ZFSStrctrdBlck3DRans::Muscl(){
  (this->*reconstructSurfaceData)();
}

void ZFSStrctrdBlck3DRans::AusmDV(ZFSFloat* QLeft, ZFSFloat*  QRight, const ZFSId dim, const ZFSId I){
  ZFSFloat pFactor[3]={F0,F0,F0};
  static const ZFSFloat gamma = m_block->m_gamma;
  static const ZFSFloat gammaMinusOne = gamma - 1.0;
  static const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_cells->surfaceMetrics[I]);

  //left side
  const ZFSFloat RHOL = QLeft[PV->RHO];
  const ZFSFloat FRHOL = F1/RHOL;
  ZFSFloat UL = QLeft[PV->U];
  ZFSFloat VL = QLeft[PV->V];
  ZFSFloat WL = QLeft[PV->W];
  const ZFSFloat PL = QLeft[PV->P];
  const ZFSFloat NUTILDEL = QLeft[PV->RANS_VAR[0]];

  //right side
  const ZFSFloat RHOR = QRight[PV->RHO];
  const ZFSFloat FRHOR = F1/RHOR;
  ZFSFloat UR = QRight[PV->U];
  ZFSFloat VR = QRight[PV->V];
  ZFSFloat WR = QRight[PV->W];
  const ZFSFloat PR = QRight[PV->P];
  const ZFSFloat NUTILDER = QRight[PV->RANS_VAR[0]];

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
                         WL * surf[ dim*3+2 ])) * FDGRAD;


  const ZFSFloat UUR = ((UR * surf[ dim*3+0 ] +
                         VR * surf[ dim*3+1 ] +
                         WR * surf[ dim*3+2 ])) * FDGRAD;

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

  static const ZFSId IJK[3]={1,m_nCells[2],m_nCells[1]*m_nCells[2]};
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
  ZFSFloat *const RESTRICT flux = ALIGNED_MF(m_cells->flux);

  flux[I+noCells*CV->RHO_U] = RHOU2 * ( UL + UR) + ARHOU2 * ( UL - UR) + PLR*pFactor[0] + RHOUL*UUL2+RHOUR*UUR2;
  flux[I+noCells*CV->RHO_V] = RHOU2 * ( VL + VR) + ARHOU2 * ( VL - VR) + PLR*pFactor[1] + RHOUL*UUL3+RHOUR*UUR3;
  flux[I+noCells*CV->RHO_W] = RHOU2 * ( WL + WR) + ARHOU2 * ( WL - WR) + PLR*pFactor[2] + RHOUL*UUL4+RHOUR*UUR4;
  flux[I+noCells*CV->RHO_E] = RHOU2 * ( e0 + e1) + ARHOU2 * ( e0 - e1);
  flux[I+noCells*CV->RHO]   = RHOU;
  flux[I+noCells*CV->RANS_VAR[0]] = RHOU2*(NUTILDEL+NUTILDER) + ARHOU2*(NUTILDEL-NUTILDER);
}


void ZFSStrctrdBlck3DRans::viscousFluxRANS(){
  (this->*viscFluxMethod)();
}

void ZFSStrctrdBlck3DRans::viscousFlux_SA(){
  computeTurbViscosity();

  //OTHER variables required to calculate the laminar viscous fluxes
  const ZFSInt noCells = m_noStrctrdCells;
  const ZFSFloat rPrL = F1/m_block->m_Pr;
  const ZFSFloat rPrT = F1/0.9;
  const ZFSFloat rRe = F1 / m_block->m_Re0;
  const ZFSFloat gammaMinusOne = m_block->m_gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  ZFSFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  ZFSFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  ZFSFloat* const RESTRICT w = &m_cells->pvariables[PV->W][0];
  ZFSFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  ZFSFloat* const RESTRICT nuTilde=&m_cells->pvariables[PV->RANS_VAR[0]][0];
  ZFSFloat* const RESTRICT T = &m_cells->temperature[0];

  ZFSFloat* const RESTRICT mut = &m_cells->fq[FQ->MU_T][0];
  ZFSFloat* const RESTRICT eflux=&m_cells->eFlux[0];
  ZFSFloat* const RESTRICT fflux=&m_cells->fFlux[0];
  ZFSFloat* const RESTRICT gflux=&m_cells->gFlux[0];
  ZFSFloat* const RESTRICT vflux=&m_cells->viscousFlux[0];
  ZFSFloat* const RESTRICT sa_1flux=&m_cells->saFlux1[0];
  ZFSFloat* const RESTRICT sa_2flux=&m_cells->saFlux2[0];


  //junoh
  const ZFSFloat lamXCoordinates =-9.02;
  for (ZFSId cellId=0; cellId< m_noStrctrdCells; cellId++) {    
    if(m_cells->coordinates[0][cellId]<lamXCoordinates){
      nuTilde[cellId]=0.0;    
    }
  }


  for (ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers+1; k++) {
    for (ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers+1; j++) {
      for (ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers+1; i++) {
        //get the adjacent cells;
        const ZFSId IJK   = cellIndex(i,j,k);
        const ZFSId IPJK  = cellIndex((i+1),j, k);
        const ZFSId IPJPK = cellIndex((i+1),(j+1), k);
        const ZFSId IJPK  = cellIndex(i,(j+1), k);
        const ZFSId IJKP  = cellIndex(i,j,(k+1));
        const ZFSId IPJKP = cellIndex((i+1),j,(k+1));
        const ZFSId IPJPKP= cellIndex((i+1),(j+1),(k+1));
        const ZFSId IJPKP = cellIndex(i,(j+1),(k+1));

        const ZFSId IMJK  = cellIndex((i-1),j, k);
        const ZFSId IJMK  = cellIndex(i,(j-1), k);
        const ZFSId IJKM  = cellIndex(i,j,(k-1));

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


        const ZFSFloat dnutldxi=F1B4*(nuTilde[IPJPKP]+nuTilde[IPJPK]+nuTilde[IPJKP]+nuTilde[IPJK]-nuTilde[IJPKP]-nuTilde[IJPK]-nuTilde[IJKP]-nuTilde[IJK]);
        const ZFSFloat dnutldet=F1B4*(nuTilde[IPJPKP]+nuTilde[IJPKP]+nuTilde[IPJPK]+nuTilde[IJPK]-nuTilde[IPJKP]-nuTilde[IJKP]-nuTilde[IPJK]-nuTilde[IJK]);
        const ZFSFloat dnutldze=F1B4*(nuTilde[IPJPKP]+nuTilde[IJPKP]+nuTilde[IPJKP]+nuTilde[IJKP]-nuTilde[IPJPK]-nuTilde[IJPK]-nuTilde[IPJK]-nuTilde[IJK]);

        const ZFSFloat uAvg=F1B8*(u[IPJPKP]+u[IJPKP]+u[IJPK]+u[IPJPK]+u[IPJKP]+u[IJKP]+u[IJK]+u[IPJK]);
        const ZFSFloat vAvg=F1B8*(v[IPJPKP]+v[IJPKP]+v[IJPK]+v[IPJPK]+v[IPJKP]+v[IJKP]+v[IJK]+v[IPJK]);
        const ZFSFloat wAvg=F1B8*(w[IPJPKP]+w[IJPKP]+w[IJPK]+w[IPJPK]+w[IPJKP]+w[IJKP]+w[IJK]+w[IPJK]);

        const ZFSFloat nutldAvg=F1B8*(nuTilde[IPJPKP]+nuTilde[IJPKP]+nuTilde[IJPK]+nuTilde[IPJPK]+nuTilde[IPJKP]+nuTilde[IJKP]+nuTilde[IJK]+nuTilde[IPJK]);

        const ZFSFloat nuAvg=F1B8*((zfsSUTHERLANDLAW(abs(T[IPJPKP]))/rho[IPJPKP])
                                   +(zfsSUTHERLANDLAW(abs(T[IJPKP]))/rho[IJPKP])
                                   +(zfsSUTHERLANDLAW(abs(T[IJPK]))/rho[IJPK])
                                   +(zfsSUTHERLANDLAW(abs(T[IPJPK]))/rho[IPJPK])
                                   +(zfsSUTHERLANDLAW(abs(T[IPJKP]))/rho[IPJKP])
                                   +(zfsSUTHERLANDLAW(abs(T[IJKP]))/rho[IJKP])
                                   +(zfsSUTHERLANDLAW(abs(T[IJK]))/rho[IJK])
                                   +(zfsSUTHERLANDLAW(abs(T[IPJK]))/rho[IPJK]));

        const ZFSFloat muAvg=F1B8*(zfsSUTHERLANDLAW(T[IPJPKP])+zfsSUTHERLANDLAW(T[IJPKP])+zfsSUTHERLANDLAW(T[IJPK])+zfsSUTHERLANDLAW(T[IPJPK])+zfsSUTHERLANDLAW(T[IPJKP])+zfsSUTHERLANDLAW(T[IJKP])+zfsSUTHERLANDLAW(T[IJK])+zfsSUTHERLANDLAW(T[IPJK]));
        const ZFSFloat mutAvg = F1B8*(mut[IPJPKP]+mut[IJPKP]+mut[IJPK]+mut[IPJPK]+mut[IPJKP]+mut[IJKP]+mut[IJK]+mut[IPJK]);


        const ZFSFloat dnutldx = dnutldxi * cornerMetrics[ xsd * 3 + xsd ] +
          dnutldet * cornerMetrics[ ysd * 3 + xsd ] +
          dnutldze * cornerMetrics[ zsd * 3 + xsd ];

        const ZFSFloat dnutldy = dnutldxi * cornerMetrics[ xsd * 3 + ysd ] +
          dnutldet * cornerMetrics[ ysd * 3 + ysd ] +
          dnutldze * cornerMetrics[ zsd * 3 + ysd ];

        const ZFSFloat dnutldz = dnutldxi * cornerMetrics[ xsd * 3 + zsd ] +
          dnutldet * cornerMetrics[ ysd * 3 + zsd ] +
          dnutldze * cornerMetrics[ zsd * 3 + zsd ];

        const ZFSFloat Frj = rRe/m_cells->cornerJac[IJK];


        const ZFSFloat sax1 = Frj*(nuAvg+(1.0+m_sa.cb2)*nutldAvg)*
          (dnutldx * cornerMetrics[ xsd * 3 + xsd ]+
           dnutldy * cornerMetrics[ xsd * 3 + ysd ]+
           dnutldz * cornerMetrics[ xsd * 3 + zsd ]);

        const ZFSFloat say1 = Frj*(nuAvg+(1.0+m_sa.cb2)*nutldAvg)*
          (dnutldx * cornerMetrics[ ysd * 3 + xsd ]+
           dnutldy * cornerMetrics[ ysd * 3 + ysd ]+
           dnutldz * cornerMetrics[ ysd * 3 + zsd ]);

        const ZFSFloat saz1 = Frj*(nuAvg+(1.0+m_sa.cb2)*nutldAvg)*
          (dnutldx * cornerMetrics[ zsd * 3 + xsd ]+
           dnutldy * cornerMetrics[ zsd * 3 + ysd ]+
           dnutldz * cornerMetrics[ zsd * 3 + zsd ]);

        const ZFSFloat sax2 = -Frj*m_sa.cb2*
          (dnutldx * cornerMetrics[ xsd * 3 + xsd ]+
           dnutldy * cornerMetrics[ xsd * 3 + ysd ]+
           dnutldz * cornerMetrics[ xsd * 3 + zsd ]);
        const ZFSFloat say2 = -Frj*m_sa.cb2*
          (dnutldx * cornerMetrics[ ysd * 3 + xsd ]+
           dnutldy * cornerMetrics[ ysd * 3 + ysd ]+
           dnutldz * cornerMetrics[ ysd * 3 + zsd ]);

        const ZFSFloat saz2=  -Frj*m_sa.cb2*
          (dnutldx * cornerMetrics[ zsd * 3 + xsd ]+
           dnutldy * cornerMetrics[ zsd * 3 + ysd ]+
           dnutldz * cornerMetrics[ zsd * 3 + zsd ]);


        //for dwdy
        const ZFSFloat dw1 = w[IPJK+2]- w[IMJK+2];
        const ZFSFloat dw2 = w[IJPK+2]- w[IJMK+2];
        const ZFSFloat dw3 = w[IJKP+2]- w[IJKM+2];

        //for dvdz
        const ZFSFloat dv1 = v[IPJK+1]- v[IMJK+1];
        const ZFSFloat dv2 = v[IJPK+1]- v[IJMK+1];
        const ZFSFloat dv3 = v[IJKP+1]- v[IJKM+1];

        //for dudz
        const ZFSFloat du1 = u[IPJK]- u[IMJK];
        const ZFSFloat du2 = u[IJPK]- u[IJMK];
        const ZFSFloat du3 = u[IJKP]- u[IJKM];

        const ZFSFloat vorti = (m_cells->cellMetrics[IJK][ xsd * 3 + ysd ]*dw1)+
          (m_cells->cellMetrics[IJK][ ysd * 3 + ysd ]*dw2)+
          (m_cells->cellMetrics[IJK][ zsd * 3 + ysd ]*dw3)-

          (m_cells->cellMetrics[IJK][ xsd * 3 + zsd ]*dv1)-
          (m_cells->cellMetrics[IJK][ ysd * 3 + zsd ]*dv2)-
          (m_cells->cellMetrics[IJK][ zsd * 3 + zsd ]*dv3);

        const ZFSFloat vortj = (m_cells->cellMetrics[IJK][ xsd * 3 + zsd ]*du1)+
          (m_cells->cellMetrics[IJK][ ysd * 3 + zsd ]*du2)+
          (m_cells->cellMetrics[IJK][ zsd * 3 + zsd ]*du3)-

          (m_cells->cellMetrics[IJK][ xsd * 3 + xsd ]*dw1)-
          (m_cells->cellMetrics[IJK][ ysd * 3 + xsd ]*dw2)-
          (m_cells->cellMetrics[IJK][ zsd * 3 + xsd ]*dw3);

        const ZFSFloat vortk = (m_cells->cellMetrics[IJK][ xsd * 3 + xsd ]*dv1)+
          (m_cells->cellMetrics[IJK][ ysd * 3 + xsd ]*dv2)+
          (m_cells->cellMetrics[IJK][ zsd * 3 + xsd ]*dv3)-

          (m_cells->cellMetrics[IJK][ xsd * 3 + ysd ]*du1)-
          (m_cells->cellMetrics[IJK][ ysd * 3 + ysd ]*du2)-
          (m_cells->cellMetrics[IJK][ zsd * 3 + ysd ]*du3);

        ZFSFloat s = (vorti*vorti)+(vortj*vortj)+(vortk*vortk);
        s = F1B2*sqrt(s)/m_cells->cellJac[IJK];

        // assuming wall distance function
        const ZFSFloat distance = m_cells->fq[FQ->WALLDISTANCE][IJK];
        const ZFSFloat Fdist2 = 1.0/(distance*distance);
        const ZFSFloat chi = nuTilde[IJK]*rho[IJK]/(zfsSUTHERLANDLAW(T[IJK])/rho[IJK]);
        const ZFSFloat chip3 = chi*chi*chi;
        const ZFSFloat Fv1 =  chip3/(chip3+m_sa.cv1to3);
        const ZFSFloat Fv2 = F1 - (chi/(F1+chi*Fv1));

        const ZFSFloat term = nuTilde[IJK]*Fdist2*m_sa.Fkap2;
        const ZFSFloat stilde = s+term*Fv2*rRe;
        const ZFSFloat r = min(10.0, rRe*term/stilde);

        const ZFSFloat g = r+m_sa.cw2*(pow(r,6)-r);
        const ZFSFloat Fwterm = (1+m_sa.cw3to6)/(pow(g,6)+m_sa.cw3to6);
        const ZFSFloat Fw = g*pow(Fwterm,(1.0/6.0));
        const ZFSFloat prodValue = rho[IJK]*m_sa.cb1*(F1-m_sa.Ft2)*stilde*nuTilde[IJK];
        const ZFSFloat destValue = rRe* rho[IJK]*(m_sa.cw1*Fw-m_sa.cb1*m_sa.Fkap2*m_sa.Ft2)*pow( nuTilde[IJK],2.0)*Fdist2;

        m_cells->prodDest[IJK] = (prodValue-destValue)*m_cells->cellJac[IJK];


        //mue=zfsSUTHERLANDLAW(TAvg);

        // compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

        // tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
        //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
        //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
        ZFSFloat tau1 = F4B3 * ( dudxi * cornerMetrics[ xsd * 3 + xsd ] +
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
        ZFSFloat tau2 = dudxi * cornerMetrics[ xsd * 3 + ysd ] +
          dudet * cornerMetrics[ ysd * 3 + ysd ] +
          dudze * cornerMetrics[ zsd * 3 + ysd ] +

          dvdxi * cornerMetrics[ xsd * 3 + xsd ] +
          dvdet * cornerMetrics[ ysd * 3 + xsd ] +
          dvdze * cornerMetrics[ zsd * 3 + xsd ];

        // compute tau3 = du/dz + dw/dx

        // tau_xz = du/dxi * dxi/dz + du/deta * deta/dz + du/dzeta * dzeta/dz
        //        + dw/dxi * dxi/dx + dw/deta * deta/dx + dw/dzeta * dzeta/dx
        ZFSFloat tau3 = dudxi * cornerMetrics[ xsd * 3 + zsd ] +
          dudet * cornerMetrics[ ysd * 3 + zsd ] +
          dudze * cornerMetrics[ zsd * 3 + zsd ] +

          dwdxi * cornerMetrics[ xsd * 3 + xsd ] +
          dwdet * cornerMetrics[ ysd * 3 + xsd ] +
          dwdze * cornerMetrics[ zsd * 3 + xsd ];

        // compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

        // tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
        //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
        //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
        ZFSFloat tau4 = F4B3 * ( dvdxi * cornerMetrics[ xsd * 3 + ysd ] +
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
        ZFSFloat tau5 = dvdxi * cornerMetrics[ xsd * 3 + zsd ] +
          dvdet * cornerMetrics[ ysd * 3 + zsd ] +
          dvdze * cornerMetrics[ zsd * 3 + zsd ] +

          dwdxi * cornerMetrics[ xsd * 3 + ysd ] +
          dwdet * cornerMetrics[ ysd * 3 + ysd ] +
          dwdze * cornerMetrics[ zsd * 3 + ysd ];

        // compute tau6 = 2 dw/dz - 2/3 ( du/dx + dv/dy + dw/dz )

        // tau_zz = 4/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
        //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
        //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
        ZFSFloat tau6 = F4B3 * ( dwdxi * cornerMetrics[ xsd * 3 + zsd ] +
                        dwdet * cornerMetrics[ ysd * 3 + zsd ] +
                        dwdze * cornerMetrics[ zsd * 3 + zsd ] ) -

          F2B3 * ( dudxi * cornerMetrics[ xsd * 3 + xsd ] +
                   dudet * cornerMetrics[ ysd * 3 + xsd ] +
                   dudze * cornerMetrics[ zsd * 3 + xsd ] ) -

          F2B3 * ( dvdxi * cornerMetrics[ xsd * 3 + ysd ] +
                   dvdet * cornerMetrics[ ysd * 3 + ysd ] +
                   dvdze * cornerMetrics[ zsd * 3 + ysd ] );

        ZFSFloat dTdx = dTdxi * cornerMetrics[ xsd * 3 + xsd ] +
          dTdet * cornerMetrics[ ysd * 3 + xsd ] +
          dTdze * cornerMetrics[ zsd * 3 + xsd ];

        ZFSFloat dTdy = dTdxi * cornerMetrics[ xsd * 3 + ysd ] +
          dTdet * cornerMetrics[ ysd * 3 + ysd ] +
          dTdze * cornerMetrics[ zsd * 3 + ysd ];

        ZFSFloat dTdz = dTdxi * cornerMetrics[ xsd * 3 + zsd ] +
          dTdet * cornerMetrics[ ysd * 3 + zsd ] +
          dTdze * cornerMetrics[ zsd * 3 + zsd ];

        const ZFSFloat muOverRe = (muAvg+mutAvg) * rRe / m_cells->cornerJac[IJK]; // divide by Jacobian
        tau1*=muOverRe;
        tau2*=muOverRe;
        tau3*=muOverRe;
        tau4*=muOverRe;
        tau5*=muOverRe;
        tau6*=muOverRe;

        const ZFSFloat muH = rRe*FgammaMinusOne/ m_cells->cornerJac[IJK]*((muAvg*rPrL)+(mutAvg*rPrT));

        const ZFSFloat qx=muH*dTdx+uAvg*tau1+vAvg*tau2+wAvg*tau3;
        const ZFSFloat qy=muH*dTdy+uAvg*tau2+vAvg*tau4+wAvg*tau5;
        const ZFSFloat qz=muH*dTdz+uAvg*tau3+vAvg*tau5+wAvg*tau6;

        //efluxes
        eflux[ 6* IJK     ]    = tau1 * cornerMetrics[ xsd * 3 + xsd ] +
          tau2 * cornerMetrics[ xsd * 3 + ysd ] +
          tau3 * cornerMetrics[ xsd * 3 + zsd ];

        eflux[ 6* IJK + 1 ] = tau2 * cornerMetrics[ xsd * 3 + xsd ] +
          tau4 * cornerMetrics[ xsd * 3 + ysd ] +
          tau5 * cornerMetrics[ xsd * 3 + zsd ];

        eflux[ 6* IJK + 2 ] = tau3 * cornerMetrics[ xsd * 3 + xsd ] +
          tau5 * cornerMetrics[ xsd * 3 + ysd ] +
          tau6 * cornerMetrics[ xsd * 3 + zsd ];

        eflux[ 6* IJK + 3 ] = qx * cornerMetrics[ xsd * 3 + xsd ] +
          qy * cornerMetrics[ xsd * 3 + ysd ] +
          qz * cornerMetrics[ xsd * 3 + zsd ];
        eflux[ 6* IJK + 4 ]=  sax1;// diffusion of nutilde for every cell
        eflux[ 6* IJK + 5 ]=  sax2;// diffusion of nutilde for every cell


        //ffluxes
        fflux[ 6 * IJK ]    = tau1 * cornerMetrics[ ysd * 3 + xsd ] +
          tau2 * cornerMetrics[ ysd * 3 + ysd ] +
          tau3 * cornerMetrics[ ysd * 3 + zsd ];

        fflux[ 6* IJK + 1 ] = tau2 * cornerMetrics[ ysd * 3 + xsd ] +
          tau4 * cornerMetrics[ ysd * 3 + ysd ] +
          tau5 * cornerMetrics[ ysd * 3 + zsd ];

        fflux[ 6* IJK + 2 ] = tau3 * cornerMetrics[ ysd * 3 + xsd ] +
          tau5 * cornerMetrics[ ysd * 3 + ysd ] +
          tau6 * cornerMetrics[ ysd * 3 + zsd ];

        fflux[ 6* IJK + 3 ] = qx * cornerMetrics[ ysd * 3 + xsd ] +
          qy * cornerMetrics[ ysd * 3 + ysd ] +
          qz * cornerMetrics[ ysd * 3 + zsd ];
        fflux[ 6* IJK + 4 ]=say1;// diffusion of nutilde for every cell
        fflux[ 6* IJK + 5 ]=say2;// diffusion of nutilde for every cell

        //gfluxes
        gflux[ 6 * IJK ]    = tau1 * cornerMetrics[ zsd * 3 + xsd ] +
          tau2 * cornerMetrics[ zsd * 3 + ysd ] +
          tau3 * cornerMetrics[ zsd * 3 + zsd ];

        gflux[ 6* IJK + 1 ] = tau2 * cornerMetrics[ zsd * 3 + xsd ] +
          tau4 * cornerMetrics[ zsd * 3 + ysd ] +
          tau5 * cornerMetrics[ zsd * 3 + zsd ];

        gflux[ 6* IJK + 2 ] = tau3 * cornerMetrics[ zsd * 3 + xsd ] +
          tau5 * cornerMetrics[ zsd * 3 + ysd ] +
          tau6 * cornerMetrics[ zsd * 3 + zsd ];

        gflux[ 6* IJK + 3 ] = qx * cornerMetrics[ zsd * 3 + xsd ] +
          qy * cornerMetrics[ zsd * 3 + ysd ] +
          qz * cornerMetrics[ zsd * 3 + zsd ];
        gflux[ 6* IJK + 4 ]=saz1;// diffusion of nutilde for every cell
        gflux[ 6* IJK + 5 ]=saz2;// diffusion of nutilde for every cell
      }
    }
  }

  for (ZFSId var=0; var<CV->noVariables-2; var++) {
    for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
      for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
        for (ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers; i++) {
          const ZFSId IJK    = cellIndex(i,j,k);
          const ZFSId IJMK   = cellIndex(i,(j-1),k);
          const ZFSId IJKM   = cellIndex(i,j,(k-1));
          const ZFSId IJMKM  = cellIndex(i,(j-1),(k-1));

          vflux[3*IJK]=F1B4*(eflux[6*IJK+var]+eflux[6*IJKM+var]+eflux[6*IJMK+var]+eflux[6*IJMKM+var]);//*m_cells->area[0][IJK];
        }
      }
    }


    for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
      for (ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers; j++) {
        for (ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
          const ZFSId IJK   = cellIndex(i,j,k);
          const ZFSId IMJK  = cellIndex((i-1),j,k);
          const ZFSId IJKM  = cellIndex(i,j,(k-1));
          const ZFSId IMJKM = cellIndex((i-1),j,(k-1));

          vflux[3*IJK+1]=F1B4*(fflux[6*IJK+var]+fflux[6*IJKM+var]+fflux[6*IMJK+var]+fflux[6*IMJKM+var]);//*m_cells->area[1][IJK];
        }
      }
    }

    for (ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers; k++) {
      for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
        for (ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
          const ZFSId IJK   = cellIndex(i,j,k);
          const ZFSId IMJK  = cellIndex((i-1),j,k);
          const ZFSId IJMK  = cellIndex(i,(j-1),k);
          const ZFSId IMJMK = cellIndex((i-1),(j-1),k);

          vflux[3*IJK+2]=F1B4*(gflux[6*IJK+var]+gflux[6*IMJK+var]+gflux[6*IJMK+var]+gflux[6*IMJMK+var]);
        }
      }
    }

    for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
      for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
        for (ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++){
          const ZFSId IJK=cellIndex(i,j,k);
          const ZFSId IMJK=cellIndex(i-1,j,k);
          const ZFSId IJMK=cellIndex(i,j-1,k);
          const ZFSId IJKM=cellIndex(i,j,k-1);
          m_cells->rightHandSide[var][IJK]+= vflux[3*IJK]-vflux[3*IMJK]+vflux[3*IJK+1]-vflux[3*IJMK+1]+vflux[3*IJK+2]-vflux[3*IJKM+2];
        }
      }
    }
  }

  for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
    for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
      for (ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers; i++) {
        const ZFSId IJK    = cellIndex(i,j,k);
        const ZFSId IJMK   = cellIndex(i,(j-1),k);
        const ZFSId IJKM   = cellIndex(i,j,(k-1));
        const ZFSId IJMKM  = cellIndex(i,(j-1),(k-1));

        sa_1flux[3*IJK]=F1B4*(eflux[6*IJK+4]+eflux[6*IJKM+4]+eflux[6*IJMK+4]+eflux[6*IJMKM+4]);//*m_cells->area[0][IJK];
        sa_2flux[3*IJK]=F1B4*(eflux[6*IJK+5]+eflux[6*IJKM+5]+eflux[6*IJMK+5]+eflux[6*IJMKM+5]);//*m_cells->area[0][IJK];
      }
    }
  }


  for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
    for (ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers; j++) {
      for (ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
        const ZFSId IJK   = cellIndex(i,j,k);
        const ZFSId IMJK  = cellIndex((i-1),j,k);
        const ZFSId IJKM  = cellIndex(i,j,(k-1));
        const ZFSId IMJKM = cellIndex((i-1),j,(k-1));

        sa_1flux[3*IJK+1]=F1B4*(fflux[6*IJK+4]+fflux[6*IJKM+4]+fflux[6*IMJK+4]+fflux[6*IMJKM+4]);
        sa_2flux[3*IJK+1]=F1B4*(fflux[6*IJK+5]+fflux[6*IJKM+5]+fflux[6*IMJK+5]+fflux[6*IMJKM+5]);
      }
    }
  }

  for (ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers; k++) {
    for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
      for (ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
        const ZFSId IJK   = cellIndex(i,j,k);
        const ZFSId IMJK  = cellIndex((i-1),j,k);
        const ZFSId IJMK  = cellIndex(i,(j-1),k);
        const ZFSId IMJMK = cellIndex((i-1),(j-1),k);

        sa_1flux[3*IJK+2]=F1B4*(gflux[6*IJK+4]+gflux[6*IMJK+4]+gflux[6*IJMK+4]+gflux[6*IMJMK+4]);
        sa_2flux[3*IJK+2]=F1B4*(gflux[6*IJK+5]+gflux[6*IMJK+5]+gflux[6*IJMK+5]+gflux[6*IMJMK+5]);
      }
    }
  }

  // separate loop for adding the prodn nad destrn terms for tur kin viscosity transport variable
  for (ZFSId k=m_noGhostLayers; k<m_nCells[0]-m_noGhostLayers; k++) {
    for (ZFSId j=m_noGhostLayers; j<m_nCells[1]-m_noGhostLayers; j++) {
      for (ZFSId i=m_noGhostLayers; i<m_nCells[2]-m_noGhostLayers; i++) {
        const ZFSId IJK=cellIndex(i,j,k);
        const ZFSId IMJK=cellIndex(i-1,j,k);
        const ZFSId IJMK=cellIndex(i,j-1,k);
        const ZFSId IJKM=cellIndex(i,j,k-1);
        const ZFSFloat dissipation_term = ((sa_1flux[3*IJK]  -sa_1flux[3*IMJK])  +((sa_2flux[3*IJK]  -sa_2flux[3*IMJK])*nuTilde[IJK]))*rho[IJK]*m_sa.Fsigma +
                                          ((sa_1flux[3*IJK+1]-sa_1flux[3*IJMK+1])+((sa_2flux[3*IJK+1]-sa_2flux[3*IJMK+1])*nuTilde[IJK]))*rho[IJK]*m_sa.Fsigma +
                                          ((sa_1flux[3*IJK+2]-sa_1flux[3*IJKM+2])+((sa_2flux[3*IJK+2]-sa_2flux[3*IJKM+2])*nuTilde[IJK]))*rho[IJK]*m_sa.Fsigma ;

        m_cells->rightHandSide[CV->RANS_FIRST][IJK]+=dissipation_term;
        m_cells->rightHandSide[CV->RANS_FIRST][IJK]+=m_cells->prodDest[IJK];
      }
    }
  }
}

void ZFSStrctrdBlck3DRans::viscousFlux_FS(){
  cout << "back";
}
void ZFSStrctrdBlck3DRans::viscousFlux_KOmega(){
  cout << "back";
}

void ZFSStrctrdBlck3DRans::computeTurbViscosity()
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

template<ZFSStrctrdBlck3DRans::fluxmethod ausm,ZFSId noVars>
void ZFSStrctrdBlck3DRans::Muscl_(){
  TRACE();

  //stencil identifier
  const ZFSId IJK[3]={1,m_nCells[2], m_nCells[1]*m_nCells[2]};

  ZFSFloat* RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  ZFSFloat* RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  ZFSFloat* RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);

  // const ZFSUint noVars = CV->noVariables;
  const ZFSUint noCells = m_noStrctrdCells;
  ZFSFloat* RESTRICT cellVariables= ALIGNED_F(m_cells->pvariables[0]);
  ZFSFloat* RESTRICT cellRhs= ALIGNED_MF(m_cells->rightHandSide[0]);
  ZFSFloat* RESTRICT qleft= ALIGNED_MF(m_block->m_QLeft);
  ZFSFloat* RESTRICT qright= ALIGNED_MF(m_block->m_QRight);
  ZFSFloat* RESTRICT flux = ALIGNED_F(m_cells->flux);

  /////////IMPORTANT PARAMETER
  //ZFSFloat epsi=F1;
  //ZFSFloat kappa=F1B3;
  /////////END IMPORTANT PARAMETER

  for(ZFSId dim=0; dim<nDim; dim++){
    for(ZFSId k=m_noGhostLayers-1; k<m_nCells[0]-m_noGhostLayers; k++) {
      for(ZFSId j=m_noGhostLayers-1; j<m_nCells[1]-m_noGhostLayers; j++) {
        for(ZFSId i=m_noGhostLayers-1; i<m_nCells[2]-m_noGhostLayers; i++) {
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

#ifdef ZFS_EXTRA_DEBUG
            ZFSFloat Flux = (m_cells->flux[IM1+noCells*v]-m_cells->flux[I+v*noCells]);
            m_block-> convFluxOut[dim][v*noCells+I]=Flux;
#endif
          }
        }
      }
    }
  }
}

template void ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV,5>();
template void ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV,6>();
template void ZFSStrctrdBlck3DRans::Muscl_<&ZFSStrctrdBlck3DRans::AusmDV,7>();


inline ZFSId ZFSStrctrdBlck3DRans::cellIndex(ZFSInt i, ZFSInt j, ZFSInt k)
{
  return i+(j+k*m_nCells[1])*m_nCells[2];
}

inline ZFSId ZFSStrctrdBlck3DRans::getCellIdfromCell( ZFSId origin, ZFSInt incI, ZFSInt incJ, ZFSInt incK )
{
  return origin + incI + incJ * m_nCells[2] + incK * m_nCells[2] * m_nCells[1];
}

inline ZFSFloat ZFSStrctrdBlck3DRans::getPSI(ZFSId I, ZFSId dim) {
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
