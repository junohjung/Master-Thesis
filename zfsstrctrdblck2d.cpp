#include "zfsstrctrdblck2d.h"
#include "zfsglobals.h"
#include "zfsconstants.h"
#include "zfsiolib.h"


/**
 *Constructor
 */

ZFSStrctrdBlck2D::ZFSStrctrdBlck2D( ZFSId blockId, bool* propertiesGroups, const MPI_Comm comm): ZFSStrctrdBlck<2>( blockId, propertiesGroups, comm)
{
  TRACE();
  const ZFSLong oldAllocatedBytes = allocatedBytes();

  //count the no of necessary FQ fields and allocate
  initializeFQField();

  // compute the cell center coordinates from point coordinates
  computeCellCentreCoordinates();

  if(m_rans) {
    m_strctrdBndryCnd = new ZFSStrctrdBndryCnd2D<true>( this, m_noSpecies );
  } else {
    m_strctrdBndryCnd = new ZFSStrctrdBndryCnd2D<false>( this, m_noSpecies );
  }

  // assign coordinates to all ghost points
  addGhostPointCoordinateValues();

  // allocate memory for aux data maps (cf,cp)
  //allocateAuxDataMaps();

  //if we are Rans we should allocate a new RANS block
  if(m_rans==true){
    m_ransBlck = new ZFSStrctrdBlck2DRans(this);
  }

  // allocate and compute metrics and jacobian
  allocateMetrics();
  allocateJacobian();
  computeCellCentreCoordinates();
  computeMetrics();
  computeJacobian();

  initFluxMethod();

  m_convergence=false;

  //Assign handlers to the correct boundary conditions
  assignBndryCells();

  if(m_rans) {
    if(domainId()==0) {
      cout << "Computing wall distance" << endl;
    }
    m_strctrdBndryCnd->computeWallDistances();
  }

  printAllocatedMemory( oldAllocatedBytes, "ZFSStrctrdBlck2D", m_zfsStrctrdComm );
} 

void ZFSStrctrdBlck2D::initFluxMethod()
{
  //set the MUSCL-scheme to the right function
  if(m_rans==true){
    reconstructSurfaceData=&ZFSStrctrdBlck2D::MusclRANS;
    viscFluxMethod=&ZFSStrctrdBlck2D::viscousFluxRANS;
  }else{
    viscFluxMethod=&ZFSStrctrdBlck2D::viscousFluxLES;
    if(m_limiter) {
      switch(string2enum(m_limiterMethod)) {
      case ALBADA: {
        zfs_log << "Using VAN ALBADA limiter!" << endl;
        reconstructSurfaceData=&ZFSStrctrdBlck2D::MusclAlbada;
        break;
      }
      default: {
        stringstream errorMessage;
        errorMessage << "Limiter function " << m_limiterMethod << " not implemented!" << endl;
        zfsTerm(1, __CALLING_FUNCTION__, errorMessage.str());
      }
      }
    } else {
      if(m_musclScheme=="Standard"){
        zfs_log << "Using unlimited MUSCL! (standard Formulation)" << endl;
        if(m_ausmScheme=="Standard"){
          switch(CV->noVariables){
          case 4:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES,4>;break;}
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES,5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES,6>;break;}
          default:{
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            zfsTerm(1, __CALLING_FUNCTION__);
          }
          }
        }else if(m_ausmScheme=="PTHRC"){
          switch(CV->noVariables){
          case 4:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,4>;break;}
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,6>;break;}
          default:{
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            zfsTerm(1, __CALLING_FUNCTION__);
          }
          }
        }else if(m_ausmScheme=="AUSMDV"){
          switch(CV->noVariables){
          case 4:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmDV,4>;break;}
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmDV,5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmDV,6>;break;}
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
        if(m_ausmScheme=="Standard"){
          switch(CV->noVariables){
          case 4:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES,4>;break;}
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES,5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES,6>;break;}
          default:{
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            zfsTerm(1, __CALLING_FUNCTION__);
          }
          }
        }else if(m_ausmScheme=="PTHRC"){
          switch(CV->noVariables){
          case 4:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,4>;break;}
          case 5:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,5>;break;}
          case 6:{ reconstructSurfaceData=&ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,6>;break;}
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
}

ZFSStrctrdBlck2D::~ZFSStrctrdBlck2D()
{
  TRACE();

  delete m_strctrdBndryCnd;
}


void ZFSStrctrdBlck2D::computeCellLength(){
  //this function can be moved into the MusclSchemeStreched later but for testing it is easier
  // REMEMBER: FOR MOVINg GRIDS THIS NEEDS TO BE CALLED EACH TIME

  for(ZFSId j=0; j<m_nCells[0]; j++){
    for(ZFSId i=0; i<m_nCells[1];i++){
      const ZFSId cellId=cellIndex(i,j);
      const ZFSId P1 = getPointIdFromCell(i,j);
      const ZFSId P2 = getPointIdFromPoint( P1, 1, 0);
      const ZFSId P3 = getPointIdFromPoint( P1, 1, 1);
      const ZFSId P4 = getPointIdFromPoint( P1, 0, 1);
      //----------Idirection
      //face 1
      const ZFSFloat f1x=F1B2*(m_coordinates[0][P1]+m_coordinates[0][P4]);
      const ZFSFloat f1y=F1B2*(m_coordinates[1][P1]+m_coordinates[1][P4]);
      //face 2
      const ZFSFloat f2x=F1B2*(m_coordinates[0][P2]+m_coordinates[0][P3]);
      const ZFSFloat f2y=F1B2*(m_coordinates[1][P2]+m_coordinates[1][P3]);
      m_cells->cellLength[0][cellId]=sqrt(POW2(f2x-f1x)+POW2(f2y-f1y));
      //----------Jdirection
      //face 3
      const ZFSFloat f3x=F1B2*(m_coordinates[0][P1]+m_coordinates[0][P2]);
      const ZFSFloat f3y=F1B2*(m_coordinates[1][P1]+m_coordinates[1][P2]);
      //face 4
      const ZFSFloat f4x=F1B4*(m_coordinates[0][P3]+m_coordinates[0][P4]);
      const ZFSFloat f4y=F1B4*(m_coordinates[1][P3]+m_coordinates[1][P4]);
      m_cells->cellLength[1][cellId]=sqrt(POW2(f4x-f3x)+POW2(f4y-f3y));
    }
  }
}


/** initalize the solution step
 * 
 * @author:Pascal Meysonnat
 * @date: 01.01.1010
 *
 */

void ZFSStrctrdBlck2D::initSolutionStep()
{
  TRACE();
  //Compute infinity values from property file
  //and (if no restart) fill cells according
  //to the initialCondition property
  initialCondition();

  if( m_restart ) {
    loadRestartFile();
  }
  //writeGridPointsWithGhostPoints();
  //saveOutput(true);
  //saveOutputPartitions<2>();
  //zfsTerm(-1, __CALLING_FUNCTION__, "Pascal killed it");
  setTimeStep();

  //Get the correct values
  //in the exchange ghostcells
  exchange();
  
  //Call the init function of each BC
  initBndryCnds();

  //Apply boundary conditions
  //and fill the non-exchange ghostcells
  applyBoundaryCondition();
  
  //Check for NaNs
  checkNans();

  computeConservativeVariables();
}

/**
 * Computation of infinity values for the conservative and primitive variables
 * Initialization ot the entire flow field
 * structure is based on Daniel Hartmann's function zfsfvblock3d.cpp initialCondition
 */

void ZFSStrctrdBlck2D::initialCondition()
{
  TRACE();

  const ZFSFloat gammaMinusOne = m_gamma - 1.0;
  ZFSFloat UT;
  //ZFSFloat pressure=F0;
  //ZFSFloat Frho;         switched off cause of compiler

  PV->TInfinity = 1.0 / ( 1.0 + F1B2 * gammaMinusOne * POW2(m_Ma));
  UT = m_Ma * sqrt( PV->TInfinity );
  PV->UInfinity = UT * cos( m_angle[ 0 ] ) * cos( m_angle[ 1 ] );
  PV->VInfinity = UT * sin( m_angle[ 0 ] ) * cos( m_angle[ 1 ] );
  PV->VVInfinity[ 0 ] = PV->UInfinity;
  PV->VVInfinity[ 1 ] = PV->VInfinity;
  PV->PInfinity = pow( PV->TInfinity, (m_gamma / gammaMinusOne)) / m_gamma;

  // compute conservative variables
  CV->rhoInfinity = pow( PV->TInfinity, ( 1.0 / gammaMinusOne ) );
  CV->rhoUInfinity = CV->rhoInfinity * PV->UInfinity;
  CV->rhoVInfinity = CV->rhoInfinity * PV->VInfinity;
  CV->rhoVVInfinity[ 0 ] = CV->rhoUInfinity;
  CV->rhoVVInfinity[ 1 ] = CV->rhoVInfinity;
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

  zfs_log << "**************************" << endl;
  zfs_log << "Initial Condition summary" << endl;
  zfs_log << "**************************" << endl;
  zfs_log << "Re = " << m_Re << endl;
  zfs_log << "Re0 = " << m_Re0 << endl;
  zfs_log << "Ma = " << m_Ma << endl;
  zfs_log << "TInfinity = " << PV->TInfinity << endl;
  zfs_log << "UInfinity = " << PV->UInfinity << endl;
  zfs_log << "VInfinity = " << PV->VInfinity << endl;
  zfs_log << "PInfinity = " << PV->PInfinity << endl;
  zfs_log << "rhoInfinity = " << CV->rhoInfinity << endl;
  zfs_log << "rhoEInfinity = " << CV->rhoEInfinity << endl;
  zfs_log << "referenceTime = " << m_timeRef << endl;

  if(domainId() == 0 )
    {
      cout << "**************************" << endl;
      cout << "Initial Condition summary" << endl;
      cout << "**************************" << endl;
      cout << "Re = " << m_Re << endl;
      cout << "Re0 = " << m_Re0 << endl;
      cout << "Ma = " << m_Ma << endl;
      cout << "TInfinity = " << PV->TInfinity << endl;
      cout << "UInfinity = " << PV->UInfinity << endl;
      cout << "VInfinity = " << PV->VInfinity << endl;
      cout << "PInfinity = " << PV->PInfinity << endl;
      cout << "rhoInfinity = " << CV->rhoInfinity << endl;
      cout << "rhoEInfinity = " << CV->rhoEInfinity << endl;
      cout << "referenceTime = " << m_timeRef << endl;
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
        for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
          //go through every cell
          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          for( ZFSId i = 0; i < nDim; i++ ) {
            m_cells->pvariables[PV->VV[i]][cellId] = PV->VVInfinity[i];
          }

          m_cells->pvariables[PV->P][cellId]= PV->PInfinity;

          if(m_rans) {
            m_cells->pvariables[PV->RANS_VAR[0]][cellId] = PV->ransInfinity[0];
          }
        }
        break;
      }
      case 43 :
        {
          //parallel inflow field
          for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
              //go through every cell
              m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
              for( ZFSId i = 0; i < nDim; i++ ) {
                m_cells->pvariables[PV->VV[i]][cellId] = F0;
              }

              m_cells->pvariables[PV->P][cellId]= PV->PInfinity;

               ZFSFloat radius = sqrt(POW2(m_cells->coordinates[0][cellId] - 0.5) + POW2(m_cells->coordinates[1][cellId] - 0.5));
               //impose pressure peak in the middle of the domain
               if(radius <= 0.05) {
                 ZFSFloat pAmp = 0.005;
                 ZFSFloat pressureSignal = sin(radius/0.05*PI)*pAmp + PV->PInfinity;
                 m_cells->pvariables[PV->P][cellId]= pressureSignal;
              }
            }
          break;
        }
      case 333 :
      {
        //parallel inflow field
        for(ZFSId cellId=0; cellId<m_noStrctrdCells; cellId++) {
          //go through every cell
          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          for( ZFSId i = 0; i < nDim; i++ ) {
            m_cells->pvariables[PV->VV[i]][cellId] = PV->VVInfinity[i];
          }

          m_cells->pvariables[PV->P][cellId]= PV->PInfinity;

          //impose pressure peak in the middle of the domain
          if(m_cells->coordinates[0][cellId] > 0.4 && m_cells->coordinates[0][cellId] < 0.5) {
            ZFSFloat pAmp = 0.005;
            ZFSFloat xCoordinate = m_cells->coordinates[0][cellId] - 0.4;
            ZFSFloat pressureSignal = sin(xCoordinate/0.1*PI)*pAmp + PV->PInfinity;
            m_cells->pvariables[PV->P][cellId]= pressureSignal/(m_gamma-F1);
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

          for(ZFSId j=0; j<m_nCells[0]; j++) {
            for(ZFSId i=0; i<m_nCells[1]; i++) {
              const ZFSId cellId=cellIndex(i,j);
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
                                                                    + 2*PI1/K*(3*eta*eta - 2*eta*eta*eta)), PV->UInfinity);
              }

              m_cells->pvariables[PV->V][cellId] =PV->VInfinity;
              m_cells->pvariables[PV->RHO][cellId] =CV->rhoInfinity;
              m_cells->pvariables[PV->P][cellId] = PV->PInfinity;

              if(m_rans) {
                m_cells->pvariables[PV->RANS_VAR[0]][cellId] = PV->ransInfinity[0];
              }
            }
          }

          break;
      }

      default:
      {
        //put the parallel flow field input in here
        //force output that no specific initial condition was chosen
        zfs_log << "No (correct) initial Condition is given! Used initial Condtion of parallel inflow!!!!!" << endl;
        break;
      }
      }
  }
}

void ZFSStrctrdBlck2D::initMovingGrid()
{
  TRACE();

  ZFSId pointId=0;

  //First approach: save whole mesh in m_mgInitCoordinates (for analytical channel with indentation)

  for(ZFSId j=0; j<m_nPoints[0]; ++j) {
    for(ZFSId i=0; i<m_nPoints[1]; ++i) {
      pointId= pointIndex(i,j);
      for( ZFSInt isd = xsd; isd < nDim; ++isd) {
        m_mgInitCoordinates[isd][pointId] = m_coordinates[isd][pointId];
      }
    }
  }

}

void ZFSStrctrdBlck2D::assignBndryCells()
{
  TRACE();
  m_strctrdBndryCnd->assignBndryCnds();
}

void ZFSStrctrdBlck2D::initBndryCnds()
{
  TRACE();
  m_strctrdBndryCnd->correctBndryCndIndices();
}

void ZFSStrctrdBlck2D::applyBoundaryCondition()
{
  TRACE();
  // treat Dirichlet and Neumann BC in one go!!!
  m_strctrdBndryCnd->applyDirichletNeumannBC();
}

void ZFSStrctrdBlck2D::computeCellCentreCoordinates()
{
  //function to compute the coordinates at cell centre
  //calculated over I, J loop but changed to one array
  for(ZFSId j=0; j < m_nCells[0]; ++j){
    for(ZFSId i=0; i < m_nCells[1]; ++i){
      const ZFSId IJ = getPointIdFromCell(i, j);
      const ZFSId IP1J= getPointIdFromPoint(IJ,1,0);
      const ZFSId IJP1= getPointIdFromPoint(IJ,0,1);
      const ZFSId IP1JP1= getPointIdFromPoint(IJ,1,1);
      const ZFSId cellId = cellIndex(i, j);

      for(ZFSId dim = 0; dim < nDim; dim++){
        //average the coordinates for cell centre data
        m_cells->coordinates[dim][cellId] = F1B4*(m_coordinates[dim][IJ]+m_coordinates[dim][IP1J]+m_coordinates[dim][IJP1]+m_coordinates[dim][IP1JP1]);
      }
    }
  }
}

bool ZFSStrctrdBlck2D::maxResidual()
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
  //ZFSInt accumCounter=0;
  for(ZFSId dim=0; dim<nDim; dim++){
    maxResIndex[dim]=F0;
  }

  if(!m_localTimeStep) {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++) {
      for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; i++) {
        cellId=cellIndex(i,j);
        tmpResidual = m_timeStep / (m_cfl *m_cells->cellJac[cellId] )* fabs( m_cells->rightHandSide[CV->RHO][cellId]);
        m_avrgResidual += tmpResidual;
        if( tmpResidual > maxResidual1 ) {
          maxResIndex[0]=i-m_noGhostLayers;
          maxResIndex[1]=j-m_noGhostLayers;
          maxResidual1=tmpResidual;
        }
      }
    }
  } else {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++) {
      for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; i++) {
        cellId=cellIndex(i,j);
        tmpResidual = m_cells->localTimeStep[cellId] / (m_cfl *m_cells->cellJac[cellId] )* fabs( m_cells->rightHandSide[CV->RHO][cellId]);
        m_avrgResidual += tmpResidual;
        if( tmpResidual > maxResidual1 ) {
          maxResIndex[0]=i-m_noGhostLayers;
          maxResIndex[1]=j-m_noGhostLayers;
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

  MPI_Allreduce(&localAvrgResidual, &accumAvrgResidual, 1, MPI_DOUBLE, MPI_SUM, m_zfsStrctrdComm);
  MPI_Allreduce(&localMaxResidual, &globalMaxResidual, 1, MPI_DOUBLE, MPI_MAX, m_zfsStrctrdComm);
  m_avrgResidual = accumAvrgResidual;//m_residualRcv.avrgRes;
  maxResidualOrg=globalMaxResidual;
  m_avrgResidual=m_avrgResidual/m_totalGridCells;

  // write first residuals;
  if( fabs(m_firstMaxResidual) < epsilon ) {
    m_firstMaxResidual = zfsMAX(epsilon,globalMaxResidual);
    m_firstAvrgResidual =zfsMAX(epsilon,m_avrgResidual);
    if(approx(localMaxResidual, maxResidualOrg, m_eps)){
      //write out values into residual file
      FILE* f_residual;
      f_residual = fopen("./Residual", "a+");
      fprintf(f_residual, "#MaxRes_1: %1.10e \n", m_firstMaxResidual);
      fprintf(f_residual, "#MaxAvgRes_1: %1.10e \n",m_firstAvrgResidual );
      fprintf(f_residual, "#iter, physTime, time, dT, wLoad, avrgRes, maxRes, blockId, i, j");
      fclose(f_residual);
    }
  }

  // normalize residuals
  globalMaxResidual = globalMaxResidual / m_firstMaxResidual;
  m_avrgResidual = (m_avrgResidual / m_firstAvrgResidual);

  if(std::isnan(m_avrgResidual)) {
    cerr << "Solution diverged, average residual is nan " << endl;
    zfs_log << "Solution diverged, average residual is nan " << endl;
    saveOutput(true);
    savePartitions();
    zfsTerm(1,__CALLING_FUNCTION__,"Solution diverged, average residual is nan ");
  }

  //convergence Check
  m_convergence=false;
  if( maxResidual1 < m_convergenceCriterion ) {
    m_convergence = true;
  }

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
    fprintf(f_residual, " %d", m_nOffsetCells[1]+maxResIndex[0]);//i
    fprintf(f_residual, " %d", m_nOffsetCells[0]+maxResIndex[1]);//j
    fprintf(f_residual, "\n");
    fclose(f_residual);
  }

  if( maxResidual1 < m_convergenceCriterion ) {
    return true;
  } else {
    return false;
  }
}

inline ZFSFloat ZFSStrctrdBlck2D::crossProduct( ZFSFloat vec1[2], ZFSFloat vec2[2])
{
  ZFSFloat result = vec1[xsd] * vec2[ysd] - vec1[ysd] * vec2[xsd];
  return result;
}


inline ZFSId ZFSStrctrdBlck2D::getPointIdFromCell( ZFSInt i, ZFSInt j)
{
  return i + (j * (m_nCells[1] + 1));
}

inline ZFSId ZFSStrctrdBlck2D::getPointIdFromPoint( ZFSId origin, ZFSInt incI,
                                                    ZFSInt incJ )
{
  return origin + incI + incJ * m_nPoints[1];

}

inline ZFSId ZFSStrctrdBlck2D::getCellIdFromCell( ZFSId origin, ZFSInt incI, ZFSInt incJ )
{
  return origin + incI + incJ * m_nCells[1];
}

inline ZFSId ZFSStrctrdBlck2D::cellIndex(ZFSInt i, ZFSInt j)
{
  return i + (j * m_nCells[1]);
}

inline ZFSId ZFSStrctrdBlck2D::pointIndex(ZFSInt i, ZFSInt j)
{
  return i + (j * m_nPoints[1]);
}

ZFSFloat ZFSStrctrdBlck2D::pressure(ZFSId cellId){
  return m_cells->pvariables[PV->P][cellId];
}

void ZFSStrctrdBlck2D::addGhostPointCoordinateValues()
{
  TRACE();

  if(m_debugOutput){
    for(ZFSId j=0; j<(m_nCells[0]); j++ ){
      for(ZFSId i=0; i<(m_nCells[1]); i++){
        ZFSId cellId = cellIndex(i,j);
        m_cells->fq[FQ->CELLID][cellId]=cellId;
        m_cells->fq[FQ->BLOCKID][cellId]=domainId();
      }
    }
  }

  //1) extrapolate GhostPoints
  extrapolateGhostPointCoordinates();
  //2) communicate GhostPoints

  if(noDomains()>1) exchangePoints();

  extrapolateGhostPointCoordinatesBC();

  //3) write the totalGridFile with GhostPoints
  if(m_savePartitionOutput) {
    writeGridPointsWithGhostPoints();
  }
}


void ZFSStrctrdBlck2D::extrapolateGhostPointCoordinates()
{
  //This function mirrors the grid points on the faces
  //Can this be written in a more efficient way????
  TRACE();


  //i-direction
  ZFSId pointId, FixPointId, MirrorPointId;

  for(ZFSId j = m_noGhostLayers; j < ( m_nPoints[0] - m_noGhostLayers ); ++j )
    {
      for(ZFSId i = 0; i < m_noGhostLayers; ++i)
        {
          pointId = ( m_noGhostLayers - 1 -i ) + ( j * m_nPoints[1] ); //pointId in Array
          FixPointId = ( m_noGhostLayers - i ) + ( j * m_nPoints[1] ); //point about which everything is mirrored
          MirrorPointId = ( m_noGhostLayers + 1 - i ) + ( j * m_nPoints[1] );

          for(ZFSId dim =0; dim < nDim; ++dim)
            {
              m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
            }

          // if(j == 2) {
          //    cout << "dim: " << nDim << " xPointId: " << m_coordinates[0][pointId] << " xFixPoint: " << m_coordinates[0][FixPointId] << " xMirror: " << m_coordinates[0][MirrorPointId] << " pointId: " << pointId << " fixPointId: " << FixPointId << " mirrorPointId: " << MirrorPointId << endl;
          // }
          //coordinates at the other end!!

          pointId = (m_nPoints[1]-i-1)+(j*m_nPoints[1]);
          FixPointId = (m_nPoints[1]-m_noGhostLayers-1)+(j*m_nPoints[1]);
          MirrorPointId = (m_nPoints[1]-1-(2*m_noGhostLayers - i))+(j*m_nPoints[1]);
          for(ZFSId dim =0; dim < nDim; ++dim)
            {
              m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
            }
        }
    }

  //j-direction

  for(ZFSId j = 0; j < m_noGhostLayers; ++j )
    {
      for(ZFSId i = m_noGhostLayers; i < ( m_nPoints[1] - m_noGhostLayers ); ++i)
        {
          pointId = i+(m_noGhostLayers-j-1)*m_nPoints[1]; //pointId in Array
          FixPointId = i+(m_noGhostLayers-j)*m_nPoints[1]; //point about which everything is mirrored
          MirrorPointId = i+(m_noGhostLayers+1-j)*m_nPoints[1];
          for(ZFSId dim =0; dim < nDim; ++dim)
            {
              m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
            }
          //coordinates at the other end!!
          pointId = i+(m_nPoints[0]-j-1)*m_nPoints[1];
          FixPointId =i+ (m_nPoints[0]-m_noGhostLayers-1)*m_nPoints[1];
          MirrorPointId =i+ (m_nPoints[0]-1-(2*m_noGhostLayers - j))*m_nPoints[1];
          for(ZFSId dim =0; dim < nDim; ++dim)
            {
              m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
            }
        }
    }


  //corner points missing yet!! They only need to be calculated if a Visualisation tool is used to
  //show the ghost points and the grid


  //in i-direction

  for(ZFSId j=0; j < ( m_nPoints[0] ); ++j )
    {
      for(ZFSId i=0; i < m_noGhostLayers; ++i)
        {

          pointId = (m_noGhostLayers -1 -i)+(j*m_nPoints[1]); //pointId in Array
          FixPointId = (m_noGhostLayers-i)+(j*m_nPoints[1]); //point about which everything is mirrored
          MirrorPointId = (m_noGhostLayers+1 - i)+(j*m_nPoints[1]); //m_noGhostLayers+(m_noGhostLayers-i)
          for(ZFSId dim =0; dim < nDim; ++dim)
            {
              m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
            }
          //coordinates at the other end!!

          pointId = (m_nPoints[1]-i-1)+(j*m_nPoints[1]);
          FixPointId = (m_nPoints[1]-m_noGhostLayers-1)+(j*m_nPoints[1]);
          MirrorPointId = (m_nPoints[1]-1-(2*m_noGhostLayers - i))+(j*m_nPoints[1]);
          for(ZFSId dim =0; dim < nDim; ++dim)
            {
              m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
            }

        }
    }

  //in j-direction

  for(ZFSId j=0; j < m_noGhostLayers; ++j )
    {
      for(ZFSId i=0; i< ( m_nPoints[1] ); ++i)
        {

          pointId = i + ( m_noGhostLayers - j - 1 ) * m_nPoints[1]; //pointId in Array
          FixPointId = i + ( m_noGhostLayers - j ) * m_nPoints[1]; //point about which everything is mirrored
          MirrorPointId = i + ( m_noGhostLayers + 1-j ) * m_nPoints[1];
          for(ZFSId dim =0; dim < nDim; ++dim)
            {
              m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
            }
          //coordinates at the other end!!
          pointId = i+(m_nPoints[0]-j-1)*m_nPoints[1];
          FixPointId =i+ (m_nPoints[0]-m_noGhostLayers-1)*m_nPoints[1];
          MirrorPointId =i+ (m_nPoints[0]-1-(2*m_noGhostLayers - j))*m_nPoints[1];
          for(ZFSId dim =0; dim < nDim; ++dim)
            {
              m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
            }
        }
    }
}

void ZFSStrctrdBlck2D::extrapolateGhostPointCoordinatesBC()
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
    ZFSInt pos[2],fix[2],mirror[2],ij[2],extendij[2];
    ZFSInt pointId,FixPointId,MirrorPointId;

    extendij[0]=1;extendij[1]=1;;
    extendij[index]=0;

      for(ij[1]=start[1]; ij[1]<end[1]+extendij[1]; ++ij[1]) {
        for(ij[0]=start[0]; ij[0]<end[0]+extendij[0]; ++ij[0]) {
          for(ZFSInt m=0;m<2;++m) {
            if(index==m) {
              if(step==1) {
                pos[m]=ij[m]+1;
                fix[m]=start[m];
                mirror[m]=2*fix[m]-pos[m];
              } else {
                pos[m]=ij[m];
                fix[m]=end[m];
                mirror[m]=2*fix[m]-pos[m];
              }
            } else {
              pos[m]=ij[m];
              fix[m]=ij[m];
              mirror[m]=ij[m];
            }
          }//m

          pointId       =pointIndex(pos[0],pos[1]);
          FixPointId    =pointIndex(fix[0],fix[1]);
          MirrorPointId =pointIndex(mirror[0],mirror[1]);

          for(ZFSId dim =0; dim < nDim; dim++) {
            m_coordinates[dim][pointId]=(2*m_coordinates[dim][FixPointId]-m_coordinates[dim][MirrorPointId]);
          }
        }//ij
      }
  }//bcid
}


void ZFSStrctrdBlck2D::computeSurfaceMetrics(){
  TRACE();
  for (ZFSInt j = 0; j < this->m_nCells[0]-1; j++) {
    for(ZFSInt i = 0 ; i< this->m_nCells[1]-1; i++) {
      //determine global cellID
      ZFSId cellId = this -> cellIndex(i, j);
      //determine global point ID for local cell IDs
      ZFSId IJ = getPointIdFromCell (i,j);
      ZFSId IPJ = getPointIdFromPoint (IJ , 1, 0);
      ZFSId IPJP = getPointIdFromPoint (IJ , 1 , 1);
      ZFSId IJP = getPointIdFromPoint (IJ , 0 , 1);

      //auxiliary variables
      ZFSFloat DcoordDxi [2];

      // Face I //
      for(ZFSInt isd = xsd; isd < nDim; ++isd) {
        DcoordDxi[isd]  = m_coordinates[isd][IPJP]- m_coordinates[isd][IPJ];
      }

      //compute Dxi
      m_cells->surfaceMetrics[cellId][0] = DcoordDxi[1];
      m_cells->surfaceMetrics[cellId][1] = -DcoordDxi[0];
      //store Dxi

      // Face I //
      for(ZFSInt isd = xsd; isd < nDim; ++isd) {
        DcoordDxi[isd]  = m_coordinates[isd][IPJP]- m_coordinates[isd][IJP];
      }

      m_cells->surfaceMetrics[cellId][2] = -DcoordDxi[1];
      m_cells->surfaceMetrics[cellId][3] = DcoordDxi[0];
    }
  }
}

void ZFSStrctrdBlck2D::computeModSurfaceMetrics(){
  TRACE();

  computeSurfaceMetrics();
}

void ZFSStrctrdBlck2D::computeCornerMetrics(){
  TRACE();
  const ZFSInt noCells = m_noStrctrdCells;

  for( ZFSInt j = m_noGhostLayers - 1; j < this->m_nCells[0] - m_noGhostLayers; ++j )
    {
      for( ZFSInt i = m_noGhostLayers -1; i < this->m_nCells[1] - m_noGhostLayers; ++i )
        {
          // determine global cell ID
          ZFSId cellId = this->cellIndex(i, j);

          // auxilliary variables
          ZFSFloat   DcoordDxi[2];
          ZFSFloat  DcoordDeta[2];

          for( ZFSInt isd = xsd; isd < nDim; ++isd )
            {
              // looks complicated, but what happens is that we always catch the point Id of ipjp
              // from the neighboring cell and build the centered difference

              //compute d(x,y,z)/dxi
              DcoordDxi[isd] = F1B2 *
                ( m_coordinates[isd][ getPointIdFromPoint( getPointIdFromCell(i + 1, j ), 1, 1 ) ] -
                  m_coordinates[isd][ getPointIdFromPoint( getPointIdFromCell(i - 1, j ), 1, 1 ) ] );

              //compute d(x,y,z)/deta
              DcoordDeta[isd] = F1B2 *
                ( m_coordinates[isd][ getPointIdFromPoint( getPointIdFromCell(i, j + 1 ), 1, 1 ) ] -
                  m_coordinates[isd][ getPointIdFromPoint( getPointIdFromCell(i, j - 1 ), 1, 1 ) ] );

            }
          m_cells->cornerMetrics[0*noCells+cellId ] =   DcoordDeta[1]; //DxiDx
          m_cells->cornerMetrics[1*noCells+cellId ] = - DcoordDeta[0]; //DxiDy
          m_cells->cornerMetrics[2*noCells+cellId ] = - DcoordDxi[1]; //DetaDx
          m_cells->cornerMetrics[3*noCells+cellId ] =   DcoordDxi[0]; //DetaDy
        }
    }
}

void ZFSStrctrdBlck2D::computeModCornerMetrics(){
  TRACE();

  //has to be implemented
  computeCornerMetrics();
}

void ZFSStrctrdBlck2D::computeCellMetrics(){
  TRACE();

  for( ZFSInt j = m_noGhostLayers - 1; j < this->m_nCells[0] - 1; ++j )
    {
      for( ZFSInt i = m_noGhostLayers - 1; i< this->m_nCells[1] - 1; ++i )
        {
          ZFSId cellId = this->cellIndex( i, j );
          // auxilliary variables
          ZFSFloat   DcoordDxi[ 2 ];
          ZFSFloat  DcoordDeta[ 2 ];

          for( ZFSInt isd = xsd; isd < nDim; isd++)
            {
              DcoordDxi[isd]   = ( m_cells->coordinates[ isd ][ cellIndex( i + 1, j ) ] -
                                   m_cells->coordinates[ isd ][ cellIndex( i - 1, j ) ] ) * F1B2;

              DcoordDeta[isd]  = ( m_cells->coordinates[ isd ][ cellIndex( i, j + 1 ) ] -
                                   m_cells->coordinates[ isd ][ cellIndex( i, j - 1 ) ] ) * F1B2;
            }

          m_cells->cellMetrics[cellId][ 0 ] =   DcoordDeta[1]; //DxiDx
          m_cells->cellMetrics[cellId][ 1 ] = - DcoordDeta[0]; //DxiDy
          m_cells->cellMetrics[cellId][ 2 ] = - DcoordDxi[1]; //DetaDx
          m_cells->cellMetrics[cellId][ 3 ] =   DcoordDxi[0]; //DetaDy
        }
    }
}

// jacobian for the viscous fluxes
void ZFSStrctrdBlck2D::computeCornerJacobian(){
  TRACE();
  //jacobian in the physical space is the inverse of the jacobian in computational space
  const ZFSInt noCells = m_noStrctrdCells;

  for(ZFSInt j = m_noGhostLayers -1 ; j < this->m_nCells[0] - m_noGhostLayers; ++j)
    {
      for(ZFSInt i = m_noGhostLayers -1 ; i < this->m_nCells[1] - m_noGhostLayers; ++i)
        {
          ZFSId cellId = cellIndex(i,j);

          ZFSFloat cornerJac =
            m_cells->cornerMetrics[ (xsd * 2 + xsd)*noCells + cellId ] *
            m_cells->cornerMetrics[ (ysd * 2 + ysd)*noCells + cellId ] -
            m_cells->cornerMetrics[ (ysd * 2 + xsd)*noCells + cellId ] *
            m_cells->cornerMetrics[ (xsd * 2 + ysd)*noCells + cellId ];

          // since metric terms are with omitted jacobian
          // there is factor of J^3; multiplied with J^-1 (invJac) we get J^2
          // --> take square root to get J
          this->m_cells->cornerJac[cellId] = cornerJac;
        }
    }
}

// same as TFS routine JACGP, more exact calculation of control volume
// for corner point than with corner metrics
void ZFSStrctrdBlck2D::computeModCornerJacobian(){
  TRACE();


  ZFSFloatScratchSpace subJ(this->m_noStrctrdCells, 4, __CALLING_FUNCTION__, "subJ");
  ZFSFloatScratchSpace subJtmp(this->m_noStrctrdCells, 4, __CALLING_FUNCTION__, "subJtmp");

  ZFSFloat** __restrict coords = m_coordinates;

  subJ.fill(1234.56);
  subJtmp.fill(5678.9);

  for(ZFSInt j = m_noGhostLayers - 2 ; j < this->m_nCells[0] - m_noGhostLayers; ++j ){

    for(ZFSInt i = m_noGhostLayers -2 ; i < this->m_nCells[1] - m_noGhostLayers; ++i ){

      const  ZFSId cellId = cellIndex( i, j ) ;
      const  ZFSId centCellId = cellIndex( i+1 , j+1 );
      const  ZFSId tmpId = getPointIdFromCell( i , j );

      const  ZFSId ij =   getPointIdFromPoint( tmpId, 1 , 1);
      const  ZFSId ipj =  getPointIdFromPoint( ij, 1, 0 ) ;
      const  ZFSId ijp =  getPointIdFromPoint( ij, 0, 1 ) ;
      const  ZFSId ipjp = getPointIdFromPoint( ij, 1, 1 ) ;

      // auxilliary variables for surface values
      ZFSFloat S1[2];
      ZFSFloat S2[2];
      ZFSFloat S1P[2];
      ZFSFloat S2P[2];


      for( ZFSInt isd = xsd; isd < nDim; isd++ ){

        S1[isd] = F1B2 * ( coords[isd][ijp] + coords[isd][ij] );
        S2[isd] = F1B2 * ( coords[isd][ipj] + coords[isd][ij] );
        S1P[isd] = F1B2 * ( coords[isd][ipjp] + coords[isd][ipj] );
        S2P[isd] = F1B2 * ( coords[isd][ipjp] + coords[isd][ijp] );
      }


      // vectors for metrics
      ZFSFloat tmpX1[2];
      ZFSFloat tmpX2[2];

      ///////////////////////////////////////////
      ////////// subjacobian 1 //////////////////
      ///////////////////////////////////////////

      for( ZFSInt isd = xsd; isd < nDim; isd++ )
        {
          //averaging corner points is in 2D not necessary, cause of calculated surfaces mid points above

          //setting up vectors for new metric terms
          tmpX1[isd] = (m_cells->coordinates[isd][centCellId] - coords[isd][ij])/ sqrt(2) ;
          tmpX2[isd] = (S1[isd] - S2[isd])/ sqrt(2);
        }

      subJ(cellId,0) = F0;
      subJ(cellId,0)= crossProduct( tmpX1, tmpX2 );

      ///////////////////////////////////////////
      ////////// subjacobian 2 //////////////////
      ///////////////////////////////////////////

      for( ZFSInt isd = xsd; isd < nDim; isd++ ){

        tmpX1[isd] = ( S1P[isd] - S2[isd] ) / sqrt(2);
        tmpX2[isd] = ( m_cells->coordinates[isd][centCellId] - m_coordinates[isd][ipj] )/ sqrt(2);
      }

      subJ(cellId,1) = F0;
      subJ(cellId,1)= crossProduct( tmpX1, tmpX2 );

      ///////////////////////////////////////////
      ////////// subjacobian 3 //////////////////
      ///////////////////////////////////////////

      for( ZFSInt isd = xsd; isd < nDim; isd++ ){

        tmpX1[isd] = ( S2P[isd] - S1[isd] )/ sqrt(2);
        tmpX2[isd] = (  m_coordinates[isd][ijp] - m_cells->coordinates[isd][centCellId] )/ sqrt(2);
      }

      subJ(cellId,2)= F0;
      subJ(cellId,2)= crossProduct( tmpX1, tmpX2 );

      ///////////////////////////////////////////
      ////////// subjacobian 4 //////////////////
      ///////////////////////////////////////////

      for( ZFSInt isd = xsd; isd < nDim; isd++ ){

        tmpX1[isd] = ( m_coordinates[isd][ipjp] - m_cells->coordinates[isd][centCellId] ) / sqrt(2);
        tmpX2[isd] = ( S2P[isd] - S1P[isd] ) / sqrt(2);
      }

      subJ(cellId,3) = F0;
      subJ(cellId,3)= crossProduct( tmpX1, tmpX2 );



    }
  }



  //////////////////////////////////////////////
  ///// assemble subjacobians //////////////////
  //////////////////////////////////////////////

  //copy into dummy array
  for( ZFSInt i = 0; i < m_noStrctrdCells; i++ )
    {
      for( ZFSInt j = 0; j < 4; j++ )
        {
          subJtmp(i,j) = subJ(i,j);
        }
    }

  //shift subjacobians

  for( ZFSInt j = m_noGhostLayers-1; j < this->m_nCells[0] - m_noGhostLayers; j++ )
    {
      for( ZFSInt i = m_noGhostLayers-1; i < this->m_nCells[1] - m_noGhostLayers; i++ )
        {
          ZFSId cellId = cellIndex( i, j );

          subJ(cellId,0) = subJtmp( cellIndex(i-1, j-1) , 3 );
          subJ(cellId,1) = subJtmp( cellIndex(i  , j-1) , 2 );
          subJ(cellId,2) = subJtmp( cellIndex(i-1, j  ) , 1 );
          subJ(cellId,3) = subJtmp( cellIndex(i  , j  ) , 0 );
        }
    }


  //finally jacobian at corner point!
  for( ZFSInt j = m_noGhostLayers-1; j < this->m_nCells[0] - m_noGhostLayers; j++ )
    {
      for( ZFSInt i = m_noGhostLayers-1; i < this->m_nCells[1] - m_noGhostLayers; i++ )
        {
          ZFSId cellId = cellIndex( i, j );

          m_cells->cornerJac[cellId] = F0;
          for( ZFSInt jacId = 0; jacId < 4; jacId ++)
            {
              m_cells->cornerJac[cellId] += subJ(cellId,jacId);
            }


          m_cells->cornerJac[cellId] =  m_cells->cornerJac[cellId];
        }
    }
}

void ZFSStrctrdBlck2D::computeCellJacobian(){
  TRACE();
  // Jacobian: Dxi/Dx * Deta/Dy - Deta/Dx * Dxi/Dy

  for(ZFSInt j = m_noGhostLayers -1 ; j < this->m_nCells[0] - 1 ; ++j)
    {
      for(ZFSInt i = m_noGhostLayers -1 ; i < this->m_nCells[1] - 1; ++i)
        {
          ZFSId cellId = cellIndex(i,j);

          ZFSFloat cellJac =
            m_cells->cellMetrics[ cellId ][ xsd * 2 + xsd ] *
            m_cells->cellMetrics[ cellId ][ ysd * 2 + ysd ] -
            m_cells->cellMetrics[ cellId ][ ysd * 2 + xsd ] *
            m_cells->cellMetrics[ cellId ][ xsd * 2 + ysd ];

          this->m_cells->cellJac[cellId] = cellJac;
        }
    }
}

void ZFSStrctrdBlck2D::computeModCellJacobian(){
  TRACE();

  //has to be implemented
  computeCellJacobian();
}

void ZFSStrctrdBlck2D::Muscl(ZFSInt zfsNotUsed(timerId))
{
  TRACE();

  //when this movingGrid part is implemented, put this method in the upper class ZFSStrdtrdBlck, because in 2D and 2D the methods
  //are the same.


  // if (m_movingGrid)
  //   {
  //     if(m_RKStep == 0)
  //    {
  //      saveGrid();
  //    }

  //     //moving the grid
  //     moveGrid(m_gridMovingMethod);

  //     //compute the volume fluxes
  //     computeDxt();
  //   }

  (this->*reconstructSurfaceData)();
}



//Muscl reconstruction with Albada limiter
void ZFSStrctrdBlck2D::MusclAlbada(){
  TRACE();

  const ZFSUint noCells = m_noStrctrdCells;
  ZFSFloat* RESTRICT flux = ALIGNED_F(m_cells->flux);
  ZFSId cellId=0, IP1=0, IM1=0, IP2=0;
  //grid stretching factors
  ZFSFloat DS=F0, DSP1=F0, DSM1=F0, DSP=F0, DSM=F0;
  //stencil identifier
  ZFSId IJ[2] = {1, m_nCells[1]};

  //flow variables differences
  //ZFSFloat DQ=F0, DQP1=F0, DQM1=F0;
  ZFSFloatScratchSpace DQ(CV->noVariables, __CALLING_FUNCTION__, "DQ" );
  ZFSFloatScratchSpace DQP1(CV->noVariables, __CALLING_FUNCTION__, "DQP1" );
  ZFSFloatScratchSpace DQM1(CV->noVariables, __CALLING_FUNCTION__, "DQM1" );

  //left and right state
  ZFSFloat epsLim=m_eps;
  ZFSFloat smps=F0;
  ZFSFloat dummy=F0, dummy1=F0;

  ZFSFloat pIM2=F0, pIM1=F0,  pIP1=F0,  pIP2=F0;
  for(ZFSId i=0; i<CV->noVariables; i++)
    {
      DQ[i]=F0;
      DQP1[i]=F0;
      DQM1[i]=F0;
    }
  //reduce to onedimensional arrays
  ZFSFloat* __restrict x = &m_cells->coordinates[0][0];
  ZFSFloat* __restrict y = &m_cells->coordinates[1][0];

  ZFSFloat phi=F0, psi=F0, vel=F0;

  /////////IMPORTANT PARAMETER
  //ZFSFloat epsi=F1;
  //ZFSFloat kappa=F1B3;
  /////////END IMPORTANT PARAMETER
  for(ZFSId dim=0; dim<nDim; ++dim)
    {
      for(ZFSId j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers; ++j)
        {
          for(ZFSId i=m_noGhostLayers-1; i<m_nCells[1]-m_noGhostLayers; ++i)
            {
              //cell ids
              cellId=cellIndex(i, j);

              IP1=cellId+IJ[dim];
              IM1=cellId-IJ[dim];
              IP2=cellId+2*IJ[dim];

              //distances q_i+1 - q_i
              DS=sqrt(POW2(x[IP1]-x[cellId]) + POW2(y[IP1]-y[cellId]));
              //distances q_i - q_i-1
              DSM1=sqrt(POW2(x[cellId]-x[IM1]) + POW2(y[cellId]-y[IM1]));
              DSP1=sqrt(POW2(x[IP2]-x[IP1]) + POW2(y[IP2]-y[IP1]));
              //account for grid stretching
              //like tfs comment
              //DSP=F2*DS/POW2(DSP1+DS);
              //DSM=F2*DS/POW2(DSM1+DS);
              //like tfs
              DSP=DS/POW2(DSP1+DS);
              DSM=DS/POW2(DSM1+DS);


              for(ZFSId var=0; var < CV->noVariables; ++var){
                DQ[var]=m_cells->variables[var][IP1]-m_cells->variables[var][cellId];

                DQP1[var]=m_cells->variables[var][IP2]-m_cells->variables[var][IP1];

                DQM1[var]=m_cells->variables[var][cellId]-m_cells->variables[var][IM1];
                //limiter

              }
              vel=F0;
              for(ZFSId dim1=0; dim1<nDim; ++dim1){
                vel += POW2( m_cells->variables[ CV->RHO_VV[dim1] ][ IM1 ] / m_cells->variables[ CV->RHO ][ IM1 ]);
              }
              pIM2 = m_cells->variables[ CV->RHO_E ][ IM1 ] - F1B2 * m_cells->variables[ CV->RHO ][ IM1 ] * vel;
              vel=F0;
              for(ZFSId dim1=0; dim1 < nDim; ++dim1){
                vel += POW2( m_cells->variables[ CV->RHO_VV[dim1] ][ cellId ] / m_cells->variables[ CV->RHO ][ cellId ]);
              }
              pIM1 = m_cells->variables[ CV->RHO_E ][ cellId ] - F1B2 * m_cells->variables[ CV->RHO ][ cellId ] * vel;

              vel=F0;
              for(ZFSId dim1=0; dim1 < nDim; ++dim1){
                vel += POW2( m_cells->variables[ CV->RHO_VV[dim1] ][ IP2 ] / m_cells->variables[ CV->RHO ][ IP2 ]);
              }
              pIP2 = m_cells->variables[ CV->RHO_E ][ IP2 ] - F1B2 * m_cells->variables[ CV->RHO ][ IP2 ] * vel;
              vel=F0;
              for(ZFSId dim1=0; dim1 < nDim; ++dim1){
                vel += POW2( m_cells->variables[ CV->RHO_VV[dim1] ][ IP1 ] / m_cells->variables[ CV->RHO ][ IP1 ]);
              }
              pIP1 = m_cells->variables[ CV->RHO_E ][ IP1 ] - F1B2 * m_cells->variables[ CV->RHO ][ IP1 ] * vel;

              smps = DS*DSP1;

              dummy = fabs( pIM2 - F2 * pIM1 + pIP1 ) / ( pIM2 +F2 * pIM1 + pIP1 );
              dummy1 = fabs( pIM1 - F2 * pIP1 + pIP2 ) / ( pIM1 + F2 * pIP1 + pIP2 );

              psi = zfsMIN( F1 , F6 * zfsMAX( dummy , dummy1 ) );
              epsLim = zfsMAX( m_eps , pow( F1B2 * smps, F5) );


              for(ZFSId var=0; var< CV->noVariables; ++var){
                phi=F1B2-(F1B2-zfsMAX(F0,(DQP1[var]*DQM1[var]*smps+F1B2*epsLim)/(POW2(DQP1[var]*DS)+POW2(DQM1[var]*DSP1)+epsLim)))*psi;

                m_QLeft[var]= m_cells->variables[var][cellId]+DSM*(DSM1*DQ[var]+DS*DQM1[var])*phi;
                m_QRight[var]= m_cells->variables[var][IP1]-DSP*(DS*DQP1[var]+DSP1*DQ[var])*phi;

                //PHI(IP,IM,ID)=F2-(F2-MAX(F0,(DQ(IP,ID)*DQ(IM,ID)*SMSP+EPSMP2)
                //     &     /((DQ(IP,ID)*SM)**2+(DQ(IM,ID)*SP)**2+EPSMP)))*PSI
              }


              AusmLES(m_QLeft, m_QRight, dim, cellId); //Flux balance in AUSM


            }
        }


      //FLUX BALANCE
      for(ZFSId v=0; v < CV->noVariables; ++v)
        {
          for(ZFSId j=m_noGhostLayers; j < m_nCells[0]-m_noGhostLayers; ++j)
            {
              for(ZFSId i=m_noGhostLayers; i < m_nCells[1]-m_noGhostLayers; ++i)
                {
                  const ZFSId I=cellIndex(i,j);
                  IM1=I-IJ[dim];
                  m_cells->rightHandSide[v][I]+=flux[IM1+noCells*v]-flux[I+noCells*v];
                }
            }
        }

    }
}

void ZFSStrctrdBlck2D::MusclRANS(){
  m_ransBlck->Muscl();
}

template<ZFSStrctrdBlck2D::fluxmethod ausm,ZFSId noVars>
void ZFSStrctrdBlck2D::Muscl_(){
  TRACE();
  //stencil identifier
  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSId IJK[2]={1,m_nCells[1]};

  const ZFSFloat *const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const ZFSFloat *const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);

  const ZFSFloat *const RESTRICT cellVariables = ALIGNED_F(m_cells->pvariables[0]);
  ZFSFloat *const RESTRICT cellRhs= ALIGNED_MF(m_cells->rightHandSide[0]);
  ZFSFloat *const RESTRICT qleft= ALIGNED_MF(m_QLeft);
  ZFSFloat *const RESTRICT qright= ALIGNED_MF(m_QRight);
  ZFSFloat *const RESTRICT flux = ALIGNED_F(m_cells->flux);

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

        (this->*ausm)(m_QLeft, m_QRight, dim, I); //Flux balance in AUSM
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

//standard Ausm
template void ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES,5>();
template void ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES,6>();
template void ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES,7>();
//pthrc Ausm
template void ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,5>();
template void ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,6>();
template void ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,7>();
//ausm dv
template void ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmDV,5>();
template void ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmDV,6>();
template void ZFSStrctrdBlck2D::Muscl_<&ZFSStrctrdBlck2D::AusmDV,7>();


template<ZFSStrctrdBlck2D::fluxmethod ausm,ZFSId noVars>
void ZFSStrctrdBlck2D::MusclStretched_(){
TRACE();

  //stencil identifier
  const ZFSUint noCells = m_noStrctrdCells;
  const ZFSId IJK[2]={1,m_nCells[1]};
  const ZFSFloat *const RESTRICT cellVariables= ALIGNED_F(m_cells->pvariables[0]);
  const ZFSFloat *const RESTRICT cellLength= ALIGNED_F(m_cells->cellLength[0]);
  ZFSFloat *const RESTRICT cellRhs= ALIGNED_MF(m_cells->rightHandSide[0]);
  ZFSFloat *const RESTRICT qleft= ALIGNED_MF(m_QLeft);
  ZFSFloat *const RESTRICT qright= ALIGNED_MF(m_QRight);
  ZFSFloat *const RESTRICT flux = ALIGNED_F(m_cells->flux);
  /////////IMPORTANT PARAMETER
  const ZFSFloat phi =F1;
  const ZFSFloat kappa=F0;//F1B3;
  /////////END IMPORTANT PARAMETER
  for(ZFSId dim=0; dim<nDim; dim++) {
    const ZFSUint dimOffset = dim*m_noStrctrdCells;
    const ZFSFloat *const RESTRICT length = ALIGNED_F(cellLength+dimOffset);

    for(ZFSId j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers; j++) {
      for(ZFSId i=m_noGhostLayers-1; i<m_nCells[1]-m_noGhostLayers; i++) {
        const ZFSId I=cellIndex(i,j);
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

        (this->*ausm)(m_QLeft, m_QRight, dim, I); //Flux balance in AUSM
      }
    }

  //FLUX BALANCE
    for(ZFSUint v=0; v<noVars; v++) {
      for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++) {
        for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; i++) {
          const ZFSId I=cellIndex(i,j);
          const ZFSId IM1=I-IJK[dim];
          const ZFSUint offset = v*m_noStrctrdCells;
          ZFSFloat *const RESTRICT rhs = ALIGNED_F(cellRhs+offset);
          rhs[I]+= flux[IM1+noCells*v]-flux[I+noCells*v];

#ifdef ZFS_EXTRA_DEBUG
          ZFSFloat Flux = (m_cells->flux[IM1+v*noCells]-m_cells->flux[I+v*noCells]);
          convFluxOut[dim][v*m_noStrctrdCells+I]=Flux;
#endif
        }
      }
    }
  }
}
//standard Ausm
template void ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES,5>();
template void ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES,6>();
template void ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES,7>();
//pthrc
template void ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,5>();
template void ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,6>();
template void ZFSStrctrdBlck2D::MusclStretched_<&ZFSStrctrdBlck2D::AusmLES_PTHRC,7>();


void ZFSStrctrdBlck2D::Ausm()
{
  //Ausm routines have been moved and are called from inside Muscl (better performance)
}


/** \brief AUSM CENTRAL as in TFS
 *  can be used for moving grids, dxt term is included
 */
//inline void ZFSStrctrdBlck2D::AusmNew(ZFSFloat* QLeft, ZFSFloat* QRight, const ZFSId dim, const ZFSId I)
inline void ZFSStrctrdBlck2D::AusmLES(ZFSFloat* QLeft, ZFSFloat* QRight, const ZFSId dim, const ZFSId I)
{
  ZFSFloat pFactor[2]={F0,F0};
  const ZFSFloat gamma = m_gamma;
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

  const ZFSUint noCells = m_noStrctrdCells;
  ZFSFloat *const RESTRICT flux = ALIGNED_MF(m_cells->flux);

  flux[I+noCells*CV->RHO_U] = RHOU2 * ( UL + UR ) + AbsRHO_U2 * ( UL - UR ) + PLR * pFactor[ 0 ];
  flux[I+noCells*CV->RHO_V] = RHOU2 * ( VL + VR ) + AbsRHO_U2 * ( VL - VR ) + PLR * pFactor[ 1 ];
  flux[I+noCells*CV->RHO_E] = RHOU2 * ( e0 + e1)  + AbsRHO_U2 * ( e0 - e1 ) + PLR * dxdtau;
  flux[I+noCells*CV->RHO]   = RHOU;
}


/**
 *  Same AUSM scheme as AusmNew with additional damping controlled
 *  by the 4th order pressure derivative. Pressure needs to computed
 *  beforehand.
 *
 */
inline void ZFSStrctrdBlck2D::AusmLES_PTHRC(ZFSFloat* QLeft, ZFSFloat* QRight, ZFSId dim, ZFSId I)
{
 ZFSFloat pFactor[2]={F0,F0};
  const ZFSFloat gamma = m_gamma;
  const ZFSFloat gammaMinusOne = gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_cells->surfaceMetrics[ I ]);
  const ZFSFloat *const RESTRICT p = ALIGNED_F(m_cells->pvariables[ PV->P ]);

  const ZFSFloat dxdtau = m_cells->dxt[dim][I];

  //calculate pressure
  const ZFSFloat PL= QLeft[ PV->P];
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

  //compute splitting pressure
  const ZFSId IPJK  = getCellIdFromCell(I,  1,0);
  const ZFSId IMJK  = getCellIdFromCell(I, -1,0);
  const ZFSId IP2JK = getCellIdFromCell(I,  2,0);
  const ZFSId IM2JK = getCellIdFromCell(I, -2,0);

  const ZFSId IJPK = getCellIdFromCell(I, 0,1);
  const ZFSId IJMK = getCellIdFromCell(I, 0,-1);
  const ZFSId IJP2K = getCellIdFromCell(I, 0,2);
  const ZFSId IJM2K = getCellIdFromCell(I, 0,-2);

  const ZFSFloat p4I4 = F4 * (p[IPJK] + p[IMJK]) -F6 * (p[I]) - p[IP2JK] - p[IM2JK];
  const ZFSFloat p4J4 = F4 * (p[IJPK] + p[IJMK]) -F6 * (p[I]) - p[IJP2K] - p[IJM2K];

  const ZFSFloat cfac = 1.0/1.3;
  const ZFSFloat pfac = fabs(p4I4) + fabs(p4J4);
  ZFSFloat fac = cfac*pfac;
  fac = min(1/64.0,fac*5.0);

  const ZFSFloat PLR= PL*(F1B2+ fac * MAL) + PR*(F1B2- fac * MAR);

  const ZFSFloat RHO_AL = RHOL*AL;
  const ZFSFloat RHO_AR = RHOR*AR;

  const ZFSFloat PLfRHOL = PL/RHOL;
  const ZFSFloat PRfRHOR = PR/RHOR;

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
}

void ZFSStrctrdBlck2D::AusmDV(ZFSFloat* QLeft, ZFSFloat*  QRight, const ZFSId dim, const ZFSId I)
{
  ZFSFloat pFactor[2]={F0,F0};
  const ZFSFloat gamma = m_gamma;
  const ZFSFloat gammaMinusOne = gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  const ZFSFloat *const RESTRICT surf = ALIGNED_F(m_cells->surfaceMetrics[I]);

  //left side
  const ZFSFloat RHOL = QLeft[PV->RHO];
  const ZFSFloat FRHOL = F1/RHOL;
  ZFSFloat UL = QLeft[PV->U];
  ZFSFloat VL = QLeft[PV->V];
  const ZFSFloat PL = QLeft[PV->P];

  //right side
  const ZFSFloat RHOR = QRight[PV->RHO];
  const ZFSFloat FRHOR = F1/RHOR;
  ZFSFloat UR = QRight[PV->U];
  ZFSFloat VR = QRight[PV->V];
  const ZFSFloat PR = QRight[PV->P];

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
}

void ZFSStrctrdBlck2D::computeTimeStep()
{
  TRACE();
  m_timeStep=1000.0;
  const ZFSFloat *const RESTRICT dxtx = ALIGNED_F(m_cells->dxt[0]);
  const ZFSFloat *const RESTRICT dxty = ALIGNED_F(m_cells->dxt[1]);

  for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++)
    {
      for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers;i++)
        {
          const ZFSId cellId = cellIndex(i,j);
          const ZFSFloat *const RESTRICT metric = ALIGNED_F(m_cells->cellMetrics[ cellId ]);
          const ZFSFloat Frho = F1 / m_cells->pvariables[ PV->RHO ][cellId];

          // compute the speed of sound
          const ZFSFloat speedOfSound = sqrt (m_gamma * pressure(cellId) * Frho);

          // no need for simplified metrics, since information is already contained
          // in cell metrics
          const ZFSFloat lenXi = sqrt( POW2( metric[0] ) +
                                       POW2( metric[1] ) );

          const ZFSFloat lenEt = sqrt( POW2( metric[2] ) +
                                       POW2( metric[3] ) );

          // contravariant velocities
          ZFSFloat U_c = F0;
          ZFSFloat V_c = F0;

          for( ZFSInt isd = xsd; isd < nDim; isd ++ )
            {
              U_c += m_cells->pvariables[ PV->VV[isd] ][ cellId ] * metric[ xsd * nDim + isd ];
              V_c += m_cells->pvariables[ PV->VV[isd] ][ cellId ] * metric[ ysd * nDim + isd ];

            }

          // subtract grid velocity
          U_c -= dxtx[cellId];
          V_c -= dxty[cellId];

          U_c = fabs(U_c);
          V_c = fabs(V_c);

          // has area information in it due to metric terms
          const ZFSFloat eigenvalue = U_c + V_c + speedOfSound * ( lenXi + lenEt );

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

void ZFSStrctrdBlck2D::updateSpongeLayer()
{
  TRACE();
  if(m_useSponge) m_strctrdBndryCnd->updateSpongeLayer();
}

bool ZFSStrctrdBlck2D::rungeKuttaStep()
{
  TRACE();
  const ZFSId noVars = CV->noVariables;
  ZFSId cellId;
  ZFSFloat factor;

#ifdef ZFS_EXTRA_DEBUG
  savePartitions(true); //testing only
#endif

  // set old variables
  if( m_RKStep == 0 ) {
    m_workload += (ZFSFloat)m_noStrctrdCells*m_workloadIncrement;

    if( approx(m_time,F0,m_eps)) {
      m_workloadIncrement = F1 / m_workload;
      m_workload = F1;
    }

    for(ZFSId varId=0; varId< noVars; varId++) {
      for(cellId =0; cellId<m_noStrctrdCells; cellId++) {
        m_cells->oldVariables[varId][cellId]=m_cells->variables[varId][cellId];
      }
    }

    //only save old Jacobian for moving grid
    // if (m_movingGrid) {
    //   for(cellId =0; cellId<m_noStrctrdCells; cellId++) {
    //          this->m_cells->oldCellJac[cellId] = this->m_cells->cellJac[cellId];
    //   }
    // }
  }

  //for moving grid all geometrical variables
  //need to be computed at every RK step
  // if (m_movingGrid) {
  //   computeCellCentreCoordinates();
  //   computeMetrics();
  //   computeJacobian();
  // }

  switch( m_rungeKuttaOrder )
    {
    case 2:
      {
        //for moving grids we take the old Jacobian into account
        // if (m_movingGrid) {
        //   for(ZFSId varId=0; varId<noVars; varId++) {

        //       for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++) {
        //      for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers;i++) {
        //        cellId = cellIndex(i,j;
        //        m_cells->variables[varId][cellId]=(m_cells->oldVariables[varId][cellId]*m_cells->oldCellJac[cellId] + m_RKalpha[m_RKStep]*m_timeStep * m_cells->rightHandSide[varId][cellId])/m_cells->cellJac[cellId];
        //      }
        //     }
        //   }
        // } else {
        if(m_localTimeStep) {
          for(ZFSId varId=0; varId<noVars; varId++) {
            for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++) {
              for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers;i++) {
                cellId = cellIndex(i,j);
                factor=(m_RKalpha[m_RKStep]*m_cells->localTimeStep[cellId]) / m_cells->cellJac[cellId];
                m_cells->variables[varId][cellId]=m_cells->oldVariables[varId][cellId] + factor * m_cells->rightHandSide[varId][cellId];
              }
            }
          }
        } else {
          for(ZFSId varId=0; varId<noVars; varId++) {
            for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++) {
              for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers;i++) {
                cellId = cellIndex(i,j);
                factor=(m_RKalpha[m_RKStep]*m_timeStep) / m_cells->cellJac[cellId];
                m_cells->variables[varId][cellId]=m_cells->oldVariables[varId][cellId] + factor * m_cells->rightHandSide[varId][cellId];
              }
            }
          }
        }
        break;
      }
    case 3:
      {
        for(ZFSId varId=0; varId<noVars; varId++) {
          for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++) {
            for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers;i++){
              cellId = cellIndex(i,j);
              factor=(m_RKalpha[m_RKStep]*m_timeStep) / m_cells->cellJac[cellId];
              m_cells->variables[varId][cellId] = m_RKalpha[ m_RKStep ] * m_cells->variables[varId][cellId] + (F1 - m_RKalpha[ m_RKStep ]) *
                m_cells->oldVariables[varId][cellId] - factor * m_cells->rightHandSide[varId][cellId];
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

void ZFSStrctrdBlck2D::viscousFlux(){
  (this->*viscFluxMethod)();
}

void ZFSStrctrdBlck2D::viscousFluxRANS(){
  m_ransBlck->viscousFluxRANS();
}


void ZFSStrctrdBlck2D::viscousFluxLES()
{
  TRACE();
  const ZFSInt noCells = m_noStrctrdCells;
  const ZFSFloat rPr = F1/m_Pr;
  const ZFSFloat rRe = F1/m_Re0;
  const ZFSFloat gammaMinusOne = m_gamma - 1.0;
  const ZFSFloat FgammaMinusOne = F1 / gammaMinusOne;

  ZFSFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  ZFSFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  ZFSFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  ZFSFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  ZFSFloat* const RESTRICT T = &m_cells->temperature[0];

  ZFSFloat *const RESTRICT eflux= ALIGNED_MF(m_cells->eFlux);
  ZFSFloat *const RESTRICT fflux= ALIGNED_MF(m_cells->fFlux);
  ZFSFloat *const RESTRICT vflux= ALIGNED_MF(m_cells->viscousFlux);

  for(ZFSId j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers+1; j++) {
    for(ZFSId ii=m_noGhostLayers-1; ii<m_nCells[1]-m_noGhostLayers+1; ii++) {
      const ZFSId I=cellIndex(ii,j);
      T[I] = m_gamma*p[I]/rho[I];
    }
  }

  ZFSFloat tau1, tau2, tau4;
  ZFSFloat dTdx, dTdy;

  for(ZFSId j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers+1; j++) {
    for(ZFSId i=m_noGhostLayers-1; i<m_nCells[1]-m_noGhostLayers+1; i++) {
      //get the adjacent cells;
      const ZFSId IJ   = cellIndex(i,j);
      const ZFSId IPJ  = cellIndex((i+1),j);
      const ZFSId IPJP = cellIndex((i+1),(j+1));
      const ZFSId IJP  = cellIndex(i,(j+1));

      const ZFSFloat dudxi=F1B2*(u[IPJP]+u[IPJ]-u[IJP]-u[IJ]);
      const ZFSFloat dudet=F1B2*(u[IPJP]+u[IJP]-u[IPJ]-u[IJ]);

      const ZFSFloat dvdxi=F1B2*(v[IPJP]+v[IPJ]-v[IJP]-v[IJ]);
      const ZFSFloat dvdet=F1B2*(v[IPJP]+v[IJP]-v[IPJ]-v[IJ]);

      const ZFSFloat dTdxi=F1B2*(T[IPJP]+T[IPJ]-T[IJP]-T[IJ]);
      const ZFSFloat dTdet=F1B2*(T[IPJP]+T[IJP]-T[IPJ]-T[IJ]);

      const ZFSFloat uAvg=F1B4*(u[IJP]+u[IPJP]+u[IJ]+u[IPJ]);
      const ZFSFloat vAvg=F1B4*(v[IJP]+v[IPJP]+v[IJ]+v[IPJ]);

      const ZFSFloat mue=F1B4*(zfsSUTHERLANDLAW(T[IJP])+zfsSUTHERLANDLAW(T[IPJP])+zfsSUTHERLANDLAW(T[IJ])+zfsSUTHERLANDLAW(T[IPJ]));

      const ZFSFloat cornerMetrics[9] = {m_cells->cornerMetrics[0*noCells+IJ],
                                         m_cells->cornerMetrics[1*noCells+IJ],
                                         m_cells->cornerMetrics[2*noCells+IJ],
                                         m_cells->cornerMetrics[3*noCells+IJ]};

      // compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

      // tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
      //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
      //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
      tau1 = F4B3 * ( dudxi * cornerMetrics[ xsd * 2 + xsd ]   +
                      dudet * cornerMetrics[ ysd * 2 + xsd ] ) -

        F2B3 * ( dvdxi * cornerMetrics[ xsd * 2 + ysd ]   +
                 dvdet * cornerMetrics[ ysd * 2 + ysd ] );

      // compute tau2 = du/dy + dv/dx

      // tau_xy = du/dxi * dxi/dy + du/deta * deta/dy + du/dzeta * dzeta/dy
      //        + dv/dxi * dxi/dx + dv/deta * deta/dx + dv/dzeta * dzeta/dx
      tau2 = dudxi * cornerMetrics[ xsd * 2 + ysd ] +
        dudet * cornerMetrics[ ysd * 2 + ysd ] +

        dvdxi * cornerMetrics[ xsd * 2 + xsd ] +
        dvdet * cornerMetrics[ ysd * 2 + xsd ];

      // compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

      // tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
      //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
      //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
      tau4 = F4B3 * ( dvdxi * cornerMetrics[ xsd * 2 + ysd ] +
                      dvdet * cornerMetrics[ ysd * 2 + ysd ] ) -

        F2B3 * ( dudxi * cornerMetrics[ xsd * 2 + xsd ] +
                 dudet * cornerMetrics[ ysd * 2 + xsd ]);


      dTdx = dTdxi * cornerMetrics[ xsd * 2 + xsd ] +
        dTdet * cornerMetrics[ ysd * 2 + xsd ];

      dTdy = dTdxi * cornerMetrics[ xsd * 2 + ysd ] +
        dTdet * cornerMetrics[ ysd * 2 + ysd ];

      const ZFSFloat mueOverRe = mue * rRe / m_cells->cornerJac[IJ]; // divide by Jacobian
      tau1*=mueOverRe;
      tau2*=mueOverRe;
      tau4*=mueOverRe;

      const ZFSFloat mueH=FgammaMinusOne*mueOverRe*rPr;
      const ZFSFloat qx=mueH*dTdx+uAvg*tau1+vAvg*tau2;
      const ZFSFloat qy=mueH*dTdy+uAvg*tau2+vAvg*tau4;

      //efluxes
      eflux[ 3 * IJ ]     = tau1 * cornerMetrics[ xsd * 2 + xsd ] +
        tau2 * cornerMetrics[ xsd * 2 + ysd ];

      eflux[ 3 * IJ + 1 ] = tau2 * cornerMetrics[ xsd * 2 + xsd ] +
        tau4 * cornerMetrics[ xsd * 2 + ysd ];

      eflux[ 3 * IJ + 2 ] = qx * cornerMetrics[ xsd * 2 + xsd ] +
        qy * cornerMetrics[ xsd * 2 + ysd ];

      //ffluxes
      fflux[ 3 * IJ ]     = tau1 * cornerMetrics[ ysd * 2 + xsd ] +
        tau2 * cornerMetrics[ ysd * 2 + ysd ];

      fflux[ 3 * IJ + 1 ] = tau2 * cornerMetrics[ ysd * 2 + xsd ] +
        tau4 * cornerMetrics[ ysd * 2 + ysd ];

      fflux[ 3 * IJ + 2 ] = qx * cornerMetrics[ ysd * 2 + xsd ] +
        qy * cornerMetrics[ ysd * 2 + ysd ];

    }
  }

  for(ZFSId var = 0; var < ( CV->noVariables ) - 1; ++var) {
    for(ZFSId j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; ++j) {
      for(ZFSId i=m_noGhostLayers-1; i<m_nCells[1]-m_noGhostLayers; ++i) {
        const ZFSId IJ    = cellIndex(i,j);
        const ZFSId IJM   = cellIndex(i,(j-1));

        vflux[ 2 * IJ ] = F1B2 * ( eflux[ 3 * IJ + var ] + eflux[ 3 * IJM + var ]);

#ifdef ZFS_EXTRA_DEBUG
        viscFluxOut[0][var*m_noStrctrdCells+IJ]= vflux[3*IJ];
#endif
      }
    }

    for(ZFSId j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers; ++j) {
      for(ZFSId i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; ++i) {
        const ZFSId IJ    = cellIndex(i,j);
        const ZFSId IMJ   = cellIndex((i-1),j);

        vflux[ 2 * IJ + 1 ] = F1B2 * ( fflux[ 3 * IJ + var ] + fflux[ 3 * IMJ + var ]);

#ifdef ZFS_EXTRA_DEBUG
        viscFluxOut[1][var*m_noStrctrdCells+IJ]= vflux[3*IJ+1];
#endif
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
}



void ZFSStrctrdBlck2D::loadRestartBC2600() {
  if(m_bc2600IsActive && !m_bc2600InitialStartup) {
    if(domainId()==0) {
      cout << "Loading BC2600 values..." << endl;
    }

    ZFSInt bcCells[2] = {m_nInputBlockCells[0],m_noGhostLayers};
    ZFSInt noCellsBC = bcCells[0]*bcCells[1];
    ZFSInt bcOffset[2] = {0,0};
    ZFSFloatScratchSpace tmpRestartVars(noCellsBC*CV->noVariables, __CALLING_FUNCTION__, "m_tmpRestartVars2600");

    if(domainId()==0) {
      stringstream restartFileName;
      ZFSInt restartFileId =-1;
      ZFSString restartFile = *(ZFSContext::getProperty("restartVariablesFileName", m_blockId, __CALLING_FUNCTION__, (ZFSString*) NULL )->asString(0));
      restartFileName << outputDir() << restartFile;

      restartFileId = io_openfile("hdf5" ,(restartFileName.str()).c_str(),"collective", MPI_COMM_SELF);
      stringstream pathStr;
      pathStr << "/block" << m_inputBlockId << "/bc2600" << endl;
      const char* path = (pathStr.str()).c_str();

      for(ZFSId var=0; var<CV->noVariables; var++) {
        cout << "Loading " << m_variableNames[var] << " offset: " << var*noCellsBC << endl;
        io_read_ddataset_part1d1(restartFileId, path, m_variableNames[var].c_str(), nDim, bcOffset, bcCells, &tmpRestartVars[var*noCellsBC]);
      }

      io_closefile(restartFileId);
    }

    MPI_Bcast(&tmpRestartVars[0], noCellsBC*CV->noVariables, MPI_DOUBLE, 0, m_zfsStrctrdComm);

    if(domainId()==0) {
      cout << "Loading BC2600 values... SUCCESSFUL!" << endl;
    }

    if(m_bc2600) {
      ZFSId startGC[2] = {0,0};
      ZFSId endGC[2] = {0,0};

      if(m_bc2600noOffsetCells[0] == 0) { startGC[0] = m_noGhostLayers; }
      if(m_bc2600noOffsetCells[0]+m_bc2600noActiveCells[0] == bcCells[0]) { endGC[0] = m_noGhostLayers; }

      for(ZFSInt i = 0; i< m_noGhostLayers; i++) {
        for(ZFSInt j = startGC[0]; j<m_bc2600noCells[0]-endGC[0]; j++) {
          ZFSId cellId = cellIndex(i,j);
          ZFSId globalI = i;
          ZFSId globalJ = m_bc2600noOffsetCells[0]-m_noGhostLayers+j;
          ZFSId cellIdBC = globalI+globalJ*bcCells[1];

          //load values from restart field
          for(ZFSId var=0; var<CV->noVariables; var++) {
            m_cells->variables[var][cellId] = tmpRestartVars[var*noCellsBC+cellIdBC];
          }
        }
      }


      //Fix diagonal cells at end of domain
      if(m_bc2600noOffsetCells[0] + m_bc2600noActiveCells[0] == m_nInputBlockCells[0]) {
        for(ZFSInt i = 0; i<m_noGhostLayers; i++) {
          const ZFSId cellIdA2 = cellIndex(i,m_noGhostLayers+m_bc2600noActiveCells[0]-2);
          const ZFSId cellIdA1 = cellIndex(i,m_noGhostLayers+m_bc2600noActiveCells[0]-1);
          const ZFSId cellIdG1 = cellIndex(i,m_noGhostLayers+m_bc2600noActiveCells[0]);
          for(ZFSId var=0; var<CV->noVariables; var++) {
            const ZFSFloat distA1A2 = sqrt(POW2(m_cells->coordinates[0][cellIdA1]-m_cells->coordinates[0][cellIdA2])+
                                           POW2(m_cells->coordinates[2][cellIdA1]-m_cells->coordinates[1][cellIdA2]));
            const ZFSFloat slope = (m_cells->variables[var][cellIdA1]-m_cells->variables[var][cellIdA2])/distA1A2;
            const ZFSFloat distG1A1 = sqrt(POW2(m_cells->coordinates[0][cellIdG1]-m_cells->coordinates[0][cellIdA1])+
                                           POW2(m_cells->coordinates[2][cellIdG1]-m_cells->coordinates[1][cellIdA1]));
            m_cells->variables[var][cellIdG1] = m_cells->variables[var][cellIdA1] + distG1A1*slope;
          }
        }
      }
    }
  }
}

void ZFSStrctrdBlck2D::computePrimitiveVariables() {
  const ZFSFloat gammaMinusOne = m_gamma - 1.0;

  ZFSFloat** const RESTRICT cvars = m_cells->variables;
  ZFSFloat** const RESTRICT pvars = m_cells->pvariables;

  for(ZFSId j=m_noGhostLayers; j < m_nCells[0]-m_noGhostLayers; ++j) {
    for(ZFSId i=m_noGhostLayers; i < m_nCells[1]-m_noGhostLayers; ++i) {
      const ZFSId cellId = cellIndex(i,j);
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


/** brief parallel: coordinates the communication (exchange)
 *
 *
 */
void ZFSStrctrdBlck2D::exchange()
{
  TRACE();
  RECORD_TIMER_START(m_tcomm);
  RECORD_TIMER_START(m_texchange);

  if(noDomains()>1) {
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
  }

  RECORD_TIMER_STOP(m_texchange);
  RECORD_TIMER_STOP(m_tcomm);
}


void ZFSStrctrdBlck2D::gather(){
  TRACE();
  ZFSId cellId;
  for(ZFSInt nghbr=0; nghbr<m_cmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt* startInfo=m_cmnctnFlag->startInfoSNDcells[nghbr];
    ZFSInt* endInfo= m_cmnctnFlag->endInfoSNDcells[nghbr];
    ZFSFloat* bufferSnd = m_cmnctnFlag->m_bufferCellsSnd[nghbr];
    ZFSInt pos=0;

    for(ZFSId var=0; var<PV->noVariables; var++) {
      for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
        for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
          cellId = i +(j*m_nCells[1]);
          bufferSnd[pos]=m_cells->pvariables[var][cellId];
          pos++;
        }
      }
    }
  }
}

void ZFSStrctrdBlck2D::send(){
  TRACE();
  for(ZFSId nghbr=0; nghbr<m_cmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt tag = domainId()+(m_cmnctnFlag->m_tagHelperSND[nghbr])*noDomains();
    ZFSInt err = MPI_Isend((void*)&m_cmnctnFlag->m_bufferCellsSnd[nghbr][0], m_cmnctnFlag->m_noNghbrDomainCellBufferSizeSnd[nghbr], MPI_DOUBLE, m_cmnctnFlag->m_sndNghbrId[nghbr], tag, m_zfsStrctrdComm, &m_cmnctnFlag->mpi_sndRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck2D::receive(){
  TRACE();
  for(ZFSId nghbr=0; nghbr<m_cmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt tag = m_cmnctnFlag->m_rcvNghbrId[nghbr]+(m_cmnctnFlag->m_tagHelperRCV[nghbr])*noDomains();
    ZFSInt err = MPI_Irecv((void*)&m_cmnctnFlag->m_bufferCellsRcv[nghbr][0],m_cmnctnFlag->m_noNghbrDomainCellBufferSizeRcv[nghbr], MPI_DOUBLE,m_cmnctnFlag->m_rcvNghbrId[nghbr] ,tag, m_zfsStrctrdComm, &m_cmnctnFlag->mpi_rcvRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck2D::scatter(){
  TRACE();

  ZFSId cellId;
  //the ordering of the grid points can be different from
  //sending instance ==> reorder it and copy it to the
  //right place

  for(ZFSInt nghbr=0; nghbr< m_cmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt j2, i2, id2;
    ZFSInt* step2 = m_cmnctnFlag->stepInfoRCV[nghbr];
    ZFSInt* order = m_cmnctnFlag->orderInfo[nghbr];
    ZFSInt start2[2];
    ZFSInt end2[2];
    ZFSInt len2[2];
    ZFSInt totalCells=1;
    ZFSInt len1[2];

    for(ZFSInt j=0; j<nDim; j++) {
      len1[j]=m_cmnctnFlag->endInfoRCVcells[nghbr][j] - m_cmnctnFlag->startInfoRCVcells[nghbr][j];
      if(len1[j]!=0){
        totalCells*=len1[j];
      }
    }

    for(ZFSInt j=0; j<nDim; j++) {
      start2[j]=0;
      end2[j]=len1[order[j]]-1;
      len2[j]=len1[order[j]];
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
      j2=start2[1];
      for(ZFSInt j=startInfo[1]; j<endInfo[1]; j++) {
        i2=start2[0];
        for(ZFSInt i=startInfo[0]; i<endInfo[0]; i++) {
          id2=var*totalCells+i2+j2*len2[0];
          cellId = i +j*m_nCells[1];
          m_cells->pvariables[var][cellId]= bufferRcv[id2];

          i2+=step2[0];
          pos++;
        }
        j2+=step2[1];
      }
    }
  }
}

void ZFSStrctrdBlck2D::gatherPoints()
{
  TRACE();
  ZFSId pointId;
  for(ZFSInt nghbr=0; nghbr< m_cmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt* startInfo=m_cmnctnFlag->startInfoSNDpoints[nghbr];
    ZFSInt* endInfo= m_cmnctnFlag->endInfoSNDpoints[nghbr];
    ZFSFloat* bufferSnd = m_cmnctnFlag->m_bufferPointsSnd[nghbr];
    ZFSInt pos=0;
    for(ZFSId dim=0; dim<nDim; dim++) {
      for(ZFSInt j=startInfo[1]; j<endInfo[1]+1; j++) {
        for(ZFSInt i=startInfo[0]; i<endInfo[0]+1; i++) {
          pointId = i +(j*m_nPoints[1]);
          bufferSnd[pos]=m_coordinates[dim][pointId];
          pos++;
        }
      }
    }
  }
}

void ZFSStrctrdBlck2D::sendPoints()
{
  TRACE();
  for(ZFSId nghbr=0; nghbr<m_cmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt tag = domainId()+(m_cmnctnFlag->m_tagHelperSND[nghbr])*noDomains();
    ZFSInt err = MPI_Isend((void*)&m_cmnctnFlag->m_bufferPointsSnd[nghbr][0], m_cmnctnFlag->m_noNghbrDomainPointBufferSizeSnd[nghbr], MPI_DOUBLE, m_cmnctnFlag->m_sndNghbrId[nghbr], tag, m_zfsStrctrdComm, &m_cmnctnFlag->mpi_sndRequest[nghbr]);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void ZFSStrctrdBlck2D::receivePoints()
{
  TRACE();
  for(ZFSId nghbr=0; nghbr<m_cmnctnFlag->noNghbrDomainsNormal; nghbr++) {
    ZFSInt tag = m_cmnctnFlag->m_rcvNghbrId[nghbr]+(m_cmnctnFlag->m_tagHelperRCV[nghbr])*noDomains();

    MPI_Irecv((void*)&m_cmnctnFlag->m_bufferPointsRcv[nghbr][0],m_cmnctnFlag->m_noNghbrDomainPointBufferSizeRcv[nghbr], MPI_DOUBLE,m_cmnctnFlag->m_rcvNghbrId[nghbr] ,tag, m_zfsStrctrdComm, &m_cmnctnFlag->mpi_rcvRequest[nghbr]);
  }
}

void ZFSStrctrdBlck2D::scatterPoints()
{
  TRACE();
  ZFSId pointId;
  //the ordering of the grid points can be different from
  //sending instance ==> reorder it and copy it to the
  //right place
  for(ZFSInt nghbr=0; nghbr< m_cmnctnFlag->noNghbrDomainsNormal; nghbr++) {
      ZFSInt j2, i2, id2;
      ZFSInt* step2 = m_cmnctnFlag->stepInfoRCV[nghbr];
      ZFSInt* order = m_cmnctnFlag->orderInfo[nghbr];
      ZFSInt start2[2];
      ZFSInt end2[2];
      ZFSInt len2[2];
      ZFSInt totalPoints=1;
      ZFSInt len1[2];
      for(ZFSInt j=0; j<nDim; j++) {
        len1[j]=m_cmnctnFlag->endInfoRCVpoints[nghbr][j] - m_cmnctnFlag->startInfoRCVpoints[nghbr][j]+1;
        totalPoints*=len1[j];
      }

      for(ZFSInt j=0; j<nDim; j++) {
        start2[j]=0;
        end2[j]=len1[order[j]]-1;
        len2[j]=len1[order[j]];
        if(step2[j]<0) {
          ZFSInt dummy=start2[j];
          start2[j]=end2[j];
          end2[j]=dummy;
        }
      }

      ZFSInt* startInfo=m_cmnctnFlag->startInfoRCVpoints[nghbr];
      ZFSInt* endInfo= m_cmnctnFlag->endInfoRCVpoints[nghbr];

      ZFSFloat* bufferRcv = m_cmnctnFlag->m_bufferPointsRcv[nghbr];
      for(ZFSId dim=0; dim<nDim; dim++) {
        j2=start2[1];
        for(ZFSInt j=startInfo[1]; j<endInfo[1]+1; j++) {
          i2=start2[0];
          for(ZFSInt i=startInfo[0]; i<endInfo[0]+1; i++) {
            id2=dim*totalPoints+i2+(j2*len2[0]);
            pointId = i +(j*m_nPoints[1]);
            m_coordinates[dim][pointId]= bufferRcv[id2];
            i2+=step2[0];
          }
          j2+=step2[1];
        }
      }
  }
}


inline ZFSFloat ZFSStrctrdBlck2D::getPSI(ZFSId I, ZFSId dim) {
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
