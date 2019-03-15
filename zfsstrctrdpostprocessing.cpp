#include "zfsstrctrdpostprocessing.h"

#include "zfsglobals.h"
#include "zfsstrctrdcell.h"
#include "zfsstrctrdblck.h"
#include "zfsiolib.h"

using namespace std;

/** \brief Constructor for the postprocessing block
 *
 * \author Andreas Lintermann
 * \date 12.09.2012
 *
 * Reads in options for postprocessing by calling initProcessingBlock().
 *
 *
 **/
template <ZFSInt nDim, class Block>
ZFSStrctrdPostprocessing<nDim, Block>::ZFSStrctrdPostprocessing()
  : m_postprocessing(0) {
  TRACE();
}

/** \brief Destructor for the postprocessing block
 *
 * \author Andreas Lintermann
 * \date 26.08.2012
 *
 * \tparam[in] T celltype
 *
 **/
template <ZFSInt nDim, class Block>
ZFSStrctrdPostprocessing<nDim, Block>::~ZFSStrctrdPostprocessing()
{
  TRACE();

  if(m_postprocessing)
  {
    if(m_noPostprocessingOps>0)
      for (ZFSInt op = 0; op < m_noPostprocessingOps; op++)
      {
        switch(string2enum(m_postprocessingOps[op]))
        {
        case ZFSPP_AVERAGE_PRE:
        case ZFSPP_AVERAGE_IN:
        case ZFSPP_AVERAGE_POST:
        case ZFSPP_TAUW_PRE:
        case ZFSPP_LOAD_AVERAGED_SOLUTION_PRE:
        {
          zfsDeallocate( m_summedVars );
          zfsDeallocate( m_square );
          if( m_useKahan ){
            zfsDeallocate( m_cSum );
            zfsDeallocate( m_ySum );
            zfsDeallocate( m_tSum );

            zfsDeallocate( m_cSquare );
            zfsDeallocate( m_ySquare );
            zfsDeallocate( m_tSquare );

            if( m_kurtosis ){
              zfsDeallocate( m_cCube );
              zfsDeallocate( m_yCube );
              zfsDeallocate( m_tCube );

              zfsDeallocate( m_cFourth );
              zfsDeallocate( m_yFourth );
              zfsDeallocate( m_tFourth );
            }
            else if( m_skewness ){
              zfsDeallocate( m_cCube );
              zfsDeallocate( m_yCube );
              zfsDeallocate( m_tCube );
            }
          }
          if( m_kurtosis ){
            zfsDeallocate( m_cube );
            zfsDeallocate( m_fourth );
          }
          else if( m_skewness ){
            zfsDeallocate( m_cube );
          }
          break;
        }
        case ZFSPP_COMPUTE_PRODUCTION_TERMS_PRE:
        {
          zfsDeallocate( m_production );
          break;
        }
        case ZFSPP_MOVING_AVERAGE_IN:
        case ZFSPP_MOVING_AVERAGE_PRE:
        case ZFSPP_MOVING_AVERAGE_POST:
          {
            break;
          }
        default:
        {
          zfsTerm(1,__CALLING_FUNCTION__,"Unknown postprocessing operation");
        }
        }
      }
    delete [] m_postprocessingOps;
  }
}

template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::initStrctrdPostprocessing(){
  TRACE();


  m_movingGrid = 0;
  m_movingGrid = *(ZFSContext::getProperty( "movingGrid", ppblock()->blockId(), __CALLING_FUNCTION__, &m_movingGrid))->asInt(0);
  //m_movingGrid = ppblock()->m_movingGrid; //testen


  /*! \page propertyPage1
    \section postprocessing
    <code>ZFSInt ZFSStrctrdPostprocessing::postprocessing</code>\n
    default = <code>0</code>\n\n
    This property determines the postrprocessing.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_postprocessing = 0;
  m_postprocessing = *(ZFSContext::getProperty( "postprocessing", ppblock()->blockId(), __CALLING_FUNCTION__, &m_postprocessing))->asInt(0);

  /*! \page propertyPage1
    \section pp_skewness
    <code>ZFSInt ZFSPostprocesssingBlock::m_skewness</code>\n
    default = <code>0</code>\n\n
    This propertpy determines if skewness is computed when ZFSPP_AVERAGE_PRE/IN/POST is activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */

  m_averageVorticity = 0;
  m_averageVorticity = *(ZFSContext::getProperty( "pp_averageVorticity", ppblock()->blockId(), __CALLING_FUNCTION__, &m_averageVorticity))->asInt(0);

  m_skewness = 0;
  m_skewness = *(ZFSContext::getProperty( "pp_skewness", ppblock()->blockId(), __CALLING_FUNCTION__, &m_skewness))->asInt(0);

  /*! \page propertyPage1
    \section pp_kurtosis
    <code>ZFSInt ZFSPostprocesssingBlock::m_kurtosis</code>\n
    default = <code>0</code>\n\n
    This property determines if kurtosis (and skewness) is computed when ZFSPP_AVERAGE_PRE/IN/POST is activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_kurtosis = 0;
  m_kurtosis = *(ZFSContext::getProperty( "pp_kurtosis", ppblock()->blockId(), __CALLING_FUNCTION__, &m_kurtosis))->asInt(0);
  if( m_kurtosis ) m_skewness = 1;

  /*! \page propertyPage1
    \section pp_turbulentProduction
    <code>ZFSInt ZFSPostprocesssingBlock::m_computeProductionTerms</code>\n
    default = <code>0</code>\n\n
    Determines if the turbulent production terms should be computed after the normal averaging
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_computeProductionTerms = 0;
  m_computeProductionTerms = *(ZFSContext::getProperty( "pp_turbulentProduction", ppblock()->blockId(), __CALLING_FUNCTION__, &m_computeProductionTerms))->asInt(0);

  //  if(m_kurtosis==1 || m_skewness==1)
  //    zfsTerm(1,__CALLING_FUNCTION__,"check zfsfvblock and all computations because pressure ampl was introduced at position 11");

  /*! \page propertyPage1
    \section pp_twoPass
    <code>ZFSInt ZFSPostprocesssingBlock::m_twoPass</code>\n
    default = <code>0</code>\n\n
    This property determines if two-pass averaging is performed in ZFSPP_AVERAGE_PRE/POST.\n
    Either m_twoPass or m_useKahan should be activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_twoPass = 0;
  m_twoPass = *(ZFSContext::getProperty( "pp_twoPass", ppblock()->blockId(), __CALLING_FUNCTION__, &m_twoPass))->asInt(0);

  /*! \page propertyPage1
    \section pp_useKahan
    <code>ZFSInt ZFSPostprocesssingBlock::m_useKahan</code>\n
    default = <code>0</code>\n\n
    This property determines if kahan summation is performed in ZFSPP_AVERAGE_PRE/IN/POST. \n
    Either m_twoPass or m_useKahan should be activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_useKahan = 0;
  m_useKahan = *(ZFSContext::getProperty( "pp_useKahan", ppblock()->blockId(), __CALLING_FUNCTION__, &m_useKahan))->asInt(0);

  /*! \page propertyPage1
    \section pp_fileName
    <code>ZFSInt ZFSPostprocesssingBlock::m_postprocessFileName</code>\n
    default = <code>""</code>\n\n
    This property determines a filename for averaging.\n
    <ul>
    <li><code>filename</code></li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */


  m_postprocessFileName = "";
  m_postprocessFileName = *(ZFSContext::getProperty( "pp_fileName", ppblock()->blockId(), FUN_, &m_postprocessFileName))->asString(0);

  //Init vars
  m_averageInterval = 0;
  m_averageRestart = 0;
  m_averageRestartInterval = 0;
  m_averageStartTimestep = 0;
  m_averageStopTimestep = 0;
  m_noPostprocessingOps = 0;
  m_noSamples = 0;

  m_noVariables = ppblock()->noVariables();

  if(m_postprocessing)
    {
      for(ZFSInt i=0; i<3; i++)
        {
          tvecpost tmp;
          vector <ZFSString> st;
          m_postprocessingMethods.push_back(tmp);
          m_postprocessingMethodsDesc.push_back(st);
        }

      /*! \page propertyPage1
        \section postprocessingOps
        <code>ZFSString* ZFSBlock::m_postprocessingOps</code>\n
        default = <code>empty</code>\n\n
        This property is a list of postprocessing operations to be performed
        <ul>
        <li><code>ZFSPP_AVERAGE_PRE</code> </li>
        <li><code>ZFSPP_AVERAGE_POST</code> </li>
        <li><code>ZFSPP_AVERAGE_IN</code> </li>
        <li><code>ZFSPP_MOVING_AVERAGE_PRE</code> </li>
        <li><code>ZFSPP_MOVING_AVERAGE_POST</code> </li>
        <li><code>ZFSPP_MOVING_AVERAGE_IN</code> </li>
        </ul>\n
        Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
      */
      m_noPostprocessingOps = ZFSContext::getProperty("postprocessingOps", ppblock()->blockId(), AT_, (ZFSString*) NULL, 1 )->count();
      m_postprocessingOps = NULL;
      if(m_noPostprocessingOps>0)
        {
          //zfsAlloc( m_postprocessingOps, m_noPostprocessingOps, "m_postprocessingOps",  __CALLING_FUNCTION__ );
          m_postprocessingOps = new ZFSString[m_noPostprocessingOps];
          for (ZFSInt op = 0; op < m_noPostprocessingOps; op++)
            {
              m_postprocessingOps[op] = *(ZFSContext::getProperty("postprocessingOps", ppblock()->blockId(), __CALLING_FUNCTION__, &m_postprocessingOps[op] ,m_noPostprocessingOps )->asString(op));

              switch(string2enum(m_postprocessingOps[op]))
                {
                case ZFSPP_AVERAGE_IN:
                  {
                    m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
                    m_postprocessingMethods[1].push_back(&ZFSStrctrdPostprocessing::averageSolutionsInSolve);
                    initAverageIn();
                    initAverageVariables();
                    break;
                  }
                case ZFSPP_AVERAGE_PRE:
                  {
                    m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
                    m_postprocessingMethods[0].push_back(&ZFSStrctrdPostprocessing::averageSolutions);
                    initTimeStepProperties();
                    initAverageVariables();
                    break;
                  }
                case ZFSPP_AVERAGE_POST:
                  {
                    m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
                    m_postprocessingMethods[2].push_back(&ZFSStrctrdPostprocessing::averageSolutions);
                    initTimeStepProperties();
                    initAverageVariables();
                    break;
                  }
                case ZFSPP_LOAD_AVERAGED_SOLUTION_PRE:
                  {
                    m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
                    m_postprocessingMethods[0].push_back(&ZFSStrctrdPostprocessing::loadAveragedSolution);
                    initTimeStepProperties();
                    initAverageVariables();
                    break;
                  }
                case ZFSPP_COMPUTE_PRODUCTION_TERMS_PRE:
                  {
                    m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
                    m_postprocessingMethods[0].push_back(&ZFSStrctrdPostprocessing::computeProductionTerms);
                    initProductionVariables();
                    break;
                  }
                case ZFSPP_TAUW_PRE:
                  {
                    m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
                    m_postprocessingMethods[0].push_back(&ZFSStrctrdPostprocessing::computeAverageSkinFriction);
                    initTimeStepProperties();
                    initAverageVariables();
                    break;
                  }
                case ZFSPP_MOVING_AVERAGE_IN:
                  {
                    m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
                    m_postprocessingMethods[1].push_back(&ZFSStrctrdPostprocessing::movingAverage);
                    initTimeStepProperties();
                    initMovingAverage();
                    break;
                  }
                case ZFSPP_MOVING_AVERAGE_PRE:
                  {
                    m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
                    m_postprocessingMethods[0].push_back(&ZFSStrctrdPostprocessing::movingAveragePost);
                    initTimeStepProperties();
                    initMovingAverage();
                    break;
                  }
                case ZFSPP_MOVING_AVERAGE_POST:
                  {
                    m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
                    m_postprocessingMethods[2].push_back(&ZFSStrctrdPostprocessing::movingAveragePost);
                    initTimeStepProperties();
                    initMovingAverage();
                    break;
                  }
                default:
                  {
                    zfsTerm(1,__CALLING_FUNCTION__,"Unknown postprocessing operation");
                  }
                }
            }
        }
    }
}

/**
 *
 * @author Frederik Temme, Dez 22, 2015
 * modified 22.12.2015
 */
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::postprocessPreInit(){
  TRACE();

  if(m_averageRestart) {
    m_restartTimeStep = ppblock()->m_restartTimeStep;
  }

  if(m_noPostprocessingOps!=0 && m_postprocessing == 1)
    {
      for (ZFSInt op = 0; op < m_noPostprocessingOps; op++)
        {

          switch(string2enum(m_postprocessingOps[op]))
            {
            case ZFSPP_AVERAGE_POST:
            case ZFSPP_AVERAGE_PRE:
            case ZFSPP_AVERAGE_IN:
              {
                // load the averaging restart
                if(m_restartTimeStep > m_averageStartTimestep && m_restartTimeStep <= m_averageStopTimestep)
                  {
                    ZFSString name = ppblock()->outputDir() + "PostprocessingRestart_";
                    ZFSChar buf1[10];
                    ZFSChar buf2[10];
                    sprintf( buf1, "%d", m_averageStartTimestep );
                    sprintf( buf2, "%d", m_restartTimeStep );
                    name.append(buf1);
                    name += "-";
                    name.append(buf2);
                    name.append(ppblock()->m_outputFormat);

                    zfs_log << "\n\n"
                            << "    ^^^^^^^^^^^^^^^ Entering postprocessing mode PreInit ^^^^^^^^^^^^^^^^ \n"
                            << "    ^   - Loading restart for mean flow calculation: " << name << "\n"
                            << "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n\n";

                    ppblock()->loadAverageRestartFile( name.c_str(), m_summedVars, m_square, m_cube, m_fourth );
                  }
                break;
              }
            case ZFSPP_TAUW_PRE:
              {
                cout << "Loading postprocessing file " << m_postprocessFileName << endl;
                if(m_averageRestart) {
                  ppblock()->loadAverageRestartFile( m_postprocessFileName.c_str(), m_summedVars, m_square, m_cube, m_fourth );
                  computeAveragedSolution();
                } else {
                  ppblock()->loadAveragedVariables(m_postprocessFileName.c_str());
                }
                cout << "Loading postprocessing file " << m_postprocessFileName << " ... SUCCESSFULL!" << endl;

                break;
              }
            default:
              {
                // has to be filled
              }
            }
        }
    }
}

template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::postprocessPreSolve()
{
  TRACE();

  if(m_noPostprocessingOps!=0 && m_postprocessing == 1)
    {
      zfs_log << "\n\n"
              << "    ^^^^^^^^^^^^^^^ Entering postprocessing mode PreSolve ^^^^^^^^^^^^^^^ \n"
              << "    ^   - Activated operations are:\n";
      for (ZFSInt op = 0; op < m_noPostprocessingOps; op++)
        zfs_log << "    ^      + " << m_postprocessingOps[op] << "\n";
      zfs_log << "    ^   - Running:\n";

      for (ZFSInt op = 0; op < (signed)m_postprocessingMethods[0].size(); op++)
        {
          zfs_log << "    ^      + " << m_postprocessingMethodsDesc[0][op] << "\n";
          (this->*(m_postprocessingMethods[0][op]))();
        }
      zfs_log << "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n" << endl;
    }
}

template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::postprocessInSolve()
{
  TRACE();

  if(m_noPostprocessingOps!=0 && m_postprocessing == 1)
    {
      for (ZFSInt op = 0; op < (signed)m_postprocessingMethods[1].size(); op++)
        {
          (this->*(m_postprocessingMethods[1][op]))();
        }
    }
}

template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::postprocessPostSolve()
{
  TRACE();

  if(m_noPostprocessingOps!=0 && m_postprocessing == 1)
    {

      zfs_log << "\n\n"
              << "    ^^^^^^^^^^^^^^^ Entering postprocessing mode PostSolve ^^^^^^^^^^^^^^^ \n"
              << "    ^   - Activated operations are:\n\n";
      for (ZFSInt op = 0; op < m_noPostprocessingOps; op++)
        zfs_log << "    ^      + " << m_postprocessingOps[op] << "\n";
      zfs_log << "    ^   - Running:\n";

      for (ZFSInt op = 0; op < (signed)m_postprocessingMethods[2].size(); op++)
        {
          zfs_log << "    ^      + " << m_postprocessingMethodsDesc[2][op] << "\n";
          (this->*(m_postprocessingMethods[2][op]))();
        }
      zfs_log << "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n" << endl;
    }
}


/**
 *
 * @author Frederik Temme, Nov 15, 2015
 * modified 17.12.2015
 *
 * modified 5.10.2016 by Thomas Luerkens
 *
 */
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::averageSolutionsInSolve(){
  TRACE();
  //Average only at the right timestep
  if (globalTimeStep >= m_averageStartTimestep && globalTimeStep <= m_averageStopTimestep) {
    if(ppblock()->getGridMovingMethod()!=9) {
      if((globalTimeStep - m_averageStartTimestep) % m_averageInterval == 0) {
        addAveragingSample();
      }
    } else {
      if((globalTimeStep - m_averageStartTimestep) % m_averageInterval == 0) {
        m_noSamples++;

        ppblock()->spanwiseWaveReorder();
        addTempWaveSample();
        if(ppblock()->domainId()==0) {cout << ">>>>> wave sample with interval  " << m_averageInterval << " time steps at time step: " << globalTimeStep << " has been added  <<<<<" << endl;}
      }
    }
  }

  //Compute the averaged solution and write to file
  if(globalTimeStep == m_averageStopTimestep){
    saveAveragedSolution(globalTimeStep);
  }
}

/** \brief Adds for the travelling wave setups
 *
 * This method is similar to addAveragingSample but
 * was adapted for the phase-averaging of the travelling
 * wave setups
 *
 * \author Marian Albers (original by Ansgar Niemoeller)
 * \date 15.12.2016
 *
 *
 **/
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::addTempWaveSample(){

  TRACE();
  const ZFSId noCells = ppblock()->getNoCells();
  const ZFSId noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);

  for(ZFSInt cellId=0; cellId<noCells; cellId++) {
    ZFSInt offset = 0;
    ZFSInt offsetSquare = 0;
    // Primitive variables
    for(ZFSId varId = 0; varId < m_noVariables; varId++) {
      m_summedVars[varId+offset][cellId] += m_tempWaveSample[varId][cellId];
    }
    offset += m_noVariables;

    // Vorticities
    if(m_averageVorticity) {
      for(ZFSId varId = 0; varId < noAveragedVorticities; varId++) {
        m_summedVars[varId+offset][cellId] += m_tempWaveSample[varId+offset][cellId];
      }
      offset += noAveragedVorticities;
    }

    // squares of velocity components ( uu, vv, ww )
    for(ZFSId varId = 0; varId < nDim; varId++) {
      m_square[varId+offsetSquare][cellId] += m_tempWaveSample[varId][cellId]*m_tempWaveSample[varId][cellId];
    }
    offsetSquare += nDim;

    //product of different velocity componets ( uv, vw, wu )
    for(ZFSId varId = 0; varId < 2*nDim-3; varId++) {
      m_square[varId+offsetSquare][cellId] += m_tempWaveSample[varId%nDim][cellId]*m_tempWaveSample[(varId+1)%nDim][cellId];
    }
    offsetSquare += 2*nDim-3;

    //square of pressure (pp)
    m_square[offsetSquare][cellId] += m_tempWaveSample[nDim+1][cellId]*m_tempWaveSample[nDim+1][cellId];
    offsetSquare += 1;

    //squares of the vorticity
    if(m_averageVorticity) {
      for( ZFSId varId=0; varId < noAveragedVorticities; varId++) {
        m_square[offsetSquare+varId][cellId] += m_tempWaveSample[varId+m_noVariables][cellId] * m_tempWaveSample[varId+m_noVariables][cellId];
      }
      offsetSquare += noAveragedVorticities;
    }

    //third and fouth powers of velocity components (skewness and kurtosis)
    if(m_kurtosis) {
      for(ZFSId varId =0; varId < nDim; varId++) {
        m_cube[varId][cellId] += pow( m_tempWaveSample[varId][cellId], 3);
        m_fourth[varId][cellId] += pow( m_tempWaveSample[varId][cellId], 4);
      }
    }
    //third power if velocity components (skewness)
    if(m_skewness) {
      for(ZFSId varId =0; varId < nDim; varId++) {
        m_cube[varId][cellId] += pow( m_tempWaveSample[varId][cellId], 3);
      }
    }
  }
}






/**
 *
 * @author Frederik Temme, Jan 14, 2016
 * modified 19.1.2016
 *
 */
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::averageSolutions(){
  TRACE();

  /***********************************************************
  in properties: solutionInterval has to be equal to pp_averageInterval
                 if using this method.
  ***********************************************************/

  zfs_log <<"    ^          * Averaging solutions ";


  //set the summationstart for averaging
  ZFSFloat avgStart = m_restartTimeStep;

  for(ZFSInt avgTimestep=avgStart; avgTimestep<=m_averageStopTimestep; avgTimestep+=m_averageInterval){

    // load samples
    stringstream filename;
    filename << ppblock()->outputDir() << avgTimestep << ppblock()->m_outputFormat;
    ppblock()->loadSampleFile(filename.str());
    ppblock()->exchange();
    ppblock()->applyBoundaryCondition();
    addAveragingSample();

    //Write average restart file
    if (m_averageRestartInterval!=0 && (avgTimestep>=m_averageStartTimestep && avgTimestep % m_averageRestartInterval == 0 && avgTimestep<=m_averageStopTimestep)){
      saveAverageRestart();
    }
  }

  saveAveragedSolution(m_averageStopTimestep);
}

template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block >::saveAveragedSolution(ZFSInt endTimeStep) {
  TRACE();

  computeAveragedSolution();

  //output filename
  ZFSString name = ppblock()->outputDir() + "Mean_";
  ZFSChar buf1[10];
  ZFSChar buf2[10];
  sprintf( buf1, "%d", m_averageStartTimestep );
  sprintf( buf2, "%d", endTimeStep );
  name.append(buf1);
  name += "-";
  name.append(buf2);

  zfs_log << "         ^   saving averaged variables " << name << endl;

  ppblock()->saveAveragedVariables(name, getNoPPVars(), m_summedVars);
}


/** \brief Computes the mean variables from summed vars
 *
 * Computes the correct averaged solution from all added samples
 * Also computes the RMS components of the velocities, the pressure
 * and the vorticities (if desired)
 *
 * \author Marian Albers (original by Ansgar Niemoeller)
 * \date 01.02.2016
 *
 *
 **/
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::computeAveragedSolution() {
  TRACE();

  const ZFSId noCells = ppblock()->getNoCells();
  const ZFSId noAveragedVorticities
      = (m_averageVorticity != 0) * (nDim * 2 - 3);

  const ZFSFloat weight = F1/m_noSamples;//F1/(((m_averageStopTimestep-m_averageStartTimestep)/m_averageInterval)+1);
  ZFSInt offset =0;
  ZFSInt offsetSquare =0;

  //mean of summed primitive variables
  for( ZFSId cellId=0; cellId < noCells; cellId++){
    for( ZFSId varId=0; varId < m_noVariables; varId++){
      m_summedVars[varId+offset][cellId] *= weight;
    }
  }
  offset+=m_noVariables;

  //Weighting of summed vorticity variables -> mean
  if (m_averageVorticity) {
    for (ZFSId cellId = 0; cellId < noCells; cellId++) {
      for (ZFSId varId = 0; varId < noAveragedVorticities; varId++) {
        m_summedVars[varId + offset][cellId] *= weight;
      }
    }

    offset += noAveragedVorticities;
  }

  //compute mean(u'u'),mean(v'v'),mean(w'w') ( e.g.mean(u'u')=mean(u^2)-(u_mean))^2 )
  for( ZFSId cellId=0; cellId < noCells; cellId++){
    for( ZFSId varId=0; varId < nDim; varId++){
      m_summedVars[varId+offset][cellId] = weight*m_square[varId][cellId] - m_summedVars[varId][cellId]*m_summedVars[varId][cellId];
    }
  }
  offset+=nDim;
  offsetSquare+=nDim;


  //compute mean(u'v'),mean(v'w'),mean(w'u') ( e.g. mean(u'v')=mean(u*v)-u_mean*v_mean )
  for( ZFSId cellId=0; cellId < noCells; cellId++){
    for( ZFSId varId=0; varId < 2*nDim-3; varId++){
      m_summedVars[varId+offset][cellId] = weight*m_square[varId+offsetSquare][cellId] - m_summedVars[varId%nDim][cellId]*m_summedVars[(varId+1)%nDim][cellId];

    }
  }
  offset+= 2*nDim-3;
  offsetSquare+= 2*nDim-3;

  if( m_kurtosis ){
    //compute skewness and kurtosis of velocity components
    //e.g. skewness(u) = mean(u^3) - 3*u_mean*mean(u^2) + 2*u_mean^3
    //e.g. kurtosis(u) = mean(u^4) - 4*u_mean*mean(u^3) + 6*u_mean^2*mean(u^2) - 3*u_mean^4
    for( ZFSId cellId=0; cellId < noCells; cellId++){
      for( ZFSId varId=0; varId < nDim; varId++){
        m_summedVars[varId+offset][cellId] = weight*m_cube[varId][cellId] - 3*weight*m_summedVars[varId][cellId]*m_square[varId][cellId] + 2*pow(m_summedVars[varId][cellId],3);

        m_summedVars[varId+offset+nDim][cellId] = weight*m_fourth[varId][cellId] - 4*weight*m_cube[varId][cellId]*m_summedVars[varId][cellId] + 6*weight*m_square[varId][cellId]*m_summedVars[varId][cellId]*m_summedVars[varId][cellId] - 3*pow(m_summedVars[varId][cellId],4);

      }
    }
    offset+= 2*nDim;

  } else if( m_skewness ){
    //compute skewness of velocity components
    for( ZFSId cellId=0; cellId < noCells; cellId++ ){
      for( ZFSId varId=0; varId < nDim; varId++ ){
        m_summedVars[varId+offset][cellId] = weight*m_cube[varId][cellId] - 3*weight*m_summedVars[varId][cellId]*m_square[varId][cellId] + 2*pow(m_summedVars[varId][cellId],3);
      }
    }
    offset+= nDim;
  }

  // compute p'*p'
  for( ZFSId cellId=0; cellId < noCells; cellId++ ){
    m_summedVars[offset][cellId] = weight*m_square[offsetSquare][cellId] - m_summedVars[nDim+1][cellId]*m_summedVars[nDim+1][cellId]; // pressure
  }
  offset+= 1;
  offsetSquare+=1;

  if(m_averageVorticity) {
    //compute vorticity symmetric rms
    for( ZFSId cellId=0; cellId < noCells; cellId++){
      for( ZFSId varId=0; varId < nDim; varId++){
        m_summedVars[varId+offset][cellId] = weight*m_square[varId+offsetSquare][cellId] - m_summedVars[varId+m_noVariables][cellId]*m_summedVars[varId+m_noVariables][cellId];
      }
    }

    offset += noAveragedVorticities;
    offsetSquare += noAveragedVorticities;
  }


  // add aditional variables here -> otherwise the code for kurtosis will overwrite the summed vars of the new introduced variables

}

/** \brief Adds one sample to the summedVars
 *
 * Adds a sample of the current time step to the
 * summedVars (and other) fields
 *
 * \author Marian Albers (original by Ansgar Niemoeller)
 * \date 12.12.2016
 *
 *
 **/
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::addAveragingSample() {
  TRACE();

  zfs_log << "     ^        * Averaging timestep " << globalTimeStep << "\n";
  m_noSamples++;

  const ZFSId noCells = ppblock()->getNoCells();
  ZFSFloatScratchSpace cellVars(m_noVariables, __CALLING_FUNCTION__, "cellVars");

  const ZFSId noAveragedVorticities
      = (m_averageVorticity != 0) * (nDim * 2 - 3);

  if (m_averageVorticity) {
    ppblock()->computeVorticity();
  }

  for(ZFSInt cellId =0; cellId <noCells; cellId++){

    /* List of Variables in m_summedVars after following computation
       0=mean(u)      1=mean(v)      2=mean(w)
       3=mean(rho)    4=mean(p)
       5=mean(u'u')   6=mean(v'v')  7=mean(w'w')
       8=mean(u'v')   9=mean(v'w') 10=mean(w'u')
       11=skew u      12=skew v     13=skew w
       14=kurt u      15=kurt v     16=kurt w
    */

    // Calculation of primitive variables
    ppblock()->getSampleVariables(cellId, cellVars.begin());

    if(m_useKahan){// Kahan summation

      /* Kahan summation pseudocode
         sum=0; c=0;
         for i=0 to input.length
         y = input[i] -c;
         t = sum + y;
         c = (t-sum) - y;
         sum = t;
      */

      ZFSInt offsetSquare =0;
      for( ZFSId varId=0; varId < m_noVariables; varId++ ){ //sum up all primitive variables
        m_ySum[varId][cellId] = cellVars[varId] - m_cSum[varId][cellId];
        m_tSum[varId][cellId] = m_summedVars[varId][cellId] + m_ySum[varId][cellId];
        m_cSum[varId][cellId] = ( m_tSum[varId][cellId] - m_summedVars[varId][cellId] ) - m_ySum[varId][cellId];
        m_summedVars[varId][cellId] = m_tSum[varId][cellId];
      }
      for( ZFSId varId=0; varId < nDim; varId++){ //squares of velocity components (u*u,v*v,w*w)
        m_ySquare[varId][cellId] = ( cellVars[varId]*cellVars[varId] ) - m_cSquare[varId][cellId];
        m_tSquare[varId][cellId] = m_square[varId][cellId] + m_ySquare[varId][cellId];
        m_cSquare[varId][cellId] = ( m_tSquare[varId][cellId] - m_square[varId][cellId] ) - m_ySquare[varId][cellId];
        m_square[varId][cellId] = m_tSquare[varId][cellId];
      }
      offsetSquare+=3;
      for( ZFSId varId=0; varId < nDim; varId++){ //products of different velocity components in order (u*v,v*w,w*)
        m_ySquare[varId+offsetSquare][cellId] = (cellVars[varId%3]*cellVars[(varId+1)%3]) - m_cSquare[varId+offsetSquare][cellId];
        m_tSquare[varId+offsetSquare][cellId] = m_square[varId+offsetSquare][cellId] + m_ySquare[varId+offsetSquare][cellId];
        m_cSquare[varId+offsetSquare][cellId] = (m_tSquare[varId+offsetSquare][cellId] - m_square[varId+offsetSquare][cellId] ) - m_ySquare[varId+offsetSquare][cellId];
        m_square[varId+offsetSquare][cellId] = m_tSquare[varId+offsetSquare][cellId];
      }
      if( m_kurtosis ){ //compute third and fourth power of velocity components for skewness and kurtosis
        for( ZFSId varId=0; varId < nDim; varId++ ){
          m_yCube[varId][cellId] = pow(cellVars[varId],3) - m_cCube[varId][cellId];
          m_tCube[varId][cellId] = m_cube[varId][cellId] + m_yCube[varId][cellId];
          m_cCube[varId][cellId] = ( m_tCube[varId][cellId] - m_cube[varId][cellId] ) - m_yCube[varId][cellId];
          m_cube[varId][cellId] = m_tCube[varId][cellId];

          m_yFourth[varId][cellId] = pow(cellVars[varId],4) - m_cFourth[varId][cellId];
          m_tFourth[varId][cellId] = m_fourth[varId][cellId] + m_yFourth[varId][cellId];
          m_cFourth[varId][cellId] = ( m_tFourth[varId][cellId] - m_fourth[varId][cellId] ) - m_yFourth[varId][cellId];
          m_fourth[varId][cellId] = m_tFourth[varId][cellId];
        }
      }
      else if( m_skewness ){ //compute only third power of velocity components for skewness
        for( ZFSId varId=0; varId < nDim; varId++ ){
          m_yCube[varId][cellId] = pow(cellVars[varId],3) - m_cCube[varId][cellId];
          m_tCube[varId][cellId] = m_cube[varId][cellId] + m_yCube[varId][cellId];
          m_cCube[varId][cellId] = ( m_tCube[varId][cellId] - m_cube[varId][cellId] ) - m_yCube[varId][cellId];
          m_cube[varId][cellId] = m_tCube[varId][cellId];
        }
      }

    } else { //normal summation
      // Reset offsets
      ZFSInt offset       = 0;
      ZFSInt offsetSquare = 0;

      // Primitive variables
      for (ZFSId varId = 0; varId < m_noVariables; varId++) {
        m_summedVars[varId + offset][cellId] += cellVars[varId];
      }
      offset += m_noVariables;

      // Vorticities
      if (m_averageVorticity) {
        for( ZFSId varId=0; varId < noAveragedVorticities; varId++ ){
          m_summedVars[varId + offset][cellId]
            += ppblock()->getSampleVorticity(cellId,varId);
        }
        offset += noAveragedVorticities;
      }

      for( ZFSId varId=0; varId < nDim; varId++){ //squares of velocity components (u*u,v*v(,w*w))
        m_square[varId][cellId] += cellVars[ varId ] * cellVars[varId];
      }
      offsetSquare+=nDim;
      for( ZFSId varId=0; varId < 2*nDim-3; varId++){ //products of different velocity components (u*v(,v*w,w*u))
        m_square[offsetSquare+varId][cellId] += (cellVars[varId%nDim])*(cellVars[(varId+1)%nDim]);
      }
      offsetSquare+=2*nDim-3;
      m_square[offsetSquare][cellId] += cellVars[ nDim+1 ] * cellVars[ nDim+1 ]; // squares of pressure  p*p
      offsetSquare+=1;

      // add aditional variables here -> otherwise the code for kurtosis will overwrite the summed vars of the new introduced variables
      //squares of the vorticity
      if(m_averageVorticity) {
        for( ZFSId varId=0; varId < noAveragedVorticities; varId++) {
          m_square[offsetSquare+varId][cellId] += ppblock()->getSampleVorticity(cellId,varId) * ppblock()->getSampleVorticity(cellId,varId);
        }
        offsetSquare += noAveragedVorticities;
      }


      if( m_kurtosis ){ //third and fourth powers of velocity components (skewness and kurtosis)
        for( ZFSId varId=0; varId < nDim; varId++ ){
          m_cube[varId][cellId] += pow( cellVars[ varId ], 3 );
          m_fourth[varId][cellId] += pow( cellVars[ varId ], 4 );
        }
      }
      else if( m_skewness ){ //third powers of velocity components (skewness)
        for( ZFSId varId=0; varId < nDim; varId++ ){
          m_cube[varId][cellId] += pow( cellVars[ varId ], 3 );
        }
      }
    }
  }
}

/** \brief Loads an averaged file again
 *
 * Loads an Postprocessing Mean file into the corresponding field
 * such that further pp actions can be performed
 *
 * Specify file name in 'pp_fileName' in the property file and
 * set 'pp_averageRestart' if this is an PostprocessingRestart and not
 * a Mean file
 *
 * \author Marian Albers
 * \date 01.02.2016
 *
 *
 **/
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::loadAveragedSolution() {
  cout << "Loading postprocessing file " << m_postprocessFileName << endl;
  if(m_averageRestart) {
    ppblock()->loadAverageRestartFile( m_postprocessFileName.c_str(), m_summedVars, m_square, m_cube, m_fourth );
    computeAveragedSolution();
  } else {
    ppblock()->loadAveragedVariables(m_postprocessFileName.c_str());
  }

  cout << "Filling ghost-cells..." << endl;
  ppblock()->ppFillGhostCells();
  cout << "Filling ghost-cells... FINISHED!" << endl;

  const ZFSId noCells = ppblock()->getNoCells();
  for(ZFSId cellId=0; cellId<noCells; cellId++) {
    for(ZFSId var=0; var<m_noVariables; var++) {
      ppblock()->saveVarToPrimitive(cellId, var, m_summedVars[var][cellId]);
    }
  }

  cout << "Computing conservative variables" << endl;
  ppblock()->computeConservativeVariables();

  if(ppblock()->isMovingGrid()){
    cout << "Moving grid to correct position!" << endl;
    ppblock()->moveGrid(false,true);
  }
}

/** \brief Computes skin friction of an averaged field
 *
 * \author Marian Albers
 * \date 12.12.2016
 *
 *
 **/
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::computeAverageSkinFriction() {
  cout << "Saving averaged variables to primitive variables" << endl;
  const ZFSId noCells = ppblock()->getNoCells();
  for(ZFSId cellId=0; cellId<noCells; cellId++) {
    for(ZFSId var=0; var<m_noVariables; var++) {
      ppblock()->saveVarToPrimitive(cellId, var, m_summedVars[var][cellId]);
    }
  }

  cout << "Computing conservative variables" << endl;
  ppblock()->computeConservativeVariables();

  if(ppblock()->isMovingGrid()){
    cout << "Moving grid to correct position!" << endl;
    ppblock()->moveGrid(false,true);
  }
  ppblock()->saveAuxData();
  zfsTerm(1,__CALLING_FUNCTION__, "My work here is done!");
}

/** \brief Computes the production terms from an averaged field 
 *
 * \author Marian Albers
 * \date 10.01.2017
 *
 *
 **/
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::computeProductionTerms() {
  TRACE();

  const ZFSId noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);
  const ZFSId offset = m_noVariables + noAveragedVorticities;

  ZFSFloat* ubar = &m_summedVars[0][0];
  ZFSFloat* vbar = &m_summedVars[1][0];
  ZFSFloat* wbar = &m_summedVars[2][0];

  ZFSFloat* uu = &m_summedVars[offset+0][0];
  ZFSFloat* vv = &m_summedVars[offset+1][0];
  ZFSFloat* ww = &m_summedVars[offset+2][0];
  ZFSFloat* uv = &m_summedVars[offset+3][0];
  ZFSFloat* vw = &m_summedVars[offset+4][0];
  ZFSFloat* uw = &m_summedVars[offset+5][0];

  for(ZFSInt cellId =0; cellId < ppblock()->getNoCells(); cellId++){
    m_production[0][cellId] = uu[cellId]*ppblock()->dvardxyz(cellId,0,ubar) +
                              uv[cellId]*ppblock()->dvardxyz(cellId,1,ubar) +
                              uw[cellId]*ppblock()->dvardxyz(cellId,2,ubar);
    m_production[1][cellId] = uv[cellId]*ppblock()->dvardxyz(cellId,0,vbar) +
                              vv[cellId]*ppblock()->dvardxyz(cellId,1,vbar) +
                              vw[cellId]*ppblock()->dvardxyz(cellId,2,vbar);
    m_production[2][cellId] = uw[cellId]*ppblock()->dvardxyz(cellId,0,wbar) +
                              vw[cellId]*ppblock()->dvardxyz(cellId,1,wbar) +
                              ww[cellId]*ppblock()->dvardxyz(cellId,2,wbar);
  }

  ppblock()->saveProductionTerms(m_postprocessFileName.c_str(), m_production);

  
  zfsTerm(1,__CALLING_FUNCTION__, "My work here is done!");
}

template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::movingAverage(){
  TRACE();
}

template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::movingAveragePost(){
  TRACE();
}

/**
 * \fn void ZFSStrctrdPostprocessing<nDim, Block>::initAverageIn()
 * \brief Initializes properties for averaging during solver run
 *
 * \author Andreas Lintermann (last modified Ansgar Niemoeller 07/14)
 * \date 14.09.2012
 *
 * \tparam[in] Block blocktype
 *
 **/
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::initAverageIn()
{
  TRACE();

  initTimeStepProperties();

  /*! \page propertyPage1
    \section pp_averageRestart
    <code>ZFSString* ZFSBlock::m_averageRestart</code>\n
    default = <code>0</code>\n\n
    This property determines if we should restart from our last averaging.
    <ul>
    <li><code>0</code> turned off </li>
    <li><code>1</code> turned on</li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageRestart = *(ZFSContext::getProperty( "pp_averageRestart", ppblock()->blockId(), __CALLING_FUNCTION__, &m_averageRestart))->asInt(0);

  /*! \page propertyPage1
    \section pp_averageRestartInterval
    <code>ZFSString* ZFSBlock::m_averageRestartInterval</code>\n
    default = <code>0</code>\n\n
    This property determines the interval to write averaging restart files. Has to be a multiple of averageInterval and of restartInterval.
    <ul>
    <li><code>interval</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageRestartInterval = *(ZFSContext::getProperty( "pp_averageRestartInterval", ppblock()->blockId(), __CALLING_FUNCTION__, &m_averageRestartInterval))->asInt(0);

  if(m_averageRestartInterval % ppblock()->restartInterval() != 0) {
      zfsTerm(1,__CALLING_FUNCTION__, "The property 'averageRestartInterval' has to be a multiple of the property 'restartInterval'...");
    }
}

/** \brief allocates memory for averageSolutions() and averageSolutionsInSolve()
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 * \param[in] noInternalCells the number of internal cells on which averaging is performed
 *
 **/
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::initAverageVariables()
{
  TRACE();
  const ZFSId noCells = ppblock()->getNoCells();
  const ZFSId noVars  = getNoPPVars();
  const ZFSId noSquareVars = getNoPPSquareVars();

  zfsAlloc( m_summedVars, noVars, noCells , "m_summedVars", F0, FUN_ ); // +1 for pressure ampl
  zfsAlloc( m_square, noSquareVars, noCells, "m_square", F0, FUN_ ); // + 1 for pressure ampl

  if( m_kurtosis )
  {
    zfs_log << "Allocating cube and fourth field for kurtosis computation for " << noCells << " cells" << endl;
    zfsAlloc( m_cube, nDim, noCells, "m_cube", F0, FUN_ );
    zfsAlloc( m_fourth, nDim, noCells, "m_fourth", F0, FUN_ );
  }
  else if( m_skewness /*&& !m_twoPass*/ )
  {
    zfsAlloc( m_cube, nDim, noCells, "m_cube", F0, FUN_ );
  }

  if( m_useKahan ) // allocate memory for kahan summation
  {
    zfs_log << "m_useKahan is activated" << endl;
    zfsAlloc( m_cSum, m_noVariables+3*(nDim-1), noCells, "m_cSum", F0, FUN_ );
    zfsAlloc( m_tSum, m_noVariables+3*(nDim-1), noCells, "m_tSum", F0, FUN_ );
    zfsAlloc( m_ySum, m_noVariables+3*(nDim-1), noCells, "m_ySum", F0, FUN_ );
    zfsAlloc( m_cSquare, 3*(nDim-1), noCells, "m_cSquare", F0, FUN_ );
    zfsAlloc( m_tSquare, 3*(nDim-1), noCells, "m_tSquare", F0, FUN_ );
    zfsAlloc( m_ySquare, 3*(nDim-1), noCells, "m_ySquare", F0, FUN_ );
    if( m_kurtosis )
    {
      zfsAlloc( m_cCube,  nDim, noCells, "m_cCube", F0, FUN_ );
      zfsAlloc( m_tCube, nDim, noCells, "m_tCube", F0, FUN_ );
      zfsAlloc( m_yCube, nDim, noCells, "m_yCube", F0, FUN_ );
      zfsAlloc( m_cFourth, nDim, noCells, "m_cFourth", F0, FUN_ );
      zfsAlloc( m_tFourth, nDim, noCells, "m_tFourth", F0, FUN_ );
      zfsAlloc( m_yFourth, nDim, noCells, "m_yFourth", F0, FUN_ );
    }
    else if( m_skewness ){
      zfsAlloc( m_cCube, nDim, noCells, "m_cCube", F0, FUN_ );
      zfsAlloc( m_tCube, nDim, noCells, "m_tCube", F0, FUN_ );
      zfsAlloc( m_yCube, nDim, noCells, "m_yCube", F0, FUN_ );
    }
  }

  zfsAlloc(m_avgVariableNames, noVars, "m_avgVariableNames", __CALLING_FUNCTION__);

  //Mean values
  m_avgVariableNames[0]="um";
  m_avgVariableNames[1]="vm";
  if(nDim == 3) {
    m_avgVariableNames[2]="wm";
  }
  m_avgVariableNames[nDim]="rhom";
  m_avgVariableNames[nDim+1]="pm";
  ZFSInt offset = m_noVariables;

  // Vorticities
  if (m_averageVorticity) {
    if (nDim == 3) {
      m_avgVariableNames[offset + 0] = "vortxm";
      m_avgVariableNames[offset + 1] = "vortym";
      m_avgVariableNames[offset + 2] = "vortzm";
      offset += 3;
    } else {
      m_avgVariableNames[offset + 0] = "vortzm";
      offset += 1;
    }
  }

  //reynolds stress components
  m_avgVariableNames[offset+0]="uu";
  m_avgVariableNames[offset+1]="vv";
  offset += 2;
  if( nDim==3 ){
    m_avgVariableNames[offset + 0] = "ww";
    m_avgVariableNames[offset + 1] = "uv";
    m_avgVariableNames[offset + 2] = "vw";
    m_avgVariableNames[offset + 3] = "uw";
    offset += 4;
  } else {
    m_avgVariableNames[offset + 0] = "uv";
    offset += 1;
  }

  // Skewness variables
  if(noVars > m_noVariables+3*(nDim-1) + 1  && m_skewness ){
    if(nDim==3){
      m_avgVariableNames[offset+0]="uuu";
      m_avgVariableNames[offset+1]="vvv";
      m_avgVariableNames[offset+2]="www";
      offset+=3;
    }
    else {
      m_avgVariableNames[offset+0]="uuu";
      m_avgVariableNames[offset+1]="vvv";
      offset+=2;
    }
  }

  // Kurtosis variables
  if( (noVars > m_noVariables+3*(nDim-1)+nDim + 1) && m_kurtosis){
    if(nDim==3){
      m_avgVariableNames[offset+0]="uuuu";
      m_avgVariableNames[offset+1]="vvvv";
      m_avgVariableNames[offset+2]="wwww";
      offset+=3;
    }
    else {
      m_avgVariableNames[offset+0]="uuuu";
      m_avgVariableNames[offset+1]="vvvv";
      offset+=2;
    }
  }

  //pressure fluctuation
  m_avgVariableNames[offset+0] = "pp";
  offset+=1;

  //rms of the vorticities
  if (m_averageVorticity) {
    if(nDim==3) {
      m_avgVariableNames[offset + 0] = "vortrmsx";
      m_avgVariableNames[offset + 1] = "vortrmsy";
      m_avgVariableNames[offset + 2] = "vortrmsz";
      offset += 3;
    } else {
      m_avgVariableNames[offset + 0] = "vortrmsz";
      offset += 1;
    }
  }
}

template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::initProductionVariables()
{
  const ZFSId noCells = ppblock()->getNoCells();
  zfsAlloc( m_production, nDim, noCells, "m_production", F0, FUN_);
}

/**
 * \fn void ZFSStrctrdPostprocessing<nDim, Block>::initTimeStepProperties()
 * \brief Initializes timestep properties for postprocessing
 *
 * \author Andreas Lintermann (last modified Ansgar Niemoeller, 07/14)
 * \date 14.09.2012
 * 
 * reads properties pp_averageStartTimestep, pp_averageStopTimestep and pp_averageInterval
 *
 * \tparam[in] Block blocktype
 *
 **/
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::initTimeStepProperties()
{
  TRACE();

  /*! \page propertyPage1
    \section pp_averageInterval
    <code>ZFSString* ZFSBlock::m_averageInterval</code>\n
    default = <code>0</code>\n\n
    This property determines the interval of the solutions used for averaging.
    <ul>
    <li><code>interval</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageInterval = *(ZFSContext::getProperty( "pp_averageInterval", ppblock()->blockId(), __CALLING_FUNCTION__, &m_averageInterval))->asInt(0);


  if(m_averageInterval == 0)
    {
      zfsTerm(1,__CALLING_FUNCTION__, "Please specify the property 'averageInterval' ...");
    }

  /*! \page propertyPage1
    \section pp_averageStartTimestep
    <code>ZFSString* ZFSBlock::m_averageStartTimestep</code>\n
    default = <code>0</code>\n\n
    This property determines the start timestep used for averaging.
    <ul>
    <li><code>timestep</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageStartTimestep = *(ZFSContext::getProperty( "pp_averageStartTimestep", ppblock()->blockId(), __CALLING_FUNCTION__, &m_averageStartTimestep))->asInt(0);

  if(m_averageStartTimestep == 0)
    {
      zfsTerm(1,__CALLING_FUNCTION__, "Please specify the property 'averageStartTimestep' ...");
    }

  /*! \page propertyPage1
    \section pp_averageStopTimestep
    <code>ZFSString* ZFSBlock::m_averageStopTimestep</code>\n
    default = <code>0</code>\n\n
    This property determines the stop timestep used for averaging.
    <ul>
    <li><code>timestep</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageStopTimestep = *(ZFSContext::getProperty( "pp_averageStopTimestep", ppblock()->blockId(), __CALLING_FUNCTION__, &m_averageStopTimestep))->asInt(0);

  if(m_averageStopTimestep == 0)
    {
      zfsTerm(1,__CALLING_FUNCTION__, "Please specify the property 'averageStopTimestep' ...");
    }

  /*! \page propertyPage1
    \section pp_averageInterval
    <code>ZFSString* ZFSBlock::m_averageRestart</code>\n
    default = <code>0</code>\n\n
    This property determines if we should restart from our last averaging.
    <ul>
    <li><code>0</code> turned off </li>
    <li><code>1</code> turned on</li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageRestart = *(ZFSContext::getProperty( "pp_averageRestart", ppblock()->blockId(), __CALLING_FUNCTION__, &m_averageRestart))->asInt(0);


  /*! \page propertyPage1
    \section pp_averageRestartInterval
    <code>ZFSString* ZFSBlock::m_averageRestartInterval</code>\n
    default = <code>0</code>\n\n
    This property determines the interval to write averaging restart files. Has to be a multiple of averageInterval and of restartInterval.
    <ul>
    <li><code>interval</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageRestartInterval = *(ZFSContext::getProperty( "pp_averageRestartInterval", ppblock()->blockId(), __CALLING_FUNCTION__, &m_averageRestartInterval))->asInt(0);
  if(m_averageRestartInterval % m_averageInterval != 0)
    {   }

/**
 * \fn void ZFSStrctrdPostprocessing<nDim, Block>::initMovingAverage()
 * \brief initializes properties and allocates memory for moving averaging
 *
 * \author Ansgar Niemoeller
 * \date 09.07.2014
 *
 * reads properties pp_movingAvgInterval, pp_movingAvgDataPoints and pp_averageVorticity
 *
 * \param[in] grid pointer to the grid
 *
 **/
}
template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::initMovingAverage()
{
  TRACE();

  const ZFSId noCells = ppblock()->getNoCells();

  /*! \page propertyPage1
    \section pp_movingAvgInterval
    <code>ZFSInt ZFSPostprocesssingBlock::m_movingAvgInterval</code>\n
    default = <code>1</code>\n\n
    This property determines the interval between timesteps considered for moving average,
    e.g. if set to 2 at an averaging timestep n (see pp_averagStartTimestep, pp_averageStopTimestep, pp_averageInterval)
    the timesteps n, n-2, n-4, ... are used to compute the moving average\n
    Note: this has to be a factor of m_averageInterval\n
    see also pp_movingAvgDataPoints, pp_averageVorticity
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_movingAvgInterval = 1;
  m_movingAvgInterval = *(ZFSContext::getProperty( "pp_movingAvgInterval", ppblock()->blockId(), __CALLING_FUNCTION__, &m_movingAvgInterval))->asInt(0);
  if( m_movingAvgInterval<1 ){
    TERMM( 1, "m_movingAvgInterval has to be >=1");
  }
  if( m_averageInterval%m_movingAvgInterval!=0 || m_movingAvgInterval>m_averageInterval ){
    TERMM( 1, "m_movingAvgInterval has to be a factor of m_averageInterval");
  }

  /*! \page propertyPage1
    \section pp_movingAvgDataPoints
    <code>ZFSInt ZFSPostprocesssingBlock::m_movingAvgDataPoints</code>\n
    This property determines the number of timesteps (data points) used for moving average computation\n
    see also pp_movingAvgInterval, pp_averageVorticity
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_movingAvgDataPoints = *(ZFSContext::getProperty( "pp_movingAvgDataPoints", ppblock()->blockId(), __CALLING_FUNCTION__, &m_movingAvgDataPoints))->asInt(0);
  if( m_movingAvgDataPoints < 2 ){
    TERMM( 1, "m_movingAvgDataPoints has to be at least 2");
  }

  /*! \page propertyPage1
    \section pp_averageVorticity
    <code>ZFSInt ZFSPostprocesssingBlock::m_pp.m_averageVorticity</code>\n
    default = <code>0</code>\n\n
    This property determines if the vorticity vector is considered in computation of averages
    <ul>c
    <li><code>0</code> disabled</li>
    <li><code>1</code> enabled</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_averageVorticity = 0;
  m_averageVorticity = *(ZFSContext::getProperty( "pp_averageVorticity", ppblock()->blockId(), __CALLING_FUNCTION__, &m_averageVorticity))->asInt(0);

  m_movingAvgCounter = 0;

  m_movAvgNoVariables = m_noVariables;
  if( m_averageVorticity==1 ){
    m_movAvgNoVariables += ( 2*nDim - 3 );
  }
  zfsAlloc( m_movAvgVariables, noCells, m_movAvgNoVariables*m_movingAvgDataPoints, "m_movAvgVariables", F0, FUN_ );
  zfsAlloc( m_movAvgVarNames, 2*m_movAvgNoVariables, "m_movAvgVarNames", ZFSString("default"), FUN_ );
}


template <ZFSInt nDim, class Block>
void ZFSStrctrdPostprocessing<nDim, Block>::saveAverageRestart(){
  TRACE();

  if (m_postprocessing && (globalTimeStep>=m_averageStartTimestep && globalTimeStep<=m_averageStopTimestep)){
    ZFSString name = ppblock()->outputDir() + "PostprocessingRestart_";
    ZFSChar buf1[10];
    ZFSChar buf2[10];
    sprintf( buf1, "%d", m_averageStartTimestep );
    sprintf( buf2, "%d", globalTimeStep );
    name.append(buf1);
    name += "-";
    name.append(buf2);
    zfs_log << "       ^     saving average restart " << name << endl;

    ppblock()->saveAverageRestart(name, getNoPPVars(), m_summedVars, m_square, m_cube, m_fourth);
  }
}

/** \brief Returns number of postprocessing variables
 *
 * \author Marian Albers
 * \date 01.08.2016
 *
 *
 **/
template <ZFSInt nDim, class Block>
ZFSInt ZFSStrctrdPostprocessing<nDim, Block>::getNoPPVars() {
  TRACE();
  ZFSId c = 0;
  if(m_kurtosis) c=2;
  else if(m_skewness) c=1;


  // Determine number of averaged variables
  // Mean vorticities and symmetric rms components (6 for 3D, 4 for 2D)
  const ZFSId noAveragedVorticities
    = (m_averageVorticity != 0) * (nDim * 2 - 3);
  const ZFSId noVars = m_noVariables            // primitive variables
                       + noAveragedVorticities  // mean vorticities
                       + 3 * (nDim - 1)         // Reynolds stress components
                       + nDim * c               // skewness/kurtosis
                       + 1                      // pressure amplitude p'
                       + noAveragedVorticities; // rms vorticities

  return noVars;
}

/** \brief Returns number of pp Square variables
 *
 * \author Marian Albers
 * \date 01.08.2016
 *
 *
 **/
template <ZFSInt nDim, class Block>
ZFSInt ZFSStrctrdPostprocessing<nDim, Block>::getNoPPSquareVars() {
  TRACE();

  // Determine number of averaged variables
  const ZFSId noAveragedVorticities
      = (m_averageVorticity != 0) * (nDim * 2 - 3);
  const ZFSId noVars = 3*(nDim-1)               // uu,vv,ww,uv,uw,vw
                       + 1                      // pressure amplitude p'
                       + noAveragedVorticities; // vorticity rms
  return noVars;
}

// Explicit instantiations for 2D and 3D
template class ZFSStrctrdPostprocessing<2, ZFSStrctrdBlck<2>>;
template class ZFSStrctrdPostprocessing<3, ZFSStrctrdBlck<3>>;
