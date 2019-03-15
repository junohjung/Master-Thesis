#ifndef ZFSSTRCTRDPOSTPROCESSING_H
#define ZFSSTRCTRDPOSTPROCESSING_H

#include <vector>
#include <map>
#include "zfstypes.h"

#include "zfsconstants.h"
#include "zfsfunctions.h"

// Forard declarations
class ZFSStrctrdCell;

template <ZFSInt nDim, class Block>
class ZFSStrctrdPostprocessing
{
  typedef ZFSStrctrdCell PPCell;

public:
  ZFSStrctrdPostprocessing();
  ~ZFSStrctrdPostprocessing();
  
  void postprocessPreInit();
  void postprocessPreSolve();
  void postprocessPostSolve();
  void postprocessInSolve();

protected:

  void initStrctrdPostprocessing();
  void initAverageIn();
  void initAverageVariables();
  void initTimeStepProperties();
  void initMovingAverage();
  void initProductionVariables();
  void averageSolutionsInSolve();
  void averageSolutions();
  void addAveragingSample();
  void addTempWaveSample();
  void saveAveragedSolution(ZFSInt);
  void computeAveragedSolution();
  void computeAverageSkinFriction();
  void movingAverage();
  void movingAveragePost();
  void computeProductionTerms();
  void loadAveragedSolution();

  void saveAverageRestart();

  void loadMeanFile( const char* fileName );
  void getSampleVariables( ZFSId cellId, const ZFSFloat*& vars );
  ZFSInt getNoPPVars();
  ZFSInt getNoPPSquareVars();

  ZFSInt m_postprocessing;
  ZFSInt m_noPostprocessingOps;
  ZFSString* m_postprocessingOps;

  // this vector holds all the functions to be called
  // it is initialized with three entries (each position of postprocessing)
  typedef void(ZFSStrctrdPostprocessing::* tpost)();
  typedef std::vector < tpost > tvecpost;
  std::vector < tvecpost > m_postprocessingMethods;
  std::vector <std::vector <ZFSString> >  m_postprocessingMethodsDesc;

public:
  ZFSInt m_restartTimeStep;

protected:
  // Number of variables of the corresponding block
  ZFSId m_noVariables;

  //Averaging variables
  ZFSFloat** m_summedVars;
  ZFSFloat** m_square;
  ZFSFloat** m_cube;
  ZFSFloat** m_fourth;
  ZFSFloat** m_tempWaveSample;

  //Kahan summation
  ZFSInt m_useKahan;
  ZFSFloat** m_cSum;
  ZFSFloat** m_ySum;
  ZFSFloat** m_tSum;
  ZFSFloat** m_cSquare;
  ZFSFloat** m_ySquare;
  ZFSFloat** m_tSquare;
  ZFSFloat** m_cCube;
  ZFSFloat** m_yCube;
  ZFSFloat** m_tCube;
  ZFSFloat** m_cFourth;
  ZFSFloat** m_yFourth;
  ZFSFloat** m_tFourth;

  ZFSId m_twoPass;
  ZFSId m_skewness;
  ZFSId m_kurtosis;
  ZFSId m_averageVorticity = 0;
  ZFSFloat** m_production;
  ZFSString* m_avgVariableNames;

  ZFSId m_movingAvgInterval;
  ZFSInt m_movingAvgDataPoints;
  ZFSInt m_movingAvgCounter;
  ZFSInt m_movAvgNoVariables;
  ZFSFloat** m_movAvgVariables;
  ZFSString* m_movAvgVarNames;

  //averaging
  ZFSInt m_averageInterval;
  ZFSInt m_averageStartTimestep;
  ZFSInt m_averageStopTimestep;
  ZFSInt m_averageRestartInterval;
  ZFSInt m_averageRestart;
  ZFSInt m_noSamples;

  // moving Grid
  ZFSId m_movingGrid;


  static const ZFSInt xsd = 0;
  static const ZFSInt ysd = 1;
  static const ZFSInt zsd = 2;

  ZFSInt m_computeProductionTerms;

  ZFSString m_postprocessFileName;

 private:
  ZFSString m_blockType;
  Block* ppblock() { return static_cast<Block*>(this); } ///< CRTP
};
#endif
