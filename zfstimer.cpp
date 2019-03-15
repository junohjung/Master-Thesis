#include "zfstimer.h"

#include <algorithm>

using namespace std;
//--------------------------------------------------------------------------------

ZFSTimers& timers() {
  static ZFSTimers timers;
  return timers;
}

vector<ZFSFunctionTiming> ZFSProfile::s_functionTimings;


ZFSBool operator<(const ZFSFunctionTiming& a, const ZFSFunctionTiming& b) {
  if ( a.getDeltaCpuTime() == b.getDeltaCpuTime() ) return ( a.getInitCpuTime() < b.getInitCpuTime() );
  return ( a.getDeltaCpuTime() < b.getDeltaCpuTime() );
}

ZFSProfile::~ZFSProfile() {
  const clock_t exitCpuTime = cpuTime();
  const ZFSFloat exitWallTime = wallTime();
  const ZFSFloat thresholdPercentage = 0.5;
  stringstream sstream;
  sstream << "    CPU      WALL   FUNCTION                    >> profile: '" << m_name <<"' <<";
  const string header = sstream.str();
  for ( std::size_t i = 0; i < header.size(); i++ ) { zfs_log << "_"; } zfs_log << endl;
  zfs_log << header << endl;
  for ( std::size_t i = 0; i < header.size(); i++ ) { zfs_log << "-"; } zfs_log << endl;
  ZFSId counter = 0;
  ZFSId supCounter = 0;
  if ( s_functionTimings.size() > 0 ) {
    sort( s_functionTimings.begin(), s_functionTimings.end() );
    reverse( s_functionTimings.begin(), s_functionTimings.end() );
    for ( vector<ZFSFunctionTiming>::size_type i = 0; i < s_functionTimings.size(); i++ ) {
      if ( s_functionTimings[i].getInitCpuTime() < m_initCpuTime ) continue;
      const ZFSFloat relCpuTime = 100.0*getCpuTimeSecs(s_functionTimings[i].getDeltaCpuTime())/max(1e-15,getCpuTimeSecs(exitCpuTime-m_initCpuTime));
      const ZFSFloat relWallTime = 100.0*s_functionTimings[i].getDeltaWallTime()/max(1e-15,(exitWallTime-m_initWallTime));
      if ( relCpuTime < thresholdPercentage ) { supCounter++; continue; }
      char buffer[7];
      sprintf(buffer,"%6.2f", relCpuTime);
      char buffer2[7];
      sprintf(buffer2,"%6.2f", relWallTime );
      zfs_log << buffer << "%   " << buffer2 << "%   " << s_functionTimings[i].getName() << endl;
      counter++;
    }
    if ( supCounter > 0 ) {
      zfs_log << "  .....     .....   (" << supCounter << " shorter timings with CPU<" << thresholdPercentage << "% were suppressed)" << endl;
    }
  }
  if ( counter == 0 ) {
    zfs_log << "No timings recorded for timer '" << m_name << "'." << endl; 
  }
  for ( std::size_t i = 0; i < header.size(); i++ ) { zfs_log << "-"; } zfs_log << endl;
  zfs_log << "Total cpu time:  " << printTime(getCpuTimeSecs(exitCpuTime-m_initCpuTime)) << endl;
  zfs_log << "Total wall time: " << printTime(exitWallTime-m_initWallTime) << endl;
  for ( std::size_t i = 0; i < header.size(); i++ ) { zfs_log << "_"; } zfs_log << endl;  
}

ZFSId ZFSProfile::getTimingId( string name ) { 
  ZFSId tId = -1;
  if ( s_functionTimings.size() > 0 ) {
    for ( vector<ZFSFunctionTiming>::size_type i = 0; i < s_functionTimings.size(); i++ ) {
      if ( s_functionTimings[i].getName() == name ) tId = i;
    }
  }
  if ( tId < 0 ) {
    tId = static_cast<ZFSId>(s_functionTimings.size()); 
    s_functionTimings.push_back( ZFSFunctionTiming(name) ); 
  }
  ASSERT( tId > -1, "Non-existing timer" );
  return tId; 
}

ZFSString ZFSProfile::printTime( ZFSFloat secs ) {
  stringstream time;
  time.str("");
  ZFSFloat rem = secs;
  if ( rem > 86400.0 ) {
    const ZFSFloat div = floor(rem/86400.0);
    time << ((ZFSId)div) << " days, ";
    rem -= div*86400.0;
  }
  if ( rem > 3600.0 ) {
    const ZFSFloat div = floor(rem/3600.0);
    time << ((ZFSId)div) << " hours, ";
    rem -= div*3600.0;
  }
  if ( rem > 60.0 ) {
    const ZFSFloat div = floor(rem/60.0);
    time << ((ZFSId)div) << " mins, ";
    rem -= div*60.0;
  }
  time << rem << " secs";
  const ZFSString ret = time.str();
  return ret;
}


ZFSFunctionTiming::ZFSFunctionTiming( string name ) : 
  m_initCpuTime(cpuTime()), m_deltaCpuTime(0), m_tmpCpuTime(0), 
  m_initWallTime(wallTime()), m_deltaWallTime(0.0), m_tmpWallTime(-1.0), m_name(name) {}

ZFSFunctionTiming::~ZFSFunctionTiming() { m_name = "<deleted>"; }

ZFSFunctionTiming& ZFSFunctionTiming::operator=(const ZFSFunctionTiming& t) {
  m_initCpuTime = t.m_initCpuTime; m_deltaCpuTime = t.m_deltaCpuTime; m_tmpCpuTime = t.m_tmpCpuTime;
  m_initWallTime = t.m_initWallTime; m_deltaWallTime = t.m_deltaWallTime; m_tmpWallTime = t.m_tmpWallTime;
  m_name = t.m_name;
  return *this;
}

void ZFSFunctionTiming::in() { 
  m_tmpCpuTime = cpuTime(); 
  m_tmpWallTime = wallTime(); 
}
 
void ZFSFunctionTiming::out() { 
  if ( m_tmpCpuTime > 0 ) m_deltaCpuTime += (cpuTime()-m_tmpCpuTime); 
  if ( m_tmpWallTime > 0.0 ) m_deltaWallTime += (wallTime()-m_tmpWallTime);
  m_tmpCpuTime = -1;
  m_tmpWallTime = -1.0;
}

