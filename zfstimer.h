#ifndef ZFSTIMER_H
#define ZFSTIMER_H

#include <iomanip>
#include <string>
#include <time.h>
#include "zfsconstants.h"
#include "zfsfunctions.h"
#include "zfsglobalvariables.h"
#include "zfsinfoout.h"
#include "zfscompiler.h"
#include "zfslikwid.h"

#if defined(ZFS_MS_COMPILER)
  #include <Windows.h>
#else
  #include <sys/time.h>
  #include <sys/times.h>
  #include <unistd.h>
#endif

#if defined(ZFS_MS_COMPILER)
  static const unsigned __int64 epoch = ((unsigned __int64)116444736000000000ULL);
#endif

#ifdef ZFS_TIMER_FUNCTION
#define NEW_TIMER_GROUP(id, groupName) \
    const ZFSId id = timers().newGroup(groupName)
#define NEW_TIMER(id, timerName, groupId) \
    const ZFSId id = timers().newGroupTimer(timerName, groupId)
#define NEW_SUB_TIMER(id, timerName, timerId) \
    const ZFSId id = timers().newSubTimer(timerName, timerId)
#define NEW_TIMER_GROUP_STATIC(id, groupName) \
    static const ZFSId id = timers().newGroup(groupName)
#define NEW_TIMER_STATIC(id, timerName, groupId) \
    static const ZFSId id = timers().newGroupTimer(timerName, groupId)
#define NEW_SUB_TIMER_STATIC(id, timerName, timerId) \
    static const ZFSId id = timers().newSubTimer(timerName, timerId)
#define NEW_TIMER_GROUP_NOCREATE(id, groupName) \
    id = timers().newGroup(groupName)
#define NEW_TIMER_NOCREATE(id, timerName, groupId) \
    id = timers().newGroupTimer(timerName, groupId)
#define NEW_SUB_TIMER_NOCREATE(id, timerName, timerId) \
    id = timers().newSubTimer(timerName, timerId)
#define START_TIMER(timerId)              timers().startTimer(timerId)
#define RECORD_TIMER_START(timerId)       timers().recordTimerStart(timerId)
#define RECORD_TIMER_STOP(timerId)        timers().recordTimerStop(timerId)
#define RETURN_TIMER(timerId)             timers().returnTimer(timerId)
#define RETURN_TIMER_TIME(timerId)        timers().returnTimerTime(timerId)
#define STOP_TIMER(timerId)               timers().stopTimer(timerId)
#define STOP_ALL_TIMERS()                 timers().stopAllTimers()
#define RECORD_TIMER(timerId)             timers().recordTimer(timerId)
#define RECORD_ALL_TIMER()                timers().recordAllTimer()
#define DISPLAY_TIMER(timerId)            timers().displayTimer(timerId)
#define DISPLAY_TIMER_INTERM(timerId)     timers().displayTimerNoToggleDisplayed(timerId)
#define DISPLAY_TIMER_OFFSET(timerId,ivl) if(globalTimeStep%ivl==0) {timers().recordTimerStop(timerId); timers().displayTimerNoToggleDisplayed(timerId); timers().recordTimerStart(timerId);}
#define DISPLAY_ALL_GROUP_TIMERS(groupId) timers().displayAllTimers(groupId)
#define DISPLAY_ALL_TIMERS()              timers().displayAllTimers()
#define RESET_TIMER(timerId)              timers().resetTimer(timerId)
#define RESET_TIMERS()                    timers().resetTimers()
#define RESET_RECORD(timerId)             timers().resetRecord(timerId)
#define RESET_ALL_RECORDS()               timers().resetRecords()
#define SET_TIMER(timeValue,timerId)      timers().setTimer(timeValue,timerId)
#else
#define NEW_TIMER_GROUP(id, groupName)                  do { } while (false) 
#define NEW_TIMER(id, timerName, groupId)               do { } while (false) 
#define NEW_SUB_TIMER(id, timerName, timerId)           do { } while (false) 
#define NEW_TIMER_GROUP_STATIC(id, groupName)           do { } while (false) 
#define NEW_TIMER_STATIC(id, timerName, groupId)        do { } while (false) 
#define NEW_SUB_TIMER_STATIC(id, timerName, timerId)    do { } while (false) 
#define NEW_TIMER_GROUP_NOCREATE(id, groupName)         do { } while (false) 
#define NEW_TIMER_NOCREATE(id, timerName, groupId)      do { } while (false) 
#define NEW_SUB_TIMER_NOCREATE(id, timerName, timerId)  do { } while (false) 
#define START_TIMER(timerId)                            do { } while (false) 
#define RECORD_TIMER_START(timerId)                     do { } while (false) 
#define RECORD_TIMER_STOP(timerId)                      do { } while (false) 
#define RETURN_TIMER(timerId)                           do { } while (false) 
#define STOP_TIMER(timerId)                             do { } while (false) 
#define STOP_ALL_TIMERS()                               do { } while (false) 
#define RECORD_TIMER(timerId)                           do { } while (false) 
#define RECORD_ALL_TIMER()                              do { } while (false) 
#define DISPLAY_TIMER(timerId)                          do { } while (false) 
#define DISPLAY_ALL_GROUP_TIMERS(groupId)               do { } while (false) 
#define DISPLAY_ALL_TIMERS()                            do { } while (false) 
#define RESET_TIMER(timerId)                            do { } while (false) 
#define RESET_TIMERS()                                  do { } while (false) 
#define RESET_RECORD(timerId)                           do { } while (false) 
#define RESET_ALL_RECORDS()                             do { } while (false)
#define #define SET_TIMER(timeValue,timerId)            do { } while (false)
#endif

/// \brief ZFSTimers manages all ZFS Timers and allows primitive profiling.
///
/// Usage:
/// - NEW_TIMER(string name) creates a new timer with name "name" and
/// returns its index (a ZFSId that you can use to access it).
/// - RESET_TIMER(ZFSId timerId): resets timerId.
/// - START_TIMER(timerId)/STOP_TIMER(timerId) work as expected.
/// - DISPLAY_TIMER(timerId) writes the timerId name and time to the log.
/// - DISPLAY_ALL_TIMERS: displays all timers.
///
/// Example:
///
/// RESET_TIMER;
/// ZFSId myTimer;
/// myTimer = NEW_TIMER("My timer");
///
/// START_TIMER(myTimer);
/// f1(); // Function will be timed.
/// STOP_TIMER(myTimer);
/// f2(); // Function will not be timed.
/// START_TIMER(myTimer);
/// f3(); // Function will be timed.
/// STOP_TIMER(myTimer);
///
/// DISPLAY_TIMER(myTimer);
///
/// NOTE: ZFS_TIMER_FUNCTION macro enable/disables timing ZFS.
///
class ZFSTimers {
  friend ZFSTimers& timers();
 public:
  inline ZFSId newGroup(const std::string groupName);
  inline ZFSId newGroupTimer(const std::string timerName, const ZFSId groupId);
  inline ZFSId newSubTimer(const std::string timerName, const ZFSId timerId);
  inline ZFSFloat returnTimer(const ZFSId timerId);
  inline ZFSFloat returnTimerTime(const ZFSId timerId);
  inline void startTimer(const ZFSId timerId); // start
  inline void stopTimer(const ZFSId timerId); // stop
  inline void resetTimer(const ZFSId timerId); // reset
  inline void recordTimer(const ZFSId timerId); // record
  inline void recordTimerStart(const ZFSId timerId); // reset + start
  inline void recordTimerStop(const ZFSId timerId); // stop + record
  inline void recordTimers();
  inline void stopAllTimers();
  inline void displayTimer(const ZFSId timerId);
  inline void displayTimerNoToggleDisplayed(const ZFSId timerId);
  inline void displayAllTimers();
  inline void displayAllTimers(const ZFSId groupId);
  inline void resetTimers();
  inline void resetRecord(const ZFSId timerId);
  inline void resetRecords();
  inline void setTimer(const ZFSFloat timeValue, const ZFSId timerId);
  
 private:
  ZFSTimers() { }
  ~ZFSTimers() {}
  // delete: copy construction, and copy assignment
  ZFSTimers(ZFSTimers&);
  ZFSTimers& operator=(const ZFSTimers&);

  struct Timer {
    Timer(const std::string n, const ZFSId g)
        : name(n), group(g), cpuTime(0), oldCpuTime(0), recordedTime(0),
          status(Timer::Uninitialized), subTimers(0), displayed(false) {}
    std::string name; ///< Timer Name
    ZFSId group;      ///< Group Id
    ZFSFloat cpuTime;       ///< CPU time
    ZFSFloat oldCpuTime;    ///< Old CPU time (for timer restart)
    ZFSFloat recordedTime;  ///< Time recorded on the timer.
    ZFSId status;           ///< Timer's status, see enum:
    enum { Uninitialized = 0,
         Running = 1,
         Stopped = 0
    };
    std::vector<ZFSId> subTimers;
    bool displayed;
  };

  std::vector<std::string> m_groups;
  std::vector<Timer> m_timers;
  typedef std::vector<Timer>::iterator TIt;
  

  inline ZFSFloat time();
  inline void displayTimer_(const ZFSId timerId, const ZFSBool toggleDisplayed=true, const ZFSId tIndent = 0, const ZFSFloat superTime = -F1);  
  inline void displayTimerHeader_();
  inline void displayTimerGroupHeader_(const ZFSId groupId);
  inline ZFSId indent(const ZFSId pIndent) const { return pIndent + 2; };
};

ZFSTimers& timers();

/// Returns Wall-Clock time in seconds
inline ZFSFloat ZFSTimers::time() {
  /// Timers are not portable.
  /// We should use the C++ <chrono> library instead.
  /// The factors convert the cpu time to seconds
#ifdef ZFS_MPI_TIMER
  return MPI_Wtime();
#else
 #if _POSIX_TIMERS > 0
  timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  return static_cast<ZFSFloat>(t.tv_sec)
      + static_cast<ZFSFloat>(t.tv_nsec / 1000000000.0);
 #else
  struct timeval t;
  gettimeofday(&t, NULL);
  return static_cast<ZFSFloat>(t.tv_sec)
      + static_cast<ZFSFloat>(t.tv_usec / 1000000.0);
 #endif
#endif
}

inline void ZFSTimers::setTimer(const ZFSFloat timeValue, const ZFSId timerId) {
  m_timers[timerId].cpuTime = timeValue;
}

inline void ZFSTimers::resetRecord(const ZFSId timerId) {
  m_timers[timerId].recordedTime = 0.0;
}

inline void ZFSTimers::resetRecords() {
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    resetRecord((ZFSInt) timerId);
  }
}

inline void ZFSTimers::recordTimerStop(const ZFSId timerId) {
  if(timerId < 0){
    return;
  }
#ifdef ZFS_SYNCHRONIZE_TIMERS
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  stopTimer(timerId);
  recordTimer(timerId);
}

inline void ZFSTimers::recordTimer(const ZFSId timerId) {
  m_timers[timerId].recordedTime += returnTimer(timerId);
}
inline void ZFSTimers::recordTimers() {
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
  recordTimer(timerId);
  }
}

inline void ZFSTimers::resetTimer(const ZFSId timerId) {
  m_timers[timerId].cpuTime = 0.0;
  m_timers[timerId].status = Timer::Stopped;
}

inline void ZFSTimers::resetTimers() {
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    resetTimer(timerId);
  }
}

/// Creates a new timer group and returns its groupId.
inline ZFSId ZFSTimers::newGroup(const std::string name) {
    m_groups.push_back( name );
  return m_groups.size() - 1;
}

/// Creates a new timer and returns its timerId.
inline ZFSId ZFSTimers::newGroupTimer(const std::string name, const ZFSId groupId) {
  ASSERT(static_cast<std::size_t>(groupId) < m_groups.size() && groupId > -1,
         "groupId: " << groupId << " does not exists | name: " << name);
  m_timers.push_back(Timer(name, groupId));
  const ZFSId newTimerId = m_timers.size() - 1;
  return newTimerId;
}

/// Creates a new timer and returns its timerId.
inline ZFSId ZFSTimers::newSubTimer(const std::string name, const ZFSId timerId) {
  if(timerId < 0) return -1;

  ASSERT(static_cast<std::size_t>(timerId) < m_timers.size(),
         "timerId " << timerId << " does not exist when trying to create subtimer with name " << name);
         
  const ZFSId groupId = m_timers[ timerId ].group;
  m_timers.push_back(Timer(name, groupId));
  const ZFSId newTimerId = m_timers.size() - 1;
  m_timers[ timerId ].subTimers.push_back(newTimerId);
  return newTimerId;
}

inline void ZFSTimers::recordTimerStart(const ZFSId timerId) {
  if(timerId < 0) {
    return;
  }
#ifdef ZFS_SYNCHRONIZE_TIMERS
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  resetTimer(timerId);
  startTimer(timerId);
}

inline void ZFSTimers::startTimer(const ZFSId timerId) {
  ASSERT(m_timers[timerId].status != Timer::Running,
         "The timer " << m_timers[timerId].name << " with id: " << timerId
         << " can't be started because it is already running!");
  const ZFSFloat t = time();
  m_timers[timerId].oldCpuTime = m_timers[timerId].cpuTime;
  m_timers[timerId].cpuTime = t;
  m_timers[timerId].status = Timer::Running;

  // Enable likwid counter if likwid is enabled
#ifdef WITH_LIKWID
  LIKWID_MARKER_START(std::to_string(timerId).c_str());
#endif
}


/// Returns the timer Value.
inline double ZFSTimers::returnTimer(const ZFSId timerId) {
#if ZFS_TIMERS_AVERAGE_OVER_DOMAINS
    const ZFSFloat t = m_timers[ timerId ].cpuTime;
    ZFSFloat tmp_rcv = 0.0;
    MPI_Reduce(&t, &tmp_rcv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return tmp_rcv/globalNoDomains();
#else
    return m_timers[ timerId ].cpuTime;
#endif
}

//Returns the recorded time
inline double ZFSTimers::returnTimerTime(const ZFSId timerId) {
  return m_timers[timerId].recordedTime;
}

/// Stops the timer and sets its final value.
inline void ZFSTimers::stopTimer(const ZFSId timerId) {
  if(m_timers[timerId].status == Timer::Running) {
    const ZFSFloat t = time();
    m_timers[timerId].cpuTime = t - m_timers[timerId].cpuTime + m_timers[timerId].oldCpuTime;
    m_timers[timerId].status = Timer::Stopped;

    // Stop likwid counter if likwid is enabled
#ifdef WITH_LIKWID
    LIKWID_MARKER_STOP(std::to_string(timerId).c_str());
#endif
  }
}


// Stops all timers.
inline void ZFSTimers::stopAllTimers() {
  for(std::size_t i = 0, e = m_timers.size(); i != e; ++i) {
    if(m_timers[i].status == Timer::Running) {
      stopTimer(i);
    }
  }
}

inline void ZFSTimers::displayTimer_(const ZFSId timerId, const ZFSBool toggleDisplayed, const ZFSId tIndent, const ZFSFloat superTime) {
  bool running = false;
    if (m_timers[timerId].displayed) { return; }

    if (m_timers[timerId].status == Timer::Running) {
      running = true;
      stopTimer(timerId);
    }
    zfs_log.width(50);
    zfs_log.setf(std::ios::left);
    std::stringstream indentedName;

    // Calculate time relative to the parent timer
    ZFSFloat percentage;
    if (superTime < F0) {
      // If the parent time is less than zero, that means that there is no parent timer
      // and the percentage should be 100%
      percentage = 100*F1;
    } else if (approx(superTime, F0, ZFSFloatEps)) {
      // If the parent time is approximately zero, that probably means that the timer was never
      // run - therefore the percentage is set to 0%
      percentage = F0;
    } else {
      // Otherwise calculate the percentage as the fraction of this timer vs. the parent timer times 100%
      percentage = 100*m_timers[timerId].recordedTime/superTime;
    }

    indentedName << std::string(tIndent,' ');
    indentedName << "[" << std::fixed << std::setprecision(1) << std::setw(4) << std::setfill('0') 
                 << std::right << percentage << std::left << "%] ";
    indentedName << m_timers[timerId].name ;
    zfs_log << indentedName.str() << std::right;
    zfs_log.precision(6);
    zfs_log.width(20);
    zfs_log << m_timers[timerId].recordedTime << std::left << " [sec]";
    // Show output of likwid performance counters if likwid is enabled
#ifdef WITH_LIKWID
    // If the timer wasn't called set the MFLops to 0.00
    zfs_log << "   ";
    if (approx(m_timers[timerId].recordedTime, F0, ZFSFloatEps)) {
      zfs_log << "0.00";
    } else {
      zfs_log << "${timer_" << timerId << "}";
    }
    zfs_log << " [DP MFlops/s]";
#endif
    if ( toggleDisplayed ) m_timers[timerId].displayed = true;
    zfs_log  << std::endl;
    for( std::size_t sub = 0, last = m_timers[timerId].subTimers.size(); sub < last; ++sub) {
      const ZFSId new_indent = indent(tIndent);
      displayTimer_(m_timers[timerId].subTimers[sub], toggleDisplayed, new_indent, m_timers[timerId].recordedTime);
    }
    if(running) {
      startTimer(timerId);
    }

}

inline void ZFSTimers::displayTimerHeader_() {
}

inline void ZFSTimers::displayTimerGroupHeader_(const ZFSId groupId) {
  zfs_log <<
      "--------------------------------------------------------------------------------"
   << std::endl;
  zfs_log.width(50);
  zfs_log.precision(12);
  zfs_log.setf(std::ios::left);
  zfs_log << "Group";
  zfs_log.width(40);
  zfs_log << m_groups[groupId] << std::endl;
 
}

inline void ZFSTimers::displayAllTimers() {
  ASSERT(m_timers.size() > 0, "ERROR: no timers have been created!");
  for(std::size_t groupId = 0, e = m_groups.size(); groupId != e; ++groupId) {
    displayAllTimers(groupId);
  }
}

inline void ZFSTimers::displayAllTimers(const ZFSId groupId) {
  ASSERT(m_timers.size() > 0, "ERROR: no timers have been created!");
  ASSERT(static_cast<std::size_t>(groupId) < m_groups.size() && groupId > -1, "ERROR: groupId does not exists");
  for (std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    m_timers[timerId].displayed = false;
  }
  displayTimerGroupHeader_(groupId);
  displayTimerHeader_();
  for (std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    if(m_timers[timerId].group == groupId) {
      displayTimer_(timerId);
    }
  }
  for (std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    m_timers[timerId].displayed = false;
  }
}

inline void ZFSTimers::displayTimer(const ZFSId timerId) {
  ASSERT(static_cast<std::size_t>(timerId) < m_timers.size(), "ERROR: timer timerId does not exist");
  displayTimerHeader_();
  displayTimer_(timerId);
}

inline void ZFSTimers::displayTimerNoToggleDisplayed(const ZFSId timerId) {
  ASSERT(static_cast<std::size_t>(timerId) < m_timers.size(), "ERROR: timer timerId does not exist");
  displayTimerHeader_();
  displayTimer_(timerId,false);
}

//------------------------------------------------------------------------------

class ZFSProfile;
class ZFSFunctionTiming;


inline clock_t cpuTime() {     
  return clock(); 
}

inline ZFSFloat wallTime() { 
  return MPI_Wtime();  
}

/** 
 * \brief This class counts the static execution time of a function
 * \author Lennart Schneiders
 * \date 14.02.2013
 */
class ZFSFunctionTiming
{
  public:
    explicit ZFSFunctionTiming( std::string name );
    ~ZFSFunctionTiming();
    ZFSFunctionTiming& operator=(const ZFSFunctionTiming& t);
    void in();
    void out();
    clock_t getInitCpuTime() const { return m_initCpuTime; }
    clock_t getDeltaCpuTime() const { return m_deltaCpuTime; }
    ZFSFloat getInitWallTime() const { return m_initWallTime; }
    ZFSFloat getDeltaWallTime() const { return m_deltaWallTime; }
    std::string getName() const { return m_name; }
    
  private:
    clock_t m_initCpuTime;
    clock_t m_deltaCpuTime;
    clock_t m_tmpCpuTime;
    ZFSFloat m_initWallTime;
    ZFSFloat m_deltaWallTime;
    ZFSFloat m_tmpWallTime;
    std::string m_name;
};

ZFSBool operator<(const ZFSFunctionTiming& a, const ZFSFunctionTiming& b);

/** 
 * \brief This class collects all function timings and produces a profiling for certain areas of the code
 * \author Lennart Schneiders
 * \date 14.02.2013
 */
class ZFSProfile
{ 
public:
  explicit ZFSProfile( const std::string& name ) : m_initCpuTime(cpuTime()), m_initWallTime(wallTime()), m_name(name) {}
  ~ZFSProfile(); 
  static ZFSId getTimingId( std::string name );
  static ZFSFloat getCpuTimeSecs( clock_t cput ) { return ( static_cast<ZFSFloat>(cput)/static_cast<ZFSFloat>(CLOCKS_PER_SEC) ); }
  static ZFSString printTime( ZFSFloat secs );
  static std::vector<ZFSFunctionTiming> s_functionTimings;

private:
  const clock_t m_initCpuTime;
  const ZFSFloat m_initWallTime;
  const std::string m_name;
  
};

//-------------------------------------------------------


#endif // ZFSTIMER_H
