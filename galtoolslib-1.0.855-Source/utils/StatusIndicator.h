#ifndef _utl_StatusIndicator_h
#define _utl_StatusIndicator_h

#include <time.h>
#include <sys/time.h>
#include <sstream>

namespace utl {
  /** \brief Nicely formed status indicator, showing percent completed in both numbers and as a bar.
   *
   * This class is OMP thread safe.
   */
  class StatusIndicator {
  public:
    /** \brief Construct a status indicator with a string and number of steps
     *
     * \param marker is the string to print in front of the status indicator
     * \param steps the total number of steps expected.  If steps = 0 then the
     * number of steps is unknown and only accumulated time is displayed
     */
    StatusIndicator(const std::string& marker, size_t steps);
    ~StatusIndicator();
    
    /** \brief Increase the number of finished steps by one and reprint. */
    void refresh(const std::string &text=std::string());
    
    //! Check the status of the indicators
    static bool isOff();
    
    //! Turn off output from ALL status indicators
    static void turnOff();
    
    //! Turn on output from ALL status indicators
    static void turnOn();
    
  private:
    static unsigned long counter;
    static unsigned long last;
    unsigned long fcounter;
    static bool quiet;
    //#pragma omp threadprivate(utl::StatusIndicator::quiet)
    
    std::string fmarker; //!< Store the marker
    size_t fsteps, fcurrentStep; //!< Number of steps and current step.
    long ftInit; //!< Initial time when counter started, used for time estimates
    double ftEstimate; //!< Estimated time needed to complete the project
  };
  
  template <typename T>
    bool from_string(T &t, const std::string &s){
    std::istringstream iss(s);
    return !(iss >> t).fail();
  }
  
}

#endif
