#include "StatusIndicator.h"

#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <unistd.h>

//Check that we are using a terminal
const static bool terminal = isatty(2);

//Initialize the static counter of statusIndicators
unsigned long utl::StatusIndicator::counter = 1;
unsigned long utl::StatusIndicator::last = 0;
bool utl::StatusIndicator::quiet = false;

//Status indicator
utl::StatusIndicator::StatusIndicator(const std::string & marker, size_t steps) : 
  fmarker(marker), 
  fsteps(steps), 
  fcurrentStep(0), 
  ftEstimate(0) 
{
#pragma omp critical (statusIndicator)
  {
    fcounter = ++counter;
  }
  
  ftInit = time(NULL);
  refresh();
}


utl::StatusIndicator::~StatusIndicator(){
  //Reset last if we were the last indicator to give output
#pragma omp critical (statusIndicator)
  {
    if (! quiet ) std::cerr<<std::endl;
    if (last == fcounter){
      last = 0;
    }
  }
}

bool utl::StatusIndicator::isOff() {
  bool out;
#pragma omp critical (statusIndicator)
  {
    out = quiet;
  }
  return out;
}

void utl::StatusIndicator::turnOff(){
#pragma omp critical (statusIndicator)
  {
    quiet = true;
  }
}


void utl::StatusIndicator::turnOn(){
#pragma omp critical (statusIndicator)
  {
    quiet = false;
  }
}


void utl::StatusIndicator::refresh(const std::string &text){
  
  //No output but increase the step size
  bool quit=false;
#pragma omp critical (statusIndicator)
  if (quiet) {
    ++fcurrentStep;
    quit=true;
  }
  
  if (quit)
    return;
  
  
  // Print everything if we are on a terminal, otherwise print only every 10th line
  // for steps and every 50th for open ended statuses.  Should work adequately for gardian.
  
  if ( ! terminal && fsteps != 0 && fcurrentStep % 10 != 0 && fcurrentStep != fsteps ) {
#pragma omp critical (statusIndicator)
    ++fcurrentStep;
    return;
  }
  
  if ( ! terminal && fsteps == 0 && fcurrentStep % 50 != 0 ) {
#pragma omp critical (statusIndicator)
    ++fcurrentStep;
    return;
  }
  
  
  //rotary at the end needs these signs
  static const std::string rot = "\\|/-"; 
  
  int totalTime = time(NULL) - ftInit;
  
#pragma omp critical (statusIndicator)
  {
    //Do not step on the toes of other indicators
    if ( terminal && last != 0 && last != fcounter )
      std::cerr<<std::endl;
    
    //Move to front of line
    if ( terminal ) 
      std::cerr<<"\r";
    else if ( fcurrentStep > 0 || (last != 0 && last != fcounter) )
      std::cerr<<std::endl;
    
    std::cerr<<fmarker<<" ";
    
    if ( fsteps != 0 ) {
      
      //The fraction finished
      double fraction = double(fcurrentStep)/double(fsteps);
      
      //Calculate the percent (integer)
      int percent = int(fraction*100);
      
      //Build the string
      std::string output = "";
      for (int i = 0; i < 20; ++i){
	if (i < percent/5) {
	  output += "#";
	}else{
	  output += "-";
	}
      }
      std::cerr<<output<<" ("<<std::setw(3)<<percent<<"%) ";
      
      //Calculate the estimate time left, basically time from start divided by
      //the fraction finished
      if (fraction > 0) {
	
	if ( ftEstimate == 0 ) {
	  ftEstimate = double(totalTime)/fraction;
	} else {
	  ftEstimate = (ftEstimate + 3*double(totalTime)/fraction)/4.;
	}
	
	int timeHours = int(ftEstimate)/3600;
	int timeMinutes = (int(ftEstimate)-timeHours*3600)/60;
	int timeSeconds = int(ftEstimate)-timeHours*3600 -timeMinutes*60;
	int timeLeftSeconds = int(ftEstimate - totalTime);
	int timeLeftHours = timeLeftSeconds/3600;
	int timeLeftMinutes = (timeLeftSeconds - timeLeftHours*3600)/60;
	timeLeftSeconds = (timeLeftSeconds - timeLeftHours*3600 - timeLeftMinutes*60);
	
	std::cerr<<"[Estimated time left: ";
	if (timeLeftHours){
	  std::cerr<<std::setw(3)<<timeLeftHours<<"h";
	} else {
	  std::cerr<<std::setw(4)<<"";
	}
	
	if (timeLeftMinutes || timeLeftHours){
	  std::cerr<<std::setw(2)<<timeLeftMinutes<<"m";
	} else {
	  std::cerr<<std::setw(3)<<"";
	}
	std::cerr<<std::setw(2)<<timeLeftSeconds<<"s";
	std::cerr<<" / ";
	
	if (timeHours){
	  std::cerr<<std::setw(3)<<timeHours<<"h";
	}
	
	if (timeMinutes || timeHours){
	  std::cerr<<std::setw(2)<<timeMinutes<<"m";
	}
	std::cerr<<std::setw(2)<<timeSeconds<<"s] ";
      }
      
    } else {
      
      //Unlimited steps, only tima accumulated printed
      int timeHours = totalTime/3600;
      int timeMinutes = (totalTime-timeHours*3600)/60;
      int timeSeconds = totalTime-timeHours*3600 -timeMinutes*60;
      
      std::cerr<<"[Time elapsed: ";
      if (timeHours){
	std::cerr<<std::setw(3)<<timeHours<<"h";
      }
      if (timeMinutes || timeHours){
	std::cerr<<std::setw(2)<<timeMinutes<<"m";
      }
      std::cerr<<std::setw(2)<<timeSeconds<<"s] ";
      
    }
    
    if (text != "")
      std::cerr<<text<<" ";
    
    if ( terminal )
      std::cerr<<rot[fcurrentStep % 4];
    
    std::cerr<<std::flush;
    
    //Set the last value so we know what indicator is active
    last = fcounter;
    
    //Increase the step
    ++fcurrentStep;
    
  } //Critical section

}

