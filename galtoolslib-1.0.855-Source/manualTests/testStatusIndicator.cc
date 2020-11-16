#include <StatusIndicator.h>
#include <iostream>
#include <sstream>
#include <unistd.h>

int main() {

   std::cerr<<"Test program for status indicators, requires manual inspection."<<std::endl;
   std::cerr<<"Output depends on if it is shown on a terminal (detected by isatty)."<<std::endl;
   std::cerr<<"If it is shown on a terminal we should have a line updating itself at each step."<<std::endl;
   std::cerr<<"Otherwise it should update every 10th time when we have known number of steps"<<std::endl;
   std::cerr<<"or every 50th time when number of steps is unknows."<<std::endl;

   std::cerr<<std::endl;
   std::cerr<<"We'll start with a limited step status indicator of 100 steps."<<std::endl;
   std::cerr<<std::endl;

   utl::StatusIndicator *status = new utl::StatusIndicator("100 steps", 100);

   for (size_t i(0); i<100; ++i) {
      status->refresh();
      usleep(50000); //Sleep .05 seconds
   }

   delete status;

   std::cerr<<std::endl;
   std::cerr<<"Then an unlimited step status indicator."<<std::endl;
   std::cerr<<std::endl;

   status = new utl::StatusIndicator("Unlimited",0);

   for (size_t i(0); i<300; ++i) {
      status->refresh();
      usleep(20000); //Sleep .02 seconds
   }

   delete status;

   std::cerr<<std::endl;
   std::cerr<<"Now we are going to turn the status off in the middle"<<std::endl;
   std::cerr<<std::endl;

   status = new utl::StatusIndicator("Turning off status with 50 steps", 50);

   for (size_t i(0); i<15; ++i) {
      status->refresh();
      usleep(200000); //Sleep .2 seconds
   }

   std::cerr<<std::endl<<std::endl;
   std::cerr<<"Nothing should be output between this line"<<std::endl;

   utl::StatusIndicator::turnOff();
   if (! utl::StatusIndicator::isOff() )
      std::cerr<<"The isOff function is returning bogus information"<<std::endl;

   for (size_t i(0); i<20; ++i) {
      status->refresh();
      usleep(200000); //Sleep .2 seconds
   }
   utl::StatusIndicator::turnOn();

   std::cerr<<"and this line but the counter should have advanced 20 steps"<<std::endl;
   std::cerr<<std::endl;

   for (size_t i(0); i<15; ++i) {
      status->refresh();
      usleep(200000); //Sleep .2 seconds
   }

   delete status;

   std::cerr<<std::endl;
   std::cerr<<"There should always be a single empty line between these messages and the status indicators."<<std::endl;
   std::cerr<<"Now we will have two indicators that should be interleaved with no empty lines between them."<<std::endl;
   std::cerr<<std::endl;

   status = new utl::StatusIndicator("Status 1", 15);
   utl::StatusIndicator *status2 = new utl::StatusIndicator("Status 2", 15);

   for (size_t i(0); i < 15; ++i) {
      status->refresh();
      status2->refresh();
      usleep(100000);
   }

   delete status;
   delete status2;

   std::cerr<<std::endl;
   std::cerr<<"Finally we'll try the additional text output using an unlimited status indicator that counts the number of steps."<<std::endl;
   std::cerr<<std::endl;

   status = new utl::StatusIndicator("Unlimited counting",0);

   for (size_t i(0); i < 200; ++i) {
      std::ostringstream os;
      os << "Step number " << i+1;
      status->refresh(os.str());
      usleep(50000);
   }

   delete status;

   return 0;
}
