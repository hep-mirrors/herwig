#include "Herwig++/Utilities/HerwigRun.h"
#include "ThePEG/Utilities/Timer.h"
// Any headers needed for analysis go here

int main(int argc, char * argv[]) {
  try {
    // create the main timer object
    Herwig::HerwigRun hw(argc,argv);
    ThePEG::MainTimer timer(".timer");
    if (hw.isRunMode() && hw.preparedToRun()) {
      int step = hw.getN() / 100;
      for(int i = 0; i<hw.getN(); i++) {
	hw.generateEvent();
	// Add analysis code here
	// To retrieve the particles at the end of the event use
	// tPVector particles = hw.getFinalState();
	// To get the particles at the end of a step, e.g. after showering
	// use (for step s)
	// tPVector particles = hw.getFinalState(s);
	// Then do analysis
	if ((i+1) % step == 0)
	  std::cout << "Generated event: " << i+1 
		    << " of " << hw.getN() << "\r" << std::flush;
      }
      hw.eventGenerator()->finalize();
    } else if(hw.isRunMode()) {
      std::cerr << "Error: Expecting a run but there is no EventGenerator!\n";
    }

  }
  catch ( std::exception & e ) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  std::cout << std::endl;
  return 0;
}

