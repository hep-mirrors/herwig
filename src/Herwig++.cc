#include "Herwig++/Utilities/HerwigRun.h"
// Any headers needed for analysis go here



int main(int argc, char * argv[]) {
  try {
    Herwig::HerwigRun hw(argc,argv);

    if (hw.isRunMode() && hw.preparedToRun()) {
      for(int i = 0; i<hw.getN(); i++) {
	hw.generateEvent();
	// Add analysis code here
	// To retrieve the particles at the end of the event use
	// tPVector particles = hw.getFinalState();
	// To get the particles at the end of a step, e.g. after showering
	// use (for step s)
	// tPVector particles = hw.getFinalState(s);
	// Then do analysis
      }
    } else if(hw.isRunMode()) {
      std::cout << "Error: Expecting a run but there is no EventGenerator!\n"
		<< std::endl;
    }

  }
  catch ( std::exception & e ) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  catch ( ... ) {
    std::cerr << "Unknown Exception\n";
    return 2;
  }

  return 0;
}

