#include "Herwig++/Utilities/HerwigRun.h"
#include "versionstring.h"

int main(int argc, char * argv[]) {
  setVersionString();
  try {
    Herwig::HerwigRun hw(argc,argv);
    if (hw.isRunMode() && hw.preparedToRun()) {
      int step = hw.getN() >= 100 ? hw.getN() / 100 : 1 ;
      for(int i = 0; i<hw.getN(); i++) {
	hw.generateEvent();
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
    return EXIT_FAILURE;
  }
  std::cout << std::endl;
  return EXIT_SUCCESS;
}

