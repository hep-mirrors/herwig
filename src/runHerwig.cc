#include "Pythia7/Repository/EventGenerator.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/PDT/StandardMatchers.h"
#include "Pythia7/PDT/PYDECYDummy.h"
#include "Pythia7/Utilities/Debug.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Utilities/Timer.h"
#include "Pythia7/Utilities/DynamicLoader.h"
#include "Pythia7/Misc/Exception.h"

int main(int argc, char * argv[]) {
  using namespace Pythia7;

  string run;
  long N = -1;
  long seed = 0;

  for ( int iarg = 1; iarg < argc; ++iarg ) {
    string arg = argv[iarg];
    if ( arg == "-r" ) run = argv[++iarg];
    else if ( arg == "-l" ) DynamicLoader::appendPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-l" )
      DynamicLoader::appendPath(arg.substr(2));
    else if ( arg == "-L" ) DynamicLoader::prependPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-L" )
      DynamicLoader::prependPath(arg.substr(2));
    else if ( arg == "-d" ) Debug::level = atoi(argv[++iarg]);
    else if ( arg.substr(0,2) == "-d" && arg.substr(0,4) != "-dHw" )
	Debug::level = atoi(arg.substr(2).c_str());
    else if ( arg == "-dHw" ) Herwig::HwDebug::level = atoi(argv[++iarg]);
    else if ( arg.substr(0,4) == "-dHw" )
	Herwig::HwDebug::level = atoi(arg.substr(4).c_str());
    else if ( arg == "-N" ) N = atoi(argv[++iarg]);
    else if ( arg.substr(0,2) == "-N" ) N = atoi(arg.substr(2).c_str());
    else if ( arg == "-seed" ) seed = atoi(argv[++iarg]);
    else if ( arg == "-h" ) {
    cerr << "Usage: " << argv[0]
	 << " [-N num-events] [-seed random-generator-seed] [-d debug-level] [-dHw herwig-debug-level] [-l load-path] [-L first-load-path] run-file"
	 << endl;
      return 3;
    }
    else
      run = arg;
  }

  if ( run.empty() ) {
    cerr << "No run-file specified." << endl;
    return 1;
  }

  try {

    PersistentIStream is(run);
    EGPtr eg;
    is >> eg;

    MainTimer timer(".runHerwig.timer");

    breakPythia7();

    if (eg) {
      if ( seed > 0 ) eg->random().setSeed(seed);
      eg->go(1, N);
    } else std::cout<<"eg = nil"<<endl;

  }
  catch ( std::exception & e ) {
    cerr << e.what() << endl;
    return 1;
  }
  catch ( ... ) {
    cerr << "Unknown Exception\n";
    return 2;
  }

  return 0;
}

