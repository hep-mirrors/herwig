#include "Pythia7/Repository/Repository.h"
#include "Pythia7/PDT/EnumParticles.h"
#include "Pythia7/Utilities/Debug.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Utilities/DynamicLoader.h"
#include <fstream>


int main(int argc, char * argv[]) {

  using namespace Pythia7;
  using namespace Herwig;                               

  Debug::level   = 1;
  HwDebug::level = 1;
  string repo =  "HerwigDefaults.rpo";                  
  string file;
  bool init = false;

  for ( int iarg = 1; iarg < argc; ++iarg ) {
    string arg = argv[iarg];
    if ( arg == "-d" ) Debug::level = atoi(argv[++iarg]);
    else if ( arg.substr(0,2) == "-d" && arg.substr(0,4) != "-dHw" )
	Debug::level = atoi(arg.substr(2).c_str());
    else if ( arg == "-dHw" ) HwDebug::level = atoi(argv[++iarg]);
    else if ( arg.substr(0,4) == "-dHw" )
	HwDebug::level = atoi(arg.substr(4).c_str());
    else if ( arg == "-r" ) repo = argv[++iarg];
    else if ( arg == "-init" ) {
      init = true;
      Debug::level = 0;
    }
    else if ( arg == "-l" ) DynamicLoader::appendPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-l" )
      DynamicLoader::appendPath(arg.substr(2));
    else if ( arg == "-L" ) DynamicLoader::prependPath(argv[++iarg]);
    else if ( arg.substr(0,2) == "-L" )
      DynamicLoader::prependPath(arg.substr(2));
    else if ( arg == "-h" ) {
      cerr << "Usage: " << argv[0]
	   << " {cmdfile} [-d debug-level] [-dHw herwig-debug-level]"
	   << " [-r input-repository-file]"
	   << " [-l load-path] [-L first-load-path]" << endl;
      return 3;
    }
    else
      file = arg;
  }

  try {

		if ( init ) {
    	breakPythia7();
      {
				HoldFlag<> setup(InterfaceBase::NoReadOnly);
				if ( file.empty() ) file = "HerwigDefaults.in";
				ifstream is(file.c_str());
				Repository::read(is, cout);
				Repository::update();
      }
      Repository::save(repo);
    } else {
      Repository::load(repo);
      breakPythia7();
      if ( file.size() && file != "-" ) {
				ifstream is(file.c_str());
				Repository::read(is, cout);
      } else {
				Repository::read(cin, cout, "Herwig++> ");
      }
    }
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



