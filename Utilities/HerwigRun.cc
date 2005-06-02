#include "HerwigRun.h"

using namespace Herwig;
using namespace ThePEG;

const string usage = " init|read|run "
                   " [-N num-events] [-seed random-generator-seed] "
                   "[-d debug-level] [-dHw herwig-debug-level] [-l load-path] "
                   "[-L first-load-path] [-r repo-file] "
                   "[-i initialization file] [run-file]\n";

HerwigRun::HerwigRun() {}

HerwigRun::~HerwigRun() {}

HerwigRun::HerwigRun(int argc, char **argv) 
: N(-1),
  ngen(0),
  repo("HerwigDefaults.rpo"),
  egCreated(false), 
  Status(UNKNOWN),
  timer(".HerwigRun.timer"), 
  isInitialized(false),
  errorFlag(false) 
{
  std::string runType;
  if(argc > 1) runType = argv[1];
  else runType = "";

  if(runType == "init") Status = INIT;
  else if(runType == "read") Status = READ;
  else if(runType == "run") Status = RUN;
  else if(runType == "-h" || runType == "--help") {
    std::cerr << "Usage: " << argv[0] << usage;
    HerwigRun::printHelp(std::cerr);
    errorFlag = true;
    return;
  } else {
    std::cerr << "Usage: " << argv[0] << usage;
    errorFlag = true; 
    return;
  }

  string spaths = SystemUtils::getenv("HERWIG_USER_MODULES");
  vector<string> vpaths = StringUtils::split(spaths, ":");
  for(int i = 0; i<vpaths.size(); i++) DynamicLoader::appendPath(vpaths[i]);
  //DynamicLoader::appendPath(ThePEG_path);
  //DynamicLoader::appendPath(Herwig_path);

  for ( int iarg = 2; iarg < argc; ++iarg ) {
    std::string arg = argv[iarg];
    if ( arg == "-r" ) repo = argv[++iarg];
    else if(arg == "-i") repoin = argv[++iarg];
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
      std::cerr << "Usage: " << argv[0] << usage;
      HerwigRun::printHelp(std::cerr);
      errorFlag = true;
      return;
    }
    else
      run = arg;
  }

  if ( Status == INIT ) {
    breakThePEG();
    {
      HoldFlag<> setup(InterfaceBase::NoReadOnly);
      if ( repoin.empty() ) repoin = "HerwigDefaults.in";
      ifstream is(repoin.c_str());
      Repository::read(is, cout);
      Repository::update();
    }
    Repository::save(repo);
  } else if(Status == READ) {
    Repository::load(repo);
    breakThePEG();
    if ( run.size() && run != "-" ) {
      ifstream is(run.c_str());
      Repository::read(is, std::cout);
    } else {
      Repository::read(std::cin, std::cout, "Herwig++> ");
    }
  } else if ( run.empty() ) {
    std::cerr << "No run-file specified." << endl;
    errorFlag = true;
    return;
  }
}
  
EGPtr HerwigRun::eventGenerator() {
  if(errorFlag || Status != RUN) return EGPtr();
  if(!egCreated) {
    PersistentIStream is(run);
    is >> eg;
    breakThePEG();
    egCreated = true;
    if ( seed > 0 ) eg->random().setSeed(seed);
    return eg;
  } else { return eg; }
}

EventPtr HerwigRun::generateEvent() {
  lastEvent = EventPtr();
  if(Status != RUN || errorFlag) return lastEvent;
  if(!isInitialized) {
    eg->initialize();
    isInitialized = true;
  }
  if(!egCreated) eventGenerator();
  if(eg) {
    if(ngen < N) { ngen++; lastEvent = eg->shoot(); return lastEvent; }
    else {
      std::cerr << "Can only generate " << N << " events.\n";
      return lastEvent;
    }
  } else return lastEvent;
}

int HerwigRun::getN() const { return N; }
int HerwigRun::getNGen() const { return ngen; }
std::string HerwigRun::runName() const { return run; }
std::string HerwigRun::repositoryFile() const { return repo; }
std::string HerwigRun::repositoryInput() const { return repoin; }
HerwigRun::RunStatus HerwigRun::status() const { return Status; }
bool HerwigRun::isRunMode() const { return (status() == RUN); }
bool HerwigRun::isReadMode() const { return (status() == READ); }
bool HerwigRun::isInitMode() const { return (status() == INIT); }

void HerwigRun::printHelp(std::ostream &out) {
  out << endl
      << "One of the following options is required.\n"
      << "==============================================================\n"
      << "init    - Reread default file and create .rpo\n"
      << "read    - Read input file and create run file\n"
      << "run     - Read run file and run\n"
      << "==============================================================\n"
      << "These are optional commands\n"
      << "==============================================================\n"
      << "-N      - Set number of events in the run\n"
      << "-d      - Sets the ThePEG debug level (see ThePEG::Debug)\n"
      << "-dHw    - Sets the Herwig debug level (see Herwig::Debug)\n"
      << "-r      - Changes the repository file from HerwigDefaults.rpo\n"
      << "-i      - Changes the repo input file from HerwigDefaults.in\n"
      << "-l      - Adds path to dynamically load library (to end)\n"
      << "-L      - Adds path to dynamically load library (to beginning)\n"
      << "-seed   - Sets the random seed on initialization\n"
      << "-h      - Displays this help message\n"
      << "==============================================================\n";
}

StepVector HerwigRun::getSteps(EventPtr e) {
  if(!e) e = lastEvent;
  if(!e) {
    std::cerr << "Invalid request: HerwigRun::getSteps on a null event\n";
    return StepVector();
  }
  return e->primaryCollision()->steps();
}

tPVector HerwigRun::getFinalState(int step, EventPtr e) {
  if(!e) e = lastEvent;
  if(!e) {
    std::cerr << "Invalid request: HerwigRun::getFinalState on a null event\n";
    return tPVector();
  }
  if(step < 0) return e->getFinalState();
  else return e->primaryCollision()->step(step)->getFinalState();
}

ParticleSet HerwigRun::getAllParticles(int step, EventPtr e) {
  if(!e) e = lastEvent;
  if(!e) {
    std::cerr << "Invalid request: HerwigRun::getFinalState on a null event\n";
    return ParticleSet();
  }
  return e->primaryCollision()->step(step)->all();
}

ParticleSet HerwigRun::getIntermediates(int step, EventPtr e) {
  if(!e) e = lastEvent;
  if(!e) {
    std::cerr << "Invalid request: HerwigRun::getFinalState on a null event\n";
    return ParticleSet();
  }
  return e->primaryCollision()->step(step)->intermediates();
}

ParticleSet HerwigRun::getOutgoing(int step, EventPtr e) {
  if(!e) e = lastEvent;
  if(!e) {
    std::cerr << "Invalid request: HerwigRun::getFinalState on a null event\n";
    return ParticleSet();
  }
  return e->primaryCollision()->step(step)->particles();
}

bool HerwigRun::preparedToRun() {
  if(eventGenerator()) return true;
  else return false;
}
