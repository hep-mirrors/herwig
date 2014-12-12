#include <fstream>
#include "ParallelRunAnalysis.h"
#include "ThePEG/Handlers/SamplerBase.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <boost/asio.hpp>

using namespace Herwig;

ParallelRunAnalysis::ParallelRunAnalysis() {}

void ParallelRunAnalysis::doinitrun() {
  string logfilename = generator()->runName() + ".parallel";
  ofstream log(logfilename.c_str(),ofstream::app);
  log << "hostname> " << boost::asio::ip::host_name() << "\n" << flush;
  log.close();
}

void ParallelRunAnalysis::dofinish() {
  AnalysisHandler::dofinish();
}

void ParallelRunAnalysis::analyze(tEventPtr, long currev, int, int) {
  long totev = generator()->N();
  long i = currev;
  bool skip = currev%(max(totev/100, 1L));
  if ( i > totev/2 ) i = totev-i;
  while ( skip && i >= 10 && !(i%10) ) i /= 10;
  if ( i == 1 || i == 2 || i == 5 ) skip = false;
  if ( skip && currev%5000!=0) return;

  tEHPtr curEvtHandler = generator()->currentEventHandler();
  long attempts = (dynamic_ptr_cast<StdEHPtr>(curEvtHandler))->sampler()->attempts();
  char str[128];
  sprintf(str,"event> %lu/%lu/%lu xs = %.10E pb +/- %.10E pb\n",
          currev,attempts,totev,
          double(curEvtHandler->integratedXSec()/picobarn),
          double(curEvtHandler->integratedXSecErr()/picobarn));
  string logfilename = generator()->runName() + ".parallel";
  ofstream log(logfilename.c_str(),ofstream::app);
  log << str << flush;
  log.close();
}

NoPIOClassDescription<ParallelRunAnalysis> ParallelRunAnalysis::initParallelRunAnalysis;
// Definition of the static class description member.

void ParallelRunAnalysis::Init() {

  static ClassDocumentation<ParallelRunAnalysis> documentation
    ("Analysis for combining parallel runs with Herwig.");

}

