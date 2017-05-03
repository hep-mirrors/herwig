// -*- C++ -*-
//
// ProcessHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ProcessHandler class.
//

#include "ProcessHandler.h"
#include "MPISampler.h"

#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Cuts/Cuts.h"

#include "Herwig/Utilities/GaussianIntegrator.h"

#include "gsl/gsl_sf_bessel.h"


#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ProcessHandler.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ProcessHandler::ProcessHandler()
  : theBinStrategy(2) {}

ProcessHandler::ProcessHandler(const ProcessHandler & x)
  : Interfaced(x), LastXCombInfo<>(x), 
    theSampler(x.theSampler), theHandler(x.theHandler), 
    theCuts(x.theCuts), theSubProcess(x.theSubProcess),
    theXCombs(x.theXCombs), theXSecs(x.theXSecs),
    theBinStrategy(x.theBinStrategy), theMEXMap(x.theMEXMap),
    theMaxDims(x.theMaxDims) {}


ProcessHandler::~ProcessHandler() {}

void ProcessHandler::initialize(tSubHdlPtr sub, tCutsPtr cut, tEHPtr eh) {
  /*
    This method should be called during the "read" phase! However due to
    the unpredictable order in which objects are set up by ThePEG it is not 
    guaranteed that the EventHandler is present. So we move that to 
    doinitrun of the MPIHandler!
  */

  if( !eh || !sub || !cut )
    throw Exception() << "ProcessHandler was created without specifying a SubProcess "
		      << "or a Cut Object or the EventHandler."
		      << Exception::runerror;
  theHandler = eh;
  theCuts = cut;
  theSubProcess = sub;

  theSampler = new_ptr(MPISampler());

  Energy maxEnergy = lumiFn().maximumCMEnergy();

  xCombs().clear();
  xSecs().clear();

  cuts()->initialize(sqr(maxEnergy), lumiFn().Y());

  CutsPtr kincuts = (*subProcess()).cuts()? (*subProcess()).cuts(): cuts();
  if ( (*subProcess()).cuts() ) kincuts->initialize(sqr(maxEnergy), lumiFn().Y());
  PExtrPtr pextract = (*subProcess()).pExtractor();

  // Use an empty ckkw handler for the additional interactions:
  tCascHdlPtr ckkw = tCascHdlPtr();

  PartonPairVec vpc = pextract->getPartons(maxEnergy, incoming(), *kincuts);

  // The last parton bin pair was in fact the bins corresponding to
  // the incoming particles, so we remove them, but save them to
  // keep them referenced.
  PBPair orig = vpc.back();
  vpc.pop_back();

  for ( PartonPairVec::iterator ppit = vpc.begin();
	ppit != vpc.end(); ++ppit )
    for ( MEVector::const_iterator meit = (*subProcess()).MEs().begin();
	  meit != (*subProcess()).MEs().end(); ++meit )
      addME(maxEnergy, subProcess(), pextract, kincuts, ckkw, *meit, *ppit);


  xSecs().resize(xCombs().size());

  theMaxDims.clear();
  switch ( binStrategy() ) {
  case 0: {
    theMaxDims.push_back(0);
    for ( int i = 0, N = xCombs().size(); i < N; ++i )
      theMaxDims[0] = max(theMaxDims[0], xCombs()[i]->nDim());
    break;
  }
  case 1: {
    for ( int i = 0, N = xCombs().size(); i < N; ++i )
      theMEXMap[xCombs()[i]->matrixElement()].push_back(xCombs()[i]);
    MEXMap::const_iterator mei = theMEXMap.begin();
    for ( int i = 0, N = theMEXMap.size(); i < N; ++i, ++mei) {
      theMaxDims.push_back(0);
      for ( int j = 0, M = mei->second.size(); j < M; ++j )
        theMaxDims[i] = max(theMaxDims[i], mei->second[j]->nDim());
    }
    break;
  }
  case 2: {
    for ( int i = 0, N = xCombs().size(); i < N; ++i )
      theMaxDims.push_back(xCombs()[i]->nDim());
    break;
  }
  }

  tMPISamplerPtr smplr = dynamic_ptr_cast<tMPISamplerPtr>(sampler());
  smplr->setProcessHandler(this);
  //sample 2 PSpoints to check for zero xsec in advance
  smplr->initialize();
}

void ProcessHandler::
addME(Energy maxEnergy, tSubHdlPtr sub, tPExtrPtr extractor, tCutsPtr cuts,
      tCascHdlPtr ckkw, tMEPtr me, const PBPair & pBins) {

  typedef MEBase::DiagramVector DiagramVector;
  typedef map<string,DiagramVector> DiagramMap;
  cPDPair pin(pBins.first->parton(), pBins.second->parton());
  DiagramVector diag = me->diagrams();
  DiagramMap tdiag;
  DiagramMap tmdiag;
  for ( int i = 0, N = diag.size(); i < N; ++i ) {
    if ( diag[i]->partons()[0] == pin.first &&
         diag[i]->partons()[1] == pin.second )
      tdiag[diag[i]->getTag()].push_back(diag[i]);
    if ( diag[i]->partons()[0] == pin.second &&
         diag[i]->partons()[1] == pin.first )
      tmdiag[diag[i]->getTag()].push_back(diag[i]);
  }

  bool mirror = false;
  if ( ( mirror = tdiag.empty() ) ) tdiag = tmdiag;
  for ( DiagramMap::iterator dit = tdiag.begin(); dit != tdiag.end(); ++dit ) {

    //todo: hope that it is no problem that I take the EventHandler here and not the ProcessHandler:
    StdXCombPtr xcomb =
      new_ptr(StandardXComb(maxEnergy, incoming(), eventHandler(), 
			    sub, extractor, ckkw, pBins, cuts, me, dit->second, mirror));

    if ( xcomb->checkInit() ) xCombs().push_back(xcomb);

    else generator()->logWarning( 
      InitError() << "The matrix element '"
      << xcomb->matrixElement()->name() << "' cannot generate the diagram '"
      << dit->first << "' when used together with the parton extractor '"
      << xcomb->pExtractor()->name()
      << "'. The corresponding diagram is switched off." << Exception::warning);
  }
}

tStdXCombPtr ProcessHandler::select(int bin, double weight) {
  int i = upper_bound(xSecs().begin(), xSecs().end(), UseRandom::rnd()*xSecs().back())
    - xSecs().begin();
  tStdXCombPtr lastXC;
  switch ( binStrategy() ) {
  case 0:
    lastXC = xCombs()[i];
    break;
  case 1: {
    MEXMap::iterator mei = theMEXMap.begin();
    for ( int j = 0; j < bin; ++j) ++mei;
    lastXC = mei->second[i];
    break;
  }
  case 2:
    lastXC = xCombs()[bin];
    break;
  }
  // clean up the old XComb object before switching to a new one
  if ( theLastXComb && theLastXComb != lastXC ) theLastXComb->clean();
  theLastXComb = lastXC;

  lastXC->select(weight);
  lastXC->accept();
  lastXC->matrixElement()->setXComb(lastXC);
  return lastXC;
}


CrossSection ProcessHandler::
dSigDR(const pair<double,double> ll, Energy2 maxS,
       int ibin, int nr, const double * r) {

  PPair inc = make_pair(incoming().first->produceParticle(),
                        incoming().second->produceParticle());
  SimplePhaseSpace::CMS(inc, maxS);

  XVector xv;
  switch ( binStrategy() ) {
  case 0:
    xv = xCombs();
    break;
  case 1: {
    MEXMap::iterator mei = theMEXMap.begin();
    for ( int i = 0; i < ibin; ++i) ++mei;
    xv = mei->second;
    break;
  }
  case 2:
    xv = XVector(1, xCombs()[ibin]);
    break;
  }

  xSecs().resize(xv.size());
  for ( int i = 0, N = xv.size(); i < N; ++i ) xv[i]->prepare(inc);
  CrossSection sum = 0.0*nanobarn;
  for ( int i = 0, N = xv.size(); i < N; ++i )
    xSecs()[i] = ( sum += xv[i]->dSigDR(ll, nr, r) );

  return sum;
}


CrossSection ProcessHandler::dSigDR(const vector<double> & r) {
  double jac = 1.0;
  pair<double,double> ll = lumiFn().generateLL(&r[0], jac);
  Energy2 maxS = sqr(lumiFn().maximumCMEnergy())/exp(ll.first + ll.second);
  int bin = sampler()->lastBin();
  CrossSection x = jac*lumiFn().value(incoming(), ll.first, ll.second)
    *dSigDR(ll, maxS, bin, nDim(bin) - lumiDim(), &r[lumiDim()]);
  return x;
}

int ProcessHandler::nBins() const {
  switch ( binStrategy() ) {
  case 0: return 1;
  case 1: return theMEXMap.size();
  case 2: return xCombs().size();
  }
  return -1;
}

void ProcessHandler::doinitrun() {
  Interfaced::doinitrun();

  sampler()->initrun();

  for ( int i = 0, N = xCombs().size(); i < N; ++i )
    xCombs()[i]->reset();

  double weight(0);
  //sample N PSpoints to get an estimate of the xsec
  for(unsigned int i=0; i<100000; i++){
    weight = sampler()->generate();
    tStdXCombPtr lastXC = select(sampler()->lastBin(), weight);
  }
}

CrossSection ProcessHandler::integratedXSec() const {
  if ( sampler()->integratedXSec() == 0.0*nanobarn )
    return sampler()->maxXSec();

  Stat tot;
  for ( int i = 0, N = xCombs().size(); i < N; ++i ) {
    const StandardXComb & x = *xCombs()[i];
    Stat s;
    s = Stat(x.stats().attempts(), x.stats().accepted(),
             x.stats().sumWeights(), sampler()->integratedXSec(),
             sampler()->sumWeights());
    tot += s;
  }
  return tot.xSec();
} 

void ProcessHandler::statistics(ostream & os, Stat & tot) const {

  if ( statLevel() == 0 ) return;
  map<cPDPair, Stat> partonMap;
  map<MEPtr, Stat> meMap;
  map<PExtrPtr, Stat> extractMap;
  //  Stat tot;

  for ( int i = 0, N = xCombs().size(); i < N; ++i ) {
    const StandardXComb & x = *xCombs()[i];
    Stat s;
    s = Stat(x.stats().attempts(), x.stats().accepted(),
             x.stats().sumWeights(), sampler()->integratedXSec(),
             sampler()->sumWeights());
    partonMap[x.partons()] += s;
    meMap[x.matrixElement()] += s;
    extractMap[x.pExtractor()] += s;
    tot += s;
  }

  string line = "======================================="
    "=======================================\n";

  if ( tot.accepted <= 0 ) {
    os << line << "No events generated by event handler '" << name() << "'."
       << endl;
    return;
  }

  os //<< line << "Statistics for event handler \'" << name() << "\':\n"
     << "                                       "
     << "generated    number of    Cross-section\n"
     << "                                       "
     << "   events     attempts             (nb)\n";

  os << line << "Total:" << setw(42) << tot.accepted << setw(13)
     << tot.attempted << setw(17) << tot.xSec()/nanobarn << endl
     << line;

  if ( statLevel() == 1 ) return;

  os << "Per matrix element breakdown:\n";
  for ( map<MEPtr, Stat>::iterator i = meMap.begin();
        i != meMap.end(); ++i ) {
    string n = i->first->name();
    n.resize(37, ' ');
    os << n << setw(11) << i->second.accepted << setw(13)
       << i->second.attempted << setw(17) << i->second.xSec()/nanobarn << endl;
  }
  os << line;

  if ( statLevel() == 2 ) return;

  os << "Per parton extractor breakdown:\n";
  for ( map<PExtrPtr, Stat>::iterator i = extractMap.begin();
        i != extractMap.end(); ++i ) {
    string n = i->first->name();
    n.resize(37, ' ');
    os << n << setw(11) << i->second.accepted << setw(13)
       << i->second.attempted << setw(17) << i->second.xSec()/millibarn << endl;
  }
  os << line;

  os << "Per incoming partons breakdown:\n";
  for ( map<cPDPair, Stat>::iterator i = partonMap.begin();
        i != partonMap.end(); ++i ) {
    string n = i->first.first->PDGName() + " " + i->first.second->PDGName();
    n.resize(37, ' ');
    os << n << setw(11) << i->second.accepted << setw(13)
       << i->second.attempted << setw(17) << i->second.xSec()/millibarn << endl;
  }
  os << line;

  if ( statLevel() == 3 ) return;

  os << "Detailed breakdown:\n";
  double xsectot = sampler()->integratedXSec()/
    (sampler()->sumWeights()*nanobarn);
  for ( int i = 0, N = xCombs().size(); i < N; ++i ) {
    const StandardXComb & x = *xCombs()[i];
    os << "(" << x.pExtractor()->name() << ") "
       << x.partons().first->PDGName() << " "
       << x.partons().second->PDGName()

       << " (" << x.matrixElement()->name() << " "
       << x.lastDiagram()->getTag() << ") " << endl
       << setw(48) << x.stats().accepted() << setw(13) << x.stats().attempts()
       << setw(17) << x.stats().sumWeights()*xsectot << endl;
  }

  os << line;
}

void ProcessHandler::persistentOutput(PersistentOStream & os) const {
  os << theBinStrategy << theSubProcess << theCuts << theLastXComb
     << theXCombs << ounit(theXSecs, nanobarn)
     << theMaxDims << theMEXMap
     << theSampler << theHandler;
}

void ProcessHandler::persistentInput(PersistentIStream & is, int) {
  is >> theBinStrategy >> theSubProcess >> theCuts >> theLastXComb
     >> theXCombs >> iunit(theXSecs, nanobarn)
     >> theMaxDims >> theMEXMap
     >> theSampler >> theHandler;
}

ClassDescription<ProcessHandler> ProcessHandler::initProcessHandler;
// Definition of the static class description member.

void ProcessHandler::Init() {

  static ClassDocumentation<ProcessHandler> documentation
    ("There is soon documentation for the ProcessHandler class");


  /*
   * Object will be created outside of *.in files 
   *
  static Switch<ProcessHandler,int> interfaceBinStrategy
    ("BinStrategy",
     "The strategy to be used when sampling different ThePEG::XComb "
     "objects. An ThePEG::XComb objet represents a pair of incoming "
     "parton types as defined by a ThePEG::PartonExtractor and a "
     "matrix element.",
     &ProcessHandler::theBinStrategy, 2, false, false);

  static SwitchOption interfaceBinStrategy0
    (interfaceBinStrategy,
     "AllAtOnce",
     "All bins are sampled together.",
     0);

  static SwitchOption interfaceBinStrategy1
    (interfaceBinStrategy,
     "PerME",
     "All bins which have the same matrix element object are sampled together.",
     1);

  static SwitchOption interfaceBinStrategy2
    (interfaceBinStrategy,
     "Individual",
     "All bins are sampled individually.",
     2);
  */

}
