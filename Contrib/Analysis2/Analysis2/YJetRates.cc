// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#include "YJetRates.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"

#include "ThePEG/Utilities/StringUtils.h"

#include "ThePEG/Utilities/Throw.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "YJetRates.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Analysis2;

YJetRates::~YJetRates() {}

void YJetRates::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _nMax << _options << _yn;
}

void YJetRates::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _nMax >> _options >> _yn;
}

ClassDescription<YJetRates> YJetRates::initYJetRates;
// Definition of the static class description member.

void YJetRates::Init() {

  static ClassDocumentation<YJetRates> documentation
    ("Analysis for jet rates");

  static Parameter<YJetRates,unsigned int> interfaceNMax
    ("NMax",
     "Maximum n up to which y(n,n+1) should be booked",
     &YJetRates::_nMax, 5, 2, 0,
     false, false, Interface::lowerlim);


  static Command<YJetRates> interfaceY
    ("R",
     "Insert options for given n.",
     &YJetRates::R, false);

}

string YJetRates::R (string in) {

  istringstream m (StringUtils::car(in));
  unsigned int mult; m >> mult;
  if (mult < 2 || mult > _nMax) {
    Throw<InitException> () << "YJetRates : n must be between 2 and nMax";
  }

  if (_options.size() < _nMax+1) {
    _options.resize(_nMax+1,"");
  }

  _options[mult] = StringUtils::cdr(in);

  return "";

}

void YJetRates::doinit() throw(InitException) {

  Analysis2Base::doinit();

  _yn.clear();

  if (!jetFinder())
    Throw<InitException>() << "YJetRates : No JetFinder has been set, giving up.";
  if (_options.size() < 3)
    Throw<InitException>() << "YJetRates : Need at least options for R_2";
  if (_options.size() < _nMax+1)
    Throw<InitException>() << "YJetRates : It seems that no options where set.";
  
  int plotFlags = HistogramOutput::Frame | HistogramOutput::Errorbars;
  Histogram2Options options (plotFlags);
  options.differential = false;

  _yn.resize(_nMax+1,"");

  for (unsigned int i=2; i<=_nMax; ++i) {
    ostringstream name("");
    name << "R_" << i;
    _yn[i] = name.str();
    insert(name.str(),_options[i],options);
  }
}

void YJetRates::analyze(const tPVector &) {

  pair<vector<Lorentz5Momentum>,double> ev;

  while (*eventExtractor() >> ev) {

    jetFinder()->use(ev.first);

    vector<double> mergings (_nMax+1,0.);

    for (unsigned int i=2; i<=_nMax; ++i)
      mergings[i] = jetFinder()->getYMerge(i);

    tcHistogram2Ptr refhisto = histogram(_yn[2]);
    assert(refhisto);

    _sumWeights += ev.second;

    for (vector<pair<double,double> >::const_iterator b = refhisto->binning().begin();
	 b != refhisto->binning().end(); ++b) {

      double y = (b->first+b->second)/2.;
      
      unsigned int mult = _nMax;

      while (mult > 1) {

	if (mergings[mult] > y)
	  break;

	--mult;

      }

      if (++mult > _nMax)
	continue;

      book(y,_yn[mult],ev.second);
      
    }

  }

}

void YJetRates::dofinish() {

  for (unsigned int i =2; i<=_nMax; ++i) {
    histogram(_yn[i])->rescale(1./_sumWeights);
    finish(_yn[i]);
  }

  Analysis2Base::dofinish();

}
