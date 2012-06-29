// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#include "YMerge.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"

#include "ThePEG/Utilities/StringUtils.h"

#include "ThePEG/Utilities/Throw.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "YMerge.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Analysis2;

YMerge::~YMerge() {}

void YMerge::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _nMax << _options << _yn;
}

void YMerge::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _nMax >> _options >> _yn;
}

ClassDescription<YMerge> YMerge::initYMerge;
// Definition of the static class description member.

void YMerge::Init() {

  static ClassDocumentation<YMerge> documentation
    ("Analysis for jet merging scales");

  static Parameter<YMerge,unsigned int> interfaceNMax
    ("NMax",
     "Maximum n up to which y(n,n+1) should be booked",
     &YMerge::_nMax, 5, 2, 0,
     false, false, Interface::lowerlim);


  static Command<YMerge> interfaceY
    ("Y",
     "Insert options for given n.",
     &YMerge::Y, false);

}

string YMerge::Y (string in) {

  istringstream m (StringUtils::car(in));
  unsigned int mult; m >> mult;
  if (mult < 2 || mult > _nMax) {
    Throw<InitException> () << "YMerge : n must be between 2 and nMax";
  }

  if (_options.size() < _nMax+1) {
    _options.resize(_nMax+1,"");
  }

  _options[mult] = StringUtils::cdr(in);

  return "";

}

void YMerge::doinit() throw(InitException) {

  Analysis2Base::doinit();

  _yn.clear();

  if (!jetFinder())
    Throw<InitException>() << "YMerge : No JetFinder has been set, giving up.";
  if (_options.size() < 3)
    Throw<InitException>() << "YMerge : Need at least options for y_{2,3}";
  if (_options.size() < _nMax+1)
    Throw<InitException>() << "YMerge : It seems that no options where set.";
  
  int plotFlags = HistogramOutput::Xlog | HistogramOutput::Ylog
    | HistogramOutput::Frame | HistogramOutput::Errorbars;
  Histogram2Options options (plotFlags);

  _yn.resize(_nMax+1,"");

  for (unsigned int i=2; i<=_nMax; ++i) {
    ostringstream name("");
    name << "y_" << i << "_" << (i+1);
    _yn[i] = name.str();
    insert(name.str(),_options[i],options);
  }
}

void YMerge::analyze(const tPVector &) {

  pair<vector<Lorentz5Momentum>,double> ev;

  while (*eventExtractor() >> ev) {

    jetFinder()->use(ev.first);
    for (unsigned int i =2; i<=_nMax; ++i) {
      book(jetFinder()->getYMerge(i),_yn[i],ev.second);
    }

  }

}

void YMerge::dofinish() {

  for (unsigned int i =2; i<=_nMax; ++i) {
    finish(_yn[i]);
  }
  Analysis2Base::dofinish();

}
