// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#include "Histogram2Output.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Histogram2Output.tcc"
#endif

using namespace Analysis2;

Histogram2Output::~Histogram2Output() {}

void Histogram2Output::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _prefix << _mctitle;
}

void Histogram2Output::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _prefix >> _mctitle;
}

ClassDescription<Histogram2Output> Histogram2Output::initHistogram2Output;
// Definition of the static class description member.

void Histogram2Output::Init() {

  static ClassDocumentation<Histogram2Output> documentation
    ("Histogram2Output is the base class for all output options"
     " for Histogram2 objects.");


  static Parameter<Histogram2Output,string> interfacePrefix
    ("Prefix",
     "Set the output file or prefix.",
     &Histogram2Output::_prefix, "",
     false, false);

  static Parameter<Histogram2Output,string> interfaceMcTitle
    ("MCModel",
     "Set the name for the MC model used.",
     &Histogram2Output::_mctitle, "",
     false, false);

}


void Histogram2Output::initialize (const string& name) {
  if (_out.is_open()) _out.close();
  _out.open((prefix()+name+".dat").c_str());
}


void Histogram2Output::put (Histogram2Ptr histo, const Histogram2Options& options, const string&) {

  initialize(histo->name());

  currentOStream () << "# The following histograms have been produced by\n"
		    << "# by Analysis2 " << "\n"
		    << "# unless a channel is explictly marked as data.\n"
		    << "# A model used is refered to " << _mctitle << "\n"
		    << endl;

  vector<string> channels = histo->channels();
  
  for (vector<string>::iterator c = channels.begin();
       c != channels.end(); ++c) {
    histo->output(currentOStream(),*c, options.channelFlags);
    currentOStream() << endl << endl;
  }

}

ostream& Histogram2Output::currentOStream () {
  return _out;
}

