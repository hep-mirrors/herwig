// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Histogram class.
//

#include "Histogram.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Histogram.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

Histogram::~Histogram() {}

void Histogram::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void Histogram::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<Histogram> Histogram::initHistogram;
// Definition of the static class description member.

void Histogram::Init() {

  static ClassDocumentation<Histogram> documentation
    ("There is no documentation for the Histogram class");

}

void Histogram::topdrawOutput(ofstream & out,
			      bool frame,
			      bool error,
			      bool xlog, bool ylog,
			      string colour,
			      string title,  string titlecase,
			      string left,   string leftcase,
			      string bottom, string bottomcase)
{
  // output the title info if needed
  if(frame)
    {
      out << "NEW FRAME\n";
      out << "SET FONT DUPLEX\n";
      out << "TITLE TOP \""    << title     << "\"\n";
      out << "CASE      \""    << titlecase << "\"\n";
      out << "TITLE LEFT \""   << left      << "\"\n";
      out << "CASE       \""   << leftcase  << "\"\n";
      out << "TITLE BOTTOM \"" << bottom     << "\"\n";
      out << "CASE        \""  << bottomcase << "\"\n";
      out << "SET ORDER X Y DX DY \n";
    }
  // scales
  if(xlog) out << "SET SCALE X LOG " << endl;
  if(ylog) out << "SET SCALE Y LOG " << endl;
  // set the x limits
  out << "SET LIMITS X " << _binlimits[0] << " " << _binlimits[_nbin-2] << endl;
  // work out the y points
  vector<double> yout;
  double ymax=-1e100,delta;
  for(unsigned int ix=0;ix<_nbin-2;++ix)
    {
      delta = 0.5*(_binlimits[ix+1]-_binlimits[ix]);
      yout.push_back(0.5*_bincontents[ix+1]->total()/(delta*numberOfPoints()));
      ymax=max(ymax,yout.back());
    }
  for(unsigned int ix=1;ix<_data.size();++ix)
    {ymax=max(ymax,_data[ix]+_error[ix]);}
  out << "SET LIMITS Y 0. " << ymax << endl;
  // the histogram from the event generator
  for(unsigned int ix=0;ix<_nbin-2;++ix)
    {
      delta = 0.5*(_binlimits[ix+1]-_binlimits[ix]);
      out << _binlimits[ix]+delta << "\t"
	  << yout[ix] << "\t"
	  << delta << "\t"
	  << 0.5*_bincontents[ix+1]->stdDev()/(delta*numberOfPoints()) << "\n";
    }
  out << "HIST " << colour << endl;
  // the real experimental data
  for(unsigned int ix=0;ix<_nbin-2;++ix)
    {
      delta = 0.5*(_binlimits[ix+1]-_binlimits[ix]);
      out << _binlimits[ix]+delta << "\t"
	  << _data[ix]            << "\t"
	  << delta                << "\t"
	  << _error[ix]           << "\n";
    }
  out << "PLOT " << endl;
}
