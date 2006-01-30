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
      if (error) out << "SET ORDER X Y DX DY \n";
      else out << "SET ORDER X Y DX\n";
    }
  // scales
  if(xlog) out << "SET SCALE X LOG " << endl;
  if(ylog) out << "SET SCALE Y LOG " << endl;
  // set the x limits

  const unsigned int lastDataBinIndx = _bins.size()-2;

  out << "SET LIMITS X " << _bins[1].limit << " " 
      << _bins[lastDataBinIndx+1].limit << endl;
  // work out the y points
  vector<double> yout;
  double ymax=-1e100,ymin=1e100;
  unsigned int numPoints = _globalStats.numberOfPoints();
  if (numPoints == 0) ++numPoints;

  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix)
    {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      double value = 0.5*_prefactor*_bins[ix].contents / (delta*numPoints);
      yout.push_back(value);
      ymax=max(ymax, max(value, _bins[ix].data+_bins[ix].error) );
      if(yout.back()>0.) ymin=min(ymin,value);
      if(_bins[ix].data>0) ymin=min(ymin,_bins[ix].data);
    }

  out << "SET LIMITS Y " << (ylog ? ymin:0. ) << " " << ymax << endl;

  // the histogram from the event generator
  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix)
    {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      out << _bins[ix].limit+delta << '\t' << yout[ix-1] << '\t' << delta;
      if (error) {
	out << '\t' << 0.5*sqrt(_bins[ix].contents)/(delta*numPoints);
      }
      out << '\n';
    }
  out << "HIST " << colour << endl;

  if (_havedata) {
    // the real experimental data
    for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix)
      {
	double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
	out << _bins[ix].limit+delta << '\t' << _bins[ix].data << '\t' << delta;
	if (error) {
	  out  << '\t' << _bins[ix].error;
	}
	out << '\n';
      }
    out << "PLOT " << endl;
  }
}

void Histogram::normaliseToData()
{
  double numer(0.),denom(0.);
  unsigned int numPoints = _globalStats.numberOfPoints();
  for(unsigned int ix=1;ix<_bins.size()-1;++ix)
    {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      double value = 0.5*_bins[ix].contents / (delta*numPoints);
      if(_bins[ix].error>0.)
	{
	  double var=sqr(_bins[ix].error);
	  numer += _bins[ix].data*value/var;
	  denom += sqr(value)/var;
	}
    }
  _prefactor=numer/denom;
}

void Histogram::chiSquared(double & chisq, unsigned int & ndegrees, double minfrac)
{
  chisq =0.;
  ndegrees=_bins.size()-2;
  unsigned int numPoints = _globalStats.numberOfPoints();
  for(unsigned int ix=1;ix<_bins.size()-1;++ix)
    {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      double value = 0.5*_prefactor*_bins[ix].contents / (delta*numPoints);
      double error=_bins[ix].error;
      if(error>0.)
	{
	  if(error/_bins[ix].data<minfrac) error=minfrac*_bins[ix].data;
	  double var=sqr(error)
	    + _bins[ix].contents*sqr(0.5*_prefactor / (delta*numPoints));
	  chisq += sqr(_bins[ix].data-value)/var;
	}
    }
}

