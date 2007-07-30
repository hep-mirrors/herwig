// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Histogram class.
//

#include "Histogram.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"

using namespace Herwig;
NoPIOClassDescription<Histogram> Histogram::initHistogram;
// Definition of the static class description member.
void Histogram::Init() {

  static ClassDocumentation<Histogram> documentation
    ("The Histogram class implements a simple histogram include data"
     " points for comparision with experimental results.");

}

string Histogram::versionstring = "";

void Histogram::topdrawOutput(ostream & out,
			      bool frame,
			      bool errorbars,
			      bool xlog, bool ylog,
			      string colour,
			      string title,  string titlecase,
			      string left,   string leftcase,
			      string bottom, string bottomcase,
			      bool smooth) const
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
      if (versionstring != "") {
	out << "TITLE RIGHT \"" << versionstring << "\"\n";
	out << "CASE        \"\"\n";
      }
      if (errorbars) out << "SET ORDER X Y DX DY \n";
      else out << "SET ORDER X Y DX\n";
    }
  // scales
  if(xlog && frame) out << "SET SCALE X LOG " << endl;
  if(ylog && frame) out << "SET SCALE Y LOG " << endl;
  // set the x limits

  const unsigned int lastDataBinIndx = _bins.size()-2;
  if (xlog && frame) {
    out << "SET LIMITS X " << _bins[1].limit << " " 
	<< _bins[lastDataBinIndx+1].limit << endl;
  }
  // work out the y points
  vector<double> yout;
  double ymax=-9.8765e34,ymin=9.8765e34;
  unsigned int numPoints = _globalStats.numberOfPoints();
  if (numPoints == 0) ++numPoints;

  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix)
    {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      double value = 0.5*_prefactor*_bins[ix].contents / (delta*numPoints);
      yout.push_back(value);
      ymax=max(ymax, max(value, _bins[ix].data+_bins[ix].dataerror) );
      if(yout.back()>0.) ymin=min(ymin,value);
      if(_bins[ix].data>0) ymin=min(ymin,_bins[ix].data);
    }
  if (ymin > 1e34)  ymin = 1e-34;
  if (ymax < 1e-33) ymax = 1e-33;
  if (ymax < 10*ymin) ymin = 0.1*ymax;

  if (ylog && frame) {
    out << "SET LIMITS Y " << ymin << " " << ymax << endl;
  }
  // the histogram from the event generator
  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix)
    {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      out << _bins[ix].limit+delta << '\t' << yout[ix-1] << '\t' << delta;
      if (errorbars) {
	out << '\t' << 0.5*sqrt(_bins[ix].contentsSq)/(delta*numPoints);
      }
      out << '\n';
    }
  // N.B. in td smoothing only works for histograms with uniform binning.
  if(!smooth) {
      out << "HIST " << colour << endl;
  } else {
      out << "SMOOTH Y LEVEL 2 " << endl;
      out << "JOIN " << colour   << endl;
  }
  if (_havedata) {
    // the real experimental data
    for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix)
      {
	double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
	out << _bins[ix].limit+delta << '\t' << _bins[ix].data << '\t' << delta;
	if (errorbars) {
	  out  << '\t' << _bins[ix].dataerror;
	}
	out << '\n';
      }
    out << "PLOT " << endl;
  }
}

void Histogram::simpleOutput(ostream & out, bool errorbars) const {
  // simple ascii output (eg for gnuplot)
  // work out the y points
  vector<double> yout;
  unsigned int numPoints = _globalStats.numberOfPoints();  
  if (numPoints == 0) ++numPoints;

  const unsigned int lastDataBinIndx = _bins.size()-2;
  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
    double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
    double value = 0.5*_prefactor*_bins[ix].contents / (delta*numPoints);
    yout.push_back(value);
  }

  out << "# " << numPoints << " entries, mean +- sigma = "
      << _globalStats.mean() << " +- " 
      << _globalStats.stdDev() << "\n"
      << "# xlo xhi ynorm " 
      << (errorbars ? "ynorm_err " : "")
      << (_havedata ? "data " : "")
      << (_havedata && errorbars ? "dataerr " : "")
      << "y_entr\n";

  // the histogram from the event generator
  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
    double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
    out << _bins[ix].limit << " "
	<< _bins[ix+1].limit << " " 
	<< yout[ix-1];
    if (errorbars) {
      out << " " << 0.5*sqrt(_bins[ix].contentsSq)/(delta*numPoints);
    }
    if (_havedata) {
      out << " " << _bins[ix].data;
      if (errorbars)
	out << " " << _bins[ix].dataerror;
    }
    out << " " << _bins[ix].contents << '\n';
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
      if(_bins[ix].dataerror>0.)
	{
	  double var = sqr(_bins[ix].dataerror);
	  numer += _bins[ix].data * value/var;
	  denom += sqr(value)/var;
	}
    }
  _prefactor = numer/denom;
}

void Histogram::chiSquared(double & chisq, 
			   unsigned int & ndegrees, double minfrac) const
{
  chisq =0.;
  ndegrees=_bins.size()-2;
  unsigned int numPoints = _globalStats.numberOfPoints();
  for(unsigned int ix=1;ix<_bins.size()-1;++ix)
    {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      double value = 0.5*_prefactor*_bins[ix].contents / (delta*numPoints);
      double error = _bins[ix].dataerror;
      if(error>0.)
	{
	  if(error/_bins[ix].data < minfrac) 
	    error = minfrac*_bins[ix].data;
	  double var=sqr(error)
	    + _bins[ix].contentsSq * sqr(0.5*_prefactor / (delta*numPoints));
	  chisq += sqr(_bins[ix].data - value) / var;
	}
    }
}

void Histogram::normaliseToCrossSection() {
  unsigned int numPoints = _globalStats.numberOfPoints();
  if (numPoints == 0) ++numPoints;
  _prefactor=CurrentGenerator::current().eventHandler()->histogramScale()*
    numPoints/nanobarn;
}

void Histogram::topdrawOutputAverage(ostream & out,
			      bool frame,
			      bool errorbars,
			      bool xlog, bool ylog,
			      string colour,
			      string title,  string titlecase,
			      string left,   string leftcase,
			      string bottom, string bottomcase) const {
  // output the title info if needed
  if(frame) {
    out << "NEW FRAME\n";
    out << "SET FONT DUPLEX\n";
    out << "TITLE TOP \""    << title     << "\"\n";
    out << "CASE      \""    << titlecase << "\"\n";
    out << "TITLE LEFT \""   << left      << "\"\n";
    out << "CASE       \""   << leftcase  << "\"\n";
    out << "TITLE BOTTOM \"" << bottom     << "\"\n";
    out << "CASE        \""  << bottomcase << "\"\n";
    if (errorbars) out << "SET ORDER X Y DX DY \n";
    else out << "SET ORDER X Y DX\n";
  }
  // scales
  if(xlog && frame) out << "SET SCALE X LOG " << endl;
  if(ylog && frame) out << "SET SCALE Y LOG " << endl;
  // set the x limits
  
  const unsigned int lastDataBinIndx = _bins.size()-2;
  if (xlog && frame) {
    out << "SET LIMITS X " << _bins[1].limit << " " 
	<< _bins[lastDataBinIndx+1].limit << endl;
  }
  // work out the y points
  vector<double> yout;
  double ymax=-9.8765e34,ymin=9.8765e34;
  unsigned int numPoints = _globalStats.numberOfPoints();
  if (numPoints == 0) ++numPoints;

  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
    double value = _prefactor*_bins[ix].contents / _bins[ix].points;
    yout.push_back(value);
    ymax=max(ymax, max(value, _bins[ix].data+_bins[ix].dataerror) );
    if(yout.back()>0.) ymin=min(ymin,value);
    if(_bins[ix].data>0) ymin=min(ymin,_bins[ix].data);
  }
  if (ymin > 1e34)  ymin = 1e-34;
  if (ymax < 1e-33) ymax = 1e-33;
  if (ymax < 10*ymin) ymin = 0.1*ymax;
  
  if (ylog && frame) {
    out << "SET LIMITS Y " << ymin << " " << ymax << endl;
  }
  // the histogram from the event generator
  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
    double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
    out << _bins[ix].limit+delta << '\t' << yout[ix-1] << '\t' << delta;
    if (errorbars) {
      out << '\t' << 0.5*sqrt(_bins[ix].contentsSq)/(delta*numPoints);
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
	if (errorbars) {
	  out  << '\t' << _bins[ix].dataerror;
	}
	out << '\n';
      }
    out << "PLOT " << endl;
  }
}

