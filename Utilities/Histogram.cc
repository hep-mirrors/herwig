// -*- C++ -*-
//
// Histogram.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Histogram class.
//

#include "Histogram.h"
#include "HerwigStrategy.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"

using namespace Herwig;

IBPtr Histogram::clone() const {
  return new_ptr(*this);
}

IBPtr Histogram::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<Histogram> Histogram::initHistogram;
// Definition of the static class description member.
void Histogram::Init() {

  static ClassDocumentation<Histogram> documentation
    ("The Histogram class implements a simple histogram include data"
     " points for comparision with experimental results.");

}

void Histogram::rivetOutput(ostream & out,
                            string histogramname,
                            string analysisname,
                            string title,
                            string xlabel,
                            string ylabel,
                            bool rawcount,
                            double multiplicator) const {

  // total number of histogram entries
  double numPoints = static_cast<double>(_total);
  if (numPoints == 0) numPoints += 1.0;


  // collect y values and errors
  vector<double> yvalues;
  vector<double> yerrors;
  const unsigned int lastBinIndex = _bins.size() - 2;
  for (size_t ix = 1; ix <= lastBinIndex; ++ix) {
    const double delta = _bins[ix+1].limit - _bins[ix].limit;
    double factor = rawcount ? _prefactor * multiplicator
                             : _prefactor * multiplicator / (numPoints * delta);
    const double value = factor * _bins[ix].contents;
    const double error = factor * sqrt(_bins[ix].contentsSq);
    yvalues.push_back(value);
    yerrors.push_back(error);
  }

  // file header
  out << "## " << numPoints << " entries, mean +- sigma = "
      << _globalStats.mean() << " +- " 
      << _globalStats.stdDev() << "\n";
  out << "## xlo xhi y yerr\n" ;
  out << "##\n";
  out << "# BEGIN HISTOGRAM /" << analysisname << "/" << histogramname << "\n";
  out << "AidaPath=/" << analysisname << "/" << histogramname << "\n";
  if ( title  != string() ) out << "Title=" << title << "\n";
  if ( xlabel != string() ) out << "XLabel=" << xlabel << "\n";
  if ( ylabel != string() ) out << "YLabel=" << ylabel << "\n";

  // histogram data
  for (size_t i = 0; i < yvalues.size(); ++i) {
    out << _bins[i+1].limit << "\t"
	<< _bins[i+2].limit << "\t" 
	<< yvalues[i] << "\t"
	<< yerrors[i] << "\n";
  }

  // footer
  out << "# END HISTOGRAM\n";
  out << "\n";

}


void Histogram::topdrawOutput(ostream & out,
			      unsigned int flags,
			      string colour,
			      string title,  string titlecase,
			      string left,   string leftcase,
			      string bottom, string bottomcase
			      ) const {
  using namespace HistogramOptions;
  bool frame     = ( flags & Frame )     == Frame;
  bool errorbars = ( flags & Errorbars ) == Errorbars;
  bool xlog      = ( flags & Xlog )      == Xlog;
  bool ylog      = ( flags & Ylog )      == Ylog;
  bool smooth    = ( flags & Smooth )    == Smooth;
  bool rawcount  = ( flags & Rawcount )  == Rawcount;

  // output the title info if needed
  if(frame) {
    out << "NEW FRAME\n";
    if(_havedata) out << "SET WINDOW X 1.6 8 Y 3.5 9\n"; 
    else          out << "SET WINDOW X 1.6 8 Y 1.6 9\n";
    out << "SET FONT DUPLEX\n";
    out << "TITLE TOP \""    << title     << "\"\n";
    out << "CASE      \""    << titlecase << "\"\n";
    out << "TITLE LEFT \""   << left      << "\"\n";
    out << "CASE       \""   << leftcase  << "\"\n";
    out << (errorbars ? "SET ORDER X Y DX DY \n" : "SET ORDER X Y DX\n");
    if (HerwigStrategy::version != "") {
      out << "TITLE RIGHT \"" << HerwigStrategy::version << "\"\n";
      out << "CASE        \"\"\n";
    }
    if(_havedata) out << "SET AXIS BOTTOM OFF\n";
    else {
      out << "TITLE BOTTOM \"" << bottom     << "\"\n";
      out << "CASE        \""  << bottomcase << "\"\n";
    }
  }
  // scales
  if(xlog && frame) out << "SET SCALE X LOG " << endl;
  if(ylog && frame) out << "SET SCALE Y LOG " << endl;
  // set the x limits
  const unsigned int lastDataBinIndx = _bins.size()-2;
  if (frame) {
    out << "SET LIMITS X " << _bins[1].limit << " " 
	<< _bins[lastDataBinIndx+1].limit << endl;
  }
  // work out the y points
  vector<double> yout;
  double ymax=-9.8765e34,ymin=9.8765e34;
  double numPoints = _total;
  if (numPoints == 0) numPoints += 1.;

  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
    double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
    
    double factor = rawcount ? _prefactor : 0.5 * _prefactor / (numPoints * delta);

    double value = factor*_bins[ix].contents;
    yout.push_back(value);
    ymax=max(ymax, max(value, _bins[ix].data+_bins[ix].dataerror) );
    if(yout.back()>0.) ymin=min(ymin,value);
    if(_bins[ix].data>0) ymin=min(ymin,_bins[ix].data);
  }
  if (ymin > 1e34)  ymin = 1e-34;
  if (ymax < 1e-33) ymax = 1e-33;
  if (ymax < 10*ymin) ymin = 0.1*ymax;
  // make the y range slightly larger
  double fac=pow(ymax/ymin,0.1);
  ymax *= fac;
  ymin /= fac;

  if (ylog && frame) {
    out << "SET LIMITS Y " << ymin << " " << ymax << endl;
  }
  // the histogram from the event generator
  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
    double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
    double factor = rawcount ? _prefactor : 0.5 * _prefactor / (numPoints * delta);
    out << _bins[ix].limit+delta << '\t' << yout[ix-1] << '\t' << delta;
    if (errorbars) {
      out << '\t' << factor*sqrt(_bins[ix].contentsSq);
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
    for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      out << _bins[ix].limit+delta << '\t' << _bins[ix].data << '\t' << delta;
      if (errorbars) out  << '\t' << _bins[ix].dataerror;
      out << '\n';
    }
    out << "PLOT " << endl;
    out << "SET WINDOW X 1.6 8 Y 2.5 3.5\n";
    out << "SET LIMITS X " << _bins[1].limit << " " 
	<< _bins[lastDataBinIndx+1].limit << "\n";
    double ymax=0.;
    out << _bins[1].limit << "\t" << _bins[1].dataerror/_bins[1].data << "\n";
    for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      if(_bins[ix].data!=0.) {
	if(_bins[ix].dataerror/_bins[ix].data>ymax) 
	  ymax=_bins[ix].dataerror/_bins[ix].data;
	out << _bins[ix].limit+delta << '\t' 
	    <<  _bins[ix].dataerror/_bins[ix].data << '\n';
      }
      else {
	out << _bins[ix].limit+delta << '\t' 
	    <<  1. << '\n';
      }
    }
    if(_bins[lastDataBinIndx].data!=0.) {
      out << _bins[lastDataBinIndx+1].limit << "\t" 
	  << _bins[lastDataBinIndx].dataerror/_bins[lastDataBinIndx].data << "\n";
      out << _bins[lastDataBinIndx+1].limit << "\t" 
	  <<-_bins[lastDataBinIndx].dataerror/_bins[lastDataBinIndx].data << "\n";
    }
    else {
      out << _bins[lastDataBinIndx+1].limit << "\t" <<  1. << "\n";
      out << _bins[lastDataBinIndx+1].limit << "\t" << -1. << "\n";
    }
    for(unsigned int ix=lastDataBinIndx;ix>=1;--ix) {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      if(_bins[ix].data!=0.) {
	out << _bins[ix].limit+delta << '\t' 
	    <<  -_bins[ix].dataerror/_bins[ix].data << '\n';
      }
      else {
	out << _bins[ix].limit+delta << '\t' 
	    <<  -1. << '\n';
      }
    }
    if(_bins[1].data!=0.) {
      out << _bins[1].limit << "\t" << -_bins[1].dataerror/_bins[1].data << "\n";
    }
    else {
      out << _bins[1].limit << "\t" << -1. << "\n";
    }
    out << "set scale y lin\n";
    out << "set limits y " << -ymax << " " << ymax << "\n";
    out << "set fill full\n";
    out << "join yellow fill yellow\n";
    for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) { 
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      if(_bins[ix].data!=0.) {
	out << _bins[ix].limit+delta << "\t" 
	    << (yout[ix-1]-_bins[ix].data)/_bins[ix].data << "\n";
      }
      else if(_bins[ix].dataerror!=0.) {
	out << _bins[ix].limit+delta << "\t" 
	    << (yout[ix-1]-_bins[ix].data)/_bins[ix].dataerror << "\n";
      }
      else {
	out << _bins[ix].limit+delta << "\t" << 0. << "\n";
      }
    } 
    out << "join\n";
    out << "SET WINDOW X 1.6 8 Y 1.6 2.5\n";
    out << "SET LIMITS X " << _bins[1].limit << " " 
	<< _bins[lastDataBinIndx+1].limit << "\n";
    out << "SET AXIS BOTTOM ON\n";
    out << "TITLE BOTTOM \"" << bottom     << "\"\n";
    out << "CASE        \""  << bottomcase << "\"\n";
    ymax =0.;
    double ymin=0.;
    for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
      double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
      double error =  sqrt(sqr(0.5*sqrt(_bins[ix].contentsSq)/(delta*numPoints))+
			   sqr(_bins[ix].dataerror));
      double point=(yout[ix-1]-_bins[ix].data)/error;
      if(point<ymin) ymin=point;
      if(point>ymax) ymax=point;
      out << _bins[ix].limit+delta << '\t' 
	  << point  << '\n'; 
    }
    out << "set limits y " << ymin << " " << ymax << "\n";
    out << "JOIN" << endl;
  }
}

double Histogram::dataNorm() const {
  double norm(0.0);
  if (_havedata) {
    const unsigned int lastDataBinIndx = _bins.size()-2;
    for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
      double delta = _bins[ix+1].limit-_bins[ix].limit;
      double value = _bins[ix].data;
      norm += delta*value;
    }
  } else {
    norm = -1.0;
  }
  return norm;
}

unsigned int Histogram::visibleEntries() const {
  unsigned int numPoints(0);
  const unsigned int lastDataBinIndx = _bins.size()-2;
  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
    numPoints += static_cast<unsigned int>( _bins[ix].contents );
  }
  return numPoints;
}

void Histogram::simpleOutput(ostream & out, bool errorbars, 
			     bool normdata) {
  // simple ascii output (eg for gnuplot)
  // work out the y points
  vector<double> yout;
  unsigned int numPoints = visibleEntries();
  if (numPoints == 0) ++numPoints;
  double datanorm(1.0); 
  double chisq(0.0), minfrac(0.05);
  unsigned int ndof(0);
  if (_havedata) {
    if (normdata) datanorm = dataNorm();
    normaliseToData();
    chiSquared(chisq, ndof, minfrac);    
  }
  prefactor(1.0);

  const unsigned int lastDataBinIndx = _bins.size()-2;
  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
    double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
    double value = 0.5*_prefactor*_bins[ix].contents / (delta*numPoints);
    yout.push_back(value);
  }

  out << "# " << numPoints << " entries, mean +- sigma = "
      << _globalStats.mean() << " +- " 
      << _globalStats.stdDev() << "\n";
  if (_havedata) {
    out << "# chi^2/dof = " << chisq << "/" << ndof << " = " 
	<< chisq/double(ndof) << " (min err = " << minfrac << ")\n";
    if (datanorm) {
      out << "# data normalised by factor " << datanorm << "\n"; 
    }
  }
  out << "# xlo xhi ynorm " 
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
      out << " " << _bins[ix].data/datanorm;
      if (errorbars)
	out << " " << _bins[ix].dataerror/datanorm;
    }
    out << " " << _bins[ix].contents << '\n';
  }
}

vector<double> Histogram::dumpBins() const {
  vector<double> bincontents(_bins.size());
  for (size_t i=0; i < _bins.size(); ++i)
    bincontents[i] = _bins[i].contents;
  return bincontents;
}

Histogram Histogram::ratioWith(const Histogram & h2) const {
  const size_t numBins = _bins.size();
  assert( numBins > 2 && numBins == h2._bins.size());
  Histogram ratio(*this);
  for (size_t i=0; i < numBins; ++i) {
    assert(_bins[i].limit == h2._bins[i].limit);
    if (h2._bins[i].contents > 0.0)
      ratio._bins[i].contents /= h2._bins[i].contents;
    else
      ratio._bins[i].contents = 0.0;
  }
  return ratio;
}


void Histogram::normaliseToData() {
  double numer(0.),denom(0.);
  double numPoints = _total;
  for(unsigned int ix=1;ix<_bins.size()-1;++ix) {
    double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
    double value = 0.5*_bins[ix].contents / (delta*numPoints);
    if(_bins[ix].dataerror>0.) {
      double var = sqr(_bins[ix].dataerror);
      numer += _bins[ix].data * value/var;
      denom += sqr(value)/var;
    }
  }
  _prefactor = numer/denom;
}

void Histogram::chiSquared(double & chisq, 
			   unsigned int & ndegrees, double minfrac) const {
  chisq =0.;
  ndegrees=0;
  double numPoints = _total;
  for(unsigned int ix=1;ix<_bins.size()-1;++ix) {
    double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);
    double value = 0.5*_prefactor*_bins[ix].contents / (delta*numPoints);
    double error = _bins[ix].dataerror;
    if(error>0.) {
      if(abs(error/_bins[ix].data) < minfrac) error = minfrac*_bins[ix].data;
      double var=sqr(error)
	+ _bins[ix].contentsSq * sqr(0.5*_prefactor / (delta*numPoints));
      chisq += sqr(_bins[ix].data - value) / var;
      ++ndegrees;
    }
  }
}

void Histogram::normaliseToCrossSection() {
  double numPoints = _total;
  if (numPoints == 0) numPoints += 1.;
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


vector<double> Histogram::LogBins(double xmin, unsigned nbins, double base) {
  vector<double> limits;
  for (unsigned e = 0; e <= nbins; e++)
    limits.push_back(xmin * pow(base, (int)e));
  return limits;
}


void Histogram::topdrawMCatNLO(ostream & out,
			       unsigned int flags,
			       string,
			       string title
			       ) const {
  
  using namespace HistogramOptions;
  //  bool frame     = ( flags & Frame )     == Frame;
  bool errorbars = ( flags & Errorbars ) == Errorbars;
  //  bool xlog      = ( flags & Xlog )      == Xlog;
  bool ylog      = ( flags & Ylog )      == Ylog;
  //  bool smooth    = ( flags & Smooth )    == Smooth;
  //  bool rawcount  = ( flags & Rawcount )  == Rawcount;

  double myFactor = _prefactor / _total * 1000.;


  // output the title info if needed
  out << "     ( 22-Apr-10 18:28\n\n";
  out << "   NEW PLOT\n\n\n";
  out << " ( SET FONT DUPLEX\n";
  out << "  SET TITLE SIZE 2\n";
  out << " TITLE 12.8 9 ANGLE -90 \" MLM   22-Apr-10 18:28\"\n";

  out << "  ( SET FONT DUPLEX\n";
  out << "  SET TITLE SIZE  -1.2247\n";
  out << "  SET LABEL SIZE  -1.2247\n";
  out << "  SET TICKS TOP OFF SIZE   0.0245\n";
  out << "  SET WINDOW X   1.5000 TO   12.000\n";
  out << "  SET WINDOW Y   1.0000 TO   8.7917\n";
  out << "  TITLE   1.5000   8.9617 \" "<<title<<"\"\n";
  out << "  TITLE   9.8719   8.6217 \" INT= "<<visibleEntries()*myFactor<<"\"\n";
  out << "  TITLE   9.8719   8.4517 \" ENT= "<<visibleEntries()<<"\"\n";
  out << "  TITLE   9.8719   8.2817 \" OFL= 2.258E+01\"\n";
  out << "  SET ORD X Y \n";
  out << "  9.8719   8.1117\n";
  out << "  12.000   8.1117\n";
  out << "  JOIN TEXT\n";
  out << "    9.8719   8.1117\n";
  out << "    9.8719   8.7917\n";
  out << "  JOIN TEXT\n";
  out << "  SET TITLE SIZE  -1.8371\n";
  out << " TITLE BOTTOM \""<<title<<"\"\n";
  out << "  TITLE    0.42188   7.37500 ANGLE 90 \" \"\n";
  if (ylog) {
    out << "  SET SCALE Y LOG\n"; }
  else {
    out << "  SET SCALE Y LIN\n"; }
  out << "  SET TICKS TOP OFF\n";

  // set the x limits
  const unsigned int lastDataBinIndx = _bins.size()-2;

  out << "  SET LIMITS X " << _bins[1].limit << " " 
      << _bins[lastDataBinIndx+1].limit << endl;
  // work out the y points
  vector<double> yout;
  double ymax=-9.8765e34,ymin=9.8765e34;
  double numPoints = _total;
  if (numPoints == 0) numPoints += 1.;

  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
    // no differential dist & in pb
    double factor = _prefactor / numPoints *1000.;
    double value = factor*_bins[ix].contents;

    yout.push_back(value);
    ymax=max(ymax, max(value, _bins[ix].data+_bins[ix].dataerror) );
    if(yout.back()>0.) ymin=min(ymin,value);
    if(_bins[ix].data>0) ymin=min(ymin,_bins[ix].data);
  }
  if (ymin > 1e34)  ymin = 1e-34;
  if (ymax < 1e-33) ymax = 1e-33;
  if (ymax < 10*ymin) ymin = 0.1*ymax;
  // make the y range slightly larger
  double fac=pow(ymax/ymin,0.1);
  ymax *= fac;
  ymin /= fac;
  //if (ylog) {
  //  out << "  SET LIMITS Y " << ymin << " " << ymax << endl;
  //}  
  out << "  SET LIMITS Y " << ymin << " " << ymax << endl;
  out << "  SET ORDER X Y DY\n";
  out << " (  "<<title<<"\n";

  out << " ( INT= "<<visibleEntries()*myFactor<<"  ENTRIES=  "<<visibleEntries()*myFactor<<"\n";

  // the histogram from the event generator
  for(unsigned int ix=1; ix<=lastDataBinIndx; ++ix) {
    double delta = 0.5*(_bins[ix+1].limit-_bins[ix].limit);

    //double factor = rawcount ? _prefactor : 0.5 * _prefactor / (numPoints * delta);

    // no differential distributions and in pb
    double factor = _prefactor / numPoints * 1000.;
    
    out << "    "<<_bins[ix].limit+delta << "    " << yout[ix-1] << "    " << delta;
    if (errorbars) {
      out << "    " << factor*sqrt(_bins[ix].contentsSq);
    }
    out << '\n';
  }
  out << "  HIST SOLID\n";
  out << "   PLOT\n"; 
}
