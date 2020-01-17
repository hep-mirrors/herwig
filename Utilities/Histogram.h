// -*- C++ -*-
//
// Histogram.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Histogram_H
#define HERWIG_Histogram_H
//
// This is the declaration of the Histogram class.
//
#include "Histogram.fh"
#include "ThePEG/Interface/Interfaced.h"
#include "Statistic.h"
#include <string>

namespace Herwig {

using namespace ThePEG;

  /**
   * Options for histogram output. 
   * They can be combined using the '|' operator, e.g. 'Frame | Ylog'
   */
  namespace HistogramOptions {
    const unsigned int None      = 0;      /**< No options */
    const unsigned int Frame     = 1;      /**< Plot on new frame */
    const unsigned int Errorbars = 1 << 1; /**< Plot error bars */
    const unsigned int Xlog      = 1 << 2; /**< log scale for x-axis */
    const unsigned int Ylog      = 1 << 3; /**< log scale for y-axis */
    const unsigned int Smooth    = 1 << 4; /**< smooth the line */
    const unsigned int Rawcount  = 1 << 5; /**< don't normalize to unit area */
  }

/**
 * The Histogram class is a simple histogram for the Analysis handlers.
 *
 * @see \ref HistogramInterfaces "The interfaces"
 * defined for Histogram.
 */
class Histogram: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   * @param lower The lower limit of the histogram
   * @param upper The upper limit of the histogram
   * @param nbin  Number of bins
   */
  Histogram(double lower=0., double upper=0., unsigned int nbin=0)
  : _globalStats(), _havedata(false), _bins(nbin+2),_prefactor(1.),_total(0.) {
    if (upper<lower) swap(upper,lower);
    _bins[0].limit=-1.e100;
    double limit(lower);
    double width((upper-lower)/nbin);
    for(unsigned int ix=1; ix <= nbin; ++ix) {
      _bins[ix].limit=limit;
      limit += width;
    }
    _bins.back().limit=limit;
  }

  /**
   * Constructor for variable width bins
   * @param limits The lower limits for the bins followed by the upper limit of the last bin
   */
  Histogram(vector<double> limits) 
    : _globalStats(), _havedata(false), _bins(limits.size()+1), _prefactor(1.),_total(0.) {
    _bins[0].limit=-1.e100;
    for (size_t i=1; i<=limits.size(); ++i)
      _bins[i].limit=limits[i-1];
  }

  /**
   * Constructor with data included
   * @param limits The lower limits for the bins followed by the upper limit of the last bin
   * @param data The data
   * @param dataerror The errors on the data
   */
  Histogram(vector<double> limits, vector<double> data, vector<double> dataerror)
    : _globalStats(), _havedata(true), _bins(limits.size()+1), _prefactor(1.),_total(0.) {
    _bins[0].limit=-1.e100;
    for (size_t i=1; i<=limits.size(); ++i)
      _bins[i].limit=limits[i-1];
    
    // no data goes into _bins[0] or _bins.back()!
    for (size_t i=1; i<=min(limits.size()-1,data.size()); ++i)
      _bins[i].data=data[i-1];
    
    for (size_t i=1; i<=min(limits.size()-1,dataerror.size()); ++i)
      _bins[i].dataerror=dataerror[i-1];
  }
  
  //@}

public:

  /**
   *  Operator to add a point to the histogrma with unit weight
   */
  void operator += (double input) {
    addWeighted(input,1.0);
  }

  /**
   *  Function to add a weighted point to the histogram
   */
  void addWeighted(double input, double weight) {
    if(std::isnan(input)) return;
    unsigned int ibin;
    for(ibin=1; ibin<_bins.size(); ++ibin) {
      if(input<_bins[ibin].limit)
	break;
    }
    _bins[ibin-1].contents   += weight;
    _bins[ibin-1].contentsSq += sqr(weight);
    _globalStats += weight * input;
    _total += weight;
  }

  /**
   *  Number of bins (not counting the overflow)
   */
  unsigned int numberOfBins() const { 
    return _bins.size()-2;
  }


  /**
   *  Get the prefactor
   */
  double prefactor() const {
    return _prefactor;
  }

  /**
   *  Set the prefactor
   */
  void   prefactor(double in ) {
    _prefactor=in;
  }

  /**
   *  Access to the statistics on the total entry of the histogram
   */
  const Statistic & globalStatistics() const {
    return _globalStats;
  }

  /**
   *  Normalise the distributions to the data
   */
  void normaliseToData();

  /**
   *  Normalise the distributions to the total cross section
   */
  void normaliseToCrossSection();

  /**
   *  Return the chi squared
   * @param chisq The chi squared
   * @param ndegrees The number of points
   * @param minfrac The minimum fractional error on the data point
   */
  void chiSquared(double & chisq, 
		  unsigned int & ndegrees, double minfrac=0.) const;

  /**
   * @brief Output as file ready for usage with flat2aida and other Rivet tools
   * @param out           The output stream
   * @param histogramname The histogram name identifying the histogram. Required
   *                      for comparisons (e.g. with rivet-mkhtml or with
   *                      compare-histos)
   * @param analysisname  The analysis name
   * @param title         The title for the top of the plot in LaTeX format
   * @param xlabel        The x label in LaTeX format
   * @param ylabel        The y label in LaTeX format
   * @param rawcount      Don't normalise to unit area.
   * @param multiplicator Factor the histogram is multiplied with.
   * N.B.  Experimental data is not output.
   */
  void rivetOutput(ostream & out,
                   string histogramname = string("default"),
                   string analysisname = string("default"),
                   string title  = string(),
                   string xlabel = string(),
                   string ylabel = string(),
                   bool rawcount = false,
                   double multiplicator = 1.0) const;

  /**
   *  Output as a topdrawer file. The histogram is normalised to unit area
   * @param out The output stream
   * @param flags A bitmask of flags from HistogramOptions, e.g. Frame|Ylog
   * @param colour The colour for the line
   * @param title  The title for the top of the plot
   * @param titlecase topdraw format for the title
   * @param left   Left axis lable
   * @param leftcase topdraw format for left axis label
   * @param bottom  Bottom axis lable
   * @param bottomcase Bottom axis lable ofr topdraw
   * N.B. in td smoothing only works for histograms with uniform binning.
   */
  void topdrawOutput(ostream & out,
		     unsigned int flags = 0,
		     string colour = string("BLACK"),
		     string title = string(),
		     string titlecase = string(),
		     string left = string(),
		     string leftcase = string(),
		     string bottom = string(),
		     string bottomcase = string()
		     ) const;

  void topdrawMCatNLO(ostream & out,
		      unsigned int flags =0 ,
		      string colour = string("BLACK"),
		      string title = string()
		      ) const;

  /**
   *  Output as a topdrawer file. A bin by bin average is taken.
   * @param out The output stream
   * @param frame output on a new graph
   * @param errorbars output data points with error bars
   * @param xlog  log scale on x axis
   * @param ylog  log scale on y axis
   * @param colour The colour for the line
   * @param title  The title for the top of the plot
   * @param titlecase topdraw format for the title
   * @param left   Left axis lable
   * @param leftcase topdraw format for left axis label
   * @param bottom  Bottom axis lable
   * @param bottomcase Bottom axis lable ofr topdraw
   */
  void topdrawOutputAverage(ostream & out,
			    bool frame,
			    bool errorbars,
			    bool xlog, bool ylog,
			    string colour=string("BLACK"),
			    string title=string(),
			    string titlecase =string(),
			    string left=string(),
			    string leftcase =string(),
			    string bottom=string(),
			    string bottomcase =string()) const;

  /**
   * get the number of visible entries (all entries without those in the
   * under- and overflow bins) in the histogram.  This assumes integer
   * entries, ie it gives wrong results for weighted histograms.
   */
  unsigned int visibleEntries() const;

  /**
   * Compute the normalisation of the data. 
   */
  double dataNorm() const;

  /**
   * Output into a simple ascii file, easily readable by gnuplot.
   */
  void simpleOutput(ostream & out, bool errorbars, bool normdata=false);

  /**
   * Dump bin data into a vector
   */
  vector<double> dumpBins() const;

  /**
   * Returns a new histogram containing bin-by-bin ratios of two histograms
   */
  Histogram ratioWith(const Histogram & h2) const;


  /**
   * @brief Returns limits for bins with exponentially increasing widths.
   *        For usage with the variable-bin-width Histogram constructor.
   * @param xmin  Lower limit of the first bin, needs to be > 0
   * @param nbins Number of bins
   * @param base  The base, needs to be > 1
   */
  static
  vector<double> LogBins(double xmin, unsigned nbins, double base = 10.0);


public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Histogram & operator=(const Histogram &) = delete;

private:

 /**
   *  Global statistics of all data that went into the histogram.
   */
  Statistic _globalStats;

 /**
   * Set to true if there is experimental data available
   */
  bool _havedata;

 /**
   *  One bin of the histogram. limit is the _lower_ bound of the bin.
   */
  struct Bin {
    /**
     *  Default constructor
     */
    Bin() : contents(0.0), contentsSq(0.0), 
	    limit(0.0), data(0.0), dataerror(0.0), points(0) {}
    /**
     *  Contents of the bin
     */
    double contents;

    /**
     *  Contents squared for the error
     */
    double contentsSq;

    /**
     * The limit for the bin
     */
    double limit;

    /**
     *  The experimental value for the bin
     */
    double data;

    /**
     *  The error on the experimental value for the bin
     */
    double dataerror;

    /**
     *  The number of points in the bin
     */
    long points;
  };

  /**
   *  The histogram bins. _bins[0] is the underflow, _bins.back() the overflow
   */
  vector<Bin> _bins;

  /**
   *  Prefactors to multiply the output by
   */
  double _prefactor;

  /**
   *  Total entry
   */
  double _total;


public:

  /**
   * The vector of bins
   */
  vector<Bin> bins() const { return _bins; }

};

}

#endif /* HERWIG_Histogram_H */
