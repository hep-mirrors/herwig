// -*- C++ -*-
//
// Histogram.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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

// workaround for OS X bug where isnan() and isinf() are hidden
// when <iostream> is included
extern "C" int isnan(double) throw();

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
  inline Histogram(double lower=0., double upper=0., unsigned int nbin=0);

  /**
   * Constructor for variable width bins
   * @param limits The lower limits for the bins followed by the upper limit of the last bin
   */
  inline Histogram(vector<double> limits);

  /**
   * Constructor with data included
   * @param limits The lower limits for the bins followed by the upper limit of the last bin
   * @param data The data
   * @param dataerror The errors on the data
   */
  inline Histogram(vector<double> limits, vector<double> data, vector<double> dataerror);
  //@}

public:

  /**
   *  Operator to add a point to the histogrma with unit weight
   */
  inline void operator+=(double);

  /**
   *  Function to add a weighted point to the histogram
   */
  inline void addWeighted(double data, double weight);

  /**
   *  Number of bins (not counting the overflow)
   */
  inline unsigned int numberOfBins() const;

  /**
   *  Get the prefactor
   */
  inline double prefactor() const;

  /**
   *  Set the prefactor
   */
  inline void   prefactor(double );

  /**
   *  Access to the statistics on the total entry of the histogram
   */
  inline const Statistic & globalStatistics() const; 

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<Histogram> initHistogram;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Histogram & operator=(const Histogram &);

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
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Histogram. */
template <>
struct BaseClassTrait<Herwig::Histogram,1> {
  /** Typedef of the first base class of Histogram. */
  typedef Herwig::Statistic NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Histogram class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Histogram>
  : public ClassTraitsBase<Herwig::Histogram> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::Histogram"; }
};

/** @endcond */

}

#include "Histogram.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Histogram.tcc"
#endif


// void SampleHistogram::printMoments(char* name, double Nmax, double dN, 
// 				   double x0, double x1) {
//   ofstream out(name);
//   if (!out) {
//     cerr << "SampleHistoGram::printMoments: ERROR! Can't open file" << endl;
//   }

//   time_t now_t;
//   now_t = time(0);
//   out << "# created " << ctime(&now_t)
//       << "# by SampleHistogram::printMoments(..., "
//       << Nmax << ", " << dN << ")" << endl
//       << "# " << this->samples() << " entries, mean +/- sigma = " 
//       << this->mean() << " +/- " << this->stdDev() << endl;

//   double x0N, x1N, delta, hi;
//   for (double N=dN; N < Nmax; N += dN) {
//     double fN = 0.0;
//     for(int i = 0; i < howManyBuckets-1; i++) {
//       x0N = pow(bucketLimit[i], N);
//       x1N = pow(bucketLimit[i+1], N);
//       delta = (bucketLimit[i+1] - bucketLimit[i]);
//       if (delta > 0 && this->samples() > 0 
// 	  && bucketLimit[i] >= x0 && bucketLimit[i] <= x1
// 	  && bucketLimit[i+1] >= x0 && bucketLimit[i+1] <= x1) {
// 	hi = double(bucketCount[i+1]/(delta*(this->samples())));
// 	fN += hi*(x1N-x0N)/N;
//       }
//     }
//     out << N 
// 	<< "  " 
// 	<< fN << endl;
//   }
//   out.close();
// }


#endif /* HERWIG_Histogram_H */
