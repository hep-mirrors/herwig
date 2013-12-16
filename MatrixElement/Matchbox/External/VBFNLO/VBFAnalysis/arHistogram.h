// -*- C++ -*-
//
// arHistogram.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ARNOLD_arHistogram_H
#define ARNOLD_arHistogram_H
//
// This is the declaration of the arHistogram class.
//
#include "arHistogram.fh"
#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Utilities/Statistic.h"
#include <string>

// workaround for OS X bug where isnan() and isinf() are hidden
// when <iostream> is included
extern "C" int isnan(double) throw();

namespace arnold {

using namespace Herwig;

  /**
   * Options for histogram output. 
   * They can be combined using the '|' operator, e.g. 'Frame | Ylog'
   */
  namespace arHistogramOptions {
    const unsigned int None      = 0;      /**< No options */
    const unsigned int Frame     = 1;      /**< Plot on new frame */
    const unsigned int Errorbars = 1 << 1; /**< Plot error bars */
    const unsigned int Xlog      = 1 << 2; /**< log scale for x-axis */
    const unsigned int Ylog      = 1 << 3; /**< log scale for y-axis */
    const unsigned int Smooth    = 1 << 4; /**< smooth the line */
    const unsigned int Rawcount  = 1 << 5; /**< don't normalize to unit area */
  }

/**
 * The arHistogram class is a simple histogram for the Analysis handlers.
 *
 * @see \ref arHistogramInterfaces "The interfaces"
 * defined for arHistogram.
 */
class arHistogram: public Interfaced {
  friend class ar2DHistogram;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   * @param lower The lower limit of the histogram
   * @param upper The upper limit of the histogram
   * @param nbin  Number of bins
   * @param nbin  Smear over several bins
   * @param nbin  Set if variable is cyclic (e.g. phi_jj). Only relevant for smearing.
   */
  inline arHistogram(double lower=0., double upper=0., unsigned int nbin=0, bool smearing=false, bool cyclic=false);

  /**
   * Constructor for variable width bins
   * @param limits The lower limits for the bins followed by the upper limit of the last bin
   */
  inline arHistogram(vector<double> limits);

  /**
   * Constructor with data included
   * @param limits The lower limits for the bins followed by the upper limit of the last bin
   * @param data The data
   * @param dataerror The errors on the data
   */
  inline arHistogram(vector<double> limits, vector<double> data, vector<double> dataerror);
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
   * @param flags A bitmask of flags from arHistogramOptions, e.g. Frame|Ylog
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
  void simpleOutput(ostream & out, bool errorbars, bool normdata=false, bool print_dsigma_over_dx = false);

  /**
   * Dump bin data into a vector
   */
  vector<double> dumpBins() const;

  /**
   * Returns a new histogram containing bin-by-bin ratios of two histograms
   */
  arHistogram ratioWith(const arHistogram & h2) const;


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
  static NoPIOClassDescription<arHistogram> initarHistogram;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  arHistogram & operator=(const arHistogram &);

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
	    limit(0.0), data(0.0), dataerror(0.0) {}
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

  /**
   *  Smear input over two bins
   */
  bool doSmearing;

  /**
   *  For activated smearing: Is the histogram depicting a cyclic variable?
   */
  bool isCyclic;
};

/**
 * The ar2D Histogram class is a simple 2D histogram for the Analysis handlers.
 *
 * @see \ref ar2DHistogramInterfaces "The interfaces"
 * defined for ar2DHistogram.
 */
class ar2DHistogram: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   * @param lower The lower limit of the histogram
   * @param upper The upper limit of the histogram
   * @param nbin  Number of bins
   * @param nbin  Smear over several bins
   * @param nbin  Set if variable is cyclic (e.g. phi_jj). Only relevant for smearing.
   */
  inline ar2DHistogram(double lowerx=0., double upperx=0., double lowery=0., double uppery=0., unsigned int nbinx=0, unsigned int nbiny=0, bool smearing=false, bool cyclic=false);

  /**
   *  Function to add a weighted point to the histogram
   */
  inline void addWeighted(double datax, double datay, double weight);

  /**
   * Output into a simple ascii file, easily readable by gnuplot.
   */
  void simpleOutput(ostream & out, bool errorbars, bool normdata=false, bool print_dsigma_over_dx = false);

  /**
   *  Normalise the distributions to the data
   */
  void normaliseToData();

  /**
   *  Return the chi squared
   * @param chisq The chi squared
   * @param ndegrees The number of points
   * @param minfrac The minimum fractional error on the data point
   */
  void chiSquared(double & chisq, 
		  unsigned int & ndegrees, double minfrac=0.) const;
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

  vector<arHistogramPtr> histo;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<ar2DHistogram> initar2DHistogram;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ar2DHistogram & operator=(const ar2DHistogram &);

  double theLowerX, theUpperX, theLowerY, theUpperY;
  unsigned int thNBinX, theNBinY;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of arHistogram. */
template <>
struct BaseClassTrait<arnold::arHistogram,1> {
  /** Typedef of the first base class of arHistogram. */
  typedef Herwig::Statistic NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the arHistogram class and the shared object where it is defined. */
template <>
struct ClassTraits<arnold::arHistogram>
  : public ClassTraitsBase<arnold::arHistogram> {
  /** Return a platform-independent class name */
  static string className() { return "arnold::arHistogram"; }
};

/** This template specialization informs ThePEG about the
 *  base classes of ar2DHistogram. */
template <>
struct BaseClassTrait<arnold::ar2DHistogram,1> {
  /** Typedef of the first base class of ar2DHistogram. */
  typedef Herwig::Statistic NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ar2DHistogram class and the shared object where it is defined. */
template <>
struct ClassTraits<arnold::ar2DHistogram>
  : public ClassTraitsBase<arnold::ar2DHistogram> {
  /** Return a platform-independent class name */
  static string className() { return "arnold::ar2DHistogram"; }
};

/** @endcond */

}

#include "arHistogram.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "arHistogram.tcc"
#endif


#endif /* ARNOLD_arHistogram_H */
