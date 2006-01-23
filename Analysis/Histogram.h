// -*- C++ -*-
#ifndef HERWIG_Histogram_H
#define HERWIG_Histogram_H
//
// This is the declaration of the Histogram class.
//

#include "Statistic.h"
#include "Histogram.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The Histogram class is a simple histogram for the Analysis handlers.
 *
 * @see \ref HistogramInterfaces "The interfaces"
 * defined for Histogram.
 */
class Histogram: public Statistic {

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
   * @param limits The limits for the bins.
   */
  inline Histogram(vector<double> limits);

  /**
   * Constructor with data included
   * @param limits The limits for the bins.
   * @param data The data
   * @param error The errors on the data
   */
  inline Histogram(vector<double> limits, vector<double> data, vector<double> error);

  /**
   * The copy constructor.
   */
  inline Histogram(const Histogram &);

  /**
   * The destructor.
   */
  virtual ~Histogram();
  //@}

public:

  /**
   *  Operator to add a point to the histogrma with unit weight
   */
  inline void operator+=(double);

  /**
   *  Number of bins (not counting the overflow)
   */
  inline unsigned int numberOfBins();

  /**
   *  Output as a topdrawe file
   * @param out The output stream
   * @param frame output on a new graph
   * @param error output data points with error bars
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
  void topdrawOutput(ofstream & out,
		     bool frame,
		     bool error,
		     bool xlog, bool ylog,
		     string colour=string("BLACK"),
		     string title=string(),
		     string titlecase =string(),
		     string left=string(),
		     string leftcase =string(),
		     string bottom=string(),
		     string bottomcase =string());

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<Histogram> initHistogram;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Histogram & operator=(const Histogram &);

private:

  /**
   *  Number of bins
   */
  unsigned int _nbin;

  /**
   *  The contents of the bin
   */
  vector<StatisticPtr> _bincontents;

  /**
   *  The limits for the bins
   */
  vector<double> _binlimits;

  /**
   *  The data for the bins
   */
  vector<double> _data;

  /**
   *  THe error on the data
   */
  vector<double> _error;
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
  static string className() { return "Herwig++::Histogram"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the Histogram class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwKtJet.so HwAnalysis.so"; }
};

/** @endcond */

}

#include "Histogram.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Histogram.tcc"
#endif

#endif /* HERWIG_Histogram_H */

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


