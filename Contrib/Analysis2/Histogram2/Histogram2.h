// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#ifndef Analysis2_Histogram2_H
#define Analysis2_Histogram2_H
//
// This is the declaration of the Histogram2 class.
//

#include <cassert>

#include "ThePEG/Interface/Interfaced.h"
#include "Histogram2.fh"

#include "ThePEG/Utilities/StringUtils.h"

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 * A channel in a histogram.
 *
 * @author Simon Plaetzer
 */
class HistogramChannel {

public:

  /**@name Constructors */
  //@{

  /**
   * Default constructor
   */
  inline HistogramChannel ();

  /**
   * Construct giving the number of bins.
   */
  inline explicit HistogramChannel (unsigned int bins, bool counting = true);

  /**
   * Construct giving contents of bins.
   */
  inline explicit HistogramChannel (const vector<pair<double,double> >&,
				    double underflow = 0., double overflow = 0.);
  //@}

public:

  /**@name Methods to book events */
  //@{
  /**
   * Book the weight into the bin with
   * the given index.
   */
  inline void book (unsigned int bin, double weight = 1.);

  /**
   * Book the weight into the underflow
   */
  inline void bookUnderflow (double weight = 1.);

  /**
   * Book the weight into the overflow
   */
  inline void bookOverflow (double weight = 1.);

  /**
   * Signal a numerical indefinite event
   */
  inline void nanEvent ();
  //@}

public:

  /**@name Channel algebra. It takes care of properly computing uncertainties. */
  //@{
  /**
   * Add another channel.
   */
  HistogramChannel& operator += (const HistogramChannel&);

  /**
   * Subtract another channel.
   */
  HistogramChannel& operator -= (const HistogramChannel&);

  /**
   * Multiply by another channel.
   */
  HistogramChannel& operator *= (const HistogramChannel&);

  /**
   * Multiply the channel by a factor.
   */
  HistogramChannel& operator *= (double);

  /**
   * Add a constant to each bin.
   */
  HistogramChannel& operator += (double);

  /**
   * Multiply the channel by a factor
   * including uncertainties.
   */
  HistogramChannel& operator *= (pair<double,double>);

  /**
   * Divide two channels.
   */
  HistogramChannel& operator /= (const HistogramChannel&);

  /**
   * Divide the channel by a factor.
   */
  inline HistogramChannel& operator /= (double);

  /**
   * Divide the channel by a factor
   * including uncertainties.
   */
  HistogramChannel& operator /= (pair<double,double>);
  //@}

public:

  /**@name Access channel statistics. */
  //@{
  /**
   * Return true, if this is a counting channel.
   */
  inline bool isCountingChannel () const;

  /**
   * Return the bins.
   */
  inline vector<pair<double,double> > bins () const;

  /**
   * Return the bin content for the
   * given index. The first entry is
   * the sum of weights, the second
   * the sum of squared weights.
   */
  inline pair<double,double> bin (unsigned int) const;

  /**
   * Return the bin entries.
   */
  inline vector<unsigned long> binEntries () const;

  /**
   * Return the bin entries for the
   * given index.
   */
  inline unsigned long binEntries (unsigned int) const;

  /**
   * Explicitly set the content of the
   * bin with given index.
   */
  inline void bin (unsigned int, pair<double,double>, unsigned long en = 0);

  /**
   * Return under- and overflow
   */
  inline pair<double,double> outOfRange () const;

  /**
   * Return the visible entries
   */
  inline unsigned long visible () const;

  /**
   * Return the total entries
   */
  inline unsigned long total () const;

  /**
   * Return the histogram of numerical
   * indefinite weights
   */
  inline vector<unsigned long> nanWeights () const;

  /**
   * Return the total number of events with
   * indefinite weights.
   */
  unsigned long nanWeightEvents () const;

  /**
   * Return the number of numerical indefinite events
   */
  inline unsigned long nanEvents () const;
  //@}

public:

  /**@name Normalization and statistics */
  //@{

  /**
   * Finish this channel. Provided for future
   * use.
   */
  inline void finish ();

  /**
   * Return the variance of the bin
   * content for the given bin.
   */
  inline double binVariance (unsigned int) const;

  /**
   * Return the mean of weights
   * in the given bin.
   */
  inline double weightMean (unsigned int) const;

  /**
   * Return the variance of weights
   * in the given bin.
   */
  inline double weightVariance (unsigned int) const;

  /**
   * This channel is a differential distribution:
   * divide each bin by its width.
   */
  void differential (const vector<pair<double,double> >&);

  /**
   * Return the summed visible bin content
   */
  pair<double,double> binSum () const;

  /**
   * Return the average visible bin content
   */
  inline pair<double,double> binAverage () const;

  /**
   * Return the visible integral of the channel
   * using the given binning.
   */
  pair<double,double> integrate (const vector<pair<double,double> >&) const;

  /**
   * Return the weighted average visible bin content
   */
  pair<double,double> average (const vector<pair<double,double> >&) const;

  /**
   * Return a channel containing this/channel - 1
   */
  HistogramChannel delta (const HistogramChannel& channel) const;

  /**
   * Return a channel containing the \f$\chi^2\f$ taking
   * this as hypothesis and the given channel as data.
   * Rescale the data error to data*minfrac, if
   * data error/data < minfrac.
   */
  HistogramChannel chi2 (const HistogramChannel& channel, double minfrac = .0) const;

  /**
   * Generate a channel containing
   * the profile histogram obtained from this
   * channel, i.e. a channel containing
   * the mean of weights and its variance
   * as content/content square values
   */
  HistogramChannel profile () const;

  //@}

public:

  /**@name Utilities for store/load */
  //@{

  /**
   * Write out the channel
   */
  void write (ostream&, const string&);

  /**
   * Read in the channel
   */
  string read (istream&);

  //@}

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
   */
  void persistentInput(PersistentIStream & is);
  //@}

private:

  /**
   * Wether or not the number of entries
   * is meaningfull.
   */
  bool _isCountingChannel;

  /**
   * The contents and uncertainties squared,
   * to be accessed by the bin id.
   */
  vector<pair<double,double> > _bins;

  /**
   * The number of entries in each bin,
   * to be accessed by the bin id.
   */
  vector<unsigned long> _binEntries;

  /**
   * The contents of under- and overflow
   */
  pair<double,double> _outOfRange;

  /**
   * The number of visible entries
   */
  unsigned long _visible;

  /**
   * The number of total entries
   */
  unsigned long _total;

  /**
   * The number of numerical indefinite
   * events.
   */
  unsigned long _nanEvents;

  /**
   * A histogram of the number of
   * indefinite weights.
   */
  vector<unsigned long> _nanWeights;

  /**
   * Has finish() already been called?
   */
  bool _finished;

};

/**\ingroup Analysis2
 * Persistent output of histogram channel
 */
inline PersistentOStream& operator << (PersistentOStream& os, const HistogramChannel&);

/**\ingroup Analysis2
 * Persistent input of histogram channel
 */
inline PersistentIStream& operator >> (PersistentIStream& is, HistogramChannel&);

/**
 * Add two histogram channels
 */
inline HistogramChannel operator + (const HistogramChannel& a, const HistogramChannel& b);

/**
 * Subtract two histogram channels
 */
inline HistogramChannel operator - (const HistogramChannel& a, const HistogramChannel& b);

/**
 * Multiply two histogram channels
 */
inline HistogramChannel operator * (const HistogramChannel& a, const HistogramChannel& b);

/**
 * Divide two histogram channels
 */
inline HistogramChannel operator / (const HistogramChannel& a, const HistogramChannel& b);

/**
 * Get the next xml tag from a stream,
 * i.e. everything not being a blank line and starting
 * (module white-spaces) with '<' followed by
 * an alphabetical character or /
 */
inline string getNextTag (istream&);

/**
 * Convert a string to the template
 * type by reading is in from an istringstream.
 */
template<class T>
inline void fromString (const string& str, T& result) {
  istringstream buffer (str);
  buffer >> result;
}

/**\ingroup Analysis2
 * Options for channel output
 */
namespace ChannelOutput {

  /**
   * Default output: bin_low bin_high content error
   */
  const unsigned int Default     = 0;

  /**
   * Instead of bin_low _bin_high put a single
   * column (bin_low+bin_high)/2
   */
  const unsigned int Bincenters  = 1;

  /**
   * Do not output uncertainties
   */
  const unsigned int NoErrorbars = 1 << 1;

  /**
   * Output the histogram of events
   * with numerical indefinite weights
   */
  const unsigned int NanEvents   = 1 << 2;

  /**
   * Output global statistics in a comment
   */
  const unsigned int Statistics  = 1 << 3;

}

/**\ingroup Analysis2
 * A more sophisticated histogram. It is almost completely
 * compatible to Histogram, but will use a different
 * way to actually dump the histogram for use with
 * the favorite plotting program.
 *
 * @see \ref Histogram2Interfaces "The interfaces"
 * defined for Histogram2.
 */
class Histogram2: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{

  /**
   * The default constructor.
   */
  inline Histogram2 ();

  /**
   * The destructor.
   */
  virtual ~Histogram2();
  //@}

public:

  /**@name Constructors */
  //@{
  /**
   * Construct giving the range and a number of
   * bins and insert a counting channel with name.
   * If name is empty, no channel is inserted.
   * Default name is given for backward compatibility.
   */
  explicit Histogram2 (double low, double high,
		       unsigned int bins,
		       const string& name = "mc");

  /**
   * Construct giving a binning and insert a counting
   * channel with name.
   * If name is empty, no channel is inserted.
   */
  explicit Histogram2 (const vector<pair<double,double> >&, const string& name);

  /**
   * Construct giving a name of a data file
   * and a name for the resulting channel.
   * The data file is assumed to be in the format
   * bin_low bin_high data error
   */
  explicit Histogram2 (const string& dataFile, const string& dataName);

  /**
   * Set a name for the histogram
   */
  inline void setName (const string&);

  //@}

public:

  /**
   * Book an event into the given channel
   */
  void book (const string& name, double event, double weight = 1.);

  /**@name Channel administration */
  //@{
  /**
   * Return true, if a channel by that name
   * exists.
   */
  inline bool haveChannel (const string&) const;

  /**
   * Access a channel
   */
  inline HistogramChannel& channel (const string&);

  /**
   * Access a channel
   */
  inline HistogramChannel channel (const string&) const;

  /**
   * Insert a new channel
   */
  inline void insertChannel (const string&, const HistogramChannel&);

  /**
   * Insert a new counting channel
   */
  inline void insertChannel (const string&);

  /**
   * Remove a channel
   */
  inline HistogramChannel removeChannel (const string&);
  //@}

public:

  /**@name Access the binning */
  //@{
  /**
   * Return the binning
   */
  inline const vector<pair<double,double> >& binning () const;

  /**
   * Return the bin hash map
   */
  inline const map<double,unsigned int>& binhash () const;

  /**
   * Return the range
   */
  inline pair<double,double> range () const;
  //@}

public:

  /**@name Normalisation and statistical tests */
  //@{

  /**
   * Finish the given channel
   */
  inline void finish (const string&);

  /**
   * The given channel is a differential distribution,
   * dvidie each bin by its bin width.
   */
  inline void differential (const string&);

  /**
   * Return the integral of the given channel
   */
  inline pair<double,double> integrate (const string&) const;

  /**
   * Normalise the given channel to unity
   */
  inline void normalise (const string&);

  /**
   * Normalise the given channel to the
   * integral of another one.
   */
  inline void normalise (const string&, const string&);

  /**
   * Rescale the given channel
   */
  inline void rescale (const string&, double);

  /**
   * Rescale the all channels
   */
  inline void rescale (double);

  /**
   * Normalise a channel to the total cross section.
   * Default argument for backwards compatibility.
   */
  inline void normaliseToCrossSection(const string& name = "mc");

  /**
   * Return the \f$\chi^2\f$ per DOF taking the first
   * channel as hypotheses and the second one as data/
   */
  inline pair<double,double> chi2perDOF (const string&, const string&, double minfrac = .0) const;
  //@}

public:

  /**@name Support for histogram output */
  //@{
  /**
   * Return the list of all channels
   * in this histogram
   */
  vector<string> channels () const;

  /**
   * Output a channel to the given
   * ostream as a ascii table.
   */
  void output (ostream&, const string& name,
	       unsigned int flags = 0,
	       char comment = '#') const;
  //@}

public:

  /*@name Support for parallel runs */
  //@{
  /**
   * Write out the histogram to a file
   * name.h2
   */
  void store (const string& name);

  /**
   * Read in the histogram from a file
   * name.h2 . Return false, if something
   * did go wrong.
   */
  bool load (const string& name);

  /**
   * Read in the histogram from a file
   * name.h2 into a new histogram. It returns
   * null, if reading failed.
   */
  Histogram2Ptr loadToHistogram (const string& name) const;

  /**
   * Combine histograms. Each run is assumed
   * to be in a subdirectory prefix.runId,
   * having written the histogram to a file
   * prefix.runId/name.h2
   *
   * runId is assumed to range from 0 to numRuns-1
   *
   * dataChannel indicates a channel which is
   * to be included only once.
   *
   * All other channels are added into this
   * histogram. If a mcChannel is given,
   * this is used to compute the total cross section
   * for the combined histogram.
   */
  void combine (const string& prefix, const string& name,
		unsigned int numRuns, const string& dataChannel = "",
		const string& mcChannel = "MC");

  /**
   * Return the cross section accumulated in this histogram.
   */
  inline CrossSection xSec () const;

  /**
   * Set the cross section accumulated in this histogram.
   */
  inline void xSec (CrossSection);

  //@}

public:

  /**@name Backward compatibility to Histogram. Note that prefactor is
   * is not supported (instead use multiplication operator and/or
   * the normalise methods for the channel(s) to be considered).
   */
  //@{

  /**
   * Construct from vector of lower limits and
   * upper limit of last bin and insert an initial counting
   * channel mc.
   */
  explicit Histogram2 (vector<double> limits);

  /**
   * Construct giving the lower limits, followed by the
   * upper limit of the last bin. Inserts channels data
   * and mc
   */
  explicit Histogram2 (vector<double> limits,
		       vector<double> data,
		       vector<double> dataerror);

  /**
   * Operator to add a point to the histogrma with unit weight.
   */
  inline void operator+=(double);

  /**
   * Function to add a weighted point to the histogram.
   */
  inline void addWeighted(double data, double weight);

  /**
   *  Number of bins (not counting the overflow)
   */
  inline unsigned int numberOfBins() const;

  /**
   * Normalise the distributions to the data.
   */
  inline void normaliseToData();

  /**
   * Return the chi squared
   * It assumes two channels to be present:
   * mc and data. This is not checked!
   */
  inline void chiSquared(double & chisq, 
			 unsigned int & ndegrees,
			 double minfrac=0.) const;

  /**
   * Returns a new histogram containing bin-by-bin
   * ratios of two histograms
   */
  Histogram2 ratioWith(const Histogram2& h2) const;

  //@}

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

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The binning. This identifies bin ids as the
   * vector indices.
   */
  vector<pair<double,double> > _binning;

  /**
   * Seperatly keep the lowest bins lower and the highest
   * bins uper bound to identify events ind under-/overflow.
   */
  pair<double,double> _range;

  /**
   * Bin hashing: Map each bins upper bound
   * to the bin id.
   */
  map<double, unsigned int> _binhash;

  /**
   * Map channels by their name
   */
  map<string,HistogramChannel> _channels;

  /**
   * The cross section this histogram accumulated
   * according to the event generator.
   */
  CrossSection _xSec;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static ClassDescription<Histogram2> initHistogram2;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Histogram2 & operator=(const Histogram2 &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Histogram2. */
template <>
struct BaseClassTrait<Analysis2::Histogram2,1> {
  /** Typedef of the first base class of Histogram2. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Histogram2 class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::Histogram2>
  : public ClassTraitsBase<Analysis2::Histogram2> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::Histogram2"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Analysis2Base is implemented. It may also include several, space-separated,
   * libraries if the class Analysis2Base depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "Analysis2.so"; }
};

/** @endcond */

}

#include "Histogram2.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Histogram2.tcc"
#endif

#endif /* Analysis2_Histogram2_H */
