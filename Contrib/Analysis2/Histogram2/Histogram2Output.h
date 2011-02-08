// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#ifndef Analysis2_Histogram2Output_H
#define Analysis2_Histogram2Output_H
//
// This is the declaration of the Histogram2Output class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Histogram2Output.fh"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Histogram2.h"

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 *
 * Options for histogram output. 
 * They can be combined using the '|' operator, e.g. 'Frame | Ylog'
 */
namespace HistogramOutput {
  /** Default behaviour */
  const unsigned int Default   = 0;
  /** Frame */
  const unsigned int Frame     = 1;
  /** Errorbars */
  const unsigned int Errorbars = 1 << 1;
  /** Log x axis */
  const unsigned int Xlog      = 1 << 2;
  /** Log y axis */
  const unsigned int Ylog      = 1 << 3;
  /** Smooth */
  const unsigned int Smooth    = 1 << 4;
  /** Rawcount */
  const unsigned int Rawcount  = 1 << 5;
}

/**\ingroup Analysis2
 * Simple struct to contain options for plotting
 * output of an observable.
 */
struct Histogram2Options {

  /**
   * Default constructor
   */  
  inline Histogram2Options ()
    : plotFlags(0), channelFlags(0), title(""), datatitle(""),
      xlabel(""), ylabel (""), differential(true) {}
  
  /**
   * Constructor giving initial values
   */
  inline Histogram2Options (int pFlags,
			    int cFlags = 0,
			    string t = "",
			    string dt = "",
			    string x = "",
			    string y = "",
			    bool diff = true)
    : plotFlags(pFlags), channelFlags(cFlags), title(t),
      datatitle(dt), xlabel(x), ylabel (y), differential(diff) {}
  
  /** Flags for plotting */
  int plotFlags;

  /** Flags for channel output */
  int channelFlags;

  /** Title of the histogram */  
  string title;

  /** Title of data, if present */
  string datatitle;
  
  /** X label */
  string xlabel;

  /** Y label */
  string ylabel;

  /** wether or not this is a differential observable */
  bool differential;
  
};

/**\ingroup Analysis2
 * Base class for Histogram2 output.
 *
 * @author Simon Plaetzer
 *
 * @see \ref Histogram2OutputInterfaces "The interfaces"
 * defined for Histogram2Output.
 */
class Histogram2Output: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Histogram2Output();

  /**
   * The copy constructor.
   */
  inline Histogram2Output(const Histogram2Output &);

  /**
   * The destructor.
   */
  virtual ~Histogram2Output();
  //@}

public:

  /**
   * Initialize the curret ostream to point
   * to name. This is a helper only, it has
   * to be called by either the implemented AnalysisHandler
   * or the out method.
   */
  void initialize (const string& name);

  /**
   * Output the given histogram.
   * The default just dumps to an ascii file.
   */
  virtual void put (Histogram2Ptr, const Histogram2Options&, const string& dataChannel);

  /**
   * Initialize. Similar to doinitrun,
   * in fact called by doinitrun, but
   * should _not_ call generator(), as it
   * may be used 'offline' to combine
   * parallel runs.
   */
  virtual inline void initialize ();

  /**
   * Finalize. Similar to dofinish,
   * in fact called by dofinish, but
   * should _not_ call generator(), as it
   * may be used 'offline' to combine
   * parallel runs.
   */
  virtual inline void finalize ();

protected:

  /**
   * Return the currently used ostream
   */
  virtual ostream& currentOStream ();

  /**
   * Get the output prefix or filename
   */
  inline string prefix () const;

  /**
   * Get the model name
   */
  inline string mctitle () const;

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


protected:

  /** @name Standard Interfaced functions. */
  //@{
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
  //@}


private:

  /**
   * The prefix
   */
  string _prefix;

  /**
   * The MC model used
   */
  string _mctitle;

  /**
   * Pointer to currently used ostream
   */
  ofstream _out;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<Histogram2Output> initHistogram2Output;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Histogram2Output & operator=(const Histogram2Output &);

};

/**\ingroup Analysis2
 * Persistent output of histogram options
 */
inline PersistentOStream& operator << (PersistentOStream& os, const Histogram2Options& options);

/**\ingroup Analysis2
 * Persistent input of histogram options
 */
inline PersistentIStream& operator >> (PersistentIStream& is, Histogram2Options& options);

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Histogram2Output. */
template <>
struct BaseClassTrait<Analysis2::Histogram2Output,1> {
  /** Typedef of the first base class of Histogram2Output. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Histogram2Output class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::Histogram2Output>
  : public ClassTraitsBase<Analysis2::Histogram2Output> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::Histogram2Output"; }
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

#include "Histogram2Output.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Histogram2Output.tcc"
#endif

#endif /* Analysis2_Histogram2Output_H */
