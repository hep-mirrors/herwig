// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#ifndef Analysis2_GnuplotOutput_H
#define Analysis2_GnuplotOutput_H
//
// This is the declaration of the GnuplotOutput class.
//

#include "Analysis2/Histogram2/Histogram2Output.h"
#include "GnuplotOutput.fh"

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 * 
 * Gnuplot output of Histogram2
 *
 * @see \ref GnuplotOutputInterfaces "The interfaces"
 * defined for GnuplotOutput.
 */
class GnuplotOutput: public Histogram2Output {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline GnuplotOutput();

  /**
   * The copy constructor.
   */
  inline GnuplotOutput(const GnuplotOutput&);

  /**
   * The destructor.
   */
  virtual ~GnuplotOutput();
  //@}

public:

  /**
   * Output the given histogram.
   * The default just dumps to an ascii file.
   */
  virtual void put (Histogram2Ptr, const Histogram2Options&, const string&);

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

public:

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

private:

  /**
   * The range of subprocess multiplicities
   * to plot, separated by a blank. If empty,
   * no subprocess multiplicity channels are
   * plotted.
   */
  string _subproMult;

  /**
   * Count the histograms already put.
   */
  unsigned int _allHistos;

  /**
   * Wether or not to plot the ratio
   * (if available)
   */
  bool _ratio;

  /**
   * Wether or not to plot the chi2
   * (if available)
   */
  bool _chi2;

  /**
   * The gnuplot program
   */
  string _gnuplot;

  /**
   * The TeX file containing the plots
   */
  ofstream _TeXfile;

  /**
   * The Makefile
   */
  ofstream _makefile;

  /**
   * Wether or not initialilzation has been done.
   */
  bool _doneInitialize;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<GnuplotOutput> initGnuplotOutput;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GnuplotOutput & operator=(const GnuplotOutput &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GnuplotOutput. */
template <>
struct BaseClassTrait<Analysis2::GnuplotOutput,1> {
  /** Typedef of the first base class of GnuplotOutput. */
  typedef Analysis2::Histogram2Output NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GnuplotOutput class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::GnuplotOutput>
  : public ClassTraitsBase<Analysis2::GnuplotOutput> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::GnuplotOutput"; }
  /**
   * The name of a file containing the dynamic library where the class
   * GnuplotOutput is implemented. It may also include several, space-separated,
   * libraries if the class GnuplotOutput depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "Analysis2.so"; }
};

/** @endcond */

}

#include "GnuplotOutput.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GnuplotOutput.tcc"
#endif

#endif /* Analysis2_GnuplotOutput_H */
