// -*- C++ -*-
#ifndef HERWIG_TopDalitzAnalysis_H
#define HERWIG_TopDalitzAnalysis_H
//
// This is the declaration of the TopDalitzAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "TopDalitzAnalysis.fh"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "Herwig++/Interfaces/KtJetInterface.h"
#include "../Utilities/Histogram.h"
#include "KtJet/KtJet.h"
#include "KtJet/KtLorentzVector.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the TopDalitzAnalysis class.
 *
 * @see \ref TopDalitzAnalysisInterfaces "The interfaces"
 * defined for TopDalitzAnalysis.
 */
class TopDalitzAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline TopDalitzAnalysis();

  /**
   * The copy constructor.
   */
  inline TopDalitzAnalysis(const TopDalitzAnalysis &);
  //@}

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);

  /**
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   */
  virtual void analyze(tPPtr particle);

  /**
   *  Identifies which step(2) final state particles originate
   *  from the top/antitop, which originate from the b/bbar...
   */
  tPVector particleID(PPtr,tPVector);

  /**
   *  Function to cluster each decay to 2 jets and calculate the 
   *  corresponding point on the Dalitz plot.
   */
  void     dalitz(tPVector);

  /**
   *  Function to cluster to 3 jets and calculate delta(R).
   */
  void     threeJetAnalysis(Energy2,tPVector,tPVector);

  /**
   *  Histogram for delta R.
   */
  Histogram _deltaR;

  /**
   *  Histogram for log(y3).
   */
  Histogram _logy3;

  /**
   *  Histogram for the b-quark energy spectrum.
   */
  Histogram _xb_bquark;

  /**
   *  Histogram for the b-quark energy spectrum around the Sudakov peak.
   */
  Histogram _xb_bquark_peak;

  /**
   *  Histogram for the B-hadron energy spectrum.
   */
  Histogram _xB_Bhad;
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
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<TopDalitzAnalysis> initTopDalitzAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TopDalitzAnalysis & operator=(const TopDalitzAnalysis &);

private:

  /**
   *  Output stream
   */
  ofstream _output[6];

  /**
   *  Number of outputs
   */
  unsigned int _nout;

  /**
   *  The interface between Herwig++ and KtJet
   */
  Herwig::KtJetInterface _kint;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TopDalitzAnalysis. */
template <>
struct BaseClassTrait<Herwig::TopDalitzAnalysis,1> {
  /** Typedef of the first base class of TopDalitzAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TopDalitzAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TopDalitzAnalysis>
  : public ClassTraitsBase<Herwig::TopDalitzAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TopDalitzAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * TopDalitzAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class TopDalitzAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwKtJet.so HwLEPJetAnalysis.so"; }
};

/** @endcond */

}

#include "TopDalitzAnalysis.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TopDalitzAnalysis.tcc"
#endif

#endif /* HERWIG_TopDalitzAnalysis_H */
