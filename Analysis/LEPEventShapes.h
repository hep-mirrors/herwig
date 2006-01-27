// -*- C++ -*-
#ifndef HERWIG_LEPEventShapes_H
#define HERWIG_LEPEventShapes_H
//
// This is the declaration of the LEPEventShapes class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "EventShapes.h"
#include "LEPEventShapes.fh"
#include "Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The LEPEventShapes class performs the analysis of global event shapes and
 * compares with LEP data. This handler is solely intended as a slave
 * handler for the EventShapesMasterAnalysis class.
 *
 * @see \ref LEPEventShapesInterfaces "The interfaces"
 * defined for LEPEventShapes
 */
class LEPEventShapes: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline LEPEventShapes();

  /**
   * The copy constructor.
   */
  inline LEPEventShapes(const LEPEventShapes &);

  /**
   * The destructor.
   */
  virtual ~LEPEventShapes();
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

protected:

  /**
   *  Methods to compute the four jet angles, assumes the jets are energy ordered
   */
  //@{
  /**
   *  Compute \f$\cos\chi_{BZ}\f$
   */
  inline double cosChiBZ(vector<Lorentz5Momentum>);

  /**
   *  Compute \f$\cos\Phi_{KSW}\f$.
   */ 
  inline double cosPhiKSW(vector<Lorentz5Momentum>);
  
  /**
   *  Compute \f$\cos\Theta_{NR}\f$
   */
  inline double cosThetaNR(vector<Lorentz5Momentum>); 

  /**
   *  Compute \f$\cos\alpha_{34}\f$
   */
  inline double cosAlpha34(vector<Lorentz5Momentum>); 
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<LEPEventShapes> initLEPEventShapes;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LEPEventShapes & operator=(const LEPEventShapes &);

private:

  /**
   *  Histogram for \f$1-T\f$ distribution.
   */
  HistogramPtr _omthr;

  /**
   * Histogram for the major distribution
   */
  HistogramPtr _maj;

  /**
   * Histogram for the minor distribution
   */
  HistogramPtr _min;

  /**
   * Histogram for the oblateness distribution
   */
  HistogramPtr _obl;

  /**
   * Histogram for the sphericity distribution
   */
  HistogramPtr _sph;

  /**
   * Histogram for the aplanarity distribution
   */
  HistogramPtr _apl;

  /**
   * Histogram for the planarity distribution
   */
  HistogramPtr _pla;

  /**
   * Histogram for the C distribution
   */
  HistogramPtr _c;

  /**
   * Histogram for the D distribution
   */
  HistogramPtr _d;

  /**
   * Histogram for the \f$M_{\rm high}\f$ distribution
   */
  HistogramPtr _mhi;

  /**
   * Histogram for the \f$M_{\rm low}\f$ distribution
   */
  HistogramPtr _mlo;

  /**
   * Histogram for the \f$M_{\rm high}-M_{\rm low}\f$ distribution
   */
  HistogramPtr _mdiff;

  /**
   * Histogram for the \f$B_{\rm max}\f$ distribution
   */
  HistogramPtr _bmax;

  /**
   * Histogram for the \f$B_{\rm min}\f$ distribution
   */
  HistogramPtr _bmin;

  /**
   * Histogram for the \f$B_{\rm max}+B_{\rm min}\f$ distribution
   */
  HistogramPtr _bsum;

  /**
   * Histogram for the \f$B_{\rm max}-B_{\rm min}\f$ distribution
   */
  HistogramPtr _bdiff;

  /**
   * Pointer to the object which calculates the event shapes
   */
  EventShapesPtr _shapes;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LEPEventShapes. */
template <>
struct BaseClassTrait<Herwig::LEPEventShapes,1> {
  /** Typedef of the first base class of LEPEventShapes. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LEPEventShapes class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LEPEventShapes>
  : public ClassTraitsBase<Herwig::LEPEventShapes> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::LEPEventShapes"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the LEPEventShapes class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwKtJet.so HwAnalysis.so"; }
};

/** @endcond */

}

#include "LEPEventShapes.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LEPEventShapes.tcc"
#endif

#endif /* HERWIG_LEPEventShapes_H */
