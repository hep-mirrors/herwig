// -*- C++ -*-
#ifndef HERWIG_GlobalParameters_H
#define HERWIG_GlobalParameters_H
//
// This is the declaration of the GlobalParameters class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Config/Herwig.h"
#include <ThePEG/CLHEPWrap/SystemOfUnits.h>
#include <ThePEG/CLHEPWrap/PhysicalConstants.h>
#include "GlobalParameters.fh"

namespace Herwig {

using namespace ThePEG;

  /** \ingroup Utilities
   * 
   *  This class provides access to interfaced parameters which are
   *  used across the various packages of Herwig++. 
   *
   * @see \ref GlobalParametersInterfaces "The interfaces"
   * defined for GlobalParameters.
   */
class GlobalParameters: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline GlobalParameters();

  /**
   * The copy constructor.
   */
  inline GlobalParameters(const GlobalParameters &);

  /**
   * The destructor.
   */
  virtual ~GlobalParameters();
  //@}

public:

  /**
   *  Access to the various parameters
   */
  //@{
  /**
   * It returns the effective gluon mass, which is necessary for
   * the cluster hadronization handler; however, it is used also
   * at the end of the showering, because setting the physical
   * massless shell gluons on this effective mass shell requires
   * kinematical reshuffling in order to conserve energy-momentum,
   * and this procedure is the same as done in the showering to
   * compensate the recoil of the emission.
   */
  inline Energy effectiveGluonMass() const;

  /** 
   * It returns roughly the hadronization scale, that is the energy
   * scale such that if a particle has a width above this value,
   * then the particle decays before hadronizing (if it is coloured)
   * and in any case (even if it is not coloured, in order to proper
   * handle non-QCD radiation) its decay must be treated inside the 
   * Shower. Notice that this parameter is not used anywhere in the 
   * Hadronization, but only by the Shower when the multi-scale option
   * is on (but we prefer anywhere to put it here, rather than 
   * in one of Shower classes, because it could be useful somewhere 
   * else in the future, like in the multiparton or soft model).
   * Notice that if you want multi-scale showering but without
   * considering the widths of W, Z, top as a scale in the showering,
   * it is enough to set this parameter above such widths.        
   */
  inline Energy hadronizationScale() const;

  /**
   * It returns true/false according if the ThePEG string fragmentation
   * model is switched on/off. In the case is off, the usual
   * Herwig++ cluster hadronization model will be used for the
   * hadronization.
   */   
  inline bool isThePEGStringFragmentationON() const;

  /**
   * It returns true/false according if the soft underlying model
   * is switched on/off. 
   */
  inline bool isSoftUnderlyingEventON() const;

  /**
   * It returns minimum virtuality^2 of partons to use in calculating 
   * distances. It is used both in the Showering and Hadronization.
   */
  inline Energy2 minVirtuality2() const;

  /**
   * It returns the maximum displacement that is allowed for a particle
   * (used to determine the position of a cluster with two components).
   */
  inline Length maxDisplacement() const;

  /**
   * It returns the conversion factor from GeV to Millimeter.
   */
  inline double conversionFactorGeVtoMillimeter() const;
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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<GlobalParameters> initGlobalParameters;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GlobalParameters & operator=(const GlobalParameters &);

private:

  /**
   *  The effective gloun mass for the cluster model.
   */
  Energy _effectiveGluonMass;

  /**
   *  Hadronization scale
   */
  Energy _hadronizationScale;

  /**
   *  Are we using string fragmentation
   */
  int _stringFragmentationMode;

  /**
   *  Is the soft underlying event on/off
   */
  int _softUnderlyingEventMode;

  /**
   * The minimum virtuality^2 of partons to use in calculating 
   * distances.
   */
  Energy2 _minVirtuality2;

  /**
   * The maximum displacement that is allowed for a particle
   * (used to determine the position of a cluster with two components).
   */
  Length _maxDisplacement;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GlobalParameters. */
template <>
struct BaseClassTrait<Herwig::GlobalParameters,1> {
  /** Typedef of the first base class of GlobalParameters. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GlobalParameters class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GlobalParameters>
  : public ClassTraitsBase<Herwig::GlobalParameters> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::GlobalParameters"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the GlobalParameters class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "libHwUtils.so"; }
};

/** @endcond */

}

#include "GlobalParameters.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GlobalParameters.tcc"
#endif

#endif /* HERWIG_GlobalParameters_H */
