// -*- C++ -*-
#ifndef HERWIG_GlobalParameters_H
#define HERWIG_GlobalParameters_H
//
// This is the declaration of the GlobalParameters class.

#include <ThePEG/Interface/Interfaced.h>
#include "Herwig++/Config/Herwig.h"
#include <ThePEG/CLHEPWrap/SystemOfUnits.h>
#include <ThePEG/CLHEPWrap/PhysicalConstants.h>


namespace Herwig {

using namespace ThePEG;

  /** \ingroup Utilities
   * 
   *  This class provides access to interfaced parameters which are
   *  used across the various packages (subdirectories) of Herwig++. 
   *
   *  @see Herwig
   */
class GlobalParameters: public ThePEG::Interfaced {

public:

  /**
   * Standard ctors and dtor.
   */
  inline GlobalParameters();
  inline GlobalParameters(const GlobalParameters &);
  virtual ~GlobalParameters();

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

public:

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

protected:

  /**
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();

  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<GlobalParameters> initGlobalParameters;

  /**
   * Private and non-existent assignment operator.
   */
  GlobalParameters & operator=(const GlobalParameters &);

  Energy _effectiveGluonMass;
  Energy _hadronizationScale;
  int _stringFragmentationMode;
  int _softUnderlyingEventMode;
  Energy2 _minVirtuality2;
  Length _maxDisplacement;

};

}


namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of GlobalParameters.
 */
template <>
struct BaseClassTrait<Herwig::GlobalParameters,1> {
  typedef ThePEG::Interfaced NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::GlobalParameters>: public ClassTraitsBase<Herwig::GlobalParameters> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/GlobalParameters"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "GlobalParameters.so"; }

};

}

#include "GlobalParameters.icc"

#endif /* HERWIG_GlobalParameters_H */
