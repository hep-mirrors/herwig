// -*- C++ -*-
#ifndef HERWIG_DrellYanHardGenerator_H
#define HERWIG_DrellYanHardGenerator_H
//
// This is the declaration of the DrellYanHardGenerator class.
//

#include "Herwig++/Shower/Nason/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "DrellYanHardGenerator.fh"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig++/Utilities/Histogram.h"



namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DrellYanHardGenerator class.
 *
 * @see \ref DrellYanHardGeneratorInterfaces "The interfaces"
 * defined for DrellYanHardGenerator.
 */
class DrellYanHardGenerator: public HardestEmissionGenerator {

  /**
   * Typedef for the BeamParticleData object
   */
  typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

public:

  /**
   * The default constructor.
   */
  inline DrellYanHardGenerator();

  /**
   *  Members which must be overridden in the inheriting classes
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual NasonTreePtr generateHardest(ShowerTreePtr,EvolverPtr);

  /**
   *  Member to decide if the inheriting class can handle this process
   */
  virtual bool canHandle(ShowerTreePtr);
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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

  //METHODS TULLY ADDED

 /**
   * Returns the result of the (nason matrix element) splitting function
   * for the current (yb,yj,p)
   */
  virtual double getResult();

  /**
   *  generates the hardest emission (yj,p)
   */
  void getEvent(bool & quarkfirst,unsigned int & process);

private:

 
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DrellYanHardGenerator> initDrellYanHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DrellYanHardGenerator & operator=(const DrellYanHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr _alphaS;

  /**
   *  The prefactor for the overestimate of the true distribution
   */
  double _prefactor;

  /**
   *  The power for the overestimate of the true distribution
   */
  double _power;

  /**
   *  The rapidity of the gauge boson
   */
  double _yb;

  /**
   *  The mass of the gauge boson
   */
  Energy _mass;

  /**
   *  the rapidity of the jet
   */
  double  _yj;

  /**
   *  The transverse momentum of the jet
   */
  Energy _pt;

  /**
   *  Pointers to the BeamParticleData objects
   */
  vector<tcBeamPtr> _beams;

  /**
   *  Pointers to the ParticleDataObjects for the partons
   */
  vector<tcPDPtr> _partons;

  /**
   *  Whether the quark is in the + or - z direction
   */
  bool _quarkplus;

  /**
   *  The momenta of the new particles
   */
  vector<Lorentz5Momentum> _pnew;

  HistogramPtr  _hyb;
  HistogramPtr  _hplow;
  HistogramPtr  _hphigh;
  HistogramPtr  _hyj;  
  HistogramPtr  _weighta,_weightb,_weightc;
  /**
   *  vector of points for scatter plot
   */
  vector<Energy> _ptplot;
  vector<double> _yjplot;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DrellYanHardGenerator. */
template <>
struct BaseClassTrait<Herwig::DrellYanHardGenerator,1> {
  /** Typedef of the first base class of DrellYanHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DrellYanHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DrellYanHardGenerator>
  : public ClassTraitsBase<Herwig::DrellYanHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::DrellYanHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DrellYanHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class DrellYanHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "DrellYanHardGenerator.so"; }
};

/** @endcond */

}

#include "DrellYanHardGenerator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DrellYanHardGenerator.tcc"
#endif

#endif /* HERWIG_DrellYanHardGenerator_H */
