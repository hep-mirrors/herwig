// -*- C++ -*-
#ifndef HERWIG_DrellYanMECorrection_H
#define HERWIG_DrellYanMECorrection_H
//
// This is the declaration of the DrellYanMECorrection class.
//

#include "MECorrectionBase.h"
#include "DrellYanMECorrection.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DrellYanMECorrection class.
 *
 * @see \ref DrellYanMECorrectionInterfaces "The interfaces"
 * defined for DrellYanMECorrection.
 */
class DrellYanMECorrection: public MECorrectionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline DrellYanMECorrection();

  /**
   * The copy constructor.
   */
  inline DrellYanMECorrection(const DrellYanMECorrection &);

  /**
   * The destructor.
   */
  virtual ~DrellYanMECorrection();
  //@}

public:

  /**
   *  Members to override those in the base class and implemented 
   *  the matrix element correction
   */
  //@{
  /**
   *  Can the matrix element correction handle a given hard process or decay
   */
  virtual bool canHandle(ShowerTreePtr);

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual void applyHardMatrixElementCorrection(ShowerTreePtr);

  /**
   * Apply the soft matrix element correction
   * @param initial The particle from the hard process which started the 
   * shower
   * @param The initial particle in the current branching
   * @param br The branching struct
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr initial,
				     ShowerParticlePtr parent,Branching br);
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

  /**
   *  Return the momenta and type of hard matrix element correction
   * @param quarks The original incoming particles.
   * @param boson The original outgoing gauge boson
   * @param iemit Whether the first (0) or second (1) particle emitted
   * the radiation
   * @param itype The type of radiated particle (0 is gluon, 1 is quark 
   *              and 2 is antiquark)
   * @param pnew The momenta of the new particles
   * @return Whether or not the matrix element correction needs to be applied
   */
  bool applyHard(ShowerParticleVector quarks, PPtr boson,unsigned int & iemit,
		 unsigned int & itype,vector<Lorentz5Momentum> & pnew);

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
  inline virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DrellYanMECorrection> initDrellYanMECorrection;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DrellYanMECorrection & operator=(const DrellYanMECorrection &);

private:

  /**
   *  Option for the generation of the phase-space
   */
  unsigned int _phasespaceopt;

  /**
   *  Relative weight for the \f$q\bar{q}\f$ and \f$qg\f$ channels
   */
  double _channelwgt;

  /**
   *  Weights for the channels as a vector
   */
  vector<double> _channelweights;
  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DrellYanMECorrection. */
template <>
struct BaseClassTrait<Herwig::DrellYanMECorrection,1> {
  /** Typedef of the first base class of DrellYanMECorrection. */
  typedef Herwig::MECorrectionBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DrellYanMECorrection class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DrellYanMECorrection>
  : public ClassTraitsBase<Herwig::DrellYanMECorrection> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::DrellYanMECorrection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DrellYanMECorrection is implemented. It may also include several, space-separated,
   * libraries if the class DrellYanMECorrection depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNewShower.so"; }
};

/** @endcond */

}

#include "DrellYanMECorrection.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DrellYanMECorrection.tcc"
#endif

#endif /* HERWIG_DrellYanMECorrection_H */
