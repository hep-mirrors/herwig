// -*- C++ -*-
#ifndef HERWIG_SMWFermionsPOWHEGDecayer_H
#define HERWIG_SMWFermionsPOWHEGDecayer_H
//
// This is the declaration of the SMWFermionsPOWHEGDecayer class.
//

#include "SMWDecayer.h"
#include "Herwig/Utilities/Maths.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SMWFermionsPOWHEGDecayer class.
 *
 * @see \ref SMWFermionsPOWHEGDecayerInterfaces "The interfaces"
 * defined for SMWFermionsPOWHEGDecayer.
 */
class SMWFermionsPOWHEGDecayer: public SMWDecayer {

public:

  /**
   * The default constructor.
   */
  SMWFermionsPOWHEGDecayer();

  /**
   *  Virtual members to be overridden by inheriting classes
   *  which implement hard corrections 
   */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return FSR;}


  /**
   *  Apply the POWHEG style correction
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr);
  //@}

  virtual double me2(const int ichan, const Particle & part,
		     const ParticleVector & decay, MEOption meopt) const;

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMWFermionsPOWHEGDecayer> initSMWFermionsPOWHEGDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMWFermionsPOWHEGDecayer & operator=(const SMWFermionsPOWHEGDecayer &);

  /**
   *  Pointer to the fermion-antifermion W vertex
   */
  AbstractFFVVertexPtr FFWVertex() const {return FFWVertex_;}

  /**
   *  Pointer to the fermion-antifermion G vertex
   */
  AbstractFFVVertexPtr FFGVertex() const {return FFGVertex_;}

  /**
   *  Real emission term, for use in generating the hardest emission
   */
  double calculateRealEmission(double x1, double x2, 
			       vector<PPtr> hardProcess,
			       double phi, double muj, double muk,
			       int iemit, bool subtract) const;

  /**
   *  Check the sign of the momentum in the \f$z\f$-direction is correct.
   */
  bool checkZMomenta(double x1, double x2, double x3, double y, Energy pT,
		     double muj, double muk) const;

  /**
   *  Calculate the Jacobian
   */
  InvEnergy calculateJacobian(double x1, double x2, Energy pT,
			      double muj, double muk) const;

  /**
   *  Calculate the ratio between NLO & LO ME
   */
  double meRatio(vector<cPDPtr> partons, 
		 vector<Lorentz5Momentum> momenta,
		 unsigned int iemitter,bool subtract) const;
  /**
   *  Calculate the LO ME
   */
  double loME(const vector<cPDPtr> & partons, 
	      const vector<Lorentz5Momentum> & momenta) const;

  /**
   *  Calculate the NLO real emission piece of ME
   */
  InvEnergy2 realME(const vector<cPDPtr> & partons, 
		  const vector<Lorentz5Momentum> & momenta) const;

  /**
   *  Generate a real emission event
   */
  bool getEvent(vector<PPtr> hardProcess);

private:

  /**
   *  Pointer to the fermion-antifermion W vertex
   */
  AbstractFFVVertexPtr FFWVertex_;

  /**
   *  Pointer to the fermion-antifermion G vertex
   */
  AbstractFFVVertexPtr FFGVertex_;

  /**
   *  The colour factor 
   */
  double CF_;

  /**
   *  The W mass
   */
  mutable Energy mW_;


  // TODO: delete this
  mutable double mu_;

  /**
   *  The reduced mass of particle 1
   */
  mutable double mu1_;
  /**
   *  The reduced mass of particle 1 squared
   */
  mutable double mu12_;

  /**
   *  The reduceed mass of particle 2
   */
  mutable double mu2_;

  /**
   *  The reduceed mass of particle 2 squared
   */
  mutable double mu22_;


  /**
   *  The strong coupling
   */
  mutable double aS_;

  /**
   * The scale
   */
  mutable Energy2 scale_;

  /**
   *  Stuff for the POWHEG correction
   */
  //@{
  /**
   *  ParticleData object for the gluon
   */
  tcPDPtr gluon_;

  /**
   *  The cut off on pt, assuming massless quarks.
   */
  Energy pTmin_;

  //  radiative variables (pt,y)
  Energy pT_;

  /**
   *  The ParticleData objects for the fermions
   */
  vector<tcPDPtr> partons_;

  /**
   * The fermion momenta
   */
  vector<Lorentz5Momentum> quark_;

  /**
   *  The momentum of the radiated gauge boson
   */
  Lorentz5Momentum gauge_;

  /**
   *  The W boson
   */
  PPtr wboson_;

  /**
   *  W mass squared
   */
  Energy2 mw2_;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMWFermionsPOWHEGDecayer. */
template <>
struct BaseClassTrait<Herwig::SMWFermionsPOWHEGDecayer,1> {
  /** Typedef of the first base class of SMWFermionsPOWHEGDecayer. */
  typedef Herwig::SMWDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMWFermionsPOWHEGDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMWFermionsPOWHEGDecayer>
  : public ClassTraitsBase<Herwig::SMWFermionsPOWHEGDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMWFermionsPOWHEGDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SMWFermionsPOWHEGDecayer is implemented. It may also include several, space-separated,
   * libraries if the class SMWFermionsPOWHEGDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPerturbativeDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SMWFermionsPOWHEGDecayer_H */
