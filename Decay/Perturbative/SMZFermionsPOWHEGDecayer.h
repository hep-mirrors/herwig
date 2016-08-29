// -*- C++ -*-
#ifndef HERWIG_SMZFermionsPOWHEGDecayer_H
#define HERWIG_SMZFermionsPOWHEGDecayer_H
//
// This is the declaration of the SMZFermionsPOWHEGDecayer class.
//

#include "SMZDecayer.h"
#include "Herwig/Utilities/Maths.h"

namespace Herwig {
  
using namespace ThePEG;


/**
 * Here is the documentation of the SMZFermionsPOWHEGDecayer class.
 *
 * @see \ref SMZFermionsPOWHEGDecayerInterfaces "The interfaces"
 * defined for SMZFermionsPOWHEGDecayer.
 */
  class SMZFermionsPOWHEGDecayer: public SMZDecayer{

public:

  /**
   * The default constructor.
   */
    SMZFermionsPOWHEGDecayer();
    
  /**
   *  Has a POWHEG style correction
   */
    virtual POWHEGType hasPOWHEGCorrection() {return FSR;}

  /**
   *  Apply the POWHEG style correction
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr);
  
  /**
   *  Virtual members to be overridden by inheriting classes
   *  which implement hard corrections 
   */
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
  static ClassDescription<SMZFermionsPOWHEGDecayer> initSMZFermionsPOWHEGDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMZFermionsPOWHEGDecayer & operator=(const SMZFermionsPOWHEGDecayer &);

  /**
   *  Pointer to the fermion-antifermion Z vertex
   */
  AbstractFFVVertexPtr FFZVertex() const {return FFZVertex_;}

  /**
   *  Pointer to the fermion-antifermion Z vertex
   */
  AbstractFFVVertexPtr FFGVertex() const {return FFGVertex_;}

  /**
   *  Real emission term, for use in generating the hardest emission
   */
  double calculateRealEmission(double x1, double x2, 
			       vector<PPtr> hardProcess,
			       double phi,
			       bool subtract,
			       int emitter) const;

  /**
   *  Real emission term, for use in generating the hardest emission
   */
  double calculateRealEmission(double x1, double x2, 
			       vector<PPtr> hardProcess,
			       double phi,
			       bool subtract) const;

  /**
   *  Check the sign of the momentum in the \f$z\f$-direction is correct.
   */
  bool checkZMomenta(double x1, double x2, double x3, double y, Energy pT) const;

  /**
   *  Calculate the Jacobian
   */
  InvEnergy calculateJacobian(double x1, double x2, Energy pT) const;

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
   *  Pointer to the fermion-antifermion Z vertex
   */
  AbstractFFVVertexPtr FFZVertex_;

  /**
   *  Pointer to the fermion-antifermion G vertex
   */
  AbstractFFVVertexPtr FFGVertex_;

  /**
   *  The colour factor 
   */
  double CF_;

  /**
   *  The Z mass
   */
  mutable Energy mZ_;

  /**
   *  The reduced mass
   */
  mutable double mu_;

  /**
   *  The square of the reduced mass
   */
  mutable double mu2_;

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
   *  The Z boson
   */
  PPtr zboson_;

  /**
   *  Higgs mass squared
   */
  Energy2 mz2_;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMZFermionsPOWHEGDecayer. */
template <>
struct BaseClassTrait<Herwig::SMZFermionsPOWHEGDecayer,1> {
  /** Typedef of the first base class of SMZFermionsPOWHEGDecayer. */
  typedef Herwig::SMZDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMZFermionsPOWHEGDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMZFermionsPOWHEGDecayer>
  : public ClassTraitsBase<Herwig::SMZFermionsPOWHEGDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMZFermionsPOWHEGDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SMZFermionsPOWHEGDecayer is implemented. It may also include several, space-separated,
   * libraries if the class SMZFermionsPOWHEGDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPerturbativeDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SMZFermionsPOWHEGDecayer_H */
