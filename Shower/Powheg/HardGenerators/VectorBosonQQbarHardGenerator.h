// -*- C++ -*-
#ifndef HERWIG_VectorBosonQQbarHardGenerator_H
#define HERWIG_VectorBosonQQbarHardGenerator_H
//
// This is the declaration of the VectorBosonQQbarHardGenerator class.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Vectors/LorentzRotation.h"

namespace Herwig {

using namespace ThePEG;
using namespace std;
/**
 * Here is the documentation of the VectorBosonQQbarHardGenerator class.
 *
 * @see \ref VectorBosonQQbarHardGeneratorInterfaces "The interfaces"
 * defined for VectorBosonQQbarHardGenerator.
 */
class VectorBosonQQbarHardGenerator: public HardestEmissionGenerator {

public:

  /**
   * The default constructor.
   */
  VectorBosonQQbarHardGenerator() : _Ptmin(1.*GeV) {}

  /**
   *  Members which must be overridden in the inheriting classes
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr);

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
  virtual IBPtr clone() const {return new_ptr(*this);}


  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinit();
  //@}

private:

  bool getEvent();
 
  /**
   * Constructs the post-emission momenta of q, qbar, g
   */
  bool constructVectors();

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<VectorBosonQQbarHardGenerator> 
  initVectorBosonQQbarHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VectorBosonQQbarHardGenerator & operator=(const VectorBosonQQbarHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr _alphaS;

  /**
   *  Pointer to the object calculating the EM coupling
   */
  ShowerAlphaPtr _alphaEM;

  /**
   *  ParticleData object for the gluon
   */
  tcPDPtr _gluon;

  /**
   *  ParticleData object for the photon
   */
  tcPDPtr _gamma;

  /**
   *  The cut off on pt, assuming massless quarks.
   */
  Energy _Ptmin;

  /**
   *  The ParticleData objects for the fermions
   */
  vector<tcPDPtr> _partons;

  /**
   * The fermion momenta
   */
  vector<Lorentz5Momentum> _quark;

  /**
   *  The momentum of the radiated gauge boson
   */
  Lorentz5Momentum _gauge;

  /**
   *  The type of interaction
   */
  ShowerInteraction::Type _inter;



  /**
   *  The dalitz variables (xq,xqb,xg=2-xq-xqb). These are
   *  the COM energies of the q,qb,g divided by 0.5*sqrt(_s).
   */
  double _xq;
  double _xqb;
  double _xg;

  //  radiative variables (pt,y)
  double _y;
  Energy _pt;

  //the (mass dependent) shower transverse momentum
  Energy2 _showerPt;

  //com energy
  Energy2 _s;
 
  // The phi angle of the rotation of emitter-gluon plane 
  double _phi;

  // iemit = 0 quark emission: =1 antiquark emission
  int _iemitter;
  int _ispectator;

  /**
   *  The EW gauge boson
   */
  PPtr _boson;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VectorBosonQQbarHardGenerator. */
template <>
struct BaseClassTrait<Herwig::VectorBosonQQbarHardGenerator,1> {
  /** Typedef of the first base class of VectorBosonQQbarHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VectorBosonQQbarHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VectorBosonQQbarHardGenerator>
  : public ClassTraitsBase<Herwig::VectorBosonQQbarHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VectorBosonQQbarHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VectorBosonQQbarHardGenerator is implemented. 
   * It may also include several, space-separated,
   * libraries if the class VectorBosonQQbarHardGenerator depends
   * on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VectorBosonQQbarHardGenerator_H */
