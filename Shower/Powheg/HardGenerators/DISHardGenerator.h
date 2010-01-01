// -*- C++ -*-
#ifndef HERWIG_DISHardGenerator_H
#define HERWIG_DISHardGenerator_H
//
// This is the declaration of the DISHardGenerator class.
//

#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "ThePEG/PDF/BeamParticleData.h"

namespace Herwig {

using namespace ThePEG;

/**
 *  Typedef for BeamParticleData pointers
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

/**
 * The DISHardGenerator class implements the generation of the hardest emission 
 * in the POWHEG approach for the DIS process.
 *
 * @see \ref DISHardGeneratorInterfaces "The interfaces"
 * defined for DISHardGenerator.
 */
class DISHardGenerator: public HardestEmissionGenerator {

public:

  /**
   * The default constructor.
   */
  DISHardGenerator();

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

protected:

  /**
   *  Calculate the coefficient A for the correlations
   */
  inline double A(tcPDPtr qin, tcPDPtr qout, tcPDPtr lin, tcPDPtr lout);

  /**
   *  Generate a Compton process
   */
  void generateCompton();

  /**
   *  Matrix element piece for the Compton process
   */
  double comptonME(double xT,double xp, double zp, double phi);

  /**
   *  Generate a BGF process
   */
  void generateBGF();

  /**
   *  Matrix element piece for the Compton process
   */
  double BGFME(double xT,double xp, double zp, double phi);


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DISHardGenerator> initDISHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DISHardGenerator & operator=(const DISHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaS_;

  /**
   *  Weight for the compton channel
   */
  double comptonWeight_;

  /**
   *  Weight for the BGF channel
   */
  double BGFWeight_;

  /**
   *  Minimum value of \f$p_T\f$
   */
  Energy pTmin_;

  /**
   *  Beam particle
   */
  tcBeamPtr beam_;

  /**
   *  Partons
   */
  tcPDPtr partons_[2];

  /**
   *  PDF object
   */
  tcPDFPtr pdf_;

  /**
   *  Rotation to the Breit frame
   */
  LorentzRotation rot_;

  /**
   *  Electroweak parameters
   */
  //@{
  /**
   *  \f$\sin\theta_W\f$
   */
  double sinW_;

  /**
   *  \f$\cos\theta_W\f$
   */
  double cosW_;

  /**
   *  The square of the Z mass
   */
  Energy2 mz2_;

  /**
   *  The coefficient for the correlations
   */
  double acoeff_;
  //@}

  /**
   *  Gluon particle data object
   */
  PDPtr gluon_;

  /**
   *  \f$Q^2\f$
   */
  Energy2 q2_;

  /**
   *  
   */
  double xB_;

  /**
   *
   */
  double l_;

  /**
   *  Lepton momenta
   */
  Lorentz5Momentum pl_[2];

  /**
   *  Quark momenta
   */
  Lorentz5Momentum pq_[2];

  /**
   *  q
   */
  Lorentz5Momentum q_;

  /**
   *  Compton parameters
   */
  Energy pTCompton_;
  bool ComptonISFS_;
  vector<Lorentz5Momentum> ComptonMomenta_;

  /**
   *  BGF parameters
   */
  Energy pTBGF_;
  vector<Lorentz5Momentum> BGFMomenta_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DISHardGenerator. */
template <>
struct BaseClassTrait<Herwig::DISHardGenerator,1> {
  /** Typedef of the first base class of DISHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DISHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DISHardGenerator>
  : public ClassTraitsBase<Herwig::DISHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DISHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DISHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class DISHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DISHardGenerator_H */
