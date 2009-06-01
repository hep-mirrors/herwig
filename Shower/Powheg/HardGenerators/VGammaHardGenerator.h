// -*- C++ -*-
#ifndef HERWIG_VGammaHardGenerator_H
#define HERWIG_VGammaHardGenerator_H
//
// This is the declaration of the VGammaHardGenerator class.
//

#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The VGammaHardGenerator class implements the generation of the
 * hard QCD radiation in \f$q\bar q \to W^\pm/Z^0\gamma$ events.
 *
 * @see \ref VGammaHardGeneratorInterfaces "The interfaces"
 * defined for VGammaHardGenerator.
 */
class VGammaHardGenerator: public HardestEmissionGenerator {

  /**
   * Typedef for the BeamParticleData object
   */
  typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

public:

  /**
   * The default constructor.
   */
  VGammaHardGenerator();

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
   *  Generate a \f$q \bar q \to V \gamma \f$ configuration
   */
  void generateQQbarG();

  /**
   *   The ratio of leading to NLO matrix elements
   */ 
  double QQbarGratio();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<VGammaHardGenerator> initVGammaHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VGammaHardGenerator & operator=(const VGammaHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaS_;

  /**
   *  ParticleData object of the gluon
   */
  tcPDPtr gluon_;

  /**
   *  ParticleData object of the gluon
   */
  tcPDPtr photon_;

  /**
   *  ParticleData object of the 
   */
  tcPDPtr boson_;

  /**
   *  The transverse momentum of the jet
   */
  Energy pTmin_;
private:

  /**
   *  Vertices
   */
  //@{
  /**
   *   FFPVertex
   */
  AbstractFFVVertexPtr FFPvertex_;

  /**
   *   FFWVertex
   */
  AbstractFFVVertexPtr FFWvertex_;

  /**
   *   FFZVertex
   */
  AbstractFFVVertexPtr FFZvertex_;

  /**
   *  WWW Vertex
   */ 
  AbstractVVVVertexPtr WWWvertex_;

  /**
   *  FFG Vertex
   */ 
  AbstractFFVVertexPtr FFGvertex_;
  //@}

  /**
   *  Properties of the incoming particles
   */
  //@{
  /**
   *  Pointers to the BeamParticleData objects
   */
  vector<tcBeamPtr> beams_;
  
  /**
   *  Pointers to the ParticleDataObjects for the partons
   */
  vector<tcPDPtr> partons_;
  //@}

  /**
   *  Whether the quark is in the + or - z direction
   */
  bool quarkplus_;

  /**
   *  Prefactor for the overestimate for \f$q\bar q\to V \gamma\f$
   */
  double qqgFactor_;

  /**
   * The power, \f$n\f$, for the sampling
   */
  double power_;

  /**
   *  Born variables
   */
  //@{
  /**
   *  Rapidity of the photon
   */
  double photonRapidity_;

  /**
   *  Rapidity of the gauge boson
   */
  double bosonRapidity_;

  /**
   *  \f$p_T\f$ of the photon
   */
  Energy photonpT_;

  /**
   *  Azimuth of the photon
   */
  double photonAzimuth_;

  /**
   * gauge boson mass
   */
  Energy bosonMass_;

  /**
   * Mass of the boson/photon system
   */
  Energy systemMass_;

  /**
   *  Momentum fractions for the LO process
   */ 
  double x_[2];

  /**
   * CMS energy squared of the hadron collision
   */
  Energy2  s_;

  /**
   * CMS energy of the hadron collision
   */
  Energy  rs_;
  //@}

  /**
   *  Momenta etc for the \f$q\bar q \to V \gamma g \f$ process
   */
  //@{
  /**
   *  Momentum of the vector boson
   */
  Lorentz5Momentum pVqqbar_;

  /**
   *  Momentum of the photon
   */
  Lorentz5Momentum pGammaqqbar_;

  /**
   *  Momentum of the gluon
   */
  Lorentz5Momentum pGqqbar_;

  /**
   *  Momentum of the incoming quark
   */
  Lorentz5Momentum pQqqbar_;

  /**
   *  Momentum of the incoming antiquark
   */
  Lorentz5Momentum pQbarqqbar_;

  /**
   *  The transverse momentum
   */
  Energy pTqqbar_;
  //@}

  /**
   *  Momenta etc for the \f$qg  \to V \gamma q \f$ process
   */
  //@{
  /**
   *  Momentum of the vector boson
   */
  Lorentz5Momentum pVqg_;

  /**
   *  Momentum of the photon
   */
  Lorentz5Momentum pGammaqg_;

  /**
   *  Momentum of the gluon
   */
  Lorentz5Momentum pGqg_;

  /**
   *  Momentum of the incoming quark
   */
  Lorentz5Momentum pQinqg_;

  /**
   *  Momentum of the incoming antiquark
   */
  Lorentz5Momentum pQoutqg_;

  /**
   *  The transverse momentum
   */
  Energy pTqg_;
  //@}

  /**
   *  Momenta etc for the \f$g\bar q \to V \gamma \bar q \f$ process
   */
  //@{
  /**
   *  Momentum of the vector boson
   */
  Lorentz5Momentum pVgqbar_;

  /**
   *  Momentum of the photon
   */
  Lorentz5Momentum pGammagqbar_;

  /**
   *  Momentum of the gluon
   */
  Lorentz5Momentum pGgqbar_;

  /**
   *  Momentum of the incoming quark
   */
  Lorentz5Momentum pQingqbar_;

  /**
   *  Momentum of the incoming antiquark
   */
  Lorentz5Momentum pQoutgqbar_;

  /**
   *  The transverse momentum
   */
  Energy pTgqbar_;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VGammaHardGenerator. */
template <>
struct BaseClassTrait<Herwig::VGammaHardGenerator,1> {
  /** Typedef of the first base class of VGammaHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VGammaHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VGammaHardGenerator>
  : public ClassTraitsBase<Herwig::VGammaHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VGammaHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VGammaHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class VGammaHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VGammaHardGenerator_H */
