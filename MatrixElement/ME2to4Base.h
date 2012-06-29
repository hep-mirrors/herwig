// -*- C++ -*-
#ifndef HERWIG_ME2to4Base_H
#define HERWIG_ME2to4Base_H
//
// This is the declaration of the ME2to4Base class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the ME2to4Base class.
 *
 * @see \ref ME2to4BaseInterfaces "The interfaces"
 * defined for ME2to4Base.
 */
class ME2to4Base: public MEBase {

public:

  /**
   * The default constructor.
   */
  ME2to4Base();

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;
  //@}

protected:

  /**
   *  Generation of the masses of the off-shell bosons
   */
  bool bosonMass(tcPDPtr boson,Energy & mb, double r,
		 Energy minMass, Energy maxMass);

  /**
   * Used internally by generateKinematics, after calculating the
   * limits on cos(theta).
   */
  virtual double getCosTheta(Energy2 sHat, Energy2 m12, Energy2 m22,
			     double cthmin, double cthmax, double r);

  /**
   *  Access to ParticleData objects
   */
  /**
   *  W plus
   */
  inline tcPDPtr WPlus() const {return _wPlus;}

  /**
   *  W minus
   */
  inline tcPDPtr WMinus() const {return _wMinus;}

  /**
   *  Z0
   */
  inline tcPDPtr Z0() const {return _z0;}

  /**
   *  Photon
   */
  inline tcPDPtr gamma() const {return _gamma;}
  //@}

  /**
   *  Access to the vertices
   */
  //@{
  /**
   *   FFPVertex
   */
  inline AbstractFFVVertexPtr vertexFFP() const {return _vertexFFP;}

  /**
   *   FFPVertex
   */
  inline AbstractFFVVertexPtr vertexFFZ() const {return _vertexFFZ;}

  /**
   *   FFPVertex
   */
  inline AbstractFFVVertexPtr vertexFFW() const {return _vertexFFW;}

  /**
   *  WWW Vertex
   */ 
  inline AbstractVVVVertexPtr vertexWWW() const {return _vertexWWW;}
  //@}

  /**
   *  set the option for the generation of the polar angle
   */
  inline void setWeightOption(unsigned int in) {_weightOpt=in;}

  /**
   *  Set the probability for generating according to \f$1/m^2\f$
   *  rather than a Breit-Wigner
   */
  inline void setSamplingProbability(double in) {_prob=in;}

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
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<ME2to4Base> initME2to4Base;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ME2to4Base & operator=(const ME2to4Base &);

private:

  /**
   *  Parameter for the sampling of the gauge boson mass
   */
  double _prob;

  /**
   *  Whether or not to generate the bosons on-shell
   */
  bool _onShell;

  /**
   *  Particle Data pointers
   */
  //@{
  /**
   *  W plus
   */
  cPDPtr _wPlus;

  /**
   *  W minus
   */
  cPDPtr _wMinus;

  /**
   *  Z0
   */
  cPDPtr _z0;

  /**
   *  Photon
   */
  cPDPtr _gamma;
  //@}

  /**
   *  Vertices
   */
  //@{
  /**
   *   FFPVertex
   */
  AbstractFFVVertexPtr _vertexFFP;

  /**
   *   FFWVertex
   */
  AbstractFFVVertexPtr _vertexFFW;

  /**
   *   FFZVertex
   */
  AbstractFFVVertexPtr _vertexFFZ;

  /**
   *  WWW Vertex
   */ 
  AbstractVVVVertexPtr _vertexWWW;
  //@}

  /**
   *  Option for the generation of the azimuth
   */
  unsigned int _weightOpt;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ME2to4Base. */
template <>
struct BaseClassTrait<Herwig::ME2to4Base,1> {
  /** Typedef of the first base class of ME2to4Base. */
  typedef MEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ME2to4Base class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ME2to4Base>
  : public ClassTraitsBase<Herwig::ME2to4Base> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ME2to4Base"; }
};

/** @endcond */

}

#endif /* HERWIG_ME2to4Base_H */
