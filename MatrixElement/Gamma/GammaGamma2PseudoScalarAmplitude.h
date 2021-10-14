// -*- C++ -*-
#ifndef Herwig_GammaGamma2PseudoScalarAmplitude_H
#define Herwig_GammaGamma2PseudoScalarAmplitude_H
//
// This is the declaration of the GammaGamma2PseudoScalarAmplitude class.
//

#include "GammaGammaAmplitude.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the GammaGamma2PseudoScalarAmplitude class.
 *
 * @see \ref GammaGamma2PseudoScalarAmplitudeInterfaces "The interfaces"
 * defined for GammaGamma2PseudoScalarAmplitude.
 */
class GammaGamma2PseudoScalarAmplitude: public GammaGammaAmplitude {

public:

  /**
   * The default constructor.
   */
  GammaGamma2PseudoScalarAmplitude() : F00_({0.274/GeV,0.274/GeV,0.344/GeV}),
				       LambdaP2_({0.6*GeV2,0.6*GeV2,0.6*GeV2})
  {}

public:

  /** @name Virtual functions required by GammaGammaAmplitude class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const {
    return 0;
  }

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const {
    return 2;
  }

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim(unsigned int) const {
    return 0;
  }

  /**
   *  The Feynman diagrams (iopt=0 gamma gamma, iopt=1, e+e-)
   */
  virtual vector<DiagPtr> getDiagrams(unsigned int iopt) const;

  /**
   *  The matrix element
   */
  virtual double me2(const vector<VectorWaveFunction> & v1,
		     const vector<VectorWaveFunction> & v2,
		     const Energy2 & t1, const Energy2 & t2,
		     const Energy2 & scale, 
		     const vector<Lorentz5Momentum> & momenta,
		     const cPDVector & partons,
		     DVector & dweights ) const; 
  
  /**
   * Matrix element for spin correlations
   */
  virtual ProductionMatrixElement me(const vector<VectorWaveFunction> & v1,
				     const vector<VectorWaveFunction> & v2,
				     tParticleVector & particles) const {
    ScalarWaveFunction(particles[0],outgoing,true);
    double output(0);
    return helicityAmplitude(v1,v2,particles[0]->mass(),output);
  }

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

  /**
   *  Calculation of the helicity amplitudes for the process
   */
  ProductionMatrixElement helicityAmplitude(const vector<VectorWaveFunction> & v1,
					    const vector<VectorWaveFunction> & v2,
					    const Energy & M, double & output) const;
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GammaGamma2PseudoScalarAmplitude & operator=(const GammaGamma2PseudoScalarAmplitude &) = delete;

private:

  /**
   *   Parameters for the form-factors
   */
  //@{
  /**
   *    Form factors at $Q_1^2=Q_2^2=0$
   */
  vector<InvEnergy> F00_;

  /**
   * Pole mass squared parameter for the form factors
   */
  vector<Energy2> LambdaP2_;
  //@}

};

}

#endif /* Herwig_GammaGamma2PseudoScalarAmplitude_H */
