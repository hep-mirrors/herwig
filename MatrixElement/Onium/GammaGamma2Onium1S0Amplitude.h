// -*- C++ -*-
#ifndef Herwig_GammaGamma2Onium1S0Amplitude_H
#define Herwig_GammaGamma2Onium1S0Amplitude_H
//
// This is the declaration of the GammaGamma2Onium1S0Amplitude class.
//

#include "Herwig/MatrixElement/Gamma/GammaGammaAmplitude.h"
#include "OniumParameters.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The GammaGamma2Onium1S0Amplitude class implements the matrix element for \f$\gamma\gamma\to^{1}\!\!S_0\f$ quarkonium.
 * states.
 *
 * @see \ref GammaGamma2Onium1S0AmplitudeInterfaces "The interfaces"
 * defined for GammaGamma2Onium1S0Amplitude.
 */
class GammaGamma2Onium1S0Amplitude: public GammaGammaAmplitude {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  GammaGamma2Onium1S0Amplitude() : O1_(ZERO), state_(ccbar), n_(1),
				   Lambda2_(sqr(3.0969*GeV)), mOpt_(0)
  {}
  //@}

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
    return mOpt_;
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
		     const cPDVector & , DVector & ) const {
    Energy M  = momenta.back().mass();
    double output(0.);
    helicityAmplitude(v1,v2,M,output);
    // coupling factors
    double eQ = state_==ccbar ? 2./3. : -1./3.;
    double alpha = generator()->standardModel()->alphaEM();
    return 128.*output*O1_/M/scale*sqr(Constants::pi*alpha*sqr(eQ)/(1.-t1/Lambda2_)/(1.-t2/Lambda2_));
  }
  
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

  /**
   * Generate the mass of the \f$\gamma\gamma\f$ system
   */
  virtual Energy generateW(double r, const tcPDVector & partons, Energy Wmin,
			   Energy Wmax, Energy2 & jacW, Energy2 scale);

  /**
   * Generate internal degrees of freedom given 'nDim()' uniform
   * random numbers in the interval ]0,1[. To help the phase space
   * generator, the 'dSigHatDR()' should be a smooth function of these
   * numbers, although this is not strictly necessary. Return
   * false if the chosen points failed the kinematical cuts.
   */
  virtual double generateKinematics(const double * r,
				    const Energy2 & scale, 
				    vector<Lorentz5Momentum> & momenta,
				    const tcPDVector & partons);

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
  GammaGamma2Onium1S0Amplitude & operator=(const GammaGamma2Onium1S0Amplitude &) = delete;

private:

  /**
   *   Parameters for the form-factors
   */
  //@{
  /**
   *  Access to the parameters for the quarkonium states
   */
  OniumParametersPtr params_;
  
  /**
   *  The \f$O_1\f$ colour-singlet coefficient
   */
  Energy3 O1_;

  /**
   *  Type of state
   */
  OniumState state_;

  /**
   *  Principal quantum number
   */
  unsigned int n_;
  /**
   * Pole mass squared parameter for the form factors
   */
  Energy2 Lambda2_;

  /**
   *  Option for the mass generation
   */
  unsigned int mOpt_;

  /**
   *  The mass generator for the onium state
   */
  GenericMassGeneratorPtr massGen_;
  //@}

};

}

#endif /* Herwig_GammaGamma2Onium1S0Amplitude_H */
