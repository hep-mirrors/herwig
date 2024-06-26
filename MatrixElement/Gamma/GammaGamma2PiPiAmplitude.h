// -*- C++ -*-
#ifndef Herwig_GammaGamma2PiPiAmplitude_H
#define Herwig_GammaGamma2PiPiAmplitude_H
//
// This is the declaration of the GammaGamma2PiPiAmplitude class.
//

#include "GammaGammaAmplitude.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The GammaGamma2PiPiAmplitude class implements the amplitude for \f$\gamma\gamma\to\pi^+\pi^-\f$.
 *
 * @see \ref GammaGamma2PiPiAmplitudeInterfaces "The interfaces"
 * defined for GammaGamma2PiPiAmplitude.
 */
class GammaGamma2PiPiAmplitude: public GammaGammaAmplitude {

public:
  
  /**
   * The default constructor.
   */
  GammaGamma2PiPiAmplitude() : mode_(0)
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
  virtual int nDim(unsigned int iopt) const {
    return iopt+2;
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
    ScalarWaveFunction(particles[1],outgoing,true);
    vector<Lorentz5Momentum> pout = {particles[0]->momentum(),
				     particles[1]->momentum()};
    double output(0);
    DVector dwgt;
    return helicityAmplitude(v1,v2,pout,output,dwgt);
  }

  /**
   * Generate the mass of the \f$\gamma\gamma\f$ system
   */
  virtual Energy generateW(double , const tcPDVector & partons, Energy Wmin, Energy Wmax,Energy2 & jacW, Energy2 scale);
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

  /**
   *  Calculation of the helicity amplitudes for the process
   */
  ProductionMatrixElement helicityAmplitude(const vector<VectorWaveFunction> & v1,
					    const vector<VectorWaveFunction> & v2,
					    const vector<Lorentz5Momentum> & out,
					    double & output, DVector & dweights) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GammaGamma2PiPiAmplitude & operator=(const GammaGamma2PiPiAmplitude &) = delete;

private:

  /**
   *   Which particles to produce
   */
  unsigned int mode_;

};

}

#endif /* Herwig_GammaGamma2PiPiAmplitude_H */
