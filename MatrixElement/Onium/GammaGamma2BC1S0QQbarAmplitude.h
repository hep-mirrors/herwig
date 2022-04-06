// -*- C++ -*-
#ifndef Herwig_GammaGamma2BC1S0QQbarAmplitude_H
#define Herwig_GammaGamma2BC1S0QQbarAmplitude_H
//
// This is the declaration of the GammaGamma2BC1S0QQbarAmplitude class.
//

#include "Herwig/MatrixElement/Gamma/GammaGammaAmplitude.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "OniumParameters.h"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * Here is the documentation of the GammaGamma2BC1S0QQbarAmplitude class.
 *
 * @see \ref GammaGamma2BC1S0QQbarAmplitudeInterfaces "The interfaces"
 * defined for GammaGamma2BC1S0QQbarAmplitude.
 */
class GammaGamma2BC1S0QQbarAmplitude: public GammaGammaAmplitude {

public:

  /**
   * The default constructor.
   */
  GammaGamma2BC1S0QQbarAmplitude() :O1_(ZERO), n_(1)
  {}

public:

  /** @name Virtual functions required by GammaGammaAmplitude class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const {
    return 2;
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
    return 5+iopt;
  }

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(unsigned int iopt,  const cPDVector & partons, tcDiagPtr diag) const;

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

  /**
   *  The Feynman diagrams (iopt=0 gamma gamma, iopt=1, e+e-)
   */
  virtual vector<DiagPtr> getDiagrams(unsigned int iopt) const;

  /**
   *  The matrix element
   */
  virtual double me2(const vector<VectorWaveFunction> & v1,
		     const vector<VectorWaveFunction> & v2,
		     const Energy2 & , const Energy2 & ,
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
    vector<SpinorWaveFunction> v4;
    vector<SpinorBarWaveFunction>  ubar5;
    assert(abs(particles[0]->id())==541);
    assert(particles[1]->id()<0);
    assert(particles[2]->id()>0);
    SpinorWaveFunction   (v4   ,particles[1],outgoing,true );
    SpinorBarWaveFunction(ubar5,particles[2],outgoing,true );
    vector<Lorentz5Momentum> pout = {particles[0]->momentum(),
				     particles[1]->momentum(),
				     particles[2]->momentum()};
    Energy2 scale = (v1[0].momentum()+v2[0].momentum()).m2();
    DVector dwgt;
    double output(0);
    return helicityAmplitude(scale,v1,v2,pout,v4,ubar5,output,dwgt);
  }
  
  /**
   * Generate the mass of the \f$\gamma\gamma\f$ system
   */
  virtual Energy generateW(double r, const tcPDVector & partons,
			   Energy Wmax, Energy2 & jacW, Energy2 scale);

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

protected:

  /**
   *  Calculation of the helicity amplitudes for the process
   */
  ProductionMatrixElement helicityAmplitude(const Energy2 & scale,
					    const vector<VectorWaveFunction> & v1,
  					    const vector<VectorWaveFunction> & v2,
  					    const vector<Lorentz5Momentum> & out,
					    const vector<SpinorWaveFunction> v4,
					    const vector<SpinorBarWaveFunction> ubar5,
  					    double & output, DVector & dweights) const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GammaGamma2BC1S0QQbarAmplitude & operator=(const GammaGamma2BC1S0QQbarAmplitude &) = delete;

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
   *  Principal quantum number
   */
  unsigned int n_;
  //@}

};

}

#endif /* Herwig_GammaGamma2BC1S0QQbarAmplitude_H */
