// -*- C++ -*-
#ifndef Herwig_GammaGammaAmplitude_H
#define Herwig_GammaGammaAmplitude_H
//
// This is the declaration of the GammaGammaAmplitude class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "GammaGammaAmplitude.fh"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {
  
using namespace ThePEG;
using Helicity::VectorWaveFunction;
  
/**
 * The GammaGammaAmplitude class provides a base class for the implementation
 * of amplitudes for $\gamma\gamma\to X$
 *
 * @see \ref GammaGammaAmplitudeInterfaces "The interfaces"
 * defined for GammaGammaAmplitude.
 */
class GammaGammaAmplitude: public Interfaced {

public:

  /**
   * The default constructor.
   */
  GammaGammaAmplitude()
  {}
  
public:

  /** @name Virtual functions */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const = 0;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const = 0;

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim(unsigned int iopt) const = 0;

  /**
   *  The Feynman diagrams (iopt=0 gamma gamma, iopt=1, e+e-)
   */
  virtual vector<DiagPtr> getDiagrams(unsigned int iopt) const = 0;

  /**
   *  The matrix element
   */
  virtual double me2(const vector<VectorWaveFunction> & v1,
		     const vector<VectorWaveFunction> & v2,
		     const Energy2 & t1, const Energy2 & t2,
		     const Energy2 & scale, 
		     const vector<Lorentz5Momentum> & momenta,
		     const cPDVector & partons,
		     DVector & dweights) const = 0;
  
  /**
   *  The matrix element for spin correlations
   */
  virtual ProductionMatrixElement me(const vector<VectorWaveFunction> & v1,
				     const vector<VectorWaveFunction> & v2,
				     tParticleVector & particles) const = 0;

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
   * Generate the mass of the \f$\gamma\gamma\f$ system
   * @param r The random number to use for the generation of W
   * @param partons The outgoing partons
   * @param The minimum value of W
   * @param Wmax The maximum for W
   * @param jacW The jacobian for the gneeration of W
   * @param scale The scale to be used to make result dimensiionless
   * @return The generated value of W
   */
  virtual Energy generateW(double r, const tcPDVector & partons, Energy Wmin, Energy Wmax, Energy2 & jacW,Energy2 scale) = 0;
  //@}

protected:

  /**
   *   book the matrix element and indices for the loops
   */
  ProductionMatrixElement bookME(vector<unsigned int> &ihMax,
				 unsigned int ih1, unsigned ih2,
				 const vector<PDT::Spin> &spin) const;

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GammaGammaAmplitude & operator=(const GammaGammaAmplitude &) = delete;

};

}

#endif /* Herwig_GammaGammaAmplitude_H */
