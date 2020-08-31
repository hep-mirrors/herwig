// -*- C++ -*-
#ifndef Herwig_MEGammaGamma2PseudoScalar_H
#define Herwig_MEGammaGamma2PseudoScalar_H
//
// This is the declaration of the MEGammaGamma2PseudoScalar class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::VectorWaveFunction;

/**
 * The MEGammaGamma2PseudoScalar class implements the production of
 * $\pi^0$, $\eta$ and $\eta^\prime$ in photon-photon collisions
 *
 * @see \ref MEGammaGamma2PseudoScalarInterfaces "The interfaces"
 * defined for MEGammaGamma2PseudoScalar.
 */
class MEGammaGamma2PseudoScalar: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEGammaGamma2PseudoScalar() : F00_({0.274/GeV,0.274/GeV,0.344/GeV}),
				LambdaP2_({0.6*GeV2,0.6*GeV2,0.6*GeV2})
  {}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const {
    return sHat();
  }

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics();

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const {
    return 1;
  }

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

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   *  Construct the vertex of spin correlations.
   */
  virtual void constructVertex(tSubProPtr);
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

protected:

  /**
   * Matrix element for \f$\gamma\gamma\to q\bar{q}\f$
   * @param p1   The wavefunctions for the first  incoming photon
   * @param p2   The wavefunctions for the second incoming photon
   * @param momenta The momenta
   * @param calc Whether or not to calculate the matrix element
   */
  double helicityME(vector<VectorWaveFunction> &p1,vector<VectorWaveFunction> &p2,
		    const vector<Lorentz5Momentum> & momenta, int iloc, bool calc) const;

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEGammaGamma2PseudoScalar & operator=(const MEGammaGamma2PseudoScalar &) = delete;

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

private:

  /**
   *  Matrix element
   */
  mutable ProductionMatrixElement me_;

};

}

#endif /* Herwig_MEGammaGamma2PseudoScalar_H */
