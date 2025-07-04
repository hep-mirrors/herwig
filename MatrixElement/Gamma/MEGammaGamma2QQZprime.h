// -*- C++ -*-
#ifndef HERWIG_MEGammaGamma2QQZprime_H
#define HERWIG_MEGammaGamma2QQZprime_H

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "Herwig/Models/HiddenValley/HiddenValleyModel.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "ThePEG/PDT/ParticleData.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Internal matrix element for gamma gamma -> j j ZPrime using the Hidden Valley model.
 *
 * Allows flavor-summed jet production from two incoming photons with
 * spin correlations and helicity wavefunction support.
 */
class MEGammaGamma2QQZprime : public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEGammaGamma2QQZprime();

  /**
   * Return matrix element squared
   */
  virtual double me2() const;

  /**
   * Return scale (used for alpha_s or alpha_em)
   */
  virtual Energy2 scale() const;

  /**
   * Generate momenta for gamma gamma -> j j ZPrime
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * Set internal kinematics (momenta, sHat, etc.)
   */
  virtual void setKinematics();

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * No specific diagram structure is required.
   */
  virtual void getDiagrams() const;

  /**
   * Construct spin-correlation-aware vertex
   */
  virtual void constructVertex(tSubProPtr sub);

  /**
   * Initialize class description and interfaces
   */
  static void Init();

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
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *> colourGeometries(tcDiagPtr) const;

  /**
   * Add all possible diagrams with the add() function.
   */
  Selector<MEBase::DiagramIndex> diagrams(const DiagramVector & diags) const;

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

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans)
   ;

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}


protected:

  /**
   * Clone this object
   */
  virtual IBPtr clone() const;

  /**
   * Full clone
   */
  virtual IBPtr fullclone() const;

  /**
   * Initialization routine
   */
  virtual void doinit();

  /**
   *  Members to calculate the matrix elements
   */
  //@{
  /**
   * Matrix element for \f$gamma gamma -> q qbar ZPrime\f$
   * @param g1        The wavefunctions for the first  incoming gluon
   * @param g2        The wavefunctions for the second incoming gluon
   * @param q         The wavefunction  for the outgoing quark
   * @param qbar      The wavefunction  for the outgoing antiquark
   * @param zprime    The wavefunction for the outgoing ZPrime
   * @param flow      The colour flow
   */
   double gammagammaME(vector<VectorWaveFunction> &g1, vector<VectorWaveFunction> &g2,
     vector<SpinorBarWaveFunction> &qbar, vector<SpinorWaveFunction> &q,
     VectorWaveFunction &zprime, unsigned int iflow) const;

private:

  // Prevent copying
  MEGammaGamma2QQZprime & operator=(const MEGammaGamma2QQZprime &) = delete;

private:

  /** @name Configuration parameters */
  //@{
  /**
   * Minimum and maximum quark flavors (e.g., 1 = d, 6 = t)
   */
  int minFlavor_, maxFlavor_;

  /**
  * Scale option for gamma gamma -> j j ZPrime matrix element
  */
  unsigned int scaleOption_;

  /**
   * Pointer to the ZPrime particle (should match PDG code 32 by default)
   */
  PDPtr Zprime_;

  /**
   *  On-shell mass for the ZPrime particle
   */
  Energy mzp_;

  /**
   *  On-shell width for the ZPrime particle
   */
  Energy wzp_;

  /**
   * Pointer to HiddenValleyModel instance
   */
  tcHiddenValleyPtr model_;

  /**
   * Pointer to the ZPrime vertex (FFV coupling to fermions)
   */
  AbstractFFVVertexPtr ZprimeVertex_;
  //@}

  /** @name Matrix element state */
  //@{
  /**
   * Pointer to selected quark flavors (PDPtr objects)
   * used for persistency or interface hooks (may remain unset)
   */
  PDPtr quark1_, quark2_, photon_;

  /**
   * External momenta for each particle
   */
  Lorentz5Momentum p1_, p2_, p3_, p4_, p5_;

  /**
   * Mandelstam variables
   */
  Energy2 sHat_, tHat_, uHat_;

  /**
   *  Colour flow
   */
  mutable unsigned int flow_;

  /**
   *  Diagram
   */
  mutable unsigned int diagram_;

  /**
   * ME object for spin-correlation-aware vertex (optional)
   */
  mutable ProductionMatrixElement me_;
  //@}
};

}

#endif // HERWIG_MEGammaGamma2QQZprime_H
