// -*- C++ -*-
#ifndef HERWIG_MEGammaGamma2WW_H
#define HERWIG_MEGammaGamma2WW_H
//
// This is the declaration of the MEGammaGamma2WW class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVVVVVertex.fh"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEGammaGamma2WW class provides the matrix elements for
 * \f$\gamma\gamma\to f \bar{f}\f$.
 *
 * @see \ref MEGammaGamma2WWInterfaces "The interfaces"
 * defined for MEGammaGamma2WW.
 */
class MEGammaGamma2WW : public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEGammaGamma2WW();

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
  virtual Energy2 scale() const;

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

  /**
   * Matrix element for \f$\gamma\gamma\to q\bar{q}\f$
   * @param p1   The wavefunctions for the first  incoming photon
   * @param p2   The wavefunctions for the second incoming photon
   * @param w1   The wavefunctions for the first  outgoing W
   * @param w2   The wavefunctions for the second outgoing W
   * @param calc Whether or not to calculate the matrix element
   */
  double helicityME(vector<VectorWaveFunction> & p1,
		    vector<VectorWaveFunction> & p2,
		    vector<VectorWaveFunction> & w1,
		    vector<VectorWaveFunction> & w2, bool calc) const;

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEGammaGamma2WW & operator=(const MEGammaGamma2WW &) = delete;

private:

  /**
   *  Treatment of the the W mass
   */
  unsigned int massOption_;

  /**
   *  Pointer to the gammaWW vertex
   */
  AbstractVVVVertexPtr WWWVertex_;

  /**
   *  Pointer to the gammagammaWW vertex
   */
  AbstractVVVVVertexPtr WWWWVertex_;

  /**
   *  Matrix element
   */
  mutable ProductionMatrixElement me_;

};

}

#endif /* HERWIG_MEGammaGamma2WW_H */
