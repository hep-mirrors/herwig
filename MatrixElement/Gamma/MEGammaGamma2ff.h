// -*- C++ -*-
#ifndef HERWIG_MEGammaGamma2ff_H
#define HERWIG_MEGammaGamma2ff_H
//
// This is the declaration of the MEGammaGamma2ff class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEGammaGamma2ff class provides the matrix elements for
 * \f$\gamma\gamma\to f \bar{f}\f$.
 *
 * @see \ref MEGammaGamma2ffInterfaces "The interfaces"
 * defined for MEGammaGamma2ff.
 */
class MEGammaGamma2ff: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEGammaGamma2ff();

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
   * @param f    The wavefunction  for the outgoing fermion
   * @param fbar The wavefunction  for the outgoing antifermion
   * @param calc Whether or not to calculate the matrix element
   */
  double helicityME(vector<VectorWaveFunction> &p1,vector<VectorWaveFunction> &p2,
		    vector<SpinorBarWaveFunction> & f,
		    vector<SpinorWaveFunction> & fbar, bool calc) const;
protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEGammaGamma2ff> initMEGammaGamma2ff;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEGammaGamma2ff & operator=(const MEGammaGamma2ff &);

private:
  
  /**
   *  Which processes to include
   */
  int process_;

  /**
   *  Pointer to the photon vertex
   */
  AbstractFFVVertexPtr vertex_;

  /**
   *  Matrix element
   */
  mutable ProductionMatrixElement me_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEGammaGamma2ff. */
template <>
struct BaseClassTrait<Herwig::MEGammaGamma2ff,1> {
  /** Typedef of the first base class of MEGammaGamma2ff. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEGammaGamma2ff class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEGammaGamma2ff>
  : public ClassTraitsBase<Herwig::MEGammaGamma2ff> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEGammaGamma2ff"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEGammaGamma2ff is implemented. It may also include several, space-separated,
   * libraries if the class MEGammaGamma2ff depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEGammaGamma.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEGammaGamma2ff_H */
