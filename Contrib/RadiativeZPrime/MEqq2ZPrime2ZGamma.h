// -*- C++ -*-
#ifndef RADIATIVEZPRIME_MEqq2ZPrime2ZGammma_H
#define RADIATIVEZPRIME_MEqq2ZPrime2ZGamma_H
//
// This is the declaration of the MEqq2ZPrime2ZGamma class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.fh"

namespace RadiativeZPrime {

using namespace ThePEG;
using namespace Herwig;

/**
 * The MEqq2ZPrime2ZGamma class implements the matrix element for
 * \f$q\bar{q}\to Z'\to Z\gamma\f$.
 *
 * @see \ref MEqq2ZPrime2ZGammaInterfaces "The interfaces"
 * defined for MEqq2ZPrime2ZGamma.
 */
class MEqq2ZPrime2ZGamma: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEqq2ZPrime2ZGamma();

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

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const { return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const { return new_ptr(*this);}
  //@}

protected:

  /**
   * Matrix element for \f$q\bar{q}\to \gamma/Z \to f\bar{f}\f$.
   * @param fin  Spinors for incoming quark
   * @param ain  Spinors for incoming antiquark
   * @param Zout Polarization vectors for outgoing Z
   * @param gammaout Polarization vectors for the the outgoing photon
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double qqME(vector<SpinorWaveFunction>    & fin ,
		 vector<SpinorBarWaveFunction> & ain ,
		 vector<VectorWaveFunction> & Zout,
		 vector<VectorWaveFunction> & gammaout,
		 bool me) const;

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
  static ClassDescription<MEqq2ZPrime2ZGamma> initMEqq2ZPrime2ZGamma;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEqq2ZPrime2ZGamma & operator=(const MEqq2ZPrime2ZGamma &);

private:

  /**
   *  Pointer to the vertices for the helicity calculations
   */
  //@{
  /**
   *  Pointer to the \f$Z'\f$ vertex
   */
  AbstractFFVVertexPtr _theFFZPrimeVertex;

  /**
   *  Pointer to the \f$\gamma Z' Z vertex\f$
   */
  AbstractVVVVertexPtr _theGammaZPrimeZVertex;
  //@}

  /**
   *  Pointer to the \f$Z'\f$ ParticleData object
   */
  tcPDPtr _zPrime;

  /**
   *  Switches to control the particles in the hard process
   */
  //@{
  /**
   *  Allowed flavours for the incoming quarks
   */
  unsigned int _maxflavour;
  //@}

  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement _me;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEqq2ZPrime2ZGamma. */
template <>
struct BaseClassTrait<RadiativeZPrime::MEqq2ZPrime2ZGamma,1> {
  /** Typedef of the first base class of MEqq2ZPrime2ZGamma. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEqq2ZPrime2ZGamma class and the shared object where it is defined. */
template <>
struct ClassTraits<RadiativeZPrime::MEqq2ZPrime2ZGamma>
  : public ClassTraitsBase<RadiativeZPrime::MEqq2ZPrime2ZGamma> {
  /** Return a platform-independent class name */
  static string className() { return "RadiativeZPrime::MEqq2ZPrime2ZGamma"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEqq2ZPrime2ZGamma is implemented. It may also include several, space-separated,
   * libraries if the class MEqq2ZPrime2ZGamma depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "RadiativeZPrime.so"; }
};

/** @endcond */

}

#endif /* RADIATIVEZPRIME_MEqq2ZPrime2ZGamma_H */
