// -*- C++ -*-
#ifndef RADIATIVEZPRIME_MEqq2ZPrime2ZGamma2ffGamma_H
#define RADIATIVEZPRIME_MEqq2ZPrime2ZGamma2ffGamma_H
//
// This is the declaration of the MEqq2ZPrime2ZGamma2ffGamma class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.fh"

namespace RadiativeZPrime {

using namespace ThePEG;
using namespace Herwig;

/**
 * The MEqq2ZPrime2ZGamma2ffGamma class simulates \f$q\bar{q}\to Z'\to Z\gamma\f$
 * including the decay of the Z boson.
 *
 * @see \ref MEqq2ZPrime2ZGamma2ffGammaInterfaces "The interfaces"
 * defined for MEqq2ZPrime2ZGamma2ffGamma.
 */
class MEqq2ZPrime2ZGamma2ffGamma: public MEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline MEqq2ZPrime2ZGamma2ffGamma()  :_maxflavour(5) {}

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
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

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

  /**
   * Used internally by generateKinematics, after calculating the
   * limits on cos(theta).
   */
  double getCosTheta(double cthmin, double cthmax, const double r);

  /**
   *  Matrix elements for the different subprocesses
   */
  //@{
  /**
   * Matrix element for \f$q\bar{q}\to Z/\gamma g\f$.
   * @param fin   Spinors for incoming quark
   * @param ain   Spinors for incoming antiquark
   * @param gout  Polarization vectors for the outgoing photon
   * @param lm    Spinors for outgoing lepton
   * @param lp    Spinors for outgoing antilepton
   * @param me    Whether or not to calculate the matrix element for spin correlations
   **/
  InvEnergy2 qqbarME(vector<SpinorWaveFunction> & fin,
		     vector<SpinorBarWaveFunction> & ain,
		     vector<VectorWaveFunction> & gout,
		     vector<SpinorBarWaveFunction> & lm, vector<SpinorWaveFunction> & lp,
		     bool me=false) const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
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
  static ClassDescription<MEqq2ZPrime2ZGamma2ffGamma> initMEqq2ZPrime2ZGamma2ffGamma;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEqq2ZPrime2ZGamma2ffGamma & operator=(const MEqq2ZPrime2ZGamma2ffGamma &) = delete;

private:

  /**
   *  Vertices for the helicity amplitude calculation
   */
  //@{
  /**
   *  Pointer to the Z vertex
   */
  AbstractFFVVertexPtr _theFFZVertex;
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
   *  @name Pointers to the \f$Z^0\f$ and \f$Z'\f$ ParticleData objects
   */
  //@{
  /**
   *  Pointer to the Z ParticleData object
   */
  tcPDPtr _z0;

  /**
   *  Pointer to the \f$Z'\f$ ParticleData object
   */
  tcPDPtr _zPrime;
  //@}

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
  
  /**
   *  Storage of the off-shell Z mass to avoid the need to recalculate
   */
  Energy2 _mz2;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEqq2ZPrime2ZGamma2ffGamma. */
template <>
struct BaseClassTrait<RadiativeZPrime::MEqq2ZPrime2ZGamma2ffGamma,1> {
  /** Typedef of the first base class of MEqq2ZPrime2ZGamma2ffGamma. */
  typedef MEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEqq2ZPrime2ZGamma2ffGamma class and the shared object where it is defined. */
template <>
struct ClassTraits<RadiativeZPrime::MEqq2ZPrime2ZGamma2ffGamma>
  : public ClassTraitsBase<RadiativeZPrime::MEqq2ZPrime2ZGamma2ffGamma> {
  /** Return a platform-independent class name */
  static string className() { return "RadiativeZPrime::MEqq2ZPrime2ZGamma2ffGamma"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEqq2ZPrime2ZGamma2ffGamma is implemented. It may also include several, space-separated,
   * libraries if the class MEqq2ZPrime2ZGamma2ffGamma depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "RadiativeZPrime.so"; }
};

/** @endcond */

}

#endif /* RADIATIVEZPRIME_MEqq2ZPrime2ZGamma2ffGamma_H */
