// -*- C++ -*-
#ifndef HERWIG_MENeutralCurrentDIS_H
#define HERWIG_MENeutralCurrentDIS_H
//
// This is the declaration of the MENeutralCurrentDIS class.
//

#include "DISBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MENeutralCurrentDIS class provides the matrix elements for
 * neutral current DIS.
 *
 *  For consistency both the incoming and outgoing quarks are assumed to be massless.
 *
 * @see \ref MENeutralCurrentDISInterfaces "The interfaces"
 * defined for MENeutralCurrentDIS.
 */
class MENeutralCurrentDIS: public DISBase {

public:

  /**
   * The default constructor.
   */
  MENeutralCurrentDIS();

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
   * Matrix element for \f$\ell q\to \gamma/Z \to \ell q\f$.
   * @param f1 Fermion on lepton line
   * @param a1 Anti-fermion on lepton line
   * @param f2 Fermion on quark line
   * @param a2 Anti-fermion on quark line
   * @param lorder The order of particles on the lepton line
   * @param qorder The order of particles on the quark line
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double helicityME(vector<SpinorWaveFunction>    & f1 ,
		    vector<SpinorWaveFunction>    & f2,
		    vector<SpinorBarWaveFunction> & a1 ,
		    vector<SpinorBarWaveFunction> & a2,
		    bool lorder, bool qorder,
		    bool me) const;


  /**
   *  Option for treatment of \f$\gamma/Z\f$ terms
   */
  inline unsigned int gammaZOption() const {return _gammaZ;}

  /**
   *  Calculate the coefficient A for the correlations in the hard
   *  radiation
   */
  virtual double A(tcPDPtr lin, tcPDPtr lout, tcPDPtr qin, tcPDPtr qout,
		   Energy2 scale) const;

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MENeutralCurrentDIS> initMENeutralCurrentDIS;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MENeutralCurrentDIS & operator=(const MENeutralCurrentDIS &) = delete;

private:

  /**
   *  Pointer to the vertices for the helicity calculations
   */
  //@{
  /**
   *  Pointer to the Z vertex
   */
  AbstractFFVVertexPtr _theFFZVertex;

  /**
   *  Pointer to the photon vertex
   */
  AbstractFFVVertexPtr _theFFPVertex;
  //@}

  /**
   *  Pointers to the intermediate resonances
   */
  //@{
  /**
   *  Pointer to the Z ParticleData object
   */
  tcPDPtr _z0;

  /**
   *  Pointer to the photon ParticleData object
   */
  tcPDPtr _gamma;
  //@}

  /**
   *  Switches to control the particles in the hard process
   */
  //@{
  /**
   *  Minimumflavour of the incoming quarks
   */
  int _minflavour;

  /**
   *  Maximum flavour of the incoming quarks
   */
  int _maxflavour;

  /**
   *  Whether to include both \f$Z^0\f$ and \f$\gamma\f$ or only one
   */
  unsigned int _gammaZ;
  //@}

  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement _me;

  /**
   *  Electroweak parameters
   */
  //@{
  /**
   *  \f$\sin\theta_W\f$
   */
  double _sinW;

  /**
   *  \f$\cos\theta_W\f$
   */
  double _cosW;

  /**
   *  The square of the Z mass
   */
  Energy2 _mz2;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MENeutralCurrentDIS. */
template <>
struct BaseClassTrait<Herwig::MENeutralCurrentDIS,1> {
  /** Typedef of the first base class of MENeutralCurrentDIS. */
  typedef Herwig::DISBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MENeutralCurrentDIS class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MENeutralCurrentDIS>
  : public ClassTraitsBase<Herwig::MENeutralCurrentDIS> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MENeutralCurrentDIS"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MENeutralCurrentDIS is implemented. It may also include several, space-separated,
   * libraries if the class MENeutralCurrentDIS depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEDIS.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MENeutralCurrentDIS_H */
