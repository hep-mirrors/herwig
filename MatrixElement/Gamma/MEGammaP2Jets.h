// -*- C++ -*-
#ifndef HERWIG_MEGammaP2Jets_H
#define HERWIG_MEGammaP2Jets_H
//
// This is the declaration of the MEGammaP2Jets class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEGammaP2Jets class implements the matrix elements for 
 * pointlike gamma+hadron -> jets.
 *
 * @see \ref MEGammaP2JetsInterfaces "The interfaces"
 * defined for MEGammaP2Jets.
 */
class MEGammaP2Jets: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEGammaP2Jets();

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
  //@}

protected:

  /**
   *  Members to calculate the matrix elements
   */
  //@{
  /**
   * Matrix element for \f$\gamma g\to q \bar{q}\f$.
   * @param gmin Polarization vectors for the incoming photon 
   * @param glin Polarization vectors for the incoming gluon
   * @param fout Spinors for the outgoing quark
   * @param aout Spinors for the outgoing antiquark
   * @param calc  Whether or not to calculate the matrix element for spin correlations
   */
  double gammagluonME(vector<VectorWaveFunction> & gmin,
		      vector<VectorWaveFunction> & glin,
		      vector<SpinorBarWaveFunction> & fout, 
		      vector<SpinorWaveFunction> & aout,
		      bool calc) const;

  /**
   * Matrix element for \f$\gamma q\to g q\f$.
   * @param gmin Polarization vectors for the incoming photon
   * @param fin  Spinors for the incoming quark
   * @param gout Polarization vectors for the outgong gluon
   * @param fout Spinors for the outgoing quark
   * @param calc  Whether or not to calculate the matrix element for spin correlations
   */
  double gammaquarkME(vector<VectorWaveFunction> & gmin,
		      vector<SpinorWaveFunction> & fin,
		      vector<VectorWaveFunction> & gout,
		      vector<SpinorBarWaveFunction> & fout,
		      bool calc) const;

  /**
   * Matrix element for \f$\gamma q\to g q\f$.
   * @param gmin Polarization vectors for the incoming photon
   * @param fin  Spinors for the incoming antiquark
   * @param gout Polarization vectors for the outgong gluon
   * @param fout Spinors for the outgoing antiquark
   * @param calc  Whether or not to calculate the matrix element for spin correlations
   */
  double gammaantiquarkME(vector<VectorWaveFunction> & gmin,
			  vector<SpinorBarWaveFunction> & fin,
			  vector<VectorWaveFunction> & gout,
			  vector<SpinorWaveFunction> & fout,
			  bool calc) const;
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
  static ClassDescription<MEGammaP2Jets> initMEGammaP2Jets;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEGammaP2Jets & operator=(const MEGammaP2Jets &);

private:

  /**
   *  Pointer to the quark-antiquark-gluon vertex
   */
  AbstractFFVVertexPtr _gluonvertex;

  /**
   *  Pointer to the quark-antiquark-photon vertex
   */
  AbstractFFVVertexPtr _photonvertex;

  /**
   *  Allowed processes
   */
  unsigned int _process;

  /**
   *  Minimum flavour
   */
  int _minflavour;

  /**
   *  Maximum flavour
   */
  int _maxflavour;
  
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
 *  base classes of MEGammaP2Jets. */
template <>
struct BaseClassTrait<Herwig::MEGammaP2Jets,1> {
  /** Typedef of the first base class of MEGammaP2Jets. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEGammaP2Jets class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEGammaP2Jets>
  : public ClassTraitsBase<Herwig::MEGammaP2Jets> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEGammaP2Jets"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEGammaP2Jets is implemented. It may also include several, space-separated,
   * libraries if the class MEGammaP2Jets depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEGammaHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEGammaP2Jets_H */
