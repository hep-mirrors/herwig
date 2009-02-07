// -*- C++ -*-
#ifndef HERWIG_MEPP2VV_H
#define HERWIG_MEPP2VV_H
//
// This is the declaration of the MEPP2VV class.
//

#include "Herwig++/MatrixElement/HwME2to2Base.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "Herwig++/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Models/StandardModel/StandardCKM.h"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * Here is the documentation of the MEPP2VV class.
 *
 * @see \ref MEPP2VVInterfaces "The interfaces"
 * defined for MEPP2VV.
 */
class MEPP2VV: public HwME2to2Base {

public:

  /**
   * The default constructor.
   */
  MEPP2VV();

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
  virtual Energy2 scale() const;

  /**
   * Return the process being run (WW/ZZ/WZ).
   */
  virtual int process() const { return _process; }

  /**
   * Return the factorisation scale.
   */
  virtual Energy mu_F() const { return mu_F_; }

  /**
   * Return the factorisation scale.
   */
  virtual Energy mu_UV() const { return mu_UV_; }

  /**
   * Return the maximum number of incoming flavours.
   */
  virtual int maxflavour() const { return _maxflavour; }

  /**
   * Return the CKM matrix elements.
   */
  Complex CKM(int ix,int iy) const { return _ckm[ix][iy]; }

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

  /**
   * Used internally by generateKinematics, after calculating the
   * limits on cos(theta).
   */
  virtual double getCosTheta(double ctmin, double ctmax, const double * r);
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
   * Matrix element for \f$f\bar{f}\toW^+W^-\to f\bar{f} f\bar{f}\f$.
   * @param f1  Spinors for the incoming fermion
   * @param f2  Spinors for the incoming antifermion
   * @param a1  Spinors for first  outgoing fermion
   * @param a2  Spinors for second outgoing fermion
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double helicityME(vector<SpinorWaveFunction>    & f1,
		    vector<SpinorBarWaveFunction> & a1,
		    vector<VectorWaveFunction>    & v1,
		    vector<VectorWaveFunction>    & v2,
		    bool me) const;

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
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2VV> initMEPP2VV;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2VV & operator=(const MEPP2VV &);

private:

  /**
   *  Vertices
   */
  //@{
  /**
   *   FFPVertex
   */
  AbstractFFVVertexPtr _vertexFFP;

  /**
   *   FFWVertex
   */
  AbstractFFVVertexPtr _vertexFFW;

  /**
   *   FFZVertex
   */
  AbstractFFVVertexPtr _vertexFFZ;

  /**
   *  WWW Vertex
   */ 
  AbstractVVVVertexPtr _vertexWWW;
  //@}

  /**
   * The ckm matrix elements (unsquared, to allow interference)
   */
  Complex _ckm[3][3];

  /**
   *  Processes
   */
  unsigned int _process;

  /**
   *  Allowed flavours of the incoming quarks
   */
  int _maxflavour;

  /**
   *  Processes
   */
  bool _mixingInWW;

  /**
   *  The matrix element
   */
  ProductionMatrixElement _me;

  /**
   * Selects a dynamic (sHat) or fixed factorization scale
   */
  unsigned int scaleopt_;

  /**
   * The factorization and renormalization scale respectively
   */
  Energy mu_F_, mu_UV_;

  /**
   *  Prefactor if variable scale used
   */
  double scaleFact_;

  /**
   *  Interfaced flag to turn on / off spin correlations for vector bosons.
   */
  unsigned int spinCorrelations_;

  /**
   *  Interfaced flag to invoke debugging (comparison with MCFM).
   */
  unsigned int debugMCFM_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2VV. */
template <>
struct BaseClassTrait<Herwig::MEPP2VV,1> {
  /** Typedef of the first base class of MEPP2VV. */
  typedef Herwig::HwME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2VV class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2VV>
  : public ClassTraitsBase<Herwig::MEPP2VV> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2VV"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2VV is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2VV depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2VV_H */
