// -*- C++ -*-
#ifndef HERWIG_MEee2ZZ_H
#define HERWIG_MEee2ZZ_H
//
// This is the declaration of the MEee2ZZ class.
//

#include "Herwig++/MatrixElement/ME2to4Base.h"
#include "Herwig++/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEee2ZZ class simulates \f$e^+e^-\to Z^0Z^0\f$ including
 * the decay products of the \f$Z^0\f$ bosons.
 *
 * @see \ref MEee2ZZInterfaces "The interfaces"
 * defined for MEee2ZZ.
 */
class MEee2ZZ: public ME2to4Base {

public:

  /**
   * The default constructor.
   */
  MEee2ZZ();

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

public:

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
   * @param a2  Spinors for first  outgoing fermion
   * @param a4  Spinors for second outgoing fermion
   * @param a1  Spinors for first  outgoing antifermion
   * @param a2  Spinors for second outgoing antifermion
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double helicityME(vector<SpinorWaveFunction> &    f1 ,
		    vector<SpinorBarWaveFunction> & a1 ,
		    vector<SpinorWaveFunction> &    f2 ,
		    vector<SpinorBarWaveFunction> & a2 ,
		    vector<SpinorWaveFunction> &    f3 ,
		    vector<SpinorBarWaveFunction> & a4 ,
		    bool me) const;

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
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<MEee2ZZ> initMEee2ZZ;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2ZZ & operator=(const MEee2ZZ &);

private:

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
 *  base classes of MEee2ZZ. */
template <>
struct BaseClassTrait<Herwig::MEee2ZZ,1> {
  /** Typedef of the first base class of MEee2ZZ. */
  typedef Herwig::ME2to4Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEee2ZZ class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEee2ZZ>
  : public ClassTraitsBase<Herwig::MEee2ZZ> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEee2ZZ"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEee2ZZ is implemented. It may also include several, space-separated,
   * libraries if the class MEee2ZZ depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMELepton.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEee2ZZ_H */
