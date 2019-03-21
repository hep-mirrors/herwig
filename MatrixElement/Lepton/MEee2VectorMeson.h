// -*- C++ -*-
//
// MEee2VectorMeson.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_MEee2VectorMeson_H
#define THEPEG_MEee2VectorMeson_H
//
// This is the declaration of the MEee2VectorMeson class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig/PDT/GenericMassGenerator.fh"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {

using namespace ThePEG;

using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::VectorWaveFunction;

/**
 * The MEee2VectorMeson class is designed to produce neutral vector mesons
 * in \f$e^+e^-\f$ collisions and is primarily intended to test the hadronic
 * decay package.
 *
 * @see \ref MEee2VectorMesonInterfaces "The interfaces"
 * defined for MEee2VectorMeson.
 */
class MEee2VectorMeson: public MEBase {

public:

  /**
   * The default constructor.
   */
  inline MEee2VectorMeson() :_coupling(0.0012), _lineshape(false) 
  {}

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
   *  Set up the spin correlations 
   */
  virtual void constructVertex(tSubProPtr sub);
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
   *  Member to return the helicity amplitudes
   */
  ProductionMatrixElement HelicityME(vector<SpinorWaveFunction> fin,
				     vector<SpinorBarWaveFunction> ain,
				     vector<VectorWaveFunction> vout,double& me) const;
private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEee2VectorMeson> initMEee2VectorMeson;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2VectorMeson & operator=(const MEee2VectorMeson &) = delete;

private:

  /**
   *  The vector meson being produced
   */
  PDPtr _vector;

  /**
   *  The coupling
   */
  double _coupling;

  /**
   *  Use the mass generator for the line shape
   */
  bool _lineshape;

  /**
   *  Pointer to the mass generator for the Higgs
   */
  GenericMassGeneratorPtr _massgen;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEee2VectorMeson. */
template <>
struct BaseClassTrait<Herwig::MEee2VectorMeson,1> {
  /** Typedef of the first base class of MEee2VectorMeson. */
  typedef MEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEee2VectorMeson class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEee2VectorMeson>
  : public ClassTraitsBase<Herwig::MEee2VectorMeson> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEee2VectorMeson"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEee2VectorMeson is implemented. It may also include several, space-separated,
   * libraries if the class MEee2VectorMeson depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMELepton.so"; }
};

/** @endcond */

}

#endif /* THEPEG_MEee2VectorMeson_H */
