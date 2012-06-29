// -*- C++ -*-
//
// MEee2gZ2qq.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEee2gZ2qq_H
#define HERWIG_MEee2gZ2qq_H
//
// This is the declaration of the MEee2gZ2qq class.
//

#include "Herwig++/MatrixElement/HwME2to2Base.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "Herwig++/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEee2gZ2qq class implements the matrix element
 * for \f$e^+e^-\to Z/\gamma \to q\bar{q}\f$ including spin correlations.
 * The class includes greater control over the type of quark produced than is available
 * in the corresponding matrix element from ThePEG, in addition to spin correlations.
 *
 * @see \ref MEee2gZ2qqInterfaces "The interfaces"
 * defined for MEee2gZ2qq.
 */
class MEee2gZ2qq: public HwME2to2Base {

public:

  /**
   * The default constructor.
   */
  inline MEee2gZ2qq() : _minflav(1), _maxflav(5), _massopt(1)  
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

private:

  /**
   * Member to calculate the matrix element
   * @param fin  Spinors for incoming fermion
   * @param ain  Spinors for incoming antifermion
   * @param fout Spinors for outgoing fermion
   * @param aout Spinors for outgong antifermion
   * @param me   Spin summed Matrix element
   * @param cont The continuum piece of the matrix element
   * @param BW   The Z piece of the matrix element
   */
  ProductionMatrixElement HelicityME(vector<SpinorWaveFunction>    & fin,
				     vector<SpinorBarWaveFunction> & ain,
				     vector<SpinorBarWaveFunction> & fout,
				     vector<SpinorWaveFunction>    & aout,
				     double & me,
				     double & cont,
				     double & BW ) const;
  
private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEee2gZ2qq> initMEee2gZ2qq;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2gZ2qq & operator=(const MEee2gZ2qq &);

private:

  /**
   *  Pointer to the fermion-antifermion Z vertex
   */
  AbstractFFVVertexPtr _theFFZVertex;
  
  /**
   *  Pointer to the fermion-antifermion photon vertex
   */
  AbstractFFVVertexPtr _theFFPVertex;
  
  /**
   *  Pointer to the particle data object for the Z
   */
  PDPtr _Z0;

  /**
   *  Pointer to the particle data object for the photon
   */
  PDPtr _gamma;

  /**
   *  The minimum PDG of the quarks to be produced
   */
   unsigned int _minflav;

  /**
   *  The maximum PDG of the quarks to be produced
   */
   unsigned int _maxflav;

  /**
   *  Option for the treatment of the top quark mass
   */
  unsigned int _massopt;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEee2gZ2qq. */
template <>
struct BaseClassTrait<Herwig::MEee2gZ2qq,1> {
  /** Typedef of the first base class of MEee2gZ2qq. */
  typedef Herwig::HwME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEee2gZ2qq class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEee2gZ2qq>
  : public ClassTraitsBase<Herwig::MEee2gZ2qq> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEee2gZ2qq"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MEee2gZ2qq class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwMELepton.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEee2gZ2qq_H */
