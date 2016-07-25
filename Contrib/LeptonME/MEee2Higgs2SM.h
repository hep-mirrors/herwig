// -*- C++ -*-
//
// MEee2Higgs2SM.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEee2Higgs2SM_H
#define HERWIG_MEee2Higgs2SM_H
//
// This is the declaration of the MEee2Higgs2SM class.
//

#include "ThePEG/MatrixElement/ME2to2Base.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEee2Higgs2SM class implements the production of an \f$s\f$-channel
 * Higgs in \f$e^+e^-\f$ collisions in order to allow easy tests of Higgs
 * decays. It should not be used for physics studies.
 *
 * @see \ref MEee2Higgs2SMInterfaces "The interfaces"
 * defined for MEee2Higgs2SM.
 */
class MEee2Higgs2SM: public ME2to2Base {

public:

  /**
   * The default constructor.
   */
  inline MEee2Higgs2SM() : allowed_(0) {}

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
   * set up the spin correlations
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
   *  The matrix element
   * @param fin The incoming spinor wavefunction
   * @param ain The incoming spinorbar wavefunction
   * @param fout The outgoing spinor bar wavefunction
   * @param aout The outgoing spinor wavefunction
   * @param me The spin averaged matrix element
   */
  ProductionMatrixElement HelicityME(vector<SpinorWaveFunction> fin,
				     vector<SpinorBarWaveFunction> ain,
				     vector<SpinorBarWaveFunction> fout,
				     vector<SpinorWaveFunction> aout,double& me) const;
  /**
   *  \f$H\to gg\f$ matrix element
   * @param fin The incoming spinor wavefunction
   * @param ain The incoming spinorbar wavefunction
   * @param g1 Outgoing gluon wavefunction
   * @param g2 Outgoing gluon wavefunction
   * @param me The spin averaged matrix element
   */
  ProductionMatrixElement ggME(vector<SpinorWaveFunction> fin,
			       vector<SpinorBarWaveFunction> ain,
			       vector<VectorWaveFunction> g1,
			       vector<VectorWaveFunction> g2, 
			       double & me) const;
  
private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEee2Higgs2SM> initMEee2Higgs2SM;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2Higgs2SM & operator=(const MEee2Higgs2SM &);

private:

  /**
   *  Pointer to the Higgs fermion-antifermion vertex
   */
  AbstractFFSVertexPtr FFHVertex_;

  /**
   *  Pointer to Higgs-gluon-gluon vertex
   */
  AbstractVVSVertexPtr HGGVertex_;

  /**
   * Allowed outgoing particles
   */
  int allowed_;

  /**
   *  Pointer to the Higgs ParticleData object
   */
  PDPtr h0_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEee2Higgs2SM. */
template <>
struct BaseClassTrait<Herwig::MEee2Higgs2SM,1> {
  /** Typedef of the first base class of MEee2Higgs2SM. */
  typedef ME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEee2Higgs2SM class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEee2Higgs2SM>
  : public ClassTraitsBase<Herwig::MEee2Higgs2SM> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEee2Higgs2SM"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEee2Higgs2SM is implemented. It may also include several, space-separated,
   * libraries if the class MEee2Higgs2SM depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "LeptonME.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEee2Higgs2SM_H */
