// -*- C++ -*-
//
// MEQCD2to2.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEQCD2to2_H
#define HERWIG_MEQCD2to2_H
//
// This is the declaration of the MEQCD2to2 class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVVertex.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEQCD2to2 class provides the matrix elements for \f$2\to2\f$ 
 * QCD scattering processes in hadron-hadron collisions.
 *
 * @see \ref MEQCD2to2Interfaces "The interfaces"
 * defined for MEQCD2to2.
 */
class MEQCD2to2: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEQCD2to2();

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
   *  Members to calculate the matrix elements
   */
  //@{
  /**
   * Matrix element for \f$gg\to gg\f$.
   * @param g1   The wavefunctions for the first  incoming gluon
   * @param g2   The wavefunctions for the second incoming gluon
   * @param g3   The wavefunctions for the first  outgoing gluon
   * @param g4   The wavefunctions for the second outgoing gluon
   * @param flow The colour flow
   */
  double gg2ggME(vector<VectorWaveFunction> &g1,vector<VectorWaveFunction> &g2,
		 vector<VectorWaveFunction> &g3,vector<VectorWaveFunction> &g4,
		 unsigned int flow) const;

  /**
   * Matrix element for \f$gg\to q\bar{q}\f$
   * @param g1   The wavefunctions for the first  incoming gluon
   * @param g2   The wavefunctions for the second incoming gluon
   * @param q    The wavefunction  for the outgoing quark
   * @param qbar The wavefunction  for the outgoing antiquark
   * @param flow The colour flow
   */
  double gg2qqbarME(vector<VectorWaveFunction> &g1,vector<VectorWaveFunction> &g2,
		    vector<SpinorBarWaveFunction> & q,vector<SpinorWaveFunction> & qbar,
		    unsigned int flow) const;

  /**
   * Matrix element for \f$q\bar{q}\to gg\f$
   * @param q    The wavefunction  for the incoming quark
   * @param qbar The wavefunction  for the incoming antiquark
   * @param g1   The wavefunctions for the first  outgoing gluon
   * @param g2   The wavefunctions for the second outgoing gluon
   * @param flow The colour flow
   */
  double qqbar2ggME(vector<SpinorWaveFunction> & q,vector<SpinorBarWaveFunction> & qbar,
		    vector<VectorWaveFunction> &g1,vector<VectorWaveFunction> &g2,
		    unsigned int flow) const;

  /**
   * Matrix element for \f$qg\to qg\f$
   * @param qin  The wavefunction for the incoming quark
   * @param g2   The wavefunction for the incoming gluon
   * @param qout The wavefunction for the outgoing quark
   * @param g4   The wavefunction for the outgoing gluon
   * @param flow The colour flow
   */
  double qg2qgME(vector<SpinorWaveFunction> & qin,vector<VectorWaveFunction> &g2,
		 vector<SpinorBarWaveFunction> & qout,vector<VectorWaveFunction> &g4,
		 unsigned int flow) const;

  /**
   * Matrix elements for \f$\bar{q}g\to \bar{q}g\f$.
   * @param qin  The wavefunction for the incoming antiquark
   * @param g2   The wavefunction for the incoming gluon
   * @param qout The wavefunction for the outgoing antiquark
   * @param g4   The wavefunction for the outgoing gluon
   * @param flow The colour flow
   */
  double qbarg2qbargME(vector<SpinorBarWaveFunction> & qin,
		       vector<VectorWaveFunction> &g2,
		       vector<SpinorWaveFunction> & qout,vector<VectorWaveFunction> &g4,
		       unsigned int flow) const;

  /**
   * Matrix element for \f$qq\to qq\f$
   * @param q1 The wavefunction  for the first  incoming quark
   * @param q2 The wavefunction  for the second incoming quark
   * @param q3 The wavefunction  for the first  outgoing quark
   * @param q4 The wavefunction  for the second outgoing quark
   * @param flow The colour flow
   */
  double qq2qqME(vector<SpinorWaveFunction> & q1, vector<SpinorWaveFunction> & q2,
		 vector<SpinorBarWaveFunction> & q3, vector<SpinorBarWaveFunction> & q4,
		 unsigned int flow) const;

  /**
   * Matrix element for \f$\bar{q}\bar{q}\to \bar{q}\bar{q}\f$
   * @param q1 The wavefunction  for the first  incoming antiquark
   * @param q2 The wavefunction  for the second incoming antiquark
   * @param q3 The wavefunction  for the first  outgoing antiquark
   * @param q4 The wavefunction  for the second outgoing antiquark
   * @param flow The colour flow
   */
  double qbarqbar2qbarqbarME(vector<SpinorBarWaveFunction> & q1,
			     vector<SpinorBarWaveFunction> & q2,
			     vector<SpinorWaveFunction>    & q3,
			     vector<SpinorWaveFunction>    & q4,
			     unsigned int flow) const;

  /**
   * Matrix element for \f$q\bar{q}\to q\bar{q}\f$
   * @param q1 The wavefunction  for the incoming quark
   * @param q2 The wavefunction  for the incoming antiquark
   * @param q3 The wavefunction  for the outgoing quark
   * @param q4 The wavefunction  for the outgoing antiquark
   * @param flow The colour flow
   */
  double qqbar2qqbarME(vector<SpinorWaveFunction> & q1,
		       vector<SpinorBarWaveFunction> & q2,
		       vector<SpinorBarWaveFunction>    & q3,
		       vector<SpinorWaveFunction>    & q4,
		       unsigned int flow) const;
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEQCD2to2> initMEQCD2to2;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEQCD2to2 & operator=(const MEQCD2to2 &) = delete;

private:

  /**
   *  Vertices needed to compute the diagrams
   */
  //@{
  /**
   *  \f$gggg\f$ vertex
   */
  AbstractVVVVVertexPtr _ggggvertex;

  /**
   *  \f$ggg\f$ vertex
   */
  AbstractVVVVertexPtr _gggvertex;

  /**
   *  \f$q\bar{q}g\f$ vertex
   */
  AbstractFFVVertexPtr _qqgvertex;
  //@}

  /**
   *  Maximum numbere of quark flavours to include
   */
  unsigned int _maxflavour;

  /**
   *  Processes to include
   */
  unsigned int _process;

  /**
   *  Colour flow
   */
  mutable unsigned int _flow;

  /**
   *  Diagram
   */
  mutable unsigned int _diagram;

  /**
   *  Matrix element
   */
  mutable ProductionMatrixElement _me;

  /**
   *  ParticleData objects of the partons
   */  
  //@{
  /**
   *  The gluon
   */
  PDPtr _gluon;
  
  /**
   *  the quarks
   */
  vector<PDPtr> _quark;

  /**
   *  the antiquarks
   */
  vector<PDPtr> _antiquark;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEQCD2to2. */
template <>
struct BaseClassTrait<Herwig::MEQCD2to2,1> {
  /** Typedef of the first base class of MEQCD2to2. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEQCD2to2 class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEQCD2to2>
  : public ClassTraitsBase<Herwig::MEQCD2to2> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEQCD2to2"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEQCD2to2 is implemented. It may also include several, space-separated,
   * libraries if the class MEQCD2to2 depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEQCD2to2_H */
