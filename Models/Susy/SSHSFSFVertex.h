// -*- C++ -*-
//
// SSHSFSFVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSHSFSFVertex_H
#define HERWIG_SSHSFSFVertex_H
//
// This is the declaration of the SSHSFSFVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/SSSVertex.h"
#include "MSSM.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This is the implementation of the coupling for a Higgs in the MSSM
 * to a pair of sfermions.
 *
 * @see \ref SSHSFSFVertexInterfaces "The interfaces"
 * defined for SSHSFSFVertex.
 */
class SSHSFSFVertex: public SSSVertex {

  /** A vector of MixingMatrix pointers. */
  typedef vector<MixingMatrixPtr> MMPVector;

public:

  /**
   * The default constructor.
   */
  SSHSFSFVertex();

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

  
  /**
   * Calculate the coupling at the given scale.
   * @param q2 The scale at which to evaluate the coupling
   * @param particle1 The first particle at the vertex
   * @param particle2 The second particle at the vertex
   * @param particle3 The third particle at the vertex
   */
  void setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2, 
		   tcPDPtr particle3);

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
  static ClassDescription<SSHSFSFVertex> initSSHSFSFVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSHSFSFVertex & operator=(const SSHSFSFVertex &);

private:

  /** @name Functions to calculate the coupling based on the sfermion type. */
  //@{
  /**
   * Calculate the coupling for the first higgs 
   * @param q2 scale
   * @param higgs The ID of the higgs
   * @param smID The ID of the SM particle to which it is a partner.
   * @param alpha The mass eigenstate of an sfermion
   * @param beta  The mass eigenstate of the other sfermion
   */
  void downSF(Energy2 q2, long higgs, long smID, unsigned int alpha, unsigned int beta);
  
  /**
   * Calculate the coupling for the second higgs 
   * @param q2 scale
   * @param higgs The ID of the higgs
   * @param smID The ID of the SM particle to which it is a partner.
   * @param alpha The mass eigenstate of an sfermion
   * @param beta  The mass eigenstate of the other sfermion
   */
  void upSF(Energy2 q2, long higgs, long smID, unsigned int alpha, unsigned int beta);
  
  /**
   * Calculate the coupling for the third higgs 
   * @param q2 scale
   * @param higgs The ID of the higgs
   * @param smID The ID of the SM particle to which it is a partner.
   * @param alpha The mass eigenstate of an sfermion
   * @param beta  The mass eigenstate of the other sfermion
   */
  void leptonSF(Energy2 q2, long higgs, long smID, unsigned int alpha, unsigned int beta);
  
  /**
   *  Calculate the coupling for the charged higgs 
   * @param q2 scale
   * @param id1 The ID of the first sfermion
   * @param id2 The ID of the second sfermion
   */
  void chargedHiggs(Energy2 q2, long id1, long id2);
  
  //@}

private:
  
  /**
   * A vector containing pointers to the mixing matrices, 0 = stop, 
   * 1 = sbottom, 2 = stau 
   */
  MMPVector theMix;

  /**
   * A vector containing the trilinear couplings, quarks then leptons
   */
  vector<complex<Energy> > theTriC;
    
  /**
   * The value of \f$\sin\alpha\f$.
   */
  double theSinA;

  /**
   * The value of \f$\cos\alpha\f$.
   */
  double theCosA;

  /**
   * The value of \f$\sin\beta\f$.
   */
  double theSinB;

  /**
   * The value of \f$\cos\beta\f$.
   */
  double theCosB;

  /**
   * The value of \f$\tan\beta\f$. 
   */
  double theTanB;

  /**
   * The value of \f$\sin(\alpha + \beta)\f$.
   */
  double theSinAB;

  /**
   * The value of \f$\cos(\alpha + \beta)\f$.
   */
  double theCosAB;

  /**
   * The mass of the \f$W\f$. 
   */
  Energy theMw;

  /**
   * The mass of the \f$Z\f$. 
   */
  Energy theMz;
  
  /**
   * The \f$\mu\f$ parameter. 
   */
  Energy theMu;

  /**
   * The value of \f$\sin\theta_W\f$
   */
  double theSw;

  /**
   * The value of \f$\cos\theta_W\f$
   */
  double theCw;

  /**
   * The value of the coupling when it was last evaluated
   */
  complex<Energy> theCoupLast;
  
  /**
   * The scale at which the coupling was last evaluated
   */
  Energy2 theq2Last;
  
  /**
   * The value of g coupling when it was last evaluated
   */
  double thegLast;
  
  /**
   * The ID of the higgs when the vertex was last evaluated 
   */
  long theHLast;
  
  /**
   * The ID of the first sfermion when the vertex was last evaluated 
   */
  long theSF1Last;
  
  /**
   * The ID of the second sfermion when the vertex was last evaluated 
   */
  long theSF2Last;

  /**
   *  The Model
   */
  tMSSMPtr theMSSM;
};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSHSFSFVertex. */
template <>
struct BaseClassTrait<Herwig::SSHSFSFVertex,1> {
  /** Typedef of the first base class of SSHSFSFVertex. */
  typedef ThePEG::Helicity::SSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSHSFSFVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSHSFSFVertex>
  : public ClassTraitsBase<Herwig::SSHSFSFVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SSHSFSFVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSHSFSFVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSHSFSFVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SSHSFSFVertex_H */
