// -*- C++ -*-
//
// SSHPPVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSHPPVertex_H
#define HERWIG_SSHPPVertex_H
//
// This is the declaration of the SSHPPVertex class.
//

#include "Herwig/Models/General/VVSLoopVertex.h"
#include "MSSM.h"

namespace Herwig {

  /**
   * This class implements the effective vertex coupling a higgs
   * to a pair of photons in the MSSM. The loops include third generation
   * sfermions and fermions, \f$W^\pm\f$ and \f$H^\pm\f$
   *
   * @see \ref SSHPPVertexInterfaces "The interfaces"
   * defined for SSHPPVertex.

   */
class SSHPPVertex: public VVSLoopVertex {
  
public:
  
  /**
   * The default constructor.
   */
  SSHPPVertex();
  
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
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param particle1 ParticleData pointer to first particle
   *@param particle2 ParticleData pointer to second particle
   *@param particle3 ParticleData pointer to third particle
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
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
  
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  
private:
  
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SSHPPVertex> initSSHPPVertex;
  
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSHPPVertex & operator=(const SSHPPVertex &);
  
private:

  /**
   * Switch to turn off squark trilinear couplings via A terms for testing
   */
  bool theIncludeTriLinear;

  /**
   *  Treat the pseudoscalar as scalar for comparison with ISAJET
   */
  bool thePseudoScalarTreatment;

  /**
   * A pointer to the MSSM object
   */
  tMSSMPtr theMSSM;

  /**
   * Value of \f$\sin\theta_W\f$
   */
  double theSw;
  
  /**
   * The mass of the \f$W\f$ boson.
   */
  Energy theMw;

  /**
   * The factor \f$\frac{M_Z}{\cos\theta_W}\f$ 
   */
  Energy theZfact;  
  
  /**
   * The mixing matrix factor \f$Q^{2i}_{11}Q^{2i}_{11}\f$ 
   * for the \f$\tilde{t}\f$
   */
  Complex theQt1L;
  
  /**
   * The mixing matrix factor \f$Q^{2i}_{12}Q^{2i}_{12}\f$ 
   * for the \f$\tilde{t}\f$
   */
  Complex theQt1R;
  
  /**
   * The mixing matrix factor \f$Q^{2i}_{12}Q^{2i}_{12}\f$ 
   * for the \f$\tilde{t}\f$
   */
  Complex theQt1LR;

  /**
   * The mixing matrix factor \f$Q^{2i}_{21}Q^{2i}_{21}\f$ 
   * for the \f$\tilde{t}\f$
   */
  Complex theQt2L;

  /**
   * The mixing matrix factor \f$Q^{2i}_{22}Q^{2i}_{22}\f$ 
   * for the \f$\tilde{t}\f$
   */
  Complex theQt2R;

  /**
   * The mixing matrix factor \f$Q^{2i}_{22}Q^{2i}_{22}\f$ 
   * for the \f$\tilde{t}\f$
   */
  Complex theQt2LR;

 /**
   * The mixing matrix factor \f$Q^{2i-1}_{11}Q^{2i-1}_{11}\f$ 
   * for the \f$\tilde{b}\f$
   */
  Complex theQb1L;

 /**
   * The mixing matrix factor \f$Q^{2i-1}_{12}Q^{2i-1}_{12}\f$ 
   * for the \f$\tilde{b}\f$
   */
  Complex theQb1R;

 /**
   * The mixing matrix factor \f$Q^{2i-1}_{12}Q^{2i-1}_{12}\f$ 
   * for the \f$\tilde{b}\f$
   */
  Complex theQb1LR;

  /**
   * The mixing matrix factor \f$Q^{2i-1}_{21}Q^{2i-1}_{21}\f$ 
   * for the \f$\tilde{b}\f$
   */
  Complex theQb2L;

  /**
   * The mixing matrix factor \f$Q^{2i-1}_{22}Q^{2i-1}_{22}\f$ 
   * for the \f$\tilde{b}\f$
   */
  Complex theQb2R; 

  /**
   * The mixing matrix factor \f$Q^{2i-1}_{22}Q^{2i-1}_{22}\f$ 
   * for the \f$\tilde{b}\f$
   */
  Complex theQb2LR; 
  
  /**
   * The mixing matrix factor \f$L^{2i}_{11}L^{2i}_{11}\f$ 
   * for the \f$\tilde{\tau}\f$
   */
  Complex theLt1L;
  
  /**
   * The mixing matrix factor \f$L^{2i}_{12}L^{2i}_{12}\f$ 
   * for the \f$\tilde{\tau}\f$
   */
  Complex theLt1R;
  
  /**
   * The mixing matrix factor \f$L^{2i}_{12}L^{2i}_{12}\f$ 
   * for the \f$\tilde{\tau}\f$
   */
  Complex theLt1LR;

  /**
   * The mixing matrix factor \f$L^{2i}_{21}L^{2i}_{21}\f$ 
   * for the \f$\tilde{\tau}\f$
   */
  Complex theLt2L;

  /**
   * The mixing matrix factor \f$L^{2i}_{22}L^{2i}_{22}\f$ 
   * for the \f$\tilde{\tau}\f$
   */
  Complex theLt2R;

  /**
   * The mixing matrix factor \f$L^{2i}_{22}L^{2i}_{22}\f$ 
   * for the \f$\tilde{\tau}\f$
   */
  Complex theLt2LR;
  
  /**
   * A pointer to the top quark ParticleData object 
   */
  tPDPtr thetop;

  /**
   * A pointer to the bottom quark ParticleData object 
   */
  tPDPtr thebot;

  /**
   * A pointer to the tau lepton ParticleData object 
   */
  tPDPtr thetau;

  /**
   * The sfermion masses 
   */
  vector<Energy> theSfmass;
  
  /**
   * The value of \f$ \tan\beta \f$ 
   */
  double theTanB;

  /**
   * The value of \f$ \sin\alpha \f$ 
   */
  double theSinA;

  /**
   * The value of \f$ \cos\alpha \f$ 
   */
  double theCosA;

  /**
   * The value of \f$ \sin\beta \f$ 
   */
  double theSinB;

  /**
   * The value of \f$ \cos\beta \f$ 
   */
  double theCosB;

  /**
   * The value of \f$ \sin(\alpha + \beta) \f$ 
   */
  double theSinApB;

  /**
   * The value of \f$ \cos(\alpha + \beta) \f$ 
   */
  double theCosApB;

  /**
   * The value of \f$ \sin(\beta-\alpha) \f$ 
   */
  double theSinBmA;

  /**
   * The value of \f$ \cos(\beta-\alpha) \f$ 
   */
  double theCosBmA;

  /**
   * The U mixing matrix
   */
  tMixingMatrixPtr theU;

  /**
   * The V mixing matrix
   */
  tMixingMatrixPtr theV;  

  /**
   * Last value of the coupling calculated
   */
  Complex theCouplast;
  
  /**
   * The scale \f$q^2\f$ at which coupling was last evaluated
   */
  Energy2 theq2last;
  
  /**
   * Whether we have calculated the tensor coefficents yet 
   */
  bool theHaveCoeff;

  /**
   *  ID of the higgs
   */
  long theLastID;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSHPPVertex. */
template <>
struct BaseClassTrait<Herwig::SSHPPVertex,1> {
  /** Typedef of the first base class of SSHPPVertex. */
  typedef Herwig::VVSLoopVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSHPPVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSHPPVertex>
  : public ClassTraitsBase<Herwig::SSHPPVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SSHPPVertex"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SSHPPVertex class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SSHPPVertex_H */
