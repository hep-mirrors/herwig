// -*- C++ -*-
#ifndef Herwig_GtoG3S1OctetSplitFn_H
#define Herwig_GtoG3S1OctetSplitFn_H
//
// This is the declaration of the GtoG3S1OctetSplitFn class.
//

#include "Herwig/Shower/QTilde/SplittingFunctions/SudakovFormFactor.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "Herwig/MatrixElement/Onium/OniumParameters.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The GtoG3S1OctetSplitFn class implements the \f$1\to1\f$ octet splitting \f$g\to J\Psi,\Upsilon\f$
 *
 * @see \ref GtoG3S1OctetSplitFnInterfaces "The interfaces"
 * defined for GtoG3S1OctetSplitFn.
 */
class GtoG3S1OctetSplitFn: public SudakovFormFactor {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  GtoG3S1OctetSplitFn() : O8_(0.015*GeV2*GeV), state_(ccbar), n_(1),
			  m_(1.2*GeV)
  {}
  //@}

public:
  
  /**
   *  Determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const {
    if(ids.size()!=2) return false;
    if(ids[0]->id()!=ParticleID::g) return false;
    int iq=4+state_;
    long idtest = iq*110+3 + (n_-1)*100000;
    if(ids[1]->id() != idtest) return false;
    return true;
  }

  /**
   *  Method to return the evolution scale given the
   *  transverse momentum, \f$p_T\f$ and \f$z\f$.
   */
  virtual Energy calculateScale(double , Energy , IdList ids,unsigned int ) {
    return 2.*ids[1]->mass();
  }


  /**
   *  Members to generate the scale of the next branching
   */
  //@{
  /**
   * Return the scale of the next time-like branching. If there is no 
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param enhance The radiation enhancement factor
   * defined.
   */
  virtual ShoKinPtr generateNextTimeBranching(const Energy startingScale,
					      const IdList &ids,
					      const RhoDMatrix & rho,
					      double enhance, double detuning);

  /**
   * Return the scale of the next space-like decay branching. If there is no 
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param stoppingScale stopping scale for the evolution
   * @param minmass The minimum mass allowed for the spake-like particle.
   * @param ids The PDG codes of the particles in the splitting
   * defined.
   * @param enhance The radiation enhancement factor
   */
  virtual ShoKinPtr generateNextDecayBranching(const Energy ,
						 const Energy ,
						 const Energy ,
						 const IdList &,
						 const RhoDMatrix & ,
						 double ,
						 double )
  {assert(false);}

  /**
   * Return the scale of the next space-like branching. If there is no 
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param x The fraction of the beam momentum
   * defined.
   * @param beam The beam particle
   * @param enhance The radiation enhancement factor
   */
   virtual ShoKinPtr generateNextSpaceBranching(const Energy ,
						 const IdList &,double ,
						 const RhoDMatrix & ,
						 double ,
						 tcBeamPtr ,
						 double )
  {assert(false);}
  //@}
  
public:
  
  /**
   *  Azimuthal angle generation
   */
  //@{
  /**
   * Generate the azimuthal angle of the branching for forward evolution
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiForward(ShowerParticle & ,const IdList &,
				    ShoKinPtr , const RhoDMatrix & )
  {assert(false);}

  /**
   *  Generate the azimuthal angle of the branching for backward evolution
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiBackward(ShowerParticle & ,const IdList & ,
				     ShoKinPtr, const RhoDMatrix & )
  {assert(false);}

  /**
   *  Generate the azimuthal angle of the branching for ISR in decays
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiDecay(ShowerParticle & ,const IdList &,
				  ShoKinPtr , const RhoDMatrix & )
  {assert(false);}
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

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GtoG3S1OctetSplitFn & operator=(const GtoG3S1OctetSplitFn &) = delete;

private:
  
  /**
   *  Access to the parameters for the quarkonium states
   */
  OniumParametersPtr params_;

  /**
   *  The \f$O_1\f$ colour-singlet coefficient
   */
  Energy3 O8_;

  /**
   *  Type of state
   */
  OniumState state_;

  /**
   *  Principal quantum number
   */
  unsigned int n_;

  /**
   *  The quark mass
   */
  Energy m_;

  /**
   *  Fixed value of \f$\alpha_S\f$
   */
  double fixedAlphaS_;
};

}

#endif /* Herwig_GtoG3S1OctetSplitFn_H */
