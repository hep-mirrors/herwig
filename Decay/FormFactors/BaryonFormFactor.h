// -*- C++ -*-
#ifndef HERWIG_BaryonFormFactor_H
#define HERWIG_BaryonFormFactor_H
//
// This is the declaration of the BaryonFormFactor class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "BaryonFormFactor.fh"
#include "ThePEG/Config/Complex.h"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   * The BaryonFormFactor class is the base class for the implementation
   * of the form-factors for the weak decay of a baryon. 
   *
   *  This is designed so that the form factors can be used for both the semi-leptonic
   *  decays and using factorization for hadronic weak decays.
   *
   *  The form factors are given below for the decay \f$X(p_0)\to Y(p_1)\f$ with 
   *  \f$q_\mu=(p_0-p_1)_\mu\f$.
   *
   *
   *  The form-factors are defined to be 
   *
   *   \f[\bar{u}(p_1) \left[  \gamma^\mu \left(F^V_1+F^A_1 \gamma_5\right)
   *  +\frac{i}{(m_0+m_1)}\sigma_{\mu\nu}q^\nu\left(F^V_2+F^A_2\gamma_5\right)
   *  +\frac1{(m_0+m_1)}q^\mu\left(F^V_3+F^A_3\gamma_5\right)\right] u(p_0) \f]
   *
   *   for the \f$\frac12\to\frac12\f$ transition and
   *
   *  \f[\bar{u}^\alpha(p_1) \left[ g_{\alpha\mu}\left(G^V_1+G^A_1 \gamma_5\right)
   *   +\frac1{(m_0+m_1)}p_{0\alpha}\gamma_\mu\left(G^V_2+G^A_2\gamma_5\right)
   *   +\frac1{(m_0+m_1)^2}p_{0\alpha}p_{1\mu}\left(G^V_3+G^A_3\gamma_5\right)\right.\f]
   * \f[ \left.
   *   +\frac1{(m_0+m_1)^2}p_{0\alpha}q_\mu\left(G^V_4+G^A_4\gamma_5\right)\right]
   * \gamma_5 u(p_0) 
   * \f]
   *   for the \f$\frac12\to\frac32\f$ transition.
   *
   *  These definitions differ from those in the liturature because we have divided some
   *  terms by the sum of the baryon masses to ensure that the form-factors are all
   *  dimensionless.
   *
   *  In many applications, particularly for the decay of baryons containing a heavy
   *  quark, an alternative version of the form factors in terms of the velocities
   *  of the baryons is used. This form is
   *
   *   \f[\bar{u}(v') \left[  \gamma^\mu \left(F_1-G_1 \gamma_5\right)
   *                               v^\mu \left(F_2-G_2 \gamma_5\right)
   *                            {v'}^\mu \left(F_3-G_3\gamma_5\right)\right] u(v) \f]
   *
   *   for the \f$\frac12\to\frac12\f$ transition and
   *
   *  \f[\bar{u}^\alpha(v') \left[ v_\alpha\gamma_\mu\left(N_1-K_1 \gamma_5\right)
   *                              +v_\alpha v^\mu    \left(N_2-K_2 \gamma_5\right)
   *                              +v_\alpha {v'}^\mu \left(N_3-K_3 \gamma_5\right)
   *                              +g_\alpha^\mu      \left(N_4-K_4 \gamma_5\right)\right]
   * \gamma_5 u(v) 
   * \f]
   *   for the \f$\frac12\to\frac32\f$ transition.
   *
   *  In terms of these form factors the form factors we use are
   *
   * \f[F^V_1= F_1+\frac12(m_0+m_1)\left(\frac{F_2}{m_0}+\frac{F_3}{m_1}\right)\f]
   * \f[F^V_2=\frac12(m_0+m_1)\left(\frac{F_2}{m_0}+\frac{F_3}{m_1}\right)\f]
   * \f[F^V_3=\frac12(m_0+m_1)\left(\frac{F_2}{m_0}-\frac{F_3}{m_1}\right)\f]
   * \f[F^A_1=-G_1+\frac12(m_0-m_1)\left(\frac{G_2}{m_0}+\frac{G_3}{m_1}\right)\f]
   * \f[F^A_2=-\frac12(m_0+m_1)\left(\frac{G_2}{m_0}+\frac{G_3}{m_1}\right)\f]
   * \f[F^A_3=\frac12(m_0+m_1)\left(-\frac{G_2}{m_0}+\frac{G_3}{m_1}\right)\f]
   *
   *   for the \f$\frac12\to\frac12\f$ transition and
   *
   * \f[G_1^V = N_4\f]
   * \f[G_2^V = N_1\frac{(m_0+m_1)}{m_0}\f]
   * \f[G_3^V = \frac{(m_0+m_1)^2}{m_0}\left(\frac{N_2}{m_0}+\frac{N_3}{m_1}\right)\f]
   * \f[G_4^V = N_2\frac{(m_0+m_1)^2}{m^2_0}\f]
   * \f[G_1^A =-K_4\f]
   * \f[G_2^A =-K_1\frac{(m_0+m_1)}{m_0}\f]
   * \f[G_3^A =-\frac{(m_0+m_1)^2}{m_0}\left(\frac{K_2}{m_0}+\frac{K_3}{m_1}\right)\f]
   * \f[G_4^A =-K_2\frac{(m_0+m_1)^2}{m^2_0}\f]
   *
   *   for the \f$\frac12\to\frac32\f$ transition.
   *
   * 
   * @see ScalarFormFactor
   * 
   */

class BaryonFormFactor: public Interfaced {

public:

  /**
   * Default constructor
   */
  BaryonFormFactor() : _numbermodes() {}

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

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

public:

  /** @name Functions to give information about the form factors available. */
  //@{

  /**
   * Find the location for a given pair of particle. 
   * \param in   PDG code for the incoming baryon.
   * \param out  PDG code for the outgoing baryon.
   * \param cc  particles or charge conjugates stored in form factor.
   * @return The location in the vectors storing the data. 
   */
  int formFactorNumber(int in,int out,bool & cc) const {
    int output(-1);
    unsigned int ix(0);
    if(_incomingid.size()==0){return output;}
    do {
      if(_incomingid[ix]== in && _outgoingid[ix]== out) {
	cc=false;
	output=ix;
      }
      else if (_incomingid[ix]==-in && _outgoingid[ix]==-out) {
	cc=true;
	output=ix;
      }
      ++ix;
    }
    while(ix<_incomingid.size()&&output<0);
    return output;
  }

  /**
   * Get the particle ids for an entry.
   * @param iloc The location in the list.
   * @param id0 The PDG code for the incoming baryon.
   * @param id1 The PDG code for the outgoing baryon.
   */
  void particleID(int iloc,int& id0,int& id1) {
    id0=_incomingid[iloc];
    id1=_outgoingid[iloc];
  }

  /**
   * Information on the form factor.
   * @param iloc The location in the list.
   * @param ispin The spin of the incoming baryon.
   * @param ospin The spin of the outgoing baryon.
   * @param spect1 The PDG code of the first spectator quark.
   * @param spect2 The PDG code of the second spectator quark.
   * @param inquark The PDG code for decaying incoming quark.
   * @param outquark The PDG code for the outgoing quark produced in the decay.
   */
  void formFactorInfo(int iloc,int & ispin,int & ospin,int & spect1,
		      int & spect2, int & inquark,int & outquark) {
    ispin    = _incomingJ[iloc];
    ospin    = _outgoingJ[iloc];
    spect1   = _spectator1[iloc];
    spect2   = _spectator2[iloc];
    inquark  = _inquark[iloc];
    outquark = _outquark[iloc]; 
  }

  /**
   * Information on the form factor.
   * @param in The PDG code of the incoming baryon.
   * @param out The PDG code of the outgoing baryon.
   * @param ispin The spin of the incoming baryon.
   * @param ospin The spin of the outgoing baryon.
   * @param spect1 The PDG code of the first spectator quark.
   * @param spect2 The PDG code of the second spectator quark.
   * @param inquark The PDG code for decaying incoming quark.
   * @param outquark The PDG code for the outgoing quark produced in the decay.
   */
  void formFactorInfo(int in,int out,int & ispin,int & ospin,int & spect1,
		      int & spect2, int & inquark,int & outquark) {
    bool dummy;
    unsigned int ix=formFactorNumber(in,out,dummy);
    formFactorInfo(ix,ispin,ospin,spect1,spect2,inquark,outquark);
  }

  /**
   * number of form factors
   */
  unsigned int numberOfFactors() const {return _incomingid.size();}
  //@}

public:

  /** @name Form Factors */
  //@{
  /**
   * The form factor for the weak decay of a spin \f$\frac12\f$ baryon to a 
   * spin \f$\frac12\f$ baryon.  
   * This method is virtual and must be implementented in classes
   * inheriting from this which include spin\f$\frac12\f$  to spin \f$\frac12\f$
   * form factors.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming baryon.
   * @param id1 The PDG code of the outgoing baryon.
   * @param m0 The mass of the incoming baryon.
   * @param m1 The mass of the outgoing baryon.
   * @param f1v The form factor \f$F^V_1\f$.
   * @param f2v The form factor \f$F^V_2\f$.
   * @param f3v The form factor \f$F^V_3\f$.
   * @param f1a The form factor \f$F^A_1\f$.
   * @param f2a The form factor \f$F^A_2\f$.
   * @param f3a The form factor \f$F^A_3\f$.
   */
  virtual void SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc, int id0, int id1,
					  Energy m0, Energy m1,
					  Complex & f1v,Complex & f2v,Complex & f3v,
					  Complex & f1a,Complex & f2a,Complex & f3a);

  /**
   * The form factor for the weak decay of a spin \f$\frac12\f$ baryon to a 
   * spin \f$\frac32\f$ baryon.  
   * This method is virtual and must be implementented in classes
   * inheriting from this which include spin\f$\frac12\f$  to spin \f$\frac32\f$
   * form factors.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming baryon.
   * @param id1 The PDG code of the outgoing baryon.
   * @param m0 The mass of the incoming baryon.
   * @param m1 The mass of the outgoing baryon.
   * @param g1v The form factor \f$G^V_1\f$.
   * @param g2v The form factor \f$G^V_2\f$.
   * @param g3v The form factor \f$G^V_3\f$.
   * @param g4v The form factor \f$G^V_4\f$.
   * @param g1a The form factor \f$G^A_1\f$.
   * @param g2a The form factor \f$G^A_2\f$.
   * @param g3a The form factor \f$G^A_3\f$.
   * @param g4a The form factor \f$G^A_4\f$.
   */
  virtual void SpinHalfSpinThreeHalfFormFactor(Energy2 q2,int iloc, int id0, int id1,
					       Energy m0, Energy m1,
					       Complex & g1v,Complex & g2v,Complex & g3v,
					       Complex & g4v,Complex & g1a,Complex & g2a,
					       Complex & g3a,Complex & g4a);
  //@}

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

protected:  

  /**
   * Add a form factor to the list.
   * @param in The PDG code of the incoming baryon.
   * @param out The PDG code of the outgoing baryon.
   * @param inspin The spin of the incoming baryon.
   * @param outspin The spin of the outgoing baryon.
   * @param spect1 The PDG code of the first  spectator quark.
   * @param spect2 The PDG code of the second spectator quark.
   * @param inquark The PDG code for decaying incoming quark.
   * @param outquark The PDG code for the outgoing quark produced in the decay.
   */
  void addFormFactor(int in,int out,int inspin, int outspin, int spect1,
		     int spect2, int inquark,int outquark) {
    _incomingid.push_back(in);
    _outgoingid.push_back(out);
    _incomingJ.push_back(inspin);
    _outgoingJ.push_back(outspin);
    _spectator1.push_back(spect1);
    _spectator2.push_back(spect2);
    _inquark.push_back(inquark);
    _outquark.push_back(outquark);
  }

  /**
   *  Set initial number of modes
   * @param nmodes The number of modes.
   */
  void initialModes(unsigned int nmodes) {_numbermodes=nmodes;}

  /**
   * Get the initial number of modes
   */
  unsigned int initialModes() const {return _numbermodes;}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<BaryonFormFactor> initBaryonFormFactor;

  /**
   * Private and non-existent assignment operator.
   */
  BaryonFormFactor & operator=(const BaryonFormFactor &) = delete;

  private:

  /**
   * the id's of the incoming particles
   */
  vector<int> _incomingid;

  /**
   * the id's of the  outgoing particles
   */
  vector<int> _outgoingid;

  /**
   * spin of the incoming particle
   */
  vector<int> _incomingJ;

  /**
   * spin of the outgoing particle
   */
  vector<int> _outgoingJ;

  /**
   * the id of the first spectator quark
   */
  vector<int> _spectator1;

  /**
   * the id of the second spectator quark
   */
  vector<int> _spectator2;

  /**
   * the id of the decaying quark
   */
  vector<int> _inquark;

  /**
   * the id of the outgoing quark
   */
  vector<int> _outquark;

  /**
   * The initial number of modes
   */
  unsigned int _numbermodes;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * BaryonFormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::BaryonFormFactor,1> {
  /** Typedef of the base class of BaryonFormFactor. */
  typedef Interfaced NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * BaryonFormFactor class.
 */
template <>
 struct ClassTraits<Herwig::BaryonFormFactor>
  : public ClassTraitsBase<Herwig::BaryonFormFactor> {
  /** Return the class name. */
  static string className() { return "Herwig::BaryonFormFactor"; }
};

/** @endcond */

}

#endif /* HERWIG_BaryonFormFactor_H */
