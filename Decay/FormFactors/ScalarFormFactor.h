// -*- C++ -*-
//
// ScalarFormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ScalarFormFactor_H
#define HERWIG_ScalarFormFactor_H
//
// This is the declaration of the ScalarFormFactor class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "ScalarFormFactor.fh"
#include "ThePEG/Config/Complex.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 * 
 *  The ScalarFormFactor class is the base class for the form factors
 *  for the weak decay of scalar to scalar, vector and tensor mesons.
 *
 *  This is designed so that the form factors can be used for both the semi-leptonic
 *  decays and using factorization for hadronic weak decays. It also includes the
 *  necessary form factors for the weak-radiative decays \f$B\to K\gamma\f$
 *   and \f$B\to K \ell^+\ell^-\f$.
 *
 *  The form factors are given below for the decay \f$X(p_X)\to Y(p_Y)\f$ with 
 *  \f$q_\mu=(p_X-p_Y)_\mu\f$.
 *
 *  The scalar-scalar form factor is defined by
 *
 *  \f[ \langle Y(p_Y)|(V-A)_\mu|X(p_X)\rangle = 
 *      \left\{(p_X+p_Y)_\mu-\frac{m_X^2-m_Y^2}{q^2}q_\mu\right\}f_+(q^2)
 *     +\left\{\frac{m_X^2-m_Y^2}{q^2}q_\mu\right\}f_0(q^2).  \f]
 *
 *  This is the form factor for the standard weak current. For the weak radiative
 *  decays we also need the penguin current.
 *  
 *  \f[ \langle Y(p_Y)|J^\sigma_\mu|X(p_X)\rangle = \frac{i}{m_X+m_Y}
 *   \left\{q^2(p_X+p_Y)_\mu-(m_X^2-m_Y^2)q_\mu\right\}f_T(q^2)\f]
 *
 *
 *  The scalar-vector form factors are defined so that
 *
 * \f[ \langle Y(p_Y)|(V-A)_\mu|X(p_X)\rangle = -i\epsilon^*_\mu(m_X+m_Y)A_1(q^2)
 * +i(p_X+p_Y)_\mu\epsilon^*\cdot q \frac{A_2(q^2)}{m_X+m_Y}\f]
 *\f[\phantom{\langle Y(p_Y)|(V-A)_\mu|X(p_X)\rangle =}
 * +iq_\mu\epsilon^*\cdot q \frac{2m_Y}{q^2}\left(A_3(q^2)-A_0(q^2)\right)
 * +\epsilon_{\mu\nu\rho\sigma}\epsilon^{*\nu}p_X^\rho p_Y^\sigma \frac{2V(q^2)}{m_X+m_Y},
 * \f]
 *
 *  where the form factor \f$A_3(q^2)\f$ can be defined in terms of \f$A_1\f$
 *  and \f$A_2\f$ using
 *
 * \f[ A_3(q^2) = \frac{m_X+m_Y}{2m_Y}A_1(q^2)-\frac{m_X-m_Y}{2m_Y}A_2(q^2)\f]
 *
 *  and \f$A_0(0)=A_3(0)\f$.
 *
 *  As with the scalar to scalar currents this is the form factor for the standard
 *  weak current. We also need the penguin current
 *
 * \f[ \langle Y(p_Y)|J^\sigma_\mu|X(p_X)\rangle = 
 * i\epsilon_{\mu\nu\rho\sigma}\epsilon^{*\nu}p_X^\rho p_Y^\sigma2T_1(q^2)
 * +T_2(q^2)\left\{\epsilon^*_\mu(m_X^2-m^2_Y)-\epsilon^*\cdot q (p_X+p_Y)_\mu\right\}
 * \f]\f[\phantom{\langle Y(p_Y)|J^\sigma_\mu|X(p_X)\rangle = }
 * +T_3(q^2)\epsilon^*\cdot q\left\{q_\mu-\frac{q^2}{m^2_X-m^2_Y}(p_X+p_Y)_\mu\right\},
 * \f]
 * 
 * with \f$T_1(0)=T_2(0)\f$.
 *
 *  The scalar-tensor form factors are defined as
 *
 * \f[ \langle Y(p_Y)|(V-A)_\mu|X(p_x)\rangle = 
 *  i h(q^2) \epsilon_{\mu\nu\lambda\rho} \epsilon^{*\nu\alpha} p_{Y\alpha}
 *    (p_X+p_Y)^\lambda(p_X-p_Y)^\rho
 *    -k(q^2)\epsilon^*_{\mu\nu}p_Y^\nu 
 * \f]\f[\phantom{\langle Y(p_Y)|(V-A)_mu|X(p_x)\rangle =}
 *    -b_+(q^2)\epsilon^*_{\alpha\beta}p_X^\alpha p_X^\beta(p_X+p_Y)_\mu
 *    -b_-(q^2)\epsilon^*_{\alpha\beta}p_X^\alpha p_X^\beta(p_X-p_Y)_\mu.
 * \f]
 *
 *  This is the base class and contains virtual methods which should return
 *  the form factors described above. This class stores information on the
 *  incoming and outgoing particles for a given form factor, the spin of the
 *  outgoing particle and the id's of the quarks.
 *
 *  Classes inheriting from this class should specify which combinations of
 *  particles etc are allowed using the addFormFactor member.
 *
 * @see BaryonFormFactor
 */

class ScalarFormFactor: public Interfaced {

public:

  /**
   * Default constructor
   */
  ScalarFormFactor() : _numbermodes(0) {}

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
   * \param in   PDG code for the incoming meson.
   * \param out  PDG code for the outgoing meson.
   * \param cc  particles or charge conjugates stored in form factor.
   * @return The location in the vectors storing the data. 
   */
  int formFactorNumber(int in,int out,bool & cc) const  {
    if(_incomingid.size()==0) return -1;
    int output(-1);unsigned int ix(0);
    do {
      if(_incomingid[ix]== in && _outgoingid[ix]== out) {
	cc=false;
	output=ix;
      }
      else if (_incomingid[ix]==-in && _outgoingid[ix]==-out) {
	cc=true;
	output=ix;
      }
      else if(_incomingid[ix]==-in && _outgoingid[ix]==out &&
	      (abs(_outgoingid[ix])/100)%10==(abs(_outgoingid[ix])/10)%10) {
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
   * @param id0 The PDG code for the incoming meson.
   * @param id1 The PDG code for the outgoing meson.
   */
  void particleID(unsigned int iloc,int& id0,int& id1) const {
    id0=_incomingid[iloc];
    id1=_outgoingid[iloc];
  }

  /**
   * Information on the form factor.
   * @param iloc The location in the list.
   * @param ispin The spin of the outgoing meson.
   * @param spect The PDG code of the spectator quark.
   * @param inquark The PDG code for decaying incoming quark.
   * @param outquark The PDG code for the outgoing quark produced in the decay.
   */
  void formFactorInfo(unsigned int & iloc,int & ispin,int & spect,
		      int & inquark, int & outquark) const {
    ispin    = _outgoingJ[iloc];
    spect    = _spectator[iloc];
    inquark  = _inquark[iloc];
    outquark = _outquark[iloc];
  }

  /**
   * Information on the form factor.
   * @param in The PDG code of the incoming meson.
   * @param out The PDG code of the outgoing meson.
   * @param ispin The spin of the outgoing meson.
   * @param spect The PDG code of the spectator quark.
   * @param inquark The PDG code for decaying incoming quark.
   * @param outquark The PDG code for the outgoing quark produced in the decay.
   */
  void formFactorInfo(int in,int out,int & ispin,
		      int & spect,int & inquark, int & outquark) const {
    bool dummy;
    unsigned int ix(formFactorNumber(in,out,dummy));
    formFactorInfo(ix,ispin,spect,inquark,outquark);
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
   * The form factor for the weak decay of a scalar to a scalar.  
   * This method is virtual and must be implementented in classes
   * inheriting from this which include scalar to scalar form factors.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param f0 The form factor \f$f_0\f$. 
   * @param fp The form factor \f$f_+\f$.
   */
  virtual void ScalarScalarFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
				      Energy m0,Energy m1,Complex & f0,
				      Complex & fp) const;

  /**
   * The form factor for the weak decay of a scalar to a vector. This method is virtual
   * and must be implemented in classes inheriting from this which include scalar to
   * vector form factors.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param A0 The form factor \f$A_0\f$
   * @param A1 The form factor \f$A_1\f$
   * @param A2 The form factor \f$A_2\f$
   * @param V  The form factor \f$V\f$
   */
  virtual void ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0, int id1,
				      Energy m0, Energy m1,Complex & A0,
				      Complex & A1,Complex & A2, Complex & V) const;

  /**
   * The form factor for the weak decay of a scalar to a tensor. This method is virtual
   * and must be implemented in classes inheriting from this which include scalar to
   * tensor form factors.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param h  The form factor \f$h\f$.
   * @param k  The form factor \f$k\f$. 
   * @param bp The form factor \f$b_+\f$.
   * @param bm The form factor \f$b_-\f$.
   */
  virtual void ScalarTensorFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
				      Energy m0, Energy m1, complex<InvEnergy2> & h,
				      Complex & k, complex<InvEnergy2> & bp,
				      complex<InvEnergy2> & bm) const;

  /**
   * The form factor for the weak penguin decay of a scalar meson to a scalar meson.
   * This method is virtual
   * and must be implemented in classes inheriting from this which include scalar to
   * scalar penguin form factors. 
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param fT The form factor \f$f_T\f$.
   */
  virtual void ScalarScalarSigmaFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
					   Energy m0, Energy m1,Complex & fT) const;

  /**
   * The form factor for the weak penguin decay of a scalar meson to a vector meson. 
   * This method is virtual
   * and must be implemented in classes inheriting from this which include scalar to
   * vector penguin form factors.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param T1 The form factor \f$T_1\f$.
   * @param T2 The form factor \f$T_2\f$.
   * @param T3 The form factor \f$T_3\f$.
   */
  virtual void ScalarVectorSigmaFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
					   Energy m0, Energy m1, Complex & T1,
					   Complex & T2, Complex & T3) const;

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
   * @param in The PDG code of the incoming meson.
   * @param out The PDG code of the outgoing meson.
   * @param spin The spin of the outgoing meson.
   * @param spect The PDG code of the spectator quark.
   * @param inquark The PDG code for decaying incoming quark.
   * @param outquark The PDG code for the outgoing quark produced in the decay.
   */
  void addFormFactor(int in,int out,int spin, int spect,
		     int inquark, int outquark) {
    _incomingid.push_back(in);
    _outgoingid.push_back(out);
    _outgoingJ.push_back(spin);
    _spectator.push_back(spect);
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
  static AbstractClassDescription<ScalarFormFactor> initScalarFormFactor;

  /**
   * Private and non-existent assignment operator.
   */
  ScalarFormFactor & operator=(const ScalarFormFactor &);

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
   * spin of the outgoing particle
   */
  vector<int> _outgoingJ;

  /**
   * the id of the spectator quark
   */
  vector<int> _spectator;

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
 * ScalarFormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::ScalarFormFactor,1> {
  /** Typedef of the base class of ScalarFormFactor. */
  typedef Interfaced NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * ScalarFormFactor class.
 */
template <>
struct ClassTraits<Herwig::ScalarFormFactor>
  : public ClassTraitsBase<Herwig::ScalarFormFactor> {
  /** Return the class name. */
  static string className() { return "Herwig::ScalarFormFactor"; }
};

/** @endcond */

}

#endif /* HERWIG_ScalarFormFactor_H */
