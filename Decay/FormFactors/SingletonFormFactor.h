// -*- C++ -*-
#ifndef HERWIG_SingletonFormFactor_H
#define HERWIG_SingletonFormFactor_H
//
// This is the declaration of the SingletonFormFactor class.
//

#include "BaryonFormFactor.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/ParticleData.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The SingletonFormFactor class implements the form-factors from
 *  PRD43, 2939 for the decay of spin-1/2 baryons containing bottom and charm
 *  quarks.
 *
 * @see BaryonFormFactor
 * 
 */
class SingletonFormFactor: public BaryonFormFactor {

public:

  /**
   * Default constructor
   */
  SingletonFormFactor();

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

public:

  /** @name Form Factors */
  //@{
  /**
   * The form factor for the weak decay of a spin \f$\frac12\f$ baryon to a 
   * spin \f$\frac12\f$ baryon.  
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
  //@}

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * Private and non-existent assignment operator.
   */
  SingletonFormFactor & operator=(const SingletonFormFactor &) = delete;

private:

  /**
   * The charm quark mass.
   */
  Energy _mcharm;

  /**
   * The strange quark mass.
   */
  Energy _mstrange;

  /**
   *  The mixing angle for the \f$\Lambda\f$.
   */
  double _thetalambda;

  /**
   *  The mixing angle for the \f$\Sigma\f$.
   */
  double _thetasigma;

  /**
   *  The mixing angle for the \f$\Xi\f$.
   */
  double _thetaxi;

  /**
   *  The mixing angle for the \f$\Xi'\f$.
   */
  double _thetaxip;

  /**
   *  The pole masses for the \f$q^2\f$ dependence of the form factors.
   */
  vector<Energy> _polemass;

  /**
   *  The \f$\xi\f$ parameter for the form factors.
   */
  vector<double> _xi;

  /**
   *  The normalisation factor, \f$N_{mM}\f$, for the form factors.
   */
  vector<double> _nmM;

  /**
   *  The mass of the quark for the form factor.
   */
  vector<Energy> _mquark;

};

}

#endif /* HERWIG_SingletonFormFactor_H */
