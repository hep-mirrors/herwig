// -*- C++ -*-
#ifndef HERWIG_ChengHeavyBaryonFormFactor_H
#define HERWIG_ChengHeavyBaryonFormFactor_H
//
// This is the declaration of the ChengHeavyBaryonFormFactor class.
//

#include "BaryonFormFactor.h"
#include "ThePEG/PDT/ParticleData.h"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   *
   *  The ChengHeavyBaryonFormFactor class implements the form-factors
   *  from PRD53, 1457 and PRD56, 2799 for the semi-leptonic decay of bottom
   *  and charm hadrons. It is also intended for use with the results in PRD56, 2799
   *  for non-leptonic decays using the factorization approximation.
   *
   * @see BaryonFormFactor
   */

class ChengHeavyBaryonFormFactor: public BaryonFormFactor {

public:

  /**
   * Default constructor
   */
  ChengHeavyBaryonFormFactor();

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

  /**
   * The form factor for the weak decay of a spin \f$\frac12\f$ baryon to a 
   * spin \f$\frac32\f$ baryon.
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
  ChengHeavyBaryonFormFactor & operator=(const ChengHeavyBaryonFormFactor &) = delete;

private:

  /** @name The quark masses. */
  //@{
  /**
   *  The mass of the up quark.
   */
  Energy _mu;

  /**
   *  The mass of the down quark.
   */
  Energy _md;

  /**
   *  The mass of the strange quark.
   */
  Energy _ms;

  /**
   *  The mass of the charm quark.
   */
  Energy _mc;

  /**
   *  The mass of the bottom quark.
   */
  Energy _mb;
  //@}

  /**
   *  The normalization parameter \f$N_{fi}\f$ is the flavour and spin overlap
   *  from PRD56, 2799.
   */
  vector<double> _Nfi;

  /**
   * The ratio of the flavour and spin overlap for the outgoing over incoming
   * baryon, \f$\eta\f$ from PRD56, 2799.
   */
  vector<double> _eta;

  /**
   *  the form factor \f$f_1\f$ evaluated at the maximum value of \f$q^2\f$.
   */
  vector<double> _f1;

  /**
   *  the form factor \f$f_2\f$ evaluated at the maximum value of \f$q^2\f$.
   */
  vector<double> _f2;

  /**
   *  the form factor \f$f_3\f$ evaluated at the maximum value of \f$q^2\f$.
   */
  vector<double> _f3;

  /**
   *  the form factor \f$g_1\f$ evaluated at the maximum value of \f$q^2\f$.
   */
  vector<double> _g1;

  /**
   *  the form factor \f$g_2\f$ evaluated at the maximum value of \f$q^2\f$.
   */
  vector<double> _g2;

  /**
   *  the form factor \f$g_3\f$ evaluated at the maximum value of \f$q^2\f$.
   */
  vector<double> _g3;

  /**
   *  The pole mass for the vector form factors for the \f$b\to c\f$ transitiion.
   */
  Energy _mVbc;

  /**
   *  The pole mass for the vector form factors for the \f$b\to s\f$ transitiion.
   */
  Energy _mVbs;

  /**
   *  The pole mass for the vector form factors for the \f$c\to s\f$ transitiion.
   */
  Energy _mVcs;

  /**
   *  The pole mass for the vector form factors for the \f$b\to d\f$ transitiion.
   */
  Energy _mVbd;

  /**
   *  The pole mass for the vector form factors for the \f$c\to u\f$ transitiion.
   */
  Energy _mVcu;

  /**
   *  The pole mass for the axial-vector form factors for the \f$b\to c\f$ transitiion.
   */
  Energy _mAbc;

  /**
   *  The pole mass for the axial-vector form factors for the \f$b\to s\f$ transitiion.
   */
  Energy _mAbs;

  /**
   *  The pole mass for the axial-vector form factors for the \f$c\to s\f$ transitiion.
   */
  Energy _mAcs;

  /**
   *  The pole mass for the axial-vector form factors for the \f$b\to d\f$ transitiion.
   */
  Energy _mAbd;

  /**
   *  The pole mass for the axial-vector form factors for the \f$c\to u\f$ transitiion.
   */
  Energy _mAcu;

};

}

#endif /* HERWIG_ChengHeavyBaryonFormFactor_H */
