// -*- C++ -*-
#ifndef HERWIG_NonLeptonicHyperonDecayer_H
#define HERWIG_NonLeptonicHyperonDecayer_H
//
// This is the declaration of the NonLeptonicHyperonDecayer class.
//
#include "Baryon1MesonDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  This is a general class for the non-leptonic decay of hyperons. The
 *  decays are given in terms of the invariant amplitudes
 *  \f[\bar{u}_{B_j} \left\{A+B\gamma_5\right\}u_{B_i}\f]
 *  where \f$B_j\f$ is the outgoing baryon and \f$B_i\f$ is the incoming baryon.
 *
 *  The default amplitudes are taken from the fit in hep-ph/9902351, 
 *  N.B. due to the sign conventions in hep-ph/9902351 the B amplitudes
 *  have the opposite sign.
 *
 * @see Baryon1MesonDecayerBase
 * 
 */
class NonLeptonicHyperonDecayer: public Baryon1MesonDecayerBase {

public:

  /**
   * Default constructor.
   */
  NonLeptonicHyperonDecayer();

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;


  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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

  /**
   *  Coupling Members.
   */
  //@{
  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac12\f$ and a scalar.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A The coupling \f$A\f$ described above.
   * @param B The coupling \f$B\f$ described above.
   */
  virtual void halfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
				      Complex& A,Complex& B) const;

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

  /**
   * Initialize this object to the begining of the run phase.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   * Private and non-existent assignment operator.
   */
  NonLeptonicHyperonDecayer & operator=(const NonLeptonicHyperonDecayer &) = delete;

private:

  /**
   * PDG code for the incoming baryon.
   */
  vector<long> _incomingB;

  /**
   * PDG code for the outgoing baryon.
   */
  vector<long> _outgoingB;

  /**
   * PDG code for the outgoing meson
   */
  vector<long> _outgoingM;

  /**
   * The \f$A\f$ coefficient.
   */
  vector<double> _a;

  /**
   * The \f$B\f$ coefficient.
   */
  vector<double> _b;

  /**
   * the maximum weights for the decays
   */
  vector<double> _maxweight;

  /**
   *  initial size fo the vectors
   */
  unsigned int _initsize;
};

}


#endif /* HERWIG_NonLeptonicHyperonDecayer_H */
