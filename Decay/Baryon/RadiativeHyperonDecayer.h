// -*- C++ -*-
#ifndef HERWIG_RadiativeHyperonDecayer_H
#define HERWIG_RadiativeHyperonDecayer_H
//
// This is the declaration of the RadiativeHyperonDecayer class.
//

#include "Baryon1MesonDecayerBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The RadiativeHyperonDecayer class provides the matrix elements
 * for the decay of hyperons in which a photon is produced using the 
 * model of Phys. Rev. D 59 (1999) 054019 [arXiv:hep-ph/9902431].
 *
 * @see \ref RadiativeHyperonDecayerInterfaces "The interfaces"
 * defined for RadiativeHyperonDecayer.
 */
class RadiativeHyperonDecayer: public Baryon1MesonDecayerBase {

public:

  /**
   * The default constructor.
   */
  RadiativeHyperonDecayer();

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();
protected:

  /**
   *  Coupling Members.
   */
  //@{

  /**
   * Couplings for spin-\f$\frac12\f$ to spin-\f$\frac12\f$ and a vector.
   * This method must be implemented in any class inheriting from this one
   * which includes \f$\frac12\to\frac12+1\f$ decays. 
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A1 The coupling \f$A_1\f$ described above.
   * @param A2 The coupling \f$A_2\f$ described above.
   * @param B1 The coupling \f$B_1\f$ described above.
   * @param B2 The coupling \f$B_2\f$ described above.
   */
  virtual void halfHalfVectorCoupling(int imode, Energy m0, Energy m1, Energy m2,
				      Complex& A1,Complex& A2,
				      Complex& B1,Complex& B2) const;
  //@}

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RadiativeHyperonDecayer & operator=(const RadiativeHyperonDecayer &) = delete;

private:

  /**
   * PDG code for the incoming baryon.
   */
  vector<long> incomingB_;

  /**
   * PDG code for the outgoing baryon.
   */
  vector<long> outgoingB_;

  /**
   * The \f$A\f$ coefficient.
   */
  vector<InvEnergy> A_;

  /**
   * The \f$B\f$ coefficient.
   */
  vector<InvEnergy> B_;

  /**
   * the maximum weights for the decays
   */
  vector<double> maxweight_;

  /**
   *  initial size fo the vectors
   */
  unsigned int initsize_;
};

}

#endif /* HERWIG_RadiativeHyperonDecayer_H */
