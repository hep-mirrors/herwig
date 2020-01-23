// -*- C++ -*-
#ifndef Herwig_DtoKPiPiFOCUS_H
#define Herwig_DtoKPiPiFOCUS_H
//
// This is the declaration of the DtoKPiPiFOCUS class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/FormFactors/KMatrix.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DtoKPiPiFOCUS class.
 *
 * @see \ref DtoKPiPiFOCUSInterfaces "The interfaces"
 * defined for DtoKPiPiFOCUS.
 */
class DtoKPiPiFOCUS: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  DtoKPiPiFOCUS();
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const tPDVector & outgoing,
	     const vector<Lorentz5Momentum> & momenta,
	     MEOption meopt) const;

  /**
   *   Construct the SpinInfos for the particles produced in the decay
   */
  virtual void constructSpinInfo(const Particle & part,
				 ParticleVector outgoing) const;

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DtoKPiPiFOCUS & operator=(const DtoKPiPiFOCUS &);

private:

  /**
   *  Parameters for the Blatt-Weisskopf form-factors
   */
  //@{
  /**
   *  Radial size for the \f$D^0\f$
   */
  InvEnergy rD0_;

  /**
   *  Radial size for the light resonances
   */
  InvEnergy rRes_;
  //@}

  /**
   *  The K-matrices for the s-wave component
   */
  //@{
  /**
   *  The \f$I=\frac12\f$ component
   */
  KMatrixPtr KHalf_;
  
  /**
   *  The \f$I=\frac32\f$ component
   */
  KMatrixPtr KThreeHalf_;
  //@}
  
  /**
   * Parameters for the f$pf$-vector
   */
  //@{
  /**
   *  Pole couplings
   */
  Energy g1_,g2_;

  /**
   * \f$\beta\f$
   */
  Energy beta_;
  
  /**
   * \f$\theta\f$
   */
  double theta_;

  /**
   *  \f$\gamma\f$
   */
  vector<double> gamma_;

  /**
   *  \f$c_{1i}\f$
   */
  vector<double> c1_;

  /**
   *  \f$c_{2i}\f$
   */
  vector<double> c2_;

  /**
   *  \f$c_{3i}\f$
   */
  vector<double> c3_;
  //@}

  /**
   *  Parameters for the phase-space integration
   */
  //@{
  /**
   *  Maximum weights for the various modes
   */
  double maxWgt_;

  /**
   *  Weights for the different integration channels
   */
  vector<double> weights_;
  //@}

  /**
   *  Masses and widths of the resonances
   */
  //@{
  /**
   *  Masses
   */
  vector<Energy> mRes_;

  /**
   * Widths
   */
  vector<Energy> wRes_;
  //@}

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

};

}

#endif /* Herwig_DtoKPiPiFOCUS_H */
