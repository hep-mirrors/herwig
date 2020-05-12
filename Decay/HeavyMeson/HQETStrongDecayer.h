// -*- C++ -*-
#ifndef Herwig_HQETStrongDecayer_H
#define Herwig_HQETStrongDecayer_H
//
// This is the declaration of the HQETStrongDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 * The HQETStrongDecayer class implements the strong decays of excited heavy, i.e. bottom and charm mesons, using
 * heavy quark effective theory results
 *
 * @see \ref HQETStrongDecayerInterfaces "The interfaces"
 * defined for HQETStrongDecayer.
 */
class HQETStrongDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  HQETStrongDecayer();

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
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 7.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  bool twoBodyMEcode(const DecayMode & dm, int & mecode, double & coupling) const;

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
  HQETStrongDecayer & operator=(const HQETStrongDecayer &) = delete;

private:

  /**
   *  Coupings for the decays
   */
  //@{
  /**
   *  Pion decay constant
   */
  Energy fPi_;

  /**
   *  Coupling for decays within the \f$(0^-,1^-)\f$ multiplet
   */
  double g_;

  /**
   *  Coupling for decays within the \f$(1^+ ,2^+)\f$ multiplet
   */
  double h_;

  /**
   *  Coupling for decays within the \f$(0^+ ,1^+)\f$ multiplet
   */
  double f_;

  /**
   *  D_1 mixing angle (up and down)
   */
  double psiL_;

  /**
   *  D_1 mixing angle (strange)
   */
  double psiS_;

  /**
   *   \f$\eta-\pi^0 \f$ mixing for isospin violating decays
   */
  double deltaEta_;
  //@}

  /**
   *   A momentum scale characterising the convergence of the
   *   derivative expansion. We expect Lambda_ ~ 1 GeV.
   */
  Energy Lambda_;
  //@}

  /**
   *  Particles etc for thye different modes
   */
  //@{
  /**
   * the PDG codes for the incoming particles
   */
  vector<int> incoming_;

  /**
   * the PDG codes for outgoing heavy meson
   */
  vector<int> outgoingH_;

  /**
   * the PDG codes for outgoing light meson
   */
  vector<int> outgoingL_;

  /**
   *  Type of decay
   */
  vector<int> type_;

  /**
   * the maximum weight for the decay
   */
  vector<double> maxWeight_;
  //@}

  /**
   *  Storage of wavefunctions etx
   */
  //@{
  /**
   *   Storage of the \f$\rho\f$ matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Storage of polarization vectors of the decaying particle
   */
  mutable vector<Helicity::LorentzPolarizationVector> vecIn_;

  /**
   *  Storage of polarization tensors of the decaying particle
   */
  mutable vector<Helicity::LorentzTensor<double> > tensorIn_;

  /**
   *  Storage of polarization vectors of the decay product
   */
  mutable vector<Helicity::LorentzPolarizationVector> vecOut_;
  //@}

};

}

#endif /* Herwig_HQETStrongDecayer_H */
