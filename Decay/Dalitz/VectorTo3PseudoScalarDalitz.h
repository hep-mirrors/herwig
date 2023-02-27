// -*- C++ -*-
#ifndef Herwig_VectorTo3PseudoScalarDalitz_H
#define Herwig_VectorTo3PseudoScalarDalitz_H
//
// This is the declaration of the VectorTo3PseudoScalarDalitz class.
//

#include "DalitzBase.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The VectorTo3PseudoScalarDalitz class provides a base class for the decay of vector mesons to 3 pseudoscalar mesons
 *
 * @see \ref VectorTo3PseudoScalarDalitzInterfaces "The interfaces"
 * defined for VectorTo3PseudoScalarDalitz.
 */
class VectorTo3PseudoScalarDalitz: public DalitzBase {

public:

  /**
   * The default constructor.
   */
  VectorTo3PseudoScalarDalitz() : coupling_(1./GeV)
  {}

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan,const Particle & part,
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

  /**
   * Method to return an object to calculate the 3 body partial width.
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;

  /**
   * The matrix element to be integrated for the three-body decays as a function
   * of the invariant masses of pairs of the outgoing particles.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element
   */
  virtual double threeBodyMatrixElement(const int imode, const Energy2 q2,
					const  Energy2 s3, const Energy2 s2, const 
					Energy2 s1, const Energy m1,
					const Energy m2, const Energy m3) const;

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
   *  Calculate the amplitude
   */
  virtual complex<InvEnergy2> amplitude(int ichan) const {
    complex<InvEnergy2> amp(ZERO);
    int iloc=-1;
    for(int ix=0;ix<int(resonances().size());++ix) {
      ++iloc;
      if(channel1()>=0) {
	if(ix!=channel1() && ix!=channel2()) continue;
      }
      if(ichan>=0&&ichan!=iloc) continue;
      amp += resAmp(ix);
    }
    return amp;
  }

  /**
   * Calculate the amplitude for the ith resonance
   */
  complex<InvEnergy2> resAmp(unsigned int i) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VectorTo3PseudoScalarDalitz & operator=(const VectorTo3PseudoScalarDalitz &) = delete;

private:

  /**
   *  Coupling for the normalisation of the mode
   */
  InvEnergy coupling_;

  /**
   *   Storage of the kinematics
   */
  //@{
  /**
   *   Mass of the parent
   */
  mutable Energy mD_;

  /**
   *   Masses of the children
   */
  mutable Energy mOut_[3];

  /**
   *   Masses of the children
   */
  mutable Energy m2_[3][3];
  //@}
  
  /**
   *  Storage of polarization vectors to try and increase
   *  speed
   */
  mutable vector<Helicity::LorentzPolarizationVector> vectors_;
  
  /**
   *   Storage of the \f$\rho\f$ matrix
   */
  mutable RhoDMatrix rho_;

};

}

#endif /* Herwig_VectorTo3PseudoScalarDalitz_H */
