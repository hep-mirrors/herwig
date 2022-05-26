// -*- C++ -*-
#ifndef Herwig_ScalarTo3ScalarDalitz_H
#define Herwig_ScalarTo3ScalarDalitz_H
//
// This is the declaration of the ScalarTo3ScalarDalitz class.
//

#include "DalitzBase.h"

namespace Herwig {
using namespace ThePEG;
  
  
/**
 * The ScalarTo3ScalarDalitz class provides a base class for the implementation
 * of weak three-body decays of bottom and charm mesons
 *
 * @see \ref ScalarTo3ScalarDalitzInterfaces "The interfaces"
 * defined for ScalarTo3ScalarDalitz.
 */
class ScalarTo3ScalarDalitz: public DalitzBase {
  
public:

  /**
   * The default constructor.
   */
  ScalarTo3ScalarDalitz() : f0gpi_(0.09), f0gK_(0.02), 
			    useResonanceMass_(false)
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
  virtual Complex amplitude(int ichan) const {
    Complex amp(0.);
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
  Complex resAmp(unsigned int i) const;

  /**
   *  Access to the invariants
   */
  const Energy & mInv(unsigned int i,unsigned int j) const {
    return m2_[i][j];
  }

  /**
   *  Masses of the children
   */
  const Energy & mOut(unsigned int i) const {
    return mOut_[i];
  }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ScalarTo3ScalarDalitz & operator=(const ScalarTo3ScalarDalitz &) = delete;

private:
  
  /**
   *  Parameters for the \f$f_0(980)\f$
   */
  //@{
  /**
   * \f$g_\pi\f$ coupling for the \f$f_0(980)\f$ width
   */
  double f0gpi_;

  /**
   * \f$g_K\f$ coupling for the \f$f_0(980)\f$ width
   */
  double f0gK_;
  //@}

  /**
   *  Choice of the mass to use in the denominator of the expressions
   */
  bool useResonanceMass_;
  //@}
  
private :

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
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;
};

}

#endif /* Herwig_ScalarTo3ScalarDalitz_H */
