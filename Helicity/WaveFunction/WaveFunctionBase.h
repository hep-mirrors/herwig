// -*- C++ -*-
#ifndef HERWIG_WaveFunctionBase_H
#define HERWIG_WaveFunctionBase_H
//
// This is the declaration of the WaveFunctionBase class.

#include <ThePEG/CLHEPWrap/Lorentz5Vector.h>
#include <ThePEG/CLHEPWrap/LorentzVector.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/Config/Complex.h>

namespace Herwig {
namespace Helicity {
enum Direction {incoming,outgoing,intermediate};
using namespace ThePEG;

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 * This class is the base class for all wavefunctions for use in helicity amplitude
 * calculations in the Herwig++. The general approach is to use a similar philosophy 
 * to the FORTRAN HELAS code but with additional structure.
 *
 * This class contains the storage of the particle type and 5-momentum 
 * and methods to set/access this information.
 *
 * The methods for the wavefunction itself will we implemented in the classes
 * derived from this one for the specific spin type, for example scalar, spinor,
 * vector and tensor. 
 *
 *  @see ScalarWaveFunction
 *  @see SpinorWaveFunction
 *  @see SpinorBarWaveFunction
 *  @see VectorWaveFunction
 *  @see TensorWaveFunction
 */
class WaveFunctionBase{

public:

  /**
   * Default constructor.
   */
  inline WaveFunctionBase();

  /**
   * Virtual destructor to keep compiler happy.
   */
  ~WaveFunctionBase();

  /**
   * Subscript operator to access momentum.
   */
  inline Energy operator [](int) const;

  /**
   * Get the x component of the momentum.
   */
  inline Energy px() const;

  /**
   * Get the y component of the momentum.
   */
  inline Energy py() const;

  /**
   * Get the z component of the momentum.
   */
  inline Energy pz() const;

  /**
   * Get the energy.
   */
  inline Energy e() const;

  /**
   * Get the mass.
   */
  inline Energy mass() const;

  /**
   * Get offshell mass squared.
   */
  inline Energy2 m2() const;

  /**
   * Set components of the momentum.
   */

  /**
   * Set the x component of the momentum.
   */
  inline void setPx(Energy);

  /**
   * Set the y component of the momentum.
   */
  inline void setPy(Energy);

  /**
   * Set the z component of the momentum.
   */
  inline void setPz(Energy);

  /**
   * Set the energy.
   */
  inline void setE(Energy);

  /**
   * Set the mass.
   */
  inline void setMass(Energy);

  /**
   * Set 5 momentum.
   */
  inline void setMomentum(const Lorentz5Momentum &);

  /**
   * Set all components of momentum.
   */
  inline void setMomentum(Energy,Energy,Energy,Energy,Energy);

  /** 
   * Set 4-momentum components.
   */
  inline void setMomentum(Energy,Energy,Energy,Energy);

  /** 
   * Set 4-momentum. 
   */
  inline void setMomentum(LorentzVector);

  /**
   * Set mass zero momentum.
   */
  inline void setMomentum(Energy);

  /**
   * Set 4 momentum and mass.
   */
  inline void setMomentum(LorentzVector,Energy);

  /**
   * Zero the 4 momentum.
   */
  inline void setMomentum();

  /** 
   * Get the particle id.
   */
  inline int id();

  /** 
   * Get 2s+1 for the particle.
   */
  inline int iSpin();

  /**
   * Get the particle pointer.
   */
  inline const tcPDPtr & getParticle() const;

  /** 
   * Direction of particle.
   */
  inline Direction direction();
  inline void direction(Direction);
  inline const Lorentz5Momentum & getMomentum() const ;

protected:

  /**
   * Set the particle pointer.
   */
  inline void setParticle(const tcPDPtr &);

private:

  /**
   * Check particle type and set pointer.
   */
  void checkParticle(const tcPDPtr &);

private:

  /**
   * Constant pointer to the particle info.
   */
  tcPDPtr _particle;

  /**
   * Lorentz 5 momentum.
   */
  Lorentz5Momentum _momentum;

  /**
   * Incoming or outgoing.
   */
  Direction _dir;
};
}
}

#include "WaveFunctionBase.icc"

#endif
