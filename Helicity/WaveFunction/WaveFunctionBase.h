// -*- C++ -*-
#ifndef HERWIG_WaveFunctionBase_H
#define HERWIG_WaveFunctionBase_H
//
// This is the declaration of the WaveFunctionBase class.

#include <ThePEG/CLHEPWrap/Lorentz5Vector.h>
#include <ThePEG/CLHEPWrap/LorentzVector.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/Config/Complex.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *  Definition of the enumerated values used for the direction of the 
 *  particles in the calculation of the wavefunction.
 */
enum Direction 
{
  incoming, /**< An incoming particle. */
  outgoing, /**< An outgoing particle. */
  intermediate /**< An intermediate particle. */
};

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
 * The methods for the wavefunction itself will be implemented in the classes
 * derived from this one for the specific spin type, for example scalar, spinor,
 * vector and tensor. 
 *
 *  @see ScalarWaveFunction
 *  @see SpinorWaveFunction
 *  @see SpinorBarWaveFunction
 *  @see VectorWaveFunction
 *  @see RSSpinorWaveFunction
 *  @see RSSpinorBarWaveFunction
 *  @see TensorWaveFunction
 */
class WaveFunctionBase{

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline WaveFunctionBase();

  /**
   * Destructor.
   */
  ~WaveFunctionBase();
  //@}

  /**
   * Access to the momentum components and mass
   */
  //@{
  /**
   * Subscript operator to access momentum.
   * @param iloc The required component 0,1,2 is \f$p_x\f$, \f$p_y\f$, \f$p_z\f$,
   * 3 is \f$E\f$ and 4 is the mass.
   */
  //inline Energy operator [](int iloc) const;

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
   * Get off-shell mass squared.
   */
  inline Energy2 m2() const;

  /**
   *  Access to the 5-momentum
   */
  inline const Lorentz5Momentum & getMomentum() const ;
  //@}

  /**
   * Set the components of the momentum and the mass
   */
  //@{
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
   * @param px The x-component of the momentum.
   * @param py The x-component of the momentum.
   * @param pz The x-component of the momentum.
   * @param E  The energy.
   * @param m  The mass.
   */
  inline void setMomentum(Energy px,Energy py,Energy pz,Energy E,Energy m);

  /** 
   * Set 4-momentum components.
   * @param px The x-component of the momentum.
   * @param py The x-component of the momentum.
   * @param pz The x-component of the momentum.
   * @param E  The energy.
   */
  inline void setMomentum(Energy px,Energy py,Energy pz,Energy E);

  /** 
   * Set 4-momentum using a vector.
   * @param p The momentum.
   */
  inline void setMomentum(LorentzMomentum p);

  /**
   * Set mass and zero momentum.
   * @param m The mass
   */
  inline void setMomentum(Energy m);

  /**
   * Set 4 momentum and mass.
   * @param p The momentum.
   * @param m The mass
   */
  inline void setMomentum(LorentzMomentum p,Energy m);

  /**
   * Zero the 4 momentum and mass.
   */
  inline void setMomentum();
  //@}

  /**
   *  Access to the particle properties
   */
  //@{
  /** 
   * Get the particle id.
   */
  inline int id();

  /** 
   * Get 2s+1 for the particle.
   */
  inline PDT::Spin iSpin();

  /**
   * Get the particle pointer.
   */
  inline const tcPDPtr & getParticle() const;

  /** 
   * Get the direction of particle.
   */
  inline Herwig::Helicity::Direction direction() const;

  /**
   * Set the direction of the particle
   */
  inline void direction(Herwig::Helicity::Direction);
  //@}

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
