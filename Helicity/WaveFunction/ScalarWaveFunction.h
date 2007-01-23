// -*- C++ -*-
#ifndef HERWIG_ScalarWaveFunction_H
#define HERWIG_ScalarWaveFunction_H
//
// This is the declaration of the ScalarWaveFunction class.

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/ScalarSpinInfo.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/Helicity/RhoDMatrix.h>

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
using ThePEG::Helicity::tScalarSpinPtr;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::RhoDMatrix;

/** \ingroup Helicity
 *  \author Peter Richardson
 * 
 *  This class is the base class for scalar wavefunctions for use in 
 *  helicity amplitude calculations in the Herwig++. The general approach 
 *  is to use a similar philosophy to the FORTRAN HELAS code but with 
 *  additional structure.
 *
 *  This class stores the scalar wavefunction as a complex number and inherits
 *  from the WaveFunctionBase class for the storage of the particles
 *  momentum and type.
 * 
 *  @see WaveFunctionBase
 */
class ScalarWaveFunction : public WaveFunctionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor, set the momentum and Wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param wave The wavefunction.
   */
  inline ScalarWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,
			    Complex wave);

  /**
   * Constructor, set the momentum, direction and Wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param wave The wavefunction.
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,
			    Complex wave,Direction dir);

  /**
   * Constructor, set all components of momentum, mass, direction
   * and Wavefunction.
   * @param px The x-component of the momentum.
   * @param py The x-component of the momentum.
   * @param pz The x-component of the momentum.
   * @param E  The energy.
   * @param m  The mass.
   * @param part The ParticleData pointer
   * @param wave The wavefunction.
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(Energy px,Energy py,Energy pz,Energy E,Energy m,
			    const tcPDPtr & part,Complex wave,Direction dir);

  /**
   * Constructor, set all components of momentum, direction and Wavefunction.
   * @param px The x-component of the momentum.
   * @param py The x-component of the momentum.
   * @param pz The x-component of the momentum.
   * @param E  The energy.
   * @param part The ParticleData pointer
   * @param wave The wavefunction.
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(Energy px,Energy py,Energy pz,Energy E,const tcPDPtr & part,
			    Complex wave,Direction dir);

  /**
   * Constructor, set the momentum, direction and Wavefunction.
   * @param p The 4-momentum.
   * @param part The ParticleData pointer
   * @param wave The wavefunction.
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(LorentzVector p,const tcPDPtr & part,Complex wave,
			    Direction dir);

  /**
   * Constructor, set the mass, direction and Wavefunction.
   * @param m The mass.
   * @param part The ParticleData pointer
   * @param wave The wavefunction.
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(Energy m,const tcPDPtr & part,Complex wave,Direction dir);

  /**
   * Constructor, set the momentum, mass, direction and Wavefunction.
   * @param p The 4-momentum.
   * @param m The mass.
   * @param part The ParticleData pointer
   * @param wave The wavefunction.
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(LorentzVector p,Energy m,const tcPDPtr & part,Complex wave,
			    Direction dir);

  /**
   * Constructor,set the 5-momentum and zero the wavefunction.
   * @param p The 5-momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(Lorentz5Momentum p,const tcPDPtr & part,Direction dir); 

  /** 
   * Constructor, set all components of momentum, mass, direction and zero the
   * wavefunction.
   * @param px The x-component of the momentum.
   * @param py The x-component of the momentum.
   * @param pz The x-component of the momentum.
   * @param E  The energy.
   * @param m  The mass. 
   * @param part The ParticleData pointer
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(Energy px,Energy py,Energy pz,Energy E,Energy m,
			    const tcPDPtr & part, Direction dir);

  /**
   * Constructor, set all components of momentum, direction and zero the 
   * wavefunction.
   * @param px The x-component of the momentum.
   * @param py The x-component of the momentum.
   * @param pz The x-component of the momentum.
   * @param E  The energy.
   * @param part The ParticleData pointer
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(Energy px,Energy py,Energy pz,Energy E, 
			    const tcPDPtr & part,Direction dir);

  /**
   * Constructor, set the momentum, direction and zero the wavefunction.
   * @param p The 4-momentum.
   * @param part The ParticleData pointer
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(LorentzVector p,const tcPDPtr & part,Direction dir);

  /**
   * Constructor, set the mass, direction and zero the wavefunction.
   * @param m The mass.
   * @param part The ParticleData pointer
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(Energy m,const tcPDPtr & part,Direction dir);

  /**
   * Constructor, set the momentum, mass, direction and zero the wavefunction
   * @param p The 4-momentum.
   * @param m The mass.
   * @param part The ParticleData pointer
   * @param dir The direction of the particle.
   */
  inline ScalarWaveFunction(LorentzVector p,Energy m,const tcPDPtr & part,Direction dir);

  /**
   * Special constructor which set's up a particle's SpinInfo.
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the ScalarSpinInfo object
   */
  ScalarWaveFunction(tPPtr part,Direction dir,bool time, bool vertex);

  /**
   * Special constructor which set's up a particle's SpinInfo.
   * @param part The particle to setup
   * @param rho \f$\rho\f$ the rho matrix for the particle.
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the ScalarSpinInfo object
   */
  ScalarWaveFunction(tPPtr part,RhoDMatrix& rho,Direction dir,bool time,bool vertex);

  /**
   * Default constructor.
   */
  inline ScalarWaveFunction();

  /**
   * Subscript operator for the wavefunction.
   * This is provided for consistency with the other wavefunctions and returns
   * the wavefunction regardless of the index.
   */
  inline Complex operator ()(int ) const;

  /**
   * Set components by index.
   * This is provided for consistency with the other wavefunctions and sets
   * the wavefunction regardless of the index.
   */
  inline Complex & operator () (int);

  /**
   * Assignment. 
   */
  inline ScalarWaveFunction & operator = (const ScalarWaveFunction &);

  /**
   * Return the wavefunction.
   */
  inline const Complex & wave() const;

  /**
   * Functions to reset the wavefunction and momentum (to speed the code up).
   */
  //@{
  /**
   * Reset the momentum, particle type and direction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline void reset(const Lorentz5Momentum & p, const tcPDPtr & part, Direction dir);

  /**
   * Reset the momentum and particle type.
   * @param p The momentum.
   * @param dir The direction.
   */
  inline void reset(const Lorentz5Momentum & p,Direction dir);

  /**
   * Reset the momentum.
   * @param p The momentum.
   */
  inline void reset(const Lorentz5Momentum & p);

  /**
   * Reset the wavefunction.
   * @param wave The wavefunction
   */
  inline void reset(Complex wave);

  /**
   * Reset the particle type and direction.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline void reset(const tcPDPtr & part,Direction dir);

  /**
   * Reset the particle type.
   * @param part The ParticleData pointer.
   */
  inline void reset(const tcPDPtr & part);
  //@}  

private:
  
  /**
   * Check the particle type.
   * @param part The ParticleData pointer.
   */
  inline void checkParticle(const tcPDPtr & part);

private:

  /**
   * Complex number to store the wavefunction.
   */
  Complex _wf;

};
}
}

#include "ScalarWaveFunction.icc"

#endif
