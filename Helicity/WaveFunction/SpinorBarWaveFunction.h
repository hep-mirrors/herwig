// -*- C++ -*-
#ifndef HERWIG_SpinorBarWaveFunction_H
#define HERWIG_SpinorBarWaveFunction_H
//
// This is the declaration of the SpinorBarWaveFunction class.

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzSpinorBar.h>
#include <ThePEG/Helicity/FermionSpinInfo.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/Helicity/RhoDMatrix.h>

namespace Herwig {

using ThePEG::Helicity::LorentzSpinorBar;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::defaultDRep;
using ThePEG::Helicity::tFermionSpinPtr;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::RhoDMatrix;

namespace Helicity {

using namespace ThePEG;

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  The SpinorBarWaveFunction class is designed to store the wavefunction
 *  of a barred spinor in a form suitable for use in helicity amplitude 
 *  calculations of the matrix element using a similar philosophy to the 
 *  FORTRAN HELAS code.
 *
 *  In addition to storing the spinor using the LorentzSpinorBar class
 *  it inherits from the WaveFunctionBase class to provide storage of
 *  the momentum and particleData for the fermion.
 *
 *  This class also contains the code which does the actually calculation 
 *  of the barred spinor for an external particle using either of the 
 *  Dirac matrix representations currently supported in the 
 *  HelicityDefinitions class.
 *
 *  When calculating the wavefunction the direction of the particle is used,
 *
 *  \e i.e. 
 *  - incoming calculates a \f$\bar{v}\f$ spinor.
 *  - outgoing calculates a \f$\bar{u}\f$ spinor.
 *
 *  N.B. In our convention 0 is the \f$-\frac12\f$ helicity state and 
 *        1 is the \f$+\frac12\f$ helicity state
 *
 *  @see WaveFunctionBase
 *  @see LorentzSpinorBar
 *  @see HelicityDefinitions
 */
class SpinorBarWaveFunction : public WaveFunctionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor, set the momentum and the components of the spinor and Dirac
   * matrix representation/
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param s1 The first component
   * @param s2 The second component
   * @param s3 The third component
   * @param s4 The fourth component
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,
			       Complex s1,Complex s2,Complex s3,Complex s4,
			       DiracRep drep=defaultDRep);
  
  
  /**
   * Constructor, set the momentum and the wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param wave The wavefunction.
   * @param dir The direction of the particle
   */
  inline SpinorBarWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,
			       const LorentzSpinorBar & wave,Direction dir=intermediate);
  
  /**
   * Constructor, set the momentum, helicity, direction and Dirac representation.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,
			       unsigned int ihel,Direction dir,
			       DiracRep drep=defaultDRep);
  
  /**
   * Constructor, set the momentum components and mass, helicity, direction and
   * Dirac representation.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param m  The mass.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(Energy px,Energy py,Energy pz,Energy E,Energy m,
			       const tcPDPtr & part,unsigned int ihel,Direction dir,
			       DiracRep drep=defaultDRep);
  
  /**
   * Constructor, set the momentum components, helicity, direction and
   * Dirac representation.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(Energy px,Energy py,Energy pz,Energy E,
			       const tcPDPtr & part,unsigned int ihel,Direction dir,
			       DiracRep drep=defaultDRep);
  
  /**
   * Constructor, set the 4-momentum, helicity, direction and
   * Dirac representation.
   * @param p the 4-momentum
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(LorentzVector p,const tcPDPtr & part,unsigned int ihel,
			       Direction dir,DiracRep drep=defaultDRep);
  
  /**
   * Constructor, set the mass and zero the momentum, set the helicity, direction and
   * Dirac representation.
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(Energy m,const tcPDPtr & part,unsigned int ihel,
			       Direction dir, DiracRep drep=defaultDRep);
  
  /**
   * Constructor, set the 4-momentum, mass, helicity, direction and
   * Dirac representation.
   * @param p the 4-momentum
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(LorentzVector p,Energy m,const tcPDPtr & part,
			       unsigned int ihel,Direction dir,
			       DiracRep drep=defaultDRep);

  /**
   * Constructor, set the momentum, direction and Diracrepresentation, zero the 
   * wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(Lorentz5Momentum p,const tcPDPtr & part,Direction dir,
			       DiracRep drep=defaultDRep); 

  /**
   * Constructor, set the momentum components, mass, direction and
   * Dirac representation, zero the wavefunction.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param m  The mass.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(Energy px,Energy py,Energy pz,Energy E,Energy m,
			       const tcPDPtr & part,Direction dir,
			       DiracRep drep=defaultDRep);

  /**
   * Constructor, set the momentum components, direction and
   * Dirac representation, zero the wavefunction.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(Energy px,Energy py,Energy pz,Energy E,
			       const tcPDPtr & part,Direction dir,
			       DiracRep drep=defaultDRep);
  
  /**
   * Constructor set the 4-momentum, direction and
   * Dirac representation, zero the wavefunction.
   * @param p The 4-momentum
   * @param part The ParticleData pointer.
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(LorentzVector p,const tcPDPtr & part,Direction dir,
			       DiracRep drep=defaultDRep);
  
  /**
   * Constructor set the mass, direction and
   * Dirac representation, zero the momentum and wavefunction.
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(Energy m,const tcPDPtr & part,Direction dir,
			       DiracRep drep=defaultDRep);
  
  /**
   * Constructor set the 4-momentum, mass, direction and
   * Dirac representation, zero the wavefunction.
   * @param p The 4-momentum
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(LorentzVector p,Energy m,const tcPDPtr & part,
			       Direction dir,DiracRep drep=defaultDRep);

  /**
   * Special constructor which calculates all the helicities, as LorentzSpinorBar
   * objects, and sets up a particle's SpinInfo.
   * @param wave The spinors for the different helicities.
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the FermionSpinInfo object 
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(vector<LorentzSpinorBar>& wave, tPPtr part,Direction dir,
			       bool time, bool vertex, DiracRep drep=defaultDRep);

  /**
   * Special constructor which calculates all the helicities, as LorentzSpinorBar
   * objects, and sets up a particle's SpinInfo.
   * @param wave The spinors for the different helicities.
   * @param rho  The \f$\rho\f$ matrix for the particle
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the FermionSpinInfo object 
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(vector<LorentzSpinorBar>& wave,RhoDMatrix& rho,
			       tPPtr part,Direction dir,bool time, bool vertex,
			       DiracRep drep=defaultDRep);
  
  /**
   * Special constructor which calculates all the helicities, as SpinorBarWaveFunction
   * objects, and sets up a particle's SpinInfo.
   * @param wave The spinors for the different helicities.
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the FermionSpinInfo object 
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(vector<SpinorBarWaveFunction>& wave, tPPtr part,
			       Direction dir,bool time, bool vertex,
			       DiracRep drep=defaultDRep);

  /**
   * Special constructor which calculates all the helicities, as SpinorBarWaveFunction
   * objects, and sets up a particle's SpinInfo.
   * @param wave The spinors for the different helicities.
   * @param rho  The \f$\rho\f$ matrix for the particle
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the FermionSpinInfo object 
   * @param drep The Dirac representation.
   */
  inline SpinorBarWaveFunction(vector<SpinorBarWaveFunction>& wave,RhoDMatrix& rho,
			       tPPtr part,Direction dir,bool time, bool vertex,
			       DiracRep drep=defaultDRep);

  /**
   * Default constructor.
   */
  inline SpinorBarWaveFunction(DiracRep=defaultDRep);
  
  /**
   * Destructor.
   */
  inline ~SpinorBarWaveFunction();
  //@}

  /**
   * Assignment. 
   */
  inline SpinorBarWaveFunction & operator = (const SpinorBarWaveFunction &);

  /**
   *  Access to the wavefunction and its components.
   */
  //@{
  /**
   * Subscript operator for the wavefunction.
   */
  inline Complex operator ()(int ) const;

  /**
   * Set components by index.
   */
  inline Complex & operator () (int);

  /**
   * Return wavefunction as LorentzSpinor.
   */
  inline LorentzSpinorBar Wave() const;

  /**
   * Get the first spin component component.
   */
  inline Complex s1() const;

  /**
   * Get the second spin component component.
   */
  inline Complex s2() const;

  /**
   * Get the third spin component component.
   */
  inline Complex s3() const;

  /**
   * Get the fourth spin component component.
   */
  inline Complex s4() const;

  /**
   * Set first spin component.
   */
  inline void setS1(Complex);

  /**
   * Set second spin component.
   */
  inline void setS2(Complex);

  /**
   * Set third spin component.
   */
  inline void setS3(Complex);

  /**
   * Set fourth spin component.
   */
  inline void setS4(Complex);
  //@}

  /**
   * Reset functions.
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
   * Reset the momentum and direction
   * @param p The momentum.
   * @param dir The direction
   */
  inline void reset(const Lorentz5Momentum & p,Direction dir);

  /**
   * Reset the momentum.
   * @param p The momentum.
   */
  inline void reset(const Lorentz5Momentum & p);

  /**
   * Reset the helicity (calculates the new spinor).
   * @param ihel The helicity (0,1 as described above.)
   * @param drep The Dirac matrix representation.
   */
  inline void reset(unsigned int ihel,DiracRep drep=defaultDRep);

  /**
   * Reset particle type and direction.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline void reset(const tcPDPtr & part,Direction dir);

  /**
   * Reset particle type.
   * @param part The ParticleData pointer.
   */
  inline void reset(const tcPDPtr & part);
  //@}

  /**
   * Calculate the spinors, as LorentzSpinorBar objects,
   * for all helicities, create and set up the SpinInfo object
   * @param wave The spinors for the different helicities.
   * @param part The particle to setup
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the FermionSpinInfo object 
   */
  inline void constructSpinInfo(vector<LorentzSpinorBar>& wave,tPPtr part,bool time,
				bool vertex=true);

  /**
   * Calculate the spinors, as LorentzSpinorBar objects,
   * for all helicities, create and set up the SpinInfo object
   * @param wave The spinors for the different helicities.
   * @rho The \f$\rho\f$ matrix for the particle
   * @param part The particle to setup
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the FermionSpinInfo object 
   */
  inline void constructSpinInfo(vector<LorentzSpinorBar>& wave,RhoDMatrix& rho,
				tPPtr part,bool time,bool vertex=true);

  /**
   * Calculate the spinors, as SpinorBarWaveFunction objects,
   *  for all helicities, create and set up the SpinInfo object
   * @param wave The spinors for the different helicities.
   * @param part The particle to setup
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the FermionSpinInfo object 
   */
  inline void constructSpinInfo(vector<SpinorBarWaveFunction>& wave,tPPtr part,bool time,
				bool vertex=true);

  /**
   * Calculate the spinors, as SpinorBar objects,
   *  for all helicities, create and set up the SpinInfo object
   * @param wave The spinors for the different helicities.
   * @rho The \f$\rho\f$ matrix for the particle
   * @param part The particle to setup
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the FermionSpinInfo object 
   */
  inline void constructSpinInfo(vector<SpinorBarWaveFunction>& wave,RhoDMatrix& rho,
				tPPtr part,bool time,bool vertex=true);
private:

  /**
   * Zero the wavefunction.
   */
  inline void zeroWaveFunction(DiracRep=defaultDRep);

  /**
   * Calcuate the wavefunction.
   * @param ihel The helicity (0,1 as described above.)
   * @param drep The Dirac matrix representation.
   */
  void calculateWaveFunction(unsigned int ihel,DiracRep drep=defaultDRep);

  /**
   * Check particle spin and set pointer.
   * @param part The ParticleData pointer.
   */
  inline void checkParticle(const tcPDPtr & part);

  /**
   * Calculate the spinors, as LorentzSpinorBar objects,
   * for all the helicities and set up the SpinInfo object.
   * @param wave The spinors for the different helicities
   * @param spin Pointer to the FermionSpinInfo object
   * @param vertex Whether or not to set up the FermionSpinInfo object 
   */
  inline void constructSpinInfo(vector<LorentzSpinorBar>& wave,tFermionSpinPtr spin,
				bool vertex=true);
  
  /**
   * Calculate the spinors, as SpinorBarWaveFunction objects,
   *  for all the helicities and set up the SpinInfo object.
   * @param wave The spinors for the different helicities
   * @param spin Pointer to the FermionSpinInfo object
   * @param vertex Whether or not to set up the FermionSpinInfo object 
   */
  inline void constructSpinInfo(vector<SpinorBarWaveFunction>& wave,tFermionSpinPtr spin,
				bool vertex=true);
private:
  
  /**
   * Storage of the Lorentz SpinorBar wavefunction.
   */
  LorentzSpinorBar _wf;
  
};
}
}

#include "SpinorBarWaveFunction.icc"

#endif




