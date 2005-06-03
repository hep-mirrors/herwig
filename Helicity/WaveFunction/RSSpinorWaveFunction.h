// -*- C++ -*-
#ifndef HERWIG_RSSpinorWaveFunction_H
#define HERWIG_RSSpinorWaveFunction_H
// This is the declaration of the RSSpinorWaveFunction class.

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzRSSpinor.h>
#include <ThePEG/Helicity/RSFermionSpinInfo.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/Helicity/RhoDMatrix.h>

namespace Herwig {

using ThePEG::Helicity::LorentzRSSpinor;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::defaultDRep;
using ThePEG::Helicity::tRSFermionSpinPtr;
using ThePEG::Helicity::RSFermionSpinInfo;
using ThePEG::Helicity::RhoDMatrix;

namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  The RSSpinorWaveFunction class is designed to store the wavefunction
 *  of a spin-3/2 particle in a form suitable for use in helicity amplitude
 *  calculations of the matrix element using a similar philosophy to the 
 *  FORTRAN HELAS code.
 *
 *  In addition to storing the spinor using the LorentzRSSpinor class
 *  it inherits from the <code>WaveFunctionBase</code> class to provide storage of
 *  the momentum and particleData for the fermion.
 *
 *  This class also contains the code which does the actually calculation of the
 *  spinor for an external particle using either of the Dirac matrix representations
 *  currently supported in the <code>HelicityDefinitions</code> class.
 *
 *  When calculating the wavefunction the direction of the particle is used,
 *
 *  \e i.e. 
 *  - incoming calculates a \f$u\f$ spinor.
 *  - outgoing calculates a \f$v\f$ spinor.
 *
 *  The spinors are calculated using a Clebsch-Gordon decomposition in the rest-frame
 *  for a massive particle and boosted to the lab-frame. For massless particles the
 *  calculation is performed in the lab-frame (N.B. there are only two helicities
 *  \f$\pm\frac32\f$ in this case.)
 *
 *  N.B. In our convention 0 is the \f$-\frac32\f$ helicity state,
 *        1 is the \f$-\frac12\f$ helicity state,
 *        2 is the \f$+\frac12\f$ helicity state
 *        3 is the \f$+\frac32\f$ helicity state and 
 *
 * @see WaveFunctionBase
 * @see LorentzRSSpinor
 * @see HelicityDefinitions
 * 
 */
class RSSpinorWaveFunction: public WaveFunctionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor, set the momentum and the components of the spinor and Dirac
   * matrix representation.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param xs1 The first  spinor component of the \f$x\f$ vector.
   * @param xs2 The second spinor component of the \f$x\f$ vector.
   * @param xs3 The third  spinor component of the \f$x\f$ vector.
   * @param xs4 The fourth spinor component of the \f$x\f$ vector.
   * @param ys1 The first  spinor component of the \f$y\f$ vector.
   * @param ys2 The second spinor component of the \f$y\f$ vector.
   * @param ys3 The third  spinor component of the \f$y\f$ vector.
   * @param ys4 The fourth spinor component of the \f$y\f$ vector.
   * @param zs1 The first  spinor component of the \f$z\f$ vector.
   * @param zs2 The second spinor component of the \f$z\f$ vector.
   * @param zs3 The third  spinor component of the \f$z\f$ vector.
   * @param zs4 The fourth spinor component of the \f$z\f$ vector.
   * @param ts1 The first  spinor component of the \f$t\f$ vector.
   * @param ts2 The second spinor component of the \f$t\f$ vector.
   * @param ts3 The third  spinor component of the \f$t\f$ vector.
   * @param ts4 The fourth spinor component of the \f$t\f$ vector.
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,
			      Complex xs1, Complex xs2, Complex xs3, Complex xs4,
			      Complex ys1, Complex ys2, Complex ys3, Complex ys4,
			      Complex zs1, Complex zs2, Complex zs3, Complex zs4,
			      Complex ts1, Complex ts2, Complex ts3, Complex ts4,
			      DiracRep drep=defaultDRep);

  /**
   * Constructor, set the momentum and the wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param wave The wavefunction.
   */
  inline RSSpinorWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,
			      LorentzRSSpinor & wave);

  /**
   * Constructor, set the momentum, helicity, direction and Dirac representation.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2,3 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,
			      unsigned int ihel,
			      Direction dir,DiracRep drep=defaultDRep);

  /**
   * Constructor, set the momentum components and mass, helicity, direction and
   * Dirac representation.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param m  The mass.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2,3 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(Energy px,Energy py,Energy pz,Energy E,Energy m,
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
   * @param ihel The helicity (0,1,2,3 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(Energy px,Energy py,Energy pz,Energy E,
			      const tcPDPtr & part,unsigned int ihel,Direction dir,
			      DiracRep drep=defaultDRep);

  /**
   * Constructor, set the 4-momentum, helicity, direction and
   * Dirac representation.
   * @param p the 4-momentum
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2,3 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(LorentzVector p,const tcPDPtr & part,unsigned int ihel,
			      Direction dir,DiracRep drep=defaultDRep);
  
  /**
   * Constructor, set the mass and zero the momentum, set the helicity, direction and
   * Dirac representation.
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2,3 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(Energy m,const tcPDPtr & part,unsigned int ihel,
			      Direction dir, DiracRep drep=defaultDRep);

  /**
   * Constructor, set the 4-momentum, mass, helicity, direction and
   * Dirac representation.
   * @param p the 4-momentum
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2,3 as described above.)
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(LorentzVector p,Energy m,const tcPDPtr & part,
			      unsigned int ihel,
			      Direction dir,DiracRep drep=defaultDRep);

  /**
   * Constructor, set the momentum, direction and Diracrepresentation, zero the 
   * wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(Lorentz5Momentum p,const tcPDPtr & part,Direction dir,
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
  inline RSSpinorWaveFunction(Energy px,Energy py,Energy pz,Energy E,Energy m,
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
  inline RSSpinorWaveFunction(Energy px,Energy py,Energy pz,Energy E,
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
  inline RSSpinorWaveFunction(LorentzVector p,const tcPDPtr & part,Direction dir,
			      DiracRep drep=defaultDRep);

  /**
   * Constructor set the mass, direction and
   * Dirac representation, zero the momentum and wavefunction.
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(Energy m,const tcPDPtr & part,Direction dir,
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
  inline RSSpinorWaveFunction(LorentzVector p,Energy m,const tcPDPtr & part,
			      Direction dir,DiracRep drep=defaultDRep);

  /**
   * Special constructor which calculates all the helicities and sets up a particle's
   * SpinInfo.
   * @param wave The spinors for the different helicities.
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the RSFermionSpinInfo object 
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(vector<LorentzRSSpinor>& wave, tPPtr part,Direction dir,
			      bool time, bool vertex, DiracRep drep=defaultDRep);

  /**
   * Special constructor which calculates all the helicities and sets up a particle's
   * SpinInfo.
   * @param wave The spinors for the different helicities.
   * @param rho The \f$\rho\f$ matrix for the particle
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the RSFermionSpinInfo object 
   * @param drep The Dirac representation.
   */
  inline RSSpinorWaveFunction(vector<LorentzRSSpinor>& wave, RhoDMatrix& rho,tPPtr part,
			      Direction dir,bool time, bool vertex,
			      DiracRep drep=defaultDRep);

  /**
   * Default constructor
   */
  inline RSSpinorWaveFunction(DiracRep=defaultDRep);

  /**
   * Destructor 
   */
  inline ~RSSpinorWaveFunction();
  //@}

  /**
   * Assignment. 
   */
  inline RSSpinorWaveFunction & operator = (const RSSpinorWaveFunction &);

  /**
   *  Access to the wavefunction and its components.
   */
  //@{
  /**
   * subscript operator for the wavefunction
   * Set components by index.
   */
  inline Complex operator ()(int,int ) const;
  /**
   * subscript operator for the wavefunction
   * Set components by index.
   */
  inline Complex & operator () (int,int);

  /**
   * return wavefunction as LorentzRSSpinor
   */
  inline LorentzRSSpinor Wave() const;

  /**
   * Get first spinor component for the x vector
   */
  inline Complex xs1() const;

  /**
   * Get second spinor component for the x vector
   */
  inline Complex xs2() const;

  /**
   * Get third  spinor component for the x vector
   */
  inline Complex xs3() const;

  /**
   * Get fourth  spinor component for the x vector
   */
  inline Complex xs4() const;

  /**
   * Get first spinor component for the y vector
   */
  inline Complex ys1() const;

  /**
   * Get second spinor component for the y vector
   */
  inline Complex ys2() const;
  
  /**
   * Get third spinor component for the y vector
   */
  inline Complex ys3() const;
  
  /**
   * Get fourth spinor component for the y vector
   */
  inline Complex ys4() const;
  
  /**
   * Get first spinor component for the z vector
   */
  inline Complex zs1() const;
  
  /**
   * Get second spinor component for the z vector
   */
  inline Complex zs2() const;
  
  /**
   * Get third spinor component for the z vector
   */
  inline Complex zs3() const;
  
  /**
   * Get fourth spinor component for the z vector
   */
  inline Complex zs4() const;
  
  /**
   * Get first spinor component for the t vector
   */
  inline Complex ts1() const;
  
  /**
   * Get second spinor component for the t vector
   */
  inline Complex ts2() const;
  
  /**
   * Get third spinor component for the t vector
   */
  inline Complex ts3() const;
  
  /**
   * Get fourth spinor component for the t vector
   */
  inline Complex ts4() const;
  
  /**
   * Set first spinor component for the x vector
   */
  inline void setXS1(Complex);
  
  /**
   * Set second spinor component for the x vector
   */
  inline void setXS2(Complex);
  
  /**
   * Set third spinor component for the x vector
   */
  inline void setXS3(Complex);
  
  /**
   * Set fourth spinor component for the x vector
   */
  inline void setXS4(Complex);
  
  /**
   * Set first spinor component for the y vector
   */
  inline void setYS1(Complex);
  
  /**
   * Set second spinor component for the y vector
   */
  inline void setYS2(Complex);
  
  /**
   * Set third spinor component for the y vector
   */
  inline void setYS3(Complex);
  
  /**
   * Set fourth spinor component for the y vector
   */
  inline void setYS4(Complex);
  
  /**
   * Set first spinor component for the z vector
   */
  inline void setZS1(Complex);
  
  /**
   * Set second spinor component for the z vector
   */
  inline void setZS2(Complex);
  
  /**
   * Set third spinor component for the z vector
   */
  inline void setZS3(Complex);
  
  /**
   * Set fourth spinor component for the z vector
   */
  inline void setZS4(Complex);
  
  /**
   * Set first spinor component for the t vector
   */
  inline void setTS1(Complex);
  
  /**
   * Set second spinor component for the t vector
   */
  inline void setTS2(Complex);
  
  /**
   * Set third spinor component for the t vector
   */
  inline void setTS3(Complex);
  
  /**
   * Set fourth spinor component for the t vector
   */
  inline void setTS4(Complex);
  //@}

  /**
   * reset functions
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
   * @param ihel The helicity (0,1,2,3 as described above.)
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
   * Calculate the spinors for all helicities, create and set up the SpinInfo object
   * @param wave The spinors for the different helicities.
   * @param part The particle to setup
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the RSFermionSpinInfo object 
   */
  inline void constructSpinInfo(vector<LorentzRSSpinor>& wave,tPPtr part,bool time,
				bool vertex=true);

  /**
   * Calculate the spinors for all helicities, create and set up the SpinInfo object
   * @param wave The spinors for the different helicities.
   * @param rho The \f$\rho\f$ matrix for the decaying particle.
   * @param part The particle to setup
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param vertex Whether or not to create the RSFermionSpinInfo object 
   */
  inline void constructSpinInfo(vector<LorentzRSSpinor>& wave,RhoDMatrix& rho,tPPtr part,
				bool time,bool vertex=true);

private:

  /**
   * Zero the wavefunction.
   */
  inline void zeroWaveFunction(DiracRep=defaultDRep);

  /**
   * Calcuate the wavefunction.
   * @param ihel The helicity (0,1,2,3 as described above.)
   * @param drep The Dirac matrix representation.
   */
  void calculateWaveFunction(unsigned int ihel,DiracRep drep=defaultDRep);

  /**
   * Check particle spin and set pointer.
   * @param part The ParticleData pointer.
   */
  inline void checkParticle(const tcPDPtr & part);

  /**
   * Calculate the spinors for all the helicities and set up the SpinInfo object.
   * @param wave The spinors for the different helicities
   * @param spin Pointer to the RSFermionSpinInfo object
   * @param vertex Whether or not to set up the RSFermionSpinInfo object 
   */
  inline void constructSpinInfo(vector<LorentzRSSpinor>& wave,tRSFermionSpinPtr spin,
				bool vertex=true);

private:

  /**
   * storage of the Lorentz RSSpinor
   */
  LorentzRSSpinor _wf;

};

}
}

#include "RSSpinorWaveFunction.icc"

#endif /* HERWIG_RSSpinorWaveFunction_H */
