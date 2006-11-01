// -*- C++ -*-
#ifndef HERWIG_VectorWaveFunction_H
#define HERWIG_VectorWaveFunction_H
//
// This is the declaration of the VectorWaveFunction class.
//
#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzPolarizationVector.h>
#include <ThePEG/Helicity/VectorSpinInfo.h>
#include <ThePEG/Helicity/RhoDMatrix.h>
#include <ThePEG/EventRecord/Particle.h>

namespace Herwig {
namespace Helicity {

/** \ingroup Helicity
 * Definition of the enumerated values of the phase to include in the 
 * calculation of the polarization vector.
 */
enum VectorPhase 
{
  vector_phase, /**< Include the phase factor.*/
  vector_nophase, /**< No phase-factor. */
  default_vector_phase=vector_nophase /**< Default option.*/
};

using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::tVectorSpinPtr;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::RhoDMatrix;

/** \ingroup Helicity
 *
 *  \author Peter Richardson
 *
 *  The VectorWaveFunction class is designed to store the wavefunction
 *  of a vector in a form suitable for use in helicity amplitude calculations 
 *  of the matrix element using a similar philosophy to the FORTRAN HELAS code.
 *
 *  In addition to storing the vector using the LorentzPolarizationVector class
 *  it inherits from the WaveFunctionBase class to provide storage of the 
 *  momentum and particleData for the vector boson.
 *
 *  This class also contains the code which does the actually calculation of the
 *  vector wavefunction.
 *
 *  There are two choices available for the calculation of the wavefunction.
 *  These are set using the VectorPhase enumeration which specifies a default choice.
 *  The first choice, vector_phase, includes a phase factor \f$\exp(\pm i \phi)\f$
 *  for the \f$\pm\f$ helicity states while the second, vector_nophase, does not.
 *
 *  N.B. In our convention 0 is the \f$-1\f$ helicity state and 
 *        1 is the \f$0\f$ helicity state
 *        2 is the \f$+1\f$ helicity state
 *
 *  @see WaveFunctionBase
 *  @see LorentzPolarizationVector
 */
class VectorWaveFunction : public WaveFunctionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor, set the momentum and Wavefunction, the direction can also
   * be specified. 
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param wave The wavefunction, \e i.e. the polarization vector.
   * @param dir The direction of the particle.
   */
  inline VectorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
			    const LorentzPolarizationVector & wave,
			    Direction  dir=intermediate);

  /**
   * Constructor, set the momentum and components of the wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param x The x component of the polarization vector
   * @param y The y component of the polarization vector
   * @param z The z component of the polarization vector
   * @param t The t component of the polarization vector
   */
  inline VectorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,const Complex & x,
			    const Complex & y,const Complex & z, const Complex & t);
  
  /**
   * Constructor, set the momentum, helicity and direction, optionally the choice
   * of the phase.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2 as described above.)
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline VectorWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,
			    unsigned int ihel,Direction dir,
			    VectorPhase phase=default_vector_phase);

  /**
   * Constructor, set the momentum components and mass, helicity and direction,
   * optionally the choice of the phase.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param m  The mass.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2 as described above.)
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline VectorWaveFunction(Energy px,Energy py,Energy pz,Energy E,Energy m,
			    const tcPDPtr & part,unsigned int ihel,Direction dir,
			    VectorPhase phase=default_vector_phase);

  /**
   * Constructor, set the momentum components, helicity and direction,
   * optionally the choice of the phase.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2 as described above.)
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline VectorWaveFunction(Energy px,Energy py,Energy pz,Energy E,const tcPDPtr & part,
			    unsigned int ihel,Direction dir,
			    VectorPhase phase=default_vector_phase);

  /**
   * Constructor, set the 4-momentum, helicity and direction,
   * optionally the choice of the phase.
   * @param p The 4-momentum.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2 as described above.)
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline VectorWaveFunction(LorentzVector p,const tcPDPtr & part,
			    unsigned int ihel,Direction dir,
			    VectorPhase phase=default_vector_phase);

  /**
   * Constructor, set the mass, zero the momentum and set the helicity and direction,
   * optionally the choice of the phase.
   * @param m The mass. 
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2 as described above.)
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline VectorWaveFunction(Energy m,const tcPDPtr & part,
			    unsigned int ihel,Direction dir,
			    VectorPhase phase=default_vector_phase);

  /**
   * Constructor, set the 4-momentum, mass, helicity and direction,
   * optionally the choice of the phase.
   * @param p The 4-momentum.
   * @param m The mass. 
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2 as described above.)
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline VectorWaveFunction(LorentzVector p,Energy m,const tcPDPtr & part,
			    unsigned int ihel,
			    Direction dir,VectorPhase phase=default_vector_phase);

  /**
   * Constructor, set the 5-momentum and direction, zero the wavefunction.
   * @param p The 5-momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline VectorWaveFunction(Lorentz5Momentum p,const tcPDPtr & part,Direction dir); 

  /**
   * Constructor, set the momentum components, mass and direction,
   * zero the wavefunction.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param m  The mass.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline VectorWaveFunction(Energy px,Energy py,Energy pz,Energy E,Energy m,
			    const tcPDPtr & part,Direction dir);

  /**
   * Constructor, set the momentum components and direction,
   * zero the wavefunction.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline VectorWaveFunction(Energy px,Energy py,Energy pz,Energy E,const tcPDPtr & part,
			    Direction dir);

  /**
   * Constructor, set the 4-momentum  and direction,
   * zero the wavefunction.
   * @param p The 4-momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline VectorWaveFunction(LorentzVector p,const tcPDPtr & part,Direction dir);

  /**
   * Constructor, set the mass   and direction,
   * zero the wavefunction and momentum.
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline VectorWaveFunction(Energy m,const tcPDPtr & part,Direction dir);

  /**  
   * Constructor, set the 4-momentum, mass  and direction,
   * zero the wavefunction.
   * @param p The 4-momentum.
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline VectorWaveFunction(LorentzVector p,Energy m,const tcPDPtr & part,Direction dir);

  /**
   * Special constructor which calculates all the polarization vectors,
   * as LorentzPolarizationVector objects, for all helicities and sets up a particle's
   * SpinInfo.
   * @param wave The polarization vectors for the different helicities.
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param massless Whether or not the particle is massless
   * @param vertex Whether or not to create the VectorSpinInfo object 
   * @param phase The phase choice.
   */
  inline VectorWaveFunction(vector<LorentzPolarizationVector>& wave, tPPtr part,
			    Direction dir, bool time, bool massless, bool vertex,
			    VectorPhase phase=default_vector_phase);

  /**
   * Special constructor which calculates all the polarization vectors,
   * as LorentzPolarizationVector objects, for all helicities and sets up a particle's
   * SpinInfo.
   * @param wave The polarization vectors for the different helicities.
   * @param rho The \f$\rho\f$ matrix for the particle.
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param massless Whether or not the particle is massless
   * @param vertex Whether or not to create the VectorSpinInfo object 
   * @param phase The phase choice.
   */

  inline VectorWaveFunction(vector<LorentzPolarizationVector>& wave, RhoDMatrix& rho,
			    tPPtr part,Direction dir, bool time, bool massless,
			    bool vertex,VectorPhase phase=default_vector_phase);

  /**
   * Special constructor which calculates all the polarization vectors,
   * as VectorWaveFunction objects, for all helicities and sets up a particle's
   * SpinInfo.
   * @param wave The polarization vectors for the different helicities.
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param massless Whether or not the particle is massless
   * @param vertex Whether or not to create the VectorSpinInfo object 
   * @param phase The phase choice.
   */
  inline VectorWaveFunction(vector<VectorWaveFunction>& wave, tPPtr part,
			    Direction dir, bool time, bool massless, bool vertex,
			    VectorPhase phase=default_vector_phase);

  /**
   * Special constructor which calculates all the polarization vectors,
   * as LorentzPolarizationVector objectsm for all helicities and sets up a particle's
   * SpinInfo.
   * @param wave The polarization vectors for the different helicities.
   * @param rho The \f$\rho\f$ matrix for the particle.
   * @param part The particle to setup
   * @param dir The direction.
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param massless Whether or not the particle is massless
   * @param vertex Whether or not to create the VectorSpinInfo object 
   * @param phase The phase choice.
   */
  inline VectorWaveFunction(vector<VectorWaveFunction>& wave, RhoDMatrix& rho,
			    tPPtr part,Direction dir, bool time, bool massless,
			    bool vertex,VectorPhase phase=default_vector_phase);

  /**
   * Default constructor.
   */
  inline VectorWaveFunction();

  /**
   * Destructor.
   */
  inline ~VectorWaveFunction();
  //@}


  /**
   * Assignment. 
   */
  inline VectorWaveFunction & operator = (const VectorWaveFunction &);

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
   * Return wavefunction as polarization vector. 
   */
  inline const LorentzPolarizationVector & wave() const;
  
  /**
   * Get x component.
   */
  inline Complex x() const;
  
  /**
   * Get y component.
   */
  inline Complex y() const;
  
  /**
   * Get z component.
   */
  inline Complex z() const;
  
  /**
   * Get t component.
   */
  inline Complex t() const;
  
  /**
   * Set x component
   */
  inline void setX(const Complex&);
  
  /**
   * Set y component
   */
  inline void setY(const Complex&);
  
  /**
   * Set z component
   */
  inline void setZ(const Complex&);
  
  /**
   * Set t component
   */
  inline void setT(const Complex&);
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
   * Reset the helicity (recalculation the polarization vector).
   * @param ihel The new helicity (0,1,2 as described above.)
   * @param phase The phase choice.
   */
  inline void reset(unsigned int ihel,VectorPhase phase=default_vector_phase);

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

  /**
   * Calculate the polarization vectors, as LorentzPolarizationVector objects,
   * for all helicities, create and set up the SpinInfo object
   * @param wave The polarization vectors for the different helicities.
   * @param part The particle to setup
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param massless Whether or not the particle is massless
   * @param phase The phase choice.
   * @param vertex Whether or not to create the VectorSpinInfo object 
   */
  inline void constructSpinInfo(vector<LorentzPolarizationVector>& wave,tPPtr part,
				bool time, bool massless,
				VectorPhase phase=default_vector_phase,bool vertex=true);

  /**
   * Calculate the polarization vectors, as LorentzPolarizationVector objects,
   * for all helicities, create and set up the SpinInfo object
   * @param wave The polarization vectors for the different helicities.
   * @param rho The \f$\rho\f$ matrix for the particle
   * @param part The particle to setup
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param massless Whether or not the particle is massless
   * @param phase The phase choice.
   * @param vertex Whether or not to create the VectorSpinInfo object 
   */
  inline void constructSpinInfo(vector<LorentzPolarizationVector>& wave,RhoDMatrix& rho,
				tPPtr part,bool time, bool massless,
				VectorPhase phase=default_vector_phase,bool vertex=true);

  /**
   * Calculate the polarization vectors, as VectorWaveFunction objects,
   * for all helicities, create and set up the SpinInfo object
   * @param wave The polarization vectors for the different helicities.
   * @param part The particle to setup
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param massless Whether or not the particle is massless
   * @param phase The phase choice.
   * @param vertex Whether or not to create the VectorSpinInfo object 
   */
  inline void constructSpinInfo(vector<VectorWaveFunction>& wave,tPPtr part,
				bool time, bool massless,
				VectorPhase phase=default_vector_phase,bool vertex=true);

  /**
   * Calculate the polarization vectors, as VectorWaveFunction objects,
   * for all helicities, create and set up the SpinInfo object
   * @param wave The polarization vectors for the different helicities.
   * @param rho The \f$\rho\f$ matrix for the particle
   * @param part The particle to setup
   * @param time Is this is timelike (true) or spacelike (false ) particle?
   * @param massless Whether or not the particle is massless
   * @param phase The phase choice.
   * @param vertex Whether or not to create the VectorSpinInfo object 
   */
  inline void constructSpinInfo(vector<VectorWaveFunction>& wave,RhoDMatrix& rho,
				tPPtr part,bool time, bool massless,
				VectorPhase phase=default_vector_phase,bool vertex=true);

private:

  /** 
   * Zero the wavefunction.
   */
  inline void zeroWaveFunction();

  
  /**
   * Calculate the wavefunction
   * @param ihel The helicity  (0,1,2 as described above.)
   * @param phase The phase choice.
   */
  void calculateWaveFunction(unsigned int ihel,VectorPhase phase=default_vector_phase);

  /**
   * Check the particle type.
   * @param part The ParticleData pointer.
   */
  inline void checkParticle(const tcPDPtr & part);

  /**
   * Calculate the polarization vectors, as LorentzPolarizationVector objects,
   * for all the helicities and set up the
   * SpinInfo object.
   * @param wave The polarization vectors for the different helicities
   * @param spin Pointer to the VectorSpinInfo object
   * @param massless Whether or not the particle is massless
   * @param phase The phase choice.
   * @param vertex Whether or not to set up the VectorSpinInfo object 
   */
  inline void constructSpinInfo(vector<LorentzPolarizationVector>& wave,
				tVectorSpinPtr spin, bool massless,
				VectorPhase phase=default_vector_phase,bool vertex=true);

  /**
   * Calculate the polarization vectors, as VectorWaveFunction objects,
   * for all the helicities and set up the SpinInfo object.
   * @param wave The polarization vectors for the different helicities
   * @param spin Pointer to the VectorSpinInfo object
   * @param massless Whether or not the particle is massless
   * @param phase The phase choice.
   * @param vertex Whether or not to set up the VectorSpinInfo object 
   */
  inline void constructSpinInfo(vector<VectorWaveFunction>& wave,
				tVectorSpinPtr spin, bool massless,
				VectorPhase phase=default_vector_phase,bool vertex=true);

private:
  
  /**
   * Storage of the wavefunction as a Lorentz Vector.
   */
  LorentzPolarizationVector _wf;
  
};

}
}

#include "VectorWaveFunction.icc"

#endif
