// -*- C++ -*-
#ifndef HERWIG_ScalarWaveFunction_H
#define HERWIG_ScalarWaveFunction_H
//
// This is the declaration of the ScalarWaveFunction class.

#include "WaveFunctionBase.h"
namespace Herwig {
namespace Helicity {
using namespace ThePEG;

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

  /**
   * Default constructors (set the momentum and Wavefunction).
   */
  inline ScalarWaveFunction(const Lorentz5Momentum &,
			    const tcPDPtr &,Complex);

  /**
   * Use a 5-momentum.
   */
  inline ScalarWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,
			    Complex,Direction);

  /**
   * Set all components of momentum.
   */
  inline ScalarWaveFunction(Energy,Energy,Energy,Energy,Energy,
			    const tcPDPtr &,Complex,Direction);

  /**
   * Set 4-momentum components.
   */
  inline ScalarWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,
			    Complex,Direction);

  /**
   * Set 4-momentum.
   */
  inline ScalarWaveFunction(LorentzVector,const tcPDPtr &,Complex,Direction);

  /**
   * Set mass zero momentum.
   */
  inline ScalarWaveFunction(Energy,const tcPDPtr &,Complex,Direction);

  /**
   * Set 4 momentum and mass.
   */
  inline ScalarWaveFunction(LorentzVector,Energy,
			    const tcPDPtr &,Complex,Direction);

  /**
   * Default constructors (set the momentum and zero the Wavefunction)
   * use 5 momentum.
   */
  inline ScalarWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction); 

  /** 
   * Set all components of momentum.
   */
  inline ScalarWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			    Direction);

  /**
   * Set 4-momentum components.
   */
  inline ScalarWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction);

  /**
   * Set 4-momentum.
   */
  inline ScalarWaveFunction(LorentzVector,const tcPDPtr &,Direction);

  /**
   * Set mass zero momentum.
   */
  inline ScalarWaveFunction(Energy,const tcPDPtr &,Direction);

  /**
   * Set 4 momentum and mass.
   */
  inline ScalarWaveFunction(LorentzVector,Energy,const tcPDPtr &,Direction);

  /**
   * Default constructor.
   */
  inline ScalarWaveFunction();

  /**
   * Destructor.
   */
  inline ~ScalarWaveFunction();

  /**
   * Subscript operator for the wavefunction.
   */
  inline Complex operator ()(int ) const;

  /**
   * Set components by index..
   */
  inline Complex & operator () (int);

  /**
   * Assignment. 
   */
  inline ScalarWaveFunction & operator = (const ScalarWaveFunction &);

  /**
   * Return the wavefunction.
   */
  inline Complex Wave() const;

  /**
   * Functions to reset the wavefunction and momentum (to speed the code up).
   */

  /**
   * Reset the momentum, particle type and direction.
   */
  inline void reset(const Lorentz5Momentum &, const tcPDPtr &, Direction);

  /**
   * Reset the momentum and particle type.
   */
  inline void reset(const Lorentz5Momentum &,Direction);

  /**
   * Reset the momentum.
   */
  inline void reset(const Lorentz5Momentum &);

  /**
   * Reset the wavefunction.
   */
  inline void reset(Complex);

  /**
   * Reset the particle type and direction.
   */
  inline void reset(const tcPDPtr &,Direction);

  /**
   * Reset the particle type.
   */
  inline void reset(const tcPDPtr &);
  
private:
  
  /**
   * Check the particle type.
   */
  inline void checkParticle(const tcPDPtr &);

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
