// -*- C++ -*-
#ifndef HERWIG_SpinorBarWaveFunction_H
#define HERWIG_SpinorBarWaveFunction_H
//
// This is the declaration of the SpinorBarWaveFunction class.

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzSpinorBar.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>

namespace Herwig {

using ThePEG::Helicity::LorentzSpinorBar;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::defaultDRep;

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
 *  i.e. ipart=-1 (incoming) calculates a v-bar spinor
 *       ipart=+1 (outgoing) calculates a u-bar spinor
 *
 *  @see WaveFunctionBase
 *  @see LorentzSpinorBar
 *  @see HelicityDefinitions
 */
class SpinorBarWaveFunction : public WaveFunctionBase {

public:

  /**
   * Default constructors (set the momentum and Wavefunction).
   */
  
  /**
   * Use a 5-momentum and specify all components (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,Complex,
			       Complex,Complex,Complex,DiracRep=defaultDRep);
  
  /**
   * Use a 5-momentum and a LorentzSpinorBar.
   */
  inline SpinorBarWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,
			       LorentzSpinorBar &);
  
  /**
   * Use a 5-momentum (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,
			       int,Direction,DiracRep=defaultDRep);
  
  /**
   * Set all components of momentum (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,Energy,
			       const tcPDPtr &,int,Direction,DiracRep=defaultDRep);
  
  /**
   * Set 4-momentum components (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,
			       int,Direction,DiracRep=defaultDRep);
  
  /**
   * Set 4-momentum (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(LorentzVector,const tcPDPtr &,int,Direction,
			       DiracRep=defaultDRep);
  
  /**
   * Set mass zero momentum (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(Energy,const tcPDPtr &,int,Direction,
			       DiracRep=defaultDRep);
  
  /**
   * Set 4 momentum and mass (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(LorentzVector,Energy,const tcPDPtr &,int,Direction,
			       DiracRep=defaultDRep);
  
  /**
   * Default constructors (set the momentum and zero the Wavefunction).
   */

  /**
   * Use 5 momentum (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction,
			       DiracRep=defaultDRep); 

  /**
   * Set all components of momentum (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,Energy,
			       const tcPDPtr &,Direction,DiracRep=defaultDRep);

  /**
   * Set 4-momentum components (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction,
			       DiracRep=defaultDRep);
  
  /**
   * Set 4-momentum (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(LorentzVector,const tcPDPtr &,Direction,
			       DiracRep=defaultDRep);
  
  /**
   * Set mass zero momentum (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(Energy,const tcPDPtr &,Direction,DiracRep=defaultDRep);
  
  /**
   * Set 4 momentum and mass (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(LorentzVector,Energy,const tcPDPtr &,Direction,
			       DiracRep=defaultDRep);
  
  /**
   * Default constructor (specify Dirac representation).
   */
  inline SpinorBarWaveFunction(DiracRep=defaultDRep);
  
  /**
   * Destructor.
   */
  inline ~SpinorBarWaveFunction();

  /**
   * Subscript operator for the wavefunction.
   */
  inline Complex operator ()(int ) const;

  /**
   * Set components by index.
   */
  inline Complex & operator () (int);

  /**
   * Assignment. 
   */
  inline SpinorBarWaveFunction & operator = (const SpinorBarWaveFunction &);

  /**
   * Return wavefunction as LorentzSpinor.
   */
  inline LorentzSpinorBar Wave() const;

  /**
   * Get components.
   */
  inline Complex s1() const;
  inline Complex s2() const;
  inline Complex s3() const;
  inline Complex s4() const;

  /**
   * Set components.
   */
  inline void setS1(Complex);
  inline void setS2(Complex);
  inline void setS3(Complex);
  inline void setS4(Complex);

  /**
   * Reset functions.
   */

  /**
   * Reset momentum, particle type and direction.
   */
  inline void reset(const Lorentz5Momentum &, const tcPDPtr &, Direction);

  /**
   * Reset momentum and direction.
   */
  inline void reset(const Lorentz5Momentum &,Direction);

  /** 
   * Reset momentum.
   */ 
  inline void reset(const Lorentz5Momentum &);

  /**
   * Reset the helicity (calculates the new spinor).
   * (specify Dirac representation).
   */
  inline void reset(int,DiracRep=defaultDRep);

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
   * Zero the wavefunction (specify Dirac representation).
   */
  inline void zeroWaveFunction(DiracRep=defaultDRep);
  
  /**
   * Calculate the wavefunction (specify Dirac representation).
   */
  void calculateWaveFunction(int,DiracRep=defaultDRep);
  
  /**
   * Check particle spin and set pointer.
   */
  inline void checkParticle(const tcPDPtr &);
  
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




