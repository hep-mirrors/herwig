// -*- C++ -*-
#ifndef HERWIG_SpinorWaveFunction_H
#define HERWIG_SpinorWaveFunction_H
//
// This is the declaration of the SpinorWaveFunction class.

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzSpinor.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>

namespace Herwig {

using ThePEG::Helicity::LorentzSpinor;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::defaultDRep;

namespace Helicity {

using namespace ThePEG;

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  The SpinorWaveFunction class is designed to store the wavefunction
 *  of a spinor in a form suitable for use in helicity amplitude calculations 
 *  of the matrix element using a similar philosophy to the FORTRAN HELAS code.
 *
 *  In addition to storing the spinor using the LorentzSpinor class
 *  it inherits from the WaveFunctionBase class to provide storage of
 *  the momentum and particleData for the fermion.
 *
 *  This class also contains the code which does the actually calculation 
 *  of the spinor for an external particle using either of the Dirac matrix 
 *  representations currently supported in the HelicityDefinitions class.
 *
 *  When calculating the wavefunction the direction of the particle is used,
 *
 *  i.e. ipart=-1 (incoming) calculates a u spinor
 *       ipart=+1 (outgoing) calculates a v spinor
 *
 *  @see WaveFunctionBase
 *  @see LorentzSpinor
 *  @see HelicityDefinitions
 */
class SpinorWaveFunction : public WaveFunctionBase {

public:

  /**
   * Default constructors (set the momentum and Wavefunction).
   */

  /**
   * Use a 5-momentum and specify all components.
   */
  inline SpinorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,Complex,
			    Complex,Complex,Complex,DiracRep=defaultDRep);

  /**
   * Use a 5-momentum and a LorentzSpinor.
   */
  inline SpinorWaveFunction(const Lorentz5Momentum &, const tcPDPtr &,LorentzSpinor &);

  /**
   * Use a 5-momentum.
   */
  inline SpinorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,int,Direction,
			    DiracRep=defaultDRep);

  /**
   * Set all components of momentum.
   */
  inline SpinorWaveFunction(Energy,Energy,Energy,Energy,Energy,
			    const tcPDPtr &,int,Direction,DiracRep=defaultDRep);

  /**
   * Set 4-momentum components.
   */
  inline SpinorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int,
			    Direction,DiracRep=defaultDRep);

  /**
   * Set 4-momentum.
   */
  inline SpinorWaveFunction(LorentzVector,const tcPDPtr &,int,
			    Direction,DiracRep=defaultDRep);

  /**
   * Set mass zero momentum.
   */
  inline SpinorWaveFunction(Energy,const tcPDPtr &,int,Direction,
			    DiracRep=defaultDRep);

  /**
   * Set 4 momentum and mass.
   */
  inline SpinorWaveFunction(LorentzVector,Energy,const tcPDPtr &,int,Direction,
			    DiracRep=defaultDRep);

  /**
   * Default constructors (set the momentum and zero the Wavefunction).
   */

  /**
   * Use 5 momentum.
   */
  inline SpinorWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep); 

  /**
   * Set all components of momentum.
   */
  inline SpinorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			    Direction,DiracRep=defaultDRep);

  /**
   * Set 4-momentum components (default Dirac representation).
   */
  inline SpinorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep);

  /**
   * Set 4-momentum.
   */
  inline SpinorWaveFunction(LorentzVector,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep);

  /**
   * Set mass zero momentum.
   */
  inline SpinorWaveFunction(Energy,const tcPDPtr &,Direction,DiracRep=defaultDRep);

  /**
   * Set 4 momentum and mass.
   */
  inline SpinorWaveFunction(LorentzVector,Energy,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep);

  /**
   * Default constructor.
   */
  inline SpinorWaveFunction(DiracRep=defaultDRep);

  /**
   * Destructor.
   */
  inline ~SpinorWaveFunction();

  /**
   * Uubscript operator for the wavefunction.
   */
  inline Complex operator ()(int ) const;

  /**
   * Set components by index.
   */
  inline Complex & operator () (int);

  /**
   * Assignment. 
   */
  inline SpinorWaveFunction & operator = (const SpinorWaveFunction &);

  /**
   * Return wavefunction as LorentzSpinor.
   */
  inline LorentzSpinor Wave() const;

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
   * Reset momentum and particle type.
   */
  inline void reset(const Lorentz5Momentum &,Direction);

  /**
   * Reset the momentum.
   */
  inline void reset(const Lorentz5Momentum &);

  /**
   * Reset the helicity (calculates the new spinor).
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
   * Zero the wavefunction.
   */
  inline void zeroWaveFunction(DiracRep=defaultDRep);

  /**
   * Calcuate the wavefunction.
   */
  void calculateWaveFunction(int,DiracRep=defaultDRep);

  /**
   * Check particle spin and set pointer.
   */
  inline void checkParticle(const tcPDPtr &);

private:

  /**
   * Storage of the Lorentz Spinor.
   */
  LorentzSpinor _wf;

};

}

}

#include "SpinorWaveFunction.icc"

#endif




