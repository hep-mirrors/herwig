// -*- C++ -*-
#ifndef HERWIG_RSSpinorWaveFunction_H
#define HERWIG_RSSpinorWaveFunction_H
// This is the declaration of the RSSpinorWaveFunction class.

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzRSSpinor.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>

namespace Herwig {

using ThePEG::Helicity::LorentzRSSpinor;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::defaultDRep;

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
 * - ipart=-1 (incoming) calculates a u spinor
 * - ipart=+1 (outgoing) calculates a v spinor
 *
 *  The spinors are calculated using a Clebsch-Gordon decomposition in the rest-frame
 *  for a massive particle and boosted to the lab-frame. For massless particles the
 *  calculation is performed in the lab-frame (N.B. there are only two helicities +/-3/2
 *  in this case.)
 *
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
   * Use a 5-momentum and specify all components
   */
  inline RSSpinorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,
			      Complex,Complex,Complex,Complex,
			      Complex,Complex,Complex,Complex,
			      Complex,Complex,Complex,Complex,
			      Complex,Complex,Complex,Complex,DiracRep=defaultDRep);

  /**
   * Use a 5-momentum and a LorentzRSSpinor
   */
  inline RSSpinorWaveFunction(const Lorentz5Momentum &,
			      const tcPDPtr &,LorentzRSSpinor &);

  /**
   * use a 5-momentum 
   */
  inline RSSpinorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,int,Direction,
			      DiracRep=defaultDRep);

  /**
   * set all components of momentum 
   */
  inline RSSpinorWaveFunction(Energy,Energy,Energy,Energy,Energy,
			      const tcPDPtr &,int,Direction,DiracRep=defaultDRep);

  /**
   * set 4-momentum components 
   */
  inline RSSpinorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int,
			      Direction,DiracRep=defaultDRep);

  /**
   * set 4-momentum
   */
  inline RSSpinorWaveFunction(LorentzVector,const tcPDPtr &,int,
			      Direction,DiracRep=defaultDRep);
  
  /**
   * set mass zero momentum 
   */
  inline RSSpinorWaveFunction(Energy,const tcPDPtr &,int,Direction,
			      DiracRep=defaultDRep);

  /**
   * set 4 momentum and mass
   */
  inline RSSpinorWaveFunction(LorentzVector,Energy,const tcPDPtr &,int,Direction,
			      DiracRep=defaultDRep);

  /**
   * default constructors (set the momentum and zero the Wavefunction)
   * use 5 momentum 
   */
  inline RSSpinorWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep); 


  /**
   * set all components of momentum 
   */
  inline RSSpinorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			    Direction,DiracRep=defaultDRep);

  /**
   * set 4-momentum components (default Dirac representation)
   */
  inline RSSpinorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep);

  /**
   * set 4-momentum
   */
  inline RSSpinorWaveFunction(LorentzVector,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep);

  /**
   * set mass zero momentum
   */
  inline RSSpinorWaveFunction(Energy,const tcPDPtr &,Direction,DiracRep=defaultDRep);

  /**
   * set 4 momentum and mass 
   */
  inline RSSpinorWaveFunction(LorentzVector,Energy,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep);

  /**
   * default constructor
   */
  inline RSSpinorWaveFunction(DiracRep=defaultDRep);

  /**
   * destructor 
   */
  inline ~RSSpinorWaveFunction();
  //@}

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
   * Assignment. 
   */
  inline RSSpinorWaveFunction & operator = (const RSSpinorWaveFunction &);

  /**
   * return wavefunction as LorentzRSSpinor
   */
  inline LorentzRSSpinor Wave() const;

  /**
   * Get components
   */
  inline Complex xs1() const;
  /**
   * Get components
   */
  inline Complex xs2() const;
  /**
   * Get components
   */
  inline Complex xs3() const;
  /**
   * Get components
   */
  inline Complex xs4() const;
  /**
   * Get components
   */
  inline Complex ys1() const;
  /**
   * Get components
   */
  inline Complex ys2() const;
  /**
   * Get components
   */
  inline Complex ys3() const;
  /**
   * Get components
   */
  inline Complex ys4() const;
  /**
   * Get components
   */
  inline Complex zs1() const;
  /**
   * Get components
   */
  inline Complex zs2() const;
  /**
   * Get components
   */
  inline Complex zs3() const;
  /**
   * Get components
   */
  inline Complex zs4() const;
  /**
   * Get components
   */
  inline Complex ts1() const;
  /**
   * Get components
   */
  inline Complex ts2() const;
  /**
   * Get components
   */
  inline Complex ts3() const;
  /**
   * Get components
   */
  inline Complex ts4() const;

  /**
   * Set components
   */
  inline void setXS1(Complex);
  /**
   * Set components
   */
  inline void setXS2(Complex);
  /**
   * Set components
   */
  inline void setXS3(Complex);
  /**
   * Set components
   */
  inline void setXS4(Complex);
  /**
   * Set components
   */
  inline void setYS1(Complex);
  /**
   * Set components
   */
  inline void setYS2(Complex);
  /**
   * Set components
   */
  inline void setYS3(Complex);
  /**
   * Set components
   */
  inline void setYS4(Complex);
  /**
   * Set components
   */
  inline void setZS1(Complex);
  /**
   * Set components
   */
  inline void setZS2(Complex);
  /**
   * Set components
   */
  inline void setZS3(Complex);
  /**
   * Set components
   */
  inline void setZS4(Complex);
  /**
   * Set components
   */
  inline void setTS1(Complex);
  /**
   * Set components
   */
  inline void setTS2(Complex);
  /**
   * Set components
   */
  inline void setTS3(Complex);
  /**
   * Set components
   */
  inline void setTS4(Complex);

  /**
   * reset functions
   */

  /**
   * reset momentum, particle type and direction
   */
  inline void reset(const Lorentz5Momentum &, const tcPDPtr &, Direction);

  /**
   * reset momentum and particle type
   */
  inline void reset(const Lorentz5Momentum &,Direction);

  /**
   * reset the momentum
   */
  inline void reset(const Lorentz5Momentum &);

  /**
   * reset the helicity (calculates the new spinor)
   */
  inline void reset(int,DiracRep=defaultDRep);

  /**
   * reset the particle type and direction
   */
  inline void reset(const tcPDPtr &,Direction);

  /**
   * reset the particle type
   */
  inline void reset(const tcPDPtr &);

private:

  /**
   * zero the wavefunction
   */
  inline void zeroWaveFunction(DiracRep=defaultDRep);

  /**
   * calcuate the wavefunction
   */
  void calculateWaveFunction(int,DiracRep=defaultDRep);

  /**
   * check particle spin and set pointer
   */
  inline void checkParticle(const tcPDPtr &);

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
