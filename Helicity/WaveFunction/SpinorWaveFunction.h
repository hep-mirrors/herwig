// -*- C++ -*-
#ifndef HERWIG_SpinorWaveFunction_H
#define HERWIG_SpinorWaveFunction_H
//
// This is the declaration of the <!id>SpinorWaveFunction<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <id>SpinorWaveFunction<!!id> class is designed to store the wavefunction
//  of a spinor in a form suitable for use in helicity amplitude calculations of
//  the matrix element using a similar philosophy to the FORTRAN HELAS code.
//
//  In addition to storing the spinor using the <id>LorentzSpinor<!id> class
//  it inherits from the <!id>WaveFunctionBase<!!id> class to provide storage of
//  the momentum and particleData for the fermion.
//
//  This class also contains the code which does the actually calculation of the
//  spinor for an external particle using either of the Dirac matrix representations
//  currently supported in the <!id>HelicityDefinitions<!!id> class.
//
//  When calculating the wavefunction the direction of the particle is used,
//
//  i.e. ipart=-1 (incoming) calculates a u spinor
//       ipart=+1 (outgoing) calculates a v spinor
//
// CLASSDOC SUBSECTION See also:
//
// <a href="WaveFunctionBase.html">WaveFunctionBase.h</a>,
// <a href="LorentzSpinor.html">LorentzSpinor.h</a>,
// <a href="HelicityDefinitions.html">HelicityDefinitions.h</a>.
//
// Author: Peter Richardson
//
#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzSpinor.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>

namespace Herwig {

using ThePEG::Helicity::LorentzSpinor;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::defaultDRep;

namespace Helicity {

using namespace ThePEG;

class SpinorWaveFunction : public WaveFunctionBase {

public:

  // default constructors (set the momentum and Wavefunction)

  // use a 5-momentum and specify all components
  inline SpinorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,Complex,
			    Complex,Complex,Complex,DiracRep=defaultDRep);

  // use a 5-momentum and a LorentzSpinor
  inline SpinorWaveFunction(const Lorentz5Momentum &, const tcPDPtr &,LorentzSpinor &);

  // use a 5-momentum 
  inline SpinorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,int,Direction,
			    DiracRep=defaultDRep);

  // set all components of momentum 
  inline SpinorWaveFunction(Energy,Energy,Energy,Energy,Energy,
			    const tcPDPtr &,int,Direction,DiracRep=defaultDRep);

 // set 4-momentum components 
  inline SpinorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int,
			    Direction,DiracRep=defaultDRep);

  // set 4-momentum
  inline SpinorWaveFunction(LorentzVector,const tcPDPtr &,int,
			    Direction,DiracRep=defaultDRep);

  // set mass zero momentum 
  inline SpinorWaveFunction(Energy,const tcPDPtr &,int,Direction,
			    DiracRep=defaultDRep);

  // set 4 momentum and mass
  inline SpinorWaveFunction(LorentzVector,Energy,const tcPDPtr &,int,Direction,
			    DiracRep=defaultDRep);

  // default constructors (set the momentum and zero the Wavefunction)

  // use 5 momentum 
  inline SpinorWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep); 

  // set all components of momentum 
  inline SpinorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			    Direction,DiracRep=defaultDRep);

  // set 4-momentum components (default Dirac representation)
  inline SpinorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep);

  // set 4-momentum
  inline SpinorWaveFunction(LorentzVector,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep);

  // set mass zero momentum
  inline SpinorWaveFunction(Energy,const tcPDPtr &,Direction,DiracRep=defaultDRep);

  // set 4 momentum and mass 
  inline SpinorWaveFunction(LorentzVector,Energy,const tcPDPtr &,Direction,
			    DiracRep=defaultDRep);

  // default constructor
  inline SpinorWaveFunction(DiracRep=defaultDRep);

  // destructor 
  inline ~SpinorWaveFunction();

  // subscript operator for the wavefunction
  inline Complex operator ()(int ) const;
  // Set components by index.
  inline Complex & operator () (int);

  // Assignment. 
  inline SpinorWaveFunction & operator = (const SpinorWaveFunction &);

  // return wavefunction as LorentzSpinor
  inline LorentzSpinor Wave() const;

  // Get components
  inline Complex s1() const;
  inline Complex s2() const;
  inline Complex s3() const;
  inline Complex s4() const;

  // Set components
  inline void setS1(Complex);
  inline void setS2(Complex);
  inline void setS3(Complex);
  inline void setS4(Complex);

  // reset functions

  // reset momentum, particle type and direction
  inline void reset(const Lorentz5Momentum &, const tcPDPtr &, Direction);

  // reset momentum and particle type
  inline void reset(const Lorentz5Momentum &,Direction);

  // reset the momentum
  inline void reset(const Lorentz5Momentum &);

  // reset the helicity (calculates the new spinor)
  inline void reset(int,DiracRep=defaultDRep);

  // reset the particle type and direction
  inline void reset(const tcPDPtr &,Direction);

  // reset the particle type
  inline void reset(const tcPDPtr &);

private:

  // zero the wavefunction
  inline void zeroWaveFunction(DiracRep=defaultDRep);



  // calcuate the wavefunction
  void calculateWaveFunction(int,DiracRep=defaultDRep);

  // check particle spin and set pointer
  inline void checkParticle(const tcPDPtr &);

private:

  // storage of the Lorentz Spinor
  LorentzSpinor _wf;
};

}
}

#include "SpinorWaveFunction.icc"

#endif




