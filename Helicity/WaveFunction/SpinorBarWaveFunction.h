// -*- C++ -*-
#ifndef HERWIG_SpinorBarWaveFunction_H
#define HERWIG_SpinorBarWaveFunction_H
//
// This is the declaration of the <!id>SpinorBarWaveFunction<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <id>SpinorBarWaveFunction<!!id> class is designed to store the wavefunction
//  of a barred spinor in a form suitable for use in helicity amplitude calculations of
//  the matrix element using a similar philosophy to the FORTRAN HELAS code.
//
//  In addition to storing the spinor using the <id>LorentzSpinorBar<!id> class
//  it inherits from the <!id>WaveFunctionBase<!!id> class to provide storage of
//  the momentum and particleData for the fermion.
//
//  This class also contains the code which does the actually calculation of the barred
//  spinor for an external particle using either of the Dirac matrix representations
//  currently supported in the <!id>HelicityDefinitions<!!id> class.
//
//  When calculating the wavefunction the direction of the particle is used,
//
//  i.e. ipart=-1 (incoming) calculates a v-bar spinor
//       ipart=+1 (outgoing) calculates a u-bar spinor
//
// CLASSDOC SUBSECTION See also:
//
// <a href="WaveFunctionBase.html">WaveFunctionBase.h</a>,
// <a href="LorentzSpinorBar.html">LorentzSpinor.h</a>,
// <a href="HelicityDefinitions.html">HelicityDefinitions.h</a>.
//
// Author: Peter Richardson
//
#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzSpinorBar.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>

namespace Herwig {

using ThePEG::Helicity::LorentzSpinorBar;
using ThePEG::Helicity::HelicityDefinitions;
using ThePEG::Helicity::DiracRep;

namespace Helicity {

using namespace ThePEG;

class SpinorBarWaveFunction : public WaveFunctionBase {

public:

  // default constructors (set the momentum and Wavefunction)
  
  // use a 5-momentum and specify all components (default Dirac representation)
  inline SpinorBarWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,Complex,
			       Complex,Complex,Complex);
  // use a 5-momentum and specify all components (specify Dirac representation)
  inline SpinorBarWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,Complex,
			       Complex,Complex,Complex,DiracRep);
  
  // use a 5-momentum and a LorentzSpinorBar 
  inline SpinorBarWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,
			       LorentzSpinorBar &);
  
  // use a 5-momentum (default Dirac representation)
  inline SpinorBarWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,int,Direction);
  // use a 5-momentum (specify Dirac representation)
  inline SpinorBarWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,
			       int,Direction,DiracRep);
  
  // set all components of momentum (default Dirac representation)
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,Energy,
			       const tcPDPtr &,int,Direction);
  // set all components of momentum (specify Dirac representation)
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,Energy,
			       const tcPDPtr &,int,Direction,DiracRep);
  
  // set 4-momentum components (default Dirac representation)
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int,
			       Direction);
  // set 4-momentum components (specify Dirac representation)
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,
			       int,Direction,DiracRep);
  
  // set 4-momentum (default Dirac representation)
  inline SpinorBarWaveFunction(LorentzVector,const tcPDPtr &,int,Direction);
  // set 4-momentum (specify Dirac representation)
  inline SpinorBarWaveFunction(LorentzVector,const tcPDPtr &,int,Direction,DiracRep);
  
  // set mass zero momentum (default Dirac representation)
  inline SpinorBarWaveFunction(Energy,const tcPDPtr &,int,Direction);
  // set mass zero momentum (specify Dirac representation)
  inline SpinorBarWaveFunction(Energy,const tcPDPtr &,int,Direction,DiracRep);
  
  // set 4 momentum and mass (default Dirac representation)
  inline SpinorBarWaveFunction(LorentzVector,Energy,const tcPDPtr &,int,Direction);
  // set 4 momentum and mass (specify Dirac representation)
  inline SpinorBarWaveFunction(LorentzVector,Energy,const tcPDPtr &,int,Direction,
			       DiracRep);
  
  // default constructors (set the momentum and zero the Wavefunction)
  
  // use 5 momentum (default Dirac representation)
  inline SpinorBarWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction); 
  // use 5 momentum (specify Dirac representation)
  inline SpinorBarWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction,DiracRep); 
  
  // set all components of momentum (default Dirac representation)
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,Energy,
			       const tcPDPtr &,Direction);
  // set all components of momentum (specify Dirac representation)
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,Energy,
			       const tcPDPtr &,Direction,DiracRep);
  
  // set 4-momentum components (default Dirac representation)
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction);
  // set 4-momentum components (specify Dirac representation)
  inline SpinorBarWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction,
			       DiracRep);
  
  // set 4-momentum (default Dirac representation)
  inline SpinorBarWaveFunction(LorentzVector,const tcPDPtr &,Direction);
  // set 4-momentum (specify Dirac representation)
  inline SpinorBarWaveFunction(LorentzVector,const tcPDPtr &,Direction,DiracRep);
  
  // set mass zero momentum (default Dirac representation)
  inline SpinorBarWaveFunction(Energy,const tcPDPtr &,Direction);
  // set mass zero momentum (specify Dirac representation)
  inline SpinorBarWaveFunction(Energy,const tcPDPtr &,Direction,DiracRep);
  
  // set 4 momentum and mass (default Dirac representation)
  inline SpinorBarWaveFunction(LorentzVector,Energy,const tcPDPtr &,Direction);
  // set 4 momentum and mass (specify Dirac representation)
  inline SpinorBarWaveFunction(LorentzVector,Energy,const tcPDPtr &,Direction,DiracRep);
  
  // default constructor (default Dirac representation)
  inline SpinorBarWaveFunction();
  // default constructor (specify Dirac representation)
  inline SpinorBarWaveFunction(DiracRep);
  
  // destructor
  inline ~SpinorBarWaveFunction();

  // subscript operator for the wavefunction
  inline Complex operator ()(int ) const;

  // Set components by index.
  inline Complex & operator () (int);

  // Assignment. 
  inline SpinorBarWaveFunction & operator = (const SpinorBarWaveFunction &);

  // return wavefunction as LorentzSpinor
  inline LorentzSpinorBar Wave() const;

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

  // reset momentum and direction
  inline void reset(const Lorentz5Momentum &,Direction);

  // reset momentum 
  inline void reset(const Lorentz5Momentum &);

  // reset helicity (recalculates the spinor)
  inline void reset(int );
  // reset the helicity (calculates the new spinor) (specify Dirac representation)
  inline void reset(int,DiracRep);


  // reset the particle type and direction
  inline void reset(const tcPDPtr &,Direction);

  // reset the particle type
  inline void reset(const tcPDPtr &);

private:

  // zero the wavefunction (default Dirac representation)
  inline void zeroWaveFunction();
  // zero the wavefunction (specify Dirac representation)
  inline void zeroWaveFunction(DiracRep);
  
  // calculate the wavefunction (default Dirac representation)
  inline void calculateWaveFunction(int);
  // calculate the wavefunction (specify Dirac representation)
  void calculateWaveFunction(int,DiracRep);
  
  // check particle spin and set pointer
  inline void checkParticle(const tcPDPtr &);
  
private:
  
  // storage of the Lorentz SpinorBar wavefunction
  LorentzSpinorBar _wf;
  
};
}
}

#include "SpinorBarWaveFunction.icc"

#endif




