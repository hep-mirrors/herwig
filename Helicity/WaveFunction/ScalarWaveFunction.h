// -*- C++ -*-
#ifndef HERWIG_ScalarWaveFunction_H
#define HERWIG_ScalarWaveFunction_H
//
// This is the declaration of the <!id>ScalarWaveFunction<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is the base class for scalar wavefunctions for use in helicity amplitude
// calculations in the Herwig++. The general approach is to use a similar philosophy 
// to the FORTRAN HELAS code but with additional structure.
//
// This class stores the scalar wavefunction as a complex number and inherits
// from the <!id>WaveFunctionBase<!!id> class for the storage of the particles
// momentum and type.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="WaveFunctionBase.html">WaveFunctionBase.h</a>.
// 
// Author: Peter Richardson
//
#include "WaveFunctionBase.h"
namespace Herwig {
namespace Helicity {
using namespace ThePEG;

class ScalarWaveFunction : public WaveFunctionBase {

public:

  // default constructors (set the momentum and Wavefunction)
  inline ScalarWaveFunction(const Lorentz5Momentum &,
			    const tcPDPtr &,Complex);

  // use a 5-momentum
  inline ScalarWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,
			    Complex,Direction);

  // set all components of momentum
  inline ScalarWaveFunction(Energy,Energy,Energy,Energy,Energy,
			    const tcPDPtr &,Complex,Direction);

  // set 4-momentum components
  inline ScalarWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,
			    Complex,Direction);

  // set 4-momentum
  inline ScalarWaveFunction(LorentzVector,const tcPDPtr &,Complex,Direction);

  // set mass zero momentum
  inline ScalarWaveFunction(Energy,const tcPDPtr &,Complex,Direction);

  // set 4 momentum and mass
  inline ScalarWaveFunction(LorentzVector,Energy,
			    const tcPDPtr &,Complex,Direction);

  // default constructors (set the momentum and zero the Wavefunction)
  // use 5 momentum
  inline ScalarWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction); 

  // set all components of momentum
  inline ScalarWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			    Direction);

  // set 4-momentum components 
  inline ScalarWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction);

  // set 4-momentum 
  inline ScalarWaveFunction(LorentzVector,const tcPDPtr &,Direction);

  // set mass zero momentum
  inline ScalarWaveFunction(Energy,const tcPDPtr &,Direction);

  // set 4 momentum and mass
  inline ScalarWaveFunction(LorentzVector,Energy,const tcPDPtr &,Direction);

  // default constructor
  inline ScalarWaveFunction();

  // destructor
  inline ~ScalarWaveFunction();

  // subscript operator for the wavefunction
  inline Complex operator ()(int ) const;

  // Set components by index.
  inline Complex & operator () (int);

  // Assignment. 
  inline ScalarWaveFunction & operator = (const ScalarWaveFunction &);

  // return the wavefunction
  inline Complex Wave() const;

  // functions to reset the wavefunction and momentum (to speed the code up)  
  // reset functions
  // reset the momentum, particle type and direction
  inline void reset(const Lorentz5Momentum &, const tcPDPtr &, Direction);

  // reset the momentum and particle type
  inline void reset(const Lorentz5Momentum &,Direction);

  // reset the momentum
  inline void reset(const Lorentz5Momentum &);

  // reset the wavefunction
  inline void reset(Complex);

  // reset the particle type and direction
  inline void reset(const tcPDPtr &,Direction);

  // reset the particle type
  inline void reset(const tcPDPtr &);
  
private:
  
  // check the particle type
  inline void checkParticle(const tcPDPtr &);

private:

  // complex number to store the wavefunction
  Complex _wf;

};
}
}

#include "ScalarWaveFunction.icc"

#endif
