// -*- C++ -*-
#ifndef HERWIG_TensorWaveFunction_H
#define HERWIG_TensorWaveFunction_H
//
// This is the declaration of the <!id>TensorWaveFunction<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <id>TensorWaveFunction<!!id> class is designed to store the wavefunction
//  of a tensor in a form suitable for use in helicity amplitude calculations of
//  the matrix element using a similar philosophy to the FORTRAN HELAS code.
//
//  In addition to storing the tensor using the <id>LorentzTensor<!id> class
//  it inherits from the <!id>WaveFunctionBase<!!id> class to provide storage of
//  the momentum and particleData for the vector boson.
//
//  This class also contains the code which does the actually calculation of the
//  tensor wavefunction.
//
//  There are two choices available for the calculation of the wavefunction.
//  These are set using the TensorPhase enumeration which specifies a default choice.
//  The first choice, tensor_phase, includes a phase factor exp(+/- i phi) for the 
//  +/- helicity states while the second, tensor_nophase, does not.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="WaveFunctionBase.html">WaveFunctionBase.h</a>,
// <a href="LorentzTensor.html">LorentzTensor.h</a>,
// <a href="VectorWaveFunction.html">VectorWaveFunction.h</a>.
//
// Author: Peter Richardson
//

#include "WaveFunctionBase.h"
#include "VectorWaveFunction.h"
#include <ThePEG/Helicity/LorentzTensor.h>

namespace Herwig {
namespace Helicity {
using ThePEG::Helicity::LorentzTensor;

enum TensorPhase {tensor_phase, tensor_nophase, default_tensor_phase=tensor_nophase};

class TensorWaveFunction : public WaveFunctionBase {

public:

  // default constructors (set the momentum and Wavefunction)
  inline TensorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,
			    Complex,Complex,Complex,Complex,Complex,Complex,
			    Complex,Complex,Complex,Complex,Complex,Complex,
			    Complex,Complex,Complex,Complex);

  // use a 5-momentum (specify phase choice)
  inline TensorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,int,Direction,
			    TensorPhase=default_tensor_phase);


  // set all components of momentum (specify phase choice)
  inline TensorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			    int,Direction,TensorPhase=default_tensor_phase);

  // set 4-momentum components (specify phase choice)
  inline TensorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int,Direction,
			    TensorPhase=default_tensor_phase);

  // set 4-momentum (specify phase choice) 
  inline TensorWaveFunction(LorentzMomentum,const tcPDPtr &,int,Direction,
			    TensorPhase=default_tensor_phase);

  // set mass zero momentum (specify phase choice)
  inline TensorWaveFunction(Energy,const tcPDPtr &,int,Direction,
			    TensorPhase=default_tensor_phase);

  // set 4 momentum and mass (specify phase choice)
  inline TensorWaveFunction(LorentzMomentum,Energy,const tcPDPtr &,int,Direction,
			    TensorPhase=default_tensor_phase);

  // default constructors (set the momentum and zero the Wavefunction)

  // use 5 momentum
  inline TensorWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction); 

  // set all components of momentum
  inline TensorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			   Direction);

  // set 4-momentum components 
  inline TensorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction);

  // set 4-momentum 
  inline TensorWaveFunction(LorentzMomentum,const tcPDPtr &,Direction);

  // set mass zero momentum
  inline TensorWaveFunction(Energy,const tcPDPtr &,Direction);

  // set 4 momentum and mass
  inline TensorWaveFunction(LorentzMomentum,Energy,const tcPDPtr &,Direction);

  // default constructor
  inline TensorWaveFunction();

  // destructor
  inline ~TensorWaveFunction();

  // subscript operator for the wavefunction
  inline Complex operator ()(int,int ) const;

  // Set components by index.
  inline Complex & operator () (int,int);

  // Assignment. 
  inline TensorWaveFunction & operator = (const TensorWaveFunction &);

  // return wavefunction as polarization vector
  inline LorentzTensor Wave() const;

  // Get components
  inline Complex xx() const;
  inline Complex yx() const;
  inline Complex zx() const;
  inline Complex tx() const;
  inline Complex xy() const;
  inline Complex yy() const;
  inline Complex zy() const;
  inline Complex ty() const;
  inline Complex xz() const;
  inline Complex yz() const;
  inline Complex zz() const;
  inline Complex tz() const;
  inline Complex xt() const;
  inline Complex yt() const;
  inline Complex zt() const;
  inline Complex tt() const;

  // Set Components
  inline void setXX(Complex);
  inline void setYX(Complex);
  inline void setZX(Complex);
  inline void setTX(Complex);
  inline void setXY(Complex);
  inline void setYY(Complex);
  inline void setZY(Complex);
  inline void setTY(Complex);
  inline void setXZ(Complex);
  inline void setYZ(Complex);
  inline void setZZ(Complex);
  inline void setTZ(Complex);
  inline void setXT(Complex);
  inline void setYT(Complex);
  inline void setZT(Complex);
  inline void setTT(Complex);

  // reset functions

  // reset momentum, particle type and direction
  inline void reset(const Lorentz5Momentum &, const tcPDPtr &, Direction);

  // reset momentum and direction
  inline void reset(const Lorentz5Momentum &,Direction);

  // reset momentum
  inline void reset(const Lorentz5Momentum &);

  // reset helicity (recalculate the tensor )
  inline void reset(int,TensorPhase=default_tensor_phase);

  // reset particle type and direction
  inline void reset(const tcPDPtr &,Direction);

  // reset particle type
  inline void reset(const tcPDPtr &);

private:

  // zero the wavefunction
  inline void zeroWaveFunction();

  // calculate the wavefunction (specify phase choice)
  void calculateWaveFunction(int,TensorPhase=default_tensor_phase);

  // check particle spin and set pointer
  inline void checkParticle(const tcPDPtr &);

private:

  // storage of the wavefunction as a Lorentz Tensor
  LorentzTensor _wf;

};
}
}

#include "TensorWaveFunction.icc"

#endif
