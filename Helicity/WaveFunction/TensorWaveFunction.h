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
//  There are two choices available for the calculation of the wavefunction,
//  the first, ivector=1, includes a phase factor exp(+/- i phi) for the 
//  +/- helicity states while the second ivector=2 does not.
//
//  The static variable _ivector_default controls which of these is used if the user
//  does not specify a choice.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="WaveFunctionBase.html">WaveFunctionBase.h</a>,
// <a href="LorentzTensor.html">LorentzTensor.h</a>.
// <a href="VectorWaveFunction.html">VectorWaveFunction.h</a>.
//
// Author: Peter Richardson
//

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzTensor.h>

namespace Herwig {
namespace Helicity {
using ThePEG::Helicity::LorentzTensor;

class TensorWaveFunction : public WaveFunctionBase {

public:

  // default constructors (set the momentum and Wavefunction)
  inline TensorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,
			    Complex,Complex,Complex,Complex,Complex,Complex,
			    Complex,Complex,Complex,Complex,Complex,Complex,
			    Complex,Complex,Complex,Complex);

  // use a 5-momentum (default phase choice)
  inline TensorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,int,int);
  // use a 5-momentum (specify phase choice)
  inline TensorWaveFunction(int,const Lorentz5Momentum &,const tcPDPtr &,int,int);

  // set all components of momentum (default phase choice)
  inline TensorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,int,int);
  // set all components of momentum (specify phase choice)
  inline TensorWaveFunction(int,Energy,Energy,Energy,Energy,Energy,
			    const tcPDPtr &,int,int);

  // set 4-momentum components (default phase choice)
  inline TensorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int,int);
  // set 4-momentum components (specify phase choice)
  inline TensorWaveFunction(int,Energy,Energy,Energy,Energy,const tcPDPtr &,int,int);

  // set 4-momentum (default phase choice)
  inline TensorWaveFunction(LorentzVector,const tcPDPtr &,int,int);
  // set 4-momentum (specify phase choice) 
  inline TensorWaveFunction(int,LorentzVector,const tcPDPtr &,int,int);

  // set mass zero momentum (default phase choice)
  inline TensorWaveFunction(Energy,const tcPDPtr &,int,int);
  // set mass zero momentum (specify phase choice)
  inline TensorWaveFunction(int,Energy,const tcPDPtr &,int,int);

  // set 4 momentum and mass (default phase choice)
  inline TensorWaveFunction(LorentzVector,Energy,const tcPDPtr &,int,int);
  // set 4 momentum and mass (specify phase choice)
  inline TensorWaveFunction(int,LorentzVector,Energy,const tcPDPtr &,int,int);

  // default constructors (set the momentum and zero the Wavefunction)

  // use 5 momentum
  inline TensorWaveFunction(Lorentz5Momentum,const tcPDPtr &,int); 

  // set all components of momentum
  inline TensorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,int);

  // set 4-momentum components 
  inline TensorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int);

  // set 4-momentum 
  inline TensorWaveFunction(LorentzVector,const tcPDPtr &,int);

  // set mass zero momentum
  inline TensorWaveFunction(Energy,const tcPDPtr &,int);

  // set 4 momentum and mass
  inline TensorWaveFunction(LorentzVector,Energy,const tcPDPtr &,int);

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
  inline void reset(const Lorentz5Momentum &, const tcPDPtr &, int);

  // reset momentum and direction
  inline void reset(const Lorentz5Momentum &,int);

  // reset momentum
  inline void reset(const Lorentz5Momentum &);

  // reset helicity (recalculate the tensor )
  // default phase choice
  inline void reset(int);
  // specfiy phase choice
  inline void reset(int,int);

  // reset particle type and direction
  inline void reset(const tcPDPtr &,int);

  // reset particle type
  inline void reset(const tcPDPtr &);

private:

  // zero the wavefunction
  inline void zeroWaveFunction();

  // calcuate the wavefunction (default phase choice)
  inline void calculateWaveFunction(int);

  // calculate the wavefunction (specify phase choice)
  void calculateWaveFunction(int,int);

  // check particle spin and set pointer
  inline void checkParticle(const tcPDPtr &);

private:

  // storage of the wavefunction as a Lorentz Tensor
  LorentzTensor _wf;

  // default definition of the phase 
  static const int _ivector_default=2;

};
}
}

#include "TensorWaveFunction.icc"

#endif
