// -*- C++ -*-
#ifndef HERWIG_VectorWaveFunction_H
#define HERWIG_VectorWaveFunction_H
//
// This is the declaration of the <!id>VectorWaveFunction<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <id>VectorWaveFunction<!!id> class is designed to store the wavefunction
//  of a vector in a form suitable for use in helicity amplitude calculations of
//  the matrix element using a similar philosophy to the FORTRAN HELAS code.
//
//  In addition to storing the vector using the <id>LorentzPolarizationVector<!id> class
//  it inherits from the <!id>WaveFunctionBase<!!id> class to provide storage of
//  the momentum and particleData for the vector boson.
//
//  This class also contains the code which does the actually calculation of the
//  vector wavefunction.
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
// <a href="LorentzPolarizationVector.html">LorentzPolarizationVector.h</a>.
//
// Author: Peter Richardson
//

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzPolarizationVector.h>

namespace Herwig {
namespace Helicity {

using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

class VectorWaveFunction : public WaveFunctionBase {

public:
  
  // default constructors (set the momentum and Wavefunction)
  // use a 5-momentum and specify all components
  inline VectorWaveFunction(const Lorentz5Momentum &,tcPDPtr,const Complex &,
			    const Complex &,const Complex &, const Complex &);
  
  // use a 5-momentum (default phase choice)
  inline VectorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,int,int);
  // use a 5-momentum (specify phase choice)
  inline VectorWaveFunction(int, const Lorentz5Momentum &,const tcPDPtr &,int,int);

  // set all components of momentum (default phase choice)
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,int,int);
  // set all components of momentum (specify phase choice)
  inline VectorWaveFunction(int,Energy,Energy,Energy,Energy,Energy,
			    const tcPDPtr &,int,int);

  // set 4-momentum components (default phase choice)
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int,int);
  // set 4-momentum components (specify phase choice)
  inline VectorWaveFunction(int,Energy,Energy,Energy,Energy,const tcPDPtr &,int,int);

  // set 4-momentum (default phase choice)
  inline VectorWaveFunction(LorentzVector,const tcPDPtr &,int,int);
  // set 4-momentum (specify phase choice)
  inline VectorWaveFunction(int,LorentzVector,const tcPDPtr &,int,int);

  // set mass zero momentum (default phase choice)
  inline VectorWaveFunction(Energy,const tcPDPtr &,int,int);
  // set mass zero momentum (specify phase choice)
  inline VectorWaveFunction(int,Energy,const tcPDPtr &,int,int);

  // set 4 momentum and mass (default phase choice)
  inline VectorWaveFunction(LorentzVector,Energy,const tcPDPtr &,int,int);
  // set 4 momentum and mass (specify phase choice)
  inline VectorWaveFunction(int, LorentzVector,Energy,const tcPDPtr &,int,int);

  // default constructors (set the momentum and zero the Wavefunction)

  // use 5 momentum
  inline VectorWaveFunction(Lorentz5Momentum,const tcPDPtr &,int); 

  // set all components of momentum
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,int);

  // set 4-momentum components 
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int);

  // set 4-momentum 
  inline VectorWaveFunction(LorentzVector,const tcPDPtr &,int);

  // set mass zero momentum
  inline VectorWaveFunction(Energy,const tcPDPtr &,int);

  // set 4 momentum and mass
  inline VectorWaveFunction(LorentzVector,Energy,const tcPDPtr &,int);

  // default constructor
  inline VectorWaveFunction();

  // destructor
  inline ~VectorWaveFunction();

  // subscript operator for the wavefunction
  inline Complex operator ()(int ) const;

  // Set components by index.
  inline Complex & operator () (int);

  // Assignment. 
  inline VectorWaveFunction & operator = (const VectorWaveFunction &);

  // return wavefunction as polarization vector
  inline LorentzPolarizationVector Wave() const;
  
  // Get position and time.
  inline Complex x() const;
  inline Complex y() const;
  inline Complex z() const;
  inline Complex t() const;
  
  // Set position and time.
  inline void setX(const Complex&);
  inline void setY(const Complex&);
  inline void setZ(const Complex&);
  inline void setT(const Complex&);

  // reset functions

  // reset the momentum, particle type and direction
  inline void reset(const Lorentz5Momentum &, const tcPDPtr &, int);

  // reset the momentum and direction
  inline void reset(const Lorentz5Momentum &,int);

  // reset the momentum
  inline void reset(const Lorentz5Momentum &);

  // reset the helicity (recalculation the polarization vector)
  // default phase choice
  inline void reset(int);
  // specify phase choice
  inline void reset(int,int);

  // reset the particle type and direction
  inline void reset(const tcPDPtr &,int);

  // reset the particle type
  inline void reset(const tcPDPtr &);

private:

  // zero the wavefunction
  inline void zeroWaveFunction();

  // calcuate the wavefunction
  // default phase choice
  inline void calculateWaveFunction(int);
  // specify phase choice
  void calculateWaveFunction(int,int);

  // check particle spin and set pointer
  inline void checkParticle(const tcPDPtr &);

private:
  
  // storage of the wavefunction as a Lorentz Vector
  LorentzPolarizationVector _wf;
  
  // default definition of the phase 
  static const int _ivector_default=2;

};
}
}

#include "VectorWaveFunction.icc"

#endif
