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
//  There are two choices available for the calculation of the wavefunction.
//  These are set using the VectorPhase enumeration which specifies a default choice.
//  The first choice, vector_phase, includes a phase factor exp(+/- i phi) for the 
//  +/- helicity states while the second, vector_nophase, does not.
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

enum VectorPhase {vector_phase, vector_nophase, default_vector_phase=vector_nophase};

using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

class VectorWaveFunction : public WaveFunctionBase {

public:
  
  // default constructors (set the momentum and Wavefunction)
  // use a 5-momentum and specify all components
  inline VectorWaveFunction(const Lorentz5Momentum &,tcPDPtr,const Complex &,
			    const Complex &,const Complex &, const Complex &);
  
  // use a 5-momentum (specify phase choice)
  inline VectorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,int,Direction,
			    VectorPhase=default_vector_phase);

  // set all components of momentum (specify phase choice)
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			    int,Direction,VectorPhase=default_vector_phase);

  // set 4-momentum components (specify phase choice)
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int,Direction,
			    VectorPhase=default_vector_phase);

  // set 4-momentum (specify phase choice)
  inline VectorWaveFunction(LorentzVector,const tcPDPtr &,int,Direction,
			    VectorPhase=default_vector_phase);

  // set mass zero momentum (specify phase choice)
  inline VectorWaveFunction(Energy,const tcPDPtr &,int,Direction,
			    VectorPhase=default_vector_phase);

  // set 4 momentum and mass (specify phase choice)
  inline VectorWaveFunction(LorentzVector,Energy,const tcPDPtr &,int,Direction,
			    VectorPhase=default_vector_phase);

  // default constructors (set the momentum and zero the Wavefunction)

  // use 5 momentum
  inline VectorWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction); 

  // set all components of momentum
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			    Direction);

  // set 4-momentum components 
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction);

  // set 4-momentum 
  inline VectorWaveFunction(LorentzVector,const tcPDPtr &,Direction);

  // set mass zero momentum
  inline VectorWaveFunction(Energy,const tcPDPtr &,Direction);

  // set 4 momentum and mass
  inline VectorWaveFunction(LorentzVector,Energy,const tcPDPtr &,Direction);

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
  inline void reset(const Lorentz5Momentum &, const tcPDPtr &, Direction);

  // reset the momentum and direction
  inline void reset(const Lorentz5Momentum &,Direction);

  // reset the momentum
  inline void reset(const Lorentz5Momentum &);

  // reset the helicity (recalculation the polarization vector)
  inline void reset(int,VectorPhase=default_vector_phase);

  // reset the particle type and direction
  inline void reset(const tcPDPtr &,Direction);

  // reset the particle type
  inline void reset(const tcPDPtr &);

private:

  // zero the wavefunction
  inline void zeroWaveFunction();

  // calcuate the wavefunction
  // default phase choice
  inline void calculateWaveFunction(int);
  // specify phase choice
  void calculateWaveFunction(int,VectorPhase=default_vector_phase);

  // check particle spin and set pointer
  inline void checkParticle(const tcPDPtr &);

private:
  
  // storage of the wavefunction as a Lorentz Vector
  LorentzPolarizationVector _wf;
  
};

}
}

#include "VectorWaveFunction.icc"

#endif
