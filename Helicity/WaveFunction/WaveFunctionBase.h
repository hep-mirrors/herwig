// -*- C++ -*-
#ifndef HERWIG_WaveFunctionBase_H
#define HERWIG_WaveFunctionBase_H
//
// This is the declaration of the <!id>WaveFunctionBase<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is the base class for all wavefunctions for use in helicity amplitude
// calculations in the Herwig++. The general approach is to use a similar philosophy 
// to the FORTRAN HELAS code but with additional structure.
//
// This class contains the storage of the particle type and 5-momentum 
// and methods to set/access this information.
//
// The methods for the wavefunction itself will we implemented in the classes
// derived from this one for the specific spin type, for example scalar, spinor,
// vector and tensor. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="ScalarWaveFunction.html">.h</a>,
// <a href="SpinorWaveFunction.html">.h</a>,
// <a href="SpinorBarWaveFunction.html">.h</a>,
// <a href="VectorWaveFunction.html">.h</a>,
// <a href="TensorWaveFunction.html">.h</a>.
// 
// Author: Peter Richardson
//

#include <ThePEG/CLHEPWrap/Lorentz5Vector.h>
#include <ThePEG/CLHEPWrap/LorentzVector.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/Config/Complex.h>

namespace Herwig {
namespace Helicity {
enum Direction {incoming,outgoing,intermediate};
using namespace ThePEG;

class WaveFunctionBase{

public:

  // default constructor
  inline WaveFunctionBase();
  // virtual destructor to keep compiler happy
  ~WaveFunctionBase();
  // subscript operator to access momentum
  inline Energy operator [](int) const;
  // get the x component of the momentum
  inline Energy px() const;
  // get the y component of the momentum
  inline Energy py() const;
  // get the z component of the momentum
  inline Energy pz() const;
  // get the energy
  inline Energy e() const;
  // get the mass
  inline Energy mass() const;
  // get offshell mass squared
  inline Energy2 m2() const;
  // Set components of the momentum
  // set the x component of the momentum
  inline void setPx(Energy);
  // set the y component of the momentum
  inline void setPy(Energy);
  // set the z component of the momentum
  inline void setPz(Energy);
  // set the energy
  inline void setE(Energy);
  // set the mass
  inline void setMass(Energy);
  // Set 5 momentum
  inline void setMomentum(const Lorentz5Momentum &);
  // set all components of momentum
  inline void setMomentum(Energy,Energy,Energy,Energy,Energy);
  // set 4-momentum components
  inline void setMomentum(Energy,Energy,Energy,Energy);
  // set 4-momentum 
  inline void setMomentum(LorentzVector);
  // set mass zero momentum
  inline void setMomentum(Energy);
  // set 4 momentum and mass
  inline void setMomentum(LorentzVector,Energy);
  // zero the 4 momentum
  inline void setMomentum();
  // get the particle id
  inline int id();
  // get 2s+1 for the particle
  inline int iSpin();
  // get the particle pointer
  inline const tcPDPtr & getParticle() const;
  // direction of particle
  inline Direction direction();
  inline void direction(Direction);
  inline const Lorentz5Momentum & getMomentum() const ;

protected:

  // set the particle pointer
  inline void setParticle(const tcPDPtr &);

private:

  // check particle type and set pointer
  void checkParticle(const tcPDPtr &);

private:

  // constant pointer to the particle info
  tcPDPtr _particle;
  // lorentz 5 momentum
  Lorentz5Momentum _momentum;
  // incoming or outgoing
  Direction _dir;
};
}
}

#include "WaveFunctionBase.icc"

#endif
