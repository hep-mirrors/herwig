// -*- C++ -*-
#ifndef HERWIG_EpsFunction_H
#define HERWIG_EpsFunction_H
//
// This is the declaration of the <!id>EpsFunction<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class is desgined to combine 5-momenta and polarization vectors together
//  with the result being the product with the eps function. The class is purely static
//  and contains no data
//
// CLASSDOC SUBSECTION See also:
//
// <a href="LorentzPolarizationVector.html">LorentzPolarizationVector.h</a>,
// <a href="Lorentz5Vector.html">Lorentz5Vector.h</a>.
// 
//  Author: Peter Richardson
//

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

// #include "EpsFunction.fh"
// #include "EpsFunction.xh"

namespace Herwig {
namespace Helicity {

using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

class EpsFunction {
  
public:
  
  // various functions to return the contraction of the vectors with the epsilon
  // function 
  static inline LorentzPolarizationVector product(Lorentz5Momentum,Lorentz5Momentum,
						  Lorentz5Momentum);
  
private:
  
  EpsFunction();
  EpsFunction(const EpsFunction & x);
  EpsFunction & operator=(const EpsFunction & x);
  
};

}
}

// CLASSDOC OFF


#include "EpsFunction.icc"

#endif /* HERWIG_EpsFunction_H */
