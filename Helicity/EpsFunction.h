// -*- C++ -*-
#ifndef HERWIG_EpsFunction_H
#define HERWIG_EpsFunction_H
//
// This is the declaration of the EpsFunction class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

// #include "EpsFunction.fh"
// #include "EpsFunction.xh"

namespace Herwig {
namespace Helicity {

using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  This class is designed to combine 5-momenta and polarization 
 *  vectors together with the result being the product with the 
 *  eps function. The class is purely static and contains no data.
 *
 *  @see LorentzPolarizationVector
 *  @see Lorentz5Vector
 */
class EpsFunction {
  
public:
  
  /** 
   * Various functions to return the contraction of the vectors with 
   * the epsilon function.
   */ 
  static inline LorentzPolarizationVector product(const Lorentz5Momentum &,
						  const Lorentz5Momentum &,
						  const Lorentz5Momentum &);

  static inline LorentzPolarizationVector product(const LorentzPolarizationVector &,
						  const Lorentz5Momentum &,
						  const Lorentz5Momentum &);

  static inline LorentzPolarizationVector product(const Lorentz5Momentum &,
						  const LorentzPolarizationVector &,
						  const Lorentz5Momentum &);

  static inline LorentzPolarizationVector product(const Lorentz5Momentum &,
						  const Lorentz5Momentum &,
						  const LorentzPolarizationVector &);

  static inline LorentzPolarizationVector product(const LorentzPolarizationVector &,
						  const LorentzPolarizationVector &,
						  const Lorentz5Momentum &);

  static inline LorentzPolarizationVector product(const Lorentz5Momentum &,
						  const LorentzPolarizationVector &,
						  const LorentzPolarizationVector &);

  static inline LorentzPolarizationVector product(const LorentzPolarizationVector &,
						  const Lorentz5Momentum &,
						  const LorentzPolarizationVector &);

  static inline LorentzPolarizationVector product(const LorentzPolarizationVector &,
						  const LorentzPolarizationVector &,
						  const LorentzPolarizationVector &);

private:
  
  EpsFunction();
  EpsFunction(const EpsFunction & x);
  EpsFunction & operator=(const EpsFunction & x);
  
};

}
}


#include "EpsFunction.icc"

#endif /* HERWIG_EpsFunction_H */
