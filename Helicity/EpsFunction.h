// -*- C++ -*-
#ifndef HERWIG_EpsFunction_H
#define HERWIG_EpsFunction_H
//
// This is the declaration of the EpsFunction class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/LorentzTensor.h"

namespace Herwig {
namespace Helicity {

using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::LorentzTensor;

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
   *  Return the product 
   *  \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   * @param v1 The first  vector \f$v_{1\alpha}\f$.
   * @param v2 The second vector \f$v_{2\alpha}\f$.
   * @param v3 The third  vector \f$v_{3\alpha}\f$.
   * @return The product 
   * \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   */
  static inline LorentzPolarizationVector product(const Lorentz5Momentum & v1,
						  const Lorentz5Momentum & v2,
						  const Lorentz5Momentum & v3);

  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   * @param v1 The first  vector \f$v_{1\alpha}\f$.
   * @param v2 The second vector \f$v_{2\alpha}\f$.
   * @param v3 The third  vector \f$v_{3\alpha}\f$.
   * @return The product 
   * \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   */
  static inline LorentzPolarizationVector product(const LorentzPolarizationVector & v1,
						  const Lorentz5Momentum & v2,
						  const Lorentz5Momentum & v3);

  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   * @param v1 The first  vector \f$v_{1\alpha}\f$.
   * @param v2 The second vector \f$v_{2\alpha}\f$.
   * @param v3 The third  vector \f$v_{3\alpha}\f$.
   * @return The product 
   * \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   */
  static inline LorentzPolarizationVector product(const Lorentz5Momentum & v1,
						  const LorentzPolarizationVector & v2,
						  const Lorentz5Momentum & v3);

  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   * @param v1 The first  vector \f$v_{1\alpha}\f$.
   * @param v2 The second vector \f$v_{2\alpha}\f$.
   * @param v3 The third  vector \f$v_{3\alpha}\f$.
   * @return The product 
   * \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   */
  static inline LorentzPolarizationVector product(const Lorentz5Momentum & v1,
						  const Lorentz5Momentum & v2,
						  const LorentzPolarizationVector & v3);

  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   * @param v1 The first  vector \f$v_{1\alpha}\f$.
   * @param v2 The second vector \f$v_{2\alpha}\f$.
   * @param v3 The third  vector \f$v_{3\alpha}\f$.
   * @return The product 
   * \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   */
  static inline LorentzPolarizationVector product(const LorentzPolarizationVector & v1,
						  const LorentzPolarizationVector & v2,
						  const Lorentz5Momentum & v3);

  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   * @param v1 The first  vector \f$v_{1\alpha}\f$.
   * @param v2 The second vector \f$v_{2\alpha}\f$.
   * @param v3 The third  vector \f$v_{3\alpha}\f$.
   * @return The product 
   * \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   */
  static inline LorentzPolarizationVector product(const Lorentz5Momentum & v1,
						  const LorentzPolarizationVector & v2,
						  const LorentzPolarizationVector & v3);

  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   * @param v1 The first  vector \f$v_{1\alpha}\f$.
   * @param v2 The second vector \f$v_{2\alpha}\f$.
   * @param v3 The third  vector \f$v_{3\alpha}\f$.
   * @return The product 
   * \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   */
  static inline LorentzPolarizationVector product(const LorentzPolarizationVector & v1,
						  const Lorentz5Momentum & v2,
						  const LorentzPolarizationVector & v3);

  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   * @param v1 The first  vector \f$v_{1\alpha}\f$.
   * @param v2 The second vector \f$v_{2\alpha}\f$.
   * @param v3 The third  vector \f$v_{3\alpha}\f$.
   * @return The product 
   * \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   */
  static inline LorentzPolarizationVector product(const LorentzPolarizationVector & v1,
						  const LorentzPolarizationVector & v2,
						  const LorentzPolarizationVector & v3);

  /**
   *  Return the product \f$\epsilon^{\mu\nu\alpha\beta}v_{1\alpha}v_{2\beta}\f$.
   */
  static inline LorentzTensor product(const LorentzPolarizationVector & v1,
				      const LorentzPolarizationVector & v2);

private:

  /**
   * This is a static class so the constructor is private.
   */  
  EpsFunction();

  /**
   * This is a static class so the copy-constructor is private. 
   */
  EpsFunction(const EpsFunction & x);

  /**
   * This is a statoc class so the assignment operator is private.
   */
  EpsFunction & operator=(const EpsFunction & x);
  
};

}
}


#include "EpsFunction.icc"

#endif /* HERWIG_EpsFunction_H */
