// -*- C++ -*-
#ifndef HERWIG_ThreeBodyOnShellME_H
#define HERWIG_ThreeBodyOnShellME_H
// This is the declaration of the ThreeBodyOnShellME class.

#include <ThePEG/PDT/Decayer.h>
#include "CLHEP/GenericFunctions/AbsFunction.hh"
#include "Herwig++/Decay/DecayIntegrator.fh"

namespace Herwig {
using namespace Genfun;
using namespace ThePEG; 

/** \ingroup PDT
 *
 *  This is a simple class to allow the integration of the threeBodyMatrixElement
 *  member of the DecayIntegrator class using the ThreeBodyIntegrator class
 *
 * @see DecayIntegrator
 * @see ThreeBodyAllOnCalculator
 *
 */
class ThreeBodyOnShellME: public Genfun::AbsFunction {


  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;
  
  /**
   * Clone method
   */
   ThreeBodyOnShellME *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;
  
public:
  
  /**
   * Constructor
   * @param decayer Pointer to the decayer.
   * @param mode The mode we are integrating
   */
  ThreeBodyOnShellME(DecayIntegratorPtr decayer, int mode);
  
  /**
   * Destructor
   */
  virtual ~ThreeBodyOnShellME();

  /**
   * The number of variables this is a function of (i.e. 7)
   */
  virtual unsigned int dimensionality() const ;     
  
  /**
   * Copy constructor
   */
  ThreeBodyOnShellME(const ThreeBodyOnShellME &right);
  
  /**
   * Retreive function value
   */
  virtual double operator ()(double) const {return 0.;}

  /**
   * Retreive function value
   */
  virtual double operator ()(const Argument & a) const ;
  
private:
  
  /**
   * It is illegal to assign a function
   */
  const ThreeBodyOnShellME & 
  operator=(const ThreeBodyOnShellME &right);
  
private:

  /**
   * Pointer to the decayer.
   */
  DecayIntegratorPtr _decayer;

  /**
   * The mode we are integrating
   */
  int _mode;
};

}

#include "ThreeBodyOnShellME.icc"

#endif /* HERWIG_ThreeBodyOnShellME_H */
