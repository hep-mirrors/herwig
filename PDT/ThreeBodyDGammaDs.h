// -*- C++ -*-
#ifndef HERWIG_ThreeBodyDGammaDs_H
#define HERWIG_ThreeBodyDGammaDs_H
// This is the declaration of the ThreeBodyDGammaDs class.

#include "CLHEP/GenericFunctions/AbsFunction.hh"

namespace Herwig {
using namespace ThePEG;
using namespace Genfun;

/** \ingroup PDT
 *
 * This class encpasulates the dGammads member of a DecayIntegrator as
 * a function to allow it to be integrated to give the width.
 *
 * @see DecayIntegrator
 * 
 */
class ThreeBodyDGammaDs: public Genfun::AbsFunction {

  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;

  /**
   * Clone method
   */
   ThreeBodyDGammaDs *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;
  
public:

  /**
   * Constructor
   * @param decayer A pointer to the DecayIntegrator class
   * @param mode The mode in the decayer being integrated.
   */
  ThreeBodyDGammaDs(DecayIntegratorPtr decayer,int mode);

  /**
   * Copy constructor
   */
  ThreeBodyDGammaDs(const ThreeBodyDGammaDs & x) 
    : Genfun::AbsFunction(), _decayer(x._decayer), _mode(x._mode) {}
  // this one's required because CLHEP::Genfun::AbsFunction 
  // doesn't implement its copy constructor. Aaaargh!

  /**
   *  The number of variables for the function (in this case 5)
   */  
  virtual unsigned int dimensionality() const ;     
  
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
  const ThreeBodyDGammaDs & 
  operator=(const ThreeBodyDGammaDs &right);

private:

  /**
   * Pointer to the decayer
   */
  DecayIntegratorPtr _decayer;

  /**
   *  The mode from the decayer.
   */
  int _mode;

};

}

#include "ThreeBodyDGammaDs.icc"

#endif /* HERWIG_ThreeBodyDGammaDs_H */
