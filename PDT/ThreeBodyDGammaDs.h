// -*- C++ -*-
#ifndef HERWIG_ThreeBodyDGammaDs_H
#define HERWIG_ThreeBodyDGammaDs_H
// This is the declaration of the ThreeBodyDGammaDs class.

#include "CLHEP/GenericFunctions/AbsFunction.hh"

namespace Herwig {
using namespace ThePEG;
using namespace Genfun;

/**
 * This class encpasulates the dGammads member of a DecayIntegrator as
 * a function to allow it to be integrated to give the width.
 *
 * @see DecayIntegrator
 * 
 */
class ThreeBodyDGammaDs: public Genfun::AbsFunction {

  FUNCTION_OBJECT_DEF(ThreeBodyDGammaDs)

public:

  /**
   * Constructor
   */
  ThreeBodyDGammaDs(DecayIntegratorPtr,int);
  
  /**
   * Destructor
   */
  virtual ~ThreeBodyDGammaDs();
  
  virtual unsigned int dimensionality() const ;     
  
  /**
   * Copy constructor
   */
  ThreeBodyDGammaDs(const ThreeBodyDGammaDs &right);
  
  /**
   * Retreive function value
   */
  virtual double operator ()(double) const {return 0.;}
  virtual double operator ()(const Argument & a) const ;
  
private:
  
  /**
   * It is illegal to assign a function
   */
  const ThreeBodyDGammaDs & 
  operator=(const ThreeBodyDGammaDs &right);

private:

  DecayIntegratorPtr _decayer;
  int _mode;

};

}

#include "ThreeBodyDGammaDs.icc"

#endif /* HERWIG_ThreeBodyDGammaDs_H */
