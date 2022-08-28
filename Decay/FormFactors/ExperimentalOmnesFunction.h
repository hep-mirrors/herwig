// -*- C++ -*-
#ifndef Herwig_ExperimentalOmnesFunction_H
#define Herwig_ExperimentalOmnesFunction_H
//
// This is the declaration of the ExperimentalOmnesFunction class.
//

#include "OmnesFunction.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the ExperimentalOmnesFunction class.
 *
 * @see \ref ExperimentalOmnesFunctionInterfaces "The interfaces"
 * defined for ExperimentalOmnesFunction.
 */
class ExperimentalOmnesFunction: public OmnesFunction {

public:
  
  /**
   * The default constructor.
   */
  ExperimentalOmnesFunction();

  /**
   *  Method to return the function value
   */
  virtual Complex D(Energy2 s) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

public:

  /**
   *  Integrand for the Omnes function
   */
  InvEnergy4 operator ()(Energy2 xpoint) const {
    InvEnergy4 output = InvEnergy4();
    Energy q(sqrt(xpoint));
    if(abs(xpoint-s_)>sqr(epsCut_)) 
      output= (*interpolator_)(q)/xpoint/(xpoint-s_);
    return output;
  }
  /** Return type for the GaussianIntegrator */
  typedef InvEnergy4 ValType;
  /** Argument type for the GaussianIntegrator */
  typedef Energy2 ArgType;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}
  
protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ExperimentalOmnesFunction & operator=(const ExperimentalOmnesFunction &) = delete;

private:

  /**
   *  Energy values for the experimental data on the phase shift
   */
  vector<Energy> energy_;

  /**
   *  Experimental values of the phase shift
   */
  vector<double> phase_;

  /**
   *  Energy values for the interpolation table for the Omnes function.
   */
  vector<Energy> omnesEnergy_;

  /**
   * Real part of the Omnes function for the interpolation table
   */
  vector<double> omnesFunctionRe_;

  /**
   * Imaginary part of the Omnes function for the interpolation table
   */
  vector<double> omnesFunctionIm_;

  /**
   * set up of the interpolation table
   */
  bool initialize_;

  /**
   * Number of points for the intepolation of the experimental Omnes function
   */
  unsigned int nPoints_;

  /**
   * Interpolators for the experimental Omnes function.
   */
  //@{
  /**
   *  The interpolator for the real part
   */
  mutable Interpolator<double,Energy>::Ptr oRe_;

  /**
   *  The interpolator for the imaginary part
   */
  mutable Interpolator<double,Energy>::Ptr oIm_;
  //@}

  /**
   *  Cut-off parameter for the integral of the experimental function
   */ 
  Energy epsCut_;

  /**
   *  Size of the vectors for the experimental data 
   */
  unsigned int nsizea_;

  /**
   * Size of the vectors for the interpolation tables
   */
  unsigned int nsizeb_;

  /**
   *  Interpolator and scale for the Omes function integral
   */
  //@
  /**
   *  The interpolator
   */
  Interpolator<double,Energy>::Ptr interpolator_;

  /**
   *  The scale
   */
  Energy2 s_;
  //@}

};

}

#endif /* Herwig_ExperimentalOmnesFunction_H */
