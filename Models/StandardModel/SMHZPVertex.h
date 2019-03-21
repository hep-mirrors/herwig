// -*- C++ -*-
#ifndef Herwig_SMHZPVertex_H
#define Herwig_SMHZPVertex_H
//
// This is the declaration of the SMHZPVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/GeneralVVSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SMHZPVertex class.
 *
 * @see \ref SMHZPVertexInterfaces "The interfaces"
 * defined for SMHZPVertex.
 */
class SMHZPVertex: public GeneralVVSVertex {

public:
  
  /**
   * The default constructor.
   */
  SMHZPVertex();

  /**
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param part1 ParticleData pointer to first particle
   *@param part2 ParticleData pointer to second particle
   *@param part3 ParticleData pointer to third particle
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
			   tcPDPtr part3);

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
  SMHZPVertex & operator=(const SMHZPVertex &) = delete;

private:

  /**
   *  Functions for the loops from PLB 276 350
   */
  //@{
  /**
   *  The \f$I_1(\tau,\lambda)\f$ function from PLB 276 350
   */
  Complex I1(double tau,double lambda) const;

  /**
   *  The \f$I_2(\tau,\lambda)\f$ function from PLB 276 350
   */
  Complex I2(double tau,double lambda) const;

  /**
   *  The \f$f(\tau)\f$ function from  PLB 276 350
   */
  Complex f(double tau) const;

  /**
   *  The \f$g(\tau)\f$ function from  PLB 276 350
   */
  Complex g(double tau) const;
  //@}
  
private:

  /**
   *Storage of couplings
   */
  //@{
  /**
   * Last value of the coupling calculated
   */
  Complex _couplast;

  /**
   * The scale \f$q^2\f$ at which coupling was last evaluated
   */
  Energy2 _q2last;
  //@}

  /**
   * Pointer to Standard Model object
   */
  tcHwSMPtr _theSM;

  /**
   * The mass of the \f$W\f$ boson.
   */
  Energy _mw;

  /**
   * The mass of the \f$Z^0\f$ boson.
   */
  Energy _mz;

  /**
   * define quark mass scheme (fixed/running)
   */
  unsigned int _massopt;

  /**
   * The minimum flavour number in quark loops
   */
  int _minloop;

  /**
   * The maximum flavour number in quark loops
   */
  int _maxloop;
  
};

}

#endif /* Herwig_SMHZPVertex_H */
