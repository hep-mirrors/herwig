// -*- C++ -*-
#ifndef Herwig_FOCUSDptoKmPipPip_H
#define Herwig_FOCUSDptoKmPipPip_H
//
// This is the declaration of the FOCUSDptoKmPipPip class.
//

#include "Herwig/Decay/FormFactors/KMatrix.h"
#include "WeakDalitzDecay.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The FOCUSDptoKmPipPip class implements the fits of FOCUS,
 * Phys.Lett. B653 (2007) 1-11 for the decay \f$D^+\to K^-\pi^+\pi^-\f$
 *
 * @see \ref FOCUSDptoKmPipPipInterfaces "The interfaces"
 * defined for FOCUSDptoKmPipPip.
 */
class FOCUSDptoKmPipPip: public WeakDalitzDecay {

public:

  /**
   * The default constructor.
   */
  FOCUSDptoKmPipPip();
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

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

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   *  Calculate the amplitude
   */
  virtual Complex amplitude(int ichan) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FOCUSDptoKmPipPip & operator=(const FOCUSDptoKmPipPip &) = delete;

private:

  /**
   *  The K-matrices for the s-wave component
   */
  //@{
  /**
   *  The \f$I=\frac12\f$ component
   */
  KMatrixPtr KHalf_;
  
  /**
   *  The \f$I=\frac32\f$ component
   */
  KMatrixPtr KThreeHalf_;
  //@}
  
  /**
   * Parameters for the f$pf$-vector
   */
  //@{
  /**
   *  Pole couplings
   */
  Energy g1_,g2_;

  /**
   * \f$\beta\f$
   */
  Energy beta_;
  
  /**
   * \f$\theta\f$
   */
  double theta_;

  /**
   *  \f$\gamma\f$
   */
  vector<double> gamma_;

  /**
   *  \f$c_{1i}\f$
   */
  vector<double> c1_;

  /**
   *  \f$c_{2i}\f$
   */
  vector<double> c2_;

  /**
   *  \f$c_{3i}\f$
   */
  vector<double> c3_;
  //@}
};

}

#endif /* Herwig_FOCUSDptoKmPipPip_H */
