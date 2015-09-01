// -*- C++ -*-
#ifndef Herwig_RPVWSSVertex_H
#define Herwig_RPVWSSVertex_H
//
// This is the declaration of the RPVWSSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "Herwig/Models/Susy/MixingMatrix.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the RPVWSSVertex class.
 *
 * @see \ref RPVWSSVertexInterfaces "The interfaces"
 * defined for RPVWSSVertex.
 */
class RPVWSSVertex: public Helicity::VSSVertex {

public:

  /**
   * The default constructor.
   */
  RPVWSSVertex();

 /**
  * Calculate the couplings.
  * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
  * @param part1 The ParticleData pointer for the first  particle.
  * @param part2 The ParticleData pointer for the second particle.
  * @param part3 The ParticleData pointer for the third  particle.
  */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
                           tcPDPtr part2,tcPDPtr part3);

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
  RPVWSSVertex & operator=(const RPVWSSVertex &);

private:

  /**
   * Value of \f$sin(\theta_w)\f$
   */
  double _sw;

  /**
   * Value of \f$cos(\theta_w)\f$
   */
  double _cw;

  /**
   * Stau mixing matrix
   */
  tMixingMatrixPtr _stau;
  
  /**
   * Stop mixing matrix
   */
  tMixingMatrixPtr _stop;
  
  /**
   * Sbottom mixing matrix
   */
  tMixingMatrixPtr _sbottom;

  /**
   * The value of \f$\sin2\theta_W\f$ 
   */
  double _s2w;

  /**
   * The value of \f$\cos2\theta_W\f$ 
   */
  double _c2w;

  /**
   *  Coupling of Z to scalar and pseudoscalar
   */
  vector<vector<Complex> > Cijeo_;

  /**
   *  Coupling of W to scalar charged
   */
  vector<vector<Complex> > Cijec_;

  /**
   *  Coupling of W to pseudoscalar charged
   */
  vector<vector<Complex> > Cijco_;
  
  /**
   *  Coupling of Z to charged Higgs
   */
  vector<vector<Complex> > Cijc_;

  /**
   *  Which interactions to include 
   */
  unsigned int _interactions;

  /**
   * Scale at which the coupling was last evaluated
   */
  Energy2 _q2last;
  
   /**
    * The up type sfermion present when the vertex was evaluated. 
    */
   long _ulast;
 
   /**
    * The down type sfermion present when the vertex was evaluated. 
    */
   long _dlast;
 
   /**
    * The gauge boson present when the vertex was last evaluated. 
    */
   long _gblast;
  
  /**
   * The value of the mixing matrix dependent part when the vertex was
   * last evaluated
   */
  Complex _factlast;

  /**
   * Value of coupling when last evaluated
   */
  Complex _couplast;
};

}

#endif /* Herwig_RPVWSSVertex_H */
