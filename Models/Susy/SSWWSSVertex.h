// -*- C++ -*-
#ifndef Herwig_SSWWSSVertex_H
#define Herwig_SSWWSSVertex_H
//
// This is the declaration of the SSWWSSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VVSSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SSWWSSVertex class.
 *
 * @see \ref SSWWSSVertexInterfaces "The interfaces"
 * defined for SSWWSSVertex.
 */
class SSWWSSVertex: public VVSSVertex {

public:
  
  /**
   * The default constructor.
   */
  SSWWSSVertex();

  /**
   *  Calculate the coupling
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
			   tcPDPtr part3,tcPDPtr part4);
  

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
  SSWWSSVertex & operator=(const SSWWSSVertex &) = delete;

private:

  /**
   * Value of \f$sin(\theta_w)\f$
   */
  double sw_;

  /**
   * Value of \f$cos(\theta_w)\f$
   */
  double cw_;
  
  /**
   * Scale at which the coupling was last evaluated
   */
  Energy2 q2last_;

  /**
   * Value of coupling when last evaluated
   */
  Complex couplast_;

  /**
   * Stau mixing matrix
   */
  tMixingMatrixPtr stau_;
  
  /**
   * Stop mixing matrix
   */
  tMixingMatrixPtr stop_;

  /**
   * Sbottom mixing matrix
   */
  tMixingMatrixPtr sbottom_;
  
  /**
   * The up type sfermion present when the vertex was evaluated. 
   */
  long ulast_;

  /**
   * The down type sfermion present when the vertex was evaluated. 
   */
  long dlast_;

  /**
   * The first gauge boson present when the vertex was last evaluated. 
   */
  long gblast1_;

  /**
   * The second gauge boson present when the vertex was last evaluated. 
   */
  long gblast2_;
  
  /**
   * The value of the mixing matrix dependent part when the vertex was
   * last evaluated
   */
  Complex factlast_;
};

}

#endif /* Herwig_SSWWSSVertex_H */
