// -*- C++ -*-
#ifndef THEPEG_SextetPVVVertex_H
#define THEPEG_SextetPVVVertex_H
//
// This is the declaration of the SextetPVVVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SextetPVVVertex class.
 *
 * @see \ref SextetPVVVertexInterfaces "The interfaces"
 * defined for SextetPVVVertex.
 */
  class SextetPVVVertex: public Helicity::VVVVertex {

public:

  /**
   * The default constructor.
   */
  SextetPVVVertex() {
    colourStructure(ColourStructure::DELTA);
  }

  /** Calculate the coupling
   *@param q2 The scale at which to evaluate the coupling
   *@param part1 The first interacting particle 
   *@param part2 The second interacting particle 
   *@param part3 The third interacting particle 
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
			   tcPDPtr part2,tcPDPtr part3);

public:

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
  SextetPVVVertex & operator=(const SextetPVVVertex &) = delete;

private:
  
  /**
   * Store the value of the coupling when last evaluated
   */
  Complex coupLast_;

  /**
   * Store the scale at which coupling was last evaluated
   */
  Energy2 q2Last_;

};

}

#endif /* THEPEG_SextetPVVVertex_H */
