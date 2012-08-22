// -*- C++ -*-
#ifndef Herwig_RPVHGGVertex_H
#define Herwig_RPVHGGVertex_H
//
// This is the declaration of the RPVHGGVertex class.
//

#include "Herwig++/Models/General/VVSLoopVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the RPVHGGVertex class.
 *
 * @see \ref RPVHGGVertexInterfaces "The interfaces"
 * defined for RPVHGGVertex.
 */
class RPVHGGVertex: public VVSLoopVertex {

public:

  /**
   * The default constructor.
   */
  RPVHGGVertex();
  
  /** 
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param particle1 ParticleData pointer to first particle
   *@param particle2 ParticleData pointer to second particle
   *@param particle3 ParticleData pointer to third particle
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			   tcPDPtr particle3);

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RPVHGGVertex & operator=(const RPVHGGVertex &);

};

}

#endif /* Herwig_RPVHGGVertex_H */
