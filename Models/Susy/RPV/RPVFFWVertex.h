// -*- C++ -*-
#ifndef Herwig_RPVFFWVertex_H
#define Herwig_RPVFFWVertex_H
//
// This is the declaration of the RPVFFWVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig/Models/Susy/SusyBase.h"
#include "Herwig/Models/Susy/MixingMatrix.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the RPVFFWVertex class.
 *
 * @see \ref RPVFFWVertexInterfaces "The interfaces"
 * defined for RPVFFWVertex.
 */
class RPVFFWVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  RPVFFWVertex();

  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr part1,
                           tcPDPtr part2, tcPDPtr part3);

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
  RPVFFWVertex & operator=(const RPVFFWVertex &) = delete;

private:

  /** @name Quark flavour parameters */
  //@{ 
  /**
   * True, if a diagonal CKM matrix should be assumed. This ignores
   * the CKM object of the StandardModel.
   */
  bool _diagonal;
  
  /**
   *  The elements of the CKM matrix.
   */
  vector<vector<Complex> > _ckm;
  //@}
  
  /**
   * Store \f$sin(\theta_W)\f$
   */
  double _sw;

  /**
   * Store the neutralino mixing matrix
   */
  tMixingMatrixPtr _theN;

  /**
   * Store the U-type chargino mixing matrix
   */
  tMixingMatrixPtr _theU;

  /**
   * Store the V-type chargino mixing matrix
   */
  tMixingMatrixPtr _theV;

 /**
   * The value of the coupling when it was last evaluated
   */
  Complex _couplast;
  
  /**
   * The scale at which the coupling was last evaluated
   */
  Energy2 _q2last;

  /**
   * The id of the first chargino the last time the vertex was evaluated
   */
  long _id1last;

  /**
   * The id of the second chargino the last time the vertex was evaluated
   */
  long _id2last;

  /**
   * The value of the left coupling when it was last evaluated
   */
  Complex _leftlast;

  /**
   * The value of the right coupling when it was last evaluated
   */
  Complex _rightlast;

  /**
   *  Which interactions to include
   */
  unsigned int _interactions;
};

}

#endif /* Herwig_RPVFFWVertex_H */
