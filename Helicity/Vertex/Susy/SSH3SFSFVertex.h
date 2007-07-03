// -*- C++ -*-
#ifndef HERWIG_SSH3SFSFVertex_H
#define HERWIG_SSH3SFSFVertex_H
//
// This is the declaration of the SSH3SFSFVertex class.
//

#include "Herwig++/Helicity/Vertex/Scalar/SSSVertex.h"
#include "Herwig++/Models/Susy/MSSM.h"
#include "SSH3SFSFVertex.fh"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/**
 * This is the implementation of the coupling of the \f$H_3^0\f$ higgs of
 * the MSSM to sfermion pairs. It inherits from SSSVertex and implements
 * the setCoupling member.
 * 
 * @see SSSVertex
 */
class SSH3SFSFVertex: public SSSVertex {

public:

  /**
   * The default constructor.
   */
  SSH3SFSFVertex();

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

  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
                           tcPDPtr part2,tcPDPtr part3);

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SSH3SFSFVertex> initSSH3SFSFVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSH3SFSFVertex & operator=(const SSH3SFSFVertex &);

  /**
   * Pointer to the MSSM object
   */
  tMSSMPtr _theSS;

  /**
   * \f$M_w\f$
   */
  Energy _mw;

  /**
   * \f$\sin(\theta_w)\f$
   */
  double _sw;
  
  /**
   * \f$\tan(\beta)\f$
   */
  double _tb;

  /**
   * \f$\mu\f$ paramter
   */
  Energy _mu;

  /**
   * vector holding the trilinear couplings quarks then leptons
   */
  vector<Energy> _trilinear;

  /**
   * The scale at which the coupling was last evaluated  
   */
  Energy2 _q2last;

  /**
   * The value of the EW coupling when last evaluated
   */
  Complex _glast;

  /**
   * Value of fermion dependent part when last evaluated
   */
  complex<Energy> _hfact;
   /**
   * Id of first scalar when coupling was last evaluated
   */
  long _id1last;

  /**
   * Id of first scalar when coupling was last evaluated
   */
  long _id2last;
};
}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSH3SFSFVertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::SSH3SFSFVertex,1> {
  /** Typedef of the first base class of SSH3SFSFVertex. */
  typedef Herwig::Helicity::SSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSH3SFSFVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::SSH3SFSFVertex>
  : public ClassTraitsBase<Herwig::Helicity::SSH3SFSFVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SSH3SFSFVertex"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SSH3SFSFVertex class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSusy.so HwSusyVertex.so"; }
};

/** @endcond */

}

#include "SSH3SFSFVertex.icc"

#endif /* HERWIG_SSH3SFSFVertex_H */
