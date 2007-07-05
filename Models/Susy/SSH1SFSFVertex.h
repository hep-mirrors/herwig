// -*- C++ -*-
#ifndef HERWIG_SSH1SFSFVertex_H
#define HERWIG_SSH1SFSFVertex_H
//
// This is the declaration of the SSH1SFSFVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/SSSVertex.h"
#include "Herwig++/Models/Susy/MSSM.h"
#include "SSH1SFSFVertex.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * The SSH1SFSFVertex class implements the coupling of the \f$H_1^0\f$ higgs  
 * of the MSSM to the sfermion pairs. It inherits from SSSVertex and 
 * implements the setCoupling method 
 * 
 * @see SSSVertex
 */
class SSH1SFSFVertex: public SSSVertex {

public:

  /**
   * The default constructor.
   */
  SSH1SFSFVertex();

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
  static ClassDescription<SSH1SFSFVertex> initSSH1SFSFVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSH1SFSFVertex & operator=(const SSH1SFSFVertex &);
  
  /**
   * Pointer to the MSSM object
   */
  tMSSMPtr _theSS;

  /**
   * Higgs mixing angle \f$\sin(\alpha)\f$
   */
  double _sinAlpha;

  /**
   * Higgs mixing angle \f$\cos(\alpha)\f$
   */
  double _cosAlpha;
 
  /**
   * \f$\sin(\beta)\f$
   */
  double _sb;

  /**
   * \f$\cos(\beta)\f$
   */
  double _cb;
  
  /**
   * \f$M_z\f$
   */
  Energy _mz;
  
  /**
   * \f$M_w\f$
   */
  Energy _mw;

  /**
   * \f$\sin(\theta_w)\f$
   */
  double _sw;

  /**
   * \f$\cos(\theta_w)\f$
   */
  double _cw;

  /**
   * \f$\mu\f$
   */
  Energy _mu;

  /**
   * vector holding the trilinear couplings quarks then leptons
   */
  vector<Energy> _trilin;

  /**
   * \f$\sin(\alpha + \beta)\f$
   */
  double _sinAB;
  
  /**
   * stop mixing matrix pointer
   */
  tMixingMatrixPtr _stop;
  
  /**
   * sbottom mixing matrix pointer
   */
  tMixingMatrixPtr _sbottom;
  
  /**
   * stau mixing matrix
   */
  tMixingMatrixPtr _stau;

  /**
   * Value of EW coupling when last evaluated
   */
  Complex _glast;

  /**
   * Scale at which the coupling was last evaluated
   */
  Energy2 _q2last;

  /**
   * Value of mixing matrix dependent part when last evaluated
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
 *  base classes of SSH1SFSFVertex. */
template <>
struct BaseClassTrait<Herwig::SSH1SFSFVertex,1> {
  /** Typedef of the first base class of SSH1SFSFVertex. */
  typedef ThePEG::Helicity::SSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSH1SFSFVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSH1SFSFVertex>
  : public ClassTraitsBase<Herwig::SSH1SFSFVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SSH1SFSFVertex"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SSH1SFSFVertex class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#include "SSH1SFSFVertex.icc"

#endif /* HERWIG_SSH1SFSFVertex_H */
