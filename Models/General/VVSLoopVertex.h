// -*- C++ -*-
#ifndef HERWIG_VVSLoopVertex_H
#define HERWIG_VVSLoopVertex_H
// // This is the declaration of the VVSLoopVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/GeneralVVSVertex.h"
#include "VVSLoopVertex.fh"
 
namespace Herwig {
using namespace ThePEG;

/**
 * The <code>VVSLoopVertex</code> is designed to
 * calculate the coefficents for the terms in the
 * Passarino-Veltman tensor reduction scheme. A vertex
 * class should inherit from this and implement it's own 
 * setCoupling member from which the VVSLoopVertex  
 * setCoupling member is called.
 */
class VVSLoopVertex: public Helicity::GeneralVVSVertex {

public:

  /**
   * The default constructor.
   */
  VVSLoopVertex() : masses(0), type(0), couplings(0), Npart_(0), loopToolsInit_(false) {
    kinematics(true);
  }

  /** 
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param part1 ParticleData pointer to first particle
   *@param part2 ParticleData pointer to first particle
   *@param part3 ParticleData pointer to first particle
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1, tcPDPtr part2,
			   tcPDPtr part3);

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
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

protected:
  
  /**
   * Vector of loop masses
   */
  vector<Energy> masses;
  
  /**
   * Vector of loop types
   */
  vector<PDT::Spin> type;

  /**
   * The left and right couplings for a fermion loop 
   */
  vector<pair<Complex, Complex> > couplings;

  /**
   * Set the number of particles in the loop 
   */
  void setNParticles(unsigned int npart) { Npart_ = npart; }

  /**
   *  Is loopTools initialized
   */
  bool loopToolsInitialized() { return loopToolsInit_; }

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<VVSLoopVertex> initVVSLoopVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VVSLoopVertex & operator=(const VVSLoopVertex &);

private:
  
  /**
   * The number of particles in the loop 
   */
  unsigned int Npart_;

  /**
   *  Loop tools initialised ?
   */
  bool loopToolsInit_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VVSLoopVertex. */
template <>
struct BaseClassTrait<Herwig::VVSLoopVertex,1> {
  /** Typedef of the first base class of VVSLoopVertex. */
  typedef Helicity::GeneralVVSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VVSLoopVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VVSLoopVertex>
  : public ClassTraitsBase<Herwig::VVSLoopVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VVSLoopVertex"; }
};

/** @endcond */

}

#endif /* HERWIG_VVSLoopVertex_H */
