// -*- C++ -*-
#ifndef HERWIG_SVVLoopVertex_H
#define HERWIG_SVVLoopVertex_H
//
// This is the declaration of the SVVLoopVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/GeneralSVVVertex.h"
#include "ThePEG/PDT/PDT.h"
#include "SVVLoopVertex.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * The <code>SVVLoopVertex</code> is designed to
 * calculate the coefficents for the terms in the
 * Passarino-Veltman tensor reduction scheme. A vertex
 * class should inherit from this and implement it's own 
 * setCoupling member from which the SVVLoopVertex  
 * setCoupling member is called.
 */
class SVVLoopVertex: public Helicity::GeneralSVVVertex {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SVVLoopVertex();
  //@}
  
  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();
  
  /** 
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param part1 ParticleData pointer to first particle
   *@param part2 ParticleData pointer to first particle
   *@param part3 ParticleData pointer to first particle
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1, tcPDPtr part2,
			   tcPDPtr part3);
  
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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

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
  inline void setNParticles(unsigned int npart);
  
private:
  
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SVVLoopVertex> initSVVLoopVertex;
  
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SVVLoopVertex & operator=(const SVVLoopVertex &);

private:
  
  /**
   * The number of particles in the loop 
   */
  unsigned int theNpart;
};
  
}

#include "SVVLoopVertex.icc"

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SVVLoopVertex. */
template <>
struct BaseClassTrait<Herwig::SVVLoopVertex,1> {
  /** Typedef of the first base class of SVVLoopVertex. */
  typedef Helicity::GeneralSVVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SVVLoopVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SVVLoopVertex>
  : public ClassTraitsBase<Herwig::SVVLoopVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SVVLoopVertex"; }
};

/** @endcond */

}
 
#endif /* HERWIG_SVVLoopVertex_H */
