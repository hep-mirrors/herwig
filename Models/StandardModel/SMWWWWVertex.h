// -*- C++ -*-
//
// SMWWWWVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMWWWWVertex_H
#define HERWIG_SMWWWWVertex_H
//
// This is the declaration of the SMWWWWVertex class.
//
#include "ThePEG/Helicity/Vertex/Vector/VVVVVertex.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {
using namespace ThePEG;
    
/** \ingroup Helicity
 *
 *  The SMWWWWVertex class is the implementation of the
 *  coupling of four electroweak gauge bosons in the SM. 
 *  It inherits from VVVVVertex nad implements the setCoupling member.
 *
 *  @see VVVVVVertex
 *  @see VertexBase
 */
class SMWWWWVertex: public VVVVVertex {
  
public:
  
  /**
   * Default constructor.
   */
  SMWWWWVertex();
  
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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param part4 The ParticleData pointer for the fourth particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3,
			   tcPDPtr part4);
  
protected:
  
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}
  
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<SMWWWWVertex> initSMWWWWVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMWWWWVertex & operator=(const SMWWWWVertex &) = delete;

  /**
   *  Intermediate particles
   */
  //@{
  /**
   * The ParticleData pointer for the photon
   */
  tcPDPtr _gamma;

  /**
   * The ParticleData pointer for the \f$Z^0\f$ boson.
   */
  tcPDPtr _Z0;

  /**
   * The ParticleData pointer for the \f$W^+\f$ boson.
   */
  tcPDPtr _wplus;

  /**
   * The ParticleData pointer for the \f$W^-\f$ boson.
   */
  tcPDPtr _wminus;
  //@}

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex _couplast;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 _q2last;

  /**
   *  The factors for the different bosons.
   */
  vector<double> _vfact;

  /**
   *  \f$\sin^2\theta_W\f$.
   */
  double _sw2;

  /**
   *  \f$\cos^2\theta_W\f$.
   */
  double _cw2;
  //@}
};
 
}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */
  
/**
 * The following template specialization informs ThePEG about the
 * base class of SMWWWWVertex.
 */
template <>
struct BaseClassTrait<Herwig::SMWWWWVertex,1> {
  /** Typedef of the base class of SMWWWWVertex. */
  typedef ThePEG::Helicity::VVVVVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::SMWWWWVertex>
  : public ClassTraitsBase<Herwig::SMWWWWVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::SMWWWWVertex"; }
  
};

/** @endcond */
  
}


#endif /* HERWIG_SMWWWWVertex_H */
