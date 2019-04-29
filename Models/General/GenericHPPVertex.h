// -*- C++ -*-
#ifndef Herwig_GenericHPPVertex_H
#define Herwig_GenericHPPVertex_H
//
// This is the declaration of the GenericHPPVertex class.
//

#include "VVSLoopVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/SSSVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The GenericHPPVertex class implements the coupling of the Higgs boson to gluons
 * in a generic model using the other vertices coupling the Higgs boson
 * to coloured particles.
 *
 * @see \ref GenericHPPVertexInterfaces "The interfaces"
 * defined for GenericHPPVertex.
 */
class GenericHPPVertex: public VVSLoopVertex {

public:

  /**
   *  Struct to store the stuff needed for the vertices
   */
  struct Interaction {

    /**
     *  Constructor
     */
    Interaction(PDPtr in_particle=PDPtr(),SSSVertexPtr in_scalar=SSSVertexPtr(),
		FFSVertexPtr in_fermion=FFSVertexPtr(),
		VVSVertexPtr in_vector=VVSVertexPtr()) 
      : particle(in_particle), scalar(in_scalar), fermion(in_fermion), vector(in_vector)
    {}

    /**
     *  The particle
     */
    PDPtr particle;

    /**
     *  The scalar vertex
     */
    SSSVertexPtr scalar;

    /**
     *  The fermion vertex
     */
    FFSVertexPtr fermion;

    /**
     *  The vector vertex
     */
    VVSVertexPtr vector;
  };

public:

  /**
   * The default constructor.
   */
  GenericHPPVertex();

  /** 
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param part1 ParticleData pointer to first particle
   *@param part2 ParticleData pointer to first particle
   *@param part3 ParticleData pointer to first particle
   */
  virtual void setCoupling (Energy2 q2, tcPDPtr part1, tcPDPtr part2, tcPDPtr part3);

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

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GenericHPPVertex & operator=(const GenericHPPVertex &) = delete;

private:

  /**
   *  Initialize the vertex
   */
  void initializeVertex();

private:

  /**
   *  The bosons to consider as particles
   */
  vector<PDPtr> bosons_;

  /**
   *  Has it been set up?
   */
  bool setup_;

  /**
   *  The vertices to use
   */
  map<cPDPtr,vector<Interaction> > vertices_;
  
  /**
   * The scale \f$q^2\f$ at which coupling was last evaluated
   */
  Energy2 q2Last_;

  /**
   * Last value of the coupling calculated
   */
  Complex coupLast_;

  /**
   *  Lasst id of the boson
   */
  long idLast_;

  /**
   *  The model for running masses
   */
  tcHwSMPtr model_;
};

/**
 *  Persistent output of the Interaction struct
 */
inline PersistentOStream & operator<<(PersistentOStream & os, 
				      const GenericHPPVertex::Interaction  & x) {
  os << x.particle << x.scalar << x.fermion << x.vector;
  return os;
}

/**
 *  Persistent input of the Interaction struct
 */
inline PersistentIStream & operator>>(PersistentIStream & is,
				      GenericHPPVertex::Interaction & x) {
  is >> x.particle >> x.scalar >> x.fermion >> x.vector;
  return is;
}

}

#endif /* Herwig_GenericHPPVertex_H */
