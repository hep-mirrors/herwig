// -*- C++ -*-
//
// StandardModel.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_StandardModel_H
#define HERWIG_StandardModel_H
//
// This is the declaration of the StandardModel class.

#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/Models/StandardModel/RunningMassBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVVertex.h"
#include "Herwig++/Models/General/ModelGenerator.fh"
#include "StandardModel.fh"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/** \ingroup Models
 *  
 *  This is the Herwig++ StandardModel class which inherits from ThePEG 
 *  Standard Model class and implements additional Standard Model couplings, 
 *  access to vertices for helicity amplitude calculations etc.
 *
 *  @see StandardModelBase
 */
class StandardModel: public StandardModelBase {
  
  /**
   * Some typedefs for the pointers.
   */
  //@{

  /**
   * Pointer to the RunningMassBase object 
   */
  typedef Ptr<Herwig::RunningMassBase>::pointer runPtr;

  /**
   * Transient pointer to the RunningMassBase object 
   */
  typedef Ptr<Herwig::RunningMassBase>::transient_pointer trunPtr;
  //@}

public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  StandardModel();

  /**
   * Copy-constructor.
   */
  StandardModel(const StandardModel &);

  /**
   * Destructor
   */
  virtual ~StandardModel();
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
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

public:
  
  /**
   * The left and right couplings of the Z^0 including sin and cos theta_W.  
   */
  //@{
  /**
   *  The left-handed coupling of a neutrino
   */
  inline double lnu() const;

  /**
   *  The left-handed coupling of a charged lepton.
   */
  inline double le() const;

  /**
   *  The left-handed coupling of an up type quark.
   */
  inline double lu() const;

  /**
   *  The left-handed coupling of a down type quark.
   */
  inline double ld() const;

  /**
   *  The right-handed coupling of a neutrino
   */
  inline double rnu() const;

  /**
   *  The right-handed coupling of a charged lepton.
   */
  inline double re() const;

  /**
   *  The right-handed coupling of an up type quark.
   */
  inline double ru() const;

  /**
   *  The right-handed coupling of a down type quark.
   */
  inline double rd() const;
  //@}

  /**
   * Pointers to the objects handling the vertices.
   */
  //@{
  /**
   * Pointer to the fermion-fermion-Z vertex
   */
  inline tAbstractFFVVertexPtr  vertexFFZ() const;

  /**
   * Pointer to the fermion-fermion-photon vertex
   */
  inline tAbstractFFVVertexPtr  vertexFFP() const;

  /**
   * Pointer to the fermion-fermion-gluon vertex
   */
  inline tAbstractFFVVertexPtr  vertexFFG() const;

  /**
   * Pointer to the fermion-fermion-W vertex
   */
  inline tAbstractFFVVertexPtr  vertexFFW() const;

  /**
   * Pointer to the fermion-fermion-Higgs vertex
   */
  virtual inline tAbstractFFSVertexPtr  vertexFFH() const;

  /**
   * Pointer to the triple gluon vertex
   */
  inline tAbstractVVVVertexPtr  vertexGGG() const;

  /**
   * Pointer to the triple electroweak gauge boson vertex.
   */
  inline tAbstractVVVVertexPtr  vertexWWW() const;

  /**
   * Pointer to the two electroweak gauge boson Higgs vertex.
   */
  virtual inline tAbstractVVSVertexPtr  vertexWWH() const;

  /**
   * Pointer to the quartic electroweak gauge boson vertex.
   */
  inline tAbstractVVVVVertexPtr vertexWWWW() const;

  /**
   * Pointer to the quartic gluon vertex
   */
  inline tAbstractVVVVVertexPtr vertexGGGG() const;

 /**
   * Pointer to the quartic gluon vertex
   */
  virtual inline tAbstractVVSVertexPtr vertexHGG() const;

 /**
   * Pointer to the quartic gluon vertex
   */
  inline tAbstractVVSVertexPtr vertexHPP() const;

  /**
   *  Total number of vertices
   */
  inline unsigned int numberOfVertices() const;

  /**
   * Access to a vertex from the list
   */
  inline tVertexBasePtr vertex(unsigned int); 
  //@}  

  /**
   * Return the running mass for a given scale \f$q^2\f$ and particle type.
   * @param scale The scale \f$q^2\f$.
   * @param part The ParticleData object for the particle
   */
  inline Energy mass(Energy2 scale,tcPDPtr part) const;
  
  /**
   * Return a pointer to the object handling the running mass.
   */
  inline trunPtr massPtr() const;
  
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

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

protected:

  /**
   *  Add a vertex to the list
   */
  inline void addVertex(VertexBasePtr);

private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<StandardModel> initStandardModel;
  
  /** 
   * Private and non-existent assignment operator.
   */
  StandardModel & operator=(const StandardModel &);
  
private:

  /**
   * Pointers to the vertices for Standard Model helicity amplitude
   * calculations.
   */
  //@{
  /**
   * Pointer to the fermion-fermion-Z vertex
   */
  AbstractFFVVertexPtr _theFFZVertex;

  /**
   * Pointer to the fermion-fermion-photon vertex
   */
  AbstractFFVVertexPtr _theFFPVertex;

  /**
   * Pointer to the fermion-fermion-gluon vertex
   */
  AbstractFFVVertexPtr _theFFGVertex;

  /**
   * Pointer to the fermion-fermion-W vertex
   */
  AbstractFFVVertexPtr _theFFWVertex;

  /**
   * Pointer to the fermion-fermion-Higgs vertex
   */
  AbstractFFSVertexPtr _theFFHVertex;

  /**
   * Pointer to the two electroweak gauge boson Higgs vertex.
   */
  AbstractVVSVertexPtr _theWWHVertex;

  /**
   * Pointer to the triple gluon vertex
   */
  AbstractVVVVertexPtr _theGGGVertex;

  /**
   * Pointer to the triple electroweak gauge boson vertex.
   */
  AbstractVVVVertexPtr _theWWWVertex;

  /**
   * Pointer to the quartic gluon vertex
   */
  AbstractVVVVVertexPtr _theGGGGVertex;

  /**
   * Pointer to the quartic electroweak gauge boson vertex.
   */
  AbstractVVVVVertexPtr _theWWWWVertex;

  /**
   * Pointer to higgs-gluon-gluon vertex
   */
  AbstractVVSVertexPtr _theHGGVertex;

  /**
   * Pointer to higgs-gamma-gamma vertex
   */
  AbstractVVSVertexPtr _theHPPVertex; 
  
  /**
   *  Full list of vertices as a vector to allow searching
   */
  vector<VertexBasePtr> _vertexlist;
  //@}

  /**
   * The running mass.
   */
  runPtr _theRunningMass;

  /**
   * Pointer to ModelGenerator Class
   */
  ModelGeneratorPtr _theModelGenerator;
};

}  

#include "StandardModel.icc"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of StandardModel.
 */
template <>
struct BaseClassTrait<Herwig::StandardModel,1> {
    /** Typedef of the base class of StandardModel. */
  typedef StandardModelBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::StandardModel>
  : public ClassTraitsBase<Herwig::StandardModel> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::StandardModel"; }
};

/** @endcond */

}


#endif /* HERWIG_StandardModel_H */
