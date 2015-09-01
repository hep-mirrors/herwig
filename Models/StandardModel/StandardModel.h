// -*- C++ -*-
//
// StandardModel.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_StandardModel_H
#define HERWIG_StandardModel_H
//
// This is the declaration of the StandardModel class.

#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Models/StandardModel/RunningMassBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractSSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractSSSSVertex.h"
#include "Herwig/Models/General/ModelGenerator.fh"
#include "StandardModel.fh"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/** \ingroup Models
 *  
 *  This is the Herwig StandardModel class which inherits from ThePEG 
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

protected:
  /**
   * Should the default vertices be considered for generic diagrams
   */
  virtual bool registerDefaultVertices() const { return true; }

public:
  
  /**
   * The left and right couplings of the Z^0 including sin and cos theta_W.  
   */
  //@{
  /**
   *  The left-handed coupling of a neutrino
   */
  double lnu() const {
    return 0.25/sqrt(sin2ThetaW()*(1.-sin2ThetaW()))*(vnu()+anu());
  }

  /**
   *  The left-handed coupling of a charged lepton.
   */
  double le() const {
    return 0.25/sqrt(sin2ThetaW()*(1.-sin2ThetaW()))*(ve()+ae());
  }

  /**
   *  The left-handed coupling of an up type quark.
   */
  double lu() const {
    return 0.25/sqrt(sin2ThetaW()*(1.-sin2ThetaW()))*(vu()+au());
  }

  /**
   *  The left-handed coupling of a down type quark.
   */
  double ld() const {
    return 0.25/sqrt(sin2ThetaW()*(1.-sin2ThetaW()))*(vd()+ad());
  }

  /**
   *  The right-handed coupling of a neutrino
   */
  double rnu() const {
    return 0.25/sqrt(sin2ThetaW()*(1.-sin2ThetaW()))*(vnu()-anu());
  }

  /**
   *  The right-handed coupling of a charged lepton.
   */
  double re() const {
    return 0.25/sqrt(sin2ThetaW()*(1.-sin2ThetaW()))*(ve()-ae());
  }

  /**
   *  The right-handed coupling of an up type quark.
   */
  double ru() const {
    return 0.25/sqrt(sin2ThetaW()*(1.-sin2ThetaW()))*(vu()-au());
  }

  /**
   *  The right-handed coupling of a down type quark.
   */
  double rd() const {
    return 0.25/sqrt(sin2ThetaW()*(1.-sin2ThetaW()))*(vd()-ad());
  }
  //@}

  /**
   * Pointers to the objects handling the vertices.
   */
  //@{
  /**
   * Pointer to the fermion-fermion-Z vertex
   */
  virtual tAbstractFFVVertexPtr  vertexFFZ() const {
    return FFZVertex_;
  }

  /**
   * Pointer to the fermion-fermion-photon vertex
   */
  virtual tAbstractFFVVertexPtr  vertexFFP() const {
    return FFPVertex_;
  }

  /**
   * Pointer to the fermion-fermion-gluon vertex
   */
  virtual tAbstractFFVVertexPtr  vertexFFG() const {
    return FFGVertex_;
  }
  
  /**
   * Pointer to the fermion-fermion-W vertex
   */
  virtual tAbstractFFVVertexPtr  vertexFFW() const {
    return FFWVertex_;
  }

  /**
   * Pointer to the fermion-fermion-Higgs vertex
   */
  virtual tAbstractFFSVertexPtr  vertexFFH() const {
    return FFHVertex_;
  }

  /**
   * Pointer to the triple gluon vertex
   */
  virtual tAbstractVVVVertexPtr  vertexGGG() const {
    return GGGVertex_;
  }
  
  /**
   * Pointer to the triple electroweak gauge boson vertex.
   */
  virtual tAbstractVVVVertexPtr  vertexWWW() const {
    return WWWVertex_;
  }

  /**
   * Pointer to the two electroweak gauge boson Higgs vertex.
   */
  virtual tAbstractVVSVertexPtr  vertexWWH() const {
    return WWHVertex_;
  }

  /**
   * Pointer to the quartic electroweak gauge boson vertex.
   */
  virtual tAbstractVVVVVertexPtr vertexWWWW() const {
    return WWWWVertex_;
  }

  /**
   * Pointer to the quartic gluon vertex
   */
  virtual tAbstractVVVVVertexPtr vertexGGGG() const {
    return GGGGVertex_;
  }
  
  /**
   * Pointer to the quartic gluon vertex
   */
  virtual tAbstractVVSVertexPtr vertexHGG() const {
    return HGGVertex_;
  }

  /**
   * Pointer to the quartic gluon vertex
   */
  virtual tAbstractVVSVertexPtr vertexHPP() const {
    return HPPVertex_;
  }

  /**
   * Pointer to the triple Higgs vertex
   */
  virtual tAbstractSSSVertexPtr vertexHHH() const {
    return HHHVertex_;
  }

  /**
   * Pointer to the WWHH vertex
   */
  virtual tAbstractVVSSVertexPtr vertexWWHH() const {
    return WWHHVertex_;
  }

  /**
   *  Total number of vertices
   */
  unsigned int numberOfVertices() const {
    return vertexList_.size() + extraVertices_.size();
  }

  /**
   * Access to a vertex from the list
   */
  tVertexBasePtr vertex(size_t ix) const {
    const size_t S = vertexList_.size();
    if ( ix < S )
      return vertexList_[ix];
    else
      return extraVertices_[ix - S];
  }
  //@}  

  /**
   * Return the running mass for a given scale \f$q^2\f$ and particle type.
   * @param scale The scale \f$q^2\f$.
   * @param part The ParticleData object for the particle
   */
  Energy mass(Energy2 scale,tcPDPtr part) const {
    return runningMass_->value(scale,part);
  }
  
  /**
   * Return a pointer to the object handling the running mass.
   */
  trunPtr massPtr() const {
    return runningMass_;
  }
  
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
  void addVertex(VertexBasePtr in) {
    if ( in )
      vertexList_.push_back(in);
  }

  /**
   *  Helper function to reset the mass of a ParticleData object
   *  in BSM models.
   */
  void resetMass(long id, Energy mass,tPDPtr particle=tPDPtr());

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
  AbstractFFVVertexPtr FFZVertex_;

  /**
   * Pointer to the fermion-fermion-photon vertex
   */
  AbstractFFVVertexPtr FFPVertex_;

  /**
   * Pointer to the fermion-fermion-gluon vertex
   */
  AbstractFFVVertexPtr FFGVertex_;

  /**
   * Pointer to the fermion-fermion-W vertex
   */
  AbstractFFVVertexPtr FFWVertex_;

  /**
   * Pointer to the fermion-fermion-Higgs vertex
   */
  AbstractFFSVertexPtr FFHVertex_;

  /**
   * Pointer to the two electroweak gauge boson Higgs vertex.
   */
  AbstractVVSVertexPtr WWHVertex_;

  /**
   * Pointer to the triple gluon vertex
   */
  AbstractVVVVertexPtr GGGVertex_;

  /**
   * Pointer to the triple electroweak gauge boson vertex.
   */
  AbstractVVVVertexPtr WWWVertex_;

  /**
   * Pointer to the quartic gluon vertex
   */
  AbstractVVVVVertexPtr GGGGVertex_;

  /**
   * Pointer to the quartic electroweak gauge boson vertex.
   */
  AbstractVVVVVertexPtr WWWWVertex_;

  /**
   * Pointer to higgs-gluon-gluon vertex
   */
  AbstractVVSVertexPtr HGGVertex_;

  /**
   * Pointer to higgs-gamma-gamma vertex
   */
  AbstractVVSVertexPtr HPPVertex_; 

  /**
   * Pointer to triple Higgs vertex
   */
  AbstractSSSVertexPtr HHHVertex_; 

  /**
   * Pointer to  WWHH vertex
   */
  AbstractVVSSVertexPtr WWHHVertex_; 
  
  /**
   *  Full list of vertices as a vector to allow searching
   */
  vector<VertexBasePtr> vertexList_;

  /**
   * Additional vertices to be considered in automatic ME construction
   */
  vector<VertexBasePtr> extraVertices_;
  //@}

  /**
   * The running mass.
   */
  runPtr runningMass_;

  /**
   * Pointer to ModelGenerator Class
   */
  ModelGeneratorPtr modelGenerator_;

};

}

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
