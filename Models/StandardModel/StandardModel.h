// -*- C++ -*-
#ifndef HERWIG_StandardModel_H
#define HERWIG_StandardModel_H
//
// This is the declaration of the StandardModel class.

#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/Models/StandardModel/RunningMassBase.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/Helicity/Vertex/Vector/VVVVVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/GeneralSVVVertex.h"
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
  inline tFFVVertexPtr  vertexFFZ() const;

  /**
   * Pointer to the fermion-fermion-photon vertex
   */
  inline tFFVVertexPtr  vertexFFP() const;

  /**
   * Pointer to the fermion-fermion-gluon vertex
   */
  inline tFFVVertexPtr  vertexFFG() const;

  /**
   * Pointer to the fermion-fermion-W vertex
   */
  inline tFFVVertexPtr  vertexFFW() const;

  /**
   * Pointer to the fermion-fermion-Higgs vertex
   */
  inline tFFSVertexPtr  vertexFFH() const;

  /**
   * Pointer to the triple gluon vertex
   */
  inline tVVVVertexPtr  vertexGGG() const;

  /**
   * Pointer to the triple electroweak gauge boson vertex.
   */
  inline tVVVVertexPtr  vertexWWW() const;

  /**
   * Pointer to the two electroweak gauge boson Higgs vertex.
   */
  inline tVVSVertexPtr  vertexWWH() const;

  /**
   * Pointer to the quartic electroweak gauge boson vertex.
   */
  inline tVVVVVertexPtr vertexWWWW() const;

  /**
   * Pointer to the quartic gluon vertex
   */
  inline tVVVVVertexPtr vertexGGGG() const;

 /**
   * Pointer to the quartic gluon vertex
   */
  inline tGeneralSVVVertexPtr vertexHGG() const;

 /**
   * Pointer to the quartic gluon vertex
   */
  inline tGeneralSVVVertexPtr vertexHPP() const;

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
  virtual void doinit() throw(InitException);
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
  FFVVertexPtr _theFFZVertex;

  /**
   * Pointer to the fermion-fermion-photon vertex
   */
  FFVVertexPtr _theFFPVertex;

  /**
   * Pointer to the fermion-fermion-gluon vertex
   */
  FFVVertexPtr _theFFGVertex;

  /**
   * Pointer to the fermion-fermion-W vertex
   */
  FFVVertexPtr _theFFWVertex;

  /**
   * Pointer to the fermion-fermion-Higgs vertex
   */
  FFSVertexPtr _theFFHVertex;

  /**
   * Pointer to the two electroweak gauge boson Higgs vertex.
   */
  VVSVertexPtr _theWWHVertex;

  /**
   * Pointer to the triple gluon vertex
   */
  VVVVertexPtr _theGGGVertex;

  /**
   * Pointer to the triple electroweak gauge boson vertex.
   */
  VVVVertexPtr _theWWWVertex;

  /**
   * Pointer to the quartic gluon vertex
   */
  VVVVVertexPtr _theGGGGVertex;

  /**
   * Pointer to the quartic electroweak gauge boson vertex.
   */
  VVVVVertexPtr _theWWWWVertex;

  /**
   * Pointer to higgs-gluon-gluon vertex
   */
  GeneralSVVVertexPtr _theHGGVertex;

  /**
   * Pointer to higgs-gamma-gamma vertex
   */
  GeneralSVVVertexPtr _theHPPVertex; 
  
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
