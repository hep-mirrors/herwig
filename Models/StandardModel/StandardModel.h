// -*- C++ -*-
#ifndef HERWIG_StandardModel_H
#define HERWIG_StandardModel_H
//
// This is the declaration of the StandardModel class.

#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/Models/StandardModel/RunningMassBase.h"
#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Helicity/Vertex/Vector/VVVVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/VVSVertex.h"
#include "Herwig++/Helicity/Vertex/Vector/VVVVVertex.h"

namespace Herwig {

using namespace ThePEG;

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
   * Some typedefs for the pointers to vertices.
   */
  typedef Ptr<Herwig::Helicity::FFVVertex>::pointer FFVPtr;
  typedef Ptr<Herwig::Helicity::FFVVertex>::transient_pointer tFFVPtr;
  typedef Ptr<Herwig::Helicity::VVVVertex>::pointer VVVPtr;
  typedef Ptr<Herwig::Helicity::VVVVertex>::transient_pointer tVVVPtr;
  typedef Ptr<Herwig::Helicity::FFSVertex>::pointer FFSPtr;
  typedef Ptr<Herwig::Helicity::FFSVertex>::transient_pointer tFFSPtr;
  typedef Ptr<Herwig::Helicity::VVSVertex>::pointer VVSPtr;
  typedef Ptr<Herwig::Helicity::VVSVertex>::transient_pointer tVVSPtr;
  typedef Ptr<Herwig::Helicity::VVVVVertex>::pointer VVVVPtr;
  typedef Ptr<Herwig::Helicity::VVVVVertex>::transient_pointer tVVVVPtr;
  typedef Ptr<Herwig::RunningMassBase>::pointer runPtr;
  typedef Ptr<Herwig::RunningMassBase>::transient_pointer trunPtr;

public:
  
  /**
   * Standard ctors and dtor.
   */
  inline StandardModel();
  inline StandardModel(const StandardModel &);
  virtual ~StandardModel();
  
public:
  
  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
  /**
   * The left and right couplings of the Z^0 including sin and cos theta_W.  
   */
  inline double lnu() const;
  inline double le() const;
  inline double lu() const;
  inline double ld() const;
  inline double rnu() const;
  inline double re() const;
  inline double ru() const;
  inline double rd() const;
  
  /**
   * Pointers to the objects handling the vertices.
   */
  inline tFFVPtr  vertexFFZ() const;
  inline tFFVPtr  vertexFFP() const;
  inline tFFVPtr  vertexFFG() const;
  inline tFFVPtr  vertexFFW() const;
  inline tFFSPtr  vertexFFH() const;
  inline tVVVPtr  vertexGGG() const;
  inline tVVVPtr  vertexWWW() const;
  inline tVVSPtr  vertexWWH() const;
  inline tVVVVPtr vertexWWWW() const;
  inline tVVVVPtr vertexGGGG() const;
  
  /**
   * Return the running mass for a given scale and particle type.
   */
  inline double mass(Energy2 scale,tcPDPtr) const;
  
  /**
   * Return a pointer to the object handling the running mass.
   */
  inline trunPtr massPtr() const;
  
protected:
  
  /**
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  
protected:
  
  /** 
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  
  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  
  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();
  
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<StandardModel> initStandardModel;
  
  /** 
   * Private and non-existent assignment operator.
   */
  StandardModel & operator=(const StandardModel &);
  
  /**
   * Pointers to the vertices standard model fermion-antifermion 
   * gauge boson couplings.
   */
  FFVPtr _theFFZVertex;
  FFVPtr _theFFPVertex;
  FFVPtr _theFFGVertex;
  FFVPtr _theFFWVertex;

  /**
   * Standard model fermion-antifermion scalar couplings.
   */
  FFSPtr _theFFHVertex;

  /**
   * Standard model vector-vector-scalar couplings.
   */
  VVSPtr _theWWHVertex;

  /**
   * Standard model vector-vector-vector couplings.
   */
  VVVPtr _theGGGVertex;
  VVVPtr _theWWWVertex;

  /**
   * Standard model vector-vector-vector couplings.
   */
  VVVVPtr _theGGGGVertex;
  VVVVPtr _theWWWWVertex;

  /**
   * The running mass.
   */
  runPtr _theRunningMass;
  
};

}  

#include "StandardModel.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of StandardModel.
 */
template <>
struct BaseClassTrait<Herwig::StandardModel,1> {
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
  static string className() { return "/Herwig++/StandardModel"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwStandardModel.so"; }

};

}


#endif /* HERWIG_StandardModel_H */
