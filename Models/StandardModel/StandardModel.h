// -*- C++ -*-
#ifndef HERWIG_StandardModel_H
#define HERWIG_StandardModel_H
//
// This is the declaration of the <!id>StandardModel<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This is the Herwig++ <!id>StandardModel<!!id> class which inherits from ThePEG 
//  Standard Model class and implements additional Standard Model couplings, access
//  to vertices for helicity amplitude calculations etc.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="StandardModelBase.html">StandardModelBase.h</a>.
// 

#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/Models/StandardModel/RunningMassBase.h"
#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Helicity/Vertex/Vector/VVVVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Helicity/Vertex/Scalar/VVSVertex.h"
#include "Herwig++/Helicity/Vertex/Vector/VVVVVertex.h"

namespace Herwig {
using namespace ThePEG;
class StandardModel: public StandardModelBase {
  
// some typedefs for the pointers to vertices
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
  
  inline StandardModel();
  inline StandardModel(const StandardModel &);
  virtual ~StandardModel();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
  inline double lnu() const;
  inline double le() const;
  inline double lu() const;
  inline double ld() const;
  inline double rnu() const;
  inline double re() const;
  inline double ru() const;
  inline double rd() const;
  // The left and right couplings of the Z^0 including sin and cos theta_W
  
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
  // pointers to the objects handling the vertices
  
  inline double mass(Energy2 scale,tcPDPtr) const;
  // Return the running mass for a given scale and particle type
  
  inline trunPtr massPtr() const;
  // return a pointer to the object handling the running mass
  
protected:
  
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.
  
protected:
  
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.
  
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.
  
  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.
  
private:
  
  static ClassDescription<StandardModel> initStandardModel;
  // Describe a concrete class with persistent data.
  
  StandardModel & operator=(const StandardModel &);
  // Private and non-existent assignment operator.
  
  // pointers to the vertices
  // standard model fermion-antifermion gauge boson couplings
  FFVPtr _theFFZVertex;
  FFVPtr _theFFPVertex;
  FFVPtr _theFFGVertex;
  FFVPtr _theFFWVertex;
  // standard model fermion-antifermion scalar couplings
  FFSPtr _theFFHVertex;
  // standard model vector-vector-scalar couplings
  VVSPtr _theWWHVertex;
  // standard model vector-vector-vector couplings
  VVVPtr _theGGGVertex;
  VVVPtr _theWWWVertex;
  // standard model vector-vector-vector couplings
  VVVVPtr _theGGGGVertex;
  VVVVPtr _theWWWWVertex;
  // the running mass
  runPtr _theRunningMass;
  
};
}  
#include "StandardModel.icc"

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of StandardModel.
template <>
struct BaseClassTrait<Herwig::StandardModel,1> {
  typedef StandardModelBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::StandardModel>
  : public ClassTraitsBase<Herwig::StandardModel> {
  static string className() { return "/Herwig++/StandardModel"; }
  // Return the class name.
  static string library() { return "libHwStandardModel.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}


#endif /* HERWIG_StandardModel_H */
