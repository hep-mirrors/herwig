// -*- C++ -*-
#ifndef HERWIG_ShowerHandler_H
#define HERWIG_ShowerHandler_H
//
// This is the declaration of the <!id>ShowerHandler<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is the main driver of the shower: it is responsible for <BR>
// the proper handling of all other specific collaborating classes <BR>
// and for the storing of the produced particles in the event record.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:CascadeHandler.html">CascadeHandler.h</a>, <BR>
// <a href="http:InsideRangeShowerEvolver.html">InsideRangeShowerEvolver.h</a>, <BR>
// <a href="http:MECorrections.html">MECorrections.h</a>, <BR>
// <a href="http:ShowerConstrainer.html">ShowerConstrainer.h</a>, <BR>
// <a href="http:ShowerParticle.html">ShowerParticle.h</a>
// 

#include "Pythia7/Handlers/CascadeHandler.h"
#include "Herwig++/Config/GlobalParameters.h"
#include "ShowerParticle.h"
#include "MECorrections.h"
#include "ShowerConstrainer.h"
#include "InsideRangeShowerEvolver.h"

namespace Herwig {


using namespace Pythia7;

class ShowerHandler: public Pythia7::CascadeHandler {

public:

  inline ShowerHandler();
  inline ShowerHandler(const ShowerHandler &);
  virtual ~ShowerHandler();
  // Standard ctors and dtor.

public:

  virtual void cascade();
  // The main method which manages the all showering.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<ShowerHandler> initShowerHandler;
  // Describe a concrete class with persistent data.

  ShowerHandler & operator=(const ShowerHandler &);
  //  Private and non-existent assignment operator.

  void createShowerParticlesFromP7Particles( const tPartCollHdlPtr ch, 
					    CollecShoParPtr & hardProcessParticles );
  // From the Pythia7 particles entering the hard subprocess, create
  // the corresponding starting <!id>ShowerParticle<!!id> objects and 
  // put them in the vector <!id>hardProcessParticles<!!id>. 
  // Notice that the transformation from Pythia7 <!id>ColourLine<!!id> 
  // objects into <!id>ShowerColourLine<!!id> ones must be properly handled.

  void fillPositions();
  // It fills the positions information for all the ShowerParticle 
  // objects in _particles, at the end of the showering.

  void debuggingInfo();
  // Print debugging information.

  void fillEvenRecord( const tPartCollHdlPtr ch );
  // At the end of the Showering, transform ShowerParticle objects
  // into Pythia7 particles and fill the event record with them.
  // Notice that the parent/child relationships and the 
  // transformation from ShowerColourLine objects into Pythia7
  // ColourLine ones must be properly handled.


  void eventShape(const tCollecShoParPtr & p, 
		  vector<double> & lam, vector<Vector3> & n);

  Ptr<GlobalParameters>::pointer _pointerGlobalParameters; 
  Ptr<MECorrections>::pointer _pointerMECorrections;
  Ptr<ShowerConstrainer>::pointer _pointerShowerConstrainer;
  Ptr<InsideRangeShowerEvolver>::pointer _pointerInsideRangeShowerEvolver;

  CollecShoParPtr _particles;   
 
};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of ShowerHandler.
template <>
struct BaseClassTrait<Herwig::ShowerHandler,1> {
  typedef Pythia7::CascadeHandler NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ShowerHandler>: public ClassTraitsBase<Herwig::ShowerHandler> {
  static string className() { return "/Herwig++/ShowerHandler"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "ShowerHandler.icc"

#endif /* HERWIG_ShowerHandler_H */
