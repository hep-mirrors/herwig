// -*- C++ -*-
#ifndef HERWIG_SplitFun_H
#define HERWIG_SplitFun_H
//
// This is the declaration of the <!id>SplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This is abstract class from which all splitting function classes, <BR>
// whatever their interaction type (QCD, QED, EWK,...) and multiplicity <BR>
// <I>1-&GT;N</I>, are inheriting from.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun1to2.html">SplitFun1to2.h</a>, <BR>
// <a href="http:SplitFun1to3.html">SplitFun1to3.h</a>.

#include "ShowerConfig.h"
#include "ThePEG/Pointer/Ptr.h"
#include "ThePEG/Pointer/ReferenceCounted.h"
#include "ThePEG/Pointer/PtrTraits.h"
#include "ThePEG/Pointer/RCPtr.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ShowerIndex.h"


namespace Herwig {

using namespace ThePEG;

class SplitFun: public ReferenceCounted  {

public:

  inline SplitFun();
  inline SplitFun(const SplitFun &);
  virtual ~SplitFun();
  // Standard ctors and dtor.

  inline SplitFun( const ShowerIndex::InteractionType interaction,
		   const int inputNumBranchingProducts,
                   const long inputIdEmitter, const Energy inputMassEmitter);
  // Specifies the interaction type of the vertex <I>0-&GT;1+2+...+N</I> 
  // the number of branching products (<I>N</I>), and the PDG id and mass of the emitter.
  // (Notice that the id and masses of branching products are of 
  // responsability of the classes that inherit from this class).

  inline ShowerIndex::InteractionType interactionType() const;
  // Type (QCD, QED, EWK,...) of the emission (branching) vertex.
  
  inline int numBranchingProducts() const;
  // The number of branching products.

  inline long idEmitter() const;
  inline Energy massEmitter() const;
  // PDG id and mass of the emitting (showering) particle.

private:

  SplitFun & operator=(const SplitFun &);
  //  Private and non-existent assignment operator.
  
  ShowerIndex::InteractionType _interaction;
  int _numProducts;
  long _idEmitter;
  Energy _mEmitter;

};

}

#include "SplitFun.icc"

#endif /* HERWIG_SplitFun_H */
