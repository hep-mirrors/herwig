// -*- C++ -*-
#ifndef HERWIG_SplitFun1to2_H
#define HERWIG_SplitFun1to2_H
//
// This is the declaration of the <!id>SplitFun1to2<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This is an abstract class which defines the common interface <BR>
// for all <I>1-&GT;2</I> splitting functions, both for Initial State <BR>
// and Final State radiation. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun.html">SplitFun.h</a>, <BR>
// <a href="http:SplitFun1to3.html">SplitFun1to3.h</a>, <BR>
// <a href="http:QtoQGSplitFun.html">QtoQGSplitFun.h</a>, <BR>
// <a href="http:QtoQGammaSplitFun.html">QtoQGammaSplitFun.h</a>, <BR>
// <a href="http:GtoGGSplitFun.html">GtoGGSplitFun.h</a>, <BR>
// <a href="http:GtoQQbarSplitFun.html">GtoQQbarSplitFun.h</a>.
// 

#include "SplitFun.h"
#include "Pythia7/PDT/ParticleData.h"
#include "Pythia7/PDT/EnumParticles.h"
#include "ShowerColourLine.h"


namespace Herwig {

using namespace Pythia7;

class SplitFun1to2: public SplitFun {

public:

  inline SplitFun1to2();
  inline SplitFun1to2(const SplitFun1to2 &);
  virtual ~SplitFun1to2();
  // Standard ctors and dtor.

  inline SplitFun1to2( const ShowerIndex::InteractionType interaction,
		       const long inputIdEmitter,
		       const long inputIdFirstProduct, const long inputIdSecondProduct );
  inline SplitFun1to2( const ShowerIndex::InteractionType interaction,
		       const long inputIdEmitter, const Energy inputMassEmitter,
		       const long inputIdFirstProduct, const Energy inputMassFirstProduct,
		       const long inputIdSecondProduct, const Energy inputMassSecondProduct );
  // It Specifies the interaction type (QCD,QED,EWK,...) of the vertex <I>A -&GT; B + C</I>, 
  // the PDG id for the emitter (<I>A</I>),
  // the ids of the two emission products (<I>B</I> and <I>C</I>).
  // If the mass is not specified, then the nominal mass is assumed
  // (therefore for Susy sparticles with large widths their masses
  //  should be given because they can be different from the nominal ones).
  
  inline long idFirstProduct() const;
  inline Energy massFirstProduct() const;
  inline long idSecondProduct() const;
  inline Energy massSecondProduct() const;
  // PDG ids and masses of the two emission products.

  virtual Complex fullFun( const double z, const Energy2 qtilde2, 
			   const double phi ) = 0;
  virtual Complex integratedFun( const double z, const Energy2 qtilde2 ) = 0;
  virtual Complex fullFunWithHelicities( const double z, 
					 const Energy2 qtilde2, 
					 const double phi, 
					 const int h0, const int h1, 
					 const int h2 ) = 0;
  virtual Complex integratedFunWithHelicities( const double z,
					       const Energy2 qtilde2, 
					       const int h0, const int h1, 
					       const int h2 ) = 0;
  // Pure virtual methods, that must be provided by the derived
  // classes, which return the exact values of the splitting function
  // <I>0 -&GT; 1 + 2</I> evaluated in terms of some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, helicities
  // of the three particles.  Notice that, to be more general, the
  // returned value is complex, although for standard QCD it is indeed
  // real.  

  virtual Complex overestimateFullFun( const double z, const double phi );
  virtual Complex overestimateIntegratedFun( const double z );
  virtual Complex overestimateFullFunWithHelicities( const double z, const double phi,
						     const int h0, const int h1, const int h2 );
  virtual Complex overestimateIntegratedFunWithHelicities( const double z,
							   const int h0, const int h1, const int h2 );  
  // Virtual methods that could be overrided in a derived class in the
  // case it is known some approximate expressions for the exact
  // expressions of the splitting function (the ones provided by the
  // pure virtual methods above), and such that the former are always
  // above than the latter: <!id>| overestimateXXX(...) | &GT;= |
  // XXX(...) |. <!!id> The default definition of overestimateXXX(...) 
  // simply returns XXX(...).  Notice that we left the overestimated
  // functions and hence their integrals and inverses not depend on
  // qtilde2 since these are supposed to be as simple as possible.

  virtual Complex integOverIntegratedFun(const double z); 
  virtual Complex invIntegOverIntegratedFun(const double r); 
  // The indefinite integral of the overestimated splitting function
  // <!id>overestimateIntegratedFun(z)<!!id> and its inverse. 

  virtual void colourConnection( const ShoColinePair & parentShoColinePair,
				 ShoColinePair & firstProductShoColinePair,
				 ShoColinePair & secondProductShoColinePair );
  // The default definition of this method does nothing.
  // However, QCD 1-&GT;2 splitting functions which derive from this
  // class should override this method, in order to provide
  // the proper colour connection between the emitting parent
  // and the branching products. More precisely, the first argument
  // is a pair of pointers to <!id>ShowerColourLine<!!id> objects, 
  // which are associated with, respectively, the colour (first element
  // of the pair) and anticolour (second element of the pair)
  // of the emitting particle; one (but not both) of the these 
  // two pointers can be null (this is the case for quarks/antiquarks,
  // whereas gluons have both not null).
  // The other two arguments of the method are also pair of pointers
  // to <!id>ShowerColourLine<!!id> objects, for respectively the two 
  // branching products. Again, for each pair, the first element
  // is associated with the colour line and the second element
  // is associated with the anticolur line; one (but not both)
  // of these two elements can be null. The <!id>ShowerColourLine<!!id>
  // objects pointed by any of the two elements of both pairs
  // can be either one of object pointerd by the radiating parent,
  // or can be a new one (and because of that, the pointers must be
  // reference counted, not transient).
  // We prefer to use <!id>ShowerColourLine<!!id> in this method
  // rather than <!id>ShowerParticles<!!id> for two reason: first,
  // there is no coupling between <!id>SplitFun1to2<!!id> class
  // (and its derived classes) and the <!id>ShowerParticle<!!id> class; 
  // and second, we avoid to repeat for each type of splitting
  // function the same operation of setting the color lines to
  // proper particles. 

public:

  static void Init();
  // Standard Init function used to initialize the interfaces


private:

  static AbstractClassDescription<SplitFun1to2> initSplitFun1to2;
  // Describe an abstract base class with persistent data.

  SplitFun1to2 & operator=(const SplitFun1to2 &);
  //  Private and non-existent assignment operator.

  long _id1;
  Energy _m1;
  long _id2;
  Energy _m2;

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of SplitFun1to2.
template <>
struct BaseClassTrait<Herwig::SplitFun1to2,1> {
  typedef Herwig::SplitFun NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::SplitFun1to2>: public ClassTraitsBase<Herwig::SplitFun1to2> {
  static string className() { return "/Herwig++/SplitFun1to2"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "SplitFun1to2.icc"

#endif /* HERWIG_SplitFun1to2_H */
