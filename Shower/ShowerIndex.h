// -*- C++ -*-
#ifndef HERWIG_ShowerIndex_H
#define HERWIG_ShowerIndex_H
//
// This is the declaration of the <!id>ShowerIndex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// Simple struct that is used as key for the multimap of the <BR>
// possible Sudakov form factors. In fact, given: <BR>
// <OL>
//   <LI> a certain particle, specified by its <!id>id<!!id>; <BR> 
//        (it must be always defined)
//   <LI> the interaction type (<!id>QCD<!!id> or <!id>QED<!!id> or <!id>EWK<!!id>, <BR>
//        or others if you want to add other extra interactions) of the <BR>
//        radiation it has emitted; in the case it didn't emit any radiation, <BR>
//        and as default value, is set to <!id>UNDEFINED<!!id>;
//   <LI> the time-order of the particle: <!id>IS<!!id> if initial state (space-like), <BR> 
//        or <!id>FS<!!id> if final state (time-like); the latter is used as default;
// </OL>
// there can be more than one Sudakov form factor, for example: <BR>
// <I>(gluon,QCD,FS)</I> has the following Sudakov form factors: <BR>
// <I> g -&GT; u ubar , g -&GT; d dbar , g -&GT; s sbar , g -&GT; c cbar, g -&GT; b bar </I><BR>
// (we can also include <I> g -&GT; g g g </I>). <BR>
// Notice that in practice we can even distinguish between <BR>
// particle (id>0) and antiparticle (id&LT;0), but in practice <BR>
// we assume CP-conserving interactions, and therefore only <BR>
// the particle (id&GT;0) will be considered. <BR>
// Before a ShowerParticle object emits eventually a radiation <BR>
// (that means forever in the case it does not emit any radiation) <BR>
// its interaction type is set to <!id>UNDEFINED<!!id>. <BR>
//
// We use a simple plain struct, rather than a proper class <BR> 
// with encapsulated private members and accessory get/set methods, <BR>
// because we feel that this is more appropriate for this very simple and <BR>
// straightforward data structure which consists of just a triplet. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:ShowerParticle.html">ShowerParticle.h</a>, <BR>
// <a href="http:SplittingGenerator.html">SplittingGenerator.h</a>.
// 
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


namespace Herwig {

  using namespace ThePEG;

  struct ShowerIndex {

    friend PersistentOStream & operator<< (PersistentOStream & os, const ShowerIndex & x);
    friend PersistentIStream & operator>> (PersistentIStream & is, ShowerIndex & x);
    // This operator overloading is necessary in order to store persistently
    // the multimap of &LT; key = showerIndex object, value = pointer to Sudakov object &GT;
    // (in class <!class>SplittingGenerator<!!class>).

    enum { NumInteractionTypes = 3, NumTimeOrderType = 2 };
    enum InteractionType { UNDEFINED=-1, QCD, QED, EWK };  
    enum TimeOrderType { UNINITIALIZED=-1, IS, FS }; // InitialState, FinaleState particle

    ShowerIndex() : id(0), interaction(UNDEFINED), timeFlag(UNINITIALIZED) {}
    // Default constructor.

    bool operator< (const ShowerIndex & rhs) const; 
    // This operator overloading is necessary in order to use <!id>ShowerIndex<!!id>
    // as a key of a multimap (used in class <!class>SplittingGenerator<!!class>).

    static InteractionType int2Interaction(const int position);
    static TimeOrderType   int2TimeOrder(const int position);
    // These conversion static methods are necessary for overloading
    // the input operator&LT;&LT; (see above), because it is a compiling 
    // error if you try to convert an int to a enum const (whereas
    // the opposite conversion, which is used in the output operator&GT;&GT;
    // is done automatically by the compiler). Furthermore, these
    // methods, together with the the constants <!id>NumInteractionTypes<!!id>
    // and <!id>NumTimeOrderType<!!id> , are useful to loop over an enum,
    // which is not possible directly.

    long id;
    InteractionType interaction;
    TimeOrderType timeFlag;

  };
  
} // end Herwig namespace


#endif // HERWIG_ShowerIndex_H
