// -*- C++ -*-
#ifndef HERWIG_ShowerIndex_H
#define HERWIG_ShowerIndex_H
//
// This is the declaration of the ShowerIndex class.

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


namespace Herwig {

  using namespace ThePEG;

  /** \ingroup Shower
   *
   *  Simple struct that is used as key for the multimap of the 
   *  possible Sudakov form factors. In fact, given: <BR>
   *  <OL>
   *    <LI> a certain particle, specified by its id; <BR>  
   *         (it must be always defined)
   *    <LI> the interaction type (QCD or QED or EWK,
   *         or others if you want to add other extra interactions) of the 
   *         radiation it has emitted; in the case it didn't emit any radiation,
   *         and as default value, is set to UNDEFINED;
   *    <LI> the time-order of the particle: IS if initial state (space-like), 
   *         or FS if final state (time-like); the latter is used as default;
   *  </OL>
   *  there can be more than one Sudakov form factor, for example:
   *  <I>(gluon,QCD,FS)</I> has the following Sudakov form factors: <BR>
   *  <I> g -&GT; u ubar , g -&GT; d dbar , g -&GT; s sbar , g -&GT; c cbar, 
   *      g -&GT; b bar </I><BR>
   *      (we can also include <I> g -&GT; g g g </I>). <BR>
   *  Notice that in practice we can even distinguish between
   *  particle (id>0) and antiparticle (id&LT;0), but in practice
   *  we assume CP-conserving interactions, and therefore only
   *  the particle (id&GT;0) will be considered.
   *  Before a ShowerParticle object emits eventually a radiation
   *  (that means forever in the case it does not emit any radiation)
   *  its interaction type is set to UNDEFINED. <BR>
   *
   *  We use a simple plain struct, rather than a proper class
   *  with encapsulated private members and accessory get/set methods,
   *  because we feel that this is more appropriate for this very simple and
   *  straightforward data structure which consists of just a triplet. 
   * 
   *  @see ShowerParticle
   *  @see SplittingGenerator
   */
  struct ShowerIndex {

    /**
     * This operator overloading is necessary in order to store persistently
     * the multimap of &LT; key = showerIndex object, value = pointer to 
     * Sudakov object &GT; (in class SplittingGenerator).
     */
    friend PersistentOStream & operator<<(PersistentOStream & os, const ShowerIndex & x);
    friend PersistentIStream & operator>>(PersistentIStream & is, ShowerIndex & x);

    enum { NumInteractionTypes = 3, NumTimeOrderType = 2 };
    enum InteractionType { UNDEFINED=-1, QCD, QED, EWK };  
    enum TimeOrderType { UNINITIALIZED=-1, IS, FS }; // InitialState, FinaleState particle

    /**
     * Default constructor.
     */
    ShowerIndex() : id(0), interaction(UNDEFINED), timeFlag(UNINITIALIZED) {}

    /**
     * This operator overloading is necessary in order to use ShowerIndex
     * as a key of a multimap (used in class SplittingGenerator).
     */
    bool operator< (const ShowerIndex & rhs) const; 

    /** 
     * These conversion static methods are necessary for overloading
     * the input operator&LT;&LT; (see above), because it is a compiling 
     * error if you try to convert an int to a enum const (whereas
     * the opposite conversion, which is used in the output operator&GT;&GT;
     * is done automatically by the compiler). Furthermore, these
     * methods, together with the the constants NumInteractionTypes
     * and NumTimeOrderType , are useful to loop over an enum,
     * which is not possible directly.
     */
    //static InteractionType int2Interaction(const int position);
    //static TimeOrderType   int2TimeOrder(const int position);

    long id;
    InteractionType interaction;
    TimeOrderType timeFlag;

  };
  
} // end Herwig namespace


#endif // HERWIG_ShowerIndex_H
