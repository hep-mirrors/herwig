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
   *  ShowerIndex is a simple struct which is used as key for the multimap of the 
   *  possible Sudakov form factors. In fact, given:
   *  -   a certain particle, specified by its id;  
   *  -   the interaction type (QCD or QED or EWK,
   *         or others if you want to add other extra interactions) of the 
   *         radiation it has emitted, in the case it didn't emit any radiation,
   *         and as default value, is set to UNDEFINED;
   *  - the time-order of the particle: IS if initial state (space-like), 
   *    or FS if final state (time-like), the latter is used as default;
   *  there can be more than one Sudakov form factor, for example
   *  a final-state gluon has the following Sudakov form factors: 
   *  \f$g\to d\bar{d}\f$, \f$g\to u\bar{u}\f$, \f$g\to s\bar{s}\f$,\f$g\to c\bar{c}\f$,
   *  \f$g\to b\bar{b}\f$, \f$g\to gg\f$
   *
   *  Notice that in principle we can even distinguish between
   *  particle (id>0) and antiparticle (id<0), but in practice
   *  we assume CP-conserving interactions, and therefore only
   *  the particle (id>0) will be considered.
   *  Before a ShowerParticle object emits eventually a radiation
   *  (that means forever in the case it does not emit any radiation)
   *  its interaction type is set to UNDEFINED.
   *
   *  We use a simple plain struct, rather than a proper class
   *  with encapsulated private members and accessory get/set methods,
   *  because we feel that this is more appropriate for this very simple and
   *  straightforward data structure which just consists of a triplet. 
   * 
   *  @see ShowerParticle
   *  @see SplittingGenerator
   */
  struct ShowerIndex {

    /**
     * This operator overloading is necessary in order to store persistently
     * the multimap of key = showerIndex object, value = pointer to 
     * Sudakov object (in class SplittingGenerator).
     */
    //@{
    /**
     * Output operator
     * @param os The output stream
     * @param x  The ShowerIndex being outputted.
     */
    friend PersistentOStream & operator<<(PersistentOStream & os, const ShowerIndex & x);

    /**
     * Input operator
     * @param is The input stream
     * @param x  The ShowerIndex being inputted.
     */
    friend PersistentIStream & operator>>(PersistentIStream & is, ShowerIndex & x);
    //@}

    /**
     *  Enumeration storing the number of possible interactions and directions
     */
    enum { NumInteractionTypes = 3, NumTimeOrderType = 2 };

    /**
     *  Enumeration for the interaction type
     */
    enum InteractionType { UNDEFINED=-1, QCD, QED, EWK };

    /**
     *  Enumeration storing the type of radiation, i.e. initial or final state.
     */  
    enum TimeOrderType { UNINITIALIZED=-1, IS, FS }; 

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
     *  The PDG code of the particle
     */
    long id;

    /**
     *  The interaction type
     */
    InteractionType interaction;

    /**
     *  The type of radiation
     */
    TimeOrderType timeFlag;

  };
  
} // end Herwig namespace


#endif // HERWIG_ShowerIndex_H
