// -*- C++ -*-
#ifndef HERWIG_PartnerFinder_H
#define HERWIG_PartnerFinder_H
//
// This is the declaration of the <!id>PartnerFinder<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is responsible of two related tasks:
// <OL>
//  <LI> it finds the partners (pairs of ShowerParticles objects) <BR> 
//       for each interaction types (QCD,QED,EWK...) and within the <BR>
//       scale range specified by the <!class>ShowerVariables<!!class> object;
//  <LI> for each pair of partners (and interaction therefore) <BR>
//       it sets the initial evolution scales of both of them.
// </OL>
// Notice that a given particle has, in general, a different partner <BR>
// for each different interaction; however, given a partner, its <BR>
// initial evolution scale, <!id>Q<!!id>, is purely a kinematical relationship <BR> 
// between the pair, without dependence on the dynamics (i.e. type of interaction). <BR>
// Therefore, there must be different methods for finding partners of different <BR>
// interaction types, but an unique common method to calculate the initial evolving <BR>
// showering scale. <BR>
// Notice that:
// <UL>
//  <LI> to be more general, one should define an abstract class 
//       <!id>AbsPartnerFinder<!!id> <BR>
//       which has, exactly like the present <!class>PartnerFinder<!!class> class, <BR>
//       the definition of the methods: <BR>
//         <!id>setQCDInitialEvolutionScales<!!id> <BR>
//         <!id>setQEDInitialEvolutionScales<!!id> <BR>
//         <!id>setEWKInitialEvolutionScales<!!id> <BR>
//     which are quite general, whereas the other one is declared pure virtual, <BR>
//     without definition: <BR>
//       <!id> virtual ...  calculateInitialEvolutionScales(...) = 0; <!!id> <BR>
//     and then having a concrete <!id>PartonFinder<!!id> class which inherits <BR>
//     from <!id>AbsPartnerFinder<!!id> and provides a definition for this virtual method. <BR>
//     In fact, it is only in this method that a specific choice of ordering variable <BR>
//     must be made. Therefore, if we wanted a different choice, we could define another <BR>
//     class, <!id>AlternativePartonFinder<!!id>, 
//     which also inherits from <!id>AbsPartnerFinder<!!id>, <BR>
//     but provides a different definition of the method 
//     <!id>calculateInitialEvolutionScales<!!id>. <BR>
// </UL>
//
// ***LOOKHERE*** 
//              - If the method  calculateInitialEvolutionScales
//                were not used anywhere else a part in this class,
//                then it should be moved in the private part of the class.
//              - It could be easy and straightforward to add to 
//                     setQCDInitialEvolutionScales  
//                (and, maybe, but less likely, also for the other ones)
//                either a pointer to the Collision object (which is
//                already available in the class InsideRangeShowerEvolver
//                that calls this class Partner Finder) --- in the
//                case some information from the hard subprocess is
//                necessary to find colour partners (expecially in the
//                case of baryon number nonconserving processes in 
//                Rp violating Susy) --- and/or a boolean flag that
//                tells whether or not the input particles are
//                associates to a decay --- and therefore finding
//                the colour partners, at least for the scale of
//                the decay, should not involved searches through
//                the particle history.
// ***endLOOKHERE***
//                

#include "ThePEG/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"


namespace Herwig {

using namespace ThePEG;

typedef pair<tShowerParticlePtr,tShowerParticlePtr> ShowerPPair;

class PartnerFinder: public ThePEG::HandlerBase {

public:

  inline PartnerFinder();
  inline PartnerFinder(const PartnerFinder &);
  virtual ~PartnerFinder();
  // Standard ctors and dtor.

  bool setQCDInitialEvolutionScales(const tShowerVarsPtr showerVariables,
				    const ShowerParticleVector particles,
                                    const bool isDecayCase = false);
  bool setQEDInitialEvolutionScales(const tShowerVarsPtr showerVariables,
				    const ShowerParticleVector particles,
                                    const bool isDecayCase = false);
  bool setEWKInitialEvolutionScales(const tShowerVarsPtr showerVariables,
				    const ShowerParticleVector particles,
                                    const bool isDecayCase = false);
  // Given in input a collection of particles (<!id>ShowerParticle<!!id> objects),
  // each of these methods set the initial evolution scales of those particles, 
  // between the ones given in input, that do not have yet their evolution scale set, 
  // according to a given interaction type (QCD, QED, EWK,...), and within the 
  // constraints specified in the <!id>showerVariables<!!id> object. 
  // The input collection of particles can be either the full collection of 
  // showering particles (kept in the main class <!class>ShowerHandler<!!class>),
  // in the case isDecayCase is false, or simply, in the case isDecayCase is true,
  // the decaying particle and its decay products.    
  // The methods returns true, unless something wrong (inconsistencies,
  // or undefined values) happens.

  pair<Energy,Energy> calculateInitialEvolutionScales(const ShowerPPair &particlePair, 
		                                      const tShowerVarsPtr showerVariables);
  // Given a pair of particles, supposedly partners w.r.t. an interaction,
  // this method returns their initial evolution scales as a pair.
  // If something wrong happens, it returns the null ( Energy() , Energy() ) pair. 
  // This method is used by the above <!id>setXXXInitialEvolutionScales<!!id> methods.

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

  static ClassDescription<PartnerFinder> initPartnerFinder;
  // Describe a concrete class with persistent data.

  PartnerFinder & operator=(const PartnerFinder &);
  //  Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of PartnerFinder.
template <>
struct BaseClassTrait<Herwig::PartnerFinder,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::PartnerFinder>: public ClassTraitsBase<Herwig::PartnerFinder> {
  static string className() { return "/Herwig++/PartnerFinder"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "PartnerFinder.icc"

#endif /* HERWIG_PartnerFinder_H */
