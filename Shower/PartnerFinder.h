// -*- C++ -*-
#ifndef HERWIG_PartnerFinder_H
#define HERWIG_PartnerFinder_H
//
// This is the declaration of the PartnerFinder class.

#include "ThePEG/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"


namespace Herwig {

using namespace ThePEG;

typedef pair<tShowerParticlePtr,tShowerParticlePtr> ShowerPPair;

/** \ingroup Shower
 *
 *  This class is responsible of two related tasks:
 *  <OL>
 *  <LI> it finds the partners (pairs of ShowerParticles objects)
 *       for each interaction types (QCD,QED,EWK...) and within the
 *       scale range specified by the ShowerVariables object;
 *  <LI> for each pair of partners (and interaction therefore)
 *       it sets the initial evolution scales of both of them.
 *  </OL>
 *  Notice that a given particle has, in general, a different partner
 *  for each different interaction; however, given a partner, its 
 *  initial evolution scale, Q, is purely a kinematical relationship 
 *  between the pair, without dependence on the dynamics (i.e. type of interaction).
 *  Therefore, there must be different methods for finding partners of different
 *  interaction types, but an unique common method to calculate the initial evolving
 *  showering scale.
 *  Notice that:
 *  <UL>
 *  <LI> to be more general, one should define an abstract class 
 *       AbsPartnerFinder which has, exactly like the present PartnerFinder
 *       class, the definition of the methods: <BR>
 *          setQCDInitialEvolutionScales <BR>
 *          setQEDInitialEvolutionScales <BR>
 *          setEWKInitialEvolutionScales <BR>
 *     which are quite general, whereas the other one is declared pure virtual, <BR>
 *     without definition: <BR>
 *       virtual ...  calculateInitialEvolutionScales(...) = 0;  <BR>
 *     and then having a concrete PartonFinder class which inherits
 *     from AbsPartnerFinder and provides a definition for this virtual method. 
 *     In fact, it is only in this method that a specific choice of ordering 
 *     variable must be made. Therefore, if we wanted a different choice, 
 *     we could define another class, AlternativePartonFinder, which also 
 *     inherits from AbsPartnerFinder, but provides a different definition 
 *     of the method calculateInitialEvolutionScales.
 *  </UL>
 *
 *  ***LOOKHERE*** 
 *              - If the method  calculateInitialEvolutionScales
 *                were not used anywhere else a part in this class,
 *                then it should be moved in the private part of the class.
 *              - It could be easy and straightforward to add to 
 *                     setQCDInitialEvolutionScales  
 *                (and, maybe, but less likely, also for the other ones)
 *                either a pointer to the Collision object (which is
 *                already available in the class InsideRangeShowerEvolver
 *                that calls this class Partner Finder) --- in the
 *                case some information from the hard subprocess is
 *                necessary to find colour partners (expecially in the
 *                case of baryon number nonconserving processes in 
 *                Rp violating Susy) --- and/or a boolean flag that
 *                tells whether or not the input particles are
 *                associates to a decay --- and therefore finding
 *                the colour partners, at least for the scale of
 *                the decay, should not involved searches through
 *                the particle history.
 *  ***endLOOKHERE***
 */                
class PartnerFinder: public ThePEG::HandlerBase {

public:

  /**
   * Standard ctors and dtor.
   */
  inline PartnerFinder();
  inline PartnerFinder(const PartnerFinder &);
  virtual ~PartnerFinder();

  /**
   * Given in input a collection of particles (ShowerParticle objects),
   * each of these methods set the initial evolution scales of those particles, 
   * between the ones given in input, that do not have yet their evolution scale set, 
   * according to a given interaction type (QCD, QED, EWK,...), and within the 
   * constraints specified in the showerVariables object. 
   * The input collection of particles can be either the full collection of 
   * showering particles (kept in the main class ShowerHandler,
   * in the case isDecayCase is false, or simply, in the case isDecayCase 
   * is true, the decaying particle and its decay products.    
   * The methods returns true, unless something wrong (inconsistencies,
   * or undefined values) happens.
   */
  bool setQCDInitialEvolutionScales(const tShowerVarsPtr showerVariables,
				    const ShowerParticleVector &particles,
                                    const bool isDecayCase = false);
  bool setQEDInitialEvolutionScales(const tShowerVarsPtr showerVariables,
				    const ShowerParticleVector &particles,
                                    const bool isDecayCase = false);
  bool setEWKInitialEvolutionScales(const tShowerVarsPtr showerVariables,
				    const ShowerParticleVector &particles,
                                    const bool isDecayCase = false);

  /**
   * Given a pair of particles, supposedly partners w.r.t. an interaction,
   * this method returns their initial evolution scales as a pair.
   * If something wrong happens, it returns the null (Energy(),Energy()) pair. 
   * This method is used by the above setXXXInitialEvolutionScales 
   * methods.
   */
  pair<Energy,Energy> calculateInitialEvolutionScales(const ShowerPPair &, 
		                                      const tShowerVarsPtr);
  pair<Energy,Energy> calculateFinalFinalScales(const ShowerPPair &,
		                                const tShowerVarsPtr);
  pair<Energy,Energy> calculateInitialInitialScales(const ShowerPPair &,
		                                    const tShowerVarsPtr);
  pair<Energy,Energy> calculateInitialFinalScales(const ShowerPPair &,
		                                  const tShowerVarsPtr);

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
  static ClassDescription<PartnerFinder> initPartnerFinder;

  /**
   * Private and non-existent assignment operator.
   */
  PartnerFinder & operator=(const PartnerFinder &);


  /**
   *  Approach to use for setting the colour partners
   */
  int _approach;

};

}

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of PartnerFinder.
 */
template <>
struct BaseClassTrait<Herwig::PartnerFinder,1> {
  typedef ThePEG::HandlerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::PartnerFinder>: 
    public ClassTraitsBase<Herwig::PartnerFinder> {
  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/PartnerFinder"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }

};

}

#include "PartnerFinder.icc"

#endif /* HERWIG_PartnerFinder_H */
