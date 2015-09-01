// -*- C++ -*-
#ifndef HERWIG_AcerDet_H
#define HERWIG_AcerDet_H
//
// This is the declaration of the AcerDet class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_HEPEVT.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The AcerDet class is designed to interface Herwig with the
 * AcerDet fast detector simulation
 *
 * @see \ref AcerDetInterfaces "The interfaces"
 * defined for AcerDet.
 */
class AcerDet: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  inline AcerDet() {}

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);
  //@}

public:

  /**
   *  Members to analyse the information from AcerDet
   */
  //@{
  /**
   *  Number of photons
   */
  inline unsigned int numberOfPhotons() const {return _nphoton;}

  /**
   *  Momenta of the photons
   */
  inline const vector<LorentzMomentum> & photonMomentum() const {return _photonMomentum;}

  /**
   *  Number of leptons
   */
  inline unsigned int numberOfLeptons() const {return _nlepton;}

  /**
   *  Momenta of the leptons
   */
  inline const vector<LorentzMomentum> & leptonMomentum() const {return _leptonMomentum;}

  /**
   *  PDG codes for the leptons
   */
  inline const vector<int> & leptonID() const {return _leptonID;}

  /**
   *  Number of jets
   */
  inline unsigned int numberOfJets() const {return _njet;}

  /**
   *  Momenta of the jets
   */
  inline const vector<LorentzMomentum> & jetMomentum() const {return _jetMomentum;}

  /**
   * PDG codes for the jets (0 for light 5 for bottom)
   * probably not very relible
   */
  inline const vector<int> & jetID() const {return _jetID;}

  /**
   *  The missing \f$E_T\f$ from the calorimeter
   * @return \f$p_x\f$, \f$p_y\f$
   */
  inline const pair<Energy,Energy> & missingETCalorimeter() const {_etcalo;}

  /**
   *  The missing \f$E_T\f$ from the neutrinos
   * @return \f$p_x\f$, \f$p_y\f$
   */
  inline const pair<Energy,Energy> & missingETNeutrino() const {_etneutrino;}

  /**
   *  The missing \f$E_T\f$ from the user identified stable neutral particles
   * @return \f$p_x\f$, \f$p_y\f$
   */
  inline const pair<Energy,Energy> & missingETStable() const {_etstable;}
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<AcerDet> initAcerDet;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AcerDet & operator=(const AcerDet &);

private:

  /**
   *  Converter from HepMC to HEPEVT
   */
  HepMC::IO_HEPEVT * _converter;

  /**
   *  Number of photons
   */
  unsigned int _nphoton;

  /**
   *  Momenta of the photons
   */
  vector<LorentzMomentum> _photonMomentum;

  /**
   *  Number of leptons
   */
  unsigned int _nlepton;

  /**
   *  Momenta of the leptons
   */
  vector<LorentzMomentum> _leptonMomentum;

  /**
   * PDG codes for the leptons
   */
  vector<int> _leptonID;

  /**
   *  Number of jets
   */
  unsigned int _njet;

  /**
   *  Momenta of the jets
   */
  vector<LorentzMomentum> _jetMomentum;

  /**
   * PDG codes for the jets
   */
  vector<int> _jetID;

  /**
   *  The missing \f$E_T\f$ from the calorimeter
   */
  pair<Energy,Energy> _etcalo;

  /**
   *  The missing \f$E_T\f$ from the neutrinos
   */
  pair<Energy,Energy> _etneutrino;

  /**
   *  The missing \f$E_T\f$ from the user identified stable neutral particles
   */
  pair<Energy,Energy> _etstable;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AcerDet. */
template <>
struct BaseClassTrait<Herwig::AcerDet,1> {
  /** Typedef of the first base class of AcerDet. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AcerDet class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::AcerDet>
  : public ClassTraitsBase<Herwig::AcerDet> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::AcerDet"; }
  /**
   * The name of a file containing the dynamic library where the class
   * AcerDet is implemented. It may also include several, space-separated,
   * libraries if the class AcerDet depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwAcerDet.so"; }
};

/** @endcond */

}

#endif /* HERWIG_AcerDet_H */
