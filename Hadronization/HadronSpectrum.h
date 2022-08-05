// -*- C++ -*-
#ifndef Herwig_HadronSpectrum_H
#define Herwig_HadronSpectrum_H
//
// This is the declaration of the HadronSpectrum class.
//

#include "ThePEG/Interface/Interfaced.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HadronSpectrum class.
 *
 * @see \ref HadronSpectrumInterfaces "The interfaces"
 * defined for HadronSpectrum.
 */
class HadronSpectrum: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HadronSpectrum();

  /**
   * The destructor.
   */
  virtual ~HadronSpectrum();
  //@}

public:

  /** @name Partonic content */
  //@{

  /**
   * Return the id of the gluon
   */
  virtual long gluonId() const = 0;

  /**
   * Return the ids of all hadronizing quarks
   */
  virtual const vector<long>& hadronizingQuarks() const = 0;

  /**
   * The light hadronizing quarks
   */
  virtual const vector<long>& lightHadronizingQuarks() const = 0;

  /**
   * The heavy hadronizing quarks
   */
  virtual const vector<long>& heavyHadronizingQuarks() const = 0;

  /**
   * Return true if any of the possible three input particles contains
   * the indicated heavy quark.  false otherwise. In the case that
   * only the first particle is specified, it can be: an (anti-)quark,
   * an (anti-)diquark an (anti-)meson, an (anti-)baryon; in the other
   * cases, each pointer is assumed to be either (anti-)quark or
   * (anti-)diquark.
   */
  virtual bool hasHeavy(long id, tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr()) const = 0;

  /**
   * Return true, if any of the possible input particle pointer is an
   * exotic quark, e.g. Susy quark; false otherwise.
   */
  virtual bool isExotic(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr()) const = 0;

  //@}

  /**
   * Return true if the two or three particles in input can be the components 
   * of a hadron; false otherwise.
   */
  virtual bool canBeHadron(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3 = PDPtr()) const = 0;

    /**
   *  Check if can't make a diquark from the partons
   */
  virtual bool canMakeDiQuark(tcPPtr p1, tcPPtr p2) const = 0;

  /**
   * Return the particle data of the diquark (anti-diquark) made by the two 
   * quarks (antiquarks) par1, par2.
   * @param par1 (anti-)quark data pointer
   * @param par2 (anti-)quark data pointer
   */
  virtual PDPtr makeDiquark(tcPDPtr par1, tcPDPtr par2) const = 0;

  /**
   * Return the quark flavour which should be considered to set the
   * minimum mass for a minimal cluster splitting.
   */
  virtual long minimalSplitQuark() const = 0;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HadronSpectrum & operator=(const HadronSpectrum &) = delete;

};

}

#endif /* Herwig_HadronSpectrum_H */
