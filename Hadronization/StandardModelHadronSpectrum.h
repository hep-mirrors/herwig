// -*- C++ -*-
#ifndef Herwig_StandardModelHadronSpectrum_H
#define Herwig_StandardModelHadronSpectrum_H
//
// This is the declaration of the StandardModelHadronSpectrum class.
//

#include "Herwig/Hadronization/HadronSpectrum.h"
#include "CheckId.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the StandardModelHadronSpectrum class.
 *
 * @see \ref StandardModelHadronSpectrumInterfaces "The interfaces"
 * defined for StandardModelHadronSpectrum.
 */
class StandardModelHadronSpectrum: public HadronSpectrum {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  StandardModelHadronSpectrum();

  /**
   * The destructor.
   */
  virtual ~StandardModelHadronSpectrum();
  //@}

public:

  /** @name Partonic content */
  //@{

  /**
   * Return the id of the gluon
   */
  virtual long gluonId() const { return ParticleID::g; }

  /**
   * Return the ids of all hadronizing quarks
   */
  virtual const vector<long>& hadronizingQuarks() const {
    static vector<long> hadronizing =
      { ParticleID::d, ParticleID::u, ParticleID::s, ParticleID::c, ParticleID::b };
    return hadronizing;
  }

  /**
   * The light hadronizing quarks
   */
  virtual const vector<long>& lightHadronizingQuarks() const {
    static vector<long> light =
      { ParticleID::d, ParticleID::u, ParticleID::s };
    return light;
  }

  /**
   * The heavy hadronizing quarks
   */
  virtual const vector<long>& heavyHadronizingQuarks() const {
    static vector<long> heavy =
      { ParticleID::c, ParticleID::b };
    return heavy;
  }

  /**
   * Return true if any of the possible three input particles contains
   * the indicated heavy quark.  false otherwise. In the case that
   * only the first particle is specified, it can be: an (anti-)quark,
   * an (anti-)diquark an (anti-)meson, an (anti-)baryon; in the other
   * cases, each pointer is assumed to be either (anti-)quark or
   * (anti-)diquark.
   */
  virtual bool hasHeavy(long id, tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr()) const {
    if ( abs(id) == ParticleID::c )
      return CheckId::hasCharm(par1,par2,par3);
    if ( abs(id) == ParticleID::b )
      return CheckId::hasBottom(par1,par2,par3);
    return false;
  }

  /**
   * Return true, if any of the possible input particle pointer is an
   * exotic quark, e.g. Susy quark; false otherwise.
   */
  virtual bool isExotic(tcPDPtr par1, tcPDPtr par2 = PDPtr(), tcPDPtr par3 = PDPtr()) const {
    return CheckId::isExotic(par1,par2,par3);
  }

  //@}

  /**
   * Return true if the two or three particles in input can be the components 
   * of a hadron; false otherwise.
   */
  virtual bool canBeHadron(tcPDPtr par1, tcPDPtr par2 , tcPDPtr par3 = PDPtr()) const {
    return CheckId::canBeHadron(par1,par2,par3);
  }

  /**
   * Check if can't make a diquark from the partons
   */
  virtual bool canMakeDiQuark(tcPPtr p1, tcPPtr p2) const {
    long id1 = p1->id(), id2 = p2->id();
    return QuarkMatcher::Check(id1) && QuarkMatcher::Check(id2) && id1*id2>0;
  }

  /**
   * Return the particle data of the diquark (anti-diquark) made by the two 
   * quarks (antiquarks) par1, par2.
   * @param par1 (anti-)quark data pointer
   * @param par2 (anti-)quark data pointer
   */
  virtual PDPtr makeDiquark(tcPDPtr par1, tcPDPtr par2) const {
    return CheckId::makeDiquark(par1,par2);
  }

  /**
   * Return the quark flavour which should be considered to set the
   * minimum mass for a minimal cluster splitting.
   */
  virtual long minimalSplitQuark() const {
    //cout << "ParticleID::d= " << ParticleID::d << endl;
    //cout << getParticleData(ParticleID::d)->id() << endl;
    return ParticleID::d;
  }

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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  StandardModelHadronSpectrum & operator=(const StandardModelHadronSpectrum &)  = delete;

};

}

#endif /* Herwig_StandardModelHadronSpectrum_H */
