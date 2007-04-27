#ifndef _AMEGIC_INTERFACE_H_
#define _AMEGIC_INTERFACE_H_

// This is the declaration of the AmegicInterface class.

#include "ThePEG/Interface/Interfaced.h"
#include "Amegic.H"
#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Config/Herwig.h"
#include "ThePEG/Repository/Strategy.fh"
#include <fstream>

namespace Herwig {

using namespace ThePEG;
using namespace AMEGIC;

/** \ingroup Interfaces
 * 
 *  AmegicInterface is a class used to interface the functions in 
 *  Amegic with Herwig++. In particular, the total cross section and 
 *  momentum generation functions are interfaced.
 *
 */
class AmegicInterface : public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline AmegicInterface();

  /**
   * The copy constructor.
   */
  inline AmegicInterface(const AmegicInterface &);

  /**
   * The destructor.
   */
  ~AmegicInterface();
  //@}

public:

  /**
   *  Particles from an event
   * @param anti ??
   */
  ParticleVector OneEvent(bool anti = false);

  /**
   *  Incoming particles
   */
  inline PDVector parents();

  /**
   *  Outgoing particles 
   */
  inline PDVector children();

  /**
   *  Convert particle between Amegic and Herwig++
   */
  static Flavour convertParticle(PDPtr);

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Change all pointers to Interfaced objects corresponding clones.
   */
  inline virtual void rebind(const TranslationMap &trans) throw(RebindException);
 
  /**
   * Return pointers to all Interfaced objects referred to by this.
   */
  inline virtual IVector getReferences();

  /**
   * Read setup info from the standard stream.
   */
  virtual void readSetup(istream &is) throw(SetupException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<AmegicInterface> initAmegicInterface;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AmegicInterface & operator=(const AmegicInterface &);

private:

  /**
   *  Initialize Amegic
   */
  void initializeAmegic();

  /**
   *  Set the processs
   */
  void setProcess(string);

  /**
   *  Create the particles
   */
  PDVector createParticles(string);

  /**
   *  the Amegic process
   */
  Amegic *process;

  /**
   *  Name of the setup directory
   */ 
  string setupDirectory;

  /**
   *  String for the process
   */
  string processString;

  /**
   *  Incoming particles
   */
  PDVector inParticles;

  /**
   *  Outgoing particles
   */
  PDVector outParticles;

  /**
   *  Whether initialized or not
   */
  bool initialized;

  /**
   *  Cross section
   */
  double crossSection;

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AmegicInterface. */
template<>
struct BaseClassTrait<Herwig::AmegicInterface,1> {
  /** Typedef of the first base class of AmegicInterface. */
   typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AmegicInterface class and the shared object where it is defined. */
template<>
struct ClassTraits<Herwig::AmegicInterface> 
  : public ClassTraitsBase<Herwig::AmegicInterface> {
  /** Return a platform-independent class name */
   static string className() { return "Herwig++::AmegicInterface"; }
  /** Return the name of the shared library be loaded to get
   *  access to the AmegicInterface class and every other class it uses
   *  (except the base class). */
   static string libName() { return "libHwInterfaces.so"; }
};

/** @endcond */

}

#include "AmegicInterface.icc"

#endif // _AMEGIC_INTERFACE_H_
