// -*- C++ -*-
#ifndef HERWIG_BaryonWidthGenerator_H
#define HERWIG_BaryonWidthGenerator_H
//
// This is the declaration of the BaryonWidthGenerator class.
//

#include "GenericWidthGenerator.h"
#include "BaryonWidthGenerator.fh"
#include "Herwig/Decay/Baryon/Baryon1MesonDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 * The BaryonWidthGenerator class is designed to automatically calculate the running
 * width for a given particle using information from the decayModes and the 
 * Baryon1MesonDecayer to construct the running width.
 *
 * It inherits from the GenericWidthGenerator.
 *
 * @see GenericWidthGenerator
 * @see GenericMassGenerator
 *
 */
class BaryonWidthGenerator: public GenericWidthGenerator {

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

public:

  /**
   * Output the initialisation info for the database
   * @param output The stream to output the information to
   * @param header output the header.
   **/
  virtual void dataBaseOutput(ofstream & output,bool header=true);


  /**
   * The \f$1\to2\f$ width for outgoing particles which can be off-shell.
   * @param iloc The location of the mode in the list.
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the first outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @return The partial width.
   */
  virtual Energy partial2BodyWidth(int iloc,Energy m0,Energy m1,Energy m2) const;

protected: 

  /**
   * Perform the set up for a mode, this is called by the base class
   * @param mode The decay mode
   * @param decayer The decayer for the mode.
   * @param imode The number of the mode.
   */
  virtual void setupMode(tcDMPtr mode, tDecayIntegratorPtr decayer, unsigned int imode);

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<BaryonWidthGenerator> initBaryonWidthGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BaryonWidthGenerator & operator=(const BaryonWidthGenerator &);

private:

  /**
   *  vector of pointers to the decayers 
   */
  vector<Baryon1MesonDecayerBasePtr> _baryondecayers;

  /**
   * Location of the decay mode in the decayer
   */
  vector<int> _modeloc;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BaryonWidthGenerator. */
template <>
struct BaseClassTrait<Herwig::BaryonWidthGenerator,1> {
  /** Typedef of the first base class of BaryonWidthGenerator. */
  typedef Herwig::GenericWidthGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BaryonWidthGenerator class and the shared object where it is defined. */
template <>
 struct ClassTraits<Herwig::BaryonWidthGenerator>
  : public ClassTraitsBase<Herwig::BaryonWidthGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BaryonWidthGenerator"; }
  /** Return the name of the shared library be loaded to get
   *  access to the BaryonWidthGenerator class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwBaryonDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_BaryonWidthGenerator_H */
