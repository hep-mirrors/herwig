// -*- C++ -*-
#ifndef HERWIG_BaryonWidthGenerator_H
#define HERWIG_BaryonWidthGenerator_H
//
// This is the declaration of the BaryonWidthGenerator class.
//

#include "GenericWidthGenerator.h"
#include "BaryonWidthGenerator.fh"
#include "Herwig++/Decay/Baryon/Baryon1MesonDecayerBase.h"

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

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline BaryonWidthGenerator();

  /**
   * The copy constructor.
   */
  inline BaryonWidthGenerator(const BaryonWidthGenerator &);

  /**
   * The destructor.
   */
  virtual ~BaryonWidthGenerator();
  //@}

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
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
  /**
   *  First coupling for the A terms
   */
  vector<double> _Afact1;

  /**
   *  Second coupling for the A terms
   */
  vector<double> _Afact2;

  /**
   *  Third coupling for the A terms
   */
  vector<double> _Afact3;

  /**
   *  Fourth coupling for the A terms
   */
  vector<double> _Afact4;

  /**
   *  Fifth coupling for the A terms
   */
  vector<double> _Afact5;

  /**
   *  Sixth coupling for the A terms
   */
  vector<double> _Afact6;

  /**
   *  First coupling for the B terms
   */
  vector<double> _Bfact1;

  /**
   *  Second coupling for the B terms
   */
  vector<double> _Bfact2;

  /**
   *  Third coupling for the B terms
   */
  vector<double> _Bfact3;

  /**
   *  Fourth coupling for the B terms
   */
  vector<double> _Bfact4;

  /**
   *  Fifth coupling for the B terms
   */
  vector<double> _Bfact5;

  /**
   *  Sixth coupling for the B terms
   */
  vector<double> _Bfact6;

};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

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
  static string className() { return "Herwig++::BaryonWidthGenerator"; }
  /** Return the name of the shared library be loaded to get
   *  access to the BaryonWidthGenerator class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwPDT.so"; }
};

}

#include "BaryonWidthGenerator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BaryonWidthGenerator.tcc"
#endif

#endif /* HERWIG_BaryonWidthGenerator_H */
