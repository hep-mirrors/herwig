// -*- C++ -*-
#ifndef Herwig_HeavyMesonWidthGenerator_H
#define Herwig_HeavyMesonWidthGenerator_H
//
// This is the declaration of the HeavyMesonWidthGenerator class.
//

#include "GenericWidthGenerator.h"
#include "Herwig/Decay/HeavyMeson/HQETStrongDecayer.h"
#include "Herwig/Decay/HeavyMeson/HQETRadiativeDecayer.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The HeavyMesonWidthGenerator class calculates the running width
 * for heavy mesons.
 *
 * @see \ref HeavyMesonWidthGeneratorInterfaces "The interfaces"
 * defined for HeavyMesonWidthGenerator.
 */
class HeavyMesonWidthGenerator: public GenericWidthGenerator {

public:

  /**
   * The default constructor.
   */
  HeavyMesonWidthGenerator();

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HeavyMesonWidthGenerator & operator=(const HeavyMesonWidthGenerator &) = delete;

private:

  /**
   *  Coupings for the decays
   */
  //@{
  /**
   *  Pion decay constant
   */
  Energy fPi_;

  /**
   *  Coupling for decays within the \f$(0^-,1^-)\f$ multiplet
   */
  double g_;

  /**
   *  Coupling for decays within the \f$(0^+,1^+)\f$ multiplet
   */
  double gp_;

  /**
   *  Coupling for decays from the \f$(0^+ ,1^+)\f$ multiplet
   */
  double h_;

  /**
   *  Coupling for decays from the \f$(1^+ ,2^+)\f$ multiplet
   */
  double hp_;

  /**
   *  Coupling for decays from the \f$(2^- ,3^-)\f$ multiplet
   */
  double k_;

  /**
   *  Coupling for decays from the \f$(1^- ,2^-)\f$ multiplet
   */
  double kp_;

  /**
   *  Coupling for decays from the 2S \f$(0^-,1^-)\f$ multiplet
   */
  double gtilde_;
  
  /**
   *  D_1 mixing angle (up and down)
   */
  double psiL_;

  /**
   *  D_1 mixing angle (strange)
   */
  double psiS_;
  //@}

  /**
   *   A momentum scale characterising the convergence of the
   *   derivative expansion. We expect Lambda_ ~ 1 GeV.
   */
  Energy Lambda_;
  //@}

  /**
   *  Check if couplings set
   */
  bool couplingsSet_;
};

}

#endif /* Herwig_HeavyMesonWidthGenerator_H */
