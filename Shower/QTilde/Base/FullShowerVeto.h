// -*- C++ -*-
#ifndef Herwig_FullShowerVeto_H
#define Herwig_FullShowerVeto_H
//
// This is the declaration of the FullShowerVeto class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "FullShowerVeto.fh"
#include "Herwig/Shower/Core/Base/ShowerTree.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the FullShowerVeto class.
 *
 * @see \ref FullShowerVetoInterfaces "The interfaces"
 * defined for FullShowerVeto.
 */
class FullShowerVeto: public Interfaced {

public:

  /**
   * The default constructor.
   */
  FullShowerVeto() : type_(1), behaviour_(0) {}

  /**
   *   Apply the veto
   */
  int applyVeto(ShowerTreePtr); 

  /**
   *  Which type of processes to consider
   */
  unsigned int type() const {return type_;}

  /**
   *  What to do if the event is vetoed
   */
  unsigned int behaviour() const {return behaviour_;}

protected:

  /**
   *  Determine whether to not to veto the shower, to be implemented in inheriting classes
   */
  virtual bool vetoShower() = 0;

  /**
   *  Incoming particles to the hard process
   */
  const vector<tPPtr> & incoming() {return incoming_;}

  /**
   *  Outgoing particles from the hard process
   */
  const vector<tPPtr> & outgoing() {return outgoing_;}

  /**
   *  The final-state particles at the end of the shower
   */
  const vector<tPPtr> & finalState();


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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FullShowerVeto & operator=(const FullShowerVeto &);

private:

  /**
   *  Switches
   */
  //@{
  /**
   *  Which type of processes to consider
   */
  unsigned int type_;

  /**
   *  What to do if the event is vetoed
   */
  unsigned int behaviour_;
  //}

  /**
   *  Temporary storage
   */
  //@{
  /**
   *  Incoming to hard process
   */
  vector<tPPtr> incoming_;

  /**
   *  Outgoing from the hard process
   */
  vector<tPPtr> outgoing_;

  /**
   *  Final State particles
   */
  vector<tPPtr> finalState_;
  //@}
};

}

#endif /* Herwig_FullShowerVeto_H */
