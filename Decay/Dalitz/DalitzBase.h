// -*- C++ -*-
#ifndef Herwig_DalitzBase_H
#define Herwig_DalitzBase_H
//
// This is the declaration of the DalitzBase class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "DalitzResonance.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The DalitzBase class provides a base class for the implementation of three-body Dalitz decays.
 *
 * @see \ref DalitzBaseInterfaces "The interfaces"
 * defined for DalitzBase.
 */
class DalitzBase: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  DalitzBase() : rParent_(5./GeV), useAllK0_(false),
		 maxWgt_(1.), channel1_(-1), channel2_(-1),
		 incoming_(0), outgoing_({0,0,0}) {
    // intermediates
    generateIntermediates(true);
  }

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

  /**
   *   Set the parameters for a decay mode
   */
  string addChannel(string arg);

  /**
   *   Set the parameters for a decay mode
   */
  string setExternal(string arg);

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

  /**
   *  Add a new resonance
   */
  void addResonance(const DalitzResonancePtr & R) {resonances_.push_back(R);}

  /**
   *  Set up the phase-space
   */
  void createMode(tPDPtr in, tPDVector out);

  /**
   *  Access to the resonances
   */
  const vector<DalitzResonancePtr> & resonances() const {return resonances_;}

  /**
   *   Access to the first channel
   */
  const int & channel1() const {
    return channel1_;
  }
    
  /**
   *   Access to the second channel
   */
  const int & channel2() const {
    return channel2_;
  }

  /**
   *  Radius of the parent
   */
  const InvEnergy & parentRadius() const {
    return rParent_;
  }

  /**
   *  Weights for the phase-space channels
   */
  const vector<double> & weights() const {return weights_;}
      
protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DalitzBase & operator=(const DalitzBase &) = delete;

private:

  /**
   *   The radii for the Blatt-Weisskopf form factors
   */
  //@{
  /**
   *   For the decaying particles
   */
  InvEnergy rParent_;
  //@}

  /**
   *  Take all \f$K_0\f$ mesnos to be the same
   */
  bool useAllK0_;

  /**
   *  Vector containing the intermediate resonances
   */
  vector<DalitzResonancePtr> resonances_;
  
  /**
   *   Parameters for the phase-space sampling
   */
  //@{
  /**
   *  Maximum weight for the decay
   */
  double maxWgt_;

  /**
   *  Weights for the phase-space channels
   */
  vector<double> weights_;
  
private:

  /**
   *  Control over channels to check fit fractions
   */
  int channel1_, channel2_;

  /**
   *  The incoming particle
   */
  long incoming_;

  /**
   *  The outgoing pairtcles
   */
  array<long,3> outgoing_;

};

}

#endif /* Herwig_DalitzBase_H */
