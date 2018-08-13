// -*- C++ -*-
#ifndef Herwig_MEMultiChannel_H
#define Herwig_MEMultiChannel_H
//
// This is the declaration of the MEMultiChannel class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig/Decay/PhaseSpaceMode.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEMultiChannel class.
 *
 * @see \ref MEMultiChannelInterfaces "The interfaces"
 * defined for MEMultiChannel.
 */
class MEMultiChannel: public MEBase {

public:

  /**
   * The default constructor.
   */
  MEMultiChannel() {};

  /**
   * The destructor.
   */
  virtual ~MEMultiChannel();
  //@}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  double me2() const {return me2(-1);}

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  //@}

protected :

  /**
   *   Add a phase-space mode
   */
  void addMode(PhaseSpaceModePtr mode) {
    modes_.push_back(mode);
  }

  /**
   *  Access to the modes
   */
  tPhaseSpaceModePtr mode(unsigned int imode) const {
    return modes_[imode];
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

protected:

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan) const = 0;

  /**
   *   Access the mode currently being used
   */
  unsigned int iMode() {return iMode_;}

  /**
   *  Set the mode currently begining used
   */
  void iMode(unsigned int in) const {iMode_=in;}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEMultiChannel & operator=(const MEMultiChannel &);

private:

  /**
   *   The phase-space modes
   */
  vector<PhaseSpaceModePtr> modes_;

  /**
   *   Map from the phase space channels to the modes
   */
  mutable vector<map <int,int> > channelMap_;

  /**
   *   The mode currently begining used
   */
  mutable unsigned int iMode_;
};

}

#endif /* Herwig_MEMultiChannel_H */
