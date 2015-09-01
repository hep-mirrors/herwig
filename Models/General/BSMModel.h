// -*- C++ -*-
#ifndef Herwig_BSMModel_H
#define Herwig_BSMModel_H
//
// This is the declaration of the BSMModel class.
//

#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Utilities/CFileLineReader.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the BSMModel class.
 *
 * @see \ref BSMModelInterfaces "The interfaces"
 * defined for BSMModel.
 */
class BSMModel: public Herwig::StandardModel {

public:

  /**
   * The default constructor.
   */
  BSMModel();

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
   * Get name of SLHA decay file
   */
  const string & decayFile() const {return decayFile_;}

  /**
   * Set name of SLHA decay file
   */
  void decayFile(string in) {decayFile_ = in;}

  /**
   *  Read the decays
   */
  void decayRead();

  /**
   * Read decaymodes from LHA file
   * @param ifs input stream containg data
   * @param decay string containing name of parent and value of total width
   */
  void readDecay(CFileLineReader & ifs, string decay) const;

  /**
   * Create a DecayMode object in the repository
   * @param tag The tag identifying the decay mode including the prefix
   * 'decaymode'
   * @param brat Branching ratio of this mode 
   */
  void createDecayMode(string tag, double brat) const;

  /**
   * Create a DecayMode object in the repository
   * @param tag The tag identifying the decay mode including the prefix
   * 'decaymode'
   * @param brat Branching ratio of this mode 
   */
  vector<pair<double,string> > createWZDecayModes(string tag, double brat,
						  tcPDPtr boson,
						  Energy maxMass) const;

  /**
   *  read the decays
   */
  bool readDecays() const {return readDecays_;}

  /**
   *  set the reading of the decays
   */ 
  void readDecays(bool in) {readDecays_=in;}

  /**
   *  Map of PDG ids from file to those used internally
   */
  map<long,long> & idMap() {return idMap_;}

  /**
   *  Get ParticleData object with Id mapping
   */
  PDPtr getBSMParticleData(PID id) const {
    map<long,long>::const_iterator it = idMap_.find(id);
    if(it==idMap_.end())
      return getParticleData(id);
    else
      return getParticleData(it->second);
  }

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

  /**
   * Overloaded function from Interfaced
   */
  virtual bool preInitialize() const {
    return true;
  }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BSMModel & operator=(const BSMModel &);

private:

  /**
   *  Name of the decay file
   */
  string decayFile_;

  /**
   *  Read the decays from the file
   */
  bool readDecays_;

  /**
   * Whether or not to replace the top decay modes with those from
   * the SLHA files
   */
  bool topModesFromFile_;

  /**
   *  Tolerance for branching ratios
   */
  double tolerance_;

  /**
   *  Map of ids from files to those used by Herwig
   */
  map<long,long> idMap_;

};

}

#endif /* Herwig_BSMModel_H */
