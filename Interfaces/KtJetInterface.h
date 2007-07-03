#ifndef _KTJET_INTERFACE_H_
#define _KTJET_INTERFACE_H_

// This is the declaration of the KtJetInterface class.

#include "KtJet/KtLorentzVector.h"
#include "KtJet/KtEvent.h"
#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Config/Herwig.h"
#include "ThePEG/Repository/Strategy.fh"
#include <fstream>
#include <map>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Interfaces
 * 
 *  Interface to allow the KtJET jet clustering library to be used.
 */
class KtJetInterface {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  KtJetInterface() {}

  /**
   * The copy constructor.
   */
  ~KtJetInterface() {}

public:

  /**
   * Clear map between ThePEG particles and KtJet vectors
   */
  inline void clearMap() { Kt2PythiaMap.clear(); }

public:

  /**
   *  Convert ThePEG particles to KtJet vectors
   */
  vector<KtJet::KtLorentzVector> convert(const tPVector &);

  /**
   *  Convert  KtJet Lorentz vectors to ThePEG momenta
   */
  static vector<LorentzMomentum> convert(const vector<KtJet::KtLorentzVector> &);

  /**
   *  Get the PDG code for a KtJet vector
   */ 
  int getThePEGID(KtJet::KtLorentzVector &);

  /**
   * Convert KtJet vector back to ThePEG
   */
  static LorentzMomentum convert(const KtJet::KtLorentzVector & kt);

 private:

  /**
   *  Convert one particle to KtJet
   */
  static KtJet::KtLorentzVector convert(tcPPtr);

  /**
   *  Map between Herwig++ and KtJet
   */
  std::map<int, int> Kt2PythiaMap;

};

}

#endif // _KTJET_INTERFACE_H_
