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
using namespace KtJet;

/** \ingroup Interfaces
 * 
 *  Some comment should be provided!
 */
class KtJetInterface {

 public:

  KtJetInterface() {}
  ~KtJetInterface() {}

  inline void clearMap() { Kt2PythiaMap.clear(); }

 public:

  vector<KtLorentzVector> convertToKtVectorList(const tPVector &); 
  int getThePEGID(KtLorentzVector &);

 private:

  KtLorentzVector convertToKtVector(const PPtr &);
  std::map<int, int> Kt2PythiaMap;

};

}

#endif // _KTJET_INTERFACE_H_
