#ifndef _KTJET_INTERFACE_H_
#define _KTJET_INTERFACE_H_

//
// This is the declaration of the <!id>KtJetInterface<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// <!id>AmegicInterface<!!id> is a class used to interface the
// functions in Amegic with Herwig++. In particular, the total
// cross section and momentum generation functions are interfaced.
//
// CLASSDOC SUBSECTION See also:
//

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

class KtJetInterface {
 public:
  KtJetInterface() {}
  ~KtJetInterface() {}
  // Standard ctors and dtor

 public:
  vector<KtLorentzVector> convertToKtVectorList(tPVector &); 
  int getThePEGID(KtLorentzVector &);
 private:
  KtLorentzVector convertToKtVector(PPtr);
  std::map<int, int> Kt2PythiaMap;
};

}

#endif // _KTJET_INTERFACE_H_
