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
#include "Pythia7/Config/Pythia7.h"
#include "Herwig++/Config/Herwig.h"
#include "Pythia7/Repository/Strategy.fh"
#include <fstream>
#include <map>

namespace Herwig {

using namespace Pythia7;
using namespace KtJet;

class KtJetInterface {
 public:
  KtJetInterface() {}
  ~KtJetInterface() {}
  // Standard ctors and dtor

 public:
  vector<KtLorentzVector> convertToKtVectorList(tPVector &); 
  int getPythia7ID(KtLorentzVector &);
 private:
  KtLorentzVector convertToKtVector(PPtr);
  std::map<int, int> Kt2PythiaMap;
};

}

#endif // _KTJET_INTERFACE_H_
