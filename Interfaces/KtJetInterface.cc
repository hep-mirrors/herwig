#include "Herwig++/Interfaces/KtJetInterface.h"
#include "Pythia7/EventRecord/Particle.h"

using namespace Herwig;
using namespace Pythia7;

vector<KtJet::KtLorentzVector> KtJetInterface::convertToKtVectorList(tPVector &pv) {
  vector<KtLorentzVector> rval;
  for(tPVector::iterator it = pv.begin(); it != pv.end(); it++) {
    rval.push_back(KtJetInterface::convertToKtVector(*it));
    Kt2PythiaMap[rval.back().getID()] = (*it)->number();
  }
  return rval;
}
   
KtLorentzVector KtJetInterface::convertToKtVector(PPtr p) {
  return KtJet::KtLorentzVector(p->momentum());
}

int KtJetInterface::getPythia7ID(KtJet::KtLorentzVector &kv) {
  return Kt2PythiaMap[kv.getID()];
}
