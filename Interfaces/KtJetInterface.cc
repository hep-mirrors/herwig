#include "Herwig++/Interfaces/KtJetInterface.h"
#include "ThePEG/EventRecord/Particle.h"

using namespace Herwig;
using namespace ThePEG;

vector<KtJet::KtLorentzVector> KtJetInterface::convertToKtVectorList(const tPVector &pv) {
  vector<KtLorentzVector> rval;
  for(tPVector::const_iterator it = pv.begin(); it != pv.end(); it++) {
    rval.push_back(KtJet::KtLorentzVector((*it)->momentum()));
    Kt2PythiaMap[rval.back().getID()] = (*it)->number();
  }
  return rval;
}
   
KtLorentzVector KtJetInterface::convertToKtVector(const PPtr &p) {
  return KtJet::KtLorentzVector(p->momentum());
}

int KtJetInterface::getThePEGID(KtJet::KtLorentzVector &kv) {
  return Kt2PythiaMap[kv.getID()];
}
