#include "Herwig++/Interfaces/KtJetInterface.h"
#include "Herwig++/Utilities/HerwigRun.h"

int main(int argc, char * argv[]) {
  using namespace Pythia7;
  tPVector fs;
  //KtJet::KtEvent ev;
  Herwig::KtJetInterface kint;

  /*** set KtEvent flags */
  int type = 0; // e+e-
  int angle = 1; // Angular
  int recom = 2; // E
  double ycut = 0.03;

  try {
    Herwig::HerwigRun hw(argc,argv);

    if(hw.isRunMode() && hw.preparedToRun()) {
      for(int i = 0; i<hw.getN(); i++) {
        fs = hw.getFinalState(1);
        KtJet::KtEvent ev = KtJet::KtEvent(kint.convertToKtVectorList(fs), type, angle, recom);
        ev.findJetsY(ycut);
        cout << "KtJet says there are " << ev.getNJets() << endl;
      }
    } else if(hw.isRunMode()) {
      std::cout << "Error: Expected a run but there is no EventGenerator object!\n" << endl;
    }
  }
  catch ( std::exception & e ) {
    cerr << e.what() << endl;
    return 1;
  }
  catch ( ... ) {
    cerr << "Unknown Exception\n";
    return 2;
  }

  return 0;
}

