// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptonDalitzAnalysis class.
//

#include "LeptonDalitzAnalysis.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LeptonDalitzAnalysis.tcc"
#endif


using namespace Herwig;

LeptonDalitzAnalysis::~LeptonDalitzAnalysis() {}

void LeptonDalitzAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  if(_nout[0]<50000)
    {
      while(true)
	{
	  tPVector pert=event->primaryCollision()->step(1)->getFinalState();
	  tPPtr quark,anti,gluon;
	  for(unsigned int ix=0;ix<pert.size();++ix)
	    {
	      if(pert[ix]->id()>-6&&pert[ix]->id()<0)
		{
		  if(!anti) anti=pert[ix];
		  else break;
		}
	      else if(pert[ix]->id()>0&&pert[ix]->id()<6)
		{
		  if(!quark) quark=pert[ix];
		  else break;
		}
	      else if(pert[ix]->id()==21)
		{
		  if(!gluon) gluon=pert[ix];
		  else break;
		}
	    }
	  if(!quark||!anti||!gluon) break;
	  Energy total=quark->momentum().e()+anti->momentum().e()+gluon->momentum().e();
	  double x1=quark->momentum().e()/total*2.;
	  double x2=anti->momentum().e()/total*2.;
	  ++_nout[0];
	  _output[0] << x1 << " " << x2 << endl;
	  break;
	};
    }
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
}

LorentzRotation LeptonDalitzAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void LeptonDalitzAnalysis::analyze(const tPVector & particles) {
  if(_nout[1]>50000) return;
  _kint.clearMap();
  KtJet::KtEvent ev = KtJet::KtEvent(_kint.convertToKtVectorList(particles), 1, 1, 1);
  ev.findJetsN(3);
  vector<KtJet::KtLorentzVector> ktjets;
  vector<Lorentz5Momentum> jets;
  ktjets = ev.getJetsPt();
  if (ktjets.size() == 3) 
    {
      for (int j=0; j<3; j++){jets.push_back(ktjets[j]);}
      double x1=0.,x2=0.,xtemp;
      Energy total=jets[0].e()+jets[1].e()+jets[2].e();
      for(unsigned int ix=0;ix<3;++ix)
	{
	  xtemp=2.*jets[ix].rho()/total;
	  if(xtemp>x1){x2=x1;x1=xtemp;}
	  else if(xtemp>x2){x2=xtemp;}
	}
      //x1=2.*jets[0].e()/total;
      //x2=2.*jets[1].e()/total;
      //if(UseRandom::rndbool()) x2=2.-x1-x2;
      if(UseRandom::rndbool())
	_output[1] << x1 << " " << x2 << "\n";
      else
	_output[1] << x2 << " " << x1 << "\n";
      ++_nout[1];
    }
}

void LeptonDalitzAnalysis::analyze(tPPtr) {}

NoPIOClassDescription<LeptonDalitzAnalysis> LeptonDalitzAnalysis::initLeptonDalitzAnalysis;
// Definition of the static class description member.

void LeptonDalitzAnalysis::Init() {

  static ClassDocumentation<LeptonDalitzAnalysis> documentation
    ("There is no documentation for the LeptonDalitzAnalysis class");

}

