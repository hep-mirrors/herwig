// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TauTo4MesonAnalysis class.
//

#include "TauTo4MesonAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

void TauTo4MesonAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tPVector hadrons=event->getFinalState();
  map<tPPtr,ParticleVector> taus;
  for(unsigned int ix=0;ix<hadrons.size();++ix) {
    PPtr mother=hadrons[ix];
    do {
      if(!mother->parents().empty()) mother=mother->parents()[0];
      else                           mother=tPPtr();
    }
    while(mother&&abs(mother->id())!=ParticleID::tauminus);
    if(mother&&abs(mother->id())==ParticleID::tauminus) {
      if(taus.find(mother)==taus.end()) {
	taus.insert(make_pair(mother,ParticleVector()));
      }
      taus[mother].push_back(hadrons[ix]);
    }
  }
  map<tPPtr,ParticleVector>::const_iterator tit;
  for(tit=taus.begin();tit!=taus.end();++tit) {
    if(tit->second.size()!=5) continue;
    ParticleVector decay=tit->second;
    int tsign=tit->first->id()/abs(tit->first->id());
    vector<Lorentz5Momentum> ppi0,ppip,ppim,pkp,pkm,pk0,pk0bar,peta,pgamma;
    Lorentz5Momentum ptotal;
    for(unsigned int ix=0;ix<decay.size();++ix) {
      long id = decay[ix]->id()*tsign;
      if(abs(id)!=ParticleID::nu_tau) ptotal+=decay[ix]->momentum();
      if(         id ==ParticleID::piplus)  ppip  .push_back(decay[ix]->momentum());
      else if(    id ==ParticleID::piminus) ppim  .push_back(decay[ix]->momentum());
      else if(abs(id)==ParticleID::pi0)     ppi0  .push_back(decay[ix]->momentum());
      else if(    id ==ParticleID::Kplus)   pkp   .push_back(decay[ix]->momentum());
      else if(    id ==ParticleID::Kminus)  pkm   .push_back(decay[ix]->momentum());
      else if(    id ==ParticleID::K0)      pk0   .push_back(decay[ix]->momentum());
      else if(    id ==ParticleID::Kbar0)   pk0bar.push_back(decay[ix]->momentum());
      else if(abs(id)==ParticleID::eta)     peta  .push_back(decay[ix]->momentum());
      else if(abs(id)==ParticleID::gamma)   pgamma.push_back(decay[ix]->momentum());
    }
    if(ppi0.size()==3&&ppim.size()==1) {
      *_mpipi[0]     += (ppi0[0]+ppim[0]).m()/GeV;
      *_mpipi[0]     += (ppi0[1]+ppim[0]).m()/GeV;
      *_mpipi[0]     += (ppi0[2]+ppim[0]).m()/GeV;
      *_mpipi[1]     += (ppi0[0]+ppi0[1]).m()/GeV;
      *_mpipi[1]     += (ppi0[0]+ppi0[2]).m()/GeV;
      *_mpipi[1]     += (ppi0[1]+ppi0[2]).m()/GeV;
      *_mpipipi[0]   += (ppi0[0]+ppi0[1]+ppi0[2]).m()/GeV;
      *_mpipipi[1]   += (ppi0[0]+ppi0[1]+ppim[0]).m()/GeV;
      *_mpipipi[1]   += (ppi0[0]+ppi0[2]+ppim[0]).m()/GeV;
      *_mpipipi[1]   += (ppi0[1]+ppi0[2]+ppim[0]).m()/GeV;
      *_mpipipipi[0] += (ppi0[0]+ppi0[1]+ppi0[2]+ppim[0]).m()/GeV;
    }
    else if(ppi0.size()==1&&ppip.size()==1&&ppim.size()==2) {
      *_mpipi[2] +=(ppi0[0]+ppip[0]).m()/GeV;
      *_mpipi[3] +=(ppi0[0]+ppim[0]).m()/GeV;
      *_mpipi[3] +=(ppi0[0]+ppim[1]).m()/GeV;
      *_mpipi[4] +=(ppip[0]+ppim[0]).m()/GeV;
      *_mpipi[4] +=(ppip[0]+ppim[1]).m()/GeV;
      *_mpipi[5] +=(ppim[0]+ppim[1]).m()/GeV;
      *_mpipipi[2]   += (ppi0[0]+ppip[0]+ppim[0]).m()/GeV;
      *_mpipipi[2]   += (ppi0[0]+ppip[0]+ppim[1]).m()/GeV;
      *_mpipipi[3]   += (ppip[0]+ppim[0]+ppim[1]).m()/GeV;
      *_mpipipi[4]   += (ppi0[0]+ppim[0]+ppim[1]).m()/GeV;
      *_mpipipipi[1] += (ppi0[0]+ppip[0]+ppim[0]+ppim[1]).m()/GeV;
    }
  }
}

NoPIOClassDescription<TauTo4MesonAnalysis> TauTo4MesonAnalysis::initTauTo4MesonAnalysis;
// Definition of the static class description member.

void TauTo4MesonAnalysis::Init() {

  static ClassDocumentation<TauTo4MesonAnalysis> documentation
    ("There is no documentation for the TauTo4MesonAnalysis class");

}

inline void TauTo4MesonAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _mpipi[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P203P2-3 mass in T2-3RN0T1P203P203P203P2-3",
			   "GX XGX X         GX XWGXGXGX XGX XGX XGX X",
			   "1/SdS/dm0P203P2-31/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P203P2-31/GeV",
			   " XGX XGX XX    ");
  _mpipi[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P203P203 mass in T2-3RN0T1P203P203P203P2-3",
			   "GX XGX X         GX XWGXGXGX XGX XGX XGX X",
			   "1/SdS/dm0P203P2031/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P203P2031/GeV",
			   " XGX XGX XX    ");
  _mpipipi[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",
			     "P203P203P203 mass in T2-3RN0T1P203P203P203P2-3",
			     "GX XGX XGX X         GX XWGXGXGX XGX XGX XGX X",
			     "1/SdS/dm0P203P203P2031/GeV2-13",
			     "  G G   XGX XGX XGX XX    X  X",
			     "m0P203P203P2031/GeV",
			     " XGX XGX XGX XX    ");
  _mpipipi[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",
			     "P203P203P2-3 mass in T2-3RN0T1P203P203P203P2-3",
			     "GX XGX XGX X         GX XWGXGXGX XGX XGX XGX X",
			     "1/SdS/dm0P203P203P2-31/GeV2-13",
			     "  G G   XGX XGX XGX XX    X  X",
			     "m0P203P203P2-31/GeV",
			     " XGX XGX XGX XX    ");
  _mpipipipi[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "P203P203P203P2-3 mass in T2-3RN0T1P203P203P203P2-3",
			       "GX XGX XGX XGX X         GX XWGXGXGX XGX XGX XGX X",
			       "1/SdS/dm0P203P203P203P2-31/GeV2-13",
			       "  G G   XGX XGX XGX XGX XX    X  X",
			       "m0P203P203P203P2-31/GeV",
			       " XGX XGX XGX XGX XX    ");
  _mpipi[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P203P2+3 mass in T2-3RN0T1P203P2+3P2-3P2-3",
			   "GX XGX X         GX XWGXGXGX XGX XGX XGX X",
			   "1/SdS/dm0P203P2+31/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P203P2+31/GeV",
			   " XGX XGX XX    ");
  _mpipi[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P203P2-3 mass in T2-3RN0T1P203P2+3P2-3P2-3",
			   "GX XGX X         GX XWGXGXGX XGX XGX XGX X",
			   "1/SdS/dm0P203P2-31/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P203P2-31/GeV",
			   " XGX XGX XX    ");
  _mpipi[4]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2+3P2-3 mass in T2-3RN0T1P203P2+3P2-3P2-3",
			   "GX XGX X         GX XWGXGXGX XGX XGX XGX X",
			   "1/SdS/dm0P2+3P2-31/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P2+3P2-31/GeV",
			   " XGX XGX XX    ");
  _mpipi[5]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P2-3 mass in T2-3RN0T1P203P2+3P2-3P2-3",
			   "GX XGX X         GX XWGXGXGX XGX XGX XGX X",
			   "1/SdS/dm0P2-3P2-31/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P2-3P2-31/GeV",
			   " XGX XGX XX    ");
  _mpipipi[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P203P2+3P2-3 mass in T2-3RN0T1P203P2+3P2-3P2-3",
			   "GX XGX XGX X         GX XWGXGXGX XGX XGX XGX X",
			   "1/SdS/dm0P203P2+3P2-31/GeV2-13",
			   "  G G   XGX XGX XGX XX    X  X",
			   "m0P203P2+3P2-31/GeV",
			   " XGX XGX XGX XX    ");
  _mpipipi[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2+3P2-3P2-3 mass in T2-3RN0T1P203P2+3P2-3P2-3",
			   "GX XGX XGX X         GX XWGXGXGX XGX XGX XGX X",
			   "1/SdS/dm0P2+3P2-3P2-31/GeV2-13",
			   "  G G   XGX XGX XGX XX    X  X",
			   "m0P2+3P2-3P2-31/GeV",
			   " XGX XGX XGX XX    ");
  _mpipipi[4]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P203P2-3P2-3 mass in T2-3RN0T1P203P2+3P2-3P2-3",
			   "GX XGX XGX X         GX XWGXGXGX XGX XGX XGX X",
			   "1/SdS/dm0P203P2-3P2-31/GeV2-13",
			   "  G G   XGX XGX XGX XX    X  X",
			   "m0P203P2-3P2-31/GeV",
			   " XGX XGX XGX XX    ");
  _mpipipipi[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "P203P2+3P2-3P2-3 mass in T2-3RN0T1P203P2+3P2-3P2-3",
			       "GX XGX XGX XGX X         GX XWGXGXGX XGX XGX XGX X",
			       "1/SdS/dm0P203P2+3P2-3P2-31/GeV2-13",
			       "  G G   XGX XGX XGX XGX XX    X  X",
			       "m0P203P2+3P2-3P2-31/GeV",
			       " XGX XGX XGX XGX XX    ");
}

inline void TauTo4MesonAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  for(unsigned int ix=0;ix<5;++ix) {
    _mpipi  .push_back(new_ptr(Histogram(0.,1.8,200)));
    _mpipipi.push_back(new_ptr(Histogram(0.,1.8,200)));
  }
  _mpipi  .push_back(new_ptr(Histogram(0.,1.8,200)));
  for(unsigned int ix=0;ix<2;++ix) {
    _mpipipipi.push_back(new_ptr(Histogram(0.,1.8,200)));
  }
}
