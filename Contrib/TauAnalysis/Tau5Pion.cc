// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tau5Pion class.
//

#include "Tau5Pion.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"


using namespace Herwig;

void Tau5Pion::analyze(tEventPtr event, long ieve, int loop, int state) {
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
    if(tit->second.size()!=6) continue;
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
    if(ppi0.size()==2&&ppim.size()==2&&ppip.size()==1) {
      *_pipi1[0]+=(ppim[0]+ppim[1]).m()/GeV;
      *_pipi1[1]+=(ppim[0]+ppip[0]).m()/GeV;
      *_pipi1[2]+=(ppim[0]+ppi0[0]).m()/GeV;
      *_pipi1[2]+=(ppim[0]+ppi0[1]).m()/GeV;
      *_pipi1[1]+=(ppim[1]+ppip[0]).m()/GeV;
      *_pipi1[2]+=(ppim[1]+ppi0[0]).m()/GeV;
      *_pipi1[2]+=(ppim[1]+ppi0[1]).m()/GeV;
      *_pipi1[3]+=(ppip[0]+ppi0[0]).m()/GeV;
      *_pipi1[3]+=(ppip[0]+ppi0[1]).m()/GeV;
      *_pipi1[4]+=(ppi0[0]+ppi0[1]).m()/GeV;
      *_pipipi1[0]+=(ppim[0]+ppim[1]-ptotal).m()/GeV;
      *_pipipi1[1]+=(ppim[0]+ppip[0]-ptotal).m()/GeV;
      *_pipipi1[2]+=(ppim[0]+ppi0[0]-ptotal).m()/GeV;
      *_pipipi1[2]+=(ppim[0]+ppi0[1]-ptotal).m()/GeV;
      *_pipipi1[1]+=(ppim[1]+ppip[0]-ptotal).m()/GeV;
      *_pipipi1[2]+=(ppim[1]+ppi0[0]-ptotal).m()/GeV;
      *_pipipi1[2]+=(ppim[1]+ppi0[1]-ptotal).m()/GeV;
      *_pipipi1[3]+=(ppip[0]+ppi0[0]-ptotal).m()/GeV;
      *_pipipi1[3]+=(ppip[0]+ppi0[1]-ptotal).m()/GeV;
      *_pipipi1[4]+=(ppi0[0]+ppi0[1]-ptotal).m()/GeV;
      *_pipipipi1[0]+=(ptotal-ppim[0]).m()/GeV;
      *_pipipipi1[0]+=(ptotal-ppim[1]).m()/GeV;
      *_pipipipi1[1]+=(ptotal-ppip[0]).m()/GeV;
      *_pipipipi1[2]+=(ptotal-ppi0[0]).m()/GeV;
      *_pipipipi1[2]+=(ptotal-ppi0[1]).m()/GeV;
      *_q1+=ptotal.m()/GeV;
    }
    else if(ppi0.size()==4&&ppim.size()==1) {
      *_pipi2[0]+=(ppim[0]+ppi0[0]).m()/GeV;
      *_pipi2[0]+=(ppim[0]+ppi0[1]).m()/GeV;
      *_pipi2[0]+=(ppim[0]+ppi0[2]).m()/GeV;
      *_pipi2[0]+=(ppim[0]+ppi0[3]).m()/GeV;
      *_pipi2[1]+=(ppi0[0]+ppi0[1]).m()/GeV;
      *_pipi2[1]+=(ppi0[0]+ppi0[2]).m()/GeV;
      *_pipi2[1]+=(ppi0[0]+ppi0[3]).m()/GeV;
      *_pipi2[1]+=(ppi0[1]+ppi0[2]).m()/GeV;
      *_pipi2[1]+=(ppi0[1]+ppi0[3]).m()/GeV;
      *_pipi2[1]+=(ppi0[2]+ppi0[3]).m()/GeV;
      *_pipipi2[0]+=(ppim[0]+ppi0[0]-ptotal).m()/GeV;
      *_pipipi2[0]+=(ppim[0]+ppi0[1]-ptotal).m()/GeV;
      *_pipipi2[0]+=(ppim[0]+ppi0[2]-ptotal).m()/GeV;
      *_pipipi2[0]+=(ppim[0]+ppi0[3]-ptotal).m()/GeV;
      *_pipipi2[1]+=(ppi0[0]+ppi0[1]-ptotal).m()/GeV;
      *_pipipi2[1]+=(ppi0[0]+ppi0[2]-ptotal).m()/GeV;
      *_pipipi2[1]+=(ppi0[0]+ppi0[3]-ptotal).m()/GeV;
      *_pipipi2[1]+=(ppi0[1]+ppi0[2]-ptotal).m()/GeV;
      *_pipipi2[1]+=(ppi0[1]+ppi0[3]-ptotal).m()/GeV;
      *_pipipi2[1]+=(ppi0[2]+ppi0[3]-ptotal).m()/GeV;
      *_pipipipi2[0]+=(ptotal-ppim[0]).m()/GeV;
      *_pipipipi2[1]+=(ptotal-ppi0[0]).m()/GeV;
      *_pipipipi2[1]+=(ptotal-ppi0[1]).m()/GeV;
      *_pipipipi2[1]+=(ptotal-ppi0[2]).m()/GeV;
      *_pipipipi2[1]+=(ptotal-ppi0[3]).m()/GeV;
      *_q2+=ptotal.m()/GeV;
    }
    else if(ppim.size()==3&&ppip.size()==2) {
      *_pipi3[0]+=(ppip[0]+ppip[1]).m()/GeV;
      *_pipi3[1]+=(ppim[0]+ppip[0]).m()/GeV;
      *_pipi3[1]+=(ppim[0]+ppip[1]).m()/GeV;
      *_pipi3[1]+=(ppim[1]+ppip[0]).m()/GeV;
      *_pipi3[1]+=(ppim[1]+ppip[1]).m()/GeV;
      *_pipi3[1]+=(ppim[2]+ppip[0]).m()/GeV;
      *_pipi3[1]+=(ppim[2]+ppip[1]).m()/GeV;
      *_pipi3[2]+=(ppim[0]+ppim[1]).m()/GeV;
      *_pipi3[2]+=(ppim[0]+ppim[2]).m()/GeV;
      *_pipi3[2]+=(ppim[1]+ppim[2]).m()/GeV;
      *_pipipi3[0]+=(ppip[0]+ppip[1]-ptotal).m()/GeV;
      *_pipipi3[1]+=(ppim[0]+ppip[0]-ptotal).m()/GeV;
      *_pipipi3[1]+=(ppim[0]+ppip[1]-ptotal).m()/GeV;
      *_pipipi3[1]+=(ppim[1]+ppip[0]-ptotal).m()/GeV;
      *_pipipi3[1]+=(ppim[1]+ppip[1]-ptotal).m()/GeV;
      *_pipipi3[1]+=(ppim[2]+ppip[0]-ptotal).m()/GeV;
      *_pipipi3[1]+=(ppim[2]+ppip[1]-ptotal).m()/GeV;
      *_pipipi3[2]+=(ppim[0]+ppim[1]-ptotal).m()/GeV;
      *_pipipi3[2]+=(ppim[0]+ppim[2]-ptotal).m()/GeV;
      *_pipipi3[2]+=(ppim[1]+ppim[2]-ptotal).m()/GeV;
      *_pipipipi3[0]+=(ptotal-ppim[0]).m()/GeV;
      *_pipipipi3[0]+=(ptotal-ppim[1]).m()/GeV;
      *_pipipipi3[0]+=(ptotal-ppim[2]).m()/GeV;
      *_pipipipi3[1]+=(ptotal-ppip[0]).m()/GeV;
      *_pipipipi3[1]+=(ptotal-ppip[1]).m()/GeV;
      *_q3+=ptotal.m()/GeV;
    }
  }
}


NoPIOClassDescription<Tau5Pion> Tau5Pion::initTau5Pion;
// Definition of the static class description member.

void Tau5Pion::Init() {

  static ClassDocumentation<Tau5Pion> documentation
    ("There is no documentation for the Tau5Pion class");

}

void Tau5Pion::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _pipi1[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P2-3 mass in T2-3RN0T12P2-3P2+32P203",
			   "GX XGX X         GX XWGXGX GX XGX X GX X",
			   "1/SdS/dm0P2-3P2-31/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P2-3P2-31/GeV",
			   " XGX XGX XX    ");
  _pipi1[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P2+3 mass in T2-3RN0T12P2-3P2+32P203",
			   "GX XGX X         GX XWGXGX GX XGX X GX X",
			   "1/SdS/dm0P2-3P2+31/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P2-3P2+31/GeV",
			   " XGX XGX XX    ");
  _pipi1[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P203 mass in T2-3RN0T12P2-3P2+32P203",
			   "GX XGX X         GX XWGXGX GX XGX X GX X",
			   "1/SdS/dm0P2-3P2031/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P2-3P2031/GeV",
			   " XGX XGX XX    ");
  _pipi1[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2+3P203 mass in T2-3RN0T12P2-3P2+32P203",
			   "GX XGX X         GX XWGXGX GX XGX X GX X",
			   "1/SdS/dm0P2+3P2031/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P2+3P2031/GeV",
			   " XGX XGX XX    ");
  _pipi1[4]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P203P203 mass in T2-3RN0T12P2-3P2+32P203",
			   "GX XGX X         GX XWGXGX GX XGX X GX X",
			   "1/SdS/dm0P203P2031/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P203P2031/GeV",
			   " XGX XGX XX    ");
  _pipipi1[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",
			     "P2+3P203P203 mass in T2-3RN0T12P2-3P2+32P203",
			     "GX XGX XGX X         GX XWGXGX GX XGX X GX X",
			     "1/SdS/dm0P2+3P203P2031/GeV2-13",
			     "  G G   XGX XGX XGX XX    X  X",
			     "m0P2+3P203P2031/GeV",
			     " XGX XGX XGX XX    ");
  _pipipi1[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",
			     "P2-3P203P203 mass in T2-3RN0T12P2-3P2+32P203",
			     "GX XGX XGX X         GX XWGXGX GX XGX X GX X",
			     "1/SdS/dm0P2-3P203P2031/GeV2-13",
			     "  G G   XGX XGX XGX XX    X  X",
			     "m0P2-3P203P2031/GeV",
			     " XGX XGX XGX XX    ");
  _pipipi1[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",
			     "P2+3P2-3P203 mass in T2-3RN0T12P2-3P2+32P203",
			     "GX XGX XGX X         GX XWGXGX GX XGX X GX X",
			     "1/SdS/dm0P2+3P2-3P2031/GeV2-13",
			     "  G G   XGX XGX XGX XX    X  X",
			     "m0P2+3P2-3P2031/GeV",
			     " XGX XGX XGX XX    ");
  _pipipi1[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",
			     "P2-3P2-3P203 mass in T2-3RN0T12P2-3P2+32P203",
			     "GX XGX XGX X         GX XWGXGX GX XGX X GX X",
			     "1/SdS/dm0P2-3P2-3P2031/GeV2-13",
			     "  G G   XGX XGX XGX XX    X  X",
			     "m0P2-3P2-3P2031/GeV",
			     " XGX XGX XGX XX    ");
  _pipipi1[4]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",
			     "P2+3P2-3P2-3 mass in T2-3RN0T12P2-3P2+32P203",
			     "GX XGX XGX X         GX XWGXGX GX XGX X GX X",
			     "1/SdS/dm0P2+3P2-3P2-31/GeV2-13",
			     "  G G   XGX XGX XGX XX    X  X",
			     "m0P2+3P2-3P2-31/GeV",
			     " XGX XGX XGX XX    ");
  _pipipipi1[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
		     "RED",
		     "P2-3P2+32P203 mass in T2-3RN0T12P2-3P2+32P203",
		     "GX XGX X GX X         GX XWGXGX GX XGX X GX X",
		     "1/SdS/dm0P2-3P2+32P2031/GeV2-13",
		     "  G G   XGX XGX X GX XX    X  X",
		     "m0P2-3P2+32P2031/GeV",
		     " XGX XGX X GX XX    ");
  _pipipipi1[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
		     "RED",
		     "2P2-32P203 mass in T2-3RN0T12P2-3P2+32P203",
		     " GX X GX X         GX XWGXGX GX XGX X GX X",
		     "1/SdS/dm02P2-32P2031/GeV2-13",
		     "  G G   X GX X GX XX    X  X",
		     "m02P2-32P2031/GeV",
		     " X GX X GX XX    ");
  _pipipipi1[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
		     "RED",
		     "2P2-3P2+3P203 mass in T2-3RN0T12P2-3P2+32P203",
		     " GX XGX XGX X         GX XWGXGX GX XGX X GX X",
		     "1/SdS/dm02P2-3P2+3P2031/GeV2-13",
		     "  G G   X GX XGX XGX XX    X  X",
		     "m02P2-3P2+3P2031/GeV",
		     " X GX XGX XGX XX    ");
  _q1->topdrawOutput(output,Frame|Errorbars|Ylog,
		     "RED",
		     "2P2-3P2+32P203 mass in T2-3RN0T12P2-3P2+32P203",
		     " GX XGX X GX X         GX XWGXGX GX XGX X GX X",
		     "1/SdS/dm02P2-3P2+32P2031/GeV2-13",
		     "  G G   X GX XGX X GX XX    X  X",
		     "m02P2-3P2+32P2031/GeV",
		     " X GX XGX X GX XX    ");
  _pipi2[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P203 mass in T2-3RN0T1P2-34P203",
			   "GX XGX X         GX XWGXGXGX X GX X",
			   "1/SdS/dm0P2-3P2031/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P2-3P2031/GeV",
			   " XGX XGX XX    ");
  _pipi2[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P203P203 mass in T2-3RN0T1P2-34P203",
			   "GX XGX X         GX XWGXGXGX X GX X",
			   "1/SdS/dm0P203P2031/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P203P2031/GeV",
			   " XGX XGX XX    ");
  _pipipi2[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",
			     "P203P203P203 mass in T2-3RN0T1P2-34P203",
			     "GX XGX XGX X         GX XWGXGXGX X GX X",
			     "1/SdS/dm0P203P203P2031/GeV2-13",
			     "  G G   XGX XGX XGX XX    X  X",
			     "m0P203P203P2031/GeV",
			     " XGX XGX XGX XX    ");
  _pipipi2[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			     "RED",
			     "P2-3P203P203 mass in T2-3RN0T1P2-34P203",
			     "GX XGX XGX X         GX XWGXGXGX X GX X",
			     "1/SdS/dm0P2-3P203P2031/GeV2-13",
			     "  G G   XGX XGX XGX XX    X  X",
			     "m0P2-3P203P2031/GeV",
			     " XGX XGX XGX XX    ");
  _pipipipi2[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "P203P203P203P203 mass in T2-3RN0T1P2-34P203",
			       "GX XGX XGX XGX X         GX XWGXGXGX X GX X",
			       "1/SdS/dm0P203P203P203P2031/GeV2-13",
			       "  G G   XGX XGX XGX XGX XX    X  X",
			       "m0P203P203P203P2031/GeV",
			       " XGX XGX XGX XGX XX    ");
  _pipipipi2[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "P2-3P203P203P203 mass in T2-3RN0T1P2-34P203",
			       "GX XGX XGX XGX X         GX XWGXGXGX X GX X",
			       "1/SdS/dm0P2-3P203P203P2031/GeV2-13",
			       "  G G   XGX XGX XGX XGX XX    X  X",
			       "m0P2-3P203P203P2031/GeV",
			       " XGX XGX XGX XGX XX    ");
  _q2->topdrawOutput(output,Frame|Errorbars|Ylog,
		     "RED",
		     " mass in T2-3RN0T1P2-34P203",
		     "         GX XWGXGXGX X GX X",
		     "1/SdS/dm0P2-34P2031/GeV2-13",
		     "  G G   XGX X GX XX    X  X",
		     "m0P2-34P2031/GeV",
		     " XGX X GX XX    ");



  _pipi3[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2+3P2+3 mass in T2-3RN0T13P2-32P2+3",
			   "GX XGX X         GX XWGXGX GX X GX X",
			   "1/SdS/dm0P2+3P2+31/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P2+3P2+31/GeV",
			   " XGX XGX XX    ");
  _pipi3[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2+3P2-3 mass in T2-3RN0T13P2-32P2+3",
			   "GX XGX X         GX XWGXGX GX X GX X",
			   "1/SdS/dm0P2+3P2-31/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P2+3P2-31/GeV",
			   " XGX XGX XX    ");
  _pipi3[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P2-3 mass in T2-3RN0T13P2-32P2+3",
			   "GX XGX X         GX XWGXGX GX X GX X",
			   "1/SdS/dm0P2-3P2-31/GeV2-13",
			   "  G G   XGX XGX XX    X  X",
			   "m0P2-3P2-31/GeV",
			   " XGX XGX XX    ");
  _pipipi3[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P2-3P2-3 mass in T2-3RN0T13P2-32P2+3",
			   "GX XGX XGX X         GX XWGXGX GX X GX X",
			   "1/SdS/dm0P2-3P2-3P2-31/GeV2-13",
			   "  G G   XGX XGX XGX XX    X  X",
			   "m0P2-3P2-3P2-31/GeV",
			   " XGX XGX XGX XX    ");
  _pipipi3[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P2-3P2+3 mass in T2-3RN0T13P2-32P2+3",
			   "GX XGX XGX X         GX XWGXGX GX X GX X",
			   "1/SdS/dm0P2-3P2-3P2+31/GeV2-13",
			   "  G G   XGX XGX XGX XX    X  X",
			   "m0P2-3P2-3P2+31/GeV",
			   " XGX XGX XGX XX    ");
  _pipipi3[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P2+3P2+3 mass in T2-3RN0T13P2-32P2+3",
			   "GX XGX XGX X         GX XWGXGX GX X GX X",
			   "1/SdS/dm0P2-3P2+3P2+31/GeV2-13",
			   "  G G   XGX XGX XGX XX    X  X",
			   "m0P2-3P2+3P2+31/GeV",
			   " XGX XGX XGX XX    ");
  _pipipipi3[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P2-3P2+3P2+3 mass in T2-3RN0T13P2-32P2+3",
			   "GX XGX XGX XGX X         GX XWGXGX GX X GX X",
			   "1/SdS/dm0P2-3P2-3P2+3P2+31/GeV2-13",
			   "  G G   XGX XGX XGX XGX XX    X  X",
			   "m0P2-3P2-3P2+3P2+31/GeV",
			   " XGX XGX XGX XGX XX    ");
  _pipipipi3[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			   "RED",
			   "P2-3P2-3P2-3P2+3 mass in T2-3RN0T13P2-32P2+3",
			   "GX XGX XGX XGX X         GX XWGXGX GX X GX X",
			   "1/SdS/dm0P2-3P2-3P2-3P2+31/GeV2-13",
			   "  G G   XGX XGX XGX XGX XX    X  X",
			   "m0P2-3P2-3P2-3P2+31/GeV",
			   " XGX XGX XGX XGX XX    ");
  _q3->topdrawOutput(output,Frame|Errorbars|Ylog,
		     "RED",
		     " mass in T2-3RN0T13P2-32P2+3",
		     "         GX XWGXGX GX X GX X",
		     "1/SdS/dm03P2-32P2+1/GeV2-13",
		     "  G G   X GX X GX X    X  X",
		     "m03P2-32P2+1/GeV",
		     " X GX X GX X    ");

}

void Tau5Pion::doinitrun() {
  AnalysisHandler::doinitrun();
  for(unsigned int ix=0;ix<5;++ix) {
    _pipi1  .push_back(new_ptr(Histogram(0.,1.8,200)));
    _pipipi1.push_back(new_ptr(Histogram(0.,1.8,200)));
  }
  for(unsigned int ix=0;ix<3;++ix) {
    _pipipipi1.push_back(new_ptr(Histogram(0.,1.8,200)));
    _pipi3.push_back(new_ptr(Histogram(0.,1.8,200)));
    _pipipi3.push_back(new_ptr(Histogram(0.,1.8,200)));
  }
  for(unsigned int ix=0;ix<2;++ix) {
    _pipi2.push_back(new_ptr(Histogram(0.,1.8,200)));
    _pipipi2.push_back(new_ptr(Histogram(0.,1.8,200)));
    _pipipipi2.push_back(new_ptr(Histogram(0.,1.8,200)));
    _pipipipi3.push_back(new_ptr(Histogram(0.,1.8,200)));
  }
  _q1=new_ptr(Histogram(0.,1.8,200));
  _q2=new_ptr(Histogram(0.,1.8,200));
  _q3=new_ptr(Histogram(0.,1.8,200));
}
