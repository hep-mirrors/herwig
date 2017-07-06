// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VHTest class.
//

#include "VHTest.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Herwig;

void VHTest::analyze(tEventPtr event, long, int , int ) {
  // Rotate to CMS, extract final state particles and call analyze(particles).
  set<tPPtr> particles;
  event->selectFinalState(inserter(particles));
  tPPtr h0;
  tParticleVector part2;
  tParticleVector leptons;
  for(set<tPPtr>::const_iterator it=particles.begin();
      it!=particles.end();++it) {
    if((**it).id()==ParticleID::h0) h0=*it;
    else {
      tPPtr parent=*it;
      do {
	if(abs(parent->id())==ParticleID::Wplus||
	   parent->id()==ParticleID::Z0) break;
	parent = parent->parents()[0];
      }
      while(!parent->parents().empty());
      if(!parent) {
	part2.push_back(*it);
	continue;
      }
      else if(abs(parent->id())==ParticleID::Wplus||
	      parent->id()==ParticleID::Z0) {
	leptons.push_back(*it);
      }
      else {
	part2.push_back(*it);
      }
    }
  }
  if(!h0) return;
  if(leptons.size()!=2) return;
  Lorentz5Momentum pv = leptons[0]->momentum()+leptons[1]->momentum();
  Lorentz5Momentum pvh = pv+h0->momentum();
  *_higgspt += h0->momentum().perp()/GeV;
  *_vpt += pv.perp()/GeV;
  *_vhpt += pvh.perp()/GeV;
  // callFastjet using R-parameter of 1 to get inclusive jets
  vector<fastjet::PseudoJet> fastjet_particles;
  for (unsigned int j=0; j<part2.size(); j++) {
    fastjet::PseudoJet p(part2[j]->momentum().x()/GeV, 
			 part2[j]->momentum().y()/GeV, 
			 part2[j]->momentum().z()/GeV, 
			 part2[j]->momentum().e()/GeV);
    p.set_user_index(j);
    fastjet_particles.push_back(p);
  }
  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, 1.,
				 recomb_scheme, strategy);
  fastjet::ClusterSequence cs(fastjet_particles, jet_def);
  vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
  double yh = h0->momentum().rapidity();
  double yv = pv.rapidity();
  double yvh = pvh.rapidity();
  double yjet = jets[0].rapidity();
  double ptj=jets[0].perp();
  int njet[3]={0,0,0};
  for(unsigned int ix=0;ix<jets.size();++ix) {
    if(jets[ix].perp()>10.) ++njet[0];
    if(jets[ix].perp()>40.) ++njet[1];
    if(jets[ix].perp()>80.) ++njet[2];
  }
  *_jetpt += ptj;
  if(ptj>10.) {
    *_yj[0] += yjet;
    *_yjyh[0] += yjet-yh;
    *_yjyv[0] += yjet-yv;
    *_yjyhv[0] += yjet-yvh;
  }
  if(ptj>40.) {
    *_yj[1] += yjet;
    *_yjyh[1] += yjet-yh;
    *_yjyv[1] += yjet-yv;
    *_yjyhv[1] += yjet-yvh;
  }
  if(ptj>80.) {
    *_yj[2] += yjet;
    *_yjyh[2] += yjet-yh;
    *_yjyv[2] += yjet-yv;
    *_yjyhv[2] += yjet-yvh;
  }
  for(unsigned int ix=0;ix<3;++ix) *_njet[ix] +=njet[ix];
}

IBPtr VHTest::clone() const {
  return new_ptr(*this);
}

IBPtr VHTest::fullclone() const {
  return new_ptr(*this);
}

void VHTest::persistentOutput(PersistentOStream & os) const {
}

void VHTest::persistentInput(PersistentIStream & is, int) {
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VHTest,AnalysisHandler>
describeHerwigVHTest("Herwig::VHTest", "libfastjet.so HadronJetTest.so");

void VHTest::Init() {

  static ClassDocumentation<VHTest> documentation
    ("There is no documentation for the VHTest class");

}

void VHTest::dofinish() {
  AnalysisHandler::dofinish();
  ofstream file;
  string fname = generator()->filename() + string("-") + name() + string(".top");
  file.open(fname.c_str());
  using namespace HistogramOptions;
  _higgspt->topdrawOutput(file,Frame|Ylog,"BLACK","Higgs Pt","",
			  "1/SdS/dp0T1/GeV2-13",
			  "  G G   X X    X  X",
			  "p0T1/GeV",
			  " X X    ");
  _higgspt->normaliseToCrossSection();
  _higgspt->topdrawOutput(file,Frame|Ylog,"BLACK","Higgs Pt","",
			  "dS/dp0T1/nbGeV2-13",
			  " G   X X      X  X",
			  "p0T1/GeV",
			  " X X    ");
  _vpt->topdrawOutput(file,Frame|Ylog,"BLACK","V Pt","",
			  "1/SdS/dp0T1/GeV2-13",
			  "  G G   X X    X  X",
			  "p0T1/GeV",
			  " X X    ");
  _vpt->normaliseToCrossSection();
  _vpt->topdrawOutput(file,Frame|Ylog,"BLACK","V Pt","",
			  "dS/dp0T1/nbGeV2-13",
			  " G   X X      X  X",
			  "p0T1/GeV",
			  " X X    ");
  _vhpt->topdrawOutput(file,Frame|Ylog,"BLACK","VH Pt","",
			  "1/SdS/dp0T1/GeV2-13",
			  "  G G   X X    X  X",
			  "p0T1/GeV",
			  " X X    ");
  _vhpt->normaliseToCrossSection();
  _vhpt->topdrawOutput(file,Frame|Ylog,"BLACK","VH Pt","",
			  "dS/dp0T1/nbGeV2-13",
			  " G   X X      X  X",
			  "p0T1/GeV",
			  " X X    ");
  _jetpt->topdrawOutput(file,Frame|Ylog,"BLACK","Hardest Jet Pt","",
			  "1/SdS/dp0T1/GeV2-13",
			  "  G G   X X    X  X",
			  "p0T1/GeV",
			  " X X    ");
  _jetpt->normaliseToCrossSection();
  _jetpt->topdrawOutput(file,Frame|Ylog,"BLACK","Hardest Jet Pt","",
			  "dS/dp0T1/nbGeV2-13",
			  " G   X X      X  X",
			  "p0T1/GeV",
			  " X X    ");
  string title;
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==0) title ="Hardest jet rapidity p0T1>10 GeV";
    else if(ix==1) title ="Hardest jet rapidity p0T1>40 GeV";
    else if(ix==2) title ="Hardest jet rapidity p0T1>80 GeV";
    _yj[ix]->topdrawOutput(file,Frame,"BLACK",title,
			   "                      X X       ",
			   "1/SdS/dy0j1",
			   "  G G   X X","y0j1"," X X");
    _yj[ix]->normaliseToCrossSection();
    _yj[ix]->topdrawOutput(file,Frame,"BLACK",title,
			   "                      X X       ",
			   "1/SdS/dy0j1/nb",
			   "  G G   X X   ","y0j1"," X X");
  }
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==0) title ="Hardest jet rapidity - higgs rapidity p0T1>10 GeV";
    else if(ix==1) title ="Hardest jet rapidity - higgs rapidity p0T1>40 GeV";
    else if(ix==2) title ="Hardest jet rapidity - higgs rapidity p0T1>80 GeV";
    _yjyh[ix]->topdrawOutput(file,Frame,"BLACK",title,
			     "                                       X X       ",
			     "1/SdS/d(y0j1-y0h1)",
			     "  G G    X X  X X    ","y0j1-y0h1"," X X  X X");
    _yjyh[ix]->normaliseToCrossSection();
    _yjyh[ix]->topdrawOutput(file,Frame,"BLACK",title,
			     "                                       X X       ",
			     "1/SdS/d(y0j1-y0h1)/nb",
			     "  G G    X X  X X    ","y0j1-y0h1"," X X  X X");
    
  }
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==0) title ="Hardest jet rapidity - V rapidity p0T1>10 GeV";
    else if(ix==1) title ="Hardest jet rapidity - V rapidity p0T1>40 GeV";
    else if(ix==2) title ="Hardest jet rapidity - V rapidity p0T1>80 GeV";
    _yjyv[ix]->topdrawOutput(file,Frame,"BLACK",title,
			     "                                   X X       ",
			     "1/SdS/d(y0j1-y0V1)",
			     "  G G    X X  X X    ","y0j1-y0V1"," X X  X X");
    _yjyv[ix]->normaliseToCrossSection();
    _yjyv[ix]->topdrawOutput(file,Frame,"BLACK",title,
			     "                                   X X       ",
			     "1/SdS/d(y0j1-y0V1)/nb",
			     "  G G    X X  X X    ","y0j1-y0V1"," X X  X X");
    
  }
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==0) title ="Hardest jet rapidity - HV rapidity p0T1>10 GeV";
    else if(ix==1) title ="Hardest jet rapidity - HV rapidity p0T1>40 GeV";
    else if(ix==2) title ="Hardest jet rapidity - HV rapidity p0T1>80 GeV";
    _yjyhv[ix]->topdrawOutput(file,Frame,"BLACK",title,
			      "                                    X X       ",
			      "1/SdS/d(y0j1-y0hV1)",
			      "  G G    X X  X  X    ","y0j1-y0hV1"," X X  X  X");
    _yjyhv[ix]->normaliseToCrossSection();
    _yjyhv[ix]->topdrawOutput(file,Frame,"BLACK",title,
			      "                                    X X       ",
			      "1/SdS/d(y0j1-y0hV1)/nb",
			      "  G G    X X  X  X    ","y0j1-y0hV1"," X X  X  X");
    
  }
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==0)      title ="Number of jets p0T1>10 GeV";
    else if(ix==1) title ="Number of jets p0T1>40 GeV";
    else if(ix==2) title ="Number of jets p0T1>80 GeV";
    _njet[ix]->topdrawOutput(file,Frame,"BLACK",title,
			     "                X X       ",
			     "1/SdS/dN0jet1",
			     "  G G   X   X","N0jet1"," X   X");
    _njet[ix]->normaliseToCrossSection();
    _njet[ix]->topdrawOutput(file,Frame,"BLACK",title,
			     "                X X       ",
			     "1/SdS/dN0jet1/nb",
			     "  G G   X   X","N0jet1"," X   X");
  }
}

void VHTest::doinitrun() {
  AnalysisHandler::doinitrun();
  _higgspt = new_ptr(Histogram(0.,1000.,1000));
  _jetpt   = new_ptr(Histogram(0.,1000.,1000));
  _vpt = new_ptr(Histogram(0.,1000.,1000));
  _vhpt = new_ptr(Histogram(0.,1000.,1000));
  for(unsigned int ix=0;ix<3;++ix) {
    _yj  [ix] = new_ptr(Histogram(-10.,10.,200));
    _yjyh[ix] = new_ptr(Histogram(-10.,10.,200));
    _yjyv[ix] = new_ptr(Histogram(-10.,10.,200));
    _yjyhv[ix] = new_ptr(Histogram(-10.,10.,200));
    _njet[ix] = new_ptr(Histogram(-0.5,10.5,11));
  }
}
