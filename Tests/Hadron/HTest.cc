// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HTest class.
//

#include "HTest.h"
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

void HTest::analyze(tEventPtr event, long, int , int) {
  // Rotate to CMS, extract final state particles and call analyze(particles).
  set<tPPtr> particles;
  event->selectFinalState(inserter(particles));
  tPPtr h0;
  tParticleVector part2;
  for(set<tPPtr>::const_iterator it=particles.begin();
      it!=particles.end();++it) {
    if((**it).id()==ParticleID::h0) h0=*it;
    else part2.push_back(*it);
  }
  if(!h0) return;
  *_higgspt += h0->momentum().perp()/GeV;
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
  }
  if(ptj>40.) {
    *_yj[1] += yjet;
    *_yjyh[1] += yjet-yh;
  }
  if(ptj>80.) {
    *_yj[2] += yjet;
    *_yjyh[2] += yjet-yh;
  }
  for(unsigned int ix=0;ix<3;++ix) *_njet[ix] +=njet[ix];
}

IBPtr HTest::clone() const {
  return new_ptr(*this);
}

IBPtr HTest::fullclone() const {
  return new_ptr(*this);
}

void HTest::persistentOutput(PersistentOStream & ) const {
}

void HTest::persistentInput(PersistentIStream & , int) {
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<HTest,AnalysisHandler>
describeHerwigHTest("Herwig::HTest", "libfastjet.so HadronJetTest.so");

void HTest::Init() {

  static ClassDocumentation<HTest> documentation
    ("There is no documentation for the HTest class");

}

void HTest::dofinish() {
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

void HTest::doinitrun() {
  AnalysisHandler::doinitrun();
  _higgspt = new_ptr(Histogram(0.,1000.,1000));
  _jetpt   = new_ptr(Histogram(0.,1000.,1000));
  for(unsigned int ix=0;ix<3;++ix) {
    _yj  [ix] = new_ptr(Histogram(-10.,10.,200));
    _yjyh[ix] = new_ptr(Histogram(-10.,10.,200));
    _njet[ix] = new_ptr(Histogram(-0.5,10.5,11));
  }
}
