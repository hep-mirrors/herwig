// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaMETest class.
//

#include "GammaMETest.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

void GammaMETest::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  tPVector particles = event->getFinalState();
  for(tPVector::const_iterator it=particles.begin();it!=particles.end();++it) {
    if(_cos.find((**it).id())!=_cos.end()) {
      *_cos[(**it).id()] += (**it).momentum().cosTheta();
      *_phi[(**it).id()] += (**it).momentum().phi()+Constants::pi;
      *_y  [(**it).id()] += (**it).momentum().rapidity();
      *_pt [(**it).id()] += (**it).momentum().perp()/GeV;
    }
    else {
      HistogramPtr ncos = new_ptr(Histogram(-1.,1.,200));
      HistogramPtr nphi = new_ptr(Histogram(0.,2.*Constants::pi,200));
      HistogramPtr ny   = new_ptr(Histogram(-10.0,10.0,200));
      HistogramPtr npt  = new_ptr(Histogram(0.,400.,200));
      *ncos += (**it).momentum().cosTheta();
      *nphi += (**it).momentum().phi()+Constants::pi;
      _cos.insert(make_pair((**it).id(),ncos));
      _phi.insert(make_pair((**it).id(),nphi));
      _y  .insert(make_pair((**it).id(),ny  ));
      _pt .insert(make_pair((**it).id(),npt ));
    }
  }
}

IBPtr GammaMETest::clone() const {
  return new_ptr(*this);
}

IBPtr GammaMETest::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<GammaMETest> GammaMETest::initGammaMETest;
// Definition of the static class description member.

void GammaMETest::Init() {

  static ClassDocumentation<GammaMETest> documentation
    ("There is no documentation for the GammaMETest class");

}

void GammaMETest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  for(map<int,HistogramPtr>::const_iterator it=_cos.begin();it!=_cos.end();++it) {
    string title = "cos of polar angle for " + getParticleData(it->first)->PDGName();
    it->second->topdrawOutput(outfile,Frame,"BLACK",title);
  }
  for(map<int,HistogramPtr>::const_iterator it=_y.begin();it!=_y.end();++it) {
    string title = "rapidity for " + getParticleData(it->first)->PDGName();
    it->second->topdrawOutput(outfile,Frame,"BLACK",title);
  }
  for(map<int,HistogramPtr>::const_iterator it=_pt.begin();it!=_pt.end();++it) {
    string title = "pT for " + getParticleData(it->first)->PDGName();
    it->second->topdrawOutput(outfile,Frame,"BLACK",title);
  }
  for(map<int,HistogramPtr>::const_iterator it=_phi.begin();it!=_phi.end();++it) {
    string title = "azimuthal angle for " + getParticleData(it->first)->PDGName();
    it->second->topdrawOutput(outfile,Frame,"BLACK",title);
  }
}
