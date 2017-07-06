// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FermionTest class.
//

#include "FermionTest.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

void FermionTest::analyze(tEventPtr event, long ieve, int loop, int state) {
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

IBPtr FermionTest::clone() const {
  return new_ptr(*this);
}

IBPtr FermionTest::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<FermionTest,AnalysisHandler>
describeHerwigFermionTest("Herwig::FermionTest", "LeptonTest.so");

void FermionTest::Init() {

  static ClassDocumentation<FermionTest> documentation
    ("There is no documentation for the FermionTest class");

}

void FermionTest::dofinish() {
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
