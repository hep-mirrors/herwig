// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2ZZ class.
//

#include "MEee2ZZ.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;

MEee2ZZ::MEee2ZZ() {
  setWeightOption(3);
  setSamplingProbability(0.);
}

NoPIOClassDescription<MEee2ZZ> MEee2ZZ::initMEee2ZZ;
// Definition of the static class description member.

void MEee2ZZ::Init() {

  static ClassDocumentation<MEee2ZZ> documentation
    ("The MEee2ZZ class implements the matrix element for"
     " e+e->Z0Z0 including the Z0 decays");

}

unsigned int MEee2ZZ::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2ZZ::orderInAlphaEW() const {
  return 2;
}

void MEee2ZZ::getDiagrams() const {
  // define pointers to the particle data objects that are required
  typedef Selector<tDMPtr> DecaySelector;
  tcPDPtr em = getParticleData(ParticleID::eminus);
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  // decay modes to include
  DecaySelector zdec = Z0()->decaySelector();
  vector<PDPair> zdecays;
  for(DecaySelector::const_iterator cit=zdec.begin();cit!=zdec.end();++cit) {
    if(cit->second->orderedProducts().size()!=2) continue;
    if(cit->second->orderedProducts()[0]->id()>0)
      zdecays.push_back(make_pair(cit->second->orderedProducts()[0],
				  cit->second->orderedProducts()[1]));
    else
      zdecays.push_back(make_pair(cit->second->orderedProducts()[1],
				  cit->second->orderedProducts()[0]));
    if(zdecays.back().first->id()<0) swap(zdecays.back().first,
					  zdecays.back().second);
  }
  // construct the diagrams
  for(unsigned int z1=0;z1<zdecays.size();++z1) {
    for(unsigned int z2=z1;z2<zdecays.size();++z2) {
      add(new_ptr((Tree2toNDiagram(3), em, em, ep, 1, Z0(), 2, Z0(), 
		   4, zdecays[z1].first,4, zdecays[z1].second,
		   5, zdecays[z2].first,5, zdecays[z2].second, -1)));
      add(new_ptr((Tree2toNDiagram(3), em, em, ep, 1, Z0(), 2, Z0(), 
		   5, zdecays[z1].first,5, zdecays[z1].second,
		   4, zdecays[z2].first,4, zdecays[z2].second, -2)));
    }
  }
}

Energy2 MEee2ZZ::scale() const {
  return sHat();
}

double MEee2ZZ::me2() const {
  // setup momenta and particle data for the external wavefunctions 
  // incoming
  SpinorWaveFunction    em_in( meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction ep_in( meMomenta()[1],mePartonData()[1],incoming);
  // outgoing
  SpinorWaveFunction    ff2_out(meMomenta()[2],mePartonData()[2],outgoing);
  SpinorBarWaveFunction af3_out(meMomenta()[3],mePartonData()[3],outgoing);
  SpinorWaveFunction    ff4_out(meMomenta()[4],mePartonData()[4],outgoing);
  SpinorBarWaveFunction af5_out(meMomenta()[5],mePartonData()[5],outgoing);
  vector<SpinorWaveFunction> f1,f2,f3;
  vector<SpinorBarWaveFunction> a1,a2,a3;
  // calculate the wavefunctions
  for(unsigned int ix=0;ix<2;++ix) {
    em_in.reset(ix);
    f1.push_back(em_in);
    ep_in.reset(ix);
    a1.push_back(ep_in);
    ff2_out.reset(ix);
    f2.push_back(ff2_out);
    af3_out.reset(ix);
    a2.push_back(af3_out);
    ff4_out.reset(ix);
    f3.push_back(ff4_out);
    af5_out.reset(ix);
    a3.push_back(af5_out);
  }
  return helicityME(f1,a1,f2,a2,f3,a3,false)*sqr(sHat()*UnitRemoval::InvE2);
}

double MEee2ZZ::helicityME(vector<SpinorWaveFunction> &    f1 ,
			   vector<SpinorBarWaveFunction> & a1 ,
			   vector<SpinorWaveFunction> &    f2 ,
			   vector<SpinorBarWaveFunction> & a2 ,
			   vector<SpinorWaveFunction> &    f3 ,
			   vector<SpinorBarWaveFunction> & a3 ,
			   bool calc) const {
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // particle data for the t-channel intermediate
  tcPDPtr em = getParticleData(ParticleID::eminus);
  // store intermediate wavefunctions separately for each
  // of the four incoming helicity combinations
  vector<VectorWaveFunction> interZ1(4), interZ2(4);
  // intermediate Z0's
  size_t index = 0;
  // loop over spinor helicity states
  for(unsigned int ohel2 = 0; ohel2 < 2; ++ohel2) {
    for(unsigned int ohel3 = 0; ohel3 < 2; ++ohel3, ++index) {
      // evaluate vertex to get wavefunction of intermediate first Z
      interZ1.at(index) = vertexFFZ()->
	evaluate(scale(),1,Z0() ,f2[ohel2],a2[ohel3]);
      // evaluate vertex to get wavefunction of intermediate second Z
      interZ2.at(index) = vertexFFZ()->
	evaluate(scale(),1,Z0() ,f3[ohel2],a3[ohel3]);
    }
  }
  // need to watch out for coloured products, increases ME by factor 3
  double prefactor = mePartonData()[2]==mePartonData()[4] ? 0.5 : 1.;
  if (mePartonData()[2]->coloured()) prefactor *= 3.0;
  if (mePartonData()[4]->coloured()) prefactor *= 3.0;

  // t-channel intermediate
  SpinorWaveFunction inter;
  // loop over the helicities
  double full_me = 0.0;
  vector<double> me(2,0.0);
  vector<Complex> diag(2,0.0);
  // first  Z0
  for (size_t i = 0; i<4; ++i) { 
    // second Z0
    for (size_t j = 0; j<4; ++j) {
      for(unsigned int ihel1 = 0; ihel1<2;++ihel1) {
	for(unsigned int ihel2 = 0; ihel2<2;++ihel2) {
	  // first diagram
	  inter   = vertexFFZ()->evaluate(scale(),1,em,f1[ihel1],interZ1[j]);
	  diag[0] = vertexFFZ()->evaluate(scale(),inter,a1[ihel2],interZ2[i]);
	  inter   = vertexFFZ()->evaluate(scale(),1,em,f1[ihel1],interZ2[i]);
	  diag[1] = vertexFFZ()->evaluate(scale(),inter,a1[ihel2],interZ1[j]);
	  // individual diagrams
	  for (size_t ii=0; ii<2; ++ii) {
	    me[ii] += prefactor * std::norm(diag[ii]);
	  }
	  // full matrix element
	  diag[0] += diag[1];
	  full_me += prefactor * std::norm(diag[0]);
	  // storage of the matrix element for spin correlations
	  if(calc) menew(ihel1,ihel2,int(j/2),int(j%2),int(i/2),int(i%2)) = diag[0];
	}
      }
    }
  }
  // save individual diagram matrix elements in XComb()
  // and return total ME  (0.25 for spin average)
  DVector save(2);
  for (size_t i = 0; i < 2; ++i) {
    save[i] = 0.25 * me[i];
  }
  meInfo(save);
  return 0.25 * full_me;
}

Selector<MEBase::DiagramIndex>
MEee2ZZ::diagrams(const DiagramVector & diags) const {
  double lastT1(0.5),lastT2(0.5);
  if ( lastXCombPtr() ) {
    lastT1 = meInfo()[0];
    lastT2 = meInfo()[1];
  }
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 )       sel.insert(lastT1 , i);
    else if ( diags[i]->id() == -2 )  sel.insert(lastT2 , i);
  }
  return sel;
}

Selector<const ColourLines *>
MEee2ZZ::colourGeometries(tcDiagPtr) const {
  static ColourLines c1("");
  static ColourLines c2("6 -7");
  static ColourLines c3("8 -9");
  static ColourLines c4("6 -7, 8 -9");
  Selector<const ColourLines *> sel;
  if(mePartonData()[2]->coloured()&&mePartonData()[4]->coloured())
    sel.insert(1.,&c4);
  else if(mePartonData()[2]->coloured())
    sel.insert(1.,&c2);
  else if(mePartonData()[4]->coloured())
    sel.insert(1.,&c3);
  else
    sel.insert(1.,&c1);
  return sel;
}
