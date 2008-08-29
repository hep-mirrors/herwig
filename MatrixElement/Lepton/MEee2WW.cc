// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2WW class.
//

#include "MEee2WW.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;

MEee2WW::MEee2WW() {
  setWeightOption(1);
  setSamplingProbability(0.);
}

NoPIOClassDescription<MEee2WW> MEee2WW::initMEee2WW;
// Definition of the static class description member.

void MEee2WW::Init() {

  static ClassDocumentation<MEee2WW> documentation
    ("The MEee2WW class implements the matrix element for"
     " e+e- -> W+W- including the W decays");

}

unsigned int MEee2WW::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2WW::orderInAlphaEW() const {
  return 2;
}

void MEee2WW::getDiagrams() const {
  // define pointers to the particle data objects that are required
  typedef Selector<tDMPtr> DecaySelector;
  tcPDPtr em = getParticleData(ParticleID::eminus);
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  tcPDPtr nu_e = getParticleData(ParticleID::nu_e);
  // decay modes to include
  // for W+
  DecaySelector wpdec = WPlus()->decaySelector();
  vector<PDPair> wpdecays;
  for(DecaySelector::const_iterator cit=wpdec.begin();cit!=wpdec.end();++cit) {
    if(cit->second->orderedProducts().size()!=2) continue;
    if(cit->second->orderedProducts()[0]->id()>0)
      wpdecays.push_back(make_pair(cit->second->orderedProducts()[0],
				   cit->second->orderedProducts()[1]));
    else
      wpdecays.push_back(make_pair(cit->second->orderedProducts()[1],
				   cit->second->orderedProducts()[0]));
    if(wpdecays.back().first->id()<0) swap(wpdecays.back().first,
					   wpdecays.back().second);
  }
  // for W-
  DecaySelector wmdec = WMinus()->decaySelector();
  vector<PDPair> wmdecays;
  for(DecaySelector::const_iterator cit=wmdec.begin();cit!=wmdec.end();++cit) {
    if(cit->second->orderedProducts().size()!=2) continue;
    if(cit->second->orderedProducts()[0]->id()>0)
      wmdecays.push_back(make_pair(cit->second->orderedProducts()[0],
				   cit->second->orderedProducts()[1]));
    else
      wmdecays.push_back(make_pair(cit->second->orderedProducts()[1],
				   cit->second->orderedProducts()[0]));
    if(wmdecays.back().first->id()<0) swap(wmdecays.back().first,
					   wmdecays.back().second);
  }
  // construct the diagrams
  for(unsigned int wp=0;wp<wpdecays.size();++wp) {
    for(unsigned int wm=0;wm<wmdecays.size();++wm) {
      // s-channel Z0 for W+W- production
      add(new_ptr((Tree2toNDiagram(2), em, ep, 1,    Z0(), 3, WMinus(), 3, WPlus(), 
		   4, wmdecays[wm].first,4, wmdecays[wm].second,
		   5, wpdecays[wp].first,5, wpdecays[wp].second, -1)));
      // s-channel photon for W+W- production
      add(new_ptr((Tree2toNDiagram(2), em, ep, 1, gamma(), 3, WMinus(), 3, WPlus(), 
		   4, wmdecays[wm].first,4, wmdecays[wm].second,
		   5, wpdecays[wp].first,5, wpdecays[wp].second, -2)));
	  // t channel for W+W- production
      add(new_ptr((Tree2toNDiagram(3), em, nu_e, ep, 1, WMinus(), 2, WPlus(), 
		   4, wmdecays[wm].first,4, wmdecays[wm].second,
		   5, wpdecays[wp].first,5, wpdecays[wp].second, -3)));
    }
  }
}

Energy2 MEee2WW::scale() const {
  return sHat();
}

double MEee2WW::me2() const {
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

double MEee2WW::helicityME(vector<SpinorWaveFunction> &    f1 ,
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
  tcPDPtr nu_e = getParticleData(ParticleID::nu_e);
  // store intermediate wavefunctions separately for each
  // of the four incoming helicity combinations
  vector<VectorWaveFunction> interP(4), interZ(4), interWp(4), interWm(4);
  size_t index = 0;
  // s-channel diagrams
  // loop over spinor helicity states
  for(unsigned int inhel0 = 0; inhel0 < 2; ++inhel0) {
    for(unsigned int inhel1 = 0; inhel1 < 2 ; ++inhel1, ++index) {
      // evaluate first vertex to get wavefunction of intermediate
      interP.at(index) = vertexFFP()->
	evaluate(sHat(),1,gamma(),f1[inhel0],a1[inhel1]);
      interZ.at(index) = vertexFFZ()->
	evaluate(sHat(),1,Z0()   ,f1[inhel0],a1[inhel1]);
    }
  }
  // intermediate W+/-
  index = 0;
  // loop over spinor helicity states
  for(unsigned int ohel2 = 0; ohel2 < 2; ++ohel2) {
    for(unsigned int ohel3 = 0; ohel3 < 2; ++ohel3, ++index) {
      // evaluate vertex to get wavefunction of intermediate W-
      interWm.at(index) = vertexFFW()->
	evaluate(scale(),1,WMinus(),f2[ohel2],a2[ohel3]);
      // evaluate vertex to get wavefunction of intermediate W+
      interWp.at(index) = vertexFFW()->
	evaluate(scale(),1,WPlus() ,f3[ohel2],a3[ohel3]);
    }
  }
  // need to watch out for coloured products, increases ME by factor 3
  double prefactor = 1.0;
  if (mePartonData()[2]->coloured()) prefactor *= 3.0;
  if (mePartonData()[4]->coloured()) prefactor *= 3.0;
  // t-channel intermediate
  SpinorWaveFunction inter_nu_e;
  // loop over the helicities
  double full_me = 0.0;
  vector<double> me(3,0.0);
  vector<Complex> diag(3,0.0);
  // Wplus
  for (size_t i = 0; i<4; ++i) { 
    // Wminus
    for (size_t j = 0; j<4; ++j) {
      // gamma/Z0/e_nu
      for (size_t k = 0; k<4; ++k) {
	// s-channel photon
	diag[0] = vertexWWW()->evaluate(scale(),interP[k],interWp[i],interWm[j]);
	// s-channel Z0
	diag[1] = vertexWWW()->evaluate(scale(),interZ[k],interWp[i],interWm[j]);
	// t-channel neutrino
	int ihel1 = k/2, ihel2 = k%2;
	inter_nu_e = vertexFFW()->evaluate(scale(),1,nu_e,f1[ihel1],interWm[j]);
	diag[2] = vertexFFW()->evaluate(scale(),inter_nu_e,a1[ihel2],interWp[i]);
	// individual diagrams
        for (size_t ii=0; ii<3; ++ii) {
          me[ii] += prefactor * std::norm(diag[ii]);
        }
	// full matrix element
	diag[0] += diag[1]+diag[2];
        full_me += prefactor * std::norm(diag[0]);
	// storage of the matrix element for spin correlations
	if(calc) menew(ihel1,ihel2,int(j/2),int(j%2),int(i/2),int(i%2)) = diag[0];
      }
    }
  }
  // save individual diagram matrix elements in XComb()
  // and return total ME  (0.25 for spin average)
  DVector save(3);
  for (size_t i = 0; i < 3; ++i) {
    save[i] = 0.25 * me[i];
  }
  meInfo(save);
  return 0.25 * full_me;
}

Selector<MEBase::DiagramIndex>
MEee2WW::diagrams(const DiagramVector & diags) const {
  double lastGamma(1./3.),lastZ(1./3.),lastT(1./3.);
  if ( lastXCombPtr() ) {
    lastGamma = meInfo()[0];
    lastZ     = meInfo()[1];
    lastT     = meInfo()[2];
  }
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 )       sel.insert(lastZ    , i);
    else if ( diags[i]->id() == -2 )  sel.insert(lastGamma, i);
    else if ( diags[i]->id() == -3 )  sel.insert(lastT    , i);
  }
  return sel;
}

Selector<const ColourLines *>
MEee2WW::colourGeometries(tcDiagPtr) const {
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
