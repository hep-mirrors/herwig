// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2Baryons class.
//

#include "MEee2Baryons.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"
#include "ThePEG/Helicity/HelicityFunctions.h"


typedef LorentzVector<complex<InvEnergy> > LorentzPolarizationVectorInvE;

using namespace Herwig;

MEee2Baryons::MEee2Baryons() {}


void MEee2Baryons::getDiagrams() const {
  int idiag=0;
  tPDPtr em    = getParticleData(ParticleID::eminus);
  tPDPtr ep    = getParticleData(ParticleID::eplus);
  tPDPtr gamma = getParticleData(ParticleID::gamma);
  for(unsigned int iform=0;iform<formFactor_->numberOfFactors();++iform) {
    // particles from the form factor
    int id0,id1;
    formFactor_->particleID (iform,id0,id1);
    tcPDPtr p1 = getParticleData( id0);
    tcPDPtr p2 = getParticleData(-id1);
    idiag+=1;
    channelMap_[idiag]=iform;
    add(new_ptr((Tree2toNDiagram(2), em, ep, 1, gamma, 3, p1, 3, p2, -idiag)));
  }
}

Energy2 MEee2Baryons::scale() const {
  return sHat();
}

double MEee2Baryons::me2() const {
  bool cc = false;
  int in  =     mePartonData()[2]->id() ;
  int out = abs(mePartonData()[3]->id());
  Energy m1 = meMomenta()[2].mass();
  Energy m2 = meMomenta()[3].mass();
  unsigned int imode = formFactor_->formFactorNumber(in,out,cc);
  // set up the matrix element
  assert(mePartonData()[2]->iSpin()==PDT::Spin1Half &&
	 mePartonData()[3]->iSpin()==PDT::Spin1Half );
  me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
				    mePartonData()[2]->iSpin(),mePartonData()[3]->iSpin()));
  // compute the lepton current
  SpinorWaveFunction    em_in( meMomenta()[0],mePartonData()[0],incoming);
  SpinorBarWaveFunction ep_in( meMomenta()[1],mePartonData()[1],incoming);
  vector<SpinorWaveFunction> f1;
  vector<SpinorBarWaveFunction> a1;
  for(unsigned int ix=0;ix<2;++ix) {
    em_in.reset(ix);
    f1.push_back(em_in);
    ep_in.reset(ix);
    a1.push_back(ep_in);
  }
  // compute the leptonic current
  LorentzPolarizationVectorInvE lepton[2][2];
  InvEnergy2 pre = SM().alphaEM(sHat())*4.*Constants::pi/sHat();
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      lepton[ix][iy]= pre*f1[ix].dimensionedWave().vectorCurrent(a1[iy].dimensionedWave());
    }
  }
  // baryonic current
  vector<LorentzSpinor   <SqrtEnergy> > f2(2);
  vector<LorentzSpinorBar<SqrtEnergy> > a2(2);
  for(unsigned int ix=0;ix<2;++ix) {
    a2[ix] = HelicityFunctions::dimensionedSpinorBar(-meMomenta()[2],ix,Helicity::outgoing);
    f2[ix] = HelicityFunctions::dimensionedSpinor   (-meMomenta()[3],ix,Helicity::outgoing);
  }
  Lorentz5Momentum q(meMomenta()[2]+meMomenta()[3]);
  q.rescaleMass();
  Energy2 q2 = q.mass2();
  // get the form factors
  Complex f1v(0.),f2v(0.),f3v(0.),f1a(0.),f2a(0.),f3a(0.);
  formFactor_->SpinHalfSpinHalfFormFactor(q2,imode,in,out,m1,m2,
					  f1v,f2v,f3v,f1a,f2a,f3a,
					  BaryonFormFactor::TimeLike);
  Complex left  = f1v - f1a + f2v -double((m1-m2)/(m1+m2))*f2a;
  Complex right = f1v + f1a + f2v +double((m1-m2)/(m1+m2))*f2a;
  Lorentz5Momentum diff = meMomenta()[2]-meMomenta()[3];
  double output(0.);
  for(unsigned int ohel1=0;ohel1<2;++ohel1) {
    for(unsigned int ohel2=0;ohel2<2;++ohel2) {
      LorentzPolarizationVectorE 
	vtemp = f2[ohel2].generalCurrent(a2[ohel1],left,right);
      complex<Energy> vspin=f2[ohel2].scalar      (a2[ohel1]);      
      complex<Energy> aspin=f2[ohel2].pseudoScalar(a2[ohel1]);
      vtemp-= (f2v*vspin+f2a*aspin)/(m1+m2)*diff;
      vtemp+= (f3v*vspin+f3a*aspin)/(m1+m2)*q;
      for(unsigned int ihel1=0;ihel1<2;++ihel1) {
	for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	  Complex amp =  lepton[ihel1][ihel2].dot(vtemp);
	  output += std::norm(amp);
	  me_(ihel1,ihel2,ohel1,ohel2) = amp;
	}
      }
    }
  }
  // prefactors
  output *= 0.25;
  // polarization stuff
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
			 beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    output = me_.average(rho[0],rho[1]);
  }
  return output;
}

unsigned int MEee2Baryons::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2Baryons::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEee2Baryons::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *>
MEee2Baryons::colourGeometries(tcDiagPtr) const {
  static ColourLines none("");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &none);
  return sel;
}

IBPtr MEee2Baryons::clone() const {
  return new_ptr(*this);
}

IBPtr MEee2Baryons::fullclone() const {
  return new_ptr(*this);
}

void MEee2Baryons::doinit() {
  HwMEBase::doinit();
  massOption(vector<unsigned int>(2,1));
}

void MEee2Baryons::persistentOutput(PersistentOStream & os) const {
  os << formFactor_ << channelMap_;
}

void MEee2Baryons::persistentInput(PersistentIStream & is, int) {
  is >> formFactor_ >> channelMap_;
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<MEee2Baryons,HwMEBase>
describeHerwigMEee2Baryons("Herwig::MEee2Baryons",
			   "HwMELeptonLowEnergy.so");

void MEee2Baryons::Init() {

  static ClassDocumentation<MEee2Baryons> documentation
    ("The MEee2Baryons class implements the matrix element for "
     "e+e- > baryon antibaryon using the baryon formfactor");
  
  static Reference<MEee2Baryons,BaryonFormFactor> interfaceFormFactor
    ("FormFactor",
     "The form factor.",
     &MEee2Baryons::formFactor_, false, false, true, false, false);

}

void MEee2Baryons::constructVertex(tSubProPtr) {
}
