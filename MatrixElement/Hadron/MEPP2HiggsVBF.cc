// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsVBF class.
//

#include "MEPP2HiggsVBF.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include <numeric>
#include "Herwig/Utilities/Maths.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "Herwig/Shower/Base/ShowerProgenitor.h"
#include "Herwig/Shower/Base/ShowerTree.h"
#include "Herwig/Shower/Base/Branching.h"
#include "Herwig/Shower/Base/HardTree.h"

using namespace Herwig;

// namespace {
// using namespace Herwig;
// using namespace ThePEG;
// using namespace ThePEG::Helicity;

// void debuggingMatrixElement(bool BGF,
// 			    tcPDPtr partons1, tcPDPtr partons2,
// 			    tcPDPtr partons3, tcPDPtr partons4,
// 			    const Lorentz5Momentum & psystem0,
// 			    const Lorentz5Momentum & psystem1,
// 			    const Lorentz5Momentum & pother0,
// 			    const Lorentz5Momentum & pother1,
// 			    const Lorentz5Momentum & p0,
// 			    const Lorentz5Momentum & p1,
// 			    const Lorentz5Momentum & p2,
// 			    const Lorentz5Momentum & phiggs,
// 			    Energy2 Q12, Energy2 scale,
// 			    double old) {
//   // get the vertex and the boson
//   tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>
//     (CurrentGenerator::current().standardModel());
//   assert(hwsm);
//   tcPDPtr boson;
//   AbstractFFVVertexPtr weakVertex;
//   AbstractFFVVertexPtr strongVertex = hwsm->vertexFFG();
//   AbstractVVSVertexPtr higgsVertex  = hwsm->vertexWWH();
//   if(partons1->id()==partons2->id()) {
//     weakVertex = hwsm->vertexFFZ();
//     boson = hwsm->getParticleData(ParticleID::Z0);
//   }
//   else {
//     weakVertex = hwsm->vertexFFW();
//     boson = hwsm->getParticleData(ParticleID::Wplus);
//   }
//   tcPDPtr hdata = hwsm->getParticleData(ParticleID::h0);
//   tcPDPtr gluon = hwsm->getParticleData(ParticleID::g);
//   SpinorWaveFunction q1,q2;
//   SpinorBarWaveFunction qbar1,qbar2;
//   if(partons1->id()>0) {
//     q1     = SpinorWaveFunction   (psystem0,partons1,incoming);
//     qbar1  = SpinorBarWaveFunction(psystem1,partons2,outgoing);
//   }
//   else {
//     q1     = SpinorWaveFunction   (psystem1,partons2,outgoing);
//     qbar1  = SpinorBarWaveFunction(psystem0,partons1,incoming);
//   }
//   if(partons3->id()>0) {
//     q2    = SpinorWaveFunction   (pother0,partons3,incoming);
//     qbar2 = SpinorBarWaveFunction(pother1,partons4,outgoing);
//   }
//   else {
//     q2    = SpinorWaveFunction   (pother1,partons4,outgoing);
//     qbar2 = SpinorBarWaveFunction(pother0,partons3,incoming);
//   }
//   ScalarWaveFunction higgs(phiggs,hdata,outgoing);
//   if(!BGF) {
//     SpinorWaveFunction q1p;
//     SpinorBarWaveFunction qbar1p;
//     if(partons1->id()>0) {
//       q1p    = SpinorWaveFunction   (p0         ,partons1,incoming);
//       qbar1p = SpinorBarWaveFunction(p1         ,partons2,outgoing);
//     }
//     else {
//       q1p    = SpinorWaveFunction   (p1         ,partons2,outgoing);
//       qbar1p = SpinorBarWaveFunction(p0         ,partons1,incoming);
//     }
//     VectorWaveFunction    gl(p2,gluon,outgoing);
//     double lome(0.),realme(0.);
//     for(unsigned int lhel1=0;lhel1<2;++lhel1) {
//       q2.reset(lhel1);
//       for(unsigned int lhel2=0;lhel2<2;++lhel2) { 
// 	qbar2.reset(lhel2);
// 	VectorWaveFunction off1 
// 	  = weakVertex->evaluate(scale,3,boson,q2,qbar2);
// 	VectorWaveFunction off2
// 	  = higgsVertex->evaluate(scale,3,boson,off1,higgs);
// 	for(unsigned int qhel1=0;qhel1<2;++qhel1) {
// 	  q1.reset(qhel1);
// 	  q1p.reset(qhel1);
// 	  for(unsigned int qhel2=0;qhel2<2;++qhel2) {
// 	    qbar1.reset(qhel2);
// 	    qbar1p.reset(qhel2);
// 	    Complex diag = weakVertex->evaluate(scale,q1,qbar1,off2);
// 	    lome += norm(diag);
// 	    for(unsigned int ghel=0;ghel<2;++ghel) {
// 	      gl.reset(2*ghel);
// 	      SpinorWaveFunction inter1 = 
// 		strongVertex->evaluate(Q12,5,q1p.particle(),q1p,gl);
// 	      Complex diag1 = weakVertex->evaluate(scale,inter1,qbar1p,off2);
// 	      SpinorBarWaveFunction inter2 = 
// 		strongVertex->evaluate(Q12,5,qbar1p.particle(),qbar1p,gl);
// 	      Complex diag2 = weakVertex->evaluate(scale,q1p,inter2,off2);
// 	      realme += norm(diag1+diag2);
// 	    }
// 	  }
// 	}
//       }
//     }
//     double test1 = realme/lome/hwsm->alphaS(Q12)*Q12*UnitRemoval::InvE2;
//     cerr << "testing ratio A " << old/test1 << "\n";
//   }
//   else {
//     SpinorWaveFunction q1p;
//     SpinorBarWaveFunction qbar1p;
//     if(partons1->id()>0) {
//       q1p    = SpinorWaveFunction   (p2,partons1->CC(),outgoing);
//       qbar1p = SpinorBarWaveFunction(p1,partons2      ,outgoing);
//     }
//     else {
//       q1p    = SpinorWaveFunction   (p1,partons2      ,outgoing);
//       qbar1p = SpinorBarWaveFunction(p2,partons1->CC(),outgoing);
//     }
//     VectorWaveFunction    gl(p0,gluon,incoming);
//     double lome(0.),realme(0.);
//     for(unsigned int lhel1=0;lhel1<2;++lhel1) {
//       q2.reset(lhel1);
//       for(unsigned int lhel2=0;lhel2<2;++lhel2) { 
// 	qbar2.reset(lhel2);
// 	VectorWaveFunction off1 
// 	  = weakVertex->evaluate(scale,3,boson,q2,qbar2);
// 	VectorWaveFunction off2
// 	  = higgsVertex->evaluate(scale,3,boson,off1,higgs);
// 	for(unsigned int qhel1=0;qhel1<2;++qhel1) {
// 	  q1.reset(qhel1);
// 	  q1p.reset(qhel1);
// 	  for(unsigned int qhel2=0;qhel2<2;++qhel2) {
// 	    qbar1.reset(qhel2);
// 	    qbar1p.reset(qhel2);
// 	    Complex diag = weakVertex->evaluate(scale,q1,qbar1,off2);
// 	    lome += norm(diag);
// 	    for(unsigned int ghel=0;ghel<2;++ghel) {
// 	      gl.reset(2*ghel);
// 	      SpinorWaveFunction inter1 = 
// 		strongVertex->evaluate(Q12,5,q1p.particle(),q1p,gl);
// 	      Complex diag1 = weakVertex->evaluate(scale,inter1,qbar1p,off2);
// 	      SpinorBarWaveFunction inter2 = 
// 		strongVertex->evaluate(Q12,5,qbar1p.particle(),qbar1p,gl);
// 	      Complex diag2 = weakVertex->evaluate(scale,q1p,inter2,off2);
// 	      realme += norm(diag1+diag2);
// 	    }
// 	  }
// 	}
//       }
//     }
//     double test1 = realme/lome/hwsm->alphaS(Q12)*Q12*UnitRemoval::InvE2;
//     cerr << "testing ratio B " << old/test1 << "\n";
//   }
// }
// }

MEPP2HiggsVBF::MEPP2HiggsVBF() : comptonWeight_(8.), BGFWeight_(30.), 
				 pTmin_(1.*GeV),initial_(10.),final_(8.),
				 procProb_(0.5), comptonInt_(0.), bgfInt_(0.),
				 nover_(0),maxwgt_(make_pair(0.,0.))
{}

void MEPP2HiggsVBF::doinit() {
  gluon_ = getParticleData(ParticleID::g);
  // integrals of me over phase space
  double r5=sqrt(5.),darg((r5-1.)/(r5+1.)),ath(0.5*log((1.+1./r5)/(1.-1./r5)));
  comptonInt_ = 2.*(-21./20.-6./(5.*r5)*ath+sqr(Constants::pi)/3.
		    -2.*Math::ReLi2(1.-darg)-2.*Math::ReLi2(1.-1./darg));
  bgfInt_ = 121./9.-56./r5*ath;
  // get the vertex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm)
    throw InitException() << "Wrong type of StandardModel object in "
			  << "MEPP2HiggsVBF::doinit() the Herwig"
			  << " version must be used" 
			  << Exception::runerror;
  // set the vertex
  setWWHVertex(hwsm->vertexWWH());
  higgs(getParticleData(ParticleID::h0));
  MEfftoffH::doinit();
}

void MEPP2HiggsVBF::dofinish() {
  MEfftoffH::dofinish();
  if(nover_==0) return;
  generator()->log() << "VBFMECorrection when applying the hard correction " 
		     << nover_ << " weights larger than one were generated of which"
		     << " the largest was " << maxwgt_.first << " for the QCD compton"
		     << " processes and "   << maxwgt_.second << " for the BGF process\n";
}

void MEPP2HiggsVBF::getDiagrams() const {
  // get the quark particle data objects as we'll be using them
  tcPDPtr q[6],qbar[6];
  for ( int ix=0; ix<5; ++ix ) {
    q   [ix] = getParticleData(ix+1);
    qbar[ix] = q[ix]->CC();
  }
  // WW processes
  if(process()==0||process()==1) {
    std::vector<pair<tcPDPtr,tcPDPtr> > parentpair;
    parentpair.reserve(6);
    // don't even think of putting 'break' in here!
    switch(maxFlavour()) {
    case 5:
      if (minFlavour()<=4)
      parentpair.push_back(make_pair(getParticleData(ParticleID::b),
				     getParticleData(ParticleID::c)));
      if (minFlavour()<=2)
      parentpair.push_back(make_pair(getParticleData(ParticleID::b),
				     getParticleData(ParticleID::u)));
    case 4:
      if (minFlavour()<=3)
      parentpair.push_back(make_pair(getParticleData(ParticleID::s),
				     getParticleData(ParticleID::c)));
      if (minFlavour()<=1)
      parentpair.push_back(make_pair(getParticleData(ParticleID::d),
				     getParticleData(ParticleID::c)));
    case 3:
      if (minFlavour()<=2)
      parentpair.push_back(make_pair(getParticleData(ParticleID::s),
				     getParticleData(ParticleID::u)));
    case 2:
      if (minFlavour()<=1)
      parentpair.push_back(make_pair(getParticleData(ParticleID::d),
				     getParticleData(ParticleID::u)));
    default:
      ;
    }
    for(unsigned int ix=0;ix<parentpair.size();++ix) {
      for(unsigned int iy=0;iy<parentpair.size();++iy) {
 	// q1 q2 -> q1' q2' h
	if(parentpair[ix].first->id()<parentpair[iy].second->id()) {
	  add(new_ptr((Tree2toNDiagram(4), parentpair[ix].first, WMinus(), WPlus(), 
		       parentpair[iy].second, 1, parentpair[ix].second, 3, 
		       parentpair[iy].first, 2, higgs(),-1)));
	}
	else {
	  add(new_ptr((Tree2toNDiagram(4), parentpair[iy].second, WPlus(), WMinus(), 
		       parentpair[ix].first, 1, parentpair[iy].first, 3,
		       parentpair[ix].second, 2, higgs(),-1)));
	}
	// q1 qbar2 -> q1' qbar2' h
	add(new_ptr((Tree2toNDiagram(4), parentpair[ix].first, WMinus(), WPlus(), 
		     parentpair[iy].first->CC(), 1,
		     parentpair[ix].second, 3, parentpair[iy].second->CC(),
		     2, higgs(),-1)));
	add(new_ptr((Tree2toNDiagram(4),parentpair[iy].second, WPlus(), WMinus(),
		     parentpair[ix].second->CC(), 1, parentpair[iy].first,
		     3, parentpair[ix].first->CC(), 
		     2, higgs(),-1)));
	// qbar1 qbar2 -> qbar1' qbar2' h
	if(parentpair[ix].first->id()<parentpair[ix].second->id()) {
	  add(new_ptr((Tree2toNDiagram(4), parentpair[ix].first->CC(), WPlus(), WMinus(), 
		       parentpair[iy].second->CC(), 1,
		       parentpair[ix].second->CC(), 3, parentpair[iy].first->CC(),
		       2, higgs(),-1))); 
	}
	else {
	  add(new_ptr((Tree2toNDiagram(4), parentpair[iy].second->CC(), WMinus(), WPlus(),
		       parentpair[ix].first->CC(), 1, 
		       parentpair[iy].first->CC(), 3, parentpair[ix].second->CC(),
		       2, higgs(),-1))); 
	}
      }
    }
  }
  // ZZ processes
  if(process()==0||process()==2) {
    for(unsigned int ix=minFlavour()-1;ix<maxFlavour();++ix) {
      for(unsigned int iy=ix;iy<maxFlavour();++iy) {
	// q    q    -> q    q    H
	add(new_ptr((Tree2toNDiagram(4), q[ix], Z0(), Z0(), q[iy], 
		     1, q[ix], 3, q[iy], 2, higgs(),-2))); 
	// qbar qbar -> qbar qbar H
	add(new_ptr((Tree2toNDiagram(4), qbar[ix], Z0(), Z0(), qbar[iy], 
		     1, qbar[ix], 3, qbar[iy], 2, higgs(),-2)));
      }
      // q    qbar -> q    qbar H
      for(unsigned int iy=minFlavour()-1;iy<maxFlavour();++iy) {
	add(new_ptr((Tree2toNDiagram(4), q[ix], Z0(), Z0(), qbar[iy], 
		     1, q[ix], 3, qbar[iy], 2, higgs(),-2))); 
      }
    }
  }
}

void MEPP2HiggsVBF::persistentOutput(PersistentOStream & os) const {
  os << initial_ << final_
     << alpha_ << ounit(pTmin_,GeV) << comptonWeight_ << BGFWeight_ << gluon_
     << comptonInt_ << bgfInt_ << procProb_;
}

void MEPP2HiggsVBF::persistentInput(PersistentIStream & is, int) {
  is >> initial_ >> final_
     >> alpha_ >> iunit(pTmin_,GeV) >> comptonWeight_ >> BGFWeight_ >> gluon_
     >> comptonInt_ >> bgfInt_ >> procProb_;
}

ClassDescription<MEPP2HiggsVBF> MEPP2HiggsVBF::initMEPP2HiggsVBF;
// Definition of the static class description member.

void MEPP2HiggsVBF::Init() {

  static ClassDocumentation<MEPP2HiggsVBF> documentation
    ("The MEPP2HiggsVBF class implements Higgs production via vector-boson fusion");
  static Reference<MEPP2HiggsVBF,ShowerAlpha> interfaceShowerAlphaQCD
    ("ShowerAlphaQCD",
     "The object calculating the strong coupling constant",
     &MEPP2HiggsVBF::alpha_, false, false, true, false, false);

  static Parameter<MEPP2HiggsVBF,Energy> interfacepTMin
    ("pTMin",
     "The minimum pT",
     &MEPP2HiggsVBF::pTmin_, GeV, 1.*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<MEPP2HiggsVBF,double> interfaceComptonWeight
    ("ComptonWeight",
     "Weight for the overestimate ofthe compton channel",
     &MEPP2HiggsVBF::comptonWeight_, 50.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<MEPP2HiggsVBF,double> interfaceBGFWeight
    ("BGFWeight",
     "Weight for the overestimate of the BGF channel",
     &MEPP2HiggsVBF::BGFWeight_, 100.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEPP2HiggsVBF,double> interfaceProcessProbability
    ("ProcessProbability",
     "The probabilty of the QCD compton process for the process selection",
     &MEPP2HiggsVBF::procProb_, 0.3, 0.0, 1.,
     false, false, Interface::limited);

}

HardTreePtr MEPP2HiggsVBF::generateHardest(ShowerTreePtr tree,
					   vector<ShowerInteraction::Type> inter) {
  bool found = false;
  // check if generating QCD radiation
  for(unsigned int ix=0;ix<inter.size();++ix) {
    found |= inter[ix]==ShowerInteraction::QCD;
  }
  if(!found) return HardTreePtr();
  pair<    tShowerParticlePtr,    tShowerParticlePtr> first,second;
  pair<tcBeamPtr,tcBeamPtr> beams;
  pair<tPPtr,tPPtr> hadrons;
  // get the incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(!first.first) {
      first.first    = cit->first->progenitor();
      beams.first    = cit->first->beam();
      hadrons.first  = cit->first->original()->parents()[0];
    }
    else {
      second.first   = cit->first->progenitor();
      beams.second   = cit->first->beam();
      hadrons.second = cit->first->original()->parents()[0];
    }
  }
  // and the outgoing
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
 	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(cjt->first->progenitor()->id()==ParticleID::h0) {
      higgs_ = cjt->first->progenitor();
    }
    else {
      if(abs(cjt->first->progenitor()->id())>5) continue;
      if(cjt->first->progenitor()->colourLine()&&
	 cjt->first->progenitor()->colourLine()==first.first->colourLine()) {
	first.second = cjt->first->progenitor();
	continue;
      }
      if(cjt->first->progenitor()->antiColourLine()&&
	 cjt->first->progenitor()->antiColourLine()==first.first->antiColourLine()) {
	first.second = cjt->first->progenitor();
	continue;
      }
      if(cjt->first->progenitor()->colourLine()&&
	 cjt->first->progenitor()->colourLine()==second.first->colourLine()) {
	second.second = cjt->first->progenitor();
	continue;
      }
      if(cjt->first->progenitor()->antiColourLine()&&
	 cjt->first->progenitor()->antiColourLine()==second.first->antiColourLine()) {
	second.second = cjt->first->progenitor();
	continue;
      }
    }
  }
  // loop over the two possible emitting systems
  q_ [0] = first .second->momentum()-first .first->momentum();
  q2_[0] = -q_[0].m2();
  q_ [1] = second.second->momentum()-second.first->momentum();
  q2_[1] = -q_[1].m2();
  for(unsigned int ix=0;ix<2;++ix) {
    if(ix==1) {
      swap(first,second);
      swap(beams.first,beams.second);
    }
    // check beam, all particles
    assert(beams.first  && higgs_ &&
	   first .first &&  first.second && 
	   second.first && second.second);
    // beam and pdf
    beam_[ix] = beams.first;
    pdf_ [ix] = beam_[ix]->pdf();
    assert(beam_[ix] && pdf_[ix] );
    // Particle data objects
    partons_[ix][0] =  first. first->dataPtr();
    partons_[ix][1] =  first.second->dataPtr();
    partons_[ix][2] = second. first->dataPtr();
    partons_[ix][3] = second.second->dataPtr();
    // extract the born variables
    xB_[ix] = first.first->x();
    Lorentz5Momentum pb     = first.first->momentum();
    Axis axis(q_[ix].vect().unit());
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    rot_[ix] = LorentzRotation();
    if(axis.perp2()>1e-20) {
      rot_[ix].setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      rot_[ix].rotateX(Constants::pi);
    }
    if(abs(1.-q_[ix].e()/q_[ix].vect().mag())>1e-6) 
      rot_[ix].boostZ( q_[ix].e()/q_[ix].vect().mag());
    pb *= rot_[ix];
    if(pb.perp2()/GeV2>1e-20) {
      Boost trans = -1./pb.e()*pb.vect();
      trans.setZ(0.);
      rot_[ix].boost(trans);
    }
    // momenta of the particles
    phiggs_ [ix]    = rot_[ix]*higgs_->momentum();
    pother_ [ix][0] = rot_[ix]*second. first->momentum();
    pother_ [ix][1] = rot_[ix]*second.second->momentum();
    psystem_[ix][0] = rot_[ix]* first. first->momentum();
    psystem_[ix][1] = rot_[ix]* first.second->momentum();
    q_[ix] *= rot_[ix];
    pTCompton_[ix] = pTBGF_[ix] = ZERO;
    // generate a compton point
    generateCompton(ix);
    // generate a BGF point
    generateBGF(ix);
  }
  // no valid emissions, return
  if(pTCompton_[0]<ZERO && pTCompton_[1]<ZERO&&
     pTBGF_    [0]<ZERO && pTBGF_    [1]<ZERO) {
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(QuarkMatcher::Check(cit->first->progenitor()->data()))
	cit->first->maximumpT(pTmin_,ShowerInteraction::QCD);
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(QuarkMatcher::Check(cit->first->progenitor()->data()))
	cit->first->maximumpT(pTmin_,ShowerInteraction::QCD);
    }
    return HardTreePtr();
  }
  // find the maximum pT emission
  unsigned int system = 0;
  bool isCompton = false;
  Energy pTmax = -GeV;
  for(unsigned int ix=0;ix<2;++ix) {
    if(pTCompton_[ix]>pTmax) {
      pTmax = pTCompton_[ix];
      isCompton = true;
      system = ix;
    }
    if(pTBGF_[ix]>pTmax) {
      pTmax = pTBGF_[ix];
      isCompton = false;
      system = ix;
    }
  }
  if(system==0) {
    swap(first,second);
    swap(beams.first,beams.second);
  }
  // add the non emitting particles
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  spaceBranchings.push_back(new_ptr(HardBranching(second.first,SudakovPtr(),
						  HardBranchingPtr(),
						  HardBranching::Incoming)));
  allBranchings.push_back(spaceBranchings.back());
  allBranchings.push_back(new_ptr(HardBranching(second.second,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  allBranchings.push_back(new_ptr(HardBranching(higgs_,SudakovPtr(),
						HardBranchingPtr(),
						HardBranching::Outgoing)));
  allBranchings[0]->colourPartner(allBranchings[1]);
  allBranchings[1]->colourPartner(allBranchings[0]);
  rot_[system].invert();
  // compton hardest
  if(isCompton) {
    for(unsigned int ix=0;ix<ComptonMomenta_[system].size();++ix) {
      ComptonMomenta_[system][ix].transform(rot_[system]);
    }
    ShowerParticlePtr newqout (new_ptr(ShowerParticle(partons_[system][1],true)));
    newqout->set5Momentum(ComptonMomenta_[system][1]);
    ShowerParticlePtr newg(new_ptr(ShowerParticle(gluon_,true)));
    newg->set5Momentum(ComptonMomenta_[system][2]);
    ShowerParticlePtr newqin   (new_ptr(ShowerParticle(partons_[system][0],false )));
    newqin->set5Momentum(ComptonMomenta_[system][0]);
    if(ComptonISFS_[system]) {
      ShowerParticlePtr newspace(new_ptr(ShowerParticle(partons_[system][0],false)));
      newspace->set5Momentum(ComptonMomenta_[system][0]-ComptonMomenta_[system][2]);
      HardBranchingPtr spaceBranch(new_ptr(HardBranching(newqin,SudakovPtr(),
							 HardBranchingPtr(),
							 HardBranching::Incoming)));
      HardBranchingPtr offBranch(new_ptr(HardBranching(newspace,SudakovPtr(),
						       spaceBranch,
						       HardBranching::Incoming)));
      spaceBranch->addChild(offBranch);
      HardBranchingPtr g(new_ptr(HardBranching(newg,SudakovPtr(),spaceBranch,
					       HardBranching::Outgoing)));
      spaceBranch->addChild(g);
      spaceBranch->type(offBranch->branchingParticle()->id()>0 ? 
			ShowerPartnerType::QCDColourLine : ShowerPartnerType::QCDAntiColourLine);
      HardBranchingPtr outBranch(new_ptr(HardBranching(newqout,SudakovPtr(),
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
      spaceBranchings.push_back(spaceBranch);
      allBranchings.push_back(offBranch);
      allBranchings.push_back(outBranch);
      ColinePtr newin(new_ptr(ColourLine())),newout(new_ptr(ColourLine()));
      newin ->addColoured(newqin  ,newspace->dataPtr()->iColour()!=PDT::Colour3);
      newin ->addColoured(newg    ,newspace->dataPtr()->iColour()!=PDT::Colour3);
      newout->addColoured(newspace,newspace->dataPtr()->iColour()!=PDT::Colour3);
      newout->addColoured(newqout ,newspace->dataPtr()->iColour()!=PDT::Colour3);
      newout->addColoured(newg    ,newspace->dataPtr()->iColour()==PDT::Colour3);
    }
    else {
      ShowerParticlePtr newtime(new_ptr(ShowerParticle(partons_[system][1],true)));
      newtime->set5Momentum(ComptonMomenta_[system][1]+ComptonMomenta_[system][2]);
      HardBranchingPtr spaceBranch(new_ptr(HardBranching(newqin,SudakovPtr(),
							 HardBranchingPtr(),
							 HardBranching::Incoming)));
      HardBranchingPtr offBranch(new_ptr(HardBranching(newtime,SudakovPtr(),
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
      HardBranchingPtr g(new_ptr(HardBranching(newg,SudakovPtr(),offBranch,
					       HardBranching::Outgoing)));
      HardBranchingPtr outBranch(new_ptr(HardBranching(newqout,SudakovPtr(),offBranch,
						       HardBranching::Outgoing)));
      offBranch->addChild(outBranch);
      offBranch->addChild(g);
      offBranch->type(offBranch->branchingParticle()->id()>0 ? 
		       ShowerPartnerType::QCDColourLine : ShowerPartnerType::QCDAntiColourLine);
      spaceBranchings.push_back(spaceBranch);
      allBranchings.push_back(spaceBranch);
      allBranchings.push_back(offBranch);

      ColinePtr newin(new_ptr(ColourLine())),newout(new_ptr(ColourLine()));
      newin ->addColoured(newqin  ,newqin->dataPtr()->iColour()!=PDT::Colour3);
      newin ->addColoured(newtime ,newqin->dataPtr()->iColour()!=PDT::Colour3);
      newin ->addColoured(newg    ,newqin->dataPtr()->iColour()!=PDT::Colour3);
      newout->addColoured(newg    ,newqin->dataPtr()->iColour()==PDT::Colour3);
      newout->addColoured(newqout ,newqin->dataPtr()->iColour()!=PDT::Colour3);
    }
  }
  // BGF hardest
  else {
    for(unsigned int ix=0;ix<BGFMomenta_[system].size();++ix) {
      BGFMomenta_[system][ix].transform(rot_[system]);
    }
    ShowerParticlePtr newq   (new_ptr(ShowerParticle(partons_[system][1],true)));
    newq->set5Momentum(BGFMomenta_[system][1]);
    ShowerParticlePtr newqbar(new_ptr(ShowerParticle(partons_[system][0]->CC(),true)));
    newqbar->set5Momentum(BGFMomenta_[system][2]);
    ShowerParticlePtr newg   (new_ptr(ShowerParticle(gluon_,false)));
    newg->set5Momentum(BGFMomenta_[system][0]);
    ShowerParticlePtr newspace(new_ptr(ShowerParticle(partons_[system][0],false)));
    newspace->set5Momentum(BGFMomenta_[system][0]-BGFMomenta_[system][2]);
    HardBranchingPtr spaceBranch(new_ptr(HardBranching(newg,SudakovPtr(),HardBranchingPtr(),
						       HardBranching::Incoming)));
    HardBranchingPtr offBranch(new_ptr(HardBranching(newspace,SudakovPtr(),spaceBranch,
						     HardBranching::Incoming)));
    HardBranchingPtr qbar(new_ptr(HardBranching(newqbar,SudakovPtr(),spaceBranch,
						HardBranching::Outgoing)));
    spaceBranch->addChild(offBranch);
    spaceBranch->addChild(qbar);
    spaceBranch->type(offBranch->branchingParticle()->id()>0 ? 
		     ShowerPartnerType::QCDColourLine : ShowerPartnerType::QCDAntiColourLine);
    HardBranchingPtr outBranch(new_ptr(HardBranching(newq,SudakovPtr(),
						     HardBranchingPtr(),
						     HardBranching::Outgoing)));
    spaceBranchings.push_back(spaceBranch);
    allBranchings.push_back(offBranch);
    allBranchings.push_back(outBranch);

    ColinePtr newin(new_ptr(ColourLine())),newout(new_ptr(ColourLine()));
    newout->addColoured(newspace,newspace->dataPtr()->iColour()!=PDT::Colour3);
    newout->addColoured(newq    ,newspace->dataPtr()->iColour()!=PDT::Colour3);
    newout->addColoured(newg    ,newspace->dataPtr()->iColour()!=PDT::Colour3);
    newin ->addColoured(newg    ,newspace->dataPtr()->iColour()==PDT::Colour3);
    newin ->addColoured(newqbar ,newspace->dataPtr()->iColour()==PDT::Colour3);
  }
  allBranchings[3]->colourPartner(allBranchings[4]);
  allBranchings[4]->colourPartner(allBranchings[3]);
  HardTreePtr newTree(new_ptr(HardTree(allBranchings,spaceBranchings,
				       ShowerInteraction::QCD)));
  // Set the maximum pt for all other emissions and connect hard and shower tree
  Energy pT = isCompton ? pTCompton_[system] : pTBGF_[system];
  // incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    // set maximum pT
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      cit->first->maximumpT(pT,ShowerInteraction::QCD);
    set<HardBranchingPtr>::iterator cjt=newTree->branchings().begin();
    if(cit->first->progenitor()==first.first) {
      ++cjt;++cjt;++cjt;
    }
    newTree->connect(cit->first->progenitor(),*cjt);
    tPPtr beam =cit->first->original();
    if(!beam->parents().empty()) beam=beam->parents()[0];
    (*cjt)->beam(beam);
    HardBranchingPtr parent=(*cjt)->parent();
    while(parent) {
      parent->beam(beam);
      parent=parent->parent();
    };
  }
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    // set maximum pT
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      cit->first->maximumpT(pT,ShowerInteraction::QCD);
    for(set<HardBranchingPtr>::iterator cjt=newTree->branchings().begin();
	cjt!=newTree->branchings().end();++cjt) {
      if((*cjt)->branchingParticle()->isFinalState()&&
	 (*cjt)->branchingParticle()->id()==cit->first->progenitor()->id()) {
	newTree->connect(cit->first->progenitor(),*cjt);
      }
    }
  }
  // set the evolution partners and scales
  ShowerParticleVector particles;
  for(set<HardBranchingPtr>::iterator cit=newTree->branchings().begin();
      cit!=newTree->branchings().end();++cit) {
    particles.push_back((*cit)->branchingParticle());
  }
  for(set<HardBranchingPtr>::iterator cjt=newTree->branchings().begin();
      cjt!=newTree->branchings().end();++cjt) {
    if(cjt==newTree->branchings().begin()) {
      (**cjt).showerMomentum((**cjt).branchingParticle()->momentum());
      ++cjt;
      (**cjt).showerMomentum((**cjt).branchingParticle()->momentum());
      ++cjt;
      (**cjt).showerMomentum((**cjt).branchingParticle()->momentum());
      ++cjt;
    }
  }
  return newTree;
}

void MEPP2HiggsVBF::generateCompton(unsigned int system) {
  // calculate the A coefficient for the correlations
  acoeff_   = A(partons_[system][0],partons_[system][1],
		partons_[system][2],partons_[system][3]);
  // maximum value of the xT
  double xT = sqrt((1.-xB_[system])/xB_[system]);
  double xTMin = 2.*pTmin_/sqrt(q2_[system]);
  double zp;
  // prefactor
  double a = alpha_->overestimateValue()*comptonWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.);
  l_ = 2.*pother_[system][0]/sqrt(q2_[system]);
  m_ = 2.*pother_[system][1]/sqrt(q2_[system]);
  vector<double> azicoeff;
  do {
    wgt = 0.;
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // dz
    zp = UseRandom::rnd();
    xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    // check allowed
    if(xp<xB_[system]||xp>1.) continue;
    // phase-space piece of the weight
    wgt = 8.*(1.-xp)*zp/comptonWeight_;
    // PDF piece of the weight
    Energy2 mu2 = q2_[system]*((1.-xp)*(1-zp)*zp/xp+1.);
    double pdf =  pdf_[system]->xfx(beam_[system],partons_[system][0],
				    mu2    ,xB_[system]/xp)/
                  pdf_[system]->xfx(beam_[system],partons_[system][0],
				    scale(),xB_[system]   );
    wgt *= max(pdf,0.);
    // me piece of the weight
    // double me = comptonME(system,xT,xp,zp,phi);
    double x2 = 1.-(1.-zp)/xp;
    azicoeff = ComptonME(xp,x2,xT,l_,m_);
    double me = 4./3.*alpha_->ratio(0.25*q2_[system]*sqr(xT))*
      (azicoeff[0]+0.5*azicoeff[2]+0.5*azicoeff[4]);
    wgt *= me;
    if(wgt>1.||wgt<0.) {
      ostringstream wstring;
      wstring << "MEPP2HiggsVBF::generateCompton() "
	      << "Weight greater than one or less than zero"
	      << "wgt = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
  }
  while(xT>xTMin&&UseRandom::rnd()>wgt);
  if(xT<=xTMin) {
    pTCompton_[system]=-GeV;
    return;
  }
  // generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi)),sphi(sin(phi));
    phiwgt =  azicoeff[0]+azicoeff[5]*sphi*cphi
      +azicoeff[1]*cphi+azicoeff[2]*sqr(cphi)
      +azicoeff[3]*sphi+azicoeff[4]*sqr(sphi);
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in MEPP2HiggsVBF"
				  << "::generateCompton() to"
				  << " generate phi" << Exception::eventerror;
  // momenta for the configuration
  Energy Q(sqrt(q2_[system]));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		      -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  pTCompton_[system] = 0.5*Q*xT;
  ComptonMomenta_[system].resize(3);
  ComptonMomenta_[system][0] = p0;
  ComptonMomenta_[system][1] = p1;
  ComptonMomenta_[system][2] = p2;
  ComptonISFS_[system] = zp>xp;
}

double MEPP2HiggsVBF::comptonME(unsigned int system, double xT, 
				double xp, double zp, double phi) {
  // scale and prefactors
  double CFfact = 4./3.*alpha_->ratio(0.25*q2_[system]*sqr(xT));
  Energy Q(sqrt(q2_[system]));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  //set NLO momenta
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		      -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  Lorentz5Momentum qnlo = p2+p1-p0; 
  // Breit frame variables
  Lorentz5Momentum r1 = -p0/x1;
  Lorentz5Momentum r2 =  p1/x2;
  // electroweak parameters
  double c0L,c1L,c0R,c1R;
  // W
  if(partons_[system][0]->id()!=partons_[system][1]->id()) {
    c0L = sqrt(0.5);
    c0R = 0;
    c1L = sqrt(0.5);
    c1R = 0;
  }
  // Z
  else {
    if(abs(partons_[system][0]->id())%2==0) {
      c0L = 
	generator()->standardModel()->vu()+
	generator()->standardModel()->au();
      c0R =
	generator()->standardModel()->vu()-
	generator()->standardModel()->au();
    }
    else {
      c0L = 
	generator()->standardModel()->vd()+
	generator()->standardModel()->ad();
      c0R =
	generator()->standardModel()->vd()-
	generator()->standardModel()->ad();
    }
    if(abs(partons_[system][2]->id())%2==0) {
      c1L = 
	generator()->standardModel()->vu()+
	generator()->standardModel()->au();
      c1R =
	generator()->standardModel()->vu()-
	generator()->standardModel()->au();
    }
    else {
      c1L = 
	generator()->standardModel()->vd()+
	generator()->standardModel()->ad();
      c1R =
	generator()->standardModel()->vd()-
	generator()->standardModel()->ad();
    }
    c0L *= 0.25;
    c0R *= 0.25;
    c1L *= 0.25;
    c1R *= 0.25;
  }
  // Matrix element variables
  double G1 = sqr(c0L*c1L)+sqr(c0R*c1R);
  double G2 = sqr(c0L*c1R)+sqr(c0R*c1L);
  Energy4 term1,term2,loME;
  if(partons_[system][0]->id()>0) {
    if(partons_[system][2]->id()>0) {
      term1 = loMatrixElement(r1                 ,pother_[system][0],
			      qnlo+r1            ,pother_[system][1],G1,G2);
      term2 = loMatrixElement(r2-qnlo            ,pother_[system][0],
			      r2                 ,pother_[system][1],G1,G2);
      loME  = loMatrixElement(psystem_[system][0],pother_[system][0],
			      psystem_[system][1],pother_[system][1],G1,G2);
    }
    else {
      term1 = loMatrixElement(r1                 ,pother_[system][1],
			      qnlo+r1            ,pother_[system][0],G1,G2);
      term2 = loMatrixElement(r2-qnlo            ,pother_[system][1],
			      r2                 ,pother_[system][0],G1,G2);
      loME  = loMatrixElement(psystem_[system][0],pother_[system][1],
			      psystem_[system][1],pother_[system][0],G1,G2);
    }
  }
  else {
    if(partons_[system][2]->id()>0) {
      term1 = loMatrixElement(qnlo+r1            ,pother_[system][0],
			      r1                 ,pother_[system][1],G1,G2);
      term2 = loMatrixElement(r2                 ,pother_[system][0],
			      r2-qnlo            ,pother_[system][1],G1,G2);
      loME  = loMatrixElement(psystem_[system][1],pother_[system][0],
			      psystem_[system][0],pother_[system][1],G1,G2);
    }
    else {
      term1 = loMatrixElement(qnlo+r1,pother_[system][1],r1     ,
			      pother_[system][0],G1,G2);
      term2 = loMatrixElement(r2     ,pother_[system][1],r2-qnlo,
			      pother_[system][0],G1,G2);
      loME  = loMatrixElement(psystem_[system][1],pother_[system][1],
			      psystem_[system][0],pother_[system][0],G1,G2);
    }
  }
  double R1 = term1/loME;
  double R2 = sqr(x2)/(sqr(x2)+sqr(xT))*(term2/loME);
//   debuggingMatrixElement(false,
// 			 partons_[system][0],partons_[system][1],
// 			 partons_[system][2],partons_[system][3],
// 			 psystem_[system][0],psystem_[system][1],
// 			 pother_ [system][0],pother_ [system][1],
// 			 p0,p1,p2,phiggs_[system],q2_[system],scale(),
// 			 8.*Constants::pi/(1.-xp)/(1.-zp)*(R1+sqr(xp)*(sqr(x2)+sqr(xT))*R2));
//   cerr << "testing pieces A " << R1 << " " << sqr(xp)*(sqr(x2)+sqr(xT)) << " " << R2 << "\n";
  return CFfact*(R1+sqr(xp)*(sqr(x2)+sqr(xT))*R2);   
}

void MEPP2HiggsVBF::generateBGF(unsigned int system) {
  // maximum value of the xT
  double xT = (1.-xB_[system])/xB_[system];
  double xTMin = 2.*pTmin_/sqrt(q2_[system]);
  double zp;
  // prefactor
  double a = alpha_->overestimateValue()*BGFWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.);
  l_ = 2.*pother_[system][0]/sqrt(q2_[system]);
  m_ = 2.*pother_[system][1]/sqrt(q2_[system]);
  vector<double> azicoeff;
  do {
    wgt = 0.;
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // dzp
    zp = UseRandom::rnd();
    xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    // check allowed
    if(xp<xB_[system]||xp>1.) continue;
    // phase-space piece of the weight
    wgt = 8.*sqr(1.-xp)*zp/BGFWeight_;
    // PDF piece of the weight
    Energy2 mu2 = q2_[system]*((1.-xp)*(1-zp)*zp/xp+1.);
    wgt *= pdf_[system]->xfx(beam_[system],gluon_             ,
			     mu2    ,xB_[system]/xp)/
           pdf_[system]->xfx(beam_[system],partons_[system][0],
			     scale(),xB_[system]);
    // me piece of the weight
    //double me = BGFME(system,xT,xp,zp,phi);
    double x1 = -1./xp;
    double x2 = 1.-(1.-zp)/xp;
    double x3 = 2.+x1-x2;
    azicoeff = BGFME(xp,x2,x3,xT,l_,m_);
    double me = 0.5*alpha_->ratio(0.25*q2_[system]*sqr(xT))*
      (azicoeff[0]+0.5*azicoeff[2]+0.5*azicoeff[4]);
    wgt *= me;
    if(wgt>1.||wgt<0.) {
      ostringstream wstring;
      wstring << "MEPP2HiggsVBF::generateBGF() "
	      << "Weight greater than one or less than zero"
	      << "wgt = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
  }
  while(xT>xTMin&&UseRandom::rnd()>wgt);
  if(xT<=xTMin) {
    pTBGF_[system] = -GeV;
    return;
  }
  // generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi)),sphi(sin(phi));
    phiwgt =  azicoeff[0]+azicoeff[5]*sphi*cphi
      +azicoeff[1]*cphi+azicoeff[2]*sqr(cphi)
      +azicoeff[3]*sphi+azicoeff[4]*sqr(sphi);
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in MEPP2HiggsVBF"
				  << "::generateBGF() to"
				  << " generate phi" << Exception::eventerror;
  // momenta for the configuration
  Energy Q(sqrt(q2_[system]));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		      -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  pTBGF_[system] = 0.5*Q*xT;
  BGFMomenta_[system].resize(3);
  BGFMomenta_[system][0] = p0;
  BGFMomenta_[system][1] = p1;
  BGFMomenta_[system][2] = p2;
}

double MEPP2HiggsVBF::BGFME(unsigned int system, double xT,
			    double xp, double zp, double phi) {
  // scale and prefactors
  double TRfact = 0.5*alpha_->ratio(0.25*q2_[system]*sqr(xT));
  Energy Q(sqrt(q2_[system]));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  // Set NLO momenta
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		      -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  Lorentz5Momentum qnlo = p2+p1-p0; 
  // Breit frame variables
  Lorentz5Momentum r2 =  p1/x2;
  Lorentz5Momentum r3 = -p2/x3;
  // electroweak parameters
  double c0L,c1L,c0R,c1R;
  // W
  if(partons_[system][0]->id()!=partons_[system][1]->id()) {
    c0L = sqrt(0.5);
    c0R = 0;
    c1L = sqrt(0.5);
    c1R = 0;
  }
  // Z
  else {
    if(abs(partons_[system][0]->id())%2==0) {
      c0L = 
	generator()->standardModel()->vu()+
	generator()->standardModel()->au();
      c0R =
	generator()->standardModel()->vu()-
	generator()->standardModel()->au();
    }
    else {
      c0L = 
	generator()->standardModel()->vd()+
	generator()->standardModel()->ad();
      c0R =
	generator()->standardModel()->vd()-
	generator()->standardModel()->ad();
    }
    if(abs(partons_[system][2]->id())%2==0) {
      c1L = 
	generator()->standardModel()->vu()+
	generator()->standardModel()->au();
      c1R =
	generator()->standardModel()->vu()-
	generator()->standardModel()->au();
    }
    else {
      c1L = 
	generator()->standardModel()->vd()+
	generator()->standardModel()->ad();
      c1R =
	generator()->standardModel()->vd()-
	generator()->standardModel()->ad();
    }
    c0L *= 0.25;
    c0R *= 0.25;
    c1L *= 0.25;
    c1R *= 0.25;
  }
  // Matrix element variables
  double G1 = sqr(c0L*c1L)+sqr(c0R*c1R);
  double G2 = sqr(c0L*c1R)+sqr(c0R*c1L);
  Energy4 term2,term3,loME;
  if(partons_[system][0]->id()>0) {
    if(partons_[system][2]->id()>0) {
      term2 = loMatrixElement(r2-qnlo,pother_[system][0],
			      r2     ,pother_[system][1],G1,G2);
      term3 = loMatrixElement(r3     ,pother_[system][0],
			      qnlo+r3,pother_[system][1],G1,G2);
      loME  = loMatrixElement(psystem_[system][0],pother_[system][0],
			      psystem_[system][1],pother_[system][1],G1,G2);
    }
    else {
      term2 = loMatrixElement(r2-qnlo,pother_[system][1],
			      r2     ,pother_[system][0],G1,G2);
      term3 = loMatrixElement(r3     ,pother_[system][1],
			      qnlo+r3,pother_[system][0],G1,G2);
      loME  = loMatrixElement(psystem_[system][0],pother_[system][1],
			      psystem_[system][1],pother_[system][0],G1,G2);
    }
  }
  else {
    if(partons_[system][2]->id()>0) {
      term2 = loMatrixElement(r2     ,pother_[system][0],
			      r2-qnlo,pother_[system][1],G1,G2);
      term3 = loMatrixElement(qnlo+r3,pother_[system][0],
			      r3     ,pother_[system][1],G1,G2);
      loME  = loMatrixElement(psystem_[system][1],pother_[system][0],
			      psystem_[system][0],pother_[system][1],G1,G2);
    }
    else {
      term2 = loMatrixElement(r2     ,pother_[system][1],
			      r2-qnlo,pother_[system][0],G1,G2);
      term3 = loMatrixElement(qnlo+r3,pother_[system][1],
			      r3     ,pother_[system][0],G1,G2);
      loME  = loMatrixElement(psystem_[system][1],pother_[system][1],
			      psystem_[system][0],pother_[system][0],G1,G2);
    }
  }
  double R3 = sqr(x3)/(sqr(x3)+sqr(xT))*(term3/loME);
  double R2 = sqr(x2)/(sqr(x2)+sqr(xT))*(term2/loME);
//   debuggingMatrixElement(true,
// 			 partons_[system][0],partons_[system][1],
// 			 partons_[system][2],partons_[system][3],
// 			 psystem_[system][0],psystem_[system][1],
// 			 pother_ [system][0],pother_ [system][1],
// 			 p0,p1,p2,phiggs_[system],q2_[system],
// 			 8.*Constants::pi/zp/(1.-zp)*(sqr(xp)*(sqr(x3)+sqr(xT))*R3+
// 						      sqr(xp)*(sqr(x2)+sqr(xT))*R2));
  return TRfact*
         (sqr(xp)*(sqr(x3)+sqr(xT))*R3+
          sqr(xp)*(sqr(x2)+sqr(xT))*R2);
   
}

Energy4 MEPP2HiggsVBF::loMatrixElement(const Lorentz5Momentum &p1,
				       const Lorentz5Momentum &p2,
				       const Lorentz5Momentum &q1,
				       const Lorentz5Momentum &q2,
				       double G1, double G2) const {
  return G1*(p1*p2)*(q1*q2) + G2*(p1*q2)*(q1*p2);
}
  
void MEPP2HiggsVBF::initializeMECorrection(ShowerTreePtr tree, double & initial,
					   double & final) {
  systems_.clear();
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data())) {
      systems_.push_back(tChannelPair());
      systems_.back().hadron   = cit->first->original()->parents()[0];
      systems_.back().beam     = cit->first->beam();
      systems_.back().incoming = cit->first->progenitor();
      systems_.back().pdf      = systems_.back().beam->pdf();
    }
  }
  vector<ShowerParticlePtr> outgoing;
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(cjt->first->progenitor()->id()==ParticleID::h0) 
      higgs_ = cjt->first->progenitor();
    else if(QuarkMatcher::Check(cjt->first->progenitor()->data()))
      outgoing.push_back(cjt->first->progenitor());
  }
  assert(outgoing.size()==2&&higgs_);
  // match up the quarks
  for(unsigned int ix=0;ix<systems_.size();++ix) {
    if(systems_[ix].incoming->colourLine()) {
      for(unsigned int iy=0;iy<outgoing.size();++iy) {
	if(outgoing[iy]->colourLine()==systems_[ix].incoming->colourLine()) {
	  systems_[ix].outgoing=outgoing[iy];
	  break;
	}
      }
    }
    else {
      for(unsigned int iy=0;iy<outgoing.size();++iy) {
	if(outgoing[iy]->antiColourLine()==systems_[ix].incoming->antiColourLine()) {
	  systems_[ix].outgoing=outgoing[iy];
	  break;
	}
      }
    }
  }
  assert(systems_[0].outgoing&&systems_[1].outgoing);
  assert(systems_.size()==2);
  initial = initial_;
  final   = final_;
}

void MEPP2HiggsVBF::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  static const double eps = 1e-6;
  // select emitting line
  if(UseRandom::rndbool()) swap(systems_[0],systems_[1]);
  // extract the born variables
  q_[0] = systems_[0].outgoing->momentum()-systems_[0].incoming->momentum();
  q2_[0] = -q_[0].m2();
  Energy Q = sqrt(q2_[0]);
  xB_[0] = systems_[0].incoming->x();
  // construct lorentz transform from lab to breit frame
  Lorentz5Momentum phadron =  systems_[0].hadron->momentum();
  phadron.setMass(0.*GeV);
  phadron.rescaleEnergy();
  Lorentz5Momentum pcmf = phadron+0.5/xB_[0]*q_[0];
  pcmf.rescaleMass();
  LorentzRotation rot(-pcmf.boostVector());
  Lorentz5Momentum pbeam = rot*phadron;
  Axis axis(pbeam.vect().unit());
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  rot.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  Lorentz5Momentum pout = rot*(systems_[1].outgoing->momentum()+higgs_->momentum());
  rot.rotateZ(-atan2(pout.y(),pout.x()));
  // calculate the A coefficient for the correlations
  acoeff_   = A(systems_[0].incoming->dataPtr(),systems_[0].outgoing->dataPtr(),
		systems_[1].incoming->dataPtr(),systems_[1].outgoing->dataPtr());
  vector<double> azicoeff;
  // select the type of process
  bool BGF = UseRandom::rnd()>procProb_;
  double wgt,xp,zp,x1,x2,x3,xperp;
  l_ = 2.*(rot*systems_[1].incoming->momentum())/Q;
  m_ = 2.*(rot*systems_[1].outgoing->momentum())/Q;
  // compton process
  if(!BGF) {
    wgt = generateComptonPoint(xp,zp);
    if(xp<eps) return;
    // common pieces
    Energy2 mu2 = q2_[0]*((1.-xp)*(1-zp)*zp/xp+1);
    wgt *= 2./3./Constants::pi*alpha_->value(mu2)/procProb_;
    // PDF piece
    wgt *= systems_[0].pdf->xfx(systems_[0].beam,
				systems_[0].incoming->dataPtr(),mu2    ,xB_[0]/xp)/
           systems_[0].pdf->xfx(systems_[0].beam,
				systems_[0].incoming->dataPtr(),scale(),xB_[0]   );
    // numerator factors
    wgt /= (1.-xp)*(1.-zp);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = ComptonME(xp,x2,xperp,l_,m_);
  }
  else {
    wgt = generateBGFPoint(xp,zp);
    if(xp<1e-6) return;
    // common pieces 
    Energy2 mu2 = q2_[0]*((1.-xp)*(1-zp)*zp/xp+1);
    wgt *= 0.25/Constants::pi*alpha_->value(mu2)/(1.-procProb_);
    // PDF piece
    wgt *= systems_[0].pdf->xfx(systems_[0].beam,
				gluon_                         ,mu2    ,xB_[0]/xp)/
           systems_[0].pdf->xfx(systems_[0].beam,
				systems_[0].incoming->dataPtr(),scale(),xB_[0]   );
    // numerator factors
    wgt /= (1.-zp);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = BGFME(xp,x2,x3,xperp,l_,m_);
  }
  // compute the azimuthal average of the weight
  wgt *= azicoeff[0]+0.5*(azicoeff[2]+azicoeff[4]);
  // finally factor as picked one line
  wgt *= 2.;
  // decide whether or not to accept the weight
  if(UseRandom::rnd()>wgt) return;
  // if accepted generate generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi)),sphi(sin(phi));
    phiwgt =  azicoeff[0]+azicoeff[5]*sphi*cphi
      +azicoeff[1]*cphi+azicoeff[2]*sqr(cphi)
      +azicoeff[3]*sphi+azicoeff[4]*sqr(sphi);
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in VBFMECorrection"
				  << "::applyHardMatrixElementCorrection() to"
				  << " generate phi" << Exception::eventerror;
  // compute the new incoming and outgoing momenta
  Lorentz5Momentum p1 = Lorentz5Momentum( 0.5*Q*xperp*cos(phi), 0.5*Q*xperp*sin(phi),
					  -0.5*Q*x2,0.*GeV,0.*GeV);
  p1.rescaleEnergy();
  Lorentz5Momentum p2 = Lorentz5Momentum(-0.5*Q*xperp*cos(phi),-0.5*Q*xperp*sin(phi),
					 -0.5*Q*x3,0.*GeV,0.*GeV);
  p2.rescaleEnergy();
  Lorentz5Momentum pin(0.*GeV,0.*GeV,-0.5*x1*Q,-0.5*x1*Q,0.*GeV);
  // debugging code to test vs helicity amplitude expression for matrix elements
//   double cphi(cos(phi)),sphi(sin(phi));
//   double old = (azicoeff[0]+azicoeff[5]*sphi*cphi
// 		+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi)
// 		+azicoeff[3]*sphi+azicoeff[4]*sqr(sphi));
//   if(!BGF) {
//     old *= 8.*Constants::pi/(1.-xp)/(1.-zp);
//   }
//   else {
//     old *= 8.*Constants::pi/zp/(1.-zp);
//   }
//   debuggingMatrixElement(BGF,
// 			 systems_[0].incoming->dataPtr(),
// 			 systems_[0].outgoing->dataPtr(),
// 			 systems_[1].incoming->dataPtr(),
// 			 systems_[1].outgoing->dataPtr(),
// 			 rot*systems_[0].incoming->momentum(),
// 			 rot*systems_[0].outgoing->momentum(),
// 			 rot*systems_[1].incoming->momentum(),
// 			 rot*systems_[1].outgoing->momentum(),
// 			 pin,p1,p2,rot*higgs_->momentum(),
// 			 q2_[0],scale(),old);
  // we need inverse of the rotation, i.e back to lab from breit
  rot.invert();
  // transform the momenta to lab frame
  pin *= rot;
  p1  *= rot;
  p2  *= rot;
  // test to ensure outgoing particles can be put on-shell
  if(!BGF) {
    if(p1.e()<systems_[0].outgoing->dataPtr()->constituentMass()) return;
    if(p2.e()<gluon_                         ->constituentMass()) return;
  }
  else {
    if(p1.e()<systems_[0].outgoing->dataPtr()      ->constituentMass()) return;
    if(p2.e()<systems_[0].incoming->dataPtr()->CC()->constituentMass()) return;
  }
  // stats for weights > 1
  if(wgt>1.) {
    ++nover_;
    if(!BGF) maxwgt_.first  = max(maxwgt_.first ,wgt);
    else     maxwgt_.second = max(maxwgt_.second,wgt);
  }
  // create the new particles and add to ShowerTree
  bool isquark = systems_[0].incoming->colourLine();
  if(!BGF) {
    PPtr newin  = new_ptr(Particle(*systems_[0].incoming));
    newin->set5Momentum(pin);
    PPtr newg   = gluon_                         ->produceParticle(p2 );
    PPtr newout = systems_[0].outgoing->dataPtr()->produceParticle(p1 ); 
    ColinePtr col=isquark ? 
      systems_[0].incoming->colourLine() : systems_[0].incoming->antiColourLine();
    ColinePtr newline=new_ptr(ColourLine());
    // final-state emission
    if(xp>zp) {
      col->removeColoured(newout,!isquark);
      col->addColoured(newin,!isquark);
      col->addColoured(newg,!isquark);
      newline->addColoured(newg,isquark);
      newline->addColoured(newout,!isquark);
    }
    // initial-state emission
    else {
      col->removeColoured(newin ,!isquark);
      col->addColoured(newout,!isquark);
      col->addColoured(newg,isquark);
      newline->addColoured(newg,!isquark);
      newline->addColoured(newin,!isquark);
    }
    PPtr orig;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()!=systems_[0].incoming) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newin);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
      cit->first->progenitor(sp);
      tree->incomingLines()[cit->first]=sp;
      sp->x(xB_[0]/xp);
      cit->first->perturbative(xp>zp);
      if(xp<=zp) orig=cit->first->original();
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(cit->first->progenitor()!=systems_[0].outgoing) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newout);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newout,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(xp<=zp);
      if(xp>zp) orig=cit->first->original();
    }
    assert(orig);
    // add the gluon
    ShowerParticlePtr sg=new_ptr(ShowerParticle(*newg,1,true));
    ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(orig,newg,sg));
    gluon->perturbative(false);
    tree->outgoingLines().insert(make_pair(gluon,sg));
    tree->hardMatrixElementCorrection(true);
  }
  else {
    PPtr newin   = gluon_                   ->produceParticle(pin);
    PPtr newqbar = systems_[0].incoming->dataPtr()->CC()->produceParticle(p2 );
    PPtr newout  = systems_[0].outgoing->dataPtr()      ->produceParticle(p1 );
    ColinePtr col=isquark ? systems_[0].incoming->colourLine() : systems_[0].incoming->antiColourLine();
    ColinePtr newline=new_ptr(ColourLine()); 
    col    ->addColoured(newin  ,!isquark);
    newline->addColoured(newin  , isquark);
    col    ->addColoured(newout ,!isquark);
    newline->addColoured(newqbar, isquark);
    PPtr orig;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()!=systems_[0].incoming) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newin);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
      cit->first->progenitor(sp);
      tree->incomingLines()[cit->first]=sp;
      sp->x(xB_[0]/xp);
      cit->first->perturbative(false);
      orig=cit->first->original();
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(cit->first->progenitor()!=systems_[0].outgoing) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newout);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newout,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(true);
    }
    assert(orig);
    // add the (anti)quark
    ShowerParticlePtr sqbar=new_ptr(ShowerParticle(*newqbar,1,true));
    ShowerProgenitorPtr qbar=new_ptr(ShowerProgenitor(orig,newqbar,sqbar));
    qbar->perturbative(false);
    tree->outgoingLines().insert(make_pair(qbar,sqbar));
    tree->hardMatrixElementCorrection(true);
  }
}

double MEPP2HiggsVBF::A(tcPDPtr qin1, tcPDPtr qout1,
			tcPDPtr qin2, tcPDPtr ) {
  double output;
  // charged current
  if(qin1->id()!=qout1->id()) {
    output = 2;
  }
  // neutral current
  else {
    double cvl,cal,cvq,caq;
    if(abs(qin2->id())%2==0) {
      cvl = generator()->standardModel()->vu();
      cal = generator()->standardModel()->au();
    }
    else {
      cvl = generator()->standardModel()->vd();
      cal = generator()->standardModel()->ad();
    }
    if(abs(qin1->id())%2==0) {
      cvq = generator()->standardModel()->vu();
      caq = generator()->standardModel()->au();
    }
    else {
      cvq = generator()->standardModel()->vd();
      caq = generator()->standardModel()->ad();
    }
    output = 8.*cvl*cal*cvq*caq/(sqr(cvl)+sqr(cal))/(sqr(cvq)+sqr(caq));
  }
  if(qin1->id()<0) output *= -1.;
  if(qin2->id()<0) output *= -1;
  return output;
}

double MEPP2HiggsVBF::generateComptonPoint(double &xp, double & zp) {
  static const double maxwgt = 50.;
  double wgt,xperp2,x2;
  do {
    xp  = UseRandom::rnd();
    double zpmin = xp, zpmax = 1./(1.+xp*(1.-xp));
    zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
    wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
    if(UseRandom::rndbool()) swap(xp,zp);
    xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp;
    x2     = 1.-(1.-zp)/xp;
    wgt *= 2.*(1.+sqr(xp)*(sqr(x2)+1.5*xperp2))/(1.-xp)/(1.-zp);
    if(wgt>maxwgt) 
    if(wgt>maxwgt) {
      ostringstream wstring;
      wstring << "MEPP2HiggsVBF::generateComptonPoint() "
	      << "Weight greater than maximum"
	      << "wgt = " << wgt << " maxwgt = " << maxwgt << "\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
  }
  while(wgt<UseRandom::rnd()*maxwgt);
  return comptonInt_/((1.+sqr(xp)*(sqr(x2)+1.5*xperp2))/(1.-xp)/(1.-zp));
}

double MEPP2HiggsVBF::generateBGFPoint(double &xp, double & zp) {
  static const double maxwgt = 25.;
  double wgt;
  double x2,x3,xperp2;
  do {
    xp = UseRandom::rnd();
    double zpmax = 1./(1.+xp*(1.-xp)), zpmin = 1.-zpmax;
    zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
    wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
    double x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp;
    wgt *= sqr(xp)/(1.-zp)*(sqr(x3)+sqr(x2)+3.*xperp2);
    if(wgt>maxwgt) {
      ostringstream wstring;
      wstring << "DISBase::generateBGFPoint "
	      << "Weight greater than maximum "
	      << "wgt = " << wgt << " maxwgt = 1\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
  }
  while(wgt<UseRandom::rnd()*maxwgt);
  return bgfInt_/sqr(xp)*(1.-zp)/(sqr(x3)+sqr(x2)+3.*xperp2);
}

bool MEPP2HiggsVBF::softMatrixElementVeto(ShowerProgenitorPtr initial,
					  ShowerParticlePtr parent,Branching br) {
  bool veto = !UseRandom::rndbool(parent->isFinalState() ? 1./final_ : 1./initial_);
  // check if me correction should be applied
  long id[2]={initial->id(),parent->id()};
  if(id[0]!=id[1]||id[1]==ParticleID::g) return veto;
  // if not from the right side
  if(initial->progenitor()!=systems_[0].incoming &&
     initial->progenitor()!=systems_[0].outgoing) return veto;
  // get the pT
  Energy pT=br.kinematics->pT();
  // check if hardest so far
  if(pT<initial->highestpT()) return veto;
  double kappa(sqr(br.kinematics->scale())/q2_[0]),z(br.kinematics->z());
  double zk((1.-z)*kappa);
  // final-state
  double wgt(0.);
  if(parent->isFinalState()) {
    double zp=z,xp=1./(1.+z*zk);
    double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    double x2 = 1.-(1.-zp)/xp;
    vector<double> azicoeff = ComptonME(xp,x2,xperp,l_,m_);
    wgt = (azicoeff[0]+0.5*azicoeff[2]+0.5*azicoeff[4])*
      xp/(1.+sqr(z))/final_;
    if(wgt<.0||wgt>1.) {
      ostringstream wstring;
      wstring << "Soft ME correction weight too large or "
	      << "negative for FSR in MEPP2HiggsVBF::"
	      << "softMatrixElementVeto() soft weight " 
	      << " xp = " << xp << " zp = " << zp
	      << " weight = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(), 
					 Exception::warning) );
    }
  }
  else {
    double xp = 2.*z/(1.+zk+sqrt(sqr(1.+zk)-4.*z*zk));
    double zp = 0.5* (1.-zk+sqrt(sqr(1.+zk)-4.*z*zk));
    double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    double x1 = -1./xp, x2 = 1.-(1.-zp)/xp, x3 = 2.+x1-x2;
    // compton
    if(br.ids[0]!=ParticleID::g) {
      vector<double> azicoeff = ComptonME(xp,x2,xperp,l_,m_);
      wgt = (azicoeff[0]+0.5*azicoeff[2]+0.5*azicoeff[4])*
	xp*(1.-z)/(1.-xp)/(1.+sqr(z))/(1.-zp+xp-2.*xp*(1.-zp));
    }
    // BGF
    else {
      vector<double> azicoeff = BGFME(xp,x2,x3,xperp,l_,m_);
      wgt = (azicoeff[0]+0.5*azicoeff[2]+0.5*azicoeff[4])*
	xp/(1.-zp+xp-2.*xp*(1.-zp))/(sqr(z)+sqr(1.-z));
    }
    wgt /=initial_;
    if(wgt<.0||wgt>1.) {
      ostringstream wstring;
      wstring << "Soft ME correction weight too large or "
	      << "negative for ISR in MEPP2HiggsVBF::"
	      << "softMatrixElementVeto() soft weight " 
	      << " xp = " << xp << " zp = " << zp
	      << " weight = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(), 
					 Exception::warning) );
    }
  }
  // if not vetoed
  if(UseRandom::rndbool(wgt)) return false;
  // otherwise
  parent->vetoEmission(br.type,br.kinematics->scale());
  return true;
}

vector<double> MEPP2HiggsVBF::ComptonME(double xp, double x2, double xperp,
					LorentzVector<double> l,
					LorentzVector<double> m) {
  vector<double> output(6,0.);
  double cos2 =   x2 /sqrt(sqr(x2)+sqr(xperp));
  double sin2 = xperp/sqrt(sqr(x2)+sqr(xperp));
  // no phi dependence
  output[0] = l.t()*m.t()-l.z()*m.z()*sqr(cos2)+0.5*acoeff_*cos2*(l.t()*m.z()-l.z()*m.t());
  // cos(phi)
  output[1] = sin2*(-l.x()*m.t()-l.t()*m.x() 
		    + 0.5*acoeff_*cos2*(l.z()*m.x()-m.z()*l.x()));
  // cos(phi)^2
  output[2] = +sqr(sin2)*l.x()*m.x();
  // sin(phi)
  output[3] = sin2*(-l.t()*m.y()-l.y()*m.t()
		    + 0.5*acoeff_*cos2*(l.z()*m.y()-m.z()*l.y()));
  // sin(phi)^2
  output[4] = +sqr(sin2)*l.y()*m.y();
  // sin(phi)cos(phi)
  output[5] = +sqr(sin2)*(m.y()*l.x()+m.x()*l.y());
  // additional factors
  double denom = -l.z()*m.z()+l.t()*m.t()+0.5*acoeff_*(l.t()*m.z()-l.z()*m.t());
  double fact = sqr(xp)*(sqr(x2)+sqr(xperp))/denom;
  for(unsigned int ix=0;ix<output.size();++ix) output[ix] *=fact;
  output[0] += 1.;
  return output;
}

vector<double> MEPP2HiggsVBF::BGFME(double xp, double x2, double x3,
				    double xperp, 
				    LorentzVector<double> l,
				    LorentzVector<double> m) {
  vector<double> output(6,0.);
  double denom = -l.z()*m.z()+l.t()*m.t()+0.5*acoeff_*(l.t()*m.z()-l.z()*m.t());
  double cos2  =   x2 /sqrt(sqr(x2)+sqr(xperp));
  double sin2  = xperp/sqrt(sqr(x2)+sqr(xperp));
  double fact2 = sqr(xp)*(sqr(x2)+sqr(xperp))/denom;
  double cos3  =   x3 /sqrt(sqr(x3)+sqr(xperp));
  double sin3  = xperp/sqrt(sqr(x3)+sqr(xperp));
  double fact3 = sqr(xp)*(sqr(x3)+sqr(xperp))/denom;
  // no phi dependence
  output[0] = 
    fact2*(l.t()*m.t()-l.z()*m.z()*sqr(cos2)
	   + 0.5*acoeff_*cos2*(l.t()*m.z()-l.z()*m.t())) + 
    fact3*(l.t()*m.t()-l.z()*m.z()*sqr(cos3)
	   - 0.5*acoeff_*cos3*(l.t()*m.z()-l.z()*m.t()));
  // cos(phi)
  output[1] = 
    fact2*sin2*( - l.x()*m.t()-l.t()*m.x() 
		 + 0.5*acoeff_*cos2*(l.z()*m.x()-m.z()*l.x())) -
    fact3*sin3*( - l.x()*m.t()-l.t()*m.x() 
		 - 0.5*acoeff_*cos3*(l.z()*m.x()-m.z()*l.x())) ;
  // cos(phi)^2
  output[2] = (fact2*sqr(sin2)+fact3*sqr(sin3))*l.x()*m.x();
  // sin(phi)
  output[3] = 
    fact2*sin2*( - l.t()*m.y()-l.y()*m.t()
		 + 0.5*acoeff_*cos2*(l.z()*m.y()-m.z()*l.y())) -
    fact3*sin3*( - l.t()*m.y()-l.y()*m.t()
		 - 0.5*acoeff_*cos3*(l.z()*m.y()-m.z()*l.y()));
  // sin(phi)^2
  output[4] = (fact2*sqr(sin2)+fact3*sqr(sin3))*l.y()*m.y();
  // sin(phi)cos(phi)
  output[5] = (fact2*sqr(sin2)+fact3*sqr(sin3))*(m.y()*l.x()+m.x()*l.y());
  // return the answer
  return output;
}
