// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPPto1S0Jet class.
//

#include "MEPPto1S0Jet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

void MEPPto1S0Jet::doinit() {
  HwMEBase::doinit();
  // get the non-perturbative ME
  O1_ = params_->singletMEProduction<0>(state_,n_,0,0);
  // set the mass option
  massOption(vector<unsigned int>({mOpt_+1,0}));
}

IBPtr MEPPto1S0Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPPto1S0Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPPto1S0Jet::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2) << oenum(state_) << n_ << process_ << mOpt_;
}

void MEPPto1S0Jet::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2) >> ienum(state_) >> n_ >> process_ >> mOpt_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPPto1S0Jet,HwMEBase>
describeHerwigMEPPto1S0Jet("Herwig::MEPPto1S0Jet",
			   "HwOniumParameters.so HwMEHadronOnium.so");

void MEPPto1S0Jet::Init() {

  static ClassDocumentation<MEPPto1S0Jet> documentation
    ("The MEPPto1S0Jet class implements the q qbar -> 1S0 g, g q to 1S0 q"
     " and g g to 1S0 g processes");

  static Reference<MEPPto1S0Jet,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &MEPPto1S0Jet::params_, false, false, true, false, false);
  
  static Switch<MEPPto1S0Jet,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &MEPPto1S0Jet::state_, ccbar, false, false);
  static SwitchOption interfaceStateccbar
    (interfaceState,
     "ccbar",
     "Charmonium state",
     ccbar);
  static SwitchOption interfaceStatebbbar
    (interfaceState,
     "bbbar",
     "Bottomonium state",
     bbbar);
  
  static Parameter<MEPPto1S0Jet,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &MEPPto1S0Jet::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Switch<MEPPto1S0Jet,unsigned int> interfaceProcess
    ("Process",
     "Which processes to generate",
     &MEPPto1S0Jet::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Generate all the processes",
     0);
  static SwitchOption interfaceProcessGQto1S0Q
    (interfaceProcess,
     "GQto1S0Q",
     "The g q -> 1S0 q process",
     1);
  static SwitchOption interfaceProcessGQbarto1S0Qbar
    (interfaceProcess,
     "GQbarto1S0Qbar",
     "The g qbar -> 1S0 qbar process",
     2);
  static SwitchOption interfaceProcessQQbarto1S0G
    (interfaceProcess,
     "QQbarto1S0G",
     "The q qbar -> 1S0 g process",
     3);
  static SwitchOption interfaceProcessGGto1S0G
    (interfaceProcess,
     "GGto1S0G",
     "The g g -> 1S0 g process",
     4);

  static Switch<MEPPto1S0Jet,unsigned int> interfaceMassOption
    ("MassOption",
     "Mass of the treatment of the 1S0 mass",
     &MEPPto1S0Jet::mOpt_, 0, false, false);
  static SwitchOption interfaceMassOptionOnShell
    (interfaceMassOption,
     "OnShell",
     "Use the on-shell mass",
     0);
  static SwitchOption interfaceMassOptionOffShell
    (interfaceMassOption,
     "OffShell",
     "Use an off-shell mass generated by the MassGenerator object for the 1S0 state.",
     1);

}

void MEPPto1S0Jet::getDiagrams() const {
  // construct the meson PDG code from quark ids
  unsigned int iq = 4+state_;
  tcPDPtr ps = getParticleData(long(iq*110+1 + (n_-1)*100000));
  tcPDPtr g = getParticleData(ParticleID::g);
  // processes involving quarks
  for ( int i = 1; i <= 3; ++i ) {
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();
    if(process_ == 0 || process_ == 1)
      add(new_ptr((Tree2toNDiagram(3), g, g, q , 1, ps, 2, q , -1)));
    if(process_ == 0 || process_ == 2)
      add(new_ptr((Tree2toNDiagram(3), g, g, qb, 1, ps, 2, qb, -2)));
    if(process_ == 0 || process_ == 3)
      add(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, ps, 3, g, -3)));
  }
  // g g -> 1S0 g (s,t,u 4-point)
  if(process_ == 0 || process_ == 4) {
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, ps, 3, g , -4)));
    add(new_ptr((Tree2toNDiagram(3), g, g, g, 1, ps, 2, g , -5)));
    add(new_ptr((Tree2toNDiagram(3), g, g, g, 2, ps, 1, g , -6)));
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, ps, 1, g , -7)));
  }
}

Selector<MEBase::DiagramIndex>
MEPPto1S0Jet::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ||
	 diags[i]->id() == -2 ||
	 diags[i]->id() == -3 ) sel.insert(1.0, i);
    else
      sel.insert(meInfo()[abs(diags[i]->id())-4],i);
  return sel;
}

Selector<const ColourLines *>
MEPPto1S0Jet::colourGeometries(tcDiagPtr diag) const {
  // g q -> 1S0 q
  static ColourLines cgq   ("1 2  5, -1 -2  3");
  // g qbar -> 1S0 qbar
  static ColourLines cgqbar("1 2 -3, -1 -2 -5");
  // q qbar -> 1S0 g
  static ColourLines cqqbar("1 3  5, -2 -3 -5");
  // g g -> 1S0 g
  static ColourLines cs[2]={ColourLines("1 3 5, -1 2, -2 -3 -5"),
			    ColourLines("1 -2, -1 -3 -5, 2 3 5")};
  static ColourLines ct[2]={ColourLines("1 2 5, -1 -2 3, -3 -5"),
			    ColourLines("1 2 -3, -1 -2 -5, 3 5")};
  static ColourLines cu[2]={ColourLines("1 5, -1 -2 3, -3 2 -5"),
			    ColourLines("1 2 -3, -1 -5, 3 -2 5")};
  // 4 point
  static ColourLines c4[2]={ColourLines("1 -2,  2  4, -1 -4"),
			    ColourLines("1  4, -4 -2,  2 -1")};
  // create the selector
  Selector<const ColourLines *> sel;
  if      (diag->id() == -1) sel.insert(1.0, &cgq   );
  else if (diag->id() == -2) sel.insert(1.0, &cgqbar);
  else if (diag->id() == -3) sel.insert(1.0, &cqqbar);
  else if (diag->id() == -4) {
    sel.insert(0.5, &cs[0]);
    sel.insert(0.5, &cs[1]);
  }
  else if (diag->id() == -5) {
    sel.insert(0.5, &ct[0]);
    sel.insert(0.5, &ct[1]);
  }
  else if (diag->id() == -6) {
    sel.insert(0.5, &cu[0]);
    sel.insert(0.5, &cu[1]);
  }
  else if (diag->id() == -7) {
    sel.insert(0.5, &c4[0]);
    sel.insert(0.5, &c4[1]);
  }
  return sel;
}

Energy2 MEPPto1S0Jet::scale() const {
  return sHat();
}

double MEPPto1S0Jet::me2() const {
  // return value
  double output(0.);
  // mass of the 1S0 state
  Energy  M  = meMomenta()[2].mass();
  Energy2 M2 = sqr(meMomenta()[2].mass());
  if(mePartonData()[0]->id()==ParticleID::g) {
    // g qbar -> 1S0 qbar
    if(mePartonData()[1]->id()==ParticleID::g) {
      // weights for the different diagrams
      DVector save(4,0.);
      save[0] = tHat()*uHat()/(2.*sqr(tHat()+uHat()));
      save[1] = (tHat()*(sqr(M2)*(sqr(sHat()) +  sHat()*tHat() + sqr(tHat())) - 
			 M2*sHat()*(2*sqr(sHat()) + sHat()*tHat() + sqr(tHat())) + 
			 sqr(sHat())*(2*sqr(sHat()) + 2*sHat()*tHat() + sqr(tHat()))))/
	(2.*sqr(M2 - uHat())*uHat()*sqr(tHat()+uHat()));
      save[2] = (uHat()*(sqr(M2)*(sqr(sHat()) + sHat()*uHat() + sqr(uHat())) - 
			 M2*sHat()*(2*sqr(sHat()) + sHat()*uHat() + sqr(uHat())) + 
			 sqr(sHat())*(2*sqr(sHat()) + 2*sHat()*uHat() + sqr(uHat()))))/
	(2.*sqr(M2 - tHat())*tHat()*sqr(tHat()+uHat()));
      save[3] = (tHat()*uHat()*(pow<4,1>(M2) + sqr(tHat())*sqr(uHat()) - 
				2*pow<3,1>(M2)*(tHat() + uHat()) - 
				M2*tHat()*uHat()*(tHat() + uHat()) + 
				sqr(M2)*(sqr(tHat()) + 
					 3*tHat()*uHat() + sqr(uHat()))))/
	(sqr(M2 - tHat())*sqr(M2 - uHat())*sqr(tHat()+uHat()));
      meInfo(save);
      output = 64./3.*O1_*pow(Constants::pi*standardModel()->alphaS(scale()),3)/M*
	(sqr(-sqr(M2) + M2*sHat() + sqr(tHat()) + tHat()*uHat() + sqr(uHat()))*
	 (pow<4,1>(M2) + pow<4,1>(sHat()) + pow<4,1>(tHat()) + pow<4,1>(uHat())))/
	(4.*sHat()*tHat()*uHat()*sqr((M2-sHat())*(M2-tHat())*(M2-uHat())));
      // test vs NPB 291 731
      // Energy3 R02 = params_->radialWaveFunctionSquared(state_,n_);
      // Energy6 Q(sHat()*tHat()*uHat());
      // Energy4 P(sHat()*tHat()+tHat()*uHat()+uHat()*sHat());
      // double test = 16.*Constants::pi*sqr(sHat())*
      // 	Constants::pi*pow(standardModel()->alphaS(scale()),3)*R02/M/sqr(sHat())/Q/sqr(Q-M2*P)*sqr(P)*
      // 	(pow<4,1>(M2)-2*sqr(M2)*P+sqr(P)+2.*M2*Q);
      // cerr << "testing matrix element " << output << " " << test << " "
      // 	   << (output-test)/(output+test) << " " << output/test << "\n";
    }
    else if(mePartonData()[1]->id()<0) {
      // spin sum version
      double total = 2.*(sqr(sHat())+sqr(uHat()))/pow<4,1>(M);
      // final factors
      output = -32.*O1_*pow<3,1>(M*Constants::pi*standardModel()->alphaS(scale()))/(27.*tHat()*sqr(tHat()-M2))*total;
      // analytic test
      // Energy3 R02 = params_->radialWaveFunctionSquared(state_,n_);
      // double test = -32.*sqr(Constants::pi)*pow(standardModel()->alphaS(scale()),3)*R02*(sqr(sHat())+sqr(uHat()))/(9.*M*tHat()*sqr(tHat()-M2));
      // cerr << "testing matrix element " << output << " " << test << " " << (output-test)/(output+test) << " " << output/test << "\n";
    }
    // g q -> 1S0 q
    else if(mePartonData()[1]->id()<6) {
      // spin sum version
      double total = 2.*(sqr(sHat())+sqr(uHat()))/pow<4,1>(M);
      // final factors
      output = -32.*O1_*pow<3,1>(M*Constants::pi*standardModel()->alphaS(scale()))/(27.*tHat()*sqr(tHat()-M2))*total;
      // analytic test
      // Energy3 R02 = params_->radialWaveFunctionSquared(state_,n_);
      // double test = -32.*sqr(Constants::pi)*pow(standardModel()->alphaS(scale()),3)*R02*(sqr(sHat())+sqr(uHat()))/(9.*M*tHat()*sqr(tHat()-M2));
      // cerr << "testing matrix element " << output << " " << test << " " << (output-test)/(output+test) << " " << output/test << "\n";
    }
    else assert(false);
  }
  // q qbar -> 1S0 g
  else if(mePartonData()[0]->id()==-mePartonData()[1]->id()) {
    // spin sum version
    double total = 2.*(sqr(tHat())+sqr(uHat()))/pow<4,1>(M);
    // final factors
    output = 256.*O1_*pow<3,1>(M*Constants::pi*standardModel()->alphaS(scale()))/(81.*sHat()*sqr(sHat()-M2))*total;
    // analytic test
    // Energy3 R02 = params_->radialWaveFunctionSquared(state_,n_);
    // double test = 256.*sqr(Constants::pi)*pow(standardModel()->alphaS(scale()),3)*
    //   R02*(sqr(tHat())+sqr(uHat()))/(27.*M*sHat()*sqr(sHat()-M2));
    // cerr << "testing matrix element " << output << " " << test << " " << (output-test)/(output+test) << " " << output/test << "\n";
  }
  else
    assert(false);
  return output;
}

void MEPPto1S0Jet::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.reserve(4);
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // get them in the right order
  bool swapped(false);
  if(hard[0]->id()==-hard[1]->id()) {
    if(hard[0]->id()<0) swapped = true;
  }
  else if(hard[0]->id()!=ParticleID::g) {
    swapped=true;
  }
  if(swapped) {
    swap(hard[0],hard[1]);
    swap(hard[2],hard[3]);
  }
  // boost to partonic CMS
  Lorentz5Momentum pcms = hard[0]->momentum()+hard[1]->momentum();
  LorentzRotation boost(-pcms.boostVector());
  for(PPtr part : hard) part->transform(boost);
  // extract kinematic variables
  Energy  M  = hard[2]->mass();
  Energy2 M2 = sqr(M);
  double phi = hard[2]->momentum().phi();
  Energy2 sh = (hard[0]->momentum()+hard[1]->momentum()).m2();
  Energy2 th = (hard[0]->momentum()-hard[2]->momentum()).m2();
  Energy2 uh = (hard[0]->momentum()-hard[3]->momentum()).m2();
  // set basis states and compute the matrix element
  ProductionMatrixElement me;
  ScalarWaveFunction(      hard[2],outgoing,true);
  if(hard[0]->id()==ParticleID::g) {
    // g g -> 1S0 g
    if(hard[1]->id()==ParticleID::g) {
      vector<VectorWaveFunction> g1,g2,g4;
      VectorWaveFunction( g1,hard[0],incoming,false, true,true,vector_phase);
      VectorWaveFunction( g2,hard[1],incoming,false, true,true,vector_phase);
      VectorWaveFunction( g4,hard[3],outgoing,true , true,true,vector_phase);
      ProductionMatrixElement me(PDT::Spin1,PDT::Spin1,PDT::Spin0,PDT::Spin1);
      Complex phase = exp(Complex(0.,phi));
      me(0,0,0,0) =  phase;
      me(0,0,0,2) = -sqr(M2)/(phase*sqr(sh));
      me(0,2,0,0) = sqr(th)/(phase*sqr(sh));
      me(0,2,0,2) = sqr(uh)/(pow(phase,3)*sqr(sh));
      me(2,0,0,0) = (pow(phase,3)*sqr(uh))/sqr(sh);
      me(2,0,0,2) = (phase*sqr(th))/sqr(sh);
      me(2,2,0,0) = -((sqr(M2)*phase)/sqr(sh));
      me(2,2,0,2) = 1./phase;
      // test the average result
      // double aver = me.average();
      // double test = 2.*(pow<4,1>(M2)+pow<4,1>(sh)+pow<4,1>(th)+pow<4,1>(uh))/pow<4,1>(sh);
      // cerr << "testing spin correlations " << test << " " << me.average() << " "
      // 	   << abs(test-aver)/(test+aver) << "\n";
    }
    // g qbar -> 1S0 qbar
    else if(hard[1]->id()<0) {
      vector<VectorWaveFunction>    g1;
      vector<SpinorBarWaveFunction> q2;
      vector<SpinorWaveFunction>    q4;
      VectorWaveFunction(   g1,hard[0],incoming,false,true,true,vector_phase);
      SpinorBarWaveFunction(q2,hard[1],incoming,false,true);
      SpinorWaveFunction(   q4,hard[3],outgoing,true ,true);
      g1[1]=g1[2];
      // matrix element
      me = ProductionMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin0,PDT::Spin1Half);
      if(!swapped) {
      	me(0,0,0,0) = sh/M2;
      	me(0,1,0,1) =-exp(Complex(0.,-2.*phi))*uh/M2;
      	me(2,1,0,1) = sh/M2;
      	me(2,0,0,0) =-exp(Complex(0., 2.*phi))*uh/M2;
      }
      else {
      	me(0,0,0,0) = -exp(Complex(0., phi))*sh/M2;
      	me(0,1,0,1) =  exp(Complex(0., phi))*uh/M2;
      	me(2,1,0,1) = -exp(Complex(0.,-phi))*sh/M2;
      	me(2,0,0,0) =  exp(Complex(0.,-phi))*uh/M2;
      }
      // Helicity code version
      // complex<InvEnergy3> fact = sqrt(-2/th)/M2*Complex(0.,1.);
      // for(unsigned int ih2=0;ih2<2;++ih2) {
      //   for(unsigned int ih4=0;ih4<2;++ih4) {
      //     LorentzPolarizationVectorE fCurrent = q4[ih4].dimensionedWave().vectorCurrent(q2[ih2].dimensionedWave());
      //     LorentzPolarizationVector eps = fact*Helicity::epsilon(fCurrent,hard[2]->momentum(),hard[0]->momentum());
      //     for(unsigned int ih1=0;ih1<2;++ih1) {
      // 	Complex amp = eps*g1[ih1].wave();
      // 	if(norm(me(2*ih1,ih2,0,ih4))>1e-10)  cerr << "testing in hel loop B " << ih1 << " " << ih2 << " " << ih4 << " "
      // 						  << amp << " " << me(2*ih1,ih2,0,ih4) << " " << amp/me(2*ih1,ih2,0,ih4)
      // 						  << " " << norm(amp/me(2*ih1,ih2,0,ih4)) << "\n";
      //     }
      //   }
      // }
    }
    // g q -> 1S0 q
    else if(hard[1]->id()<6) {
      vector<VectorWaveFunction> g1;
      vector<SpinorWaveFunction> q2;
      vector<SpinorBarWaveFunction> q4;
      VectorWaveFunction(   g1,hard[0],incoming,false,true,true,vector_phase);
      SpinorWaveFunction(   q2,hard[1],incoming,false,true);
      SpinorBarWaveFunction(q4,hard[3],outgoing,true ,true);
      g1[1]=g1[2];
      // matrix element
      me = ProductionMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin0,PDT::Spin1Half);
      if(!swapped) {
	me(0,0,0,0) = sh/M2;
	me(0,1,0,1) =-exp(Complex(0.,-2.*phi))*uh/M2;
	me(2,1,0,1) = sh/M2;
	me(2,0,0,0) =-exp(Complex(0., 2.*phi))*uh/M2;
      }
      else {
	me(0,0,0,0) = -exp(Complex(0., phi))*sh/M2;
	me(0,1,0,1) =  exp(Complex(0., phi))*uh/M2;
	me(2,1,0,1) = -exp(Complex(0.,-phi))*sh/M2;
	me(2,0,0,0) =  exp(Complex(0.,-phi))*uh/M2;
      }
      // Helicity code version
      // complex<InvEnergy3> fact = sqrt(-2/th)/M2*Complex(0.,1.);
      // for(unsigned int ih2=0;ih2<2;++ih2) {
      //   for(unsigned int ih4=0;ih4<2;++ih4) {
      //     LorentzPolarizationVectorE fCurrent = q2[ih2].dimensionedWave().vectorCurrent(q4[ih4].dimensionedWave());
      //     LorentzPolarizationVector eps = fact*Helicity::epsilon(fCurrent,hard[2]->momentum(),hard[0]->momentum());
      //     for(unsigned int ih1=0;ih1<2;++ih1) {
      // 	Complex amp = eps*g1[ih1].wave();
      // 	if(norm(me(2*ih1,ih2,0,ih4))>1e-10)  cerr << "testing in hel loop B " << ih1 << " " << ih2 << " " << ih4 << " "
      // 						  << amp << " " << me(2*ih1,ih2,0,ih4) << " " << amp/me(2*ih1,ih2,0,ih4)
      // 						  << " " << norm(amp/me(2*ih1,ih2,0,ih4)) << "\n";
      //     }
      //   }
      // }
    }
    else
      assert(false);
  }
  else if(hard[0]->id()==-hard[1]->id()) {
    vector<SpinorWaveFunction>    q1;
    vector<SpinorBarWaveFunction> q2;
    vector<VectorWaveFunction>    g4;
    SpinorWaveFunction(   q1,hard[0],incoming,false,true);
    SpinorBarWaveFunction(q2,hard[1],incoming,false,true);
    VectorWaveFunction(   g4,hard[3],outgoing,true,true,true,vector_phase);
    g4[1]=g4[2];
    // matrix element
    me = ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0,PDT::Spin1);
    if(!swapped) {
      me(0,1,0,0) = -th/M2;
      me(0,1,0,2) =  exp(Complex(0.,-2.*phi))*uh/M2;
      me(1,0,0,0) =  exp(Complex(0., 2.*phi))*uh/M2;
      me(1,0,0,2) = -th/M2;
    }
    else {
      me(0,1,0,0) =  exp(Complex(0., 2.*phi))*th/M2;
      me(0,1,0,2) = -uh/M2;
      me(1,0,0,0) = -uh/M2;
      me(1,0,0,2) =  exp(Complex(0.,-2.*phi))*th/M2;
    }
    // Helicity code version
    // complex<InvEnergy3> fact = sqrt(2./sh)/M2*Complex(0.,1.);
    // for(unsigned int ih1=0;ih1<2;++ih1) {
    //   for(unsigned int ih2=0;ih2<2;++ih2) {
    // 	LorentzPolarizationVectorE fCurrent = q1[ih1].dimensionedWave().vectorCurrent(q2[ih2].dimensionedWave());
    // 	LorentzPolarizationVector eps = fact*Helicity::epsilon(fCurrent,hard[2]->momentum(),hard[3]->momentum());
    // 	for(unsigned int ih4=0;ih4<2;++ih4) {
    // 	  Complex amp = eps*g4[ih4].wave();
    // 	  //if(norm(me(ih1,ih2,0,2*ih4))>1e-10)
    // 	    cerr << "testing in hel loop B " << ih1 << " " << ih2 << " " << ih4 << " "
    // 		 << amp << " " << me(ih1,ih2,0,2*ih4) << " " << amp/me(ih1,ih2,0,2*ih4)
    // 		 << " " << norm(amp/me(ih1,ih2,0,2*ih4)) << "\n";
    // 	}
    //   }
    // }
  }
  else
    assert(false);
  // construct the vertex
  HardVertexPtr hardvertex = new_ptr(HardVertex());
  // // set the matrix element for the vertex
  hardvertex->ME(me);
  // set the pointers and to and from the vertex
  for(unsigned int i = 0; i < hard.size(); ++i)
    hard[i]->spinInfo()->productionVertex(hardvertex);
  // boost back to lab
  boost = LorentzRotation(pcms.boostVector());
  for(PPtr part : hard)
    part->transform(boost);
}