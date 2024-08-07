// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPPto3D1Jet class.
//

#include "MEPPto3D1Jet.h"
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
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

void MEPPto3D1Jet::doinit() {
  HwMEBase::doinit();
  // get the non-perturbative ME
  O1_ = params_->singletMEProduction<2>(state_,n_,1,1);
  // set the mass option
  massOption(vector<unsigned int>({mOpt_+1,0}));
}

IBPtr MEPPto3D1Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPPto3D1Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPPto3D1Jet::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2*GeV2*GeV2) << oenum(state_) << n_ << mOpt_;
}

void MEPPto3D1Jet::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2*GeV2*GeV2) >> ienum(state_) >> n_ >> mOpt_;
}

//The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPPto3D1Jet,HwMEBase>
describeHerwigMEPPto3D1Jet("Herwig::MEPPto3D1Jet",
			     "HwOniumParameters.so HwMEHadronOnium.so");

void MEPPto3D1Jet::Init() {

  static ClassDocumentation<MEPPto3D1Jet> documentation
    ("The MEPPto3D1Jet class implements the g g to 3D1 g processes");

  static Reference<MEPPto3D1Jet,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &MEPPto3D1Jet::params_, false, false, true, false, false);
  
  static Switch<MEPPto3D1Jet,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &MEPPto3D1Jet::state_, ccbar, false, false);
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
  
  static Parameter<MEPPto3D1Jet,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &MEPPto3D1Jet::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Switch<MEPPto3D1Jet,unsigned int> interfaceMassOption
    ("MassOption",
     "Mass of the treatment of the 3D1 mass",
     &MEPPto3D1Jet::mOpt_, 0, false, false);
  static SwitchOption interfaceMassOptionOnShell
    (interfaceMassOption,
     "OnShell",
     "Use the on-shell mass",
     0);
  static SwitchOption interfaceMassOptionOffShell
    (interfaceMassOption,
     "OffShell",
     "Use an off-shell mass generated by the MassGenerator object for the 3D1 state.",
     1);
}

void MEPPto3D1Jet::getDiagrams() const {
  // construct the meson PDG code from quark ids
  unsigned int iq = 4+state_;
  tcPDPtr ps = getParticleData(long(iq*110 + 30003 + (n_-1)*100000));
  tcPDPtr g = getParticleData(ParticleID::g);
  add(new_ptr((Tree2toNDiagram(2), g, g, 1, ps, 1, g , -1)));
}

Selector<MEBase::DiagramIndex>
MEPPto3D1Jet::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ) sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *>
MEPPto3D1Jet::colourGeometries(tcDiagPtr ) const {
  static ColourLines c1("1 -2,  2  4, -1 -4");
  static ColourLines c2("1  4, -4 -2,  2 -1");
  Selector<const ColourLines *> sel;
  sel.insert(0.5, &c1);
  sel.insert(0.5, &c2);
  return sel;
}

Energy2 MEPPto3D1Jet::scale() const {
  return sHat();
}

double MEPPto3D1Jet::me2() const {
  Energy M = meMomenta()[2].mass();
  Energy2 M2=sqr(M);
  Energy4 um(sqr(uHat()-M2)),tm(sqr(tHat()-M2)),sm(sqr(sHat()-M2));
  double output = 256*O1_*pow(Constants::pi*standardModel()->alphaS(scale()),3)/
    (243.*pow<3,1>(M)*pow<5,1>(M2-tHat())*pow<5,1>(M2-uHat())*pow<5,1>(tHat()+uHat()))*
    ((102*pow<10,1>(M2)*pow<3,1>(tHat())-286*pow<9,1>(M2)*pow<4,1>(tHat())+
      275*pow<8,1>(M2)*pow<5,1>(tHat())-227*pow<7,1>(M2)*pow<6,1>(tHat())+
      410*pow<6,1>(M2)*pow<7,1>(tHat())-470*pow<5,1>(M2)*pow<8,1>(tHat())+
      245*pow<4,1>(M2)*pow<9,1>(tHat())-49*pow<3,1>(M2)*pow<10,1>(tHat())+
      302*pow<10,1>(M2)*sqr(tHat())*uHat()-1732*pow<9,1>(M2)*pow<3,1>(tHat())*uHat()+
      3840*pow<8,1>(M2)*pow<4,1>(tHat())*uHat()-5004*pow<7,1>(M2)*pow<5,1>(tHat())*uHat()+
      5137*pow<6,1>(M2)*pow<6,1>(tHat())*uHat()-4220*pow<5,1>(M2)*pow<7,1>(tHat())*uHat()+
      2190*pow<4,1>(M2)*pow<8,1>(tHat())*uHat()-580*pow<3,1>(M2)*pow<9,1>(tHat())*uHat()+
      67*sqr(M2)*pow<10,1>(tHat())*uHat()+302*pow<10,1>(M2)*tHat()*sqr(uHat())-
      2844*pow<9,1>(M2)*sqr(tHat())*sqr(uHat())+10289*pow<8,1>(M2)*pow<3,1>(tHat())*sqr(uHat())-
      19569*pow<7,1>(M2)*pow<4,1>(tHat())*sqr(uHat())+23585*pow<6,1>(M2)*pow<5,1>(tHat())*sqr(uHat())-
      19534*pow<5,1>(M2)*pow<6,1>(tHat())*sqr(uHat())+10358*pow<4,1>(M2)*pow<7,1>(tHat())*sqr(uHat())-
      2822*pow<3,1>(M2)*pow<8,1>(tHat())*sqr(uHat())+210*sqr(M2)*pow<9,1>(tHat())*sqr(uHat())+
      25*M2*pow<10,1>(tHat())*sqr(uHat())+102*pow<10,1>(M2)*pow<3,1>(uHat())-
      1732*pow<9,1>(M2)*tHat()*pow<3,1>(uHat())+10289*pow<8,1>(M2)*sqr(tHat())*pow<3,1>(uHat())-
      29536*pow<7,1>(M2)*pow<3,1>(tHat())*pow<3,1>(uHat())+
      47908*pow<6,1>(M2)*pow<4,1>(tHat())*pow<3,1>(uHat())-
      47528*pow<5,1>(M2)*pow<5,1>(tHat())*pow<3,1>(uHat())+
      28602*pow<4,1>(M2)*pow<6,1>(tHat())*pow<3,1>(uHat())-
      8984*pow<3,1>(M2)*pow<7,1>(tHat())*pow<3,1>(uHat())+774*sqr(M2)*pow<8,1>(tHat())*pow<3,1>(uHat())+
      100*M2*pow<9,1>(tHat())*pow<3,1>(uHat())+5*pow<10,1>(tHat())*pow<3,1>(uHat())-
      286*pow<9,1>(M2)*pow<4,1>(uHat())+3840*pow<8,1>(M2)*tHat()*pow<4,1>(uHat())-
      19569*pow<7,1>(M2)*sqr(tHat())*pow<4,1>(uHat())+47908*pow<6,1>(M2)*pow<3,1>(tHat())*pow<4,1>(uHat())-
      63536*pow<5,1>(M2)*pow<4,1>(tHat())*pow<4,1>(uHat())+
      47093*pow<4,1>(M2)*pow<5,1>(tHat())*pow<4,1>(uHat())-
      17653*pow<3,1>(M2)*pow<6,1>(tHat())*pow<4,1>(uHat())+2006*sqr(M2)*pow<7,1>(tHat())*pow<4,1>(uHat())+
      220*M2*pow<8,1>(tHat())*pow<4,1>(uHat())+25*pow<9,1>(tHat())*pow<4,1>(uHat())+
      275*pow<8,1>(M2)*pow<5,1>(uHat())-5004*pow<7,1>(M2)*tHat()*pow<5,1>(uHat())+
      23585*pow<6,1>(M2)*sqr(tHat())*pow<5,1>(uHat())-47528*pow<5,1>(M2)*pow<3,1>(tHat())*pow<5,1>(uHat())+
      47093*pow<4,1>(M2)*pow<4,1>(tHat())*pow<5,1>(uHat())-
      21968*pow<3,1>(M2)*pow<5,1>(tHat())*pow<5,1>(uHat())+3147*sqr(M2)*pow<6,1>(tHat())*pow<5,1>(uHat())+
      340*M2*pow<7,1>(tHat())*pow<5,1>(uHat())+60*pow<8,1>(tHat())*pow<5,1>(uHat())-
      227*pow<7,1>(M2)*pow<6,1>(uHat())+5137*pow<6,1>(M2)*tHat()*pow<6,1>(uHat())-
      19534*pow<5,1>(M2)*sqr(tHat())*pow<6,1>(uHat())+28602*pow<4,1>(M2)*pow<3,1>(tHat())*pow<6,1>(uHat())-
      17653*pow<3,1>(M2)*pow<4,1>(tHat())*pow<6,1>(uHat())+3147*sqr(M2)*pow<5,1>(tHat())*pow<6,1>(uHat())+
      390*M2*pow<6,1>(tHat())*pow<6,1>(uHat())+90*pow<7,1>(tHat())*pow<6,1>(uHat())+
      410*pow<6,1>(M2)*pow<7,1>(uHat())-4220*pow<5,1>(M2)*tHat()*pow<7,1>(uHat())+
      10358*pow<4,1>(M2)*sqr(tHat())*pow<7,1>(uHat())-8984*pow<3,1>(M2)*pow<3,1>(tHat())*pow<7,1>(uHat())+
      2006*sqr(M2)*pow<4,1>(tHat())*pow<7,1>(uHat())+340*M2*pow<5,1>(tHat())*pow<7,1>(uHat())+
      90*pow<6,1>(tHat())*pow<7,1>(uHat())-470*pow<5,1>(M2)*pow<8,1>(uHat())+
      2190*pow<4,1>(M2)*tHat()*pow<8,1>(uHat())-2822*pow<3,1>(M2)*sqr(tHat())*pow<8,1>(uHat())+
      774*sqr(M2)*pow<3,1>(tHat())*pow<8,1>(uHat())+220*M2*pow<4,1>(tHat())*pow<8,1>(uHat())+
      60*pow<5,1>(tHat())*pow<8,1>(uHat())+245*pow<4,1>(M2)*pow<9,1>(uHat())-
      580*pow<3,1>(M2)*tHat()*pow<9,1>(uHat())+210*sqr(M2)*sqr(tHat())*pow<9,1>(uHat())+
      100*M2*pow<3,1>(tHat())*pow<9,1>(uHat())+25*pow<4,1>(tHat())*pow<9,1>(uHat())-
      49*pow<3,1>(M2)*pow<10,1>(uHat())+67*sqr(M2)*tHat()*pow<10,1>(uHat())+
      25*M2*sqr(tHat())*pow<10,1>(uHat())+5*pow<3,1>(tHat())*pow<10,1>(uHat())));
  // // test vsPRD 45, 116 PR45, 116
  // Energy7 R02 = params_->secondDerivativeRadialWaveFunctionSquared(state_,n_);
  // Energy6 Q(sHat()*tHat()*uHat());
  // Energy4 P(sHat()*tHat()+tHat()*uHat()+uHat()*sHat());
  // double test = 16.*Constants::pi*20.*Constants::pi*pow(standardModel()->alphaS(scale()),3)*R02/
  //   (9.*pow<3,1>(M)*pow<5,1>(Q-M2*P))*
  //   (-49*pow<3,1>(M2)*pow<5,1>(P) - 
  //     sqr(M2)*pow<4,1>(P)*(20*pow<3,1>(M2) - 67*Q) - 
  //    52*pow<7,1>(M2)*sqr(Q) + 
  //     894*pow<4,1>(M2)*pow<3,1>(Q) - 
  //     5*M2*pow<4,1>(Q) + 
  //     sqr(M2)*P*Q*(4*pow<6,1>(M2) - 
  // 		   1742*pow<3,1>(M2)*Q - 361*sqr(Q)) - 
  //    M2*pow<3,1>(P)*
  //     (102*pow<6,1>(M2) + 295*pow<3,1>(M2)*Q - 
  //      25*sqr(Q)) + 
  //     sqr(P)*Q*(902*pow<6,1>(M2) + 729*pow<3,1>(M2)*Q + 
  // 		    5*sqr(Q)));
  // cerr << "testing matrix element " << output << " " << test << " "
  //      << (output-test)/(output+test) << " " << output/test << "\n";
  return output;
}

void MEPPto3D1Jet::constructVertex(tSubProPtr sub) {
  using namespace ThePEG::Helicity;
  // extract the particles in the hard process
  // only one order
  ParticleVector hard;
  hard.reserve(4);
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // boost to partonic CMS
  Lorentz5Momentum pcms = hard[0]->momentum()+hard[1]->momentum();
  LorentzRotation boost(-pcms.boostVector());
  for(PPtr part : hard) part->transform(boost);
  // set the wavefunctions
  vector<VectorWaveFunction> g1,g2,psi,g4;
  VectorWaveFunction( g1,hard[0],incoming,false, true,true,vector_phase);
  VectorWaveFunction( g2,hard[1],incoming,false, true,true,vector_phase);
  VectorWaveFunction(psi,hard[2],outgoing,true ,false,true,vector_phase);
  VectorWaveFunction( g4,hard[3],outgoing,true , true,true,vector_phase);
  // extract kinematic variables
  Energy M = hard[2]->mass();
  Energy2 M2 = sqr(hard[2]->mass());
  double phi = hard[2]->momentum().phi();
  Energy2 sh = (hard[0]->momentum()+hard[1]->momentum()).m2();
  Energy2 th = (hard[0]->momentum()-hard[2]->momentum()).m2();
  Energy2 uh = (hard[0]->momentum()-hard[3]->momentum()).m2();
  Energy2 um(uh-M2),tm(th-M2),sm(sh-M2);
  Complex phase = exp(Complex(0.,phi));
  // Energy rstu = sqrt(th*uh/sh);
  // calculate the matrix element
  ProductionMatrixElement me(PDT::Spin1,PDT::Spin1,PDT::Spin1,PDT::Spin1);
  me(0,0,0,0)=(2*M2*sqr(phase)*sh*(10*pow<5,1>(M2)*sm+6*M2*pow<3,1>(sm)*th*uh+sqr(sm)*sqr(th)*sqr(uh)-6*pow<3,1>(M2)*sm*(sqr(th)-13*th*uh+sqr(uh))+pow<4,1>(M2)*(11*sqr(th)+54*th*uh+11*sqr(uh))+sqr(M2)*(-7*pow<4,1>(th)+34*pow<3,1>(th)*uh+66*sqr(th)*sqr(uh)+34*th*pow<3,1>(uh)-7*pow<4,1>(uh))))/(sqrt(15)*sqr(sm)*pow<3,1>(tm)*pow<3,1>(um));
  me(0,0,0,2)=(-4*sqrt(1.6666666666666667)*pow<3,1>(M2)*(M2+sh)*th*uh)/(sqr(sm)*sqr(tm)*sqr(um));
  me(0,0,1,0)=(-4*sqrt(0.13333333333333333)*pow<3,1>(M)*phase*sh*sqrt(sh)*(th-uh)*sqrt(th*uh)*(3*sqr(M2)*(sh-uh)+sh*uh*(sh+uh)+M2*(5*sqr(sh)+2*sh*uh+3*sqr(uh))))/(sqr(sm)*pow<3,1>(tm)*pow<3,1>(um));
  me(0,0,1,2)=(4*sqrt(3.3333333333333335)*pow<7,1>(M)*sqrt(sh)*(th-uh)*sqrt(th*uh))/(phase*sqr(sm)*sqr(tm)*sqr(um));
  me(0,0,2,0)=(-4*M2*sqr(sh)*th*uh*(16*pow<3,1>(M2)+15*sqr(M2)*sm-sm*th*uh+M2*(5*sqr(th)+2*th*uh+5*sqr(uh))))/(sqrt(15)*sqr(sm)*pow<3,1>(tm)*pow<3,1>(um));
  me(0,0,2,2)=(4*sqrt(1.6666666666666667)*pow<3,1>(M2)*sh*(M2*sm+2*th*uh))/(sqr(phase)*sqr(sm)*sqr(tm)*sqr(um));
  me(0,2,0,0)=(-2*M2*th*(2*pow<4,1>(M2)*(sh+th)+2*sqr(sh)*sqr(th)*(sh+th)+3*pow<3,1>(M2)*(5*sqr(sh)+6*sh*th+5*sqr(th))+M2*sh*th*(7*sqr(sh)+2*sh*th+7*sqr(th))-sqr(M2)*(sh+th)*(15*sqr(sh)+10*sh*th+7*sqr(th)))*uh)/(sqrt(15)*pow<3,1>(sm)*sqr(tm)*pow<3,1>(um));
  me(0,2,0,2)=(-2*M2*sh*sqr(uh)*(2*pow<4,1>(M2)-sm*sqr(th)*uh+M2*th*(th-uh)*(5*th+6*uh)+pow<3,1>(M2)*(13*th+15*uh)-sqr(M2)*(8*sqr(th)+17*th*uh+7*sqr(uh))))/(sqrt(15)*sqr(phase)*pow<3,1>(sm)*pow<3,1>(tm)*sqr(um));
  me(0,2,1,0)=(2*sqrt(0.13333333333333333)*pow<3,1>(M)*sqrt(sh)*th*(pow<3,1>(M2)*(8*sh-4*th)-sqr(M2)*(3*sh-11*th)*(sh+th)-sh*th*(sh+th)*(3*sh+th)+M2*(-3*pow<3,1>(sh)+sqr(sh)*th-5*sh*sqr(th)+3*pow<3,1>(th)))*sqrt(th*uh))/(phase*pow<3,1>(sm)*sqr(tm)*pow<3,1>(um));
  me(0,2,1,2)=(2*sqrt(0.13333333333333333)*pow<3,1>(M)*sqrt(sh)*uh*sqrt(th*uh)*(pow<3,1>(M2)*(8*sh-4*uh)-sqr(M2)*(3*sh-11*uh)*(sh+uh)-sh*uh*(sh+uh)*(3*sh+uh)+M2*(-3*pow<3,1>(sh)+sqr(sh)*uh-5*sh*sqr(uh)+3*pow<3,1>(uh))))/(pow(phase,3)*pow<3,1>(sm)*pow<3,1>(tm)*sqr(um));
  me(0,2,2,0)=(-2*M2*sh*sqr(th)*(2*pow<4,1>(M2)-sm*th*sqr(uh)+pow<3,1>(M2)*(15*th+13*uh)+M2*uh*(-6*sqr(th)+th*uh+5*sqr(uh))-sqr(M2)*(7*sqr(th)+17*th*uh+8*sqr(uh))))/(sqrt(15)*sqr(phase)*pow<3,1>(sm)*sqr(tm)*pow<3,1>(um));
  me(0,2,2,2)=(-2*M2*th*uh*(2*pow<4,1>(M2)*(sh+uh)+2*sqr(sh)*sqr(uh)*(sh+uh)+3*pow<3,1>(M2)*(5*sqr(sh)+6*sh*uh+5*sqr(uh))+M2*sh*uh*(7*sqr(sh)+2*sh*uh+7*sqr(uh))-sqr(M2)*(sh+uh)*(15*sqr(sh)+10*sh*uh+7*sqr(uh))))/(sqrt(15)*pow(phase,4)*pow<3,1>(sm)*pow<3,1>(tm)*sqr(um));
  me(2,0,0,0)=(-2*M2*pow(phase,4)*th*uh*(2*pow<4,1>(M2)*(sh+uh)+2*sqr(sh)*sqr(uh)*(sh+uh)+3*pow<3,1>(M2)*(5*sqr(sh)+6*sh*uh+5*sqr(uh))+M2*sh*uh*(7*sqr(sh)+2*sh*uh+7*sqr(uh))-sqr(M2)*(sh+uh)*(15*sqr(sh)+10*sh*uh+7*sqr(uh))))/(sqrt(15)*pow<3,1>(sm)*pow<3,1>(tm)*sqr(um));
  me(2,0,0,2)=(-2*M2*sqr(phase)*sh*sqr(th)*(2*pow<4,1>(M2)-sm*th*sqr(uh)+pow<3,1>(M2)*(15*th+13*uh)+M2*uh*(-6*sqr(th)+th*uh+5*sqr(uh))-sqr(M2)*(7*sqr(th)+17*th*uh+8*sqr(uh))))/(sqrt(15)*pow<3,1>(sm)*sqr(tm)*pow<3,1>(um));
  me(2,0,1,0)=(-2*sqrt(0.13333333333333333)*pow<3,1>(M)*pow(phase,3)*sqrt(sh)*uh*sqrt(th*uh)*(pow<3,1>(M2)*(8*sh-4*uh)-sqr(M2)*(3*sh-11*uh)*(sh+uh)-sh*uh*(sh+uh)*(3*sh+uh)+M2*(-3*pow<3,1>(sh)+sqr(sh)*uh-5*sh*sqr(uh)+3*pow<3,1>(uh))))/(pow<3,1>(sm)*pow<3,1>(tm)*sqr(um));
  me(2,0,1,2)=(-2*sqrt(0.13333333333333333)*pow<3,1>(M)*phase*sqrt(sh)*th*(pow<3,1>(M2)*(8*sh-4*th)-sqr(M2)*(3*sh-11*th)*(sh+th)-sh*th*(sh+th)*(3*sh+th)+M2*(-3*pow<3,1>(sh)+sqr(sh)*th-5*sh*sqr(th)+3*pow<3,1>(th)))*sqrt(th*uh))/(pow<3,1>(sm)*sqr(tm)*pow<3,1>(um));
  me(2,0,2,0)=(-2*M2*sqr(phase)*sh*sqr(uh)*(2*pow<4,1>(M2)-sm*sqr(th)*uh+M2*th*(th-uh)*(5*th+6*uh)+pow<3,1>(M2)*(13*th+15*uh)-sqr(M2)*(8*sqr(th)+17*th*uh+7*sqr(uh))))/(sqrt(15)*pow<3,1>(sm)*pow<3,1>(tm)*sqr(um));
  me(2,0,2,2)=(-2*M2*th*(2*pow<4,1>(M2)*(sh+th)+2*sqr(sh)*sqr(th)*(sh+th)+3*pow<3,1>(M2)*(5*sqr(sh)+6*sh*th+5*sqr(th))+M2*sh*th*(7*sqr(sh)+2*sh*th+7*sqr(th))-sqr(M2)*(sh+th)*(15*sqr(sh)+10*sh*th+7*sqr(th)))*uh)/(sqrt(15)*pow<3,1>(sm)*sqr(tm)*pow<3,1>(um));
  me(2,2,0,0)=(4*sqrt(1.6666666666666667)*pow<3,1>(M2)*sqr(phase)*sh*(M2*sm+2*th*uh))/(sqr(sm)*sqr(tm)*sqr(um));
  me(2,2,0,2)=(-4*M2*sqr(sh)*th*uh*(16*pow<3,1>(M2)+15*sqr(M2)*sm-sm*th*uh+M2*(5*sqr(th)+2*th*uh+5*sqr(uh))))/(sqrt(15)*sqr(sm)*pow<3,1>(tm)*pow<3,1>(um));
  me(2,2,1,0)=(4*sqrt(3.3333333333333335)*pow<7,1>(M)*phase*sqrt(sh)*sqrt(th*uh)*(-th+uh))/(sqr(sm)*sqr(tm)*sqr(um));
  me(2,2,1,2)=(4*sqrt(0.13333333333333333)*pow<3,1>(M)*sh*sqrt(sh)*(th-uh)*sqrt(th*uh)*(3*sqr(M2)*(sh-uh)+sh*uh*(sh+uh)+M2*(5*sqr(sh)+2*sh*uh+3*sqr(uh))))/(phase*sqr(sm)*pow<3,1>(tm)*pow<3,1>(um));
  me(2,2,2,0)=(-4*sqrt(1.6666666666666667)*pow<3,1>(M2)*(M2+sh)*th*uh)/(sqr(sm)*sqr(tm)*sqr(um));
  me(2,2,2,2)=(2*M2*sh*(10*pow<5,1>(M2)*sm+6*M2*pow<3,1>(sm)*th*uh+sqr(sm)*sqr(th)*sqr(uh)-6*pow<3,1>(M2)*sm*(sqr(th)-13*th*uh+sqr(uh))+pow<4,1>(M2)*(11*sqr(th)+54*th*uh+11*sqr(uh))+sqr(M2)*(-7*pow<4,1>(th)+34*pow<3,1>(th)*uh+66*sqr(th)*sqr(uh)+34*th*pow<3,1>(uh)-7*pow<4,1>(uh))))/(sqrt(15)*sqr(phase)*sqr(sm)*pow<3,1>(tm)*pow<3,1>(um));
  // test the spin averaged result
  // Energy6 Q(sh*th*uh);
  // Energy4 P(sh*th+th*uh+uh*sh);
  // double test = (16*sqr(M2)*(-49*pow<3,1>(M2)*pow<5,1>(P)-sqr(M2)*pow<4,1>(P)*(20*pow<3,1>(M2)-67*Q)-52*pow<7,1>(M2)*sqr(Q)+894*pow<4,1>(M2)*pow<3,1>(Q)-5*M2*pow<4,1>(Q)+sqr(M2)*P*Q*(4*pow<6,1>(M2)-1742*pow<3,1>(M2)*Q-361*sqr(Q))-M2*pow<3,1>(P)*(102*pow<6,1>(M2)+295*pow<3,1>(M2)*Q-25*sqr(Q))+sqr(P)*Q*(902*pow<6,1>(M2)+729*pow<3,1>(M2)*Q+5*sqr(Q))))/(15.*pow<5,1>(Q-M2*P));
  
  // double aver = me.average();
  // cerr << "testing spin correlations " << test << " " << me.average() << " "
  //      << abs(test-aver)/(test+aver) << "\n";
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
