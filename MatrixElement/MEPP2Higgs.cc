// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2Higgs class.
//

#include "MEPP2Higgs.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "GeneralHardME.h"
#include "HardVertex.h"

using namespace Herwig;

ClassDescription<MEPP2Higgs> MEPP2Higgs::initMEPP2Higgs;
// Definition of the static class description member.

void MEPP2Higgs::persistentOutput(PersistentOStream & os) const {
  os << hggvertex << ffhvertex << theSM << branchingopt 
     << shapeopt << widthopt << processopt;
}

void MEPP2Higgs::persistentInput(PersistentIStream & is, int) {
  is >> hggvertex >> ffhvertex >> theSM >> branchingopt 
     >> shapeopt >> widthopt >> processopt;
}
void MEPP2Higgs::Init() {

  static ClassDocumentation<MEPP2Higgs> documentation
    ("The MEPP2Higgs class implements the matrix elements for"
     " Higgs production (with decay H->W-W+) in hadron-hadron collisions.");

//////////////////////////////////////////////////////////////////////////////////////////
  static Switch<MEPP2Higgs,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &MEPP2Higgs::shapeopt, 1, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonanse",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "RunningBreitWigner",
     "Breit-Wigner Higgs resonanse with internal running width",
     2);
  static SwitchOption interfaceImprovedShape
    (interfaceShapeOption,
     "ImprovedBreitWigner",
     "Improved Higgs s-channel resonanse (hep-ph/9505211)",
     6);

//////////////////////////////////////////////////////////////////////////////////////////
  static Switch<MEPP2Higgs,unsigned int> interfaceWidthOption
    ("WidthScheme",
     "Option for the treatment of the Higss Width calculation",
     &MEPP2Higgs::widthopt, 1, false, false);
  static SwitchOption interfaceFixedWidth
    (interfaceWidthOption,
     "FixedHiggsWidth",
     "Fixed Higgs width, taken from ThePEGParticles.in",
     1);
  static SwitchOption interfaceRunningWidth
    (interfaceWidthOption,
     "RunningHiggsWidth",
     "Running Higgs width, according to Breit-Wigner (does not work...)",
     2);
  static SwitchOption interfaceNLLWidth
    (interfaceWidthOption,
     "NLLcorrectedHiggsWidth",
     "NLL corrected Higgs width (a-la FORTRAN HERWIG)",
     3);
  static SwitchOption interfaceLOWidthOption
    (interfaceWidthOption,
     "LOHiggsWidth",
     "LO Higgs width (formula taken from The \"Higgs Hunter's Guide\")",
     4);
  static SwitchOption interfaceWidthTestOption
    (interfaceWidthOption,
     "HiggsWidthTests",
     "Test option",
     5);

//////////////////////////////////////////////////////////////////////////////////////////
  static Switch<MEPP2Higgs,unsigned int> interfaceBranchingOption
    ("BranchAllowed",
     "Option to switch on/off branchings in the total Higgs width calculation",
     &MEPP2Higgs::branchingopt, 3, false, false);
  static SwitchOption interfaceFermionOnlyBranchings
    (interfaceBranchingOption,
     "FermionBR",
     "Fermion branchings in the full Higgs width",
     1);
  static SwitchOption interfacePlusWWZZBranchings
    (interfaceBranchingOption,
     "FermionWZBR",
     "Fermion and WW/ZZ branchings in the full Higgs width",
     2);
  static SwitchOption interfacePlusGammaBranching
    (interfaceBranchingOption,
     "FermionWZGammaBR",
     "Fermion 2gamma, and WW/ZZ branchings in the full Higgs width",
     3);
  static SwitchOption interfacePlausGluonBranching
    (interfaceBranchingOption,
     "FermionWZGammaGluonBR",
     "Fermion 2gamma, 2glions, and WW/ZZ branchings in the full Higgs width",
     4);

//////////////////////////////////////////////////////////////////////////////////////////
  static Switch<MEPP2Higgs,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2Higgs::processopt, 1, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     1);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "qqbar",
     "Only include the incoming q qbar subprocess",
     2);
  static SwitchOption interfaceProcessgg
    (interfaceProcess,
     "gg",
     "Only include the incoming gg subprocess",
     3);

//////////////////////////////////////////////////////////////////////////////////////////
  static Parameter<MEPP2Higgs,unsigned int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The minimum flavour of the incoming quarks in the hard process",
     &MEPP2Higgs::minflavouropt, 3, 3, 5,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the incoming quarks in the hard process",
     &MEPP2Higgs::maxflavouropt, 5, 3, 5,
     false, false, Interface::limited);
}

void MEPP2Higgs::doinit() throw(InitException) {
  MEBase::doinit();
  // get the vertex pointers from the SM object
  theSM = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!theSM) {
    throw InitException() << "Wrong type of StandardModel object in MEPP2Higgs::doinit(),"
                          << " the Herwig++ version must be used" 
                          << Exception::runerror;
  }
  hggvertex = dynamic_ptr_cast<SVVLoopVertexPtr>(theSM->vertexHGG());
//  hggvertex = dynamic_ptr_cast<GeneralSVVVertexPtr>(theSM->vertexHGG());
  ffhvertex = theSM->vertexFFH();
  PDPtr h0 = getParticleData(ParticleID::h0);
}

unsigned int MEPP2Higgs::orderInAlphaS() const {
  return 2;
}

unsigned int MEPP2Higgs::orderInAlphaEW() const {
  return 1;
}

Energy2 MEPP2Higgs::scale() const {
  return sHat();
}

int MEPP2Higgs::nDim() const {
  return 0;
}


bool MEPP2Higgs::generateKinematics(const double *) {
  Lorentz5Momentum pout = meMomenta()[0] + meMomenta()[1];
  pout.rescaleMass();
  meMomenta()[2].setMass(pout.mass());
  meMomenta()[2] = LorentzMomentum(pout.x(),pout.y(),pout.z(),pout.t());
  jacobian(1.0);

  Energy2 s(sHat());
  PDPtr h0 = getParticleData(ParticleID::h0);
  if (s < 0.*GeV2 && widthopt != 1) {
    throw MEException() << "Warning! Shat < 0.0 (" << s/GeV2 
                        << "), so I can not use widthopt != 1" 
                        << Exception::warning;
  } else {
    switch (widthopt) {
    case 1:
      break;
    case 2: {
      Energy runwidth = h0->generateWidth(sqrt(s));
      h0->width(runwidth);
      break;}
    case 3: {
      Energy runwidth = calcNLLRunningWidth(sqrt(s));
      h0->width(runwidth);
      break;
    }
    case 4: {
      Energy runwidth = calcLORunningWidth(sqrt(s));
      h0->width(runwidth);
      break;
    }
    case 5: {
      Energy fixwidth1 = h0->width();
      Energy runwidth1 = h0->generateWidth(sqrt(s));
      Energy runwidth2 = calcNLLRunningWidth(sqrt(s));
      cout << "Test: MEPP2Higgs::generateKinematics --> mass = " << sqrt(s)/GeV 
                                                       << ": fw = " << fixwidth1/GeV 
                                                       << ", rw1 = " << runwidth1/GeV 
                                                       << ", rw2 = " << runwidth2/GeV << endl;
      break;
    }
    default: 
      throw MEException() << "Unknown width option in MEPP2Higgs::generateKinematics" 
                          << Exception::runerror;
      break;
    }
  }
//  cout << "Higgs width = " << ounit(h0->width(),GeV) << endl;

  // check whether it passes all the cuts: returns true if it does
  vector<LorentzMomentum> out(1,meMomenta()[2]);
  tcPDVector tout(1,mePartonData()[2]);
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}


void MEPP2Higgs::getDiagrams() const {
  tcPDPtr h0=getParticleData(ParticleID::h0);
  // gg -> H process
  if(processopt==1||processopt==3) {
    tcPDPtr g=getParticleData(ParticleID::g);
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, h0, -1)));
  }
  // q qbar -> H processes
  if(processopt==1||processopt==2) {
    for (unsigned int i = minflavouropt; i <= maxflavouropt; ++i) {
      tcPDPtr q = getParticleData(i);
      tcPDPtr qb = q->CC();
    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, h0, -2)));
    }
  }
}


CrossSection MEPP2Higgs::dSigHatDR() const {
  tcPDPtr h0 = getParticleData(ParticleID::h0);
  double cs = me2() * jacobian() * h0->width()/sqrt(sHat());
  return UnitRemoval::InvE2 * sqr(hbarc) * cs;
}


double MEPP2Higgs::me2() const {
  double output(0.0);
  useMe();
  ScalarWaveFunction hout(meMomenta()[2],mePartonData()[2],outgoing);

// Safety code to garantee reliable behaviour of Higgs shape limits (important for heavy and broad Higgs resonance).
  Energy hmass = meMomenta()[2].m();
  PDPtr h0 = getParticleData(ParticleID::h0);
  Energy mass = h0->mass();
  Energy halfmass = .5*mass;
  if (.0*GeV > hmass) return 0.0;
// stricly speaking the condition is applicable if h0->widthUpCut() == h0->widthLoCut()...
  if (h0->widthLoCut() > halfmass) {
    if ((mass + h0->widthUpCut() < hmass || mass - h0->widthLoCut() > hmass)) return 0.0;
  } else {
    if (mass + halfmass < hmass || halfmass > hmass) return 0.0;
  }

  if (mePartonData()[0]->id() == ParticleID::g && mePartonData()[1]->id() == ParticleID::g) {
    VectorWaveFunction gin1(meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction gin2(meMomenta()[1],mePartonData()[1],incoming);

    vector<VectorWaveFunction> g1,g2;
    for(unsigned int i = 0; i < 2; ++i) {
      gin1.reset(2*i);
      g1.push_back(gin1);
      gin2.reset(2*i);
      g2.push_back(gin2);
    }
    output = ggME(g1,g2,hout,false);
  } else {
    if (mePartonData()[0]->id() == -mePartonData()[1]->id()) {
      SpinorWaveFunction    qin (meMomenta()[0],mePartonData()[0],incoming);
      SpinorBarWaveFunction qbin(meMomenta()[1],mePartonData()[1],incoming);

      vector<SpinorWaveFunction> fin;
      vector<SpinorBarWaveFunction> ain;
      for (unsigned int i = 0; i < 2; ++i) {
        qin.reset(i);
        fin.push_back(qin);
        qbin.reset(i);
        ain.push_back(qbin);
      }
      output = qqME(fin,ain,hout,false);
    }
    else {
    throw MEException() << "Unknown subprocess in MEPP2Higgs::me2()" 
                        << Exception::runerror;
    }
  }
  return output;
}

Selector<MEBase::DiagramIndex> MEPP2Higgs::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for (DiagramIndex i = 0; i < diags.size(); ++i)
    sel.insert(1.0, i);
  return sel;
}


Selector<const ColourLines *> MEPP2Higgs::colourGeometries(tcDiagPtr diag) const {
  // colour lines
  static const ColourLines line1("1 -2,2 -1");
  static const ColourLines line2("1 -2");
  // select the colour flow
  Selector<const ColourLines *> sel;
  if (diag->id() == -1) {
    sel.insert(1.0, &line1);
  } else {
    sel.insert(1.0, &line2);
  }
  // return the answer
  return sel;
}


void MEPP2Higgs::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  if(hard[0]->id() < hard[1]->id()) {
    swap(hard[0],hard[1]);
  }
  // identify the process and calculate the matrix element
  if(hard[0]->id() == ParticleID::g && hard[1]->id() == ParticleID::g) {
    vector<VectorWaveFunction> g1,g2;
    vector<SpinorBarWaveFunction> q;
    vector<SpinorWaveFunction> qbar;
    VectorWaveFunction (g1,hard[0],incoming,false,true,true);
    VectorWaveFunction (g2,hard[1],incoming,false,true,true);
    ScalarWaveFunction hout(hard[2],outgoing,true,true);
    g1[1] = g1[2];
    g2[1] = g2[2];
    ggME(g1,g2,hout,true);
  } else {
    vector<SpinorWaveFunction>    q1;
    vector<SpinorBarWaveFunction> q2;
    SpinorWaveFunction    (q1,hard[0],incoming,false,true);
    SpinorBarWaveFunction (q2,hard[1],incoming,false,true);
    ScalarWaveFunction     hout(hard[2],outgoing,true,true);
    qqME(q1,q2,hout,true);
  }
  // construct the vertex
  HardVertexPtr hardvertex = new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int i = 0; i < 3; ++i) {
    dynamic_ptr_cast<SpinfoPtr>(hard[i]->spinInfo())->setProductionVertex(hardvertex);
  }
}


double MEPP2Higgs::ggME(vector<VectorWaveFunction> g1, 
                          vector<VectorWaveFunction> g2, 
                          ScalarWaveFunction &, 
                          bool calc) const {
  PDPtr h0 = getParticleData(ParticleID::h0);
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  Energy2 s(sHat());
  double me2(0.0);
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 2; ++j) {
      ScalarWaveFunction higgsWF = hggvertex->evaluate(s,shapeopt,h0,g1[i],g2[j]);
      Complex diag = higgsWF.wave();
      me2 += real(diag*conj(diag));
      if(calc) newme(2*i, 2*j, 0) = diag;
    }
  }
  if(calc) _me.reset(newme);
  // initial colour and spin factors: colour -> (8/64) and spin -> (1/4)
  return me2/(32.);
}


double MEPP2Higgs::qqME(vector<SpinorWaveFunction> & fin, 
                          vector<SpinorBarWaveFunction> & ain, 
                          ScalarWaveFunction &, 
                          bool calc) const {
  PDPtr h0 = getParticleData(ParticleID::h0);
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0);
  Energy2 s(scale());
  double me2(0.0);
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 2; ++j) {
      ScalarWaveFunction higgsWF = ffhvertex->evaluate(s,shapeopt,h0,fin[i],ain[j]);
      Complex diag = higgsWF.wave();
      me2+=real(diag*conj(diag));
      if(calc) newme(i, j, 0) = diag;
    }
  }
  if(calc) _me.reset(newme);
  // final colour/spin factors
  return 3.*me2/9./4.;
}


// Taken from HERWIG 6510 with some simplifications.
Energy MEPP2Higgs::calcNLLRunningWidth(Energy Mh) const {
  Energy2 q2 = sqr(Mh);
  Energy QCDLambda = theSM->LambdaQCD(q2);
//  double alphaEM   = theSM->alphaEM(q2);
  double alphaEM   = theSM->alphaEM();
  double alphaS    = theSM->alphaS(q2);
  double sw2       = theSM->sin2ThetaW();
  Energy partialw[14] = {0.0*GeV};

  tPDPtr w = getParticleData(ParticleID::Wplus);
  tPDPtr z = getParticleData(ParticleID::Z0);
  Energy wmnom = w->mass();
  Energy zmnom = z->mass();
  Energy wgnom = w->width ();
  Energy zgnom = z->width ();

  double ncolour = theSM->Nc();
  double Ca = ncolour;
  double Cf = (sqr(ncolour)-1.0)/(2.0*Ca);
  const unsigned int minflavour = 1;
  const unsigned int maxflavour = 6;
  double pi = Constants::pi;

// All calculation are being done for Monte-Carlo QCD Lambda, except Higgs width...
  double bcoeff4=(11.*Ca-10.)/(12.*pi);
  double kfac=Ca*(67./18.-sqr(pi)/6.)-25./9.;
  QCDLambda /= exp(kfac/(4.*pi*bcoeff4))/sqrt(2.);

// H->fermion pair
  Energy qmass[7];
  double nflavour=minflavour-1;
  for (unsigned int i = minflavour; i <= maxflavour; ++i) {
    tcPDPtr quark = getParticleData(i);
    qmass[i] = quark->mass();
    if (2.0*qmass[i] < Mh) nflavour+=1.;
  }
  double k1 = 5./sqr(pi);
  double k0 = 3./(4.*sqr(pi));
  double beta0 = (11.*Ca-2.0*nflavour)/3.;
  double beta1 = (34.*sqr(Ca)-(10.*Ca+6.*Cf)*nflavour)/3.;
  double gam0 = -8.;
  double gam1 = -404./3.+40.*nflavour/9.;
  double SClog = log(sqr(Mh/QCDLambda));
  double Cd = 1.+(k1/k0-2.*gam0+gam0*beta1/sqr(beta0)*log(SClog)+(gam0*beta1-gam1*beta0)/sqr(beta0))/(beta0*SClog);
  Energy2 GFermiINV = 8.*sw2*sqr(wmnom)/alphaEM;

// quarks: partialw[1-6]
  for (unsigned int i = minflavour; i <= maxflavour; ++i) {
    Energy mf = qmass[i];
    double xf = sqr(mf/Mh);
    if (mf > QCDLambda) mf *= pow(log(Mh/QCDLambda)/log(mf/QCDLambda),gam0/(2.0*beta0));
    if (xf < 0.25) partialw[i] = Ca*Mh*sqr(mf)*pow(1.0-4.0*xf,1.5)*Cd/GFermiINV;
  }

// leptons: partialw[7-9]
  for (unsigned int i = 0; i < 3; ++i) {
    tcPDPtr lepton = getParticleData(11+2*i);
    Energy mf = lepton->mass();
    double xf = sqr(mf/Mh);
    if (xf < 0.25) partialw[7+i] = Mh*sqr(mf)*pow(1.0-4.0*xf,1.5)/GFermiINV;
  }

// H->W*W*/Z*Z*: partialw[10,11]
  if (branchingopt > 1) {
    double xw = sqr(wmnom/Mh);
    double xz = sqr(zmnom/Mh);
    double yw = wmnom*wgnom/sqr(Mh);
    double yz = zmnom*zgnom/sqr(Mh);
    partialw[10] = Mh*Mh*Mh*HwDoubleBW(xw,yw)/2./GFermiINV;
    partialw[11] = Mh*Mh*Mh*HwDoubleBW(xz,yz)/4./GFermiINV;
  }

// H->gamma,gamma: partialw[12]
  if (branchingopt > 2) {
    double taut = sqr(2.0*qmass[ParticleID::t]/Mh);
    double tauw = sqr(2.0*wmnom/Mh);
    std::pair<double,double> ftaut = HwW2 (taut);
    std::pair<double,double> ftauw = HwW2 (tauw);
    double re = 4.0/3.0*(-2.0*taut*(1.0+(1.0-taut)*ftaut.first))+(2.0+3.0*tauw*(1+(2.0-tauw)*ftauw.first));
    double im = 4.0/3.0*(-2.0*taut*(    (1.0-taut)*ftaut.second))+(   3.0*tauw*(  (2.0-tauw)*ftauw.second));
    partialw[12] = sqr(alphaEM/pi)*Mh*Mh*Mh*(sqr(re)+ sqr(im))/32./GFermiINV;
  }

// H->gluon,gluon: partialw[13]
  if (branchingopt > 3) {
    double tau = sqr(2.0*qmass[ParticleID::t]/Mh);
    std::pair<double,double> ftau = HwW2 (tau);
    double re = 1+(1.0-tau)*ftau.first;
    double im =   (1.0-tau)*ftau.second;
    partialw[13] = sqr(alphaS/pi)*Mh*Mh*Mh*sqr(tau)*(sqr(re)+ sqr(im))/4./GFermiINV;
  }

  Energy higgswidth = Energy();
  for (unsigned int i = 1; i < 14; ++i) {
    higgswidth += partialw[i];
  }
  return higgswidth;
}

// Taken from HERWIG 6510.
double MEPP2Higgs::HwDoubleBW(double x, double y) const {
// Calculate the Double Breit-Wigner Integral
//  x=(mw/mh)**2, y=mw*gw/mh**2
  double limit = 0.425;
  double nbin = 25;
  double itgerl = 0.0;
  if (y < 0.0) return itgerl;
  if (x > limit) {
// Direct Integration
    double fac1 = 0.25/nbin;
    for (unsigned int i = 0; i < nbin; ++i) {
      double x1 = (i+0.5)*fac1;
      double fac2 = (sqr(1-sqrt(x1))-x1)/nbin;
      double sq = 1.0;
      int j = 0;
      while (j < nbin && 0.0 < sq) {
        double x2 = (j+0.5)*fac2+x1;
        sq = 1.0+x1*x1+x2*x2-2*(x1+x2+x1*x2);
        itgerl += 2.0*(sqr(1-x1-x2)+8.0*x1*x2)*sqrt(sq)/(sqr(x1-x)+y*y)*y/(sqr(x2-x)+y*y)*y*fac1*fac2;
        ++j;
      }
    }
  } else {
// Integration using tan theta substitution
    double th1low = atan ((0.0-x)/y);
    double th1high = atan ((1.0-x)/y);
    double fac1 = (th1high-th1low)/nbin;
    for (unsigned int i = 0; i < nbin; ++i) {
      double th1 = (i+0.5)*fac1+th1low;
      double x1 = y*tan(th1)+x;
      double x2max = min(x1,sqr(1-sqrt(x1)));
      double th2low = atan ((0-x)/y);
      double th2high = atan ((x2max-x)/y);
      double fac2 = (th2high-th2low)/nbin;
      double sq = 1.0;
      int j = 0;
      while (j < nbin && 0.0 < sq) {
        double th2 = (j+0.5)*fac2+th2low;
        double x2 = y*tan(th2)+x;
        double sq = 1.0+x1*x1+x2*x2-2*(x1+x2+x1*x2);
        itgerl += 2.0*(sqr(1-x1-x2)+8*x1*x2)*sqrt(sq)*fac1*fac2;
        ++j;
      }
    }
  }
  itgerl *= 1/sqr(Constants::pi);
  return itgerl;
}

// Taken from HERWIG 6510.
std::pair<double,double> MEPP2Higgs::HwW2(double tau) const {
  std::pair<double,double> ac;
  double pi = Constants::pi;
  if (tau > 1.0) {
    ac.first = sqr(asin(1.0/sqrt(tau)));
  } else if (tau < 1.0) {
    double FNsqr = sqrt(1-tau);
    double FNlog = log((1+FNsqr)/(1-FNsqr));
    ac.first = -0.25 * (sqr(FNlog)-pi*pi);
    ac.second = 0.5*pi*FNlog;
  } else {
    ac.first = sqr(0.5*pi);
  }
  return ac;
}

// Taken from HERWIG 6510 with simplifications.
Energy MEPP2Higgs::calcLORunningWidth(Energy Mh) const {
  Energy2 q2 = sqr(Mh);
  Energy QCDLambda = theSM->LambdaQCD(q2);
//  double alphaEM   = theSM->alphaEM(q2);
  double alphaEM   = theSM->alphaEM();
  double alphaS    = theSM->alphaS(q2);
  double sw2       = theSM->sin2ThetaW();
  Energy partialw[13] = {Energy()};

  tPDPtr w = getParticleData(ParticleID::Wplus);
  tPDPtr z = getParticleData(ParticleID::Z0);
  Energy wmnom = w->mass();
  Energy zmnom = z->mass();

  double ncolour = theSM->Nc();
  double Ca = ncolour;
  const unsigned int minflavour = 1;
  const unsigned int maxflavour = 6;
  double pi = Constants::pi;

// All calculation are being done for Monte-Carlo QCD Lambda, except Higgs width...
  double bcoeff4=(11.*Ca-10.)/(12.*pi);
  double kfac(Ca*(67./18.-sqr(pi)/6.)-25./9.);
  QCDLambda /= exp(kfac/(4.*pi*bcoeff4))/sqrt(2.);

// H->fermion pair
  Energy qmass[7];
  double nflavour = minflavour-1;
  for (unsigned int i = minflavour; i <= maxflavour; ++i) {
    tcPDPtr quark = getParticleData(i);
    qmass[i] = quark->mass();
    if (2.0*qmass[i] < Mh) nflavour+=1.;
  }
  Energy2 GFermiINV = 8.*sw2*sqr(wmnom)/alphaEM;

// quarks: partialw[1-6]
  for (unsigned int i = minflavour; i <= maxflavour; ++i) {
    Energy mf = qmass[i];
    double xf = sqr(mf/Mh);
    if (xf < 0.25) partialw[i] = Ca*Mh*sqr(mf)*pow(1.0-4.0*xf,1.5)/GFermiINV;
  }

// leptons: partialw[7-9]
  for (unsigned int i = 0; i < 3; ++i) {
    tcPDPtr lepton = getParticleData(11+2*i);
    Energy mf = lepton->mass();
    double xf = sqr(mf/Mh);
    if (xf < 0.25) partialw[7+i] = Mh*sqr(mf)*pow(1.0-4.0*xf,1.5)/GFermiINV;
  }

// H->W*W*/Z*Z*: partialw[10,11]
  if (branchingopt > 1) {
    double xfw = sqr(wmnom/Mh);
    double xfz = sqr(zmnom/Mh);
    if (2.0*wmnom < Mh) partialw[10] = Mh*Mh*Mh*sqrt(1-4.0*xfw)*(1-xfw+0.75*sqr(xfw))/2./GFermiINV;
    if (2.0*zmnom < Mh) partialw[11] = Mh*Mh*Mh*sqrt(1-4.0*xfz)*(1-xfz+0.75*sqr(xfz))/4./GFermiINV;
  }

// H->gamma,gamma: partialw[12]
  if (branchingopt > 2) {
    double taut = sqr(2.0*qmass[ParticleID::t]/Mh);
    double tauw = sqr(2.0*wmnom/Mh);
    std::pair<double,double> ftaut = HwW2 (taut);
    std::pair<double,double> ftauw = HwW2 (tauw);
    double re = 4.0/3.0*(-2.0*taut*(1.0+(1.0-taut)*ftaut.first))+(2.0+3.0*tauw*(1+(2.0-tauw)*ftauw.first));
    double im = 4.0/3.0*(-2.0*taut*(    (1.0-taut)*ftaut.second))+(   3.0*tauw*(  (2.0-tauw)*ftauw.second));
    partialw[12] = sqr(alphaEM/pi)*Mh*Mh*Mh*(sqr(re)+ sqr(im))/32./GFermiINV;
  }

// H->gluon,gluon: partialw[13]
  if (branchingopt > 3) {
    double tau = sqr(2.0*qmass[ParticleID::t]/Mh);
    std::pair<double,double> ftau = HwW2 (tau);
    double re = 1+(1.0-tau)*ftau.first;
    double im =   (1.0-tau)*ftau.second;
    partialw[13] = sqr(alphaS/pi)*Mh*Mh*Mh*sqr(tau)*(sqr(re)+ sqr(im))/4./GFermiINV;
  }

  Energy higgswidth = Energy();
  for (unsigned int i = 1; i < 13; ++i) {
    higgswidth += partialw[i];
  }
  return higgswidth;
}
