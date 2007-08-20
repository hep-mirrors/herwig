// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsQQ class.
//

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "HardVertex.h"
#include "MEPP2HiggsQQ.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MEPP2HiggsQQ.tcc"
#endif

using namespace Herwig;

MEPP2HiggsQQ::~MEPP2HiggsQQ() {}

unsigned int MEPP2HiggsQQ::orderInAlphaS() const {
  return 3;
}

unsigned int MEPP2HiggsQQ::orderInAlphaEW() const {
  return 1;
}

void MEPP2HiggsQQ::doinit() throw(InitException) {
  // get the vertex pointers from the SM object
  theSM = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!theSM) {
    throw InitException() << "Wrong type of StandardModel object in MEPP2HiggsQQ::doinit(),"
                          << " the Herwig++ version must be used" 
                          << Exception::runerror;
  }
  hggvertex = dynamic_ptr_cast<SVVLoopVertexPtr>(theSM->vertexHGG());
  ffhvertex = theSM->vertexFFH();
  // call the base class
  ME2to2Base::doinit();
}

void MEPP2HiggsQQ::persistentOutput(PersistentOStream & os) const {
  os << hggvertex << ffhvertex 
     << shapeopt << processopt << widthopt << branchingopt 
     << minflavouropt << maxflavouropt << quarkopt 
     << theSM;
}

void MEPP2HiggsQQ::persistentInput(PersistentIStream & is, int) {
  is >> hggvertex >> ffhvertex 
     >> shapeopt >> processopt >> widthopt >> branchingopt 
     >> minflavouropt >> maxflavouropt >> quarkopt 
     >> theSM;
}

ClassDescription<MEPP2HiggsQQ> MEPP2HiggsQQ::initMEPP2HiggsQQ;
// Definition of the static class description member.

void MEPP2HiggsQQ::Init() {

  static ClassDocumentation<MEPP2HiggsQQ> documentation
    ("The MEPP2HiggsQQ class implements the matrix elements for"
     " Higgs production (with decay H->QQ) in hadron-hadron collisions.");

//////////////////////////////////////////////////////////////////////////////////////////
  static Switch<MEPP2HiggsQQ,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the masses in the loop diagrams",
     &MEPP2HiggsQQ::shapeopt, 1, false, false);
  static SwitchOption interfaceStandardResonanse
    (interfaceShapeOption,
     "BreitWigner",
     "Breit-Wigner Higgs s-channel resonanse",
     1);
  static SwitchOption interfaceImprovedResonanse
    (interfaceShapeOption,
     "ImprovedBreitWigner",
     "Improved Higgs s-channel resonanse",
     6);

//////////////////////////////////////////////////////////////////////////////////////////
  static Switch<MEPP2HiggsQQ,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2HiggsQQ::processopt, 1, false, false);
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
  static Switch<MEPP2HiggsQQ,unsigned int> interfaceQuark
    ("DecayMode",
     "Which decay mode to include",
     &MEPP2HiggsQQ::quarkopt, 1, false, false);
  static SwitchOption interfaceBosonTT
    (interfaceQuark,
     "tt",
     "H->ttbar decay channel only",
     1);
  static SwitchOption interfaceBosonBB
    (interfaceQuark,
     "bb",
     "H->bbbar decay channel only",
     2);
  static SwitchOption interfaceBosonQQ
    (interfaceQuark,
     "All",
     "Both H->tt/bb decay channels",
     3);

//////////////////////////////////////////////////////////////////////////////////////////
  static Switch<MEPP2HiggsQQ,unsigned int> interfaceBranchingOption
    ("BranchAllowed",
     "Option to switch on/off branchings in the total Higgs width calculation",
     &MEPP2HiggsQQ::branchingopt, 3, false, false);
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
  static Switch<MEPP2HiggsQQ,unsigned int> interfaceWidthOption
    ("WidthScheme",
     "Option for the treatment of the Higss Width calculation",
     &MEPP2HiggsQQ::widthopt, 1, false, false);
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
//  static SwitchOption interfaceNLLWidth
//    (interfaceWidthOption,
//     "NLLcorrectedHiggsWidth",
//     "NLL corrected Higgs width (a-la FORTRAN HERWIG)",
//     3);
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
  static Parameter<MEPP2HiggsQQ,unsigned int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The minimum flavour of the incoming quarks in the hard process",
     &MEPP2HiggsQQ::minflavouropt, 3, 3, 5,
     false, false, Interface::limited);
  static Parameter<MEPP2HiggsQQ,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the incoming quarks in the hard process",
     &MEPP2HiggsQQ::maxflavouropt, 5, 3, 5,
     false, false, Interface::limited);

}

bool MEPP2HiggsQQ::generateKinematics(const double * r) {
  bool cutres = ME2to2Base::generateKinematics(r);
  if (cutres) {
    PDPtr h0 = getParticleData(ParticleID::h0);
    Energy2 s(sHat());
    if (s < 0.*GeV2 && widthopt != 1) {
      throw Exception() << "Warning! Shat < 0.0 (" << s/GeV2 
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
//        Energy runwidth = calcNLLRunningWidth(sqrt(s));
//        h0->width(runwidth);
        break;
      }
      case 4: {
        Energy runwidth = calcLORunningWidth(sqrt(s));
        h0->width(runwidth);
        break;
      }
      default: 
        throw Exception() << "Unknown width option in MEPP2Higgs::generateKinematics" 
                          << Exception::runerror;
        break;
      }
    }
//  cout << "Higgs width = " << ounit(h0->width(),GeV) << endl;
  }
  return cutres;
}

void MEPP2HiggsQQ::getDiagrams() const {
  tcPDPtr qout,qbout;
  tcPDPtr g  = getParticleData(ParticleID::g);
  tcPDPtr h0 = getParticleData(ParticleID::h0);
  switch (quarkopt) {
    case 1: {
      qout = getParticleData(ParticleID::t);
      qbout = qout->CC();
      break;
    }
    case 2: {
      qout = getParticleData(ParticleID::b);
      qbout = qout->CC();
      break;
    }
    case 3:
    case 4:
    default:
      std::cerr << "Invalid decay channel option! tt/bb/All are available only!" << std::endl;
      break;
   }

  // gg->H->qq
  if (1 == processopt || 3 == processopt) {
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, h0, 3, qout, 3, qbout, -1)));
  }
  // qq->H->qq
  if (1 == processopt || 2 == processopt) {
    for (unsigned int i = minflavouropt; i <= maxflavouropt; ++i ) {
      tcPDPtr qin = getParticleData(i);
      tcPDPtr qbin = qin->CC();
      add(new_ptr((Tree2toNDiagram(2), qin, qbin, 1, h0, 3, qout, 3, qbout, -2)));
    }
  }
}

Energy2 MEPP2HiggsQQ::scale() const {
  return sHat();
}

double MEPP2HiggsQQ::me2() const {
  double me(0.0);
  useMe();

  // gg->H->qq
  if (mePartonData()[0]->id()==ParticleID::g && mePartonData()[1]->id()==ParticleID::g) {
    VectorWaveFunction     gin1(meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction     gin2(meMomenta()[1],mePartonData()[1],incoming);
    SpinorWaveFunction     qout(meMomenta()[2],mePartonData()[2],outgoing);
    SpinorBarWaveFunction qbout(meMomenta()[3],mePartonData()[3],outgoing);

    vector<VectorWaveFunction> g1,g2;
    vector<SpinorWaveFunction> fout;
    vector<SpinorBarWaveFunction> aout;
    for(unsigned int ix = 0; ix < 2; ++ix) {
      gin1.reset(2*ix);
      g1.push_back(gin1);
      gin2.reset(2*ix);
      g2.push_back(gin2);
      qout.reset(ix); 
      fout.push_back(qout);
      qbout.reset(ix); 
      aout.push_back(qbout);
    }
    me = ggME (g1,g2,fout,aout,false);
  }
  // qq->H->qq
  else if(mePartonData()[0]->id()==-mePartonData()[1]->id()) {
    SpinorWaveFunction      qin(meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction  qbin(meMomenta()[1],mePartonData()[1],incoming);
    SpinorWaveFunction     qout(meMomenta()[2],mePartonData()[2],outgoing);
    SpinorBarWaveFunction qbout(meMomenta()[3],mePartonData()[3],outgoing);

    vector<SpinorWaveFunction> fin,fout;
    vector<SpinorBarWaveFunction> ain,aout;
    for(unsigned int ix = 0; ix < 2; ++ix) {
      qin.reset(ix);
      fin.push_back(qin);
      qbin.reset(ix);
      ain.push_back(qbin);
      qout.reset(ix); 
      fout.push_back(qout);
      qbout.reset(ix); 
      aout.push_back(qbout);
    }
    me = qqME(fin,ain,fout,aout,false); 
  } else {
    throw Exception() << "Unknown subprocess in MEPP2HiggsQQ::me2()" 
                      << Exception::runerror;
  }
  return me;
}

double MEPP2HiggsQQ::ggME(vector<VectorWaveFunction> & g1,
                          vector<VectorWaveFunction> & g2,
                          vector<SpinorWaveFunction> & q1,
                          vector<SpinorBarWaveFunction> & q2,
                          bool calc) const {
  tcPDPtr h0 = getParticleData(ParticleID::h0);
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,PDT::Spin1,PDT::Spin1);

  double me2(0.0);
  Energy2 s(sHat());
  for(int in1 = 0; in1 < 2; ++in1) {
    for(int in2 = 0; in2 < 2; ++in2) {
      for(int out1 = 0; out1 < 2; ++out1) {
        for(int out2 = 0; out2 < 2; ++out2) {
          ScalarWaveFunction higgsWF = hggvertex->evaluate(s,shapeopt,h0,g1[in1],g2[in2]);
          Complex me = ffhvertex->evaluate(s,q1[out1],q2[out2],higgsWF);
          me2 += real(me*conj(me));
          if(calc) newme(2*in1,2*in2,out1,out2) = me;
        }
      }
    }
  }
  if(calc) _me.reset(newme);
  // initial colour and spin factors: colour -> (8/64) and spin -> (1/4)
  // final colour : colour -> (3)
  return 3.*me2/(32.);
}

double MEPP2HiggsQQ::qqME(vector<SpinorWaveFunction>    & fin,
                          vector<SpinorBarWaveFunction> & ain,
                          vector<SpinorWaveFunction>    & fout,
                          vector<SpinorBarWaveFunction> & aout,
                          bool calc) const {
  tcPDPtr h0 = getParticleData(ParticleID::h0);
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1,PDT::Spin1);

  double me2(0.);
  Energy2 s(sHat());
  for(int in1 = 0; in1 < 2; ++in1) {
    for(int in2 = 0; in2 < 2; ++in2) {
      for(int out1 = 0; out1 < 2; ++out1) {
        for(int out2 = 0; out2 < 2; ++out2) {
          ScalarWaveFunction higgsWF = ffhvertex->evaluate(s,shapeopt,h0,fin[in1],ain[in2]);
          Complex me = ffhvertex->evaluate(s,fout[out1],aout[out2],higgsWF);
          me2+=real(me*conj(me));
          if(calc) newme(in1,in2,out1,out2) = me;
        }
      }
    }
  }
  if(calc) _me.reset(newme);
  // final colour/spin factors
  // final colour and spin factors: 
  // colour (initial) -> (3/9) 
  // colour (final) -> (3) 
  // spin -> (1/4)
  return me2*3.*3./9./4.;
}

Selector<MEBase::DiagramIndex> 
MEPP2HiggsQQ::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for (DiagramIndex i = 0; i < diags.size(); ++i) {
    if(abs(diags[i]->id())<4) sel.insert(1.0, i);
    else sel.insert(diagwgt[abs(diags[i]->id())-4], i);
  }
  return sel;
}

Selector<const ColourLines *>
MEPP2HiggsQQ::colourGeometries(tcDiagPtr diag) const {
  // colour lines for q,qbar -> H -> Z,Z
  static const ColourLines gluon("1 -2, 2 -1,4 -5");
  static const ColourLines iquark("1 -2,4 -5");
  // select the colour flow
  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 ) {
    sel.insert(1.0, &gluon);
  } else {
    sel.insert(1.0, &iquark);
  }
  return sel;
}

void MEPP2HiggsQQ::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);

  if(hard[0]->id()<hard[1]->id()) swap(hard[0],hard[1]);
  // identify the process and calculate the matrix element
  if(hard[0]->id()==ParticleID::g && hard[1]->id()==ParticleID::g) {
    vector<VectorWaveFunction> g1,g2;
    vector<SpinorWaveFunction> q1;
    vector<SpinorBarWaveFunction> q2;

    VectorWaveFunction    (g1,hard[0],incoming,false,true,true);
    VectorWaveFunction    (g2,hard[1],incoming,false,true,true);
    SpinorWaveFunction    (q1,hard[2],outgoing,true,true);
    SpinorBarWaveFunction (q2,hard[3],outgoing,true,true);
    g1[1] = g1[2];
    g2[1] = g2[2];
    ggME(g1,g2,q1,q2,true);
  } else {
    vector<SpinorWaveFunction>    q1,q3;
    vector<SpinorBarWaveFunction> q2,q4;
    unsigned int order[2]={2,3};
//    if(hard[0]->id() < 0) {order[0]=3; order[1]=2;}

    SpinorWaveFunction    (q1,hard[      0 ],incoming,false,true);
    SpinorBarWaveFunction (q2,hard[      1 ],incoming,false,true);
    SpinorWaveFunction    (q3,hard[order[0]],outgoing,true,true);
    SpinorBarWaveFunction (q4,hard[order[1]],outgoing,true,true);
    qqME(q1,q2,q3,q4,true);
  }

  // get the spin info objects
  SpinfoPtr spin[4];
  for (unsigned int i = 0; i < 4; ++i) {
    spin[i] = dynamic_ptr_cast<SpinfoPtr>(hard[i]->spinInfo());
  }
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int i = 0; i < 4; ++i) {
    spin[i]->setProductionVertex(hardvertex);
  }
}

// Taken from HERWIG 6510.
std::pair<double,double> MEPP2HiggsQQ::HwW2(double tau) const {
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
Energy MEPP2HiggsQQ::calcLORunningWidth(Energy Mh) const {
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
