// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2HiggsVV class.
//
#include "MEPP2HiggsVV.h"
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
#include "GeneralHardME.h"
#include "HardVertex.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MEPP2HiggsVV.tcc"
#endif

using namespace Herwig;

MEPP2HiggsVV::~MEPP2HiggsVV() {}

unsigned int MEPP2HiggsVV::orderInAlphaS() const {
  return 3;
}

unsigned int MEPP2HiggsVV::orderInAlphaEW() const {
  return 1;
}

void MEPP2HiggsVV::doinit() throw(InitException) {
  // get the vertex pointers from the SM object
  theSM = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!theSM) {
    throw InitException() << "Wrong type of StandardModel object in MEPP2HiggsQQ::doinit(),"
                          << " the Herwig++ version must be used" 
                          << Exception::runerror;
  }
  vvhvertex = theSM->vertexWWH();
  hggvertex = dynamic_ptr_cast<SVVLoopVertexPtr>(theSM->vertexHGG());
  hppvertex = dynamic_ptr_cast<SVVLoopVertexPtr>(theSM->vertexHPP());
  ffhvertex = theSM->vertexFFH();

  // call the base class
  ME2to2Base::doinit();
}

void MEPP2HiggsVV::persistentOutput(PersistentOStream & os) const {
  os << hggvertex << vvhvertex << ffhvertex << hppvertex
     << shapeopt << processopt 
     << minflavouropt << maxflavouropt << bosonopt
     << theSM;
}

void MEPP2HiggsVV::persistentInput(PersistentIStream & is, int) {
  is >> hggvertex >> vvhvertex >> ffhvertex >> hppvertex
     >> shapeopt >> processopt 
     >> minflavouropt >> maxflavouropt >> bosonopt
     >> theSM;
}

ClassDescription<MEPP2HiggsVV> MEPP2HiggsVV::initMEPP2HiggsVV;
// Definition of the static class description member.

void MEPP2HiggsVV::Init() {

  static ClassDocumentation<MEPP2HiggsVV> documentation
    ("The MEPP2HiggsVV class implements the matrix elements for"
     " Higgs production (with decay H->WW/ZZ) in hadron-hadron collisions.");

//////////////////////////////////////////////////////////////////////////////////////////
  static Switch<MEPP2HiggsVV,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the masses in the loop diagrams",
     &MEPP2HiggsVV::shapeopt, 1, false, false);
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
  static Switch<MEPP2HiggsVV,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2HiggsVV::processopt, 1, false, false);
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
  static Switch<MEPP2HiggsVV,unsigned int> interfaceBoson
    ("DecayMode",
     "Which decay mode to include",
     &MEPP2HiggsVV::bosonopt, 1, false, false);
  static SwitchOption interfaceBosonWW
    (interfaceBoson,
     "WW",
     "H->WW decay channel only",
     1);
  static SwitchOption interfaceBosonZZ
    (interfaceBoson,
     "ZZ",
     "H->ZZ decay channel only",
     2);
  static SwitchOption interfaceBosonVV
    (interfaceBoson,
     "GammaGamma",
     "H->gamma,gamma decay channels",
     3);

//////////////////////////////////////////////////////////////////////////////////////////
  static Parameter<MEPP2HiggsVV,unsigned int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The minimum flavour of the incoming quarks in the hard process",
     &MEPP2HiggsVV::minflavouropt, 3, 3, 5,
     false, false, Interface::limited);
  static Parameter<MEPP2HiggsVV,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the incoming quarks in the hard process",
     &MEPP2HiggsVV::maxflavouropt, 5, 3, 5,
     false, false, Interface::limited);

}

void MEPP2HiggsVV::getDiagrams() const {
  tcPDPtr g  = getParticleData(ParticleID::g);
  tcPDPtr h0 = getParticleData(ParticleID::h0);
  tcPDPtr v1,v2;
  switch (bosonopt) {
    case 1: {
      v1 = getParticleData(ParticleID::Wplus);
      v2 = getParticleData(ParticleID::Wminus);
      break;}
    case 2: {
      v1 = getParticleData(ParticleID::Z0);
      v2 = getParticleData(ParticleID::Z0);
      break;}
    case 3: {
      v1 = getParticleData(ParticleID::gamma);
      v2 = getParticleData(ParticleID::gamma);
      break;}
    default: {
      throw MEException() << "Invalid decay channel option! W/Z/Gamma modes are available only!" 
                          << Exception::runerror;
      break;}
   }

  // gg->H->VV
  if (1 == processopt || 3 == processopt) {
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, h0, 3, v1, 3, v2, -1)));
  }
  // qqbar->H->VV
  if (1 == processopt || 2 == processopt) {
    for (unsigned int i = minflavouropt; i <= maxflavouropt; ++i ) {
      tcPDPtr q = getParticleData(i);
      tcPDPtr qb = q->CC();
    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, h0, 3, v1, 3, v2, -2)));
    }
  }
}

bool MEPP2HiggsVV::generateKinematics(const double * r) {
  return ME2to2Base::generateKinematics(r);
}

Energy2 MEPP2HiggsVV::scale() const {
  return sHat();
}

double MEPP2HiggsVV::me2() const {
  double me(0.0);
  useMe();

  // gg->H->VV
  if (mePartonData()[0]->id()==ParticleID::g && mePartonData()[1]->id()==ParticleID::g) {
    VectorWaveFunction gin1(meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction gin2(meMomenta()[1],mePartonData()[1],incoming);
    VectorWaveFunction vout1(meMomenta()[2],mePartonData()[2],outgoing);
    VectorWaveFunction vout2(meMomenta()[3],mePartonData()[3],outgoing);

    vector<VectorWaveFunction> g1,g2;
    for(unsigned int i = 0; i < 2; ++i) {
      gin1.reset(2*i);
      g1.push_back(gin1);
      gin2.reset(2*i);
      g2.push_back(gin2);
    }
    vector<VectorWaveFunction> v1,v2;
    switch (bosonopt) {
      case 1:
      case 2: {
        for(unsigned int i = 0; i < 3; ++i) {
          vout1.reset(i); 
          v1.push_back(vout1);
          vout2.reset(i); 
          v2.push_back(vout2);
        }
        break;}
      case 3: {
        for(unsigned int i = 0; i < 2; ++i) {
          vout1.reset(2*i);
          v1.push_back(vout1);
          vout2.reset(2*i);
          v2.push_back(vout2);
        }
        break;}
      default: {
        throw MEException() << "Invalid decay channel option! W/Z/Gamma modes are available only!" 
                            << Exception::runerror;
        break;}
    }
    me = ggME (g1,g2,v1,v2,false);
  }
  // qqbar->H->VV
  else if(mePartonData()[0]->id()==-mePartonData()[1]->id()) {
    SpinorWaveFunction     qin(meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction qbin(meMomenta()[1],mePartonData()[1],incoming);
    VectorWaveFunction   vout1(meMomenta()[2],mePartonData()[2],outgoing);
    VectorWaveFunction   vout2(meMomenta()[3],mePartonData()[3],outgoing);

    vector<SpinorWaveFunction> fin;
    vector<SpinorBarWaveFunction> ain;
    for(unsigned int i = 0; i < 2; ++i) {
      qin.reset(i);
      fin.push_back(qin);
      qbin.reset(i);
      ain.push_back(qbin);
    }
    vector<VectorWaveFunction> v1,v2;
    switch (bosonopt) {
      case 1:
      case 2: {
        for(unsigned int i = 0; i < 3; ++i) {
          vout1.reset(i); 
          v1.push_back(vout1);
          vout2.reset(i); 
          v2.push_back(vout2);
        }
        break;}
      case 3: {
        for(unsigned int i = 0; i < 2; ++i) {
          vout1.reset(2*i);
          v1.push_back(vout1);
          vout2.reset(2*i);
          v2.push_back(vout2);
        }
        break;}
    }
    // calculate the matrix element
    me = qqME(fin,ain,v1,v2,false); 
  } else {
    throw MEException() << "Unknown subprocess in MEPP2HiggsVV::me2()" 
                      << Exception::runerror;
  }
  return me;
}

double MEPP2HiggsVV::ggME(vector<VectorWaveFunction> g1,
                          vector<VectorWaveFunction> g2,
                          vector<VectorWaveFunction> v1,
                          vector<VectorWaveFunction> v2, bool calc) const {
  tcPDPtr h0 = getParticleData(ParticleID::h0);
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,PDT::Spin1,PDT::Spin1);

  double me2(0.0);
  Energy2 s(sHat());
  for(int in1 = 0; in1 < 2; ++in1) {
    for(int in2 = 0; in2 < 2; ++in2) {
      for(int out1 = 0; out1 < 3; ++out1) {
        for(int out2 = 0; out2 < 3; ++out2) {
          ScalarWaveFunction higgsWF = hggvertex->evaluate(s,shapeopt,h0,g1[in1],g2[in2]);
          Complex me(0.0);
          switch (bosonopt) {
            case 1:
            case 2: {me = vvhvertex->evaluate(s,v1[out1],v2[out2],higgsWF); break;}
            case 3: {me = hppvertex->evaluate(s,higgsWF,v1[out1],v2[out2]); break;}
          }
          me2 += real(me*conj(me));
          if(calc) newme(2*in1,2*in2,out1,out2) = me;
        }
      }
    }
  }
  if(calc) _me.reset(newme);
  // final colour and spin factors: 
  // colour -> (8/64) 
  // spin -> (1/4)
  double colourspinfactor = 32.;
  if (2 == bosonopt) colourspinfactor *= 2.;  // extra factor for Z

  return me2/colourspinfactor;
}

double MEPP2HiggsVV::qqME(vector<SpinorWaveFunction>    fin,
                          vector<SpinorBarWaveFunction> ain,
                          vector<VectorWaveFunction> v1,
                          vector<VectorWaveFunction> v2, bool calc) const {
  tcPDPtr h0 = getParticleData(ParticleID::h0);
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1,PDT::Spin1);

  double me2(0.);
  Energy2 s(sHat());
  for(int in1 = 0; in1 < 2; ++in1) {
    for(int in2 = 0; in2 < 2; ++in2) {
      for(int out1 = 0; out1 < 3; ++out1) {
        for(int out2 = 0; out2 < 3; ++out2) {
          ScalarWaveFunction higgsWF = ffhvertex->evaluate(s,shapeopt,h0,fin[in1],ain[in2]);
          Complex me(0.0);
          switch (bosonopt) {
            case 1:
            case 2: {me = vvhvertex->evaluate(s,v1[out1],v2[out2],higgsWF); break;}
            case 3: {me = hppvertex->evaluate(s,higgsWF,v1[out1],v2[out2]); break;}
          }
          me2+=real(me*conj(me));
          if(calc) newme(in1,in2,out1,out2) = me;
        }
      }
    }
  }
  if(calc) _me.reset(newme);
  // final colour and spin factors: 
  // colour -> (3/9) 
  // spin -> (1/4)
  double colourspinfactor = 3.*4.;
  if (2 == bosonopt) colourspinfactor *= 2.;  // extra factor for Z

  return me2/colourspinfactor;
}

Selector<MEBase::DiagramIndex> 
MEPP2HiggsVV::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for (DiagramIndex i = 0; i < diags.size(); ++i) {
    if(abs(diags[i]->id())<4) sel.insert(1.0, i);
    else sel.insert(diagwgt[abs(diags[i]->id())-4], i);
  }
  return sel;
}

Selector<const ColourLines *>
MEPP2HiggsVV::colourGeometries(tcDiagPtr diag) const {
  // colour lines for q,qbar -> H -> W-,W+
  static const ColourLines gluon("1 -2, 2 -1");
  static const ColourLines quark("1 -2");
  // select the colour flow
  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 ) {
    sel.insert(1.0, &gluon);
  } else {
    sel.insert(1.0, &quark);
  }
  return sel;
}

void MEPP2HiggsVV::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);

  // identify the process and calculate the matrix element
  if(hard[0]->id()<hard[1]->id()) swap(hard[0],hard[1]);
  if(hard[0]->id()==ParticleID::g && hard[1]->id()==ParticleID::g) {
    vector<VectorWaveFunction> g1,g2;
    vector<VectorWaveFunction> v1,v2;
    VectorWaveFunction (g1,hard[0],incoming,false,true,true);
    VectorWaveFunction (g2,hard[1],incoming,false,true,true);
    g1[1] = g1[2];
    g2[1] = g2[2];
    VectorWaveFunction (v1,hard[2],outgoing,true,true,true);
    VectorWaveFunction (v2,hard[3],outgoing,true,true,true);
    ggME(g1,g2,v1,v2,true);
  } else {
    vector<SpinorWaveFunction>    q1;
    vector<SpinorBarWaveFunction> q2;
    vector<VectorWaveFunction> v1,v2;
    SpinorWaveFunction    (q1,hard[0],incoming,false,true);
    SpinorBarWaveFunction (q2,hard[1],incoming,false,true);
    VectorWaveFunction (v1,hard[2],outgoing,true,true,true);
    VectorWaveFunction (v2,hard[3],outgoing,true,true,true);
    qqME(q1,q2,v1,v2,true);
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
