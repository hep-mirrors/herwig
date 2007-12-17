// -*- C++ -*-
//
// MEPP2Higgs.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
#include "Herwig++/MatrixElement/General/HardVertex.h"

using namespace Herwig;

ClassDescription<MEPP2Higgs> MEPP2Higgs::initMEPP2Higgs;
// Definition of the static class description member.

void MEPP2Higgs::persistentOutput(PersistentOStream & os) const {
  os << hggvertex << ffhvertex << theSM << shapeopt << processopt 
     << minflavouropt << maxflavouropt << _hmass << ounit(_mh,GeV) << ounit(_wh,GeV);
}

void MEPP2Higgs::persistentInput(PersistentIStream & is, int) {
  is >> hggvertex >> ffhvertex >> theSM >> shapeopt >> processopt 
     >> minflavouropt >> maxflavouropt >> _hmass >> iunit(_mh,GeV) >> iunit(_wh,GeV);
}

void MEPP2Higgs::Init() {

  static ClassDocumentation<MEPP2Higgs> documentation
    ("The MEPP2Higgs class implements the matrix elements for"
     " Higgs production (with decay H->W-W+) in hadron-hadron collisions.");

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
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);

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

  static Parameter<MEPP2Higgs,unsigned int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The minimum flavour of the incoming quarks in the hard process",
     &MEPP2Higgs::minflavouropt, 4, 3, 5,
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
  ffhvertex = theSM->vertexFFH();
  // get the mass generator for the higgs
  PDPtr h0 = getParticleData(ParticleID::h0);
  _mh = h0->mass();
  _wh = h0->generateWidth(_mh);
  if(h0->massGenerator()) {
    _hmass=dynamic_ptr_cast<SMHiggsMassGeneratorPtr>(h0->massGenerator());
  }
  if(shapeopt==2&&!_hmass) throw InitException()
    << "If using the mass generator for the line shape in MEPP2Higgs::doinit()"
    << "the mass generator must be an instance of the SMHiggsMassGenerator class"
    << Exception::runerror;
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
  using Constants::pi;
  InvEnergy2 bwfact;
  if(shapeopt==1) {
    bwfact = mePartonData()[2]->generateWidth(sqrt(sHat()))*sqrt(sHat())/pi/
      (sqr(sHat()-sqr(_mh))+sqr(_mh*_wh));
  }
  else {
    bwfact = _hmass->BreitWignerWeight(sqrt(sHat()),0);
  }
  double cs = me2() * jacobian() * pi * double(UnitRemoval::E4 * bwfact/sHat());
  return UnitRemoval::InvE2 * sqr(hbarc) * cs;
}

double MEPP2Higgs::me2() const {
  double output(0.0);
  useMe();
  ScalarWaveFunction hout(meMomenta()[2],mePartonData()[2],outgoing);

// Safety code to garantee the reliable behaviour of Higgs shape limits 
// (important for heavy and broad Higgs resonance).
  Energy hmass = meMomenta()[2].m();
  tcPDPtr h0 = mePartonData()[2];
  Energy mass = h0->mass();
  Energy halfmass = .5*mass;
  if (.0*GeV > hmass) return 0.0;
  // stricly speaking the condition is applicable if h0->widthUpCut() == h0->widthLoCut()...
  if (h0->widthLoCut() > halfmass) {
    if ((mass + h0->widthUpCut() < hmass || mass - h0->widthLoCut() > hmass)) return 0.0;
  } else {
    if (mass + halfmass < hmass || halfmass > hmass) return 0.0;
  }

  if (mePartonData()[0]->id() == ParticleID::g && 
      mePartonData()[1]->id() == ParticleID::g) {
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
    throw Exception() << "Unknown subprocess in MEPP2Higgs::me2()" 
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
                          ScalarWaveFunction & in, 
                          bool calc) const {
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  Energy2 s(sHat());
  double me2(0.0);
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 2; ++j) {
      Complex diag = hggvertex->evaluate(s,in,g1[i],g2[j]);
      me2 += norm(diag);
      if(calc) newme(2*i, 2*j, 0) = diag;
    }
  }
  if(calc) _me.reset(newme);
  // initial colour and spin factors: colour -> (8/64) and spin -> (1/4)
  return me2/32.;
}


double MEPP2Higgs::qqME(vector<SpinorWaveFunction> & fin, 
                          vector<SpinorBarWaveFunction> & ain, 
                          ScalarWaveFunction & in, 
                          bool calc) const {
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0);
  Energy2 s(scale());
  double me2(0.0);
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 2; ++j) {
      Complex diag = ffhvertex->evaluate(s,fin[i],ain[j],in);
      me2+=norm(diag);
      if(calc) newme(i, j, 0) = diag;
    }
  }
  if(calc) _me.reset(newme);
  // final colour/spin factors
  return me2/12.;
}
