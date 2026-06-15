// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MESimpleGammaJet class.
//

#include "MESimpleGammaJet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig/MatrixElement/HardVertex.h"

using namespace Herwig;

static const double fact = 1.0;

MESimpleGammaJet::MESimpleGammaJet()
  : _process(0), _maxflavour(5) {
  massOption(std::vector<unsigned int>(2,1));
}

void MESimpleGammaJet::doinit() {
  HwMEBase::doinit();
  _gamma = getParticleData(ThePEG::ParticleID::gamma);
  ThePEG::Ptr<Herwig::StandardModel>::transient_const_pointer
    hwsm = ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardModel>
    ::transient_const_pointer>(standardModel());
  if(!hwsm)
    throw InitException() << "Must be Herwig::StandardModel in MESimpleGammaJet::doinit()"
                          << Exception::runerror;
  _theFFGammaVertex = hwsm->vertexFFP();
  _theQQGVertex     = hwsm->vertexFFG();
}

void MESimpleGammaJet::getDiagrams() const {
  typedef std::vector<std::pair<long,long> > Pairvector;
  Pairvector parentpair;
  parentpair.reserve(6);
  switch(_maxflavour) {
  case 5: parentpair.push_back(std::make_pair(ParticleID::b, ParticleID::bbar)); [[fallthrough]];
  case 4: parentpair.push_back(std::make_pair(ParticleID::c, ParticleID::cbar)); [[fallthrough]];
  case 3: parentpair.push_back(std::make_pair(ParticleID::s, ParticleID::sbar)); [[fallthrough]];
  case 2: parentpair.push_back(std::make_pair(ParticleID::d, ParticleID::dbar)); [[fallthrough]];
  case 1: parentpair.push_back(std::make_pair(ParticleID::u, ParticleID::ubar)); [[fallthrough]];
  default: ;
  }
  tcPDPtr g = getParticleData(ParticleID::g);
  for (std::pair<long,long> parent : parentpair) {
    tcPDPtr qNeg1 = getParticleData(parent.first);
    tcPDPtr qNeg2 = getParticleData(parent.second);
    tcPDPtr qPos1 = qNeg2->CC();
    tcPDPtr qPos2 = qNeg1->CC();
    if(_process==0||_process==1) {
      add(new_ptr((Tree2toNDiagram(3), qNeg1, qNeg2, qNeg2, 2, g, 1, _gamma, -1)));
      add(new_ptr((Tree2toNDiagram(3), qNeg1, qNeg1, qNeg2, 1, g, 2, _gamma, -2)));
    }
    if(_process==0||_process==2) {
      add(new_ptr((Tree2toNDiagram(3), qNeg1, qPos1, g, 2, qPos1, 1, _gamma, -5)));
      add(new_ptr((Tree2toNDiagram(2), qNeg1, g, 1, qNeg1, 3, qPos1, 3, _gamma, -6)));
    }
    if(_process==0||_process==3) {
      add(new_ptr((Tree2toNDiagram(3), qPos2, qNeg2, g, 2, qNeg2, 1, _gamma, -11)));
      add(new_ptr((Tree2toNDiagram(2), qPos2, g, 1, qPos2, 3, qNeg2, 3, _gamma, -12)));
    }
  }
}

Energy2 MESimpleGammaJet::scale() const {
  return sqr(91.188)*GeV2;
}

double MESimpleGammaJet::me2() const {
  double output(ZERO);

  // Photon: only transverse helicities 0 and 2
  std::vector<VectorWaveFunction> A;
  {
    VectorWaveFunction Aout(meMomenta()[3], mePartonData()[3], outgoing);
    for (unsigned int ix = 0; ix < 2; ++ix) {
      Aout.reset(2*ix); // 0,2
      A.push_back(Aout);
    }
  }

  // q g -> q γ
  if (mePartonData()[0]->id() <= 6 && mePartonData()[0]->id() > 0 &&
      mePartonData()[1]->id() == ParticleID::g) {

    std::vector<SpinorWaveFunction>     fin;
    std::vector<VectorWaveFunction>     gin;
    std::vector<SpinorBarWaveFunction>  fout;

    SpinorWaveFunction     qin (meMomenta()[0], mePartonData()[0], incoming);
    VectorWaveFunction     glin(meMomenta()[1], mePartonData()[1], incoming);
    SpinorBarWaveFunction  qout(meMomenta()[2], mePartonData()[2], outgoing);

    for (unsigned int ix = 0; ix < 2; ++ix) {
      qin.reset(ix);         fin.push_back(qin);
      glin.reset(2*ix);      gin.push_back(glin);  // 0,2
      qout.reset(ix);        fout.push_back(qout);
    }
    output = qgME(fin, gin, fout, A);
  }
  // qbar g -> qbar γ
  else if (mePartonData()[0]->id() >= -6 && mePartonData()[0]->id() < 0 &&
           mePartonData()[1]->id() == ParticleID::g) {

    std::vector<SpinorBarWaveFunction>  ain;
    std::vector<VectorWaveFunction>     gin;
    std::vector<SpinorWaveFunction>     aout;

    SpinorBarWaveFunction qbin(meMomenta()[0], mePartonData()[0], incoming);
    VectorWaveFunction    glin(meMomenta()[1], mePartonData()[1], incoming);
    SpinorWaveFunction    qbout(meMomenta()[2], mePartonData()[2], outgoing);

    for (unsigned int ix = 0; ix < 2; ++ix) {
      qbin.reset(ix);        ain.push_back(qbin);
      glin.reset(2*ix);      gin.push_back(glin);  // 0,2
      qbout.reset(ix);       aout.push_back(qbout);
    }
    output = qbargME(ain, gin, aout, A);
  }
  // q qbar -> γ g
  else {
    std::vector<SpinorWaveFunction>     fin;
    std::vector<SpinorBarWaveFunction>  ain;
    std::vector<VectorWaveFunction>     gout;

    SpinorWaveFunction     qin (meMomenta()[0], mePartonData()[0], incoming);
    SpinorBarWaveFunction  qbin(meMomenta()[1], mePartonData()[1], incoming);
    VectorWaveFunction     glout(meMomenta()[2], mePartonData()[2], outgoing);

    for (unsigned int ix = 0; ix < 2; ++ix) {
      qin.reset(ix);         fin.push_back(qin);
      qbin.reset(ix);        ain.push_back(qbin);
      glout.reset(2*ix);     gout.push_back(glout); // 0,2
    }
    output = qqbarME(fin, ain, gout, A);
  }

  return output * sqr(fact);
}

unsigned int MESimpleGammaJet::orderInAlphaS() const { return 1; }
unsigned int MESimpleGammaJet::orderInAlphaEW() const { return 1; }

Selector<MEBase::DiagramIndex>
MESimpleGammaJet::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    int id = abs(diags[i]->id());
    if     (id <= 2 ) sel.insert(meInfo()[id- 1],i);
    else if(id <= 4 ) sel.insert(meInfo()[id- 3],i);
    else if(id <= 6 ) sel.insert(meInfo()[id- 5],i);
    else if(id <= 8 ) sel.insert(meInfo()[id- 7],i);
    else if(id <= 10) sel.insert(meInfo()[id- 9],i);
    else if(id <= 12) sel.insert(meInfo()[id-11],i);
  }
  return sel;
}

Selector<const ColourLines *>
MESimpleGammaJet::colourGeometries(tcDiagPtr diag) const {
  static const ColourLines cqqbar[2]={ColourLines("1 -2 4,-3 -4"),
                                      ColourLines("1 4, -4 2 -3")};
  static const ColourLines cqg   [2]={ColourLines("1 2 -3,3 4"),
                                      ColourLines("1 -2,2 3 4")};
  static const ColourLines cqbarg[2]={ColourLines("-1 -2 3,-3 -4"),
                                      ColourLines("-1 2,-2 -3 -4")};
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
  case  1 : case  3: sel.insert(1.0, &cqqbar[0]); break;
  case  2 : case  4: sel.insert(1.0, &cqqbar[1]); break;
  case  5 : case  7: sel.insert(1.0, &cqg[0]);    break;
  case  6 : case  8: sel.insert(1.0, &cqg[1]);    break;
  case  9 : case 11: sel.insert(1.0, &cqbarg[0]); break;
  case 10 : case 12: sel.insert(1.0, &cqbarg[1]); break;
  }
  return sel;
}

IBPtr MESimpleGammaJet::clone() const     { return new_ptr(*this); }
IBPtr MESimpleGammaJet::fullclone() const { return new_ptr(*this); }

void MESimpleGammaJet::persistentOutput(PersistentOStream & os) const {
  os << _theFFGammaVertex << _theQQGVertex << _gamma << _process << _maxflavour;
}

void MESimpleGammaJet::persistentInput(PersistentIStream & is, int) {
  is >> _theFFGammaVertex >> _theQQGVertex >> _gamma >> _process >> _maxflavour;
}

// The following static variable is needed for the type description system.
DescribeClass<MESimpleGammaJet,HwMEBase>
  describeHerwigMESimpleGammaJet("Herwig::MESimpleGammaJet", "HwMEHadron.so");

void MESimpleGammaJet::Init() {

  static ClassDocumentation<MESimpleGammaJet> documentation
    ("The MESimpleGammaJet class implements the matrix element for gamma+jet production");

  static Parameter<MESimpleGammaJet,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MESimpleGammaJet::_maxflavour, 5, 0, 8, false, false, true);

  static Switch<MESimpleGammaJet,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MESimpleGammaJet::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess, "All",   "Include all subprocesses", 0);
  static SwitchOption interfaceProcessqqbar
    (interfaceProcess, "qqbar", "Only include q qbar -> gamma g process", 1);
  static SwitchOption interfaceProcessqg
    (interfaceProcess, "qg",    "Only include the q g -> gamma q process", 2);
  static SwitchOption interfaceProcessqbarg
    (interfaceProcess, "qbarg", "Only include the qbar g -> gamma qbar process", 3);
}

double MESimpleGammaJet::qqbarME(std::vector<SpinorWaveFunction> & fin,
                                 std::vector<SpinorBarWaveFunction> & ain,
                                 std::vector<VectorWaveFunction> & gout,
                                 std::vector<VectorWaveFunction> & Aout,
                                 bool calc) const {
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
                                             PDT::Spin1,PDT::Spin1));
  unsigned int ihel1,ihel2,ohel1,ohel2;
  double me[3]={0.,0.,0.};
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
        inters=_theQQGVertex->evaluate(scale(),5,mePartonData()[0],
                                       fin[ihel1],gout[ohel1]);
        interb=_theQQGVertex->evaluate(scale(),5,mePartonData()[1],
                                       ain[ihel2],gout[ohel1]);
        for(ohel2=0;ohel2<2;++ohel2) {
          Complex d0 = _theFFGammaVertex->evaluate(scale(),fin[ihel1],interb,Aout[ohel2]);
          Complex d1 = _theFFGammaVertex->evaluate(scale(),inters,ain[ihel2],Aout[ohel2]);
          me[1] += norm(d0);
          me[2] += norm(d1);
          Complex dt = d0 + d1;
          me[0] += norm(dt);
          if(calc) _me(ihel1,ihel2,2*ohel1,2*ohel2) = dt;
        }
      }
    }
  }
  double colspin = 1./9./4.;
  colspin *= 4.;
  DVector save;
  for(unsigned int ix=0;ix<3;++ix) {
    me[ix]*=colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0];
}

double MESimpleGammaJet::qgME(std::vector<SpinorWaveFunction> & fin,
                              std::vector<VectorWaveFunction> & gin,
                              std::vector<SpinorBarWaveFunction> & fout,
                              std::vector<VectorWaveFunction> & Aout,
                              bool calc) const {
  if (calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half, PDT::Spin1,
                                              PDT::Spin1Half, PDT::Spin1));

  unsigned int ihel1, ihel2, ohel1, ohel2;
  double me[3] = {0., 0., 0.};
  Complex diag[2];

  SpinorWaveFunction     inters;
  SpinorBarWaveFunction  interb;

  for (ihel1 = 0; ihel1 < 2; ++ihel1) {
    for (ihel2 = 0; ihel2 < 2; ++ihel2) {
      for (ohel1 = 0; ohel1 < 2; ++ohel1) {

        // q g -> q* -> q gamma
        //
        // The intermediate object is on the incoming-quark line.
        // Use the incoming quark ParticleData, not the outgoing CC.
        inters = _theQQGVertex->evaluate(scale(), 5, mePartonData()[0],
                                         fin[ihel1], gin[ihel2]);

        // q g -> q gamma with the gluon attached to the outgoing-quark side.
        //
        // fout is a SpinorBarWaveFunction for the outgoing quark, so the FFG
        // vertex needs the charge-conjugate ParticleData here.
        //
        // Crucially: do NOT use mePartonData()[1], since that is the gluon.
        interb = _theQQGVertex->evaluate(scale(), 5, mePartonData()[2]->CC(),
                                         fout[ohel1], gin[ihel2]);

        for (ohel2 = 0; ohel2 < 2; ++ohel2) { // photon helicities 0,2

          diag[0] = _theFFGammaVertex->evaluate(scale(),
                                                inters, fout[ohel1], Aout[ohel2]);

          diag[1] = _theFFGammaVertex->evaluate(scale(),
                                                fin[ihel1], interb, Aout[ohel2]);

          me[1] += norm(diag[0]);
          me[2] += norm(diag[1]);

          diag[0] += diag[1];
          me[0] += norm(diag[0]);

          if (calc) _me(ihel1, 2*ihel2, ohel1, 2*ohel2) = diag[0];
        }
      }
    }
  }

  // Initial-state spin and colour average:
  // q g has colour average 1/(3*8), spin average 1/(2*2).
  double colspin = 1./24./4.;

  // C_F N_c colour factor
  colspin *= 4.;

  DVector save;
  for (unsigned int ix = 0; ix < 3; ++ix) {
    me[ix] *= colspin;
    if (ix > 0) save.push_back(me[ix]);
  }

  meInfo(save);
  return me[0];
}

double MESimpleGammaJet::qbargME(std::vector<SpinorBarWaveFunction> & fin,
                                 std::vector<VectorWaveFunction> & gin,
                                 std::vector<SpinorWaveFunction> & fout,
                                 std::vector<VectorWaveFunction> & Aout,
                                 bool calc) const {
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
                                             PDT::Spin1Half,PDT::Spin1));
  unsigned int ihel1,ihel2,ohel1,ohel2;
  double me[3]={0.,0.,0.};
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
        inters=_theQQGVertex->evaluate(scale(),5,mePartonData()[2]->CC(),
                                       fout[ohel1],gin[ihel2]);
        interb=_theQQGVertex->evaluate(scale(),5,mePartonData()[0],
                                       fin[ihel1],gin[ihel2]);
        for(ohel2=0;ohel2<2;++ohel2) { // photon: 0,2
          Complex d0 = _theFFGammaVertex->evaluate(scale(),inters,fin[ihel1],Aout[ohel2]);
          Complex d1 = _theFFGammaVertex->evaluate(scale(),fout[ohel1],interb,Aout[ohel2]);
          me[1] += norm(d0);
          me[2] += norm(d1);
          Complex dt = d0 + d1;
          me[0] += norm(dt);
          if(calc) _me(ihel1,2*ihel2,ohel1,2*ohel2) = dt;
        }
      }
    }
  }
  double colspin = 1./24./4.;
  colspin *= 4.;
  DVector save;
  for(unsigned int ix=0;ix<3;++ix) {
    me[ix]*=colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0];
}

void MESimpleGammaJet::constructVertex(tSubProPtr sub) {
  ParticleVector hard(4);

  // incoming: keep quark first
  hard[0] = sub->incoming().first;
  hard[1] = sub->incoming().second;
  if ((hard[0]->id() < 0 && hard[1]->id() <= 6) || hard[0]->id() == ParticleID::g)
    swap(hard[0], hard[1]);

  // outgoing: photon in hard[3], the other leg in hard[2]
  bool haveGamma = false;
  for (unsigned int ix = 0; ix < 2; ++ix) {
    if (sub->outgoing()[ix]->id() == ParticleID::gamma) {
      hard[3] = sub->outgoing()[ix];
      haveGamma = true;
    } else {
      hard[2] = sub->outgoing()[ix];
    }
  }
  if (!haveGamma || !hard[2]) return;

  // Photon transverse polarizations (0,2)
  std::vector<VectorWaveFunction> Aout;
  {
    VectorWaveFunction A(hard[3]->momentum(), hard[3]->dataPtr(), outgoing);
    for (unsigned int ix = 0; ix < 2; ++ix) {
      A.reset(2*ix);
      Aout.push_back(A);
    }
  }

  // Identify hard process and compute ME (use only 0/2 for the gluon as well)
  if (hard[0]->id() <= 6 && hard[0]->id() > 0 && hard[1]->id() == ParticleID::g) {
    // q g -> q γ
    std::vector<SpinorWaveFunction>     fin;
    std::vector<VectorWaveFunction>     gin;
    std::vector<SpinorBarWaveFunction>  fout;

    SpinorWaveFunction     qin (hard[0]->momentum(), hard[0]->dataPtr(), incoming);
    VectorWaveFunction     glin(hard[1]->momentum(), hard[1]->dataPtr(), incoming);
    SpinorBarWaveFunction  qout(hard[2]->momentum(), hard[2]->dataPtr(), outgoing);

    for (unsigned int ix = 0; ix < 2; ++ix) {
      qin.reset(ix);         fin.push_back(qin);
      glin.reset(2*ix);      gin.push_back(glin);
      qout.reset(ix);        fout.push_back(qout);
    }
    qgME(fin, gin, fout, Aout, true);
  }
  else if (hard[0]->id() >= -6 && hard[0]->id() < 0 && hard[1]->id() == ParticleID::g) {
    // qbar g -> qbar γ
    std::vector<SpinorBarWaveFunction>  ain;
    std::vector<VectorWaveFunction>     gin;
    std::vector<SpinorWaveFunction>     aout;

    SpinorBarWaveFunction qbin(hard[0]->momentum(), hard[0]->dataPtr(), incoming);
    VectorWaveFunction    glin(hard[1]->momentum(), hard[1]->dataPtr(), incoming);
    SpinorWaveFunction    qbout(hard[2]->momentum(), hard[2]->dataPtr(), outgoing);

    for (unsigned int ix = 0; ix < 2; ++ix) {
      qbin.reset(ix);        ain.push_back(qbin);
      glin.reset(2*ix);      gin.push_back(glin);
      qbout.reset(ix);       aout.push_back(qbout);
    }
    qbargME(ain, gin, aout, Aout, true);
  }
  else {
    // q qbar -> γ g
    std::vector<SpinorWaveFunction>     fin;
    std::vector<SpinorBarWaveFunction>  ain;
    std::vector<VectorWaveFunction>     gout;

    SpinorWaveFunction     qin (hard[0]->momentum(), hard[0]->dataPtr(), incoming);
    SpinorBarWaveFunction  qbin(hard[1]->momentum(), hard[1]->dataPtr(), incoming);
    VectorWaveFunction     glout(hard[2]->momentum(), hard[2]->dataPtr(), outgoing);

    for (unsigned int ix = 0; ix < 2; ++ix) {
      qin.reset(ix);         fin.push_back(qin);
      qbin.reset(ix);        ain.push_back(qbin);
      glout.reset(2*ix);     gout.push_back(glout);
    }
    qqbarME(fin, ain, gout, Aout, true);
  }

  // Build the hard vertex (guard against missing spinInfo)
  HardVertexPtr hv = new_ptr(HardVertex());
  hv->ME(_me);
  for (unsigned int ix = 0; ix < 4; ++ix) {
    if (hard[ix] && hard[ix]->spinInfo())
      hard[ix]->spinInfo()->productionVertex(hv);
  }
}
