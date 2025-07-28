// -*- C++ -*-  MEGammaGamma2QQZprime.cc
// ---------------------------------------------------------------
// gamma gamma -> ZPrime j j  (internal 2 -> 3 ME for Herwig 7.x)
// ---------------------------------------------------------------
#include <complex>
#include <cmath>
#include "MEGammaGamma2QQZprime.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Models/HiddenValley/HiddenValleyFFZPrimeVertex.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/Cuts/Cuts.h"

using namespace Herwig;

MEGammaGamma2QQZprime::MEGammaGamma2QQZprime() :
  minFlavor_(1), maxFlavor_(6), scaleOption_(0), mzp_(), wzp_()
{}

DescribeClass<MEGammaGamma2QQZprime, HwMEBase>
describeHerwigMEGammaGamma2QQZprime("Herwig::MEGammaGamma2QQZprime",
                                    "HwMEGammaGamma.so");

void MEGammaGamma2QQZprime::Init() {
  static ClassDocumentation<MEGammaGamma2QQZprime>
      documentation("Internal ME for gamma gamma -> j j ZPrime");

  static Parameter<MEGammaGamma2QQZprime,int>
      pMin("MinFlavour",
           "Smallest quark flavour to include (1=d … 6=t).",
           &MEGammaGamma2QQZprime::minFlavor_, 1, 1, 6,
           false, false, Interface::limited);

  static Parameter<MEGammaGamma2QQZprime,int>
      pMax("MaxFlavour",
           "Largest quark flavour to include.",
           &MEGammaGamma2QQZprime::maxFlavor_, 6, 1, 6,
           false, false, Interface::limited);

   static Switch<MEGammaGamma2QQZprime,unsigned int> interfaceScaleOption
     ("Scale",
      "Scale option for gamma gamma -> j j Z' matrix element",
      &MEGammaGamma2QQZprime::scaleOption_, 0, false, false);

   static SwitchOption interfaceScaleOptionSHat
     (interfaceScaleOption,
      "SHat",
      "Fixed scale SHat",
      0);

   static SwitchOption interfaceScaleOptionMHat
     (interfaceScaleOption,
      "InvMass",
      "Invariant mass scale",
      1);
}

Energy2 MEGammaGamma2QQZprime::scale() const {
  if(scaleOption_==0)
    return sHat();
  else if(scaleOption_==1)
    return sqr(mePartonData()[2]->mass()+mePartonData()[3]->mass()+mePartonData()[4]->mass());
  else
  assert(false);
}

int MEGammaGamma2QQZprime::nDim() const {
  return 6; // // 5 for a 2->3 map +1 for the flavour choice
}

unsigned int MEGammaGamma2QQZprime::orderInAlphaS() const {
  return 0;
}

unsigned int MEGammaGamma2QQZprime::orderInAlphaEW() const {
  return 3;
}

IBPtr MEGammaGamma2QQZprime::clone() const {
  return new_ptr(*this);
}

IBPtr MEGammaGamma2QQZprime::fullclone() const {
  return new_ptr(*this);
}

void MEGammaGamma2QQZprime::setKinematics() {
  HwMEBase::setKinematics();
}

void MEGammaGamma2QQZprime::persistentOutput(PersistentOStream &os) const {
  os << minFlavor_ << maxFlavor_ << ZprimeVertex_ << Zprime_
     << quark1_ << quark2_ << photon_ << model_ << ounit(mzp_,GeV) << ounit(wzp_,GeV);
}

void MEGammaGamma2QQZprime::persistentInput(PersistentIStream &is, int) {
  is >> minFlavor_ >> maxFlavor_ >> ZprimeVertex_ >> Zprime_
     >> quark1_ >> quark2_ >> photon_ >> model_ >> iunit(mzp_,GeV) >> iunit(wzp_,GeV);
}

void MEGammaGamma2QQZprime::doinit() {
  HwMEBase::doinit();

  if (minFlavor_ > maxFlavor_)
    throw InitException() << "MinFlavour > MaxFlavour" << Exception::runerror;

  photon_  = getParticleData(ParticleID::gamma);
  Zprime_  = getParticleData(32); // PDG 32 for ZPrime

  model_ = dynamic_ptr_cast<tcHiddenValleyPtr>(generator()->standardModel());
  if (!model_)
    throw InitException() << "HiddenValleyModel required." << Exception::runerror;

  ZprimeVertex_ = model_->FFZPVertex();
  if (!ZprimeVertex_)
    throw InitException() << "ZPrime vertex missing." << Exception::runerror;
  ZprimeVertex_->init();

  mzp_ = Zprime_->mass();
  wzp_ = Zprime_->width();
}

void MEGammaGamma2QQZprime::rebind(const TranslationMap & trans) {
  ZprimeVertex_ = trans.translate(ZprimeVertex_);
  Zprime_       = trans.translate(Zprime_);
  photon_       = trans.translate(photon_);
  quark1_       = trans.translate(quark1_);
  quark2_       = trans.translate(quark2_);
  HwMEBase::rebind(trans);
}

IVector MEGammaGamma2QQZprime::getReferences() {
  IVector ret = HwMEBase::getReferences();
  ret.push_back(ZprimeVertex_);
  ret.push_back(Zprime_);
  ret.push_back(photon_);
  ret.push_back(quark1_);
  ret.push_back(quark2_);
  return ret;
}

bool MEGammaGamma2QQZprime::generateKinematics(const double * r) {
  jacobian(1.0);

  int nf = minFlavor_ + int((maxFlavor_ - minFlavor_ + 1) * r[5]);
  if (nf < 1 || nf > 6) return false;
  tPDPtr q     = getParticleData(nf);
  tPDPtr qbar  = q->CC();
  Energy mq    = q->mass();

  Energy sqrtS = sqrt(sHat());
  if (sqrtS < 2.0 * mq + mzp_) return false;

  meMomenta()[2].setMass(mq);
  meMomenta()[3].setMass(mq);
  meMomenta()[4].setMass(mzp_);

  Energy2 mmin = sqr(mq + mzp_);
  Energy2 mmax = sqr(sqrtS - mq);
  if (mmin >= mmax) return false;

  double rho = r[1];
  Energy2 m132 = mmin + rho * (mmax - mmin);
  Energy m13 = sqrt(m132);

  Energy p1 = ZERO;
  try {
    p1 = SimplePhaseSpace::getMagnitude(sHat(), m13, mq);
  } catch (ImpossibleKinematics & e) {
    return false;
  }

  Energy p2 = ZERO;
  try {
    p2 = SimplePhaseSpace::getMagnitude(m132, mq, mzp_);
  } catch (ImpossibleKinematics & e) {
    return false;
  }

  Energy ptmin = lastCuts().minKT(qbar);
  double ctmin = -1.0, ctmax = 1.0;
  if (ptmin > ZERO) {
    double ctm = 1.0 - sqr(ptmin / p1);
    if (ctm <= 0.0) return false;
    ctmin = max(ctmin, -sqrt(ctm));
    ctmax = min(ctmax,  sqrt(ctm));
  }

  double cos1 = getCosTheta(ctmin, ctmax, r[0]);
  double sin1 = sqrt(1.0 - sqr(cos1));
  double phi1 = Constants::twopi * r[2];

  Lorentz5Momentum p13(sin1 * p1 * cos(phi1),
                       sin1 * p1 * sin(phi1),
                       cos1 * p1,
                       sqrt(sqr(p1) + m132),
                       m13);

  meMomenta()[3].setVect(-p13.vect());
  meMomenta()[3].setMass(mq);
  meMomenta()[3].rescaleEnergy();

  bool ok = Kinematics::twoBodyDecay(p13, mq, mzp_, -1.0 + 2.0 * r[3], Constants::twopi * r[4], meMomenta()[2], meMomenta()[4]);
  if (!ok) return false;

  vector<LorentzMomentum> out = {meMomenta()[2],meMomenta()[3],meMomenta()[4]};
  tcPDVector tout = { q, qbar, Zprime_ };
  if (!lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]))
    return false;
  jacobian(jacobian() * p1 * p2 * Constants::pi / sHat() * (mmax - mmin)/GeV2);
  return true;
}


CrossSection MEGammaGamma2QQZprime::dSigHatDR() const {
  return me2()*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
}

void MEGammaGamma2QQZprime::getDiagrams() const {
  for (int f = minFlavor_; f <= maxFlavor_; ++f) {
    tcPDPtr q  = getParticleData(f);
    tcPDPtr qb = q->CC();

    /* t-channel, central, u-channel + gamma1/gamma2 interchanged, 6 total */
    add(new_ptr((Tree2toNDiagram(3), photon_, qb, photon_,
                 1,q,4,q,2,qb,4,Zprime_,-1)));
    add(new_ptr((Tree2toNDiagram(4), photon_, qb, qb, photon_,
                 1,q,3,qb,2,Zprime_,-2)));
    add(new_ptr((Tree2toNDiagram(3), photon_, qb, photon_,
                 1,q,2,qb,5,qb,5,Zprime_,-3)));
    add(new_ptr((Tree2toNDiagram(3), photon_, q, photon_,
                 2,q,4,q,1,qb,4,Zprime_,-4)));
    add(new_ptr((Tree2toNDiagram(4), photon_, q, q, photon_,
                 3,q,1,qb,2,Zprime_,-5)));
    add(new_ptr((Tree2toNDiagram(3), photon_, q, photon_,
                 2,q,1,qb,5,qb,5,Zprime_,-6)));
  }
}

double MEGammaGamma2QQZprime::me2() const {
  double me = 0.0;

  if (mePartonData()[0]->id() == ParticleID::gamma) {
    VectorWaveFunction g1w(meMomenta()[0], mePartonData()[0], incoming);
    VectorWaveFunction g2w(meMomenta()[1], mePartonData()[1], incoming);
    SpinorBarWaveFunction qw(meMomenta()[2], mePartonData()[2], outgoing);
    SpinorWaveFunction qbarw(meMomenta()[3], mePartonData()[3], outgoing);

    // ZPrime is a vector, not scalar!
    VectorWaveFunction zprime(meMomenta()[4], mePartonData()[4], outgoing);

    vector<VectorWaveFunction> g1, g2;
    vector<SpinorBarWaveFunction> q;
    vector<SpinorWaveFunction> qbar;

    for (unsigned int ix = 0; ix < 2; ++ix) {
      g1w.reset(2 * ix); g1.push_back(g1w);
      g2w.reset(2 * ix); g2.push_back(g2w);
      qw.reset(ix);      q.push_back(qw);
      qbarw.reset(ix);   qbar.push_back(qbarw);
    }

    // Use the updated ME function
    me = gammagammaME(g1, g2, q, qbar, zprime, 0);
  }

  return me * sHat() * UnitRemoval::InvE2;
}


double MEGammaGamma2QQZprime::gammagammaME(
  vector<VectorWaveFunction> &g1,
  vector<VectorWaveFunction> &g2,
  vector<SpinorBarWaveFunction> &qbar,
  vector<SpinorWaveFunction> &q,
  VectorWaveFunction &zprime,
  unsigned int iflow) const {

  Energy2 mt = scale();
  Energy mass = q[0].mass();

  if (iflow != 0) {
    me_.reset(ProductionMatrixElement(PDT::Spin1, PDT::Spin1,
                                      PDT::Spin1Half, PDT::Spin1Half,
                                      PDT::Spin1));
  }

  double total = 0.0;
  double sumflow[1] = {0.0};
  double sumdiag[6] = {0.0};

  Complex diag[6];

  for (unsigned int ihel1 = 0; ihel1 < 2; ++ihel1) {
    for (unsigned int ihel2 = 0; ihel2 < 2; ++ihel2) {
      for (unsigned int ohel1 = 0; ohel1 < 2; ++ohel1) {
        for (unsigned int ohel2 = 0; ohel2 < 2; ++ohel2) {

          // Compute each of the 6 diagrams separately
          diag[0] = ZprimeVertex_->evaluate(mt,
                       ZprimeVertex_->evaluate(mt, 1, q[ohel1].particle()->CC(),
                                               q[ohel1], g1[ihel1], mass),
                       qbar[ohel2], g2[ihel2]);

          diag[1] = ZprimeVertex_->evaluate(mt,
                       ZprimeVertex_->evaluate(mt, 1, q[ohel1].particle()->CC(),
                                               q[ohel1], zprime, mass),
                       qbar[ohel2], g2[ihel2]);

          diag[2] = ZprimeVertex_->evaluate(mt,
                       ZprimeVertex_->evaluate(mt, 1, q[ohel1].particle()->CC(),
                                               q[ohel1], g2[ihel2], mass),
                       qbar[ohel2], g1[ihel1]);

          diag[3] = ZprimeVertex_->evaluate(mt,
                       ZprimeVertex_->evaluate(mt, 1, q[ohel1].particle()->CC(),
                                               q[ohel1], g2[ihel2], mass),
                       qbar[ohel2], g1[ihel1]);

          diag[4] = ZprimeVertex_->evaluate(mt,
                       ZprimeVertex_->evaluate(mt, 1, q[ohel1].particle()->CC(),
                                               q[ohel1], zprime, mass),
                       qbar[ohel2], g1[ihel1]);

          diag[5] = ZprimeVertex_->evaluate(mt,
                       ZprimeVertex_->evaluate(mt, 1, q[ohel1].particle()->CC(),
                                               q[ohel1], g1[ihel1], mass),
                       qbar[ohel2], g2[ihel2]);

          Complex amp = Complex(0.0, 0.0);
          for (int i = 0; i < 6; ++i) {
            amp += diag[i];
            sumdiag[i] += norm(diag[i]);
          }

          sumflow[0] += norm(amp);
          total += norm(amp);

          if (iflow != 0)
            me_(2*ihel1, 2*ihel2, ohel1, ohel2, 0) = amp;
        }
      }
    }
  }

  // Choose flow and diagram
  flow_ = 1;

  // Choose diagram probabilistically from sumdiag[]
  double prob = UseRandom::rnd();
  double totalNorm = 0.0;
  for (int i = 0; i < 6; ++i) totalNorm += sumdiag[i];

  for (int i = 0; i < 6; ++i) {
    if (prob <= sumdiag[i]/totalNorm) {
      diagram_ = i + 1;
      break;
    }
    prob -= sumdiag[i]/totalNorm;
  }

  return total / 4.0;
}

Selector<const ColourLines *>
MEGammaGamma2QQZprime::colourGeometries(tcDiagPtr diag) const {
  static const ColourLines cLines[6] = {
    ColourLines("5 4 -2 -6"),   // diagram -1
    ColourLines("5 -2 -3 -6"),  // diagram -2
    ColourLines("4 -2 -5 -6"),  // diagram -3
    ColourLines("5 4 2 -6"),    // diagram -4
    ColourLines("5 3 2 -6"),    // diagram -5
    ColourLines("4 2 -5 -6")    // diagram -6
  };

  int id = -diag->id();  // diagram IDs are negative
  assert(id >= 1 && id <= 6);

  Selector<const ColourLines *> sel;
  sel.insert(1.0, &cLines[id - 1]);
  return sel;
}


Selector<MEBase::DiagramIndex>
MEGammaGamma2QQZprime::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for (DiagramIndex i = 0; i < diags.size(); ++i) {
    if (diags[i]->id() == -int(diagram_))
      sel.insert(1.0, i);
    else
      sel.insert(0.0, i);
  }
  return sel;
}

void MEGammaGamma2QQZprime::constructVertex(tSubProPtr sub) {
  // Extract the particles
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  for (unsigned int ix = 0; ix < 3; ++ix)
    hard.push_back(sub->outgoing()[ix]);

  // Ensure incoming photon ordering
  if (hard[0]->id() < 0)
    std::swap(hard[0], hard[1]);

  // Identify ZPrime in final state and move to slot 4
  if (hard[2]->id() == 32) std::swap(hard[2], hard[4]);
  if (hard[3]->id() == 32) std::swap(hard[3], hard[4]);

  // Ensure quark is a particle, antiquark is an antiparticle
  if (hard[2]->id() < 0)
    std::swap(hard[2], hard[3]);

  // Photon-initiated gamma gamma -> q qbar ZPrime process
  if (hard[0]->id() == ParticleID::gamma) {
    std::vector<VectorWaveFunction> g1, g2;
    std::vector<SpinorBarWaveFunction> q;
    std::vector<SpinorWaveFunction> qbar;

    // Construct wavefunctions
    VectorWaveFunction(g1, hard[0], incoming, false, true, true);
    VectorWaveFunction(g2, hard[1], incoming, false, true, true);
    SpinorBarWaveFunction(q,    hard[2], outgoing, true, true);
    SpinorWaveFunction(  qbar, hard[3], outgoing, true, true);
    std::vector<VectorWaveFunction> zp;
    VectorWaveFunction(zp, hard[4], outgoing, true, true, true);

    g1[1] = g1[2];
    g2[1] = g2[2];

    gammagammaME(g1, g2, q, qbar, zp[1], flow_);
  }
  else {
    assert(false);  // Unexpected incoming particles
  }

  // Build the vertex
  HardVertexPtr hardvertex = new_ptr(HardVertex());
  hardvertex->ME(me_);

  // Attach all particles to the vertex
  for (unsigned int ix = 0; ix < 5; ++ix)
    hard[ix]->spinInfo()->productionVertex(hardvertex);
}
