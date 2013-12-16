// -*- C++ -*-
//
// HEJMEBase.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HEJMEBase class.
//

#include "HEJMEBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/ColourBasis.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/DiagramDrawer.h"

using namespace Herwig;

HEJMEBase::HEJMEBase()
  : theJets(0) {}

HEJMEBase::~HEJMEBase() {}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void HEJMEBase::factory(Ptr<HEJFactory>::ptr f) {
  theFactory = f;
}

void HEJMEBase::persistentOutput(PersistentOStream & os) const {
  os << theFirstIncoming << theSecondIncoming << theNGluons << theFactory;
}

void HEJMEBase::persistentInput(PersistentIStream & is, int) {
  is >> theFirstIncoming >> theSecondIncoming >> theNGluons >> theFactory;
}

Selector<MEBase::DiagramIndex> 
HEJMEBase::diagrams(const DiagramVector & v) const {
  assert(v.size() == 1);
  Selector<MEBase::DiagramIndex> sel;
  sel.insert(1.,0);
  return sel;
}

void HEJMEBase::setXComb(tStdXCombPtr xc) {
  if ( !theJets )
    theFactory->makeCMultijet(xc);
  if ( !xc->hasMeta(HEJMetaKeys::Jets) )
    xc->meta(HEJMetaKeys::Jets,*theJets);
  MatchboxMEBase::setXComb(xc);
}

double HEJMEBase::me2() const {

  assert(lastXCombPtr()->hasMeta(HEJMetaKeys::Jets));
  CMultijet& jets = lastXCombPtr()->meta<CMultijet>(HEJMetaKeys::Jets);

  CLHEPConverter convert;

  convert(meMomenta()[1],pa);
  convert(meMomenta()[0],pb);

  pOut.resize(meMomenta().size()-2);
  size_t n = meMomenta().size();
  for ( size_t k = 2; k < n; ++k ) {
    convert(meMomenta()[n+1-k],pOut[k-2]);
  }

  /*
  cerr << name() << "::me2() evaluating ...\n";
  cerr << "using PS point\n"
       << "HEJ conventions:\n";
  cerr << pa.x() << " "
       << pa.y() << " "
       << pa.z() << " "
       << pa.t() << "\n"
       << pb.x() << " "
       << pb.y() << " "
       << pb.z() << " "
       << pb.t() << "\n";
  for ( size_t k = 0; k < meMomenta().size()-2; ++k )
    cerr << pOut[k].x() << " "
	 << pOut[k].y() << " "
	 << pOut[k].z() << " "
	 << pOut[k].t() << "\n";
  cerr << "Herwig++ conventions\n";
  for ( size_t k = 0; k < meMomenta().size(); ++k )
    cerr << meMomenta()[k]/GeV << "\n";
  cerr << flush;
  */

  jets.reset();
  jets.setAPtype(mePartonData()[1]->id());
  jets.setBPtype(mePartonData()[0]->id());
  jets.setXA(lastX2());
  jets.setXB(lastX1());

  double res = jets.UberME2(pOut,pa,pb,true,true);

  return pow(lastSHat()/GeV2,(double)(meMomenta().size()-4))*res;

}

Selector<const ColourLines *>
HEJMEBase::colourGeometries(tcDiagPtr diag) const {

  Ptr<Tree2toNDiagram>::tcptr diagptr = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::tcptr>(diag);

  assert(lastXCombPtr()->hasMeta(HEJMetaKeys::Jets));
  CMultijet& jets = lastXCombPtr()->meta<CMultijet>(HEJMetaKeys::Jets);
  jets.reset();
  jets.setAPtype(mePartonData()[1]->id());
  jets.setBPtype(mePartonData()[0]->id());
  jets.setXA(lastX2());
  jets.setXB(lastX1());
  jets.SetupLHInformation(jets.LHParticlesPtr(),&pOut,&pa,&pb);

  size_t n = meMomenta().size();

  /*
  cerr << name() << " got LesHouchesParticles (Herwig++ conventions)\n";
  cerr << jets.LesHouchesParticles()[1].PID << " "
       << jets.LesHouchesParticles()[1].COL1 << " "
       << jets.LesHouchesParticles()[1].COL2 << " "
       << jets.LesHouchesParticles()[1].p.x() << " "
       << jets.LesHouchesParticles()[1].p.y() << " "
       << jets.LesHouchesParticles()[1].p.z() << " "
       << jets.LesHouchesParticles()[1].p.t() << "\n"
       << jets.LesHouchesParticles()[0].PID << " "
       << jets.LesHouchesParticles()[0].COL1 << " "
       << jets.LesHouchesParticles()[0].COL2 << " "
       << jets.LesHouchesParticles()[0].p.x() << " "
       << jets.LesHouchesParticles()[0].p.y() << " "
       << jets.LesHouchesParticles()[0].p.z() << " "
       << jets.LesHouchesParticles()[0].p.t() << "\n";
  for ( size_t k = n-1; k > 1; --k )
    cerr << jets.LesHouchesParticles()[k].PID << " "
	 << jets.LesHouchesParticles()[k].COL1 << " "
	 << jets.LesHouchesParticles()[k].COL2 << " "
	 << jets.LesHouchesParticles()[k].p.x() << " "
	 << jets.LesHouchesParticles()[k].p.y() << " "
	 << jets.LesHouchesParticles()[k].p.z() << " "
	 << jets.LesHouchesParticles()[k].p.t() << " "
	 << jets.LesHouchesParticles()[k].p.rapidity()
	 << "\n";
  cerr << flush;
  */

  multimap<int,int> colours;

  if ( jets.LesHouchesParticles()[1].COL1 )
    colours.insert(make_pair(jets.LesHouchesParticles()[1].COL1,0));
  if ( jets.LesHouchesParticles()[1].COL2 )
    colours.insert(make_pair(-jets.LesHouchesParticles()[1].COL2,0));
  if ( jets.LesHouchesParticles()[0].COL1 )
    colours.insert(make_pair(jets.LesHouchesParticles()[0].COL1,1));
  if ( jets.LesHouchesParticles()[0].COL2 )
    colours.insert(make_pair(-jets.LesHouchesParticles()[0].COL2,1));
  for ( size_t k = n-1; k > 1; --k ) {
    if ( jets.LesHouchesParticles()[k].COL1 )
      colours.insert(make_pair(jets.LesHouchesParticles()[k].COL1,n-k+1));
    if ( jets.LesHouchesParticles()[k].COL2 )
      colours.insert(make_pair(-jets.LesHouchesParticles()[k].COL2,n-k+1));
  }

  list<list<pair<int,bool> > > flow;

  while ( !colours.empty() ) {

    pair<int,bool> source(-1,false);
    pair<int,bool> sink(-1,false);

    source.first = colours.begin()->second;
    source.second = colours.begin()->first > 0;
    int id = colours.begin()->first;
    colours.erase(colours.begin());

    map<int,int>::iterator partner =
      colours.find(-id);
    if ( partner == colours.end() )
      partner = colours.find(id);
    assert(partner != colours.end());

    sink.first = partner->second;
    sink.second = partner->first > 0;
    colours.erase(partner);

    list<pair<int,bool> > line =
      ColourBasis::colouredPath(source,sink,diagptr);
    assert(!line.empty());

    flow.push_back(line);

  }

  string res = ColourBasis::cfstring(flow);

  /*
  cerr << name() << " found colour flow: "
       << res << "\n" << "for diagram\n";
  DiagramDrawer::drawDiag(cerr,*diagptr);
  cerr << "\n" << flush;
  */

  assert(res != "");

  currentColourLines.reset(res);

  Selector<const ColourLines *> ret;
  ret.insert(1.,&currentColourLines);

  return ret;

}

Energy2 HEJMEBase::factorizationScale() const {
  assert(lastXCombPtr()->hasMeta(HEJMetaKeys::Jets));
  CMultijet& jets = lastXCombPtr()->meta<CMultijet>(HEJMetaKeys::Jets);
  return sqr(jets.renormalizationScale())*GeV2;
}

Energy2 HEJMEBase::renormalizationScale() const {
  assert(lastXCombPtr()->hasMeta(HEJMetaKeys::Jets));
  CMultijet& jets = lastXCombPtr()->meta<CMultijet>(HEJMetaKeys::Jets);
  return sqr(jets.renormalizationScale())*GeV2;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<HEJMEBase,Herwig::MatchboxMEBase>
  describeHerwigHEJMEBase("Herwig::HEJMEBase", "HwMatchbox.so HwHEJ.so");

void HEJMEBase::Init() {

  static ClassDocumentation<HEJMEBase> documentation
    ("HEJMEBase is the base class for HEJ matrix elements.");

}

