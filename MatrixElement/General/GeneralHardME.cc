// -*- C++ -*-
//
// GeneralHardME.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralHardME class.
//

#include "GeneralHardME.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"
#include "ThePEG/Utilities/Debug.h"
#include <numeric>

using namespace Herwig;

GeneralHardME::GeneralHardME() : incoming_(0, 0), outgoing_(0, 0),
				 diagrams_(0), numberOfDiagrams_(0), 
				 colour_(0), numberOfFlows_(0) , 
				 debug_(false), scaleChoice_(0),
				 scaleFactor_(1.) {
  massOption(vector<unsigned int>(2,1));
}
  
void GeneralHardME::setProcessInfo(const vector<HPDiagram> & alldiagrams,
				   ColourStructure colour,
				   bool debug, unsigned int scaleOption,
				   double scaleFactor) {
  // external particles
  incoming_ = alldiagrams.at(0).incoming;
  outgoing_ = alldiagrams.at(0).outgoing;
  diagrams_ = alldiagrams;
  numberOfDiagrams_ = alldiagrams.size();
  // debug option
  debug_ = debug;
  // scale choice
  scaleChoice_ = scaleOption; 
  scaleFactor_ = scaleFactor; 
  // OffShell options
  pair<bool, bool> offshell(make_pair(false, false));
  vector<unsigned int> mopt(2,1);
  if( getParticleData(outgoing_.first )->widthGenerator() ||
      getParticleData(outgoing_.first )-> massGenerator()) {
    offshell.first = true;
    mopt[0] = 2;
  }
  if( getParticleData(outgoing_.second)->widthGenerator() ||
      getParticleData(outgoing_.second)-> massGenerator() ) {
    offshell.second = true;
    mopt[1] = 2;
  }
  if(outgoing_.first  == incoming_.first ||
     outgoing_.first  == incoming_.second )
    mopt[0] = 0;
  if(outgoing_.second == incoming_.first ||
     outgoing_.second == incoming_.second )
    mopt[1] = 0;
  massOption(mopt);
  if( offshell.first == true &&  offshell.second == true &&
      abs(outgoing_.first) == abs(outgoing_.second)  )
    rescalingOption(3);
  // colour structure
  colourStructure_ = colour;
  switch (colour) {
  // colour neutral process
  case Colour11to11:
    colour_ = vector<DVector>(1,DVector(1,1.));
    numberOfFlows_ = 1;
    break;
  // colour neutral -> 3 3bar or swap process
  case Colour11to33bar: case Colour33barto11:
    colour_ = vector<DVector>(1,DVector(1,3.));
    numberOfFlows_ = 1;
    break;
  // colour neutral -> 8 8    process or swap
  case Colour11to88: case Colour88to11 :
    colour_ = vector<DVector>(1,DVector(1,8.));
    numberOfFlows_ = 1;
    break;
  //    colour 33    -> 33 or 3bar3bar -> 3bar3bar process
  // or colour 33bar -> 33bar
  case Colour33to33: case Colour3bar3barto3bar3bar:
  case Colour33barto33bar:
    colour_ = vector<DVector>(4, DVector(4, 0.));
    colour_[0][0] = colour_[1][1] = 2.; 
    colour_[2][2] = colour_[3][3] = 9.;
    colour_[0][1] = colour_[1][0] = -2./3.;
    colour_[0][2] = colour_[2][0] = 0.;
    colour_[0][3] = colour_[3][0] = 4.;
    colour_[1][2] = colour_[2][1] = 4.;
    colour_[1][3] = colour_[3][1] = 0.;
    colour_[2][3] = colour_[3][2] = 3.;
    numberOfFlows_ = 4;
    break;
  // colour 3 3bar -> 6 6bar
  case Colour33barto66bar: case Colour33barto6bar6:
    colour_ = vector<DVector>(10, DVector(10, 0.));
    // diagonals
    for(unsigned int ix=0;ix<4;++ix){
      colour_[ix][ix]     = 1.5;
      colour_[ix+4][ix+4] = 27./16.;
    }
    colour_[8][8]=colour_[9][9]=27./4.;
    // 0
    colour_[0][1] = colour_[1][0] = 0.5;
    colour_[0][2] = colour_[2][0] = 0.5;
    colour_[0][3] = colour_[3][0] = 0.;
    colour_[0][4] = colour_[4][0] = 3./2.;
    colour_[0][5] = colour_[5][0] = 1./2.;
    colour_[0][6] = colour_[6][0] = 1./2.;
    colour_[0][7] = colour_[7][0] = 0.;
    // 1
    colour_[1][2] = colour_[2][1] = 0.;
    colour_[1][3] = colour_[3][1] = 0.5;
    colour_[1][4] = colour_[4][1] = 1./2.;
    colour_[1][5] = colour_[5][1] = 3./2.;
    colour_[1][6] = colour_[6][1] = 0.;
    colour_[1][7] = colour_[7][1] = 1./2.;
    // 2
    colour_[2][3] = colour_[3][2] = 0.5;
    colour_[2][4] = colour_[4][2] = 1./2.;
    colour_[2][5] = colour_[5][2] = 0.;
    colour_[2][6] = colour_[6][2] = 3./2.;
    colour_[2][7] = colour_[7][2] = 1./2.;
    // 3
    colour_[3][4] = colour_[4][3] = 0.;
    colour_[3][5] = colour_[5][3] = 1./2.;
    colour_[3][6] = colour_[6][3] = 1./2.;
    colour_[3][7] = colour_[7][3] = 3./2.;
    // 4
    colour_[4][5] = colour_[5][4] = 9./16.;
    colour_[4][6] = colour_[6][4] = 9./16.;
    colour_[4][7] = colour_[7][4] = 3./16.;
    colour_[4][8] = colour_[8][4] = 9./8.;
    colour_[4][9] = colour_[9][4] = 3./8.;
    // 5
    colour_[5][6] = colour_[6][5] = 3./16.;
    colour_[5][7] = colour_[7][5] = 9./16.;
    colour_[5][8] = colour_[8][5] = 3./8.;
    colour_[5][9] = colour_[9][5] = 9./8.;
    //6
    colour_[6][7] = colour_[7][6] = 9./16.;
    colour_[6][8] = colour_[8][6] = 3./8.;
    colour_[6][9] = colour_[9][6] = 9./8.;
    // 7
    colour_[7][8] = colour_[8][7] = 9./8.;
    colour_[7][9] = colour_[9][7] = 3./8.;
    // 8
    colour_[8][9] = colour_[9][8] = 9./4.;
    numberOfFlows_ = 10;
    break;
  // colour 33bar -> 88, 88 -> 33bar
  case Colour33barto88: case Colour88to33bar:
  case Colour38to83: case Colour38to38:
  case Colour3bar8to83bar: case Colour3bar8to3bar8:
    colour_ = vector<DVector>(3, DVector(3, 0.));
    colour_[0][0] = colour_[1][1] = 16./3.;
    colour_[0][1] = colour_[1][0] = -2./3.;
    colour_[2][2]=24.;
    colour_[0][2] = colour_[2][0] = 4.;
    colour_[1][2] = colour_[2][1] = 4.;
    numberOfFlows_ = 3;
    break;
  case Colour88to88: {
    colour_ = vector<DVector>(9, DVector(9, 0.));
    double on1=19./6.,on2=64.;
    double off1=-1./3.,off2=2./3.,off3=8.,off4=16./3.,off5=-2./3.;
    // octet exchange
    colour_[0][0] = colour_[1][1] = colour_[2][2]=on1;
    colour_[3][3] = colour_[4][4] = colour_[5][5]=on1;
    colour_[0][1] = colour_[1][0] = off1;
    colour_[0][2] = colour_[2][0] = off1;
    colour_[0][3] = colour_[3][0] = off1;
    colour_[0][4] = colour_[4][0] = off1;
    colour_[0][5] = colour_[5][0] = off2;
    colour_[1][2] = colour_[2][1] = off1;
    colour_[1][3] = colour_[3][1] = off2;
    colour_[1][4] = colour_[4][1] = off1;
    colour_[1][5] = colour_[5][1] = off1;
    colour_[2][3] = colour_[3][2] = off1;
    colour_[2][4] = colour_[4][2] = off2;
    colour_[2][5] = colour_[5][2] = off1;
    colour_[3][4] = colour_[4][3] = off1;
    colour_[3][5] = colour_[5][3] = off1;
    colour_[4][5] = colour_[5][4] = off1;
    // singlet exchange
    colour_[6][6] = colour_[7][7] = colour_[8][8] = on2;
    colour_[6][7] = colour_[7][6] = off3;
    colour_[6][8] = colour_[8][6] = off3;
    colour_[7][8] = colour_[8][7] = off3;
    // octet singlet interference
    colour_[0][6] = colour_[6][0] = off4;
    colour_[1][6] = colour_[6][1] = off4;
    colour_[2][6] = colour_[6][2] = off5;
    colour_[3][6] = colour_[6][3] = off4;
    colour_[4][6] = colour_[6][4] = off5;
    colour_[5][6] = colour_[6][5] = off4;
    colour_[0][7] = colour_[7][0] = off5;
    colour_[1][7] = colour_[7][1] = off4;
    colour_[2][7] = colour_[7][2] = off4;
    colour_[3][7] = colour_[7][3] = off4;
    colour_[4][7] = colour_[7][4] = off4;
    colour_[5][7] = colour_[7][5] = off5;
    colour_[0][8] = colour_[8][0] = off4;
    colour_[1][8] = colour_[8][1] = off5;
    colour_[2][8] = colour_[8][2] = off4;
    colour_[3][8] = colour_[8][3] = off5;
    colour_[4][8] = colour_[8][4] = off4;
    colour_[5][8] = colour_[8][5] = off4;
    numberOfFlows_ = 9;
    break;
  }
  case Colour33barto18   : case Colour33barto81   : 
  case Colour38to13      : case Colour38to31      :
  case Colour3bar8to13bar: case Colour3bar8to3bar1:
    colour_ = vector<DVector>(1,DVector(1,4.));
    numberOfFlows_ = 1;
    break;
  case Colour88to18 : case Colour88to81:
    colour_ = vector<DVector>(1,DVector(1,24.));
    numberOfFlows_ = 1;
    break;
  case Colour88to66bar:
    colour_ = vector<DVector>(12, DVector(12, 0.));
    // diagonals
    for(unsigned int ix=0;ix<12;++ix) colour_[ix][ix] = 4.;
    //  1 1 block
    colour_[ 0][ 1] = colour_[ 1][ 0] =  4./3.;
    colour_[ 0][ 2] = colour_[ 2][ 0] =  4./3.;
    colour_[ 0][ 3] = colour_[ 3][ 0] =  0.5  ;
    colour_[ 1][ 2] = colour_[ 2][ 1] =  0.5  ;
    colour_[ 1][ 3] = colour_[ 3][ 1] =  4./3.;
    colour_[ 2][ 3] = colour_[ 3][2] =   4./3.;
    // 1 2 and 2 1 blocks
    colour_[ 0][ 4] = colour_[ 4][ 0] =  4./3.;
    colour_[ 1][ 5] = colour_[ 5][ 1] =  4./3.;
    colour_[ 2][ 6] = colour_[ 6][ 2] =  4./3.;
    colour_[ 3][ 7] = colour_[ 7][ 3] =  4./3.;
    colour_[ 0][ 7] = colour_[ 7][ 0] = -1./6.;
    colour_[ 1][ 6] = colour_[ 6][ 1] = -1./6.;
    colour_[ 2][ 5] = colour_[ 5][ 2] = -1./6.;
    colour_[ 3][ 4] = colour_[ 4][ 3] = -1./6.;
    // 1 3 and 3 1 blocks
    colour_[ 0][11] = colour_[11][ 0] =  4./3.;
    colour_[ 1][10] = colour_[10][ 1] =  4./3.;
    colour_[ 2][ 9] = colour_[ 9][ 2] =  4./3.;
    colour_[ 3][ 8] = colour_[ 8][ 3] =  4./3.;
    colour_[ 0][ 8] = colour_[ 8][ 0] = -1./6.;
    colour_[ 1][ 9] = colour_[ 9][ 1] = -1./6.;
    colour_[ 2][10] = colour_[10][ 2] = -1./6.;
    colour_[ 3][11] = colour_[11][ 3] = -1./6.;
    // 2 2 block
    colour_[ 4][ 5] = colour_[ 5][ 4] =  4./3.;
    colour_[ 4][ 6] = colour_[ 6][ 4] =  4./3.;
    colour_[ 4][ 7] = colour_[ 7][ 4] =  0.5  ;
    colour_[ 5][ 6] = colour_[ 6][ 5] =  0.5  ;
    colour_[ 5][ 7] = colour_[ 7][ 5] =  4./3.;
    colour_[ 6][ 7] = colour_[ 7][ 6] =  4./3.;
    // 2 3 and 3 2 blocks
    colour_[ 4][ 8] = colour_[ 8][ 4] = -0.5  ;
    colour_[ 4][ 9] = colour_[ 9][ 4] = -1./6.;
    colour_[ 4][10] = colour_[10][ 4] = -1./6.;
    colour_[ 4][11] = colour_[11][ 4] =  0.5  ;
    colour_[ 5][ 8] = colour_[ 8][ 5] = -1./6.;
    colour_[ 5][ 9] = colour_[ 9][ 5] = -0.5  ;
    colour_[ 5][10] = colour_[10][ 5] =  0.5  ;
    colour_[ 5][11] = colour_[11][ 5] = -1./6.;
    colour_[ 6][ 8] = colour_[ 8][ 6] = -1./6.;
    colour_[ 6][ 9] = colour_[ 9][ 6] =  0.5  ;
    colour_[ 6][10] = colour_[10][ 6] = -0.5  ;
    colour_[ 6][11] = colour_[11][ 6] = -1./6.;
    colour_[ 7][ 8] = colour_[ 8][ 7] =  0.5  ;
    colour_[ 7][ 9] = colour_[ 9][ 7] = -1./6.;
    colour_[ 7][10] = colour_[10][ 7] = -1./6.;
    colour_[ 7][11] = colour_[11][ 7] = -0.5  ;
    // 3 3 block
    colour_[ 8][ 9] = colour_[ 9][ 8] =  4./3.;
    colour_[ 8][10] = colour_[10][ 8] =  4./3.;
    colour_[ 8][11] = colour_[11][ 8] =  0.5  ;
    colour_[ 9][10] = colour_[10][ 9] =  0.5  ;
    colour_[ 9][11] = colour_[11][ 9] =  4./3.;
    colour_[10][11] = colour_[11][10] =  4./3.;
    numberOfFlows_ = 12;
    break;
  case Colour33to61: case Colour33to16:
  case Colour3bar3barto6bar1: case Colour3bar3barto16bar:
    colour_ = vector<DVector>(2, DVector(2, 0.));
    colour_[1][1] = colour_[0][0] = 9./4.;
    colour_[0][1] = colour_[1][0] = 3./4.;
    numberOfFlows_ = 2;
    break;
  case Colour38to3bar6: case Colour38to63bar:
    colour_ = vector<DVector>(8, DVector(8, 0.));
    // diagonals
    for(unsigned int ix=0;ix<8;++ix) colour_[ix][ix] = 3.;
    colour_[0][1] = colour_[1][0] = 1.;
    colour_[0][2] = colour_[2][0] = 1.;
    colour_[0][3] = colour_[3][0] = 0.;
    colour_[0][4] = colour_[4][0] = 1.;
    colour_[0][5] = colour_[5][0] = 3.;
    colour_[0][6] = colour_[6][0] = 1.;
    colour_[0][7] = colour_[7][0] = 0.;
    // 1
    colour_[1][2] = colour_[2][1] = 0.;
    colour_[1][3] = colour_[3][1] = 1.;
    colour_[1][4] = colour_[4][1] = 3.;
    colour_[1][5] = colour_[5][1] = 1.;
    colour_[1][6] = colour_[6][1] = 0.;
    colour_[1][7] = colour_[7][1] = 1.;
    // 2
    colour_[2][3] = colour_[3][2] = 1.;
    colour_[2][4] = colour_[4][2] = 0.;
    colour_[2][5] = colour_[5][2] = 1.;
    colour_[2][6] = colour_[6][2] = 3.;
    colour_[2][7] = colour_[7][2] = 1.;
    // 3
    colour_[3][4] = colour_[4][3] = 1.;
    colour_[3][5] = colour_[5][3] = 0.;
    colour_[3][6] = colour_[6][3] = 1.;
    colour_[3][7] = colour_[7][3] = 3.;
    // 4
    colour_[4][5] = colour_[5][4] = 1.;
    colour_[4][6] = colour_[6][4] = 0.;
    colour_[4][7] = colour_[7][4] = 1.;
    // 5
    colour_[5][6] = colour_[6][5] = 1.;
    colour_[5][7] = colour_[7][5] = 0.;
    //6
    colour_[6][7] = colour_[7][6] = 1.;
    numberOfFlows_ = 8;
    break;
  case Colour33to13bar: case Colour33to3bar1:
  case Colour3bar3barto13: case Colour3bar3barto31:
    colour_ = vector<DVector>(1, DVector(1, 6.));
    numberOfFlows_ = 1;
    break;
  case Colour33to83bar: case Colour33to3bar8:
  case Colour3bar3barto83: case Colour3bar3barto38:
  case Colour38to3bar3bar: case Colour3bar8to33:
    colour_ = vector<DVector>(3, DVector(3, 0.));
    colour_[0][0]=colour_[1][1]=colour_[2][2]=8.;
    colour_[0][1]=colour_[1][0]=-4.;
    colour_[0][2]=colour_[2][0]=-4.;
    colour_[1][2]=colour_[2][1]=-4.;
    numberOfFlows_ = 3;
    break;
  default:
    assert(false);
  }
  if(Debug::level > 1 )  {
    generator()->log() << "Set up 2->2 ME for " 
		       << getParticleData(incoming_.first )->PDGName() << " " 
		       << getParticleData(incoming_.second)->PDGName() << " -> "
		       << getParticleData(outgoing_.first )->PDGName() << " " 
		       << getParticleData(outgoing_.second)->PDGName() << "\n";
    for(unsigned int ix=0;ix<alldiagrams.size();++ix) {
      generator()->log() << "Diagram " << ix << " has "
			 << alldiagrams[ix].incoming.first  << " " 
			 << alldiagrams[ix].incoming.second << " -> "
			 << alldiagrams[ix].outgoing.first  << " " 
			 << alldiagrams[ix].outgoing.second << "\n";
      generator()->log() << "Type " << alldiagrams[ix].channelType << " " 
			 << alldiagrams[ix].ordered.first  << " " 
			 << alldiagrams[ix].ordered.second << "\n";
      if(alldiagrams[ix].intermediate)
	generator()->log() << "Intermediate " << alldiagrams[ix].intermediate->PDGName() << "\n";
      generator()->log() << "vertices " 
			 << alldiagrams[ix].vertices.first ->fullName() << " " 
			 << alldiagrams[ix].vertices.second->fullName() << "\n";
    }
  }
}

void GeneralHardME::getDiagrams() const {
  //get ParticleData pointers for external particles
  tcPDPtr ina = getParticleData(getIncoming().first);
  tcPDPtr inb = getParticleData(getIncoming().second);
  tcPDPtr outa = getParticleData(getOutgoing().first);
  tcPDPtr outb = getParticleData(getOutgoing().second);

  bool hasThreePoints = false;
  for(HPCount idx = 0; idx < numberOfDiagrams_; ++idx) {
    const HPDiagram & current = getProcessInfo()[idx];
    tcPDPtr offshell = current.intermediate;
    if(!offshell) continue;
    //t-channel
    if(current.channelType == HPDiagram::tChannel) {
      hasThreePoints = true;
      if(offshell->id() < 0) offshell = offshell->CC();
      if(current.ordered.second)
	add(new_ptr((Tree2toNDiagram(3), ina, offshell,
		     inb, 1, outa, 2, outb, -(idx+1))));
      else 
	add(new_ptr((Tree2toNDiagram(3), ina, offshell,
		     inb, 2, outa, 1, outb, -(idx+1))));
    }
    //s-channel
    else if(current.channelType == HPDiagram::sChannel) {
      hasThreePoints = true;
      add(new_ptr((Tree2toNDiagram(2), ina, inb, 1, offshell,
		   3, outa, 3, outb, -(idx+1))));
    }
    else
      throw MEException() << "getDiagrams() - Unknown diagram in matrix element "
			  << fullName() << Exception::runerror;			  
  }
  if(!hasThreePoints) {
      add(new_ptr((Tree2toNDiagram(2), ina, inb, 1, outa,
		   3, outa, 3, outb, -1)));
  }
}

unsigned int GeneralHardME::orderInAlphaS() const {
  unsigned int order(0);
  for(HPCount idx = 0; idx < numberOfDiagrams_; ++idx) {
    unsigned int tOrder = diagrams_[idx].vertices.first->orderInGs() + 
      diagrams_[idx].vertices.second->orderInGs();
    if(tOrder > order) order = tOrder;
  }
  return order;
}

unsigned int GeneralHardME::orderInAlphaEW() const {
  unsigned int order(0);
  for(HPCount idx = 0; idx < numberOfDiagrams_; ++idx) {
    unsigned int tOrder = diagrams_[idx].vertices.first->orderInGem() + 
      diagrams_[idx].vertices.second->orderInGem();
    if(tOrder > order) order = tOrder;
  }
  return order;
}

Selector<MEBase::DiagramIndex>
GeneralHardME::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if(abs(diags[i]->id()) == int(diagram_+1)) sel.insert(1., i);
  }
  return sel;
}

void GeneralHardME::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << diagrams_ << colour_ << oenum(colourStructure_)
     << numberOfDiagrams_ << numberOfFlows_ << debug_ 
     << scaleChoice_ << scaleFactor_;
}

void GeneralHardME::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> diagrams_ >> colour_ >> ienum(colourStructure_)
     >> numberOfDiagrams_ >> numberOfFlows_ >> debug_ 
     >> scaleChoice_ >> scaleFactor_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<GeneralHardME,HwMEBase>
describeHerwigGeneralHardME("Herwig::GeneralHardME", "Herwig.so");

void GeneralHardME::Init() {

  static ClassDocumentation<GeneralHardME> documentation
    ("This class is designed to be a base class for a specific spin "
     "configuration where no matrix element exists, i.e. when processes "
     "are created automaticlly for a different model.");

}

Selector<const ColourLines *>
GeneralHardME::colourGeometries(tcDiagPtr diag) const {
  // get the current diagram
  const HPDiagram & current = getProcessInfo()[abs(diag->id()) - 1];
  Selector<const ColourLines *> sel;
  switch(colourStructure_) {
  case Colour11to11:
    static ColourLines f11to11("");
    sel.insert(1.,&f11to11);
    break;
  case Colour11to33bar:
    static ColourLines f11to33bar[2]={ColourLines("4 -5"),
				      ColourLines("4 2 -5")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,&f11to33bar[1]);
    else
      sel.insert(1.,&f11to33bar[0]);
    break;
  case Colour11to88:
    static ColourLines f11to88[2]={ColourLines("4 -5, 5 -4"),
				   ColourLines("4 2 -5,5 -2 4")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,&f11to88[1]);
    else
      sel.insert(1.,&f11to88[0]);
    break;
  case Colour33to33:
    static ColourLines f33to33[4]={ColourLines("1 2 5, 3 -2 4"),
				   ColourLines("1 2 4, 3 -2 5"),
				   ColourLines("1 4, 3 5"),
				   ColourLines("1 5, 3 4")};
    static ColourLines f33to33s[4]={ColourLines("1 3:1 4, 2 3:2 5"),
				    ColourLines("1 3:2 4, 2 3:1 5"),
				    ColourLines("1 3:1 5, 2 3:2 4"),
				    ColourLines("1 3:2 5, 2 3:1 4")};
    static ColourLines f33to33t[2]={ColourLines("1 4, 2 5"),
				    ColourLines("1 5, 2 4")};
    if(current.intermediate->iColour() == PDT::Colour8)
      sel.insert(1.,current.ordered.second ? &f33to33[0] : &f33to33[1]);
    else if(current.intermediate->iColour() == PDT::Colour6) {
      sel.insert(1., &f33to33s[2*(flow_-2)+UseRandom::irnd(0,2)]);
    }
    else if(current.intermediate->iColour() == PDT::Colour3bar) {
      sel.insert(1., &f33to33t[flow_-2]);
    }
    else
      sel.insert(1.,current.ordered.second ? &f33to33[2] : &f33to33[3]);
    break;
  case Colour3bar3barto3bar3bar:
    static ColourLines 
      f3bar3barto3bar3bar[4]={ColourLines("-1 -2 -5, -3 2 -4"),
			      ColourLines("-1 -2 -4, -3 2 -5"),
			      ColourLines("-1 -4, -3 -5"),
			      ColourLines("-1 -5, -3 -4")};
    static ColourLines f3bar3barto3bar3bars[4]=
      {ColourLines("-1 -3:1 -4, -2 -3:2 -5"),
       ColourLines("-1 -3:2 -4, -2 -3:1 -5"),
       ColourLines("-1 -3:1 -5, -2 -3:2 -4"),
       ColourLines("-1 -3:2 -5, -2 -3:1 -4")};
    static ColourLines f3bar3barto3bar3bart[2]=
      {ColourLines("-1 -4, -2 -5"),
       ColourLines("-1 -5, -2 -4")};
    if(current.intermediate->iColour() == PDT::Colour8) {
      sel.insert(1.,current.ordered.second ? 
		 &f3bar3barto3bar3bar[0] : &f3bar3barto3bar3bar[1]);
    }
    else if(current.intermediate->iColour() == PDT::Colour6bar) {
      sel.insert(1., &f3bar3barto3bar3bars[2*(flow_-2)+UseRandom::irnd(0,2)]);
    }
    else if(current.intermediate->iColour() == PDT::Colour3) {
      sel.insert(1., &f3bar3barto3bar3bart[flow_-2]);
    }
    else
      sel.insert(1.,current.ordered.second ? 
		 &f3bar3barto3bar3bar[2] : &f3bar3barto3bar3bar[3]);
    break;
  case Colour33barto33bar:
    static ColourLines 
      f33barto33bar[5]={ColourLines("1 2 -3, 4 -2 -5"),
			ColourLines("1 3 4, -2 -3 -5"),
			ColourLines("1 4, -3 -5"),
			ColourLines("1 -2, 4 -5"),
			ColourLines("1 -3, 4 -5")};
    if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second) {
	sel.insert(1.,current.intermediate->iColour() == PDT::Colour8 ? 
		   &f33barto33bar[0] : &f33barto33bar[2]);
      }
      else {
	sel.insert(1.,&f33barto33bar[2*(flow_-2)+2]);
      }
    }
    else
      sel.insert(1.,current.intermediate->iColour() == PDT::Colour8 ? 
		 &f33barto33bar[1] : &f33barto33bar[3]);
    break;
  case Colour33barto11:
    static ColourLines f33barto11[2]={ColourLines("1 -2"),
				      ColourLines("1 2 -3")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,&f33barto11[1]);
    else
      sel.insert(1.,&f33barto11[0]);
    break;
  case Colour33barto88:
    static ColourLines f33barto88[5]={ColourLines("1 4, -4 2 5, -5 -3"),
				      ColourLines("1 5, -5 2 4, -4 -3"),
				      ColourLines("1 3 4, -5 -3 -2, -4 5"),
				      ColourLines("1 3 5, -4 -3 -2, -5 4"),
				      ColourLines("1 -2,4 -5, 5 -4")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,current.ordered.second ? &f33barto88[0] : &f33barto88[1]);
    else if(current.intermediate->iColour() == PDT::Colour8)
      sel.insert(1.,&f33barto88[flow_+2]);
    else
      sel.insert(1.,&f33barto88[4]);
    break;
  case Colour33barto18:
    static ColourLines f33barto18[3]={ColourLines("1 2 5, -3 -5"),
				      ColourLines("1 5, -5 2 -3"),
				      ColourLines("1 3 5,-2 -3 -5")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,current.ordered.second ? &f33barto18[0] : &f33barto18[1]);
    else
      sel.insert(1.,&f33barto18[2]);
    break;
  case Colour33barto81:
    static ColourLines f33barto81[3]={ColourLines("1 4, -4 2 -3"),
				      ColourLines("-3 -4, 1 2 4"),
				      ColourLines("1 3 4,-2 -3 -4")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,current.ordered.second ? &f33barto81[0] : &f33barto81[1]);
    else
      sel.insert(1.,&f33barto81[2]); 
    break;
  case Colour88to11:
    static ColourLines f88barto11[2]={ColourLines("1 -2, 2 -1"),
				      ColourLines("1 -2 -3, 3 2 -1")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,&f88barto11[1]);
    else
      sel.insert(1.,&f88barto11[0]);
    break;
  case Colour88to33bar:
    static ColourLines f88to33bar[5]={ColourLines("1 4, -3 -5, 3 2 -1"),
				      ColourLines("-1 -5, 1 2 -3, 3 4"),
				      ColourLines("2 -1, 1 3 4, -2 -3 -5"),
				      ColourLines("1 -2, -1 -3 -5, 2 3 4"),
				      ColourLines("1 -2, 2 -1, 4 -5")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,current.ordered.second ? &f88to33bar[0] : &f88to33bar[1]);
    else if(current.intermediate->iColour() == PDT::Colour8 )
      sel.insert(1.,&f88to33bar[flow_+2]);
    else
      sel.insert(1.,&f88to33bar[4]);
    break;
  case Colour88to88:
    static ColourLines f88to88s[9]={ColourLines("-1 2, 1 3 5, -5 4, -2 -3 -4"),
				    ColourLines("-1 2, 1 3 4, -4 5, -2 -3 -5"),
				    ColourLines(""),
				    ColourLines("1 -2, -1 -3 -4, 4 -5, 2 3 5"),
				    ColourLines(""),
				    ColourLines("1 -2, -1 -3 -5, 5 -4, 2 3 4"),
				    ColourLines("1 -2, 2 -1, 4 -5, 5 -4"),
				    ColourLines(""),
				    ColourLines("")};
    static ColourLines f88to88t[9]={ColourLines(""),
				    ColourLines("1 4, -1 -2 3, -2 -5, -3 2 5"),
				    ColourLines("-1 -4, 1 2 5, -3 -5, 3 -2 4"),
				    ColourLines("-1 -4, 1 2 -3, 3 5, 4 -2 -5"),
				    ColourLines("1 4, -1 -2 -5, 3 5, -3 2 -4"),
				    ColourLines(""),
				    ColourLines(""),
				    ColourLines("1 4, -1 -4, 3 5, -5 -3"),
				    ColourLines("")};
    static ColourLines f88to88u[9]={ColourLines("1 4, -1 -2 3, -3 -5, -4 2 5"),
				    ColourLines(""),
				    ColourLines("1 5, -1 -2 -4, 3 4, -3 2 -5"),
				    ColourLines(""),
				    ColourLines("-1 -5, 1 2 4, -3 -4, 3 -2 5"),
				    ColourLines("-1 -5, 1 2 -3, 3 4, 5 -2 -4"),
				    ColourLines(""),
				    ColourLines(""),
				    ColourLines("1 5, -1 -5, 3 4, -3 -4")};
    if(current.channelType == HPDiagram::sChannel) {
      sel.insert(1., &f88to88s[flow_]);
    }
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second) {
	sel.insert(1., &f88to88t[flow_]);
      }
      else {
	sel.insert(1., &f88to88u[flow_]);
      }
    }
    break;
  case Colour38to13:
    static ColourLines f38to13[2]={ColourLines("1 2 -3, 3 5"),
				   ColourLines("1 -2, 2 3 5")};
    if(current.channelType == HPDiagram::tChannel) {
      sel.insert(1.,&f38to13[0]);
    }
    else {
      sel.insert(1.,&f38to13[1]);
    }
    break;
  case Colour38to31:
    static ColourLines f38to31[2]={ColourLines("1 2 -3, 3 4"),
				   ColourLines("1 -2, 2 3 4 ")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,&f38to31[0]);
    else	        
      sel.insert(1.,&f38to31[1]);
    break;
  case Colour3bar8to13bar:
    static ColourLines f3bar8to13bar[2]={ColourLines("-1 2 3, -3 -5 "),
					 ColourLines("-1 2, -5 -3 -2")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,&f3bar8to13bar[0]);
    else	        
      sel.insert(1.,&f3bar8to13bar[1]);
    break;
  case Colour3bar8to3bar1:
    static ColourLines f3bar8to3bar1[2]={ColourLines("-1 2 3, -3 -4 "),
					 ColourLines("-1 2, -4 -3 -2")};
    if(current.channelType == HPDiagram::tChannel)
      sel.insert(1.,&f3bar8to3bar1[0]);
    else	        
      sel.insert(1.,&f3bar8to3bar1[1]);
    break;
  case Colour38to83:
    static ColourLines f38to83[5]={ColourLines("1 4, -4 2 -3, 3 5"),
				   ColourLines("1 -2, 2 3 4, 5 -4"),
				   ColourLines("1 2 4, -4  -3, 3 -2 5"),
				   ColourLines("1 2 -3, -4 -2  5, 3 4"),
				   ColourLines("1 5, 3 4, -3 -4")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1.,&f38to83[1]);
    else {
      if(current.intermediate->iColour() == PDT::Colour8)
	sel.insert(1.,&f38to83[flow_+2]);
      else if(current.intermediate->iColour() == PDT::Colour0)
	sel.insert(1.,&f38to83[4]);
      else
	sel.insert(1.,&f38to83[0]);
    }
    break;
  case Colour38to38:
    static ColourLines f38to38[5]={ColourLines("1 5, -5 2 -3, 3 4"),
				   ColourLines("1 -2, 2 3 5, 4 -5"),
				   ColourLines("1 2 5, -5  -3, 3 -2 4"),
				   ColourLines("1 2 -3, -5 -2  4, 3 5"),
				   ColourLines("1 4, 3 5, -3 -5")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1.,&f38to38[1]);
    else {
      if(current.intermediate->iColour() == PDT::Colour8)
	sel.insert(1.,&f38to38[flow_+2]);
      else if(current.intermediate->iColour() == PDT::Colour0)
	sel.insert(1.,&f38to38[4]);
      else
	sel.insert(1.,&f38to38[0]);
    }
    break;
  case Colour3bar8to83bar:
    static ColourLines f3bar8to83bar[5]={ColourLines("-1 -4, 3 2 4, -3 -5"),
					 ColourLines("-1 2, -4 -3 -2, 4 -5"),
					 ColourLines("-1 -2 -4,-3 2 -5,3 4"),
					 ColourLines("-1 -2 3, -5 2 4, -3 -4"),
					 ColourLines("-1 -5, -3 -4, 3 4")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1.,&f3bar8to83bar[1]);
    else {
      if(current.intermediate->iColour() == PDT::Colour8)
	sel.insert(1.,&f3bar8to83bar[flow_+2]);
      else if(current.intermediate->iColour() == PDT::Colour0)
	sel.insert(1.,&f3bar8to83bar[4]);
      else
	sel.insert(1.,&f3bar8to83bar[0]);
    }
    break;
  case Colour3bar8to3bar8:
    static ColourLines f3bar8to3bar8[4]={ColourLines("-1 -5, 3 2 5, -3 -4"),
					 ColourLines("-1 2, -5 -3 -2, 5 -4"),
					 ColourLines("-1 -2 -5,-3 2 -4,3 5"),
					 ColourLines("-1 -2 3, -4 2 5, -3 -5")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1.,&f3bar8to3bar8[1]);
    else {
      if(current.intermediate->iColour() == PDT::Colour8)
	sel.insert(1.,&f3bar8to3bar8[flow_+2]);
      else
	sel.insert(1.,&f3bar8to3bar8[0]);
    }
    break;
  case Colour88to18:
    static ColourLines f88to18[6]={ColourLines("  1  3  5, -1  2, -2 -3 -5"),
				   ColourLines(" -1 -3 -5,  1 -2,  2  3  5"),
				   ColourLines(" 1  2 -3, -1 -2 -5,  3  5"),
				   ColourLines("-1 -2  3,  1  2  5, -3 -5"),
				   ColourLines(" 1  5, -1  2  3, -3 -2 -5"),
				   ColourLines("-1 -5,  1 -2 -3,  3  2  5")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1.,&f88to18[UseRandom::irnd(0,2)]);
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second)
	sel.insert(1.,&f88to18[UseRandom::irnd(2,4)]);
      else 
	sel.insert(1.,&f88to18[UseRandom::irnd(4,6)]);
    }
    break;
  case Colour88to81:
    static ColourLines f88to81[6]={ColourLines("  1  3  4, -1  2, -2 -3 -4"),
				   ColourLines(" -1 -3 -4,  1 -2,  2  3  4"),
				   ColourLines("  1  4, -1  2  3, -3 -2 -4"),
				   ColourLines(" -1 -4,  1 -2 -3,  3  2  4"),
				   ColourLines("  1  2 -3, -1 -2 -4,  3  4"),
				   ColourLines(" -1 -2  3,  1  2  4, -3 -4")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1.,&f88to81[UseRandom::irnd(0,2)]);
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second)
	sel.insert(1.,&f88to81[UseRandom::irnd(2,4)]);
      else 
	sel.insert(1.,&f88to81[UseRandom::irnd(4,6)]);
    }
    break;
  case Colour33barto66bar:
    static ColourLines f33barto66bars[10]
      ={ColourLines("1 3 4:1, -2 -3 -5:1, 4:2 -5:2"),
	ColourLines("1 3 4:1, -2 -3 -5:2, 4:2 -5:1"),
	ColourLines("1 3 4:2, -2 -3 -5:1, 4:1 -5:2"),
	ColourLines("1 3 4:2, -2 -3 -5:2, 4:1 -5:1"),
	ColourLines(""), ColourLines(""),
	ColourLines(""), ColourLines(""),
	ColourLines("1 -2, 4:1 -5:1, 4:2 -5:2"),
	ColourLines("1 -2, 4:1 -5:2, 4:2 -5:1")};
    static ColourLines f33barto66bart[10]
      ={ColourLines(""), ColourLines(""),
	ColourLines(""), ColourLines(""),
	ColourLines("1 4:1, 4:2 2 -5:2, -3 -5:1"),
	ColourLines("1 4:1, 4:2 2 -5:1, -3 -5:2"),
	ColourLines("1 4:2, 4:1 2 -5:2, -3 -5:1"),
	ColourLines("1 4:2, 4:1 2 -5:1, -3 -5:2"),
	ColourLines(""), ColourLines("")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1.,&f33barto66bars[flow_]);
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second)
	sel.insert(1.,&f33barto66bart[flow_]);
      else
	assert(false);
    }
    break;
  case Colour33barto6bar6:
    static ColourLines f33barto6bar6s[8]
      ={ColourLines("1 3 5:1, -2 -3 -4:1, -4:2 5:2"),
	ColourLines("1 3 5:1, -2 -3 -4:2, -4:1 5:2"),
	ColourLines("1 3 5:2, -2 -3 -4:1, -4:2 5:1"),
	ColourLines("1 3 5:2, -2 -3 -4:2, -4:1 5:1"),
	ColourLines(""), ColourLines(""),
	ColourLines(""), ColourLines("")};
    static ColourLines f33barto6bar6u[8]
      ={ColourLines(""), ColourLines(""),
	ColourLines(""), ColourLines(""),
	ColourLines("1 5:1, 5:2 2 -4:1, -3 -4:2"),
    	ColourLines("1 5:1, 5:2 2 -4:2, -3 -4:1"),
    	ColourLines("1 5:2, 5:1 2 -4:1, -3 -4:2"),
    	ColourLines("1 5:2, 5:1 2 -4:2, -3 -4:1")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1.,&f33barto6bar6s[flow_]);
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second)
	assert(false);
      else
	sel.insert(1.,&f33barto6bar6u[flow_]);
    }
   break;
  case Colour88to66bar:
    static ColourLines f88to66t[12]=
      {ColourLines("1 4:2, 4:1 2:1 3, -1 2:2 -5:2, -3 -5:1"),
       ColourLines("1 4:1, 4:2 2:2 3, -1 2:1 -5:2, -3 -5:1"),
       ColourLines("1 4:2, 4:1 2:1 3, -1 2:2 -5:1, -3 -5:2"),
       ColourLines("1 4:1, 4:2 2:2 3, -1 2:1 -5:1, -3 -5:2"),
       ColourLines("1 4:2, 4:1 2:2 -5:2, -1 2:1 3, -3 -5:1"),
       ColourLines("1 4:1, 4:2 2:2 -5:2, -1 2:1 3, -3 -5:1"),
       ColourLines("1 4:2, 4:1 2:2 -5:1, -1 2:1 3, -3 -5:2"),
       ColourLines("1 4:1, 4:2 2:2 -5:1, -1 2:1 3, -3 -5:2"),
       ColourLines(""),ColourLines(""),
       ColourLines(""),ColourLines("")};
    static ColourLines f88to66u[12]=
      {ColourLines("1 2:2 4:2, 4:1 3, -1 -5:2, -3 2:1 -5:1"),
       ColourLines("1 2:1 4:1, 4:2 3, -1 -5:2, -3 2:2 -5:1"),
       ColourLines("1 2:2 4:2, 4:1 3, -1 -5:1, -3 2:1 -5:2"),
       ColourLines("1 2:1 4:1, 4:2 3, -1 -5:1, -3 2:2 -5:2"),
       ColourLines(""),ColourLines(""),
       ColourLines(""),ColourLines(""),
       ColourLines("-1 -5:1, -5:2 2:2 4:1, 1 2:1 -3, 3 4:2"),
       ColourLines("-1 -5:1, -5:2 2:2 4:2, 1 2:1 -3, 3 4:1"),
       ColourLines("-1 -5:2, -5:1 2:2 4:1, 1 2:1 -3, 3 4:2"),
       ColourLines("-1 -5:2, -5:1 2:2 4:2, 1 2:1 -3, 3 4:1")};
    static ColourLines f88to66s[12]=
      {ColourLines(""),ColourLines(""),
       ColourLines(""),ColourLines(""),
       ColourLines("1 3 4:2, 4:1 -5:2, -1 2, -2 -3 -5:1"),
       ColourLines("1 3 4:1, 4:2 -5:2, -1 2, -2 -3 -5:1"),
       ColourLines("1 3 4:2, 4:1 -5:1, -1 2, -2 -3 -5:2"),
       ColourLines("1 3 4:1, 4:2 -5:1, -1 2, -2 -3 -5:2"),
       ColourLines("-1 -3 -5:1, -5:2 4:1, 1 -2, 2 3 4:2"),
       ColourLines("-1 -3 -5:1, -5:2 4:2, 1 -2, 2 3 4:1"),
       ColourLines("-1 -3 -5:2, -5:1 4:1, 1 -2, 2 3 4:2"),
       ColourLines("-1 -3 -5:2, -5:1 4:2, 1 -2, 2 3 4:1")};
    if(current.channelType == HPDiagram::sChannel) 
      sel.insert(1., &f88to66s[flow_]);
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second) {
        sel.insert(1., &f88to66t[flow_]);
      }
      else {
        sel.insert(1., &f88to66u[flow_]);
      }
    }
    break;
  case Colour33to61:
    static ColourLines f33to61t[2]
      ={ColourLines("1 4:1, 3 2 4:2"),
        ColourLines("1 4:2, 3 2 4:1")};
    static ColourLines f33to61u[2]
      ={ColourLines("1 2 4:1, 3 4:2"),
        ColourLines("1 2 4:2, 3 4:1")};
    static ColourLines f33to61s[2]
      ={ColourLines("1 3:1 4:1, 2 3:2 4:2"),
        ColourLines("1 3:2 4:2, 2 3:1 4:1")};
    if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second) 
        sel.insert(1., &f33to61t[flow_]);
      else 
        sel.insert(1., &f33to61u[flow_]);
    }
    else {
      sel.insert(1., &f33to61s[flow_]);
    }
    break;
  case Colour33to16:
    static ColourLines f33to16t[2]
      ={ColourLines("1 2 5:1, 3 5:2"),
        ColourLines("1 2 5:2, 3 5:1")};
    static ColourLines f33to16u[2]
      ={ColourLines("1 5:1, 3 2 5:2"),
        ColourLines("1 5:2, 3 2 5:1")};
    static ColourLines f33to16s[2]
      ={ColourLines("1 3:1 5:1, 2 3:2 5:2"),
        ColourLines("1 3:2 5:2, 2 3:1 5:1")};
    if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second)
        sel.insert(1., &f33to16t[flow_]);
      else 
        sel.insert(1., &f33to16u[flow_]);
    }
    else {
      sel.insert(1., &f33to16s[flow_]);
    }
    break;
  case Colour3bar3barto6bar1:
    static ColourLines f3bar3barto6bar1t[2]
      ={ColourLines("-1 -4:1, -3 -2 -4:2"),
	ColourLines("-1 -4:2, -3 -2 -4:1")};
    static ColourLines f3bar3barto6bar1u[2]
      ={ColourLines("-1 -2 -4:1, -3 -4:2"),
        ColourLines("-1 -2 -4:2, -3 -4:1")};
    static ColourLines f3bar3barto6bar1s[2]
      ={ColourLines("-1 -3:1 -4:1, -2 -3:2 -4:2"),
        ColourLines("-1 -3:2 -4:2, -2 -3:1 -4:1")};
    if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second) 
        sel.insert(1., &f3bar3barto6bar1t[flow_]);
      else 
        sel.insert(1., &f3bar3barto6bar1u[flow_]);
    }
    else {
      sel.insert(1., &f3bar3barto6bar1s[flow_]);
    }
    break;
  case Colour3bar3barto16bar:
    static ColourLines f3bar3barto16bart[2]
      ={ColourLines("-1 -2 -5:1, -3 -5:2"),
	ColourLines("-1 -2 -5:2, -3 -5:1")};
    static ColourLines f3bar3barto16baru[2]
      ={ColourLines("-1 -5:1, -3 -2 -5:2"),
        ColourLines("-1 -5:2, -3 -2 -5:1")};
    static ColourLines f3bar3barto16bars[2]
      ={ColourLines("-1 -3:1 -5:1, -2 -3:2 -5:2"),
        ColourLines("-1 -3:2 -5:2, -2 -3:1 -5:1")};
    if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second) 
        sel.insert(1., &f3bar3barto16bart[flow_]);
      else 
        sel.insert(1., &f3bar3barto16baru[flow_]);
    }
    else {
      sel.insert(1., &f3bar3barto16bars[flow_]);
    }
    break;
  case Colour38to3bar6:
    static ColourLines f38to3bar6t[8]
      ={ColourLines("1 2:1 -3, -4 2:2 5:1, 3 5:2"),
        ColourLines("1 2:1 -3, -4 2:2 5:2, 3 5:1"),
        ColourLines("1 2:1 5:1, -4 2:2 -3, 3 5:2"),
        ColourLines("1 2:1 5:2, -4 2:2 -3, 3 5:1"),
	ColourLines(""),ColourLines(""),
	ColourLines(""),ColourLines("")};
    static ColourLines f38to3bar6u[8]
      ={ColourLines(""),ColourLines(""),
      	ColourLines(""),ColourLines(""),
	ColourLines("1 5:1, 3 2 5:2, -4 -3"),
	ColourLines("1 5:2, 3 2 5:1, -4 -3"),
	ColourLines(""),ColourLines("")};
    static ColourLines f38to3bar6s[8]
      ={ColourLines(""),ColourLines(""),
      	ColourLines(""),ColourLines(""),
      	ColourLines(""),ColourLines(""),
	ColourLines("1 -2, 2 3 5:1, -4 5:2"),
        ColourLines("1 -2, 2 3 5:2, -4 5:1")};
    if(current.channelType == HPDiagram::sChannel) 
      sel.insert(1., &f38to3bar6s[flow_]);
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second) 
        sel.insert(1., &f38to3bar6t[flow_]);
      else 
        sel.insert(1., &f38to3bar6u[flow_]);
    }
  break;
  case Colour38to63bar:
    static ColourLines f38to63baru[8]
      ={ColourLines("1 2:1 -3, -5 2:2 4:1, 3 4:2"),
        ColourLines("1 2:1 -3, -5 2:2 4:2, 3 4:1"),
        ColourLines("1 2:1 4:1, -5 2:2 -3, 3 4:2"),
        ColourLines("1 2:1 4:2, -5 2:2 -3, 3 4:1"),
	ColourLines(""),ColourLines(""),
	ColourLines(""),ColourLines("")};
    static ColourLines f38to63bart[8]
      ={ColourLines(""),ColourLines(""),
      	ColourLines(""),ColourLines(""),
	ColourLines("1 4:1, 3 2 4:2, -5 -3"),
        ColourLines("1 4:2, 3 2 4:1, -5 -3"),
	ColourLines(""),ColourLines("")};
    static ColourLines f38to63bars[8]
      ={ColourLines(""),ColourLines(""),
      	ColourLines(""),ColourLines(""),
      	ColourLines(""),ColourLines(""),
	ColourLines("1 -2, 2 3 4:1, -5 4:2"),
        ColourLines("1 -2, 2 3 4:2, -5 4:1")};
    if(current.channelType == HPDiagram::sChannel) 
      sel.insert(1., &f38to63bars[flow_]);
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second) 
        sel.insert(1., &f38to63bart[flow_]);
      else 
        sel.insert(1., &f38to63baru[flow_]);
    }
    break;
  case Colour33to13bar:
    static ColourLines f33to13bar[3]={ColourLines("1 -6, 2 -6, -5 -3 -6"),
				      ColourLines("1 2 -6, 3 -6, -5 -6"),
				      ColourLines("1 -6, 3 2 -6, -5 -6")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1., &f33to13bar[0]);
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second)
        sel.insert(1., &f33to13bar[1]);
      else
        sel.insert(1., &f33to13bar[2]);
    }
    break;
  case Colour33to3bar1:
    static ColourLines f33to3bar1[3]={ColourLines("1 -6, 2 -6, -4 -3 -6"),
				      ColourLines("1 -6, 3 2 -6, -4 -6"),
				      ColourLines("1 2 -6, 3 -6, -4 -6")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1., &f33to3bar1[0]);
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second)
        sel.insert(1., &f33to3bar1[1]);
      else
        sel.insert(1., &f33to3bar1[2]);
    }
    break;
  case Colour3bar3barto13:
    static ColourLines f3bar3barto13[3]={ColourLines("-1 6, -2 6, 5 3 6"),
					 ColourLines("-1 2 6, -3 6, 5 6"),
					 ColourLines("-1 6, -3 2 6, 5 6")};
    if(current.channelType == HPDiagram::sChannel) {
      sel.insert(1., &f3bar3barto13[0]);
    }
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second) {
        sel.insert(1., &f3bar3barto13[1]);
      }
      else {
        sel.insert(1., &f3bar3barto13[2]);
      }
    }
    break;
  case Colour3bar3barto31:
    static ColourLines f3bar3barto31[3]={ColourLines("-1 6, -2 6, 4 3 6"),
					 ColourLines("-1 6, -3 2 6, 4 6"),
					 ColourLines("-1 2 6, -3 6, 4 6")};
    if(current.channelType == HPDiagram::sChannel)
      sel.insert(1., &f3bar3barto31[0]);
    else if(current.channelType == HPDiagram::tChannel) {
      if(current.ordered.second)
        sel.insert(1., &f3bar3barto31[1]);
      else
        sel.insert(1., &f3bar3barto31[2]);
    }
    break;
  case Colour33to83bar:
    static ColourLines f33to83bar[3]={ColourLines("1 -6, 2 -6, -5 4, -4 -3 -6"),
     				      ColourLines("1 4, -4 2 -6, 3 -6, -5 -6"),
     				      ColourLines("1 -6, 3 4, -4 2 -6, -5 -6")};
    sel.insert(1., &f33to83bar[flow_]);
    break;
  case Colour33to3bar8:
    static ColourLines f33to3bar8[3]={ColourLines("1 -6, 2 -6, -4 5, -5 -3 -6"),
     				      ColourLines("1 -6, 3 5, -5 2 -6, -4 -6"),
				      ColourLines("1 5, -5 2 -6, 3 -6, -4 -6")};
    sel.insert(1., &f33to3bar8[flow_]);
    break;
  case Colour3bar3barto83:
    static ColourLines f3bar3barto83[3]={ColourLines("-1 6, -2 6, 5 -4, 4 3 6"),
					 ColourLines("-1 -4, 4 2 6, -3 6, 5 6"),
					 ColourLines("-1 6, -3 -4, 4 2 6, 5 6")};
    sel.insert(1., &f3bar3barto83[flow_]);
    break;
  case Colour3bar3barto38:
    static ColourLines f3bar3barto38[3]={ColourLines("-1 6, -2 6, 4 -5, 5 3 6"),
					 ColourLines("-1 6, -3 -5, 5 -2 6, 4 6"),
					 ColourLines("-1 -5, 5 -2 6, -3 6, 4 6")};
    sel.insert(1., &f3bar3barto38[flow_]);
    break;
  case Colour38to3bar3bar:
    static ColourLines f38to3bar3bar[3]={ColourLines("1 -2, 2 3 -6, -4 -6, -5 -6"),
					 ColourLines("1 -6, 3 2 -6, -3 -5, -4 -6"),
					 ColourLines("1 -6, 3 2 -6, -3 -4, -5 -6")};
    sel.insert(1., &f38to3bar3bar[flow_]);
    break;
  case Colour3bar8to33:
    static ColourLines f3bar8to33[3]={ColourLines("-1 2, -2 -3 6, 4 6, 5 6"),
				      ColourLines("-1 6, -3 2 6, 3 5, 4 6"),
				      ColourLines("-1 6, -3 2 6, 3 4, 5 6")};
    sel.insert(1., &f3bar8to33[flow_]);
    break;
  default:
    assert(false);
  }
  return sel;
}

double GeneralHardME::
selectColourFlow(vector<double> & flow,vector<double> & me,
		 double average) const {
  // spin average
  double output = 0.25*average;
  // special for beam polarization
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = 
      {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
       beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    for(unsigned int ix = 0;ix<numberOfFlows();++ix)
      flow[ix] = flowME_[ix].average(rho[0],rho[1]);
    for(unsigned int ix = 0;ix<numberOfDiags();++ix)
      me  [ix] = diagramME_[ix].average(rho[0],rho[1]);
    output = 0.;
    for(unsigned int ii = 0; ii < numberOfFlows(); ++ii) 
      for(unsigned int ij = 0; ij < numberOfFlows(); ++ij)
	output += real(getColourFactors()[ii][ij]*
		       flowME_[ii].average(flowME_[ij],rho[0],rho[1]));
    // correction for photons and gluons
    if(mePartonData()[0]->id()==ParticleID::g || 
       mePartonData()[0]->id()==ParticleID::gamma) output *= 1.5;
    if(mePartonData()[1]->id()==ParticleID::g ||
       mePartonData()[1]->id()==ParticleID::gamma) output *= 1.5;
  }
  // select the colour flow
  double maxWgt = UseRandom::rnd()*std::accumulate(flow.begin(),flow.end(),0.);
  flow_ = flow.size();
  for(unsigned int ix=0;ix<flow.size();++ix) {
    if(flow[ix]>=maxWgt) {
      flow_=ix;
      break;
    }
    maxWgt -= flow[ix];
  }
  assert(flow_<flow.size());
  // select the diagram
  for(unsigned int ix=0;ix<numberOfDiags();++ix) {
    const HPDiagram & current = getProcessInfo()[ix];
    bool found=false;
    for(unsigned int iy = 0; iy < current.colourFlow.size(); ++iy) {
      if(current.colourFlow[iy].first==flow_) {
	me[ix] *= sqr(current.colourFlow[iy].second);
	found = true;
      }
    }
    // set to zero if four point diagram or doesn't contribute to colour flow
    if(!found || current.channelType == HPDiagram::fourPoint) me[ix]=0.;
  }
  maxWgt = UseRandom::rnd()*std::accumulate(me.begin(),me.end(),0.);
  for(unsigned int ix=0;ix<me.size();++ix) {
    if(me[ix]>maxWgt) {
      diagram_=ix;
      break;
    }
    maxWgt -= me[ix];
  }
  // colour factors
  output /= max(1,abs(int(mePartonData()[0]->iColour())));
  output /= max(1,abs(int(mePartonData()[1]->iColour())));
  // identical particle factor
  output *=  mePartonData()[2]->id() == mePartonData()[3]->id() ? 0.5 : 1;
  // return the answer
  return output;
}

void GeneralHardME::doinitrun() {
  HwMEBase::doinitrun();
  for(unsigned int ix=0;ix<diagrams_.size();++ix) {
    diagrams_[ix].vertices.first ->initrun();
    diagrams_[ix].vertices.second->initrun();
  }
}
