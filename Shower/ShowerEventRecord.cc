// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerEventRecord class.
//

#include "ShowerEventRecord.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Base/SubtractedME.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

using namespace Herwig;

ShowerEventRecord::ShowerEventRecord() 
  : isMCatNLOSEvent_(false),isMCatNLOHEvent_(false),
    isPowhegSEvent_ (false),isPowhegHEvent_ (false)
{}

ShowerEventRecord::~ShowerEventRecord() {}


void ShowerEventRecord::updateColour(PPtr particle,
				     bool recursive) {
  // if attached to a colour line
  if(particle->colourLine()) {
    // one and only one
    if(particle->colourInfo()->colourLines().size()==1) {
      bool reset=false;
      // if colour line from hard process reconnect
      ColinePtr c1=particle->colourLine();
      if(colourLines().find(c1)!=colourLines().end()) {
	c1->removeColoured(particle);
	colourLines()[c1]->addColoured(particle);
	reset=true;
      }
      // ensure properly connected to the line
      if(!reset) {
	ColinePtr c1=particle->colourLine();
	c1->removeColoured(particle);
	c1->addColoured(particle);
      }
    }
    else {
      Ptr<MultiColour>::pointer colour = 
	dynamic_ptr_cast<Ptr<MultiColour>::pointer>(particle->colourInfo());
      vector<tcColinePtr> lines = colour->colourLines();
      for(unsigned int ix=0;ix<lines.size();++ix) {
	ColinePtr c1 = const_ptr_cast<ColinePtr>(lines[ix]);
	if(colourLines().find(c1)!=colourLines().end()) {
	  colour->colourLine(colourLines()[c1],int(ix)+1);
	  c1->removeColoured(particle);
	}
      }
    }
  }
  // if attached to an anticolour line
  if(particle->antiColourLine()) {
    // one and only one
    if(particle->colourInfo()->antiColourLines().size()==1) {
      bool reset=false;
      ColinePtr c1=particle->antiColourLine();
      // if anti colour line from hard process reconnect
      if(colourLines().find(c1)!=colourLines().end()) {
	c1->removeColoured(particle,true);
	colourLines()[c1]->addColoured(particle,true);
	reset=true;
      }
      if(!reset) {
	ColinePtr c1=particle->antiColourLine();
	c1->removeColoured(particle,true);
	c1->addColoured(particle,true);
      }
    }
    else {
      Ptr<MultiColour>::pointer colour = 
	dynamic_ptr_cast<Ptr<MultiColour>::pointer>(particle->colourInfo());
      vector<tcColinePtr> lines = colour->antiColourLines();
      for(unsigned int ix=0;ix<lines.size();++ix) {
	ColinePtr c1 = const_ptr_cast<ColinePtr>(lines[ix]);
	if(colourLines().find(c1)!=colourLines().end()) {
	  colour->antiColourLine(colourLines()[c1],int(ix)+1);
	  c1->removeColoured(particle,true);
	}
      }
    }
  }
  if(!recursive) return;
  for ( ParticleVector::const_iterator c = particle->children().begin();
	c != particle->children().end(); ++c ) {
    updateColour(*c,true);
  }
}

void ShowerEventRecord::colourIsolate(const vector<PPtr> & original,
				      const vector<PPtr> & copy) {
  // vectors must have same size
  assert(original.size()==copy.size());
  // create a temporary map with all the particles to make looping easier
  vector<PPair> particles;
  particles.reserve(original.size());
  for(unsigned int ix=0;ix<original.size();++ix)
    particles.push_back(make_pair(copy[ix],original[ix]));
  // reset the colour of the copies
  vector<PPair>::const_iterator cit;
  // make the colour connections of the copies
  for(cit=particles.begin();cit!=particles.end();++cit) {
    if((*cit).first->colourInfo()) {
      if((*cit).first->dataPtr()->iColour() != PDT::Colour6 &&
	 (*cit).first->dataPtr()->iColour() != PDT::Colour6bar)
	(*cit).first->colourInfo(new_ptr(ColourBase()));
      else 
	(*cit).first->colourInfo(new_ptr(MultiColour()));
    }
  }
  map<tColinePtr,tColinePtr> cmap;
  // make the colour connections of the copies
  // loop over the particles
  for(cit=particles.begin();cit!=particles.end();++cit) {
    // if particle has at least one colour line
    if((*cit).second->colourLine()) {
      // one and only one line
      if(int((*cit).second->colourInfo()->colourLines().size())==1) {
	// if particle has a colour line
        if(!(*cit).first->colourLine()) {
          // make new line
          tColinePtr oldline=(*cit).second->colourLine();
          ColinePtr newline=ColourLine::create((*cit).first);
          cmap[oldline]=newline;
          isolateLine(cit,particles,oldline,newline);
	}
      }
      // more than one line
      else {
        Ptr<MultiColour>::pointer colour1 = 
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cit).second->colourInfo());
        vector<tcColinePtr> lines1 = colour1->colourLines();
        Ptr<MultiColour>::pointer colour2 = 
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cit).first->colourInfo());
        vector<tcColinePtr> lines2 = colour2->colourLines();
        // loop over lines
        for(unsigned int ix=0;ix<lines1.size();++ix) {
          if( (lines2.size()>ix && !lines2[ix]) ||
              lines2.size()<=ix) {
            tColinePtr oldline = const_ptr_cast<tColinePtr>(lines1[ix]);
            ColinePtr newline = new_ptr(ColourLine());
            cmap[oldline]=newline;
            colour2->colourLine(newline, int(ix)+1);
            isolateLine(cit,particles,oldline,newline);
	  }
        }
      }
    }
    // if anticolour line
    if((*cit).second->antiColourLine()) {
      // one and only one line
      if(int((*cit).second->colourInfo()->antiColourLines().size())==1) {
	// if not already change
	if(!(*cit).first->antiColourLine()) {
	  // make new line
	  tColinePtr oldline=(*cit).second->antiColourLine();
	  ColinePtr newline=ColourLine::create((*cit).first, true);
	  cmap[oldline]=newline;
	  isolateLine(cit,particles,oldline,newline);
	}
      }
      // more than one line
      else {
	Ptr<MultiColour>::pointer colour1 = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  ((*cit).second->colourInfo());
	vector<tcColinePtr> lines1 = colour1->antiColourLines();
	Ptr<MultiColour>::pointer colour2 = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  ((*cit).first->colourInfo());
	vector<tcColinePtr> lines2 = colour2->antiColourLines();
	// loop over lines
	for(unsigned int ix=0;ix<lines1.size();++ix) {
	  if( (lines2.size()>ix && !lines2[ix]) ||
	      lines2.size()<=ix) {
	    tColinePtr oldline = const_ptr_cast<tColinePtr>(lines1[ix]);
	    ColinePtr newline = new_ptr(ColourLine());
	    cmap[oldline]=newline;
	    colour2->antiColourLine(newline, int(ix)+1);
	    isolateLine(cit,particles,oldline,newline);
	  }
	}
      }
    }    
  }
  for ( map<tColinePtr,tColinePtr>::const_iterator c = cmap.begin();
	c != cmap.end(); ++c ) {
    colourLines()[c->second] = c->first;
  }
  // sort out sinks and sources
  for(cit=particles.begin();cit!=particles.end();++cit) {
    tColinePtr cline[2];
    tColinePair cpair;
    for(unsigned int ix=0;ix<4;++ix) {
      cline[0] = ix<2 ? cit->second->colourLine() : cit->second->antiColourLine();
      cline[1] = ix<2 ? cit->first ->colourLine() : cit->first ->antiColourLine();
      if(cline[0]) {
	switch (ix) {
	case 0: case 2:
 	  cpair = cline[0]->sinkNeighbours();
	  break;
	case 1: case 3:
	  cpair = cline[0]->sourceNeighbours();
	  break;
	};
      }
      else {
	cpair = make_pair(tColinePtr(),tColinePtr());
      }
      if(cline[0]&&cpair.first) {
 	map<tColinePtr,tColinePtr>::const_iterator 
	  mit[2] = {cmap.find(cpair.first),cmap.find(cpair.second)};
	if(mit[0]!=cmap.end()&&mit[1]!=cmap.end()) {
	  if(ix==0||ix==2) {
	    cline[1]->setSinkNeighbours(mit[0]->second,mit[1]->second);
	  }
	  else {
	    cline[1]->setSourceNeighbours(mit[0]->second,mit[1]->second);
	  }
	}
      }
    }
  }
}

void ShowerEventRecord::isolateLine(vector<PPair>::const_iterator cit,
				    vector<PPair> & particles,
				    tcColinePtr oldline,
				    tColinePtr  newline) {
  // loop over particles
  for(vector<PPair>::const_iterator cjt=particles.begin();
      cjt!=particles.end();++cjt) {
    if(cjt==cit) continue;
    // if particle has colour line
    if((*cjt).second->colourLine()) {
      // if only one check if current line and reset
      if(int((*cjt).second->colourInfo()->colourLines().size())==1) {
	if((*cjt).second->colourLine()==oldline)
	  newline->addColoured((*cjt).first);
      }
      // if more than one check if each line current line and reset 
      else {
        Ptr<MultiColour>::pointer colour1 =
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cjt).second->colourInfo());
        Ptr<MultiColour>::pointer colour2 =
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cjt).first ->colourInfo());
        for(unsigned int ix=0;ix<colour1->colourLines().size();++ix) {
          if(colour1->colourLines()[ix]==oldline)
	    colour2->colourLine(newline,int(ix)+1);
	}
      }
    }  
    // if particle has anticolour line
    if((*cjt).second->antiColourLine()) {
      // if only one check if current line and reset
      if(int((*cjt).second->colourInfo()->antiColourLines().size())==1) {
	if((*cjt).second->antiColourLine()==oldline)
	  newline->addColoured((*cjt).first,true);
      }
      // if more than one check if each line current line and reset 
      else {
        Ptr<MultiColour>::pointer colour1 =
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cjt).second->colourInfo());
        Ptr<MultiColour>::pointer colour2 =
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cjt).first ->colourInfo());
        for(unsigned int ix=0;ix<colour1->antiColourLines().size();++ix) {
          if(colour1->antiColourLines()[ix]==oldline)
	    colour2->antiColourLine(newline, int(ix)+1);
        }
      }
    }
  }
}

void ShowerEventRecord::mapColour(PPtr original,
				  PPtr copy) {
  // has colour line
  if(copy->colourLine()) {
    // one and only one
    if(copy->colourInfo()->colourLines().size()==1) {
      colourLines_[copy->colourLine()] = original->colourLine();
    }
    // more than one
    else {
      Ptr<MultiColour>::pointer colour1 =
        dynamic_ptr_cast<Ptr<MultiColour>::pointer>(copy->colourInfo());
      vector<tcColinePtr> lines1 = colour1->colourLines();
      Ptr<MultiColour>::pointer colour2 =
        dynamic_ptr_cast<Ptr<MultiColour>::pointer>(original->colourInfo());
      vector<tcColinePtr> lines2 = colour2->colourLines();
      for(unsigned int ix=0;ix<lines1.size();++ix)
        colourLines_[const_ptr_cast<ColinePtr>(lines1[ix])] =
	  const_ptr_cast<ColinePtr>(lines2[ix]);
    }
  }
  // has anticolour line
  if(copy->antiColourLine()) {
    // one and only one
    if(copy->colourInfo()->antiColourLines().size()==1) {
      colourLines_[copy->antiColourLine()] = original->antiColourLine();
    }
    // more than one
    else {
      Ptr<MultiColour>::pointer colour1 =
        dynamic_ptr_cast<Ptr<MultiColour>::pointer>(copy->colourInfo());
      vector<tcColinePtr> lines1 = colour1->antiColourLines();
      Ptr<MultiColour>::pointer colour2 =
        dynamic_ptr_cast<Ptr<MultiColour>::pointer>(original->colourInfo());
      vector<tcColinePtr> lines2 = colour2->antiColourLines();
      for(unsigned int ix=0;ix<lines1.size();++ix)
        colourLines_[const_ptr_cast<ColinePtr>(lines1[ix])] =
	  const_ptr_cast<ColinePtr>(lines2[ix]);
    }
  }
}

void ShowerEventRecord::clear() {
  subProcess_ = SubProPtr();
  XComb_ = StdXCombPtr();
  incoming_ = PPair();
  outgoing_.clear();
  intermediates_.clear();
  PDFs_ = pair<PDF,PDF>();
  colourLines_.clear();
}

void ShowerEventRecord::identifyEventType() {
  isMCatNLOSEvent_ = false;
  isMCatNLOHEvent_ = false;
  isPowhegSEvent_  = false;
  isPowhegHEvent_  = false;
  Ptr<SubtractedME>::tptr subme;
  Ptr<MatchboxMEBase>::tptr me;

  Ptr<StandardXComb>::ptr sxc = dynamic_ptr_cast<Ptr<StandardXComb>::ptr>(xcombPtr());
  if ( sxc ) {
    subme = dynamic_ptr_cast<Ptr<SubtractedME>::tptr>(sxc->matrixElement());
    me = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(sxc->matrixElement());
  }
  if ( subme ) {
    if ( subme->showerApproximation() ) {
      showerApproximation_ = subme->showerApproximation();
      // separate MCatNLO and POWHEG-type corrections
      if ( !subme->showerApproximation()->needsSplittingGenerator() ) {
	if ( subme->realShowerSubtraction() )
	  isMCatNLOHEvent_ = true;
	else if ( subme->virtualShowerSubtraction() )
	  isMCatNLOSEvent_ = true;
      }
      else {
	if ( subme->realShowerSubtraction() )
	  isPowhegHEvent_ = true;
	else if ( subme->virtualShowerSubtraction() ||  subme->loopSimSubtraction() )
	  isPowhegSEvent_ = true;
      }
    }
  } else if ( me ) {
    if ( me->factory()->showerApproximation() ) {
      showerApproximation_ = me->factory()->showerApproximation();
      if ( !me->factory()->showerApproximation()->needsSplittingGenerator() ) 
	isMCatNLOSEvent_ = true;
      else
	isPowhegSEvent_ = true;
    }
  }
  // check for truncated shower
  truncatedShower_ = false;
  if (me && me->factory()->showerApproximation()) {
    if(me->factory()->showerApproximation()->needsTruncatedShower())
      truncatedShower_ = true;
  }
  else if (subme && subme->factory()->showerApproximation()) {
    if(subme->factory()->showerApproximation()->needsTruncatedShower())
      truncatedShower_ = true;
  }
}
