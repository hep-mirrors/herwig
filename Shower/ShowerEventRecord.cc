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

ShowerEventRecord::ShowerEventRecord() {}

ShowerEventRecord::~ShowerEventRecord() {}

void ShowerEventRecord::clear() {
  subProcess_ = SubProPtr();
  XComb_ = StdXCombPtr();
  incoming_ = PPair();
  outgoing_.clear();
  intermediates_.clear();
  PDFs_ = pair<PDF,PDF>();
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
