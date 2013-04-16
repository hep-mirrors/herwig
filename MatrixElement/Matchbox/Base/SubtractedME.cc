// -*- C++ -*-
//
// SubtractedME.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SubtractedME class.
//

#include "SubtractedME.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Handlers/StdXCombGroup.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/MatchboxXCombGroup.h"
#include "Herwig++/MatrixElement/Matchbox/MatchboxFactory.h"

using namespace Herwig;

SubtractedME::SubtractedME() 
  : MEGroup(), 
    theRealShowerSubtraction(false), theVirtualShowerSubtraction(false) {}

SubtractedME::~SubtractedME() {}

Ptr<MatchboxFactory>::tcptr SubtractedME::factory() const { return theFactory; }

void SubtractedME::factory(Ptr<MatchboxFactory>::tcptr f) { theFactory = f; }

bool SubtractedME::subProcessGroups() const { 
  return 
    (factory()->subProcessGroups() && !(showerApproximation() || inclusive())) ||
    factory()->subtractionData() != "";
}

bool SubtractedME::inclusive() const { return factory()->inclusive(); }

Ptr<ShowerApproximation>::tptr SubtractedME::showerApproximation() const { return factory()->showerApproximation(); }

const vector<Ptr<MatchboxMEBase>::ptr>& SubtractedME::borns() const { return factory()->bornMEs(); }

bool SubtractedME::verbose() const { return factory()->verbose(); }

IBPtr SubtractedME::clone() const {
  return new_ptr(*this);
}

IBPtr SubtractedME::fullclone() const {
  return new_ptr(*this);
}

StdXCombPtr SubtractedME::makeXComb(Energy newMaxEnergy, const cPDPair & inc,
				    tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
				    tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
				    const PBPair & newPartonBins, tCutsPtr newCuts,
				    const DiagramVector & newDiagrams, bool mir,
				    const PartonPairVec& allPBins,
				    tStdXCombPtr newHead,
				    tMEPtr newME) {

  tMEGroupPtr newMEGroup = dynamic_ptr_cast<tMEGroupPtr>(newME);
  if ( !newMEGroup )
    newMEGroup = this;
  Ptr<MatchboxXCombGroup>::ptr res =  
    new_ptr(MatchboxXCombGroup(newMaxEnergy, inc,
			       newEventHandler, newSubProcessHandler,
			       newExtractor, newCKKW,
			       newPartonBins, newCuts, newMEGroup,
			       newDiagrams, mir,
			       newHead));
  res->build(allPBins);

  theReal->prepareXComb(*res);

  if ( factory()->subtractionData() != "" ) {
    set<cPDVector> procs;
    for ( DiagramVector::const_iterator d = head()->diagrams().begin();
	  d != head()->diagrams().end(); ++d ) {
      if ( procs.find((**d).partons()) == procs.end() )
	procs.insert((**d).partons());
    }
    for ( set<cPDVector>::const_iterator p = procs.begin();
	  p != procs.end(); ++p ) {
      for ( size_t i = 0; i < (*p).size(); ++i ) {
	if ( !(*p)[i]->coloured() )
	  continue;
	if ( i > 1 && 
	     (*p)[i]->id() == ParticleID::g ) {
	  softHistograms[SoftSubtractionIndex(*p,i)] = SubtractionHistogram(0.001,10.);
	  if ( theReal->phasespace() )
	    singularLimits().insert(make_pair(i,i));
	}
	for ( size_t j = i+1; j < (*p).size(); ++j ) {
	  if ( !(*p)[j]->coloured() )
	    continue;
	  long iid = (*p)[i]->id();
	  long jid = (*p)[j]->id();
	  if ( i < 2 && j < 2 )
	    continue;
	  if ( i < 2 && j > 1 ) {
	    if ( abs(iid) < 7 && abs(jid) < 7 && iid != jid )
	      continue;
	  }
	  if ( i > 1 && j > 1 ) {
	    if ( abs(iid) < 7 && abs(jid) < 7 && iid + jid != 0 )
	      continue;
	  }
	  bool haveDipole = false;
	  for ( MEVector::const_iterator k = dependent().begin();
		k != dependent().end(); ++k ) {
	    const SubtractionDipole& dip = dynamic_cast<const SubtractionDipole&>(**k);
	    if ( ( (size_t)(dip.realEmitter()) == i && (size_t)(dip.realEmission()) == j ) ||
		 ( (size_t)(dip.realEmitter()) == j && (size_t)(dip.realEmission()) == i ) ) {
	      haveDipole = true;
	      break;
	    }
	  }
	  if ( !haveDipole )
	    continue;
	  collinearHistograms[CollinearSubtractionIndex(*p,make_pair(i,j))] = SubtractionHistogram(0.001,20.);
	  if ( theReal->phasespace() )
	    singularLimits().insert(make_pair(i,j));
	}
      }
    }
  }

  return res;

}

void SubtractedME::setXComb(tStdXCombPtr xc) {
  MEGroup::setXComb(xc);
  lastMatchboxXComb(xc);
}

MEBase::DiagramVector SubtractedME::dependentDiagrams(const cPDVector& proc,
						      tMEPtr depME) const {

  Ptr<SubtractionDipole>::tptr dipole = 
    dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(depME);

  if ( !dipole ) {
    Throw<InitException>() << "A dependent matrix element of SubtractedME "
			   << "has not been derived from SubtractionDipole. "
			   << "Please check the corresponding input file.";
  }

  return dipole->underlyingBornDiagrams(proc);

}

vector<Ptr<SubtractionDipole>::ptr> SubtractedME::dipoles() {

  if ( dependent().empty() )
    getDipoles();
  vector<Ptr<SubtractionDipole>::ptr> res;
  for ( MEVector::const_iterator k = dependent().begin();
	k != dependent().end(); ++k )
    res.push_back(dynamic_ptr_cast<Ptr<SubtractionDipole>::ptr>(*k));
  return res;
}

void SubtractedME::getDipoles() {

  if ( !dependent().empty() )
    return;

  Ptr<MatchboxMEBase>::tptr real =
    dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(head());

  if ( borns().empty() || !real )
    Throw<InitException>() << "The SubtractedME '"
			   << name() << "' could not generate "
			   << "subtraction terms for the real emission "
			   << "matrix element '" << real->name() << "'. "
			   << "Please check the corresponding input file.";

  Ptr<MatchboxMEBase>::ptr myRealEmissionME = real->cloneMe();
  ostringstream pname;
  pname << fullName() << "/" << myRealEmissionME->name();
  if ( ! (generator()->preinitRegister(myRealEmissionME,pname.str()) ) )
    throw InitException() << "Matrix element " << pname.str() << " already existing.";
  myRealEmissionME->cloneDependencies(pname.str());
  head(myRealEmissionME);
  real = myRealEmissionME;

  MEVector dipMEs;
  vector<Ptr<SubtractionDipole>::ptr> genDipoles
    = real->getDipoles(DipoleRepository::dipoles(),borns());

  if ( factory()->subtractionData() != "" ) {
    for ( vector<Ptr<SubtractionDipole>::ptr>::const_iterator d =
	    genDipoles.begin(); d != genDipoles.end(); ++d )
      (**d).doTestSubtraction();
  }

  if ( genDipoles.empty() ) {
    // probably finite real contribution, but warn
    generator()->log() << "\nWarning: No subtraction dipoles could be found for the processes:\n";
    for ( vector<PDVector>::const_iterator s = real->subProcesses().begin();
	  s != real->subProcesses().end(); ++s ) {
      generator()->log() << (*s)[0]->PDGName() << " " << (*s)[1]->PDGName()
			 << " -> ";
      for ( PDVector::const_iterator p = s->begin() + 2; p != s->end(); ++p )
	generator()->log() << (**p).PDGName() << " ";
      generator()->log() << "\n" << flush;
    }
    generator()->log() << "Assuming finite tree-level O(alphaS) correction.\n";
  }

  dipMEs.resize(genDipoles.size());
  copy(genDipoles.begin(),genDipoles.end(),dipMEs.begin());

  dependent() = dipMEs;

}

void SubtractedME::cloneRealME(const string& prefix) {

  theReal = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(head());

  if ( theReal ) {
    Ptr<MatchboxMEBase>::ptr myRealEmissionME = theReal->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myRealEmissionME->name();
    if ( ! (generator()->preinitRegister(myRealEmissionME,pname.str()) ) )
      throw InitException() << "Matrix element " << pname.str() << " already existing.";
    myRealEmissionME->cloneDependencies(pname.str());
    theReal = myRealEmissionME;
  }

  head(theReal);

}

void SubtractedME::cloneDipoles(const string& prefix) {

  MEVector dipMEs;

  for ( MEVector::const_iterator m = dependent().begin();
	m != dependent().end(); ++m ) {
    Ptr<SubtractionDipole>::tptr dip = 
      dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(*m);
    assert(dip);

    Ptr<SubtractionDipole>::ptr cloned = dip->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << cloned->name();
    if ( ! (generator()->preinitRegister(cloned,pname.str()) ) )
      throw InitException() << "Subtraction dipole " << pname.str() << " already existing.";
    cloned->cloneDependencies(pname.str());
    dipMEs.push_back(cloned);

  }

  dependent() = dipMEs;

}


vector<Ptr<SubtractionDipole>::ptr> SubtractedME::splitDipoles(const cPDVector& born) {

  vector<Ptr<SubtractionDipole>::ptr> dips = dipoles();

  vector<Ptr<SubtractionDipole>::ptr> res;

  for ( vector<Ptr<SubtractionDipole>::ptr>::iterator d = dips.begin();
	d != dips.end(); ++d ) {
    for ( DiagramVector::const_iterator p = (**d).underlyingBornME()->diagrams().begin();
	  p != (**d).underlyingBornME()->diagrams().end(); ++p )
      if ( born == (**p).partons() ) {
	res.push_back(*d);
	break;
      }
  }

  return res;

}

void SubtractedME::doRealEmissionScales() { 
  for ( MEVector::const_iterator m = dependent().begin();
	m != dependent().end(); ++m ) {
    Ptr<SubtractionDipole>::tptr dip = 
      dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(*m);
    assert(dip);
    dip->doRealEmissionScales();
  }
}

void SubtractedME::doRealShowerSubtraction() { 
  theRealShowerSubtraction = true;
  for ( MEVector::const_iterator m = dependent().begin();
	m != dependent().end(); ++m ) {
    Ptr<SubtractionDipole>::tptr dip = 
      dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(*m);
    assert(dip);
    dip->showerApproximation(showerApproximation());
    dip->doRealShowerSubtraction();
  }
}

void SubtractedME::doVirtualShowerSubtraction() { 
  theVirtualShowerSubtraction = true; 
  for ( MEVector::const_iterator m = dependent().begin();
	m != dependent().end(); ++m ) {
    Ptr<SubtractionDipole>::tptr dip = 
      dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(*m);
    assert(dip);
    dip->showerApproximation(showerApproximation());
    dip->doVirtualShowerSubtraction();
  }
}

void SubtractedME::setVetoScales(tSubProPtr) const {}

void SubtractedME::fillProjectors() {
  if ( !inclusive() && !virtualShowerSubtraction() )
    return;
  Ptr<StdXCombGroup>::tptr group = 
    dynamic_ptr_cast<Ptr<StdXCombGroup>::tptr>(lastXCombPtr());
  for ( vector<StdXCombPtr>::const_iterator d = group->dependent().begin();
	d != group->dependent().end(); ++d ) {
    if ( !(**d).kinematicsGenerated() )
      continue;
    if ( (**d).willPassCuts() )
      lastXCombPtr()->projectors().insert(1.,*d);
  }
}

double SubtractedME::reweightHead(const vector<tStdXCombPtr>& dep) {

  if ( inclusive() && !lastXComb().lastProjector() )
    return 1.;

  if ( virtualShowerSubtraction() && !lastXComb().lastProjector() ) {
    return 0.;
  }

  if ( realShowerSubtraction() ) {
    assert(showerApproximation());
    bool below = showerApproximation()->belowCutoff();
    bool haveDipole = false;
    for ( vector<tStdXCombPtr>::const_iterator d = dep.begin(); d != dep.end(); ++d ) {
      if ( !(**d).kinematicsGenerated() )
	continue;
      if ( (**d).willPassCuts() ) {
	haveDipole = true;
	break;
      }
    }
    if ( !haveDipole )
      return 1.;
    if ( below )
      return 0.;
    return 1.;
  }

  if ( virtualShowerSubtraction() || inclusive() ) {
    if ( virtualShowerSubtraction() ) {
      assert(showerApproximation());
      bool above = !showerApproximation()->belowCutoff();
      if ( above )
	return 0.;
    }
    double sum = 0.;
    size_t n = 0;
    for ( vector<tStdXCombPtr>::const_iterator d = dep.begin(); d != dep.end(); ++d ) {
      if ( !(**d).kinematicsGenerated() )
	continue;
      if ( (**d).willPassCuts() ) {
	sum += (**d).lastME2();
	++n;
      }
    }
    return
      n * lastXComb().lastProjector()->lastME2() / sum;
  }

  return 1.;

}

double SubtractedME::reweightDependent(tStdXCombPtr xc, const vector<tStdXCombPtr>& dep) {

  if ( inclusive() && !lastXComb().lastProjector() )
    return 0.;

  if ( virtualShowerSubtraction() && !lastXComb().lastProjector() ) {
    return 0.;
  }

  if ( virtualShowerSubtraction() || inclusive() ) {
    if ( xc != lastXComb().lastProjector() )
      return 0.;
    size_t n = 0;
    for ( vector<tStdXCombPtr>::const_iterator d = dep.begin(); d != dep.end(); ++d ) {
      if ( !(**d).kinematicsGenerated() )
	continue;
      if ( (**d).willPassCuts() ) {
	++n;
      }
    }
    return n;
  }

  if ( realShowerSubtraction() ) {
    bool haveDipole = false;
    for ( vector<tStdXCombPtr>::const_iterator d = dep.begin(); d != dep.end(); ++d ) {
      if ( !(**d).kinematicsGenerated() )
	continue;
      if ( (**d).willPassCuts() ) {
	haveDipole = true;
	break;
      }
    }
    if ( !haveDipole )
      return 0.;
  }

  return 1.;

}

void SubtractedME::doinit() {

  // has been deactivated by the factory
  if ( !head() ) {
    MEBase::doinit();
    return;
  }

  theReal = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(head());

  if ( theReal ) {
    getDipoles();
  }

  if ( verbose() )
    print(Repository::clog());

  MEGroup::doinit();

}

void SubtractedME::doinitrun() {

  // has been deactivated by the factory
  if ( !head() ) {
    MEBase::doinitrun();
    return;
  }

  MEGroup::doinitrun();

  theReal = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(head());

}

void SubtractedME::dofinish() {

  // has been deactivated by the factory
  if ( !head() ) {
    MEBase::dofinish();
    return;
  }

  MEGroup::dofinish();

  for ( map<CollinearSubtractionIndex,SubtractionHistogram>::
	  const_iterator b = collinearHistograms.begin();
	b != collinearHistograms.end(); ++b ) {
    b->second.dump(factory()->subtractionData(),
		   b->first.first,
		   b->first.second.first,
		   b->first.second.second);
  }

  for ( map<SoftSubtractionIndex,SubtractionHistogram>::
	  const_iterator b = softHistograms.begin();
	b != softHistograms.end(); ++b ) {
    b->second.dump(factory()->subtractionData(),
		   b->first.first,
		   b->first.second,
		   b->first.second);
  }

}

void SubtractedME::print(ostream& os) const {

  os << "--- SubtractedME setup ---------------------------------------------------------\n";

  os << " '" << name() << "' subtracting real emission\n '"
     << head()->name() << "' using the dipoles:\n";

  for ( MEVector::const_iterator d = dependent().begin();
	d != dependent().end(); ++d )
    dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(*d)->print(os);

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void SubtractedME::printLastEvent(ostream& os) const {

  os << "--- SubtractedME last event information ----------------------------------------\n";

  os << " for subtracted matrix element '" << name() << "'\n";

  os << " real emission event information:\n";
  dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(head())->printLastEvent(os);

  os << " dipoles event information:\n";
  for ( MEVector::const_iterator d = dependent().begin();
	d != dependent().end(); ++d )
    dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(*d)->printLastEvent(os);


  os << "--- end SubtractedME last event information ------------------------------------\n\n\n";

  os << flush;

}

void SubtractedME::lastEventStatistics() {
  MEGroup::lastEventStatistics();
  if ( !generator() )
    return;
  /*
  if ( verbose() )
    printLastEvent(generator()->log());
  */
  if ( !collinearHistograms.empty() )
    lastEventSubtraction();
}

SubtractedME::SubtractionHistogram::
SubtractionHistogram(double low, 
		     double up, 
		     unsigned int nbins) 
  : lower(low) {
  nbins = nbins + 1;

  double c = log10(up/low) / (nbins-1.);

  for ( unsigned int k = 1; k < nbins; ++k ) {
    bins[low*pow(10.0,k*c)] = make_pair(Constants::MaxDouble,0.);
  }

}

void SubtractedME::SubtractionHistogram::
dump(const std::string& prefix, 
     const cPDVector& proc,
     int i, int j) const {
  ostringstream fname("");
  for ( cPDVector::const_iterator p = proc.begin();
	p != proc.end(); ++p )
    fname << (**p).PDGName();
  fname << "-" << i << "-" << j;
  ofstream out((prefix+fname.str()+".dat").c_str());
  for ( map<double,pair<double,double> >::const_iterator b = bins.begin();
	b != bins.end(); ++b ) {
    map<double,pair<double,double> >::const_iterator bp = b; --bp;
    if ( b->second.first != Constants::MaxDouble &&
	 b->second.second != 0. ) {
      if ( b != bins.begin() )
	out << bp->first;
      else
	out << lower;
      out << " " << b->first
	  << " " << b->second.first
	  << " " << b->second.second
	  << "\n" << flush;
    }
  }
  ofstream gpout((prefix+fname.str()+".gp").c_str());
  gpout << "set terminal epslatex color solid;\n"
	<< "set output '" << fname.str() << "-plot.tex';\n"
	<< "set log x;\n"
	<< "set size 0.5,0.6;\n"
	<< "set yrange [0:2];\n"
	<< "set xrange [0.001:10];\n";
  if ( i != j ) {
    gpout << "set xlabel '$\\sqrt{s_{" << i << j << "}}/{\\rm GeV}$'\n";
  } else {
    gpout << "set xlabel '$E_{" << i << "}/{\\rm GeV}$'\n";
  }
  gpout << "plot 1 w lines lc rgbcolor \"#DDDDDD\" notitle, '" << fname.str() 
	<< ".dat' u (($1+$2)/2.):3:($4 < 4. ? $4 : 4.) w filledcurves lc rgbcolor \"#00AACC\" t "
	<< "'$";
  for ( size_t k = 0; k < proc.size(); k++ ) {
    if ( k == 2 )
      gpout << "\\to ";
    gpout << (proc[k]->id() < 0 ? "\\bar{" : "")
	  << (proc[k]->id() < 0 ? proc[k]->CC()->PDGName() : proc[k]->PDGName())
	  << (proc[k]->id() < 0 ? "}" : "") << " ";
  }
  gpout << "$';\n";
  gpout << "reset;\n";
}

void SubtractedME::lastEventSubtraction() {

  tStdXCombGroupPtr xc = dynamic_ptr_cast<tStdXCombGroupPtr>(lastXCombPtr());

  CrossSection xcme2 = xc->lastHeadCrossSection();
  CrossSection xcdip = ZERO;

  if ( xcme2 == ZERO )
    return;

  for ( vector<StdXCombPtr>::const_iterator d = xc->dependent().begin();
	d != xc->dependent().end(); ++d ) {
    if ( !(*d) )
      continue;
    if ( !(**d).matrixElement()->apply() )
      continue;
    xcdip += (**d).lastCrossSection();
  }

  if ( theReal->phasespace() ) {
    size_t i = lastSingularLimit()->first;
    size_t j = lastSingularLimit()->second;
    if ( i == j && 
	 softHistograms.find(SoftSubtractionIndex(head()->mePartonData(),i))
	 != softHistograms.end() ) {
      softHistograms[SoftSubtractionIndex(head()->mePartonData(),i)].
	book(meMomenta()[i].t()/GeV,abs(xcdip)/abs(xcme2));
    }
    if ( i != j &&
	 collinearHistograms.find(CollinearSubtractionIndex(head()->mePartonData(),make_pair(i,j))) 
	 != collinearHistograms.end() ) {
      double s = sqrt(2.*meMomenta()[i]*meMomenta()[j])/GeV;
      collinearHistograms[CollinearSubtractionIndex(head()->mePartonData(),make_pair(i,j))].
	book(s,abs(xcdip)/abs(xcme2));
    }
    return;
  }

  for ( size_t i = 0; i < meMomenta().size(); ++i ) {
    if ( i > 1 ) {
      if ( softHistograms.find(SoftSubtractionIndex(head()->mePartonData(),i))
	   != softHistograms.end() ) {
	softHistograms[SoftSubtractionIndex(head()->mePartonData(),i)].
	  book(meMomenta()[i].t()/GeV,abs(xcdip)/abs(xcme2));
      }
    }
    for ( size_t j = i+1; j < meMomenta().size(); ++j ) {
      if ( collinearHistograms.find(CollinearSubtractionIndex(head()->mePartonData(),make_pair(i,j))) 
	   == collinearHistograms.end() )
	continue;
      double s = sqrt(2.*meMomenta()[i]*meMomenta()[j])/GeV;
      collinearHistograms[CollinearSubtractionIndex(head()->mePartonData(),make_pair(i,j))].
	book(s,abs(xcdip)/abs(xcme2));
    }
  }

}

void SubtractedME::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb
     << theFactory << theReal 
     << theRealShowerSubtraction << theVirtualShowerSubtraction;
}

void SubtractedME::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb
     >> theFactory >> theReal 
     >> theRealShowerSubtraction >> theVirtualShowerSubtraction;
  lastMatchboxXComb(theLastXComb);
}

void SubtractedME::Init() {

  static ClassDocumentation<SubtractedME> documentation
    ("SubtractedME represents a subtracted real emission matrix element.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SubtractedME,MEGroup>
describeHerwigSubtractedME("Herwig::SubtractedME", "HwMatchbox.so");
