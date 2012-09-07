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

using namespace Herwig;

SubtractedME::SubtractedME() 
  : MEGroup(), 
    theSubtractionData(""), 
    theVerbose(false), theSubProcessGroups(false),
    theVetoScales(false) {}

SubtractedME::~SubtractedME() {}

IBPtr SubtractedME::clone() const {
  return new_ptr(*this);
}

IBPtr SubtractedME::fullclone() const {
  return new_ptr(*this);
}

void SubtractedME::setXComb(tStdXCombPtr xc) {
  MEGroup::setXComb(xc);
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

  if ( allDipoles().empty() ) {
    allDipoles() = DipoleRepository::dipoles();
  }

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

  if ( allDipoles().empty() ) {
    allDipoles() = DipoleRepository::dipoles();
  }

  Ptr<MatchboxMEBase>::tptr real =
    dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(head());

  if ( allDipoles().empty() || theBorns.empty() || !real )
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
    = real->getDipoles(allDipoles(),theBorns);

  if ( theSubtractionData != "" ) {
    theSubProcessGroups = true;
    for ( vector<Ptr<SubtractionDipole>::ptr>::const_iterator d =
	    genDipoles.begin(); d != genDipoles.end(); ++d )
      (**d).doTestSubtraction();
  }

  if ( genDipoles.empty() )
    Throw<InitException>() << "The SubtractedME '"
			   << name() << "' could not generate "
			   << "subtraction terms for the real emission "
			   << "matrix element '" << real->name() << "'. "
			   << "Please check the corresponding input file.";

  dipMEs.resize(genDipoles.size());
  copy(genDipoles.begin(),genDipoles.end(),dipMEs.begin());

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

void SubtractedME::setVetoScales(tSubProPtr subpro) const {

  if ( !vetoScales() )
    return;

  Ptr<SubtractionDipole>::tptr dipole;
  Energy pt;

  assert(head()->noMirror());

  for ( MEVector::const_iterator d = dependent().begin();
	d != dependent().end(); ++d ) {

    dipole = dynamic_ptr_cast<Ptr<SubtractionDipole>::ptr>(*d);
    assert(dipole);
    pt = dipole->lastPt();

    if ( dipole->realEmitter() == 0 ||
	 dipole->realSpectator() == 0 ) {
      if ( subpro->incoming().first->vetoScale() < 0.0*GeV2 ||
	   subpro->incoming().first->vetoScale() > sqr(pt) )
	subpro->incoming().first->vetoScale(sqr(pt));
    }

    if ( dipole->realEmitter() == 1 ||
	 dipole->realSpectator() == 1 ) {
      if ( subpro->incoming().second->vetoScale() < 0.0*GeV2 ||
	   subpro->incoming().second->vetoScale() > sqr(pt) )
	subpro->incoming().second->vetoScale(sqr(pt));
    }

    if ( dipole->realEmitter() > 1 ) {
      if ( subpro->outgoing()[dipole->realEmitter()-2]->vetoScale() < 0.0*GeV2 ||
	   subpro->outgoing()[dipole->realEmitter()-2]->vetoScale() > sqr(pt) )
	subpro->outgoing()[dipole->realEmitter()-2]->vetoScale(sqr(pt));
    }

    if ( dipole->realSpectator() > 1 ) {
      if ( subpro->outgoing()[dipole->realSpectator()-2]->vetoScale() < 0.0*GeV2 ||
	   subpro->outgoing()[dipole->realSpectator()-2]->vetoScale() > sqr(pt) )
	subpro->outgoing()[dipole->realSpectator()-2]->vetoScale(sqr(pt));
    }

    if ( subpro->outgoing()[dipole->realEmission()-2]->vetoScale() < 0.0*GeV2 ||
	 subpro->outgoing()[dipole->realEmission()-2]->vetoScale() > sqr(pt) )
      subpro->outgoing()[dipole->realEmission()-2]->vetoScale(sqr(pt));  

  }

}

void SubtractedME::doinit() {

  theReal = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(head());

  if ( theReal ) {
    getDipoles();
  }

  if ( theVerbose )
    print(Repository::clog());

  MEGroup::doinit();
}

void SubtractedME::doinitrun() {

  MEGroup::doinitrun();

  theReal = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(head());

  if ( theSubtractionData != "" ) {
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
	    theReal->phasespace()->singularLimit(i,i);
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
	    theReal->phasespace()->singularLimit(i,j);
	}
      }
    }
  }

}

void SubtractedME::dofinish() {

  MEGroup::dofinish();

  for ( map<CollinearSubtractionIndex,SubtractionHistogram>::
	  const_iterator b = collinearHistograms.begin();
	b != collinearHistograms.end(); ++b ) {
    b->second.dump(theSubtractionData,
		   b->first.first,
		   b->first.second.first,
		   b->first.second.second);
  }

  for ( map<SoftSubtractionIndex,SubtractionHistogram>::
	  const_iterator b = softHistograms.begin();
	b != softHistograms.end(); ++b ) {
    b->second.dump(theSubtractionData,
		   b->first.first,
		   b->first.second,
		   b->first.second);
  }

}

void SubtractedME::rebind(const TranslationMap & trans) {
  for ( vector<Ptr<SubtractionDipole>::ptr>::iterator d =
	  theDipoles.begin(); d != theDipoles.end(); ++d )
    *d = trans.translate(*d);
  MEGroup::rebind(trans);
}

IVector SubtractedME::getReferences() {
  IVector ret = MEGroup::getReferences();
  for ( vector<Ptr<SubtractionDipole>::ptr>::const_iterator d =
	  theDipoles.begin(); d != theDipoles.end(); ++d )
    ret.push_back(*d);
  return ret;
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
  if ( theVerbose )
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

  for ( StdDepXCVector::const_iterator d = xc->dependent().begin();
	d != xc->dependent().end(); ++d ) {
    if ( !(*d) )
      continue;
    if ( !(**d).matrixElement()->apply() )
      continue;
    xcdip += (**d).lastCrossSection();
  }

  if ( theReal->phasespace() ) {
    size_t i = theReal->phasespace()->lastSingularLimit().first;
    size_t j = theReal->phasespace()->lastSingularLimit().second;
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
  os << theDipoles << theBorns << theSubtractionData 
     << theVerbose << theSubProcessGroups << theVetoScales;
}

void SubtractedME::persistentInput(PersistentIStream & is, int) {
  is >> theDipoles >> theBorns >> theSubtractionData 
     >> theVerbose >> theSubProcessGroups >> theVetoScales;
}

void SubtractedME::Init() {

  static ClassDocumentation<SubtractedME> documentation
    ("SubtractedME represents a subtracted real emission matrix element.");

  static RefVector<SubtractedME,MatchboxMEBase> interfaceBorns
    ("Borns",
     "The underlying Born matrix elements to be considered",
     &SubtractedME::theBorns, -1, false, false, true, true, false);


  static Parameter<SubtractedME,string> interfaceSubtractionData
    ("SubtractionData",
     "File to dump subtraction check to.",
     &SubtractedME::theSubtractionData, "",
     false, false);


  static Switch<SubtractedME,bool> interfaceVerbose
    ("Verbose",
     "Print full infomation on each evaluated phase space point.",
     &SubtractedME::theVerbose, false, false, false);
  static SwitchOption interfaceVerboseOn
    (interfaceVerbose,
     "On",
     "On",
     true);
  static SwitchOption interfaceVerboseOff
    (interfaceVerbose,
     "Off",
     "Off",
     false);

  static Switch<SubtractedME,bool> interfaceSubProcessGroups
    ("SubProcessGroups",
     "Switch on or off production of sub-process groups.",
     &SubtractedME::theSubProcessGroups, false, false, false);
  static SwitchOption interfaceSubProcessGroupsOn
    (interfaceSubProcessGroups,
     "On",
     "On",
     true);
  static SwitchOption interfaceSubProcessGroupsOff
    (interfaceSubProcessGroups,
     "Off",
     "Off",
     false);

  static Switch<SubtractedME,bool> interfaceVetoScales
    ("VetoScales",
     "Switch on or off production of sub-process groups.",
     &SubtractedME::theVetoScales, false, false, false);
  static SwitchOption interfaceVetoScalesOn
    (interfaceVetoScales,
     "On",
     "On",
     true);
  static SwitchOption interfaceVetoScalesOff
    (interfaceVetoScales,
     "Off",
     "Off",
     false);


}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SubtractedME,MEGroup>
describeHerwigSubtractedME("Herwig::SubtractedME", "HwMatchbox.so");
