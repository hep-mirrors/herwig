// -*- C++ -*-
//
// OpenLoopsAmplitude.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OpenLoopsAmplitude class.
//

#include "OpenLoopsAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Utilities/DynamicLoader.h"

#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

using namespace Herwig;

#ifndef OPENLOOPSLIBS
#error Makefile.am needs to define OPENLOOPSLIBS
#endif

#ifndef OPENLOOPSPREFIX
#error Makefile.am needs to define OPENLOOPSPREFIX
#endif

OpenLoopsAmplitude::OpenLoopsAmplitude() :
  theHiggsEff(false), use_cms(true), psp_tolerance(12),
  OpenLoopsLibs_(OPENLOOPSLIBS), OpenLoopsPrefix_(OPENLOOPSPREFIX) {
}

OpenLoopsAmplitude::~OpenLoopsAmplitude() {
}

IBPtr OpenLoopsAmplitude::clone() const {
	return new_ptr(*this);
}

IBPtr OpenLoopsAmplitude::fullclone() const {
	return new_ptr(*this);
}

extern "C" void OLP_Start(const char*, int* i);
extern "C" void OLP_SetParameter(const char* ,double* ,double*,int*);
extern "C" void ol_setparameter_string(const char*, const char*);
extern "C" void OLP_PrintParameter(const char*);
extern "C" void OLP_EvalSubProcess(int*, double*, double*, double*, double*);
extern "C" void OLP_EvalSubProcess2(int*, double*, double*, double*, double*);
                                  // id  ps-point emitter polvec res
extern "C" void ol_evaluate_sc(int, double*, int, double*, double*);
extern "C" void OLP_Polvec(double*,double*,double*);

void OpenLoopsAmplitude::doinitrun() {
  MatchboxOLPME::doinitrun();
}

void OpenLoopsAmplitude::startOLP(const string& contract, int& status) {
	string tempcontract=contract;

	if ( ! (DynamicLoader::load(OpenLoopsLibs_+"/libopenloops.so") ||
		DynamicLoader::load(OpenLoopsPrefix_+"/lib/libopenloops.so") ||
		DynamicLoader::load("libopenloops.so") ||
		DynamicLoader::load(OpenLoopsLibs_+"/libopenloops.dylib") ||
		DynamicLoader::load(OpenLoopsPrefix_+"/lib/libopenloops.dylib") ||
		DynamicLoader::load("libopenloops.dylib") ) ) {
	  throw Exception() << "OpenLoopsAmplitude::startOLP(): Failed to load libopenloops.so/dylib\n"
			    << DynamicLoader::lastErrorMessage
			    << Exception::runerror;
	}

	string stabilityPrefix = factory()->runStorage() + "OpenLoops.StabilityLog";
	assert(stabilityPrefix.size() < 256);

	ol_setparameter_string("stability_logdir",stabilityPrefix.c_str());

	ol_setparameter_string("install_path",OpenLoopsPrefix_.c_str());

	int a=0;double null=0.0;double one=1.0;
	int part[10]={1,2,3,4,5,6,15,23,24,25};string stri;
	for (int i=0;i<10;i++){
	 map<long,Energy>::const_iterator it=reshuffleMasses().find(part[i]);
	 double mass;
	 if(it==reshuffleMasses().end())
	   mass = getParticleData(part[i])->hardProcessMass()/GeV;
	 else
	   mass = it->second/GeV;
	 double width=getParticleData(part[i])->hardProcessWidth()/GeV;
	 std::stringstream ss;
	 ss << part[i];
	 string str = ss.str();
	 stri="mass("+str+")";
	 OLP_SetParameter(stri.c_str(),&mass,&null,&a);
     stri="width("+str+")";
     OLP_SetParameter(stri.c_str(),&width,&null,&a);
	}
  stri="alphas";
  one=SM().alphaS();
  OLP_SetParameter( stri.c_str(),&one ,&null,&a);
  stri="alpha";
  one=SM().alphaEMMZ();
  OLP_SetParameter(stri.c_str(),&one ,&null,&a);

	OLP_Start(tempcontract.c_str(), &status);

  didStartOLP() = true;

}

void OpenLoopsAmplitude::fillOrderFile(const map<pair<Process, int>, int>& procs) {
	string orderFileName = factory()->buildStorage() + name() + ".OLPContract.lh";
	ofstream orderFile(orderFileName.c_str());
	size_t asPower = 100;
	size_t minlegs = 100;
	size_t maxlegs = 0;

	for ( map<pair<Process, int>, int>::const_iterator t = procs.begin() ; t != procs.end() ; ++t ) {
		asPower = min(asPower, static_cast<size_t>(t->first.first.orderInAlphaS));
		minlegs = min(minlegs, static_cast<size_t>(t->first.first.legs.size()));
		maxlegs = max(maxlegs, static_cast<size_t>(t->first.first.legs.size()));
	}

	orderFile << "# OLP order file created by Herwig/Matchbox for OpenLoops\n\n";
	orderFile << "CorrectionType           QCD\n";
	orderFile << "IRregularization         " << (isDR() ? "DRED" : "CDR") << "\n";
	orderFile << "extra answerfile      " << (factory()->buildStorage() + name() + ".OLPAnswer.lh") << "\n";
	orderFile << "extra psp_tolerance "<<psp_tolerance<<"\n";
	orderFile << "extra use_cms "<<(use_cms?"1":"0")<< "\n";
	if (theHiggsEff)
		orderFile << "model heft\n";
	orderFile << "\n";

	if (extraOpenLoopsPath!="")
	    orderFile << "Extra OpenLoopsPath  " << extraOpenLoopsPath << "\n";
	for ( map<pair<Process, int>, int>::const_iterator p = procs.begin() ; p != procs.end() ; ++p ) {
		std::stringstream Processstr;
		std::stringstream Typestr;
		Processstr << (*p).first.first.legs[0]->id() << " " << (*p).first.first.legs[1]->id() << " -> ";
		for ( PDVector::const_iterator o = (*p).first.first.legs.begin() + 2 ; o != (*p).first.first.legs.end() ; ++o )
			Processstr << (**o).id() << " ";
		if ( (*p).first.second == ProcessType::treeME2 ) {
			Typestr << "Tree";
		} else if ( (*p).first.second == ProcessType::colourCorrelatedME2 ) {
			Typestr << "ccTree";
		} else if ( (*p).first.second == ProcessType::spinColourCorrelatedME2 ) {
			Typestr << "sctree_polvect";
		} else if ( (*p).first.second == ProcessType::oneLoopInterference ) {
			Typestr << "Loop";
		}
		OpenLoopsProcInfo pro = OpenLoopsProcInfo((*p).second, -1, Processstr.str(), Typestr.str());
		pro.setOAs(p->first.first.orderInAlphaS);
		processmap[(*p).second] = pro;
	}

	vector < string > types;
	types.push_back("Tree");
	types.push_back("ccTree");
	types.push_back("sctree_polvect");
	types.push_back("Loop");
	for ( size_t i = asPower ; i != asPower + maxlegs - minlegs + 1 ; i++ ) {
		orderFile << "\n\nCouplingPower QCD       " << i;
		orderFile << "\n\n#AlphasPower            " << i;
		for ( vector<string>::iterator it = types.begin() ; it != types.end() ; it++ ) {
			for ( map<int, OpenLoopsProcInfo>::iterator p = processmap.begin() ; p != processmap.end() ; ++p )
			  if ( (*p).second.Tstr() == *it && i == (unsigned int) (*p).second.orderAs() ) {
					orderFile << "\nAmplitudeType " << *it << "\n";
					break;
				}
			for ( map<int, OpenLoopsProcInfo>::iterator p = processmap.begin() ; p != processmap.end() ; ++p )
			  if ( (*p).second.Tstr() == *it && i == (unsigned int) (*p).second.orderAs() ) {
					orderFile << (*p).second.Pstr() << "\n";
				}
		}
	}
	orderFile << flush;
}

bool OpenLoopsAmplitude::checkOLPContract() {
	string contractFileName = factory()->buildStorage() + name() + ".OLPAnswer.lh";
        ifstream infile(contractFileName.c_str());
	string line;
	vector < string > contractfile;
	while (std::getline(infile, line)) {
		contractfile.push_back(line);
	}
	for ( map<int, OpenLoopsProcInfo>::iterator p = processmap.begin() ; p != processmap.end() ; p++ ) {
		bool righttype = false;
		for ( vector<string>::iterator linex = contractfile.begin() ; linex != contractfile.end() ; ++linex ) {
			if ( (*linex).find("AmplitudeType ")!= std::string::npos ) {
				if ( (*linex).find(" " + (*p).second.Tstr() + " ")!= std::string::npos ) {
					righttype = true;
				} else {
					righttype = false;
				}
			}
			if ( righttype ) {
					if ( (*linex).find((*p).second.Pstr()) != std::string::npos ){
					if( (*p).second.Pstr().length() == (*linex).find("|") ) {
					string sub = (*linex).substr((*linex).find("|") + 1, (*linex).find("#") - (*linex).find("|") - 1); // | 1 23 # buggy??

					int subint;
					int subint2;
					istringstream(sub) >> subint >> subint2;
					assert(subint==1);
					(*p).second.setGID(subint2);
				}
}
			}
		}
	}
	string ids = factory()->buildStorage() + "OpenLoops.ids.dat";
	ofstream IDS(ids.c_str());

	for ( map<int, OpenLoopsProcInfo>::iterator p = processmap.begin() ; p != processmap.end() ; p++ ) {
	    idpair.insert ( std::pair<int,int>((*p).second.HID(),(*p).second.GID()) );
	    IDS << (*p).second.HID() << " " << (*p).second.GID() << "\n";
	    if ( (*p).second.GID() == -1 ) return 0;
	}
	IDS << flush;
	return 1;
}

void OpenLoopsAmplitude::getids() const{
	string line = factory()->buildStorage() + "OpenLoops.ids.dat";
	ifstream infile(line.c_str());
	int hid;
	int gid;
	while (std::getline(infile, line)) {
	   istringstream(line) >> hid>>gid;
	   idpair.insert ( std::pair<int,int>(hid,gid) );
	}
}
bool OpenLoopsAmplitude::startOLP(const map<pair<Process, int>, int>& procs) {
	string contractFileName =  factory()->buildStorage() + name() + ".OLPAnswer.lh";
	string orderFileName = factory()->buildStorage() + name() + ".OLPContract.lh";
	fillOrderFile(procs);
 	int status = -1;
	startOLP(orderFileName, status);
	if ( !checkOLPContract() ) {
	  return false;
	}
	if ( status != 1 ) return false;
	return true;
}


void OpenLoopsAmplitude::evalSubProcess() const {
  useMe();
	double units = pow(lastSHat() / GeV2, mePartonData().size() - 4.);
	fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
	double acc ;
	double scale = sqrt(mu2() / GeV2);
  if (hasRunningAlphaS()) {
    int a=0;double null=0.0;double one=1.0;
    string stri="alphas";
    one=lastAlphaS();
    OLP_SetParameter( stri.c_str(),&one ,&null,&a);
  }
  double out[7]={};

	int id = olpId()[ProcessType::oneLoopInterference] ? olpId()[ProcessType::oneLoopInterference] : olpId()[ProcessType::treeME2];

	if ( idpair.size() == 0 ) {
		getids();
		if ( Debug::level > 1 ) {
		  string parfile=factory()->runStorage() + name() + ".Parameters.dat";
		  OLP_PrintParameter(parfile.c_str());
		}
		
    }
  

	OLP_EvalSubProcess2(&((*(idpair.find(id))).second), olpMomenta(), &scale, out,&acc );
  
  
	if ( olpId()[ProcessType::oneLoopInterference] ) {
		if(calculateTreeME2())lastTreeME2(out[3] * units);
		lastOneLoopInterference((out[2])* units);
		lastOneLoopPoles(pair<double, double>(out[0] * units, out[1] * units));
	} else if ( olpId()[ProcessType::treeME2] ) {
		lastTreeME2(out[0] * units);
	}
}

void OpenLoopsAmplitude::evalColourCorrelator(pair<int, int>  ) const {
	double units = pow(lastSHat() / GeV2, mePartonData().size() - 4.);
	fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
	double acc ;
	double scale = sqrt(mu2() / GeV2);
  if (hasRunningAlphaS()) {
    int a=0;double null=0.0;double one=1.0;
    string stri="alphas";
    one=lastAlphaS();
    OLP_SetParameter( stri.c_str(),&one ,&null,&a);
  }
	int n = lastXComb().meMomenta().size();

	colourCorrelatorResults.resize(n * (n - 1) / 2);
	if ( idpair.size() == 0 ) {
		getids();
		if ( Debug::level > 1 ) {
		  string parfile=factory()->runStorage() + name() + ".Parameters.dat";
		  OLP_PrintParameter(parfile.c_str());
		}
		
	}


	int id = olpId()[ProcessType::colourCorrelatedME2];

	OLP_EvalSubProcess2(&((*(idpair.find(id))).second), olpMomenta(), &scale, &colourCorrelatorResults[0],&acc );
	for ( int i = 0 ; i < n ; ++i ){
	    for ( int j = i + 1 ; j < n ; ++j ) {
			lastColourCorrelator(make_pair(i, j), colourCorrelatorResults[i+j*(j-1)/2] * units);
	    }
	}
}

void OpenLoopsAmplitude::evalSpinColourCorrelator(pair<int , int > ) const {
        assert(false);
}


double OpenLoopsAmplitude::spinColourCorrelatedME2(pair<int,int> ij,
					 const SpinCorrelationTensor& c) const{

	double units = pow(lastSHat() / GeV2, mePartonData().size() - 4.);
	fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());

  if (hasRunningAlphaS()) {
    int a=0;double null=0.0;double one=1.0;
    string stri="alphas";
    one=lastAlphaS();
    OLP_SetParameter( stri.c_str(),&one ,&null,&a);
  }

	int emitter=ij.first+1;
	int n = lastXComb().meMomenta().size();

	if ( idpair.size() == 0 ) {
		getids();
		if ( Debug::level > 1 ) {
		  string parfile=factory()->runStorage() + name() + ".Parameters.dat";
		  OLP_PrintParameter(parfile.c_str());
		}
		
	}
	int id = (*(idpair.find(olpId()[ProcessType::spinColourCorrelatedME2]))).second;
	//double * outx =new double[n];
	spinColourCorrelatorResults.resize(n);
        double polvec[4];

	polvec[0]=c.momentum().e()/GeV;
	polvec[1]=c.momentum().x()/GeV;
	polvec[2]=c.momentum().y()/GeV;
	polvec[3]=c.momentum().z()/GeV;
        double avg= colourCorrelatedME2(ij)*(-c.diagonal());

	ol_evaluate_sc(id, olpMomenta(),emitter,polvec,&spinColourCorrelatorResults[0]);

        double corr =-1.*units * spinColourCorrelatorResults[ij.second]/c.scale()*c.momentum().dot(c.momentum());


  double Nc = generator()->standardModel()->Nc();
  double cfac = 1.;
  if ( mePartonData()[ij.first]->iColour() == PDT::Colour8 ) {
    cfac = Nc;
  } else if ( mePartonData()[ij.first]->iColour() == PDT::Colour3 ||
	      mePartonData()[ij.first]->iColour() == PDT::Colour3bar ) {
    cfac = (sqr(Nc)-1.)/(2.*Nc);
  } else assert(false);


  return
    avg + corr/cfac;

}





// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void OpenLoopsAmplitude::persistentOutput(PersistentOStream & os) const {
  os << idpair << OpenLoopsLibs_ << OpenLoopsPrefix_;
}

void OpenLoopsAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> idpair >> OpenLoopsLibs_ >> OpenLoopsPrefix_;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<OpenLoopsAmplitude, MatchboxOLPME> describeHerwigOpenLoopsAmplitude("Herwig::OpenLoopsAmplitude", "HwMatchboxOpenLoops.so");

void OpenLoopsAmplitude::Init() {

  static ClassDocumentation<OpenLoopsAmplitude> 
    documentation("OpenLoopsAmplitude implements an interface to OpenLoops.",
		  "Matrix elements have been calculated using OpenLoops \\cite{Cascioli:2011va}",
		  "%\\cite{Cascioli:2011va}\n"
		  "\\bibitem{Cascioli:2011va}\n"
		  "F.~Cascioli et al.,\n"
		  "``Scattering Amplitudes with Open Loops,''\n"
		  "arXiv:1111.5206 [hep-ph].\n"
		  "%%CITATION = ARXIV:1111.5206;%%");
  
  static Switch<OpenLoopsAmplitude,bool> interfaceHiggsEff
         ("HiggsEff",
          "Switch On/Off for effective higgs model.",
          &OpenLoopsAmplitude::theHiggsEff, false, false, false);
  static SwitchOption interfaceHiggsEffOn
         (interfaceHiggsEff,
          "On",
          "On",
          true);
  static SwitchOption interfaceHiggsEffOff
         (interfaceHiggsEff,
          "Off",
          "Off",
          false);

  static Switch<OpenLoopsAmplitude,bool> interfaceUseComplMass
  ("ComplexMassScheme",
   "Switch on or off if Complex Masses.",
   &OpenLoopsAmplitude::use_cms, true, false, false);
  static SwitchOption interfaceUseComplMassOn
  (interfaceUseComplMass,
   "True",
   "True for Complex Masses.",
   true);
  static SwitchOption interfaceUseComplMassOff
  (interfaceUseComplMass,
   "False",
   "False for no Complex Masses.",
   false);
  
  static Parameter<OpenLoopsAmplitude,int> interfacepsp_tolerance
  ("PSP_tolerance",
   "(Debug)Phase Space Tolerance. Better use e.g.: set OpenLoops:Massless 13",
   &OpenLoopsAmplitude::psp_tolerance, 12, 0, 0,
   false, false, Interface::lowerlim);
    
  static Parameter<OpenLoopsAmplitude,string> interfaceOpenLoopsLibs
    ("OpenLoopsLibs",
     "The location of OpenLoops libraries",
     &OpenLoopsAmplitude::OpenLoopsLibs_, string(OPENLOOPSLIBS),
     false, false);
    
  static Parameter<OpenLoopsAmplitude,string> interfaceOpenLoopsPrefix
    ("OpenLoopsPrefix",
     "The location of OpenLoops libraries",
     &OpenLoopsAmplitude::OpenLoopsPrefix_, string(OPENLOOPSPREFIX),
     false, false);
  
}

