// -*- C++ -*-
//
// SusyBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SusyBase class.
//

#include "SusyBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;

SusyBase::SusyBase() : readFile_(false), MPlanck_(2.4e18*GeV),
		       gravitino_(false), majoranaNeutrinos_(false),
		       tanBeta_(0), mu_(ZERO), 
		       M1_(ZERO), M2_(ZERO), M3_(ZERO),
		       mH12_(ZERO),mH22_(ZERO),
		       meL_(ZERO),mmuL_(ZERO),mtauL_(ZERO),
		       meR_(ZERO),mmuR_(ZERO),mtauR_(ZERO),
		       mq1L_(ZERO),mq2L_(ZERO),mq3L_(ZERO),
		       mdR_(ZERO),muR_(ZERO),msR_(ZERO),
		       mcR_(ZERO),mbR_(ZERO),mtR_(ZERO),
		       gluinoPhase_(1.) 
{}

IBPtr SusyBase::clone() const {
  return new_ptr(*this);
}

IBPtr SusyBase::fullclone() const {
  return new_ptr(*this);
}

void SusyBase::doinit() {
  addVertex(WSFSFVertex_);
  addVertex(NFSFVertex_);
  addVertex(GFSFVertex_);
  addVertex(HSFSFVertex_);
  addVertex(CFSFVertex_);
  addVertex(GSFSFVertex_);
  addVertex(GGSQSQVertex_);
  addVertex(WGSQSQVertex_);
  addVertex(GSGSGVertex_);
  addVertex(NNZVertex_);
  if(NNPVertex_) addVertex(NNPVertex_);
  if(GNGVertex_) addVertex(GNGVertex_);
  addVertex(CCZVertex_);
  addVertex(CNWVertex_);
  addVertex(GOGOHVertex_);
  addVertex(WHHVertex_);
  if(NCTVertex_) addVertex(NCTVertex_);
  if(gravitino_) {
    if(GVNHVertex_) addVertex(GVNHVertex_);
    if(GVNVVertex_) addVertex(GVNVVertex_);
    if(GVFSVertex_) addVertex(GVFSVertex_);
  }
  if(majoranaNeutrinos()) {
    idMap().insert(make_pair(12,17));
    idMap().insert(make_pair(14,18));
    idMap().insert(make_pair(16,19));
  }
  BSMModel::doinit();
}

void SusyBase::persistentOutput(PersistentOStream & os) const {
  os << readFile_ << gravitino_
     << NMix_ << UMix_ << VMix_ << WSFSFVertex_ 
     << NFSFVertex_ << GFSFVertex_ << HSFSFVertex_ << CFSFVertex_ 
     << GSFSFVertex_ << GGSQSQVertex_ << WGSQSQVertex_ << GSGSGVertex_ 
     << NNZVertex_ << NNPVertex_ << CCZVertex_ << CNWVertex_ 
     << GOGOHVertex_ << WHHVertex_ << GNGVertex_ << NCTVertex_
     << GVNHVertex_ << GVNVVertex_ << GVFSVertex_
     << tanBeta_ << ounit(mu_,GeV) 
     << ounit(M1_,GeV) << ounit(M2_,GeV) << ounit(M3_,GeV)
     << ounit(mH12_,GeV2) << ounit(mH22_,GeV2) 
     << ounit(meL_,GeV)  << ounit(mmuL_,GeV) << ounit(mtauL_,GeV) 
     << ounit(meR_,GeV)  << ounit(mmuR_,GeV) << ounit(mtauR_,GeV) 
     << ounit(mq1L_,GeV) << ounit(mq2L_,GeV) << ounit(mq3L_,GeV) 
     << ounit(mdR_,GeV)  << ounit(muR_,GeV)  << ounit(msR_,GeV) 
     << ounit(mcR_,GeV)  << ounit(mbR_,GeV)  << ounit(mtR_,GeV)
     << gluinoPhase_ << ounit(MPlanck_,GeV);
}

void SusyBase::persistentInput(PersistentIStream & is, int) {
  is >> readFile_  >> gravitino_
     >> NMix_ >> UMix_ >> VMix_ >> WSFSFVertex_ 
     >> NFSFVertex_ >> GFSFVertex_ >> HSFSFVertex_ >> CFSFVertex_ 
     >> GSFSFVertex_ >> GGSQSQVertex_ >> WGSQSQVertex_ >> GSGSGVertex_ 
     >> NNZVertex_ >> NNPVertex_ >> CCZVertex_ >> CNWVertex_
     >> GOGOHVertex_ >> WHHVertex_ >> GNGVertex_ >> NCTVertex_
     >> GVNHVertex_ >> GVNVVertex_ >> GVFSVertex_
     >> tanBeta_ >> iunit(mu_,GeV) 
     >> iunit(M1_,GeV) >> iunit(M2_,GeV) >> iunit(M3_,GeV)
     >> iunit(mH12_,GeV2) >> iunit(mH22_,GeV2) 
     >> iunit(meL_,GeV)  >> iunit(mmuL_,GeV) >> iunit(mtauL_,GeV) 
     >> iunit(meR_,GeV)  >> iunit(mmuR_,GeV) >> iunit(mtauR_,GeV) 
     >> iunit(mq1L_,GeV) >> iunit(mq2L_,GeV) >> iunit(mq3L_,GeV) 
     >> iunit(mdR_,GeV)  >> iunit(muR_,GeV)  >> iunit(msR_,GeV) 
     >> iunit(mcR_,GeV)  >> iunit(mbR_,GeV)  >> iunit(mtR_,GeV)
     >> gluinoPhase_ >> iunit(MPlanck_,GeV);
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SusyBase,BSMModel>
describeHerwigSusyBase("Herwig::SusyBase", "HwSusy.so");

void SusyBase::Init() {

  static ClassDocumentation<SusyBase> documentation
    ("This is the base class for any SUSY model.",
     "SUSY spectrum files follow the Les Houches accord"
     " \\cite{Skands:2003cj,Allanach:2008qq}.",
     " %\\cite{Skands:2003cj}\n"
     "\\bibitem{Skands:2003cj}\n"
     "  P.~Skands {\\it et al.},\n"
     "   ``SUSY Les Houches accord: Interfacing SUSY spectrum calculators, decay\n"
     "  %packages, and event generators,''\n"
     "  JHEP {\\bf 0407}, 036 (2004)\n"
     "  [arXiv:hep-ph/0311123].\n"
     "  %%CITATION = JHEPA,0407,036;%%\n"
     "%\\cite{Allanach:2008qq}\n"
     "\\bibitem{Allanach:2008qq}\n"
     "  B.~Allanach {\\it et al.},\n"
     "  %``SUSY Les Houches Accord 2,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 180}, 8 (2009)\n"
     "  [arXiv:0801.0045 [hep-ph]].\n"
     "  %%CITATION = CPHCB,180,8;%%\n"
     );


  static Reference<SusyBase,Helicity::AbstractVSSVertex> interfaceVertexWSS
    ("Vertex/WSFSF",
     "Reference to Susy W SF SF vertex",
     &SusyBase::WSFSFVertex_, false, false, true, false);
  
  static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexNFSF
    ("Vertex/NFSF",
     "Reference to the neutralino-fermion-sfermion vertex",
     &SusyBase::NFSFVertex_, false, false, true, false);

  static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexGFSF
    ("Vertex/GFSF",
     "Reference to the gluino-fermion-sfermion vertex",
     &SusyBase::GFSFVertex_, false, false, true, false);
  
  static Reference<SusyBase,Helicity::AbstractSSSVertex> interfaceVertexHSFSF
    ("Vertex/HSFSF",
     "Reference to the Higgs-fermion-sfermion vertex",
     &SusyBase::HSFSFVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexCFSF
   ("Vertex/CFSF",
      "Reference to the chargino-fermion-sfermion vertex",
      &SusyBase::CFSFVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVSSVertex> interfaceVertexGSFSF
   ("Vertex/GSFSF",
      "Reference to the gluon-sfermion-sfermion vertex",
      &SusyBase::GSFSFVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVVSSVertex> interfaceVertexGGSS
     ("Vertex/GGSQSQ",
      "Reference to the gluon-gluon-squark-squark vertex",
      &SusyBase::GGSQSQVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVVSSVertex> interfaceVertexWGSS
     ("Vertex/WGSQSQ",
      "Reference to the gauge boson-gluon-squark-squark vertex",
      &SusyBase::WGSQSQVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexGSGSG
     ("Vertex/GSGSG",
      "Reference to the gluon-gluino-gluino vertex",
      &SusyBase::GSGSGVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexNNZ
    ("Vertex/NNZ",
     "Reference to Z-~chi_i0-~chi_i0 vertex",
     &SusyBase::NNZVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexNNP
    ("Vertex/NNP",
     "Reference to photon-~chi_i0-~chi_i0 vertex",
     &SusyBase::NNPVertex_, false, false, true, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexGNG
    ("Vertex/GNG",
     "Reference to gluon-~chi_i0-gluino vertex",
     &SusyBase::GNGVertex_, false, false, true, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexCCZ
    ("Vertex/CCZ",
     "Reference to ~chi_i+-~chi_i-Z vertex",
     &SusyBase::CCZVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexCNW
    ("Vertex/CNW",
     "Reference to ~chi_i+-chi_i0-W vertex",
     &SusyBase::CNWVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexGOGOH
   ("Vertex/GOGOH",
    "Reference to the gaugino-gaugino-higgs vertex",
    &SusyBase::GOGOHVertex_, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVSSVertex> interfaceVertexWHH
    ("Vertex/SSWHH",
     "Reference to Susy WHHVertex",
     &SusyBase::WHHVertex_, false, false, true, false);

  static Reference<SusyBase,AbstractFFSVertex> interfaceVertexNCT
    ("Vertex/NCT",
     "Vertex for the flavour violating coupling of the top squark "
     "to the neutralino and charm quark.",
     &SusyBase::NCTVertex_, false, false, true, true, false);

  static Reference<SusyBase,AbstractRFSVertex> interfaceVertexGVNH
    ("Vertex/GVNH",
     "Vertex for the interfaction of the gravitino-neutralino"
     " and Higgs bosons",
     &SusyBase::GVNHVertex_, false, false, true, true, false);

  static Reference<SusyBase,AbstractRFVVertex> interfaceVertexGVNV
    ("Vertex/GVNV",
     "Vertex for the interfaction of the gravitino-neutralino"
     " and vector bosons",
     &SusyBase::GVNVVertex_, false, false, true, true, false);

  static Reference<SusyBase,AbstractRFSVertex> interfaceVertexGVFS
    ("Vertex/GVFS",
     "Vertex for the interfaction of the gravitino-fermion"
     " and sfermion",
     &SusyBase::GVFSVertex_, false, false, true, true, false);

  static Parameter<SusyBase,Energy> interfaceMPlanck
    ("MPlanck",
     "The Planck mass for GMSB models",
     &SusyBase::MPlanck_, GeV, 2.4e18*GeV, 1.e16*GeV, 1.e20*GeV,
     false, false, Interface::limited);

  static Switch<SusyBase,bool> interfaceMajoranaNeutrinos
    ("MajoranaNeutrinos",
     "Whether or not the neutrinos should be treated as Majorana particles",
     &SusyBase::majoranaNeutrinos_, false, false, false);
  static SwitchOption interfaceMajoranaNeutrinosYes
    (interfaceMajoranaNeutrinos,
     "Yes",
     "Neutrinos are Majorana",
     true);
  static SwitchOption interfaceMajoranaNeutrinosNo
    (interfaceMajoranaNeutrinos,
     "No",
     "Neutrinos are Dirac fermions",
     false);

}

void SusyBase::readSetup(istream & is) {
  string filename = dynamic_ptr_cast<istringstream*>(&is)->str();
  if(majoranaNeutrinos()) {
    idMap().insert(make_pair(12,17));
    idMap().insert(make_pair(14,18));
    idMap().insert(make_pair(16,19));
  }
  filename = StringUtils::stripws(filename);
  if(readFile_)
    throw SetupException() 
      << "A second SLHA file " << filename << " has been opened."
      << "This is probably unintended and as it can cause crashes"
      << " and other unpredictable behaviour it is not allowed."
      << Exception::runerror;
  CFileLineReader cfile;
  cfile.open(filename);
  if( !cfile ) throw SetupException() 
		 << "SusyBase::readSetup - An error occurred in opening the "
		 << "spectrum file \"" << filename << "\". A SUSY model cannot be "
		 << "run without this."
		 << Exception::runerror;
  useMe();
  // read first line and check if this is a Les Houches event file
  cfile.readline();
  bool lesHouches = cfile.find("<LesHouchesEvents");
  bool reading = !lesHouches;
  if(lesHouches) cfile.readline();
  //function pointer for putting all characters to lower case.
  int (*pf)(int) = tolower;
  do {
    string line = cfile.getline();
    // check for start of slha block in SLHA files
    if(lesHouches && !reading) {
      if(line.find("<slha")==0) reading = true;
      if(!cfile.readline()) break;
      continue;
    }
    // ignore comment lines
    if(line[0] == '#') {
      if(!cfile.readline()) break;
      continue;
    }
    // make everything lower case
    transform(line.begin(), line.end(), line.begin(), pf);
    // start of a block
    if(line.find("block") == 0) {
      string name = StringUtils::car(StringUtils::cdr(line), " #");
      name = StringUtils::stripws(name);
      // mixing matrix
      if((name.find("mix")  != string::npos && 
	  name.find("hmix") != 0)) {
	unsigned int row(0),col(0);
	MixingVector vals = readMatrix(cfile,row,col);
	mixings_[name] = make_pair(make_pair(row,col),vals);
      }
      else if(name.find("au") == 0 || name.find("ad") == 0 ||
	      name.find("ae") == 0 ) {
	string test = StringUtils::car(line, "#");
	while (test.find("=")!= string::npos) {
	  test = StringUtils::cdr(test, "=");
	}
	istringstream is(test);
	double scale;
	is >> scale;
	unsigned int row(0),col(0);
	MixingVector vals = readMatrix(cfile,row,col);
	if(scale>1e10) continue;
	mixings_[name] = make_pair(make_pair(row,col),vals);
      }
      else if( name == "spinfo" ) {
	readBlock(cfile,name,line,true);
      }
      else if( name.find("info") == string::npos) {
	readBlock(cfile,name,line,false);
      }
      else {
	if(!cfile.readline()) break;
      }
      continue;
    }
    else if( lesHouches && line.find("</slha") == 0 ) {
      break;
    }
    if(!cfile.readline()) break;
  }
  while(true);
  // extract the relevant parameters
  extractParameters();
  // create the mixing matrices we need
  createMixingMatrices();
  // set the masses, this has to be done after the 
  // mixing matrices have been created
  resetRepositoryMasses();
  // have now read the file
  if(decayFile()=="") decayFile(filename);
  readFile_=true;
}

void SusyBase::readBlock(CFileLineReader & cfile,string name,string linein,
			 bool stringBlock) {
  if(!cfile)
    throw SetupException() 
				    << "SusyBase::readBlock() - The input stream is in a bad state"
				    << Exception::runerror;
  // storage or the parameters
  string test = StringUtils::car(linein, "#");
  ParamMap  storeParam;
  StringMap storeString;
  bool set = true;
  // special for the alpha block
  if(name.find("alpha") == 0 ) {
    double alpha;
    cfile.readline();
    string line = cfile.getline();
    istringstream iss(line);
    iss >> alpha;
    storeParam.insert(make_pair(1,alpha));
    parameters_[name]=storeParam;
    return;
  }
  // extract the scale from the block if present
  if(test.find("=")!= string::npos) { 
    while(test.find("=")!=string::npos)
      test= StringUtils::cdr(test,"=");
    istringstream is(test);
    double scale;
    is >> scale;
    // only store the lowest scale block
    if(parameters_.find(name)!=parameters_.end()) {
      set = scale < parameters_[name][-1];
    }
    else {
      storeParam.insert(make_pair(-1,scale));
    }
  }
  while(cfile.readline()) {
    string line = cfile.getline();
    // skip comments
    if(line[0] == '#') continue;
    // reached the end
    if( line[0] == 'B' || line[0] == 'b' ||
	line[0] == 'D' || line[0] == 'd' ||
	line[0] == '<' ) {
      cfile.resetline();
      break;
    }
    istringstream is(line);
    long index;
    if(!stringBlock) {
      double value;
      if(name.find("rvlam")!= string::npos||
	 name=="rvt" || name=="rvtp" || name=="rvtpp") {
	int i,j,k;
	is >> i >> j >> k >> value;
	index = i*100+j*10+k;
      }
      else {
	is >> index >> value;
      }
      storeParam.insert(make_pair(index, value));
    }
    else {
      string value;
      is >> index >> value;
      storeString.insert(make_pair(index,value));
    }
  }
  if(set) {
    if(stringBlock) info_      [name] = storeString;
    else            parameters_[name] = storeParam ;
  }
}

const MixingVector
SusyBase::readMatrix(CFileLineReader & cfile, 
		     unsigned int & row, unsigned int & col) {
  if(!cfile)
    throw SetupException() 
      << "SusyBase::readMatrix() - The input stream is in a bad state."
      << Exception::runerror;
  unsigned int rowmax(0), colmax(0);
  MixingVector values;
  while(cfile.readline()) {
    string line = cfile.getline();
    // skip comments
    if(line[0] == '#') continue;
    // reached the end
    if( line[0] == 'B' || line[0] == 'b' ||
	line[0] == 'D' || line[0] == 'd' ||
	line[0] == '<' ) {
      cfile.resetline();
      break;
    }
    istringstream is(line);
    unsigned int index1, index2;
    double real(0.), imag(0.);   
    line = StringUtils::stripws(line);
    if(!line.empty()) {
      is >> index1 >> index2 >> real >> imag;
      values.push_back(MixingElement(index1,index2,Complex(real, imag)));
      if(index1 > rowmax) rowmax = index1;
      if(index2 > colmax) colmax = index2;
    }
  }
  col=colmax;
  row=rowmax;
  return values;
}

void SusyBase::createMixingMatrix(MixingMatrixPtr & matrix,
				  string name, const MixingVector & values,
				  MatrixSize size) {
  matrix = new_ptr(MixingMatrix(size.first,size.second));
  for(unsigned int ix=0; ix < values.size(); ++ix)
    (*matrix)(values[ix].row-1,values[ix].col-1) = values[ix].value;
  // test against stupid mixing matrices  
  for(unsigned int ix=0;ix<matrix->size().first;++ix) {
    Complex sum(0.);
    for(unsigned int iy=0;iy<matrix->size().second;++iy) {
      sum += norm((*matrix)(ix,iy));
    }
    if(abs(sum-1.)>1e-4) {
      cerr << "The sum of the mod squares of row " << ix+1
	   << " of the " << name << " block does not sum to 1. \n"
	   << "sum = " << sum.real() << ". We strongly suggest you check your SLHA file.\n";
    }
  }

  vector<long> ids;
  if(name == "nmix") {
    ids.resize(4);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035; 
  }
  else if(name == "nmnmix") {
    ids.resize(5);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035;
    ids[4] = 1000045;
  }
  else if(name == "rvnmix") {
    ids.resize(7);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035;
    if(!majoranaNeutrinos()) {
      ids[4] = 12; ids[5] = 14; ids[6] = 16;
    }
    else { 
      ids[4] = 17; ids[5] = 18; ids[6] = 19;
    }
  }
  else if(name == "rpvmix") {
    ids.resize(5);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035;
    ids[4] = 1000045;
  }
  else if(name == "umix" || name == "vmix") {
    ids.resize(2);
    ids[0] = 1000024; ids[1] = 1000037;
  }
  else if(name == "rvumix" || name == "rvvmix" ) {
    ids.resize(5);
    ids[3] = 1000024; ids[4] = 1000037; 
    ids[0] = -11; ids[1] = -13; ids[2] = -15;
  }
  else if(name == "stopmix") {
    ids.resize(2);
    ids[0] = 1000006; ids[1] = 2000006;
  }
  else if(name == "sbotmix") {
    ids.resize(2);
    ids[0] = 1000005; ids[1] = 2000005;
  }
  else if(name == "staumix") {
    ids.resize(2);
    ids[0] = 1000015; ids[1] = 2000015;
  }
  else if(name == "nmhmix") {
    ids.resize(3);
    ids[0] = 25; ids[1] = 35; ids[2] = 45;
  }
  else if(name == "nmamix") {
    ids.resize(2);
    ids[0] = 36; ids[1] = 46;
  }
  else if(name == "rvamix") {
    // some programs, i.e. spheno include the goldstone even though the standard says they shouldn't
    // obeys standard
    if(size.first==4) {
      ids.resize(4);
    }
    // includes goldstone
    else if(size.first==5) {
      ids.resize(5);
      ids[4] = 0;
    }
    else {
      throw Exception() << "SusyBase::createMixingMatrix() "
			<< "rvamix matrix must have either 4 rows (obeying standard) "
			<< "or 5 including the goldstone, not " << size.first << Exception::runerror;
    }
    ids[0] = 36     ; ids[1] = 1000017;
    ids[2] = 1000018; ids[3] = 1000019;
  }
  else if(name == "rvhmix") {
    ids.resize(5);
    ids[0] = 25; ids[1] = 35;
    ids[2] = 1000012; ids[3] = 1000014; ids[4] = 1000016;
  }
  else if(name == "rvlmix") {
    // some programs, i.e. spheno include the goldstone even though the standard says they shouldn't
    // obeys standard
    if(size.first==7) {
      ids.resize(7);
    }
    // includes goldstone
    else if(size.first==8) {
      ids.resize(8);
      ids[7] = 0;
    }
    else {
      throw Exception() << "SusyBase::createMixingMatrix() "
			<< "rvlmix matrix must have either 7 rows (obeying standard) "
			<< "or 8 including the goldstone, not " << size.first << Exception::runerror;
    }
    ids[0] = 37;
    ids[1] = -1000011;
    ids[2] = -1000013;
    ids[3] = -1000015;
    ids[4] = -2000011;
    ids[5] = -2000013;
    ids[6] = -2000015;
  }
  else if(name == "usqmix" ) {
    ids.resize(6);
    ids[0] = 1000002; ids[1] = 1000004; ids[2] = 1000006;
    ids[3] = 2000002; ids[4] = 2000004; ids[5] = 2000006;
  }
  else if(name == "dsqmix" ) {
    ids.resize(6);
    ids[0] = 1000001; ids[1] = 1000003; ids[2] = 1000005;
    ids[3] = 2000001; ids[4] = 2000003; ids[5] = 2000005;
  }
  else {
    throw SetupException() << "SusyBase::createMixingMatrix() "
			   << "cannot find correct title for mixing matrix "
			   << name << Exception::runerror;
  }
  matrix->setIds(ids);
}

void SusyBase::resetRepositoryMasses() {
  map<string,ParamMap>::const_iterator fit=parameters_.find("mass");
  if(fit==parameters_.end()) 
    throw Exception() << "BLOCK MASS not found in input file"
		      << " can't set masses of SUSY particles"
		      << Exception::runerror;
  ParamMap theMasses = fit->second;
  for(ParamMap::iterator it = theMasses.begin(); it != theMasses.end(); 
      ++it) {
    long id = it->first;
    map<long,long>::const_iterator dit = idMap().find(id);
    if(dit!=idMap().end()) id = dit->second;
    double mass = it->second;
    //a negative mass requires an adjustment to the 
    //associated mixing matrix by a factor of i
    if(mass < 0.0) adjustMixingMatrix(id);
    PDPtr part = getParticleData(id);
    if(!part) throw SetupException() 
      << "SusyBase::resetRepositoryMasses() - Particle with PDG code " << id  
      << " not found." << Exception::warning;
    if(abs(id)<=5||abs(id)==23||abs(id)==24||
       (abs(id)>=11&&abs(id)<=16))
      cerr << "SusyBase::resetRepositoryMasses() Resetting mass of " 
	   << part->PDGName() << " using SLHA "
	   << "file,\nthis can affect parts of the Standard Model simulation and"
	   << " is strongly discouraged.\n";
    // reset the masses
    resetMass(it->first,it->second*GeV,part);

    // switch on gravitino interactions?
    gravitino_ |= id== ParticleID::SUSY_Gravitino;
  }
  theMasses.clear();
}

void SusyBase::adjustMixingMatrix(long id) {
  //get correct mixing matrix
  switch(id) {
  case 1000021 :
    gluinoPhase_ = Complex(0.,1.);
    break;
  case 1000022 :
  case 1000023 :
  case 1000025 :
  case 1000035 : 
  case 1000045 :
  case 12 : case 17 : 
  case 14 : case 18 : 
  case 16 : case 19 :
    if(NMix_) {
      if(id>20||(id<=16&&NMix_->size().first>4))
	 NMix_->adjustPhase(id);
    }
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The neutralino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    break;
  case 1000024 :
  case 1000037 : 
    if(UMix_)
      UMix_->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The U-Type chargino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    if(VMix_) VMix_->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The V-Type chargino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    break;
  default : 
    throw SetupException() 
      << "SusyBase::adjustMixingMatrix - Trying to adjust mixing matrix "
      << "phase for a particle that does not have a mixing matrix "
      << "associated with it. " << id << " must have a negative mass in "
      << "the spectrum file, this should only occur for particles that mix."
      << Exception::runerror;
  }
}

void SusyBase::createMixingMatrices() {
  map<string,pair<MatrixSize, MixingVector > >::const_iterator it;
  for(it=mixings_.begin();it!=mixings_.end();++it) {
    string name=it->first;
    // create the gaugino mixing matrices
    if(name == "nmix" || name == "nmnmix" || name == "rvnmix" ) {
      createMixingMatrix(NMix_,name,it->second.second,it->second.first);
    }
    else if (name == "umix" || name == "rvumix" ) {
      createMixingMatrix(UMix_,name,it->second.second,it->second.first);
    }
    else if (name == "vmix" || name == "rvvmix" ) {
      createMixingMatrix(VMix_,name,it->second.second,it->second.first);
    }
  }
}

void SusyBase::extractParameters(bool checkmodel) {
  map<string,ParamMap>::const_iterator pit;
  ParamMap::const_iterator it;
  // try and get tan beta from hmin and extpar first
  // extract tan beta
  tanBeta_ = -1.;
  if(tanBeta_<0.) {
    pit=parameters_.find("hmix");
    if(pit!=parameters_.end()) {
      it = pit->second.find(2);
      if(it!=pit->second.end()) tanBeta_ = it->second;
    }
  }
  if(tanBeta_<0.) {
    pit=parameters_.find("extpar");
    if(pit!=parameters_.end()) {
      it = pit->second.find(25);
      if(it!=pit->second.end()) tanBeta_ = it->second;
    }
  }
  // otherwise from minpar
  if(tanBeta_<0.) {
    pit=parameters_.find("minpar");
    if(pit!=parameters_.end()) { 
      it = pit->second.find(3);
      if(it!=pit->second.end()) tanBeta_ = it->second;
    }
  }
  if(tanBeta_<0.) 
    throw Exception() << "SusyBase::extractParameters() "
		      << "Can't find tan beta in BLOCK MINPAR"
		      << " or BLOCK EXTPAR " << Exception::runerror;
  if(tanBeta_==0.)
    throw Exception() << "Tan(beta) = 0 in SusyBase::extractParameters()"
		      << Exception::runerror;
  // extract parameters from hmix
  pit=parameters_.find("hmix");
  if(pit==parameters_.end()) {
    if(generator())
      generator()->logWarning(Exception("SusyBase::extractParameters() BLOCK HMIX not found setting mu to zero\n",
					Exception::warning));
    else
      cerr << "SusyBase::extractParameters() BLOCK HMIX not found setting mu to zero\n";
     mu_=ZERO;
  }
  else {
    mu_=findValue(pit,1,"HMIX","mu")*GeV;
  }
  pit = parameters_.find("msoft");
  if( pit == parameters_.end() )
    throw Exception() << "SusyBase::extractParameters() "
		      << "BLOCK MSOFT not found in " 
		      << "SusyBase::extractParameters()"
		      << Exception::runerror;
  M1_    = findValue(pit,1 ,"MSOFT","M_1"   )*GeV;
  M2_    = findValue(pit,2 ,"MSOFT","M_2"   )*GeV;
  M3_    = findValue(pit,3 ,"MSOFT","M_3"   )*GeV;
  mH12_  = findValue(pit,21,"MSOFT","m_H1^2")*GeV2;
  mH22_  = findValue(pit,22,"MSOFT","m_H2^2")*GeV2;
  meL_   = findValue(pit,31,"MSOFT","M_eL"  )*GeV;
  mmuL_  = findValue(pit,32,"MSOFT","M_muL" )*GeV;
  mtauL_ = findValue(pit,33,"MSOFT","M_tauL")*GeV; 
  meR_   = findValue(pit,34,"MSOFT","M_eR"  )*GeV;
  mmuR_  = findValue(pit,35,"MSOFT","M_muR" )*GeV;
  mtauR_ = findValue(pit,36,"MSOFT","M_tauR")*GeV; 
  mq1L_  = findValue(pit,41,"MSOFT","M_q1L" )*GeV;
  mq2L_  = findValue(pit,42,"MSOFT","M_q2L" )*GeV;
  mq3L_  = findValue(pit,43,"MSOFT","M_q3L" )*GeV; 
  muR_   = findValue(pit,44,"MSOFT","M_uR"  )*GeV;
  mcR_   = findValue(pit,45,"MSOFT","M_cR"  )*GeV;
  mtR_   = findValue(pit,46,"MSOFT","M_tR"  )*GeV;
  mdR_   = findValue(pit,47,"MSOFT","M_dR"  )*GeV;
  msR_   = findValue(pit,48,"MSOFT","M_sR"  )*GeV;
  mbR_   = findValue(pit,49,"MSOFT","M_bR"  )*GeV;
  if(checkmodel) {
    throw Exception() << "The SusyBase class should not be used as a "
		      << "Model class, use one of the models which inherit"
		      << " from it" << Exception::runerror;
  }
}
