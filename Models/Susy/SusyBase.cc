// -*- C++ -*-
//
// SusyBase.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SusyBase class.
//

#include "SusyBase.h"
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
#include "ThePEG/PDT/MassGenerator.h"
#include "ThePEG/PDT/WidthGenerator.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;

SusyBase::SusyBase() : _readFile(false), _topModesFromFile(false),
		       _tolerance(1e-6),
		       _tanbeta(0), _mu(ZERO), 
		       theMone(ZERO), theMtwo(ZERO),
		       theMthree(ZERO) {}

IBPtr SusyBase::clone() const {
  return new_ptr(*this);
}

IBPtr SusyBase::fullclone() const {
  return new_ptr(*this);
}

void SusyBase::doinit() {
  addVertex(theWSFSFVertex);
  addVertex(theNFSFVertex);
  addVertex(theGFSFVertex);
  addVertex(vertexHSS());
  addVertex(theCFSFVertex);
  addVertex(theGSFSFVertex);
  addVertex(theGGSQSQVertex);
  addVertex(theGSGSGVertex);
  addVertex(theNNZVertex);
  addVertex(theNNPVertex);
  addVertex(theGNGVertex);
  addVertex(theCCZVertex);
  addVertex(theCNWVertex);
  addVertex(vertexGOGOH());
  addVertex(vertexWHH());
  addVertex(vertexHHH());
  StandardModel::doinit();
}

void SusyBase::persistentOutput(PersistentOStream & os) const {
  os << _readFile << _topModesFromFile 
     << theNMix << theUMix << theVMix << theWSFSFVertex 
     << theNFSFVertex << theGFSFVertex << theHSFSFVertex << theCFSFVertex 
     << theGSFSFVertex << theGGSQSQVertex << theGSGSGVertex 
     << theNNZVertex << theNNPVertex << theCCZVertex << theCNWVertex 
     << theGOGOHVertex << theWHHVertex 
     << theGNGVertex << theHHHVertex << _tanbeta << ounit(_mu,GeV) 
     << ounit(theMone,GeV) << ounit(theMtwo,GeV) << ounit(theMthree,GeV)
     << _tolerance;
}

void SusyBase::persistentInput(PersistentIStream & is, int) {
  is >> _readFile >> _topModesFromFile
     >> theNMix >> theUMix >> theVMix >> theWSFSFVertex 
     >> theNFSFVertex >> theGFSFVertex >> theHSFSFVertex >> theCFSFVertex 
     >> theGSFSFVertex >> theGGSQSQVertex >> theGSGSGVertex 
     >> theNNZVertex >> theNNPVertex >> theCCZVertex >> theCNWVertex
     >> theGOGOHVertex >> theWHHVertex
     >> theGNGVertex >> theHHHVertex >> _tanbeta >> iunit(_mu,GeV) 
     >> iunit(theMone,GeV) >> iunit(theMtwo,GeV) >> iunit(theMthree,GeV)
     >> _tolerance;
}

ClassDescription<SusyBase> SusyBase::initSusyBase;
// Definition of the static class description member.

void SusyBase::Init() {

  static ClassDocumentation<SusyBase> documentation
    ("This is the base class for any SUSY model.",
     "SUSY spectrum files follow the Les Houches accord \\cite{Skands:2003cj,Allanach:2008qq}.",
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

  static Switch<SusyBase,bool> interfaceTopModes
    ("TopModes",
     "Whether ro use the Herwig++ SM top decays or those from the SLHA file",
     &SusyBase::_topModesFromFile, false, false, false);
  static SwitchOption interfaceTopModesFile
    (interfaceTopModes,
     "File",
     "Take the modes from the files",
     true);
  static SwitchOption interfaceTopModesHerwig
    (interfaceTopModes,
     "Herwig",
     "Use the SM ones",
     false);

  static Reference<SusyBase,Helicity::AbstractVSSVertex> interfaceVertexWSS
    ("Vertex/WSFSF",
     "Reference to Susy W SF SF vertex",
     &SusyBase::theWSFSFVertex, false, false, true, false);
  
  static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexNFSF
    ("Vertex/NFSF",
     "Reference to the neutralino-fermion-sfermion vertex",
     &SusyBase::theNFSFVertex, false, false, true, false);

  static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexGFSF
    ("Vertex/GFSF",
     "Reference to the gluino-fermion-sfermion vertex",
     &SusyBase::theGFSFVertex, false, false, true, false);
  
  static Reference<SusyBase,Helicity::AbstractSSSVertex> interfaceVertexHSFSF
    ("Vertex/HSFSF",
     "Reference to the Higgs-fermion-sfermion vertex",
     &SusyBase::theHSFSFVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexCFSF
   ("Vertex/CFSF",
      "Reference to the chargino-fermion-sfermion vertex",
      &SusyBase::theCFSFVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVSSVertex> interfaceVertexGSFSF
   ("Vertex/GSFSF",
      "Reference to the gluon-sfermion-sfermion vertex",
      &SusyBase::theGSFSFVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVVSSVertex> interfaceVertexGGSS
     ("Vertex/GGSQSQ",
      "Reference to the gluon-gluon-squark-squark vertex",
      &SusyBase::theGGSQSQVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexGSGSG
     ("Vertex/GSGSG",
      "Reference to the gluon-gluino-gluino vertex",
      &SusyBase::theGSGSGVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexNNZ
    ("Vertex/NNZ",
     "Reference to Z-~chi_i0-~chi_i0 vertex",
     &SusyBase::theNNZVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexNNP
    ("Vertex/NNP",
     "Reference to photon-~chi_i0-~chi_i0 vertex",
     &SusyBase::theNNPVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexGNG
    ("Vertex/GNG",
     "Reference to gluon-~chi_i0-gluino vertex",
     &SusyBase::theGNGVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexCCZ
    ("Vertex/CCZ",
     "Reference to ~chi_i+-~chi_i-Z vertex",
     &SusyBase::theCCZVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFVVertex> interfaceVertexCNW
    ("Vertex/CNW",
     "Reference to ~chi_i+-chi_i0-W vertex",
     &SusyBase::theCNWVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractFFSVertex> interfaceVertexGOGOH
   ("Vertex/GOGOH",
    "Reference to the gaugino-gaugino-higgs vertex",
    &SusyBase::theGOGOHVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractVSSVertex> interfaceVertexWHH
    ("Vertex/SSWHH",
     "Reference to Susy WHHVertex",
     &SusyBase::theWHHVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::AbstractSSSVertex> interfaceVertexHHH
    ("Vertex/HHH",
     "Triple higgs coupling",
     &SusyBase::theHHHVertex, false, false, true, false);

  static Parameter<SusyBase,double> interfaceBRTolerance
    ("BRTolerance",
     "Tolerance for the sum of branching ratios to be difference from one.",
     &SusyBase::_tolerance, 1e-6, 1e-8, 0.01,
     false, false, Interface::limited);

}

void SusyBase::readSetup(istream & is) {
  string filename = dynamic_ptr_cast<istringstream*>(&is)->str();
  if(_readFile)
    throw SetupException() 
      << "A second SLHA file " << filename << " has been opened."
      << "This is probably unintended and as it can cause crashes"
      << " and other unpredictable behaviour it is not allowed."
      << Exception::runerror;
  ifstream file(filename.c_str());
  if( !file ) throw SetupException() 
    << "SusyBase::readSetup - An error occurred in opening the "
    << "spectrum file \"" << filename << "\". A SUSY model cannot be "
    << "run without this."
    << Exception::runerror; 

  useMe();
  //Before reading the spectrum/decay files the SM higgs 
  //decay modes, mass and width generators need to be turned off.
  PDPtr h0 = getParticleData(ParticleID::h0);
  h0->widthGenerator(WidthGeneratorPtr());
  h0->massGenerator(MassGenPtr());
  DecaySet::const_iterator dit = h0->decayModes().begin();
  DecaySet::const_iterator dend = h0->decayModes().end();
  for( ; dit != dend; ++dit ) {
    const InterfaceBase * ifb = 
      BaseRepository::FindInterface(*dit, "BranchingRatio");
    ifb->exec(**dit, "set", "0.0");
    ifb = BaseRepository::FindInterface(*dit, "OnOff");
    ifb->exec(**dit, "set", "Off");
  }
  // if taking the top modes from the file
  // delete the SM stuff
  if(_topModesFromFile) {
    PDPtr top = getParticleData(ParticleID::t);
    top->widthGenerator(WidthGeneratorPtr());
    top->massGenerator(MassGenPtr());
    DecaySet::const_iterator dit = top->decayModes().begin();
    DecaySet::const_iterator dend = top->decayModes().end();
    for( ; dit != dend; ++dit ) {
      const InterfaceBase * ifb = 
	BaseRepository::FindInterface(*dit, "BranchingRatio");
      ifb->exec(**dit, "set", "0.0");
      ifb = BaseRepository::FindInterface(*dit, "OnOff");
      ifb->exec(**dit, "set", "Off");
    }
    
  }
  string line;
  //function pointer for putting all characters to lower case.
  int (*pf)(int) = tolower;
   //Start to read in the file
   while(getline(file, line)) {
     // ignore comment lines
     if(line[0] == '#') continue;
     // make everything lower case
     transform(line.begin(), line.end(), line.begin(), pf);
     // start of a block
     if(line.find("block") == 0) {
       string name = StringUtils::car(StringUtils::cdr(line), " #");
       // mixing matrix
       if((name.find("mix")  != string::npos && 
	   name.find("hmix") != 0)) {
	 unsigned int row(0),col(0);
	 MixingVector vals = readMatrix(file,row,col);
	 _mixings[name] = make_pair(make_pair(row,col),vals);
       }
       else if(name.find("au") == 0||name.find("ad") == 0||
	       name.find("ae") == 0 ) {
	 string test = StringUtils::car(line, "#");
	 while (test.find("=")!= string::npos) {
	   test = StringUtils::cdr(test, "=");
	 }
	 istringstream is(test);
	 double scale;
	 is >> scale;
	 if(scale>1e10) continue;
	 unsigned int row(0),col(0);
	 MixingVector vals = readMatrix(file,row,col);
	 _mixings[name] = make_pair(make_pair(row,col),vals);
       }
       else if( name.find("info") == string::npos)
	 readBlock(file,name);
     }
     // decays
     else if( line.find("decay") == 0 )
       readDecay(file, line);
   }
   // extract the relevant parameters
   extractParameters();
   // create the mixing matrices we need
   createMixingMatrices();
   // set the masses, this has to be done after the 
   // mixing matrices have been created
   resetRepositoryMasses();
   // have now read the file
   _readFile=true;
}

void SusyBase::readBlock(ifstream & ifs,string name) {
  if(!ifs)
    throw SetupException() 
      << "SusyBase::readBlock() - The input stream is in a bad state"
      << Exception::runerror;
  string line;
  ParamMap store;
  // special for the alpha block
  if(name.find("alpha") == 0 ) {
    double alpha;
    getline(ifs, line);
    istringstream iss(line);
    iss >> alpha;
    store.insert(make_pair(1,alpha));
  }
  else {
    while(getline(ifs, line)) {
      if(line[0] == '#') {
	if( ifs.peek() == 'D' || ifs.peek() == 'B' ||
	    ifs.peek() == 'd' || ifs.peek() == 'b' ) break;
	else continue;
      }
      istringstream is(line);
      long index;
      double value;   
      is >> index >> value;
      store.insert(make_pair(index, value));
      if(ifs.peek() == 'B' || ifs.peek() == 'D' ||
	 ifs.peek() == 'b' || ifs.peek() == 'd' ||
	 ifs.peek() == '#') break;
    }
  }
  _parameters[name]=store;
}

void SusyBase::readDecay(ifstream & ifs, 
			 string decay) const{
  if(!ifs)
    throw SetupException() 
      <<"SusyBase::readDecay - The input stream is in a bad state";
  //if the particle is stable then next line will be the start of a 
  //new decaymode or a new block so just return if this is found
  if( ifs.peek() == 'D' || ifs.peek() == 'B' ||
      ifs.peek() == 'd' || ifs.peek() == 'b' ) return;
  istringstream iss(decay);
  string dummy;
  long parent(0);
  Energy width(ZERO);
  iss >> dummy >> parent >> iunit(width, GeV);
  PDPtr inpart = getParticleData(parent);
  if(!_topModesFromFile&&abs(parent)==ParticleID::t) return;
  if(!inpart)  {
    throw SetupException() 
    << "SusyBase::readDecay() - A ParticleData object with the PDG code "
    << parent << " does not exist. " << Exception::runerror;
    return;
  }
  inpart->width(width);
  if( width > ZERO ) inpart->cTau(hbarc/width);
  inpart->widthCut(5.*width);
  string prefix("decaymode " + inpart->name() + "->"), tag(""),line("");
  double brsum(0.);
  unsigned int nmode = 0;
  while(getline(ifs, line)) {
    if(line[0] == '#') {
      if( ifs.peek() == 'D' || ifs.peek() == 'B' ||
	  ifs.peek() == 'd' || ifs.peek() == 'b' ) break;
      else continue;
    }
    istringstream is(line);
    double brat(0.);
    unsigned int nda(0),npr(0);
    is >> brat >> nda;
    brsum += brat;
    ++nmode;
    while( true ) {
      long t;
      is >> t;
      if( is.fail() ) break; 
      if( t == abs(parent) ) {
	throw SetupException() 
	  << "An error occurred while read a decay of the " 
	  << inpart->PDGName() << ". One of its products has the same PDG code "
	  << "as the parent particle. Please check the SLHA file.\n"
	  << Exception::runerror;
      }
      tcPDPtr p = getParticleData(t);
      if( !p ) {
	throw SetupException()
	  << "SusyBase::readDecay() - An unknown PDG code has been encounterd "
	  << "while reading a decay mode. ID: " << t
	  << Exception::runerror;
      }
      else {
	++npr;
	tag += p->name() + ",";
      }
    }
    if( npr != nda ) {
      throw SetupException()
	<< "SusyBase::readDecay - While reading a decay of the " 
	<< inpart->PDGName() << " from an SLHA file, an inconsistency "
	<< "between the number of decay products and the value in "
	<< "the 'NDA' column was found. Please check if the spectrum "
	<< "file is correct.\n"
	<< Exception::warning;
    }
    if( npr > 1 ) {
      tag.replace(tag.size() - 1, 1, ";");
      createDecayMode(prefix + tag, brat);
      tag.clear();
    }
    if( ifs.peek() == 'D' || ifs.peek() == 'B' ||
	ifs.peek() == 'd' || ifs.peek() == 'b' ) break;
  }
  if( abs(brsum - 1.) > _tolerance && nmode!=0 ) {
    cerr << "Warning: The total branching ratio for " << inpart->PDGName()
	 << " from the spectrum file does not sum to 1. The branching fractions"
	 << " will be rescaled.\n";
    cerr << setprecision(13) << abs(brsum - 1.) << "\n";
  }
}

const MixingVector
SusyBase::readMatrix(ifstream & ifs, 
		     unsigned int & row, unsigned int & col) {
  if(!ifs)
    throw SetupException() 
      << "SusyBase::readMatrix() - The input stream is in a bad state."
      << Exception::runerror;
  string line;
  unsigned int rowmax(0), colmax(0);
  MixingVector values;
  while(getline(ifs, line)) {
    if(line[0] == '#') {
      if( ifs.peek() == 'D' || ifs.peek() == 'B' ||
	  ifs.peek() == 'd' || ifs.peek() == 'b' ) break;
      else continue;
    }
    istringstream is(line);
    unsigned int index1, index2;
    double real(0.), imag(0.);   
    is >> index1 >> index2 >> real >> imag;
    values.push_back(MixingElement(index1,index2,Complex(real, imag)));
    if(index1 > rowmax) rowmax = index1; 
    if(index2 > colmax) colmax = index2; 
    if( ifs.peek() == 'B' || ifs.peek() == 'D' ||
        ifs.peek() == 'b' || ifs.peek() == 'd' ) break;
  }
  col=colmax;
  row=rowmax;
  return values;
}

void SusyBase::createDecayMode(string tag, double brat) const {
  ostringstream cmd;
  cmd << tag << string(" ") 
      << setprecision(13) << brat << string(" 1 /Herwig/Decays/Mambo");
  Repository::exec(cmd.str(), cerr);
}

void SusyBase::createMixingMatrix(MixingMatrixPtr & matrix,
				  string name, const MixingVector & values,
				  MatrixSize size) {
  matrix = new_ptr(MixingMatrix(size.first,size.second));
  for(unsigned int ix=0; ix < values.size(); ++ix)
    (*matrix)(values[ix].row-1,values[ix].col-1) = values[ix].value;
  
  if(name == "nmix") {
    vector<long> ids(4);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035; 
    matrix->setIds(ids);
  }
  else if(name == "nmnmix") {
    vector<long> ids(5);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035;
    ids[4] = 1000045; 
    matrix->setIds(ids);
  }
  else if(name == "umix") {
    vector<long> ids(2);
    ids[0] = 1000024; ids[1] = 1000037; 
    matrix->setIds(ids);
  }
  else if(name == "vmix") {
    vector<long> ids(2);
    ids[0] = 1000024; ids[1] = 1000037; 
    matrix->setIds(ids);
  }
  else if(name == "stopmix") {
    vector<long> ids(2);
    ids[0] = 1000006; ids[1] = 2000006; 
    matrix->setIds(ids);
  }
  else if(name == "sbotmix") {
    vector<long> ids(2);
    ids[0] = 1000005; ids[1] = 2000005; 
    matrix->setIds(ids);
  }
  else if(name == "staumix") {
    vector<long> ids(2);
    ids[0] = 1000015; ids[1] = 2000015; 
    matrix->setIds(ids);
  }
  else if(name == "nmhmix") {
    vector<long> ids(3);
    ids[0] = 25; ids[1] = 35; ids[2] = 45;
    matrix->setIds(ids);
  }
  else if(name == "nmamix") {
    vector<long> ids(2);
    ids[0] = 36; ids[1] = 46;
    matrix->setIds(ids);
  }
  else
    throw SetupException() << "Cannot find correct title for mixing matrix "
			   << name << Exception::runerror;
}

void SusyBase::resetRepositoryMasses() {
  map<string,ParamMap>::const_iterator fit=_parameters.find("mass");
  if(fit==_parameters.end()) 
    throw Exception() << "BLOCK MASS not found in input file"
		      << " can't set masses of SUSY particles"
		      << Exception::runerror;
  ParamMap theMasses = fit->second;
  for(ParamMap::iterator it = theMasses.begin(); it != theMasses.end(); 
      ++it) {
    long id = it->first;
    double mass = it->second;
    //a negative mass requires an adjustment to the 
    //associated mixing matrix by a factor of i
    if(mass < 0.0) adjustMixingMatrix(id);
    PDPtr part = getParticleData(id);
    if(!part) throw SetupException() 
      << "SusyBase::resetRepositoryMasses() - Particle with PDG code " << id  
      << " not found." << Exception::warning;
    //Find interface nominal mass interface
    const InterfaceBase * ifb = BaseRepository::FindInterface(part, "NominalMass");
    ostringstream os;
    os << abs(it->second);
    ifb->exec(*part, "set", os.str());
  }
  theMasses.clear();
}

void SusyBase::adjustMixingMatrix(long id) {
  //get correct mixing matrix
  switch(id) {
  case 1000022 :
  case 1000023 :
  case 1000025 :
  case 1000035 : 
  case 1000045 : 
    if(theNMix)
      theNMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The neutralino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    break;
  case 1000024 :
  case 1000037 : 
    if(theUMix)
      theUMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The U-Type chargino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    if(theVMix)
      theVMix->adjustPhase(id);
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
  for(it=_mixings.begin();it!=_mixings.end();++it) {
    string name=it->first;
    // create the gaugino mixing matrices
    if(name == "nmix" || name == "nmnmix" ){
      createMixingMatrix(theNMix,name,it->second.second,it->second.first);
    }
    else if (name == "umix" ) {
      createMixingMatrix(theUMix,name,it->second.second,it->second.first);
    }
    else if (name == "vmix") {
      createMixingMatrix(theVMix,name,it->second.second,it->second.first);
    }
  }
}

void SusyBase::extractParameters(bool checkmodel) {
  map<string,ParamMap>::const_iterator pit;
  ParamMap::const_iterator it;
  // try and get tan beta from extpar first
  pit=_parameters.find("extpar");
  // extract tan beta
  _tanbeta = -1.;
  if(pit!=_parameters.end()) {
    it = pit->second.find(25);
    if(it!=pit->second.end()) _tanbeta = it->second;
  }
  // otherwise from minpar
  if(_tanbeta<0.) {
    pit=_parameters.find("minpar");
    if(pit!=_parameters.end()) { 
      it = pit->second.find(3);
      if(it!=pit->second.end()) _tanbeta = it->second;
    }
  }
  if(_tanbeta<0.) 
    throw Exception() << "Can't find tan beta in BLOCK MINPAR"
		      << " or BLOCK EXTPAR " << Exception::runerror;
  // extract parameters from hmix
  pit=_parameters.find("hmix");
  if(pit==_parameters.end()) {
    cerr << "BLOCK HMIX not found setting mu to zero\n";
    _mu=ZERO;
  }
  else {
    it = pit->second.find(1);
    if(it==pit->second.end()) {
      cerr << "mu not found in BLOCK HMIX setting to zero\n";
    }
    else {
      _mu=it->second*GeV;
    }
  }
  pit = _parameters.find("msoft");
  if( pit == _parameters.end() )
    throw Exception() << "BLOCK MSOFT not found in " 
		      << "SusyBase::extractParameters()"
		      << Exception::runerror;
  it = pit->second.find(1);
  theMone = it->second*GeV;
  it = pit->second.find(2);
  theMtwo = it->second*GeV;
  it = pit->second.find(3);
  theMthree = it->second*GeV;
  if(checkmodel) {
    throw Exception() << "The SusyBase class should not be used as a "
		      << "Model class, use one of the models which inherit"
		      << " from it" << Exception::runerror;
  }
}
