//A few standard headers
#include <iostream>
#include <sstream>
#include <cstring>
#include <fstream>
#include <cmath>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iomanip>

using namespace std;

double sqr(double x);

inline int nInt(double x) {
  int theCeiling=int(ceil(x));
  int theFloor=int(floor(x));
  if((theCeiling-theFloor!=1&&theCeiling!=theFloor)
     ||theCeiling-x>1.0||x-theFloor>1.0||x<theFloor||x>theCeiling) {
      cout << "nInt:\n"
	   << "Fatal double to integer conversion error.\n"
	   << "input double    = " << x << "\n"
	   << "integer ceiling = " << theCeiling << "\n"
	   << "integer floor   = " << theFloor << "\n"
	   << "Quitting ...";
      exit(1);
  }
  return (theCeiling-x) < (x-theFloor) ? theCeiling : theFloor;
}

int ndnsToLHAPDF(int ndns);

double parstrToparval(string varName,
		      vector<string> * parstrPtr,
		      vector<double> * parvalPtr);

void doIndividualHardProcessAssignments(int ihrd                , double * nup,
					vector<double> * idup   , vector<double> * istup,
					vector<double> * mothup1, vector<double> * mothup2,
					vector<double> * icolup1, vector<double> * icolup2,
					vector<vector<double> > * pup,
					vector<double> masses   , int itopprc);

//*****************************//
//*****************************//
//*****************************//
//                             //    
//       ATTENTION!!!!         //
//       ATTENTION!!!!         //
//       ATTENTION!!!!         //
//                             //    
// Remember to uncomment all   //
// of the currently commented  //
// commands outputting to      //
// the HW++ input file for     //
// the variables not yet       //
// declared and interfaced     //
// in the AlpGenHandler!!      // 
// These should only be the    //
// ones commented out by _4_   //
// consecutive forward slashes //
// '////' in the function      // 
// writeHWPPinFile(...)        // 
//                             //    
// These are only commented    //
// as the AlpGenHandler isn't  //
// ready for them yet!         //    
//                             //    
//*****************************//
//*****************************//
//*****************************//

void writeHWPPinFile(string prefix, int ihrd, int unwev,
		     int lhapdf, int idbmup0, int idbmup1, int idwtup,
		     double aqcdup, int nloop,
		     vector<double> * massesPtr,
		     vector<double> * parvalPtr,
		     vector<string> * parstrPtr);

string trim(string theString);

int main(int argc, char *argv[]) {

  bool usePowhegBoxConventions(true); // Control
  int  debugging(3);                  // To validate

  cout << "\n";
  cout << "------------------------------------------------------------------\n";
  cout << "  AlpGenToLH: Convert Alpgen Les Houches to Herwig++ Les Houches \n";
  cout << "                               v2.0-beta \n";
  cout << "------------------------------------------------------------------\n";
  cout << "\n";
  char* prefix;

 
  if(argv[1]) { prefix = argv[1]; } else { 
    cout << "Use: ./AlpGenToLH [input string] [number of events (optional)]\n"; exit(1);
    cout << "Note: As of version 2, the .stat and _unw.par files are "
	 << "required to read the generation parameters.\n";
  }
  int maxevents(0); int eventcount(0);
  if(argc>2) { maxevents=(atoi(argv[2])); }
  
  string lheFilename    = string(prefix) + string(".lhe");
  string unwFilename    = string(prefix) + string(".unw");
  string unwparFilename = string(prefix) + string("_unw.par");
  string statFilename   = string(prefix) + string(".stat");

  cout << "Opening input files ...\n";
  cout << "-----------------------\n";
  cout << "Unweighted events in "   << unwFilename    << ".\n"
       << "Generation settings in "
       << unwparFilename << " and " << statFilename   << ".\n\n";

  ifstream unwStream;
  ifstream unwparStream;
  ifstream statStream;
  ofstream lheStream;
  unwStream.open(unwFilename.c_str());
  unwparStream.open(unwparFilename.c_str());
  statStream.open(statFilename.c_str());
  lheStream.open(lheFilename.c_str());
  if(!unwStream) { 
    cerr << "error: Failed to open input file " << unwFilename << "\n";
    exit(1);
  }
  if(!unwparStream) { 
    cerr << "error: Failed to open input file " << unwparFilename << "\n";
    exit(1);
  }
  if(!statStream) {
    cerr << "error: Failed to open input file " << statFilename << "\n";
    exit(1);
  }

  // ******************************************************************** //
  // Dump the AlpGen *_unw.par file into the LH header (it's not so big). //
  // ******************************************************************** //

  string tmpString;
  lheStream << "<LesHouchesEvents version =\"1.0\">\n";
  lheStream << "<!--\n";
  lheStream << "File generated using AlpGen and converted with AlpGenToLH \n";
  while(unwparStream) {
    getline(unwparStream,tmpString);
    lheStream << tmpString << "\n";
  }
  lheStream << "\n\n\n" << "-->\n";
  unwparStream.close();

  // ***************************************** //
  // Read in all relevant info from *_unw.par. //
  // ***************************************** //

  int    ihrd;                       // AlpGen convention hard process code.
  double mc,mb,mt,mw,mz,mh;          // C, B, Top, W, Z & Higgs mass from *_unw.par
  double avgwgt,errwgt;              // Average weight and its error.
  int    unwev;                      // Number of unweighted events.
  double totlum;                     // Effective luminosity.
  vector<double> parval(200,-999.0); // AlpGen parameters.
  vector<string> parstr(200,"----"); // AlpGen parameter (variable) names.
  vector<double> alpgenParticleMasses;

  unwparStream.open(unwparFilename.c_str());
  while(unwparStream) {
    getline(unwparStream,tmpString);
    if(tmpString.find("hard process code") != string::npos) {
      tmpString=tmpString.substr(0,tmpString.find("!"));      
      tmpString=trim(tmpString);
      ihrd=atoi(tmpString.c_str());
    }
    if(tmpString.find("mc,mb,mt,mw,mz,mh") != string::npos) {
      tmpString=trim(tmpString.substr(0,tmpString.find("!")));
      mc=atof((tmpString.substr(0,tmpString.find_first_of(" "))).c_str());
      alpgenParticleMasses.push_back(mc);
      tmpString=trim(tmpString.substr(tmpString.find_first_of(" ")));
      mb=atof((tmpString.substr(0,tmpString.find_first_of(" "))).c_str());
      alpgenParticleMasses.push_back(mb);
      tmpString=trim(tmpString.substr(tmpString.find_first_of(" ")));
      mt=atof((tmpString.substr(0,tmpString.find_first_of(" "))).c_str());
      alpgenParticleMasses.push_back(mt);
      tmpString=trim(tmpString.substr(tmpString.find_first_of(" ")));
      mw=atof((tmpString.substr(0,tmpString.find_first_of(" "))).c_str());
      alpgenParticleMasses.push_back(mw);
      tmpString=trim(tmpString.substr(tmpString.find_first_of(" ")));
      mz=atof((tmpString.substr(0,tmpString.find_first_of(" "))).c_str());
      alpgenParticleMasses.push_back(mz);
      tmpString=trim(tmpString.substr(tmpString.find_first_of(" ")));
      mh=atof((tmpString.substr(0,tmpString.find_first_of(" "))).c_str());
      alpgenParticleMasses.push_back(mh);
    }
    if(tmpString.find("Crosssection +- error (pb)") != string::npos) {
      tmpString=trim(tmpString.substr(0,tmpString.find("!")));
      avgwgt=atof(trim(tmpString.substr(0,tmpString.find(" "))).c_str());
      errwgt=atof(trim(tmpString.substr(tmpString.find(" "))).c_str());
    }
    if(tmpString.find("unwtd events, lum (pb-1)") != string::npos) {
      tmpString=trim(tmpString.substr(0,tmpString.find("!")));
      unwev=atoi(trim(tmpString.substr(0,tmpString.find(" "))).c_str());
      totlum=atof(trim(tmpString.substr(tmpString.find(" "))).c_str());
    }
  }
  if(maxevents > unwev) { 
    cout << "-------------------------------\n";
    cout << "requested " << maxevents << " > " << unwev << " (contained in file), will use all events.\n"; maxevents = 0; }

  if(debugging>=4) {
    cout << "\nDebugging initial reading of *_unw.par:\n";
    cout << "ihrd = " << ihrd << "\n";
    cout << "mc,mb,mt,mw,mz,mh = "
	 <<  mc << "  " <<  mb << "  "
	 <<  mt << "  " <<  mw << "  "
	 <<  mz << "  " <<  mh << "\n";
    cout << "Cross section +/- error = "
	 <<  avgwgt << " +/- " << errwgt << "\n";
    cout << "Number of unweighted events = " <<  unwev << "\n";
    cout << "Effective luminosity = " <<  totlum << "\n";
  }
  unwparStream.close();

  unwparStream.open(unwparFilename.c_str());
  int index;
  while(unwparStream) {
    getline(unwparStream,tmpString);
    if(tmpString.find("!")==string::npos||
       tmpString.find("hard process code")!=string::npos||
       tmpString.find("mc,mb,mt,mw,mz,mh")!=string::npos||
       tmpString.find("Crosssection +- error (pb)")!=string::npos||
       tmpString.find("unwtd events, lum (pb-1)")!=string::npos) continue;
    tmpString=trim(tmpString);
    if(debugging>=4) cout << "\nDebugging reading paramters in *_unw.par:\n";
    if(debugging>=4) cout << "File says: " << tmpString << "\n";
    index = atoi((tmpString.substr(0,tmpString.find_first_of(" "))).c_str());
    tmpString=trim(tmpString.substr(tmpString.find_first_of(" ")));
    parval[index]=atof((tmpString.substr(0,tmpString.find_first_of(" "))).c_str());
    tmpString=trim(tmpString.substr(tmpString.find_first_of("!")+1));
    parstr[index]=tmpString;
    if(debugging>=4) cout << "We say:    "
		          << index << "  "
		          << parval[index] << "  "
		          << parstr[index] << "\n\n\n\n\n";
  }
  unwparStream.close();

  // Variables defined in parval array read from *_unw.par:

  // PDG codes for the beam particles.
  int idbmup[2]={0,0};
  if(parstrToparval("ih1",&parstr,&parval)==1)       idbmup[0] =  2212;
  else if(parstrToparval("ih1",&parstr,&parval)==-1) idbmup[0] = -2212;
  else idbmup[0] =  2212;
  if(parstrToparval("ih2",&parstr,&parval)==1)       idbmup[1] =  2212;
  else if(parstrToparval("ih2",&parstr,&parval)==-1) idbmup[1] = -2212;
  else idbmup[1] =  2212;

  // Energies of the beam particles --- implementation implicitly assumes
  // these are equal!
  double ebmup[2]={0,0};
  ebmup[0]=parstrToparval("ebeam",&parstr,&parval);
  ebmup[1]=ebmup[0];

  // LH accord pdf info variables for <init> block:
  int pdfgup[2];
  pdfgup[0]=-1;  // Simply set to 1 as in POWHEG-BOX.
  pdfgup[1]=-1;
  int pdfsup[2];
  // LHAPDF index: (note in POWHEG-BOX it is set to just -1).
  pdfsup[0]=ndnsToLHAPDF(int(parstrToparval("ndns",&parstr,&parval))); 
  pdfsup[1]=pdfsup[0];
  // LH accord flag defining weight scheme:
  // N.B. AlpGen alpsho.f UPINIT uses idwtup = 3 (this is likely better from the
  // point of view of combining events of diff multiplicity together in real life
  // i.e. in ATLAS - so we should probably use it!).
  int idwtup( 3); // As in POWHEG-BOX withnegweights 0 mode: unit wgts +1 only.
  //  int idwtup(-4); // As in POWHEG-BOX withnegweights 1 mode: +/- |xsecup| wgts only.
  // Number of processes in the file (assume all one process as
  // with alpsho.f UPINIT).
  int nprup(1);
  // Cross section, it's error, the maximum weight in the file.
  double xsecup,xerrup,xmaxup;
  xsecup = avgwgt;
  xerrup = errwgt;
  xmaxup = xsecup;
  // Process id code (to be augmented by jet multiplicity - see just below).
  int lprup(ihrd*100);

  // We augment the process code (for the LH file only) by the number of
  // (light) jets, just in case we end up connecting many different files
  // to the shower MC in parallel (otherwise it likely won't distinguish 
  // between X+0,1,2,3,...,n jet processes).
  int njets(int(parstrToparval("njets",&parstr,&parval)));
  lprup+=njets;

  // Write out some bits of info to the screen.
  cout << "No. of jets: "      << njets << "\n";
  cout << "Total xsec in pb (all processes): " << scientific
       << xsecup << " +/- " << xerrup << "\n\n";

  // NOW write out <init> block:
  lheStream << "<init>\n";
  lheStream << setw(9) << idbmup[0];
  lheStream << setw(9) << idbmup[1];
  lheStream << scientific << setprecision(5) << setw(13) << ebmup[0];
  lheStream << scientific << setprecision(5) << setw(13) << ebmup[1];
  lheStream << setw(7) << pdfgup[0];
  lheStream << setw(7) << pdfgup[1];
  if(usePowhegBoxConventions) {
    lheStream << setw(7) << -1;
    lheStream << setw(7) << -1;
  } else {
    lheStream << setw(7) << pdfsup[0];
    lheStream << setw(7) << pdfsup[1];
  }
  lheStream << setw(7) << idwtup;
  lheStream << setw(7) << nprup << "\n"; 

  lheStream << scientific << setprecision(5) << setw(13) << xsecup;
  lheStream << scientific << setprecision(5) << setw(13) << xerrup;
  if(usePowhegBoxConventions)
    lheStream << scientific << setprecision(5) << setw(13) << 1.0;
  else // else put in more info (xmaxup=xsecup).
    lheStream << scientific << setprecision(5) << setw(13) << xmaxup;
  lheStream << setw(7) << lprup << "\n";

  lheStream << "</init>\n";

  // ************************************************* //
  // All done with the init section of the LH file!    //
  // ************************************************* //

  // ***************************************************************** //
  // To write the HW++ .in file we have everything we could possibly   //
  // want except maybe the QCD coupling and no. of loops for running.  //
  // These are the only numbers we get / use from *.stat now.          //
  // ***************************************************************** //

  double aqedup(-999),aqedStat(-999);
  double aqcdup(-999),aqcdStat(-999);
  int nloop(-999);

  // Fish around for the QCD and QED alphas in .stat.
  while(statStream) {
    getline(statStream,tmpString);
    if(tmpString.find("as(MZ)") != string::npos) {
      aqcdup=atof(trim(tmpString.substr(tmpString.find_last_of("=")+1,
					tmpString.length())).c_str());
      tmpString=trim(tmpString.substr(tmpString.find_first_of("nloop=")));
      tmpString=trim(tmpString.substr(tmpString.find_first_of(" ")));
      tmpString=trim(tmpString.substr(0,tmpString.find_first_of("]")));
      nloop=atoi(tmpString.c_str());
    }
    if(tmpString.find("aem(mZ)") != string::npos) {
      tmpString=tmpString.substr(tmpString.find("aem(mZ)="));
      tmpString=tmpString.substr(tmpString.find_first_of(" "));
      tmpString=trim(tmpString);
      aqedup=atof(tmpString.c_str());
      aqedup=1/aqedup;
    }
  }
  statStream.close();
  aqedStat=aqedup;
  aqcdStat=aqcdup;

  // Write out a couple more bits of info to the screen.
  cout << "aqcdup [as(MZ)]  from stat file: " << aqcdup << "\n";
  cout << "nloop for as from stat file    : " << nloop << "\n";
  cout << "aqedup [inverse] from stat file: " << 1/aqedup << "\n";
  cout << "\n";

  writeHWPPinFile(prefix,ihrd,unwev,pdfsup[0], idbmup[0], idbmup[1], idwtup,
		  aqcdup,nloop,
		  &alpgenParticleMasses,
		  &parval,&parstr);

  // *********************************************************** //
  // All done writing the HW++ .in file!                         //
  // *********************************************************** //

  // *********************************************************** //
  // Start reading AlpGen events and writing them as LH events:  //
  // *********************************************************** //

  int nupMax(20);

  // First line of an event contains:
  double nup(0),idprup(0),xwgtup(0),scalup(0);

  // Subsequent lines for particle content contain:
  vector<double> idup,istup;
  vector<double> mothup1,mothup2;  // WARNING: not implemented ... YET!
  vector<double> icolup1,icolup2;
  vector<double> vtimup,spinup;    // WARNING: not implemented (but safe).
  vector<vector<double> > pup;
  pup.resize(5);
  for (int ixx = 0; ixx <= 4; ++ixx)
    pup[ixx].resize(nupMax);

  // Initialise pup matrix (not really necessary):
  for(unsigned int jup=0; jup < nupMax; jup++)
    for(unsigned int ixx=0; ixx<5; ixx++) pup[ixx][jup]=0;

  // Control reading AlpGen file:
  int iup(0);
  int counter(0);
  bool readInWholeEvent(false),beginNewEvent(true);
  double stringdoub;

  // Needed as input to doIndividualHardProcessAssignments only:
  int itopprc(nInt(parstrToparval("itopprc",&parstr,&parval)));

  while(unwStream && (eventcount < maxevents || maxevents==0)) { // So long as we haven't hit the EOF do ...

    if(beginNewEvent) {
      // Rest / set control variables:
      beginNewEvent=false;
      readInWholeEvent=false;
      counter=0;
      iup=0;
      // Reset variables for first line of LH event:
      nup=0; idprup=0; xwgtup=0; scalup=0; aqedup=0; aqcdup=0;
      // Reset variables for all individual particles in LH event:
      idup.clear()   ; istup.clear();
      mothup1.clear(); mothup2.clear();
      icolup1.clear(); icolup2.clear();
      vtimup.clear() ; spinup.clear();
      for(unsigned int jup=0; jup < nupMax; jup++)
	for(unsigned int ixx=0; ixx<5; ixx++) pup[ixx][jup]=0;
    }

    // Read in next thing starting from last position in
    // the file (int/real) as a double
    unwStream >> stringdoub;

    // counter counts the number of numbers read-in for a given
    // event. On starting to read a new event counter will be 0.
    counter++;

    switch (counter) {
    case  1: if(int(stringdoub)!=0) {
// 	        if(int(stringdoub)%100==0)
// 		  cout << "Processed " << fixed << setprecision(0)
// 		       << stringdoub/unwev*100
// 		       << " % of events ..." << "\r" << flush;
             }
             break;
    case  2: idprup = stringdoub;  break;
    case  3: nup    = stringdoub;  break;
    case  4: xwgtup = stringdoub;  break;
    case  5: scalup = stringdoub;  break;
    // N.B. There are no aqedup / aqcdup variables from AlpGen.

    // Initial state particle +z direction:
    case  6: idup.push_back(stringdoub)   ;
             istup.push_back(-1);  break;
    case  7: icolup1.push_back(stringdoub);
 	     mothup1.push_back(0); break; // ATTENTION: not given by AlpGen
    case  8: icolup2.push_back(stringdoub);
             mothup2.push_back(0); break; // ATTENTION: not given by AlpGen
    case  9: pup[2][iup]=stringdoub       ;
             pup[3][iup]=fabs(stringdoub) ; iup++ ; break;

    // Initial state particle -z direction:
    case 10: idup.push_back(stringdoub)   ;
             istup.push_back(-1);  break;
    case 11: icolup1.push_back(stringdoub);
             mothup1.push_back(0); break; // ATTENTION: not given by AlpGen
    case 12: icolup2.push_back(stringdoub);
             mothup2.push_back(0); break; // ATTENTION: not given by AlpGen
    case 13: pup[2][iup]=stringdoub;
             pup[3][iup]=fabs(stringdoub) ; iup++ ; break;
    }

    if(debugging<5) idprup=lprup;

    if(counter<14) continue;

    // Final state particles:
    if(counter==0+7*iup) idup.push_back(stringdoub);
    // istup gets assigned later on to just -1/+1 (I.S. / F.S.).
    if(counter==1+7*iup) { icolup1.push_back(stringdoub);
                           mothup1.push_back(1.); } // ATTENTION: not given by AlpGen
    if(counter==2+7*iup) { icolup2.push_back(stringdoub);
                           mothup2.push_back(2.); } // ATTENTION: not given by AlpGen
    if(counter==3+7*iup) pup[0][iup] = stringdoub;
    if(counter==4+7*iup) pup[1][iup] = stringdoub;
    if(counter==5+7*iup) pup[2][iup] = stringdoub;
    if(counter==6+7*iup) pup[4][iup] = stringdoub;
    if(counter==6+7*iup) istup.push_back(1.);
    if(counter==6+7*iup) iup+=1;
    if(counter==7*nup-1) readInWholeEvent = true; 

    for(int jup = 0; jup < nup; jup++) {
      pup[3][jup] = sqrt( sqr(pup[4][jup]) + sqr(pup[0][jup]) 
			+ sqr(pup[1][jup]) + sqr(pup[2][jup]) );
      vtimup.push_back(0.); // ATTENTION: not implemented - so taking 
      spinup.push_back(9.); // POWHEG-BOX default values (should be v.safe).
    }

    // ***************************************************************** //
    // Now consider assignments specific to individual hard processes:   //
    // ***************************************************************** //    

    if(readInWholeEvent)
      doIndividualHardProcessAssignments(ihrd, &nup, &idup, &istup,
					 &mothup1  , &mothup2,
					 &icolup1  , &icolup2,
					 &pup      ,
					 alpgenParticleMasses, itopprc);

    if(readInWholeEvent) {
      lheStream << "<event>\n";
      if(debugging>=5) {
	lheStream << nup    << "\t" << idprup << "\t" << xwgtup << "\t"
		  << scalup << "\t" << "0.007297352" << "\t" << "0.118\n";
      } else {
	// Bit about signing of xwgtup here is redundant as 
	// AlpGen only gives +ve weight events ...
	double signOfXwgtup = xwgtup >= 0 ? 1 : -1;
	xwgtup = idwtup==3 ? 1 : xsecup*signOfXwgtup;
	// N.B. There are no aqedup / aqcdup variables from AlpGen
	// events only the stat file has related information.
	aqedup=aqedStat;
	aqcdup=aqcdStat;
	if(usePowhegBoxConventions) aqedup=-1;
	lheStream << setw(7) << int(nup);
	lheStream << setw(7) << int(idprup);
	lheStream << scientific << setprecision(5) << setw(13) << xwgtup;
	lheStream << scientific << setprecision(5) << setw(13) << scalup;
	lheStream << scientific << setprecision(5) << setw(13) << aqedup;
	lheStream << scientific << setprecision(5) << setw(13) << aqcdup;
	lheStream << "\n";
      }
      for(int jup = 0; jup < nup; jup++) {
	if(debugging>=5) {
	  lheStream << idup[jup]     << "\t" << istup[jup]    << "\t"
		    << mothup1[jup]  << "\t" << mothup2[jup]  << "\t" 
		    << icolup1[jup]  << "\t" << icolup2[jup]  << "\t";
	  lheStream << pup[0][jup]   << "\t" << pup[1][jup]   << "\t"
		    << pup[2][jup]   << "\t" << pup[3][jup]   << "\t"
		    << pup[4][jup]   << "\t";
	  lheStream << "0"           << "\t" << "1\n";
	} else {
	  if(icolup1[jup]!=0) icolup1[jup]+=500;
	  if(icolup2[jup]!=0) icolup2[jup]+=500;
	  lheStream << setw(8) << int(idup[jup]);
	  lheStream << setw(6) << int(istup[jup])
		    << setw(6) << int(mothup1[jup])
		    << setw(6) << int(mothup2[jup])
		    << setw(6) << int(icolup1[jup])
		    << setw(6) << int(icolup2[jup]);
	  lheStream << scientific << setprecision(9) << setw(17) << pup[0][jup];
	  lheStream << scientific << setprecision(9) << setw(17) << pup[1][jup];
	  lheStream << scientific << setprecision(9) << setw(17) << pup[2][jup];
	  lheStream << scientific << setprecision(9) << setw(17) << pup[3][jup];
	  lheStream << scientific << setprecision(9) << setw(17) << pup[4][jup];
	  lheStream << scientific << setprecision(5) << setw(13) << vtimup[jup];
	  lheStream << scientific << setprecision(3) << setw(11) << vtimup[jup];
	  lheStream << "\n";
	}
      }
      lheStream << "</event>\n";
      eventcount++;
      cout << "Processed " << eventcount << " events ..." << "\r" << flush;
      beginNewEvent=true;
    }
  }
  
  cout << "\n\n";

  if( maxevents!=0 ) {
    cout << "All done (" << maxevents << " out of " <<  unwev << " events).\n";
  } else { 
    cout << "All done (" << eventcount << " events).\n";
  }

  cout << "\n\n"
       << "Wrote a LH event file " << lheFilename
       << " and a HW++ MLM merging input file "
       << prefix+string(".in") << ".\n\n";
  return 0;

}


inline double sqr(double x) { return x*x; }

int ndnsToLHAPDF(int ndns) {
  // The information in this function is based on 
  // subroutine PRNTSF from alplib/alppdf.f, LHAPDF's 
  // PDFsets.index and, finally, the .stat output that
  // results when the relevant ndns value is entered
  // in the input file.
  string Set("no PDF set found");
  double Lambda_4(0),Lambda_5_2loop(0);
  string Scheme("no PDF scheme");
  int    LHAPDFindex(-999);
  string tmpString("");

  if(ndns==1) {
    Set = "CTEQ4M"       ; Lambda_4 = 0.298 ; Lambda_5_2loop = 0.202 ; Scheme = "MS" ;
    LHAPDFindex =  19150;
  } else if(ndns==2) {
    Set = "CTEQ4L"       ; Lambda_4 = 0.298 ; Lambda_5_2loop = 0.202 ; Scheme = "MS" ;
    LHAPDFindex =  19170;
  } else if(ndns==3) {
    Set = "CTEQ4HJ"      ; Lambda_4 = 0.298 ; Lambda_5_2loop = 0.202 ; Scheme = "MS" ;
    LHAPDFindex = -99999;
  } else if(ndns==4) {
    Set = "CTEQ5M"       ; Lambda_4 = 0.326 ; Lambda_5_2loop = 0.226 ; Scheme = "MS" ;
    LHAPDFindex =  19050;
  } else if(ndns==5) {
    Set = "CTEQ5L"       ; Lambda_4 = 0.192 ; Lambda_5_2loop = 0.144 ; Scheme = "MS" ;
    LHAPDFindex =  19070;
  } else if(ndns==6) {
    Set = "CTEQ5HJ"      ; Lambda_4 = 0.326 ; Lambda_5_2loop = 0.226 ; Scheme = "MS" ;
    LHAPDFindex = -99999;
  } else if(ndns==7) {
    Set = "CTEQ6M"       ; Lambda_4 = 0.326 ; Lambda_5_2loop = 0.226 ; Scheme = "MS" ;
    LHAPDFindex=   10050;
  } else if(ndns==8) {
    Set = "CTEQ6L"       ; Lambda_4 = 0.326 ; Lambda_5_2loop = 0.226 ; Scheme = "MS" ;
    LHAPDFindex =  10041;
  } else if(ndns==9) {
    Set = "CTEQ6L1"      ; Lambda_4 = 0.215 ; Lambda_5_2loop = 0.167 ; Scheme = "MS" ;
    LHAPDFindex =  10042;
  } else if(ndns>=10&&ndns<=50) {
    Set = "CTEQ6xx"      ; Lambda_4 = 0.326 ; Lambda_5_2loop = 0.226 ; Scheme = "MS" ;
    LHAPDFindex =  10150+(ndns-10);
  } else if(ndns==101) {
    Set = "MRST99"       ; Lambda_4 = 0.321 ; Lambda_5_2loop = 0.220 ; Scheme = "MS" ;
    LHAPDFindex = -99999;
  } else if(ndns==102) {
    Set = "MRST01"       ; Lambda_4 = 0.342 ; Lambda_5_2loop = 0.239 ; Scheme = "MS" ;
    LHAPDFindex = -99999;
  } else if(ndns==103) {
    Set = "MRST01"       ; Lambda_4 = 0.310 ; Lambda_5_2loop = 0.214 ; Scheme = "MS" ;
    LHAPDFindex = -99999;
  } else if(ndns==104) {
    Set = "MRST01"       ; Lambda_4 = 0.378 ; Lambda_5_2loop = 0.267 ; Scheme = "MS" ;
    LHAPDFindex = -99999;
  } else if(ndns==105) {
    Set = "MRST01J"      ; Lambda_4 = 0.378 ; Lambda_5_2loop = 0.267 ; Scheme = "MS" ;
    LHAPDFindex = -99999;
  } else if(ndns==106) {
    Set = "MRST02LO"     ; Lambda_4 = 0.215 ; Lambda_5_2loop = 0.167 ; Scheme = "MS" ;
    LHAPDFindex = -99999;
  } else if(ndns==201) {
    Set = "MSTW2008lo"   ; Lambda_4 = 0.322 ; Lambda_5_2loop = 0.255 ; Scheme = "MS" ;
    LHAPDFindex =  21000;
  } else if(ndns==202) {
    Set = "MSTW2008nlo"  ; Lambda_4 = 0.365 ; Lambda_5_2loop = 0.255 ; Scheme = "MS" ;
    LHAPDFindex =  21100;
  } else if(ndns>=203&&ndns<=242) {
    Set = "MSTW2008lo68cl"; Lambda_4 = 0.322 ; Lambda_5_2loop = 0.255 ; Scheme = "MS" ;
    LHAPDFindex =  21000+(ndns-202);
  } else if(ndns==243) {
    Set = "MRST LO*" ; Lambda_4 = 0.365 ; Lambda_5_2loop = 0.255 ; Scheme = "MS" ;
    LHAPDFindex =  20650;
  } else if(ndns==244) {
    Set = "MRST LO**" ; Lambda_4 = 0.280 ; Lambda_5_2loop = 0.190 ; Scheme = "MS" ;
    LHAPDFindex =  20651;
  } else if(ndns==301 ) {
    Set = "CTQ6.6" ; Lambda_4 = 0.326 ; Lambda_5_2loop = 0.226 ; Scheme = "MS" ;
    LHAPDFindex =  10550;
  } else if(ndns>=302&&ndns<=345) {
    Set = "CTQ66" ; Lambda_4 = 0.326 ; Lambda_5_2loop = 0.226 ; Scheme = "MS" ;
    LHAPDFindex =  10550+(ndns-301);
  } else if(ndns==346) {
    Set = "CT09MC1" ; Lambda_4 = 0.215 ; Lambda_5_2loop = 0.167 ; Scheme = "MS" ;
    LHAPDFindex =  10771;
  } else if(ndns==347) {
    Set = "CT09MC2" ; Lambda_4 = 0.326 ; Lambda_5_2loop = 0.226 ; Scheme = "MS" ;
    LHAPDFindex =  10772;
  }

  cout << "-------------------------------\n";
  cout << "ndnsToLHAPDF found: \n";
  cout << "PDF set      = " << Set << "\n";
  cout << "ndns index   = " << ndns << "\n";
  cout << "LHAPDF index = " << LHAPDFindex << "\n";
  cout << "-------------------------------\n\n";
  return LHAPDFindex;

}


double parstrToparval(string varName,
		      vector<string> * parstrPtr,
		      vector<double> * parvalPtr) {
  for(unsigned int index=0; index<parvalPtr->size(); index++)
    if(varName==parstrPtr->at(index))
      return parvalPtr->at(index);

  return -999.0;
}

string trim(string theString) {
  int endStr    = theString.find_last_not_of(" ");
  int beginStr  = theString.find_first_not_of(" ");
  if(beginStr==0&&endStr==theString.length()-1) return theString; // No lead / trail spaces.
  theString = theString.substr(beginStr,endStr-beginStr+1);
  return theString;
}

void writeHWPPinFile(string prefix, int ihrd, int unwev,
		     int lhapdf, int idbmup0, int idbmup1, int idwtup,
		     double aqcdup, int nloop,
		     vector<double> * massesPtr,
		     vector<double> * parvalPtr,
		     vector<string> * parstrPtr) {
  ofstream hwpp;
  hwpp.open(string(prefix+".in").c_str());

  hwpp << "#############################################################\n";
  hwpp << "# Create an event generator taking the default LHCGenerator #\n";
  hwpp << "# as the starting point ...                                 #\n";
  hwpp << "#############################################################\n";
  hwpp << "cd /Herwig/Generators\n";
  hwpp << "# Copy the default LHCGenerator with its settings to a new \n";
  hwpp << "# which will be the basis of the one we use for showering: \n";
  hwpp << "cp LHCGenerator theGenerator\n";
  hwpp << "\n";
  hwpp << "#############################################################\n";
  hwpp << "# Create a LH event handler (set up & assigned below) ...   #\n";
  hwpp << "#############################################################\n";
  hwpp << "cd /Herwig/EventHandlers\n";
  hwpp << "library LesHouches.so\n";
  hwpp << "create ThePEG::LesHouchesEventHandler theLesHouchesHandler\n";
  hwpp << "\n";
  hwpp << "#############################################################\n";
  hwpp << "# Create a LH reader (set up & assigned below) ...          #\n";
  hwpp << "#############################################################\n";
  hwpp << "cd /Herwig/EventHandlers\n";
  hwpp << "library BasicLesHouchesFileReader.so\n";
  hwpp << "create Herwig::BasicLesHouchesFileReader theLHReader\n";
  hwpp << "\n";
  hwpp << "#############################################################\n";
  hwpp << "# Create an AlpGenHandler (set up & assigned below) ...     #\n";
  hwpp << "#############################################################\n";
  hwpp << "cd /Herwig/Shower\n";
  hwpp << "library AlpGenHandler.so\n";
  hwpp << "create Herwig::AlpGenHandler AlpGenHandler\n";
  hwpp << "\n";
  hwpp << "#############################################################\n";
  hwpp << "# Create an LHAPDF (set up & assigned below) ...            #\n";
  hwpp << "#############################################################\n";
  hwpp << "cd /Herwig/Partons\n";
  hwpp << "create ThePEG::LHAPDF thePDFset ThePEGLHAPDF.so\n";
  hwpp << "\n";
  hwpp << "############################################################\n";
  hwpp << "# Create a cuts object ...                                 #\n";
  hwpp << "############################################################\n";
  hwpp << "cd /Herwig/EventHandlers\n";
  hwpp << "create ThePEG::Cuts   /Herwig/Cuts/NoCuts\n";
  hwpp << "\n";
  hwpp << "#############################################################\n";
  hwpp << "# Setup the LH event handler ...                            #\n";
  hwpp << "#############################################################\n";
  hwpp << "cd /Herwig/EventHandlers\n";
  hwpp << "insert theLesHouchesHandler:LesHouchesReaders 0 theLHReader\n";
  if(idwtup==3) {
    hwpp << "set theLesHouchesHandler:WeightOption UnitWeight\n";
  } else if(idwtup==-3) {
    hwpp << "set theLesHouchesHandler:WeightOption NegUnitWeight\n";
  } else if(idwtup==4) {
    hwpp << "set theLesHouchesHandler:WeightOption VarWeight\n";
  } else {
    hwpp << "set theLesHouchesHandler:WeightOption VarNegWeight\n";
  }
  hwpp << "set theLesHouchesHandler:PartonExtractor "
       << "/Herwig/Partons/QCDExtractor\n";
  hwpp << "set theLesHouchesHandler:CascadeHandler "
       << "/Herwig/Shower/AlpGenHandler\n";
  hwpp << "set theLesHouchesHandler:HadronizationHandler "
       << "/Herwig/Hadronization/ClusterHadHandler\n";
  hwpp << "set theLesHouchesHandler:DecayHandler "
       << "/Herwig/Decays/DecayHandler\n";
  hwpp << "\n";
  hwpp << "#############################################################\n";
  hwpp << "# Set up the Evolver to veto hard emissions > scalup ...    #\n";
  hwpp << "#############################################################\n";
  hwpp << "cd /Herwig/Shower\n";
  hwpp << "# MaxTry 100 sets the maximum number of times to try \n";
  hwpp << "# showering a given shower tree to 100. \n";
  hwpp << "# HardVetoMode Yes to veto emissions with pT greater than pT_max.\n";
  hwpp << "# HardVetoScaleSource Read means pT_max comes from hepeup.SCALUP.\n";
  hwpp << "# This is what you need to set _along_with_ HardVetoMode Yes in \n";
  hwpp << "# the case of Powheg external events _AND_ mc@nlo (we know this \n";
  hwpp << "# from looking at the *MCinput file that mc@nlo generates). \n";
  hwpp << "# MeCorrMode No turns off ME corrs. IntrinsicPtGaussian  2.2*GeV \n";
  hwpp << "# is the RMS of intrinsic pT of Gaussian distribution: \n";
  hwpp << "# 2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS) \n";
  hwpp << "set Evolver:MaxTry               100\n";
  hwpp << "set Evolver:HardVetoMode         Yes\n";
  hwpp << "set Evolver:HardVetoScaleSource  Read\n";
  hwpp << "set Evolver:HardVetoReadOption   PrimaryCollision\n";
  hwpp << "set Evolver:MECorrMode           No\n";
  hwpp << "# Intrinsic pT etc should be set as part of a tune i.e. it \n";
  hwpp << "# should either be left alone (default) or set by reading in \n";
  hwpp << "# one of the tunes before theGenerator is created by copying \n";
  hwpp << "# LHCGenerator (second line of this file). The following \n";
  hwpp << "# settings were extracted from LHC-UE-EE-3-CTEQ6L1.in - In \n";
  hwpp << "# light of some bad experience with MPI we prefer to set \n";
  hwpp << "# these here manually rather than read try to read that .in. \n";
  if(idbmup0 == 2212 && idbmup1 == -2212) {
    hwpp << "set /Herwig/UnderlyingEvent/KtCut:MinKT 2.26 \n";
    hwpp << "set /Herwig/UnderlyingEvent/UECuts:MHatMin 4.52 \n";
    hwpp << "set /Herwig/Shower/Evolver:IntrinsicPtGaussian 1.9*GeV \n";
  } else {
    hwpp << "set /Herwig/UnderlyingEvent/KtCut:MinKT 2.752 \n";
    hwpp << "set /Herwig/UnderlyingEvent/UECuts:MHatMin 5.504 \n";
    hwpp << "set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.34*GeV \n";
  }
  hwpp << "# Colour reconnection (re)settings \n";
  hwpp << "set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes \n";
  hwpp << "set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.61\n";
  hwpp << "# Colour Disrupt settings \n";
  hwpp << "set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.75 \n";
  hwpp << "# Inverse hadron radius \n";
  hwpp << "set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 1.35 \n";
  hwpp << "set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes \n";
  hwpp << "set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes \n";
  hwpp << "set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2 \n";
  hwpp << "\n";
  hwpp << "#############################################################\n";
  hwpp << "# Set up kinematics reconstructor (relevant only to mc@nlo) #\n";
  hwpp << "#############################################################\n";
  hwpp << "# Options for QTildeReconstructor - not needed for Powheg but\n";
  hwpp << "# critical for mc@nlo. If using Powheg you may either leave \n";
  hwpp << "# the next two settings as they are or comment them out as \n";
  hwpp << "# you wish, for mc@nlo though you must leave them as they \n";
  hwpp << "# are! ReconstructionOption General was the old default, it \n";
  hwpp << "# ignores the colour structure for all processes - mc@nlo \n";
  hwpp << "# will give you garbage unless you set this! \n";
  hwpp << "# ReconstructionOption Colour is the new default - use the \n";
  hwpp << "# colour structure of the process to determine the \n";
  hwpp << "# reconstruction procedure. InitialInitialBoostOption \n";
  hwpp << "# determines how the boost from the system before ISR to that\n";
  hwpp << "# after ISR is applied. mc@nlo requires the old kinematics \n";
  hwpp << "# reconstruction method: \n";
  hwpp << "# InitialInitialBoostOption LongTransBoost - first apply a \n";
  hwpp << "# longitudinal and then a transverse boost. Whereas the \n";
  hwpp << "# default method now is InitialInitialBoostOption OneBoost - \n";
  hwpp << "# apply one boost from old CMS to new CMS. Both options \n";
  hwpp << "# should work with Powheg but probably it's best to use the \n";
  hwpp << "# defaults in that case by simply commenting this setting. \n";
  hwpp << "set KinematicsReconstructor:ReconstructionOption General \n";
  hwpp << "set KinematicsReconstructor:InitialInitialBoostOption LongTransBoost\n";
  hwpp << "\n";
  hwpp << "#############################################################\n";
  hwpp << "# Set up the AlpGenHandler ...                              #\n";
  hwpp << "#############################################################\n";
  hwpp << "cd /Herwig/Shower\n";
  hwpp << "set AlpGenHandler:MPIHandler  /Herwig/UnderlyingEvent/MPIHandler\n";
  hwpp << "set AlpGenHandler:RemDecayer  /Herwig/Partons/RemnantDecayer\n";
  hwpp << "set AlpGenHandler:Evolver     Evolver\n";
  hwpp << "set AlphaQCD:AlphaMZ       " << aqcdup << "\n";
  hwpp << "set AlphaQCD:NumberOfLoops " << nloop << "\n";
  hwpp << "set AlpGenHandler:ShowerAlpha  AlphaQCD\n";
  hwpp << "# Calorimeter granularity settings used by GetJet algorithm\n"; 
  hwpp << "set AlpGenHandler:NoCellsInRapidity 100\n"; 
  hwpp << "set AlpGenHandler:NoCellsInPhi       60\n";
  // AlpGen hard process code.
  hwpp << "# AlpGen hard process code.\n";
  hwpp << "set AlpGenHandler:ihrd        " << ihrd << "\n";
  // Number of (light) jets.
  int njets(int(parstrToparval("njets",parstrPtr,parvalPtr)));
  hwpp << "# No. of light jets in AlpGen process (the \"extra\" ones).\n"; 
  hwpp << "set AlpGenHandler:njets       " << njets << "\n";
  // Mimimum jet pT use for generation.
  double ptjmin(parstrToparval("ptjmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:ptjmin      " << ptjmin << "*GeV\n";
  // Mimimum parton-parton R-sep used for generation.
  double drjmin(parstrToparval("drjmin",parstrPtr,parvalPtr));
  hwpp << "# Mimimum parton-parton R-sep used for generation.\n"; 
  hwpp << "set AlpGenHandler:drjmin      " << drjmin << "\n";
  // Max |eta| for partons in generation.
  double etajmax(parstrToparval("etajmax",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:etajmax     " << etajmax << "\n";
  // Also probably want these variables fed to AlpGenHandler too ---
  // they get set in the alpsho.f AHspar routine (note the list below
  // does not include some variables from AHspar because they are already
  // included in the above eg the PDFs are already handled so I removed
  // ndns and also ptjmin drjmin are written out for the AlpGenHandler above).
  int ickkw(int(parstrToparval("ickkw",parstrPtr,parvalPtr)));
  ////  hwpp << "set AlpGenHandler:ickkw       " << ickkw << "\n";
  int ihvy(int(parstrToparval("ihvy",parstrPtr,parvalPtr)));
  hwpp << "# heavy flavour in WQQ,ZQQ,2Q etc (4=c, 5=b, 6=t):\n";
  hwpp << "set AlpGenHandler:ihvy              " << ihvy << "\n";
  int ihvy2(int(parstrToparval("ihvy2",parstrPtr,parvalPtr)));
  ////  hwpp << "set AlpGenHandler:ihvy2       " << ihvy2 << "\n";
  int itopprc(nInt(parstrToparval("itopprc",parstrPtr,parvalPtr)));
  ////  hwpp << "set AlpGenHandler:itopprc     " << itopprc << "\n";
  int nw(int(parstrToparval("nw",parstrPtr,parvalPtr)));
  if(ihrd==13) {
    nw=1;                // N.B. nw is reassigned in this way 
    if(itopprc>=3) nw=2; //      by UPEVNT (after UPINIT).
  }
  ////  hwpp << "set AlpGenHandler:nw          " << nw << "\n";
  int nz(int(parstrToparval("nz",parstrPtr,parvalPtr)));
  ////  hwpp << "set AlpGenHandler:nz          " << nz << "\n";
  int nh(int(parstrToparval("nh",parstrPtr,parvalPtr)));
  hwpp << "# Number of Higgses in the AlpGen process:\n";
  hwpp << "set AlpGenHandler:nh          " << nh << "\n";
  int nph(int(parstrToparval("nph",parstrPtr,parvalPtr)));
  hwpp << "# Number of photons in the AlpGen process:\n";
  hwpp << "set AlpGenHandler:nph         " << nph << "\n";
  double ptbmin(parstrToparval("ptbmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:ptbmin      " << ptbmin << "\n";
  double ptcmin(parstrToparval("ptcmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:ptcmin      " << ptcmin << "\n";
  double ptlmin(parstrToparval("ptlmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:ptlmin      " << ptlmin << "\n";
  double metmin(parstrToparval("metmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:metmin      " << metmin << "\n";
  double ptphmin(parstrToparval("ptphmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:ptphmin     " << ptphmin << "\n";
  double etabmax(parstrToparval("etabmax",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:etabmax     " << etabmax << "\n";
  double etacmax(parstrToparval("etacmax",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:etacmax     " << etacmax << "\n";
  double etalmax(parstrToparval("etalmax",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:etalmax     " << etalmax << "\n";
  double etaphmax(parstrToparval("etaphmax",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:etaphmax    " << etaphmax << "\n";
  double drbmin(parstrToparval("drbmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:drbmin      " << drbmin << "\n";
  double drcmin(parstrToparval("drcmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:drcmin      " << drcmin << "\n";
  double drlmin(parstrToparval("drlmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:drlmin      " << drlmin << "\n";
  double drphjmin(parstrToparval("drphjmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:drphjmin    " << drphjmin << "\n";
  double drphlmin(parstrToparval("drphlmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:drphlmin    " << drphlmin << "\n";
  double drphmin(parstrToparval("drphmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:drphmin     " << drphmin << "\n";
  double mllmin(parstrToparval("mllmin",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:mllmin      " << mllmin << "\n";
  double mllmax(parstrToparval("mllmax",parstrPtr,parvalPtr));
  ////  hwpp << "set AlpGenHandler:mllmax      " << mllmax << "\n";
  hwpp << "\n";
  hwpp << "#############################################################\n";
  hwpp << "# Set up the LH reader ...                                  #\n";
  hwpp << "#############################################################\n";
  hwpp << "cd /Herwig/EventHandlers\n";
  hwpp << "set theLHReader:WeightWarnings    false\n";
  hwpp << "# Input event file name:\n";
  hwpp << "set theLHReader:FileName               " << prefix << ".lhe\n";
  hwpp << "set theLHReader:MomentumTreatment      RescaleEnergy\n";
  hwpp << "# set theLHReader:IgnoreIDPRUP           Yes\n";
  hwpp << "set theLHReader:Cuts  /Herwig/Cuts/NoCuts\n";
  hwpp << "set theLHReader:OverSampling ForbidOverSampling\n";
  hwpp << "\n";
  hwpp << "#############################################################\n";
  hwpp << "# Set up the LHAPDF ...                                     #\n";
  hwpp << "#############################################################\n";
  hwpp << "cd /Herwig/Partons\n";
  hwpp << "# Don't try and find PDF index out from the LH file ...\n";
  hwpp << "set /Herwig/EventHandlers/theLHReader:InitPDFs false\n";
  hwpp << "# Instead set them explicitly here:\n";
  hwpp << "set thePDFset:PDFNumber       " << lhapdf << "\n";
  hwpp << "set thePDFset:RemnantHandler  HadronRemnants\n";
  hwpp << "set /Herwig/EventHandlers/theLHReader:PDFA thePDFset\n";
  hwpp << "set /Herwig/EventHandlers/theLHReader:PDFB thePDFset\n";
  hwpp << "set /Herwig/Particles/p+:PDF    thePDFset\n";
  hwpp << "set /Herwig/Particles/pbar-:PDF thePDFset\n";
  hwpp << "# The PDF for beam particles A/B - overrides particle's own "
       << "PDF above\n";
  hwpp << "set /Herwig/Shower/AlpGenHandler:PDFA    thePDFset\n";
  hwpp << "set /Herwig/Shower/AlpGenHandler:PDFB    thePDFset\n";
  hwpp << "set /Herwig/Shower/ShowerHandler:PDFA thePDFset\n";
  hwpp << "set /Herwig/Shower/ShowerHandler:PDFB thePDFset\n";
  hwpp << "\n";
  hwpp << "####################################################\n";
  hwpp << "# Set up the generator ...                         #\n";
  hwpp << "####################################################\n";
  hwpp << "cd /Herwig/Generators\n";
  hwpp << "set theGenerator:EventHandler  "
       << "/Herwig/EventHandlers/theLesHouchesHandler\n";
  hwpp << "set theGenerator:NumberOfEvents " <<  unwev << "\n";
  hwpp << "set theGenerator:RandomNumberGenerator:Seed 31122001\n";
  hwpp << "set theGenerator:PrintEvent     10\n";
  hwpp << "set theGenerator:MaxErrors      10000\n";
  hwpp << "\n";
  hwpp << "###################################################################\n";
  hwpp << "# ReDefine particle data like it is in the AlpGen parameter file. #\n";
  hwpp << "###################################################################\n";
  hwpp << "\ncd  /Herwig/Particles/ \n";
  // 'if' statements needed here to protect against mc<ms (i.e. if AlpGen has mc=0
  // for example). If this occurs then the QCD coupling can't initialise and the
  // HW++ read step fails (due to matching at flavour thresholds etc).
  if(massesPtr->at(0)>1.0) {
    hwpp << "set c:NominalMass           " << massesPtr->at(0) << "*GeV\n";
    hwpp << "set cbar:NominalMass        " << massesPtr->at(0) << "*GeV\n";
  }
  // Ditto.
  if(massesPtr->at(0)>4.0) {
    hwpp << "set b:NominalMass           " << massesPtr->at(1) << "*GeV\n";
    hwpp << "set bbar:NominalMass        " << massesPtr->at(1) << "*GeV\n";
  }
  hwpp << "set t:NominalMass           " << massesPtr->at(2) << "*GeV\n";
  hwpp << "set tbar:NominalMass        " << massesPtr->at(2) << "*GeV\n";
  hwpp << "set W+:NominalMass          " << massesPtr->at(3) << "*GeV\n";
  hwpp << "set W-:NominalMass          " << massesPtr->at(3) << "*GeV\n";
  hwpp << "set Z0:NominalMass          " << massesPtr->at(4) << "*GeV\n";
  hwpp << "set h0:NominalMass          " << massesPtr->at(5) << "*GeV\n";
  hwpp << "\n\n\n\n\n";
  hwpp << "######################################################### \n";
  hwpp << "######################################################### \n";
  hwpp << "##                                                      # \n";
  hwpp << "##   --- USER SERVICEABLE PART BELOW HERE ONLY ! ---    # \n";
  hwpp << "##                                                      # \n";
  hwpp << "######################################################### \n";
  hwpp << "######################################################### \n";
  hwpp << "\n\n\n\n\n";
  hwpp << "######################################################### \n";
  hwpp << "# Option to off shower / hadronization / decays / MPI.  # \n";
  hwpp << "######################################################### \n";
  hwpp << "cd /Herwig/EventHandlers \n";
  hwpp << "# set theLesHouchesHandler:CascadeHandler        NULL \n";
  hwpp << "set theLesHouchesHandler:HadronizationHandler  NULL \n";
  hwpp << "set theLesHouchesHandler:DecayHandler          NULL \n";
  hwpp << "# The handler for multiple parton interactions \n";
  hwpp << "set /Herwig/Shower/AlpGenHandler:MPIHandler       NULL \n";
  hwpp << "\n\n";
  hwpp << "######################################################### \n";
  hwpp << "# Recommended key MLM merging parameters below - change # \n";
  hwpp << "# for systematic error studies and / or at your peril.  # \n";
  hwpp << "######################################################### \n";
  hwpp << "cd /Herwig/Shower\n";
  hwpp << "# Is this the highest multiplicity ME in merging? \n";
  hwpp << "# 0 = no, 1 = yes . \n";
  hwpp << "set AlpGenHandler:highestMultiplicity  0 \n";
  hwpp << "# Jet ET cut to apply in jet clustering in merging.\n";
  double etclus(max(ptjmin+5,1.2*ptjmin));
  hwpp << "set AlpGenHandler:ETClus " << etclus << "*GeV\n";
  hwpp << "# Cone size used in clustering in merging.\n";
  double rclus(drjmin);
  hwpp << "set AlpGenHandler:RClus " << rclus << "\n";
  hwpp << "# Max |eta| for jets in clustering in merging.\n";
  double etaclmax(etajmax);
  hwpp << "set AlpGenHandler:EtaClusMax " << etaclmax << "\n";
  hwpp << "# Default 1.5 factor used to decide if a jet matches a parton\n";
  hwpp << "# in merging: if DR(parton,jet)<rclusfactor*rclus the parton \n";
  hwpp << "# and jet are said to have been matched.\n";
  double rclusfactor(1.5);
  hwpp << "set AlpGenHandler:RClusFactor " << rclusfactor << "\n";
  hwpp << "\n\n";
  hwpp << "################ \n";
  hwpp << "# Save the run # \n";
  hwpp << "################ \n";
  hwpp << "cd /Herwig/Generators \n";
  hwpp << "saverun " << prefix << " theGenerator\n";

}

// Now consider assignments specific to individual hard processes:
void doIndividualHardProcessAssignments(int ihrd                , double * nup,
					vector<double> * idup   , vector<double> * istup,
					vector<double> * mothup1, vector<double> * mothup2,
					vector<double> * icolup1, vector<double> * icolup2,
					vector<vector<double> > * pup,
					vector<double> masses   , int itopprc) {

  int iwch(0);
  
  // W/Z/gamma b bbar + jets ( wcjet*, ihrd=10 / wphjet*, ihrd=14 / wphqq*,
  //  ihrd=15 ), or W/Z + jets ( wjet*, ihrd=3 / zjet*, ihrd=4 ):
  // In this case we add to the list of particles the single intermediate
  // boson (at the end of the list) appropriately, and assign the relevant
  // parent-daughter and colour flow indices.
  if (ihrd<=4||ihrd==10||ihrd==14||ihrd==15) {
    iwch=0; // <--- used to determine type: W/Z/gamma
    for(int iup=int(*nup)-2; iup<int(*nup); iup++) {
      (*mothup1)[iup]=*nup+1; // Assigning, to-be-added boson, as
      (*mothup2)[iup]=0;      // the parent of each decay prod.
      if(ihrd!=2) iwch = iwch - (int((*idup)[iup])%2);
      // electron+nubar -> 11 + (-12) => -(1)+0 = -1  => W-
      // positron+nu    -> -11+ 12    => -(-1)+0 = -1 => W+
      // u dbar -> 2 -1  => 0 -(-1) = 1 => W+
      // c dbar -> 4 -1  => W+
      // etc.
    }
    // Now we start adding the intermediate boson entry:
    int iup=nInt(*nup);
    if(iwch>0)      (*idup).push_back( 24);
    else if(iwch<0) (*idup).push_back(-24);
    else            (*idup).push_back( 23);
    (*istup).push_back(2);
    (*mothup1).push_back(1);
    (*mothup2).push_back(2);
    (*pup).push_back(vector<double>(5));
    double tmp = (*pup)[3][iup-2]+(*pup)[3][iup-1]; // Vector boson energy.
    (*pup)[3][iup] = tmp;
    tmp = sqr(tmp);
    for(unsigned int ixx=0; ixx<=2; ixx++) {
      (*pup)[ixx][iup] = (*pup)[ixx][iup-2]+(*pup)[ixx][iup-1];
      tmp = tmp-sqr((*pup)[ixx][iup]);  // Vector boson mass^2 when loop ends.
    }
    (*pup)[4][iup] = sqrt(tmp);  // Set vector boson mass.
    (*icolup1).push_back(0);   // Set 1st colour line for vector boson.
    (*icolup2).push_back(0);   // Set 2nd colour line for vector boson.
    (*nup) = (*nup)+1;         // Increment number of particles to be read in the event.
  }
  // nW + mZ + kH + jets ( vbjet* / ihrd=5 ):
  else if(ihrd==5) {
    unsigned int ivstart(0),ivend(0);
    // Identify the range of the Z and W bosons in the event (AlpGen conventions).
    // Note the Z's and W's are the only things that decay - Higgs's and photons 
    // do not.
    vector<unsigned int> bosonIndex;
    for(unsigned int ixx=0; ixx<(*nup); ixx++)
      if(abs((*idup)[ixx])==24||(*idup)[ixx]==23) {
	(*istup)[ixx]=2;             //  Set the W/Z boson status to intermediate.
	bosonIndex.push_back(ixx+1); //  W/Z boson index in LH event record (1->nup).
      }
    unsigned int bosonCounter(nInt(*nup)-2*bosonIndex.size());
    for(unsigned int ixx=0; ixx<bosonIndex.size(); ixx++) {
      (*mothup1)[bosonCounter]=bosonIndex[ixx];
      (*mothup2)[bosonCounter]=0;
      bosonCounter++;
      (*mothup1)[bosonCounter]=bosonIndex[ixx];
      (*mothup2)[bosonCounter]=0;
      bosonCounter++;
    }
  }
  // t tbar + jets [ + photons ] ( 2Q*, ihrd=6 [ 2Qph*, ihrd=16 ] ):
  else if ((ihrd==6||ihrd==16)&&abs((*idup)[2])==6) {
    // Redefine the tops as intermediates in the event record.
    (*istup)[2]=2;
    (*istup)[3]=2;
    unsigned int it(3),itb(4); // Index of t & tbar in evt.record (LH index).
    if((*idup)[2]!=6) swap(it,itb);
    // Reconstruct intermediate W's from the decay products.
    for(unsigned int ixx=0; ixx<4;ixx++) {
      (*idup).push_back(0);
      (*istup).push_back(0);
      (*mothup1).push_back(0);
      (*mothup2).push_back(0);
      (*icolup1).push_back(0);
      (*icolup2).push_back(0);
      (*pup).push_back(vector<double>(5));
    }
    for(unsigned int iw=1; iw<=2; iw++) {
      int iwdec(nInt(*nup)-5+2*iw);  // First of decay products for this W (iw).
      int iwup(nInt(*nup)+iw); // Where the reco. W will go in evt.record (LH index).
      int ibup(iwup+2);    // Where the reco. b will go in evt.record (under the Ws).
      int iwch(0);
      for(unsigned int iup=iwdec; iup<=iwdec+1; iup++) {
	(*mothup1)[iup-1]=iwup;
	(*mothup2)[iup-1]=0;
	iwch=iwch-int((*idup)[iup-1])%2; // iwch = charge of W boson.
	// electron+nubar -> 11 + (-12) = -1 => W-
	// d + ubar -> 1 + (-2) = -1 => W-
	// positron+nu    -> -11+ 12    =  1 => W+
	// u + dbar -> 2 + (-1) = 1 => W+
      }
      // Make space for the b and W:
      // Fill in b and W LH record entries:
      if(iwch>0) {
	(*idup)[iwup-1]=24;
	(*idup)[ibup-1]=5;
	(*mothup1)[iwup-1]=it;
	(*mothup2)[iwup-1]=0;
	(*mothup1)[ibup-1]=it;
	(*mothup2)[ibup-1]=0;
      } else if (iwch<0) {
	(*idup)[iwup-1]=-24;
	(*idup)[ibup-1]=-5;
	(*mothup1)[iwup-1]=itb;
	(*mothup2)[iwup-1]=0;
	(*mothup1)[ibup-1]=itb;
	(*mothup2)[ibup-1]=0;
      }
      (*istup)[iwup-1]=2; // The W is an intermediate and the
      (*istup)[ibup-1]=1; // b is a final-state particle.
      // Now reconstruct W boson momentum
      double tmp=(*pup)[3][iwdec-1]+(*pup)[3][iwdec];
      (*pup)[3][iwup-1]=tmp;            // W energy.
      tmp=sqr(tmp);
      for(unsigned int ixx=0; ixx<=2; ixx++) {
	(*pup)[ixx][iwup-1] = (*pup)[ixx][iwdec-1]            // W's 3-mom.
	                    + (*pup)[ixx][iwdec];
	tmp=tmp-sqr((*pup)[ixx][iwup-1]); // Equals m^2 at the end of the loop.
      }
      (*pup)[4][iwup-1]=sqrt(tmp);      // W mass.
      // Reconstruct b momentum
      int itmp(nInt((*mothup1)[iwup-1]));
      tmp=(*pup)[3][itmp-1]-(*pup)[3][iwup-1];
      (*pup)[3][ibup-1]=tmp;            // b energy.
      tmp=sqr(tmp);
      for(unsigned int ixx=0; ixx<=2; ixx++) {
	(*pup)[ixx][ibup-1] = (*pup)[ixx][nInt((*mothup1)[iwup-1])-1] // b's 3mom.
	                    - (*pup)[ixx][iwup-1];
	tmp=tmp-sqr((*pup)[ixx][ibup-1]); // Equals m^2 at the end of the loop.
      }
      (*pup)[4][ibup-1]=sqrt(tmp);      // b mass.
      (*icolup1)[iwup-1]=0;             // W has no colour
      (*icolup2)[iwup-1]=0;             // lines.
      (*icolup1)[ibup-1]=(*icolup1)[nInt((*mothup1)[iwup-1])-1]; // b shares top
      (*icolup2)[ibup-1]=(*icolup2)[nInt((*mothup1)[iwup-1])-1]; // colour line.
    }
    (*nup)+=4;
  }
  // H t tbar + jets ( QQh*, ihrd=8 ):
  else if (ihrd==8&&abs((*idup)[3])==6) {
    // Redefine the tops as intermediates in the event record.
    (*istup)[3]=2;
    (*istup)[4]=2;
    unsigned int it(4),itb(5); // Index of t & tbar in evt.record (LH index).
    if((*idup)[3]!=6) swap(it,itb);
    // Reconstruct intermediate W's from the decay products.
    for(unsigned int ixx=0; ixx<4;ixx++) {
      (*idup).push_back(0);
      (*istup).push_back(0);
      (*mothup1).push_back(0);
      (*mothup2).push_back(0);
      (*icolup1).push_back(0);
      (*icolup2).push_back(0);
      (*pup).push_back(vector<double>(5));
    }
    for(unsigned int iw=1; iw<=2; iw++) {
      int iwdec(nInt(*nup)-5+2*iw);  // First of decay products for this W (iw).
      int iwup(nInt(*nup)+iw); // Where the reco. W will go in evt.record (LH index).
      int ibup(iwup+2);    // Where the reco. b will go in evt.record (under the Ws).
      int iwch(0);
      for(unsigned int iup=iwdec; iup<=iwdec+1; iup++) {
	(*mothup1)[iup-1]=iwup;
	(*mothup2)[iup-1]=0;
	iwch=iwch-int((*idup)[iup-1])%2; // iwch = charge of W boson.
	// electron+nubar -> 11 + (-12) = -1 => W-
	// d + ubar -> 1 + (-2) = -1 => W-
	// positron+nu    -> -11+ 12    =  1 => W+
	// u + dbar -> 2 + (-1) = 1 => W+
      }
      if(iwch>0) {
	(*idup)[iwup-1]=24;
	(*idup)[ibup-1]=5;
	(*mothup1)[iwup-1]=it;
	(*mothup2)[iwup-1]=0;
	(*mothup1)[ibup-1]=it;
	(*mothup2)[ibup-1]=0;
      } else if (iwch<0) {
	(*idup)[iwup-1]=-24;
	(*idup)[ibup-1]=-5;
	(*mothup1)[iwup-1]=itb;
	(*mothup2)[iwup-1]=0;
	(*mothup1)[ibup-1]=itb;
	(*mothup2)[ibup-1]=0;
      }
      (*istup)[iwup-1]=2; // The W is an intermediate and the
      (*istup)[ibup-1]=1; // b is a final-state particle.
      // Now reconstruct W boson momentum
      double tmp=(*pup)[3][iwdec-1]+(*pup)[3][iwdec];
      (*pup)[3][iwup-1]=tmp;            // W energy.
      tmp=sqr(tmp);
      for(unsigned int ixx=0; ixx<=2; ixx++) {
	(*pup)[ixx][iwup-1] = (*pup)[ixx][iwdec-1]            // W's 3-mom.
	                    + (*pup)[ixx][iwdec];
	tmp=tmp-sqr((*pup)[ixx][iwup-1]); // Equals m^2 at the end of the loop.
      }
      (*pup)[4][iwup-1]=sqrt(tmp);      // W mass.
      // Reconstruct b momentum
      tmp=(*pup)[3][nInt((*mothup1)[iwup-1])-1]-(*pup)[3][iwup-1];
      (*pup)[3][ibup-1]=tmp;            // b energy.
      tmp=sqr(tmp);
      for(unsigned int ixx=0; ixx<=2; ixx++) {
	(*pup)[ixx][ibup-1] = (*pup)[ixx][nInt((*mothup1)[iwup-1])-1] // b's 3mom.
	                    - (*pup)[ixx][iwup-1];
	tmp=tmp-sqr((*pup)[ixx][ibup-1]); // Equals m^2 at the end of the loop.
      }
      (*pup)[4][ibup-1]=sqrt(tmp);      // b mass.
      (*icolup1)[iwup-1]=0;             // W has no colour
      (*icolup2)[iwup-1]=0;             // lines.
      (*icolup1)[ibup-1]=(*icolup1)[nInt((*mothup1)[iwup-1])-1]; // b shares top
      (*icolup2)[ibup-1]=(*icolup2)[nInt((*mothup1)[iwup-1])-1]; // colour line.
    }
    (*nup)+=4;
  }
  // Single top production ( top*, ihrd=13):
  else if (ihrd==13) {
    int nw=1;
    if(itopprc>=3) nw=2;
    // Assign a mass to the incoming bottom quark, if there is one,
    // rescaling the energy to accommodate the mass.
    for(unsigned int ixx=1; ixx<=2; ixx++) {
      if(abs((*idup)[ixx-1])==5) {
	(*pup)[4][ixx-1]=masses[1];
	(*pup)[3][ixx-1]=sqrt(sqr((*pup)[2][ixx-1])+sqr((*pup)[4][ixx-1]));
      }
    }
    // Set the top status to that of an intermediate.
    (*istup)[2]=2;
    // Get the index of the t / tbar in the evt. record (LH convention: 1->nup).
    unsigned int it=0;
    unsigned itb=0;
    if((*idup)[2]==6)       it=3;
    else if((*idup)[2]==-6) itb=3;
    else {
      cout << "doIndividualHardProcessAssignments:\n"
	   << "wrong assumption about top position.\n"
	   << "Quitting ...";
      exit(1);
    }
    // Reconstruct intermediate W's from the decay products.
    // iwdec is the index of the first W decay product 
    unsigned int iwdec(nInt(*nup)-1);  // LH conventions: 1->nup.
    if(nw==2) iwdec=nInt(*nup)-3;
    // The W and b will go at the end of the record.
    unsigned int iwup(nInt(*nup)+1);
    unsigned int ibup(nInt(*nup)+2);
    int iwch(0);
    for(unsigned int iup=iwdec; iup<=iwdec+1; iup++) {
      (*mothup1)[iup-1]=iwup;
      (*mothup2)[iup-1]=0;
      iwch=iwch-int((*idup)[iup-1])%2; // iwch = charge of W boson.
      // electron+nubar -> 11 + (-12) = -1 => W-
      // d + ubar -> 1 + (-2) = -1 => W-
      // positron+nu    -> -11+ 12    =  1 => W+
      // u + dbar -> 2 + (-1) = 1 => W+
    }
    for(unsigned int ixx=0; ixx<2;ixx++) {
      (*idup).push_back(0);
      (*istup).push_back(0);
      (*mothup1).push_back(0);
      (*mothup2).push_back(0);
      (*icolup1).push_back(0);
      (*icolup2).push_back(0);
      (*pup).push_back(vector<double>(5));
    }
    if(iwch>0) {
      (*idup)[iwup-1]=24;
      (*idup)[ibup-1]=5;
      (*mothup1)[iwup-1]=it;
      (*mothup2)[iwup-1]=0;
      (*mothup1)[ibup-1]=it;
      (*mothup2)[ibup-1]=0;
    } else if (iwch<0) {
      (*idup)[iwup-1]=-24;
      (*idup)[ibup-1]=-5;
      (*mothup1)[iwup-1]=itb;
      (*mothup2)[iwup-1]=0;
      (*mothup1)[ibup-1]=itb;
      (*mothup2)[ibup-1]=0;
    }
    (*istup)[iwup-1]=2;
    (*istup)[ibup-1]=1;
    // Now reconstruct W boson momentum
    double tmp=(*pup)[3][iwdec-1]+(*pup)[3][iwdec];
    (*pup)[3][iwup-1]=tmp;            // W energy.
    tmp=sqr(tmp);
    for(unsigned int ixx=0; ixx<=2; ixx++) {
      (*pup)[ixx][iwup-1] = (*pup)[ixx][iwdec-1]            // W's 3-mom.
	+ (*pup)[ixx][iwdec];
      tmp=tmp-sqr((*pup)[ixx][iwup-1]); // Equals m^2 at the end of the loop.
    }
    (*pup)[4][iwup-1]=sqrt(tmp);      // W mass.
    // Reconstruct b momentum
    tmp=(*pup)[3][nInt((*mothup1)[iwup-1])-1]-(*pup)[3][iwup-1];
    (*pup)[3][ibup-1]=tmp;            // b energy.
    tmp=sqr(tmp);
    for(unsigned int ixx=0; ixx<=2; ixx++) {
      (*pup)[ixx][ibup-1] = (*pup)[ixx][nInt((*mothup1)[iwup-1])-1] // b's 3mom.
	                  - (*pup)[ixx][iwup-1];
      tmp=tmp-sqr((*pup)[ixx][ibup-1]); // Equals m^2 at the end of the loop.
    }
    (*pup)[4][ibup-1]=sqrt(tmp);      // b mass.
    (*icolup1)[iwup-1]=0;
    (*icolup2)[iwup-1]=0;
    (*icolup1)[ibup-1]=(*icolup1)[nInt((*mothup1)[iwup-1])-1];
    (*icolup2)[ibup-1]=(*icolup2)[nInt((*mothup1)[iwup-1])-1];
    (*nup)=(*nup)+2;
    if(nw==2) {
      // W DECAY
      // iwdec is the index of the first W decay product.
      iwdec=nInt(*nup)-3;
      // iwup is the index of the W in the event record (LH conventions: 1->nup).
      iwup=nInt(*nup)-6;
      iwch=0;
      int iwch(0);
      for(unsigned int iup=iwdec; iup<=iwdec+1; iup++) {
	(*mothup1)[iup-1]=iwup;
	(*mothup2)[iup-1]=0;
	iwch=iwch-int((*idup)[iup-1])%2; // iwch = charge of W boson.
	// electron+nubar -> 11 + (-12) = -1 => W-
	// d + ubar -> 1 + (-2) = -1 => W-
	// positron+nu    -> -11+ 12    =  1 => W+
	// u + dbar -> 2 + (-1) = 1 => W+
      }
      (*istup)[iwup-1]=2;
      (*icolup1)[iwup-1]=0;
      (*icolup2)[iwup-1]=0;
    }
  }

  return;
}
