#include <cstdlib> 
#include <ctime> 
#include <iostream>
#include <string>  
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;


int main (int argc, char *argv[]) { 
  string eventfile = "";
  string numevents = "1";
  if(argv[1] && argv[2]) { eventfile = argv[1]; numevents = argv[2];} else { cout << "Use: run_hpp [eventfile] [numevents]" << endl; exit(1); }
 
  cout << "Running Herwig++ with eventfile: " << eventfile << " and " << numevents << " events\n";
  //  ifstream stream;
  //ofstream outstream;
  //ofstream outorigin;
  //string allfile = "";
  // string addstr = "";
  
  // stream.open("wpnlo.in");
  // outstream.open("wpnlo.temp");
  //outorigin.open("wpnlo.orig");
  // if(!stream) { cerr << "Error: Failed to open file" << endl;}
  
  // string search_string = "eventfile.dat";
  //string replace_string = eventfile;
  //string inbuf;
 
  //  while(!stream.eof())
  //{
  //     getline(stream, inbuf);
  // int spot = inbuf.find(search_string);
      //    outorigin << inbuf << endl;
  //    if(spot >= 0)
  //    {
  //       string tmpstring = inbuf.substr(0,spot);
  //       tmpstring += replace_string;
   //      tmpstring += inbuf.substr(spot+search_string.length(), inbuf.length());
   //      inbuf = tmpstring;
  //      }
  //	 outstream << inbuf << endl;
  //
  // }
  //outstream.close();
  //stream.close();

  string runcomm1 = "mv ";
  runcomm1.append(eventfile);
  runcomm1.append(" EPEM.dat");
  system(runcomm1.c_str());
  cout << "Reading EPEM.in..." << endl;
  system("Herwig++ read EPEM.in");
  cout << "Done..." << endl;
  cout << "Running MG.run..." << endl;
  string runcomm2 = "Herwig++ run MG.run -N";
  cout << "Done!" << endl;
  runcomm2.append(numevents);
  system(runcomm2.c_str());
  string runcomm3 = "mv EPEM.dat ";
  runcomm3.append(eventfile);
  system(runcomm3.c_str());
  // system("mv wpnlo.orig wpnlo.in");
}
