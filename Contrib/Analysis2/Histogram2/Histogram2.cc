// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Histogram2 class.
//

#include "Histogram2.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Histogram2.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Utilities/Throw.h"

using namespace Analysis2;

void HistogramChannel::persistentOutput(PersistentOStream & os) const {
  os << _isCountingChannel << _bins << _binEntries << _outOfRange << _visible
     << _total << _nanEvents << _nanWeights << _finished;
}

void HistogramChannel::persistentInput(PersistentIStream & is) {
  is >> _isCountingChannel >> _bins >> _binEntries >> _outOfRange >> _visible
     >> _total >> _nanEvents >> _nanWeights >> _finished;
}

HistogramChannel& HistogramChannel::operator += (const HistogramChannel& c) {
  if (!c.isCountingChannel()) _isCountingChannel = false;
  for (unsigned int i = 0; i < _bins.size(); ++i) {
    _bins[i].first += c.bin(i).first;
    _bins[i].second += c.bin(i).second;
  }
  // uppon addition of a counting channel it remains a counting channel
  if (_isCountingChannel) {
    for (unsigned int i = 0; i < _bins.size(); ++i) {
      _binEntries[i] += c.binEntries()[i];
      _nanWeights[i] += c.nanWeights()[i];
    }
    _outOfRange.first += c.outOfRange().first;
    _outOfRange.second += c.outOfRange().second;
    _total += c.total();
    _visible += c.visible();
    _nanEvents += c.nanEvents();
  }
  return *this;
}

HistogramChannel& HistogramChannel::operator -= (const HistogramChannel& c) {
  for (unsigned int i = 0; i < _bins.size(); ++i) {
    _bins[i].first -= c.bin(i).first;
    _bins[i].second += c.bin(i).second;
  }
  _isCountingChannel = false;
  return *this;
}

HistogramChannel& HistogramChannel::operator *= (const HistogramChannel& c) {
  for (unsigned int i = 0; i < _bins.size(); ++i) {
    _bins[i].first *= c.bin(i).first;
    _bins[i].second = 
      _bins[i].second*sqr(c.bin(i).first)+c.bin(i).second*sqr(_bins[i].first);
  }
  _isCountingChannel = false;
  return *this;
}

HistogramChannel& HistogramChannel::operator *= (double factor) {
  for (vector<pair<double,double> >::iterator b = _bins.begin();
       b != _bins.end(); ++b) {
    b->first *= factor;
    b->second *= sqr(factor);
  }
  _isCountingChannel = false;
  return *this;
}

HistogramChannel& HistogramChannel::operator += (double off) {
  for (vector<pair<double,double> >::iterator b = _bins.begin();
       b != _bins.end(); ++b) {
    b->first += off;
  }
  _isCountingChannel = false;
  return *this;
}

HistogramChannel& HistogramChannel::operator *= (pair<double,double> factor) {
  for (vector<pair<double,double> >::iterator b = _bins.begin();
       b != _bins.end(); ++b) {
    b->first *= factor.first;
    b->second = b->first*factor.second + b->second*sqr(factor.first);
  }
  _isCountingChannel = false;
  return *this;
}

HistogramChannel& HistogramChannel::operator /= (pair<double,double> factor) {
  for (vector<pair<double,double> >::iterator b = _bins.begin();
       b != _bins.end(); ++b) {
    if (factor.first != 0.) {
      b->first /= factor.first;
      if (b->first != 0.)
	b->second = sqr(b->first/factor.first)*(b->second/sqr(b->first) + factor.second/sqr(factor.first));
      else
	b->second = 0.;
    } else {
      b->first = 0.;
      b->second = 0.;
    }
  }
  _isCountingChannel = false;
  return *this;
}

HistogramChannel& HistogramChannel::operator /= (const HistogramChannel& c) {
  for (unsigned int i = 0; i < _bins.size(); ++i) {
    if (c.bin(i).first != 0. && _bins[i].first != 0.) {
      _bins[i].first /= c.bin(i).first;
      _bins[i].second = sqr(_bins[i].first/c.bin(i).first)*
	(_bins[i].second/sqr(_bins[i].first)+c.bin(i).second/sqr(c.bin(i).first));
    } else {
      _bins[i].first = 0.;
      _bins[i].second = 0.;
    }
  }
  _isCountingChannel = false;
  return *this;
}

unsigned long HistogramChannel::nanWeightEvents () const {
  unsigned long all = 0;
  for (vector<unsigned long>::const_iterator n = _nanWeights.begin();
       n != _nanWeights.end(); ++n)
    all += *n;
  return all;
}

void HistogramChannel::differential (const vector<pair<double,double> >& binning) {
  for (unsigned int i=0; i<binning.size(); ++i) {
    _bins[i].first /= (binning[i].second-binning[i].first);
    _bins[i].second /= (binning[i].second-binning[i].first);
  }
}

pair<double,double> HistogramChannel::binSum () const {
  pair<double,double> s = make_pair(0.,0.);
  for (vector<pair<double,double> >::const_iterator b = _bins.begin();
       b != _bins.end(); ++b) {
    s.first += b->first;
    s.second += b->second;
  }
  return s;
}

pair<double,double> HistogramChannel::integrate (const vector<pair<double,double> >& binning) const {
  pair<double,double> integral = make_pair(0.,0.);
  for (unsigned int i = 0; i< _bins.size(); ++i) {
    integral.first += (binning[i].second-binning[i].first)*_bins[i].first;
    integral.second += sqr(binning[i].second-binning[i].first)*_bins[i].second;
  }
  return integral;
}

pair<double,double> HistogramChannel::average (const vector<pair<double,double> >& binning) const {
  pair<double,double> integral = make_pair(0.,0.);
  double volume = 0.;
  for (unsigned int i = 0; i< _bins.size(); ++i) {
    volume += (binning[i].second-binning[i].first);
    integral.first += (binning[i].second-binning[i].first)*_bins[i].first;
    integral.second += sqr(binning[i].second-binning[i].first)*_bins[i].second;
  }
  integral.first /= volume;
  integral.second /= sqr(volume);
  return integral;
}

HistogramChannel HistogramChannel::delta (const HistogramChannel& channel) const {
  HistogramChannel c(*this);
  c /= channel; c += -1.;
  return c;
}

HistogramChannel HistogramChannel::chi2 (const HistogramChannel& channel, double minfrac) const {
  HistogramChannel chi2 (*this);
  chi2 -= channel;
  chi2 *= chi2;
  for (unsigned int i = 0; i<_bins.size(); ++i) {
    double var = 0.;
    if (channel.bin(i).second/sqr(channel.bin(i).first) < minfrac)
      var = sqr(minfrac*channel.bin(i).first);
    else var = channel.bin(i).second;
    if (var != 0)
      chi2.bin(i,make_pair(chi2.bin(i).first/var,chi2.bin(i).second/var));
    else
      chi2.bin(i,make_pair(0,0));
  }
  return chi2;
}


Histogram2::Histogram2 (double low, double high,
			unsigned int bins,
			const string& name) {

  // get the bin length
  double length = (high-low)/(double)(bins);

  _range = make_pair(low,high);

  for (unsigned int i = 0; i<bins; ++i) {
    _binning.push_back(make_pair(low+i*length,low+(i+1)*length));
    _binhash.insert(make_pair(low+(i+1)*length,i));
  }

  if (!name.empty()) insertChannel(name);

}

Histogram2::Histogram2 (const vector<pair<double,double> >& binning, const string& name) {

  _range = make_pair(binning.begin()->first,binning.end()->second);
  _binning = binning;

  for (unsigned int i = 0; i<_binning.size(); ++i) {
    _binhash.insert(make_pair(_binning[i].second,i));
  }

  _range = make_pair(_binning.front().first,_binning.back().second);

  if (!name.empty()) insertChannel(name);

}

Histogram2::Histogram2 (const string& dataFile, const string& dataName) {
  ifstream data (dataFile.c_str());
  if (!data) {
    Throw<InitException>() << "Histogram2::Histogram2 : Building from datafile, but cannot open "
			   << dataFile;
  }
  vector<pair<double,double> > dataCache;
  double low, high, dataval, errstat, errsys;
  double sigma2;
  string in;
  while (getline(data,in)) {
    in = StringUtils::stripws(in);
    if (in[0] == '#') continue;
    if (in == "") continue;
    istringstream theIn (in);
    theIn >> low >> high >> dataval >> errstat >> errsys;
    _binning.push_back(make_pair(low,high));
    sigma2 = sqr(errstat) + sqr(errsys);
    dataCache.push_back(make_pair(dataval,sigma2));
  }
  for (unsigned int i = 0; i<_binning.size(); ++i) {
    _binhash.insert(make_pair(_binning[i].second,i));
  }
  _range = make_pair(_binning.front().first,_binning.back().second);
  HistogramChannel theData (dataCache);
  insertChannel(dataName,theData);

}

void Histogram2::book (const string& name, double event, double weight) {
  map<string,HistogramChannel>::iterator c = _channels.find(name);
  if (c != _channels.end()) {
    if (isnan(event) || isinf(event)) c->second.nanEvent();
    else {
      if (event < range().first) {
	c->second.bookUnderflow(weight);
	return;
      }
      if (event > range().second) {
	c->second.bookOverflow(weight);
	return;
      }
      unsigned int bin = _binhash.upper_bound(event)->second;
      c->second.book(bin,weight);
    }
  }
}

vector<string> Histogram2::channels () const {
  vector<string> all;
  for (map<string,HistogramChannel>::const_iterator c = _channels.begin();
       c != _channels.end(); ++c)
    all.push_back(c->first);
  return all;
}

void Histogram2::output (ostream& os, const string& name,
			 unsigned int flags, char comment) const {
  bool bincenters = (flags & ChannelOutput::Bincenters) == ChannelOutput::Bincenters;
  bool noerrorbars = (flags & ChannelOutput::NoErrorbars) == ChannelOutput::NoErrorbars;
  bool nanevents = (flags & ChannelOutput::NanEvents) == ChannelOutput::NanEvents;
  bool statistics = (flags & ChannelOutput::Statistics) == ChannelOutput::Statistics;

  map<string,HistogramChannel>::const_iterator c = _channels.find(name);
  if (c == _channels.end()) return;

  os << comment << " channel = " << name << endl;

  if(statistics) {
    os << comment << " total entries = " << c->second.total()
       << " , visible entries = " << c->second.visible() << endl;
    if (c->second.isCountingChannel())
      os << comment << " underflow = " << c->second.outOfRange().first
	 << " , overflow = " << c->second.outOfRange().second << endl;
    if (c->second.nanEvents() > 0)
      os << comment << " nan events = " << c->second.nanEvents();
    if (c->second.nanWeightEvents() > 0)
      os << " , nan weight events = " << c->second.nanWeightEvents() << endl;
  }

  for (unsigned int i = 0; i<_binning.size(); ++i) {
    if (bincenters) {
      os << (_binning[i].first+_binning[i].second)/2. << "\t";
    } else {
      os << _binning[i].first << "\t" << _binning[i].second << "\t";
    }
    os << c->second.bin(i).first;
    if (!noerrorbars) {
      os << "\t" << sqrt(c->second.bin(i).second);
    }
    if (nanevents) {
      os << "\t" << c->second.nanWeights()[i];
    }
    os << endl;
  }

}

HistogramChannel HistogramChannel::profile () const {
  HistogramChannel temp (_bins.size(),false);
  for (unsigned int i = 0; i< _bins.size(); ++i) {
    pair<double,double> prof;
    // mean of weights
    prof.first = weightMean(i);
    // varaince of mean
    prof.second = binEntries(i) > 0 ? weightVariance(i)/binEntries(i) : 0;
    temp.bin(i,prof,binEntries(i));
  }
  return temp;
}

void HistogramChannel::write (ostream& os, const string& name) {
  finish();
  os << "<channel"
     << " name=\"" << name << "\""
     << " counting=\"" << _isCountingChannel << "\"";

  if (_isCountingChannel) {
     os << " underflow=\"" << _outOfRange.first << "\""
	<< " overflow=\"" << _outOfRange.second << "\""
	<< " visible=\"" << _visible << "\""
	<< " total=\"" << _total << "\"";
  }

  if (_nanEvents != 0) os << " nanevents=\"" << _nanEvents << "\"";

  os << ">" << endl;

  os << "<bins>" << endl;
  for (unsigned int b = 0; b < _bins.size(); ++b) {
    os << "<bincontent sumweights=\"" << _bins[b].first
       << "\" sumsquaredweights=\"" << _bins[b].second
       << "\" entries=\"" << _binEntries[b]
       << "\"/>" << endl;
  }
  os << "</bins>" << endl;

  os << "<nanweights>" << endl;

  // only write out, if at least happend for one event
  if (nanWeightEvents() > 0)
    for (vector<unsigned long>::const_iterator n = _nanWeights.begin();
	 n != _nanWeights.end(); ++n)
      os << "<bincontent entries=\"" << *n << "\"/>" << endl;

  os << "</nanweights>" << endl;

  os << "</channel>" << endl;
}

string HistogramChannel::read (istream& is) {

  _finished = true;

  string tag;
  string name;
  map<string,string> attributes;
  map<string,string>::iterator atit;

  if (!is) return "";

  // parse the channel tag

  tag = getNextTag(is);

  if (tag.find("<channel") == string::npos) return "";

  attributes = StringUtils::xmlAttributes("channel",tag);

  atit = attributes.find("name"); if (atit == attributes.end()) return "";
  name = atit->second;

  atit = attributes.find("counting"); if (atit == attributes.end()) return "";
  fromString(atit->second,_isCountingChannel);

  atit = attributes.find("underflow"); if (atit == attributes.end() && _isCountingChannel) return "";
  if (_isCountingChannel) fromString(atit->second,_outOfRange.first);

  atit = attributes.find("overflow"); if (atit == attributes.end() && _isCountingChannel) return "";
  if (_isCountingChannel) fromString(atit->second,_outOfRange.second);

  atit = attributes.find("visible"); if (atit == attributes.end() && _isCountingChannel) return "";
  if (_isCountingChannel) fromString(atit->second,_visible);

  atit = attributes.find("total"); if (atit == attributes.end() && _isCountingChannel) return "";
  if (_isCountingChannel) fromString(atit->second,_total);

  atit = attributes.find("nanevents");
  if (atit != attributes.end()) fromString(atit->second,_nanEvents);

  // read in the bin contents

  tag = getNextTag(is);

  if (tag != "<bins>") return "";

  for (unsigned int i = 0; i < _bins.size(); ++i) {
    tag = getNextTag(is);
    if (tag.find("<bincontent") == string::npos) return "";
    attributes = StringUtils::xmlAttributes("bincontent",tag);
    atit = attributes.find("sumweights"); if (atit == attributes.end()) return "";
    fromString(atit->second,_bins[i].first);
    atit = attributes.find("sumsquaredweights"); if (atit == attributes.end()) return "";
    fromString(atit->second,_bins[i].second);
    atit = attributes.find("entries"); if (atit == attributes.end()) return "";
    fromString(atit->second,_binEntries[i]);
  }

  tag = getNextTag(is);

  if (tag != "</bins>") return "";

  // read in the nan weight histogram

  tag = getNextTag(is);

  if (tag != "<nanweights>") return "";

  bool nanweights = true;

  for (unsigned int i = 0; i < _bins.size(); ++i) {
    tag = getNextTag(is);
    if (tag == "</nanweights>") {
      nanweights = false;
      break;
    }
    if (tag.find("<bincontent") == string::npos) return "";
    attributes = StringUtils::xmlAttributes("bincontent",tag);
    atit = attributes.find("entries"); if (atit == attributes.end()) return "";
    fromString(atit->second,_nanWeights[i]);
  }

  if (nanweights) tag = getNextTag(is); if (tag != "</nanweights>") return "";

  tag = getNextTag(is); 

  if (tag != "</channel>") return "";

  return name;

}

void Histogram2::store (const string& name) {

  ofstream os ((name+".h2").c_str());
  if (!os) return;

  os << "<?xml version=\"1.0\"?>" << endl;
  os << "<Analysis2Histogram version=\"1.0\"" << " AnalysisName=\"" << Named::name() << "\">" << endl;
  os << "<!--" << endl
     << "  WARNING" << endl
     << "  Though this is valid XML, the Histogram2 class will" << endl
     << "  not be able to parse arbitraty, XML-valid changes" << endl
     << "  to this file!" << endl
     << "-->" << endl;

  os << "<xsec unit=\"nanobarn\" value=\"" << _xSec/nanobarn << "\"/>" << endl;

  os << "<binning>" << endl;

  for (vector<pair<double,double> >::const_iterator b = _binning.begin();
       b != _binning.end(); ++b) {
    os << "<bin lower=\"" << b->first << "\" upper=\"" << b->second << "\"/>" << endl;
  }

  os << "</binning>" << endl;

  os << "<channels size=\"" << _channels.size()  << "\">" << endl;

  for(map<string,HistogramChannel>::iterator c = _channels.begin();
      c != _channels.end(); ++c) {
    c->second.write(os,c->first);
  }

  os << "</channels>" << endl;

  os << "</Analysis2Histogram>" << endl;

}

bool Histogram2::load (const string& fname) {

  ifstream is ((fname+".h2").c_str());
  if (!is) return false;

  string tag = getNextTag(is);
  if (tag.find("<Analysis2Histogram") == string::npos) return false;
  string name;
  map<string,string> attributes = StringUtils::xmlAttributes("Analysis2Histogram",tag);
  map<string,string>::iterator atit = attributes.find("name"); if (atit == attributes.end()) return false;
  Named::name(atit->second);

  // get the cross section
  tag = getNextTag(is);

  if (tag.find("<xsec") == string::npos) return false;
  attributes = StringUtils::xmlAttributes("xsec",tag);
  atit = attributes.find("unit"); if (atit == attributes.end()) return false;
  if (atit->second != "nanobarn") return false; // switch units in a future version
  atit = attributes.find("value"); if (atit == attributes.end()) return false;
  double theXsec; fromString(atit->second,theXsec);
  _xSec = theXsec * nanobarn;

  // get the binning

  tag = getNextTag(is);

  if (tag != "<binning>") return false;

  tag = getNextTag(is);

  double binlower, binupper;

  while (tag != "</binning>") {
    if (tag.find("<bin") == string::npos && tag != "</binning>") return false;
    attributes = StringUtils::xmlAttributes("bin",tag);
    atit = attributes.find("lower"); if (atit == attributes.end()) return false;
    fromString(atit->second,binlower);
    atit = attributes.find("upper"); if (atit == attributes.end()) return false;
    fromString(atit->second,binupper);
    _binning.push_back(make_pair(binlower,binupper));
    tag = getNextTag(is);
  }

  if (!_binning.size()) return false;
  for (unsigned int i = 0; i<_binning.size(); ++i) {
    _binhash.insert(make_pair(_binning[i].second,i));
  }

  _range = make_pair(_binning.front().first,_binning.back().second);

  // get the channels

  tag = getNextTag(is); 

  if (tag.find("<channels") == string::npos) return false;
  attributes = StringUtils::xmlAttributes("channels",tag);
  atit = attributes.find("size"); if (atit == attributes.end()) return false;

  unsigned int numchannels;
  fromString(atit->second,numchannels);

  for (unsigned int i = 0; i<numchannels; ++i) {
    HistogramChannel ch (_binning.size());
    name = ch.read(is);
    if (name != "") _channels.insert(make_pair(name,ch));
  }

  // there has to be at least one channel
  if (!_channels.size()) return false;

  tag = getNextTag(is);
  if (tag != "</channels>") return false;

  tag = getNextTag(is);
  if (tag != "</Analysis2Histogram>") return false;

  return true;

}

Histogram2Ptr Histogram2::loadToHistogram (const string& name) const {
  Histogram2Ptr histo = new_ptr(Histogram2());
  bool ok = histo->load(name);
  if (!ok) histo = Histogram2Ptr();
  return histo;
}

void Histogram2::combine (const string& prefix, const string& name,
			  unsigned int numRuns, const string& dataChannel,
			  const string& mcChannel) {
  vector<Histogram2Ptr> inHistos;
  for (unsigned int i = 0; i<numRuns; ++i) {
    ostringstream fname ("");
    fname << prefix << "." << i << "/" << name;
    Histogram2Ptr in = loadToHistogram (fname.str());
    if (in) {
      inHistos.push_back(in);
    }
  }

  // get the total cross section
  CrossSection all = 0.*nanobarn;
  unsigned long allEvents = 0;
  for(vector<Histogram2Ptr>::iterator h = inHistos.begin();
      h != inHistos.end(); ++h) {
    if ((**h).haveChannel(mcChannel)) {
      all += (**h).xSec() * (**h).channel(mcChannel).total();
      allEvents += (**h).channel(mcChannel).total();
    }
  }
  all /= allEvents;
  xSec (all);

  if (!inHistos.size()) return;
  vector<string> channels = inHistos[0]->channels();
  _binning = inHistos[0]->binning();
  _range = inHistos[0]->range();
  _binhash = inHistos[0]->binhash();

  unsigned int numBins = _binning.size();

  for (vector<string>::iterator c = channels.begin();
       c != channels.end(); ++c) {
    if (*c == dataChannel && inHistos[0]->haveChannel(dataChannel)) {
      _channels.insert(make_pair(dataChannel,HistogramChannel(inHistos[0]->channel(dataChannel))));
      continue;
    }
    HistogramChannel ch (numBins);
    for (vector<Histogram2Ptr>::iterator h = inHistos.begin();
	 h != inHistos.end(); ++h) {
      if ((**h).haveChannel(*c)) {
	ch += (**h).channel(*c);
      }
    }
    _channels.insert(make_pair(*c,ch));
  }

}

Histogram2::Histogram2 (vector<double> limits) {
  for (unsigned int i = 0; i < limits.size()-1; ++i)
    _binning.push_back(make_pair(limits[i],limits[i+1]));

  for (unsigned int i = 0; i<_binning.size(); ++i) {
    _binhash.insert(make_pair(_binning[i].second,i));
  }

  _range = make_pair(_binning.front().first,_binning.back().second);

  insertChannel("mc");
}

Histogram2::Histogram2 (vector<double> limits,
			vector<double> data,
			vector<double> dataerror) {
  vector<pair<double,double> > dataCache;
  for (unsigned int i = 0; i < limits.size(); ++i) {
    if (i<limits.size()-1)
      _binning.push_back(make_pair(limits[i],limits[i+1]));
    dataCache.push_back(make_pair(data[i],sqr(dataerror[i])));
  }

  for (unsigned int i = 0; i<_binning.size(); ++i) {
    _binhash.insert(make_pair(_binning[i].second,i));
  }

  _range = make_pair(_binning.front().first,_binning.back().second);

  insertChannel("mc");
  insertChannel("data",HistogramChannel(dataCache));

}

Histogram2 Histogram2::ratioWith(const Histogram2& h2) const {
  Histogram2 tmp (*this);
  tmp.channel("mc") /= h2.channel("mc");
  return tmp;
}

Histogram2::~Histogram2() {}

void Histogram2::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _binning << _range << _binhash << _channels << ounit(_xSec,nanobarn);
}

void Histogram2::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _binning >> _range >> _binhash >> _channels >> iunit(_xSec,nanobarn);
}

ClassDescription<Histogram2> Histogram2::initHistogram2;
// Definition of the static class description member.

void Histogram2::Init() {

  static ClassDocumentation<Histogram2> documentation
    ("Histogram2 provides advanced histogramming.");
}

