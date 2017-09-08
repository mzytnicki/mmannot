/*
   mmannot, a multi-mapping annotation tool.
   Copyright (C) 2017 Matthias Zytnicki

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <queue>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <thread> 
#include <mutex> 
#include <atomic> 
#include <functional> 
#include <utility>
#include <regex>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <locale>
#include <libgen.h>
#include <zlib.h>
using namespace std;

class Interval;
class Read;
class Config;

static const char BAM_CIGAR_LOOKUP[] = "MIDNSHP=X";
static const char BAM_DNA_LOOKUP[] = "=ACMGRSVTWYHKDBN";

enum class Strand       { ALL, F, R };
enum class Strategy     { DEFAULT, UNIQUE, RANDOM, RATIO };
enum class Strandedness { U, F, R, FF, FR, RF };
enum class ReadsFormat  { UNKNOWN, SAM, BAM };

typedef unsigned long Position;
typedef bool(*StrandednessFunction)(bool, bool, bool);
typedef size_t(*IntervalOverlapFunction)(Interval &, Interval &);
typedef void(*PrintReadStatsFunction)(const string &, unsigned int, vector < size_t > &);
typedef bool(*RescueFunction)(vector < size_t > &);

const string defaultConfigFileName { "config.txt" };
const Position UNKNOWN = numeric_limits<Position>::max();
const unsigned int NO_CHR = numeric_limits<unsigned int>::max();
const size_t NOT_KEPT = numeric_limits<size_t>::max();
const size_t NO_ID = numeric_limits<size_t>::max();
const static string EMPTY_STRING { "" };
const static Position binSize = 16384;

bool sorted;
float overlap;
vector <Strandedness> strandednesses;
vector <StrandednessFunction> strandednessFunctions;
vector <ReadsFormat> formats;
IntervalOverlapFunction intervalOverlapFunction;
PrintReadStatsFunction printReadStatsFunction;
RescueFunction rescueFunction;
Config *config;
float rescueThreshold;
Strategy strategy;
unsigned int upstreamSize    = 1000;
unsigned int downstreamSize  = 1000;
unsigned int nInputs;
unsigned int nThreads = 1;
mutex printMutex;
bool progress = false;
bool readStats = false;
bool intervalStats = false;
ofstream readStatsFile;
ofstream intervalStatsFile;


ostream &operator<< (ostream &os, const Strand &s) {
  if      (s == Strand::F) os << "(+)";
  else if (s == Strand::R) os << "(-)";
  return os;
}


namespace std {
  template <>
    struct hash<vector<size_t>> {
      size_t operator() (const vector<size_t> &vc) const {
        size_t seed = 0;
        for(auto& i : vc) seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
      }
    };
}

class comma_numpunct : public std::numpunct<char> {
  protected:
    virtual char   do_thousands_sep () const { return ','; }
    virtual string do_grouping ()      const { return "\03"; }
};


template<typename T>
void join(vector <T> &inStrings, string &outString, const char *delim) {
  typename vector <T>::iterator it = inStrings.begin();
  stringstream ss;
  ss << *it;
  for (++it; it != inStrings.end(); it++) ss << delim << *it;
  outString = ss.str();
}

static inline void ltrim(string &s) {
  s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
}

static inline void rtrim(string &s) {
  s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
}

static inline void trim(string &s) {
  ltrim(s); rtrim(s);
}

static inline void printStats(unsigned int n, const char *s, unsigned int nHits) {
  unsigned int size = log10(static_cast<double>(nHits))+1;
  size += static_cast<unsigned int>(size / 3.0);
  cerr << "\t" << s << setw(size) << n << " (" << setw(5) << fixed << setprecision(1) << (static_cast<float>(n)/nHits*100) << "%)\n";
}

static inline void lower(string &s) {
  transform(s.begin(), s.end(), s.begin(), ::tolower);
}

static inline void split(string &inString, char delim, vector <string> &outStrings) {
  stringstream ss(inString);
  string item;
  while (getline(ss, item, delim)) {
    outStrings.push_back(item);
  }
}

static inline bool split(string &sin, char delim, string &sout1, string &sout2) {
  size_t pos = sin.find(delim);
  if (pos == string::npos) {
    return false;
  }
  sout1 = sin.substr(0, pos);
  sout2 = sin.substr(pos+1);
  trim(sout1);
  trim(sout2);
  return true;
}

static inline string rtrimTo(const string &sin, char delim) {
  size_t pos = sin.find(delim);
  if (pos == string::npos) {
    return sin;
  }
  string sout = sin.substr(0, pos);
  return sout;
}


template<typename T>
void setUnique (T &t) {
  sort(t.begin(), t.end());
  typename T::iterator it = unique(t.begin(), t.end());
  t.resize(distance(t.begin(), it)); 
}

struct RegionType {
  size_t id, subId;
  Strand strand;
  RegionType (size_t i = NO_ID, size_t si = NO_ID, Strand s = Strand::ALL): id(i), subId(si), strand(s) { }
  bool isSet () const { return (id != NO_ID); }
  friend bool operator<(const RegionType &r1, const RegionType &r2) {
    return ((r1.id < r2.id) || ((r1.id == r2.id) && (r1.subId < r2.subId)));
  }
};
ostream &operator<< (ostream &os, const RegionType &r) {
  os << r.id << ":" << r.subId << " " << r.strand;
  return os;
}
static const RegionType NO_TYPE (NO_ID, NO_ID);


struct AnnotationType {
  string source;
  regex sourceRegex;
  string type;
  Strand strand;
  AnnotationType (const string &s, const regex &r, const string &t, Strand st): source(s), sourceRegex(r), type(t), strand(st) {}
};
ostream &operator<< (ostream &os, const AnnotationType &a) {
  os << a.source << ":" << a.type << " " << a.strand;
  return os;
}

class Config {
  protected:
    vector < pair < regex, string > > synonyms;
    vector < vector < AnnotationType > > order;
    vector < tuple < string, string, size_t > > introns;
    vector < tuple < string, string, size_t, size_t > > vicinity;
    vector < RegionType > elements;
    size_t nElements = {0};

    static size_t check (const string &s1, const string &s2, vector < AnnotationType > &m) {
      for (size_t i = 0; i < m.size(); ++i) {
        AnnotationType &a = m[i];
        if ((regex_match(s1, a.sourceRegex)) && ((a.type.empty()) || (s2 == a.type))) return i;
      }
      return NO_ID;
    }

  public:
    void parse (string &fileName) {
      ifstream file(fileName.c_str());
      string   line, key, value, reg;
      size_t   pos;
      bool     inSynonyms = false;
      bool     inIntrons  = false;
      bool     inVicinity = false;
      bool     inOrder    = false;
      if (! file.good()) {
        cerr << "Error, configuration file '" << fileName << "' does not exists!" << endl;
        exit(EXIT_FAILURE);
      }
      while (getline(file, line)) {
        trim(line);
        if ((! line.empty()) && (line[0] != '#')) {
          if (line == "Synonyms:") {
            inSynonyms = true;
            inIntrons = inVicinity = inOrder = false;
          }
          else if (line == "Introns:") {
            inIntrons = true;
            inSynonyms = inVicinity = inOrder = false;
          }
          else if (line == "Vicinity:") {
            inVicinity = true;
            inSynonyms = inIntrons = inOrder = false;
          }
          else if (line == "Order:") {
            inOrder = true;
            inSynonyms = inIntrons = inVicinity = false;
          }
          else if (inSynonyms) {
            if (! split(line, ':', key, value)) {
              cerr << "Error, cannot parse line '" << line << "' in the 'Synonyms' section of the configuration file!" << endl;
              exit(EXIT_FAILURE);
            }
            if ((pos = key.find('*')) != string::npos) key.replace(pos, 1, ".*");
            try {
              synonyms.push_back(make_pair(regex(key), value));
            }
            catch (const regex_error &e) {
              cerr << "Error, cannot parse regular expression '" << key << "' in line '" << line << "' in the 'Synonyms' section of the configuration file!" << endl;
              exit(EXIT_FAILURE);
            }
          }
          else if (inIntrons) {
            if (! split(line, ':', key, value)) {
              cerr << "Error, cannot parse line '" << line << "' in the 'Introns' section of the configuration file!" << endl;
              exit(EXIT_FAILURE);
            }
            introns.push_back(make_tuple(key, value, NO_ID));
          }
          else if (inVicinity) {
            if (! split(line, ':', key, value)) {
              cerr << "Error, cannot parse line '" << line << "' in the 'Vicinity' section of the configuration file!" << endl;
              exit(EXIT_FAILURE);
            }
            vicinity.push_back(make_tuple(key, value, NO_ID, NO_ID));
          }
          else if (inOrder) {
            vector < string > fields;
            split(line, ',', fields);
            vector < AnnotationType > thisOrder;
            for (string field: fields) {
              string rest, strandStr;
              Strand strand = Strand::ALL;
              if (split(field, ' ', rest, strandStr)) {
                if      (strandStr == "+") strand = Strand::F;
                else if (strandStr == "-") strand = Strand::R;
                else {
                  cerr << "Error, cannot parse line '" << line << "' in the 'Order' section of the configuration file (last item item should be the strand: '+' or '-')!" << endl;
                  exit(EXIT_FAILURE);
                }
                field = rest;
              }
              if (split(field, ':', key, value)) {
                reg = key;
                if ((pos = reg.find('*')) != string::npos) reg.replace(pos, 1, ".*");
                try {
                  thisOrder.push_back(AnnotationType(key, regex(reg), value, strand));
                }
                catch (const regex_error &e) {
                  cerr << "Error, cannot parse regular expression '" << key << "' in line '" << line << "' in the 'Order' section of the configuration file!" << endl;
                  exit(EXIT_FAILURE);
                }
              }
              else {
                reg = field;
                if ((pos = reg.find('*')) != string::npos) reg.replace(pos, 1, ".*");
                try {
                  thisOrder.push_back(AnnotationType(field, regex(reg), "", strand));
                }
                catch (const regex_error &e) {
                  cerr << "Error, cannot parse regular expression '" << key << "' in line '" << line << "' in the 'Order' section of the configuration file!" << endl;
                  exit(EXIT_FAILURE);
                }
              }
              ++nElements;
            }
            order.push_back(thisOrder);
          }
          else {
            cerr << "Error, line '" << line << "' is not in the 'Synonyms', 'Introns', 'Vicinity', nor 'Order' section !" << endl;
            exit(EXIT_FAILURE);
          }
        }
      }
      for (auto &i: introns) {
        size_t o = getOrder(get<0>(i), "intron", Strand::ALL);
        if (o == NO_ID) {
          cerr << "Error, type '" << get<0>(i) << ":intron' (of '" << get<0>(i) << ":" << get<1>(i) << "') should be included in the 'Order:' section." << endl;
          exit(EXIT_FAILURE);
        }
        get<2>(i) = o;
      }
      for (auto &v: vicinity) {
        size_t i = getOrder(get<0>(v), "upstream", Strand::ALL);
        if (i == NO_ID) {
          cerr << "Error, type '" << get<0>(v) << ":upstream' (of '" << get<0>(v) << ":" << get<1>(v) << "') should be included in the 'Order:' section." << endl;
          exit(EXIT_FAILURE);
        }
        get<2>(v) = i;
        i = getOrder(get<0>(v), "downstream", Strand::ALL);
        if (i == NO_ID) {
          cerr << "Error, type '" << get<0>(v) << ":downstream' (of '" << get<0>(v) << ":" << get<1>(v) << "') should be included in the 'Order:' section." << endl;
          exit(EXIT_FAILURE);
        }
        get<3>(v) = i;
      }
      for (size_t i1 = 0; i1 < order.size(); ++i1) {
        for (size_t i2 = 0; i2 < order[i1].size(); ++i2) {
          elements.push_back(RegionType(i1, i2));
        }
      }
      cerr << "Order:" << endl;
      for (auto &a1: order) {
        for (auto &a2: a1) {
          cerr << a2 << "\t";
        }
        cerr << endl;
      }
    }

    string translate (const string &sin) {
      for (auto &line: synonyms) {
        if (regex_match(sin, line.first)) {
          return line.second;
        }
      }
      return sin;
    }

    size_t checkIntrons (const string &source, const string &type) {
      for (auto &p: introns) {
        if ((get<0>(p) == source) && (get<1>(p) == type)) return get<2>(p);
      }
      return NO_ID;
    }

    size_t checkUpstream (const string &source, const string &type) {
      for (auto &p: vicinity) {
        if ((get<0>(p) == source) && (get<1>(p) == type)) return get<2>(p);
      }
      return NO_ID;
    }

    size_t checkDownstream (const string &source, const string &type) {
      for (auto &p: vicinity) {
        if ((get<0>(p) == source) && (get<1>(p) == type)) return get<3>(p);
      }
      return NO_ID;
    }

    size_t getOrder(const string &s1, const string &s2, Strand s) {
      size_t n = 0;
      for (size_t i = 0; i < order.size(); ++i) {
        size_t subId;
        if ((subId = check(s1, s2, order[i])) != NO_ID) {
          return n + subId;
        }
        n += order[i].size();
      }
      return NO_ID;
    }

    size_t getNElements () {
      return nElements;
    }

    size_t getNElements (size_t i) {
      return order[i].size();
    }

    RegionType &getElement (size_t i) {
      return elements[i];
    }

    bool checkStrand (size_t i, Strand s1, bool s2) {
      RegionType &rt = getElement(i);
      if      (order[rt.id][rt.subId].strand == Strand::ALL) return true;
      else if (order[rt.id][rt.subId].strand == Strand::F)   return (((s1 == Strand::F) && (s2)) || ((s1 == Strand::R) && (! s2)));
      else                                                   return (((s1 == Strand::F) && (! s2)) || ((s1 == Strand::R) && (s2)));
    }

    string getName (size_t i) {
      return getName(getElement(i));
    }
    string getName (RegionType &regionType) {
      if (! regionType.isSet()) return "";
      AnnotationType &a = order[regionType.id][regionType.subId];
      string s = a.source;
      if (! a.type.empty()) {
        s += ":" + a.type;
      }
      if (a.strand == Strand::F) {
        s += " (+)";
      }
      else if (a.strand == Strand::R) {
        s += " (-)";
      }
      return s;
    }
    bool isUpstream (size_t id) {
      RegionType &regionType = getElement(id);
      return (order[regionType.id][regionType.subId].type == "upstream");
    }
    bool isDownstream (size_t id) {
      RegionType &regionType = getElement(id);
      return (order[regionType.id][regionType.subId].type == "downstream");
    }
};


void printReadStats (const string &name, unsigned int nHits, vector < size_t > &regions) {
  sort(regions.begin(), regions.end());
  readStatsFile << name << " \t" << nHits;
  size_t c = 0, cr = config->getNElements();
  for (size_t r: regions) {
    if (cr == r) {
      ++c;
    }
    else {
      if (cr != config->getNElements()) {
        readStatsFile << "\t" << config->getName(cr) << ": " << c;
      }
      cr = r;
      c  = 1;
    }
  }
  if (cr != config->getNElements()) readStatsFile << "\t" << config->getName(cr) << ": " << c;
  if (rescueFunction(regions)) readStatsFile << "\tRescued";
  readStatsFile << "\n";
}

void printReadStatsVoid (const string &name, unsigned int nHits, vector < size_t > &regions) {}

bool rescue (vector < size_t > &regions) {
  size_t n = regions.size();
  if (n == 1) return false;
  size_t t = ceil(n * rescueThreshold);
  size_t m = config->getNElements();
  vector < size_t > c(config->getNElements(), 0);
  for (size_t r: regions) {
    if (++c[r] >= t) {
      regions = { r };
      return true;
    }
  }
  return false;
}

bool rescueVoid (vector < size_t > &regions) {
  return false;
}

class GtfLineParser {
  protected:
    string source, type, chromosome, geneId;
    Position start, end;
    Strand strand;
    unordered_map < string, vector < string > > tags;

  public:
    GtfLineParser (string line) {
      vector <string> splittedLine;
      split(line, '\t', splittedLine);
      assert(splittedLine.size() == 9);
      type       = splittedLine[2];
      source     = splittedLine[1];
      chromosome = splittedLine[0];
      strand     = (splittedLine[6] == "+")? Strand::F: Strand::R;
      start      = stoul(splittedLine[3]);
      end        = stoul(splittedLine[4]);
      string remaining = splittedLine[8];
      trim(remaining);
      while (! remaining.empty()) {
        size_t splitPosSpace = remaining.find(" "), splitPosEq = remaining.find("="), splitPos, endVal, endTag;
        if      (splitPosEq    == string::npos) splitPos = splitPosSpace;
        else if (splitPosSpace == string::npos) splitPos = splitPosEq;
        else                                    splitPos = min<size_t>(splitPosSpace, splitPosEq);
        string tag = remaining.substr(0, splitPos), value;
        rtrim(tag);
        remaining = remaining.substr(splitPos+1);
        ltrim(remaining);
        if (remaining[0] == '"') {
          remaining = remaining.substr(1);
          endVal    = remaining.find("\"");
          value     = remaining.substr(0, endVal);
          remaining = remaining.substr(endVal+1);
        }
        else {
          endVal = remaining.find(";");
          if (endVal == string::npos) endVal = remaining.size();
          value  = remaining.substr(0, endVal);
          rtrim(value);
        }
        vector < string > values;
        split(value, ',', values);
        tags[tag] = values;
        endTag = remaining.find(";");
        if (endTag == string::npos) {
          remaining.clear();
        }
        else {
          remaining = remaining.substr(endTag+1);
          ltrim(remaining);
        }
      }
    }
    string          &getSource ()                { return source; }
    string          &getType ()                  { return type; }
    string          &getChromosome ()            { return chromosome; }
    Strand         getStrand ()                { return strand; }
    Position         getStart ()                 { return start; }
    Position         getEnd ()                   { return end; }
    bool             hasTag (const string &s)    { return (tags.find(s) != tags.end()); }
    vector <string> &getTag (const string &s)    { return tags[s];}
    string          &getGeneId ()                { return tags["geneId"].front(); }
    void             setSource (const string &s) { source = s; }
    void             setType (const string &t)   { type   = t; }
};

class XamRecord {
  protected:
    string name, chromosome;
    unsigned int flags, nHits;
    Position start;
    vector <pair <char, int>> cigar;
    size_t size;
    bool over;

  public:
    XamRecord (): start(UNKNOWN), nHits(1), over(false) { }
    void setChromosome (const string &c) { chromosome = c; }
    void setStart (Position s) { start = s; }
    void setName (const string &n) {
      name = n;
      if (isPaired()) {
        size_t pos = name.rfind('_');
        if ((pos != string::npos) && (pos < name.size()-1) && ((name[pos+1] == '1') || (name[pos+1] == '2'))) {
          name.resize(pos);
        }
      }
    }
    void setSize (size_t s) { size = s; }
    void setFlags (unsigned int f) { flags = f; }
    void setCigar (const vector < pair < char, int > > &g) { cigar = g; }
    void setNHits (unsigned int n) { nHits = n; }
    string &getChromosome() { return chromosome; }
    Position getStart() { return start; }
    string &getName () { return name; }
    vector <pair <char, int>> &getCigar () { return cigar; }
    size_t getSize () const { return size; }
    bool isMapped () const { return ((flags & 0x4) == 0);}
    bool isProperlyMapped () const { return ((flags & 0x2) == 0x2);}
    bool getStrand () const { return ((flags & 0x10) == 0);}
    bool isPaired () const { return ((flags & 0x1) == 0x1);}
    bool isFirst () const { return ((flags & 0x40) == 0x40);}
    bool isOver () const { return over; }
    void setOver () { over = true; }
    unsigned int getNHits () const { return nHits; }
};

class Reader {
  protected:
    ifstream file;
    XamRecord record;

  public:
    Reader (string &fileName): file(fileName.c_str()) {
      if (! file.good()) {
        cerr << "Error, file '" << fileName << "' does not exists!" << endl;
        exit(EXIT_FAILURE);
      }
    }
    virtual ~Reader () {}
    virtual XamRecord &getRecord() {
      getNextRecord();
      return record;
    }
    virtual void getNextRecord() = 0;
};

class SamReader: public Reader {
  public:
    SamReader (string &fileName): Reader(fileName) {
      lock_guard<mutex> lock(printMutex);
      cerr << "Reading SAM file " << fileName << endl;
    }
    virtual void getNextRecord () {
      string line;
      vector <string> splittedLine;
      vector <pair <char, int>> cigar;
      do {
        if (! getline(file, line)) {
          record.setOver();
          return;
        }
      }
      while ((line.empty()) || (line[0] == '@') || (line[0] == '#'));
      split(line, '\t', splittedLine);
      assert(splittedLine.size() >= 12);
      record.setFlags(stoi(splittedLine[1]));
      record.setChromosome(splittedLine[2]);
      record.setStart(stoul(splittedLine[3]));
      record.setName(splittedLine[0]);
      record.setSize(splittedLine[9].size());
      int value = 0;
      for (char c: splittedLine[5]) {
        if ((c >= '0') && (c <= '9')) {
          value *= 10;
          value += (c - '0');
        }
        else {
          cigar.push_back(make_pair(c, value));
          value = 0;
        }
      }
      record.setCigar(cigar);
      string key;
      for (unsigned int i = 11; i < splittedLine.size(); i++) {
        string &part = splittedLine[i];
        size_t pos   = part.find(':');
        key = part.substr(0, pos);
        if (key == "NH") {
          record.setNHits(stoul(part.substr(part.find(':', pos+1)+1)));
          return;
        }
      }
      record.setNHits(1);
    }
};

class BamReader: public Reader {
  protected:
    vector <string> chromosomes;
    gzFile file;
  public:
    BamReader (string &fileName): Reader(fileName) {
      lock_guard<mutex> lock(printMutex);
      cerr << "Reading BAM file " << fileName << endl;
      char buffer[1000000];
      int32_t tlen, nChrs, size;
      file = gzopen(fileName.c_str(), "rb");
      if (! file) {
        cerr << "Cannot open file '" << fileName << "'." << endl;
        return;
      }
      gzread(file, buffer, 4);
      buffer[4] = 0;
      gzread(file, reinterpret_cast<char*>(&tlen), 4);
      gzread(file, buffer, tlen);
      gzread(file, reinterpret_cast<char*>(&nChrs), 4);
      for (int i = 0; i < nChrs; i++) {
        gzread(file, reinterpret_cast<char*>(&size), 4);
        gzread(file, buffer, size);
        chromosomes.push_back(buffer);
        gzread(file, reinterpret_cast<char*>(&buffer), 4);
      }
      chromosomes.push_back("*");
    }
    virtual void getNextRecord () {
      if (gzeof(file)) {
        record.setOver();
        gzclose(file);
        return;
      }
      char     buffer[10000], v_c;
      int32_t  v_32, size, chrId, pos, lReadName, lSeq, flags, nCigar;
      uint32_t v_u32, binMqNl, flagNc;
      float    v_f;
      vector <pair <char, int>> cigar;
      gzread(file, reinterpret_cast<char*>(&size), 4);
      if (gzeof(file)) {
        record.setOver();
        gzclose(file);
        return;
      }
      gzread(file, reinterpret_cast<char*>(&chrId), 4);
      if (chrId == -1) record.setChromosome(chromosomes.back());
      else record.setChromosome(chromosomes[chrId]);
      gzread(file, reinterpret_cast<char*>(&pos), 4);
      record.setStart(++pos);
      gzread(file, reinterpret_cast<char*>(&binMqNl), 4);
      lReadName = binMqNl & 0xff;
      gzread(file, reinterpret_cast<char*>(&flagNc), 4);
      flags  = flagNc >> 16;
      nCigar = flagNc & 0xffff;
      record.setFlags(flags);
      gzread(file, reinterpret_cast<char*>(&lSeq), 4);
      record.setSize(lSeq);
      gzread(file, buffer, 4);
      gzread(file, buffer, 4);
      gzread(file, buffer, 4);
      gzread(file, buffer, lReadName);
      record.setName(buffer);
      cigar.reserve(nCigar);
      for (int i = 0; i < nCigar; i++) {
        gzread(file, reinterpret_cast<char*>(&v_u32), 4);
        uint32_t s = v_u32 >> 4;
        char op = BAM_CIGAR_LOOKUP[v_u32 & ((1 << 4) - 1)];
        cigar.push_back(make_pair(op, s));
      }
      record.setCigar(cigar);
      gzread(file, buffer, (lSeq+1)/2);
      gzread(file, buffer, lSeq);
      record.setNHits(1);
      string key(2, 0);
      char c;
      for (int32_t i = 33+lReadName+4*nCigar+(lSeq+1)/2+lSeq; i < size; ) {
        for (unsigned int j = 0; j < 2; j++) {
          gzread(file, &c, 1);
          key[j] = c;
        }
        gzread(file, &c, 1);
        i += 3;
        int8_t n = 1;
        v_32 = v_u32 = 0;
        switch(c) {
          case 'H':
            gzread(file, reinterpret_cast<char*>(&n), 1);
            c = 'C';
            i++;
            break;
          case 'B':
            int8_t s = 0, m = 1;
            gzread(file, &c, 1);
            n = 0;
            for (unsigned int j = 0; j < 4; j++) {
              gzread(file, reinterpret_cast<char*>(&s), 1);
              n += s * m;
              m *= 16;
            }
            i += 5;
            break;
        }
        for (int j = 0; j < n; j++) {
          switch(c) {
            case 'A':
              gzread(file, reinterpret_cast<char*>(&v_c), 1);
              i++;
              break;
            case 'c':
              gzread(file, reinterpret_cast<char*>(&v_32), 1);
              i++;
              break;
            case 'C':
              gzread(file, reinterpret_cast<char*>(&v_u32), 1);
              i++;
              break;
            case 's':
              gzread(file, reinterpret_cast<char*>(&v_32), 2);
              i += 2;
              break;
            case 'S':
              gzread(file, reinterpret_cast<char*>(&v_u32), 2);
              i += 2;
              break;
            case 'i':
              gzread(file, reinterpret_cast<char*>(&v_32), 4);
              i += 4;
              break;
            case 'I':
              gzread(file, reinterpret_cast<char*>(&v_u32), 4);
              i += 4;
              break;
            case 'f':
              gzread(file, reinterpret_cast<char*>(&v_f), 4);
              i += 4;
              break;
            case 'Z':
              while ((v_c = gzgetc(file)) != 0) {
                i++;
              }
              i++;
              break;
            default:
              cerr << "Problem with tag type '" << c << "'" << endl;
              return;
          }
        }
        if (key == "NH") record.setNHits(v_u32);
      }
    }
};

class Interval {
  protected:
    Position start, end;
  public:
    Interval (): start(UNKNOWN), end(UNKNOWN) {}
    Interval (Position s, Position e): start(s), end(e) { }
    Interval (GtfLineParser &line): Interval(line.getStart(), line.getEnd()) {}
    Position getStart () const     { return start; }
    Position getEnd   () const     { return end; }
    void     setStart (Position p) { start  = p; }
    void     setEnd   (Position p) { end    = p; }
    void     unSet    ()           { start = end = UNKNOWN; }
    bool isSet () const {
      return ((start != UNKNOWN) && (end != UNKNOWN));
    }
    Position getSize () const {
      return (end - start + 1);
    }
    size_t overlaps (const Interval &i) const {
      Position s = max<Position>(start, i.start), e = min<Position>(end, i.end);
      if (s >= e) return 0;
      return e - s;
    }
    bool includes (const Interval &i) const {
      return ((i.start >= start) && (i.end <= end));
    }
    bool isIncluded (const Interval &i) const {
      return i.includes(*this);
    }
    bool isBefore (const Interval &i) const {
      return (end < i.start);
    }
    bool isAfter (const Interval &i) const {
      return (i.isBefore(*this));
    }
    void merge (const Interval &i) {
      start = min<Position>(start, i.start);
      end   = max<Position>(end,   i.end);
    }
    void intersect (Interval &i) {
      if (! isSet()) return;
      start = max<Position>(start, i.start);
      end   = min<Position>(end,   i.end);
      if (start > end) {
        start = end = UNKNOWN;
      }
    }
    size_t getDistance (Position p) {
      if (p < start) return (start - p);
      if (p > end  ) return (p - end);
      return 0;
    }
    friend bool operator<(const Interval &i1, const Interval &i2) {
      return (i1.start < i2.start);
    }
};

ostream &operator<< (ostream &os, const Interval &i) {
  os << i.getStart() << "-" << i.getEnd();
  return os;
}


class PlacedElement {
  protected:
    Strand strand;
    unsigned int chromosomeId;
  public:
    PlacedElement (Strand st = Strand::ALL, unsigned int c = NO_CHR): strand(st), chromosomeId(c) { }
    unsigned int getChromosomeId () const { return chromosomeId; }
    Strand       getStrand ()       const { return strand; }
};


class TypedInterval: public Interval, public PlacedElement {
  protected:
    size_t type;
    string id;
  public:
    TypedInterval (Position s = UNKNOWN, Position e = UNKNOWN, size_t t = NO_ID, Strand st = Strand::ALL, unsigned int c = NO_CHR, const string &i = EMPTY_STRING): Interval(s, e), PlacedElement(st, c), type(t), id(i) { }
    TypedInterval (Interval &i, size_t t, Strand st, unsigned int c, const string &n = EMPTY_STRING): Interval(i), PlacedElement(st, c), type(t), id(n) { }
    size_t getType () const { return type; }
    string &getId () { return id; }
    friend bool operator<(const TypedInterval &i1, const TypedInterval &i2) {
      return ((i1.chromosomeId < i2.chromosomeId) || ((i1.chromosomeId == i2.chromosomeId) && (i1.start < i2.start)));
    }
};

ostream &operator<< (ostream &os, const TypedInterval &ti) {
  os << ti.getChromosomeId() << ":" << ti.getStart() << "-" << ti.getEnd() << " (" << config->getName(ti.getType()) << ")";
  return os;
}


class Transcript: public Interval {
  protected:
    vector <Interval> exons;
    vector <Interval> introns;
  public:
    Transcript (Position s = UNKNOWN, Position e = UNKNOWN): Interval(s, e) {}
    Transcript (Interval e): Interval(e) {}
    Transcript (vector < Interval > &i) {
      for (Interval &e: i) addExon(e);
    }
    Transcript (GtfLineParser &line): Transcript(line.getStart(), line.getEnd()) {}
    void addExon (const Interval &e) {
      this->Interval::merge(e);
      exons.push_back(e);
      vector < Interval > newExons;
      sort(exons.begin(), exons.end());
      Interval currentExon;
      for (Interval &exon: exons) {
        if (! currentExon.isSet()) {
          currentExon = exon;
        }
        else if (currentExon.isBefore(exon)) {
          newExons.push_back(currentExon);
          currentExon = exon;
        }
        else {
          currentExon.merge(exon);
        }
      }
      newExons.push_back(currentExon);
      exons = newExons;
    }
    void addExon (Position s, Position e) {
      addExon(Interval(s, e));
    }
    void checkStructure () {
      sort(exons.begin(), exons.end());
      if (exons.empty()) exons.push_back(*this);
      Position start = UNKNOWN;
      for (Interval &exon: exons) {
        if (start != UNKNOWN) introns.push_back(Interval(start, exon.getStart()-1));
        start = exon.getEnd()+1;
        //cout << "Exon " << exon << endl;
      }
      //cout << "\t" << (*this) << endl;
      //cout << "\tIntron size: " << introns.size() << endl;
    }
    size_t overlaps (Transcript &t) {
      if (Interval::overlaps(t) == 0) return 0;
      //cout << "passed" << endl;
      size_t o = 0;
      for (Interval &i1: exons) {
        //cout << "i1 (" << name << "): " << i1.getStart() << "-" << i1.getEnd() << endl;
        for (Interval &i2: t.exons) {
          //cout << "\ti2 (" << name << "): " << i2.getStart() << "-" << i2.getEnd() << endl;
          //cout << "\t\t" << i1.overlaps(i2) << endl;
          o += i1.overlaps(i2);
        }
      }
      return o;
    }
    bool includes (Transcript &t) {
      //cout << "include transcript: " << start << "-" << end << " vs " << t.start << "-" << t.end << endl;
      if (! Interval::includes(t)) return false;
      //cout << "\tpassed" << endl;
      for (Interval &e2: t.exons) {
        if (! any_of(exons.begin(), exons.end(), [&e2](Interval &e1){return e1.includes(e2);})) return false;
      }
      return true;
    }
    bool isIncluded (Transcript &t) {
      return t.isIncluded(*this);
    }
    void intersect (Interval &i) {
      introns.clear();
      vector <Interval> newExons;
      for (Interval &exon: exons) {
        exon.intersect(i);
        if (exon.isSet()) {
          newExons.push_back(exon);
        }
      }
      exons = newExons;
      if (exons.empty()) {
        start = end = UNKNOWN;
        return;
      }
      start = exons.front().getStart();
      end   = exons.back().getEnd();
    }
    void merge (Interval &e) {
      vector < Interval > newExons;
      for (Interval &thisExon: exons) {
        if (thisExon.overlaps(e)) {
          e.merge(thisExon);
        }
        else {
          newExons.push_back(thisExon);
        }
      }
      newExons.push_back(e);
      exons = newExons;
    }
    void merge (Transcript &t) {
      for (Interval &exon: t.getExons()) {
        merge(exon);
      }
    }
    vector <Interval> &getExons () {
      return exons;
    }
    vector <Interval> &getIntrons () {
      //cout << "\t\t# introns: " << introns.size() << endl;
      return introns;
    }
    friend ostream &operator<< (ostream &os, const Transcript &t);
};

ostream &operator<< (ostream &os, const Transcript &t) {
  for (const Interval &exon: t.exons) os << exon << " ";
  return os;
}

inline bool strandF (bool strand, bool isFirst, bool isPaired) {
  if (isPaired) {
    cerr << "Error! Strandedness 'F' should be used for single-end reads.\nExiting." << endl;
    exit(EXIT_FAILURE);
  }
  return strand;
}
inline bool strandR (bool strand, bool isFirst, bool isPaired) {
  if (isPaired) {
    cerr << "Error! Strandedness 'R' should be used for single-end reads.\nExiting." << endl;
    exit(EXIT_FAILURE);
  }
  return ! strand;
}
inline bool strandFF (bool strand, bool isFirst, bool isPaired) {
  if (! isPaired) {
    cerr << "Error! Strandedness 'FF' should be used for paired-end reads.\nExiting." << endl;
    exit(EXIT_FAILURE);
  }
  return strand;
}
inline bool strandFR (bool strand, bool isFirst, bool isPaired) {
  if (! isPaired) {
    cerr << "Error! Strandedness 'FR' should be used for paired-end reads.\nExiting." << endl;
    exit(EXIT_FAILURE);
  }
  return (strand == isFirst);
}
inline bool strandRF (bool strand, bool isFirst, bool isPaired) {
  if (! isPaired) {
    cerr << "Error! Strandedness 'RF' should be used for paired-end reads.\nExiting." << endl;
    exit(EXIT_FAILURE);
  }
  return (strand != isFirst);
}
inline bool strandU (bool strand, bool isFirst, bool isPaired) {
  return true;
}

class Read: public Transcript {
  protected:
    string name, chromosome;
    bool strand, paired, first;
    unsigned int nHits;
    size_t size;
  public:
    Read () {}
    Read (XamRecord &record, StrandednessFunction strandednessFunction): Transcript(record.getStart()), name(record.getName()), chromosome(record.getChromosome()), paired(record.isPaired()), nHits(record.getNHits()), size(0) {
      addPart(record);
      strand = strandednessFunction(record.getStrand(), record.isFirst(), record.isPaired());
    }
    void addPart(XamRecord &record) {
      Position s = record.getStart(), e = record.getStart();
      start      = min<Position>(start, record.getStart());
      first      = record.isFirst();
      for (auto &part: record.getCigar()) {
        char c = part.first;
        int  v = part.second;
        switch (c) {
          case 'M':
          case 'D':
          case '=':
          case 'X':
            e += v;
            break;
          case 'N':
            if (s != e) addExon(s, e-1);
            introns.push_back(Interval(e, e+v-1));
            //cout << "Adding intron in read " << introns.back().getStart() << "-" << introns.back().getEnd() << endl;
            e += v;
            s = e;
            break;
          case 'I':
          case 'S':
          case 'H':
          case 'P':
            break;
          default:
            cerr << "Problem in the cigar: do not understand char " << c << endl;
        }
      }
      if (s != e) addExon(s, e-1);
      --e;
      if ((end == UNKNOWN) || (e > end)) end = e;
      nHits = min<unsigned int>(nHits, record.getNHits());
      size += record.getSize();
      //cout << "Read: " << name << ", " << chromosome; for (Interval &e: exons) cout << " exon " << e.getStart() << "-" << e.getEnd(); cout << endl;
    }
    string &getName () { return name; }
    bool getStrand() { return strand; }
    string &getChromosome() { return chromosome; }
    unsigned int getNHits () { return nHits; }
    size_t getSize () { return size; }
    bool isPaired () { return paired; }
    bool hasFirst () { return first; }
    friend ostream &operator<< (ostream &os, const Read &r);
};

ostream &operator<< (ostream &os, const Read &r) {
  os << r.chromosome << ":" << static_cast<const Transcript &>(r);
  return os;
}


class Gene: public Interval, public PlacedElement {
  protected:
    string id, source, type;
    Transcript mergedTranscript, cds, utr5, utr3;
    Interval upstream, downstream;
  public:
    Gene (const string &i, const string &s, const string &t, Position st, Position en, Strand str, unsigned int ci): Interval(st, en), PlacedElement(str, ci), id(i), source(s), type(t), mergedTranscript(st, en) { }
    Gene (GtfLineParser &line, unsigned int c): Gene(line.hasTag("gene_id")? line.getTag("gene_id").front(): (line.hasTag("ID")? line.getTag("ID").front(): rtrimTo(line.getTag("Parent").front(), '.')), line.getSource(), line.getType(), line.getStart(), line.getEnd(), line.getStrand(), c) { }
    string &getId ()     { return id; }
    string &getSource () { return source; }
    string &getType ()   { return type; }
    void addExon (Interval &e) {
      Interval::merge(e);
      mergedTranscript.addExon(e);
    }
    void addCds (Interval &c) {
      addExon(c);
      if (cds.isSet()) cds.Interval::merge(c);
      else             cds = c;
    }
    void setCds () {
      if (! cds.isSet()) return;
      Interval intervalCds = cds;
      cds = mergedTranscript;
      cds.intersect(intervalCds);
    }
    void setUtr () {
      if (! cds.isSet()) return;
      Interval i5(start, cds.getStart()-1), i3(cds.getEnd()+1, end);
      utr5 = utr3 = mergedTranscript;
      utr5.intersect(i5);
      utr3.intersect(i3);
      if (strand == Strand::R) swap(utr5, utr3);
    }
    void setUpDownStream () {
      if (strand == Strand::F) {
        upstream   = Interval(((getStart() <= upstreamSize)? 1: getStart()-upstreamSize), getStart()-1);
        downstream = Interval(getEnd()+1, getEnd()+downstreamSize);
      }
      else {
        downstream = Interval(((getStart() <= downstreamSize)? 1: getStart()-downstreamSize), getStart()-1);
        upstream   = Interval(getEnd()+1, getEnd()+upstreamSize);
      }
    }
    void checkStructure () {
      mergedTranscript.checkStructure();
      start = mergedTranscript.getStart();
      end   = mergedTranscript.getEnd();
      setCds();
      setUtr();
      setUpDownStream();
    }
    size_t overlaps (Read &read) {
      return mergedTranscript.overlaps(read);
    }
    bool includes (Read &read) {
      //cout << "gene include" << endl;
      return (mergedTranscript.includes(read));
    }
    Transcript &getMergedTranscript () {
      return mergedTranscript;
    }
    Transcript &getCds () {
      return cds;
    }
    Transcript &getUtr5 () {
      return utr5;
    }
    Transcript &getUtr3 () {
      return utr3;
    }
    Interval &getUpstream () {
      return upstream;
    }
    Interval &getDownstream () {
      return downstream;
    }
    friend bool operator<(const Gene &g1, const Gene &g2) {
      return ((g1.chromosomeId < g2.chromosomeId) || ((g1.chromosomeId == g2.chromosomeId) && (g1.start < g2.start)));
    }
};

inline size_t intervalInclusion (Interval &i, Interval &r) {
  return (i.includes(r))? 1: 0;
}
inline size_t intervalOverlapPc (Interval &i, Interval &r) {
  size_t o = i.overlaps(r);
  return (r.getSize() * overlap <= o)? o: 0;
}
inline size_t intervalOverlap (Interval &i, Interval &r) {
  size_t o = i.overlaps(r);
  return (o >= overlap)? o: 0;
}

struct IntervalListPosition {
  size_t chromosomeId, intervalId;
  IntervalListPosition (): chromosomeId(0), intervalId(0) {}
  void reset () {
    chromosomeId = intervalId = 0;
  }
};

class EvaluationStructure {
  protected:
    vector < vector < pair < size_t, size_t > > > evaluation;
    vector < vector < size_t > > ids;
    unsigned int nHits;
  public:
    EvaluationStructure (Read &r): evaluation(vector<vector<pair<size_t, size_t>>>(config->getNElements(), vector<pair<size_t, size_t>>(r.getExons().size()))), ids(config->getNElements()), nHits(0) {
      for (size_t i = 0; i < config->getNElements(); ++i) {
        for (size_t j = 0; j < r.getExons().size(); ++j) {
          evaluation[i][j] = {0, 0};
        }
      }
    }
    void set (size_t id, size_t type, size_t exonId, size_t overlap, size_t distance) {
      ids[type].push_back(id);
      evaluation[type][exonId].first  = overlap;
      evaluation[type][exonId].second = distance;
      ++nHits;
    }
    void getFirst (vector <size_t> &regions) {
      if (nHits == 0) return;
      //cout << "starting" << endl;
      size_t goodId = NO_ID;
      RegionType regionType;
      vector < size_t > selected;
      size_t maxOverlap = 0;
      for (size_t i = 0; i < evaluation.size(); ++i) {
        regionType = config->getElement(i);
        //cout << "\t1: " << regionType.id << ", " << regionType.subId << endl;
        if ((goodId != NO_ID) && (regionType.id != goodId)) {
          //cout << "\t2" << endl;
          break;
        }
        bool presentInAll = true;
        size_t overlap = 0;
        for (size_t exonId = 0; (exonId < evaluation[regionType.id].size()) && (presentInAll); ++exonId) {
          if (evaluation[i][exonId].first > 0) {
            //cout << "\t\t3: " << evaluation[i][exonId].first << endl;
            overlap += evaluation[i][exonId].first;
            goodId = regionType.id;
          }
          else {
            //cout << "\t\t3: not here" << endl;
            presentInAll = false;
          }
        }
        if (overlap > 0) {
          if (overlap > maxOverlap) {
            //cout << "\t4" << endl;
            selected   = { i };
            maxOverlap = overlap;
          }
          else if (overlap == maxOverlap) {
            //cout << "\t5" << endl;
            selected.push_back(i);
          }
        }
      }
      if (selected.empty()) {
        //cout << "\t6" << endl;
        return;
      }
      if (selected.size() == 1) {
        //cout << "\t7" << endl;
        regions = selected;
        return;
      }
      size_t minDistance = numeric_limits<size_t>::max();
      for (size_t i: selected) {
        for (auto &p: evaluation[i]) {
          if (p.second < minDistance) {
            minDistance = p.second;
            regions = { i };
          }
          else if (p.second == minDistance) {
            regions.push_back(i);
          }
        }
      }
    }
    void getIds (vector <size_t> &regions, vector <size_t> &intervals) {
      for (size_t r: regions) {
        intervals.insert(intervals.end(), ids[r].begin(), ids[r].end());
      }
    }
};

class IntervalList {
  protected:
    vector <string>                         chromosomes;
    vector <string>                         unknownChromosomes;
    vector <TypedInterval>                  intervals;
    vector <unsigned int>                   chrStarts;
		unordered_map <string, vector <size_t>> bins;


  public:
    IntervalList(string &fileName) {
      ifstream file (fileName.c_str());
      unordered_map < string, unsigned int> geneHash;
      unordered_set < string > unused;
      vector < Gene > genes;
      string line, chromosome;
      unsigned long cpt;
      unsigned int chromosomeId = numeric_limits<unsigned int>::max();
      cerr << "Reading GTF file" << endl;
      for (cpt = 0; getline(file, line); cpt++) {
        if ((! line.empty()) && (line[0] != '#')) {
          //cout << line << endl;
          GtfLineParser parsedLine(line);
          parsedLine.setSource(config->translate(parsedLine.getSource()));
          parsedLine.setType(config->translate(parsedLine.getType()));
          if (parsedLine.getChromosome() != chromosome) {
            geneHash.clear();
            unused.clear();
            chromosome = parsedLine.getChromosome();
            bool seen = false;
            for (unsigned int i = 0; (i < chromosomes.size()) && (! seen); i++) {
              if (chromosomes[i] == chromosome) {
                chromosomeId = i;
                seen         = true;
              }
            }
            if (! seen) {
              chromosomeId = chromosomes.size();
              chromosomes.push_back(chromosome);
            }
          }
          if (parsedLine.getType() == "gene") {
            string geneId;
            if      (parsedLine.hasTag("ID"))      geneId = parsedLine.getTag("ID").front();
            else if (parsedLine.hasTag("gene_id")) geneId = parsedLine.getTag("gene_id").front();
            else {
              cerr << "Warning, cannot deduce gene id at line " << cpt << ": '" << line << "'." << endl;
            }
            Gene gene(parsedLine, chromosomeId);
            geneHash[geneId] = genes.size();
            genes.push_back(gene);
          }
          else if (parsedLine.getType() == "transcript") {
            string id, geneId;
            if      (parsedLine.hasTag("ID"))            id = parsedLine.getTag("ID").front();
            else if (parsedLine.hasTag("transcript_id")) id = parsedLine.getTag("transcript_id").front();
            else {
              cerr << "Warning, cannot deduce transcript id at line " << cpt << ": '" << line << "'." << endl;
            }
            if      (parsedLine.hasTag("Parent"))  geneId = parsedLine.getTag("Parent").front();
            else if (parsedLine.hasTag("gene_id")) geneId = parsedLine.getTag("gene_id").front();
            else {
              cerr << "Warning, cannot deduce transcript parent id at line " << cpt << ": '" << line << "'." << endl;
            }
            if ((unused.find(geneId) == unused.end()) && (geneHash.find(geneId) != geneHash.end())) {
              geneHash[id] = geneHash[geneId];
            }
          }
          else if (parsedLine.getType() == "exon") {
            string geneId;
            if      (parsedLine.hasTag("Parent"))  geneId = parsedLine.getTag("Parent").front();
            else if (parsedLine.hasTag("gene_id")) geneId = parsedLine.getTag("gene_id").front();
            else {
              cerr << "Warning, cannot deduce exon id at line " << cpt << ": '" << line << "'." << endl;
            }
            if (unused.find(geneId) == unused.end()) {
              auto pos = geneHash.find(geneId);
              Interval e(parsedLine);
              if (pos == geneHash.end()) {
                Gene gene(parsedLine, chromosomeId);
                gene.addExon(e);
                geneHash[geneId] = genes.size();
                genes.push_back(gene);
              }
              else {
                genes[pos->second].addExon(e);
              }
            }
          }
          else if (parsedLine.getType() == "CDS") {
            string geneId;
            if      (parsedLine.hasTag("gene_id")) geneId = parsedLine.getTag("gene_id").front();
            else if (parsedLine.hasTag("Parent"))  geneId = parsedLine.getTag("Parent").front();
            else {
              cerr << "Warning, cannot deduce CDS parent id at line " << cpt << ": '" << line << "'." << endl;
            }
            auto pos = geneHash.find(geneId);
            Interval e(parsedLine);
            if (pos == geneHash.end()) {
              Gene gene(parsedLine, chromosomeId);
              gene.addCds(e);
              geneHash[geneId] = genes.size();
              genes.push_back(gene);
            }
            else {
              genes[pos->second].addCds(e);
            }
          }
          else if (parsedLine.getType() == "5'UTR") {
            // skip this
          }
          else if (parsedLine.getType() == "3'UTR") {
            // skip this
          }
          else if (config->getOrder(parsedLine.getSource(), parsedLine.getType(), parsedLine.getStrand()) != NO_ID) {
            string id;
            if      (parsedLine.hasTag("ID"))      id = parsedLine.getTag("ID").front();
            else if (parsedLine.hasTag("gene_id")) id = parsedLine.getTag("gene_id").front();
            else {
              if (parsedLine.hasTag("Parent")) id = parsedLine.getTag("Parent").front() + "_" + parsedLine.getType();
              else {
                cerr << "Warning, cannot deduce id at line " << cpt << ": '" << line << "'." << endl;
              }
            }
            geneHash[id] = genes.size();
            genes.push_back(Gene(parsedLine, chromosomeId));
          }
          else {
            if (parsedLine.hasTag("gene_id"))       unused.insert(parsedLine.getTag("gene_id").front());
            if (parsedLine.hasTag("transcript_id")) unused.insert(parsedLine.getTag("transcript_id").front());
            if (parsedLine.hasTag("ID"))            unused.insert(parsedLine.getTag("ID").front());

          }
        }
        if (progress && (cpt % 100000 == 0)) cerr << "\t" << cpt << " lines read.\r" << flush;
      }
      cerr << "\t" << cpt << " lines read, done.  " << genes.size() << " genes found." << endl;
      for (Gene &gene: genes) {
        string      &source       = gene.getSource();
        string      &type         = gene.getType();
        Strand       strand       = gene.getStrand();
        unsigned int chromosomeId = gene.getChromosomeId();
        size_t       regionType;
        //cerr << "\t\tGene: " << gene.getSource() << ":" << gene.getType() << endl;
        gene.checkStructure();
        if ((regionType = config->getOrder(source, "CDS", strand)) != NO_ID) {
          for (Interval &exon: gene.getCds().getExons()) {
            intervals.push_back(TypedInterval(exon, regionType, strand, chromosomeId, gene.getId() + "-CDS"));
          }
        }
        if ((regionType = config->getOrder(source, "5'UTR", strand)) != NO_ID) {
          for (Interval &exon: gene.getUtr5().getExons()) {
            intervals.push_back(TypedInterval(exon, regionType, strand, chromosomeId, gene.getId() + "-5UTR"));
          }
        }
        if ((regionType = config->getOrder(source, "3'UTR", strand)) != NO_ID) {
          for (Interval &exon: gene.getUtr3().getExons()) {
            intervals.push_back(TypedInterval(exon, regionType, strand, chromosomeId, gene.getId() + "-3UTR"));
          }
        }
        if ((regionType = config->checkIntrons(source, type)) != NO_ID) {
          for (Interval &intron: gene.getMergedTranscript().getIntrons()) {
            intervals.push_back(TypedInterval(intron, regionType, strand, chromosomeId, gene.getId() + "-intron"));
          }
        }
        if ((regionType = config->checkUpstream(source, type)) != NO_ID) {
          intervals.push_back(TypedInterval(gene.getUpstream(), regionType, strand, chromosomeId, gene.getId() + "-upstream"));
        }
        if ((regionType = config->checkDownstream(source, type)) != NO_ID) {
          intervals.push_back(TypedInterval(gene.getDownstream(), regionType, strand, chromosomeId, gene.getId() + "-downstream"));
        }
        if ((regionType = config->getOrder(source, type, strand)) != NO_ID) {
          for (Interval &exon: gene.getMergedTranscript().getExons()) {
            intervals.push_back(TypedInterval(exon, regionType, strand, chromosomeId, gene.getId()));
          }
        }
      }
      sort(intervals.begin(), intervals.end());
      chrStarts = vector<unsigned int>(chromosomes.size());
      chromosomeId = numeric_limits<unsigned int>::max();
      for (unsigned int i = 0; i < intervals.size(); i++) {
        //cerr << "\t" << intervals[i] << endl;
        if (intervals[i].getChromosomeId() != chromosomeId) {
          chrStarts[chromosomeId = intervals[i].getChromosomeId()] = i;
        }
      }
			if (! sorted) {
				//cerr << "Not all sorted" << endl;
				for (unsigned int i = 0; i < intervals.size(); i++) {
					vector<size_t> &chrBins    = bins[chromosomes[intervals[i].getChromosomeId()]];
					unsigned        bin        = intervals[i].getEnd() / binSize;
					if (chrBins.size() <= bin) {
						//cerr << "Adding " << chrBins.size() << "-" << bin << ": " << i << endl;
						chrBins.insert(chrBins.end(), bin-chrBins.size()+1, i);
					}
				}
			}
      cerr << "\t" << intervals.size() << " intervals found." << endl;
    }
    void scan(Read &read, vector <size_t> &regions, vector <size_t> &selectedIntervals, IntervalListPosition &position, Strandedness strandedness) {
      if (sorted) {
        if (chromosomes[position.chromosomeId] != read.getChromosome()) {
          if (find(unknownChromosomes.begin(), unknownChromosomes.end(), read.getChromosome()) != unknownChromosomes.end()) return;
          for (position.chromosomeId = 0; (position.chromosomeId < chromosomes.size()) && (chromosomes[position.chromosomeId] != read.getChromosome()); position.chromosomeId++) ;
          if (position.chromosomeId == chromosomes.size()) {
            unknownChromosomes.push_back(read.getChromosome());
            position.reset();
            return;
          }
          position.intervalId = chrStarts[position.chromosomeId];
        }
      }
      else {
				unsigned int bin = read.getStart() / binSize;
				auto         p   = bins.find(read.getChromosome());
				if (p == bins.end()) return;
				bin = min<unsigned int>(bin, p->second.size()-1);
				position.intervalId   = p->second[bin];
				position.chromosomeId = intervals[position.intervalId].getChromosomeId();
      }
      while ((position.intervalId < intervals.size()) && (intervals[position.intervalId].getChromosomeId() == position.chromosomeId) && (intervals[position.intervalId].isBefore(read))) {
        position.intervalId++;
      }
      size_t id = position.intervalId;
      EvaluationStructure evaluation(read);
      while ((id < intervals.size()) && (intervals[id].getChromosomeId() == position.chromosomeId) && (! intervals[id].isAfter(read))) {
        TypedInterval &interval = intervals[id];
        if (config->checkStrand(interval.getType(), interval.getStrand(), read.getStrand())) {
          for (size_t exonId = 0; exonId < read.getExons().size(); ++exonId) {
            Interval &exon = read.getExons()[exonId];
            size_t o;
            if ((o = intervalOverlapFunction(interval, exon)) > 0) {
              size_t d = 0;
              if (config->isUpstream(interval.getType())) {
                d = exon.getDistance(interval.getEnd());
              }
              else if (config->isDownstream(interval.getType())) {
                d = exon.getDistance(interval.getStart());
              }
              evaluation.set(id, interval.getType(), exonId, o, d);
            }
          }
        }
        id++;
      }
      evaluation.getFirst(regions);
      //cout << "Evaluating " << regions.size() << " intervals." << endl;
      if (intervalStats) {
        evaluation.getIds(regions, selectedIntervals);
      }
    }
    TypedInterval &getInterval (size_t i) {
      return intervals[i];
    }
};

class Counter {
  protected:
    IntervalList &intervalList;
    unordered_map<string, pair <unsigned int, vector <size_t>>> readCounts;
    unordered_map<string, size_t> rawCounts;
    unordered_map<vector<size_t>, double> regionCounts;
    unordered_map<vector<size_t>, unsigned int> intervalCounts;
    unordered_map<string, vector <size_t>> readsIntervals;
    unordered_map<string, size_t> chosenId, numberSeen;
    unordered_set<string> seen;
    unsigned int nHits, nReads, nUnique, nAmbiguous, nMultiple, nUnassigned, nRescued;
    string fileName;
    void addCount(const string &read, vector <size_t> &regions, vector <size_t> &intervals, unsigned int nHits) {
      if      (regions.empty())    nUnassigned++;
      else if (regions.size() > 1) nAmbiguous++;
      else if (nHits == 1)         nUnique++;
      if ((nHits > 1) && (strategy == Strategy::DEFAULT)) {
        nMultiple++;
        auto pos = readCounts.find(read);
        if (pos == readCounts.end()) {
          readCounts[read] = make_pair(nHits-1, regions);
          rawCounts[read]  = nHits;
          ++nReads;
          if (intervalStats) readsIntervals[read] = intervals;
        }
        else {
          unordered_map<string, vector <size_t>>::iterator pri;
          if (intervalStats) pri = readsIntervals.find(read);
          pos->second.first--;
          pos->second.second.insert(pos->second.second.end(), regions.begin(), regions.end());
          if (intervalStats) readsIntervals[read].insert(readsIntervals[read].end(), intervals.begin(), intervals.end());
          if (pos->second.first == 0) {
            if (! pos->second.second.empty()) {
              printReadStatsFunction(read, nHits, pos->second.second);
              setUnique(pos->second.second);
              ++regionCounts[pos->second.second];
              if (pos->second.second.size() == 1) nRescued++;
              if (intervalStats) {
                if (! pri->second.empty()) {
                  sort(pri->second.begin(), pri->second.end());
                  ++intervalCounts[pri->second];
                }
                readsIntervals.erase(pri);
              }
            }
            readCounts.erase(pos);
            rawCounts.erase(read);
          }
        }
      }
      else {
        if (! regions.empty()) {
          bool output = false;
          if (strategy == Strategy::RANDOM) {
            if (seen.find(read) == seen.end()) {
              auto         p = chosenId.find(read);
              unsigned int i;
              if (p == chosenId.end()) {
                i                = rand() % nHits;
                chosenId[read]   = i;
                numberSeen[read] = 0;
              }
              else {
                i = p->second;
                ++numberSeen[read];
              }
              if (numberSeen[read] == i) {
                output = true;
                chosenId.erase(read);
                numberSeen.erase(read);
                seen.insert(read);
              }
            }
          }
          if ((strategy != Strategy::RANDOM) || (output)) {
            printReadStatsFunction(read, nHits, regions);
            setUnique(regions);
            regionCounts[regions] += (strategy == Strategy::RATIO)? 1.0/nHits: 1;
            if (! intervals.empty()) {
              sort(intervals.begin(), intervals.end());
              ++intervalCounts[intervals];
            }
          }
        }
        ++nReads;
      }
    }
  public:
    Counter (IntervalList &gl): intervalList(gl) { }
    void clear () {
      readCounts.clear();
      regionCounts.clear();
      rawCounts.clear();
      nHits = nReads = nUnique = nAmbiguous = nMultiple = nUnassigned = nRescued = 0;
    }
    void read (string &f, Strandedness strandedness, StrandednessFunction strandednessFunction, ReadsFormat format) {
      Reader *reader;
      fileName = f;
      if      (format == ReadsFormat::BAM) reader = new BamReader(fileName);
      else if (format == ReadsFormat::SAM) reader = new SamReader(fileName);
      else {
        if (fileName.size() < 4) {
          cerr << "Cannot deduce type from file name '" << fileName << "'.  Should be a .sam or .bam file.  Please specify it using the '-f' option." << endl;
          exit(EXIT_FAILURE);
        }
        string suffix = fileName.substr(fileName.size()-4);
        lower(suffix);
        if      (suffix == ".bam") reader = new BamReader(fileName);
        else if (suffix == ".sam") reader = new SamReader(fileName);
        else {
          cerr << "Cannot deduce type from file name '" << fileName << "'.  Should be a .sam or .bam file.  Please specify it using the '-f' option." << endl;
          exit(EXIT_FAILURE);
        }
      }
      unsigned int cpt;
      IntervalListPosition position;
      Position previousPos = 0;
      unordered_map < string, array < queue < Read >, 2 > > pendingReads;
      XamRecord &record = reader->getRecord();
      regionCounts.clear();
      for (cpt = 0; !record.isOver(); cpt++, reader->getNextRecord()) {
        if (record.isMapped()) {
          Read read;
          bool pending = true;
          if (record.isPaired()) {
            string readName = record.getName();
            size_t bucket   = record.isFirst()? 0: 1;
            auto   pos      = pendingReads.find(readName);
            if (pos != pendingReads.end()) {
              queue < Read > &reads = pos->second[1-bucket];
              if (! reads.empty()) {
                pending = false;
                read = reads.front();
                read.addPart(record);
                reads.pop();
                if ((reads.empty()) && (pos->second[bucket].empty())) pendingReads.erase(pos);
              }
            }
            if (pending) {
              read = Read(record, strandednessFunction);
              pendingReads[readName][bucket].push(read);
            }
          }
          else {
            pending = false;
            read = Read(record, strandednessFunction);
          }
          if (! pending) {
            if ((strategy != Strategy::UNIQUE) || (read.getNHits() == 1)) {
              nHits++;
              vector < size_t > regions, intervals;
              intervalList.scan(read, regions, intervals, position, strandedness);
              addCount(read.getName(), regions, intervals, read.getNHits());
            }
          }
          if (previousPos > record.getStart()) pendingReads.clear();
          previousPos = record.getStart();
        }
        if (progress && (nThreads == 1) && (cpt % 1000000 == 0)) cerr << "\t" << cpt << " lines read.\r" << flush;
      }
      cerr << "\t" << cpt << " lines read, done." << endl;
      for (auto &e: readCounts) {
        //cout << "Problem with " << e.first << ": " << e.second.first << " hits are missing." << endl;
        if (! e.second.second.empty()) {
          if ((strategy != Strategy::UNIQUE) || (rawCounts[e.first] == 1)) {
            printReadStatsFunction(e.first, rawCounts[e.first], e.second.second);
            setUnique(e.second.second);
            regionCounts[e.second.second] += (strategy == Strategy::RATIO)? 1.0/rawCounts[e.first]: 1;
            if ((rawCounts[e.first] > 1) && (e.second.second.size() == 1)) nRescued++;
          }
        }
      }
      if (intervalStats) {
        for (auto &e: readsIntervals) {
          if (! e.second.empty()) {
            sort(e.second.begin(), e.second.end());
            ++intervalCounts[e.second];
          }
        }
      }
      pendingReads.clear();
      delete reader;
    }
    unordered_map<vector<size_t>, double> &getCounts () {
      return regionCounts;
    }
    void dump () {
      cerr << "Results for " << fileName << ":" << endl;
      if (nHits == 0) {
        cerr << "\tNo hit." << endl;
      }
      else {
        cerr << "\t# reads:                       " << nReads << "\n";
        printStats(nUnique,       "# uniquely mapped reads:       ", nReads);
        printStats(nRescued,       "# multi-mapping rescued reads: ", nReads);
        cerr << "\t# hits:                        " << nHits << "\n";
        printStats(nAmbiguous,  "# ambiguous hits:              ", nHits);
        printStats(nUnassigned, "# unassigned hits:             ", nHits);
      }
      if (intervalStats) {
        vector < pair < string, unsigned int > > lines;
        for (auto &p: intervalCounts) {
          vector < string > names;
          string name;
          for (size_t i: p.first) {
            TypedInterval &interval = intervalList.getInterval(i);
            names.push_back(interval.getId() + " (" + config->getName(interval.getType()) + ")");
          }
          sort(names.begin(), names.end());
          join(names, name, " -- ");
          lines.push_back(make_pair(name, p.second));
        }
        sort(lines.begin(), lines.end());
        string currentName;
        unsigned int count;
        for (auto &line: lines) {
          if (line.first == currentName) {
            count += line.second;
          }
          else {
            if (! currentName.empty()) {
              intervalStatsFile << currentName << "\t" << count << "\n";
            }
            currentName = line.first;
            count       = line.second;
          }
        }
        if (! currentName.empty()) {
          intervalStatsFile << currentName << "\t" << count << "\n";
        }
      }
    }
};
class TableCount {
  protected:
    IntervalList &intervalList;
    unsigned int nColumns;
    unordered_map<vector<size_t>, vector<unsigned int>> regionCounts;
  public:
    TableCount(IntervalList &g): intervalList(g), nColumns(0) {}
    void addCounter(Counter &counter) {
      auto &counts = counter.getCounts();
      for (auto &count: counts) {
        auto p = regionCounts.find(count.first);
        if (p == regionCounts.end()) {
          regionCounts[count.first] = vector <unsigned int> (nInputs, 0);
          vector <unsigned int> v (nInputs, 0);
          v[nColumns] = round(count.second);
          regionCounts[count.first] = v;
        }
        else {
          p->second[nColumns] = round(count.second);
        }
      }
      ++nColumns;
    }
    void dump(ostream &outputFile, vector <string> &samples) {
      vector < vector < size_t > > lineNames;
      string lineName;
      for (auto i: regionCounts) {
        lineNames.push_back(i.first);
      }
      sort(lineNames.begin(), lineNames.end());
      outputFile << "Type";
      for (string &sample: samples) {
        outputFile << "\t" << sample;
      }
      outputFile << "\n";
      for (auto &l: lineNames) {
        vector < string > s;
        for (size_t i: l) s.push_back(config->getName(i));
        join(s, lineName, "--");
        outputFile << lineName;
        vector < unsigned int > &values = regionCounts[l];
        for (unsigned int v: values) {
          outputFile << "\t" << v;
        }
        outputFile << "\n";
      }
    }
};

inline void printUsage () {
  cerr << "Usage: mmquant [options]\n";
  cerr <<   "\tCompulsory options:\n";
  cerr <<     "\t\t-a file: annotation file in GTF format\n";
  cerr <<     "\t\t-r file1 [file2 ...]: reads in BAM/SAM format\n";
  cerr << "\tMain options:\n";
  cerr <<     "\t\t-o output: output file (default: stdout)\n";
  cerr <<     "\t\t-c config_file: configuration file (default: " << defaultConfigFileName << ")\n";
  cerr <<     "\t\t-n name1 name2...: short name for each of the reads files\n";
  cerr <<     "\t\t-s strand: string (U, F, R, FR, RF, FF, defaut: F) (use several strand types if the library strategies differ)\n";
  cerr <<     "\t\t-f format (SAM or BAM): format of the read files (default: guess from file extension)\n";
  cerr <<     "\t\t-l integer: overlap type (<0: read is included, <1: % overlap, otherwise: # nt, default: " << overlap << ")\n";
  cerr <<     "\t\t-u: reads are unsorted (default: false)\n";
  cerr <<     "\t\t-d integer: upstream region size (default: " << upstreamSize << ")\n";
  cerr <<     "\t\t-D integer: downstream region size (default: " << downstreamSize << ")\n";
  cerr <<     "\t\t-y string: quantification strategy, valid values are: default, unique, random, ratio (default: default)\n";
  cerr <<     "\t\t-e integer: attribute a read to a feature if at least N% of the hits map to the feature (default: " << static_cast<unsigned int>(rescueThreshold*100) << "%)\n";
  cerr << "\tOutput options:\n";
  cerr <<     "\t\t-p: print progress\n";
  cerr <<     "\t\t-m file: print mapping statistics for each read (slow, only work with 1 input file)\n";
  cerr <<     "\t\t-M file: print mapping statistics for each interval (slow, only work with 1 input file)\n";
  cerr <<     "\t\t-t integer: # threads (default: " << nThreads << ")\n";
  cerr <<     "\t\t-h: this help" << endl;
}

int main(int argc, char **argv) {
  overlap                   = -1.0;
  sorted                    = true;
  rescueThreshold           = 1.0;
  intervalOverlapFunction   = intervalInclusion;
  printReadStatsFunction    = printReadStatsVoid;
  rescueFunction            = rescueVoid;
  strategy                  = Strategy::DEFAULT;
  string gtfFileName, outputFileName, configFileName = defaultConfigFileName;
  vector <string> readsFileNames, names;
  if (argc == 1) {
    printUsage();
    return EXIT_SUCCESS;
  }
  for (int i = 1; i < argc; i++) {
    string s(argv[i]);
    if (! s.empty()) {
      if (s == "-a") {
        gtfFileName = string(argv[++i]);
      }
      else if (s == "-r") {
        for (++i; i < argc; ++i) {
          string t(argv[i]);
          if (t[0] == '-') {--i; break;}
          else readsFileNames.push_back(t);
        }
      }
      else if (s == "-n") {
        for (++i; i < argc; ++i) {
          string t(argv[i]);
          if (t[0] == '-') {--i; break;}
          else names.push_back(t);
        }
      }
      else if (s == "-c") {
        configFileName = argv[++i];
      }
      else if (s == "-o") {
        outputFileName = string(argv[++i]);
      }
      else if (s == "-l") {
        overlap = stof(argv[++i]);
        if      (overlap < 0.0) intervalOverlapFunction = intervalInclusion;
        else if (overlap < 1.0) intervalOverlapFunction = intervalOverlapPc;
        else                    intervalOverlapFunction = intervalOverlap;
      }
      else if (s == "-s") {
        for (++i; i < argc; ++i) {
          s = argv[i];
          if      (s == "U")  { strandednesses.push_back(Strandedness::U);  strandednessFunctions.push_back(strandU); }
          else if (s == "F")  { strandednesses.push_back(Strandedness::F);  strandednessFunctions.push_back(strandF); }
          else if (s == "R")  { strandednesses.push_back(Strandedness::R);  strandednessFunctions.push_back(strandR); }
          else if (s == "FR") { strandednesses.push_back(Strandedness::FR); strandednessFunctions.push_back(strandFR); }
          else if (s == "FF") { strandednesses.push_back(Strandedness::FF); strandednessFunctions.push_back(strandFF); }
          else if (s == "RF") { strandednesses.push_back(Strandedness::RF); strandednessFunctions.push_back(strandRF); }
          else if (s[0] == '-') {--i; break;}
          else if (s.empty())   {--i; break;}
          else {
            cerr << "Do not understand strandedness " << s << "\n" << "Exiting." << endl;
            printUsage();
            return EXIT_FAILURE;
          }
        }
      }
      else if (s == "-p") {
        progress = true;
      }
      else if (s == "-t") {
        nThreads = stoi(argv[++i]);
      }
      else if (s == "-m") {
        readStatsFile.open(string(argv[++i]));
        printReadStatsFunction = printReadStats;
        readStats              = true;
      }
      else if (s == "-M") {
        intervalStatsFile.open(string(argv[++i]));
        intervalStats = true;
      }
      else if (s == "-f") {
        for (++i; i < argc; ++i) {
          s = argv[i];
          lower(s);
          if      (s == "sam")  { formats.push_back(ReadsFormat::SAM); }
          else if (s == "bam")  { formats.push_back(ReadsFormat::BAM); }
          else if (s.empty())   {--i; break;}
          else if (s[0] == '-') {--i; break;}
          else {
            cerr << "Do not understand reads format " << s << "\n" << "Exiting." << endl;
            printUsage();
            return EXIT_FAILURE;
          }
        }
      }
      else if (s == "-e") {
        rescueThreshold = stof(argv[++i])/100.0;
        if (rescueThreshold < 1.0) rescueFunction = rescue;
      }
      else if (s == "-d") {
        upstreamSize = stoul(argv[++i]);
      }
      else if (s == "-D") {
        downstreamSize = stoul(argv[++i]);
      }
      else if (s == "-y") {
        s = argv[++i];
        lower(s);
        if      (s == "default") { strategy = Strategy::DEFAULT; }
        else if (s == "unique")  { strategy = Strategy::UNIQUE; }
        else if (s == "random")  { strategy = Strategy::RANDOM; }
        else if (s == "ratio")   { strategy = Strategy::RATIO; }
        else {
          cerr << "Do not understand strategy " << s << "\n" << "Exiting." << endl;
          printUsage();
          return EXIT_FAILURE;
        }
      }
      else if (s == "-u") {
        sorted = false;
      }
      else if (s == "-h") {
        printUsage();
        return 0;
      }
      else {
        cerr << "Error: wrong parameter '" << s << "'.\nExiting." << endl;
        printUsage();
        return EXIT_FAILURE;
      }
    }
  }
  if (gtfFileName.empty()) {
    cerr << "Missing input GTF file.\nExiting." << endl;
    printUsage();
    return EXIT_FAILURE;
  }
  if (readsFileNames.empty()) {
    cerr << "Missing input BAM file.\nExiting." << endl;
    printUsage();
    return EXIT_FAILURE;
  }
  nInputs = readsFileNames.size();
  if (names.empty()) {
    for (string &fileName: readsFileNames) {
      string n = fileName;
      size_t p = n.find_last_of("/");
      if (p != string::npos) n = n.substr(p+1); 
      p = n.find_last_of(".");
      if (p != string::npos) n = n.substr(0, p); 
      names.push_back(n);
    }
  }
  else if (names.size() != nInputs) {
    cerr << "Number of names is not equal to number of file names.\nExiting." << endl;
    printUsage();
    return EXIT_FAILURE;
  }
  if (strandednesses.size() == 0) {
    strandednesses        = vector <Strandedness>         (nInputs, Strandedness::F);
    strandednessFunctions = vector <StrandednessFunction> (nInputs, strandF);
  }
  else if (strandednesses.size() == 1) {
    if (nInputs != 1) {
      strandednesses        = vector <Strandedness>         (nInputs, strandednesses.front());
      strandednessFunctions = vector <StrandednessFunction> (nInputs, strandednessFunctions.front());
    }
  }
  else if (strandednesses.size() != nInputs) {
    cerr << "Number of strandedness is not equal to number of file names.\nExiting." << endl;
    printUsage();
    return EXIT_FAILURE;
  }
  if (formats.size() == 0) {
    formats = vector <ReadsFormat> (nInputs, ReadsFormat::UNKNOWN);
  }
  else if (formats.size() == 1) {
    if (nInputs != 1) {
      formats = vector <ReadsFormat> (nInputs, formats.front());
    }
  }
  else if (formats.size() != nInputs) {
    cerr << "Number of reads formats is not equal to number of file names.\nExiting." << endl;
    printUsage();
    return EXIT_FAILURE;
  }
  if ((readStats || intervalStats) && (nInputs != 1)) {
    cerr << "Only one reads file when providing reads or interval statistics.\nExiting." << endl;
    printUsage();
    return EXIT_FAILURE;
  }
  locale comma_locale(std::locale(), new comma_numpunct());
  cerr.imbue(comma_locale);
  ofstream of;
  streambuf *buf;
  if (outputFileName.empty()) {
    buf = cout.rdbuf();
  }
  else {
    of.open(outputFileName.c_str());
    buf = of.rdbuf();
  }
  ostream outputFile(buf);
  config = new Config();
  config->parse(configFileName);
  IntervalList intervalList (gtfFileName);
  TableCount table (intervalList);
  if (nThreads == 1) {
    Counter counter (intervalList);
    for (unsigned int i = 0; i < nInputs; i++) {
      counter.clear();
      counter.read(readsFileNames[i], strandednesses.front(), strandednessFunctions.front(), formats.front());
      counter.dump();
      table.addCounter(counter);
    }
  }
  else {
    vector < thread > threads(nThreads);
    atomic < unsigned int > i(0);
    mutex m1, m2;
    for (thread &t: threads) {
      //t = thread([&intervalList, &m, &i, &readsFileNames, &names, &table]() {
      t = thread([&]() {
          Counter counter (intervalList);
          while (i < nInputs) {
          unsigned int thisI;
          m1.lock();
          thisI = i;
          ++i;
          m1.unlock();
          counter.clear();
          counter.read(readsFileNames[thisI], strandednesses[thisI], strandednessFunctions[thisI], formats[thisI]);
          m2.lock();
          counter.dump();
          m2.unlock();
          table.addCounter(counter);
        }
      });
    }
    for (thread &t: threads) {
      t.join();
    }
  }
  table.dump(outputFile, names);
  if (readStats)     readStatsFile.close();
  if (intervalStats) intervalStatsFile.close();
  cerr << "Successfully done." << endl;
  return 0;
}
