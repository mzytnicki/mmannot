/*
   addNH, add the NH to a read-sorted SAM file.
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
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
using namespace std;

inline void printUsage () {
  cerr << "Usage: addNH [options]\n";
  cerr <<     "\t-i input file (default: stdin)\n";
  cerr <<     "\t-o output file (default: stdout)\n";
  cerr <<     "\t-h: this help" << endl;
}

void printLines(vector <string> &lines, ostream &output) {
  size_t c = lines.size();
  for (string &line: lines) {
    output << line << "\tNH:i:" << c << "\n";
  }
}

int main(int argc, char **argv) {
  ifstream inputFile;
  ofstream outputFile;
  string inputFileName, outputFileName;
  for (int i = 1; i < argc; i++) {
    string s(argv[i]);
    if (! s.empty()) {
      if (s == "-i") {
        inputFileName = string(argv[++i]);
        inputFile.open(inputFileName);
        if (! inputFile) {
          cerr << "Error: Cannot open input file '" << inputFileName << "'." << endl;
          return 1;
        }
      }
      else if (s == "-o") {
        outputFileName = string(argv[++i]);
        outputFile.open(outputFileName);
        if (! outputFile) {
          cerr << "Error: Cannot open output file '" << outputFileName << "'." << endl;
          return 1;
        }
      }
      else if (s == "-h") {
        printUsage();
        return 0;
      }
      else {
        cerr << "Error: wrong parameter '" << s << "'.\nExiting." << endl;
        printUsage();
        return 1;
      }
    }
  }
  istream &input  = (inputFileName.empty())?  cin:  inputFile;
  ostream &output = (outputFileName.empty())? cout: outputFile;
  string line, previousRead, currentRead;
  vector<string> previousLines;
  while (getline(input, line)) {
    if ((line.empty()) || (line[0] == '@')) {
      output << line << "\n";
    }
    else {
      currentRead = line.substr(0, line.find_first_of('\t'));
      if (currentRead == previousRead) {
        previousLines.push_back(line);
      }
      else {
        printLines(previousLines, output);
        previousRead  = currentRead;
        previousLines = {line};
      }
    }
  }
  printLines(previousLines, output);
  return 0;
}
