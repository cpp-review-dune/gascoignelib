#ifndef __filescanner_h
#define __filescanner_h

#include  "dataformathandler.h"
#include  "paramfile.h"

using namespace Gascoigne;

/***************************************************/

class FileScanner
{
  DataFormatHandler& DH;
  string             blocksymbol;
  bool               complain;

  void FormatToValue(const vector<string>& words);
  void print(const string& blockname) const;
  void _assert(bool b, const vector<string>& words) const;
  
public:
  
  FileScanner(DataFormatHandler& D, const ParamFile* pf, const string& b="");
  FileScanner(DataFormatHandler& D);
  void readfile(const ParamFile* pf, const string& blockname);
  void NoComplain() { complain=0; }
};

/***************************************************/

#endif
