#ifndef __filescanner_h
#define __filescanner_h

#include "dataformathandler.h"

/***************************************************/

class FileScanner
{
  DataFormatHandler& DH;
  string             blocksymbol;
  bool               complain;

  void FormatToValue(const vector<string>& words);
  void print(const string& inputname,
	     const string& blockname) const;
  void _assert(bool b, const vector<string>& words) const;
  
public:
  
  FileScanner(DataFormatHandler& D, const string& s, const string& b="");
  FileScanner(DataFormatHandler& D);
  void readfile(const string&,const string&);
  void NoComplain() { complain=0; }
};

/***************************************************/

#endif
