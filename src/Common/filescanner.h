#ifndef __filescanner_h
#define __filescanner_h

#include  "dataformathandler.h"
#include  "paramfile.h"

/***************************************************/

class FileScanner
{
  DataFormatHandler& DH;
  std::string        blocksymbol;
  bool               complain;

  void FormatToValue(const std::vector<std::string>& words);
  void print(const std::string& blockname) const;
  void _assert(bool b, const std::vector<std::string>& words) const;
  
public:
  
  FileScanner(DataFormatHandler& D, const Gascoigne::ParamFile* pf, const std::string& b="");
  FileScanner(DataFormatHandler& D);
  void readfile(const Gascoigne::ParamFile* pf, const std::string& blockname);
  void NoComplain() { complain=0; }
};

/***************************************************/

#endif
