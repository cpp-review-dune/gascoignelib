#ifndef __linescanner_h
#define __linescanner_h

#include <fstream>
#include <vector>
#include  <string>

/***************************************************/

namespace Gascoigne
{
class LineScanner
{
  std::ifstream fp;

public:

  LineScanner(const std::string& filename);
  ~LineScanner();

  int NextLine(std::vector<double>& words);
  int NextLine(std::vector<std::string>& words);
  int NextLine(std::vector<std::string>& words, const std::vector<int>& w);

  void split(std::vector<std::string>& words, const char& c) const;
  void split(std::vector<std::string>& words, const std::vector<char>& c) const;
};
}

/***************************************************/

#endif
