#ifndef  __stringutil_h
#define  __stringutil_h

#include  <string>
#include  <vector>

std::string GetBase(const char* buf, char sep='.');
std::string GetTail(const char* buf, char sep='.');
std::vector<std::string> StringSplit(const char* buf, char sep);
std::vector<std::string> StringSplit(const char* buf, char sep1, char sep2);

std::string Int2String   (int a   );
std::string Double2String(double a);

std::pair<std::string,std::vector<std::string> > SplitArgs(std::string s);

#endif
