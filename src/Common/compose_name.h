#ifndef __compose_name_h
#define __compose_name_h

#include  <string>

void compose_name(std::string&, double,std::string f="%5.2f");
void compose_name(std::string&, int);
void compose_name(std::string&, int, std::string t);
void compose_name(std::string&, int, int);
void compose_name_without_dot(std::string&, int);

#endif
