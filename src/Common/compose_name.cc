#include "compose_name.h"
#include "stdio.h"

void compose_name(std::string& s, int i)
{
  char cname[100];
  sprintf(cname,"%s.%05d",s.c_str(),i);
  s = cname;  
}

void compose_name(std::string& s, double d, std::string f)
{
  char cname[30];
  std::string format("%s");
  format += f;
  sprintf(cname,format.c_str(),s.c_str(),d);
  s = cname;  
}

void compose_name(std::string& s, int i, std::string t)
{
  char cname[30];
  std::string format("%s");
  format += t;
  //  format += "%05d";
  sprintf(cname,format.c_str(),s.c_str(),i);
  s = cname;  
}

void compose_name(std::string& s, int i, int l)
{
  char ll[1];
  sprintf(ll,"%01d",l);
  std::string format("%s.%0");
  format += ll;
  format += "d";
  char cname[30];
  sprintf(cname,format.c_str(),s.c_str(),i);
  s = cname;  
}

void compose_name_without_dot(std::string& s, int i)
{
  char cname[30];
  sprintf(cname,"%s%03d",s.c_str(),i);
  s = cname;  
}
