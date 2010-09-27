#include "compose_name.h"
#include "stdio.h"

using namespace std;

namespace Gascoigne
{
void compose_name(string& s, int i)
{
  char cname[200];
  sprintf(cname,"%s.%05d",s.c_str(),i);
  s = cname;  
}

void compose_name(string& s, double d, string f)
{
  char cname[200];
  string format("%s");
  format += f;
  sprintf(cname,format.c_str(),s.c_str(),d);
  s = cname;  
}

void compose_name(string& s, int i, string t)
{
  char cname[200];
  sprintf(cname,"%s.%05d.%s",s.c_str(),i,t.c_str());
  s = cname;  
}

void compose_name(string& s, int i, int l)
{
  char ll[1];
  sprintf(ll,"%01d",l);
  string format("%s.%0");
  format += ll;
  format += "d";
  char cname[30];
  sprintf(cname,format.c_str(),s.c_str(),i);
  s = cname;  
}

void compose_name_without_dot(string& s, int i)
{
  char cname[200];
  sprintf(cname,"%s%03d",s.c_str(),i);
  s = cname;  
}
}
