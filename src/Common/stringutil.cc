#include  "stringutil.h"
#include  "gostream.h"

#include  <iostream>
#include  <fstream>
#include  "stlio.h"

#ifdef __OLDCOMPILER__
#include  <strstream>
#include  <stdio.h>
#define STRINGSTREAM  strstream 
#define ISTRINGSTREAM istrstream 
#else
#include  <sstream>
#define STRINGSTREAM  std::stringstream
#define ISTRINGSTREAM std::istringstream
#endif

/*--------------------------------------*/

std::string Int2String   (int a   )
{
#ifdef __OLDCOMPILER__
  char c[1+(int) log(a)];
  sprintf(c,"%d",a);
  std::string r = c;
  return r;
#else
  STRINGSTREAM ss;
  ss << a;
  std::string r;
  r = ss.str();
  return r;
#endif
}

/*--------------------------------------*/

std::string Double2String(double a)
{
  STRINGSTREAM ss;
  ss << a;
  std::string r;
  r = ss.str();
  return r;
}

/*--------------------------------------*/

std::string GetBase(const char* buf, char sep)
{
  std::vector<std::string> all= StringSplit(buf,sep);
  assert(all.size());
  return all[0];
}

std::string GetTail(const char* buf, char sep)
{
  std::vector<std::string> all= StringSplit(buf,sep);
  assert(all.size());
  return all[all.size()-1];
}

/*--------------------------------------*/

std::vector<std::string> StringSplit(const char* buf, char sep)
{
  std::vector<std::string> words;
  ISTRINGSTREAM is(buf);
  while (!is.eof())
    {
      std::string t;
      getline(is,t,sep);
      if(t.size()==0) continue;
      if(t[0]!=sep) words.push_back(t);
    }
  return words;
}

std::vector<std::string> StringSplit(const char* buf, char sep1, char sep2)
{
  std::vector<std::string> words;
  ISTRINGSTREAM is(buf);
  while (!is.eof())
    {
      std::string t;
      getline(is,t,sep1);
      if(t.size()==0) continue;
      if(t[0]!=sep1) words.push_back(t);
    }
  std::vector<std::string> words2;
  for(int i=0;i<words.size();i++)
    {
      std::vector<std::string> s1 = StringSplit(words[i].c_str(),sep2);
      copy(s1.begin(),s1.end(),back_inserter(words2));
    }
      
  return words2;
}

/*-----------------------------------------*/

std::pair<std::string,std::vector<std::string> > SplitArgs(std::string s)
{
  std::vector<std::string> vs = StringSplit(s.c_str(),'_');
  if(!vs.size())
    {
      std::cerr << "SplitArgs()\t";
      std::cerr << s << " --> " <<  vs << std::endl;
      abort();
    }
  std::string name = vs[0];
  std::vector<std::string> args(vs.size()-1);
  for(int i=1;i<vs.size();i++)
    {
      args[i-1] = vs[i];
    }
  return make_pair(name,args);
}
