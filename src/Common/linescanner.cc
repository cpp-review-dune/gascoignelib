#include  "linescanner.h"
#include  "gascoignemath.h"
#include  "stlio.h"
#include  "stringutil.h"

using namespace std;

namespace Gascoigne
{

/***************************************************/

LineScanner::LineScanner(const string& filename) : 
  fp(filename.c_str()) 
{  
  if(!fp.is_open())
    {
      cout << "LineScanner::cannot open file " << filename << endl;
      abort();
    }
};

/***************************************************/

LineScanner::~LineScanner()
{  
  fp.close();
}

/***************************************************/

void LineScanner::split(vector<string>& words, const char& c) const
{
  // Splits each string of the vector "words"
  // into separate strings due to separator "c"
  // The result is stored again in "words"

  vector<string> help(words);
  words.resize(0);
  for (int i=0; i<help.size(); i++)
    {
      vector<string> s = StringSplit(help[i].c_str(), c);
      for (int j=0; j<s.size(); j++)
        {
          if (s[j]!="") words.push_back(s[j]);
        }      
    }
}

/***************************************************/

void LineScanner::split(vector<string>& words, const vector<char>& c) const
{
  for (int i=0; i<c.size(); i++)
    {
      split(words,c[i]);
    }
}

/***************************************************/

int LineScanner::NextLine(vector<double>& words)
{
  if (fp.eof()) return -1;

  string toto;
  getline(fp,toto);
  vector<string> s1 = StringSplit(toto.c_str(),' ','\t');
  
  words.resize(s1.size());
  for(int i=0;i<s1.size();i++) words[i] = atof(s1[i].c_str());

  return words.size();
}

/***************************************************/

int LineScanner::NextLine(vector<string>& words)
{
  if (fp.eof()) return -1;

  string toto;
  getline(fp,toto);

  words.resize(0);

  vector<string> s1 = StringSplit(toto.c_str(),' ');
  if ( (s1.size()) && (s1[0]!="//Block") )
    {
      if( (toto[0]=='/') && (toto[1]=='/') )
        {
          return 0;
        }
    }
  string tab = "\t";
  for (int i=0; i<s1.size(); i++)
    {
      vector<string> s2 = StringSplit(s1[i].c_str(),tab[0]);
      for (int j=0; j<s2.size(); j++)
        {
          if (s2[j]!="") words.push_back(s2[j]);
        }
    }
  
  return words.size();
}

/***************************************************/

int LineScanner::NextLine(vector<string>& words, const vector<int>& w)
{
  // usefull to read formated FORTRAN data bases
  // fills in "words" the strings of the next line
  // word[i] has fixed lenght w[i]

  if (fp.eof()) return -1;

  string toto;
  getline(fp,toto);

  int j = 0;
  for (int i=0; i<w.size(); i++)
    {
      if (j+w[i]<=toto.size())
        {
          words[i] = toto.substr(j,w[i]);
        }
      else
        {
          words[i] = "";
        }
      j += w[i];
    }
  return words.size(); 
}

/***************************************************/

}
