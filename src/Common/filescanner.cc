#include  "filescanner.h"
#include  "fixarray.h"
#include  "linescanner.h"
#include  "stlio.h"

using namespace std;

/***************************************************/

FileScanner::FileScanner(DataFormatHandler& D) : DH(D)
{
  complain = 1;
  blocksymbol = "//Block";
}

/***************************************************/

FileScanner::FileScanner(DataFormatHandler& D, const string& inputname,
			 const string& blockname) : DH(D)
{
  complain = 1;
  blocksymbol = "//Block";
  readfile(inputname,blockname);
}

/***************************************************/

void FileScanner::_assert(bool b, const vector<string>& words) const
{
  if(!b) cerr << "*** FileScanner:\tWrong number of arguments in row\t" << words << endl;
  assert(b);
}

/***************************************************/

void FileScanner::readfile(const string& inputname, const string& blockname)
{
  LineScanner LS(inputname);

  vector<string> words;
  int nwords = 0;
  bool searchblock = (blockname!="");
  bool blockfound  = !searchblock;
  bool helpfound   = 0;
  string helpname = "HELP_ME";
  
  while (nwords>=0)
    {
      nwords = LS.NextLine(words);

      if (nwords==0) continue;

      if (words[0]!=blocksymbol) continue;

      if (nwords==1)
	{
	  cout << "FileScanner::Block without name" << endl;
	}
      else if (words[1]==helpname)
	{
	  helpfound = 1;
	}
      else if (words[1]==blockname)
	{
	  blockfound = 1;
	  break;
	}
    }
  if (helpfound) print(inputname,blockname);

  if (!blockfound)
    {
      //cout << "FileScanner::missing Block " << blockname << endl;
      return;
    }
  //
  // scanning parameters in block
  //
  while (nwords>=0)
    {
      nwords = LS.NextLine(words);

      if (nwords==0) continue;

      // testing if next block begins
      //
      if (words[0]==blocksymbol) break;
      //
      // testing commentaries
      //
      if (words[0]=="/*") continue;
      if (words[0]=="//") continue;

      if (nwords==1)
	{
	  cout << "where is the parameter \"" << words[0] << "\" ?" << endl;
	  continue;
	}
      FormatToValue(words);
    }
}

/***************************************************/

void FileScanner::FormatToValue(const vector<string>& words)
{
  string fo;
  string keyword = words[0];
  DH.get(fo,keyword);

  /*----------------------------------------------*/

  if (fo=="string")
    {
      DH.setvalue(keyword,words[1]);
    }
  else if (fo=="integer")
    {
      int value = atoi(words[1].c_str());
      DH.setvalue(keyword,value);
    }
  else if (fo=="float")
    {
      float value = atof(words[1].c_str());
      DH.setvalue(keyword,value);
    }
  else if (fo=="double")
    {
      double value = atof(words[1].c_str());
      DH.setvalue(keyword,value);
    }

  /*----------------------------------------------*/

  else if (fo=="fixarray<2,double>")
    {
      fixarray<2,double> value;
      value[0] = atof(words[1].c_str());
      value[1] = atof(words[2].c_str());
      DH.setvalue(keyword,value);
    }

  /*----------------------------------------------*/

  else if (fo=="vector<double>")
    {
      int n = atoi(words[1].c_str());
      _assert(words.size()>n+1,words);
      if (n>0)
	{
	  vector<double> value(n);
	  for (int i=0; i<n; i++)
	    {
	      value[i] = atof(words[i+2].c_str());
	    }
	  DH.setvalue(keyword,value);
	}
    }
  else if (fo=="IntVector")
    {
      int n = atoi(words[1].c_str());
      _assert(words.size()>n+1,words);
      cerr << "-----" << n << " " << words.size() << endl;
      if (n>0)
	{
	  IntVector value(n);
	  for (int i=0; i<n; i++)
	    {
	      value[i] = atoi(words[i+2].c_str());
	    }
	  DH.setvalue(keyword,value);
	}
    }
  else if (fo=="vector<string>")
    {
      int n = atoi(words[1].c_str());
      _assert(words.size()>n+1,words);
      if (n>0)
	{
	  vector<string> value(n);
	  for (int i=0; i<n; i++)
	    {
	      value[i] = words[i+2];
	    }
	  DH.setvalue(keyword,value);
	}
    }

  /*----------------------------------------------*/

  else if (fo=="map<int,IntVector >")
    {
      int col = atoi(words[1].c_str());
      int n   = atoi(words[2].c_str());
      IntVector value(n);
      _assert(words.size()>n+2,words);
      for(int i=0; i<n; i++) value[i] = atoi(words[i+3].c_str());
      pair<int,IntVector > p = make_pair(col,value);
      DH.setvalue(keyword,p);
    }

  /*----------------------------------------------*/

  else if (fo=="set<vector<string> >")
    {
      int n = atoi(words[1].c_str());
      vector<string> value(n);
      _assert(words.size()>n+1,words);
      for(int i=0; i<n; i++) value[i] = words[i+2];

      DH.insertvalue(keyword,value);
    }

  else if (fo=="set<int>")
    {
      int n = atoi(words[1].c_str());
      set<int> value;
      _assert(words.size()>n+1,words);
      for(int i=0; i<n; i++) 
	{
	  int v = atoi(words[i+2].c_str());
	  value.insert(v);
	}

      DH.setvalue(keyword,value);
    }

  /*----------------------------------------------*/

  else if (fo=="StringDouble")
    {
      pair<string,double> p = make_pair(words[1], atof(words[2].c_str()));
      DH.setvalue(keyword,p);
    }
  else if ((fo=="") && complain)
    {
      cerr << " ********  FileScanner::not found " << keyword << endl;
    }
}

/***************************************************/

void FileScanner::print(const string& inputname,
			const string& blockname) const
{
  cout << "=====================" << endl;
  cout << blocksymbol << " " << blockname << endl << endl;
  DH.print(cout);
}
