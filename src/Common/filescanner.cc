#include  "filescanner.h"
#include  "fixarray.h"
#include  "linescanner.h"
#include  "stlio.h"

using namespace std;

namespace Gascoigne
{
  
/***************************************************/

FileScanner::FileScanner(DataFormatHandler& D) : DH(D)
{
  complain = 1;
  blocksymbol = "//Block";
  _i_defaultvalues_level = 0;
}

/***************************************************/

FileScanner::FileScanner(DataFormatHandler& D, const ParamFile* pf, const string& blockname) : DH(D)
{
  complain = 1;
  blocksymbol = "//Block";
  _i_defaultvalues_level = 0;
  readfile(pf,blockname);
}

/***************************************************/

void FileScanner::_assert(bool b, const vector<string>& words) const
{
  if(!b) cerr << "*** FileScanner:\tWrong number of arguments in row\t" << words << endl;
  assert(b);
}

/***************************************************/

void FileScanner::readfile(const ParamFile* pf, const string& blockname)
{
  if(pf==NULL)
    {
      return;
    }

  string inputname = pf->GetName();
  LineScanner LS(inputname);

  vector<string> words;
  int nwords = 0;
  bool searchblock = (blockname!="");
  bool blockfound  = !searchblock;
  bool helpfound   = 0;
  string helpname = "HELP_ME";

  // der folgende block macht es moeglich beliebig viele defaultwert-dateien anzugeben
  // als massnahme gegen schleifen der inklusion gibt es den defaultvalues-zaehler
  // bei jedem weiteren auslesen der defaultwerte der 'naechst-tiefen' datei wird der zaehler
  // um eins hoch gezaehlt
  //syntax:
  //   //Block DefaultValues
  //   files   3   file_a  file_b  file_c
  //   file_a  settingsa.param
  //   file_b  settingsb.param
  //   file_c  settingsc.param
  //
  if(blockname!="DefaultValues"){
     _i_defaultvalues_level++;
     if( 10 < _i_defaultvalues_level ){
       cerr<< __FILE__ << ":" << __LINE__ << ": there seems to be a 'Block DefaultValues' loop"<<endl;
       cerr<< __FILE__ << ":" << __LINE__ << ": last used file: "<< pf->GetName() <<endl;
       abort();
     }
     DataFormatHandler DFH;    
     vector<string>    vs_files;
     string            s_paramfile = pf->GetName();

     DFH.insert("files",&vs_files);  
     FileScanner FS(DFH);
     FS._i_defaultvalues_level = _i_defaultvalues_level;
     FS.NoComplain();
     FS.readfile(pf,"DefaultValues");

     for(int i=0;i<vs_files.size();i++){
       DataFormatHandler DFH2;
       string            s_filename;
       string            s_keyname = vs_files[i];

       DFH2.insert(s_keyname,&s_filename,"none");  
       FileScanner FS(DFH2);
       FS._i_defaultvalues_level = _i_defaultvalues_level;
       FS.NoComplain();
       FS.readfile(pf,"DefaultValues");

       if(s_filename!="none" && s_filename!= s_paramfile){
         ParamFile paramfile(s_filename);
         FileScanner FS2(DH);
         FS2._i_defaultvalues_level = _i_defaultvalues_level;
         FS2.NoComplain();
         FS2.readfile(&paramfile,blockname);
       }
     }
  }    
  
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
  if (helpfound) print(blockname);

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
  string keyword_type;
  string keyword = words[0];
  DH.get(keyword_type,keyword);

  /*----------------------------------------------*/

  if (keyword_type=="string")
    {
      DH.setvalue(keyword,words[1]);
    }
  else if (keyword_type=="integer")
    {
      int value = atoi(words[1].c_str());
      DH.setvalue(keyword,value);
    }
  else if (keyword_type=="bool")
    {
      if(words[1]=="true" || words[1]=="1")
        {
          DH.setvalue(keyword,true);
        }
      else
        {
          DH.setvalue(keyword,false);
        }
    }
  else if (keyword_type=="float")
    {
      float value = atof(words[1].c_str());
      DH.setvalue(keyword,value);
    }
  else if (keyword_type=="double")
    {
      double value = atof(words[1].c_str());
      DH.setvalue(keyword,value);
    }

  /*----------------------------------------------*/

  else if (keyword_type=="fixarray<2,double>")
    {
      fixarray<2,double> value;
      value[0] = atof(words[1].c_str());
      value[1] = atof(words[2].c_str());
      DH.setvalue(keyword,value);
    }

  /*----------------------------------------------*/

  else if (keyword_type=="vector<double>")
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
  else if (keyword_type=="IntVector")
    {
      int n = atoi(words[1].c_str());
      _assert(words.size()>n+1,words);
//       cerr << "-----" << n << " " << words.size() << endl;
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
  else if (keyword_type=="vector<string>")
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

  else if (keyword_type=="map<int,IntVector >")
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

  else if (keyword_type=="set<vector<string> >")
    {
      int n = atoi(words[1].c_str());
      vector<string> value(n);
      _assert(words.size()>n+1,words);
      for(int i=0; i<n; i++) value[i] = words[i+2];

      DH.insertvalue(keyword,value);
    }

  else if (keyword_type=="set<int>")
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

  else if (keyword_type=="StringDouble")
    {
      pair<string,double> p = make_pair(words[1], atof(words[2].c_str()));
      DH.setvalue(keyword,p);
    }
  else if (keyword_type=="" && complain)
    {
      cerr << " ********  FileScanner::not found " << keyword << endl;
    }
}

/***************************************************/

void FileScanner::print(const string& blockname) const
{
  cout << "=====================" << endl;
  cout << blocksymbol << " " << blockname << endl << endl;
  DH.print(cout);
}
}
