#ifndef __dataformathandler_h
#define __dataformathandler_h

#include  <set>
#include  <map>
#include  <utility>
#include  <vector>
#include  "nvector.h"

#include  <string>
#include  "fixarray.h"
#include  "gascoigne.h"

using namespace Gascoigne;
using namespace std;

/***************************************************/

class DataFormatHandler
{
private:

  typedef pair<string,string>  NameType;
  typedef pair<string,double>  StringDouble;

  typedef map<string,string*>  TypeString;
  typedef map<string,int*>     TypeInt;
  typedef map<string,bool*>     TypeBool;//neu
  typedef map<string,float*>   TypeFloat;
  typedef map<string,double*>  TypeDouble;
  typedef map<string,StringDouble*>  TypeStringDouble;

  typedef map<string,fixarray<2,double>*>     TypeFix2Double;
  typedef map<string,fixarray<3,double>*>     TypeFix3Double;

  typedef map<string,vector<double>*>        TypeVectorDouble;
  typedef map<string,IntVector*>           TypeVectorInt;
  typedef map<string,vector<string>*>         TypeVectorString;

  typedef map<string,set<int>*>           TypeSetInt;
  typedef map<string,set<vector<string> >*>  TypeSetVectorString;

  typedef map<string,map<int,IntVector >*> TypeMapIntVectorInt;

  set<NameType>        NT;
  TypeString           TS;
  TypeInt              TI;
  TypeBool             TB;//neu
  TypeFloat            TF;
  TypeDouble           TD;
  TypeFix2Double       TF2D;
  TypeFix3Double       TF3D;

  TypeVectorDouble    TND;
  TypeVectorInt       TNI;
  TypeVectorString     TVS;

  TypeSetInt           TSI;
  TypeSetVectorString  TSVS;

  TypeMapIntVectorInt TMINI;

  TypeStringDouble     TSD;

  string search(string& fo, const string& name);

public:

  // without default values
  void insert(const string&, string*);
  void insert(const string&, int*);
  void insert(const string&, bool*); //neu
  void insert(const string&, float*);
  void insert(const string&, double*);

  void insert(const string&, fixarray<2,double>*);
  void insert(const string&, fixarray<3,double>*);

  void insert(const string&, vector<double>*);
  void insert(const string&, IntVector*);
  void insert(const string&, vector<string>*);

  void insert(const string&, set<int>*);
  void insert(const string&, set<vector<string> >*);

  void insert(const string&, map<int,IntVector >*);
  void insert(const string&, map<int,string>*);

  void insert(int, StringDouble*);

  // with default values
  void insert(const string&, string*, const string&);
  void insert(const string&, int*   , int);
  void insert(const string&, bool*   , bool);//neu
  void insert(const string&, float*, float);
  void insert(const string&, double*, double);
  void insert(const string&, fixarray<2,double>*, fixarray<2,double>&);
  void insert(const string&, fixarray<3,double>*, fixarray<3,double>&);
  void insert(const string&, vector<double>*, vector<double>&);
  void insert(const string&, IntVector*, IntVector&);

  void get(string&, const string&);

  void setvalue(const string&, const string&);
  void setvalue(const string&, int);
  void setvalue(const string&, float);
  void setvalue(const string&, double);

  void setvalue(const string&, fixarray<2,double>&);
  void setvalue(const string&, fixarray<3,double>&);

  void setvalue(const string&, vector<double>&);
  void setvalue(const string&, IntVector&);
  void setvalue(const string&, vector<string>&);

  void setvalue(const string&, set<int>&);

  void setvalue(const string&, pair<int,IntVector >&);

  void setvalue(const string&, StringDouble&);

  void insertvalue(const string&, vector<string>&);

  void print(ostream&) const;
};

#endif
