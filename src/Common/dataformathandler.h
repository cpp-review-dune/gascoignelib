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

/***************************************************/

class DataFormatHandler
{
private:

  typedef std::pair<std::string,std::string>  NameType;
  typedef std::pair<std::string,double>  StringDouble;

  typedef std::map<std::string,std::string*>  TypeString;
  typedef std::map<std::string,int*>     TypeInt;
  typedef std::map<std::string,bool*>     TypeBool;//neu
  typedef std::map<std::string,float*>   TypeFloat;
  typedef std::map<std::string,double*>  TypeDouble;
  typedef std::map<std::string,StringDouble*>  TypeStringDouble;

  typedef std::map<std::string,fixarray<2,double>*>     TypeFix2Double;
  typedef std::map<std::string,fixarray<3,double>*>     TypeFix3Double;

  typedef std::map<std::string,std::vector<double>*>        TypeVectorDouble;
  typedef std::map<std::string,Gascoigne::IntVector*>           TypeVectorInt;
  typedef std::map<std::string,std::vector<std::string>*>         TypeVectorString;

  typedef std::map<std::string,std::set<int>*>           TypeSetInt;
  typedef std::map<std::string,std::set<std::vector<std::string> >*>  TypeSetVectorString;

  typedef std::map<std::string,std::map<int,Gascoigne::IntVector >*> TypeMapIntVectorInt;

  std::set<NameType>        NT;
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

  std::string search(std::string& fo, const std::string& name);

public:

  // without default values
  void insert(const std::string&, std::string*);
  void insert(const std::string&, int*);
  void insert(const std::string&, bool*); //neu
  void insert(const std::string&, float*);
  void insert(const std::string&, double*);

  void insert(const std::string&, fixarray<2,double>*);
  void insert(const std::string&, fixarray<3,double>*);

  void insert(const std::string&, std::vector<double>*);
  void insert(const std::string&, Gascoigne::IntVector*);
  void insert(const std::string&, std::vector<std::string>*);

  void insert(const std::string&, std::set<int>*);
  void insert(const std::string&, std::set<std::vector<std::string> >*);

  void insert(const std::string&, std::map<int,Gascoigne::IntVector >*);
  void insert(const std::string&, std::map<int,std::string>*);

  void insert(int, StringDouble*);

  // with default values
  void insert(const std::string&, std::string*, const std::string&);
  void insert(const std::string&, int*   , int);
  void insert(const std::string&, bool*   , bool);//neu
  void insert(const std::string&, float*, float);
  void insert(const std::string&, double*, double);
  void insert(const std::string&, fixarray<2,double>*, fixarray<2,double>&);
  void insert(const std::string&, fixarray<3,double>*, fixarray<3,double>&);
  void insert(const std::string&, std::vector<double>*, std::vector<double>&);
  void insert(const std::string&, Gascoigne::IntVector*, Gascoigne::IntVector&);

  void get(std::string&, const std::string&);

  void setvalue(const std::string&, const std::string&);
  void setvalue(const std::string&, int);
  void setvalue(const std::string&, bool);//neu
  void setvalue(const std::string&, float);
  void setvalue(const std::string&, double);

  void setvalue(const std::string&, fixarray<2,double>&);
  void setvalue(const std::string&, fixarray<3,double>&);

  void setvalue(const std::string&, std::vector<double>&);
  void setvalue(const std::string&, Gascoigne::IntVector&);
  void setvalue(const std::string&, std::vector<std::string>&);

  void setvalue(const std::string&, std::set<int>&);

  void setvalue(const std::string&, std::pair<int,Gascoigne::IntVector >&);

  void setvalue(const std::string&, StringDouble&);

  void insertvalue(const std::string&, std::vector<std::string>&);

  void print(std::ostream&) const;
};

#endif
