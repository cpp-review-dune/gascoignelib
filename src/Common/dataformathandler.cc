#include  "dataformathandler.h"
#include  "stdio.h"
#include  "stringutil.h"

using namespace std;

namespace Gascoigne
{

/***************************************************/

void DataFormatHandler::insert(const string& nm, string* pos)
{
  string type = "string";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TS.insert(make_pair(name,pos));
}

void DataFormatHandler::insert(const string& nm, int* pos)
{
  string type = "integer";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TI.insert(make_pair(name,pos));
}

void DataFormatHandler::insert(const string& nm, bool* pos)
{
  string type = "bool";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TB.insert(make_pair(name,pos));
}

void DataFormatHandler::insert(const string& nm, float* pos)
{
  string type = "float";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TF.insert(make_pair(name,pos));
}

void DataFormatHandler::insert(const string& nm, double* pos)
{
  string type = "double";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TD.insert(make_pair(name,pos));
}

/*-----------------------------------------*/

void DataFormatHandler::insert(const string& nm, fixarray<2,double>* pos)
{
  string type = "fixarray<2,double>";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TF2D.insert(make_pair(name,pos));
}

/*-----------------------------------------*/

void DataFormatHandler::insert(const string& nm, fixarray<3,double>* pos)
{
  string type = "fixarray<3,double>";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TF3D.insert(make_pair(name,pos));
}

/*-----------------------------------------*/

void DataFormatHandler::insert(const string& nm, vector<double>* pos)
{
  string type = "vector<double>";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TND.insert(make_pair(name,pos));
}

void DataFormatHandler::insert(const string& nm, IntVector* pos)
{
  string type = "IntVector";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TNI.insert(make_pair(name,pos));
}

void DataFormatHandler::insert(const string& nm, vector<string>* pos)
{
  string type = "vector<string>";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TVS.insert(make_pair(name,pos));
}

/*-----------------------------------------*/

void DataFormatHandler::insert(const string& nm, map<int,IntVector >* pos)
{
  string type = "map<int,IntVector >";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TMINI.insert(make_pair(name,pos));
}

/*-----------------------------------------*/

void DataFormatHandler::insert(const string& nm, IntSet* pos)
{
  string type = "set<int>";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TSI.insert(make_pair(name,pos));
}

/*-----------------------------------------*/

void DataFormatHandler::insert(const string& nm, set<vector<string> >* pos)
{
  string type = "set<vector<string> >";
  string name = nm;
  NameType p = make_pair(name,type);
  NT.insert(p);
  TSVS.insert(make_pair(name,pos));
}

/*-----------------------------------------*/

void DataFormatHandler::insert(int i, StringDouble* pos)
{
  char zahl[3];
  sprintf(zahl,"%01d",i);
  string name = zahl;
  string type = "StringDouble";
  NameType p = make_pair(name,type);
  NT.insert(p);
  TSD.insert(make_pair(name,pos));
}

/***************************************************/

void DataFormatHandler::insert(const string& nm, string* pos, 
			       const string& def)
{
  insert(nm,pos);
  *pos = def;
}

void DataFormatHandler::insert(const string& nm, int* pos,
			       int def)
{
  insert(nm,pos);
  *pos = def;
}

void DataFormatHandler::insert(const string& nm, bool* pos,
			       bool def)
{
  insert(nm,pos);
  *pos = def;
}

void DataFormatHandler::insert(const string& nm, float* pos,
			       float def)
{
  insert(nm,pos);
  *pos = def;
}

void DataFormatHandler::insert(const string& nm, double* pos,
			       double def)
{
  insert(nm,pos);
  *pos = def;
}

void DataFormatHandler::insert(const string& nm, fixarray<2,double>* pos,
			       fixarray<2,double>& def)
{
  insert(nm,pos);
  *pos = def;
}

void DataFormatHandler::insert(const string& nm, fixarray<3,double>* pos,
			       fixarray<3,double>& def)
{
  insert(nm,pos);
  *pos = def;
}

void DataFormatHandler::insert(const string& nm, vector<double>* pos,
			       vector<double>& def)
{
  insert(nm,pos);
  *pos = def;
}

void DataFormatHandler::insert(const string& nm, IntVector* pos,
			       IntVector& def)
{
  insert(nm,pos);
  *pos = def;
}

/***************************************************/


void DataFormatHandler::setvalue(const string& name, const string& value)
{
  TypeString::const_iterator p;
  p = TS.find(name);
  if (p!=TS.end())
    {
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, int value)
{
  TypeInt::const_iterator p;
  p = TI.find(name);
  if (p!=TI.end())
    {
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, bool value)
{
  TypeBool::const_iterator p;
  p = TB.find(name);
  if (p!=TB.end())
    {
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, float value)
{
  TypeFloat::const_iterator p;
  p = TF.find(name);
  if (p!=TF.end())
    {
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, double value)
{
  TypeDouble::const_iterator p;
  p = TD.find(name);
  if (p!=TD.end())
    {
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, fixarray<2,double>& value)
{
  TypeFix2Double::const_iterator p;
  p = TF2D.find(name);
  if (p!=TF2D.end())
    {
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, fixarray<3,double>& value)
{
  TypeFix3Double::const_iterator p;
  p = TF3D.find(name);
  if (p!=TF3D.end())
    {
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, vector<double>& value)
{
  TypeVectorDouble::const_iterator p;
  p = TND.find(name);
  if (p!=TND.end())
    {
      p->second->resize(value.size());
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, IntVector& value)
{
  TypeVectorInt::const_iterator p;
  p = TNI.find(name);
  if (p!=TNI.end())
    {
      p->second->resize(value.size());
      *(p->second) = value;
      return;
    }  
}


void DataFormatHandler::setvalue(const string& name, vector<string>& value)
{
  TypeVectorString::const_iterator p;
  p = TVS.find(name);
  if (p!=TVS.end())
    {
      p->second->resize(value.size());
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, IntSet& value)
{
  TypeSetInt::const_iterator p;
  p = TSI.find(name);
  if (p!=TSI.end())
    {
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, pair<int, IntVector >& value)
{
  TypeMapIntVectorInt::const_iterator p;
  p = TMINI.find(name);
  if (p!=TMINI.end())
    {
      p->second->insert(value);
      return;
    }  
}

void DataFormatHandler::setvalue(const string& name, pair<string,double>& value)
{
  TypeStringDouble::const_iterator p;
  p = TSD.find(name);
  if (p!=TSD.end())
    {
      *(p->second) = value;
      return;
    }  
}

void DataFormatHandler::insertvalue(const string& name, vector<string>& value)
{
  TypeSetVectorString::const_iterator p;
  p = TSVS.find(name);
  if (p!=TSVS.end())
    {
      p->second->insert(value);
      return;
    }  
}

/***************************************************/

string DataFormatHandler::search(string& fo, const string& name)
{
  set<NameType>::const_iterator p = NT.find(make_pair(name,fo));
  if (p!=NT.end())
    {
      return p->second;
    }
  string leer = "";
  return leer;
}

/***************************************************/

void DataFormatHandler::get(string& f, const string& name)
{
  vector<string> s(14);
  s[0] = "string";
  s[1] = "integer";
  s[2] = "float";
  s[3] = "double";
  s[4] = "fixarray<2,double>";
  s[5] = "vector<double>";
  s[6] = "IntVector";
  s[7] = "vector<string>";
  s[8] = "map<int,IntVector >";
  s[9] = "set<vector<string> >";
  s[10] = "StringDouble";
  s[11] = "set<int>";
  s[12] = "fixarray<3,double>";
  s[13] = "bool";
  f = "";
  for (int i=0; i<s.size(); i++)
    {
      f = search(s[i],name);
      if (f!="") break;
    }
}

/***************************************************/

void DataFormatHandler::print(ostream& s) const
{
  TypeString::const_iterator p = TS.begin();
  for (p; p!=TS.end(); p++)
    {
      s << p->first;
      s << "     ";
      s << *(p->second);
      s << "\n";
    }
  TypeInt::const_iterator p1 = TI.begin();
  for (p1; p1!=TI.end(); p1++)
    {
      s << p1->first;
      s << "     ";
      s << *(p1->second);
      s << "\n";
    }
  TypeFloat::const_iterator p2 = TF.begin();
  for (p2; p2!=TF.end(); p2++)
    {
      s << p2->first;
      s << "     ";
      s << *(p2->second);
      s << "\n";
    }
  TypeDouble::const_iterator p3 = TD.begin();
  for (p3; p3!=TD.end(); p3++)
    {
      s << p3->first;
      s << "     ";
      s << *(p3->second);
      s << "\n";
    }
  // no values
  TypeFix2Double::const_iterator p4 = TF2D.begin();
  for (p4; p4!=TF2D.end(); p4++)
    {
      s << p4->first;
      s << "\n";
    }
  TypeFix3Double::const_iterator p4a = TF3D.begin();
  for (p4a; p4a!=TF3D.end(); p4a++)
    {
      s << p4a->first;
      s << "\n";
    }
  TypeVectorDouble::const_iterator p5 = TND.begin();
  for (p5; p5!=TND.end(); p5++)
    {
      s << p5->first;
      s << "\n";
    }
  TypeVectorInt::const_iterator p6 = TNI.begin();
  for (p6; p6!=TNI.end(); p6++)
    {
      s << p6->first;
      s << "\n";
    }
  TypeVectorString::const_iterator p7 = TVS.begin();
  for (p7; p7!=TVS.end(); p7++)
    {
      s << p7->first;
      s << "\n";
    }
  TypeMapIntVectorInt::const_iterator p8 = TMINI.begin();
  for (p8; p8!=TMINI.end(); p8++)
    {
      s << p8->first;
      s << "\n";
    }
  TypeSetVectorString::const_iterator p9 = TSVS.begin();
  for (p9; p9!=TSVS.end(); p9++)
    {
      s << p9->first;
      s << "\n";
    }
  TypeStringDouble::const_iterator p10 = TSD.begin();
  for (p10; p10!=TSD.end(); p10++)
    {
      s << p10->first;
      s << "\n";
    }
}

/***************************************************/

}
