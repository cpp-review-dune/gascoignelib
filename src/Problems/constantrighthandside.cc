#include  "constantrighthandside.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

ConstantRightHandSideData::ConstantRightHandSideData
(const int ncomp, const int comp, const double d) 
  : DomainRightHandSide()
{
  _ncomp = ncomp;
  _comp  = comp;
  _d     = d;
}

ConstantRightHandSideData::ConstantRightHandSideData
(const vector<string>& args) 
  : DomainRightHandSide()
{
  if(args.size()!=3)
    {
      cerr << "ConstantRightHandSideData::ConstantRightHandSideData()\n";
      cerr << "wrong number of args (3): " << args.size() << endl;
      abort();
    }
  _ncomp = atoi(args[0].c_str());
  _comp  = atoi(args[1].c_str());
  _d     = atof(args[2].c_str());
}

/*-----------------------------------------*/

double ConstantRightHandSideData::operator()(int c, const Vertex2d& v)const 
{
  if(c==_comp) return _d;
  return 0.;
}

/*-----------------------------------------*/

double ConstantRightHandSideData::operator()(int c, const Vertex3d& v)const 
{
  if(c==_comp) return _d;
  return 0.;
}
}
