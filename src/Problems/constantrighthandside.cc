#include  "constantrighthandside.h"

/*-----------------------------------------*/

ConstantRightHandSideData::ConstantRightHandSideData
(const std::vector<std::string>& args) 
  : RightHandSideData()
{
  if(args.size()!=3)
    {
      std::cerr << "ConstantRightHandSideData::ConstantRightHandSideData()\n";
      std::cerr << "wrong number of args (3): " << args.size() << std::endl;
      abort();
    }
  ncomp = atoi(args[0].c_str());
  comp  = atoi(args[1].c_str());
  d     = atof(args[2].c_str());
}

/*-----------------------------------------*/

double ConstantRightHandSideData::operator()(int c, const Vertex2d& v)const 
{
  if(c==comp) return d;
  return 0.;
}

/*-----------------------------------------*/

double ConstantRightHandSideData::operator()(int c, const Vertex3d& v)const 
{
  if(c==comp) return d;
  return 0.;
}
