#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "constantrighthandside.h"

/*-----------------------------------------*/

ResidualFunctional::ResidualFunctional() : _DD(NULL) { _scale = 0.;}

/*-----------------------------------------*/

ResidualFunctional::ResidualFunctional(const Equation& EQ, const std::vector<std::string>& args) : _DD(NULL)
{
  Construct(EQ,args);
}

/*-----------------------------------------*/

ResidualFunctional::~ResidualFunctional()
{
  if(_DD!=NULL) {delete _DD; _DD=NULL;}
}

/*-----------------------------------------*/

void ResidualFunctional::Construct(const Equation& EQ, const std::vector<std::string>& args) 
{
  int ncomp = EQ.ncomp();

  _comp = atoi(args[0].c_str());
  _scale = atof(args[1].c_str());
  _col.clear();
  for (int ii=2; ii<args.size(); ii++) 
    {
      int color = atoi(args[ii].c_str());
      _col.insert(color);
    }
  _DD  = new DirichletDataByColor(GetComp(),GetColors(),GetScale());
}
