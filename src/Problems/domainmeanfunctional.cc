#include  "domainmeanfunctional.h"
#include  "zerodirichletdata.h"
#include  "constantrighthandside.h"

using namespace Gascoigne;

// /*-----------------------------------------*/

DomainMeanFunctional::DomainMeanFunctional(const std::vector<std::string>& args) : DomainFunctional()
{
  _x0 = 0.; _x1 = 0.; _y0 = 0.; _y1 = 0.; _z0 = 0.; _z1 = 0.;
  RightHandSideData* RH;
  _comp  = atoi(args[1].c_str());
  _ncomp = atoi(args[0].c_str());
  if (args.size()==2)
    {
      _domain = "all";

      RH = new OneComponentRightHandSideData(_ncomp,_comp);
    }
  else if (args.size()==6)
    {
      _domain = "subdomain";
      _x0 = atof(args[2].c_str());
      _x1 = atof(args[3].c_str());
      _y0 = atof(args[4].c_str());
      _y1 = atof(args[5].c_str());
      RH = new RectangleRightHandSideData(_ncomp,_comp,_x0,_x1,_y0,_y1);
    }
  else if (args.size()==8)
    {
      _domain = "subdomain";
      _x0 = atof(args[2].c_str());
      _x1 = atof(args[3].c_str());
      _y0 = atof(args[4].c_str());
      _y1 = atof(args[5].c_str());
      _z0 = atof(args[6].c_str());
      _z1 = atof(args[7].c_str());
      RH = new RectangleRightHandSideData(_ncomp,_comp,_x0,_x1,_y0,_y1,_z0,_z1);
    }
  else
    {
      std::cerr << "DomainMeanFunctional::DomainMeanFunctional()\n";
      std::cerr << "wrong number of arguments: " << args.size() << std::endl;
   }  
}

/*-----------------------------------------*/

double DomainMeanFunctional::J(const FemFunction& U, const Vertex2d& V) const
{
  if (_domain!="subdomain")
    return U[_comp].m();

  if ((V.x()>_x1) || (V.x()<_x0)) return 0.;
  if ((V.y()>_y1) || (V.y()<_y0)) return 0.;

  return U[_comp].m();
}

/*-----------------------------------------*/

double DomainMeanFunctional::J(const FemFunction& U, const Vertex3d& V) const
{
  if (_domain!="subdomain")
    return U[_comp].m();

  if ((V.x()>_x1) || (V.x()<_x0)) return 0.;
  if ((V.y()>_y1) || (V.y()<_y0)) return 0.;
  if ((V.z()>_z1) || (V.z()<_z0)) return 0.;

  return U[_comp].m();
}
