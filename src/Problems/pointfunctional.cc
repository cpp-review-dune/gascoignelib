#include  "pointfunctional.h"

/*-----------------------------------------*/

PointFunctional::~PointFunctional()
{
}

/*-----------------------------------------*/

PointFunctional::PointFunctional
(const std::vector<std::string>& args)  
  : Functional(), v2d(0), v3d(0), w(0)
{
  Init(args);
}

/*-----------------------------------------*/


void PointFunctional::Init(const std::vector<std::string>& args)
{
  type = args[0];
      
  // syntax z.B:
  // MachnenPunkt 2  point vertex_1_2_0._0.5

  if(type=="vertex2d")
    {
      int n = atoi(args[1].c_str());
      assert(n>0);
      v2d.resize(n);
      w.resize(n);
      
      //std::cerr << "n = " << n << std::endl;
      int count=2;
      assert(args.size()==2+n*4);
      for(int ii=0; ii<n; ii++)
	{
	  mycomp = int(atof(args[count++].c_str()));
	  double w0 = atof(args[count++].c_str());
	  double x  = atof(args[count++].c_str());
	  double y  = atof(args[count++].c_str());
	  v2d[ii].x() = x;  v2d[ii].y() = y;
	  w[ii] = w0;
	  //std::cerr << "point: " << w0 << "\t" << v[ii] << std::endl;
	}
    }
  else if(type=="vertex3d")
    {
      int n = atoi(args[1].c_str());
      assert(n>0);
      v3d.resize(n);
      w.resize(n);
      
      //std::cerr << "n = " << n << std::endl;
      int count=2;
      assert(args.size()==2+n*5);
      for(int ii=0; ii<n; ii++)
	{
	  mycomp = int(atof(args[count++].c_str()));
	  double w0 = atof(args[count++].c_str());
	  double x  = atof(args[count++].c_str());
	  double y  = atof(args[count++].c_str());
	  double z  = atof(args[count++].c_str());
	  v3d[ii].x() = x;  v3d[ii].y() = y; v3d[ii].z() = z;
	  w[ii] = w0;
	  //std::cerr << "point: " << w0 << "\t" << v3d[ii] << std::endl;
        }
    }
  else if(type=="id")
    {
      int n = atoi(args[1].c_str());
      ids.resize(n);
      w  .resize(n);
      
      std::cerr << "n = " << n << std::endl;
      int count=2;
      for(int ii=0;ii<n;ii++)
	{
	  w  [ii] = atof(args[count++].c_str());
	  ids[ii] = atoi(args[count++].c_str());
	}
      assert(0);
    }
  else
    {
      std::cerr << "PointFunctional::PointFunctional()\n";
      std::cerr << "unknown type: " << type << std::endl;
      abort();
    }
}

/*-----------------------------------------*/

