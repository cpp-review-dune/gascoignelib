#include  "pointfunctional.h"

/*-----------------------------------------*/

PointFunctional::~PointFunctional()
{
}

/*-----------------------------------------*/

PointFunctional::PointFunctional()  
  : Functional(), v(0), w(0) {}

/*-----------------------------------------*/

PointFunctional::PointFunctional
(const Equation& EQ, const std::vector<std::string>& args)  
  : Functional(), v(0), w(0)
{
  Init(EQ,args);
}

/*-----------------------------------------*/


void PointFunctional::Init(const Equation& EQ, const std::vector<std::string>& args)
{
  type = args[0];
      
  // syntax z.B:
  // MachnenPunkt 2  point vertex_1_2_0._0.5

  if(type=="vertex")
    {
      int n = atoi(args[1].c_str());
      assert(n>0);
      v.resize(n);
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
	  v[ii].x() = x;  v[ii].y() = y;
	  w[ii] = w0;
	  //std::cerr << "point: " << w0 << "\t" << v[ii] << std::endl;
	}
      //      ConstructByVertex(v,w);

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
      //      ConstructById(ids,w);

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

// void PointFunctional::ConstructByVertex
// (const std::vector<Vertex2d>& v0, const std::vector<double>& w0)
// {
//   v.resize(v0.size());
//   w.resize(w0.size());
//   v = v0;
//   w = w0;
// }

/*-----------------------------------------*/

// void PointFunctional::ConstructById(const std::vector<int>& ids0, const std::vector<double>& w0)
// {
//   ids.resize(ids0.size());
//   w.resize(w0.size());
//   ids = ids0;
//   w = w0;
// }
