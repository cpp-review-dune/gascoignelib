#include "curvedshapes.h"
#include "splinecontour.h"
#include "runderkreis.h"

/******************************************************/

CurvedShapes::~CurvedShapes()
{
  std::map<int,BoundaryFunction<2>* > ShapeOfColor2d;
  std::map<int,BoundaryFunction<3>* > ShapeOfColor3d;
  colors.clear();
  std::map<int,BoundaryFunction<2>* >::iterator it2d;
  std::map<int,BoundaryFunction<3>* >::iterator it3d;
  for (it2d=ShapeOfColor2d.begin();it2d!=ShapeOfColor2d.end();++it2d)
    if (it2d->second) delete it2d->second; it2d=NULL;
  for (it3d=ShapeOfColor3d.begin();it3d!=ShapeOfColor3d.end();++it3d)
    if (it3d->second) delete it3d->second; it3d=NULL;
  ShapeOfColor2d.clear();
  ShapeOfColor3d.clear();
}

/******************************************************/

CurvedShapes::CurvedShapes()
{
  colors.clear();
}

/******************************************************/

int CurvedShapes::Curved(int col) const 
{ 
  if (colors.find(col)!=colors.end()) return 1; 
  return 0;
}

/******************************************************/

void CurvedShapes::ReInit(const std::vector<BoundaryLine>& blines)
{
}

/******************************************************/

void CurvedShapes::ReInit(const std::vector<BoundaryQuad>& bquads)
{
}

/******************************************************/

void CurvedShapes::BasicInit(const std::vector<std::string>& names)
{
  if (names.size()<1) return;

  int col = atoi(names[0].c_str());

  if (names.size()<2) return;

  std::string shapename = names[1];

  if (shapename=="Kreis")
    {
      assert(names.size()==5);
      double x = atof(names[2].c_str());
      double y = atof(names[3].c_str());
      double r = atof(names[4].c_str());

      Vertex2d V(x,y);
      ShapeOfColor2d[col] = new RunderKreis<2>(V,r);

      colors.insert(col);
    }
  else if (shapename=="Oval")
    {
      assert(names.size()==6);
      double x = atof(names[2].c_str());
      double y = atof(names[3].c_str());
      double r1 = atof(names[4].c_str());
      double r2 = atof(names[5].c_str());

      Vertex2d V(x,y);
      Vertex2d R(r1,r2);
      ShapeOfColor2d[col] = new Oval<2>(V,R);

      colors.insert(col);
    }
  else if ((shapename=="Kugel") || (shapename=="FlacherKreis"))
    {
      assert(names.size()==6);
      double x = atof(names[2].c_str());
      double y = atof(names[3].c_str());
      double z = atof(names[4].c_str());
      double r = atof(names[5].c_str());
      Vertex3d V(x,y,z);
      if (shapename=="Kugel")
	{
	  ShapeOfColor3d[col] = new RunderKreis<3>(V,r);
	}
      else if (shapename=="FlacherKreis")
	{
	  ShapeOfColor3d[col] = new FlacherKreis(V,r);
	}
      colors.insert(col);
    }
  else
    {
      std::cerr << "CurvedShapes::init\n";
      std::cerr << "shapename=" << shapename << " ???" << std::endl;
      exit(1);
    }
}

