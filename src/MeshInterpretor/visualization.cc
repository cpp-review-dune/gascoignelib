#include  "visualization.h"
#include  <set>
#include  "errormacros.h"
#include  "compose_name.h"
#include  "filescanner.h"

/********************************************************************/

Visualization::Visualization() 
  : mesh(0), PointData(0), CellData(0),
    filename("none"), directory("."), GP(0)
{
  init();
}

/********************************************************************/

Visualization::Visualization(const Visualization& W) : 
  mesh(W.mesh)
{
  *this = W;
}

/********************************************************************/

Visualization::~Visualization()
{}

/********************************************************************/

Visualization& Visualization::operator=(const Visualization& W)
{
  stepname  = W.stepname;
  directory = W.directory;
  title     = W.title;
  avsa      = W.avsa;
  gmva      = W.gmva;
  vua       = W.vua;
  vigiea    = W.vigiea;
  teca      = W.teca;
  vtka      = W.vtka;
  //gnua      = W.gnua;
  pstep     = W.pstep;
  time      = W.time;
  tstep     = W.tstep;
  nexttime  = W.nexttime;
  showoutput = W.showoutput;
}

/********************************************************************/

void Visualization::init()
{
  set_directory(".");

  pstep = 1;
  avsa = gmva = vua = vigiea = gnua = teca = 0;
  time = 0.; tstep = 0.; nexttime = 0.;
  vtka = showoutput = 1;
  stepname = "solve";
}

/********************************************************************/

void Visualization::read_parameters(const std::string& dname, const std::string& fname)
{
  set_directory(dname);

  double time;

  std::vector<std::string>  planes(0);
  nvector<double> gnupos(0);

  DataFormatHandler DH;
  DH.insert("step"     ,& pstep     ,1);
  DH.insert("tstep"    ,& time       ,0.);
  DH.insert("showoutput",& showoutput,0);
  DH.insert("gnuplot"  ,& planes);
  DH.insert("gnuposition",& gnupos);
  DH.insert("vtk"      ,& vtka, 1);
  DH.insert("gmv"      ,& gmva,0);
  DH.insert("vu"       ,& vua,0);
  DH.insert("vigie"    ,& vigiea,0);
  DH.insert("gnu"      ,& gnua,0);
  DH.insert("tecplot"  ,& teca,0);
  DH.insert("avs"      ,& avsa,0);

  std::string ffname = directory;
  ffname += "/";
  ffname += fname;
  // BUBU das stimmt was mit dem pfad nicht.
  ffname=fname;

  FileScanner FS(DH,ffname,"Visualization");

  std::string stepname2 = directory;
  stepname2 += "/";
  stepname2 += stepname;
  stepname = stepname2;

  if (time>0.) set_tstep(time);

  // Gnuplot
  
  if (gnupos.size())
    {
      Vertex3d V(gnupos[0],gnupos[1],gnupos[2]);
      std::vector<GnuplotData> vgp;
      for (int i=0; i<planes.size(); i++)
	{
	  vgp.push_back(GnuplotData(planes[i],V));
	}
      if (vgp.size())
	{
	  set_gnuplotdata(vgp);
	}
    }
}

/********************************************************************/

void Visualization::set_directory(const std::string& s) 
{ 
  directory   = s;
}

/********************************************************************/

void Visualization::set_name(const std::string& s) 
{ 
  stepname = directory;
  stepname += "/";
  stepname += s;
  
  filename = stepname;
}

/********************************************************************/

void Visualization::format(const std::string& s)
{
  if      (s=="vtk")     vtka   = 1;
  else if (s=="gmv")     gmva   = 1;
  else if (s=="vu")      vua    = 1;
  else if (s=="vigie")   vigiea = 1;
  else if (s=="gnu")     gnua   = 1;
  else if (s=="tecplot") teca   = 1;
  else if (s=="avs")     avsa   = 1;
  else assert(0);
}

/********************************************************************/

int Visualization::active(int i) const
{
  if (pstep<0)
    {
      if (time<nexttime) return 0;
    }
  else 
    {
      if ( (pstep==0) || (i%pstep) ) return 0;
    }
  return 1;
}

/********************************************************************/

void Visualization::step(int i)
{
  ii=i;
  if ( (pstep==-1) ||  (!active(i)) )
    {
      if (pstep<0) nexttime = time+tstep;

      return;
    }
  filename = stepname;
  compose_name(filename,i);
}

/********************************************************************/

void Visualization::preview()
{
  if (pstep<0) return;
  
  filename = directory;
  filename += "/";
  filename +="Preview";
}

/********************************************************************/

void Visualization::write()
{
  if(mesh==0)
    {
      std::cerr << "Visualization::write()\n";
      std::cerr << "mesh pointer not set\n";
      abort();
    }
  if(filename=="none")
    {
      std::cerr << "Visualization::write()\n";
      std::cerr << "no filename set [use \"step(i)\" or \"preview()\"]\n";
      abort();
    }

  if (avsa)       avs(filename);
  if (gmva)       gmv(filename);
  if (vua)        vu(filename);
  if (gnua)       gnuplot(filename);
  if (vtka)       vtk(filename);
  if (showoutput) std::cout << "[" << filename << "]\n";
}

/********************************************************************/

int Visualization::CheckPointData() const
{
  if(PointData==NULL) return 0;

  std::set<int> comps;

  int pn = PointData->visucomp();

  // scalars
  for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
    {
      int q = p->second;

      assert(q>=0);
      assert(q<pn);

      comps.insert(q);
    }

  // vectors
  for(VisuDataInfo::viterator p=PointDataInfo->vbegin();p!=PointDataInfo->vend();++p)
    {
      for(int i=0;i<3;i++)
	{
	  int q = p->second[i];

	  assert(q>=-1);
	  assert(q<pn);

	  comps.insert(q);
	}
    }
  return comps.size();
}

/********************************************************************/

int Visualization::CheckCellData() const
{
  if(!CellData) return 0;

  std::set<int> comps;

  // scalars
  for(VisuDataInfo::siterator p=CellDataInfo->sbegin();p!=CellDataInfo->send();++p)
    {
      if( (p->second<0)||(p->second>=CellData->visucomp()))
	{
	  std::cerr << "Visualization::CheckCellData()\n";
	  std::cerr << "scalar does not exist "<<p->second<<std::endl;
	  exit(1);
	}
      comps.insert(p->second);
    }

  // vectors
  for(VisuDataInfo::viterator p=CellDataInfo->vbegin();p!=CellDataInfo->vend();++p)
    {
      for(int i=0;i<3;i++)
	{
	  if( (p->second[i]<-1)||(p->second[i]>=CellData->visucomp()))
	    {
	      std::cerr << "Visualization::CheckCellData()\n";
	      std::cerr << "vector component does not exist "<<p->second[i]<<std::endl;
	      exit(1);
	    }
	  comps.insert(p->second[i]);
	}
    }
  return comps.size();
}

/********************************************************************/

void Visualization::output_vertexs(std::ofstream& file) const
{
  if (mesh->dimension()==2)
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << mesh->vertex2d(i) << std::endl;
	}
    }
  else
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << mesh->vertex3d(i) << std::endl;
	}
    }
}

/********************************************************************/

void Visualization::output_vertexs_by_component(std::ofstream& file) const
{
  if (mesh->dimension()==2)
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << " " << mesh->vertex2d(i).x();
	}
      file << std::endl;
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << " " << mesh->vertex2d(i).y();
	}
    }
  else
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << " " << mesh->vertex3d(i).x();
	}
      file << std::endl;
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << " " << mesh->vertex3d(i).y();
	}
      file << std::endl;
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << " " << mesh->vertex3d(i).z();
	}
    }
}

/********************************************************************/

void Visualization::output_quads(std::ofstream& file, const std::string& text) const
{
  if (mesh->dimension()==2)
    {
      for (int c=0; c<mesh->ncells(); c++)
	{
	  file << text;
	  for(int i=0;i<4;i++)
	    {
	      file << mesh->vertex_of_cell(c,i)+1 << " "; 
	    }
	  file<<std::endl;      
	}
    }
}

/********************************************************************/

void Visualization::output_hexs(std::ofstream& file, const std::string& text) const
{
  if (mesh->dimension()==3)
    {
      for (int c=0; c<mesh->ncells(); c++)
	{
	  file << text;
	  for(int i=0;i<8;i++)
	    {
	      file << mesh->vertex_of_cell(c,i)+1 << " "; 
	    }
	  file<<std::endl;      
	}
    }
}

/********************************************************************/

void Visualization::output_solution(std::ofstream& file, int c) const
{
  for (int ind=0; ind<PointData->visun(); ind++)
    {
      file << PointData->visudata(ind,c) << " ";
    }
}
