#include "hangcontainer3d.h"
#include "find_in_linehang.h"

/*********************************************************************/

HangContainer3d::HangContainer3d(HangList<2>& lh2, HangList<4>& lh3) : 
  HangContainer2d(lh2),
  FaceHanging(lh3) {}

/*********************************************************************/

bool HangContainer3d::ToBeDeleted(const FaceVector& v) const
{
  return (FaceToBeDeleted.find(v) != FaceToBeDeleted.end());
}

/*********************************************************************/

bool HangContainer3d::ToBeCreated(const FaceVector& v) const
{
  return (FaceToBeCreated.find(v) != FaceToBeCreated.end());
}

/*********************************************************************/

void HangContainer3d::clear_hanging_lines()
{
  Hanging          .clear();
  VertexToBeCreated.clear();
  VertexToBeDeleted.clear();
  NotAnyMoreHanging.clear();  
}

/*********************************************************************/

void HangContainer3d::build_hanging_lines(const HangList<2>& oldhangs)
{
  clear_hanging_lines();

  for (HangList<4>::const_iterator  p = FaceHanging.begin(); 
       p!=FaceHanging.end(); p++)
    {
      const fixarray<4,int>& F = p->first;
      fixarray<2,int> edge;
      for (int i=0; i<4; i++)
	{
	  edge[0] = F[i];
	  edge[1] = F[(i+1)%4];

	  if (Hanging.find(edge)==Hanging.end())
	    {
	      // hang muss kreiiert werden
	      HangList<2>::const_iterator old = oldhangs.find(edge);
	      if (old!=oldhangs.end())
		{
		  // alter Hang
		  Hanging.insert(std::make_pair(edge,old->second));		  
		}
	      else
		{
		  int f = p->second.rneighbour();
		  int c = p->second.cneighbour();
		  Hang h(-1,f,c);
		  Hanging.insert(std::make_pair(edge,h));
		}
	    }
	}
    }
  // not any more hanging
  for (HangList<2>::const_iterator  p = oldhangs.begin(); 
       p!=oldhangs.end(); p++)
    {
      if (Hanging.find(p->first)==Hanging.end()) 
	{
	  NotAnyMoreHanging.insert(std::make_pair(p->first,p->second));
	}
    } 
}

/*********************************************************************/

void HangContainer3d::NeighbourSwapper()
{  
  for(HangList<4>::iterator p=FaceToBeDeleted.begin(); p!=FaceToBeDeleted.end(); p++)
    {
      int r = p->second.rneighbour();
      if (r<0)
	{
	  p->second.rneighbour() = p->second.cneighbour();
	  p->second.cneighbour() = r;
	}
    }
}

/*********************************************************************/

void HangContainer3d::make_consistent()
{ 
  HangContainer2d::make_consistent();
  FaceToBeCreated.make_consistent(FaceToBeDeleted);
}

/*********************************************************************/

void HangContainer3d::load_elimination(IntVector& v) const
{
  HangContainer2d::load_elimination(v);
  for(HangList<4>::const_iterator p=FaceToBeDeleted.begin(); 
      p!=FaceToBeDeleted.end(); p++)
    {
      v.push_back(p->second.hanging());
    }
}

/*********************************************************************/

void HangContainer3d::update_olds(const IntVector& v, const IntVector& c)
{
  HangContainer2d::update_olds(v,c);
  FaceToBeCreated.update(v,c);
  FaceNotAnyMore.update(v,c);      
  FaceHanging.update(v,c);      
}

/*********************************************************************/

void HangContainer3d::update_news(const IntVector& vnew, int offset)
{
  // std::cerr << "new_hangs()" << std::endl;
  // newhangs-hanging setzten fuer new-quad und linehang fuer die zukunft
  int i = offset;
  for (HangList<4>::iterator p = FaceToBeCreated.begin(); p!=FaceToBeCreated.end(); p++)
    {
      p->second.hanging() = vnew[i];
      HangList<4>::iterator Lp = FaceHanging.find(p->first);
      if(Lp!=FaceHanging.end())
	{
	  Lp->second.hanging() = vnew[i]; 
	}
      i++;
    }
  HangContainer2d::update_news(vnew,offset+FaceToBeCreated.size());
}

/*********************************************************************/

int HangContainer3d::vertex_index(const FaceVector& face) const
{
  HangList<4>::const_iterator Np = FaceToBeCreated.find(face);
  if(Np!=FaceToBeCreated.end())
    {
      return Np->second.hanging();
    }
  else
    {
      HangList<4>::const_iterator Tp = FaceNotAnyMore.find(face);
      if(Tp!=FaceNotAnyMore.end())
	{
	  return Tp->second.hanging();
	}
      else
	{
	  HangList<4>::const_iterator Lp = FaceHanging.find(face);
	  if(Lp!=FaceHanging.end())
	    {
	      return Lp->second.hanging();
	    }
	}
    }
  return -1;
}

/*********************************************************************/

void HangContainer3d::face_refine(const FaceVector& face, int f)
{
  std::pair<HangList<4>::iterator,bool> p = find_in_linehang(FaceHanging,face);
  if(p.second)
    {
      // face haengt, aber Nachbar wird nun auch verfeinert

      HangList<4>::iterator Lp = p.first;
      FaceNotAnyMore.insert(*Lp);
      FaceHanging.erase(Lp);
      
      // falls face kreiert werden soll, trage schon mal coarse
      // Nachbarn ein, der streng genommen gar nicht coarser ist.

      HangList<4>::iterator Np = FaceToBeCreated.find(face);
      if (Np!=FaceToBeCreated.end())
	{
	  Np->second.cneighbour() = f;
	}
    }
  else
    {
      // pruefe ob face vorher bestand und zum loeschen
      // vorgesehen ist
      HangList<4>::iterator dp = FaceToBeDeleted.find(face);
      if (dp!=FaceToBeDeleted.end())
	{
	  // wenn ja, dann trage ihn wieder zurueck in die hang Liste ein
	  // und loesche ihn aus der delete Liste

	  Hang& h = dp->second;
	  
	  if (h.cneighbour()==f)
	    {
	      std::swap(h.rneighbour(),h.cneighbour());
	    }
	  FaceHanging.insert(std::make_pair(dp->first,h));
	  FaceToBeDeleted.erase(dp);
	}
      else
	{
	  // face gabs vorher also nicht und trage sie deshalb
	  // in die zu kreierenden faces ein

	  Hang h(-1,f,-1);
	  FaceHanging.insert(std::make_pair(face,h));
	  FaceToBeCreated.insert(std::make_pair(face,h)); 
	}
    }
}

/*********************************************************************/

void HangContainer3d::face_coarse(const FaceVector& face, int f, int face_vertex)
{
  std::pair<HangList<4>::iterator,bool> p = find_in_linehang(FaceHanging,face);
  
  if(p.second)
    {
      // face war vorher haengend, markiere zum loeschen und nicht mehr
      // haengend in Zukunft

      HangList<4>::iterator Lp = p.first;
      Hang& h = Lp->second;

      if (h.cneighbour()==f)
	{
	  std::swap(h.rneighbour(),h.cneighbour());
	}
      FaceToBeDeleted.insert(std::make_pair(Lp->first,h));
      FaceHanging.erase(Lp);
    }
  else
    {
      Hang  h(face_vertex,-1,f);
      FaceHanging.insert(std::make_pair(face,h));
    }
}

/*********************************************************************/

void HangContainer3d::output() const
{
  std::cout << "HangContainer3d: " << std::endl;
  std::cout << "FaceHanging  " << FaceHanging.size() << std::endl;
  //  std::cout << FaceHanging <<std::endl;
  std::cout << "FaceToBeCreated  " << FaceToBeCreated.size() << std::endl;
  //std::cout << FaceToBeCreated << std::endl;
  std::cout << "FaceToBeDeleted  " << FaceToBeDeleted.size() << std::endl;
  //  std::cout << FaceToBeDeleted << std::endl;
  std::cout << "FaceNotMore  " << FaceNotMore().size() << std::endl;
  //  std::cout << FaceNotMore() << std::endl;
  std::cout << "Hanging  " << Hanging.size() << std::endl;
  //  std::cout << Hanging << std::endl;
  std::cout << "VertexToBeCreated  " << VertexToBeCreated.size() << std::endl;
  //  std::cout << VertexToBeCreated << std::endl;
  std::cout << "VertexToBeDeleted  " << VertexToBeDeleted.size() << std::endl;
  //  std::cout << VertexToBeDeleted << std::endl;
  std::cout << "NotAnyMoreHanging  " << NotAnyMoreHanging.size() << std::endl;
  //  std::cout << NotAnyMoreHanging << std::endl;
  //  std::cout << NotAnyMoreHanging << std::endl;
//   HangList<4>::const_iterator Np = FaceDeleting().begin();
//   for (;Np!=FaceDeleting().end(); Np++)
//     {
//       std::cout << Np->first << " * " << Np->second.hanging() << std::endl;
//     }
  std::cout << "-------------" << std::endl;
}

/*********************************************************************/

void HangContainer3d::line_coarse(EdgeVector& edge, int f, int edge_vertex)
{
  HangList<2>::iterator p = Hanging.find(edge);

  if (p==Hanging.end())
    {
      Hang h(edge_vertex,f,-1);
      VertexToBeDeleted.insert(std::make_pair(edge,h));
    }
  else
    {
      p->second.hanging() = edge_vertex;
    }
}

/*********************************************************************/

void HangContainer3d::line_refine(EdgeVector& edge, int f, 
				  const HangList<2>& oldhangs)
{
  if (VertexToBeCreated.find(edge)==VertexToBeCreated.end()) 
    {
      HangList<2>::const_iterator p = oldhangs.find(edge);
      
      if (p==oldhangs.end())
	{
	  Hang h(-1,f,-1);
	  VertexToBeCreated.insert(std::make_pair(edge,h));
	}
    }
}

/*********************************************************************/

// void HangContainer3d::decide_new_line_vertex(const EdgeVector& edge, 
// 					     int f)
// {
//   HangList<2>::const_iterator p = Hanging.find(edge);

//   if (p!=Hanging.end())
//     {
//       int h = p->second.hanging();
//       if (h>=0) return;  // vertex gibts schon
//     }
//   if (VertexToBeCreated.find(edge)!=VertexToBeCreated.end()) 
//     {
//       return; // vertex ist schon markiert zum kreieren
//     }
//   Hang h(-1,-1,-1);
//   VertexToBeCreated.insert(std::make_pair(edge,h));
// }

/*********************************************************************/

