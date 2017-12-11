/**
*
* Copyright (C) 2004, 2005, 2006, 2008, 2010, 2011 by the Gascoigne 3D authors
*
* This file is part of Gascoigne 3D
*
* Gascoigne 3D is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version.
*
* Gascoigne 3D is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* Please refer to the file LICENSE.TXT for further information
* on this license.
*
**/


#include  "finehierarchicalmesh3d.h"
#include  "deletecells.h"
#include  "vecalgo.h"
#include  "set2vec.h"
#include  "quad.h"

#include  "facemanager.h"
#include  "stlio.h"
#include  "regular_update.h"
#include  "coarsehierarchicalmesh3d.h"
#include  "filescanner.h"
#include  <fstream>
#include  "giota.h"

#include  "stringutil.h"
using namespace std;

namespace Gascoigne
{
typedef  triple<int,int,int>          tint;

/*------------------------------------------------------*/


/*------------------------------------------------------*/

FineHierarchicalMesh3d::FineHierarchicalMesh3d (const std::string gridname)
  : HierarchicalMesh3d(), HexLaO(hexs)
{

  if(gridname=="none")
    {
      cerr << "no \"gridname\" " << endl;
      return;
      abort();
    }
  vector<string> s = StringSplit(gridname.c_str(),'.');
  string suff = s[s.size()-1];
  if(suff=="inp")
    {
      read_inp(gridname);
    }
  else
    {
      cerr << "HierarchicalMesh::read():\tunknown suffix " << suff << endl;
      abort();
    }
}

/*------------------------------------------------------*/
/*
FineHierarchicalMesh3d::FineHierarchicalMesh3d(const FineHierarchicalMesh3d& H)
  : HierarchicalMesh(), HexLaO(hexs)
{
  *this = H;
}
*/
/*------------------------------------------------------*/
/*
FineHierarchicalMesh3d& FineHierarchicalMesh3d::operator=(const FineHierarchicalMesh3d& H)
{
  HierarchicalMesh::operator= (H);
  // copy all data
  vertexs3d = H.vertex3d();
  hexs      = H.hex();
  Bquads    = H.bquad();
  return *this;
}
*/




/*------------------------------------------------------*/






/*---------------------------------------------------*/

pair<bool,tint> FineHierarchicalMesh3d::check_inp(const string& name)
{
  //  detect some errors in input-file... and more

  ifstream file(name.c_str());
  if(!file.is_open())
    {
      cerr << "HierarchicalMesh3d::check_inp()\n";
      cerr << "cannot open file " << name << endl;
      abort();
    }

  bool first_one = 1;
  int  nv, nl, nq, nh, nt;
  int  n_unkonwn;
  file >> nv >> nt >> n_unkonwn >> n_unkonwn >> n_unkonwn;

  Vertex3d c; int ind;
  for(int i=0;i<nv;i++)
    {
      file >> ind >> c;
    }

  nh = 0; nq = 0; nl = 0;
  fixarray<8,int> ih;
  fixarray<4,int> iq;
  fixarray<2,int> il;
  for(int i=0;i<nt;i++)
    {
      string name;
      string mat;
      int ii;
      file >> ii >> mat  >> name;
      if(name=="hex")
	{
	  file >> ih;
	  nh++;
	  if( (ih[0]==0)||(ih[1]==0)||(ih[2]==0)||(ih[3]==0) ||
	      (ih[4]==0)||(ih[5]==0)||(ih[6]==0)||(ih[7]==0))
	    {
	      first_one = 0;
	    }
	}
      else if(name=="quad")
	{
	  file >> iq;
	  nq++;
	  if( (iq[0]==0)||(iq[1]==0)||(iq[2]==0)||(iq[3]==0) )
	    {
	      first_one = 0;
	    }
	}
      else if(name=="line")
	{
	  file >> il;
	  nl++;
	  if( (il[0]==0)||(il[1]==0))
	    {
	      first_one = 0;
	    }
	}
    }

  // fehlerabfragen ....
  if(nt!=(nl+nq+nh))
    {
      cerr << "wrong number of cells: " << nt << endl;
      cerr << "lines quads hexs: " << nl << " " << nq << " " << nh << endl;
      abort();
    }

  file.close();

  return make_pair(first_one,make_triple(nl,nq,nh));
}

/*---------------------------------------------------*/

void FineHierarchicalMesh3d::read_inp(const string& name)
{
  cout<<"read_inpcalled "<<endl;
  // check mesh.....
  pair<bool,tint> p = check_inp(name);
  bool first_one = p.first;
  tint n = p.second;

  int nl = n.first;
  int nq = n.second;
  int nh = n.third;

  ifstream file(name.c_str());
  if(!file.is_open())
    {
      cerr << "HierarchicalMesh3d::read_inp()\n";
      cerr << "cannot open file " << name << endl;
      abort();
    }

  int  nv,  nt, n_unkonwn;

  file >> nv >> nt >> n_unkonwn >> n_unkonwn >> n_unkonwn;

  cout << "3D Mesh:  " << nv << " nodes, ";
  cout << nl << " lines, ";
  cout << nq << " quads, ";
  cout << nh << " hexs" << endl;

  vertexs3d.reserve(nv);
  vertexs3d.resize(nv,Vertex3d());
  //hexs.     reserve(nh);
  //hexs.     resize(nh,Hex());
  Bquads.   reserve(nq);
  Bquads.   resize(nq);

  Vertex3d c;      int ind;
  for(int i=0;i<nv;i++)
    {
      file >> ind >> c;
      vertexs3d[i] = c;
    }
  fixarray<8,int> ihv;
  fixarray<4,int> iqv;
  fixarray<2,int> ilv;
  int ih = 0;
  int iq = 0;
  cout<<"read_inp nt of for:  "<<nt<<endl;
  for(int i=0;i<nt;i++)
    {
      string name;	int unknown; string matstring;
      file >> unknown >> matstring  >> name;
      if(name=="hex")
	{
	  file >> ihv;
	  //if(first_one) for(int iii=0;iii<8;iii++) ihv[iii]--;
	  //hexs[ih].vertex() = ihv;
	  //hexs[ih].material() = atoi(matstring.c_str());
	  //ih++;
	}
      else if(name=="quad")
	{
	  file >> iqv;
	  if(first_one) for(int iii=0;iii<4;iii++) iqv[iii]--;
	  
	  BoundaryQuad li;
	  li.material() = atoi(matstring.c_str());
	  li.vertex() = iqv;
	  //init_quad(li);
	  Bquads[iq++]=li;
	}
    }
  //init_edges3d();

}
/*---------------------------------------------------*/



/*---------------------------------------------------*/








}

/*---------------------------------------------------*/

