#include  "visualization.h"
#include  "errormacros.h"

/********************************************************************/

void Visualization::vtk(const std::string& bname) const
{
  std::string name = bname;
  name += ".vtk";

  std::ofstream out(name.c_str());
  FILE_ERROR(out,name);

  int nn = mesh->nnodes();

  out << "# vtk DataFile Version 3.1 "<<std::endl;
  out << "output from mmail" << std::endl;
  out << "ASCII" << std::endl;
  out << "DATASET UNSTRUCTURED_GRID" << std::endl;
  out << "POINTS " << nn << " FLOAT " << std::endl;
  
  int ne  = mesh->ncells();
  int npe = mesh->nodes_per_cell();
  int npf = 9;
  if(mesh->dimension()==2)
    { 
      for (int i=0; i<nn; i++)
	{
	  out<<  mesh->vertex2d(i) << " " << 0 << std::endl;
	}
    }
  else
    {
      npf = 12;
      for (int i=0; i<nn; i++)
	{
	  out<<  mesh->vertex3d(i) << std::endl;
	}
    }
  out << std::endl << "CELLS " << ne <<" " << (npe+1)*ne << std::endl;
  for (int c=0; c<ne; c++)
    {
      out << npe << " ";
      for(int i=0;i<npe;i++)
	{
	  out << mesh->vertex_of_cell(c,i) << " "; 
	}
      out << std::endl; 
    }     
  out << std::endl << "CELL_TYPES " << ne << std::endl;
  for(int i=0;i<ne;i++) out << npf << " ";
  out << std::endl << std::endl;

 //  Data

 if (PointData)
   {
     CheckPointData();
     for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
       {
	 if(p==PointDataInfo->sbegin()) out << "POINT_DATA " << nn << std::endl;
	 out << "SCALARS "<< p->first <<" float "<< std::endl;
	 out << "LOOKUP_TABLE default"<< std::endl;
	 for (int ind=0; ind<PointData->visun(); ind++)
	   {
	     if(mesh->dimension()==2)
	       {
		 out << PointData->visudata(ind,p->second,mesh->vertex2d(ind)) << std::endl;
	       }
	     else
	       {
		 out << PointData->visudata(ind,p->second,mesh->vertex3d(ind)) << std::endl;
	       }
	   }
	 out << std::endl<< std::endl;

       }
     for(VisuDataInfo::viterator p=PointDataInfo->vbegin();p!=PointDataInfo->vend();++p)
       {
	 out << "VECTORS "<< p->first <<" float "<< std::endl;
	 for (int ind=0; ind<PointData->visun(); ind++)
	   {
	     for(int ii=0;ii<2;ii++)
	       {
		 if(mesh->dimension()==2)
		   {
		     out << PointData->visudata(ind,p->second[ii],mesh->vertex2d(ind)) << " ";
		   }
		 else
		   {
		     out << PointData->visudata(ind,p->second[ii],mesh->vertex3d(ind)) << std::endl;
		   }
	       }
	     if(p->second[2]==-1)
	       {
		 out << 0. << " ";
	       }
	     else
	       {
		 out << PointData->visudata(ind,p->second[2]) << " ";
	       }
	     out << std::endl;
	   }
	 out << std::endl<< std::endl;
       }
   }
 
 if (CellData)
   {
     CheckCellData();
     for(VisuDataInfo::siterator p=CellDataInfo->sbegin();p!=CellDataInfo->send();++p)
       {
	 if(p==CellDataInfo->sbegin()) out << "CELL_DATA " << mesh->ncells() << std::endl;
	 out << "SCALARS "<< p->first <<" float "<< std::endl;
	 out << "LOOKUP_TABLE default"<< std::endl;

	 for (int ind=0; ind<CellData->visun(); ind++)
	   {
	     out << CellData->visudata(ind,p->second) << std::endl;
	   }
	 out << std::endl<< std::endl;
       }
     for(VisuDataInfo::viterator p=CellDataInfo->vbegin();p!=CellDataInfo->vend();++p)
       {
	 out << "VECTORS "<< p->first <<" float "<< std::endl;
	 for (int ind=0; ind<CellData->visun(); ind++)
	   {
	     for(int ii=0;ii<2;ii++)
	       {
		 out << CellData->visudata(ind,p->second[ii]) << " ";
	       }
	     if(p->second[2]==-1)
	       {
		 out << 0. << " ";
	       }
	     else
	       {
		 out << CellData->visudata(ind,p->second[2]) << " ";
	       }
	     out << std::endl;
	   }
	 out << std::endl<< std::endl;
       }
   }
 out.close();
}
