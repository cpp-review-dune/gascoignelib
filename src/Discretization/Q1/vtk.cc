/**
 *
 * Copyright (C) 2004, 2005, 2010, 2011 by the Gascoigne 3D authors
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

#include "errormacros.h"
#include "stlio.h"
#include "visualization.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne {
void
Visualization::_vtk_pointdata(ofstream& out) const
{
  if (PointData) {
    int nn = mesh->nnodes();

    CheckPointData();
    // for(VisuDataInfo::siterator
    // p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
    for (int i = 0; i < PointDataInfo->nscalars(); i++) {
      VisuDataInfo::siterator p =
        (const_cast<VisuDataInfo*>(PointDataInfo))->GetSIterator(i);
      if (i == 0)
        out << "POINT_DATA " << nn << '\n';
      out << "SCALARS " << p->first << " DOUBLE " << '\n';
      out << "LOOKUP_TABLE default" << '\n';
      for (int ind = 0; ind < PointData->visun(); ind++) {
        if (mesh->dimension() == 2) {
          assert(ind < mesh->nnodes());
          out << PointData->visudata2(ind, p->second, mesh->vertex2d(ind))
              << '\n';
        } else {
          out << PointData->visudata2(ind, p->second, mesh->vertex3d(ind))
              << '\n';
        }
      }
      out << '\n' << '\n';
    }
    out << flush;
    // for(VisuDataInfo::viterator
    // p=PointDataInfo->vbegin();p!=PointDataInfo->vend();++p) {
    for (int i = 0; i < PointDataInfo->nvectors(); i++) {
      VisuDataInfo::viterator p =
        (const_cast<VisuDataInfo*>(PointDataInfo))->GetVIterator(i);
      out << "VECTORS " << p->first << " DOUBLE " << '\n';
      for (int ind = 0; ind < PointData->visun(); ind++) {
        for (int ii = 0; ii < 2; ii++) {
          if (mesh->dimension() == 2) {
            out << PointData->visudata2(ind, p->second[ii], mesh->vertex2d(ind))
                << ' ';
          } else {
            out << PointData->visudata2(ind, p->second[ii], mesh->vertex3d(ind))
                << '\n';
          }
        }
        if (p->second[2] == -1) {
          out << 0. << ' ';
        } else {
          out << PointData->visudata(ind, p->second[2]) << ' ';
        }
        out << '\n';
      }
      out << '\n' << '\n';
    }
    out << flush;
  }
}

/* ----------------------------------------- */

void
Visualization::_vtk_celldata(ofstream& out) const
{
  if (CellData) {
    CheckCellData();

    // cout << "CellDataInfo->nscalars()" << CellDataInfo->nscalars() << '\n';
    //
    // die reihenfolge der elemente per index ist wesentlich, es reicht nicht
    // nur sie per iterator aus CellDataInfo raus zu holen
    // for(VisuDataInfo::siterator
    // p=CellDataInfo->sbegin();p!=CellDataInfo->send();++p){
    for (int i = 0; i < CellDataInfo->nscalars(); i++) {
      VisuDataInfo::siterator p =
        (const_cast<VisuDataInfo*>(CellDataInfo))->GetSIterator(i);
      if (i == 0)
        out << "CELL_DATA " << mesh->ncells() << '\n';
      // if(p==CellDataInfo->sbegin()) out << "CELL_DATA " << mesh->ncells() <<
      // '\n';
      out << "SCALARS " << p->first << " DOUBLE " << '\n';
      out << "LOOKUP_TABLE default" << '\n';

      for (int ind = 0; ind < CellData->visun(); ind++) {
        out << CellData->visudata(ind, p->second) << '\n';
      }
      out << '\n' << '\n';
    }
    out << flush;
    // cout << "CellDataInfo->nvectors()" << CellDataInfo->nvectors() << '\n';
    // for(VisuDataInfo::viterator
    // p=CellDataInfo->vbegin();p!=CellDataInfo->vend();++p){
    for (int i = 0; i < CellDataInfo->nvectors(); i++) {
      VisuDataInfo::viterator p =
        (const_cast<VisuDataInfo*>(CellDataInfo))->GetVIterator(i);
      out << "VECTORS " << p->first << " DOUBLE " << '\n';
      for (int ind = 0; ind < CellData->visun(); ind++) {
        for (int ii = 0; ii < 2; ii++) {
          out << CellData->visudata(ind, p->second[ii]) << ' ';
        }
        if (p->second[2] == -1) {
          out << 0. << ' ';
        } else {
          out << CellData->visudata(ind, p->second[2]) << ' ';
        }
        out << '\n';
      }
      out << '\n' << '\n';
    }
    out << flush;
  }
}

/* ----------------------------------------- */

void
Visualization::_vtk_points(ofstream& out) const
{
  int nn = mesh->nnodes();
  out << "POINTS " << nn << " DOUBLE " << '\n';
  assert(mesh->dimension() == 2 || mesh->dimension() == 3);
  if (mesh->dimension() == 2) {
    for (int i = 0; i < nn; i++) {
      out << mesh->vertex2d(i) << ' ' << 0 << '\n';
    }
  } else if (mesh->dimension() == 3) {
    for (int i = 0; i < nn; i++) {
      out << mesh->vertex3d(i) << '\n';
    }
  }
  out << endl;
}

/* ----------------------------------------- */

void
Visualization::_vtk_cells(ofstream& out) const
{
  int ne = mesh->ncells();

  int lenght = 0;
  for (int c = 0; c < ne; c++) {
    lenght += mesh->nodes_per_cell(c) + 1;
  }

  out << '\n' << "CELLS " << ne << ' ' << lenght << '\n';

  for (int c = 0; c < ne; c++) {
    int nle = mesh->nodes_per_cell(c);
    out << nle << ' ';
    for (int ii = 0; ii < nle; ii++) {
      out << mesh->vertex_of_cell(c, ii) << ' ';
    }
    out << '\n';
  }
  out << '\n' << "CELL_TYPES " << ne << '\n';
  for (int c = 0; c < ne; c++) {
    out << mesh->VtkType(c) << ' ';
  }
  out << endl;
}

void
Visualization::_vtk_cellmaterial(ofstream& out) const
{
  int ne = mesh->ncells();

  out << '\n'
      << "CELL_DATA " << ne << '\n'
      << "FIELD FieldData 1" << '\n'
      << "material 1 " << ne << " int" << '\n';

  for (int c = 0; c < ne; c++)
    out << mesh->material(c) << ' ';
  out << endl;
}

/* ----------------------------------------- */

void
Visualization::vtk(const string& bname) const
{
  string name = bname;
  name += ".vtk";

  ofstream out(name.c_str());
  FILE_ERROR(out, name);

  //  Header

  out << "# vtk DataFile Version 2.0 " << '\n';
  out << "output from GascoigneStd, " << title << '\n';
  out << "ASCII" << '\n';
  out << "DATASET UNSTRUCTURED_GRID" << '\n';
  out << "FIELD FieldData 1" << '\n';
  out << "TIME 1 1 double" << '\n';
  out << time << endl;

  //  Mesh
  _vtk_points(out);
  _vtk_cells(out);
  // material
  if (cellmaterial)
    _vtk_cellmaterial(out);
  //  Data
  _vtk_pointdata(out);
  _vtk_celldata(out);

  out.close();
  if (compress) {
    string command = "gzip -f " + name;
    int status = system(command.c_str());
    (void)status; // to avoid warning;
    assert(status == 0);
  }
}
} // namespace Gascoigne
