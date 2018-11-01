
/**
 *
 * Copyright (C) 2018 by the Gascoigne 3D authors
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

#include "hangingnodes.h"
/*-----------------------------------------*/

namespace Gascoigne
{
  template <>
  HangingNodes<2, 2>::HangingNodes()
  {
    wei[0] = 0.5;
    wei[1] = 0.5;

    lnoe[0][0] = 0;
    lnoe[0][1] = 1;
    lnoe[0][2] = 3;
    lnoe[1][0] = 1;
    lnoe[1][1] = 3;
    lnoe[1][2] = 2;
    lnoe[2][0] = 3;
    lnoe[2][1] = 2;
    lnoe[2][2] = 0;
    lnoe[3][0] = 2;
    lnoe[3][1] = 0;
    lnoe[3][2] = 1;

    lnop[0][0] = 0;
    lnop[0][1] = 2;
    lnop[0][2] = 1;
    lnop[1][0] = 2;
    lnop[1][1] = 8;
    lnop[1][2] = 5;
    lnop[2][0] = 8;
    lnop[2][1] = 6;
    lnop[2][2] = 7;
    lnop[3][0] = 6;
    lnop[3][1] = 0;
    lnop[3][2] = 3;
  }
  template <>
  HangingNodes<3, 2>::HangingNodes()
  {
    wei[0] = 0.375;
    wei[1] = 0.75;
    wei[2] = -0.125;

    lnoe[0][0] = 0;
    lnoe[0][1] = 1;
    lnoe[0][2] = 2;
    lnoe[1][0] = 0;
    lnoe[1][1] = 3;
    lnoe[1][2] = 6;
    lnoe[2][0] = 2;
    lnoe[2][1] = 5;
    lnoe[2][2] = 8;
    lnoe[3][0] = 6;
    lnoe[3][1] = 7;
    lnoe[3][2] = 8;
  }
  template <>
  HangingNodes<2, 3>::HangingNodes()
  {
    wei[0] = 0.5;
    wei[1] = 0.5;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        fwei[3 * i + j] = wei[i] * wei[j];

    lnoe[0][0] = 0;
    lnoe[0][1] = 2;
    lnoe[0][2] = 1;
    lnoe[1][0] = 0;
    lnoe[1][1] = 6;
    lnoe[1][2] = 3;
    lnoe[2][0] = 2;
    lnoe[2][1] = 8;
    lnoe[2][2] = 5;
    lnoe[3][0] = 6;
    lnoe[3][1] = 8;
    lnoe[3][2] = 7;

    lnop[0][0] = 0;
    lnop[0][1] = 2;
    lnop[0][2] = 6;
    lnop[0][3] = 8;
    lnop[0][4] = 4;
    lnop[1][0] = 2;
    lnop[1][1] = 20;
    lnop[1][2] = 8;
    lnop[1][3] = 26;
    lnop[1][4] = 14;
    lnop[2][0] = 6;
    lnop[2][1] = 8;
    lnop[2][2] = 24;
    lnop[2][3] = 26;
    lnop[2][4] = 16;
    lnop[3][0] = 0;
    lnop[3][1] = 18;
    lnop[3][2] = 6;
    lnop[3][3] = 24;
    lnop[3][4] = 12;
    lnop[4][0] = 0;
    lnop[4][1] = 2;
    lnop[4][2] = 18;
    lnop[4][3] = 20;
    lnop[4][4] = 10;
    lnop[5][0] = 18;
    lnop[5][1] = 20;
    lnop[5][2] = 24;
    lnop[5][3] = 26;
    lnop[5][4] = 22;
  }
  template <>
  HangingNodes<3, 3>::HangingNodes()
  {
    wei[0] = 0.375;
    wei[1] = 0.75;
    wei[2] = -0.125;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        fwei[3 * i + j] = wei[i] * wei[j];

    lnoe[0][0] = 0;
    lnoe[0][1] = 1;
    lnoe[0][2] = 2;
    lnoe[1][0] = 0;
    lnoe[1][1] = 3;
    lnoe[1][2] = 6;
    lnoe[2][0] = 2;
    lnoe[2][1] = 5;
    lnoe[2][2] = 8;
    lnoe[3][0] = 6;
    lnoe[3][1] = 7;
    lnoe[3][2] = 8;

    lnoe[4][0] = 18;
    lnoe[4][1] = 19;
    lnoe[4][2] = 20;
    lnoe[5][0] = 18;
    lnoe[5][1] = 21;
    lnoe[5][2] = 24;
    lnoe[6][0] = 20;
    lnoe[6][1] = 23;
    lnoe[6][2] = 26;
    lnoe[7][0] = 24;
    lnoe[7][1] = 25;
    lnoe[7][2] = 26;

    lnoe[8][0] = 0;
    lnoe[8][1] = 9;
    lnoe[8][2] = 18;
    lnoe[9][0] = 2;
    lnoe[9][1] = 11;
    lnoe[9][2] = 20;
    lnoe[10][0] = 6;
    lnoe[10][1] = 15;
    lnoe[10][2] = 24;
    lnoe[11][0] = 8;
    lnoe[11][1] = 17;
    lnoe[11][2] = 26;

    lnop[0][0] = 0;
    lnop[0][1] = 2;
    lnop[0][2] = 6;
    lnop[0][3] = 8;
    lnop[0][4] = 4;
    lnop[1][0] = 2;
    lnop[1][1] = 20;
    lnop[1][2] = 8;
    lnop[1][3] = 26;
    lnop[1][4] = 14;
    lnop[2][0] = 6;
    lnop[2][1] = 8;
    lnop[2][2] = 24;
    lnop[2][3] = 26;
    lnop[2][4] = 16;
    lnop[3][0] = 0;
    lnop[3][1] = 18;
    lnop[3][2] = 6;
    lnop[3][3] = 24;
    lnop[3][4] = 12;
    lnop[4][0] = 0;
    lnop[4][1] = 2;
    lnop[4][2] = 18;
    lnop[4][3] = 20;
    lnop[4][4] = 10;
    lnop[5][0] = 18;
    lnop[5][1] = 20;
    lnop[5][2] = 24;
    lnop[5][3] = 26;
    lnop[5][4] = 22;
  }


  template <>
  void HangingNodes<2, 2>::CondenseHanging(EntryMatrix &E,
                                           IntVector &indices) const
  {
    assert(indices.size() == lnoe.size());
    for (int ii = 0; ii < indices.size(); ii++)
    {
      std::array<int, 3> p = lnoe[ii];

      int &hang = indices[p[1]];

      if (!hanging(hang))
        continue;

      const std::array<int, 3> &f = regular_nodes(hang);

      if ((indices[p[2]] == f[0]) || (indices[p[2]] == f[1]))
        std::swap(p[0], p[2]);

      if (indices[p[0]] == f[0])
        hang = f[1];
      else if (indices[p[0]] == f[1])
        hang = f[0];
      else
        assert(0);

      E.add_column(p[0], p[1], wei[0]);
      E.multiply_column(p[1], wei[1]);

      E.add_row(p[0], p[1], wei[0]);
      E.multiply_row(p[1], wei[1]);
    }
  }

  template <>
  void HangingNodes<2, 3>::CondenseHanging(EntryMatrix &E,
                                           IntVector &indices) const
  {
    for (int ii = 0; ii < 4; ii++) // nur 4 kandiaten koennen haengen !!
    {
      int i = indices[2 * ii + 1];
      if (!hanging(i))
        continue;

      const std::array<int, 3> &f = regular_nodes(i);

      std::array<int, 3> p = lnoe[ii];

      if ((indices[p[0]] == f[1]) && (indices[p[1]] == f[0]))
      {
        std::swap(p[0], p[1]);
      }

      indices[p[2]] = f[2];

      E.add_column(p[0], p[2], wei[0]);
      E.add_column(p[1], p[2], wei[1]);
      E.multiply_column(p[2], wei[2]);

      E.add_row(p[0], p[2], wei[0]);
      E.add_row(p[1], p[2], wei[1]);
      E.multiply_row(p[2], wei[2]);
    }
  }

  template <>
  void HangingNodes<3, 2>::CondenseHanging(EntryMatrix &E,
                                           IntVector &indices) const
  {
    if (1)
    {
      IntVector x(0), y(0);

      for (int ii = 0; ii < 8; ii++)
      {
        int i = indices[ii];

        if (hanging(i) == 2) // 2er haengender Knoten
        {
          std::array<int, 2> Edge = GetHangingEdge(i);

          x.push_back(ii);

          for (int iii = 0; iii < 2; iii++)
          {
            int ir = Edge[iii];
            bool found = 0;
            for (int iiii = 0; (iiii < 8) && !found; iiii++)
            {
              if (ir == indices[iiii])
              {
                found = 1;
                y.push_back(iiii);
              }
            }
            if (!found)
              indices[ii] = ir;
          }
        }
      }
      assert(x.size() == y.size());
      assert(x.size() <= 3);

      for (int i = 0; i < x.size(); i++)
      {
        int i1 = x[i]; // new node !
        int i2 = y[i]; // already there

        E.multiply_column_row(i1, 0.5);
        E.add_column_row(i2, i1);
      }
    }
    if (1)
    {
      IntVector x(0), y(0);

      for (int ii = 0; ii < 8; ii++)
      {
        int j = indices[ii];

        if (hanging(j) == 4) // 4er haengender Knoten
        {
          std::array<int, 4> Face = GetHangingFace(j);

          x.push_back(ii);
          for (int i = 0; i < 4; i++)
          {
            int FaceIndex = Face[i];
            //
            // suche ob FaceIndex schon in indices sind
            //
            int jj = 0;
            bool found = 0;
            while ((jj < 8) && !found)
            {
              found = (indices[jj] == FaceIndex);
              jj++;
            }
            jj--;
            if (found)
              y.push_back(jj); // merke Kopplung in Nachbar vertex
            else
              indices[ii] = FaceIndex; // ersetze Kopplung
          }
        }
      }
      assert(y.size() == 3 * x.size());

      int counter = 0;
      for (int i = 0; i < x.size(); i++)
      {
        int i1 = x[i]; // new node !

        E.multiply_column_row(i1, 0.25);
        int last = counter + 3;

        assert(last <= y.size());

        for (; counter < last; counter++)
        {
          int i2 = y[counter]; // already there
          E.add_column_row(i2, i1);
        }
      }
    }
  }


  template <>
  void HangingNodes<3, 3>::CondenseHanging(EntryMatrix &E,
                                           IntVector &indices) const
  {
    if (1)
    {

      IntVector x(0), y(0);

      for (int ii = 0; ii < 8; ii++)
      {
        int i = indices[ii];

        if (hanging(i) == 2) // 2er haengender Knoten
        {
          auto Edge = GetHangingEdge(i);

          x.push_back(ii);

          for (int iii = 0; iii < 2; iii++)
          {
            int ir = Edge[iii];
            bool found = 0;
            for (int iiii = 0; (iiii < 8) && !found; iiii++)
            {
              if (ir == indices[iiii])
              {
                found = 1;
                y.push_back(iiii);
              }
            }
            if (!found)
              indices[ii] = ir;
          }
        }
      }
      assert(x.size() == y.size());
      assert(x.size() <= 3);

      for (int i = 0; i < x.size(); i++)
      {
        int i1 = x[i]; // new node !
        int i2 = y[i]; // already there

        E.multiply_column_row(i1, 0.5);
        E.add_column_row(i2, i1);
      }
    }

    if (1)
    {

      for (int i = 0; i < 12; i++)
      {
        std::array<int, 3> p = lnoe[i];

        int elim = p[1];
        int h = indices[elim];

        auto q = edges->find(h);

        if (q == edges->end())
          continue;

        const std::array<int, 3> &f = q->second;

        indices[elim] = f[2];

        if ((indices[p[0]] == f[1]) && (indices[p[2]] == f[0]))
        {
          std::swap(p[0], p[2]);
        }
        assert(indices[p[0]] == f[0]);
        assert(indices[p[2]] == f[1]);

        E.add_column(p[0], elim, wei[0]);
        E.add_column(p[2], elim, wei[1]);
        E.multiply_column(elim, wei[2]);

        E.add_row(p[0], elim, wei[0]);
        E.add_row(p[2], elim, wei[1]);
        E.multiply_row(elim, wei[2]);
      }

      for (int i = 0; i < 6; i++)
      {
        std::array<int, 5> lf = lnop[i];

        int elim = lf[4];
        int h = indices[elim];

        auto q = faces->find(h);

        if (q == faces->end())
          continue;

        const std::array<int, 9> &gf = q->second;

        indices[elim] = gf[8];

        std::array<int, 8> x;
        x.fill(-1);

        for (int j = 0; j < indices.size(); j++)
        {
          int k = indices[j];
          if (k == gf[2])
            x[2] = j;
          else if (k == gf[5])
            x[5] = j;
          else if (k == gf[6])
            x[6] = j;
          else if (k == gf[7])
            x[7] = j;
        }
        for (int j = 0; j < 4; j++)
        {
          int k = indices[lf[j]];
          if (k == gf[0])
            x[0] = lf[j];
          else if (k == gf[1])
            x[1] = lf[j];
          else if (k == gf[3])
            x[3] = lf[j];
          else if (k == gf[4])
            x[4] = lf[j];
          else
            assert(0);
        }
        for (int j = 0; j < 8; j++)
        {
          assert(x[j] >= 0);
          E.add_column(x[j], elim, fwei[j]);
          E.add_row(x[j], elim, fwei[j]);
        }
        E.multiply_column(elim, fwei[8]);
        E.multiply_row(elim, fwei[8]);
      }
    }
  }

}


template class Gascoigne::HangingNodes<2, 2>;
template class Gascoigne::HangingNodes<2, 3>;
template class Gascoigne::HangingNodes<3, 2>;
template class Gascoigne::HangingNodes<3, 3>;
