#include "dg.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"


#include "dgbase.h"

namespace Gascoigne
{
  /////////// INIT, Constructor, ReInit...
  template <class BASE>
  DG<BASE>::DG()
  {
    _integrator.BasicInit();
  }

  template <class BASE>
  void DG<BASE>::BasicInit(const ParamFile *pf){}

  template <class BASE>
  void DG<BASE>::ReInit(const GascoigneMesh *M)
  {
    _dofhandler.InitFromGascoigneMesh(M);
    _mesh = M;
    assert(_mesh);
  }


  ////////// STUFF Data, NodeVEctor, Cellvector ...
  template <class BASE>
  void DG<BASE>::GlobalToGlobalData(LocalParameterData &QP) const
  {
    abort();
    const GlobalParameterData &gpd = GetDataContainer().GetParameterData();

    for (auto p : gpd)
    {
      QP.insert(make_pair(p.first, *p.second));
    }
  }
  
  /* ----------------------------------------- */
  
  template<class BASE>
  void DG<BASE>::GlobalToLocalData(LocalData& QN, const DataContainer& DC,int iq) const
  {
    const auto node = DC.GetNodeData();

    QN.clear();
    
    for(auto it : node)
      {
	_dofhandler.GlobalToLocal(QN[it.first],*it.second,iq);
      }
  }


  ////////////////////////////////////////////////// Init

  template <class BASE>
  void DG<BASE>::Transformation(FemInterface::Matrix &T, int iq) const
  {
    assert(_mesh);
    int dim = _mesh->dimension();
    int ne = _mesh->nodes_per_cell(iq); // always Q1 transformation 
    T.memory(dim, ne);

    IntVector indices = _mesh->IndicesOfCell(iq);

    if (dim == 2)
    {
      for (int ii = 0; ii < ne; ii++)
      {
        Vertex2d v = _mesh->vertex2d(indices[ii]);
        T(0, ii) = v.x();
        T(1, ii) = v.y();
      }
    }
    else if (dim == 3)
    {
      for (int ii = 0; ii < ne; ii++)
      {
        Vertex3d v = _mesh->vertex3d(indices[ii]);
        T(0, ii) = v.x();
        T(1, ii) = v.y();
        T(2, ii) = v.z();
      }
    }
  }

  ////////////////////////////////////////////////// INTEGRATION


  template <class BASE>
  void DG<BASE>::Form(GlobalVector &f,
                      const GlobalVector &u,
                      const Equation &EQ,
                      double d) const
  {
    nmatrix<double> T;
    // LocalParameterData QP;
    // GlobalToGlobalData(QP);
    // EQ.SetParameterData(QP);

    LocalVector U(u.ncomp(), BASE::N);
    LocalVector F(f.ncomp(), BASE::N);

    LocalData QC;

    for (int iq = 0; iq < _dofhandler.nelements(); ++iq)
    {
      Transformation(T, iq);
      _fe.ReInit(T);

      GlobalToLocalData(__QN_master, _datacontainer,iq);
      _dofhandler.GlobalToLocal(U, u, iq);
      _integrator.Form(EQ, F, _fe, U, __QN_master, QC);
      _dofhandler.LocalToGlobal(f, F, iq, d);
    }

    EdgeForm(f, u, dynamic_cast<const DGEquation &>(EQ), d);

  }

  template <class BASE>
  void DG<BASE>::EdgeForm(GlobalVector &f,
                          const GlobalVector &u,
                          const DGEquation &EQ,
                          double d) const
  {
    nmatrix<double> T1, T2;
    // LocalParameterData QP;
    // GlobalToGlobalData(QP);
    // EQ.SetParameterData(QP);

    LocalVector U1(u.ncomp(), BASE::N);
    LocalVector U2(u.ncomp(), BASE::N);
    LocalVector F1(f.ncomp(), BASE::N);
    LocalVector F2(f.ncomp(), BASE::N);

    LocalData QC;

    for (int ie = 0; ie < _dofhandler.nedges(); ++ie)
    {
      auto &E = _dofhandler.getedge(ie);

      // master
      size_t master = E[0];
      size_t masterli = E[2];
      assert(master != -1);
      Transformation(T1, master);
      _fe.ReInit(T1);
      _dofhandler.GlobalToLocal(U1, u, master);
      GlobalToLocalData(__QN_master, _datacontainer,master);


      // slave
      size_t slave = E[1];
      size_t slaveli = E[3];
      bool internaledge = (slave != -1);
      if (internaledge)
      {
        Transformation(T2, slave);
        _feslave.ReInit(T2);
        _dofhandler.GlobalToLocal(U2, u, slave);
	GlobalToLocalData(__QN_slave, _datacontainer,slave);
      }
      
      _integrator.EdgeForm(internaledge,
                           EQ,
                           F1,
                           F2,
                           _fe,
                           _feslave,
                           masterli,
                           slaveli,
                           U1,
                           U2,
                           __QN_master,
			   __QN_slave,
                           QC);
      _dofhandler.LocalToGlobal(f, F1, master, d);
      if (internaledge)
        _dofhandler.LocalToGlobal(f, F2, slave, d);
    }
  }


  template <class BASE>
  void
  DG<BASE>::Rhs(GlobalVector &f, const DomainRightHandSide &RHS, double s) const
  {
    nmatrix<double> T;

    // LocalParameterData QP;
    // GlobalToGlobalData(QP);
    // RHS.SetParameterData(QP);

    LocalVector F(f.ncomp(), BASE::N);

    LocalData QC;

    for (int iq = 0; iq < _dofhandler.nelements(); ++iq)
    {
      Transformation(T, iq);
      _fe.ReInit(T);
      GlobalToLocalData(__QN_master, _datacontainer,iq);
      _integrator.Rhs(RHS, F, _fe, __QN_master, QC);
      _dofhandler.LocalToGlobal(f, F, iq, s);
    }
  }

  template <class BASE>
  void DG<BASE>::Matrix(MatrixInterface &A,
                        const GlobalVector &u,
                        const Equation &EQ,
                        double d) const
  {
    nmatrix<double> T;

    // LocalParameterData QP;
    // GlobalToGlobalData(QP);
    // EQ.SetParameterData(QP);

    LocalData QC;

    LocalVector U(u.ncomp(), BASE::N);

    EntryMatrix E;


    for (int iq = 0; iq < _dofhandler.nelements(); ++iq)
    {
      Transformation(T, iq);
      _fe.ReInit(T);

      _dofhandler.GlobalToLocal(U, u, iq);
      GlobalToLocalData(__QN_master, _datacontainer,iq);

      // EQ.cell(GetMesh(),iq,__U,__QN);
      _integrator.Matrix(EQ, E, _fe, U, __QN_master, QC);
      _dofhandler.LocalToGlobalMatrix(A, E, iq, d);
    }
    EdgeMatrix(A, u, dynamic_cast<const DGEquation &>(EQ), d);

  }

  template <class BASE>
  void DG<BASE>::EdgeMatrix(MatrixInterface &A,
                            const GlobalVector &u,
                            const DGEquation &EQ,
                            double d) const
  {
    nmatrix<double> T1, T2;
    // LocalParameterData QP;
    // GlobalToGlobalData(QP);
    // EQ.SetParameterData(QP);

    LocalVector U1(u.ncomp(), BASE::N);
    LocalVector U2(u.ncomp(), BASE::N);
    EntryMatrix E11, E12, E21, E22;

    LocalData QC;

    for (int ie = 0; ie < _dofhandler.nedges(); ++ie)
    {
      auto &E = _dofhandler.getedge(ie);

      // master
      size_t master = E[0];
      size_t masterli = E[2];
      assert(master != -1);
      Transformation(T1, master);
      _fe.ReInit(T1);
      _dofhandler.GlobalToLocal(U1, u, master);
      GlobalToLocalData(__QN_master, _datacontainer,master);

      // slave
      size_t slave = E[1];
      size_t slaveli = E[3];
      bool internaledge = (slave != -1);
      if (internaledge)
      {
        Transformation(T2, slave);
        _feslave.ReInit(T2);
        _dofhandler.GlobalToLocal(U2, u, slave);
	GlobalToLocalData(__QN_slave, _datacontainer,slave);
      }

      //_dofhandler.GlobalToLocalData(iq);
      _integrator.EdgeMatrix(internaledge,
                             EQ,
                             E11,
                             E12,
                             E21,
                             E22,
                             _fe,
                             _feslave,
                             masterli,
                             slaveli,
                             U1,
                             U2,
                             __QN_master,
			     __QN_slave,
                             QC);
      _dofhandler.LocalToGlobalMatrix(A, E11, master,master, d);
      if (internaledge)
      {
	_dofhandler.LocalToGlobalMatrix(A, E12, master, slave, d);
	_dofhandler.LocalToGlobalMatrix(A, E21, slave, master, d);
        _dofhandler.LocalToGlobalMatrix(A, E22, slave, slave,d);
      }
    }
  }


  template class DG<BASEQ12D>;
  template class DG<BASEQ22D>;
}

#include "finiteelement.xx"
namespace Gascoigne
{
  template class Transformation2d<BASEQ12D>;
  template class FiniteElement<2, 1, Transformation2d<BASEQ12D>, BASEQ12D>;
}
