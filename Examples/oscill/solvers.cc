#include "solvers.h"

#include "alediscretization.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
//#include "splittingilu.h"
#include "colors.h"
#include "umfilu.h"
#include "fmatrixblock.h"
#include "cuthillmckee.h"
#include "mginterpolatornested.h"


using namespace std;

namespace Gascoigne
{
/*-------------------------------------------------------*/

template <int DIM>
void FSISolver<DIM>::DeleteSolidPressure(VectorInterface& gf) const
{
    ////////////////////
    if (GetProblemDescriptor()->GetName() == "div")
    {
        DeleteSolidPressure_DIV(gf);
    }
    else
    {
        const HASHSET<int>& i_nodes = GetAleDiscretization()->GetInterfaceNodes();
        const vector<int>& s_nodes  = GetAleDiscretization()->GetSolidL2G();
        vector<int> cv;
        for (int i = 0; i < s_nodes.size(); ++i)
            if (i_nodes.find(s_nodes[i]) == i_nodes.end())
            {
                GetGV(gf)(s_nodes[i], 0) = 0.0;
            }
    }
}

template <int DIM>
void FSISolver<DIM>::SetBoundaryVectorZero(VectorInterface& gf) const
{
    StdSolver::SetBoundaryVectorZero(gf);

    DeleteSolidPressure(gf);
}

template <int DIM>
void FSISolver<DIM>::SetBoundaryVector(VectorInterface& gf) const
{
    StdSolver::SetBoundaryVector(gf);

    DeleteSolidPressure(gf);
}

template <int DIM>
void FSISolver<DIM>::AssembleMatrix(const VectorInterface& gu, double d)
{
    if (GetProblemDescriptor()->GetName() == "div")
    {
        AssembleMatrix_DIV(gu, d);
    }
    else
    {
        StdSolver::AssembleMatrix(gu, d);
        // cout << GetProblemDescriptor()->GetName()
        if ((_directsolver) || (_matrixtype == "block"))
        {
            // Modify for pressure zero in Solid-Part
            const HASHSET<int>& i_nodes = GetAleDiscretization()->GetInterfaceNodes();
            const vector<int>& s_nodes  = GetAleDiscretization()->GetSolidL2G();
            const vector<int>& f_nodes  = GetAleDiscretization()->GetFluidL2G();
            vector<int> cv;
            cv.push_back(0);
            for (int i = 0; i < s_nodes.size(); ++i)
                if (i_nodes.find(s_nodes[i]) == i_nodes.end())
                    GetMatrix()->dirichlet(s_nodes[i], cv);
        }
    }
}

template <int DIM>
void FSISolver<DIM>::DeleteSolidPressure_DIV(VectorInterface& gf) const
{
    ////////////////////
    VectorInterface old_iface("old");
    GlobalVector& old_gv         = StdSolver::GetGV(old_iface);
    // const HASHSET<int>& i_nodes = GetAleDiscretization()->GetInterfaceNodes();
    const vector<int>& s_nodes  = GetAleDiscretization()->GetSolidL2G();
    vector<int> cv;

    // u=uold everywhere
    for (auto i = 0; i < old_gv.n(); ++i)
    {
        GetGV(gf)(i, 0) = old_gv(i, 0);
    }

    //  v = vold in structure and interface
    for (auto i = 0; i < s_nodes.size(); ++i)
    {
        GetGV(gf)(s_nodes[i], 1) = old_gv(s_nodes[i], 1);
        GetGV(gf)(s_nodes[i], 2) = old_gv(s_nodes[i], 2);
    }
}

template <int DIM>
void FSISolver<DIM>::AssembleMatrix_DIV(const VectorInterface& gu, double d)
{
    StdSolver::AssembleMatrix(gu, d);
    cout << style::bb << "Solving " << GetProblemDescriptor()->GetName();
    if ((_directsolver) || (_matrixtype == "block"))
    {
        // Modify for pressure zero in Solid-Part
        const HASHSET<int>& i_nodes = GetAleDiscretization()->GetInterfaceNodes();
        const vector<int>& s_nodes  = GetAleDiscretization()->GetSolidL2G();
        const vector<int>& f_nodes  = GetAleDiscretization()->GetFluidL2G();
        vector<int> cv;
        cv.push_back(0);
        for (int i = 0; i < s_nodes.size(); ++i)
            if (i_nodes.find(s_nodes[i]) == i_nodes.end())
                GetMatrix()->dirichlet_only_row(s_nodes[i], cv);
    }
}

template <int DIM>
DiscretizationInterface* FSISolver<DIM>::NewDiscretization(int dimension, const string& discname)
{
    if (dimension == 2)
    {
        if (discname == "AleQ1")
            return new AleQ12d;
        else if (discname == "AleQ1Lps")
            return new AleQ1Lps2d;
        else if (discname == "AleQ2Lps")
            return new AleQ2Lps2d;
        else
            return StdSolver::NewDiscretization(dimension, discname);
    }
    else if (dimension == 3)
    {
        if (discname == "AleQ1Lps")
            return new AleQ1Lps3d;
        else if (discname == "AleQ2Lps")
            return new AleQ2Lps3d;
        return StdSolver::NewDiscretization(dimension, discname);
    }
    else
        abort();
}

template <>
void FSISolver<2>::reinit_element(int en, const nvector<int>& indices,
                                  HASHMAP<int, std::vector<int>>& solid_interface_cells,
                                  HASHMAP<int, std::vector<int>>& fluid_interface_cells,
                                  HASHSET<int>& fluid_cells, HASHSET<int>& solid_cells,
                                  HASHSET<int>& interface_nodes, set<int>& fluid_nodes,
                                  set<int>& solid_nodes)
{
    Chi chi;

    const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*>(GetMesh());
    assert(M);

    int nf = 0;
    int ns = 0;

    vector<int> ni;
    for (int i = 0; i < indices.size(); ++i)
    {
        int domain = chi(M->vertex2d(indices[i]));
        if (domain > 0)
        {
            ++ns;
            solid_nodes.insert(indices[i]);
        }
        if (domain < 0)
        {
            ++nf;
            fluid_nodes.insert(indices[i]);
        }
        if (domain == 0)
        {
            ni.push_back(i);
            fluid_nodes.insert(indices[i]);
            solid_nodes.insert(indices[i]);
            interface_nodes.insert(indices[i]);
        }
    }

    if ((ns > 0) && (nf > 0))
    {
        cerr << "Geht nicht, fluid & solid!" << endl;
        cerr.precision(20);

        for (int i = 0; i < indices.size(); ++i)
            cerr << M->vertex2d(indices[i]) << "\t" << chi(M->vertex2d(indices[i])) << endl;

        abort();
    }
    if (ni.size() > 0)
    {
        if (ns > 0)
        {
            solid_interface_cells[en] = ni;
            solid_cells.insert(en);
        }
        else if (nf > 0)
        {
            fluid_interface_cells[en] = ni;
            fluid_cells.insert(en);
        }
        else
        {
            solid_interface_cells[en] = ni;
            cout << "Element has interface everywhere!!!" << endl;
        }

        //	if (nf>0) cout << indices << "\t\t" << ni << endl;
    }
    else
    {
        if (ns == indices.size())
            solid_cells.insert(en);
        else if (nf == indices.size())
            fluid_cells.insert(en);
        else
            abort();
    }
}

template <>
void FSISolver<3>::reinit_element(int en, const nvector<int>& indices,
                                  HASHMAP<int, std::vector<int>>& solid_interface_cells,
                                  HASHMAP<int, std::vector<int>>& fluid_interface_cells,
                                  HASHSET<int>& fluid_cells, HASHSET<int>& solid_cells,
                                  HASHSET<int>& interface_nodes, set<int>& fluid_nodes,
                                  set<int>& solid_nodes)
{
    Chi chi;

    const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*>(GetMesh());
    assert(M);

    int nf = 0;
    int ns = 0;

    vector<int> ni;
    for (int i = 0; i < indices.size(); ++i)
    {
        int domain = chi(M->vertex3d(indices[i]));
        if (domain > 0)
        {
            ++ns;
            solid_nodes.insert(indices[i]);
        }
        if (domain < 0)
        {
            ++nf;
            fluid_nodes.insert(indices[i]);
        }
        if (domain == 0)
        {
            ni.push_back(i);
            fluid_nodes.insert(indices[i]);
            solid_nodes.insert(indices[i]);
            interface_nodes.insert(indices[i]);
        }
    }
    if ((ns > 0) && (nf > 0))
    {
        cerr << "Geht nicht, fluid & solid!" << endl;
        abort();
    }
    if (ni.size() > 0)
    {
        if (ns > 0)
        {
            solid_interface_cells[en] = ni;
            solid_cells.insert(en);
        }
        else if (nf > 0)
        {
            fluid_interface_cells[en] = ni;
            fluid_cells.insert(en);
        }
        else
        {
            solid_interface_cells[en] = ni;
            cout << "Element has interface everywhere!!!" << endl;
        }

        //	if (nf>0) cout << indices << "\t\t" << ni << endl;
    }
    else
    {
        if (ns == indices.size())
            solid_cells.insert(en);
        else if (nf == indices.size())
            fluid_cells.insert(en);
        else
        {
            cerr << ns << " " << nf << "\t" << indices.size() << endl;

            abort();
        }
    }
}

template <int DIM>
void FSISolver<DIM>::ReInitInterface(AleBaseDiscretization* ALEDISC)
{
    HASHMAP<int, std::vector<int>>& solid_interface_cells = ALEDISC->GetSolidInterfaceCells();
    HASHMAP<int, std::vector<int>>& fluid_interface_cells = ALEDISC->GetFluidInterfaceCells();
    HASHSET<int>& interface_nodes                         = ALEDISC->GetInterfaceNodes();
    HASHSET<int>& fluid_cells                             = ALEDISC->GetFluidCells();
    HASHSET<int>& solid_cells                             = ALEDISC->GetSolidCells();
    vector<int>& fluid_l2g                                = ALEDISC->GetFluidL2G();
    vector<int>& solid_l2g                                = ALEDISC->GetSolidL2G();
    HASHMAP<int, int>& fluid_g2l                          = ALEDISC->GetFluidG2L();
    HASHMAP<int, int>& solid_g2l                          = ALEDISC->GetSolidG2L();

    set<int> fluid_nodes, solid_nodes;

    solid_interface_cells.clear();
    fluid_interface_cells.clear();
    interface_nodes.clear();
    fluid_cells.clear();
    solid_cells.clear();

    int dim = GetMesh()->dimension();

    if (dim == 2)
    {
        const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*>(GetMesh());
        assert(M);
        if ((GetDiscretization()->GetName() == "Q1 Ale 2d Lps")
            || (GetDiscretization()->GetName() == "Q1 Ale 2d"))
            for (int c = 0; c < M->ncells(); ++c)
                reinit_element(c, M->IndicesOfCell(c), solid_interface_cells, fluid_interface_cells,
                               fluid_cells, solid_cells, interface_nodes, fluid_nodes, solid_nodes);

        else if ((GetDiscretization()->GetName() == "Q2 Ale 2d Lps")
                 || (GetDiscretization()->GetName() == "Q2 Ale 2d"))
            for (int c = 0; c < M->npatches(); ++c)
                reinit_element(c, *(M->IndicesOfPatch(c)), solid_interface_cells,
                               fluid_interface_cells, fluid_cells, solid_cells, interface_nodes,
                               fluid_nodes, solid_nodes);
        else
            abort();
    }
    else if (dim == 3)
    {
        const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*>(GetMesh());
        assert(M);

        if (GetDiscretization()->GetName() == "Q1 Ale 3d Lps")
            for (int c = 0; c < M->ncells(); ++c)
                reinit_element(c, M->IndicesOfCell(c), solid_interface_cells, fluid_interface_cells,
                               fluid_cells, solid_cells, interface_nodes, fluid_nodes, solid_nodes);
        else if (GetDiscretization()->GetName() == "Q2 Ale 3d Lps")
            for (int c = 0; c < M->npatches(); ++c)
                reinit_element(c, *(M->IndicesOfPatch(c)), solid_interface_cells,
                               fluid_interface_cells, fluid_cells, solid_cells, interface_nodes,
                               fluid_nodes, solid_nodes);
        else
        {
            std::cout << GetDiscretization()->GetName() << std::endl;

            abort();
        }
    }
    else
        abort();

    // Nodes Fluid & Solid,  local <-> global (fluid nodes include interface and same for solid)
    fluid_l2g.clear();
    solid_l2g.clear();
    fluid_g2l.clear();
    solid_g2l.clear();

    // l2g
    for (set<int>::const_iterator it = fluid_nodes.begin(); it != fluid_nodes.end(); ++it)
        fluid_l2g.push_back(*it);
    for (set<int>::const_iterator it = solid_nodes.begin(); it != solid_nodes.end(); ++it)
        solid_l2g.push_back(*it);

    // g2l
    for (int i = 0; i < fluid_l2g.size(); ++i)
        fluid_g2l[fluid_l2g[i]] = i;
    for (int i = 0; i < solid_l2g.size(); ++i)
        solid_g2l[solid_l2g[i]] = i;
}

// --------------------------------------------------

template <int DIM>
void FSISolver<DIM>::NewMesh(const MeshInterface* mp)
{
    StdSolver::NewMesh(mp);

    ReInitInterface(GetAleDiscretization());
}

template <>
void FSISolver<2>::PointVisu(const string& name, const GlobalVector& u, int iter) const
{
    GlobalVector U;
    U.ncomp() = u.ncomp() + 1;
    assert(2 * 2 + 2 == U.ncomp());
    U.resize(u.n());

    const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*>(GetMesh());
    assert(M);
    Chi chi;

    for (int i = 0; i < u.n(); ++i)
    {
        const Vertex2d v = M->vertex2d(i);
        int domain       = chi(v);
        for (int c = 0; c < u.ncomp(); ++c)
            U(i, c) = u(i, c);
        U(i, u.ncomp()) = domain;
    }
    StdSolver::PointVisu(name, U, iter);
}
template <>
void FSISolver<3>::PointVisu(const string& name, const GlobalVector& u, int iter) const
{
    GlobalVector U;
    U.ncomp() = u.ncomp() + 1;
    assert(2 * 3 + 2 == U.ncomp());
    U.resize(u.n());

    const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*>(GetMesh());
    assert(M);
    Chi chi;

    for (int i = 0; i < u.n(); ++i)
    {
        const Vertex3d v = M->vertex3d(i);
        int domain       = chi(v);
        for (int c = 0; c < u.ncomp(); ++c)
            U(i, c) = u(i, c);
        U(i, u.ncomp()) = domain;
    }
    StdSolver::PointVisu(name, U, iter);
}

template class FSISolver<2>;
template class FSISolver<3>;
template class FSIMultiLevelSolver<2>;
template class FSIMultiLevelSolver<3>;
}  // namespace Gascoigne
