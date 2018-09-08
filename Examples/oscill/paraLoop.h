#ifndef PARAREAL_H
#define PARAREAL_H

#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <limits>
#include <iostream>

#include "loop.h"
#include "weightedpointfunctional.h"
#include "resfunctional.h"
#include "local.h"
#include "colors.h"
#include "log.h"

namespace Gascoigne
{
enum class name_type
{
    error,
    functional
};
enum class iter_type
{
    fine_last    = 0,
    fine         = 1,
    coarse_first = 2,
    coarse       = 3,
};
enum class log_level
{
    results       = 0,
    visu_info     = 1,
    visu_detailed = 2
};

struct dir_vec
{
    std::string dir;
    VectorInterface vec;
    size_t index;
};

/*!
 * \brief This class implements the [parareal-algorithm](https://en.wikipedia.org/wiki/Parareal)
 * for parallelizing the solution of PDEs using the FEM libary Gascoigne.
 *
 * # The algorithm
 * Parareal is an algorithm parallelizing the solution of PDEs in time using
 * a fine and a coarse propagator which operate on a time discretization. This decomposes the
 * time interval into multiple subintervals, where each of these are assigned to one thread. The
 * coarse propagator \f$G\f$ is used for predicting values at the end of each subinterval which
 * the fine propagator \f$F\f$ uses to correct the solution.
 * ### Sketch of the algorithm
 * 1. Initialization
 *
 *    Initialize values on each subinterval \f$i = 0,\dots,\,N-1\f$ with the coarse propagator:
 *    \f[g_{i} \longleftarrow G(t_{i},\,t_{i-1},\,g_{i-1})\f]
 *    Set iteration counter \f$k = 1\f$
 * 2. Parallelization
 *
 *    Parallelize fine propagations over subintervals with index \f$i = k,\dots,\,N+1\f$:
 *    \f[f_i\longleftarrow F(t_i,\,t_{i-1},\,u_{i-1})\f]
 * 3. Corrector term
 *
 *    Compute corrections for all subintervals with index \f$i = k+1,\dots,\,N+1\f$:
 *    \f[f_i \longleftarrow f_i-g_i\f]
 * 4. Coarse prediction
 *
 *    On each subinterval with index \f$i = k+1,\dots,\,N+1\f$:
 *    \f[g_i \longleftarrow G(t_{i},\,t_{i-1},\,u_{i-1})\f]
 * 5. Correction
 *
 *    On each subinterval \f$i = k+1,\dots,\,N+1\f$:
 *    \f[u_i \longleftarrow g_i + f_i\f]
 *    Set \f$k\longleftarrow k+1\f$ and go back to 2.
 *
 * # The implementation
 * The algorithm depends on Gascoigne and the
 * correct definition of the problem, there is no checking. The implementation only executes the
 * bare algorithm. Checking for convergence is not available at the moment. This can be done via
 * comparing serial with parallel execution. Some properties
 * - Parallelization with OpenMP
 * - Inherits from Loop
 * - Memory consumption is higher than sequential
 * - Timestepping with theta-scheme - with different theta for fine and coarse solutions.
 *
 * Every Instance of parareal correponds to one solver on an subinterval. Thus every solver also
 * has its own map of vectors. The user
 * does not have to interact directly with them. He only needs to call runPara or
 * compare_para_serial, everything else is handled automatically by the class. For better
 * debugging one can specify the `PDBG` flag (i.e. `-PDBG=n`), where`n` can be 1 or something
 * arbitrary larger. In the case where `n=1` only some more timings are taken. Otherwise
 * additional visualization will be written.
 */
template <int DIM, log_level logging>
class parareal : public Loop<DIM>
{
public:
    /*!
     * \brief Constructor for parareal object.
     *
     * Every instance needs the index of the subinterval to identify the vectors of the fine,
     * coarse and final solution internally \param id index of parareal instance which
     * corresponds to the subinterval index.
     */
    explicit parareal(std::size_t id);

    /*!
     * \brief Output parareal specific parameters
     *
     * \param ostrm std::ostream where data should be written; cerr by default but can also be a
     * file (ofstream)
     */
    static void outputParameters(std::ostream& ostrm = std::cerr);

    /*! \brief Method setting parameters and calling the algorithm.
     * This method takes some parameters (maybe from command line) and sets them and calls
     * paraAlgo afterwards to execute parareal algorithm. If no parameters are given, the
     * parameters will be read from the paramfile. Only parameters which may be varied in some
     * tests are included.
     * \param max_iterations maximal number of iterations executed by parareal algorithm
     * \param coarse_theta \f$\theta\f$ for coarse time stepping
     * \param dtcoarse coarse time step size
     * \param theta \f$\theta\f$ for fine time stepping
     * \param dtfine fine time step size
     */

    static double runPara(const int max_iterations  = 1,
                          const double coarse_theta = get_param<double>("thetacoarse", "Equation"),
                          const double dtcoarse     = get_param<double>("dtcoarse", "Equation"),
                          const double fine_theta   = get_param<double>("theta", "Equation"),
                          const double dtfine       = get_param<double>("dt", "Equation"));

    /*! \brief Method comparing parareal with sequential execution.
     *
     * Executes the parareal algorithm with the specified coarse time step sizes and coarse
     * thetas, values for fine method is taken from the paramfile. Before execution of parareal
     * it solves the problem sequentially and calculates the norm \f$l_2\f$ vector norm between
     * the two solutions on every time node. \param max_iterations maximal number of iterations
     * executed by parareal algorithm \param dtcoarse_vec Vector with coarse time step sizes
     * that will be tested against sequential execution \param coarse_theta_vec Vector with
     * coarse theta values for step sizes specified by dtcoarse_vec
     */

    static void compare_para_serial(const int max_iterations,
                                    const std::vector<double>& dtcoarse_vec,
                                    const std::vector<double>& coarse_theta_vec);

    static ParamFile
      paramfile;  //!< Name of the file in which all problem-related parameters are saved in.
    static std::string problemlabel;  //!< Name of the problem.

protected:
    /*!
     * \brief Wrapper function for GetSolver calls via GetMultiLevelSolver returning *const*
     * pointer to Solver
     */
    const SolverInterface* GetSolver() const;

    /*!
     * \brief Wrapper function for GetMultiLevelSolver returning *const* pointer to
     * MultiLevelSolver
     */
    const MultiLevelSolverInterface* GetMultiLevelSolver() const;

private:
    /*!
     * \brief Wrapper function for GetSolver calls via GetMultiLevelSolver returning *nonconst*
     * pointer to Solver
     */
    SolverInterface* GetSolver();

    /*!
     * \brief Wrapper function for GetMultiLevelSolver returning *nonconst* pointer to
     * MultiLevelSolver
     */
    MultiLevelSolverInterface* GetMultiLevelSolver();

    /*! \brief Method executing parareal-algorithm.
     *
     * This is the method that is called after runPara has set all parameters for
     * parareal.
     */
    static double paraAlgo();

    /*!
     * \brief Time stepping method.
     *
     * This method encapsulates *n_finesteps* or *n_coarsesteps* timesteps into one call.
     * The step size used is dtfine or *dtcoarse* (depending on whether coarse or fine is
     * specified by the template parameter). This method writes a logfile if method equals
     * fine_last. \param time          current time from where the steps are executed \param u
     * VectorInterface with approximation at current time \param f             Right hand inside
     */
    template <iter_type method>
    void propagator(double time, VectorInterface& u, VectorInterface& f, const int current = 0,
                    const int n_steps = n_finesteps, const int implicit_steps = 0);

    template <iter_type method>
    void inline step(double& time, const double& dt, const int& iter, const double& theta,
                     VectorInterface& u, VectorInterface& f);

    void inline divergence_stab(VectorInterface& u, VectorInterface& f, const double& time,
                                const double& dt, const int& current);

    void visu_write(const std::string& dir, const VectorInterface& vec, size_t index) const;
    void visu(const std::string& dir, const VectorInterface& vec, size_t index) const;
    void inline visu(const dir_vec& dv) const;

    template <typename... dir_vec_pack>
    void inline visu(const dir_vec& dv, const dir_vec_pack&... dvs) const;
    /*!
     * \brief Method checking for convergence on one subinterval
     *
     * This method checks if the difference between the solution from last iteration
     * and the fine propagations from this iteration is lower than...
     */
    bool converged(GlobalVector& fine_v, GlobalVector& u_v);

    /*!
     * \brief Setting time and time step size
     *
     * This method is important for setting the different coarse and fine timesteps
     * in the equation.
     * Typically it is called inside the step method before the steps are executed.
     * \param dt is the time step size
     * \param time is for time dependent problems and zero by default
     */
    void set_time(const double dt, const double time = 0.) const;

    /*!
     * \brief Utility method for assigning a VectorInterface to a specific GlobalVector.
     *
     * Every Vector is represented by a VectorInterface (a string) which itself is just a key in
     * a map associated to a solver. If one wants to change the Vector itself it does not
     * suffice to copy the VectorInterface but one rather has to change the value in the map.
     * \param v Key in the map associated with the solver
     * \param v_ptr Value in the map associated with the solver
     */
    void setGV(VectorInterface& v, const GlobalVector& v_ptr) const;

    /*!
     * \brief Utility method for preallocating GlobalVector with correct size.
     *
     * The method calculates the correct size of a GlobalVector with parameters of the solver.
     * \param v_ptr Adress of the GlobalVector to be modified
     */
    void setGVsize(GlobalVector& v_ptr);
    void setGVsize(const VectorInterface& vec);

    /*!
     *\brief Setting up the problem.
     *
     * This method initializes:
     * - ProblemDescriptor
     * - ProblemContainer
     * - etc...
     */
    void setup_problem(const ParamFile& paramfile, const std::string& problemlabel);

    /*!
     * \brief Initialization method for setting the number of intervals.
     *
     * The method sets the number of the subintervals to an already defined number of threads
     * (with `omp_get_thread_num()`), or to the maximum number of threads available with
     * `omp_get_max_threads()`. Setting the number of subintervals to the maximum is good for
     * using the maximum potential of parallelization. If threads are defined in the paramfile,
     * they are used.
     */
    static void set_interval_threads(int& n_intervals, int& threads);

    /*!
     * \brief Method for reading parameters from the paramfile.
     *
     * Before calling this the `paramfile` has to be defined already.
     * \param paramName Name of parameter which will be searched for
     * \param block Block in which paramName is searched
     */
    template <typename T>
    static const T get_param(const std::string& paramName, const std::string& block);

    template <name_type name>
    static const std::string get_name(double dtc = dtcoarse, double tc = coarse_theta,
                                      int max_iter = max_iterations);

    std::size_t subinterval_idx;  //!< Identifies the subinterval on which the loop is acting in the
                                  //!< first iteration
    VectorInterface
      fine_sol;  //!< Vector containing the fine solution on the corresponding interval.
    VectorInterface
      coar_sol;  //!< Vector containing the coarse solution on the corresponding interval.
    VectorInterface old, tmp_vi;  //!< Vector for temporary results.
    VectorInterface
      end_sol;  //!< Vector containing the final solution on the corresponding interval.
    std::vector<DoubleVector> log_buffer;  //!< Vector for buffering log entries.

    const int c_implicit_steps = get_param<int>("implicit_coarse", "Equation");
    const int f_implicit_steps = get_param<int>("implicit_fine", "Equation");

    static double fine_theta;    //!< Theta for fine timestepping
    static double coarse_theta;  //!< Theta for coarse timestepping

    static double start_time;       //!< Start time point
    static double stop_time;        //!< End time point
    static double precompute_time;  //!< Time from which it should run parallel
    static int n_intervals;         //!< Number of intervals parareal is using. Usually this
                                    //!< is the maximum number of threads (can be higher).
    static int threads;             //!< Number of threads being used.

    static int max_iterations;      //!< Maximum number of iterations parareal should execute.
    static double interval_length;  //!< Length of one subinterval parareal uses

    static double dtfine;            //!< Time step size of the fine propagator
    static std::size_t n_finesteps;  //!< Number of fine timesteps per subinterval

    static double dtcoarse;            //!< Time step size of the coarse propagator
    static std::size_t n_coarsesteps;  //!< Number of coarse timesteps per subinterval

    static double coarse_exec_time;
    static double fine_exec_time;

    static std::vector<GlobalVector>
      u_sol_arr;                    //!< Vector containing the final solutions on every time node
    static std::ofstream func_log;  //! Ofstream object for logfile writing
};
// end class definition

/***************************************************************************************************
public methods
***************************************************************************************************/

template <int DIM, log_level logging>
parareal<DIM, logging>::parareal(size_t id)
    : subinterval_idx(id)
    , fine_sol("f" + std::to_string(id))
    , coar_sol("g" + std::to_string(id))
    , old("old" + std::to_string(id))
    , tmp_vi("tmp" + std::to_string(id))
    , end_sol("u" + std::to_string(id))
    , log_buffer(n_finesteps, DoubleVector(5, 0.))
{
    // StdLoop();
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::outputParameters(std::ostream& ostrm)
{
    // clang-format off
    ostrm << "\n========================================================="
          << "\nStart time                             " << start_time
          << "\nEnd time                               " << stop_time
          << "\nPrecompute time                        " << precompute_time
          << "\nNumber of threads is                   " << threads
          << "\nNumber of subintervals is              " << n_intervals
          << "\nInterval length                        " << interval_length
          << "\nMaximal number of iterations           " << max_iterations
          << "\nParameter for coarse method:\
              \nCoarse steps per interval:             " << n_coarsesteps
          << "\nCoarse step size:                      " << dtcoarse
          << "\nTheta for coarse steps                 " << coarse_theta
          << "\nParameter for fine method:\
              \nFine steps per interval:               " << n_finesteps
          << "\nFine step size:                        " << dtfine
          << "\nTheta for fine steps                   " << fine_theta
          << "\nUsing paramfile:                       " << paramfile
          << "\nProblemlabel is:                       " << problemlabel
          << "\n=========================================================\n\n";
    // clang-format on
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
double parareal<DIM, logging>::runPara(const int max_iterations, const double coarse_theta,
                                       const double dtcoarse, const double fine_theta,
                                       const double dtfine)
{
    assert(dtcoarse >= dtfine && "coarse step size should be bigger than fine step size");
    parareal::max_iterations = max_iterations;
    parareal::coarse_theta   = coarse_theta;
    parareal::dtcoarse       = dtcoarse;
    parareal::fine_theta     = fine_theta;
    parareal::dtfine         = dtfine;
    parareal::n_coarsesteps  = static_cast<size_t>(interval_length / dtcoarse);
    parareal::n_finesteps    = static_cast<size_t>(interval_length / dtfine);
    set_interval_threads(parareal::n_intervals, parareal::threads);
    auto ratio = static_cast<int>(dtcoarse / dtfine);
    if (ratio * n_coarsesteps != n_finesteps)
    {
        n_finesteps += n_finesteps % ratio;
        n_coarsesteps   = static_cast<int>(n_finesteps / ratio);
        interval_length = n_finesteps * dtfine;
        stop_time       = start_time + precompute_time + interval_length * n_intervals;
    }

    std::cerr << logwrite::info("Info", "Initialized Parameters");
    outputParameters();
    func_log.precision(12);
    // Checking if input is correct
    assert(isfinite(dtcoarse) && isfinite(dtfine) && "At least one timestep is infinite.");
    assert(isfinite(n_coarsesteps) && isfinite(n_finesteps) && "Amount of steps is infinite.");
    assert(isfinite(interval_length) && isfinite(start_time)
           && "interval_length or start time is infinite");
    assert(precompute_time >= 0 && "precompute time has to be greater or equal zero");
    assert(stop_time - start_time > precompute_time && "start has to be be before stop");

    assert(dtcoarse > 0 && dtfine > 0 && "At least one timestep is less or equal 0.");
    assert(n_coarsesteps != 0 && n_finesteps != 0 && "Amount of steps is 0.");
    assert(n_intervals > 1
           && "Number of intervals is 1 or 0. You should prefer sequential execution.");
    assert(interval_length > 0 && "interval_length should be greater zero!");
    assert(fine_theta > 0 && coarse_theta > 0 && "Thetas have to be positive");
    assert(dtcoarse * n_coarsesteps - dtfine * n_finesteps <= numeric_limits<double>::epsilon()
           && "Coarse and fine step sizes and number of time steps don't match");

    if (!func_log.is_open())
    {
        func_log.open(get_name<name_type::functional>());
    }
    auto exec_time = paraAlgo();

    std::cerr << sty::bb << "\n==================================================\n"
              << sty::g << "Finished parareal" << sty::b << ".\nElapsed time is\t" << exec_time
              << "\n\nTime spent in coarse method: " << coarse_exec_time
              << "\nTime spent in fine method: " << fine_exec_time << sty::n << '\n';
    func_log.close();
    return exec_time;
}

//--------------------------------------------------------------------------------------------------

#ifdef _OPENMP  // parareal can't be run without OpenMP
template <int DIM, log_level logging>
double parareal<DIM, logging>::paraAlgo()
{
    static auto equal = [&](parareal*& para_dest, VectorInterface& dest, parareal*& para_src,
                            VectorInterface& source) {
        para_dest->GetSolver()->GetGV(dest).equ(1., para_src->GetSolver()->GetGV(source));
    };
    static auto w_equ_x_plus_y_minus_z =
      [&](parareal*& para_dest, VectorInterface& dest, parareal*& src_x, VectorInterface& x,
          parareal*& src_y, VectorInterface& y, parareal*& src_z, VectorInterface& z) {
          // clang-format off
        para_dest->GetSolver()->GetGV(dest).equ(1,  src_x->GetSolver()->GetGV(x),
                                                1,  src_y->GetSolver()->GetGV(y),
                                                -1, src_z->GetSolver()->GetGV(z));
          // clang-format on
      };
    // auto add = [&](VectorInterface& dest, VectorInterface& source){
    //
    // };
    double time     = start_time;
    auto start_wall = omp_get_wtime();
    //
    // Initialization of solvers and problems and solutions from parareal
    // coar_sol   - VectorInterface for results of coarse method
    // fine_sol   - VectorInterface for results of fine method
    // end_sol    - VectorInterface for final results
    // old           - VectorInterface for temporary results
    // u             - VectorInterface used for initialization of the others
    //
    VectorInterface u("u");
    VectorInterface f("f");
    std::vector<parareal*> subinterval(threads);
    for (auto m = 0; m < threads; ++m)
    {
        subinterval[m] = new parareal<DIM, logging>(m + 1);
    }
    for (const auto& si : subinterval)
    {
        si->setup_problem(paramfile, problemlabel);

        si->GetMultiLevelSolver()->ReInit(problemlabel);
        si->GetSolver()->OutputSettings();
        cout << si->GetMeshAgent()->ncells() << " cells" << '\n';

        si->GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
        si->set_time(dtfine, start_time);

        si->GetMultiLevelSolver()->ReInitVector(f);
        si->GetMultiLevelSolver()->ReInitVector(u);
        si->GetMultiLevelSolver()->ReInitVector(si->fine_sol);
        si->GetMultiLevelSolver()->ReInitVector(si->coar_sol);
        si->GetMultiLevelSolver()->ReInitVector(si->old);
        si->GetMultiLevelSolver()->ReInitVector(si->tmp_vi);
        si->GetMultiLevelSolver()->ReInitVector(si->end_sol);

        si->setGVsize(si->GetSolver()->GetGV(si->fine_sol));
        si->setGVsize(si->GetSolver()->GetGV(si->coar_sol));
        si->setGVsize(si->GetSolver()->GetGV(si->old));
        si->setGVsize(si->GetSolver()->GetGV(si->tmp_vi));
        si->setGVsize(si->GetSolver()->GetGV(si->end_sol));

        si->InitSolution(u);
        si->InitSolution(si->fine_sol);
        si->InitSolution(si->coar_sol);
        si->InitSolution(si->end_sol);
    }

    for (auto m = 0; m < n_intervals; ++m)
    {
        subinterval[m]->setGVsize(u_sol_arr[m]);
    }
    // precompute till given time
    if (precompute_time > 0)
    {
        auto precompute_steps = static_cast<int>(precompute_time / dtfine);
        subinterval[0]->template propagator<iter_type::fine>(time, u, f, 0, precompute_steps, 0);
        subinterval[0]
          ->GetSolver()
          ->GetGV(subinterval[0]->end_sol)
          .equ(1., subinterval[0]->GetSolver()->GetGV(u));
        for (auto m = 0; m < n_intervals; ++m)
        {
            equal(subinterval[m], subinterval[m]->end_sol, subinterval[0], u);
            equal(subinterval[m], subinterval[m]->fine_sol, subinterval[0], u);
            equal(subinterval[m], subinterval[m]->coar_sol, subinterval[0], u);
            if (m > 0)
            {
                equal(subinterval[m], u, subinterval[0], u);
            }
        }
        time += precompute_time;
    }
    /// Initializing locks
    auto interval_locker = new omp_lock_t[n_intervals];
    for (auto m = 0; m < n_intervals; ++m)
    {
        omp_init_lock(&interval_locker[m]);
    }

//
// A C T U A L   A L G O R I T H M
//
#pragma omp parallel num_threads(n_intervals) firstprivate(time)  // proc_bind(close)
    {
        // 'Iteration 0'
        //
#pragma omp for ordered nowait schedule(monotonic : static)
        // Initialization
        for (auto i = 0; i < n_intervals; ++i)
        {
#pragma omp ordered
            {
                // auto im  = i % threads;
                // auto im1 = (i + 1) % threads;
                omp_set_lock(&interval_locker[i]);
                // u ⟵ G⁰(i)
                subinterval[i]->template propagator<iter_type::coarse_first>(
                  time + i * interval_length, u, f, 0, n_coarsesteps,
                  subinterval[i]->c_implicit_steps);
                // G⁰(i) ⟵ u;
                equal(subinterval[i], subinterval[i]->coar_sol, subinterval[i], u);
                // U⁰(i) ⟵ u;
                equal(subinterval[i], subinterval[i]->end_sol, subinterval[i], u);

                // Visualization for debugging
                if constexpr (logging == log_level::visu_info
                              || logging == log_level::visu_detailed)
                {
                    std::cerr << logwrite::info("Done", "Initialization on subinterval ", i);
                    subinterval[i]->visu("Results/initsol", subinterval[i]->end_sol, i);
                }

                if (i < n_intervals - 1)
                {
                    // omp_set_lock(&interval_locker[i + 1]);
                    // pass GlobalVector around
                    subinterval[i + 1]->setGV(u, subinterval[i]->GetSolver()->GetGV(u));
                    // F⁰(i) ⟵ u
                    equal(subinterval[i + 1], subinterval[i + 1]->fine_sol, subinterval[i + 1], u);
                    // omp_unset_lock(&interval_locker[i + 1]);
                }
                omp_unset_lock(&interval_locker[i]);
            }
        }  // Initialization

        // for (auto n = 0; n <= n_intervals - threads; n++)
        // {
        //
        // #pragma omp barrier  // barrier for debugging
        //
        // Do <max_iterations> parareal-iterations
        // Iterations
        for (auto k = 1; k <= max_iterations; ++k)
        {
            // fine propagations
#pragma omp for nowait schedule(static)
            for (auto m = 0; m < n_intervals - k + 1; ++m)
            {
                // #pragma omp ordered  // for debugging
                //                 {
                //omp_set_lock(&interval_locker[m]);
                // Fᵏ(m) ⟵ F(Uᵏ⁻¹(m-1))
                if (m == 0 || k == max_iterations)
                {
                    omp_set_lock(&interval_locker[m]);
                    subinterval[m]->template propagator<iter_type::fine_last>(
                      time + m * interval_length, subinterval[m]->fine_sol, f, k, n_finesteps,
                      subinterval[m]->f_implicit_steps);
                    omp_unset_lock(&interval_locker[m]);
                }
                else  // we are not in the last iteration and not exact
                {
                    omp_set_lock(&interval_locker[m]);
                    subinterval[m]->template propagator<iter_type::fine>(
                      time + m * interval_length, subinterval[m]->fine_sol, f, k, n_finesteps,
                      subinterval[m]->f_implicit_steps);
                    omp_unset_lock(&interval_locker[m]);
                }
                if constexpr (logging == log_level::visu_info
                              || logging == log_level::visu_detailed)
                {
                    std::cerr << logwrite::info("Done", "Fine solution on subinterval ", m + k - 1,
                                                " in iteration number ", k);
                    subinterval[m]->visu("Results/fineres", subinterval[m]->fine_sol,
                                               (k - 1) * 10 + m);
                }
                //omp_unset_lock(&interval_locker[m]);
                //}
            }
// barrier for debugging
// #pragma omp barrier
//
// coarse propagations and corrections
//
#pragma omp for ordered nowait schedule(monotonic : static)
            // Corrections
            for (auto m = 0; m < n_intervals - k; ++m)
            {
#pragma omp ordered
                {
                    // k-th value is exact in k-th Iteration
                    if (m == 0)
                    {
                        omp_set_lock(&interval_locker[0]);
                        // Uᵏ(k) ⟵ Fᵏ(k)
                        equal(subinterval[0], subinterval[0]->end_sol, subinterval[0],
                              subinterval[0]->fine_sol);
                        u_sol_arr[k - 1].equ(
                          1., subinterval[0]->GetSolver()->GetGV(subinterval[0]->fine_sol));
                        for (const auto& x : subinterval[0]->log_buffer)
                        {
                            func_log << x << '\n';
                        }
                        subinterval[0]->visu("Results/p", subinterval[0]->end_sol, k - 1);
                        omp_unset_lock(&interval_locker[0]);
                    }

                    // prepare predictor step
                    omp_set_lock(&interval_locker[m]);
                    equal(subinterval[m], subinterval[m]->coar_sol, subinterval[m],
                          subinterval[m]->end_sol);
                    omp_unset_lock(&interval_locker[m]);
                    // predictor step with coarse method
                    subinterval[m]->template propagator<iter_type::coarse>(
                      time + (m + 1) * interval_length, subinterval[m]->coar_sol, f, k,
                      n_coarsesteps, subinterval[m]->c_implicit_steps);

                    // omp_unset_lock(&interval_locker[m]);
                    omp_set_lock(&interval_locker[m + 1]);

                    // calculate corrector term
                    w_equ_x_plus_y_minus_z(subinterval[m + 1], subinterval[m + 1]->end_sol,
                                           subinterval[m], subinterval[m]->coar_sol,
                                           subinterval[m + 1], subinterval[m + 1]->fine_sol,
                                           subinterval[m + 1], subinterval[m + 1]->coar_sol);

                    // prepare next iteration
                    equal(subinterval[m + 1], subinterval[m + 1]->fine_sol, subinterval[m + 1],
                          subinterval[m + 1]->end_sol);
                    omp_unset_lock(&interval_locker[m + 1]);

                    // Visualization for debugging
                    if constexpr (logging == log_level::visu_info
                                  || logging == log_level::visu_detailed)
                    {
                        subinterval[m]->visu(
                          dir_vec{"Results/coarseres", subinterval[m]->coar_sol, (k - 1) * 10 + m},
                          dir_vec{"Results/iterres", subinterval[m]->end_sol, (k - 1) * 10 + m});
                        std::cerr << logwrite::info("Done", "Calculations on subinterval ", m + 1,
                                                    " Iteration number ", k, " of parareal ");
                    }
                }
            }  // Corrections

#pragma omp atomic update  // increase time
            time += interval_length;
        }  // Iterations

        //}  // compute window loop

    }  // parallel region
    auto exec_time = omp_get_wtime() - start_wall;

    for (auto m = 0; m < n_intervals - max_iterations; ++m)
    {
        // write final solution vectors
        subinterval[m]->setGVsize(u_sol_arr[m]);
        // final solutions on time nodes
        u_sol_arr[m + max_iterations].equ(
          1., subinterval[m + 1]->GetSolver()->GetGV(subinterval[m + 1]->end_sol));
        // logging
        for (const auto& x : subinterval[m + 1]->log_buffer)
        {
            func_log << x << '\n';
        }

        // Visu part
        // subinterval[m]->setGV(u, u_sol_arr[m]);
        subinterval[m + 1]->visu("Results/p", subinterval[m + 1]->end_sol,
                                       m + max_iterations);
        std::cerr << logwrite::info("Done", "Wrote subinterval ", m + max_iterations);
    }

    // Clean up
    omp_destroy_lock(interval_locker);
    delete[] interval_locker;
    for (auto i = 0; i < n_intervals; ++i)
    {
        delete subinterval[i];
    }
    return exec_time;
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::compare_para_serial(const int max_iterations,
                                                 const std::vector<double>& dtcoarse_vec,
                                                 const std::vector<double>& coarse_theta_vec)
{
    // ParamFile paramfile(paramstr);
    if constexpr (logging == log_level::results)
    {
        ofstream timings_log({"timingsI" + std::to_string(max_iterations) + ".txt"});
        timings_log << "ratio time accel\n";

        // serial part
        auto start   = omp_get_wtime();
        auto loopSeq = new parareal<DIM, logging>(0);
        loopSeq->setup_problem(paramfile, problemlabel);
        loopSeq->GetMultiLevelSolver()->ReInit(problemlabel);
        loopSeq->GetSolver()->OutputSettings();
        cout << loopSeq->GetMeshAgent()->ncells() << " cells" << '\n';
        VectorInterface f("f");
        double time = start_time;
        std::vector<VectorInterface> U_seq(n_intervals, loopSeq->old);
        for (auto i = 0; i < n_intervals; ++i)
        {
            U_seq[i].SetName({"useq" + std::to_string(i)});
            U_seq[i].SetType("node");
            loopSeq->GetSolver()->ReInitVector(U_seq[i]);
        }
        for (auto i = 0; i < n_intervals; ++i)
        {
            loopSeq->InitSolution(U_seq[i]);
        }

        loopSeq->GetMultiLevelSolver()->ReInitVector(loopSeq->old);
        loopSeq->GetMultiLevelSolver()->ReInitVector(loopSeq->end_sol);
        loopSeq->GetMultiLevelSolver()->ReInitVector(f);
        loopSeq->InitSolution(loopSeq->end_sol);
        loopSeq->InitSolution(loopSeq->old);

        if (func_log.is_open())
        {
            func_log.close();
        }
        func_log.open(get_name<name_type::functional>(dtcoarse, coarse_theta, 0));
        // sequential run
        for (auto i = 0; i < n_intervals; ++i)
        {
            loopSeq->template propagator<iter_type::fine_last>(time + i * interval_length,
                                                               loopSeq->end_sol, f);
            loopSeq->GetSolver()->Equ(U_seq[i], 1., loopSeq->end_sol);
            // time += interval_length;
            for (const auto& x : loopSeq->log_buffer)
            {
                func_log << x << '\n';
            }
            loopSeq->visu("Results/u", loopSeq->end_sol, i);
        }
        func_log.close();

        auto end         = omp_get_wtime();
        auto serial_time = end - start;
        std::cerr << sty::bb << "\n==================================================\n"
                  << sty::g << "Finished serial execution" << sty::b << ".\nElapsed time is\t"
                  << serial_time << ".\n";
        timings_log << 1 << ' ' << serial_time << ' ' << 1 << '\n';

        assert(
          dtcoarse_vec.size() == coarse_theta_vec.size()
          && "For every dtcoarse there has to be one coarse_theta given and the other way round");

        ofstream error_log;

        for (auto i = 0; i < dtcoarse_vec.size(); ++i)
        {
            func_log.open(get_name<name_type::functional>(dtcoarse_vec[i], coarse_theta_vec[i]));

            // run Parareal
            auto ratio           = static_cast<int>(dtcoarse_vec[i] / dtfine);
            double parallel_time = runPara(max_iterations, coarse_theta_vec[i], dtcoarse_vec[i]);
            timings_log << ratio << ' ' << parallel_time << ' ' << serial_time / parallel_time
                        << '\n';
            func_log.close();

            // calculate error
            error_log.open(get_name<name_type::error>(dtcoarse_vec[i], coarse_theta_vec[i]));
            error_log << "interval error\n";
            for (auto i = 0; i < n_intervals; ++i)
            {
                double norm_seq = loopSeq->GetSolver()->Norm(U_seq[i]);
                assert(u_sol_arr[i].size() == loopSeq->GetSolver()->GetGV(U_seq[i]).size()
                       && "vector size mismatch");
                u_sol_arr[i].add(-1., loopSeq->GetSolver()->GetGV(U_seq[i]));
                double error = u_sol_arr[i].norm();
                std::cerr << "Error on interval no. " << i << ": " << error / norm_seq << '\n';
                error_log << i << ' ' << error / norm_seq << '\n';
                u_sol_arr[i].zero();
            }
            error_log.close();
        }
    }
    else
    {
        std::cerr << "Comparison will only be effective with low log_level (i.e. results)" << '\n';
        std::exit(EXIT_FAILURE);
    }
}

#else  // if OpenMP is not available, exit program when methods are called

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
double parareal<DIM, logging>::paraAlgo()
{
    std::cerr << "Running parareal relies on OpenMP" << '\n';
    std::exit(EXIT_FAILURE);
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::compare_para_serial(const int max_iterations,
                                                 const std::vector<double>& dtcoarse_vec,
                                                 const std::vector<double>& coarse_theta_vec)
{
    std::cerr << "Running parareal relies on OpenMP, comparison makes no sense!" << '\n';
    std::exit(EXIT_FAILURE);
}

#endif  //_OPENMP

/***************************************************************************************************
protected methods
***************************************************************************************************/

template <int DIM, log_level logging>
const SolverInterface* parareal<DIM, logging>::GetSolver() const
{
    return StdLoop::GetMultiLevelSolver()->GetSolver();
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
SolverInterface* parareal<DIM, logging>::GetSolver()
{
    return StdLoop::GetMultiLevelSolver()->GetSolver();
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
const MultiLevelSolverInterface* parareal<DIM, logging>::GetMultiLevelSolver() const
{
    return StdLoop::GetMultiLevelSolver();
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
MultiLevelSolverInterface* parareal<DIM, logging>::GetMultiLevelSolver()
{
    return StdLoop::GetMultiLevelSolver();
}

/***************************************************************************************************
private methods
***************************************************************************************************/

template <int DIM, log_level logging>
template <iter_type method>
void parareal<DIM, logging>::propagator(double time, VectorInterface& u, VectorInterface& f,
                                        int current, int n_steps, int implicit_steps)
{
    // lambda for visualisation of steps
    auto step_visu = [&](int iter, VectorInterface& u) {
        if constexpr (method == iter_type::fine || method == iter_type::fine_last)
        {
            if (iter % 10 == 0 && iter < 310)
            {
                std::cerr << logwrite::info("Info", "Fine Visualization");
                auto step = static_cast<int>(time / dtfine);
                visu("Results_detail/subinterval" + std::to_string(subinterval_idx + current - 1)
                       + "/iter" + std::to_string(current) + "/u_fine_no",
                     u, current * static_cast<int>(pow(10, ceil(log10(step)))) + step);
            }
        }
        else if constexpr (method == iter_type::coarse || method == iter_type::coarse_first)
        {
            if (iter < 31)
            {
                std::cerr << logwrite::info("Info", "Coarse Visualization");
                auto step = static_cast<int>(time / dtcoarse);
                visu("Results_detail/subinterval" + std::to_string(subinterval_idx + current)
                       + "/iter" + std::to_string(current) + "/u_coarse_no",
                     u, current * static_cast<int>(pow(10, ceil(log10(step)))) + step);
            }
        }
    };

    double dt;
    double theta;
    auto start = omp_get_wtime();
    // Some compile-time branching...
    if constexpr (method == iter_type::fine_last || method == iter_type::fine)
    {
        dt    = dtfine;
        theta = fine_theta;
    }
    else if constexpr (method == iter_type::coarse || method == iter_type::coarse_first)
    {
        dt    = dtcoarse;
        theta = coarse_theta;
    }
    // stabilization problem if not in initialization
    // if constexpr (method != iter_type::coarse_first)
    // {
    //     divergence_stab(u, f, time, dt, current);
    // }
    // std::cout << "Current iteration: " <<  current <<'\n';
    auto fsi_eq_ptr =
      dynamic_cast<const FSI<DIM>* const>(GetSolver()->GetProblemDescriptor()->GetEquation());
    if (implicit_steps > 0)
    {
        fsi_eq_ptr->getTheta() = 1.;
        for (auto iter = 0; iter < implicit_steps; ++iter)
        {
            step<method>(time, dt, iter, 1., u, f);
            if constexpr (logging == log_level::visu_detailed)
            {
                step_visu(iter, u);
            }
        }
    }
    fsi_eq_ptr->getTheta() = theta;

    // Time stepping loop over subinterval
    for (auto iter = implicit_steps; iter < n_steps; ++iter)
    {
        // actual solving
        step<method>(time, dt, iter, theta, u, f);
        if constexpr (logging == log_level::visu_detailed)
        {
            step_visu(iter, u);
        }
    }

    auto end = omp_get_wtime();
    if constexpr (method == iter_type::fine_last || method == iter_type::fine)
    {
#pragma omp atomic update
        fine_exec_time += (end - start);
    }
    else if constexpr (method == iter_type::coarse || method == iter_type::coarse_first)
    {
#pragma omp atomic update
        coarse_exec_time += (end - start);
    }
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
template <iter_type method>
void parareal<DIM, logging>::step(double& time, const double& dt, const int& iter,
                                  const double& theta, VectorInterface& u, VectorInterface& f)
{
    if constexpr (logging == log_level::visu_info || logging == log_level::visu_detailed)
    {
        if constexpr (method == iter_type::fine || method == iter_type::fine_last)
        {
            std::cerr << logwrite::info("Fine", sty::bb, time, " -> ", time + dt, " [", theta,
                                        " ]");
        }
        else if constexpr (method == iter_type::coarse || method == iter_type::coarse_first)
        {
            std::cerr << logwrite::info("Coarse", sty::bb, time, " -> ", time + dt, " [", theta,
                                        " ]");
        }
    }

    time += dt;
    set_time(dt, time);
    GetMultiLevelSolver()->Equ(old, 1.0, u);
    GetMultiLevelSolver()->AddNodeVector("old", old);
    assert(StdLoop::Solve(u, f) == "converged");

    if constexpr (method == iter_type::fine_last)
    {
        log_buffer[iter][0] = time;
        DoubleVector juh    = StdLoop::Functionals(u, f);
        for (short i = 0; i < 4; ++i)
        {
            log_buffer[iter][i + 1] = juh[i];
        }
    }
    GetMultiLevelSolver()->DeleteNodeVector("old");
}

//--------------------------------------------------------------------------------------------------
// TODO
template <int DIM, log_level logging>
void parareal<DIM, logging>::divergence_stab(VectorInterface& u, VectorInterface& f,
                                             const double& time, const double& dt,
                                             const int& current)
{
    cerr << logwrite::info("Info", "Problem for Stabilization");
    visu("Results/before_stab", u, 10 * current + subinterval_idx);

    GetMultiLevelSolver()->SetProblem("div");
    set_time(dt, time);

    StdLoop::GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    GetMultiLevelSolver()->Equ(tmp_vi, 1.0, u);
    GetMultiLevelSolver()->AddNodeVector("old", tmp_vi);
    assert(StdLoop::Solve(u, f) == "converged");
    GetMultiLevelSolver()->DeleteNodeVector("old");
    GetSolver()->Add(tmp_vi, -1, u);
    visu(dir_vec{"Results/difference", tmp_vi, 10 * current + subinterval_idx},
         dir_vec{"Results/after_stab", u, 10 * current + subinterval_idx});

    GetMultiLevelSolver()->SetProblem("fsi");
    StdLoop::GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::visu_write(const std::string& dir, const VectorInterface& vec,
                                        size_t index) const
{
    GetSolver()->Visu(dir, vec, index);
    StdLoop::WriteMeshAndSolution(dir, vec);
}

template <int DIM, log_level logging>
void parareal<DIM, logging>::visu(const std::string& dir, const VectorInterface& vec,
                                  size_t index) const
{
    GetSolver()->Visu(dir, vec, index);
    GetSolver()->Write(vec, dir + '.' + to_string(index));
}

template <int DIM, log_level logging>
void parareal<DIM, logging>::visu(const dir_vec& dv) const
{
    GetSolver()->Visu(dv.dir, dv.vec, dv.index);
    GetSolver()->Write(dv.vec, dv.dir + '.' + to_string(dv.index));
}

template <int DIM, log_level logging>
template <typename... dir_vec_pack>
void parareal<DIM, logging>::visu(const dir_vec& dv, const dir_vec_pack&... dvs) const
{
    visu(dv);
    visu(dvs...);
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::set_time(const double dt, const double time) const
{
    for (auto i = 0; i < GetMultiLevelSolver()->nlevels(); ++i)
    {
        GetMultiLevelSolver()->GetSolver(i)->GetProblemDescriptor()->SetTime(time, dt);
    }
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::setGV(VectorInterface& v, const GlobalVector& v_ptr) const
{
    GetSolver()->GetGV(v) = v_ptr;
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::setGVsize(GlobalVector& v_ptr)
{
    auto ncomp    = GetSolver()->GetProblemDescriptor()->GetEquation()->GetNcomp();
    auto n        = GetSolver()->GetDiscretization()->n();
    v_ptr.ncomp() = ncomp;
    v_ptr.resize(n);
    // v_ptr.zero();
}

template <int DIM, log_level logging>
void parareal<DIM, logging>::setGVsize(const VectorInterface& vec)
{
    auto v_ptr = GetSolver()->GetGV(vec);
    setGVsize(v_ptr);
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::setup_problem(const ParamFile& paramfile,
                                           const std::string& problemlabel)
{
    auto Problem2d = new ProblemDescriptor2d();
    Problem2d->BasicInit(&paramfile);
    auto Div2d = new DivergenceDescriptor2d();
    Div2d->BasicInit(&paramfile);

    auto PC2d = new ProblemContainer();
    PC2d->AddProblem(problemlabel, Problem2d);
    PC2d->AddProblem("div", Div2d);

    auto FC2d = new FunctionalContainer();
    auto Px   = new WeightedPointFunctional();
    auto Py   = new WeightedPointFunctional();

    vector<Vertex2d> v2d;
    v2d.push_back(Vertex2d(0.6, 0.2));
    vector<int> cx, cy;
    cx.push_back(3);
    cy.push_back(4);
    vector<double> weigh;
    weigh.push_back(1.0);
    Px->BasicInit(v2d, cx, weigh);
    Py->BasicInit(v2d, cy, weigh);

    auto drag = new Drag();
    auto lift = new Lift();

    FC2d->AddFunctional("ux", Px);
    FC2d->AddFunctional("uy", Py);
    FC2d->AddFunctional("drag", drag);
    FC2d->AddFunctional("lift", lift);

    this->Loop<DIM>::BasicInit(&paramfile, PC2d, FC2d);
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::set_interval_threads(int& n_intervals, int& threads)
{
    if (n_intervals == 0 && threads == 0)
    {
        n_intervals = omp_get_num_threads();
        if (n_intervals < 2)
        {
            n_intervals = omp_get_max_threads();
        }
        else if (!omp_in_parallel())
        {
            std::cerr << "Number of Threads already set, so we are using this" << '\n';
        }
        else if (omp_in_parallel())
        {
            std::cerr << "At the moment calling parareal inside parallel "
                         "region is not guaranteed to work, will set it anyways."
                      << '\n';
            n_intervals = omp_get_max_threads();
        }
        assert(n_intervals > 1 && "Only 1 thread available");
        threads = n_intervals;
    }
    else
    {
        n_intervals = (n_intervals == 0) ? threads : n_intervals;
        threads     = (threads == 0) ? n_intervals : threads;
        assert(threads <= n_intervals && "More threads than subintervals make no sense");
    }
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
template <typename T>
const T parareal<DIM, logging>::get_param(const std::string& paramName, const std::string& block)
{
    T param;
    DataFormatHandler DFH;
    DFH.insert(paramName, &param, static_cast<T>(0));
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(&paramfile, block);
    return param;
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
template <name_type name>
const std::string parareal<DIM, logging>::get_name(double dtc, double tc, int max_iter)
{
    std::string dtc_str = std::to_string(dtc);
    if (dtc < 1.)
    {
        dtc_str = dtc_str.substr(dtc_str.find(".") + 1);
        dtc_str = dtc_str.substr(0, dtc_str.find_last_not_of("0") + 1);
    }
    else
    {
        dtc_str = dtc_str.substr(0, dtc_str.find_last_not_of("0") + 1);
    }
    std::string tc_str = std::to_string(tc);
    if (tc < 1.)
    {
        tc_str = tc_str.substr(tc_str.find(".") + 1);
        tc_str = tc_str.substr(0, tc_str.find_last_not_of("0") + 1);
    }
    else
    {
        tc_str = tc_str.substr(0, tc_str.find_last_not_of("0") + 1);
    }

    if constexpr (name == name_type::error)
    {
        return "para_errorC" + dtc_str + "tc" + tc_str + "I" + std::to_string(max_iter) + ".txt";
    }
    else if constexpr (name == name_type::functional)
    {
        return "para_functionalC" + dtc_str + "tc" + tc_str + "I" + std::to_string(max_iter)
               + ".txt";
    }
}

/**************************************************************************************************
Initialize static members
***************************************************************************************************
Setting up parareal:
***************************************************************************************************
     n_intervals        - Number of intervals for parareal
     max_iterations      - don't do more than max_iterations iterations
     dtfine             - fine step size
     n_finesteps        - number of fine timesteps per interval
     dtcoarse           - coarse step size
     n_coarsesteps      - number of coarse timesteps per interval
     interval_length     - length of intervals
**************************************************************************************************/

template <int DIM, log_level logging>
ParamFile parareal<DIM, logging>::paramfile = ParamFile("fsi-3.param");
template <int DIM, log_level logging>
std::string parareal<DIM, logging>::problemlabel = "fsi";

template <int DIM, log_level logging>
double parareal<DIM, logging>::fine_theta = parareal::get_param<double>("theta", "Equation");
template <int DIM, log_level logging>
double parareal<DIM, logging>::coarse_theta = parareal::get_param<double>("thetacoarse",
                                                                          "Equation");

template <int DIM, log_level logging>
double parareal<DIM, logging>::precompute_time = parareal::get_param<double>("precompute",
                                                                             "Equation");
template <int DIM, log_level logging>
int parareal<DIM, logging>::n_intervals = parareal::get_param<double>("subintervals", "Loop");
template <int DIM, log_level logging>
int parareal<DIM, logging>::threads = parareal::get_param<double>("threads", "Loop");
template <int DIM, log_level logging>
int parareal<DIM, logging>::max_iterations = static_cast<int>(1);
template <int DIM, log_level logging>
double parareal<DIM, logging>::start_time = parareal::get_param<double>("start_time", "Equation");
template <int DIM, log_level logging>
double parareal<DIM, logging>::stop_time = parareal::get_param<double>("stop_time", "Equation");

template <int DIM, log_level logging>
double parareal<DIM, logging>::interval_length =
  static_cast<double>((stop_time - start_time - precompute_time) / n_intervals);

template <int DIM, log_level logging>
double parareal<DIM, logging>::dtfine = parareal::get_param<double>("dt", "Equation");
template <int DIM, log_level logging>
size_t parareal<DIM, logging>::n_finesteps = static_cast<size_t>(interval_length / dtfine);

template <int DIM, log_level logging>
double parareal<DIM, logging>::dtcoarse = parareal::get_param<double>("dtcoarse", "Equation");
template <int DIM, log_level logging>
size_t parareal<DIM, logging>::n_coarsesteps = static_cast<size_t>(parareal::interval_length
                                                                   / parareal::dtcoarse);
template <int DIM, log_level logging>
double parareal<DIM, logging>::coarse_exec_time = 0;

template <int DIM, log_level logging>
double parareal<DIM, logging>::fine_exec_time = 0;

template <int DIM, log_level logging>
std::vector<GlobalVector> parareal<DIM, logging>::u_sol_arr = [] {
    std::vector<GlobalVector> tGridVec_ptr(parareal::n_intervals);
    tGridVec_ptr.reserve(parareal::n_intervals);
    for (auto i = 0; i < parareal::n_intervals; ++i)
    {
        tGridVec_ptr[i] = GlobalVector();
    }
    return tGridVec_ptr;
}();

template <int DIM, log_level logging>
std::ofstream parareal<DIM, logging>::func_log;
}  // namespace Gascoigne
#endif  // PARAREAL_H
