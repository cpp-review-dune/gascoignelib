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
enum name_type
{
    error,
    functional
};
enum iter_type
{
    fine_last,
    fine,
    coarse
};
enum log_level
{
    results       = 0,
    visu_info     = 1,
    visu_detailed = 2
};

namespace Gascoigne
{
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
 * compareParaSerial, everything else is handled automatically by the class. For better
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
     * \param maxIterations maximal number of iterations executed by parareal algorithm
     * \param coarse_theta \f$\theta\f$ for coarse time stepping
     * \param dtcoarse coarse time step size
     * \param theta \f$\theta\f$ for fine time stepping
     * \param dtfine fine time step size
     */
    static double runPara(const int maxIterations   = 1,
                          const double coarse_theta = getParam<double>("thetacoarse", "Equation"),
                          const double dtcoarse     = getParam<double>("dtcoarse", "Equation"),
                          const double fine_theta   = getParam<double>("theta", "Equation"),
                          const double dtfine       = getParam<double>("dt", "Equation"));

    /*! \brief Method comparing parareal with sequential execution.
     *
     * Executes the parareal algorithm with the specified coarse time step sizes and coarse
     * thetas, values for fine method is taken from the paramfile. Before execution of parareal
     * it solves the problem sequentially and calculates the norm \f$l_2\f$ vector norm between
     * the two solutions on every time node. \param maxIterations maximal number of iterations
     * executed by parareal algorithm \param dtcoarse_vec Vector with coarse time step sizes
     * that will be tested against sequential execution \param coarse_theta_vec Vector with
     * coarse theta values for step sizes specified by dtcoarse_vec
     */
    static void compareParaSerial(const unsigned int maxIterations,
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
     * \brief Wrapper function for GetSolver calls via GetMultiLevelSolver returning *nonconst*
     * pointer to Solver
     */
    SolverInterface* GetSolver();

    /*!
     * \brief Wrapper function for GetMultiLevelSolver returning *const* pointer to
     * MultiLevelSolver
     */
    const MultiLevelSolverInterface* GetMultiLevelSolver() const;

    /*!
     * \brief Wrapper function for GetMultiLevelSolver returning *nonconst* pointer to
     * MultiLevelSolver
     */
    MultiLevelSolverInterface* GetMultiLevelSolver();

private:
    /*! \brief Method executing parareal-algorithm.
     *
     * This is the method that is called after runPara has set all parameters for
     * parareal.
     */
    static double paraAlgo();

    /*!
     * \brief Time stepping method.
     *
     * This method encapsulates *noFineSteps* or *noCoarseSteps* timesteps into one call.
     * The step size used is dtfine or *dtcoarse* (depending on whether coarse or fine is
     * specified by the template parameter). This method writes a logfile if method equals
     * fine_last. \param time          current time from where the steps are executed \param u
     * VectorInterface with approximation at current time \param f             Right hand inside
     */
    template <iter_type method>
    void propagator(double time, VectorInterface& u, VectorInterface& f, const int current = 0,
                    const unsigned implicit_steps = 0);

    template <iter_type method>
    void inline step(double& time, const double& dt, const int& iter, VectorInterface& u,
                     VectorInterface& f);

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
    void SetTime(const double dt, const double time = 0.) const;

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

    /*!
     *\brief Setting up the problem.
     *
     * This method initializes:
     * - ProblemDescriptor
     * - ProblemContainer
     * - etc...
     */
    void setupProblem(const ParamFile& paramfile, const std::string& problemlabel);

    /*!
     * \brief Initialization method for setting the number of intervals.
     *
     * The method sets the number of the subintervals to an already defined number of threads
     * (with `omp_get_thread_num()`), or to the maximum number of threads available with
     * `omp_get_max_threads()`. Setting the number of subintervals to the maximum is good for
     * using the maximum potential of parallelization. If threads are defined in the paramfile,
     * they are used.
     */
    static int setNoIntervals();

    /*!
     * \brief Method for reading parameters from the paramfile.
     *
     * Before calling this the `paramfile` has to be defined already.
     * \param paramName Name of parameter which will be searched for
     * \param block Block in which paramName is searched
     */
    template <typename T>
    static const T getParam(const std::string& paramName, const std::string& block);

    template <name_type name>
    static const std::string getName(double dtc = dtcoarse, double tc = coarse_theta,
                                     unsigned maxIter = maxIterations);

    std::size_t subinterval_idx;  //!< Identifies the subinterval on which the loop is acting in the
                                  //!< first iteration
    VectorInterface
      fi_solution;  //!< Vector containing the fine solution on the corresponding interval.
    VectorInterface
      co_solution;        //!< Vector containing the coarse solution on the corresponding interval.
    VectorInterface old;  //!< Vector for temporary results.
    VectorInterface
      u_solution;  //!< Vector containing the final solution on the corresponding interval.
    std::vector<DoubleVector> log_buffer;  //!< Vector for buffering log entries.

    const int c_implicit_steps = getParam<int>("implicit coarse", "Equation");
    const int f_implicit_steps = getParam<int>("implicit fine", "Equation");

    static double fine_theta;    //!< Theta for fine timestepping
    static double coarse_theta;  //!< Theta for coarse timestepping

    static double start_time;           //!< Start time point
    static double stop_time;            //!< End time point
    static unsigned int noIntervals;    //!< Number of intervals parareal is using. Usually this
                                        //!< is the maximum number of threads
    static unsigned int maxIterations;  //!< Maximum number of iterations parareal should execute.
    static double intervalLength;       //!< Length of one subinterval parareal uses

    static double dtfine;            //!< Time step size of the fine propagator
    static std::size_t noFineSteps;  //!< Number of fine timesteps per subinterval

    static double dtcoarse;            //!< Time step size of the coarse propagator
    static std::size_t noCoarseSteps;  //!< Number of coarse timesteps per subinterval

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
    , fi_solution("f" + std::to_string(id))
    , co_solution("g" + std::to_string(id))
    , old("old" + std::to_string(id))
    , u_solution("u" + std::to_string(id))
    , log_buffer(noFineSteps, DoubleVector(5, 0.))
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
          << "\nNumber of threads is                   " << noIntervals
          << "\nInterval length                        " << intervalLength
          << "\nMaximal number of iterations           " << maxIterations
          << "\nParameter for coarse method:\
              \nCoarse steps per interval:             " << noCoarseSteps
          << "\nCoarse step size:                      " << dtcoarse
          << "\nTheta for coarse steps                 " << coarse_theta
          << "\nParameter for fine method:\
              \nFine steps per interval:               " << noFineSteps
          << "\nFine step size:                        " << dtfine
          << "\nTheta for fine steps                   " << fine_theta
          << "\nUsing paramfile:                       " << paramfile
          << "\nProblemlabel is:                       " << problemlabel
          << "\n=========================================================\n\n";
    // clang-format on
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
double parareal<DIM, logging>::runPara(const int maxIterations, const double coarse_theta,
                                       const double dtcoarse, const double fine_theta,
                                       const double dtfine)
{
    assert(dtcoarse >= dtfine && "coarse step size should be bigger than fine step size");
    parareal::maxIterations = maxIterations;
    parareal::coarse_theta  = coarse_theta;
    parareal::dtcoarse      = dtcoarse;
    parareal::fine_theta    = fine_theta;
    parareal::dtfine        = dtfine;
    parareal::noCoarseSteps = static_cast<size_t>(intervalLength / dtcoarse);
    parareal::noFineSteps   = static_cast<size_t>(intervalLength / dtfine);

    auto ratio = static_cast<int>(dtcoarse / dtfine);
    if (ratio * noCoarseSteps != noFineSteps)
    {
        noFineSteps += noFineSteps % ratio;
        noCoarseSteps  = static_cast<int>(noFineSteps / ratio);
        intervalLength = noFineSteps * dtfine;
    }

    std::cerr << style::bb << style::g << "[Info] " << style::n << "Fine Visualization\n";
    outputParameters();
    func_log.precision(12);
    // Checking if input is correct
    assert(isfinite(dtcoarse) && isfinite(dtfine) && "At least one timestep is infinite.");
    assert(isfinite(noCoarseSteps) && isfinite(noFineSteps) && "Amount of steps is infinite.");
    assert(isfinite(intervalLength) && isfinite(start_time)
           && "intervalLength or start time is infinite");
    assert(stop_time - start_time > 0 && "start has to be be before stop");

    assert(dtcoarse > 0 && dtfine > 0 && "At least one timestep is less or equal 0.");
    assert(noCoarseSteps != 0 && noFineSteps != 0 && "Amount of steps is 0.");
    assert(noIntervals > 1
           && "Number of intervals is 1 or 0. You should prefer sequential execution.");
    assert(intervalLength > 0 && "intervalLength should be greater zero!");
    assert(fine_theta > 0 && coarse_theta > 0 && "Thetas have to be positive");
    assert(dtcoarse * noCoarseSteps - dtfine * noFineSteps <= numeric_limits<double>::epsilon()
           && "Coarse and fine step sizes and number of time steps don't match");

    if (!func_log.is_open())
    {
        func_log.open(getName<functional>());
    }
    auto exec_time = paraAlgo();

    std::cerr << style::bb << "\n==================================================\n"
              << style::G << "Finished parareal" << style::b << ".\nElapsed time is\t" << exec_time
              << "\n\nTime spent in coarse method: " << coarse_exec_time
              << "\nTime spent in fine method: " << fine_exec_time << style::n << '\n';
    func_log.close();
    return exec_time;
}

//--------------------------------------------------------------------------------------------------

#ifdef _OPENMP  // parareal can't be run without OpenMP
template <int DIM, log_level logging>
double parareal<DIM, logging>::paraAlgo()
{
    double time     = start_time;
    auto start_wall = omp_get_wtime();
    //
    // Initialization of solvers and problems and solutions from parareal
    // co_solution   - VectorInterface for results of coarse method
    // fi_solution   - VectorInterface for results of fine method
    // u_solution    - VectorInterface for final results
    // old           - VectorInterface for temporary results
    // u             - VectorInterface used for initialization of the others
    //
    VectorInterface u("u");
    VectorInterface f("f");
    std::vector<parareal*> loop_interval(noIntervals);
    for (auto i = 0; i < noIntervals; ++i)
    {
        loop_interval[i] = new parareal<DIM, logging>(i + 1);
        loop_interval[i]->setupProblem(paramfile, problemlabel);

        loop_interval[i]->GetMultiLevelSolver()->ReInit(problemlabel);
        loop_interval[i]->GetSolver()->OutputSettings();
        cout << loop_interval[i]->GetMeshAgent()->ncells() << " cells" << '\n';

        loop_interval[i]->GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
        loop_interval[i]->SetTime(dtfine, start_time);
    }

    for (auto m = 0; m < noIntervals; ++m)
    {
        loop_interval[m]->GetMultiLevelSolver()->ReInitVector(f);
        loop_interval[m]->GetMultiLevelSolver()->ReInitVector(u);
        loop_interval[m]->GetMultiLevelSolver()->ReInitVector(loop_interval[m]->fi_solution);
        loop_interval[m]->GetMultiLevelSolver()->ReInitVector(loop_interval[m]->co_solution);
        loop_interval[m]->GetMultiLevelSolver()->ReInitVector(loop_interval[m]->old);
        loop_interval[m]->GetMultiLevelSolver()->ReInitVector(loop_interval[m]->u_solution);
        loop_interval[m]->setGVsize(
          loop_interval[m]->GetSolver()->GetGV(loop_interval[m]->fi_solution));
        loop_interval[m]->setGVsize(
          loop_interval[m]->GetSolver()->GetGV(loop_interval[m]->co_solution));
        loop_interval[m]->setGVsize(loop_interval[m]->GetSolver()->GetGV(loop_interval[m]->old));
        loop_interval[m]->setGVsize(
          loop_interval[m]->GetSolver()->GetGV(loop_interval[m]->u_solution));
    }

    for (auto m = 0; m < noIntervals; ++m)
    {
        loop_interval[m]->InitSolution(u);
        loop_interval[m]->InitSolution(loop_interval[m]->fi_solution);
        loop_interval[m]->InitSolution(loop_interval[m]->co_solution);
        // loop_interval[m]->InitSolution(loop_interval[m]->old);
        loop_interval[m]->InitSolution(loop_interval[m]->u_solution);
    }

    for (auto m = 0; m < noIntervals; ++m)
    {
        loop_interval[m]->setGVsize(u_sol_arr[m]);
    }

    /// Initializing locks
    auto interval_locker = new omp_lock_t[noIntervals];
    for (auto m = 0; m < noIntervals; ++m)
    {
        omp_init_lock(&interval_locker[m]);
    }

    //
    // A C T U A L   A L G O R I T H M
    //
#pragma omp parallel num_threads(noIntervals) firstprivate(time) proc_bind(close)
    {
        //
        // 'Iteration 0'
        //
#pragma omp for ordered nowait schedule(monotonic : static)
        // Initialization
        for (auto i = 0; i < noIntervals; ++i)
        {
#pragma omp ordered
            {
                omp_set_lock(&interval_locker[i]);
                // u ⟵ G⁰(i)
                loop_interval[i]->template propagator<coarse>(time + i * intervalLength, u, f, 0,
                                                              loop_interval[i]->c_implicit_steps);
                // G⁰(i) ⟵ u;
                loop_interval[i]->GetSolver()->Equ(loop_interval[i]->co_solution, 1., u);
                // U⁰(i) ⟵ u;
                loop_interval[i]->GetSolver()->Equ(loop_interval[i]->u_solution, 1., u);

                // Visualization for debugging
                if constexpr (logging > 0)
                {
                    std::cerr << style::bb << style::g << "\n[Done] " << style::n
                              << " Initialization on subinterval " << i << '\n';
                    loop_interval[i]->GetSolver()->Visu("Results/initsol",
                                                        loop_interval[i]->u_solution, i);
                    loop_interval[i]->WriteMeshAndSolution("Results/initsol",
                                                           loop_interval[i]->u_solution);
                }

                if (i < noIntervals - 1)
                {
                    // pass GlobalVector around
                    loop_interval[i + 1]->setGV(u, loop_interval[i]->GetSolver()->GetGV(u));
                    // F⁰(i) ⟵ u
                    loop_interval[i + 1]->GetSolver()->Equ(loop_interval[i + 1]->fi_solution, 1.,
                                                           u);
                }

                omp_unset_lock(&interval_locker[i]);
            }
        }  // Initialization

        //#pragma omp barrier  // barrier for debugging
        //
        // Do <maxIterations> parareal-iterations
        // Iterations
        for (auto k = 1; k <= maxIterations; ++k)
        {
            // fine propagations
#pragma omp for nowait schedule(monotonic : static)
            for (auto m = 0; m < noIntervals - k + 1; ++m)
            {
                // #pragma omp ordered  // for debugging
                //                 {
                omp_set_lock(&interval_locker[m]);
                // Fᵏ(m) ⟵ F(Uᵏ⁻¹(m-1))
                if (m == 0 || k == maxIterations)
                {
                    loop_interval[m]->template propagator<fine_last>(
                      time + m * intervalLength, loop_interval[m]->fi_solution, f, k,
                      loop_interval[m]->f_implicit_steps);
                }
                else  // we are not in the last iteration and not exact
                {
                    loop_interval[m]->template propagator<fine>(time + m * intervalLength,
                                                                loop_interval[m]->fi_solution, f, k,
                                                                loop_interval[m]->f_implicit_steps);
                }
                if constexpr (logging > 0)
                {
                    std::cerr << style::bb << style::g << "\n[Done] " << style::n
                              << " Fine solution on subinterval " << m + k - 1
                              << " in iteration number " << k << '\n';
                    loop_interval[m]->GetSolver()->Visu(
                      "Results/fineres", loop_interval[m]->fi_solution, (k - 1) * 10 + m);
                    loop_interval[m]->WriteMeshAndSolution("Results/fineres",
                                                           loop_interval[m]->fi_solution);
                }
                omp_unset_lock(&interval_locker[m]);
                //                }
            }
            // barrier for debugging
            //#pragma omp barrier
            //
            // coarse propagations and corrections
            //
#pragma omp for ordered nowait schedule(monotonic : static)
            // Corrections
            for (auto m = 0; m < noIntervals - k; ++m)
            {
#pragma omp ordered
                {
                    // k-th value is exact in k-th Iteration
                    if (m == 0)
                    {
                        omp_set_lock(&interval_locker[0]);
                        // Uᵏ(k) ⟵ Fᵏ(k)
                        loop_interval[0]
                          ->GetSolver()
                          ->GetGV(loop_interval[0]->u_solution)
                          .equ(1.,
                               loop_interval[0]->GetSolver()->GetGV(loop_interval[0]->fi_solution));
                        u_sol_arr[k - 1].equ(
                          1., loop_interval[0]->GetSolver()->GetGV(loop_interval[0]->fi_solution));
                        for (const auto& x : loop_interval[0]->log_buffer)
                        {
                            func_log << x << '\n';
                        }
                        loop_interval[0]->GetSolver()->Visu("Results/p",
                                                            loop_interval[0]->u_solution, k - 1);
                        loop_interval[0]->WriteMeshAndSolution("Results/p",
                                                               loop_interval[0]->u_solution);
                        omp_unset_lock(&interval_locker[0]);
                    }

                    // prepare predictor step
                    omp_set_lock(&interval_locker[m]);
                    loop_interval[m]
                      ->GetSolver()
                      ->GetGV(loop_interval[m]->co_solution)
                      .equ(1., loop_interval[m]->GetSolver()->GetGV(loop_interval[m]->u_solution));
                    omp_unset_lock(&interval_locker[m]);
                    // predictor step with coarse method
                    loop_interval[m]->template propagator<coarse>(
                      time + (m + 1) * intervalLength, loop_interval[m]->co_solution, f, k,
                      loop_interval[m]->c_implicit_steps);

                    // omp_unset_lock(&interval_locker[m]);
                    omp_set_lock(&interval_locker[m + 1]);

                    // calculate corrector term
                    // clang-format off
                    loop_interval[m + 1]
                      ->GetSolver()
                      ->GetGV(loop_interval[m + 1]->u_solution)
                      .equ(
                        1,  loop_interval[m]->GetSolver()->GetGV(loop_interval[m]->co_solution),
                        1,  loop_interval[m + 1]->GetSolver()->GetGV(
                                loop_interval[m + 1]->fi_solution),
                        -1, loop_interval[m + 1]->GetSolver()->GetGV(
                                loop_interval[m + 1]->co_solution));
                    // clang-format on

                    // prepare next iteration
                    loop_interval[m + 1]
                      ->GetSolver()
                      ->GetGV(loop_interval[m + 1]->fi_solution)
                      .equ(1., loop_interval[m + 1]->GetSolver()->GetGV(
                                 loop_interval[m + 1]->u_solution));
                    omp_unset_lock(&interval_locker[m + 1]);

                    // Visualization for debugging
                    if constexpr (logging > 0)
                    {
                        std::cerr << style::bb << style::g << "\n[Done] " << style::n
                                  << " Coarse solution on subinterval " << m + 1
                                  << " in iteration number " << k << '\n';
                        loop_interval[m]->GetSolver()->Visu(
                          "Results/coarseres", loop_interval[m]->co_solution, (k - 1) * 10 + m);
                        loop_interval[m]->WriteMeshAndSolution("Results/coarseres",
                                                               loop_interval[m]->co_solution);
                        std::cerr << style::bb << style::g << "\n[Done] " << style::n
                                  << " initial value on subinterval " << m + 1
                                  << " in iteration number " << k << '\n';
                        loop_interval[m]->GetSolver()->Visu(
                          "Results/iterres", loop_interval[m]->u_solution, (k - 1) * 10 + m);
                        loop_interval[m]->WriteMeshAndSolution("Results/iterres",
                                                               loop_interval[m]->u_solution);
                        std::cerr << style::bb << style::g << "[Done] " << style::b
                                  << " Iteration number " << k << " of parareal " << style::n
                                  << '\n';
                    }
                }

            }  // Corrections

#pragma omp atomic update  // increase time
            time += intervalLength;
        }  // Iterations
    }      // parallel region
    auto exec_time = omp_get_wtime() - start_wall;

    for (auto m = 0; m < noIntervals - maxIterations; ++m)
    {
        // write final solution vectors
        loop_interval[m]->setGVsize(u_sol_arr[m]);
        // final solutions on time nodes
        u_sol_arr[m + maxIterations].equ(
          1., loop_interval[m + 1]->GetSolver()->GetGV(loop_interval[m + 1]->u_solution));
        // logging
        for (const auto& x : loop_interval[m + 1]->log_buffer)
        {
            func_log << x << '\n';
        }

        // Visu part
        // loop_interval[m]->setGV(u, u_sol_arr[m]);
        loop_interval[m + 1]->GetSolver()->Visu("Results/p", loop_interval[m + 1]->u_solution,
                                                m + maxIterations);
        loop_interval[m + 1]->WriteMeshAndSolution("Results/p", loop_interval[m + 1]->u_solution);
        std::cerr << style::bb << style::g << "[Done] " << style::n << "Wrote subinterval "
                  << m + maxIterations << '\n';
    }

    // Clean up
    omp_destroy_lock(interval_locker);
    delete[] interval_locker;
    for (auto i = 0; i < noIntervals; ++i)
    {
        delete loop_interval[i];
    }
    return exec_time;
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::compareParaSerial(const unsigned int maxIterations,
                                               const std::vector<double>& dtcoarse_vec,
                                               const std::vector<double>& coarse_theta_vec)
{
    // ParamFile paramfile(paramstr);
    if constexpr (logging == results)
    {
        ofstream timings_log({"timingsI" + std::to_string(maxIterations) + ".txt"});
        timings_log << "ratio time accel\n";

        // serial part
        auto start   = omp_get_wtime();
        auto loopSeq = new parareal<DIM, logging>(0);
        loopSeq->setupProblem(paramfile, problemlabel);
        loopSeq->GetMultiLevelSolver()->ReInit(problemlabel);
        loopSeq->GetSolver()->OutputSettings();
        cout << loopSeq->GetMeshAgent()->ncells() << " cells" << '\n';
        VectorInterface f("f");
        double time = start_time;
        std::vector<VectorInterface> U_seq(noIntervals, loopSeq->old);
        for (auto i = 0; i < noIntervals; ++i)
        {
            U_seq[i].SetName({"useq" + std::to_string(i)});
            U_seq[i].SetType("node");
            loopSeq->GetSolver()->ReInitVector(U_seq[i]);
        }
        for (auto i = 0; i < noIntervals; ++i)
        {
            loopSeq->InitSolution(U_seq[i]);
        }

        loopSeq->GetMultiLevelSolver()->ReInitVector(loopSeq->old);
        loopSeq->GetMultiLevelSolver()->ReInitVector(loopSeq->u_solution);
        loopSeq->GetMultiLevelSolver()->ReInitVector(f);
        loopSeq->InitSolution(loopSeq->u_solution);
        loopSeq->InitSolution(loopSeq->old);

        if (func_log.is_open())
        {
            func_log.close();
        }
        func_log.open(getName<functional>(dtcoarse, coarse_theta, 0));
        // sequential run
        for (auto i = 0; i < noIntervals; ++i)
        {
            loopSeq->template propagator<fine_last>(time + i * intervalLength, loopSeq->u_solution,
                                                    f);
            loopSeq->GetSolver()->Equ(U_seq[i], 1., loopSeq->u_solution);
            // time += intervalLength;
            for (const auto& x : loopSeq->log_buffer)
            {
                func_log << x << '\n';
            }
            loopSeq->GetSolver()->Visu("Results/u", loopSeq->u_solution, i);
            loopSeq->WriteMeshAndSolution("Results/u", loopSeq->u_solution);
        }
        func_log.close();

        auto end         = omp_get_wtime();
        auto serial_time = end - start;
        std::cerr << "\n==================================================\n"
                     "Finished sequential execution.\nElapsed time is\t"
                  << serial_time << ".\n";
        timings_log << 1 << ' ' << serial_time << ' ' << 1 << '\n';

        assert(
          dtcoarse_vec.size() == coarse_theta_vec.size()
          && "For every dtcoarse there has to be one coarse_theta given and the other way round");

        ofstream error_log;

        for (auto i = 0; i < dtcoarse_vec.size(); ++i)
        {
            func_log.open(getName<functional>(dtcoarse_vec[i], coarse_theta_vec[i]));

            // run Parareal
            auto ratio           = static_cast<int>(dtcoarse_vec[i] / dtfine);
            double parallel_time = runPara(maxIterations, coarse_theta_vec[i], dtcoarse_vec[i]);
            timings_log << ratio << ' ' << parallel_time << ' ' << serial_time / parallel_time
                        << '\n';
            func_log.close();

            // calculate error
            error_log.open(getName<error>(dtcoarse_vec[i], coarse_theta_vec[i]));
            error_log << "interval error\n";
            for (auto i = 0; i < noIntervals; ++i)
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
void parareal<DIM, logging>::compareParaSerial(const unsigned int maxIterations,
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
                                        int current, unsigned implicit_steps)
{
    double dt;
    double theta;
    size_t noSteps;
    auto start = omp_get_wtime();
    // Some compile-time branching...
    if constexpr (method == fine_last || method == fine)
    {
        dt      = dtfine;
        theta   = fine_theta;
        noSteps = noFineSteps;
    }
    else if constexpr (method == coarse)
    {
        dt      = dtcoarse;
        theta   = coarse_theta;
        noSteps = noCoarseSteps;
    }

    auto fsi_eq_ptr =
      dynamic_cast<const FSI<DIM>* const>(GetSolver()->GetProblemDescriptor()->GetEquation());
    if (implicit_steps > 0)
    {
        fsi_eq_ptr->getTheta() = 1.;
        for (auto iter = 0; iter < implicit_steps; iter++)
        {
            step<method>(time, dt, iter, u, f);
        }
    }
    fsi_eq_ptr->getTheta() = theta;

    // Time stepping loop over subinterval
    for (auto iter = implicit_steps; iter < noSteps; ++iter)
    {
        if constexpr (logging > 0)
        {
            if constexpr (method == fine || method == fine_last)
            {
                std::cerr << style::bb << style::i << "Fine " << time << " -> " << time + dt
                          << style::n << "  [" << theta << "]\n";
            }
            else if constexpr (method == coarse)
            {
                std::cerr << style::bb << "Coarse " << time << " -> " << time + dt << style::n
                          << "  [" << theta << "]\n";
            }
        }

        // actual solving
        step<method>(time, dt, iter, u, f);

        if constexpr (logging == visu_detailed)
        {
            if constexpr (method == fine || method == fine_last)
            {
                if (iter % 10 == 0 && iter < 310)
                {
                    std::cerr << style::bb << style::g << "[Info] " << style::n
                              << "Fine Visualization\n";
                    auto step = static_cast<unsigned>(time / dtfine);
                    GetMultiLevelSolver()->GetSolver()->Visu(
                      "Results_detail/subinterval" + std::to_string(subinterval_idx + current - 1)
                        + "/iter" + std::to_string(current) + "/u_fine_no",
                      u, current * static_cast<unsigned>(pow(10, ceil(log10(step)))) + step);
                }
            }
            else if constexpr (method == coarse)
            {
                if (iter < 31)
                {
                    std::cerr << style::bb << style::g << "[Info] " << style::n
                              << "Coarse Visualization\n";
                    auto step = static_cast<unsigned>(time / dtcoarse);
                    GetMultiLevelSolver()->GetSolver()->Visu(
                      "Results_detail/subinterval" + std::to_string(subinterval_idx + current)
                        + "/iter" + std::to_string(current) + "/u_coarse_no",
                      u, current * static_cast<unsigned>(pow(10, ceil(log10(step)))) + step);
                }
            }
        }
    }

    auto end = omp_get_wtime();
    if constexpr (method == fine_last || method == fine)
    {
#pragma omp atomic update
        fine_exec_time += (end - start);
    }
    else if constexpr (method == coarse)
    {
#pragma omp atomic update
        coarse_exec_time += (end - start);
    }
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
template <iter_type method>
void parareal<DIM, logging>::step(double& time, const double& dt, const int& iter,
                                  VectorInterface& u, VectorInterface& f)
{
    time += dt;
    SetTime(dt, time);
    GetMultiLevelSolver()->Equ(old, 1.0, u);
    GetMultiLevelSolver()->AddNodeVector("old", old);
    assert(StdLoop::Solve(u, f) == "converged");

    if constexpr (method == fine_last)
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

template <int DIM, log_level logging>
void parareal<DIM, logging>::SetTime(const double dt, const double time) const
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

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
void parareal<DIM, logging>::setupProblem(const ParamFile& paramfile,
                                          const std::string& problemlabel)
{
    auto Problem2d = new ProblemDescriptor2d();
    Problem2d->BasicInit(&paramfile);

    auto PC2d = new ProblemContainer();
    PC2d->AddProblem(problemlabel, Problem2d);

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
int parareal<DIM, logging>::setNoIntervals()
{
    auto noIntervals = getParam<double>("threads", "Loop");
    if (noIntervals == 0)
    {
        noIntervals = omp_get_num_threads();
        if (noIntervals < 2)
        {
            noIntervals = omp_get_max_threads();
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
            noIntervals = omp_get_max_threads();
        }
        assert(noIntervals > 1 && "Only 1 thread available");
        return noIntervals;
    }
    else
    {
        return noIntervals;
    }
}

//--------------------------------------------------------------------------------------------------

template <int DIM, log_level logging>
template <typename T>
const T parareal<DIM, logging>::getParam(const std::string& paramName, const std::string& block)
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
const std::string parareal<DIM, logging>::getName(double dtc, double tc, unsigned maxIter)
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

    if constexpr (name == error)
    {
        return "para_errorC" + dtc_str + "tc" + tc_str + "I" + std::to_string(maxIter) + ".txt";
    }
    else if constexpr (name == functional)
    {
        return "para_functionalC" + dtc_str + "tc" + tc_str + "I" + std::to_string(maxIter)
               + ".txt";
    }
}

/**************************************************************************************************
Initialize static members
***************************************************************************************************
Setting up parareal:
***************************************************************************************************
     noIntervals        - Number of intervals for parareal
     maxIterations      - don't do more than maxIterations iterations
     dtfine             - fine step size
     noFineSteps        - number of fine timesteps per interval
     dtcoarse           - coarse step size
     noCoarseSteps      - number of coarse timesteps per interval
     intervalLength     - length of intervals
**************************************************************************************************/

template <int DIM, log_level logging>
ParamFile parareal<DIM, logging>::paramfile = ParamFile("fsi-3.param");
template <int DIM, log_level logging>
std::string parareal<DIM, logging>::problemlabel = "fsi";

template <int DIM, log_level logging>
double parareal<DIM, logging>::fine_theta = parareal::getParam<double>("theta", "Equation");
template <int DIM, log_level logging>
double parareal<DIM, logging>::coarse_theta = parareal::getParam<double>("thetacoarse", "Equation");

template <int DIM, log_level logging>
unsigned int parareal<DIM, logging>::noIntervals = static_cast<unsigned int>(setNoIntervals());
template <int DIM, log_level logging>
unsigned int parareal<DIM, logging>::maxIterations = static_cast<unsigned int>(1);
template <int DIM, log_level logging>
double parareal<DIM, logging>::start_time = parareal::getParam<double>("start_time", "Equation");
template <int DIM, log_level logging>
double parareal<DIM, logging>::stop_time = parareal::getParam<double>("stop_time", "Equation");

template <int DIM, log_level logging>
double parareal<DIM, logging>::intervalLength = static_cast<double>((stop_time - start_time)
                                                                    / noIntervals);

template <int DIM, log_level logging>
double parareal<DIM, logging>::dtfine = parareal::getParam<double>("dt", "Equation");
template <int DIM, log_level logging>
size_t parareal<DIM, logging>::noFineSteps = static_cast<size_t>(intervalLength / dtfine);

template <int DIM, log_level logging>
double parareal<DIM, logging>::dtcoarse = parareal::getParam<double>("dtcoarse", "Equation");
template <int DIM, log_level logging>
size_t parareal<DIM, logging>::noCoarseSteps = static_cast<size_t>(parareal::intervalLength
                                                                   / parareal::dtcoarse);
template <int DIM, log_level logging>
double parareal<DIM, logging>::coarse_exec_time = 0;

template <int DIM, log_level logging>
double parareal<DIM, logging>::fine_exec_time = 0;

template <int DIM, log_level logging>
std::vector<GlobalVector> parareal<DIM, logging>::u_sol_arr = [] {
    std::vector<GlobalVector> tGridVec_ptr(parareal::noIntervals);
    tGridVec_ptr.reserve(parareal::noIntervals);
    for (auto i = 0; i < parareal::noIntervals; ++i)
    {
        tGridVec_ptr[i] = GlobalVector();
    }
    return tGridVec_ptr;
}();

template <int DIM, log_level logging>
std::ofstream parareal<DIM, logging>::func_log;
}  // namespace Gascoigne
#endif  // PARAREAL_H
