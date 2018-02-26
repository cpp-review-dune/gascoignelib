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
 * correct definition of the problem, there is no checking.The implementation only executes the
 * bare algorithm. Checking for convergence is not available at the moment. This can be done via
 * comparing serial with parallel execution. Some properties
 * - Parallelization with OpenMP
 * - Inherits from MultiLevelAlgorithm
 * - Memory consumption is higher than sequential
 *
 * Every Instance of parareal correponds to one solver on an subinterval. Thus every solver also
 * has its own map of vectors, which sometimes make data-sharing a bit complicated. The user
 * does not have to interact directly with them. He only needs to call runPara, everything else
 * is handled automatically by the class.
 */
template <int DIM>
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
     *
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
                          const double theta        = getParam<double>("theta", "Equation"),
                          const double dtfine       = getParam<double>("dt", "Equation"));

    /*! \brief Method comparing parareal with sequential execution.
     *
     * Executes the parareal algorithm with the specified coarse time step sizes and coarse thetas,
     * values for fine method is taken from the paramfile.
     * Before execution of parareal it solves the problem sequentially and calculates the norm
     * \f$l_2\f$ vector norm between the two solutions on every time node.
     * \param maxIterations maximal number of iterations executed by parareal algorithm
     * \param dtcoarse_vec Vector with coarse time step sizes that will be tested against sequential
     * execution
     * \param coarse_theta_vec Vector with coarse theta values for step sizes specified by
     * dtcoarse_vec
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
     * \brief Wrapper function for GetMultiLevelSolver returning *const* pointer to MultiLevelSolver
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
     * \brief Fine time stepping method.
     *
     * This method encapsulates *noFineSteps* timesteps into one call. The step size
     * used is dtfine. This method writes a logfile.
     * \param time          current time from where the steps are executed
     * \param u             VectorInterface with approximation at current time
     * \param f             Right hand inside
     */
    template <bool last>
    void step_fine_final(double time, VectorInterface& u, VectorInterface& f);

    /*!
     * \brief Coarse time stepping method.
     *
     * This method encapsulates *noCoarseSteps* timesteps into one call. The step size
     * used is dtcoarse.
     * \param time          current time from where the steps are executed
     * \param u             VectorInterface with approximation at current time
     * \param f             Right hand inside
     */
    void step_coarse(double time, VectorInterface& u, VectorInterface& f);

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
     * using the maximum potential of parallelization.
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

    static const std::vector<GlobalVector> setTimeGridVector(std::size_t noIntervals);

    VectorInterface
      fi_solution;  //!< Vector containing the fine solution on the corresponding interval.
    VectorInterface
      co_solution;        //!< Vector containing the coarse solution on the corresponding interval.
    VectorInterface old;  //!< Vector for temporary results.
    VectorInterface
      u_solution;  //!< Vector containing the final solution on the corresponding interval.
    std::vector<DoubleVector> log_buffer;  //!< Vector for buffering log entries.

    static double theta;
    static double coarse_theta;

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

#if PDBG >= 1
    static double coarse_exec_time;
    static double fine_exec_time;
#endif

    static std::vector<GlobalVector>
      u_sol_arr;                    //!< Vector containing the final solutions on every time node
    static std::ofstream func_log;  //! Ofstream object for logfile writing
};
// end class definition

/***************************************************************************************************
public methods
***************************************************************************************************/

template <int DIM>
parareal<DIM>::parareal(size_t id)
    : fi_solution("f" + std::to_string(id))
    , co_solution("g" + std::to_string(id))
    , old("old" + std::to_string(id))
    , u_solution("u" + std::to_string(id))
    , log_buffer(noFineSteps, DoubleVector(5, 0.))
{
    // StdLoop();
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
void parareal<DIM>::outputParameters(std::ostream& ostrm)
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
          << "\nTheta for fine steps                   " << theta
          << "\nUsing paramfile:                       " << paramfile
          << "\nProblemlabel is:                       " << problemlabel
          << "\n=========================================================\n\n";
    // clang-format on
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
double parareal<DIM>::runPara(const int maxIterations, const double coarse_theta,
                              const double dtcoarse, const double theta, const double dtfine)
{
    assert(dtcoarse >= dtfine && "coarse step size should be bigger than fine step size");
    parareal::maxIterations = maxIterations;
    parareal::coarse_theta  = coarse_theta;
    parareal::dtcoarse      = dtcoarse;
    parareal::theta         = theta;
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
    assert(theta > 0 && coarse_theta > 0 && "Thetas have to be positive");
    assert(dtcoarse * noCoarseSteps - dtfine * noFineSteps <= numeric_limits<double>::epsilon()
           && "Coarse and fine step sizes and number of time steps don't match");

    auto exec_time = paraAlgo();

    // clang-format off
    std::cerr << "\n==================================================\n"
                 "Finished parareal.\nElapsed time is\t" << exec_time << '\n';
#if PDBG >= 1
    std::cerr << ".\nTime spent in coarse method: "      << coarse_exec_time
              << "\nTime spent in fine method: "         << fine_exec_time;
#endif
    // clang-format on
    return exec_time;
}

//--------------------------------------------------------------------------------------------------

#ifdef _OPENMP  // parareal can't be run without OpenMP
template <int DIM>
double parareal<DIM>::paraAlgo()
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
    for (size_t i = 0; i < noIntervals; ++i)
    {
        loop_interval[i] = new parareal<DIM>(i);
        loop_interval[i]->setupProblem(paramfile, problemlabel);

        loop_interval[i]->GetMultiLevelSolver()->ReInit(problemlabel);
        loop_interval[i]->GetSolver()->OutputSettings();
        cout << loop_interval[i]->GetMeshAgent()->ncells() << " cells" << '\n';

        loop_interval[i]->GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
        loop_interval[i]->SetTime(dtfine, start_time);
    }

    for (size_t m = 0; m < noIntervals; ++m)
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

    for (size_t m = 0; m < noIntervals; ++m)
    {
        loop_interval[m]->InitSolution(u);
        loop_interval[m]->InitSolution(loop_interval[m]->fi_solution);
        loop_interval[m]->InitSolution(loop_interval[m]->co_solution);
        // loop_interval[m]->InitSolution(loop_interval[m]->old);
        loop_interval[m]->InitSolution(loop_interval[m]->u_solution);
    }

    for (size_t m = 0; m < noIntervals; ++m)
    {
        loop_interval[m]->setGVsize(u_sol_arr[m]);
    }

    /// Initializing locks
    auto interval_locker = new omp_lock_t[noIntervals];
    for (int m = 0; m < noIntervals; ++m)
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
        for (int i = 0; i < noIntervals; ++i)
        {
#pragma omp ordered
            {
                omp_set_lock(&interval_locker[i]);
                // u ⟵ G⁰(i)
                loop_interval[i]->step_coarse(time + i * intervalLength, u, f);
                // G⁰(i) ⟵ u;
                loop_interval[i]->GetSolver()->Equ(loop_interval[i]->co_solution, 1., u);
                // U⁰(i) ⟵ u;
                loop_interval[i]->GetSolver()->Equ(loop_interval[i]->u_solution, 1., u);

// Visualization for debugging
#if PDBG > 1
                {
                    loop_interval[i]->GetSolver()->Visu("Results/initU", u, i);
                    loop_interval[i]->WriteMeshAndSolution("Results/initU", u);

                    loop_interval[i]->GetSolver()->Visu("Results/initsol",
                                                        loop_interval[i]->u_solution, i);
                    loop_interval[i]->WriteMeshAndSolution("Results/initsol",
                                                           loop_interval[i]->u_solution);

                    loop_interval[i]->GetSolver()->Visu("Results/initfine",
                                                        loop_interval[i]->fi_solution, i);
                    loop_interval[i]->WriteMeshAndSolution("Results/initfine",
                                                           loop_interval[i]->fi_solution);
                }
#endif

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
        for (int k = 1; k <= maxIterations; ++k)
        {
            // fine propagations
#pragma omp for nowait schedule(monotonic : static)
            for (int m = 0; m < noIntervals - k + 1; ++m)
            {
                // #pragma omp ordered  // for debugging
                //                 {
                omp_set_lock(&interval_locker[m]);
                // Fᵏ(m) ⟵ F(Uᵏ⁻¹(m-1))
                if (m == 0 || k == maxIterations)
                {
                    loop_interval[m]->template step_fine_final<true>(
                      time + m * intervalLength, loop_interval[m]->fi_solution, f);
                }
                else  // we are not in the last iteration and not exact
                {
                    loop_interval[m]->template step_fine_final<false>(
                      time + m * intervalLength, loop_interval[m]->fi_solution, f);
                }
#if PDBG > 1
                {
                    loop_interval[m]->GetSolver()->Visu(
                      "Results/fineres", loop_interval[m]->fi_solution, (k - 1) * 10 + m);
                    loop_interval[m]->WriteMeshAndSolution("Results/fineres",
                                                           loop_interval[m]->fi_solution);
                }
#endif
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
            for (int m = 0; m < noIntervals - k; ++m)
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
                    loop_interval[m]->step_coarse(time + (m + 1) * intervalLength,
                                                  loop_interval[m]->co_solution, f);
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
#if PDBG > 1
                    {
                        loop_interval[m]->GetSolver()->Visu(
                          "Results/coarseres", loop_interval[m]->co_solution, (k - 1) * 10 + m);
                        loop_interval[m]->WriteMeshAndSolution("Results/coarseres",
                                                               loop_interval[m]->co_solution);
                        loop_interval[m]->GetSolver()->Visu(
                          "Results/iterres", loop_interval[m]->u_solution, (k - 1) * 10 + m);
                        loop_interval[m]->WriteMeshAndSolution("Results/iterres",
                                                               loop_interval[m]->u_solution);
                    }
#endif
                }

            }  // Corrections

#pragma omp atomic update  // increase time
            time += intervalLength;
        }  // Iterations
    }      // parallel region
    auto exec_time = omp_get_wtime() - start_wall;

    for (size_t m = 0; m < noIntervals - maxIterations; ++m)
    {
        // write final solution vectors
        loop_interval[m]->setGVsize(u_sol_arr[m]);
#if PDBG > 1
        assert(
          u_sol_arr[m].size()
            == loop_interval[m + 1]->GetSolver()->GetGV(loop_interval[m + 1]->u_solution).size()
          && u_sol_arr[m].ncomp()
               == loop_interval[m + 1]->GetSolver()->GetGV(loop_interval[m + 1]->u_solution).ncomp()
          && u_sol_arr[m].n()
               == loop_interval[m + 1]->GetSolver()->GetGV(loop_interval[m + 1]->u_solution).n()
          && "Dimension mismatch");
#endif
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
        loop_interval[m + maxIterations]->GetSolver()->Visu(
          "Results/p", loop_interval[m + maxIterations]->u_solution, m + maxIterations);
        loop_interval[m + maxIterations]->WriteMeshAndSolution(
          "Results/p", loop_interval[m + maxIterations]->u_solution);
        std::cout << "wrote subinterval " << m + maxIterations << '\n';
    }

    // Clean up
    omp_destroy_lock(interval_locker);
    delete[] interval_locker;
    for (size_t i = 0; i < noIntervals; ++i)
    {
        delete loop_interval[i];
    }
    func_log.close();
    return exec_time;
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
void parareal<DIM>::compareParaSerial(const unsigned int maxIterations,
                                      const std::vector<double>& dtcoarse_vec,
                                      const std::vector<double>& coarse_theta_vec)
{
    // ParamFile paramfile(paramstr);
    ofstream timings_log({"timingsI" + to_string(maxIterations) + ".txt"});
    timings_log << "ratio time accel\n";

    // serial part
    auto start   = omp_get_wtime();
    auto loopSeq = new parareal<DIM>(0);
    loopSeq->setupProblem(paramfile, problemlabel);
    loopSeq->GetMultiLevelSolver()->ReInit(problemlabel);
    loopSeq->GetSolver()->OutputSettings();
    cout << loopSeq->GetMeshAgent()->ncells() << " cells" << '\n';
    VectorInterface f("f");
    double time = start_time;
    std::vector<VectorInterface> U_seq(noIntervals, loopSeq->old);
    for (size_t i = 0; i < noIntervals; ++i)
    {
        U_seq[i].SetName({"useq" + std::to_string(i)});
        U_seq[i].SetType("node");
        loopSeq->GetSolver()->ReInitVector(U_seq[i]);
    }
    for (size_t i = 0; i < noIntervals; ++i)
    {
        loopSeq->InitSolution(U_seq[i]);
    }

    loopSeq->GetMultiLevelSolver()->ReInitVector(loopSeq->old);
    loopSeq->GetMultiLevelSolver()->ReInitVector(loopSeq->u_solution);
    loopSeq->GetMultiLevelSolver()->ReInitVector(f);
    loopSeq->InitSolution(loopSeq->u_solution);
    loopSeq->InitSolution(loopSeq->old);

    if (func_log.is_open())
        func_log.close();
    func_log.open("functional.txt");
    // sequential run
    for (int i = 0; i < noIntervals; ++i)
    {
        loopSeq->template step_fine_final<true>(time + i * intervalLength, loopSeq->u_solution, f);
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

    assert(dtcoarse_vec.size() == coarse_theta_vec.size()
           && "For every dtcoarse there has to be one coarse_theta given and the other way round");

    string dtc = to_string(dtcoarse_vec[0]);
    ofstream error_log;

    for (size_t i = 0; i < dtcoarse_vec.size(); ++i)
    {
        dtc = to_string(dtcoarse_vec[i]);
        dtc = dtc.substr(dtc.find(".") + 1);
        func_log.open("para_functionalC" + dtc + ".txt");

        // run Parareal
        auto ratio           = static_cast<int>(dtcoarse_vec[i] / dtfine);
        double parallel_time = runPara(maxIterations, coarse_theta_vec[i], dtcoarse_vec[i]);
        timings_log << ratio << ' ' << parallel_time << ' ' << serial_time / parallel_time << '\n';
        func_log.close();

        // calculate error
        error_log.open({"error_paraC" + dtc + "I" + to_string(maxIterations) + ".txt"});
        error_log << "interval error\n";
        for (size_t i = 0; i < noIntervals; ++i)
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

#else  // if OpenMP is not available, exit program when methods are called

//--------------------------------------------------------------------------------------------------

template <int DIM>
double parareal<DIM>::paraAlgo()
{
    std::cerr << "Running parareal relies on OpenMP" << '\n';
    std::exit(EXIT_FAILURE);
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
void parareal<DIM>::compareParaSerial(const unsigned int maxIterations,
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

template <int DIM>
const SolverInterface* parareal<DIM>::GetSolver() const
{
    return StdLoop::GetMultiLevelSolver()->GetSolver();
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
SolverInterface* parareal<DIM>::GetSolver()
{
    return StdLoop::GetMultiLevelSolver()->GetSolver();
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
const MultiLevelSolverInterface* parareal<DIM>::GetMultiLevelSolver() const
{
    return StdLoop::GetMultiLevelSolver();
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
MultiLevelSolverInterface* parareal<DIM>::GetMultiLevelSolver()
{
    return StdLoop::GetMultiLevelSolver();
}

/***************************************************************************************************
private methods
***************************************************************************************************/

template <int DIM>
template <bool last>
void parareal<DIM>::step_fine_final(double time, VectorInterface& u, VectorInterface& f)
{
#if PDBG >= 1
    auto start = omp_get_wtime();
#endif
    dynamic_cast<const FSI<DIM>*>(GetSolver()->GetProblemDescriptor()->GetEquation())->getTheta() =
      theta;
    //
    SetTime(dtfine, time);
    DoubleVector juh(4, 0.);
    for (auto iter = 0; iter < noFineSteps; ++iter)
    {
        cout << "Fine " << time << " -> " << time + dtfine << '\n';
        time += dtfine;
        SetTime(dtfine, time);
        GetMultiLevelSolver()->Equ(old, 1.0, u);
        GetMultiLevelSolver()->AddNodeVector("old", old);
        assert(StdLoop::Solve(u, f) == "converged");
        if (last)
        {
            log_buffer[iter][0] = time;
            juh                             = StdLoop::Functionals(u, f);
            for (short i = 0; i < 4; i++)
            {
                log_buffer[iter][i + 1] = juh[i];
            }
        }
        GetMultiLevelSolver()->DeleteNodeVector("old");
    }

#if PDBG >= 1
    auto end = omp_get_wtime();
#pragma omp atomic update
    fine_exec_time += (end - start);
#endif
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
void parareal<DIM>::step_coarse(double time, VectorInterface& u, VectorInterface& f)
{
//
#if PDBG >= 1
    auto start = omp_get_wtime();
#endif
    dynamic_cast<const FSI<DIM>*>(GetSolver()->GetProblemDescriptor()->GetEquation())->getTheta() =
      coarse_theta;
    SetTime(dtcoarse, time);
    for (auto iter = 0; iter < noCoarseSteps; ++iter)
    {
        cout << "Coarse " << time << " -> " << time + dtfine << '\n';
        time += dtcoarse;
        SetTime(dtcoarse, time);
        GetMultiLevelSolver()->Equ(old, 1.0, u);

        GetMultiLevelSolver()->AddNodeVector("old", old);
        assert(StdLoop::Solve(u, f) == "converged");
        GetMultiLevelSolver()->DeleteNodeVector("old");
    }

#if PDBG >= 1
    auto end = omp_get_wtime();
#pragma omp atomic update
    coarse_exec_time += (end - start);
#endif
}

template <int DIM>
void parareal<DIM>::SetTime(const double dt, const double time) const
{
    for (int i = 0; i < GetMultiLevelSolver()->nlevels(); ++i)
    {
        GetMultiLevelSolver()->GetSolver(i)->GetProblemDescriptor()->SetTime(time, dt);
    }
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
void parareal<DIM>::setGV(VectorInterface& v, const GlobalVector& v_ptr) const
{
    GetSolver()->GetGV(v) = v_ptr;
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
void parareal<DIM>::setGVsize(GlobalVector& v_ptr)
{
    auto ncomp    = GetSolver()->GetProblemDescriptor()->GetEquation()->GetNcomp();
    auto n        = GetSolver()->GetDiscretization()->n();
    v_ptr.ncomp() = ncomp;
    v_ptr.resize(n);
    // v_ptr.zero();
}

//--------------------------------------------------------------------------------------------------

template <int DIM>
void parareal<DIM>::setupProblem(const ParamFile& paramfile, const std::string& problemlabel)
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

template <int DIM>
int parareal<DIM>::setNoIntervals()
{
    auto noIntervals = omp_get_num_threads();
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

//--------------------------------------------------------------------------------------------------

template <int DIM>
template <typename T>
const T parareal<DIM>::getParam(const std::string& paramName, const std::string& block)
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

template <int DIM>
const std::vector<GlobalVector> parareal<DIM>::setTimeGridVector(const size_t noIntervals)
{
    std::vector<GlobalVector> tGridVec_ptr(noIntervals);
    tGridVec_ptr.reserve(noIntervals);
    for (size_t i = 0; i < noIntervals; ++i)
    {
        tGridVec_ptr[i] = GlobalVector();
    }
    return tGridVec_ptr;
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

template <int DIM>
ParamFile parareal<DIM>::paramfile = ParamFile("fsi-3.param");
template <int DIM>
string parareal<DIM>::problemlabel = "fsi";

template <int DIM>
double parareal<DIM>::theta = parareal::getParam<double>("theta", "Equation");
template <int DIM>
double parareal<DIM>::coarse_theta = parareal::getParam<double>("thetacoarse", "Equation");

template <int DIM>
unsigned int parareal<DIM>::noIntervals = static_cast<unsigned int>(setNoIntervals());
template <int DIM>
unsigned int parareal<DIM>::maxIterations = static_cast<unsigned int>(1);
template <int DIM>
double parareal<DIM>::start_time = parareal::getParam<double>("start_time", "Equation");
template <int DIM>
double parareal<DIM>::stop_time = parareal::getParam<double>("stop_time", "Equation");

template <int DIM>
double parareal<DIM>::intervalLength = static_cast<double>((stop_time - start_time) / noIntervals);

template <int DIM>
double parareal<DIM>::dtfine = parareal::getParam<double>("dt", "Equation");
template <int DIM>
size_t parareal<DIM>::noFineSteps = static_cast<size_t>(intervalLength / dtfine);

template <int DIM>
double parareal<DIM>::dtcoarse = parareal::getParam<double>("dtcoarse", "Equation");
template <int DIM>
size_t parareal<DIM>::noCoarseSteps = static_cast<size_t>(parareal::intervalLength
                                                          / parareal::dtcoarse);
#if PDBG >= 1
template <int DIM>
double parareal<DIM>::coarse_exec_time = 0;

template <int DIM>
double parareal<DIM>::fine_exec_time = 0;
#endif

template <int DIM>
std::vector<GlobalVector> parareal<DIM>::u_sol_arr = parareal::setTimeGridVector(noIntervals);
template <int DIM>
std::ofstream parareal<DIM>::func_log("functional_para.txt");
}  // namespace Gascoigne
#endif  // PARAREAL_H
