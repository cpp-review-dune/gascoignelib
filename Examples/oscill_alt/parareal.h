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

namespace Gascoigne {
enum class name_type { error, functional };
enum class iter_type {
  fine_last   = 0,
  fine        = 1,
  coarse_first= 2,
  coarse      = 3,
};
enum class log_level { results= 0, visu_info= 1, visu_detailed= 2 };

struct dir_vec {
  std::string     dir;
  VectorInterface vec;
  size_t          index;
};

using implicit= std::integral_constant<bool, true>;
using not_implicit std::integral_constant<bool, false>;

// Step wrapper used for parareal algorithm
template <iter_type method, typename fctor, typename... fctor_types>
struct stepper {
  static inline constexpr void
    propagator(size_t n_steps, size_t implicit_steps, fctor_types&&... fctor_args) {
    for (size_t i= 0; i < implicit_steps; ++i) {
      fctor::template step<method, implicit>(std::forward<fctor_types>(fctor_args)...);
    }
    for (size_t i= implicit_steps; i < n_steps; ++i) {
      fctor::template step<method, not_implicit>(
        std::forward<fctor_types>(fctor_args)...);
    }
  }
};

template <int DIM, log_level logging, typename Loop, typename Vector>
double parareal_algorithm(const int n_intervals, const int max_iterations) {
  static auto equal= [&](parareal const* const& para_dest,
                         VectorInterface&       dest,
                         parareal const* const& para_src,
                         VectorInterface&       source) {
    auto dest_begin= para_dest->GetSolver()->GetGV(dest).begin();
    auto src_begin = para_src->GetSolver()->GetGV(source).cbegin();

    const auto length= para_src->GetSolver()->GetGV(source).cend() - src_begin;
    for (auto i= 0; i < length; ++i) {
      *dest_begin= *src_begin;
      ++src_begin;
      ++dest_begin;
    }
  };

  static auto correction_w_weigths= [&](parareal const* const& para_dest,
                                        VectorInterface&       dest,
                                        parareal const* const& src_x,
                                        VectorInterface&       x,
                                        parareal const* const& src_y,
                                        VectorInterface&       y,
                                        parareal const* const& src_z,
                                        VectorInterface&       z) {
    auto       dest_begin= para_dest->GetSolver()->GetGV(dest).begin();
    const auto weight    = para_dest->coar_weight;
    auto       x_begin   = src_x->GetSolver()->GetGV(x).cbegin();
    auto       y_begin   = src_y->GetSolver()->GetGV(y).cbegin();
    auto       z_begin   = src_z->GetSolver()->GetGV(z).cbegin();

    const auto length= src_x->GetSolver()->GetGV(x).cend() - x_begin;
    for (auto i= 0; i < length; ++i) {
      *dest_begin= weight * *x_begin + *y_begin - weight * *z_begin;
      ++x_begin;
      ++y_begin;
      ++z_begin;
      ++dest_begin;
    }
  };

  static auto set_coar_weight= [&](parareal* const& para_solver,
                                   VectorInterface& f,
                                   VectorInterface& c,
                                   double           damping) {
    const nvector<double> c_c= para_solver->GetSolver()->GetGV(c).CompNorm();
    const nvector<double> f_f= para_solver->GetSolver()->GetGV(f).CompNorm();
    nvector<double>       c_f(2 * DIM + 1, 0.0);
    para_solver->GetSolver()->GetGV(f).ScalarProductComp(
      c_f, para_solver->GetSolver()->GetGV(c));

    auto cw_iter = c_f.begin();
    auto c_c_iter= c_c.cbegin();
    auto f_f_iter= f_f.cbegin();
    for (auto i= 0; i < 2 * DIM + 1; ++i) {
      *cw_iter/= *c_c_iter * *f_f_iter;
      ++cw_iter;
      ++c_c_iter;
      ++f_f_iter;
    }
    double sc= 0.0;
    for (auto i= 0; i < DIM + 1; i++) {
      sc+= c_f[i];
    }
    sc/= DIM + 1;
    sc*= damping;
    para_solver->coar_weight= sc;
  };

  double     time      = start_time;
  const auto start_wall= omp_get_wtime();
  //
  // Initialization of solvers and problems and solutions from parareal
  // coar_sol   - VectorInterface for results of coarse method
  // fine_sol   - VectorInterface for results of fine method
  // end_sol    - VectorInterface for final results
  // old           - VectorInterface for temporary results
  // u             - VectorInterface used for initialization of the others
  //
  VectorInterface        u("u");
  VectorInterface        f("f");
  std::vector<parareal*> subinterval(n_intervals);
  std::vector<Vector>    fine_sol(n_intervals);
  std::vector<Vector>    coar_sol(n_intervals);
  std::vector<Vector>    end_sol(n_intervals);
  std::vector<Vector>    old(n_intervals);
  // std::vector<Vector> tmp_sol(n_intervals);
  double coar_weight;

  for (auto m= 0; m < n_intervals; ++m) {
    subinterval[m]= new parareal<DIM, logging>(m + 1);
    fine_sol[m]   = "f" + std::to_string(m);
    coar_sol[m]   = "g" + std::to_string(m);
    end_sol[m]    = "u" + std::to_string(m);
    old[m]        = "old" + std::to_string(m);
  }
  for (auto m= 0; m < n_intervals; ++m) {
    subinterval[m]->setup_problem(paramfile, problemlabel);

    subinterval[m]->GetMultiLevelSolver()->ReInit(problemlabel);
    subinterval[m]->GetSolver()->OutputSettings();
    cout << subinterval[m]->GetMeshAgent()->ncells() << " cells" << '\n';

    subinterval[m]->GetSolverInfos()->GetNLInfo().control().matrixmustbebuild()= 1;
    subinterval[m]->set_time(dtfine, start_time);

    subinterval[m]->GetMultiLevelSolver()->ReInitVector(f);
    subinterval[m]->GetMultiLevelSolver()->ReInitVector(u);
    subinterval[m]->GetMultiLevelSolver()->ReInitVector(fine_sol[m]);
    subinterval[m]->GetMultiLevelSolver()->ReInitVector(coar_sol[m]);
    subinterval[m]->GetMultiLevelSolver()->ReInitVector(old[m]);
    // subinterval[m]->GetMultiLevelSolver()->ReInitVector(tmp_vi[m]);
    subinterval[m]->GetMultiLevelSolver()->ReInitVector(end_sol[m]);

    subinterval[m]->setGVsize(subinterval[m]->GetSolver()->GetGV(fine_sol[m]));
    subinterval[m]->setGVsize(subinterval[m]->GetSolver()->GetGV(coar_sol[m]));
    subinterval[m]->setGVsize(subinterval[m]->GetSolver()->GetGV(old[m]));
    // subinterval[m]->setGVsize(subinterval[m]->GetSolver()->GetGV(tmp_vi[m]));
    subinterval[m]->setGVsize(subinterval[m]->GetSolver()->GetGV(end_sol[m]));

    subinterval[m]->InitSolution(u);
    subinterval[m]->InitSolution(fine_sol[m]);
    subinterval[m]->InitSolution(coar_sol[m]);
    subinterval[m]->InitSolution(end_sol[m]);
  }

  for (auto m= 0; m < n_intervals; ++m) {
    subinterval[m]->setGVsize(u_sol_arr[m]);
  }
  // precompute till given time
  if (precompute_time > 0) {
    auto precompute_steps= static_cast<int>(precompute_time / dtfine);
    subinterval[0]->template propagator<iter_type::fine>(
      time, u, f, 0, precompute_steps, 0);
    subinterval[0]
      ->GetSolver()
      ->GetGV(subinterval[0]->end_sol)
      .equ(1., subinterval[0]->GetSolver()->GetGV(u));
    for (auto m= 0; m < n_intervals; ++m) {
      equal(subinterval[m], end_sol[m], subinterval[0], u);
      equal(subinterval[m], fine_sol[m], subinterval[0], u);
      equal(subinterval[m], coar_sol[m], subinterval[0], u);
      if (m > 0) {
        equal(subinterval[m], u, subinterval[0], u);
      }
    }
    time+= precompute_time;
  }
  /// Initializing locks
  auto interval_locker= new omp_lock_t[n_intervals];
  for (auto m= 0; m < n_intervals; ++m) {
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
    for (auto i= 0; i < n_intervals; ++i) {
#pragma omp ordered
      {
        // auto im  = i % threads;
        // auto im1 = (i + 1) % threads;
        omp_set_lock(&interval_locker[i]);
        // u ⟵ G⁰(i)
        subinterval[i]->template propagator<iter_type::coarse_first>(
          time + i * interval_length,
          u,
          f,
          0,
          n_coarsesteps,
          subinterval[i]->c_implicit_steps);
        // G⁰(i) ⟵ u;
        equal(subinterval[i], subinterval[i]->coar_sol, subinterval[i], u);
        // U⁰(i) ⟵ u;
        equal(subinterval[i], subinterval[i]->end_sol, subinterval[i], u);

        // Visualization for debugging
        if constexpr (logging == log_level::visu_info
                      || logging == log_level::visu_detailed) {
          std::cerr << logwrite::info("Done", "Initialization on subinterval ", i);
          subinterval[i]->visu("Results/initsol", subinterval[i]->end_sol, i);
        }

        if (i < n_intervals - 1) {
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
    for (auto k= 1; k <= max_iterations; ++k) {
      // fine propagations
#pragma omp for nowait schedule(static)
      for (auto m= 0; m < n_intervals - k + 1; ++m) {
        // #pragma omp ordered  // for debugging
        //                 {
        // omp_set_lock(&interval_locker[m]);
        // Fᵏ(m) ⟵ F(Uᵏ⁻¹(m-1))
        if (m == 0 || k == max_iterations) {
          omp_set_lock(&interval_locker[m]);
          subinterval[m]->template propagator<iter_type::fine_last>(
            time + m * interval_length,
            fine_sol[m],
            f,
            k,
            n_finesteps,
            subinterval[m]->f_implicit_steps);
          omp_unset_lock(&interval_locker[m]);
        } else  // we are not in the last iteration and not exact
        {
          omp_set_lock(&interval_locker[m]);
          subinterval[m]->template propagator<iter_type::fine>(
            time + m * interval_length,
            fine_sol[m],
            f,
            k,
            n_finesteps,
            subinterval[m]->f_implicit_steps);
          omp_unset_lock(&interval_locker[m]);
        }
        if constexpr (logging == log_level::visu_info
                      || logging == log_level::visu_detailed) {
          std::cerr << logwrite::info("Done",
                                      "Fine solution on subinterval ",
                                      m + k - 1,
                                      " in iteration number ",
                                      k);
          subinterval[m]->visu("Results/fineres", fine_sol[m], (k - 1) * 10 + m);
        }
        // omp_unset_lock(&interval_locker[m]);
        //}
      }
      // barrier for debugging
      // #pragma omp barrier
      //
      // coarse propagations and corrections
      //

#pragma omp for ordered nowait schedule(monotonic : static)
      // Corrections
      for (auto m= 0; m < n_intervals - k; ++m) {
#pragma omp ordered
        {
          // k-th value is exact in k-th Iteration
          if (m == 0) {
            omp_set_lock(&interval_locker[0]);
            // Uᵏ(k) ⟵ Fᵏ(k)
            equal(subinterval[0],
                  subinterval[0]->end_sol,
                  subinterval[0],
                  subinterval[0]->fine_sol);
            u_sol_arr[k - 1].equ(
              1., subinterval[0]->GetSolver()->GetGV(subinterval[0]->fine_sol));
            for (const auto& x : subinterval[0]->log_buffer) {
              func_log << x << '\n';
            }

            subinterval[0]->visu("Results/p", subinterval[0]->end_sol, k - 1);
            omp_unset_lock(&interval_locker[0]);
          }

          // prepare predictor step
          omp_set_lock(&interval_locker[m]);
          equal(subinterval[m], coar_sol[m], subinterval[m], end_sol[m]);
          omp_unset_lock(&interval_locker[m]);
          // predictor step with coarse method
          subinterval[m]->template propagator<iter_type::coarse>(
            time + (m + 1) * interval_length,
            coar_sol[m],
            f,
            k,
            n_coarsesteps,
            subinterval[m]->c_implicit_steps);

          // omp_unset_lock(&interval_locker[m]);
          omp_set_lock(&interval_locker[m + 1]);
          // calculate weights
          set_coar_weight(subinterval[m + 1],
                          subinterval[m + 1]->fine_sol,
                          subinterval[m + 1]->coar_sol,
                          0.95);

          // calculate corrector term
          // clang-format off
                    correction_w_weigths(
                        subinterval[m + 1], subinterval[m + 1]->end_sol,  // source vector
                        subinterval[m], coar_sol[m],         //+ w_x * coar_sol_old
                        subinterval[m + 1], subinterval[m + 1]->fine_sol, //+ w_y * fine_sol
                        subinterval[m + 1], subinterval[m + 1]->coar_sol  //+ w_z * coar_sol_new
                    );
                    //clang-format on
                    // prepare next iteration
                    equal(subinterval[m + 1], subinterval[m + 1]->fine_sol, subinterval[m + 1],
                          subinterval[m + 1]->end_sol);
                    omp_unset_lock(&interval_locker[m + 1]);

                    // Visualization for debugging
                    if constexpr (logging == log_level::visu_info
                                  || logging == log_level::visu_detailed)
                    {
                        subinterval[m]->visu(
                          dir_vec{"Results/coarseres", coar_sol[m], (k - 1) * 10 + m},
                          dir_vec{"Results/iterres", end_sol[m], (k - 1) * 10 + m});
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
        subinterval[m + 1]->visu("Results/p", subinterval[m + 1]->end_sol, m + max_iterations);
        std::cerr << logwrite::info("Done", "Wrote subinterval ", m + max_iterations);
    }

    // Clean up
    omp_destroy_lock(interval_locker);
    return exec_time;
}



}  // namespace Gascoigne
