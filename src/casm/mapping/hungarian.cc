#include "casm/mapping/hungarian.hh"

#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace mapping {
namespace hungarian {
namespace hungarian_impl {

/// \brief For every row of cost_matrix, subtract the minimum element
///      from every entry in that row.
///
/// \param cost_matrix Cost matrix to reduce
/// \param infinity Skip rows where minimum element > infinity.
///
void reduce_cost(Eigen::MatrixXd &cost_matrix, double infinity) {
  // cost matrix dimension
  int dim = cost_matrix.rows();

  // Step 1. For every row subtract the minimum element from every
  // entry in that row.
  for (int i = 0; i < dim; i++) {
    double row_min = cost_matrix.row(i).minCoeff();
    if (row_min > infinity) continue;
    for (int j = 0; j < dim; j++) {
      cost_matrix(i, j) -= row_min;
    }
  }
}

/// \brief Update zero_marks to find zeros and "star" them
///
/// \param cost_matrix The cost matrix
/// \param zero_marks Matrix indicating the location of zeros in
///     the cost matrix and their type (prime==1, starred==-1).
/// \param tol Tolerance for finding zeros
///
void find_zeros(Eigen::MatrixXd const &cost_matrix, Eigen::MatrixXi &zero_marks,
                double tol) {
  // cost matrix dimension
  int dim = cost_matrix.rows();
  int star = -1;
  bool is_star;
  // Step 2. Find a zero in the matrix
  // find zeros and star
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      if (almost_zero(cost_matrix(i, j), tol)) {
        is_star = false;
        for (int k = 0; k < dim; k++) {
          if (zero_marks(i, k) == star) {
            is_star = true;
          }
        }
        for (int l = 0; l < dim; l++) {
          if (zero_marks(l, j) == star) {
            is_star = true;
          }
        }
        if (!is_star) {
          zero_marks(i, j) = star;
        }
      }
    }
  }
}

/// \brief Check if complete assignment has been found
///
/// \param zero_marks Matrix indicating the location of zeros in
///     the cost matrix and their type (prime==1, starred==-1).
/// \param col_covered Indicates columns that are "covered". Updated
///     so that `col_covered(j)=1` if `zero_marks(i,j)==-1` for any i.
///
/// \return true if the number of columns covered == number of cost_matrix rows
bool check_assignment(const Eigen::MatrixXi &zero_marks,
                      Eigen::VectorXi &col_covered) {
  // cover columns with starred zeros and check assignment. If N colmn
  // covered we are done else return false.
  for (int i = 0; i < zero_marks.rows(); i++) {
    for (int j = 0; j < zero_marks.cols(); j++) {
      if (zero_marks(i, j) == -1) {
        col_covered(j) = 1;
      }
    }
  }

  if (col_covered.sum() == zero_marks.rows()) {
    return true;
  } else
    return false;
}

//*******************************************************************************************

/// \brief TODO
int prime_zeros(const Eigen::MatrixXd &cost_matrix,
                Eigen::VectorXi &row_covered, Eigen::VectorXi &col_covered,
                Eigen::MatrixXi &zero_marks, double &min,
                Eigen::VectorXi &first_prime_zero, const double tol,
                double infinity) {
  int prime = 1;
  int covered = 1;
  int uncovered = 0;

  bool DONE = false;

  // In step 4 we are looking for a noncovered zero.
  // If one is found it is primed. If no star is found
  // in its row, then step 5 is called. If a star is
  // found, the row containing the prime is covered,
  // and the column containing the star is uncovered.
  // This is continued until no uncovered zeros exist.
  while (!DONE) {
    DONE = true;
    for (int i = 0; i < cost_matrix.rows(); i++) {
      for (int j = 0; j < cost_matrix.cols(); j++) {
        if (almost_zero(cost_matrix(i, j), tol) && col_covered(j) == 0 &&
            row_covered(i) == 0) {
          // find an uncovered zero and prime it
          zero_marks(i, j) = prime;
          // Determine if there is a starred zero in the row, if so, cover row i
          // and uncover the column of the starred zero, k.
          bool is_star = false;
          for (int k = 0; k < zero_marks.cols(); k++) {
            if (zero_marks(i, k) == -1 && !is_star) {
              row_covered(i) = covered;
              col_covered(k) = uncovered;
              is_star = true;
            }
          }

          // If no starred zero, save the last primed zero and proceed to
          // STEP 5.
          if (!is_star) {
            first_prime_zero(0) = i;
            first_prime_zero(1) = j;
            return 5;
          }
        }
      }
    }
    for (int i = 0; i < cost_matrix.rows(); i++) {
      for (int j = 0; j < cost_matrix.cols(); j++) {
        if (row_covered(i) == 0 && col_covered(j) == 0 &&
            almost_zero(cost_matrix(i, j), tol)) {
          DONE = false;
        }
      }
    }
  }

  // loop through all uncovered elements
  // to find the minimum value remaining.
  // This will be used in step 6.
  min = infinity;
  int return_code(-1);
  for (int i = 0; i < cost_matrix.rows(); i++) {
    for (int j = 0; j < cost_matrix.cols(); j++) {
      if (row_covered(i) == 0 && col_covered(j) == 0 &&
          cost_matrix(i, j) < min) {
        min = cost_matrix(i, j);
        return_code = 6;
      }
    }
  }
  return return_code;
}

/// \brief TODO
int alternating_path(const Eigen::MatrixXd &cost_matrix,
                     const Eigen::VectorXi &first_prime_zero,
                     Eigen::MatrixXi &zero_marks, Eigen::VectorXi &row_covered,
                     Eigen::VectorXi &col_covered) {
  bool done = false;

  // Initialize vectors to hold the locations
  // of the alternating path. The largest path
  // is at most 2*dimension of cost matrix -1.
  Eigen::VectorXi row_path(2 * cost_matrix.rows() - 1);
  Eigen::VectorXi col_path(2 * cost_matrix.cols() - 1);
  //
  // fill the paths initially with -1's so they
  // are not mistaken for points in path.
  row_path.fill(-1);
  col_path.fill(-1);
  //
  // Keep track of how long the path is.
  int path_counter = 1;

  // Set the first prime zero in the path
  // as the location specified from step 4.
  row_path(0) = first_prime_zero(0);
  col_path(0) = first_prime_zero(1);

  // These bools help test if there
  // are stars found in the prime's column.
  // Or if there are any primes found in
  // star's column.
  bool star_found = false;
  bool prime_found = false;

  // Build alternating path beginning with:
  //      uncovered primed zero (Z0)
  // ---> starred zero in same col as Z0 (Z1)
  // ---> primed zero in same row as Z1 (Z2)
  // ---> continue until a starred zero cannot be found
  //
  // The path is stored in row_path and col_path
  while (!done) {
    for (int i = 0; i < zero_marks.rows(); i++) {
      if (zero_marks(i, col_path(path_counter - 1)) == -1 && !star_found) {
        path_counter += 1;
        row_path(path_counter - 1) = i;
        col_path(path_counter - 1) = col_path(path_counter - 2);
        star_found = true;
        prime_found = false;
      }
      if (i == zero_marks.rows() - 1 && !star_found) {
        done = true;
      }
    }
    if (star_found && !done) {
      for (int j = 0; j < zero_marks.cols(); j++) {
        if (zero_marks(row_path(path_counter - 1), j) == 1 && !prime_found) {
          path_counter += 1;
          row_path(path_counter - 1) = row_path(path_counter - 2);
          col_path(path_counter - 1) = j;
          prime_found = true;
          star_found = false;
        }
      }
    }
  }

  // Unstar each starred zero and star each primed zero
  for (int i = 0; i < path_counter; i++)
    // Star primed zeros. I.e. even elements in path.
    if (i % 2 == 0 && zero_marks(row_path(i), col_path(i)) == 1) {
      zero_marks(row_path(i), col_path(i)) = -1;
    }
    // Unstar all starred zeros. I.e. odd elements.
    else if (zero_marks(row_path(i), col_path(i)) == -1) {
      zero_marks(row_path(i), col_path(i)) = 0;
    }
  for (int i = 0; i < zero_marks.rows(); i++) {
    for (int j = 0; j < zero_marks.cols(); j++) {
      if (zero_marks(i, j) == 1) {
        zero_marks(i, j) = 0;
      }
    }
  }

  // Uncover all lines
  row_covered.fill(0);
  col_covered.fill(0);
  return 3;
  // Step 3 should be run next by the while loop in main.
}

/// \brief TODO
int update_costs(Eigen::VectorXi const &row_covered,
                 Eigen::VectorXi const &col_covered, double min,
                 Eigen::MatrixXd &cost_matrix) {
  // Step 6. Add value from 4 to all covered rows, and subtract it from
  // uncovered columns. DO NOT alter stars, primes or covers. Return to 3.

  for (int i = 0; i < cost_matrix.rows(); i++) {
    for (int j = 0; j < cost_matrix.cols(); j++) {
      if (row_covered(i) == 1) {
        cost_matrix(i, j) += min;
      }
      if (col_covered(j) == 0) {
        cost_matrix(i, j) -= min;
      }
    }
  }
  return 4;
}

// *******************************************************************************************
/* Hungarian Algorithm Routines
 * Step 1: reduce the rows by smallest element
 * Step 2: star zeros
 * Step 3: cover columns with starred zeros and check assignement
 *         if K columns covered DONE. Else goto 4.
 * Step 4: Find uncovered zero and prime. if no starred zero in row
 *         goto 5. Otherwise, cover row, uncover column with star zero.
 *         Continue until all zeros are covered. Store smalles
 *         uncovered value goto 6.
 * Step 5: Build alternating prime and star zeros. Goto 3.
 * Step 6: Add value from 4 to all covered rows, and subtract it
 *         from uncovered columns. DO NOT alter stars, primes or covers.
 *         Return to 3.
 *
 */
// *******************************************************************************************
Assignment hungarian_method(Eigen::MatrixXd cost_matrix, double infinity,
                            double tol) {
  Assignment optimal_assignment;

  // cost matrix dimension
  int dim = cost_matrix.rows();

  // initialize matrix to carry stars and primes.
  Eigen::MatrixXi zero_marks = Eigen::MatrixXi::Zero(dim, dim);
  // initialize vectors to track (un)covered rows or columns
  Eigen::VectorXi row_covered(dim);
  Eigen::VectorXi col_covered(dim);

  int uncovered = 0;

  // begin uncovered
  row_covered.fill(uncovered);
  col_covered.fill(uncovered);

  // initialized by step 4
  double min = 0;

  // Vector to track the first_prime_zero found in step 4 and used in step 5.
  Eigen::VectorXi first_prime_zero(2);

  // Calling munkres step 1 to
  // reduce the rows of the cost matrix
  // by the smallest element in each row
  reduce_cost(cost_matrix, infinity);

  // Calling munkres step 2 to find
  // a zero in the reduced matrix and
  // if there is no starred zero in
  // its row or column, star it.
  // Repeat for all elements.
  //
  find_zeros(cost_matrix, zero_marks, tol);

  // Set DONE to false in order to control the while loop
  // which runs the main portion of the algorithm.
  bool DONE = false;

  // Main control loop for hungarian algorithm.
  // Can only exit the while loop if munkres step 3
  // returns true which means that there exists
  // an optimal assignment.
  //
  // Begin with step three.
  int next_step = 3;

  while (!DONE) {
    if (next_step == 3) {
      // Step 3 returns a boolean. TRUE if the assignment has been found
      // FALSE, if step 4 is required.
      if (check_assignment(zero_marks, col_covered) == true) {
        DONE = true;
      } else {
        next_step = 4;
      }
    }
    // Step 4 returns an int, either 5 or 6 (or -1 if failure detected)
    if (next_step == 4) {
      next_step = prime_zeros(cost_matrix, row_covered, col_covered, zero_marks,
                              min, first_prime_zero, tol, infinity);
    }
    // Step 5 returns an int, either 3 or 4.
    if (next_step == 5) {
      next_step = alternating_path(cost_matrix, first_prime_zero, zero_marks,
                                   row_covered, col_covered);
    }
    // Step 6 returns an int, always 4.
    if (next_step == 6) {
      next_step = update_costs(row_covered, col_covered, min, cost_matrix);
    }
    if (next_step == -1) {
      optimal_assignment.clear();
      return optimal_assignment;
    }
  }

  // Once the main control loop finishes, an
  // optimal assignment has been found. It
  // is denoted by the locations of the
  // starred zeros in the zero_marks matrix.
  // optimal_assignment vector contains the
  // results of the assignment as follows:
  // the index corresponds to the index of
  // the atom in the ideal POSCAR, and the
  // value at an index corresponds to the
  // atom of the relaxed structure.
  // i.e. optimal_assignemnt(1) = 0
  // means that the atom with index '1' in
  // the ideal POSCAR is mapped onto the atom
  // with index '0' in the relaxed CONTCAR.
  optimal_assignment.assign(cost_matrix.rows(), -1);

  for (int i = 0; i < zero_marks.rows(); i++) {
    for (int j = 0; j < zero_marks.cols(); j++) {
      if (zero_marks(i, j) == -1) {
        if (!valid_index(optimal_assignment[i])) {
          optimal_assignment[i] = j;
        } else {
          throw std::runtime_error(
              "Error in hungarian::solve: final assignment failure, multiple "
              "starred zeros");
        }
      }
    }
    if (!valid_index(optimal_assignment[i])) {
      throw std::runtime_error(
          "Error in hungarian::solve: final assignment failure, no starred "
          "zeros");
    }
  }
  return optimal_assignment;
}

}  // namespace hungarian_impl

/// \brief Find the optimal solution to the assignment problem
///
/// The assignment problem is: minimize the cost of assigning m
/// "workers" to n "tasks", where the cost of assigning "worker" i
/// to "task" j is cost_matrix(i,j).
///
/// \param cost_matrix The cost of assigning "worker" i to "task" j
///     is cost_matrix(i,j). The number of rows and columns must
///     be greater than 1. The number of rows must be equal to the
///     number of columns.
/// \param infinity Cost used for "infinity", when an assignment is
///     forced off
/// \param tol Tolerance used for comparing costs
///
/// \returns The optimal assignment solution, as a pair of
///     {cost, assignment}. The assignment vector gives
///     `j = assignment[i]`, where i is the worker (row) and j is
///     the task (column). If no solution is found (every solution
///     includes an infinity cost assignment), then the return
///     values is {infinity, {}}.
///
std::pair<double, Assignment> solve(Eigen::MatrixXd const &cost_matrix,
                                    double infinity, double tol) {
  // --- Input validation ---
  if (cost_matrix.rows() < 1) {
    throw std::runtime_error(
        "Error in hungarian::solve: cost_matrix.rows() < 1");
  }
  if (cost_matrix.cols() < 1) {
    throw std::runtime_error(
        "Error in hungarian::solve: cost_matrix.cols() < 1");
  }
  if (cost_matrix.rows() != cost_matrix.cols()) {
    throw std::runtime_error(
        "Error in hungarian::solve: cost_matrix.rows() != cost_matrix.cols()");
  }

  Assignment assignment =
      hungarian_impl::hungarian_method(cost_matrix, infinity, tol);

  // Failure, causes includes:
  // - Forced to assign task with infinity cost
  // - others?
  if (assignment.size() == 0) {
    return std::make_pair(infinity, assignment);
  }

  double cost = 0.0;
  // Find the costs associated with the optimal assignments
  for (int i = 0; i < assignment.size(); i++) {
    cost += cost_matrix(i, assignment[i]);
  }

  return std::make_pair(cost, assignment);
}

}  // namespace hungarian
}  // namespace mapping
}  // namespace CASM
