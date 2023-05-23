//
// Created by 苏立彪 on 2023/5/23.
//
#include "tcf_linear_system.h"

namespace tri_cross_field{


LinearSystem::LinearSystem() {
  solver_type_ = EigenSparseLU;
  factorize_ = false;
  preprocess_sparsity_pattern_ = false;
}

bool LinearSystem::IsAllocated() const {
  if (a_.rows() > 0)
    return true;
  else
    return false;
}

void LinearSystem::Allocate(int nbRows) {
  a_.resize(nbRows, nbRows);
  b_.resize(nbRows);
  x_.resize(nbRows);
  b_.fill(0.);
  x_.fill(0.);
}

void LinearSystem::Clear() {
  a_.setZero();
  b_.setZero();
  x_.setZero();
}

void LinearSystem::ZeroMatrix() {
  a_.setZero();
  b_.setZero();
  x_.setZero();
}

void LinearSystem::ZeroRightHandSide() { b_.fill(0.); }

void LinearSystem::ZeroSolution() { x_.fill(0.); }


int LinearSystem::SystemSolve() {
  if (solver_type_ == EigenCholeskyLLT) {
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(a_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
    x_ = solver.solve(b_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
  } else if (solver_type_ == EigenCholeskyLDLT) {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(a_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
    x_ = solver.solve(b_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
  } else if (solver_type_ == EigenSparseLU) {
    //debug(__FILE__,__LINE__,"Eigen call: solve linear system");

    if (!preprocess_sparsity_pattern_) {
      PreprocessSparsityPattern();
      //if(!preprocess_sparsity_pattern_)return -1;
    }
    if (!factorize_) {
      Factorize();
      //if(!factorize_)return -1;
    }

    x_ = solver_lu.solve(b_);
    if (solver_lu.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
  } else if (solver_type_ == EigenSparseQR) {
    /* Note: maybe another ordering method is better, see Eigen documentation */
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::NaturalOrdering<int> >
        solver;
    solver.compute(a_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
    x_ = solver.solve(b_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
  } else if (solver_type_ == EigenCG) {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
    solver.compute(a_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
    x_ = solver.solve(b_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
  } else if (solver_type_ == EigenCGLeastSquare) {
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    solver.compute(a_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
    x_ = solver.solve(b_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
  } else if (solver_type_ == EigenBiCGSTAB) {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(a_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
    x_ = solver.solve(b_);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      return -1;
    }
  }
  return 1;
}

double LinearSystem::NormInfRightHandSide() const {
  return b_.maxCoeff();
}

double LinearSystem::NormInfSolution() const {
  return x_.maxCoeff();
}

void LinearSystem::AddToMatrix(int row, int col, const double &val) {
  a_.coeffRef(row, col) += val; /* slow ! */
}

void LinearSystem::GetFromMatrix(int row, int col,
                                 double &val) const {
  val = a_.coeff(row, col);
}

void LinearSystem::AddToRightHandSide(int row, const double &val) {
  if ((int)b_.size() <= row) {
    b_.resize(row + 1);
    b_[row] = val;
  } else {
    b_[row] += val;
  }
}

bool LinearSystem::AddToRightHandSide(const std::vector<double> &b) {
  for (int i = 0; i < (int)b.size(); ++i) {
    if (b[i] != 0.) {
      b_[i] = b[i];
    }
  }
  return true;
}

void LinearSystem::GetFromRightHandSide(int row, double &val) const {
  if ((int)b_.size() <= row) val = 0.;
  val = b_[row];
}

void LinearSystem::GetFromSolution(int row, double &val) const {
  if ((int)x_.size() <= row)
    val = 0.;
  else
    val = x_[row];
}

void LinearSystem::AddToSolution(int row, const double &val) {
  if ((int)x_.size() <= row) {
    x_.resize(row + 1);
    x_[row] = val;
  } else {
    x_[row] += val;
  }
}

bool LinearSystem::AddSparseCoefficients(const std::vector<std::vector<size_t>> &columns,
                                         const std::vector<std::vector<double>> &values,
                                         bool first_time) {
  // Msg::Debug("Eigen call: add coefficients");
  std::vector<Eigen::Triplet<double, size_t> > triplets;
  triplets.reserve(values.size());
  if (first_time) {
    for (size_t i = 0; i < columns.size(); ++i) {
      for (size_t k = 0; k < columns[i].size(); ++k) {
        size_t j = columns[i][k];
        double val = values[i][k];
        triplets.emplace_back(Eigen::Triplet<double,size_t>{i, j, val});
      }
    }
    a_.setFromTriplets(triplets.begin(), triplets.end());
  } else {
    for (int i = 0; i < (int)columns.size(); ++i) {
      for (int k = 0; k < (int)columns[i].size(); ++k) {
        int j = (int)columns[i][k];
        double val = values[i][k];
        a_.coeffRef(i, j) += val; /* entry (i,j) should already exist */
      }
    }
  }
  return true;
}

bool LinearSystem::PreprocessSparsityPattern() {
  //debug(__FILE__,__LINE__,"Eigen call: analyse sparse matrix sparsity pattern");
  solver_lu.analyzePattern(a_);
  preprocess_sparsity_pattern_ = true;
  return true;
}

bool LinearSystem::Factorize() {
  //debug(__FILE__,__LINE__,"Eigen call: factorize sparse matrix");
  solver_lu.factorize(a_);
  factorize_ = true;
  return true;
}

bool LinearSystem::Solve(std::vector<double> &x) {
  if (SystemSolve() == -1)return false;
  x.clear();
  x.resize(x_.size());
  for (int i = 0; i < x_.size(); ++i) x[i] = x_[i];
  return true;
}




}