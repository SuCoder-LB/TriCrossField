//
// Created by 苏立彪 on 2023/5/19.
//

#ifndef TRICROSSFIELD__TCF_LINEAR_SYSTEM_H_
#define TRICROSSFIELD__TCF_LINEAR_SYSTEM_H_

#include <vector>

//
//class LinearSystem{
// public:
//  LinearSystem()=default;
//
//  bool IsAllocated() const;
//  void Allocate(int nbRows);
//  void Clear();
//  void ZeroMatrix();
//
//  void ZeroRightHandSide();
//  void ZeroSolution();
//  int SystemSolve();
//  double NormInfRightHandSide() const;
//  double NormInfSolution() const;
//
//  void AddToMatrix(int row, int col, const double &val);
//  void GetFromMatrix(int row, int col, double &val) const;
//  void AddToRightHandSide(int row, const double &val);
//  bool AddToRightHandSide(const std::vector<double> &b);
//  void GetFromRightHandSide(int row, double &val) const;
//  void GetFromSolution(int row, double &val) const;
//  void AddToSolution(int row, const double &val);
//
//  bool AddSparseCoefficients(
//      const std::vector<std::vector<size_t>> &columns,
//      const std::vector<std::vector<double>> &values,
//      bool first_time = false);
//
//  bool PreprocessSparsityPattern();
//  bool Factorize();
//
//  bool Solve(std::vector<double> &x);
//
//};

#include <Eigen/Sparse>

namespace tri_cross_field {


enum linearSystemEigenSolver {
  EigenCholeskyLLT,
  EigenCholeskyLDLT,
  EigenSparseLU,
  EigenSparseQR,
  EigenCG,
  EigenCGLeastSquare,
  EigenBiCGSTAB
};

class  LinearSystem {
 private:
  Eigen::VectorXd x_;
  Eigen::VectorXd b_;
  Eigen::SparseMatrix<double> a_;
  linearSystemEigenSolver solver_type_;

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver_lu;

  //for lu
  bool factorize_;
  bool preprocess_sparsity_pattern_;

 public:
  LinearSystem();

  bool IsAllocated() const;
  void Allocate(int nbRows);
  void Clear();
  void ZeroMatrix();

  void ZeroRightHandSide();
  void ZeroSolution();
  int SystemSolve();
  double NormInfRightHandSide() const;
  double NormInfSolution() const;

  void AddToMatrix(int row, int col, const double &val);
  void GetFromMatrix(int row, int col, double &val) const;
  void AddToRightHandSide(int row, const double &val);
  bool AddToRightHandSide(const std::vector<double> &b);
  void GetFromRightHandSide(int row, double &val) const;
  void GetFromSolution(int row, double &val) const;
  void AddToSolution(int row, const double &val);

  bool AddSparseCoefficients(
      const std::vector<std::vector<size_t>> &columns,
      const std::vector<std::vector<double>> &values,
      bool first_time = false);

  bool PreprocessSparsityPattern();
  bool Factorize();

  bool Solve(std::vector<double> &x);
};




}

#endif //TRICROSSFIELD__TCF_LINEAR_SYSTEM_H_
