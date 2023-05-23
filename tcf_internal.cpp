#ifndef QUAD_LAYOUT__CROSS_FIELD_H_
#define QUAD_LAYOUT__CROSS_FIELD_H_

/* System includes */
#include <vector>
#include <array>
#include <unordered_map>
#include <cstdint>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>

#include <vcg/complex/complex.h>
#include <vcg/math/base.h>
#include <vcg/space/point2.h>
#include <vcg/space/point3.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/complex/algorithms/update/topology.h>

#include "tcf_basic_types.h"
//#include "log.h"

#include "tcf_linear_system.h"

using std::vector;
using std::unordered_map;

namespace Geometry {

const double EPS = 1.e-14;

using id = uint32_t;
using sid = int64_t;
using id2 = std::array<id, 2>;
using id3 = std::array<id, 3>;
using sid3 = std::array<sid, 3>;
const id NO_ID = (id) -1;
const sid NO_SID = (sid) NO_ID;

using vec2 = std::array<double, 2>;
using vec3 = std::array<double, 3>;

/* scalar utils */
inline double clamp(double x, double lower, double upper) { return std::min(upper, std::max(x, lower)); }
inline double maxAbs(double x, double y, double z) { return std::max(std::abs(x), std::max(std::abs(y), std::abs(z))); }

/* vec3 math */
inline vec3 operator-(const vec3 &a, const vec3 &b) { return {{a[0] - b[0], a[1] - b[1], a[2] - b[2]}}; }
inline vec3 operator+(const vec3 &a, const vec3 &b) { return {{a[0] + b[0], a[1] + b[1], a[2] + b[2]}}; }
inline vec3 operator*(const double &a, const vec3 &b) { return {{a * b[0], a * b[1], a * b[2]}}; }
inline vec3 operator*(const vec3 &a, const double &b) { return {{a[0] * b, a[1] * b, a[2] * b}}; }
inline double dot(const vec3 &a, const vec3 &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline double length2(const vec3 &a) { return dot(a, a); }
inline double length(const vec3 &a) { return sqrt(length2(a)); }
inline double maxAbs(const vec3 &a) { return maxAbs(a[0], a[1], a[2]); }
inline vec3 cross(const vec3 &a, const vec3 &b) {
  return {{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}};
}
inline void normalizeFast(vec3 &a) { a = 1. / sqrt(length2(a)) * a; } /* no check, not safe, not accurate */
inline void normalize(vec3 &a) {
  double amp = maxAbs(a);
  if (amp == 0.) { return; }
  a = amp * a;
  normalizeFast(a);
}
inline double angleVectorsAlreadyNormalized(const vec3 &a, const vec3 &b) { return acos(clamp(dot(a, b), -1., 1.)); }
inline double angleVectors(vec3 a, vec3 b);

inline double angle(const vec3 &a, const vec3 &b) {
  double cosTheta = dot(a, b);
  double sinTheta = cross(a, b)[2];
  return atan2(sinTheta, cosTheta);
}

/* vec2 math */
inline vec2 operator-(const vec2 &a, const vec2 &b) { return {{a[0] - b[0], a[1] - b[1]}}; }
inline vec2 operator+(const vec2 &a, const vec2 &b) { return {{a[0] + b[0], a[1] + b[1]}}; }
inline vec2 operator*(const double &a, const vec2 &b) { return {{a * b[0], a * b[1]}}; }
inline vec2 operator*(const vec2 &a, const double &b) { return {{a[0] * b, a[1] * b}}; }
inline double dot(const vec2 &a, const vec2 &b) { return a[0] * b[0] + a[1] * b[1]; }
inline double length2(const vec2 &a) { return dot(a, a); }
inline double length(const vec2 &a) { return sqrt(length2(a)); }

/* functions */

inline double bboxDiag(const std::vector<vec3> &points);
inline vec3 triangleNormal(const vec3 &p0, const vec3 &p1, const vec3 &p2); /* normalized */
inline double triangleArea(const vec2 &a, const vec2 &b, const vec2 &c);
inline double triangleArea(const vec3 &p0, const vec3 &p1, const vec3 &p2);


/**********************/
/* inline definitions */
/**********************/

double angleVectors(vec3 a, vec3 b) {
  double ma = maxAbs(a);
  double mb = maxAbs(b);
  if (ma == 0. || mb == 0.) return DBL_MAX;
  a = ma * a;
  b = mb * b;
  normalizeFast(a);
  normalizeFast(b);
  return angleVectorsAlreadyNormalized(a, b);
}

vec3 triangleNormal(const vec3 &p0, const vec3 &p1, const vec3 &p2) {
  vec3 N = cross(p2 - p0, p1 - p0);
  if (maxAbs(N) == 0.) return {{0., 0., 0.}};
  normalize(N);
  return N;
}

double bboxDiag(const std::vector<vec3> &points) {
  vec3 mi = {DBL_MAX, DBL_MAX, DBL_MAX};
  vec3 ma = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
  for (size_t i = 0; i < points.size(); ++i)
    for (size_t d = 0; d < 3; ++d) {
      mi[d] = std::min(points[i][d], mi[d]);
      ma[d] = std::max(points[i][d], ma[d]);
    }
  return length(ma - mi);
}

double triangleArea(const vec2 &a, const vec2 &b, const vec2 &c) {
  return .5 * ((b[1] - a[1]) * (b[0] + a[0]) + (c[1] - b[1]) * (c[0] + b[0]) + (a[1] - c[1]) * (a[0] + c[0]));
}

double triangleArea(const vec3 &p0, const vec3 &p1, const vec3 &p2) {
  vec3 N = cross(p2 - p0, p1 - p0);
  return length(N) / 2.;
}

}

//template<class Mesh>
namespace tri_cross_field {

using namespace Geometry;

struct IJV {
  id2 ij;
  double val;
  bool operator<(const IJV &b) const {
    return ij[0] < b.ij[0] || (ij[0] == b.ij[0] && ij[1] < b.ij[1]);
  }
};

struct IV {
  id i;
  double val;
  bool operator<(const IV &b) const { return i < b.i; }
};

template<class T1, class T2>
T2 sort_unique_with_perm(
    const std::vector<T1> &in,
    std::vector<T1> &uniques,
    std::vector<T2> &old2new) {

  std::vector<T2> ids(in.size());
  for (T2 k = 0; k != in.size(); ++k) ids[k] = k;

  std::sort(ids.begin(), ids.end(),
            [&in](const T2 &a, const T2 &b) { return (in[a] < in[b]); }
  );

  uniques.resize(in.size());
  old2new.resize(in.size());
  for (T2 k = 0; k != in.size(); ++k) uniques[k] = in[k];

  std::sort(uniques.begin(), uniques.end());
  uniques.erase(std::unique(uniques.begin(), uniques.end()),
                uniques.end());
  T2 ic = 0; // indice current
  T2 ir = 0; // indice representant
  T2 cur_rep = 0; // indice of current representant
  while (ic < in.size()) {
    ic = ir;
    while (ic < in.size() && in[ids[ic]] == in[ids[ir]]) {
      old2new[ids[ic]] = cur_rep;
      ++ic;
    }
    ir = ic;
    ++cur_rep;
  }
  return (T2) uniques.size();
}

bool compute_triangle_adjacencies(const vector<id3> &triangles,
                                  vector<sid3> &triangle_neighbors,
                                  vector<vector<id>> &nm_triangle_neighbors,
                                  vector<id2> &uIEdges, vector<id> &old2IEdge,
                                  vector<vector<id>> &uIEdgeToOld) {
  triangle_neighbors.clear();
  triangle_neighbors.resize(triangles.size(), {NO_ID, NO_ID, NO_ID});
  nm_triangle_neighbors.clear();

  constexpr size_t NBF = 3;

  /* Store element 'faces', with duplicates for further 'equality test' */
  std::vector<id2> faces;
  for (size_t i = 0; i < triangles.size(); ++i) {
    for (size_t lf = 0; lf < NBF; ++lf) {
      id2 face = {triangles[i][lf], triangles[i][(lf + 1) % NBF]};
      std::sort(face.begin(), face.end());
      faces.push_back(face);
    }
  }

  /* Reduce 'duplicated faces' to 'unique faces', keeping the 'old2new'
   * mapping */
  size_t nbUniques = sort_unique_with_perm(faces, uIEdges, old2IEdge);

  /* Build the 'unique face to elements' mapping */
  uIEdgeToOld.resize(nbUniques);
  for (size_t i = 0; i < triangles.size(); ++i) {
    for (size_t lf = 0; lf < NBF; ++lf) {
      id facePos = NBF * i + lf;
      uIEdgeToOld[old2IEdge[facePos]].push_back(facePos);
    }
  }

  /* Loop over 'unique face to elements' and set the element adjacencies */
  const id2 NO_FACE = {NO_ID, NO_ID};
  for (size_t i = 0; i < nbUniques; ++i) {
    if (uIEdges[i] == NO_FACE) continue;
    if (uIEdges[i][0] == NO_ID) return false;
    if (uIEdges[i][1] == NO_ID) return false;
    if (uIEdgeToOld[i].size() == 1) { /* boundary */
      id eltId = uIEdgeToOld[i][0] / NBF;
      id lf = uIEdgeToOld[i][0] % NBF;
      triangle_neighbors[eltId][lf] = NO_ID;
    } else if (uIEdgeToOld[i].size() == 2) { /* regular face */
      id e1 = uIEdgeToOld[i][0] / NBF;
      id lf1 = uIEdgeToOld[i][0] % NBF;
      id e2 = uIEdgeToOld[i][1] / NBF;
      id lf2 = uIEdgeToOld[i][1] % NBF;
      triangle_neighbors[e1][lf1] = NBF * e2 + lf2;
      triangle_neighbors[e2][lf2] = NBF * e1 + lf1;
    } else if (uIEdgeToOld[i].size() > 2) { /* non manifold face */
      for (size_t j = 0; j < uIEdgeToOld[i].size(); ++j) {
        id e = uIEdgeToOld[i][j] / NBF;
        id lf = uIEdgeToOld[i][j] % NBF;
        std::vector<id> neighs;
        for (size_t k = 0; k < uIEdgeToOld[i].size(); ++k) {
          if (uIEdgeToOld[i][k] != uIEdgeToOld[i][j]) {
            neighs.push_back(uIEdgeToOld[i][k]);
          }
        }
        neighs.shrink_to_fit();
        id pos = nm_triangle_neighbors.size();
        nm_triangle_neighbors.push_back(neighs);
        triangle_neighbors[e][lf] = -((sid) pos + 1);
      }
    }
  }

  /* Reduce memory consumption */
  triangle_neighbors.shrink_to_fit();
  nm_triangle_neighbors.shrink_to_fit();

  return true;
}

/* two unknowns (x_2i, x_2i+1) for each edge */
bool stiffness_coefficient(const std::vector<vec3> &points,
                           const std::vector<std::array<id, 3> > &triangles,
                           id e, const vector<id2> &uIEdges,
                           const vector<id> &old2IEdge,
                           const vector<vector<id>> &uIEdgeToOld,
                           vector<IV> &iv, vector<IJV> &ijv) {
  if (uIEdgeToOld[e].size() != 2) {
//    error("assembly, edge {}: uIEdgeToOld[e].size() = {}", e,
//          uIEdgeToOld[e].size());
    return false;
  }
  /* edge vertices */
  id v1 = uIEdges[e][0];
  id v2 = uIEdges[e][1];
  vec3 e_x = points[v2] - points[v1];
  double lenr = length(e_x);
  if (lenr < EPS) {
//    error("edge too small: v1={}, v2={}, length = {}", v1, v2, lenr);
    return false;
  }
  e_x = (1. / lenr) * e_x;
  /* four neighbor edges */
  id bvars[4] = {NO_ID, NO_ID, NO_ID, NO_ID};
  double alpha[4] = {0., 0., 0., 0.};
  double CR_weight[4] = {-1. / 4, -1. / 4, -1. / 4, -1. / 4};
  vec3 prevN;
  for (size_t s = 0; s < 2; ++s) { /* two triangles */
    id oe = uIEdgeToOld[e][s];
    id t = oe / 3;
    id le = oe % 3;
    vec3 N = triangleNormal(points[triangles[t][0]], points[triangles[t][1]],
                            points[triangles[t][2]]);
    if (s == 1 && dot(prevN, N) < 0.) N = -1. * N;
    prevN = N;
    // N = {0,0,1};
    vec3 e_y = cross(N, e_x);
    if (maxAbs(e_y) == 0.) {
//      error("length(e_y) = {}", length(e_y));
      return false;
    }
    normalize(e_y);

    for (size_t k = 0; k < 2; ++k) { /* two other edges in triangle t */
      id aoe = 3 * t + (le + 1 + k) % 3;
      id ae = old2IEdge[aoe];
      id2 aedge = uIEdges[ae];

      bvars[2 * s + k] = ae;
      vec3 edg = points[aedge[1]] - points[aedge[0]];
      double len = length(edg);
      if (len < EPS) {
//        error("edge too small: t={}, k = {}, edge={}, length = {}", t,
//              k, aedge, len);
        return false;
      }
      edg = (1. / len) * edg;

      /* 360deg angle (used for the angle rotation) */
      double cx = dot(edg, e_x);
      double cy = dot(edg, e_y);
      alpha[2 * s + k] = atan2(cy, cx);
      if (alpha[2 * s + k] < 0.) alpha[2 * s + k] += 2. * M_PI;

      /* 180deg edge-edge angle (used for the Crouzeix Raviart) */
      double agl = 0.;
      if (aedge[0] == v1) { agl = angleVectorsAlreadyNormalized(edg, e_x); }
      else if (aedge[1] == v1) {
        agl = angleVectorsAlreadyNormalized(edg, -1. * e_x);
      } else if (aedge[0] == v2) {
        agl = angleVectorsAlreadyNormalized(-1. * edg, e_x);
      } else if (aedge[1] == v2) {
        agl = angleVectorsAlreadyNormalized(-1. * edg, -1. * e_x);
      } else {
//        error("should not happen");
        return false;
      }
      CR_weight[2 * s + k] = -2. / tan(agl);
    }
  }
  /* compute coefficients */
  double isum =
      -1. * (CR_weight[0] + CR_weight[1] + CR_weight[2] + CR_weight[3]);
  for (size_t k = 0; k < 4; ++k) CR_weight[k] = CR_weight[k] / isum;

  iv.clear();
  ijv.clear();
  id x_i = 2 * e;
  id y_i = 2 * e + 1;
  iv.push_back({x_i, 1.});
  iv.push_back({y_i, 1.});
  double Nd = 4;
  for (size_t j = 0; j < 4; ++j) {
    id x_j = 2 * bvars[j];
    id y_j = 2 * bvars[j] + 1;
    ijv.push_back({{x_i, x_j}, CR_weight[j] * cos(Nd * alpha[j])});
    ijv.push_back({{x_i, y_j}, CR_weight[j] * -1. * sin(Nd * alpha[j])});
    ijv.push_back({{y_i, x_j}, CR_weight[j] * sin(Nd * alpha[j])});
    ijv.push_back({{y_i, y_j}, CR_weight[j] * cos(Nd * alpha[j])});
  }

  return true;
}

bool prepare_system(const vector<IV> &K_diag, const vector<IJV> &K_coefs,
                    vector<vector<size_t>> &columns,
                    vector<vector<double>> &values) {
  vector<IJV> coefs = K_coefs;
  coefs.resize(coefs.size() + K_diag.size());
  for (size_t i = 0; i < K_diag.size(); ++i) {
    coefs[K_coefs.size() + i] = {{id(K_diag[i].i), id(K_diag[i].i)},
                                 K_diag[i].val};
  }
  std::sort(coefs.begin(), coefs.end());

  size_t cur_i = coefs[0].ij[0];
  size_t cur_j = coefs[0].ij[1];
  double acc_val = coefs[0].val;
  for (size_t k = 1; k < coefs.size(); ++k) {
    id i = coefs[k].ij[0];
    id j = coefs[k].ij[1];
    double v = coefs[k].val;
    if (i != cur_i) {
      if (std::abs(acc_val) > EPS) {
        columns[cur_i].push_back(cur_j);
        values[cur_i].push_back(acc_val);
      }
      cur_i = i;
      acc_val = v;
      cur_j = j;
    } else if (j != cur_j) {
      if (std::abs(acc_val) > EPS) {
        columns[cur_i].push_back(cur_j);
        values[cur_i].push_back(acc_val);
      }
      cur_i = i;
      acc_val = v;
      cur_j = j;
    } else {
      acc_val += v;
    }
  }
  if (std::abs(acc_val) > EPS) {
    columns[cur_i].push_back(cur_j);
    values[cur_i].push_back(acc_val);
  }

  return true;
}

bool compute_cross_field_with_multilevel_diffusion(
    const std::vector<std::array<double, 3> > &points,
    const std::vector<std::array<id, 2> > &lines,
    const std::vector<std::array<id, 3> > &triangles,
    std::vector<std::array<double, 3> > &triEdgeTheta,
    std::vector<std::array<double, 3>> &theta_norm,
    const CrossFieldParameter &cross_field_parameter,
    LinearSystem *linear_system) {

  int nbDiffusionLevels = cross_field_parameter.nbDiffusionLevels;

  /* Build unique edges and association with adjacent triangles,
   * including the non-manifold cases */
  vector<id2> uIEdges;
  vector<id> old2IEdge;
  vector<vector<id> > uIEdgeToOld;
  std::vector<sid3> triangle_neighbors;
  std::vector<std::vector<id> > nm_triangle_neighbors;
  bool oka = compute_triangle_adjacencies(triangles, triangle_neighbors,
                                          nm_triangle_neighbors, uIEdges,
                                          old2IEdge, uIEdgeToOld);
  if (!oka) {
//    error("failed to compute mesh adjacencies");
    return false;
  }

  /* Bbox diagonal is used later to specify the diffusion length */
  double diag = bboxDiag(points);

  if (uIEdges.size() == 0) {
//    error("no internal edges");
    return false;
  }

  /* System unknowns: cos(Nt),sin(Nt) for each edge */
  vector<double> x(2 * uIEdges.size(), 0.);

  /* Initial Dirichlet boundary conditions
   * alignment of crosses with edges (relative angle = 0)
   * theta_e = 0 => (cos4t/sin4t) = (1,0) */
  size_t nbc = 0;
  vector<bool> dirichletEdge(uIEdges.size(), false);
  vector<vec2> dirichletValue(uIEdges.size(), {1., 0.});
  for (size_t e = 0; e < uIEdges.size(); ++e)
    if (uIEdgeToOld[e].size() != 2) {
      dirichletEdge[e] = true;
      dirichletValue[e] = {1., 0.};
      nbc += 1;
    }
  for (size_t l = 0; l < lines.size(); ++l) {
    /* mark the lines as boundary conditions */
    id2 edge = {lines[l][0], lines[l][1]};
    std::sort(edge.begin(), edge.end());
    auto it = std::find(uIEdges.begin(), uIEdges.end(), edge);
    if (it != uIEdges.end()) {
      id e = id(it - uIEdges.begin());
      dirichletEdge[e] = true;
      dirichletValue[e] = {1., 0.};
      nbc += 1;
    }
  }

  vector<IV> K_diag;
  vector<IJV> K_coefs;
  vector<double> rhs(2 * uIEdges.size(), 0.);
  vector<double> Mass(2 * uIEdges.size(),
                      1.); /* diagonal for Crouzeix-Raviart */
  for (size_t e = 0; e < uIEdges.size(); ++e) {
    if (!dirichletEdge[e]) {
      vector<IV> iv;
      vector<IJV> ijv;
      bool oks = stiffness_coefficient(points, triangles, e, uIEdges,
                                       old2IEdge, uIEdgeToOld, iv, ijv);
      if (!oks) {
//        error(
//            "failed to compute stiffness matrix coefficients for e = {}", e);
        return false;
      }
      for (size_t i = 0; i < iv.size(); ++i) K_diag.push_back(iv[i]);
      for (size_t i = 0; i < ijv.size(); ++i) K_coefs.push_back(ijv[i]);
      id t1 = uIEdgeToOld[e][0] / 3;
      id t2 = uIEdgeToOld[e][1] / 3;
      double area1 =
          triangleArea(points[triangles[t1][0]], points[triangles[t1][1]],
                       points[triangles[t1][2]]);
      double area2 =
          triangleArea(points[triangles[t2][0]], points[triangles[t2][1]],
                       points[triangles[t2][2]]);
      Mass[2 * e + 0] = 1. / 3 * (area1 + area2);
      Mass[2 * e + 1] = 1. / 3 * (area1 + area2);
    } else { /* Dirichlet BC */
      K_diag.push_back({id(2 * e + 0), 1.});
      rhs[2 * e + 0] = dirichletValue[e][0];
      K_diag.push_back({id(2 * e + 1), 1.});
      rhs[2 * e + 1] = dirichletValue[e][1];
    }
  }

  double eavg = 0.;
  double emin = DBL_MAX;
  double emax = -DBL_MAX;
  for (size_t e = 0; e < uIEdges.size(); ++e) {
    double len = length(points[uIEdges[e][1]] - points[uIEdges[e][0]]);
    eavg += len;
    emin = std::min(len, emin);
    emax = std::max(len, emax);
  }
  eavg /= uIEdges.size();


  /* prepare system */
  vector<vector<size_t> > K_columns(2 * uIEdges.size());
  vector<vector<double> > K_values(2 * uIEdges.size());
  {
    bool okp = prepare_system(K_diag, K_coefs, K_columns, K_values);
    if (!okp) {
//      error("failed to prepare system");
      return false;
    }
  }

  double dtInitial = (0.1 * diag) * (0.1 * diag);
  double dtFinal = (3. * emin) * (3. * emin);

  double wti = clock();
  for (size_t e = 0; e < uIEdges.size(); ++e) {
    x[2 * e + 0] = dirichletValue[e][0];
    x[2 * e + 1] = dirichletValue[e][1];
  }

  vector<double> steps;

  if (nbDiffusionLevels > 1) {
    for (size_t i = 0; i < nbDiffusionLevels;
         ++i) { /* resolution transition */
      double dt = dtInitial + (dtFinal - dtInitial) * double(i) /
          double(nbDiffusionLevels - 1);
      steps.push_back(dt);
    }
  } else {
    steps.push_back(dtFinal);
  }

  {
    /* Initialize system (I/dt + M^-1 * L) x_(i+1) = 1/dt * x_i     (Ax = B)
     */
    vector<vector<size_t> > Acol = K_columns;
    vector<vector<double> > Aval_add = K_values; /* to get sparsity pattern */
    for (size_t i = 0; i < Aval_add.size(); ++i)
      std::fill(Aval_add[i].begin(), Aval_add[i].end(), 0.);
    vector<double> B = x;
    vector<double> norms(uIEdges.size(), 0.);
    vector<double> prevNorms = norms;

    auto &solver = *linear_system;
    solver.Allocate(2 * uIEdges.size());

//    CrossFieldLinearSystem solver(2 * uIEdges.size());

    vector<double> diag_sum(2 * uIEdges.size(), 0.);
    for (size_t i = 0; i < Acol.size(); ++i) {
      if (!dirichletEdge[i / 2]) {
        for (size_t j = 0; j < Acol[i].size(); ++j) {
          Aval_add[i][j] = 1. / Mass[i] * K_values[i][j];
        }
      } else {
        Aval_add[i] = {1.};
      }
    }
    solver.AddSparseCoefficients(Acol, Aval_add, true);
    for (size_t i = 0; i < Aval_add.size(); ++i) Aval_add[i].clear();
    for (size_t i = 0; i < Acol.size(); ++i) Acol[i].clear();
    bool okp = solver.PreprocessSparsityPattern();
    if (!okp) {
//      error("linear solver analysis failed");
      return false;
    }

    /* Loop over the changing timesteps */
    double prev_dt = DBL_MAX;
//    steps.resize(1);
//    steps[0] = DBL_MAX;
    for (int iter = 0; iter < steps.size(); ++iter) {
      if (iter > 0 && steps[iter] > prev_dt) continue;
      double dt = steps[iter];
      prev_dt = dt;

      /* Update LHS matrix with the new timestep */
      for (size_t i = 0; i < Acol.size(); ++i) {
        if (!dirichletEdge[i / 2]) {
          Acol[i] = {i};
          Aval_add[i] = {-diag_sum[i] + 1. / dt};
          diag_sum[i] = 1. / dt;
        }
      }
      bool oku = solver.AddSparseCoefficients(Acol, Aval_add, false);
      if (!oku) {
//        error("failed to update linear system");
        return false;
      }
      solver.Factorize();

      /* Loop at fixed time step */
      size_t subiter_max = cross_field_parameter.subiter;
//      constexpr size_t subiter_max = 1;
      for (size_t subiter = 0; subiter < subiter_max; ++subiter) {
        prevNorms = norms;

//        dbg(iter,subiter);
        /* Update RHS */
        B = x;
        for (size_t i = 0; i < B.size(); ++i)
          if (!dirichletEdge[i / 2]) B[i] /= dt;
        solver.AddToRightHandSide(B);

        bool oks = solver.Solve(x);
        if (!oks) {
//          error("failed to solve linear system");
          return false;
        }

//        if (iter == 9 && subiter == 24)break;

        /* Normalize cross field and gather norm stats */
        double nmi = DBL_MAX;
        double nma = -DBL_MAX;
        //su:test
        for (size_t e = 0; e < uIEdges.size(); ++e) {
          norms[e] =
              std::sqrt(x[2 * e] * x[2 * e] + x[2 * e + 1] * x[2 * e + 1]);
          nmi = std::min(nmi, norms[e]);
          nma = std::max(nma, norms[e]);
          if (!dirichletEdge[e] && norms[e] > EPS) {
            x[2 * e + 0] /= norms[e];
            x[2 * e + 1] /= norms[e];
          }
        }
        // Msg::Info("            |   norm, min = {}, max = {}", nmi, nma);
        const double EPS_NORM = 1.e-1;
        if (nma > 1. + EPS_NORM) {
          steps[iter] /= 10;
          dt = steps[iter];
          iter -= 1;
          break;
        }
        if (subiter > 0 || iter > 0) {
          double linf = 0.;
          for (size_t i = 0; i < norms.size(); ++i)
            if (!dirichletEdge[i]) {
              linf = std::max(linf, norms[i] - prevNorms[i]);
            }
          if (linf < 1.e-3) break;
        } else {

        }
      }
    }
  }
  double et = (clock() - wti) / CLOCKS_PER_SEC;

  { /* Export solution */
    triEdgeTheta.resize(triangles.size());
    theta_norm.resize(triangles.size());
    for (size_t t = 0; t < triangles.size(); ++t) {
      for (size_t le = 0; le < 3; ++le) {
        id e = old2IEdge[3 * t + le];
        double len =
            std::sqrt(x[2 * e] * x[2 * e] + x[2 * e + 1] * x[2 * e + 1]);
        double theta = (len > EPS) ?
                       1. / 4 * atan2(x[2 * e + 1] / len, x[2 * e] / len) : 0.;
        theta_norm[t][le] = len;
        triEdgeTheta[t][le] = theta;
      }
    }
  }

  return true;
}

static void Save4ROSY(vector<vec3> &dir, const std::string &path) {
  FILE *f = fopen(path.c_str(), "wt");
  fprintf(f, "%d\n", (int) dir.size());
  fprintf(f, "4\n");
  for (unsigned int i = 0; i < dir.size(); i++) {
    float dirX = (float) dir[i][0];
    float dirY = (float) dir[i][1];
    float dirZ = (float) dir[i][2];
    fprintf(f, "%f %f %f \n", dirX, dirY, dirZ);
  }
  fclose(f);
}

inline vec3 crouzeix_raviart_interpolation(vec3 lambda,
                                           vec3 edge_values[3]) {
  double x = lambda[1];
  double y = lambda[2];
  vec3 shape = {1.0 - 2.0 * y, -1.0 + 2.0 * (x + y), 1.0 - 2.0 * x};
  return shape[0] * edge_values[0] + shape[1] * edge_values[1] +
      shape[2] * edge_values[2];
}

void ComputeCrossField(PolyMesh &tri_mesh, std::vector<std::array<int, 2>> &features,
                       LinearSystem *linear_system, std::vector<std::array<double, 3>> &field,
                       const CrossFieldParameter &cfp) {
  if (!linear_system) {
    printf("没有求解器\n");
    return;
  }

  vector<vec3> dpoints;
  for (auto p : tri_mesh.vert)dpoints.push_back({p.P()[0], p.P()[1], p.P()[2]});

  std::vector<std::array<uint32_t, 2>> feature;
  feature.reserve(features.size());
  for (auto e : features) {
    feature.push_back({(uint32_t) vcg::tri::Index(tri_mesh, tri_mesh.face[e[0]].cV(e[1])),
                       (uint32_t) vcg::tri::Index(tri_mesh, tri_mesh.face[e[0]].cV((e[1] + 1) % 3))});
  }

//  for (size_t i = 0; i < tri_mesh.face.size(); i++) {
//    for (int j = 0; j < 3; j++) {
//      if (!tri_mesh.face[i].IsFaceEdgeS(j))continue;
//      feature.push_back({(uint32_t) vcg::tri::Index(tri_mesh, tri_mesh.face[i].cV(j)),
//                         (uint32_t) vcg::tri::Index(tri_mesh, tri_mesh.face[i].cV((j + 1) % 3)),});
//    }
//  }

  std::vector<std::array<uint32_t, 3>> triangles;
  for (int i = 0; i < tri_mesh.fn; ++i) {
    triangles.push_back({(uint32_t) vcg::tri::Index(tri_mesh, tri_mesh.face[i].cV(0)),
                         (uint32_t) vcg::tri::Index(tri_mesh, tri_mesh.face[i].cV(1)),
                         (uint32_t) vcg::tri::Index(tri_mesh, tri_mesh.face[i].cV(2))});
  }

  for (auto &x : feature)std::sort(x.begin(), x.end());
  std::sort(feature.begin(), feature.end());
  feature.erase(std::unique(feature.begin(), feature.end()), feature.end());

  std::vector<std::array<double, 3>> triEdgeTheta;
  std::vector<std::array<double, 3>> theta_norm;

  compute_cross_field_with_multilevel_diffusion(dpoints, feature, triangles, triEdgeTheta, theta_norm,
                                                cfp, linear_system);

  std::vector<std::array<double, 3>> vert_dir_norm(dpoints.size());
  std::vector<vec3> vertex_dir(dpoints.size());
  std::vector<std::array<double, 3>> vertex_dir_norm(dpoints.size());
  std::set<id> feature_points;
  for (auto line : feature)feature_points.insert(line[0]), feature_points.insert(line[1]);

  for (size_t f = 0; f < triangles.size(); ++f) {
    id t = f;

    vec3 N = triangleNormal(dpoints[triangles[t][0]],
                            dpoints[triangles[t][1]], dpoints[triangles[t][2]]);
    /* Compute one branch at triangle vertices */
    vec3 refCross = {0., 0., 0.};
    vec3 avgCross = {0., 0., 0.};
    vec3 lifted_dirs[3];
    for (size_t le = 0; le < 3; ++le) {
      /* Get cross angle */
      std::array<id, 2> edge = {triangles[t][le], triangles[t][(le + 1) % 3]};
      if (edge[0] > edge[1])
        std::reverse(edge.begin(), edge.end());

      double A = triEdgeTheta[f][le];

      /* Local reference frame */
      vec3 tgt = dpoints[edge[1]] - dpoints[edge[0]];
      if (length(tgt) < 1.e-16) {
//        warning("Edge (tri={},le={}), length = {}", t, le,
//                length(tgt));
        if (length(tgt) == 0.) { assert(0); }
      }
      normalize(tgt);
      vec3 tgt2 = cross(N, tgt);
      normalize(tgt2);

      vec3 cross1 = tgt * cos(A) + tgt2 * sin(A);
      normalize(cross1);

      vec3 cross2 = cross(N, cross1);
      normalize(cross2);

      if (le == 0) {
        refCross = cross1;
        avgCross = avgCross + cross1;
        lifted_dirs[le] = refCross;
      } else {
        vec3 closest = {0., 0., 0.};
        double dotmax = -DBL_MAX;
        for (int k = 0; k < 4; ++k) {
          double agl = A + double(k) / double(4) * 2. * M_PI;
          vec3 candidate = tgt * cos(agl) + tgt2 * sin(agl);
          normalize(candidate);
          if (dot(candidate, refCross) > dotmax) {
            closest = candidate;
            dotmax = dot(candidate, refCross);
          }
        }

        lifted_dirs[le] = closest;
        avgCross = avgCross + closest;
      }
    }

    for (size_t lv = 0; lv < 3; ++lv) {
      vec3 lambda = {0., 0., 0.};
      vec3 lifted_norm[3];
      lifted_norm[0][0] = theta_norm[t][0];
      lifted_norm[1][0] = theta_norm[t][1];
      lifted_norm[2][0] = theta_norm[t][2];
      lambda[lv] = 1.;
      vec3 dir = crouzeix_raviart_interpolation(lambda, lifted_dirs);
      double dir_norm = crouzeix_raviart_interpolation(lambda, lifted_norm)[0];
      vert_dir_norm[triangles[f][lv]][0] = dir_norm;
      if (length2(dir) != 0.) normalize(dir);
      for (size_t d = 0; d < 3; ++d) {
        vertex_dir[triangles[f][lv]][d] = dir[d];
        if (feature_points.count(triangles[f][lv]))vertex_dir_norm[triangles[f][lv]][0] = 1.;
      }
    }

    vec3 dir = lifted_dirs[0] + lifted_dirs[1] + lifted_dirs[2];
    normalize(dir);
    triEdgeTheta[t] = {dir[0], dir[1], dir[2]};
    theta_norm[t][0] = (theta_norm[t][0] + theta_norm[t][1] + theta_norm[t][2]) / 3;
  }
  field = triEdgeTheta;
}

};

#endif //QUAD_LAYOUT__CROSS_FIELD_H_
