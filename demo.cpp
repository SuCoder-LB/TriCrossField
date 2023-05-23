//
// Created by 苏立彪 on 2023/5/19.
//
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "tcf_linear_system.h"
#include "tcf_basic_types.h"
#include "tcf_interface.h"

#include <wrap/io_trimesh/import.h>

void SaveField(std::string &file_name, std::vector<std::array<double, 3>> &field) {
  FILE *fp = fopen(file_name.c_str(), "w");
  fprintf(fp, "%zu\n", field.size());
  for (auto x : field) {
    fprintf(fp, "%lf %lf %lf\n", x[0], x[1], x[2]);
  }
  fclose(fp);
}

int LoadFeatures(const std::string &file_name, std::vector<std::array<int, 2>> &features) {
  FILE *f = fopen(file_name.c_str(), "rt");
  if (!f) return -1;
  int Num;
  fscanf(f, "%d\n", &Num);
  for (size_t i = 0; i < (size_t) Num; i++) {
    int FIndex, EIndex;
    fscanf(f, "%d,%d\n", &FIndex, &EIndex);
    assert(FIndex >= 0);
    assert(EIndex >= 0);
    assert(EIndex < 4);
    features.push_back({FIndex, EIndex});
  }
  fclose(f);
  return 0;
}


int SaveFeatures(const std::string &file_name, std::vector<std::array<int, 2>> &features) {
  FILE *f = fopen(file_name.c_str(), "w");
  if (!f) return -1;
//  int Num;
  fprintf(f, "%zu\n", features.size());
  for (size_t i = 0; i < features.size(); i++) {
//    int FIndex, EIndex;
    fprintf(f, "%d,%d\n", features[i][0], features[i][1]);
//    assert(FIndex >= 0);
//    assert(EIndex >= 0);
//    assert(EIndex < 4);
//    features.push_back({FIndex, EIndex});
  }
  fclose(f);
  return 0;
}

int main() {

  tri_cross_field::PolyMesh poly_mesh;
  int mask;
  vcg::tri::io::ImporterOBJ<tri_cross_field::PolyMesh>::LoadMask("../gear.obj", mask);
  int err = vcg::tri::io::ImporterOBJ<tri_cross_field::PolyMesh>::Open(poly_mesh, "../gear.obj", mask);

  if (err != 0)return -1;

  std::vector<std::array<int, 2>> features;
  //根据角度求特征边
  double angle = 35.;
  {
    vcg::tri::UpdateTopology<tri_cross_field::PolyMesh>::FaceFace(poly_mesh);
    vcg::tri::UpdateFlags<tri_cross_field::PolyMesh>::VertexBorderFromFaceAdj(poly_mesh);
    vcg::tri::UpdateFlags<tri_cross_field::PolyMesh>::FaceBorderFromFF(poly_mesh);

    for (auto &f : poly_mesh.face)
      for (int j = 0; j < 3; j++)
        f.ClearFaceEdgeS(j);

    if (angle > 0)
      vcg::tri::UpdateFlags<tri_cross_field::PolyMesh>::FaceEdgeSelCrease(
          poly_mesh, (float) vcg::math::ToRad(angle));

    for (auto &f : poly_mesh.face)
      for (int j = 0; j < 3; ++j)
        if (vcg::face::IsBorder(f, j))
          f.SetFaceEdgeS(j);

    for (auto &f : poly_mesh.face)
      for (int j = 0; j < 3; ++j)
        if (!vcg::face::IsManifold(f, j))
          f.SetFaceEdgeS(j);

    for (int i = 0; i < poly_mesh.fn; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (poly_mesh.face[i].IsFaceEdgeS(j)) {
          features.push_back({i, j});
        }
      }
    }

    SaveFeatures("../gear.feature",features);
  }


  std::vector<std::array<double, 3>> field;

  auto linear_system = new tri_cross_field::LinearSystem();

  tri_cross_field::CrossFieldParameter cfp;

  tri_cross_field::ComputeCrossField(poly_mesh, features, linear_system, field, cfp);

  std::string file = "../gear.rosy";
  SaveField(file, field);

  return 0;
}