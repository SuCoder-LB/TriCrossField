//
// Created by 苏立彪 on 2023/5/11.
//

#ifndef TRICROSSFIELD__TCF_BASIC_TYPES_H_
#define TRICROSSFIELD__TCF_BASIC_TYPES_H_

#include <unordered_map>
#include <vector>
#include <array>


#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/clean.h>


namespace tri_cross_field{

struct CrossFieldParameter {
  int nbDiffusionLevels = 10;
  double thresholdNormConvergence = 1.e-2;
  int subiter = 25;
};

class PolyVertex;
class PolyFace;
class PolyEdge;

struct MyPolyTypes : public vcg::UsedTypes<
    vcg::Use<PolyVertex>::AsVertexType,
    vcg::Use<PolyEdge>::AsEdgeType,
    vcg::Use<PolyFace>::AsFaceType> {
};

class PolyVertex : public vcg::Vertex<MyPolyTypes,
                                      vcg::vertex::Coord3d,
                                      vcg::vertex::Normal3f,
                                      vcg::vertex::Color4b,
                                      vcg::vertex::Qualityf,
                                      vcg::vertex::BitFlags,
                                      vcg::vertex::VFAdj,
                                      vcg::vertex::CurvatureDirf> {
};

class PolyFace : public vcg::Face<
    MyPolyTypes,
    vcg::face::PolyInfo,
    vcg::face::VertexRef,
    vcg::face::Normal3f,
    vcg::face::Color4b,
    vcg::face::Qualityf,
    vcg::face::BitFlags,
    vcg::face::PFVAdj,
    vcg::face::PFFAdj,
    vcg::face::PVFAdj,
    vcg::face::CurvatureDirf,
    vcg::face::Mark> {
};

class PolyEdge : public vcg::Edge<
    MyPolyTypes,
    vcg::edge::VertexRef,
    vcg::edge::BitFlags> {
};

class PolyMesh : public vcg::tri::TriMesh<
    std::vector<PolyVertex>,
                 std::vector<PolyEdge>,
                 std::vector<PolyFace>> {
};

}

#endif //TRICROSSFIELD__TCF_BASIC_TYPES_H_
