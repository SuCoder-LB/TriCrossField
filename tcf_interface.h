//
// Created by 苏立彪 on 2023/5/11.
//

#ifndef TRICROSSFIELD__TCF_INTERFACE_H_
#define TRICROSSFIELD__TCF_INTERFACE_H_

#include "tcf_basic_types.h"
#include "tcf_linear_system.h"

namespace tri_cross_field {

void ComputeCrossField(PolyMesh &tri_mesh,std::vector<std::array<int,2>>&features,
                       LinearSystem*linear_system,std::vector<std::array<double,3>>&field,
                       const CrossFieldParameter &cfp);

}

#endif //TRICROSSFIELD__TCF_INTERFACE_H_
