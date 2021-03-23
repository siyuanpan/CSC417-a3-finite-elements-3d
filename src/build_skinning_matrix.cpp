#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>
#include <igl/AABB.h>
#include <igl/in_element.h>

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T,
                           Eigen::Ref<const Eigen::MatrixXd> V_skin)
{

    igl::AABB<Eigen::MatrixXd, 3> tree;
    Eigen::MatrixXd V_copy = V;
    Eigen::MatrixXi T_copy = T;
    tree.init(V_copy, T_copy);
    Eigen::VectorXi I;
    igl::in_element(Eigen::MatrixXd(V), Eigen::MatrixXi(T), Eigen::MatrixXd(V_skin), tree, I);

    typedef Eigen::Triplet<double> Tri;
    std::vector<Tri> triplets;
    triplets.reserve(I.size() * 4);

    for (int i = 0; i < I.size(); ++i)
    {
        Eigen::Vector4d phi;
        phi_linear_tetrahedron(phi, V, T.row(I(i)), V_skin.row(i).transpose());
        for (int j = 0; j < 4; ++j)
            triplets.push_back(Tri(i, T(I(i), j), phi(j)));
    }
    N.resize(V_skin.rows(), V.rows());
    N.setFromTriplets(triplets.begin(), triplets.end());
}