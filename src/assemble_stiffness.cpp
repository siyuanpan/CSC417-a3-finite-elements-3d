#include <assemble_stiffness.h>
#include <d2V_linear_tetrahedron_dq2.h>
#include <d2V_spring_particle_particle_dq2.h>
#include <igl/edge_lengths.h>
#include <igl/edges.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot,
                        Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                        double C, double D)
{
    // Eigen::MatrixXd E;
    // Eigen::VectorXd l0;
    // igl::edges(T, E);
    // igl::edge_lengths(V, E, l0);
    // typedef Eigen::Triplet<double> Tri;
    // std::vector<Tri> triplets;
    // triplets.reserve(E.rows() * 9 * 4);
    // for (int i = 0; i < E.rows(); ++i)
    // {
    //     Eigen::Matrix66d H;
    //     auto v0 = E(i, 0);
    //     auto v1 = E(i, 1);
    //     auto q0 = V.row(E(i, 0)).transpose();
    //     auto q1 = V.row(E(i, 1)).transpose();
    //     d2V_spring_particle_particle_dq2(H, q0, q1, l0(i), 1e5);
    //     for (int row = 0; row < 3; ++row)
    //     {
    //         for (int col = 0; col < 3; ++col)
    //         {
    //             auto val = -H(row, col);
    //             triplets.push_back(Tri(3 * v0 + row, 3 * v0 + col, val));
    //             triplets.push_back(Tri(3 * v0 + row, 3 * v1 + col, -val));
    //             triplets.push_back(Tri(3 * v1 + row, 3 * v0 + col, -val));
    //             triplets.push_back(Tri(3 * v1 + row, 3 * v1 + col, val));
    //         }
    //     }
    // }

    // K.resize(q.rows(), q.rows());
    // K.setFromTriplets(triplets.begin(), triplets.end());

    typedef Eigen::Triplet<double> Tri;
    std::vector<Tri> triplets;
    triplets.reserve(T.rows() * 4 * 4 * 9);
    for (int i = 0; i < T.rows(); ++i)
    {
        Eigen::Matrix1212d dV;
        d2V_linear_tetrahedron_dq2(dV, q, V, T.row(i), v0(i), C, D);
        for (int j = 0; j < 4; ++j)
        {
            for (int k = 0; k < 4; ++k)
            {
                for (int row = 0; row < 3; ++row)
                {
                    for (int col = 0; col < 3; ++col)
                    {
                        triplets.push_back(Tri(3 * T(i, j) + row, 3 * T(i, k) + col, -dV(3 * j + row, 3 * k + col)));
                    }
                }
            }
        }
    }
    K.resize(q.rows(), q.rows());
    K.setFromTriplets(triplets.begin(), triplets.end());
};
