#include <assemble_forces.h>
#include <dV_linear_tetrahedron_dq.h>
#include <dV_spring_particle_particle_dq.h>
#include <igl/edges.h>
#include <igl/edge_lengths.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D)
{
    f.resize(q.rows());
    f.setZero();

    Eigen::MatrixXd E;
    Eigen::VectorXd l0;
    igl::edges(T, E);
    igl::edge_lengths(V, E, l0);

    for (int i = 0; i < E.rows(); ++i)
    {
        auto q0 = q.segment<3>(3 * E(i, 0));
        auto q1 = q.segment<3>(3 * E(i, 1));

        Eigen::Vector6d fi;
        dV_spring_particle_particle_dq(fi, q0, q1, l0(i), 1e5);
        f.segment<3>(3 * E(i, 0)) -= fi.segment<3>(0);
        f.segment<3>(3 * E(i, 1)) -= fi.segment<3>(3);
    }

#define force(i, j) f.segment<3>(3 * T((i), (j)))
    for (int i = 0; i < T.rows(); ++i)
    {
        Eigen::Vector12d dV;
        dV_linear_tetrahedron_dq(dV, q, V, T.row(i), v0(i), C, D);
        force(i, 0) -= dV.segment<3>(0);
        force(i, 1) -= dV.segment<3>(3);
        force(i, 2) -= dV.segment<3>(6);
        force(i, 3) -= dV.segment<3>(9);
    }
#undef force
};