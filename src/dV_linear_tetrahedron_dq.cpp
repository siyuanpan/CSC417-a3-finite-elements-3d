#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q,
                              Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                              double C, double D)
{

    auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        Eigen::Matrix3d F;
        Eigen::Matrix34d t;
        for (int i = 0; i < 4; ++i)
        {
            t.col(i) = q.segment<3>(3 * element(i));
        }
        Eigen::Matrix43d dphi;
        dphi_linear_tetrahedron_dX(dphi, V, element, X);
        F = t * dphi;

        Eigen::Vector9d psi;
        dpsi_neo_hookean_dF(psi, F, C, D);

        Eigen::Matrix<double, 9, 12> B;
        for (int i = 0; i < 4; ++i)
        {
            B.block(0, 0 + 3 * i, 3, 1) = dphi.row(i).transpose();
            B.block(3, 1 + 3 * i, 3, 1) = dphi.row(i).transpose();
            B.block(6, 2 + 3 * i, 3, 1) = dphi.row(i).transpose();
        }
        dV = B.transpose() * psi;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);
}