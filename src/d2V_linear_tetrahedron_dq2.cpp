#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>

void d2V_linear_tetrahedron_dq2(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q,
                                Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                double C, double D)
{

    auto neohookean_linear_tet = [&](Eigen::Matrix1212d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        //Code to compute non-integrated hessian matrix goes here
        Eigen::Matrix3d F;
        Eigen::Matrix34d t;
        for (int i = 0; i < 4; ++i)
        {
            t.col(i) = q.segment<3>(3 * element(i));
        }
        Eigen::Matrix43d dphi;
        dphi_linear_tetrahedron_dX(dphi, V, element, X);
        F = t * dphi;

        Eigen::Matrix99d psi;
        d2psi_neo_hookean_dF2(psi, F, C, D);

        Eigen::Matrix<double, 9, 12> B;
        for (int i = 0; i < 4; ++i)
        {
            B.block(0, 0 + 3 * i, 3, 1) = dphi.row(i).transpose();
            B.block(3, 1 + 3 * i, 3, 1) = dphi.row(i).transpose();
            B.block(6, 2 + 3 * i, 3, 1) = dphi.row(i).transpose();
        }
        dV = B.transpose() * psi * B;
    };

    //integrate the non-integrated hessian across the tetrahedral element
    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);

    //DO NOT REMOVE THIS CODE This code ensures that the hessian matrix is symmetric postive definite by projecting all
    //negative eigenvalues to small, postive values.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);

    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();

    for (int i = 0; i < 12; ++i)
    {
        if (es.eigenvalues()[i] < 1e-6)
        {
            DiagEval(i, i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();
}
