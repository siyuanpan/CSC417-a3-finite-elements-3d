#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X)
{
    auto X0 = V.row(element(0)).transpose();
    auto X1 = V.row(element(1)).transpose();
    auto X2 = V.row(element(2)).transpose();
    auto X3 = V.row(element(3)).transpose();

    Eigen::Matrix3d T;
    T.col(0) = X1 - X0;
    T.col(1) = X2 - X0;
    T.col(2) = X3 - X0;

    Eigen::Matrix3d Tinv = T.inverse();
    Eigen::Vector3d NegOne = -Eigen::Vector3d::Ones();
    dphi.row(0) = NegOne.transpose() * Tinv;
    dphi.block<3, 3>(1, 0) = Tinv;
}