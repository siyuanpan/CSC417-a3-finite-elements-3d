#include <phi_linear_tetrahedron.h>
#include <iostream>

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x)
{
    auto X0 = V.row(element(0)).transpose();
    auto X1 = V.row(element(1)).transpose();
    auto X2 = V.row(element(2)).transpose();
    auto X3 = V.row(element(3)).transpose();

    Eigen::Matrix3d T;
    T.col(0) = X1 - X0;
    T.col(1) = X2 - X0;
    T.col(2) = X3 - X0;

    phi.segment<3>(1) = T.fullPivHouseholderQr().solve(x - X0);
    phi(0) = 1 - phi(1) - phi(2) - phi(3);
}