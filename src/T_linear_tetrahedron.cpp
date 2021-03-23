#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume)
{

    double mass = density * volume;
    double c0 = (1.0 / 10.0) * mass;
    double c1 = (1.0 / 20.0) * mass;

    Eigen::Matrix1212d M;
    M.block(0, 0, 3, 3) = c0 * Eigen::Matrix3d::Identity();
    M.block(0, 3, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(0, 6, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(0, 9, 3, 3) = c1 * Eigen::Matrix3d::Identity();

    M.block(3, 0, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(3, 3, 3, 3) = c0 * Eigen::Matrix3d::Identity();
    M.block(3, 6, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(3, 9, 3, 3) = c1 * Eigen::Matrix3d::Identity();

    M.block(6, 0, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(6, 3, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(6, 6, 3, 3) = c0 * Eigen::Matrix3d::Identity();
    M.block(6, 9, 3, 3) = c1 * Eigen::Matrix3d::Identity();

    M.block(9, 0, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(9, 3, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(9, 6, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(9, 9, 3, 3) = c0 * Eigen::Matrix3d::Identity();

    Eigen::Vector12d v;
    v.segment<3>(0) = qdot.segment<3>(element(0) * 3);
    v.segment<3>(3) = qdot.segment<3>(element(1) * 3);
    v.segment<3>(6) = qdot.segment<3>(element(2) * 3);
    v.segment<3>(9) = qdot.segment<3>(element(3) * 3);

    T = 0.5 * v.transpose() * M * v;
}