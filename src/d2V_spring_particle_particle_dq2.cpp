#include <d2V_spring_particle_particle_dq2.h>
#include <iostream>

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness)
{
    double l = (q0 - q1).norm();
    Eigen::Matrix3d K = stiffness * (Eigen::Matrix3d::Identity() - l0 / l * (Eigen::Matrix3d::Identity() - (q1 - q0) * (q1 - q0).transpose() / (l * l)));

    H.block<3, 3>(0, 0) = K;
    H.block<3, 3>(0, 3) = -K;
    H.block<3, 3>(3, 0) = -K;
    H.block<3, 3>(3, 3) = K;
}