#include <Eigen/Dense>
#include <EigenTypes.h>
#include <iostream>
//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template <typename Ret, typename Integrand_Func>
inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q,
                                    Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                    Integrand_Func integrand)
{
    auto circumcenter = [](const Eigen::Matrix43d &t, Eigen::Vector3d &c) {
        Eigen::Matrix3d A;
        Eigen::Vector3d b;

        const double n0 = t.row(0).squaredNorm();

        for (int k = 0; k < 3; ++k)
        {
            A.row(k) = t.row(k + 1) - t.row(0);
            b(k) = t.row(k + 1).squaredNorm() - n0;
        }

        c = 0.5 * A.fullPivHouseholderQr().solve(b);
    };

    Eigen::Matrix43d t;
    for (int i = 0; i < 4; ++i)
    {
        t.row(i) = q.segment<3>(3 * element(i)).transpose();
    }
    Eigen::Vector3d cc;
    circumcenter(t, cc);

    integrand(integrated, q, element, cc);
    integrated = volume * integrated;
}
