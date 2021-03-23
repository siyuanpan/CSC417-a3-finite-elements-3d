#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0)
{

    typedef Eigen::Triplet<double> Tri;
    std::vector<Tri> triplets;
    triplets.reserve(T.rows() * 12 * 12);
    for (int i = 0; i < T.rows(); ++i)
    {
        Eigen::Matrix1212d Mi;
        Eigen::RowVectorXi element = T.row(i);
        mass_matrix_linear_tetrahedron(Mi, qdot, element, density, v0(i));

        for (int j = 0; j < 4; ++j)
        {
            for (int k = 0; k < 4; ++k)
            {
                for (int row = 0; row < 3; ++row)
                {
                    for (int col = 0; col < 3; ++col)
                    {
                        triplets.push_back(Tri(3 * element(j) + row, 3 * element(k) + col, Mi(3 * j + row, 3 * k + col)));
                    }
                }
            }
        }
    }

    M.resize(qdot.rows(), qdot.rows());
    M.setFromTriplets(triplets.begin(), triplets.end());
}