#include <iostream>
#include <Eigen/Dense>
#include "util/eigen_helper.hpp"

int main() {

    Eigen::VectorXd v(5);
    v << 10, 20, 30, 40, 50;

    Eigen::VectorXi idx(3);
    idx << 0, 2, 4;

    Eigen::VectorXd out(3);

    eigen_fast::take(v, idx, out);

    std::cout << out << std::endl;

    return 0;
}
