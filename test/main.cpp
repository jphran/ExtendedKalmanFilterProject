//
// Created by jfrancis on 3/26/21.
//

#include "tools.h"

#include <iostream>

int main() {
    Eigen::VectorXd polarCoords(3);
    polarCoords << 5, 1, 10;

    auto boy = Tools::convertPolarToCartesian(polarCoords);
    std::cout << boy << std::endl;

    return 0;
}