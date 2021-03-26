#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the inputs:
    if(estimations.empty() || estimations.size() != ground_truth.size()) {
        std::cout << "Estimation vector is either empty or not the same size as the ground truth vector" << std::endl;
        return rmse;
    }

    VectorXd squared_residuals(estimations[0].size());
    for (int i=0; i < estimations.size(); ++i) {
        squared_residuals = squared_residuals.array() +
                            ((estimations[i] - ground_truth[i]).array() * (estimations[i] - ground_truth[i]).array());
    }

    // calculate the mean
    VectorXd mean = squared_residuals.array() / estimations.size();

    // calculate the squared root
    rmse = mean.array().sqrt();

    // return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3,4);
    // recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    double rho_magnitude = sqrt(pow(px, 2) + pow(py, 2));

    // check division by zero
    if(fabs(rho_magnitude) < 1e-3) {
        std::cout << "Error - Division by zero" << std::endl;
        return Hj;
    }

    // compute the Jacobian matrix
    Hj << px / rho_magnitude, py / rho_magnitude, 0, 0,
            -py / (pow(rho_magnitude,2)), px / (pow(rho_magnitude, 2)), 0, 0,
            (py * (vx*py - vy*px))/(pow(rho_magnitude,3)), (px*(vy*px - vx*py)) / (pow(rho_magnitude, 3)), px / rho_magnitude, py / rho_magnitude;

    return Hj;
}

Eigen::VectorXd Tools::convertPolarToCartesian(const Eigen::VectorXd& polarCoords) {
    // polar coords (rho, phi, rho-dot)
    // cart coords (x, y, vx, vy)
    Eigen::VectorXd cartCoords(4);
    cartCoords(0) = polarCoords(0) * cos(polarCoords(1));
    cartCoords(1) = polarCoords(0) * sin(polarCoords(1));
    cartCoords(2) = polarCoords(2) * cos(polarCoords(1));
    cartCoords(3) = polarCoords(2) * sin(polarCoords(1));

    return cartCoords;
}
