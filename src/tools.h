#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

  /**
  * A helper method to convert a Cartesian vector to a polar vector.
  * 
  */
  Eigen::VectorXd CartesianToPolar(const Eigen::VectorXd& x);

  /**
  * A helper method to convert a polar vector to a Cartesian vector.
  *
  */
  Eigen::VectorXd PolarToCartesian(const Eigen::VectorXd& x);

  /**
  * A helper method to normalize the angle in radians to -Pi and Pi
  */
  double normalizeAngle(double angle_rad);
};

#endif /* TOOLS_H_ */
