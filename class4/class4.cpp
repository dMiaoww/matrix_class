#include <Eigen/Eigen>
#include <cmath>
#include <iostream>
#include <ostream>

int main() {
  // 定义矩阵
  Eigen::MatrixXd a(2, 10);
  a << 0, 0, 3, 3, 1, 1, 2, 2, 1, 1, 0, 5, 5, 4, 4, 3, 3, 2, 2, 0;

  // 缩放 [[2,0],[0,1.5]]
  Eigen::Matrix2d scale;
  scale << 2, 0, 0, 1.5;

  Eigen::MatrixXd afer = scale * a;
  std::cout << "after: \n";
  std::cout << afer;

  // 绕z轴顺时针旋转
  // 旋转矩阵，[[cos, sin][-sin, cos]]
  double theta = -90.0 / 180.0 * M_PI;
  Eigen::Matrix2d rotate;
  rotate << cos(theta), sin(theta), -sin(theta), cos(theta);
  Eigen::MatrixXd aferR = rotate * a;
  std::cout << "\nafterR: \n";
  std::cout << aferR;

  // 生成4*4单位矩阵
  Eigen::Matrix4d I = Eigen::MatrixXd::Identity(4, 4);;
  std::cout << "\nI: \n";
  std::cout << I;

  // 求逆
  Eigen::Matrix2d b;
  b << 1,2,3,4;
  Eigen::Matrix2d inv = b.inverse();
  std::cout << "\ninv: \n";
  std::cout << inv;

}