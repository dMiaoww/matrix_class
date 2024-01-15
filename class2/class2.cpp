#include <Eigen/Eigen>
#include <cmath>
#include <iostream>
#include <ostream>


int main() {
  // 计算向量的模  sqrt(各元素平方之和)
  Eigen::Vector3d vectorC{0,0,3};
  std::cout << "norm of vectorC: " << vectorC.norm() << "\n";   

  // 向量归一化  转变成模为1
  Eigen::Vector3d  d = vectorC.normalized();
  std::cout << d  << std::endl; 

  // 向量点乘
  double x = vectorC.dot(d);
  std::cout << x << std::endl;

  // 向量夹角
  Eigen::Vector3d vectorE{3,0,0};
  double angle = acos(vectorC.dot(vectorE)/(vectorC.norm()*vectorE.norm()));
  std::cout << angle << std::endl;

  // 一个向量到另一个向量投影点坐标  c在e上的投影
  Eigen::Vector3d projection = (vectorC.dot(vectorE) / vectorE.squaredNorm()) * vectorE;
  std::cout << projection << "\n";
}