#include <Eigen/Eigen>
#include <iostream>
#include <ostream>


int main() {
  // 定义零矩阵
  Eigen::MatrixXd matrixA(1, 1);
  Eigen::MatrixXd matrixB(5, 5);
  std::cout << matrixA << std::endl;
  

  // 定义零向量
  Eigen::VectorXd vectorC = Eigen::VectorXd::Zero(5);
  Eigen::VectorXd vectorD = Eigen::VectorXd::Zero(5);

  vectorC << 1, 1, 1, 1, 1;

  // 访问矩阵中的元素
  matrixA(0,0) = 1;
  matrixB(4,4) = 1;
  std::cout << matrixA(0,0) << " " << matrixB(4,4) << std::endl; 

  // 加减、数乘
  auto VE = 5*(vectorC + vectorD - vectorD);

  std::cout << VE << std::endl;

  // 向量是特殊的Matrix
  // Matrix<double, row, col> , 如果row和col是固定的，就是固定大小的矩阵，比如 Matrix4f 本质上就是固定长度16的数组，存储在栈内
  // typedef Matrix<double, Dynamic, Dynamic> MatrixXd  这种是动态的矩阵，运行的时候会动态调整容量 Dynamic是-1，存储在堆上
  // 动态数组类似 std::vector 管理内存
  // 所以虽然一开始定义了a是1*1，B是5*5，但是也可以将B赋值给A
  matrixA = matrixB;
  std::cout << matrixA << std::endl;
  // 计算通过重载运算符实现计算
}