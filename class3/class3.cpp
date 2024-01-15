#include <Eigen/Eigen>
#include <iostream>
#include <ostream>


int main() {
  // 定义矩阵
  Eigen::MatrixXd b(5, 5);
  b << 1,2,3,4,5,
        1,2,3,4,5,
        1,2,3,4,5,
        1,2,3,4,5,
        1,2,3,4,5;

  // 获取矩阵的行数和列数
  std::cout << "row: " << b.rows() << "\n";
  std::cout << "col: " << b.cols() << "\n";

  // 获取矩阵元素的数量
  std::cout << "size: " << b.size() << "\n";

  // 取出某一行的行向量
  Eigen::RowVectorXd row_vector = b.row(2);
  std::cout << "row 2: " << row_vector << "\n";

  // 取出列向量
  Eigen::VectorXd col_vector = b.col(3);
  std::cout << "col 3: " << col_vector << "\n";

  // 矩阵的幂
  Eigen::MatrixXd c = b*b;
  std::cout << "c: " << c << "\n";

  // 矩阵的转置
  Eigen::MatrixXd ct = c.transpose();
  std::cout << "ct: " << ct << "\n";
}