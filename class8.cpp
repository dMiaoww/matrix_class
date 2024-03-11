
#include <Eigen/Eigen>
#include <iostream>
#include <ostream>

int gauss_rank(Eigen::MatrixXd A) {

  int rowCount = A.rows();
  int colCount = A.cols();
  double tolerance = 1e-10;

  int rank = 0;
  for (int row = 0, col = 0; row < rowCount && col < colCount; ++col) {
    int el = row; // 主元

    for (int i = row + 1; i < rowCount; ++i) // 按行搜索最大值
      if (fabs(A(i, col)) > fabs(A(el, col)))
        el = i;

    if (fabs(A(el, col)) < tolerance) {  // 跳过这一列
      continue; 
    }

    A.row(row).swap(A.row(el));

    for (int i = row + 1; i < rowCount; ++i) {
      double coffe = A(i, col) / A(row, col);
      A.row(i) -= coffe * A.row(row);
    } 
    // 处理完一个元素（该元素的下方都为0），行和列都+1
    ++row;
  }

  for (int i = 0; i < rowCount; ++i) {
    if (A.row(i).norm() > tolerance) {
      ++rank;
    }
  }

  return rank;
}





int main() {
  Eigen::MatrixXd A(4,5);
  A << 1, -1, 2, 0,3,
    -1,1,0,2,-5,
    1,-1,4,2,4,
    -2,2,-5,-1,-3;
  std::cout << "原矩阵：\n" << A << std::endl;

  int rank = A.fullPivLu().rank();
  std::cout << "eigen求秩：" << rank << std::endl; // 3
  int rank2 = gauss_rank(A);
  std::cout << "高斯求秩：" << rank2 << std::endl;  // 4
  

  return 0;
}

