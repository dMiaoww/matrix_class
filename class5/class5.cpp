#include <Eigen/Eigen>
#include <cmath>
#include <iostream>
#include <ostream>

// 实现高斯约旦消元法求解线性方程组 AX=b
bool gauss_jordan(Eigen::MatrixXd A, Eigen::VectorXd b, Eigen::VectorXd& x) {
  int n = b.size();
  for (int i = 0; i < n; i++) {
    // 搜索主元
    double maxEl = fabs(A(i, i));
    int maxRow = i;
    for (int k = i + 1; k < n; k++) {
      if (fabs(A(k, i)) > maxEl) {
        maxEl = A(k, i);
        maxRow = k; // 当前列中，最大值所在的行
      }
    }

    // 交换最大行
    A.row(maxRow).swap(A.row(i));
    std::swap(b(maxRow), b(i));

    // 把行缩放为主元为1
    for (int k = i + 1; k < n; k++) {   // 下面的行 - 倍数*当前行
      double coeff = A(k, i) / A(i, i); // 计算倍数
      for (int j = i; j < n; j++)
        A(k, j) -= coeff * A(i, j);
      b(k) -= coeff * b(i);
    }
  } // 得到高斯消元形，上三角矩阵

  for (int i = 0; i < n; i++) {
    if (A.row(i).sum() == 0 && b(i) != 0) {
      // 如果存在一行在A中所有元素都为0，而在b中对应的元素不为0，这个系统没有解
      return false;
    }
  }
  // 如果没有找到这样的行，那么这个系统至少有一个解


  // 简化到约旦形
  for (int i = n - 1; i >= 0; i--) { // 从下往上，简化成约旦形
    for (int j = i + 1; j < n; j++)
      b(i) -= A(i, j) * x(j);
    x(i) = b(i) / A(i, i);
  }
  return true;
}

int main() {
  Eigen::MatrixXd A(3, 3);
  A << 1, 2, 4, 3, 7, 2, 2, 3, 3;

  Eigen::VectorXd b(3);
  b << 7, -11, 1;

  Eigen::VectorXd x(3);
  gauss_jordan(A, b, x);
  std::cout << "Solution:\n" << x << std::endl;

  // 使用高斯约旦消元法求出行最简形式，判断有没有解
  Eigen::MatrixXd A2(4, 5);
  A2 << 1, -1, 2, 0, 3, 
  -1, 1, 0, 2, -5, 
  1, -1, 4, 2, 4, 
  -2, 2, -5, -1, -3;

  Eigen::VectorXd b2(4);
  b2 << 1, 5, 13, -1;
  Eigen::VectorXd x2(4);
  bool res = gauss_jordan(A2, b2, x2);
  std::cout << "\n" << res;

}