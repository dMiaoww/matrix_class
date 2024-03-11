
#include <Eigen/Eigen>
#include <iostream>
#include <ostream>

// 标准正交基
Eigen::MatrixXd gram_schmidt(const Eigen::MatrixXd& V) {
    Eigen::MatrixXd U(V);
    for (int i = 0; i < V.cols(); i++) {
        for (int j = 0; j < i; j++) {
            U.col(i) -= U.col(j).dot(V.col(i)) / U.col(j).squaredNorm() * U.col(j);
        }
    }
    U.colwise().normalize();
    return U;
}


// 实现QR分解，Q是正交矩阵，R是上三角矩阵，正交矩阵就是 Q.T = Q.inv
// 用来求解最小二乘问题



int main() {
  // 测试标准正交基
  Eigen::MatrixXd b1(2,2);
  b1 << 2,1,
        1,1;

  std::cout << "basis1: \n" << b1 << "\nafter: \n" << gram_schmidt(b1);

  Eigen::MatrixXd b2(2,2);
  b2 << 2,4,
        3,5;
  std::cout << "\nbasis2: \n" << b2 << "\nafter: \n" << gram_schmidt(b2) << "\n";

    Eigen::MatrixXd b3(3,3);
  b3 << 1,3,-1,
        0,1,-1,
        1,1,-1;
  std::cout << "\nbasis3: \n" << b3 << "\nafter: \n" << gram_schmidt(b3) << "\n";

      Eigen::MatrixXd b4(4,3);
  b4 << 1,-3,-1,
        1, 3,-2,
        5, 4,2,
        2,-2, 5;
  std::cout << "\nbasis4: \n" << b4 << "\nafter: \n" << gram_schmidt(b4) << "\n";



  // 测试QR分解
  Eigen::MatrixXd A1(3,3);
  A1 << 1,1,2,
        1,1,0,
        1,0,0;
  
  Eigen::MatrixXd A2(3,3);
  A2 << 2,-1,-1,
        2,0,2,
        2,-1,3;

  std::cout << "原矩阵：\n" << A1 << std::endl;

    // HouseholderQR分解
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A1);

    // 获取Q和R
    Eigen::MatrixXd Q = qr.householderQ();
    Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();

    std::cout << "Q矩阵：\n" << Q << std::endl;
    std::cout << "R矩阵：\n" << R << std::endl;

    // 验证QR分解，应该得到原矩阵
    std::cout << "QR复原：\n" << Q * R << std::endl;  


     std::cout << "原矩阵：\n" << A1 << std::endl;

    // HouseholderQR分解
    Eigen::HouseholderQR<Eigen::MatrixXd> qr2(A2);

    // 获取Q和R
    Eigen::MatrixXd Q2 = qr2.householderQ();
    Eigen::MatrixXd R2 = qr2.matrixQR().triangularView<Eigen::Upper>();

    std::cout << "Q矩阵：\n" << Q2 << std::endl;
    std::cout << "R矩阵：\n" << R2 << std::endl;

    // 验证QR分解，应该得到原矩阵
    std::cout << "QR复原：\n" << Q2 * R2 << std::endl;  


  return 0;
}

