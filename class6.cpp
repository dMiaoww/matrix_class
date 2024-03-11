// 初等矩阵和矩阵的可逆性
#include <Eigen/Eigen>
#include <iostream>
#include <ostream>

// 实现高斯消元法求矩阵的逆
Eigen::MatrixXd gaussian_inverse(Eigen::MatrixXd matrix) {
  int n = matrix.rows();
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd combined(n, 2 * n);
  combined << matrix, I; // 构造增广矩阵[A,I],把A通过初等变换变成I，I就成了A的逆矩阵
  
  for (int i = 0; i < n; ++i) {
    float maxEl = abs(combined(i, i));
    int maxRow = i;
    for (int k = i + 1; k < n; ++k) {
      if (abs(combined(k, i)) > maxEl) {
        maxEl = combined(k, i);
        maxRow = k;
      }
    }

   if (maxRow != i) {
     combined.row(maxRow).swap(combined.row(i));
   }

    for (int k = i + 1; k < n; ++k) {
      float constant = -combined(k, i) / combined(i, i);
      combined.row(k) += constant * combined.row(i);
    }
  }// 得到上三角矩阵

  for (int i = n - 1; i >= 0; --i) {
    for (int k = i - 1; k >= 0; --k) {
      float constant = -combined(k, i) / combined(i, i);
      combined.row(k) += constant * combined.row(i);
    }
  }// 从下往上，变成对角线矩阵
  for (int i = 0; i < n; ++i) {
    combined.row(i) /= combined(i, i);
  }// 变成单位矩阵

  Eigen::MatrixXd inverse = combined.block(0, n, n, n);
  return inverse;
}


// 使用lu分解方阵， L：下三角 U：上三角 A可以分解成L和U的乘积
void LUDecomposition(const Eigen::MatrixXd& matrix) {
    int n = matrix.rows();

    // Initialize L matrix as identity matrix
    Eigen::MatrixXd L = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; i++) {
        // U matrix formation
        for (int j = i; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++) {
                sum += (L(i, k) * U(k, j));
            }
            U(i, j) = matrix(i, j) - sum;
        }

        // L matrix formation
        for (int j = i; j < n; j++) {
            if (i == j)
                L(i, i) = 1;  // Diagonal values of L are 1
            else {
                double sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += (L(j, k) * U(k, i));
                }
                L(j, i) = (matrix(j, i) - sum) / U(i, i);
            }
        }
    }

    //print L and U
    std::cout << "L matrix: " << std::endl;
    std::cout << L << std::endl;
    std::cout << "U matrix: " << std::endl;
    std::cout << U << std::endl;
    std::cout << "LU:\n" << L*U << "\n";
}


int main() {
  Eigen::MatrixXd A(2,2);
  A << 1, 2, 3, 4;
  std::cout << "原矩阵：\n" << A << std::endl;
  Eigen::MatrixXd inv1 = gaussian_inverse(A);
  std::cout << "高斯法逆矩阵:\n" << inv1 << std::endl;
  Eigen::MatrixXd inv2 = A.inverse();
  std::cout << "eigen库逆矩阵:\n" << inv2 << std::endl;

  Eigen::MatrixXd AA(3,3);
  AA << 1, 2, 3, 4,5,6, 3,-3,5;
  LUDecomposition(AA); // A=LU


  // 使用eigen库plu分解 ，P是置换矩阵，相当于交换行的顺序, 结果是PA=LU
  Eigen::PartialPivLU<Eigen::MatrixXd> plu(AA);
  Eigen::MatrixXd L = Eigen::MatrixXd::Identity(AA.rows(), AA.cols());
  L += plu.matrixLU().triangularView<Eigen::StrictlyLower>();
  Eigen::MatrixXd U = plu.matrixLU().triangularView<Eigen::Upper>();  // 这样得到的包括对角线元素
  Eigen::MatrixXd P = plu.permutationP();
  std::cout << "P:" << std::endl <<  P << std::endl;
  std::cout << "L:" << L << std::endl;
  std::cout << "U:" << std::endl << U << std::endl;
  std::cout << "after plu\n" << P.transpose()*L*U << std::endl; 


  // 使用eigen库plup‘分解
  Eigen::FullPivLU<Eigen::MatrixXd> lu(AA);
  Eigen::MatrixXd L2 = Eigen::MatrixXd::Identity(AA.rows(), AA.cols());
  L2 += lu.matrixLU().triangularView<Eigen::StrictlyLower>();
  Eigen::MatrixXd U2 = lu.matrixLU().triangularView<Eigen::Upper>();
  Eigen::MatrixXd P2 = lu.permutationP();
  Eigen::MatrixXd Q = lu.permutationQ();
  std::cout << "L:" << std::endl << L2 << std::endl;
  std::cout << "U:" << std::endl << U2 << std::endl;
  std::cout << "P:" << std::endl << P2 << std::endl;
  std::cout << "Q:" << std::endl << Q << std::endl;
  std::cout << "after plup'\n" << P2.transpose()*L2*U2*Q.transpose() << std::endl; 

  
  return 0;
}

