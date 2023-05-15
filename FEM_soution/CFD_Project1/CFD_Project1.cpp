#include <iostream>
#include <vector>
#include <boost/lexical_cast/detail/converter_lexical.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/quadrature/gauss.hpp> // 高斯-勒让德数值积分
#include <Eigen/Dense>               // 使用Eigen库的稠密矩阵函数
#include <Eigen/Sparse>              // 稀疏矩阵
#include <unordered_map>             // 哈希表类型

using namespace boost::math::differentiation;
using namespace boost::math::quadrature;
using namespace Eigen;

typedef std::vector<Triplet<float>> trip;  // 三元组的存储类型

#pragma region ClassExtension
// 自定义哈希函数(这一小片是纯来自chatGPT的)
struct PairHash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

// 定义单元类
class element {
public:
    element(std::vector<int> nodeIndex, Matrix2f K_e) {
        this->nodeIndex = nodeIndex;
        this->K_e = K_e;
    };
    Matrix2f K_e;             // 单元刚度矩阵
    std::vector<int> nodeIndex; // 节点相对编号,长度2向量，分别存储下标1,2的节点全局下标
};

#pragma endregion

#pragma region Stiff_Derivation
// Phi的全局函数，并使用类模板，定义成多类型变量求导函数
template<typename W, typename X, typename Y>
promote<W,X,Y> Phi(W x, X x_i, Y h) {
    if (x <= x_i - h || x >= x_i + h) return 0;
    if (x < x_i) return 1.0 - (x_i - x) / h;
    return 1.0 - (x - x_i)/ h;
}

// 计算一个单元的刚度矩阵，返回一个2x2的矩阵
Matrix2f Elem_Stiffness(double x_i, double x_j, double h) {
    Matrix2f elem_mat;
    float xVec[]{ x_i, x_j };
    float dPhi[2]; // 存储导数的值
    // 定义求导的自变量值
    constexpr unsigned Order = 1; // 最大求解一阶导数
    auto const x = make_fvar<double, 1> (x_i + x_j)/2;
    // 对应的求导变量的阶次以及设置对应变量求导值
    constexpr unsigned Nx = 1;
    constexpr unsigned Nx_i = 0;
    constexpr unsigned Nh = 0;
    // 由于是一次函数，在(x_i + x_j)/2求导数来代替整体的导数值
    auto const variables1 = make_ftuple<float, Nx, Nx_i, Nh>((x_i + x_j) / 2, x_i, h);
    auto const variables2 = make_ftuple<float, Nx, Nx_i, Nh>((x_i + x_j) / 2, x_j, h);
    auto const& _x = std::get<0>(variables1);
    auto const& _xi = std::get<1>(variables1);
    auto const& _xj = std::get<1>(variables2);
    auto const& _h = std::get<2>(variables1);
    auto const y1 = Phi(_x, _xi, _h);
    auto const y2 = Phi(_x, _xj, _h);
    dPhi[0] = y1.derivative(Nx, Nx_i, Nh); // i, j求导阶次相同
    dPhi[1] = y2.derivative(Nx, Nx_i, Nh); 
    for (int i = 0; i < 2; i++) {
        std::vector<float> rowVec;
        for (int j = 0; j < 2; j++) {
            // 定义积分函数
            // 最后面的是代入Phi_j进行求解
            auto a = [&](const float& x) {
                // ******* 注意: 由于换标，公式也需要换标  *****
                return 2.0 * dPhi[i] * dPhi[j] + (5.0 + (double)x * x) * Phi(x, xVec[i], h) * dPhi[j];
            };
            // 建立高斯积分器,10 是积分点个数(需要注意:gauss积分器仅有double类型部分)
            auto result = gauss<float, 10>::integrate(a, x_i, x_j);
            elem_mat(i, j) = (float)result;
        }
    }
    /* Test Code: std::cout << y1.derivative(Nx, Nx_i, Nh) << " " << y2.derivative(Nx, Nx_i, Nh) << std::endl;
    result : -10, 10;*/
    return elem_mat;
}

// 组装整体刚度矩阵
SparseMatrix<float> Tot_Stiffness(int elem_num, int node_num ,std::vector<element*> elems){
    SparseMatrix<float> A(node_num, node_num);
    A.setZero();  // 将整个稀疏矩阵初始化为0矩阵
    trip triplets; // 定义三元组存储数据
    triplets.reserve(3 * elem_num + 1); // 分配非零元素的空间
    // 由于std::pair不能作为哈希表的键值使用, 因此自定义一个哈希函数
    std::unordered_map<std::pair<int, int>, float, PairHash> visited;
    std::vector<std::pair<int, int>> indexPair;
    // 为了运算速度, 每一次存储i,j是否访问过, 这个只需要使用哈希表进行存储即可，避免了查找时间消耗
    for (element* elem : elems) {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                // 由于nodeIndex是整体坐标
                int m = elem->nodeIndex[i];
                int n = elem->nodeIndex[j];
                if (visited.find(std::make_pair(m, n)) != visited.end()) {
                    visited[std::make_pair(m, n)] += elem->K_e(i, j);
                }
                else { // 没有找到
                    visited[std::make_pair(m, n)] = elem->K_e(i, j);
                    indexPair.push_back(std::make_pair(m, n));
                }
            }
        }
    }
    // 之后建立三元组放入元素,并使用三元组建立稀疏矩阵，这是一种更快的矩阵建立方式
    for (std::pair<int, int> P : indexPair) {
        triplets.push_back(Triplet<float>(P.first, P.second, visited[P]));
        // 因为三元组只能psuh_back,不能修改值,有一点时间浪费
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed(); // 压缩存储
    return A;
}

// 计算整个右端项f, 这需要传入节点对应的坐标， 以及单元向量
VectorXf fVecGenerate(int elem_num, int node_num,std::vector<element*> elems, std::vector<float> location) {
    // 对于每一个单元, 计算
    VectorXf FVec(node_num);
    FVec.setZero();
    for (int index = 0; index < elem_num; index++) {
        int i_index = elems[index]->nodeIndex[0];
        int j_index = elems[index]->nodeIndex[1];
        float x_i = location[i_index];  
        float x_j = location[j_index];
        float h = x_j - x_i;  // 使用x_j - x_i 计算出dx

        // 一个单元下两个i的相对坐标
        for (int i = 0; i < 2; i++) {
            int index_relative = elems[index]->nodeIndex[i];
            float loc_relative = location[index_relative];
            auto f = [&](const float& x) {
                return 2.0 * x * Phi(x, loc_relative , h);
            };
            float append = 4 * Phi(1, loc_relative, h);  // 边界条件处的值
            auto result = gauss<float, 10>::integrate(f, x_i, x_j) + append;
            FVec(index_relative) += result;
        }
    }
    return FVec;
}

#pragma endregion

// 一维问题求解入口函数
void Solution1d(float length,int elem_num) {
    int node_num = elem_num + 1;
    float h = length / elem_num; // the number of the element

    std::vector<element*> elems;
    std::vector<float> cord;
    for (int i = 0; i < node_num; i++) {
        cord.push_back(i * h);
    }

    // 存储相对坐标,绝对坐标和每个单元的刚度矩阵
    for (int i = 0; i < elem_num; i++) {
        element* e = new element(std::vector<int> { i, i + 1 }, Elem_Stiffness(cord[i], cord[i + 1], h));
        elems.push_back(e);
    }

    // 组装整体刚度矩阵
    SparseMatrix<float> tot_Stiff = Tot_Stiffness(elem_num, node_num, elems);
    VectorXf FVec = fVecGenerate(elem_num, node_num, elems, cord);
    
    /**********************求解方程右侧的f向量*******************************/

    // 对于本质边界条件的部分，使用对角线元素扩大方法进行求解的修正
    for (SparseMatrix<float>::InnerIterator iter(tot_Stiff, 0); iter ; ++iter) {
        std::cout << iter.row() << " " << iter.col() << " " << iter.value() << std::endl;
        if (iter.row() == 0 && iter.col() == 0) {
            iter.valueRef() *= 1e10;        
            FVec[0] = 5 * iter.value();     // 使用主对角线的元素扩大方法进行求解
            // 将右侧的数以 aK_{jj} \bar{a_j} 来代替，其中\bar{a_j}是给出的边界值
        }
    }

    std::cout << "---------- The Stiff Matrix ----------- " << std::endl <<  tot_Stiff << std::endl;
    std::cout << "---------- The FVector ----------- " << std::endl  << FVec << std::endl;

    // 创建稀疏LU分解器
    SparseLU<SparseMatrix<float>> solver;
    solver.analyzePattern(tot_Stiff);  // 分析矩阵形式
    solver.factorize(tot_Stiff);       // 进行LU分解

    // 求解线性方程组
    VectorXf x = solver.solve(FVec);
    std::cout << std::endl << x << std::endl;
}

int main()
{
    int elem_num = 10;
    float length = 1.0;
    Solution1d(length, elem_num);
    return 0;
}