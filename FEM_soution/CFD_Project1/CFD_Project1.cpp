#include <iostream>
#include <vector>
#include <boost/lexical_cast/detail/converter_lexical.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/quadrature/gauss.hpp> // ��˹-���õ���ֵ����
#include <Eigen/Dense>               // ʹ��Eigen��ĳ��ܾ�����
#include <Eigen/Sparse>              // ϡ�����
#include <unordered_map>             // ��ϣ������

using namespace boost::math::differentiation;
using namespace boost::math::quadrature;
using namespace Eigen;

typedef std::vector<Triplet<float>> trip;  // ��Ԫ��Ĵ洢����

#pragma region ClassExtension
// �Զ����ϣ����(��һСƬ�Ǵ�����chatGPT��)
struct PairHash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

// ���嵥Ԫ��
class element {
public:
    element(std::vector<int> nodeIndex, Matrix2f K_e) {
        this->nodeIndex = nodeIndex;
        this->K_e = K_e;
    };
    Matrix2f K_e;             // ��Ԫ�նȾ���
    std::vector<int> nodeIndex; // �ڵ���Ա��,����2�������ֱ�洢�±�1,2�Ľڵ�ȫ���±�
};

#pragma endregion

#pragma region Stiff_Derivation
// Phi��ȫ�ֺ�������ʹ����ģ�壬����ɶ����ͱ����󵼺���
template<typename W, typename X, typename Y>
promote<W,X,Y> Phi(W x, X x_i, Y h) {
    if (x <= x_i - h || x >= x_i + h) return 0;
    if (x < x_i) return 1.0 - (x_i - x) / h;
    return 1.0 - (x - x_i)/ h;
}

// ����һ����Ԫ�ĸնȾ��󣬷���һ��2x2�ľ���
Matrix2f Elem_Stiffness(double x_i, double x_j, double h) {
    Matrix2f elem_mat;
    float xVec[]{ x_i, x_j };
    float dPhi[2]; // �洢������ֵ
    // �����󵼵��Ա���ֵ
    constexpr unsigned Order = 1; // ������һ�׵���
    auto const x = make_fvar<double, 1> (x_i + x_j)/2;
    // ��Ӧ���󵼱����Ľ״��Լ����ö�Ӧ������ֵ
    constexpr unsigned Nx = 1;
    constexpr unsigned Nx_i = 0;
    constexpr unsigned Nh = 0;
    // ������һ�κ�������(x_i + x_j)/2��������������ĵ���ֵ
    auto const variables1 = make_ftuple<float, Nx, Nx_i, Nh>((x_i + x_j) / 2, x_i, h);
    auto const variables2 = make_ftuple<float, Nx, Nx_i, Nh>((x_i + x_j) / 2, x_j, h);
    auto const& _x = std::get<0>(variables1);
    auto const& _xi = std::get<1>(variables1);
    auto const& _xj = std::get<1>(variables2);
    auto const& _h = std::get<2>(variables1);
    auto const y1 = Phi(_x, _xi, _h);
    auto const y2 = Phi(_x, _xj, _h);
    dPhi[0] = y1.derivative(Nx, Nx_i, Nh); // i, j�󵼽״���ͬ
    dPhi[1] = y2.derivative(Nx, Nx_i, Nh); 
    for (int i = 0; i < 2; i++) {
        std::vector<float> rowVec;
        for (int j = 0; j < 2; j++) {
            // ������ֺ���
            // �������Ǵ���Phi_j�������
            auto a = [&](const float& x) {
                // ******* ע��: ���ڻ��꣬��ʽҲ��Ҫ����  *****
                return 2.0 * dPhi[i] * dPhi[j] + (5.0 + (double)x * x) * Phi(x, xVec[i], h) * dPhi[j];
            };
            // ������˹������,10 �ǻ��ֵ����(��Ҫע��:gauss����������double���Ͳ���)
            auto result = gauss<float, 10>::integrate(a, x_i, x_j);
            elem_mat(i, j) = (float)result;
        }
    }
    /* Test Code: std::cout << y1.derivative(Nx, Nx_i, Nh) << " " << y2.derivative(Nx, Nx_i, Nh) << std::endl;
    result : -10, 10;*/
    return elem_mat;
}

// ��װ����նȾ���
SparseMatrix<float> Tot_Stiffness(int elem_num, int node_num ,std::vector<element*> elems){
    SparseMatrix<float> A(node_num, node_num);
    A.setZero();  // ������ϡ������ʼ��Ϊ0����
    trip triplets; // ������Ԫ��洢����
    triplets.reserve(3 * elem_num + 1); // �������Ԫ�صĿռ�
    // ����std::pair������Ϊ��ϣ��ļ�ֵʹ��, ����Զ���һ����ϣ����
    std::unordered_map<std::pair<int, int>, float, PairHash> visited;
    std::vector<std::pair<int, int>> indexPair;
    // Ϊ�������ٶ�, ÿһ�δ洢i,j�Ƿ���ʹ�, ���ֻ��Ҫʹ�ù�ϣ����д洢���ɣ������˲���ʱ������
    for (element* elem : elems) {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                // ����nodeIndex����������
                int m = elem->nodeIndex[i];
                int n = elem->nodeIndex[j];
                if (visited.find(std::make_pair(m, n)) != visited.end()) {
                    visited[std::make_pair(m, n)] += elem->K_e(i, j);
                }
                else { // û���ҵ�
                    visited[std::make_pair(m, n)] = elem->K_e(i, j);
                    indexPair.push_back(std::make_pair(m, n));
                }
            }
        }
    }
    // ֮������Ԫ�����Ԫ��,��ʹ����Ԫ�齨��ϡ���������һ�ָ���ľ�������ʽ
    for (std::pair<int, int> P : indexPair) {
        triplets.push_back(Triplet<float>(P.first, P.second, visited[P]));
        // ��Ϊ��Ԫ��ֻ��psuh_back,�����޸�ֵ,��һ��ʱ���˷�
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed(); // ѹ���洢
    return A;
}

// ���������Ҷ���f, ����Ҫ����ڵ��Ӧ�����꣬ �Լ���Ԫ����
VectorXf fVecGenerate(int elem_num, int node_num,std::vector<element*> elems, std::vector<float> location) {
    // ����ÿһ����Ԫ, ����
    VectorXf FVec(node_num);
    FVec.setZero();
    for (int index = 0; index < elem_num; index++) {
        int i_index = elems[index]->nodeIndex[0];
        int j_index = elems[index]->nodeIndex[1];
        float x_i = location[i_index];  
        float x_j = location[j_index];
        float h = x_j - x_i;  // ʹ��x_j - x_i �����dx

        // һ����Ԫ������i���������
        for (int i = 0; i < 2; i++) {
            int index_relative = elems[index]->nodeIndex[i];
            float loc_relative = location[index_relative];
            auto f = [&](const float& x) {
                return 2.0 * x * Phi(x, loc_relative , h);
            };
            float append = 4 * Phi(1, loc_relative, h);  // �߽���������ֵ
            auto result = gauss<float, 10>::integrate(f, x_i, x_j) + append;
            FVec(index_relative) += result;
        }
    }
    return FVec;
}

#pragma endregion

// һά���������ں���
void Solution1d(float length,int elem_num) {
    int node_num = elem_num + 1;
    float h = length / elem_num; // the number of the element

    std::vector<element*> elems;
    std::vector<float> cord;
    for (int i = 0; i < node_num; i++) {
        cord.push_back(i * h);
    }

    // �洢�������,���������ÿ����Ԫ�ĸնȾ���
    for (int i = 0; i < elem_num; i++) {
        element* e = new element(std::vector<int> { i, i + 1 }, Elem_Stiffness(cord[i], cord[i + 1], h));
        elems.push_back(e);
    }

    // ��װ����նȾ���
    SparseMatrix<float> tot_Stiff = Tot_Stiffness(elem_num, node_num, elems);
    VectorXf FVec = fVecGenerate(elem_num, node_num, elems, cord);
    
    /**********************��ⷽ���Ҳ��f����*******************************/

    // ���ڱ��ʱ߽������Ĳ��֣�ʹ�öԽ���Ԫ�����󷽷�������������
    for (SparseMatrix<float>::InnerIterator iter(tot_Stiff, 0); iter ; ++iter) {
        std::cout << iter.row() << " " << iter.col() << " " << iter.value() << std::endl;
        if (iter.row() == 0 && iter.col() == 0) {
            iter.valueRef() *= 1e10;        
            FVec[0] = 5 * iter.value();     // ʹ�����Խ��ߵ�Ԫ�����󷽷��������
            // ���Ҳ������ aK_{jj} \bar{a_j} �����棬����\bar{a_j}�Ǹ����ı߽�ֵ
        }
    }

    std::cout << "---------- The Stiff Matrix ----------- " << std::endl <<  tot_Stiff << std::endl;
    std::cout << "---------- The FVector ----------- " << std::endl  << FVec << std::endl;

    // ����ϡ��LU�ֽ���
    SparseLU<SparseMatrix<float>> solver;
    solver.analyzePattern(tot_Stiff);  // ����������ʽ
    solver.factorize(tot_Stiff);       // ����LU�ֽ�

    // ������Է�����
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