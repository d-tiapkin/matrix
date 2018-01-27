#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <utility>

namespace rationalSpace {

    class Rational {
    public:
        Rational() :
                _numer(0),
                _denom(1) {}

        Rational(int numer) :
                _numer(numer),
                _denom(1) {}

        Rational(int numer, int denom) :
                _numer(numer),
                _denom(denom) { normalize(); }

        Rational &operator=(const Rational &oth) {
            _numer = oth._numer;
            _denom = oth._denom;
            normalize();

            return *this;
        }

        // Compare operations
        bool operator==(const Rational &oth) const {
            return _numer == oth._numer && _denom == oth._denom;
        }

        bool operator!=(const Rational &oth) const {
            return !(*this == oth);
        }

        bool operator<(const Rational &oth) const {
            return _numer * oth._denom < oth._numer * _denom;
        }

        bool operator>(const Rational &oth) const {
            return _numer * oth._denom > oth._numer * _denom;
        }
        // Unary operations
        Rational operator+() const {
            return Rational(*this);
        }

        Rational operator-() const {
            Rational copy(*this);
            copy._numer *= -1;
            copy.normalize();

            return copy;
        }

        // ++ and --
        Rational &operator++() {
            _numer += _denom;

            normalize();
            return *this;
        }

        Rational operator++(int) {
            Rational copy(*this);
            _numer += _denom;

            normalize();
            return copy;
        }

        Rational &operator--() {
            _numer -= _denom;

            normalize();
            return *this;
        }

        Rational operator--(int) {
            Rational copy(*this);
            _numer -= _denom;

            normalize();
            return copy;
        }

        // +=, -=, *=, /=
        Rational &operator+=(const Rational &oth) {
            _numer = _numer * oth._denom + _denom * oth._numer;
            _denom = _denom * oth._denom;

            normalize();
            return *this;
        }

        Rational &operator-=(const Rational &oth) {
            _numer = _numer * oth._denom - _denom * oth._numer;
            _denom = _denom * oth._denom;

            normalize();
            return *this;
        }

        Rational &operator*=(const Rational &oth) {
            _numer *= oth._numer;
            _denom *= oth._denom;

            normalize();
            return *this;
        }

        Rational &operator/=(const Rational &oth) {
            *this = Rational(_numer * oth._denom, _denom * oth._numer);

            normalize();
            return *this;
        }


        // Getters
        const int& numerator() const {
            return _numer;
        }

        const int& denominator() const {
            return _denom;
        }


    private:
        static int gcd(int a, int b) {
            int sign = (b > 0) ? 1 : -1;

            b = static_cast<int>(std::abs(b));
            a = static_cast<int>(std::abs(a));

            while (b > 0) {
                a %= b;
                std::swap(a, b);
            }
            return sign * a;
        }

        void normalize() {
            int cur = gcd(_numer, _denom);
            _numer /= cur;
            _denom /= cur;
        }

        int _numer, _denom;
    };

    std::ostream& operator<< (std::ostream& out, const Rational& r) {
        out << r.numerator();
        if (!(r.numerator() == 0 || r.denominator() == 1))
            out << "/" << r.denominator();
        return out;
    }

// Operator +
    Rational operator+(const Rational &lhs, const Rational &rhs) {
        Rational copy(lhs);
        return copy += rhs;
    }

    Rational operator+(int lhs, const Rational &rhs) {
        Rational copy(lhs);
        return copy += rhs;
    }

// Operator -
    Rational operator-(const Rational &lhs, const Rational &rhs) {
        Rational copy(lhs);
        return copy -= rhs;
    }

    Rational operator-(int lhs, const Rational &rhs) {
        Rational copy(lhs);
        return copy -= rhs;
    }

// Operator *
    Rational operator*(const Rational &lhs, const Rational &rhs) {
        Rational copy(lhs);
        return copy *= rhs;
    }

    Rational operator*(int lhs, const Rational &rhs) {
        Rational copy(lhs);
        return copy *= rhs;
    }

// Operator /
    Rational operator/(const Rational &lhs, const Rational &rhs) {
        Rational copy(lhs);
        return copy /= rhs;
    }

    Rational operator/(int lhs, const Rational &rhs) {
        Rational copy(lhs);
        return copy /= rhs;
    }
}  // namespace rationalSpace

namespace std {
    rationalSpace::Rational abs(const rationalSpace::Rational& a) {
        if (a.numerator() > 0) {
            return -a;
        }
        return a;
    }
}

namespace matrixSpace {

    template<class T>
    class Matrix {
    public:
        Matrix<T>(const std::vector<std::vector<T>> &m) :
                _matrix(m) {}

        template<class U>
        Matrix<T>(const Matrix<U> &mtx) {
            auto sz = mtx.size();
            std::size_t n = sz.first, m = sz.second;

            _matrix.assign(n, std::vector<T>(m));
            for (size_t i = 0; i != n; ++i) {
                for (size_t j = 0; j != m; ++j) {
                    _matrix[i][j] = static_cast<T>(mtx.at(i, j));
                }
            }
        }

        // Unary operators
        Matrix<T> &operator+=(const Matrix<T> &oth) {
            auto x = size();
            std::size_t n = x.first;
            std::size_t m = x.second;

            for (std::size_t i = 0; i != n; ++i) {
                for (std::size_t j = 0; j != m; ++j) {
                    _matrix[i][j] += oth._matrix[i][j];
                }
            }
            return *this;
        }

        template<class U>
        Matrix<T> &operator*=(const U &oth) {
            for (auto &row : _matrix) {
                for (auto &elem : row) {
                    elem *= oth;
                }
            }
            return *this;
        }

        Matrix<T> &operator*=(const Matrix<T> &oth) {
            return *this = Matrix(*this) * oth;
        }


        Matrix<T> transposed() const {
            auto x = size();
            std::size_t n = x.first;
            std::size_t m = x.second;

            std::vector<std::vector<T>> res(m, std::vector<T>(n));
            for (std::size_t i = 0; i != n; ++i) {
                for (std::size_t j = 0; j != m; ++j) {
                    res[j][i] = at(i, j);
                }
            }
            return Matrix<T>(res);
        }

        Matrix<T> &transpose() {
            return *this = transposed();
        }

        // SOLVE BY GAUSS
        void forwardGauss(std::size_t step) {
            std::size_t n = size().first;
            if (step == n) return;

            T maxElem = at(step, step);
            std::size_t maxElemPos = step;
            for (size_t i = step + 1; i != n; ++i) {
                if (std::abs(at(i, step)) > maxElem) {
                    maxElemPos = i;
                    maxElem = std::abs(at(i, step));
                }
            }
            elemRow2(step, maxElemPos);
            if (at(step, step) != 0) {
                for (size_t i = step + 1; i != n; ++i) {
                    elemRow1<T>(i, step, -at(i, step) / at(step, step));
                }
            }
            forwardGauss(step + 1);
        }

        void backGauss(std::size_t step) {
            if (at(step, step) != 0) {
                elemRow3<T>(step, 1 / at(step, step));
                for (size_t i = 0; i != step; ++i) {
                    elemRow1<T>(i, step, -at(i, step));
                }
            }
            if (step == 0) return;
            backGauss(step - 1);
        }

        // Binary operators
        Matrix<T> operator+(const Matrix<T> &oth) const {
            Matrix<T> copy(*this);
            return copy += oth;
        }

        template<class U>
        Matrix<U> operator|(const std::vector<U> &b) const {
            Matrix<U> add(*this);
            for (std::size_t i = 0; i != size().first; ++i) {
                add.atV(i).push_back(b[i]);
            }
            return add;
        }

        // Getters
        std::pair<std::size_t, std::size_t> size() const {
            if (_matrix.empty())
                return {0, 0};
            return {_matrix.size(), _matrix[0].size()};
        }

        std::vector<T> &atV(std::size_t i) {
            return _matrix[i];
        }

        const std::vector<T> &atV(std::size_t i) const {
            return _matrix[i];
        }

        T &at(std::size_t i, std::size_t j) {
            return _matrix[i][j];
        }

        const T &at(std::size_t i, std::size_t j) const {
            return _matrix[i][j];
        }

    private:
        std::vector<std::vector<T>> _matrix;

        // elementary_operations on rows
        template<class U>
        void elemRow1(std::size_t i, std::size_t j, U lambda) {
            for (size_t k = 0; k != size().second; ++k) {
                _matrix[i][k] += lambda * _matrix[j][k];
            }
        }

        void elemRow2(std::size_t i, std::size_t j) {
            std::swap(_matrix[i], _matrix[j]);
        }

        template<class U>
        void elemRow3(std::size_t i, U lambda) {
            for (size_t k = 0; k != size().second; ++k) {
                _matrix[i][k] *= lambda;
            }
        }
    };

    template<class T>
    std::ostream &operator<<(std::ostream &out, const Matrix<T> &mtx) {
        auto x = mtx.size();
        std::size_t n = x.first;
        std::size_t m = x.second;

        for (std::size_t i = 0; i != n; ++i) {
            for (std::size_t j = 0; j != m; ++j) {
                out << mtx.at(i, j);

                if (j != m - 1)
                    out << "\t\t";
            }
            if (i != n - 1)
                out << '\n';
        }
        return out;
    }

    template<class T, class U>
    Matrix<T> operator*(const Matrix<T> &mtx, const U &oth) {
        Matrix<T> copy(mtx);
        return copy *= oth;
    }

    template<class T, class U>
    Matrix<T> operator*(const U &oth, const Matrix<T> &mtx) {
        Matrix<T> copy(mtx);
        return copy *= oth;
    }

    template<class T>
    Matrix<T> operator*(const Matrix<T> &lhs, const Matrix<T> &rhs) {
        auto sz1 = lhs.size();
        auto sz2 = rhs.size();
        std::size_t
                n1 = sz1.first,
                l = sz1.second,
                m2 = sz2.second;

        std::vector<std::vector<T>> res(n1, std::vector<T>(m2));
        for (std::size_t i = 0; i != n1; ++i) {
            for (std::size_t j = 0; j != m2; ++j) {
                for (std::size_t k = 0; k != l; ++k) {
                    res[i][j] += lhs.at(i, k) * rhs.at(k, j);
                }
            }
        }
        return Matrix<T>(res);
    }
}  //namespace matrixSpace

struct Interface {
    static void gauss() {
        std::size_t n;
        std::cout << "Size of matrix:\n";
        std::cin >> n;

        std::cout << "Matrix<int> NxN:\n";
        std::vector <std::vector <rationalSpace::Rational> > inputMtx
                (n, std::vector <rationalSpace::Rational>(n));
        for (auto& row : inputMtx) {
            for (auto& elem : row) {
                int x;
                std::cin >> x;
                elem = rationalSpace::Rational(x);
            }
        }
        matrixSpace::Matrix<rationalSpace::Rational> mtx(inputMtx);

        std::cout << "\nVector<int> N:\n";
        std::vector <rationalSpace::Rational> b(n);
        for (auto& elem : b) {
            int x;
            std::cin >> x;
            elem = rationalSpace::Rational(x);
        }

        mtx = mtx | b;
        mtx.forwardGauss(0);
        std::cout << "Row echelon form:\n"
                  << mtx << "\n\n";
        mtx.backGauss(n - 1);
        std::cout << "Row canonical form:\n"
                  << mtx << "\n\n";

    }

    static void mult() {
        size_t n1, m1;
        std::cout << "Sizes of first matrix:\n";
        std::cin >> n1 >> m1;

        std::cout << "Matrix<int> N1xM1:\n";
        std::vector <std::vector <rationalSpace::Rational> > inputMtx1
                (n1, std::vector <rationalSpace::Rational>(m1));
        for (auto& row : inputMtx1) {
            for (auto& elem : row) {
                int x;
                std::cin >> x;
                elem = rationalSpace::Rational(x);
            }
        }

        size_t n2, m2;
        std::cout << "\nSizes of second matrix:\n";
        std::cin >> n2 >> m2;

        std::cout << "Matrix<int> N2xM2:\n";
        std::vector <std::vector <rationalSpace::Rational> > inputMtx2
                (n2, std::vector <rationalSpace::Rational>(m2));
        for (auto& row : inputMtx2) {
            for (auto& elem : row) {
                int x;
                std::cin >> x;
                elem = rationalSpace::Rational(x);
            }
        }

        if (n2 != m1) {
            std::cout << "\nIncorrect sizes\n";
            return;
        }
        matrixSpace::Matrix<rationalSpace::Rational> mtx1(inputMtx1);
        matrixSpace::Matrix<rationalSpace::Rational> mtx2(inputMtx2);

        std::cout << "\nResult:\n" << (mtx1 * mtx2) << "\n\n";
    }

    static void sum() {
        size_t n, m;
        std::cout << "Sizes of matrix:\n";
        std::cin >> n >> m;

        std::cout << "First matrix<int> N1xM1:\n";
        std::vector <std::vector <rationalSpace::Rational> > inputMtx1
                (n, std::vector <rationalSpace::Rational>(m));
        for (auto& row : inputMtx1) {
            for (auto& elem : row) {
                int x;
                std::cin >> x;
                elem = rationalSpace::Rational(x);
            }
        }

        std::cout << "Second matrix<int> N1xM1:\n";
        std::vector <std::vector <rationalSpace::Rational> > inputMtx2
                (n, std::vector <rationalSpace::Rational>(m));
        for (auto& row : inputMtx1) {
            for (auto& elem : row) {
                int x;
                std::cin >> x;
                elem = rationalSpace::Rational(x);
            }
        }

        matrixSpace::Matrix<rationalSpace::Rational> mtx1(inputMtx1);
        matrixSpace::Matrix<rationalSpace::Rational> mtx2(inputMtx2);

        std::cout << "\nResult:\n" << (mtx1 + mtx2) << "\n\n";
    }

    static void trans() {
        std::size_t n;
        std::cout << "Size of matrix:\n";
        std::cin >> n;

        std::cout << "Matrix<int> NxN:\n";
        std::vector <std::vector <rationalSpace::Rational> > inputMtx
                (n, std::vector <rationalSpace::Rational>(n));
        for (auto& row : inputMtx) {
            for (auto& elem : row) {
                int x;
                std::cin >> x;
                elem = rationalSpace::Rational(x);
            }
        }
        matrixSpace::Matrix<rationalSpace::Rational> mtx(inputMtx);
        std::cout << "\nTransposed matrix:\n" << mtx.transposed() << "\n\n";
    }

    static void help(){
        std::cout << "Available commands:\n"
                  << "help -- open this help menu\n\n"
                  << "exit -- close the program\n\n"
                  << "gauss -- solve system of linear equation\n"
                  <<    "\tinput(size N, matrix<int> NxN, vector N)\n"
                  <<    "\toutput(matrix<rational> NxN in \"row echelon form\",\n\t\tmatrix<rational> NxN in \"row canonical form\")\n\n"
                  << "mult -- return multiplication of two matrix\n"
                  <<    "\tinput(size N1, size M1, matrix<int> N1xM1, size N2, size M2, matrix<int> N2xM2)\n"
                  <<    "\toutput(multiplication in matrix<int> N1xM2)\n\n"
                  << "sum -- return summary of two matrix\n"
                  <<    "\tinput(size N, size M, matrix<int> NxM, matrix<int> NxM)\n"
                  <<    "\toutput(sum in matrix<int> NxM)\n\n"
                  << "trans -- return transposed matrix\n"
                  <<    "\tinput(size N, matrix<int> NxN)\n"
                  <<    "\toutput((matrix<int> NxN)^T)\n\n";
    }
};

int main() {
    std::string command;
    std::cout << "Waiting for command...\n "
              << "Type \"help\" for show help menu\n";
    std::cin >> command;
    Interface events;
    while (command != "exit") {
        if (command == "gauss") {
            events.gauss();
        } else if (command == "mult") {
            events.mult();
        } else if (command == "sum") {
            events.sum();
        } else if (command == "trans") {
            events.trans();
        } else if (command == "help") {
            events.help();
        } else {
            std::cout << "Incorrect command\n";
        }
        std::cout << "Waiting for command...\n";
        std::cin >> command;
    }

}