#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>

template <typename T>
class Polynomial {
private:
    std::vector<T> m_coefs;

public:
    Polynomial(const std::vector<T> &coefs) : m_coefs(coefs) {}

    Polynomial(const T value = T()) {
        m_coefs.push_back(value);
    }

    template<class Iter>
    Polynomial(Iter f, Iter l) {
        while (f != l) {
            m_coefs.push_back(*f);
            ++f;
        }
    }

    int Degree() const {
        int i = m_coefs.size() - 1;
        while (i > -1 && m_coefs[i] == T(0)) {
            --i;
        }
        return i;
    }

    Polynomial<T> &operator+=(const Polynomial<T> &polinom) {
        int new_degree = max(this->Degree(), polinom.Degree());
        m_coefs.resize(new_degree + 1);
        for (int i = 0; i <= polinom.Degree(); ++i) {
            m_coefs[i] += polinom[i];
        }
        return *this;
    }

    friend Polynomial<T> operator+(Polynomial<T> polinom1, const Polynomial<T> &polinom2) {
        return polinom1 += polinom2;
    }

    Polynomial<T> &operator-=(const Polynomial<T> &polinom) {
        int new_degree = max(this->Degree(), polinom.Degree());
        m_coefs.resize(new_degree + 1);
        for (int i = 0; i < polinom.Degree() + 1; ++i) {
            m_coefs[i] -= polinom.m_coefs[i];
        }
        return *this;
    }


    friend Polynomial<T> operator-(Polynomial<T> polinom1, const Polynomial<T> &polinom2) {
        return polinom1 -= polinom2;
    }

    Polynomial<T> &operator*=(const Polynomial &polinom) {
        std::vector<T> m_new(this->Degree() + polinom.Degree() + 1);
        for (int i = 0; i <= this->Degree(); ++i) {
            for (int j = 0; j <= polinom.Degree(); ++j) {
                m_new[i + j] += this->operator[](i) * polinom[j];
            }
        }
        m_coefs = m_new;
        return *this;
    }

    friend Polynomial<T> operator*(Polynomial<T> polinom1, const Polynomial<T> &polinom2) {
        return polinom1 *= polinom2;
    }

    bool operator==(const Polynomial<T> &polinom) const {
        if (this->Degree() != polinom.Degree()) {
            return false;
        }
        for (int i = 0; i <= this->Degree(); ++i) {
            if (this->m_coefs[i] != polinom.m_coefs[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const Polynomial<T> &polinom) const {
        return !(*this == polinom);
    }

    T operator[](int i) const {
        if (i > (*this).Degree()) {
            return T(0);
        } else {
            return m_coefs[i];
        }
    }

    T operator()(const T& x) const {
        T result = T(0);
        T power = T(1);
        for (int i = 0; i <= Degree(); ++i) {
            result += power * this->operator[](i);
            power *= x;
        }
        return result;
    }


    typename std::vector<T>::const_iterator begin() const {
        return m_coefs.begin();
    }

    typename std::vector<T>::const_iterator end() const {
        return m_coefs.end();
    }
/*    vector<T> get_coefs_for_test() {
        return m_coefs;
    }*/

    friend std::ostream& operator<<(std::ostream& out, Polynomial polinom) {
        if (polinom.Degree() == -1) {
            out << polinom[0];
            return out;
        }

        bool plus = false;
        for (int i = polinom.Degree(); i >= -1; --i) {
            if (polinom[i] == T(0)) {
                if (i == 0) {
                    break;
                }
                continue;
            }
            if (polinom[i] == T(-1)) {
                if (i == 0) {
                    out << -1;
                } else {
                    out << '-';
                }
            }
            if (polinom[i] == T(1)) {
                if (plus) {
                    out << '+';
                }
                if (i == 0) {
                    out << 1;
                }
            }
            if (polinom[i] != T(1) && polinom[i] != T(-1)) {
                if (plus && polinom[i] > T(0)) {
                    out << '+';
                }
                out << polinom[i];
                if (i != 0) {
                    out << '*';
                }
            }
            if (i != 0) {
                out << 'x';
            }
            if (i > 1) {
                out << '^' << i;
            }
            plus = true;
            if (i == 0) {
                break;
            }
        }
        return out;
    }

    friend Polynomial<T> operator&(const Polynomial<T> polinom1, const Polynomial<T>& polinom2) {
        Polynomial<T> cur(T(0)), result;
        for (int i = 0; i <= polinom1.Degree(); ++i) {
            cur = Polynomial<T>(polinom1[i]);
            for (int j = 0; j < i; ++j) {
                cur *= polinom2;
            }
            result += cur;
        }
        return result;
    }

    friend Polynomial<T> operator/(const Polynomial<T>& polinom1, const Polynomial<T>& polinom2) {
        if (polinom2.Degree() != -1) {
            Polynomial<T> p1 = polinom1;
            std::vector<T> coef(max(0, p1.Degree() - polinom2.Degree() + 1)),
                    m(max(0, p1.Degree() - polinom2.Degree() + 1), T(0));
            Polynomial<T> monomial(m);
            int cur_i;
            while (p1.Degree() >= polinom2.Degree()) {
                cur_i = p1.Degree() - polinom2.Degree();
                coef[cur_i] = p1[p1.Degree()] / polinom2[polinom2.Degree()];
                m[cur_i] = p1[p1.Degree()] / polinom2[polinom2.Degree()];
                monomial = Polynomial<T>(m);
                p1 -= (monomial * polinom2);
                m[cur_i] = T(0);
            }
            return Polynomial<T>(coef);
        } else {
            return T(0);
        }
    }
    friend Polynomial<T> operator%(const Polynomial<T>& polinom1, const Polynomial<T>& polinom2) {
        return polinom1 - (polinom1 / polinom2) * polinom2;
    }

    friend Polynomial<T> operator,(const Polynomial<T>& polinom1, const Polynomial<T>& polinom2) {
        if (polinom1 == T(0)) {
            if (polinom2.Degree() != -1) {
                return polinom2 / polinom2[polinom2.Degree()];
            }
            return polinom2;
        }
        return (polinom2 % polinom1, polinom2);
    }
};
