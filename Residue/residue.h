#include<iostream>
#include<cmath>
#include<algorithm>
#include<cstring>

template<bool flag>
struct myassert {
    static void f() {};
};

template<>
struct myassert<0> {
    static void f() = delete;
};

template<unsigned N, unsigned lef, unsigned range>
struct calcSquare {
    static const unsigned value = calcSquare<N, static_cast<long long>((lef + range / 2)) * (lef + range / 2) <= N ? lef + range / 2 : lef,
            static_cast<long long>((lef + range / 2)) * (lef + range / 2) <= N ? (range + 1) / 2 : range / 2>::value;
};

template<unsigned N, unsigned lef>
struct calcSquare<N, lef, 1> {
    static const unsigned value = (lef % 2 == 0 ? lef - 1 : lef);
};

template<unsigned N>
struct Square {
    static const unsigned value = calcSquare<N, 1, N>::value;
};

template<unsigned N, unsigned div>
struct checkPrime {
    static const bool value = N % div == 0 ? false : checkPrime<N, div - 2>::value;
    static const unsigned maxdivisor = value ? N : (checkPrime<N, div - 2>::value ? div :
        checkPrime<N, div - 2>::maxdivisor);
};

template<unsigned N>
struct checkPrime<N, 1> {
    static const bool value = true;
    static const unsigned maxdivisor = N;
};

template<unsigned N>
struct is_prime {
    static const bool value = (N % 2 == 0 ? false : checkPrime<N, Square<N>::value>::value);
    static const unsigned maxdivisor = (N % 2 == 0 ? 2 : checkPrime<N, Square<N>::value>::maxdivisor);
};

template<>
struct is_prime<1> {
    static const bool value = false;
    static const unsigned maxdivisor = 1;
};

template<>
struct is_prime<2> {
    static const bool value = true;
    static const unsigned maxdivisor = 2;
};

template<unsigned N>
static const bool is_prime_v = is_prime<N>::value;

template<unsigned N, unsigned p>
struct checkPrimitive {
    static const bool value = (N % 2 == 0 ? false : checkPrimitive<N % p == 0 ? N / p : 0, p>::value);
};

template<unsigned p>
struct checkPrimitive<1, p> {
    static const bool value = true;
};

template<unsigned p>
struct checkPrimitive<0, p> {
    static const bool value = false;
};

template<unsigned N, int flag>
struct checkOdd {
    static const bool value = checkPrimitive<N % 4 == 0 ? 0 : (N % 2 == 0 ? N / 2 : N), N % 4 == 0 ? 0 : (N % 2 == 0 ? is_prime<N / 2>::maxdivisor : is_prime<N>::maxdivisor)>::value;
};

template<unsigned N>
struct checkOdd<N, 0> {
    static const bool value = false;
};

template<unsigned N>
struct has_primitive_root {
    static const bool value = checkOdd<N, N % 4>::value;
};

template<>
struct has_primitive_root<1> {
    static const bool value = true;
};

template<>
struct has_primitive_root<2> {
    static const bool value = true;
};

template<>
struct has_primitive_root<4> {
    static const bool value = true;
};

template<unsigned N>
static const bool has_primitive_root_v = has_primitive_root<N>::value;

unsigned CalcPhi(unsigned N) {
    unsigned res = N;
    for (unsigned div = 2; static_cast<long long>(div) * div <= N; ++div) {
        if (N % div == 0) {
            while (N % div == 0) {
                N /= div;
            }
            res -= res / div;
        }
    }
    if (N > 1) {
        res -= res / N;
    }
    return res;
}

template<unsigned N>
class Residue {
private:
    long long res;
    static const unsigned phi;
public:
    explicit Residue(int x);
    explicit operator int() const;

    bool operator==(const Residue<N>& y) const;
    bool operator!=(const Residue<N>& y) const;

    Residue<N>& operator+=(const Residue<N>& y);
    Residue<N>& operator-=(const Residue<N>& y);
    Residue<N>& operator*=(const Residue<N>& y);
    Residue<N>& operator/=(const Residue<N>& y);

    Residue<N> operator+(const Residue<N>& y) const;
    Residue<N> operator-(const Residue<N>& y) const;
    Residue<N> operator*(const Residue<N>& y) const;
    Residue<N> operator/(const Residue<N>& y) const;

    Residue<N> getInverse() const;
    static Residue<N> getPrimitiveRoot();
    Residue<N> pow(unsigned deg) const;
    Residue<N> pow(signed deg) const = delete;
    unsigned order() const;
};

template<unsigned N>
const unsigned Residue<N>::phi = CalcPhi(N);

template<unsigned N>
Residue<N>::Residue(int x) {
    res = x % static_cast<int>(N);
    if (res < 0)
        res += N;
}

template<unsigned N>
Residue<N>::operator int() const {
    return static_cast<int>(res);
}

template<unsigned N>
bool Residue<N>::operator==(const Residue<N>& y) const {
    return (res == y.res);
}

template<unsigned N>
bool Residue<N>::operator!=(const Residue<N>& y) const {
    return !(res == y.res);
}

template<unsigned N>
Residue<N>& Residue<N>::operator+=(const Residue<N>& y) {
    res += y.res;
    if (res >= N)
        res -= N;
    return *this;
}

template<unsigned N>
Residue<N>& Residue<N>::operator-=(const Residue<N>& y) {
    res -= y.res;
    if (res < 0)
        res += N;
    return *this;
}

template<unsigned N>
Residue<N>& Residue<N>::operator*=(const Residue<N>& y) {
    res = (res * y.res) % N;
    return *this;
}

template<unsigned N>
Residue<N>& Residue<N>::operator/=(const Residue<N>& y) {
    *this *= y.getInverse();
    return *this;
}

template<unsigned N>
Residue<N> Residue<N>::operator+(const Residue<N>& y) const {
    Residue<N> z = *this;
    z += y;
    return z;
}

template<unsigned N>
Residue<N> Residue<N>::operator-(const Residue<N>& y) const {
    Residue<N> z = *this;
    z -= y;
    return z;
}

template<unsigned N>
Residue<N> Residue<N>::operator*(const Residue<N>& y) const {
    Residue<N> z = *this;
    z *= y;
    return z;
}

template<unsigned N>
Residue<N> Residue<N>::operator/(const Residue<N>& y) const {
    Residue<N> z = *this;
    z /= y;
    return z;
}

template<unsigned N>
Residue<N> Residue<N>::getInverse() const {
    myassert<is_prime_v<N>>::f();
    Residue<N> rev = pow(phi - 1);
    return rev;
}

template<unsigned N>
Residue<N> Residue<N>::pow(unsigned deg) const {
    Residue<N> x = Residue<N>(1);
    Residue<N> a = *this;
    while (deg) {
        if (deg & 1)
            x *= a;
        a *= a;
        deg >>= 1;
    }
    return x;
}

template<unsigned N>
unsigned Residue<N>::order() const {
    unsigned mn = phi;
    for (unsigned div = 1; static_cast<long long>(div) * div <= phi; ++div) {
        if (phi % div != 0)
            continue;
        if (pow(div).res == 1) {
            if (div < mn)
                mn = div;
        }
        if (pow(phi / div).res == 1) {
            if (phi / div < mn)
                mn = phi / div;
        }
    }
    return mn;
}

template<unsigned N>
Residue<N> Residue<N>::getPrimitiveRoot() {
    myassert<has_primitive_root_v<N>>::f();
    for (unsigned x = 2; x <= N; ++x) {
        Residue<N> a(x);
        bool flag = true;
        if (std::__gcd(x, N) && a.pow(phi / 2).res != 1) {
            for (unsigned st = 2; static_cast<long long>(st) * st <= phi; ++st) {
                if (phi % st != 0)
                    continue;
                if (a.pow(st).res == 1) {
                    flag = false;
                    break;
                }
                if (a.pow(phi / st).res == 1) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                return a;
            }
        }
    }
    return Residue<N>(1);
}

