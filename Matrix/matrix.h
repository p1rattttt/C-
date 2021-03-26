#include<cassert>
#include<iostream>
#include<vector>
#include<string>
#include<complex>
#include<cmath>
#include<algorithm>
#include<cstring>

const double PI = acos(-1);

class BigInteger {
private:
    friend std::istream& operator>>(std::istream& in, BigInteger& x);
    friend class Rational;
    static const long long BASE = 10000;
    static const int LEN = 4;
    bool negative = 0;
    std::vector<long long> digit = {};
    void add(const BigInteger& y, bool sign);
public:
    const BigInteger operator-() const;
    const std::string toString() const;
    BigInteger& operator=(const BigInteger& x);
    BigInteger& operator+=(const BigInteger& x);
    BigInteger& operator-=(const BigInteger& x);
    BigInteger& operator++();
    BigInteger& operator--();
    BigInteger& operator*=(const BigInteger& x);
    BigInteger& operator/=(const BigInteger& x);
    BigInteger& operator%=(const BigInteger& x);
    const BigInteger operator++(int);
    const BigInteger operator--(int);
    void swap(BigInteger& x);
    void mul(int x);
    void div(int x);

    size_t Size() const {
        return digit.size();
    }

    BigInteger(const int x) {
        int y = x;
        if (y < 0) {
            y *= -1;
            negative = 1;
        };
        do {
            digit.push_back(y % BASE);
            y /= BASE;
        } while (y != 0);
    }

    bool& sign() {
        return negative;
    }

    bool sign() const {
        return negative;
    }

    ~BigInteger() {
        negative = 0;
        digit.clear();
        digit.shrink_to_fit();
    }

    long long operator[](const int ind) const {
        return digit[ind];
    }

    long long& operator[](const int ind) {
        return digit[ind];
    }

    BigInteger(const BigInteger& x) {
        negative = x.negative;
        digit = x.digit;
    }

    BigInteger() {
        negative = 0;
        digit = {0};
    }

    void normalize() {
        while (Size() > 1 && digit.back() == 0)
            digit.pop_back();
        if (Size() == 1 && digit.back() == 0)
            negative = 0;
    }

    explicit operator bool() const {
        if (Size() > 1)
            return 1;
        return digit[0] != 0;
    }

    explicit operator int() const;
};

void BigInteger::swap(BigInteger& x) {
    digit.swap(x.digit);
    std::swap(x.negative, negative);
}

std::ostream& operator<<(std::ostream& out, const BigInteger& x) {
    out << x.toString();
    return out;
}

void BigInteger::add(const BigInteger& x, bool sign) {
    long long carry = 0;
    for (size_t i = 0; i < std::min(Size(), x.Size()) || carry; ++i) {
        if (i >= Size()) {
            digit.push_back(0);
        }
        digit[i] += (sign ? -1 : 1) * carry;
        if (i < x.Size()) {
            digit[i] += (sign ? -1 : 1) * x[i];
        }
        carry = 0;
        if (digit[i] >= BASE) {
            carry = 1;
            digit[i] -= BASE;
        }
        if (digit[i] < 0) {
            carry = 1;
            digit[i] += BASE;
        }
    }
    normalize();
}

std::istream& operator>>(std::istream& in, BigInteger& x) {
    std::string s;
    in >> s;
    x.negative = 0;
    x.digit.clear();
    bool negative = 0;
    int zero = 0;
    if (s[0] == '-') {
        negative = 1;
        zero = 1;
    }
    if (s[0] == '+') {
        zero = 1;
    }
    x.negative = negative;
    for (int i = s.size(); i > zero; i -= x.LEN) {
        if (i >= x.LEN + zero)
            x.digit.push_back(atoi(s.substr(i - x.LEN, x.LEN).c_str()));
        else
            x.digit.push_back(atoi(s.substr(zero, i - zero).c_str()));
    }
    if (!x.Size()) {
        x.digit.push_back(0);
    }
    x.normalize();
    return in;
}

void BigInteger::mul(int x) {
    int carry = 0;
    for (size_t i = 0; i < Size() || carry; ++i) {
        if (i == Size())
            digit.push_back (0);
        long long cur = carry + digit[i] * (long long)x;
        digit[i] = cur % BASE;
        carry = cur / BASE;
    }
    normalize();
}

void BigInteger::div(int x) {
    int carry = 0;
    for (int i = (int)Size() - 1; i >= 0; --i) {
        long long cur = digit[i] + (long long)carry * BASE;
        digit[i] = cur / x;
        carry = cur % x;
    }
    normalize();
}

bool operator<(const BigInteger& y, const BigInteger& x) {
    if (y.sign() && !x.sign())
        return true;
    if (!y.sign() && x.sign())
        return false;
    bool f = 0;
    bool eq = 1;
    if (y.Size() == x.Size()) {
        for (int i = static_cast<int>(x.Size()) - 1; i >= 0; --i) {
            if (y[i] < x[i]) {
                f = 1;
                eq = 0;
                break;
            }
            if (y[i] > x[i]) {
                eq = 0;
                break;
            }
        }
    } else {
        eq = 0;
        if (y.Size() < x.Size()) {
            f = 1;
        }
    }
    if (eq)
        return false;
    if (y.sign())
        f ^= 1;
    return f;
}

bool operator==(const BigInteger& y, const BigInteger& x) {
    return !(y < x) && !(x < y);
}

bool operator!=(const BigInteger& y, const BigInteger& x) {
    return !(y == x);
}

bool operator>(const BigInteger& y, const BigInteger& x) {
    return x < y;
}

bool operator<=(const BigInteger& y, const BigInteger& x) {
    return (y < x) || (y == x);
}

bool operator>=(const BigInteger& y, const BigInteger& x) {
    return (y > x) || (y == x);
}

const BigInteger BigInteger::operator-() const {
    BigInteger x(*this);
    x.sign() ^= 1;
    x.normalize();
    return x;
}

const BigInteger operator-(const BigInteger& x, const BigInteger& y);
const BigInteger operator+(const BigInteger& x, const BigInteger& y);

const std::string BigInteger::toString() const {
    std::string ans = "";
    if (negative)
        ans += '-';
    for (int i = static_cast<int>(Size()) - 1; i >= 0; --i) {
        std::string cur = std::to_string(digit[i]);
        if (i != static_cast<int>(Size()) - 1) {
            int need = LEN - cur.size();
            while(need--)
                ans += '0';
        }
        ans += cur;
    }
    return ans;
}

BigInteger& BigInteger::operator=(const BigInteger& x) {
    if (this == &x)
        return *this;
    BigInteger y = x;
    swap(y);
    return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& x) {
    if (x.negative) {
        return *this = (*this) + (-x);
    } else if ((*this).negative) {
        return *this = -((-(*this)) + x);
    } else if ((*this) < x) {
        return *this = -(x - (*this));
    }
    add(x, 1);
    return *this;
}

const BigInteger operator-(const BigInteger& x, const BigInteger& y) {
    BigInteger ans = x;
    ans -= y;
    return ans;
}

BigInteger& BigInteger::operator+=(const BigInteger& x) {
    if (negative) {
        if (x.negative) {
            return (*this) = -((-(*this)) + (-x));
        }
        return *this = x - (-(*this));
    } else if (x.negative) {
        return *this = (*this) - (-x);
    }
    if ((*this) < x) {
        return *this = x + (*this);
    }
    add(x, 0);
    return (*this);
}

const BigInteger operator+(const BigInteger& x, const BigInteger& y) {
    BigInteger ans = x;
    ans += y;
    return ans;
}

BigInteger& BigInteger::operator++() {
    *this += 1;
    return *this;
}

BigInteger& BigInteger::operator--() {
    *this -= 1;
    return *this;
}

const BigInteger BigInteger::operator++(int) {
    *this += 1;
    return (*this - 1);
}

const BigInteger BigInteger::operator--(int) {
    *this -= 1;
    return (*this + 1);
}

BigInteger operator "" _bi(unsigned long long x) {
    return BigInteger(x);
}

class Complex {
public:
    long double real = 0, image = 0;

    Complex& operator*=(const Complex& y) {
        Complex it(real * y.real - image * y.image, real * y.image + image * y.real);
        *this = it;
        return *this;
    }

    Complex& operator/=(size_t n) {
        real /= (long double)n;
        image /= (long double)n;
        return *this;
    }

    Complex operator*(const Complex& y) {
        return Complex(real * y.real - image * y.image, real * y.image + image * y.real);
    }

    Complex operator+(const Complex& y) {
        return Complex(real + y.real, image + y.image);
    }

    Complex operator-() {
        return Complex(-real, -image);
    }

    Complex operator-(const Complex& y) {
        return Complex(real - y.real, image - y.image);
    }
    Complex() : real(0), image(0) {}
    Complex(long double real, long double image): real(real), image(image) {}

    Complex(long double real): real(real), image(0) {}
};

const BigInteger operator*(const BigInteger& x, const BigInteger& y);

void fft(std::vector<Complex>& a, bool invert) {
    size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            std::swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        long double ang = 2 * PI / len * (invert ? -1 : 1);
        Complex wlen(cosl(ang), sinl(ang));
        for (size_t i = 0; i < n; i += len) {
            Complex w(1.0);
            for (size_t j = 0; j < len / 2; j++) {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (invert) {
        for (Complex& x : a) {
            x /= n;
        }
    }
}

BigInteger& BigInteger::operator*=(const BigInteger& y) {
    bool nsign = sign() ^ y.sign();
    std::vector<Complex> fx(digit.begin(), digit.end());
    std::vector<Complex> fy(y.digit.begin(), y.digit.end());
    size_t n = 1;
    while (n < std::max(fx.size(), fy.size())) {
        n <<= 1;
    }
    n <<= 1;
    fx.resize(n);
    fy.resize(n);
    fft(fx, 0);
    fft(fy, 0);
    for (size_t i = 0; i < n; ++i) {
        fx[i] *= fy[i];
    }
    fft(fx, 1);
    long long carry = 0;
    for (size_t i = 0; i < n || carry; ++i) {
        if (i >= Size()) {
            digit.push_back(0);
        }
        long long cur = 0;
        cur = carry;
        if (i < fx.size()) {
            cur += static_cast<long long>(fx[i].real + 0.5);
        }
        carry = 0;
        if (cur >= BASE) {
            carry = cur / BASE;
            cur %= BASE;
        }
        digit[i] = cur;
    }
    sign() = nsign;
    normalize();
    return *this;
}

const BigInteger operator*(const BigInteger& x, const BigInteger& y) {
    BigInteger z = x;
    z *= y;
    return z;
}

const BigInteger operator/(const BigInteger& x, const BigInteger& y) {
    BigInteger ans = x;
    ans /= y;
    return ans;
}

BigInteger& BigInteger::operator/=(const BigInteger& x) {
    bool f = 0;
    if (negative != x.negative) {
        f = 1;
    }
    negative = 0;
    BigInteger copy = x;
    copy.negative = 0;
    if (*this < copy) {
        *this = 0;
        return *this;
    }
    std::vector<int> alls(Size() * LEN);
    for (size_t i = 0, j = 0; i < alls.size(); i += LEN, ++j) {
        int cur = digit[j];
        int st = 10;
        for (int j = 0; j < LEN; ++j) {
            alls[i + j] = cur % st;
            cur /= st;
        }
    }
    while (alls.size() > 1 && alls.back() == 0) {
        alls.pop_back();
    }
    std::vector<int> all(alls.size());
    for (size_t i = 0; i < all.size(); ++i) {
        all[i] = alls[alls.size() - i - 1];
    }
    *this = 0;
    BigInteger cur;
    for (size_t i = 0; i < all.size(); ++i) {
        cur.mul(10);
        (*this).mul(10);
        cur += all[i];
        while (cur >= copy) {
            cur -= copy;
            ++(*this);
        }
    }
    negative = f;
    normalize();
    return *this;
}

const BigInteger operator%(const BigInteger& x, const BigInteger& y) {
    BigInteger mod = x;
    mod %= y;
    return mod;
}

BigInteger& BigInteger::operator%=(const BigInteger& x) {
    if (negative) {
        if (x.negative) {
            return *this = (-(*this)) % (-x);
        }
        return *this = -((-(*this)) % x);
    }
    if (x.negative) {
        return *this = -((*this) % (-x));
    }
    *this = *this - ((*this) / x) * x;
    return *this;
}

class Rational {
private:
    friend bool operator<(const Rational& y, const Rational& x);
    BigInteger num;
    BigInteger den;

    BigInteger gcd(BigInteger x, BigInteger y) {
        if (x < y) x.swap(y);
        int cnt2 = 0;
        while (y != 0) {
            if (!(x[0] & 1) && !(y[0] & 1)) {
                x.div(2);
                y.div(2);
                ++cnt2;
                continue;
            }
            if (!(x[0] & 1)) {
                x.div(2);
                if (x < y) x.swap(y);
                continue;
            }
            if (!(y[0] & 1)) {
                y.div(2);
                continue;
            }
            x -= y;
            if (x < y) x.swap(y);
        }
        for (int i = 0; i < cnt2; ++i) x.mul(2);
        return x;
    }
public:
    const Rational operator-();
    const std::string toString() const;
    void swap(Rational& x);
    Rational& operator=(const Rational& x);
    Rational& operator+=(const Rational& x);
    Rational& operator-=(const Rational& x);
    Rational& operator*=(const Rational& x);
    Rational& operator/=(const Rational& x);
    std::string asDecimal(const size_t precision = 0) const;

    Rational() {
        num = 0;
        den = 1;
    }

    Rational(const int x) {
        num = x;
        den = 1;
    }

    Rational(const BigInteger& x) {
        num = x;
        den = 1;
    }

    Rational(const Rational& x) {
        num = x.num;
        den = x.den;
    }

    void normalize() {
        BigInteger x;
        bool f = num.negative ^ den.negative;
        num.negative = f;
        den.negative = 0;
        if (num.negative) {
            x = gcd(-num, den);
        } else {
            x = gcd(num, den);
        }
        num /= x;
        den /= x;
        if (num == 0) {
            den = 1;
            num.negative = 0;
        }
    }

    explicit operator double() const {
        long double x = 0;
        x = stod(asDecimal(9));
        return x;
   }
};

const Rational Rational::operator-() {
    Rational x = *this;
    x.num.negative ^= 1;
    x.normalize();
    return x;
}

void Rational::swap(Rational& x) {
    std::swap(num, x.num);
    std::swap(den, x.den);
}

Rational& Rational::operator=(const Rational& x) {
    if (this == &x)
        return *this;
    Rational y = x;
    swap(y);
    return *this;
}

std::ostream& operator<<(std::ostream& out, const Rational& x) {
    out << x.asDecimal(9);
    return out;
}

std::istream& operator>>(std::istream& in, Rational& x) {
    int y;
    in >> y;
    x = Rational(y);
    return in;
}

bool operator<(const Rational& y, const Rational& x) {
    return y.num * x.den < y.den * x.num;
}

bool operator==(const Rational& y, const Rational& x) {
    return !(y < x) && !(x < y);
}


bool operator!=(const Rational& y, const Rational& x) {
    return !(y == x);
}

bool operator>(const Rational& y, const Rational& x) {
    return x < y;
}

bool operator<=(const Rational& y, const Rational& x) {
    return y < x || y == x;
}

bool operator>=(const Rational& y, const Rational& x) {
    return y > x || y == x;
}

const std::string Rational::toString() const {
    //assert(0);
    std::string ans = "";
    ans += num.toString();
    if (den != 1) {
        ans += '/';
        ans += den.toString();
    }
    return ans;
}

Rational& Rational::operator*=(const Rational& x) {
    num *= x.num;
    den *= x.den;
    normalize();
    return *this;
}

Rational& Rational::operator/=(const Rational& x) {
    num *= x.den;
    den *= x.num;
    normalize();
    return *this;
}

Rational& Rational::operator+=(const Rational& x) {
    num *= x.den;
    num += den * x.num;
    den *= x.den;
    normalize();
    return *this;
}

Rational& Rational::operator-=(const Rational& x) {
    num *= x.den;
    num -= den * x.num;
    den *= x.den;
    normalize();
    return *this;
}

const Rational operator*(const Rational& x, const Rational& y) {
    Rational z = x;
    z *= y;
    return z;
}

const Rational operator/(const Rational& x, const Rational& y) {
    Rational z = x;
    z /= y;
    return z;
}

const Rational operator+(const Rational& x, const Rational& y) {
    Rational z = x;
    z += y;
    return z;
}

const Rational operator-(const Rational& x, const Rational& y) {
    Rational z = x;
    z -= y;
    return z;
}

std::string Rational::asDecimal(const size_t precision) const {
    //assert(0);
    if (num == 0) {
        std::string ans = "0";
        if (precision != 0) {
            ans += ".";
        }
        for (size_t i = 0; i < precision; ++i) {
            ans += '0';
        }
        return ans;
    }
    Rational x = *this;
    BigInteger cur;
    std::string ans = "";
    if (x.num.negative) {
        ans += "-";
    }
    x.num.negative = 0;
    for (size_t i = 0; i < precision; ++i) {
        cur = x.num / x.den;
        x -= cur;
        ans += cur.toString();
        x.num.mul(10);
        if (i == 0 && precision != 0) {
            ans += '.';
        }
    }
    return ans;
}

//using Rational = double;

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
    template<unsigned M>
    friend std::ostream& operator<<(std::ostream& out, const Residue<M>& x);
    long long res;
    static const unsigned phi;
public:
    explicit Residue(int x = 0);
    explicit operator int() const;

    bool operator==(const Residue<N>& y) const;
    bool operator!=(const Residue<N>& y) const;

    Residue<N> operator-() const;
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
Residue<N> Residue<N>::operator-() const {
    return Residue<N>(N - res);
}

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
std::ostream& operator<<(std::ostream& out, const Residue<N>& x) {
    out << x.res;
    return out;
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

template<unsigned N, unsigned M, typename Field = Rational>
class Matrix {
private:
    std::vector<std::vector<Field>> arr;
    bool Gauss();
    void SwapRows(unsigned i, unsigned j); // swap of rows
    void AddRow(unsigned i, unsigned j, Field k); // add row_i * k to row_j
    void MultRow(unsigned i, Field k); // multiply row_i by k
    int curRank;
public:
    Matrix<N, M, Field>& operator+=(const Matrix<N, M, Field>& y);
    Matrix<N, M, Field>& operator-=(const Matrix<N, M, Field>& y);

    Matrix<N, M, Field>& operator*=(const Field& y);
    Matrix<N, M, Field>& operator*=(const Matrix<M, M, Field>& y);

    Field det() const;
    Matrix<M, N, Field> transposed() const;
    Matrix<N, N, Field> inverted() const;
    void invert();
    int rank() const;
    Field trace() const;

    std::vector<Field> getRow(unsigned n) const;
    std::vector<Field> getColumn(unsigned m) const;

    std::vector<Field>& operator[](unsigned ind) {
        return arr[ind];
    }

    const std::vector<Field>& operator[](unsigned ind) const {
        return arr[ind];
    }

    Matrix() {
        arr.resize(N);
        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < M; ++j) {
                arr[i].push_back(Field(0));
            }
        }
    }
    Matrix(const std::vector<std::vector<Field>>& a) : arr(a) {}
    Matrix(const std::vector<std::vector<int>>& a) {
        arr.resize(N);
        for (unsigned i = 0; i < N; ++i) {
            arr[i].resize(M);
            for (unsigned j = 0; j < M; ++j) {
                arr[i][j] = Field(a[i][j]);
            }
        }
    }
    Matrix(const std::initializer_list<std::vector<int>>& a) {
        arr.resize(N);
        unsigned ind = 0;
        for (auto i : a) {
            for (unsigned j = 0; j < M; ++j) {
                arr[ind].push_back(Field(i[j]));
            }
            ++ind;
        }
    }
};

template<unsigned N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator+=(const Matrix<N, M, Field>& y) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            arr[i][j] += y.arr[i][j];
        }
    }
    return *this;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator-=(const Matrix<N, M, Field>& y) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            arr[i][j] -= y.arr[i][j];
        }
    }
    return *this;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator+(const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
    Matrix<N, M, Field> ans = x;
    ans += y;
    return ans;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator-(const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
    Matrix<N, M, Field> ans = x;
    ans -= y;
    return ans;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=(const Field& y) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            arr[i][j] *= y;
        }
    }
    return *this;
}

template<unsigned N, unsigned M, typename Field>
Matrix<M, N, Field> Matrix<N, M, Field>::transposed() const {
    Matrix<M, N, Field> ans;
    for (unsigned i = 0; i < M; ++i) {
        for (unsigned j = 0; j < N; ++j) {
            ans[i][j] = arr[j][i];
        }
    }
    return ans;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator*(const Matrix<N, M, Field>& x, const Field& y) {
    Matrix<N, M, Field> ans = x;
    ans *= y;
    return ans;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator*(const Field& y, const Matrix<N, M, Field>& x) {
    return x * y;
}

template<unsigned N, unsigned M, unsigned L, typename Field>
Matrix<N, L, Field> operator*(const Matrix<N, M, Field>& x, const Matrix<M, L, Field>& y) {
    Matrix<N, L, Field> ans;
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < L; ++j) {
            for (unsigned z = 0; z < M; ++z) {
                ans[i][j] += x[i][z] * y[z][j];
            }
        }
    }
    return ans;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=(const Matrix<M, M, Field>& y) {
    static_assert(N == M, "bruh");
    return *this = *this * y;
}

template<unsigned N, unsigned M, typename Field>
std::vector<Field> Matrix<N, M, Field>::getRow(unsigned n) const {
    return arr[n];
}

template<unsigned N, unsigned M, typename Field>
std::vector<Field> Matrix<N, M, Field>::getColumn(unsigned m) const {
    std::vector<Field> ans(N);
    for (unsigned i = 0; i < N; ++i)
        ans[i] = arr[i][m];
    return ans;
}

template<unsigned N, unsigned M, typename Field>
Field Matrix<N, M, Field>::trace() const {
    static_assert(N == M, "Kostya loh");
    Field ans;
    for (unsigned i = 0; i < N; ++i) {
        ans += arr[i][i];
    }
    return ans;
}

template<unsigned N, unsigned M, typename Field>
Field Matrix<N, M, Field>::det() const {
    static_assert(N == M, "Kostya loh");
    Matrix<N, M, Field> copy = *this;
    bool inv = copy.Gauss();
    Field ans = Field(1);
    if (inv) ans *= static_cast<Field>(-1);
    for (unsigned i = 0; i < N; ++i) {
        ans *= copy[i][i];
    }
    return ans;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, N, Field> Matrix<N, M, Field>::inverted() const {
    static_assert(N == M, "Artur loh");
    Matrix<N, M, Field> copy = *this;
    copy.invert();
    return copy;
}

template<unsigned N, unsigned M, typename Field>
void Matrix<N, M, Field>::invert() {
    static_assert(N == M, "Artur loh");
    Matrix<N, N, Field> ans;
    for (unsigned i = 0; i < N; ++i) {
        ans[i][i] = Field(1);
    }
//    for (unsigned i = 0; i < N; ++i) {
//        for (unsigned j = 0; j < N; ++j) {
//            std::cerr << arr[i][j] << ' ';
//        }
//        std::cerr << '\n';
//    }
    for (unsigned row = 0, col = 0; row < N && col < M; ++col) {
        bool flag = false;
        if (arr[row][col] == Field(0)) {
            flag = true;
            for (unsigned next = row + 1; next < N; ++next) {
                if (arr[next][col] != Field(0)) {
                    SwapRows(row, next);
                    ans.SwapRows(row, next);
                    flag = false;
                    break;
                }
            }
        }
//        std::cout << row << ' ' << col << std::endl;
        if (flag)
            continue;
        for (unsigned next = row + 1; next < N; ++next) {
            if (arr[next][col] == Field(0))
                continue;
            Field k;
            AddRow(row, next, k = (-arr[next][col] / arr[row][col]));
            ans.AddRow(row, next, k);
        }
        ++row;
    }
//    return;
//    static_assert(row == N, "Artur loh");
    for (int i = N - 1; i >= 0; --i) {
        ans.MultRow(i, arr[i][i]);
        MultRow(i, arr[i][i]);
        for (int j = i - 1; j >= 0; --j) {
            if (arr[j][i] == Field(0))
                continue;
            Field k;
            AddRow(i, j, k = (-arr[j][i] / arr[i][i]));
            ans.AddRow(i, j, k);
        }
    }
    *this = ans;
}

template<unsigned N, unsigned M, typename Field>
int Matrix<N, M, Field>::rank() const {
    Matrix<N, M, Field> copy = *this;
    copy.Gauss();
    return copy.curRank;
}

template<unsigned N, unsigned M, typename Field>
void Matrix<N, M, Field>::SwapRows(unsigned i, unsigned j) {
    swap(arr[i], arr[j]);
}

template<unsigned N, unsigned M, typename Field>
void Matrix<N, M, Field>::AddRow(unsigned i, unsigned j, Field k) {
    for (unsigned z = 0; z < M; ++z) {
        arr[j][z] += arr[i][z] * k;
    }
}

template<unsigned N, unsigned M, typename Field>
void Matrix<N, M, Field>::MultRow(unsigned int i, Field k) {
    for (unsigned j = 0; j < M; ++j) {
        arr[i][j] /= k;
    }
}

template<unsigned N, unsigned M, typename Field>
bool Matrix<N, M, Field>::Gauss() {
    curRank = 0;
    bool inv = 0;
    for (unsigned row = 0, col = 0; row < N && col < M; ++col) {
        bool flag = 0;
        if (arr[row][col] == Field(0)) {
            flag = 1;
            for (unsigned next = row + 1; next < N; ++next) {
                if (arr[next][col] != Field(0)) {
                    SwapRows(row, next);
                    inv ^= 1;
                    flag = 0;
                    break;
                }
            }
        }
        if (flag)
            continue;
        for (unsigned next = row + 1; next < N; ++next) {
            if (arr[next][col] == Field(0))
                continue;
            AddRow(row, next, -arr[next][col] / arr[row][col]);
        }
        ++row;
        ++curRank;
    }
    return inv;
}

template<unsigned N, unsigned M, typename Field>
bool operator==(const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            if (x[i][j] != y[i][j]) {
                return false;
            }
        }
    }
    return true;
}

template<unsigned N, unsigned M, typename Field>
bool operator!=(const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
    return !(x == y);
}

