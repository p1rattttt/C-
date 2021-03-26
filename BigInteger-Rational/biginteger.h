#include<cassert>
#include<iostream>
#include<vector>
#include<string>
#include<complex>
#include<cmath>

const double PI = acos(-1);

class BigInteger {
private:
    friend std::istream& operator>>(std::istream& in, BigInteger& x);
    friend class Rational;
    static const long long BASE = 10000;
    static const int LEN = 4;
    bool negative = 0;
    std::vector<long long> digit = {};
    void add(const BigInteger& x, bool sign);
    BigInteger abs() const;
public:
    BigInteger operator-() const;
    std::string toString() const;
    BigInteger& operator=(const BigInteger& x);
    BigInteger& operator+=(const BigInteger& x);
    BigInteger& operator-=(const BigInteger& x);
    BigInteger& operator++();
    BigInteger& operator--();
    BigInteger& operator*=(const BigInteger& x);
    BigInteger& operator/=(const BigInteger& x);
    BigInteger& operator%=(const BigInteger& x);
    BigInteger operator++(int);
    BigInteger operator--(int);
    void swap(BigInteger& x);

    size_t getLen() const {
        return digit.size();
    }

    BigInteger(int x) {
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

    bool sign() const {
        return negative;
    }

    ~BigInteger() {
        negative = 0;
        digit.clear();
        digit.shrink_to_fit();
    }

    long long operator[](int ind) const {
        return digit[ind];
    }

    long long& operator[](int ind) {
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
        while (getLen() > 1 && digit.back() == 0)
            digit.pop_back();
        if (getLen() == 1 && digit.back() == 0)
            negative = 0;
    }

    explicit operator bool() const {
        if (getLen() > 1)
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

BigInteger BigInteger::abs() const {
    BigInteger ans(*this);
    ans.negative = 0;
    return ans;
}

void BigInteger::add(const BigInteger& x, bool sign) {
    long long carry = 0;
    for (size_t i = 0; i < std::min(getLen(), x.getLen()) || carry; ++i) {
        if (i >= getLen()) {
            digit.push_back(0);
        }
        digit[i] += (sign ? -1 : 1) * carry;
        if (i < x.getLen()) {
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
    if (!x.getLen()) {
        x.digit.push_back(0);
    }
    x.normalize();
    return in;
}

bool operator<(const BigInteger& y, const BigInteger& x) {
    if (y.sign() != x.sign())
        return y.sign();
    if (y.getLen() != x.getLen()) {
        return !y.sign() == (y.getLen() < x.getLen());
    }
    for (int i = y.getLen() - 1; i >= 0; --i) {
        if (y[i] != x[i]) {
            return !y.sign() == (y[i] < x[i]);
        }
    }
    return false;
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
    return !(x < y);
}

bool operator>=(const BigInteger& y, const BigInteger& x) {
    return !(y < x);
}

BigInteger BigInteger::operator-() const {
    BigInteger x(*this);
    x.negative ^= 1;
    x.normalize();
    return x;
}

BigInteger operator-(const BigInteger& x, const BigInteger& y);
BigInteger operator+(const BigInteger& x, const BigInteger& y);

std::string BigInteger::toString() const {
    std::string ans = "";
    if (negative)
        ans += '-';
    for (int i = static_cast<int>(getLen()) - 1; i >= 0; --i) {
        std::string cur = std::to_string(digit[i]);
        if (i != static_cast<int>(getLen()) - 1) {
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
    return *this = *this + (-x);
}

BigInteger operator-(const BigInteger& x, const BigInteger& y) {
    BigInteger ans = x;
    ans -= y;
    return ans;
}

BigInteger& BigInteger::operator+=(const BigInteger& x) {
    if (this->abs() < x.abs()) {
        return *this = x + *this;
    }
    add(x, negative ^ x.negative);
    return *this;
}

BigInteger operator+(const BigInteger& x, const BigInteger& y) {
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

BigInteger BigInteger::operator++(int) {
    *this += 1;
    return (*this - 1);
}

BigInteger BigInteger::operator--(int) {
    *this -= 1;
    return (*this + 1);
}

BigInteger operator "" _bi(unsigned long long x) {
    return BigInteger(x);
}

BigInteger operator*(const BigInteger& x, const BigInteger& y);

void fft(std::vector<std::complex<long double>>& a, bool invert) {
    size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        long double ang = 2 * PI / len * (invert ? -1 : 1);
        std::complex<long double> wlen(cos(ang), sin(ang));
        for (size_t i = 0; i < n; i += len) {
            std::complex<long double> w(1);
            for (size_t j = 0; j < len / 2; j++) {
                std::complex<long double> u = a[i + j];
                std::complex<long double> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (invert) {
        for (std::complex<long double>& x : a) {
            x /= n;
        }
    }
}

BigInteger& BigInteger::operator*=(const BigInteger& y) {
    if (sign()) {
        if (y.sign()) {
            return *this = (-(*this)) * (-y);
        }
        return *this = -((-(*this)) * y);
    } else if (y.sign()) {
        return *this = -((*this) * (-y));
    }
    std::vector<std::complex<long double>> fx(digit.begin(), digit.end());
    std::vector<std::complex<long double>> fy(y.digit.begin(), y.digit.end());
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
        if (i >= getLen()) {
            digit.push_back(0);
        }
        long long cur = 0;
        cur = carry;
        if (i < fx.size()) {
            cur += static_cast<long long>(fx[i].real() + 0.5);
        }
        carry = 0;
        if (cur >= BASE) {
            carry = cur / BASE;
            cur %= BASE;
        }
        digit[i] = cur;
    }
    normalize();
    return *this;
}

BigInteger operator*(const BigInteger& x, const BigInteger& y) {
    BigInteger z = x;
    z *= y;
    return z;
}

BigInteger operator/(const BigInteger& x, const BigInteger& y) {
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
    std::vector<int> alls(getLen() * LEN);
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
        cur *= 10;
        *this *= 10;
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

BigInteger operator%(const BigInteger& x, const BigInteger& y) {
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
        BigInteger t;
        while (y) {
            x %= y;
            t = x;
            x = y;
            y = t;
        }
        return x;
    }
public:
    Rational operator-();
    std::string toString() const;
    void swap(Rational& x);
    Rational& operator=(const Rational& x);
    Rational& operator+=(const Rational& x);
    Rational& operator-=(const Rational& x);
    Rational& operator*=(const Rational& x);
    Rational& operator/=(const Rational& x);
    std::string asDecimal(size_t precision = 0) const;

    Rational() {
        num = 0;
        den = 1;
    }

    Rational(int x) {
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
        x = gcd((num.negative ? -num : num), den);
        num /= x;
        den /= x;
        if (num == 0) {
            den = 1;
            num.negative = 0;
        }
    }

    explicit operator double() const {
        double x = 0;
        std::string s = asDecimal(18);
        size_t i = 0;
        for (; i < s.size() && s[i] != '.'; ++i) {
            x = x * 10 + static_cast<int>(s[i] - '0');
        }
        double cur = 1. / 10.;
        ++i;
        for (; i < s.size(); ++i) {
            x += static_cast<double>(s[i] - '0') * cur;
            cur /= 10;
        }
        return x;
    }
};

Rational Rational::operator-() {
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
    return !(x < y);
}

bool operator>=(const Rational& y, const Rational& x) {
    return !(y < x);
}

std::string Rational::toString() const {
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

Rational operator*(const Rational& x, const Rational& y) {
    Rational z = x;
    z *= y;
    return z;
}

Rational operator/(const Rational& x, const Rational& y) {
    Rational z = x;
    z /= y;
    return z;
}

Rational operator+(const Rational& x, const Rational& y) {
    Rational z = x;
    z += y;
    return z;
}

Rational operator-(const Rational& x, const Rational& y) {
    Rational z = x;
    z -= y;
    return z;
}

std::string Rational::asDecimal(size_t precision) const {
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
    BigInteger x = num;
    for (size_t i = 0; i < precision; ++i) {
        x *= 10;
    }
    x /= den;
    std::string ans = "";
    if (x.negative) {
        ans += "-";
    }
    x.negative = 0;
    std::string s = x.toString();
    if (num < den) {
        ans += "0.";
        for (size_t i = 0; i < precision - s.size(); ++i) {
            ans += "0";
        }
        ans += s;
    } else {
        ans += s.substr(0, s.size() - precision) + "." + s.substr(s.size() - precision);
    }
    return ans;
}


