#include<cmath>
#include<vector>
#include<iostream>
#include<initializer_list>

const double PI = atan2(0, -1);
const double INF = 1e9;
const double EPS = 1e-6;

struct Vector2;
class Line;

bool equal(double x, double y) {
    return fabs(x - y) < EPS;
}

struct Point {
    double x;
    double y;

    Point(double x = 0, double y = 0) : x(x), y(y) {}
    explicit Point(const Vector2& a);

    Point& operator+=(const Vector2& a);

    double distance(const Point& a) const {
        return sqrt((x - a.x) * (x - a.x) + (y - a.y) * (y - a.y));
    }
    Point& rotate(Point centre, double angle);
    Point& reflex(const Point& a);
    Point& reflex(const Line& axis);
};

bool operator==(const Point& a, const Point& b) {
    return equal(a.x, b.x) && equal(a.y, b.y);
}

bool operator!=(const Point& a, const Point& b) {
    return !(a == b);
}

struct Vector2 {
    double x;
    double y;

    Vector2() = default;
    Vector2(const Point& a) : x(a.x), y(a.y) {}
    Vector2(double x, double y) : x(x), y(y) {}
    Vector2(const Point& a, const Point& b) : x(b.x - a.x), y(b.y - a.y) {}

    Vector2& operator+=(const Vector2& b);
    Vector2& operator-=(const Vector2& b);
    Vector2& operator*=(double b);
    Vector2 operator-() const;

    Vector2& normalize();
    Vector2& rotate(double angle);
    double length() const;
};

double scalar(const Vector2& a, const Vector2& b) {
    return a.x * b.x + a.y * b.y;
}

Point& Point::operator+=(const Vector2& a) {
    x += a.x;
    y += a.y;
    return *this;
}

Vector2 operator-(const Point& a, const Point& b) {
    return Vector2(b, a);
}

Point::Point(const Vector2& a) : x(a.x), y(a.y) {}

bool operator==(const Vector2& a, const Vector2& b) {
    return equal(a.x, b.x) && equal(a.y, b.y);
}

bool operator!=(const Vector2& a, const Vector2& b) {
    return !(a == b);
}

double Vector2::length() const {
    return sqrt(x * x + y * y);
}

Point& Point::reflex(const Point& a) {
    Vector2 delta = Vector2(*this, a);
    *this += delta;
    return *this;
}

Vector2& Vector2::normalize() {
    double len = length();
    x /= len;
    y /= len;
    return *this;
}

Vector2& Vector2::rotate(double angle) {
    double nx = x * cos(angle) - y * sin(angle);
    double ny = x * sin(angle) + y * cos(angle);
    x = nx;
    y = ny;
    return *this;
}

Vector2 Vector2::operator-() const {
    return Vector2(-x, -y);
}


Vector2& Vector2::operator+=(const Vector2& b) {
    x += b.x;
    y += b.y;
    return *this;
}
Vector2& Vector2::operator-=(const Vector2& b) {
    x -= b.x;
    y -= b.y;
    return *this;
}

Vector2& Vector2::operator*=(double b) {
    x *= b;
    y *= b;
    return *this;
}

Vector2 operator+(const Vector2& a, const Vector2& b) {
    Vector2 nw = a;
    nw += b;
    return nw;
}

Vector2 operator-(const Vector2& a, const Vector2& b) {
    Vector2 nw = a;
    nw -= b;
    return nw;
}

double operator*(const Vector2& a, const Vector2& b) {
    return a.x * b.y - a.y * b.x;
}

Vector2 operator*(double a, const Vector2& b) {
    Vector2 nw = b;
    nw *= a;
    return nw;
}

Vector2 operator*(const Vector2& a, double b) {
    Vector2 nw = a;
    nw *= b;
    return nw;
}

Point& Point::rotate(Point centre, double angle) {
    Vector2 v = Vector2(centre, *this);
    v.rotate(angle);
    centre += v;
    return *this = centre;
}

Vector2 Bisect(Vector2 a, Vector2 b) {
    if (scalar(a, b) < 0)
        b = -b;
    a.normalize();
    b.normalize();
    return a + b;
}

Point median(Point a, Point b) {
    return static_cast<Point>(static_cast<Vector2>(a) + (b - a) * 0.5);
}

class Line {
private:
    friend Point intersection(const Line& a, const Line& b);
    double A = 1;
    double B = 0;
    double C = 0;
    friend bool operator==(const Line& a, const Line& b);
public:
    void normalize() {
        if (!equal(A, 0)) {
            B /= A;
            C /= A;
            A = 1.0;
        } else {
            C /= B;
            B = 1.0;
        }
    }
    Line() = default;
    Line(const Point& a, const Point& b) {
        A = a.y - b.y;
        B = b.x - a.x;
        C = a.x * b.y - b.x * a.y;
        normalize();
    }
    Line(double k, double shift) : Line(Point(0, shift), Point(1, k + shift)) {}
    Line(const Point& a, double k) : Line(a, Point(0, a.y - a.x * k)) {}

    Vector2 normal() const {
        return Vector2(A, B);
    }

    double get(Vector2 x) const {
        return x.x * A + x.y * B + C;
    }

    double distance(Point x) const {
        return std::abs(A * x.x + B * x.y + C) / sqrt(A * A + B * B);
    }
};

Point& Point::reflex(const Line& axis) {
    Vector2 n = axis.normal();
    if (axis.get(*this) > 0)
        n *= -1;
    double dist = 2.0 * axis.distance(*this);
    n *= (dist / n.length());
    *this += n;
    return *this;
}

Point intersection(const Line& a, const Line& b) {
    return Point((b.C * a.B - a.C * b.B) / (a.A * b.B - b.A * a.B), (b.A * a.C - a.A * b.C) / (a.A * b.B - b.A * a.B));
}

bool operator==(const Line& a, const Line& b) {
    return equal(a.A, b.A) && equal(a.B, b.B) && equal(a.C, b.C);
}

bool operator!=(const Line& a, const Line& b) {
    return !(a == b);
}

class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator==(const Shape& another) const = 0;
    virtual bool operator!=(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual void rotate(Point center, double angle) = 0;
    virtual void reflex(Point center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;
    virtual ~Shape() = default;
};

class Polygon : public Shape {
protected:
    std::vector<Point> vertices;
public:
    const Point centroid() const {
        Point centr;
        for (auto d : vertices) {
            centr += d;
        }
        centr.x /= static_cast<double>(verticesCount());
        centr.y /= static_cast<double>(verticesCount());
        return centr;
    }
    Polygon() = default;
    Polygon(const std::vector<Point>& vertices) : vertices(vertices) {}
    Polygon(const std::initializer_list<Point>& vertices) : vertices(vertices) {}

    size_t verticesCount() const {
        return vertices.size();
    }

    const std::vector<Point>& getVertices() const {
        return vertices;
    }

    double perimeter() const override {
        double P = 0;
        for (size_t i = 0; i < vertices.size(); ++i) {
            P += Vector2(vertices[i], vertices[i + 1]).length();
        }
        return P;
    }

    double area() const override {
        double S = 0;
        for (size_t i = 1; i + 1 < vertices.size(); ++i) {
            S += Vector2(vertices[0], vertices[i]) * Vector2(vertices[0], vertices[i + 1]);
        }
        return std::abs(S) / 2.0;
    }

    bool isConvex() const {
        bool positive = false;
        bool negative = false;
        for (size_t i = 0; i < vertices.size(); ++i) {
            if (Vector2(vertices[i], vertices[(i + 1) % verticesCount()]) * Vector2(vertices[(i + 1) % verticesCount()], vertices[(i + 2) % verticesCount()]) > 0)
                positive = 1;
            else
                negative = 1;
        }
        return positive ^ negative;
    }

    bool operator==(const Shape& another) const override {
        auto ptr = dynamic_cast<const Polygon*>(&another);
        if (ptr == nullptr)
            return false;
        const auto& b = dynamic_cast<const Polygon&>(another);
        if (b.verticesCount() != verticesCount())
            return false;
        size_t ind = 0;
        for (; ind < verticesCount(); ++ind) {
            if (b.vertices[ind] == vertices[0])
                break;
        }
        if (ind == verticesCount())
            return false;

        int sign = 0;
        size_t n = verticesCount();
        if (b.vertices[(ind + 1) % n] == vertices[1])
            sign = 1;
        else
            sign = -1;
        for (size_t i = 0; i < verticesCount(); ++i) {
            if (b.vertices[ind] != vertices[i])
                return false;
            ind = (ind + n + sign) % n;
        }
        return true;
    }

    bool operator!=(const Shape& another) const override {
        return !(*this == another);
    }

    bool isCongruentTo(const Shape& another) const override {
        auto ptr = dynamic_cast<const Polygon*>(&another);
        if (ptr == nullptr)
            return false;
        const auto& b = dynamic_cast<const Polygon&>(another);
        if (b.verticesCount() != verticesCount())
            return false;
        Point first = centroid();
        Point second = b.centroid();
        Polygon firstshape;
        Polygon secondshape;
        for (size_t i = 0; i < verticesCount(); ++i) {
            firstshape.vertices.push_back(Point(std::abs((vertices[(i + 1) % verticesCount()] - vertices[i]) *
                                                         (vertices[(i - 1 + verticesCount()) % verticesCount()] - vertices[i])), (vertices[i] - first).length()));
            secondshape.vertices.push_back(Point(std::abs((b.vertices[(i + 1) % verticesCount()] - b.vertices[i]) *
                                                          (b.vertices[(i - 1 + verticesCount()) % verticesCount()] - b.vertices[i])), (b.vertices[i] - second).length()));
        }
        bool f = firstshape == secondshape;
        return f;
    }

    bool containsPoint(Point point) const override {
        double sum = 0;
        for (size_t i = 0; i < verticesCount(); ++i) {
            Vector2 a = Vector2(point, vertices[i]);
            Vector2 b = Vector2(point, vertices[(i + 1) % verticesCount()]);
            sum += atan2(a * b, scalar(a, b));
        }
        return !equal(sum, 0);
    }

    bool isSimilarTo(const Shape& another) const override {
        auto ptr = dynamic_cast<const Polygon*>(&another);
        if (ptr == nullptr)
            return false;
        const auto& b = dynamic_cast<const Polygon&>(another);
        if (b.verticesCount() != verticesCount())
            return false;
        double mnf = 1e70;
        double mns = 1e70;
        for (size_t i = 0; i < verticesCount(); ++i) {
            mnf = std::min((vertices[(i + 1) % verticesCount()] - vertices[i]).length(), mnf);
            mns = std::min((b.vertices[(i + 1) % verticesCount()] - b.vertices[i]).length(), mns);
        }
        Polygon nw = *this;
        nw.scale(Point(0, 0), mns / mnf);
        return nw.isCongruentTo(b);
    }

    void rotate(Point center, double angle) override {
        angle = PI * angle / 180.0;
        for (auto& d : vertices) {
            d = Point(center + (d - center).rotate(angle));
        }
    }

    void reflex(Point center) override {
        for (auto& d : vertices) {
            d = Point(center - Vector2(center, d));
        }
    }

    void reflex(Line axis) override {
        for (auto& d : vertices) {
            d.reflex(axis);
        }
    }

    void scale(Point center, double coefficient) override {
        for (auto& d : vertices) {
            d = static_cast<Point>(center + static_cast<Vector2>(d - center) * coefficient);
        }
    }
};

class Ellipse: public Shape {
protected:
    std::pair<Point, Point> focus;
    double sum;
public:
    Ellipse() = default;
    Ellipse(const Point& F1, const Point& F2, double sum) : focus(std::make_pair(F1, F2)), sum(sum) {}

    std::pair<Point, Point> focuses() const {
        return focus;
    }

    Point center() const {
        return median(focus.first, focus.second);
    }

    double eccentricity() const {
        return static_cast<double>(2) * std::abs((focus.second - center()).length() / sum);
    }

    std::pair<Line, Line> directrices() const {
        Vector2 perp = ((focus.first - center()).normalize());
        Vector2 n = perp.rotate(PI / 2.0);
        double k = sum / 2.0 * eccentricity();
        return std::make_pair(Line(Point(center() + perp * k), Point(center() + perp * k + n)), Line(Point(center() - perp * k), Point(center() - perp * k + n)));
    }

    double perimeter() const override {
        double a = sum / 2.0;
        double b = sqrt(a * a - eccentricity() * eccentricity() * a * a);
        return PI * (3.0 * (a + b) - sqrt((3 * a + b) * (a + 3.0 * b)));
    }

    double area() const override {
        double a = sum / 2.0;
        double b = sqrt(a * a - eccentricity() * eccentricity() * a * a);
        return PI * a * b;
    }

    bool operator==(const Shape& another) const override {
        if (dynamic_cast<const Ellipse*>(&another) == nullptr)
            return false;
        const auto& second = dynamic_cast<const Ellipse&>(another);
        return (equal(sum, second.sum)) && ((second.focus == focus) || (second.focus.second == focus.first && second.focus.first == focus.second));
    }

    bool operator!=(const Shape& another) const override {
        return !(*this == another);
    }

    bool isCongruentTo(const Shape& another) const override {
        if (dynamic_cast<const Ellipse*>(&another) == nullptr)
            return false;
        const auto& second = dynamic_cast<const Ellipse&>(another);
        return (equal(sum, second.sum) && equal(eccentricity(), second.eccentricity()));
    }

    bool isSimilarTo(const Shape& another) const override {
        if (dynamic_cast<const Ellipse*>(&another) == nullptr)
            return false;
        const auto& second = dynamic_cast<const Ellipse&>(another);
        return equal(eccentricity(), second.eccentricity());
    }

    bool containsPoint(Point point) const override {
        return sum > point.distance(focus.first) + point.distance(focus.second) - EPS;
    }
    void rotate(Point center, double angle) override {
        angle = PI * angle / 180.0;
        focus.first.rotate(center, angle);
        focus.second.rotate(center, angle);
    }
    void reflex(Point center) override {
        focus.first.reflex(center);
        focus.second.reflex(center);
    }
    virtual void reflex(Line axis) override {
        focus.first.reflex(axis);
        focus.second.reflex(axis);
    }
    virtual void scale(Point center, double coefficient) override {
        focus.first = Point(center + (focus.first - center) * coefficient);
        focus.second = Point(center + (focus.second - center) * coefficient);
        sum *= coefficient;
    }
};

class Circle: public Ellipse {
public:
    Circle() = default;
    Circle(const Point& a, double r) : Ellipse(a, a, 2.0 * r) {}

    double radius() const {
        return sum / 2.0;
    }
};

class Triangle: public Polygon {
public:
    Triangle() = default;
    Triangle(const Point& a, const Point& b, const Point& c) : Polygon({a, b, c}) {}

    Point orthocenter() const {
        Line a = Line(vertices[0], vertices[1]);
        Line b = Line(vertices[0], vertices[2]);
        Line n1 = Line(vertices[2], Point(vertices[2] + a.normal()));
        Line n2 = Line(vertices[1], Point(vertices[1] + b.normal()));
        return intersection(n1, n2);
    }

    Circle circumscribedCircle() const {
        Point centre = Point((centroid() - orthocenter()) * 1.5 + orthocenter());
        return Circle(centre, (vertices[0] - centre).length());
    }

    Circle inscribedCircle() const {
        Vector2 v01 = vertices[1] - vertices[0];
        Vector2 v02 = vertices[2] - vertices[0];
        Vector2 v12 = vertices[2] - vertices[1];
        Vector2 v10 = -v01;
        Line bis1 = Line(vertices[0], Point(vertices[0] + Bisect(v01, v02)));
        Line bis2 = Line(vertices[1], Point(vertices[1] + Bisect(v12, v10)));
        Point centre = intersection(bis1, bis2);
        Line a = Line(vertices[0], vertices[1]);
        double r = a.distance(centre);
        return Circle(centre, r);
    }

    Line EulerLine() const {
        return Line(orthocenter(), centroid());
    }

    Circle ninePointsCircle() const {
        Point a = median(vertices[0], vertices[1]);
        Point b = median(vertices[1], vertices[2]);
        Point c = median(vertices[0], vertices[2]);
        Triangle eul(a, b, c);
        return eul.circumscribedCircle();
    }
};


class Rectangle: public Polygon {
public:
    Rectangle() = default;
    Rectangle(const Point& a, const Point& c, double k) {
        if (k < 1.0)
            k = 1.0 / k;
        Point b = Point(a + static_cast<Vector2>(((c - a) * (1.0 / sqrt(1 + k * k)))).rotate(asin(k / sqrt(1.0 + k * k))));
        Point d = Point(a + static_cast<Vector2>(((c - a) * (1.0 / sqrt(1 + k * k)))).rotate(asin(-k / sqrt(1.0 + k * k))));
        vertices = {a, b, c, d};
    }

    Point center() const {
        return median(vertices[0], vertices[2]);
    }

    std::pair<Line, Line> diagonals() const {
        return std::make_pair(Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3]));
    }
};

class Square: public Rectangle {
public:
    Square() = default;
    Square(const Point& a, const Point& c) : Rectangle(a, c, 1.0) {}

    Circle circumscribedCircle() const {
        return Circle(center(), std::abs(vertices[0].x - vertices[1].x) / sqrt(2.0));
    }

    Circle inscribedCircle() {
        return Circle(center(), std::abs(vertices[0].x - vertices[1].x) / 2.0);
    }
};

