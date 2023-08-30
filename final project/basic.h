#pragma once
using db = double;
#include <variant>
#include "bits/stdc++.h"
using namespace std;
using sparse = vector<vector<pair<int, db>>>;
class point {
public:
    db x, y;
    db norm_sqr() const {
        return x * x + y * y;
    }
    point operator+(const point &o) const {
        return {x + o.x, y + o.y};
    }
};
point operator*(db k, const point &o) {
    return {k * o.x, k * o.y};
}
class point_3d {
public:
    db x, y, z;
    db norm_sqr() const {
        return x * x + y * y + z * z;
    }
    point_3d operator+(const point_3d &o) const {
        return {x + o.x, y + o.y, z + o.z};
    }
};
point_3d operator*(db k, const point_3d &o) {
    return {k * o.x, k * o.y, k * o.z};
}
class boundary_value {
    db phi;

protected:
    boundary_value(const db &_phi) : phi(_phi) {
    }

public:
    db value() const {
        return phi;
    }
};
class Dirichlet : public boundary_value {
public:
    Dirichlet(const db &_phi) : boundary_value(_phi) {
    }
};
class Neumann : public boundary_value {
public:
    Neumann(const db &_phi) : boundary_value(_phi) {
    }
};
