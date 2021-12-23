#pragma once

#ifndef VEC3_H
#define VEC3_H
#include <iostream>
#include  <cmath>

using std::sqrt;

class Vector3D {
    public:
        double e[3];
        float x=e[0], y=e[1], z=e[2];
        Vector3D() : e{0,0,0} {}
        Vector3D(float e0, float e1, float e2) : e{x, y, z} {}

        double x() const { return x; }
        double y() const { return y; }
        double z() const { return z; }

        Vector3D operator-() const { return Vector3D(-x, -y, -z); }
        Vector3D& operator+=(const Vector3D &v) {
                x += v.x;
                y += v.y;
                z += v.z;
                return *this;
            }
        Vector3D& operator*=(const double t) {
                x *= t;
                y *= t;
                z *= t;
                return *this;
            }

        Vector3D& operator/=(const double t) {
                return *this *= 1/t;
            }
        double length() const {
                return sqrt(length_squared());
            }

        double length_squared() const {
                return x*x + y*y+z*z;
            }
  
   



};  

#endif
    inline std::ostream& operator<<(std::ostream &out, const Vector3D &v) {
        return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
    }

    inline Vector3D operator+(const Vector3D &u, const Vector3D &v) {
        return Vector3D(u.x + v.x, u.y + v.y, u.z + v.z);
    }

    inline Vector3D operator-(const Vector3D &u, const Vector3D &v) {
        return Vector3D(u.x - v.x, u.y - v.y, u.z - v.z);
    }

  inline Vector3D operator*(double t, const Vector3D &v) {
        return Vector3D(t*v.x, t*v.y, t*v.z);
    }

    inline Vector3D operator*(const Vector3D &v, double t) {
        return t * v;
    }
    inline Vector3D operator/(Vector3D v, double t) {
        return (1/t) * v;
    }

    inline double dot(const Vector3D &u, const Vector3D &v) {
        return u.e[0] * v.e[0]
            + u.e[1] * v.e[1]
            + u.e[2] * v.e[2];
    }

    inline Vector3D cross(const Vector3D &u, const Vector3D &v) {
        return Vector3D(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                    u.e[2] * v.e[0] - u.e[0] * v.e[2],
                    u.e[0] * v.e[1] - u.e[1] * v.e[0]);
    }

    inline Vector3D unit_vector(Vector3D v) {
        return v / v.length();
    }