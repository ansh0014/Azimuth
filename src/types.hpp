#ifndef AZIMUTH_TYPES_HPP
#define AZIMUTH_TYPES_HPP

#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

constexpr double GRAVITY_STANDARD = 9.80665;

struct Vector3d {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vector3d operator+(const Vector3d& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vector3d operator-(const Vector3d& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vector3d operator*(double s) const { return {x * s, y * s, z * s}; }
    Vector3d operator/(double s) const { return {x / s, y / s, z / s}; }

    double dot(const Vector3d& o) const { return x * o.x + y * o.y + z * o.z; }
    Vector3d cross(const Vector3d& o) const {
        return {
            y * o.z - z * o.y,
            z * o.x - x * o.z,
            x * o.y - y * o.x
        };
    }
    double norm_sq() const { return x * x + y * y + z * z; }
    double norm() const { return std::sqrt(norm_sq()); }

    Vector3d normalized() const {
        double n = norm();
        if (n < 1e-12) return {0.0, 0.0, 0.0};
        return *this * (1.0 / n);
    }
};

struct Quaterniond {
    double w = 1.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Quaterniond operator*(const Quaterniond& b) const {
        return {
            w * b.w - x * b.x - y * b.y - z * b.z,
            w * b.x + x * b.w + y * b.z - z * b.y,
            w * b.y - x * b.z + y * b.w + z * b.x,
            w * b.z + x * b.y - y * b.x + z * b.w
        };
    }

    Quaterniond conjugate() const {
        return {w, -x, -y, -z};
    }

    Vector3d rotate(const Vector3d& v) const {
        // v' = q * [0, v] * q^-1
        Quaterniond qv = {0.0, v.x, v.y, v.z};
        Quaterniond res = (*this) * qv * conjugate();
        return {res.x, res.y, res.z};
    }

    Quaterniond normalized() const {
        double mag = std::sqrt(w * w + x * x + y * y + z * z);
        if (mag < 1e-12) return {1.0, 0.0, 0.0, 0.0};
        double inv = 1.0 / mag;
        return {w * inv, x * inv, y * inv, z * inv};
    }

    static Quaterniond from_axis_angle(const Vector3d& axis, double angle) {
        double half = 0.5 * angle;
        double s = std::sin(half);
        Vector3d a_norm = axis.normalized();
        return {std::cos(half), a_norm.x * s, a_norm.y * s, a_norm.z * s};
    }
};

struct Pose6D {
    double timestamp = 0.0;
    Vector3d pos = {0.0, 0.0, 0.0};
    Quaterniond quat = {1.0, 0.0, 0.0, 0.0}; // w, x, y, z
};

struct RawImuSample {
    double raw_timestamp = 0.0;
    double timestamp_s = 0.0;
    double ax = 0.0, ay = 0.0, az = 0.0; // m/s^2
    double gx = 0.0, gy = 0.0, gz = 0.0; // rad/s
    bool is_duplicate = false;
};

struct CleanedImuSample {
    double timestamp_s = 0.0;
    double dt = 0.0;
    double ax = 0.0, ay = 0.0, az = 0.0; // m/s^2
    double gx = 0.0, gy = 0.0, gz = 0.0; // rad/s
    bool is_failure = false;
};

#endif // AZIMUTH_TYPES_HPP
