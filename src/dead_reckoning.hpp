#ifndef AZIMUTH_DEAD_RECKONING_HPP
#define AZIMUTH_DEAD_RECKONING_HPP

#include "types.hpp"
#include <cmath>
#include <algorithm>

// Coordinated Flight Arc & Aerodynamic Dead Reckoning Extrapolator
class DeadReckoningExtrapolator {
public:
    DeadReckoningExtrapolator()
        : initialized_(false),
          failure_start_ts_(0.0),
          last_ts_(0.0),
          lambda_v_(0.08),  // Aerodynamic drag decay
          lambda_w_(0.12)   // Rotational angular damping
    {
    }

    void start(double failure_start_ts, const Pose6D& last_pose, const Vector3d& clean_v, const Vector3d& clean_w) {
        failure_start_ts_ = failure_start_ts;
        last_ts_ = failure_start_ts;
        current_pose_ = last_pose;
        initial_w_ = clean_w;
        current_v_world_ = clean_v;

        // Transform initial velocity into body coordinates: v_body = q^-1 * v_world
        v_body_0_ = last_pose.quat.conjugate().rotate(clean_v);

        initialized_ = true;
    }

    bool is_initialized() const { return initialized_; }

    Pose6D step(double current_ts) {
        if (!initialized_) return current_pose_;

        double dt = current_ts - last_ts_;
        if (dt <= 0.0) {
            current_pose_.timestamp = current_ts;
            return current_pose_;
        }

        double dt_fail = current_ts - failure_start_ts_;
        if (dt_fail < 0.0) dt_fail = 0.0;

        // 1. Damped angular velocity: w(t) = w0 * exp(-lambda_w * dt_fail)
        Vector3d w = initial_w_ * std::exp(-lambda_w_ * dt_fail);
        double w_norm = w.norm();

        if (w_norm > 1e-8) {
            double angle = w_norm * dt;
            double half_angle = 0.5 * angle;
            double s = std::sin(half_angle) / w_norm;
            Quaterniond delta_q = {std::cos(half_angle), w.x * s, w.y * s, w.z * s};
            current_pose_.quat = (current_pose_.quat * delta_q).normalized();
        }

        // 2. Coordinated Turn: Project velocity along the evolving aircraft body heading
        double speed_decay = std::exp(-lambda_v_ * dt_fail);
        Vector3d v_body_t = v_body_0_ * speed_decay;
        Vector3d next_v_world = current_pose_.quat.rotate(v_body_t);

        // 3. Trapezoidal position integration along the flight curve
        current_pose_.pos = current_pose_.pos + (current_v_world_ + next_v_world) * (0.5 * dt);
        current_v_world_ = next_v_world;

        current_pose_.timestamp = current_ts;
        last_ts_ = current_ts;
        return current_pose_;
    }

private:
    bool initialized_;
    double failure_start_ts_;
    double last_ts_;
    Pose6D current_pose_;
    Vector3d initial_w_;
    Vector3d v_body_0_;
    Vector3d current_v_world_;
    double lambda_v_;
    double lambda_w_;
};

#endif // AZIMUTH_DEAD_RECKONING_HPP
