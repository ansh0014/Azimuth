#ifndef AZIMUTH_ESTIMATOR_HPP
#define AZIMUTH_ESTIMATOR_HPP

#include "types.hpp"
#include "dead_reckoning.hpp"
#include <cmath>
#include <algorithm>

// Pure Native C++17 Adaptive Mahony/Madgwick AHRS & Drift Controller
class NativeMadgwickAhrs {
public:
    NativeMadgwickAhrs()
        : beta_(0.12),
          sample_period_(0.005),
          q_{1.0, 0.0, 0.0, 0.0},
          gyro_bias_{0.0, 0.0, 0.0},
          bias_filter_coeff_(0.01),
          bias_threshold_(0.06),
          initialized_(false)
    {
    }

    void init(double sample_period) {
        sample_period_ = (sample_period > 0.0) ? sample_period : 0.005;
        reset();
    }

    void reset() {
        q_ = {1.0, 0.0, 0.0, 0.0};
        gyro_bias_ = {0.0, 0.0, 0.0};
        initialized_ = false;
    }

    void set_sample_period(double dt) {
        if (dt > 1e-5) sample_period_ = dt;
    }

    void set_gyro_bias(const Vector3d& bias) {
        gyro_bias_ = bias;
    }

    Quaterniond get_quaternion() const { return q_; }
    Vector3d get_gyro_bias() const { return gyro_bias_; }

    void update(const Vector3d& raw_gyro, const Vector3d& raw_accel) {
        const double dt = sample_period_;
        const double accel_norm = raw_accel.norm();
        const double gyro_norm = raw_gyro.norm();

        // 1. Maneuver Acceleration Rejection Weight (w_a)
        double accel_diff = std::fabs(accel_norm - GRAVITY_STANDARD);
        double w_a = 1.0 / (1.0 + 4.0 * accel_diff * accel_diff);

        Vector3d a_dir = {0.0, 0.0, 1.0};
        if (accel_norm > 1.0) {
            a_dir = raw_accel * (1.0 / accel_norm);
        }

        // 2. Stationary / Quasi-static Gyro Bias Estimation
        if (gyro_norm < bias_threshold_ && accel_diff < 0.20) {
            gyro_bias_.x += bias_filter_coeff_ * (raw_gyro.x - gyro_bias_.x);
            gyro_bias_.y += bias_filter_coeff_ * (raw_gyro.y - gyro_bias_.y);
            gyro_bias_.z += bias_filter_coeff_ * (raw_gyro.z - gyro_bias_.z);
        }

        Vector3d gyro_comp = raw_gyro - gyro_bias_;

        // 3. Initialize orientation from accelerometer gravity vector
        if (!initialized_ && accel_norm > 1.0) {
            double pitch = std::atan2(-a_dir.x, std::sqrt(a_dir.y * a_dir.y + a_dir.z * a_dir.z));
            double roll = std::atan2(a_dir.y, a_dir.z);

            double cr = std::cos(roll * 0.5);
            double sr = std::sin(roll * 0.5);
            double cp = std::cos(pitch * 0.5);
            double sp = std::sin(pitch * 0.5);

            q_.w = cr * cp;
            q_.x = sr * cp;
            q_.y = cr * sp;
            q_.z = -sr * sp;
            q_ = q_.normalized();
            initialized_ = true;
        }

        // 4. Accelerometer Gradient Descent Correction with Adaptive Gain
        Quaterniond grad = {0.0, 0.0, 0.0, 0.0};
        double beta_eff = beta_ * w_a;

        if (accel_norm > 1.0 && beta_eff > 1e-4) {
            double qw = q_.w, qx = q_.x, qy = q_.y, qz = q_.z;

            double f1 = 2.0 * (qx * qz - qw * qy) - a_dir.x;
            double f2 = 2.0 * (qw * qx + qy * qz) - a_dir.y;
            double f3 = 2.0 * (0.5 - qx * qx - qy * qy) - a_dir.z;

            grad.w = -2.0 * qy * f1 + 2.0 * qx * f2;
            grad.x =  2.0 * qz * f1 + 2.0 * qw * f2 - 4.0 * qx * f3;
            grad.y = -2.0 * qw * f1 + 2.0 * qz * f2 - 4.0 * qy * f3;
            grad.z =  2.0 * qx * f1 + 2.0 * qy * f2;

            double g_norm = std::sqrt(grad.w * grad.w + grad.x * grad.x + grad.y * grad.y + grad.z * grad.z);
            if (g_norm > 1e-9) {
                double inv_g = 1.0 / g_norm;
                grad.w *= inv_g; grad.x *= inv_g; grad.y *= inv_g; grad.z *= inv_g;
            }
        }

        // 5. Gyroscope quaternion derivative: 0.5 * q * [0, wx, wy, wz]
        Quaterniond q_omega = {
            0.5 * (-q_.x * gyro_comp.x - q_.y * gyro_comp.y - q_.z * gyro_comp.z),
            0.5 * ( q_.w * gyro_comp.x + q_.y * gyro_comp.z - q_.z * gyro_comp.y),
            0.5 * ( q_.w * gyro_comp.y - q_.x * gyro_comp.z + q_.z * gyro_comp.x),
            0.5 * ( q_.w * gyro_comp.z + q_.x * gyro_comp.y - q_.y * gyro_comp.x)
        };

        // 6. Integrate: q_dot = q_omega - beta_eff * grad
        q_.w += (q_omega.w - beta_eff * grad.w) * dt;
        q_.x += (q_omega.x - beta_eff * grad.x) * dt;
        q_.y += (q_omega.y - beta_eff * grad.y) * dt;
        q_.z += (q_omega.z - beta_eff * grad.z) * dt;
        q_ = q_.normalized();
    }

    Vector3d extract_world_acceleration(const Vector3d& raw_accel) const {
        // q_ maps body vector to world: a_world = q * a_body * q^-1
        Vector3d a_world = q_.rotate(raw_accel);
        a_world.z -= GRAVITY_STANDARD;
        return a_world;
    }

private:
    double beta_;
    double sample_period_;
    Quaterniond q_;
    Vector3d gyro_bias_;
    double bias_filter_coeff_;
    double bias_threshold_;
    bool initialized_;
};

class StateEstimator {
public:
    static constexpr size_t HISTORY_CAPACITY = 256;

    StateEstimator()
        : start_timestamp_(-1.0),
          protect_seconds_(2.0),
          protect_sample_count_(0),
          in_protect_phase_(true),
          dr_active_(false),
          tau_bias_(4.0),
          vel_damping_(0.6),
          history_head_(0),
          history_count_(0)
    {
        current_pose_.timestamp = 0.0;
        current_pose_.pos = {0.0, 0.0, 0.0};
        current_pose_.quat = {1.0, 0.0, 0.0, 0.0};

        current_velocity_ = {0.0, 0.0, 0.0};
        last_world_accel_ = {0.0, 0.0, 0.0};

        sum_init_world_accel_ = {0.0, 0.0, 0.0};
        sum_init_gyro_ = {0.0, 0.0, 0.0};
        world_accel_bias_ = {0.0, 0.0, 0.0};
    }

    void init(double first_timestamp, double median_dt) {
        start_timestamp_ = first_timestamp;
        ahrs_.init(median_dt);
    }

    void set_parameters(double tau_bias, double vel_damping) {
        tau_bias_ = tau_bias;
        vel_damping_ = vel_damping;
    }

    Pose6D update(const CleanedImuSample& sample) {
        if (start_timestamp_ < 0.0) {
            start_timestamp_ = sample.timestamp_s;
        }

        if (sample.is_failure) {
            if (!dr_active_) {
                dr_active_ = true;
                Vector3d clean_v = compute_clean_velocity();
                Vector3d clean_w = compute_clean_angular_velocity();
                dr_extrapolator_.start(sample.timestamp_s, current_pose_, clean_v, clean_w);
            }
            current_pose_ = dr_extrapolator_.step(sample.timestamp_s);
            return current_pose_;
        }

        double elapsed = sample.timestamp_s - start_timestamp_;
        double dt = (sample.dt > 1e-4) ? sample.dt : 0.005;

        Vector3d gyro = {sample.gx, sample.gy, sample.gz};
        Vector3d accel = {sample.ax, sample.ay, sample.az};

        // AHRS attitude estimation
        ahrs_.set_sample_period(dt);
        ahrs_.update(gyro, accel);
        current_pose_.quat = ahrs_.get_quaternion();

        // Gravity-removed Earth acceleration (raw)
        Vector3d raw_world_accel = ahrs_.extract_world_acceleration(accel);
        Vector3d world_accel = {0.0, 0.0, 0.0};

        // 1. Protect window calibration (first 2.0 seconds)
        if (in_protect_phase_) {
            if (elapsed < protect_seconds_) {
                sum_init_world_accel_ = sum_init_world_accel_ + raw_world_accel;
                sum_init_gyro_ = sum_init_gyro_ + gyro;
                protect_sample_count_++;
                current_velocity_ = {0.0, 0.0, 0.0};
                world_accel = {0.0, 0.0, 0.0};
            } else {
                in_protect_phase_ = false;
                if (protect_sample_count_ > 0) {
                    double inv_n = 1.0 / static_cast<double>(protect_sample_count_);
                    Vector3d init_gyro_bias = sum_init_gyro_ * inv_n;
                    ahrs_.set_gyro_bias(init_gyro_bias);
                    world_accel_bias_ = sum_init_world_accel_ * inv_n;
                }
            }
        }

        // 2. Post-protect active flight tracking: online leaky world_accel_bias_ high-pass filter, ZUPT, damping
        if (!in_protect_phase_) {
            double alpha = dt / (tau_bias_ + dt);
            world_accel_bias_ = world_accel_bias_ + (raw_world_accel - world_accel_bias_) * alpha;
            world_accel = raw_world_accel - world_accel_bias_;

            // Safety clamp: prevent any residual spike > 15.0 m/s^2 from blowing up integration
            double wa_norm = world_accel.norm();
            if (wa_norm > 15.0) {
                world_accel = world_accel * (15.0 / wa_norm);
            }

            // Exponential velocity damping
            current_velocity_ = current_velocity_ * std::exp(-vel_damping_ * dt);

            // Zero-Velocity Update (ZUPT) during stationary/quasi-stationary segments
            if (gyro.norm() < 0.08 && std::fabs(accel.norm() - GRAVITY_STANDARD) < 0.25) {
                current_velocity_ = current_velocity_ * std::exp(-8.0 * dt);
            }
        }

        // Trapezoidal integration
        Vector3d next_velocity = current_velocity_ + (last_world_accel_ + world_accel) * (0.5 * dt);
        current_pose_.pos = current_pose_.pos + (current_velocity_ + next_velocity) * (0.5 * dt);
        current_velocity_ = next_velocity;
        last_world_accel_ = world_accel;
        current_pose_.timestamp = sample.timestamp_s;

        // Rolling history for dead reckoning handoff
        push_history(sample.timestamp_s, current_velocity_, gyro);

        return current_pose_;
    }

private:
    NativeMadgwickAhrs ahrs_;
    DeadReckoningExtrapolator dr_extrapolator_;

    double start_timestamp_;
    double protect_seconds_;
    size_t protect_sample_count_;
    bool in_protect_phase_;
    bool dr_active_;
    double tau_bias_;
    double vel_damping_;

    Pose6D current_pose_;
    Vector3d current_velocity_;
    Vector3d last_world_accel_;

    Vector3d sum_init_world_accel_;
    Vector3d sum_init_gyro_;
    Vector3d world_accel_bias_;

    struct HistoryEntry {
        double timestamp;
        Vector3d velocity;
        Vector3d angular_velocity;
    };
    HistoryEntry history_[HISTORY_CAPACITY];
    size_t history_head_;
    size_t history_count_;

    void push_history(double t, const Vector3d& v, const Vector3d& w) {
        size_t idx = (history_head_ + history_count_) % HISTORY_CAPACITY;
        if (history_count_ == HISTORY_CAPACITY) {
            history_head_ = (history_head_ + 1) % HISTORY_CAPACITY;
        } else {
            history_count_++;
        }
        history_[idx] = {t, v, w};
    }

    Vector3d compute_clean_velocity() const {
        if (history_count_ == 0) return {0.0, 0.0, 0.0};
        size_t n = std::min(history_count_, static_cast<size_t>(150));
        Vector3d sum = {0.0, 0.0, 0.0};
        for (size_t i = 0; i < n; ++i) {
            size_t idx = (history_head_ + history_count_ - 1 - i) % HISTORY_CAPACITY;
            sum = sum + history_[idx].velocity;
        }
        return sum * (1.0 / static_cast<double>(n));
    }

    Vector3d compute_clean_angular_velocity() const {
        if (history_count_ == 0) return {0.0, 0.0, 0.0};
        size_t n = std::min(history_count_, static_cast<size_t>(50));
        Vector3d sum = {0.0, 0.0, 0.0};
        for (size_t i = 0; i < n; ++i) {
            size_t idx = (history_head_ + history_count_ - 1 - i) % HISTORY_CAPACITY;
            sum = sum + history_[idx].angular_velocity;
        }
        return sum * (1.0 / static_cast<double>(n));
    }
};

#endif // AZIMUTH_ESTIMATOR_HPP
