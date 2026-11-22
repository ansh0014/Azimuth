#ifndef AZIMUTH_DSP_CLEANER_HPP
#define AZIMUTH_DSP_CLEANER_HPP

#include "types.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>

class BiquadFilter {
public:
    BiquadFilter() : b0(1.0), b1(0.0), b2(0.0), a1(0.0), a2(0.0), x1(0.0), x2(0.0), y1(0.0), y2(0.0), initialized(false) {}

    void init(double cutoff_hz, double sample_rate_hz) {
        if (sample_rate_hz <= 0.0) sample_rate_hz = 200.0;
        if (cutoff_hz >= sample_rate_hz * 0.49) {
            cutoff_hz = sample_rate_hz * 0.49;
        }

        const double omega = 2.0 * M_PI * cutoff_hz;
        const double T = 1.0 / sample_rate_hz;
        const double K = std::tan(omega * T * 0.5);
        const double K2 = K * K;
        const double sqrt2 = 1.4142135623730950488;
        const double norm = 1.0 / (1.0 + sqrt2 * K + K2);

        b0 = K2 * norm;
        b1 = 2.0 * b0;
        b2 = b0;
        a1 = 2.0 * (K2 - 1.0) * norm;
        a2 = (1.0 - sqrt2 * K + K2) * norm;

        reset();
    }

    void reset() {
        x1 = x2 = y1 = y2 = 0.0;
        initialized = false;
    }

    double filter(double x) {
        if (!initialized) {
            x1 = x2 = x;
            y1 = y2 = x;
            initialized = true;
            return x;
        }
        double y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1;
        x1 = x;
        y2 = y1;
        y1 = y;
        return y;
    }

private:
    double b0, b1, b2;
    double a1, a2;
    double x1, x2;
    double y1, y2;
    bool initialized;
};

class DspCleaner {
public:
    static constexpr size_t HAMPEL_WINDOW = 15;
    static constexpr size_t RING_CAPACITY = 64;

    DspCleaner()
        : median_dt_(0.005),
          consecutive_zeros_(0),
          failure_latched_(false),
          last_ts_(-1.0),
          history_count_(0),
          stuck_run_len_(0)
    {
    }

    void init(double median_dt) {
        median_dt_ = (median_dt > 0.0) ? median_dt : 0.005;
        double fs = 1.0 / median_dt_;

        for (int i = 0; i < 6; ++i) {
            biquad_[i].init(18.0, fs);
        }

        consecutive_zeros_ = 0;
        failure_latched_ = false;
        last_ts_ = -1.0;
        history_count_ = 0;
        stuck_run_len_ = 0;
    }

    bool is_failure_latched() const {
        return failure_latched_;
    }

    // Process sample immediately (1-in, 1-out streaming, zero frame lag)
    CleanedImuSample clean_sample(const RawImuSample& sample) {
        CleanedImuSample out;
        out.timestamp_s = sample.timestamp_s;
        out.dt = (last_ts_ > 0.0) ? (sample.timestamp_s - last_ts_) : median_dt_;
        if (out.dt <= 0.0) out.dt = median_dt_;
        last_ts_ = sample.timestamp_s;

        // Check for 5 consecutive zeros failure tail
        double sum6 = std::fabs(sample.ax) + std::fabs(sample.ay) + std::fabs(sample.az) +
                      std::fabs(sample.gx) + std::fabs(sample.gy) + std::fabs(sample.gz);

        if (sum6 < 1e-9) {
            consecutive_zeros_++;
            if (consecutive_zeros_ >= 5) {
                failure_latched_ = true;
            }
        } else {
            consecutive_zeros_ = 0;
        }

        if (failure_latched_) {
            out.ax = 0.0; out.ay = 0.0; out.az = 0.0;
            out.gx = 0.0; out.gy = 0.0; out.gz = 0.0;
            out.is_failure = true;
            return out;
        }

        out.is_failure = false;
        double vals[6] = {sample.ax, sample.ay, sample.az, sample.gx, sample.gy, sample.gz};

        // Stuck-frame detection & interpolation
        if (history_count_ > 0) {
            const auto& prev = history_[history_count_ - 1];
            double diff = std::fabs(vals[0] - prev[0]) + std::fabs(vals[1] - prev[1]) +
                          std::fabs(vals[2] - prev[2]) + std::fabs(vals[3] - prev[3]) +
                          std::fabs(vals[4] - prev[4]) + std::fabs(vals[5] - prev[5]);

            if (diff < 1e-6) {
                stuck_run_len_++;
                if (stuck_run_len_ >= 2 && history_count_ >= 2) {
                    // Slight variation to prevent integrator freezing during stuck run
                    for (int j = 0; j < 6; ++j) {
                        vals[j] = prev[j] + (history_[history_count_ - 1][j] - history_[history_count_ - 2][j]) * 0.5;
                    }
                }
            } else {
                stuck_run_len_ = 0;
            }
        }

        // Add to Hampel history window
        if (history_count_ < HAMPEL_WINDOW) {
            for (int j = 0; j < 6; ++j) history_[history_count_][j] = vals[j];
            history_count_++;
        } else {
            for (size_t k = 0; k < HAMPEL_WINDOW - 1; ++k) {
                for (int j = 0; j < 6; ++j) history_[k][j] = history_[k + 1][j];
            }
            for (int j = 0; j < 6; ++j) history_[HAMPEL_WINDOW - 1][j] = vals[j];
        }

        // 15-sample sliding window Hampel spike filter
        if (history_count_ >= 5) {
            double window_vals[15];
            double devs[15];

            for (int axis = 0; axis < 6; ++axis) {
                for (size_t k = 0; k < history_count_; ++k) {
                    window_vals[k] = history_[k][axis];
                }
                size_t mid_idx = history_count_ / 2;
                std::nth_element(window_vals, window_vals + mid_idx, window_vals + history_count_);
                double median = window_vals[mid_idx];

                for (size_t k = 0; k < history_count_; ++k) {
                    devs[k] = std::fabs(window_vals[k] - median);
                }
                std::nth_element(devs, devs + mid_idx, devs + history_count_);
                double mad = devs[mid_idx];

                double threshold = 4.0 * 1.4826 * mad;
                if (std::fabs(vals[axis] - median) > threshold && threshold > 1e-6) {
                    vals[axis] = median;
                }
            }
        }

        // 2nd-order Butterworth low-pass filter
        out.ax = biquad_[0].filter(vals[0]);
        out.ay = biquad_[1].filter(vals[1]);
        out.az = biquad_[2].filter(vals[2]);
        out.gx = biquad_[3].filter(vals[3]);
        out.gy = biquad_[4].filter(vals[4]);
        out.gz = biquad_[5].filter(vals[5]);

        return out;
    }

private:
    double median_dt_;
    int consecutive_zeros_;
    bool failure_latched_;
    double last_ts_;
    size_t history_count_;
    int stuck_run_len_;
    double history_[HAMPEL_WINDOW][6];
    BiquadFilter biquad_[6];
};

#endif // AZIMUTH_DSP_CLEANER_HPP
