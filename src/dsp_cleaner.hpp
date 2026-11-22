#ifndef AZIMUTH_DSP_CLEANER_HPP
#define AZIMUTH_DSP_CLEANER_HPP

#include "types.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <array>

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
    static constexpr size_t RING_CAPACITY = 64;
    static constexpr size_t HAMPEL_WINDOW = 15;
    static constexpr size_t HALF_WINDOW = HAMPEL_WINDOW / 2; // 7

    DspCleaner()
        : head_(0),
          count_(0),
          median_dt_(0.005),
          consecutive_zeros_(0),
          failure_latched_(false),
          last_pushed_ts_(-1.0),
          last_emitted_ts_(-1.0)
    {
    }

    void init(double median_dt) {
        median_dt_ = (median_dt > 0.0) ? median_dt : 0.005;
        double fs = 1.0 / median_dt_;

        for (int i = 0; i < 6; ++i) {
            biquad_[i].init(18.0, fs);
        }

        head_ = 0;
        count_ = 0;
        consecutive_zeros_ = 0;
        failure_latched_ = false;
        last_pushed_ts_ = -1.0;
        last_emitted_ts_ = -1.0;
    }

    bool is_failure_latched() const {
        return failure_latched_;
    }

    bool push_sample(const RawImuSample& sample) {
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

        size_t push_idx = (head_ + count_) % RING_CAPACITY;
        buffer_[push_idx] = sample;
        count_++;
        last_pushed_ts_ = sample.timestamp_s;

        return (count_ > HALF_WINDOW);
    }

    bool has_pending() const {
        return count_ > 0;
    }

    bool pop_cleaned(CleanedImuSample& out) {
        if (count_ == 0) return false;

        if (failure_latched_) {
            const RawImuSample& raw = buffer_[head_];
            out.timestamp_s = raw.timestamp_s;
            out.dt = (last_emitted_ts_ > 0.0) ? (raw.timestamp_s - last_emitted_ts_) : median_dt_;
            if (out.dt <= 0.0) out.dt = median_dt_;
            last_emitted_ts_ = raw.timestamp_s;

            out.ax = 0.0; out.ay = 0.0; out.az = 0.0;
            out.gx = 0.0; out.gy = 0.0; out.gz = 0.0;
            out.is_failure = true;

            head_ = (head_ + 1) % RING_CAPACITY;
            count_--;
            return true;
        }

        interpolate_stuck_and_dropouts();

        size_t target_offset = (count_ > HALF_WINDOW) ? HALF_WINDOW : 0;
        size_t target_idx = (head_ + target_offset) % RING_CAPACITY;

        RawImuSample center = buffer_[target_idx];

        double cleaned_vals[6] = {center.ax, center.ay, center.az, center.gx, center.gy, center.gz};
        apply_hampel(cleaned_vals);

        out.ax = biquad_[0].filter(cleaned_vals[0]);
        out.ay = biquad_[1].filter(cleaned_vals[1]);
        out.az = biquad_[2].filter(cleaned_vals[2]);
        out.gx = biquad_[3].filter(cleaned_vals[3]);
        out.gy = biquad_[4].filter(cleaned_vals[4]);
        out.gz = biquad_[5].filter(cleaned_vals[5]);

        out.timestamp_s = center.timestamp_s;
        out.dt = (last_emitted_ts_ > 0.0) ? (center.timestamp_s - last_emitted_ts_) : median_dt_;
        if (out.dt <= 0.0) out.dt = median_dt_;
        last_emitted_ts_ = center.timestamp_s;
        out.is_failure = false;

        head_ = (head_ + 1) % RING_CAPACITY;
        count_--;
        return true;
    }

private:
    std::array<RawImuSample, RING_CAPACITY> buffer_;
    size_t head_;
    size_t count_;
    double median_dt_;
    int consecutive_zeros_;
    bool failure_latched_;
    double last_pushed_ts_;
    double last_emitted_ts_;
    BiquadFilter biquad_[6];

    void interpolate_stuck_and_dropouts() {
        if (count_ < 3) return;

        size_t i = 0;
        while (i + 1 < count_) {
            size_t idx0 = (head_ + i) % RING_CAPACITY;
            size_t run_len = 1;

            while (i + run_len < count_) {
                size_t idx1 = (head_ + i + run_len) % RING_CAPACITY;
                double diff = std::fabs(buffer_[idx0].ax - buffer_[idx1].ax) +
                              std::fabs(buffer_[idx0].ay - buffer_[idx1].ay) +
                              std::fabs(buffer_[idx0].az - buffer_[idx1].az) +
                              std::fabs(buffer_[idx0].gx - buffer_[idx1].gx) +
                              std::fabs(buffer_[idx0].gy - buffer_[idx1].gy) +
                              std::fabs(buffer_[idx0].gz - buffer_[idx1].gz);

                if (diff < 1e-6) {
                    run_len++;
                } else {
                    break;
                }
            }

            if (run_len >= 2 && (i + run_len < count_)) {
                size_t start_idx = (head_ + i) % RING_CAPACITY;
                size_t end_idx = (head_ + i + run_len) % RING_CAPACITY;

                double t0 = buffer_[start_idx].timestamp_s;
                double t1 = buffer_[end_idx].timestamp_s;
                double dt = t1 - t0;

                if (dt > 1e-9) {
                    for (size_t r = 1; r < run_len; ++r) {
                        size_t mid_idx = (head_ + i + r) % RING_CAPACITY;
                        double s = (buffer_[mid_idx].timestamp_s - t0) / dt;
                        s = std::clamp(s, 0.0, 1.0);

                        if (run_len <= 50) {
                            double h00 = (1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s);
                            double h01 = s * s * (3.0 - 2.0 * s);

                            buffer_[mid_idx].ax = h00 * buffer_[start_idx].ax + h01 * buffer_[end_idx].ax;
                            buffer_[mid_idx].ay = h00 * buffer_[start_idx].ay + h01 * buffer_[end_idx].ay;
                            buffer_[mid_idx].az = h00 * buffer_[start_idx].az + h01 * buffer_[end_idx].az;
                            buffer_[mid_idx].gx = h00 * buffer_[start_idx].gx + h01 * buffer_[end_idx].gx;
                            buffer_[mid_idx].gy = h00 * buffer_[start_idx].gy + h01 * buffer_[end_idx].gy;
                            buffer_[mid_idx].gz = h00 * buffer_[start_idx].gz + h01 * buffer_[end_idx].gz;
                        } else {
                            buffer_[mid_idx].ax = buffer_[start_idx].ax + s * (buffer_[end_idx].ax - buffer_[start_idx].ax);
                            buffer_[mid_idx].ay = buffer_[start_idx].ay + s * (buffer_[end_idx].ay - buffer_[start_idx].ay);
                            buffer_[mid_idx].az = buffer_[start_idx].az + s * (buffer_[end_idx].az - buffer_[start_idx].az);
                            buffer_[mid_idx].gx = buffer_[start_idx].gx + s * (buffer_[end_idx].gx - buffer_[start_idx].gx);
                            buffer_[mid_idx].gy = buffer_[start_idx].gy + s * (buffer_[end_idx].gy - buffer_[start_idx].gy);
                            buffer_[mid_idx].gz = buffer_[start_idx].gz + s * (buffer_[end_idx].gz - buffer_[start_idx].gz);
                        }
                    }
                }
                i += run_len;
            } else {
                i++;
            }
        }
    }

    void apply_hampel(double out_vals[6]) {
        size_t win_len = std::min(count_, HAMPEL_WINDOW);
        if (win_len < 3) return;

        double vals[15];
        double devs[15];

        for (int axis = 0; axis < 6; ++axis) {
            for (size_t k = 0; k < win_len; ++k) {
                size_t idx = (head_ + k) % RING_CAPACITY;
                const auto& s = buffer_[idx];
                switch (axis) {
                    case 0: vals[k] = s.ax; break;
                    case 1: vals[k] = s.ay; break;
                    case 2: vals[k] = s.az; break;
                    case 3: vals[k] = s.gx; break;
                    case 4: vals[k] = s.gy; break;
                    case 5: vals[k] = s.gz; break;
                }
            }

            size_t mid_idx = win_len / 2;
            std::nth_element(vals, vals + mid_idx, vals + win_len);
            double median = vals[mid_idx];

            for (size_t k = 0; k < win_len; ++k) {
                devs[k] = std::fabs(vals[k] - median);
            }
            std::nth_element(devs, devs + mid_idx, devs + win_len);
            double mad = devs[mid_idx];

            double threshold = 4.0 * 1.4826 * mad;
            double current_val = out_vals[axis];

            if (std::fabs(current_val - median) > threshold && threshold > 1e-6) {
                out_vals[axis] = median;
            }
        }
    }
};

#endif // AZIMUTH_DSP_CLEANER_HPP
