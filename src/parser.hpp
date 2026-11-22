#ifndef AZIMUTH_PARSER_HPP
#define AZIMUTH_PARSER_HPP

#include "types.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <cctype>

enum class ImuColumnOrder {
    Unknown,
    EurocGyroFirst,     // timestamp, wx, wy, wz, ax, ay, az
    DegradedAccelFirst  // timestamp, ax, ay, az, gx, gy, gz
};

class FastCsvParser {
public:
    FastCsvParser()
        : file_(nullptr),
          has_header_(false),
          order_(ImuColumnOrder::Unknown),
          ts_scale_(1.0),
          accel_scale_(1.0),
          gyro_scale_(1.0),
          median_dt_(0.005),
          last_timestamp_s_(-1.0),
          first_timestamp_s_(-1.0),
          last_raw_ts_ns_(-1)
    {
    }

    ~FastCsvParser() {
        close();
    }

    bool open(const char* filepath) {
        close();
        file_ = std::fopen(filepath, "rb");
        if (!file_) {
            return false;
        }
        std::setvbuf(file_, io_buffer_, _IOFBF, sizeof(io_buffer_));
        return inspect_and_initialize(filepath);
    }

    void close() {
        if (file_) {
            std::fclose(file_);
            file_ = nullptr;
        }
    }

    double get_median_dt() const { return median_dt_; }
    double get_first_timestamp() const { return first_timestamp_s_; }
    ImuColumnOrder get_column_order() const { return order_; }

    bool read_sample(RawImuSample& sample) {
        if (!file_ || std::feof(file_)) {
            return false;
        }

        char line[512];
        while (read_line(line, sizeof(line))) {
            if (line[0] == '\0' || line[0] == '#' || line[0] == '\r' || line[0] == '\n') {
                continue;
            }
            if (parse_line(line, sample)) {
                return true;
            }
        }
        return false;
    }

private:
    FILE* file_;
    char io_buffer_[16384];
    bool has_header_;
    ImuColumnOrder order_;
    double ts_scale_;
    double accel_scale_;
    double gyro_scale_;
    double median_dt_;
    double last_timestamp_s_;
    double first_timestamp_s_;
    int64_t last_raw_ts_ns_;

    bool read_line(char* dest, size_t max_len) {
        if (!std::fgets(dest, static_cast<int>(max_len), file_)) {
            return false;
        }
        size_t len = std::strlen(dest);
        while (len > 0 && (dest[len - 1] == '\r' || dest[len - 1] == '\n')) {
            dest[len - 1] = '\0';
            --len;
        }
        return true;
    }

    bool inspect_and_initialize(const char* filepath) {
        FILE* fp = std::fopen(filepath, "rb");
        if (!fp) return false;

        char line[512];
        double sample_ts[40];
        double sample_diffs[40];
        double norm1_list[40];
        double norm2_list[40];
        size_t n_samples = 0;
        size_t n_diffs = 0;

        bool header_checked = false;
        order_ = ImuColumnOrder::Unknown;

        while (std::fgets(line, sizeof(line), fp) && n_samples < 40) {
            size_t len = std::strlen(line);
            while (len > 0 && (line[len - 1] == '\r' || line[len - 1] == '\n')) {
                line[len - 1] = '\0';
                --len;
            }
            if (len == 0) continue;

            if (!header_checked) {
                header_checked = true;
                if (line[0] == '#' || !std::isdigit(static_cast<unsigned char>(line[0]))) {
                    has_header_ = true;
                    char header_buf[512];
                    size_t copy_len = (len < sizeof(header_buf) - 1) ? len : (sizeof(header_buf) - 1);
                    for (size_t i = 0; i < copy_len; ++i) {
                        header_buf[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(line[i])));
                    }
                    header_buf[copy_len] = '\0';

                    const char* w_ptr = std::strstr(header_buf, "w_");
                    if (!w_ptr) w_ptr = std::strstr(header_buf, "gyro");
                    const char* a_ptr = std::strstr(header_buf, "a_");
                    if (!a_ptr) a_ptr = std::strstr(header_buf, "acc");

                    if (w_ptr && a_ptr) {
                        if (w_ptr < a_ptr) {
                            order_ = ImuColumnOrder::EurocGyroFirst;
                        } else {
                            order_ = ImuColumnOrder::DegradedAccelFirst;
                        }
                    }
                    continue;
                }
            }

            double vals[7];
            int count = parse_doubles(line, vals, 7);
            if (count >= 7) {
                sample_ts[n_samples] = vals[0];
                if (n_samples > 0) {
                    double dt = sample_ts[n_samples] - sample_ts[n_samples - 1];
                    if (dt > 0.0) {
                        sample_diffs[n_diffs++] = dt;
                    }
                }
                norm1_list[n_samples] = std::sqrt(vals[1]*vals[1] + vals[2]*vals[2] + vals[3]*vals[3]);
                norm2_list[n_samples] = std::sqrt(vals[4]*vals[4] + vals[5]*vals[5] + vals[6]*vals[6]);
                n_samples++;
            }
        }
        std::fclose(fp);

        if (n_samples == 0) return false;

        double raw_median_dt = 1.0;
        if (n_diffs > 0) {
            std::sort(sample_diffs, sample_diffs + n_diffs);
            raw_median_dt = sample_diffs[n_diffs / 2];
        }

        if (raw_median_dt < 1.0) {
            ts_scale_ = 1.0;
        } else if (raw_median_dt < 1e3) {
            ts_scale_ = 1e-3;
        } else if (raw_median_dt < 1e6) {
            ts_scale_ = 1e-6;
        } else {
            ts_scale_ = 1e-9;
        }
        median_dt_ = raw_median_dt * ts_scale_;
        if (median_dt_ <= 0.0) median_dt_ = 0.005;

        // Always sort norm lists before computing medians for order or unit scale detection
        std::sort(norm1_list, norm1_list + n_samples);
        std::sort(norm2_list, norm2_list + n_samples);

        double med_norm1 = norm1_list[n_samples / 2];
        double med_norm2 = norm2_list[n_samples / 2];

        if (order_ == ImuColumnOrder::Unknown) {
            bool norm1_is_accel = (med_norm1 > 5.0 && med_norm1 < 20.0) || (med_norm1 > 0.5 && med_norm1 < 2.0);
            bool norm2_is_accel = (med_norm2 > 5.0 && med_norm2 < 20.0) || (med_norm2 > 0.5 && med_norm2 < 2.0);

            if (norm1_is_accel && !norm2_is_accel) {
                order_ = ImuColumnOrder::DegradedAccelFirst;
            } else if (norm2_is_accel && !norm1_is_accel) {
                order_ = ImuColumnOrder::EurocGyroFirst;
            } else {
                double diff1 = std::fabs(med_norm1 - GRAVITY_STANDARD);
                double diff2 = std::fabs(med_norm2 - GRAVITY_STANDARD);
                order_ = (diff1 < diff2) ? ImuColumnOrder::DegradedAccelFirst : ImuColumnOrder::EurocGyroFirst;
            }
        }
        if (order_ == ImuColumnOrder::Unknown) {
            order_ = ImuColumnOrder::DegradedAccelFirst;
        }

        double accel_norm = (order_ == ImuColumnOrder::DegradedAccelFirst) ? med_norm1 : med_norm2;
        if (accel_norm > 0.5 && accel_norm < 2.0) {
            accel_scale_ = GRAVITY_STANDARD; // g to m/s^2
        } else {
            accel_scale_ = 1.0;
        }

        double gyro_norm = (order_ == ImuColumnOrder::DegradedAccelFirst) ? med_norm2 : med_norm1;
        if (gyro_norm > 10.0) {
            gyro_scale_ = M_PI / 180.0; // deg/s to rad/s
        } else {
            gyro_scale_ = 1.0;
        }

        std::fseek(file_, 0, SEEK_SET);
        if (has_header_) {
            read_line(line, sizeof(line));
        }
        return true;
    }

    int parse_doubles(char* str, double* out, int max_vals) {
        int count = 0;
        char* cur = str;
        while (*cur && count < max_vals) {
            while (*cur == ' ' || *cur == '\t') cur++;
            if (!*cur) break;
            char* end = nullptr;
            out[count] = std::strtod(cur, &end);
            if (end == cur) break;
            count++;
            cur = end;
            while (*cur == ' ' || *cur == '\t') cur++;
            if (*cur == ',') cur++;
        }
        return count;
    }

    bool parse_line(char* line, RawImuSample& sample) {
        char* cur = line;
        while (*cur == ' ' || *cur == '\t') cur++;
        if (!*cur) return false;

        int64_t raw_ts_ns = -1;
        double ts_sec = 0.0;
        char* end = nullptr;

        if (ts_scale_ == 1e-9) {
            raw_ts_ns = std::strtoll(cur, &end, 10);
            if (end == cur) return false;
            if (*end == '.') {
                ts_sec = std::strtod(cur, &end) * ts_scale_;
                raw_ts_ns = static_cast<int64_t>(ts_sec * 1e9);
            } else {
                int64_t sec = raw_ts_ns / 1000000000LL;
                int64_t nsec = raw_ts_ns % 1000000000LL;
                ts_sec = static_cast<double>(sec) + static_cast<double>(nsec) * 1e-9;
            }
            sample.raw_timestamp = static_cast<double>(raw_ts_ns);
            sample.timestamp_s = ts_sec;
        } else {
            double raw_ts = std::strtod(cur, &end);
            if (end == cur) return false;
            sample.raw_timestamp = raw_ts;
            sample.timestamp_s = raw_ts * ts_scale_;
        }

        cur = end;
        while (*cur == ' ' || *cur == '\t') cur++;
        if (*cur == ',') cur++;

        double vals[6];
        int count = parse_doubles(cur, vals, 6);
        if (count < 6) return false;

        if (first_timestamp_s_ < 0.0) {
            first_timestamp_s_ = sample.timestamp_s;
            last_timestamp_s_ = sample.timestamp_s;
            last_raw_ts_ns_ = raw_ts_ns;
            sample.is_duplicate = false;
        } else {
            bool is_dup = false;
            if (ts_scale_ == 1e-9 && raw_ts_ns >= 0 && last_raw_ts_ns_ >= 0) {
                is_dup = (raw_ts_ns == last_raw_ts_ns_);
            } else {
                is_dup = (std::fabs(sample.timestamp_s - last_timestamp_s_) < 1e-12);
            }
            if (is_dup) {
                sample.is_duplicate = true;
            } else {
                sample.is_duplicate = false;
                last_timestamp_s_ = sample.timestamp_s;
                last_raw_ts_ns_ = raw_ts_ns;
            }
        }

        if (order_ == ImuColumnOrder::DegradedAccelFirst) {
            sample.ax = vals[0] * accel_scale_;
            sample.ay = vals[1] * accel_scale_;
            sample.az = vals[2] * accel_scale_;
            sample.gx = vals[3] * gyro_scale_;
            sample.gy = vals[4] * gyro_scale_;
            sample.gz = vals[5] * gyro_scale_;
        } else {
            sample.gx = vals[0] * gyro_scale_;
            sample.gy = vals[1] * gyro_scale_;
            sample.gz = vals[2] * gyro_scale_;
            sample.ax = vals[3] * accel_scale_;
            sample.ay = vals[4] * accel_scale_;
            sample.az = vals[5] * accel_scale_;
        }
        return true;
    }
};

#endif // AZIMUTH_PARSER_HPP
