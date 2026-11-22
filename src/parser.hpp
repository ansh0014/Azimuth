#ifndef AZIMUTH_PARSER_HPP
#define AZIMUTH_PARSER_HPP

#include "types.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
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
          first_timestamp_s_(-1.0)
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
    char io_buffer_[65536];
    bool has_header_;
    ImuColumnOrder order_;
    double ts_scale_;
    double accel_scale_;
    double gyro_scale_;
    double median_dt_;
    double last_timestamp_s_;
    double first_timestamp_s_;

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
        std::vector<double> sample_ts;
        std::vector<double> sample_diffs;
        std::vector<double> norm1_list;
        std::vector<double> norm2_list;

        bool header_checked = false;
        order_ = ImuColumnOrder::Unknown;

        while (std::fgets(line, sizeof(line), fp) && sample_ts.size() < 40) {
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
                    std::string header_str(line);
                    std::transform(header_str.begin(), header_str.end(), header_str.begin(), [](unsigned char c) {
                        return static_cast<char>(std::tolower(c));
                    });
                    size_t w_pos = header_str.find("w_");
                    if (w_pos == std::string::npos) w_pos = header_str.find("gyro");
                    size_t a_pos = header_str.find("a_");
                    if (a_pos == std::string::npos) a_pos = header_str.find("acc");

                    if (w_pos != std::string::npos && a_pos != std::string::npos) {
                        if (w_pos < a_pos) {
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
                sample_ts.push_back(vals[0]);
                if (sample_ts.size() > 1) {
                    double dt = sample_ts.back() - sample_ts[sample_ts.size() - 2];
                    if (dt > 0.0) {
                        sample_diffs.push_back(dt);
                    }
                }
                double norm1 = std::sqrt(vals[1]*vals[1] + vals[2]*vals[2] + vals[3]*vals[3]);
                double norm2 = std::sqrt(vals[4]*vals[4] + vals[5]*vals[5] + vals[6]*vals[6]);
                norm1_list.push_back(norm1);
                norm2_list.push_back(norm2);
            }
        }
        std::fclose(fp);

        if (sample_ts.empty()) return false;

        double raw_median_dt = 1.0;
        if (!sample_diffs.empty()) {
            std::sort(sample_diffs.begin(), sample_diffs.end());
            raw_median_dt = sample_diffs[sample_diffs.size() / 2];
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

        if (order_ == ImuColumnOrder::Unknown && !norm1_list.empty()) {
            std::sort(norm1_list.begin(), norm1_list.end());
            std::sort(norm2_list.begin(), norm2_list.end());
            double med_norm1 = norm1_list[norm1_list.size() / 2];
            double med_norm2 = norm2_list[norm2_list.size() / 2];

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

        double accel_norm = GRAVITY_STANDARD;
        if (!norm1_list.empty()) {
            accel_norm = (order_ == ImuColumnOrder::DegradedAccelFirst) ? norm1_list[norm1_list.size() / 2] : norm2_list[norm2_list.size() / 2];
        }
        if (accel_norm > 0.5 && accel_norm < 2.0) {
            accel_scale_ = GRAVITY_STANDARD; // g to m/s^2
        } else {
            accel_scale_ = 1.0;
        }

        double gyro_norm = 0.05;
        if (!norm1_list.empty()) {
            gyro_norm = (order_ == ImuColumnOrder::DegradedAccelFirst) ? norm2_list[norm2_list.size() / 2] : norm1_list[norm1_list.size() / 2];
        }
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
        double vals[7];
        int count = parse_doubles(line, vals, 7);
        if (count < 7) return false;

        sample.raw_timestamp = vals[0];
        sample.timestamp_s = vals[0] * ts_scale_;

        if (first_timestamp_s_ < 0.0) {
            first_timestamp_s_ = sample.timestamp_s;
            last_timestamp_s_ = sample.timestamp_s;
            sample.is_duplicate = false;
        } else {
            if (std::fabs(sample.timestamp_s - last_timestamp_s_) < 1e-12) {
                sample.is_duplicate = true;
            } else {
                sample.is_duplicate = false;
                last_timestamp_s_ = sample.timestamp_s;
            }
        }

        if (order_ == ImuColumnOrder::DegradedAccelFirst) {
            sample.ax = vals[1] * accel_scale_;
            sample.ay = vals[2] * accel_scale_;
            sample.az = vals[3] * accel_scale_;
            sample.gx = vals[4] * gyro_scale_;
            sample.gy = vals[5] * gyro_scale_;
            sample.gz = vals[6] * gyro_scale_;
        } else {
            sample.gx = vals[1] * gyro_scale_;
            sample.gy = vals[2] * gyro_scale_;
            sample.gz = vals[3] * gyro_scale_;
            sample.ax = vals[4] * accel_scale_;
            sample.ay = vals[5] * accel_scale_;
            sample.az = vals[6] * accel_scale_;
        }
        return true;
    }
};

#endif // AZIMUTH_PARSER_HPP
