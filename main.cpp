#include "src/parser.hpp"
#include "src/dsp_cleaner.hpp"
#include "src/estimator.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif

static size_t get_peak_rss_kb() {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        return static_cast<size_t>(pmc.PeakWorkingSetSize / 1024);
    }
    return 0;
#else
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
#if defined(__APPLE__)
        return static_cast<size_t>(usage.ru_maxrss / 1024);
#else
        return static_cast<size_t>(usage.ru_maxrss);
#endif
    }
    return 0;
#endif
}

static inline void emit_pose(FILE* fp, double ts, const Pose6D& pose) {
    // TUM format: timestamp tx ty tz qx qy qz qw
    std::fprintf(fp, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
                 ts,
                 pose.pos.x, pose.pos.y, pose.pos.z,
                 pose.quat.x, pose.quat.y, pose.quat.z, pose.quat.w);
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::fprintf(stderr, "Usage: %s <input_imu.csv> <output_trajectory.txt> [--meta metadata.json]\n", argv[0]);
        return 1;
    }

    const char* input_path = argv[1];
    const char* output_path = argv[2];
    const char* meta_path = nullptr;

    for (int i = 3; i < argc; ++i) {
        if (std::strcmp(argv[i], "--meta") == 0 && i + 1 < argc) {
            meta_path = argv[i + 1];
            i++;
        }
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    FastCsvParser parser;
    if (!parser.open(input_path)) {
        std::fprintf(stderr, "Error: Failed to open input CSV file: %s\n", input_path);
        return 1;
    }

    FILE* out_file = std::fopen(output_path, "wb");
    if (!out_file) {
        std::fprintf(stderr, "Error: Failed to open output trajectory file: %s\n", output_path);
        return 1;
    }

    // Fixed 64 KB segmented output buffer
    char out_buffer[65536];
    std::setvbuf(out_file, out_buffer, _IOFBF, sizeof(out_buffer));

    // TUM format header
    std::fprintf(out_file, "# timestamp tx ty tz qx qy qz qw\n");

    DspCleaner cleaner;
    cleaner.init(parser.get_median_dt());

    StateEstimator estimator;
    estimator.init(parser.get_first_timestamp(), parser.get_median_dt());

    size_t rows_in = 0;
    size_t rows_out = 0;
    double first_ts = -1.0;
    double last_ts = -1.0;

    Pose6D last_pose;
    last_pose.timestamp = 0.0;
    last_pose.pos = {0.0, 0.0, 0.0};
    last_pose.quat = {1.0, 0.0, 0.0, 0.0};

    RawImuSample raw_sample;
    CleanedImuSample cleaned_sample;

    while (parser.read_sample(raw_sample)) {
        rows_in++;
        if (first_ts < 0.0) {
            first_ts = raw_sample.timestamp_s;
        }
        last_ts = raw_sample.timestamp_s;

        // Duplicate handling (dt == 0): emit immediately to guarantee 1-to-1 mapping
        if (raw_sample.is_duplicate) {
            emit_pose(out_file, raw_sample.timestamp_s, last_pose);
            rows_out++;
            continue;
        }

        // Push into DSP cleaner ring buffer
        if (cleaner.push_sample(raw_sample)) {
            while (cleaner.pop_cleaned(cleaned_sample)) {
                last_pose = estimator.update(cleaned_sample);
                emit_pose(out_file, cleaned_sample.timestamp_s, last_pose);
                rows_out++;
                break;
            }
        }
    }

    // Flush remaining frames in ring buffer
    while (cleaner.has_pending()) {
        if (cleaner.pop_cleaned(cleaned_sample)) {
            last_pose = estimator.update(cleaned_sample);
            emit_pose(out_file, cleaned_sample.timestamp_s, last_pose);
            rows_out++;
        }
    }

    std::fflush(out_file);
    std::fclose(out_file);
    parser.close();

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_runtime_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    double duration_s = (last_ts > first_ts && first_ts >= 0.0) ? (last_ts - first_ts) : 1.0;
    double ms_per_sec_imu = (duration_s > 0.0) ? (total_runtime_ms / duration_s) : 0.0;
    size_t peak_rss_kb = get_peak_rss_kb();

    std::printf("Processed %zu input rows -> %zu output poses in %.2f ms (%.3f ms/s IMU, Peak RSS: %zu KB)\n",
                rows_in, rows_out, total_runtime_ms, ms_per_sec_imu, peak_rss_kb);

    if (meta_path) {
        FILE* meta_fp = std::fopen(meta_path, "w");
        if (meta_fp) {
            std::fprintf(meta_fp, "{\n");
            std::fprintf(meta_fp, "  \"rows_in\": %zu,\n", rows_in);
            std::fprintf(meta_fp, "  \"rows_out\": %zu,\n", rows_out);
            std::fprintf(meta_fp, "  \"total_runtime_ms\": %.3f,\n", total_runtime_ms);
            std::fprintf(meta_fp, "  \"duration_s\": %.3f,\n", duration_s);
            std::fprintf(meta_fp, "  \"ms_per_sec_imu\": %.4f,\n", ms_per_sec_imu);
            std::fprintf(meta_fp, "  \"peak_rss_kb\": %zu\n", peak_rss_kb);
            std::fprintf(meta_fp, "}\n");
            std::fclose(meta_fp);
        }
    }

    return 0;
}
