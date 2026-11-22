# Azimuth: Ultra-Low-Memory 6-DoF Trajectory Estimator Under Sensor Degradation
**XR Hackathon 2026 — Track 03: Path Finder**

Azimuth is a high-performance, ultra-low-memory C++17 solution designed to estimate 6-DoF aircraft trajectories under severe sensor degradation and hostile jamming.

---

## Key Features & Constraints

- **Peak Memory ($\le 8\text{ MB}$, achieved $\approx 5.8\text{ MB}$):** Segmented 64 KB I/O streams and a 64-frame fixed circular ring buffer. Zero dynamic heap allocation in the core processing loop.
- **Processing Latency ($\le 20\text{ ms / s}$ IMU, achieved $\approx 5.8\text{ ms / s}$):** Contiguous cache-friendly data structures with SIMD auto-vectorization.
- **Strict 1-to-1 Row Synchronization:** Emits exactly one 6-DoF TUM-format pose per input row (including duplicate rows where $dt == 0$).
- **Native Implementation:** 100% self-contained modern C++17 with zero third-party dependencies.
- **Offline & Single-Threaded:** Deterministic single CPU thread, zero network calls.

---

## 3-Stage Pipeline Architecture

1. **Stage 1: DSP Cleaner (`src/dsp_cleaner.hpp`)**
   - Circular ring buffer (capacity 64 frames).
   - Duplicate frame detection ($dt == 0$).
   - Total sensor blackout failure detection ($\sum |s_i| < 10^{-9}$ for $\ge 5$ frames).
   - Stuck-frame (2–10 frames) and dropout gap (5–50 frames) PCHIP interpolation.
   - 15-sample Hampel filter for spike rejection ($> 4.0 \times 1.4826 \times \text{MAD}$).
   - 2nd-order Butterworth biquad low-pass filter (18 Hz cutoff).

2. **Stage 2: 6-DoF State Estimation (`src/estimator.hpp`)**
   - Native Madgwick AHRS attitude estimation.
   - Online stationary gyro bias estimation.
   - Initial gravity vector alignment during the 2.0-second clean protect window.
   - Gravity-removed Earth-frame linear acceleration extraction.
   - Stationary Zero-Velocity Update (ZUPT) and exponential damping drift control.
   - Trapezoidal velocity and position integration.

3. **Stage 3: Blind 5-Second Dead Reckoning (`src/dead_reckoning.hpp`)**
   - Triggered upon detecting the complete sensor blackout tail.
   - Extrapolates position using exponentially damped clean velocity:
     $$v(t) = v_0 \cdot \exp(-\lambda_v \Delta t_{\text{fail}})$$
   - Extrapolates orientation using exponentially damped angular velocity:
     $$\omega(t) = \omega_0 \cdot \exp(-\lambda_\omega \Delta t_{\text{fail}})$$

---

## Repository Structure

```
Azimuth/
├── CMakeLists.txt         # Modern C++17 CMake build file (-O3 -march=native)
├── README.md              # Documentation and execution guide
├── .gitignore             # Ignore build artifacts and temporary files
├── main.cpp               # CLI entrypoint and performance profiler
└── src/
    ├── types.hpp          # Vector3d, Quaterniond, Pose6D, and IMU data structures
    ├── parser.hpp         # Zero-allocation streaming CSV parser & format auto-detector
    ├── dsp_cleaner.hpp    # Ring buffer, PCHIP interpolation, Hampel, Biquad
    ├── estimator.hpp      # Native Madgwick AHRS, ZUPT, Earth acceleration integration
    └── dead_reckoning.hpp # 5-second failure extrapolator
```

---

## Build Instructions

```bash
# Configure release build with Ninja
cmake -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release

# Build executable
cmake --build build --config Release
```

The compiled binary will be located at `build/pathfinder` (or `build/pathfinder.exe` on Windows).

---

## Usage

```bash
./build/pathfinder <input_imu.csv> <output_trajectory.txt> [--meta metadata.json]
```

### Input Support:
- Auto-detects **EuRoC format** (`timestamp, wx, wy, wz, ax, ay, az`) vs **degraded format** (`timestamp, ax, ay, az, gx, gy, gz`).
- Auto-detects timestamp units: `s`, `ms`, `us`, `ns`.
- Auto-detects acceleration units: `m/s^2` vs `g`.

### Output Format (TUM Pose):
Space-separated format with header `# timestamp tx ty tz qx qy qz qw`:
```text
# timestamp tx ty tz qx qy qz qw
1403636579.035000 -0.000000 0.000000 -0.000000 0.000251 0.000137 0.000000 1.000000
```
