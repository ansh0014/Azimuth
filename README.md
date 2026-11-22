# Azimuth: Ultra-Low-Memory 6-DoF Inertial Trajectory Estimator

**XR Hackathon Challenge 2026 — Track 03: Path Finder**  
**Organized by:** Samsung R&D Institute India, Bangalore  
**Team Name:** Nexus  
**Institution:** Thapar Institute of Engineering & Technology  
**Team Members:** Kavish Arora • Anshul Jagota • Dhriti Dhall  
**Presentation:** [Thapar_Nexus.pptx](Thapar_Nexus.pptx)  

Azimuth is a native, single-threaded, zero-heap-allocation system (available in both high-performance C++17 and Python) that estimates 6-DoF trajectory poses (`timestamp tx ty tz qx qy qz qw`) from degraded 6-axis IMU streams subject to progressive white noise, random-walk bias drift, scale factor drift, high-magnitude sensor spikes, stuck frames, burst dropouts, duplicate rows, and complete sensor blackout intervals.

---

## 1. System Specifications & Measured Performance

| System Metric | Architecture Specification | Measured Performance |
| :--- | :--- | :--- |
| **Language & Standard** | C++17 (`-O3 -march=native -fno-exceptions -fno-rtti -s`) | Zero external dependencies |
| **Peak Memory (RSS)** | Segmented 16 KB stream buffers + fixed stack arrays | **5136 KB (5.01 MB)** |
| **Compute Latency** | Single-pass streaming execution | **3.697 ms per 1 s of IMU data (270.5x real-time)** |
| **Execution Model** | Single CPU thread, deterministic | Zero threading or dynamic allocation overhead |
| **Frame Synchronization** | 1-to-1 pose output per input row | **11706 in -> 11706 out** (including `dt == 0` duplicates) |
| **Runtime Environment** | 100% Offline | Zero network or external runtime calls |

---

## 2. Algorithmic Architecture (3-Stage Pipeline)

### Stage 0: Zero-Allocation Streaming Parser (`src/parser.hpp`)
- **Segmented I/O:** Uses 16 KB `setvbuf` stream buffers (`char[16384]`) on standard `FILE*` descriptors without dynamic heap allocations.
- **Format & Unit Auto-Detection:** Inspects the header and the first 40 rows using fixed stack arrays (`double[40]`) to automatically detect:
  - **Column Order:** EuRoC order (`timestamp, wx, wy, wz, ax, ay, az`) vs. accel-first order (`timestamp, ax, ay, az, gx, gy, gz`) via header tokens and median 3-axis vector norms.
  - **Timestamp Units:** Auto-detects `s`, `ms`, `us`, or `ns` from the median sampling interval `median_dt`. For 19-digit nanosecond timestamps, parses via `std::strtoll` into `int64_t` and splits integer seconds and remainder nanoseconds prior to floating-point conversion to prevent IEEE-754 precision loss.
  - **Sensor Units:** Normalizes acceleration to `m/s^2` (detecting `g` vs `m/s^2`) and angular velocity to `rad/s` (detecting `deg/s` vs `rad/s`).

### Stage 1: Causal DSP Cleaning (`src/dsp_cleaner.hpp`)
Processes each incoming sample in a 1-in, 1-out zero-lag streaming pipeline:
1. **Duplicate Detection (`dt == 0`):** Compares exact 64-bit nanosecond timestamps (`raw_ts_ns == last_raw_ts_ns_`). When a duplicate timestamp is encountered, the DSP filter and state integrator are bypassed, and the current pose is emitted immediately to preserve 1-to-1 row alignment without corrupting filter state.
2. **Total Sensor Blackout Detection:** Monitors the L1 norm of all 6 sensor axes (`|ax| + |ay| + |az| + |gx| + |gy| + |gz| < 1e-9`). When 5 consecutive non-duplicate zero frames are detected, the pipeline latches into Stage 3 (Blind Dead Reckoning).
3. **Stuck-Frame Detection & Extrapolation:** Identifies runs of `>= 2` consecutive non-duplicate frames with identical 6-axis readings (`diff < 1e-6`) and extrapolates local sample gradients so the integrator does not freeze during sensor lockups.
4. **15-Sample Hampel Spike Rejection Filter:** Maintains a 15-sample sliding window per axis on stack memory. Computes the sample median and Median Absolute Deviation (`MAD`) via `std::nth_element`. Any sample exceeding `|x - median| > 4.0 * 1.4826 * MAD` is replaced with the window median, rejecting `6g` accelerometer and `1 rad/s` gyroscope spikes.
5. **2nd-Order Butterworth Biquad Low-Pass Filter:** Attenuates progressive high-frequency white noise across all 6 axes using a bilinear-transformed 2nd-order Butterworth filter with an `18.0 Hz` cutoff frequency dynamically parameterized by the detected sampling rate `fs = 1.0 / median_dt`:
   ```
   y[n] = b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] - a1 * y[n-1] - a2 * y[n-2]
   ```

### Stage 2: 6-DoF State Estimation & Drift Control (`src/estimator.hpp`)
1. **Initial 2.0-Second Stationary Calibration (`protect_seconds_ = 2.0`):**
   - Initializes roll and pitch attitude from the initial gravity vector.
   - Accumulates baseline gyroscope readings (`sum_init_gyro_`) and gravity-removed world accelerations (`sum_init_world_accel_`) across the initial stationary window to calibrate `init_gyro_bias` and `world_accel_bias_`.
2. **Adaptive-Gain Madgwick/Mahony Attitude Filter (`NativeMadgwickAhrs`):**
   - Tracks 3D orientation as a unit quaternion `q = [qw, qx, qy, qz]` mapping body-frame vectors to the world frame.
   - **Maneuver Acceleration Rejection:** Computes deviation from standard gravity `accel_diff = ||a|| - 9.80665|` and scales the accelerometer gradient-descent correction gain via a fast rational bell curve:
     ```
     w_a = 1.0 / (1.0 + 4.0 * accel_diff^2)
     beta_eff = beta * w_a
     ```
     During high-g turns or linear accelerations, `w_a -> 0`, preventing centripetal and linear accelerations from corrupting roll and pitch estimates.
   - **Quasi-Static Gyro Bias Tracking:** Continuously refines `gyro_bias_` during low-rotation, 1g flight segments (`||w|| < 0.06 rad/s` and `accel_diff < 0.20 m/s^2`).
3. **Gravity Removal & Leaky High-Pass Accelerometer Bias Tracker:**
   - Rotates body specific force into the world frame and subtracts Earth gravity:
     ```
     a_world_raw = q * a_body * q^-1 - [0, 0, 9.80665]^T
     ```
   - Tracks and removes random-walk accelerometer bias using an online exponential moving average with time constant `tau_bias = 4.0 s`:
     ```
     alpha = dt / (tau_bias + dt)
     world_accel_bias = world_accel_bias + alpha * (a_world_raw - world_accel_bias)
     a_world = a_world_raw - world_accel_bias
     ```
   - Clamps residual acceleration norm to `15.0 m/s^2`.
4. **Trapezoidal Integration, Velocity Damping & ZUPT:**
   - Integrates velocity and position via the trapezoidal rule:
     ```
     v_k = v_{k-1} + 0.5 * (a_{k-1} + a_k) * dt
     p_k = p_{k-1} + 0.5 * (v_{k-1} + v_k) * dt
     ```
   - Applies exponential velocity damping `v_k *= exp(-0.6 * dt)` to bound quadratic double-integration drift during un-aided inertial flight, and strong Zero-Velocity Update (ZUPT) damping `v_k *= exp(-8.0 * dt)` when stationary (`||w|| < 0.08` and `accel_diff < 0.25`).
   - Maintains a 256-entry circular history buffer of recent clean velocities and angular rates for seamless Stage 3 handoff.

### Stage 3: Coordinated Flight Arc Dead Reckoning (`src/dead_reckoning.hpp`)
- Triggered automatically when a total sensor blackout (`0.0` on all 6 axes for `>= 5` consecutive frames) is detected.
- Latches the last healthy pose `(p_0, q_0)`, the smoothed 150-sample average velocity `v_0`, and the 50-sample average angular velocity `w_0`.
- Transforms `v_0` into the aircraft body frame (`v_body_0 = q_0^-1 * v_0 * q_0`).
- At each blackout timestamp `t_k` (`dt_fail = t_k - t_fail_start`):
  1. Integrates orientation using exponentially damped angular velocity `w(t) = w_0 * exp(-0.12 * dt_fail)`:
     ```
     q_k = normalize(q_{k-1} * exp(0.5 * w(t) * dt))
     ```
  2. Projects decayed body velocity `v_body(t) = v_body_0 * exp(-0.08 * dt_fail)` along the evolving body heading into the world frame (`v_world(t) = q_k * v_body(t) * q_k^-1`) and integrates position via the trapezoidal rule.

---

## 3. Benchmark Results

Evaluated under progressive fault injection (`SEED = 42`, white noise, bias random walk, scale drift, `6g` spikes, stuck runs of 2-10 samples, burst dropouts of 5-50 samples, duplicate rows, and a 5.0-second total sensor blackout tail) across both a 60-second standard flight profile and a 75-second high-maneuver generalization profile:

| Benchmark Metric | Sequence 1 (60s, 11706 rows) | Sequence 2 (75s High-Maneuver, 14621 rows) |
| :--- | :--- | :--- |
| **Pose Synchronization** | `11706 in -> 11706 out` (100% Match) | `14621 in -> 14621 out` (100% Match) |
| **Peak Memory (RSS)** | **5136 KB (5.01 MB)** | **5132 KB (5.01 MB)** |
| **Total Execution Time** | **221.80 ms** | **224.13 ms** |
| **Compute Latency** | **3.697 ms/s IMU** (`270.5x` real-time) | **2.989 ms/s IMU** (`334.6x` real-time) |
| **Active Flight ATE (RMSE)** | **3.436 m** | **6.619 m** |
| **Blind 5s Dead-Reckoning ATE (RMSE)** | **4.108 m** | **9.072 m** |
| **Relative Pose Error (RPE @ 1s)** | **1.683 m/s** | **3.286 m/s** |
| **Quaternion Norm Error (`max ||q|| - 1|`)** | **6.86e-07** | **6.83e-07** |

---

## 4. Repository Structure

```text
Azimuth/
├── CMakeLists.txt         # C++17 build configuration (-O3 -march=native -fno-exceptions -fno-rtti -s)
├── README.md              # Architecture, algorithms, benchmarks, and build documentation
├── main.cpp               # C++ CLI entrypoint, 16 KB segmented stream buffers, RSS & latency profiler
└── src/
    ├── main.py            # Standalone Python entrypoint (Official Evaluator compliant)
    ├── main.cpp           # C++ entrypoint copy
    ├── CMakeLists.txt     # Build configuration
    ├── types.hpp          # Vector3d, Quaterniond, Pose6D, and IMU frame structures
    ├── parser.hpp         # Zero-heap CSV stream parser, 19-digit ns precision, format auto-detector
    ├── dsp_cleaner.hpp    # Zero-lag causal DSP cleaner: duplicate, stuck-frame, Hampel, 18 Hz Biquad
    ├── estimator.hpp      # Adaptive Madgwick AHRS, 2s calibration, leaky accel bias tracker, ZUPT
    └── dead_reckoning.hpp # Stage 3 coordinated flight arc 5-second blackout extrapolator
```

---

## 5. Build & Execution Instructions

### Option A: Python Execution (Direct / Standalone)
Requires standard Python 3.x (zero external dependencies required):
```bash
python src/main.py <input_imu.csv> <output_trajectory.txt> [--meta metadata.json]
```

### Option B: C++ Native Compilation & Execution
Prerequisites: CMake >= 3.16, C++17 Compiler (GCC / MinGW-w64, Clang, or MSVC)
```bash
cmake -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
./build/pathfinder <input_imu.csv> <output_trajectory.txt> [--meta metadata.json]
```

---

## 6. AI Disclosure

In accordance with the XR Hackathon Challenge 2026 transparency guidelines:
- **AI Tool Assistance:** Generative AI tools were utilized during development for code structuring, documentation generation, syntax validation, Markdown formatting, and scaffolding test utilities.
- **Original Engineering & Implementation:** All core algorithmic implementations—including the zero-allocation streaming parser, causal 2nd-order Butterworth biquad filter, Hampel MAD outlier rejection, adaptive-gain Madgwick attitude estimation with maneuver acceleration suppression, world-acceleration drift control, and coordinated flight arc dead-reckoning extrapolator—were designed, implemented, mathematically validated, and benchmarked by Team Nexus.
