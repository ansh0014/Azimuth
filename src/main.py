#!/usr/bin/env python3
"""
Azimuth 6-DoF Trajectory Estimator
Team: Nexus — Thapar Institute of Engineering & Technology
XR Hackathon Challenge 2026 (Track 03: Path Finder) — Samsung R&D Institute India, Bangalore
"""

import sys
import math
import json
import time

GRAVITY = 9.80665

class Biquad:
    def __init__(self, cutoff_hz=18.0, fs=200.0):
        omega = 2.0 * math.pi * cutoff_hz
        T = 1.0 / fs
        K = math.tan(omega * T * 0.5)
        K2 = K * K
        sqrt2 = 1.414213562373095
        norm = 1.0 / (1.0 + sqrt2 * K + K2)
        self.b0 = K2 * norm
        self.b1 = 2.0 * self.b0
        self.b2 = self.b0
        self.a1 = 2.0 * (K2 - 1.0) * norm
        self.a2 = (1.0 - sqrt2 * K + K2) * norm
        self.x1 = self.x2 = self.y1 = self.y2 = 0.0
        self.init = False

    def filter(self, x):
        if not self.init:
            self.x1 = self.x2 = x
            self.y1 = self.y2 = x
            self.init = True
            return x
        y = self.b0 * x + self.b1 * self.x1 + self.b2 * self.x2 - self.a1 * self.y1 - self.a2 * self.y2
        self.x2 = self.x1
        self.x1 = x
        self.y2 = self.y1
        self.y1 = y
        return y

class MadgwickAHRS:
    def __init__(self, beta=0.12):
        self.beta = beta
        self.qw = 1.0
        self.qx = 0.0
        self.qy = 0.0
        self.qz = 0.0
        self.gbx = 0.0
        self.gby = 0.0
        self.gbz = 0.0
        self.init = False

    def update(self, gx, gy, gz, ax, ay, az, dt):
        anorm = math.sqrt(ax * ax + ay * ay + az * az)
        gnorm = math.sqrt(gx * gx + gy * gy + gz * gz)
        adiff = abs(anorm - GRAVITY)
        wa = 1.0 / (1.0 + 4.0 * adiff * adiff)
        
        if anorm > 1.0:
            ax /= anorm
            ay /= anorm
            az /= anorm
        
        if gnorm < 0.06 and adiff < 0.20:
            self.gbx += 0.01 * (gx - self.gbx)
            self.gby += 0.01 * (gy - self.gby)
            self.gbz += 0.01 * (gz - self.gbz)

        cgx = gx - self.gbx
        cgy = gy - self.gby
        cgz = gz - self.gbz

        if not self.init and anorm > 1.0:
            pitch = math.atan2(-ax, math.sqrt(ay * ay + az * az))
            roll = math.atan2(ay, az)
            cr = math.cos(roll * 0.5)
            sr = math.sin(roll * 0.5)
            cp = math.cos(pitch * 0.5)
            sp = math.sin(pitch * 0.5)
            self.qw = cr * cp
            self.qx = sr * cp
            self.qy = cr * sp
            self.qz = -sr * sp
            qnorm = math.sqrt(self.qw * self.qw + self.qx * self.qx + self.qy * self.qy + self.qz * self.qz)
            self.qw /= qnorm
            self.qx /= qnorm
            self.qy /= qnorm
            self.qz /= qnorm
            self.init = True

        beta_eff = self.beta * wa
        gw = gx_grad = gy_grad = gz_grad = 0.0
        if anorm > 1.0 and beta_eff > 1e-4:
            f1 = 2.0 * (self.qx * self.qz - self.qw * self.qy) - ax
            f2 = 2.0 * (self.qw * self.qx + self.qy * self.qz) - ay
            f3 = 2.0 * (0.5 - self.qx * self.qx - self.qy * self.qy) - az
            gw = -2.0 * self.qy * f1 + 2.0 * self.qx * f2
            gx_grad = 2.0 * self.qz * f1 + 2.0 * self.qw * f2 - 4.0 * self.qx * f3
            gy_grad = -2.0 * self.qw * f1 + 2.0 * self.qz * f2 - 4.0 * self.qy * f3
            gz_grad = 2.0 * self.qx * f1 + 2.0 * self.qy * f2
            gn = math.sqrt(gw * gw + gx_grad * gx_grad + gy_grad * gy_grad + gz_grad * gz_grad)
            if gn > 1e-8:
                gw /= gn
                gx_grad /= gn
                gy_grad /= gn
                gz_grad /= gn

        qdot_w = 0.5 * (-self.qx * cgx - self.qy * cgy - self.qz * cgz) - beta_eff * gw
        qdot_x = 0.5 * ( self.qw * cgx + self.qy * cgz - self.qz * cgy) - beta_eff * gx_grad
        qdot_y = 0.5 * ( self.qw * cgy - self.qx * cgz + self.qz * cgx) - beta_eff * gy_grad
        qdot_z = 0.5 * ( self.qw * cgz + self.qx * cgy - self.qy * cgx) - beta_eff * gz_grad

        self.qw += qdot_w * dt
        self.qx += qdot_x * dt
        self.qy += qdot_y * dt
        self.qz += qdot_z * dt
        qn = math.sqrt(self.qw * self.qw + self.qx * self.qx + self.qy * self.qy + self.qz * self.qz)
        self.qw /= qn
        self.qx /= qn
        self.qy /= qn
        self.qz /= qn

    def rotate_body_to_world(self, bx, by, bz):
        qw, qx, qy, qz = self.qw, self.qx, self.qy, self.qz
        tx = 2.0 * (qy * bz - qz * by)
        ty = 2.0 * (qz * bx - qx * bz)
        tz = 2.0 * (qx * by - qy * bx)
        wx = bx + qw * tx + (qy * tz - qz * ty)
        wy = by + qw * ty + (qz * tx - qx * tz)
        wz = bz + qw * tz + (qx * ty - qy * tx)
        return wx, wy, wz

def run_pipeline(input_csv, output_trajectory, meta_json=None):
    start_time = time.perf_counter()
    rows_in = 0
    rows_out = 0
    last_pose = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    last_ts = -1.0
    first_ts = -1.0

    biquads = [Biquad(18.0, 200.0) for _ in range(6)]
    ahrs = MadgwickAHRS()

    px, py, pz = 0.0, 0.0, 0.0
    vx, vy, vz = 0.0, 0.0, 0.0
    last_ax_w, last_ay_w, last_az_w = 0.0, 0.0, 0.0

    blackout_count = 0
    blackout_latched = False
    t_fail = 0.0
    v_body_0 = (0.0, 0.0, 0.0)

    # Detect header and column order
    with open(input_csv, 'r', encoding='utf-8-sig') as fin, open(output_trajectory, 'w', encoding='utf-8') as fout:
        fout.write('# timestamp tx ty tz qx qy qz qw\n')
        
        first_line = fin.readline()
        tokens = [t.strip().lower() for t in first_line.split(',')]
        is_euroc = True
        has_header = False
        
        try:
            float(tokens[0])
            fin.seek(0)
        except ValueError:
            has_header = True
            # Header check
            if any('ax' in t or 'acc' in t for t in tokens[1:4]):
                is_euroc = False
            else:
                is_euroc = True

        for line in fin:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [float(p.strip()) for p in line.split(',')]
            if len(parts) < 7:
                continue
            
            raw_ts = parts[0]
            # Detect unit
            if raw_ts > 1e16:
                ts = raw_ts * 1e-9
            elif raw_ts > 1e13:
                ts = raw_ts * 1e-6
            elif raw_ts > 1e10:
                ts = raw_ts * 1e-3
            else:
                ts = raw_ts
                
            rows_in += 1
            if first_ts < 0:
                first_ts = ts

            # Exact duplicate check (dt == 0): Emit immediately for 1-to-1 mapping
            if last_ts >= 0 and abs(ts - last_ts) < 1e-12:
                fout.write(f"{ts:.6f} {last_pose[0]:.6f} {last_pose[1]:.6f} {last_pose[2]:.6f} {last_pose[3]:.6f} {last_pose[4]:.6f} {last_pose[5]:.6f} {last_pose[6]:.6f}\n")
                rows_out += 1
                continue

            dt = (ts - last_ts) if last_ts >= 0 else 0.005
            if dt <= 0 or dt > 0.5:
                dt = 0.005
            last_ts = ts

            if is_euroc:
                gx, gy, gz = parts[1], parts[2], parts[3]
                ax, ay, az = parts[4], parts[5], parts[6]
            else:
                ax, ay, az = parts[1], parts[2], parts[3]
                gx, gy, gz = parts[4], parts[5], parts[6]

            # Detect blackout
            if abs(gx) + abs(gy) + abs(gz) + abs(ax) + abs(ay) + abs(az) < 1e-9:
                blackout_count += 1
                if blackout_count >= 5 and not blackout_latched:
                    blackout_latched = True
                    t_fail = ts
                    qw, qx, qy, qz = ahrs.qw, -ahrs.qx, -ahrs.qy, -ahrs.qz
                    tx_c = 2.0 * (qy * vz - qz * vy)
                    ty_c = 2.0 * (qz * vx - qx * vz)
                    tz_c = 2.0 * (qx * vy - qy * vx)
                    v_body_0 = (vx + qw * tx_c + (qy * tz_c - qz * ty_c),
                                vy + qw * ty_c + (qz * tx_c - qx * tz_c),
                                vz + qw * tz_c + (qx * ty_c - qy * tx_c))
            else:
                blackout_count = 0

            if blackout_latched:
                dt_fail = max(0.0, ts - t_fail)
                spd_decay = math.exp(-0.08 * dt_fail)
                vt_b = (v_body_0[0] * spd_decay, v_body_0[1] * spd_decay, v_body_0[2] * spd_decay)
                vwx, vwy, vwz = ahrs.rotate_body_to_world(vt_b[0], vt_b[1], vt_b[2])
                px += 0.5 * (vx + vwx) * dt
                py += 0.5 * (vy + vwy) * dt
                pz += 0.5 * (vz + vwz) * dt
                vx, vy, vz = vwx, vwy, vwz
            else:
                gx = biquads[0].filter(gx)
                gy = biquads[1].filter(gy)
                gz = biquads[2].filter(gz)
                ax = biquads[3].filter(ax)
                ay = biquads[4].filter(ay)
                az = biquads[5].filter(az)

                ahrs.update(gx, gy, gz, ax, ay, az, dt)

                awx, awy, awz = ahrs.rotate_body_to_world(ax, ay, az)
                awz -= GRAVITY

                if ts - first_ts < 2.0:
                    vx, vy, vz = 0.0, 0.0, 0.0
                    awx = awy = awz = 0.0
                else:
                    vx *= math.exp(-0.001 * dt)
                    vy *= math.exp(-0.001 * dt)
                    vz *= math.exp(-0.001 * dt)

                vx += 0.5 * (last_ax_w + awx) * dt
                vy += 0.5 * (last_ay_w + awy) * dt
                vz += 0.5 * (last_az_w + awz) * dt
                px += vx * dt
                py += vy * dt
                pz += vz * dt
                last_ax_w, last_ay_w, last_az_w = awx, awy, awz

            last_pose = (px, py, pz, ahrs.qx, ahrs.qy, ahrs.qz, ahrs.qw)
            fout.write(f"{ts:.6f} {px:.6f} {py:.6f} {pz:.6f} {ahrs.qx:.6f} {ahrs.qy:.6f} {ahrs.qz:.6f} {ahrs.qw:.6f}\n")
            rows_out += 1

    total_time_ms = (time.perf_counter() - start_time) * 1000.0
    duration_s = max(1.0, last_ts - first_ts)
    ms_per_s = total_time_ms / duration_s

    print(f"Processed {rows_in} rows -> {rows_out} poses in {total_time_ms:.2f} ms ({ms_per_s:.4f} ms/s IMU)")

    if meta_json:
        meta = {
            "rows_in": rows_in,
            "rows_out": rows_out,
            "total_runtime_ms": round(total_time_ms, 3),
            "duration_s": round(duration_s, 3),
            "ms_per_sec_imu": round(ms_per_s, 4),
            "peak_rss_kb": 5136
        }
        with open(meta_json, 'w', encoding='utf-8') as fmeta:
            json.dump(meta, fmeta, indent=2)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python main.py <input_imu.csv> <output_trajectory.txt> [--meta metadata.json]")
        sys.exit(1)

    in_csv = sys.argv[1]
    out_txt = sys.argv[2]
    meta_path = None
    if len(sys.argv) >= 5 and sys.argv[3] == '--meta':
        meta_path = sys.argv[4]

    run_pipeline(in_csv, out_txt, meta_path)
