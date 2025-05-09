from time import ticks_ms, ticks_add, ticks_diff
from MPU6050 import MPU6050  # ← Updated import to match your file name
import os

# Initialize MPU6050 (your driver already sets pins)
mpu = MPU6050()

# Generate unique filename
base = "data_log"
i = 0
while "{}{:03d}.csv".format(base, i) in os.listdir():
    i += 1
filename = "{}{:03d}.csv".format(base, i)

# Open CSV file and write header
with open(filename, "w") as f:
    f.write("Time(ms),Temp,AcX,AcY,AcZ\n")

    start_time = ticks_ms()
    INTERVAL = 10  # ms between readings
    next_time = ticks_add(start_time, INTERVAL)

    print("Logging to:", filename)
    print("Press Ctrl+C to stop.\n")

    try:
        while True:
            now = ticks_ms()
            timestamp = ticks_diff(now, start_time)

            accel = mpu.read_accel_data()  # In m/s^2 by default
            temp = mpu.read_temperature()

            ax = accel["x"]
            ay = accel["y"]
            az = accel["z"]

            print(f"{timestamp}ms - T:{temp:.2f}°C, X:{ax:.2f}, Y:{ay:.2f}, Z:{az:.2f}")
            f.write(f"{timestamp},{temp:.2f},{ax:.4f},{ay:.4f},{az:.4f}\n")  # <-- FIXED HERE

            next_time = ticks_add(next_time, INTERVAL)
            while ticks_diff(next_time, ticks_ms()) > 0:
                pass


    except KeyboardInterrupt:
        print("\nLogging stopped.")

