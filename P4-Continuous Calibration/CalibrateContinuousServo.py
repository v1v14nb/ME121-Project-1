import machine
from machine import I2C, Pin, PWM
import ssd1306
import servo
import time

# Initialize I2C and OLED display
#i2c = I2C(scl=Pin(23), sda=Pin(22))
#display = ssd1306.SSD1306_I2C(128, 64, i2c)


# Draw a rectangle and greeting message on OLED
#display.text('Calibration Mode', 0, 50, 1)
#display.show()


# Initialize servos
servo1 = servo.Servo(0, freq=50)  # GPIO for servo 1
servo2 = servo.Servo(2, freq=50)  # GPIO for servo 2

# Stop the servos initially
servo1.write_us(1500) # Servos don't move
servo2.write_us(1500)

# Results matrix to fill in
results1 = [
    [0, 35],
    [500, 36.9],
    [1000, 35.5],
    [1100, 29.44],
    [1200, 23.5],  # [PWM, # revs in 30s]
    [1300, 16.28],
    [1325, 14],
    [1387, 5.66],
    [1390, 0],
    [1500, 0],	# stopped PWM, expected to have 0 movement
    [1556, 0],
    [1557, 4.74],
    [1560, 5.2],
    [1575, 7.4],
    [1600, 9.84],
    [1700, 18.6],
    [1800, 25.26],
    [1900, 30.7],
    [2000, 35.6],
    [2500, 36.24],
    [3000, 36.6],
]

# Results matrix to fill in
results2 = [
    [0, 34.58],
    [500, 34.62],
    [1000, 32.8],
    [1200, 21.45],  # [PWM, # revs in 30s]
    [1300, 14.4],
    [1400, 5.78],
    [1425, 2.88],
    [1428, 2.95],
    [1429, 1.9],
    [1430, 0],
    [1500, 0],	# stopped PWM, expected to have 0 movement
    [1541, 0],
    [1542, 2.75],
    [1550, 3.34],
    [1575, 5.48],
    [1600, 7.75],
    [1700, 15.3],
    [1800, 21.55],
    [1900, 27.4],
    [2000, 34.08],
    [2500, 33.9],
    [3000, 34.28],
]

time.sleep(1)


# Run test
servo1.write_us(0)
#servo2.write_us(0)

time.sleep(15) #run for 30 seconds

# Stop servos
servo1.write_us(1500)
servo2.write_us(1500)
