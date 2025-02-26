% Plot continuous servo motor calibration results

% data
servo1 = [
    0, 35
    500, 36.9
    1000, 35.5;
    1100, 29.44;
    1200, 23.5;
    1300, 16.28;
    1325, 14;
    1387, 5.66;
    1390, 0;
    1500, 0;
    1556, 0;
    1557, 4.74;
    1560, 5.2;
    1575, 7.4;
    1600, 9.84;
    1700, 18.6;
    1800, 25.26;
    1900, 30.7;
    2000, 35.6;
    2500, 36.24;
    3000, 36.6;
];

servo2 = [
    0, 34.58;
    500, 34.62;
    1000, 32.8;
    1200, 21.45;
    1300, 14.4;
    1400, 5.78;
    1425, 2.88;
    1428, 2.95;
    1429, 1.9;
    1430, 0;
    1500, 0;
    1541, 0;
    1542, 2.75;
    1550, 3.34;
    1575, 5.48;
    1600, 7.75;
    1700, 15.3;
    1800, 21.55;
    1900, 27.4;
    2000, 34.08;
    2500, 33.9;
    3000, 34.28;
];

% Convert revolutions per 30s to rev/s
PWM1 = servo1(:,1);
RPS1 = servo1(:,2) / 30;

PWM2 = servo2(:,1);
RPS2 = servo2(:,2) / 30;

% Plot
figure;
hold on;
plot(PWM1, RPS1, 'bo-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Servo 1 (shoulder)');
plot(PWM2, RPS2, 'ro-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Servo 2 (elbow)');

xlabel('PWM (microseconds)');
ylabel('Rate (rev/s)');
title('Continuous Servo Calibration: PWM vs Revolutions per Second');
legend('Location', 'best');
grid on;

hold off;
