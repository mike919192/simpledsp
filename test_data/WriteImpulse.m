
%this requires signal 1.4.4 because there is a bug in the previous version where zp2sos function return wrong result!
pkg load signal

clear all
close all
fs = 39e3;
fs2 = fs / 2;
n = 1000;
impulse = [1,zeros(1, n-1)];
step = fs/n;
x = 0:step:fs - step;
f0 = 200;
Q = 1.4;
N = 8;

%HP filter
[z, p, k] = butter ( N, f0 / fs2, "high" );
sos = zp2sos(z, p, k);
filtered = sosfilt(sos, impulse);
csvwrite("impulse_response/HPimpulse.csv", [2, fs, f0, Q, n, filtered]);

%LP filter
[z, p, k] = butter ( N, f0 / fs2 );
sos = zp2sos(z, p, k);
filtered = sosfilt(sos, impulse);
csvwrite("impulse_response/LPimpulse.csv", [1, fs, f0, Q, n, filtered]);

%BP filter
[f1, f2] = findIIRCutoffFreq(fs, f0, Q);
[z, p, k] = butter( N/2, [f1 / fs2, f2 / fs2] );
sos = zp2sos(z, p, k);
filtered = sosfilt(sos, impulse);
csvwrite("impulse_response/BPimpulse.csv", [3, fs, f0, Q, n, filtered]);

clear all
close all
fs = 39e3;
fs2 = fs / 2;
n = 1000;
impulse = [1,zeros(1, n-1)];
step = fs/n;
x = 0:step:fs - step;
f0 = 2000;
Q = 0.8;
N = 8;

%HP filter
[z, p, k] = butter ( N, f0 / fs2, "high" );
sos = zp2sos(z, p, k);
filtered = sosfilt(sos, impulse);
csvwrite("impulse_response/HPimpulse2.csv", [2, fs, f0, Q, n, filtered]);

%LP filter
[z, p, k] = butter ( N, f0 / fs2 );
sos = zp2sos(z, p, k);
filtered = sosfilt(sos, impulse);
csvwrite("impulse_response/LPimpulse2.csv", [1, fs, f0, Q, n, filtered]);

%BP filter
[f1, f2] = findIIRCutoffFreq(fs, f0, Q);
[z, p, k] = butter( N/2, [f1 / fs2, f2 / fs2] );
sos = zp2sos(z, p, k);
filtered = sosfilt(sos, impulse);
csvwrite("impulse_response/BPimpulse2.csv", [3, fs, f0, Q, n, filtered]);

clear all
close all
fs = 39e3;
fs2 = fs / 2;
n = 1000;
impulse = [1,zeros(1, n-1)];
step = fs/n;
x = 0:step:fs - step;
f0 = 15000;
Q = 2.0;
N = 8;

%HP filter
[z, p, k] = butter ( N, f0 / fs2, "high" );
sos = zp2sos(z, p, k);
filtered = sosfilt(sos, impulse);
csvwrite("impulse_response/HPimpulse3.csv", [2, fs, f0, Q, n, filtered]);

%LP filter
[z, p, k] = butter ( N, f0 / fs2 );
sos = zp2sos(z, p, k);
filtered = sosfilt(sos, impulse);
csvwrite("impulse_response/LPimpulse3.csv", [1, fs, f0, Q, n, filtered]);

%BP filter
[f1, f2] = findIIRCutoffFreq(fs, f0, Q);
[z, p, k] = butter( N/2, [f1 / fs2, f2 / fs2] );
sos = zp2sos(z, p, k);
filtered = sosfilt(sos, impulse);
csvwrite("impulse_response/BPimpulse3.csv", [3, fs, f0, Q, n, filtered]);

