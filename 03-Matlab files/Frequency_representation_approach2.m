addpath('f90model-main/results/')

% Step 1: Read Data
data_n0 = load('frequency_07-02-2018_moreWind_n0_N1.dat');
data_n2 = load('frequency_07-02-2018_moreWind_n2_N1.dat');
data_n5 = load('frequency_07-02-2018_moreWind_n5_N1.dat');
data_n10 = load('frequency_07-02-2018_moreWind_n10_N1.dat');
data_n20 = load('frequency_07-02-2018_moreWind_n20_N1.dat');
data_n50 = load('frequency_07-02-2018_moreWind_n50_N1.dat');
data_n100 = load('frequency_07-02-2018_moreWind_n100_N1.dat');
data_n500 = load('frequency_07-02-2018_moreWind_n500_N1.dat');

startDateTime = datetime('06-Feb-2018 23:00', 'Format', 'dd-MMM-yyyy HH:mm');
endDateTime = datetime('08-Feb-2018 01:00', 'Format', 'dd-MMM-yyyy HH:mm');
t = linspace(startDateTime, endDateTime, numel(data_n0));

subplot(4,1,1);
plot(t,data_n0);
title('n = 0')
ylabel('f (Hz)')

ylim([-0.5 0.5]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

subplot(4,1,2);
plot(t,data_n10);
title('n = 10')
ylabel('f (Hz)')

ylim([-0.5 0.5]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

subplot(4,1,3);
plot(t,data_n50);
title('n = 50')
ylabel('f (Hz)')

ylim([-0.5 0.5]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

subplot(4,1,4);
plot(t,data_n500);
title('n = 500')
ylabel('f (Hz)')

ylim([-0.5 0.5]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15
