addpath('f90model-main/results/')

% Step 1: Read Data
data_07_n0_N4 = load('frequency_07-02-2018_moreWind_n0_N4.dat');
data_07_n10_N4 = load('frequency_07-02-2018_moreWind_n10_N4.dat');
data_07_n100_N4 = load('frequency_07-02-2018_moreWind_n100_N4.dat');
data_07_n500_N4 = load('frequency_07-02-2018_moreWind_n500_N4.dat');

startDateTime = datetime('06-Feb-2018 23:00', 'Format', 'dd-MMM-yyyy HH:mm');
endDateTime = datetime('08-Feb-2018 01:00', 'Format', 'dd-MMM-yyyy HH:mm');
t = linspace(startDateTime, endDateTime, numel(data_07_N1));

subplot(4,1,1);
plot(t,data_07_n0_N4);
title('Wx4, n = 0','Interpreter','none')
ylabel('f (Hz)')

ylim([-1 1]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

subplot(4,1,2);
plot(t,data_07_n10_N4);
title('Wx4, n = 10','Interpreter','none')
ylabel('f (Hz)')

ylim([-1 1]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

subplot(4,1,3);
plot(t,data_07_n100_N4);
title('Wx4, n = 100','Interpreter','none')
ylabel('f (Hz)')

ylim([-1 1]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

subplot(4,1,4);
plot(t,data_07_n500_N4);
title('Wx4, n = 500','Interpreter','none')
ylabel('f (Hz)')

ylim([-1 1]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15