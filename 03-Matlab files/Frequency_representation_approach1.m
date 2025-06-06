addpath('f90model-main/results/')

% Step 1: Read Data
data_n0 = load('frequency_07-02-2018_mod_n0.dat');
data_n1 = load('frequency_07-02-2018_mod_n1.dat');
data_n5 = load('frequency_07-02-2018_mod_n5.dat');
data_n10 = load('frequency_07-02-2018_mod_n10.dat');
data_n20 = load('frequency_07-02-2018_mod_n20.dat');
data_n50 = load('frequency_07-02-2018_mod_n50.dat');
data_n100 = load('frequency_07-02-2018_mod_n100.dat');
data_n200 = load('frequency_07-02-2018_mod_n200.dat');
data_n500 = load('frequency_07-02-2018_mod_n500.dat');
data_n10000 = load('frequency_07-02-2018_mod_n10000.dat');

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
plot(t,data_n10000);
title('n = 10000')
ylabel('f (Hz)')

ylim([-0.5 0.5]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15
