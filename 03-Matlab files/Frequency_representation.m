addpath('f90model-main/results/')

% Step 1: Read Data
data_07 = load('frequency_07-02-2018.dat');
data_08 = load('frequency_08-02-2018.dat');

startDateTime = datetime('06-Feb-2018 23:00', 'Format', 'dd-MMM-yyyy HH:mm');
endDateTime = datetime('08-Feb-2018 01:00', 'Format', 'dd-MMM-yyyy HH:mm');
t = linspace(startDateTime, endDateTime, numel(data_07));

subplot(2,1,1);
plot(t,data_07);
title('frequency_07-02-2018.dat','Interpreter','none')
ylabel('f (Hz)')

ylim([-0.5 0.5]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

startDateTime = datetime('07-Feb-2018 23:00', 'Format', 'dd-MMM-yyyy HH:mm');
endDateTime = datetime('09-Feb-2018 01:00', 'Format', 'dd-MMM-yyyy HH:mm');
t = linspace(startDateTime, endDateTime, numel(data_08));

subplot(2,1,2);
plot(t,data_08);
title('frequency_08-02-2018.dat','Interpreter','none')
ylabel('f (Hz)')

ylim([-0.5 0.5]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

%figure();

%plot(t,data_n10);
%title('Proves','Interpreter','none')
%ylabel('f (Hz)')

%ylim([-0.5 0.5]);

%yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
%yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

