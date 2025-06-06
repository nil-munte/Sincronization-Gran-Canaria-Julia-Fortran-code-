addpath('f90model-main/results/')

% Step 1: Read Data
data_Bij1 = load('frequency_07-02-2018_moreWind_n1_N1_Bij1.dat');

startDateTime = datetime('06-Feb-2018 23:00', 'Format', 'dd-MMM-yyyy HH:mm');
endDateTime = datetime('08-Feb-2018 01:00', 'Format', 'dd-MMM-yyyy HH:mm');
t = linspace(startDateTime, endDateTime, numel(data_Bij1));

plot(t,data_Bij1);
title('B_i_j=1.0')
ylabel('f (Hz)')

%ylim([-1 1]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15