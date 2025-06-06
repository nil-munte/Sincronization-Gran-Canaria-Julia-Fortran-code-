addpath('f90model-main/results/')

% Step 1: Read Data
data_07_n500 = load('frequency_07-02-2018_moreWind_n500_N1.dat'); % Bij = 500, node_to_connect = 1
data_07_n500_N1_1 = load('frequency_07-02-2018_moreWind_n500_N1_nous_valors.dat'); % Bij = 1000, node_to_connect = 23
data_07_n500_N1_2 = load('frequency_07-02-2018_moreWind_n500_N1_nous_valors_2.dat'); % Bij = 100, node_to_connect = 8
data_07_n500_N1_3 = load('frequency_07-02-2018_moreWind_n500_N1_nous_valors_3.dat'); % Bij = 10000, node_to_connect = 2

startDateTime = datetime('06-Feb-2018 23:00', 'Format', 'dd-MMM-yyyy HH:mm');
endDateTime = datetime('08-Feb-2018 01:00', 'Format', 'dd-MMM-yyyy HH:mm');
t = linspace(startDateTime, endDateTime, numel(data_07_n500));

subplot(4,1,1);
plot(t,data_07_n500);
title('Wx1, n = 500, B_i_j = 500, node_s_t_a_r = 1')
ylabel('f (Hz)')

ylim([-1 1]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

subplot(4,1,2);
plot(t,data_07_n500_N1_1);
title('Wx1, n = 500, B_i_j = 1000, node_s_t_a_r = 23')
ylabel('f (Hz)')

ylim([-1 1]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

subplot(4,1,3);
plot(t,data_07_n500_N1_2);
title('Wx1, n = 500, B_i_j = 100, node_s_t_a_r = 8')
ylabel('f (Hz)')

ylim([-1 1]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15

subplot(4,1,4);
plot(t,data_07_n500_N1_3);
title('Wx1, n = 500, B_i_j = 10000, node_s_t_a_r = 2')
ylabel('f (Hz)')

ylim([-1 1]);

yline(0.15, '--r', 'LineWidth', 1.5);  % Dotted line at y = 0.15
yline(-0.15, '--r', 'LineWidth', 1.5); % Dotted line at y = -0.15
