addpath('f90model-main/results/')

% Step 1: Read Data
data_07_N1 = load('power_07-02-2018_n0_N1.dat');
data_07_N2 = load('power_07-02-2018_n0_N2.dat');
data_07_N4 = load('power_07-02-2018_n0_N4.dat');
data_07_N8 = load('power_07-02-2018_n0_N8.dat');
data_07_N16 = load('power_07-02-2018_n0_N16.dat');
data_07_demand = load('demand_07-02-2018.dat');

startDateTime = datetime('06-Feb-2018 23:00', 'Format', 'dd-MMM-yyyy HH:mm');
endDateTime = datetime('08-Feb-2018 01:00', 'Format', 'dd-MMM-yyyy HH:mm');
t = linspace(startDateTime, endDateTime, numel(data_07_N1));

plot(t,data_07_N1);
title('Wind generation','Interpreter','none')
ylabel('P (MW)')
xlabel('t (h)')
hold on
plot(t,data_07_N2);
hold on
plot(t,data_07_N4);
hold on
plot(t,data_07_N8);
hold on
plot(t,data_07_N16);
hold on
plot(t,data_07_demand,'-.');

legend('Wx1','Wx2','Wx4','Wx8','Wx16','10-min demand')

hold off