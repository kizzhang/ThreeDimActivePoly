traj_data = importdata('S:\Programs\Active Polymer\ABP polymer\walk_trajectory.txt',' ', 3);
par_num = str2double(cell2mat(traj_data.textdata(1)));
time = 100;
data = traj_data.data(:,2:4);
nData = size(data,1); %# number of data points
x = zeros(time,3);

for i = 1:time
    x(i,1) = data((i-1)*par_num+i+1,1);
    x(i,2) = data((i-1)*par_num+i+1,2);
    x(i,3) = data((i-1)*par_num+i+1,3);
end

numberOfDeltaT = time; %# for MSD, dt should be up to 1/4 of number of data points

msd = zeros(numberOfDeltaT,3); %# We'll store [mean, std, n]

%# calculate msd for all deltaT's

for dt = 1:numberOfDeltaT
   deltaCoords = x(1+dt:end,1:3) - x(1:end-dt,1:3);
   squaredDisplacement = sum(deltaCoords.^2,2); %# dx^2+dy^2+dz^2

   msd(dt,1) = mean(squaredDisplacement); %# averagse
   msd(dt,2) = std(squaredDisplacement); %# std
   msd(dt,3) = length(squaredDisplacement); %# n
end
figure(1);

f1= errorbar(0.0001*(1:numberOfDeltaT), msd(:,1), msd(:,2),'r');  hold on;

%% Active 
time_t = 0.0001 *(1:numberOfDeltaT);
D_t = 26.85; D_r = 100;
v_0 = 0;
D_id = D_t + v_0^2 / (6*D_r);
tau_r = 1/(2*D_r);
%plot(time_t, 6 * D_id * time_t, 'r-'); hold on;
f2=plot(time_t, 6*D_t*time_t + 2 / ((2*D_r)^2) * v_0^2 ... 
        *(time_t/tau_r + exp(-time_t/tau_r) - 1), 'k-');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel("Time")
ylabel("MSD")

legend([f1,f2], {'simulation', 'theoretical prediction'})


% %% Passive
% time_t = 0.01 *(1:numberOfDeltaT);
% D_t = 0.02; 
% 
% f2 = plot(time_t, 6 * D_t * time_t, 'r-'); hold on;
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlabel("Time")
% ylabel("MSD")
% legend([f1,f2], {'Simulation', '6Dt'})

