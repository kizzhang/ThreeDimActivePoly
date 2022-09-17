traj_data = importdata('S:\Programs\Active Polymer\ABP polymer\walk_trajectory.txt',' ', 3);
time = 100;
par_num = str2double(cell2mat(traj_data.textdata(1)));
Rg2 = zeros(time,1);
time_step = 2e-6;
% Parameters:
D = 26.85; Dr = 100; v0 = 1; b = sqrt(3*D/(4*Dr)); Pe= v0*sqrt(3/(4*Dr*D));

%% MAIN CODE %%
%%  Rg vs. N
Rg2N = zeros(par_num,1);
Rg2Ndata = traj_data.data((time-1)*100+time:time*100+time-1,2:end);

for n = 1:par_num
    R_0_vec  = sum(Rg2Ndata(1:n,:))/n;
    Rg2N(n) = sum(sum((Rg2Ndata(1:n,:) - R_0_vec) .* (Rg2Ndata(1:n,:) - R_0_vec),2))/n; 
end

Rg_vs_N = zeros(par_num,3);
for dn = 1:par_num
   deltaRg2N = Rg2N(1+dn:end) - Rg2N(1:end-dn);

   Rg_vs_N(dn,1) = mean(deltaRg2N); %# averagse
   Rg_vs_N(dn,2) = std(deltaRg2N); %# std
   Rg_vs_N(dn,3) = length(deltaRg2N); %# n
end
figure()
f1 = errorbar((1:par_num), Rg_vs_N(:,1), Rg_vs_N(:,2),'bo');  hold on;
%f1 = plot(1:par_num, Rg2N,'bo'); hold on;
f2 = plot((1:par_num), 2/9*2*b^2*(1:par_num) *4/3*(1+2/3*Pe^2),'r-');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel("N")
ylabel("Rg")

legend([f1,f2], {'simulation', 'theoretical prediction'})

%% MSD
for i = 1:time
    data = traj_data.data((i-1)*100+i:i*100+i-1,2:end);

    R0_x = sum(traj_data.data((i-1)*100+i:i*100+i-1,2))/par_num;
    R0_y = sum(traj_data.data((i-1)*100+i:i*100+i-1,3))/par_num;
    R0_z = sum(traj_data.data((i-1)*100+i:i*100+i-1,4))/par_num;
    
    Rg2_x = (data(:,1) - R0_x) .* (data(:,1) - R0_x);
    Rg2_y = (data(:,2) - R0_y) .* (data(:,2) - R0_y);
    Rg2_z = (data(:,3) - R0_z) .* (data(:,3) - R0_z);
    
    Rg2(i) = ((Rg2_x + Rg2_y + Rg2_z)' * (Rg2_x + Rg2_y + Rg2_z)) /par_num;
end 
figure();
plot((1:time)*time_step,  Rg2,  'r')
% figure(2);
% 
% p1 = errorbar(time_step*(1:time), msd(:,1), msd(:,2),'r'); 
