%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Purpose: ELEC-E8105 - Project work,...
%%%          calculates and prints data from /savedata/*
%%% Author:  Jani Arponen, 769684
%%% Date:    2020-04-03
%%% Info:    Main script: proj_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants to memory from proj_main.m
R = 0.01;
% RMSE 
rmse = @(x,hat) sqrt(mean((x - hat).^2));
%% quantizer 6
files = dir('./savedata/*_*_100_0*.mat'); % 0 = just tracking, 1 = EKF, 2 = UKF

Tab = [];
for i = 1:length(files)
    load(sprintf('%s\\%s',files(i).folder,files(i).name));
    s = split(files(i).name,'_'); % seed, estimator, dt, bitcount
    sr = split(s(5),'.'); % loss
    newrow = [
        str2double(s(1)),str2double(sr(1)),... % seed and packetloss
        rmse(X(1,:),X_ekf(1,:)),... % x pos ekf 
        rmse(X(1,:),X_ukf(1,:)),... % x pos ukf
        rmse(X(2,:),X_ekf(2,:)),... % y pos ekf
        rmse(X(2,:),X_ukf(2,:)),... % y pos ukf
        sum(U*U'.*Cr,'all'),...     % quadratic sum for LQR cost
        max(tt),...                 % final k
        str2double(s(2))            %estimator
        ];
    Tab = [Tab;newrow];
end
%% cleanup
clear files i newrow s sr;
ss = unique(Tab(:,1)); % unique seeds.
% separate data into matrices that can be plotted with surf(x,y,z)
%bb = unique(Tab(:,1)); % bitcounts
%bb = bb(bb > 0); % remove 0
pp = unique(Tab(:,2)); % packetlsoses
ee = unique(Tab(:,end)); % estimators
for s_ = 1:size(ss,1)
    for p_ = 1:size(pp,1) % inner loop all packetlossses
        for e_ = 1:size(ee,1) % inner inner loop for estimators
            % find index with unique seed-loss-estimator pair
            
            ind = and(and(Tab(:,1) == ss(s_),Tab(:,2) == pp(p_)),Tab(:,end) == ee(e_));
            
            if e_ == 1
                xe(s_,p_) = Tab(ind,3);
                ye(s_,p_) = Tab(ind,5);
                ue(s_,p_) = Tab(ind,7);
                ke(s_,p_) = Tab(ind,8);
            end
            if e_ == 2
                xu(s_,p_) = Tab(ind,4);
                yu(s_,p_) = Tab(ind,6);
                uu(s_,p_) = Tab(ind,7);
                ku(s_,p_) = Tab(ind,8);
            end
        end
    end
end
xem = mean(xe,1);
yem = mean(ye,1);
uem = mean(ue,1);
kem = mean(ke,1);
xum = mean(xu,1);
yum = mean(yu,1);
uum = mean(uu,1);
kum = mean(ku,1);
figure(1244);clf;
subplot(4,1,1)
hold on
plot(pp,xem,'-', 'Color',c.ekf,'LineWidth',2)
plot(pp,xum,'-', 'Color',c.ukf,'LineWidth',2)
for i_ = 1:size(ss,1)
    plot(pp,xe(i_,:),'o','Color',c.ekf)
    plot(pp,xu(i_,:),'o','Color',c.ukf)
end
xlabel('p_{loss} [%]')
ylabel('RMSE of x')
ylim([-0.5 3])
xlim([-1 51])
legend('EKF mean','UKF mean','EKF sample','UKF sample')
title('Performance parameters as a function of packet loss for EKF and UKF controlled systems')

subplot(4,1,2)
hold on
plot(pp,yem,'-','Color',c.ekf,'LineWidth',2)
plot(pp,yum,'-','Color',c.ukf,'LineWidth',2)
for i_ = 1:size(ss,1)
    plot(pp,ye(i_,:),'o','Color',c.ekf)
    plot(pp,yu(i_,:),'o','Color',c.ukf)
end
xlabel('p_{loss} [%]')
ylabel('RMSE of y')
ylim([-0.5 3])
xlim([-1 51])
legend('EKF mean','UKF mean','EKF sample','UKF sample')

subplot(4,1,3)
hold on
plot(pp,uem,'Color',c.ekf,'LineWidth',2)
plot(pp,uum,'Color',c.ukf,'LineWidth',2)
for i_ = 1:size(ss,1)
    plot(pp,ue(i_,:),'o','Color',c.ekf)
    plot(pp,uu(i_,:),'o','Color',c.ukf)
end
xlabel('p_{loss} [%]')
ylabel('Input cost')
ylim([300 1000])
xlim([-1 51])
legend('EKF mean','UKF mean','EKF sample','UKF sample')


subplot(4,1,4)
hold on
plot(pp,kem,'Color',c.ekf,'LineWidth',2)
plot(pp,kum,'Color',c.ukf,'LineWidth',2)
for i_ = 1:size(ss,1)
    plot(pp,ke(i_,:),'o','Color',c.ekf)
    plot(pp,ku(i_,:),'o','Color',c.ukf)
end
xlabel('p_{loss} [%]')
ylabel('final k');
ylim([200 1050])
xlim([-1 51])
legend('EKF mean','UKF mean','EKF sample','UKF sample')

