%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Purpose: ELEC-E8105 - Project work
%%% Author:  Jani Arponen, 769684
%%% Date:    2020-04-03
%%% Info:    Heavily based on Supplemental Matlab code for Exercise 4.3,...
%%%          which is under GNU GPL v. 2 or later,...
%%%          Original code by Simo Särkkä(?).
%%%          -----
%%%          Check the Run options: section to modify simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;clc;% close all;
%seed = 41;
seed = 769684; % baseline and quantizer,
%seed = 893729638;
%seed = 951287; % 65\% packet loss results
%seeds = [112,956,212,310,388,416,883,720]; % chosen from random.org 1-1000 roll, none diverge at least
%for sss = 1:size(seeds,2)
%seed = seeds(sss);
rng(seed,'twister');

%% Run options:
% model:
dt = 0.1; % sampling time
q1 = 0.05;% System noise 1
q2 = 0.05;% System noise 2
r  = 0.05;% measurement noise

% channel / quant
ploss = 0.0;   % packet loss probability [0...1]
bitcount = 0; % 0 -> no quantization

% which estimator to use for controller?
%   0: direct access to state X
%   1: EKF
%   2: UKF
estimator = 0;

% waypoints
waypoints = [0 0; 10 0; 10 10; 0 10; 0 0]'; % square

wpx = [0 10 9.5 8.5 7 5 3 1.5 0.5 0];
wpy = [0, abs(sqrt(5^2 - (wpx(2:end) -5).^2))];

%waypoints = [wpx;wpy]; % semicircle
%waypoints = [0 0; 10 10]';
wp = 1;   % index for which to follow
threshold = 1/5; % norm threshold (radius) for waypoint arrival

% params to quit simulation early
maxsteps = 1000;
os = 5; % max overshoot for x,y coordinates

% Sensor positions:
S{1}.pos = [-1;5];
S{2}.pos = [5;11];

% LQR cost matrices:
Cq = eye(4);     % state cost
Cr = 100*eye(2); % input cost

%% Plotting options:
c.true = '#000000'; % black
c.meas = '#AAAAAA'; % gray
c.lost = '#FF5733'; % red
c.ref  = '#5DADE2'; % blue
c.ctrl = '#FFC300'; % yellow
c.ekf  = '#58D68D'; % green
c.ukf  = '#7D3C98'; % purpe

%% Functions
rmse = @(x,hat) sqrt(mean(sum((x - hat).^2)));
grnd = @(m,P) P*randn(size(m,1),1); % gaussian noise
urnd = @(p) (rand > p); % times by 1 to get double instead of logical

%% Model

A =[1, 0,dt, 0;  % x
    0, 1, 0,dt;  % y
    0, 0, 1, 0;  % xdot
    0, 0, 0, 1]; % ydot

B =[0, 0;
    0, 0;
    1, 0;
    0, 1];

sysd = ss(A,B,eye(4),0,dt); % system object

% system noise

Q =[q1*dt^3/3, 0, q1*dt^2/2, 0;
    0, q2*dt^3/3, 0, q2*dt^2/2;
    q1*dt^2/2, 0, q1*dt,     0;
    0, q2*dt^2/2, 0, q2*dt];
Lq = chol(Q);

%% measurement equations:
% eq:s
h   = @(x,y,s) atan2(y - s(2), x - s(1)); % measurement from sensor s
hqu = @(x,y,s,mi,ma,b) quant_uni(h(x,y,s),mi,ma,b); % quantized measurement
Hs  = @(m) [h(m(1),m(2),S{1}.pos);h(m(1),m(2),S{2}.pos)]; % meas for UKF
%Hsq = @(m,mi,ma,b) [hqu(m(1),m(2),S{1}.pos,mi,ma,b);hqu(m(1),m(2),S{2}.pos,mi,ma,b)]; % quantized meas for UKF?
Hsq = @(m,mi,ma,b) Hs(m);
dhx = @(x,y,s) (s(2) - y) / (s(2)^2 - 2*s(2)*y + s(1)^2 - 2*s(1)*x + x^2 + y^2); % deriv of x of sensor s
dhy = @(x,y,s) (x - s(1)) / (s(2)^2 - 2*s(2)*y + s(1)^2 - 2*s(1)*x + x^2 + y^2); % deriv of y of sensor s
H   = @(x,y) [dhx(x,y,S{1}.pos) dhy(x,y,S{1}.pos) 0 0;dhx(x,y,S{2}.pos) dhy(x,y,S{2}.pos) 0 0]; % jacobian for sensors S{1} and S{2}
% noise:
R  = r^2 * eye(2);
Lr = r * eye(2);

%% Controller
Kc = lqr(sysd,Cq,Cr); % LQR gain
C =[1, 0, 0, 0;    % reference can only target x and y
    0, 1, 0, 0];
Kr = inv(C*inv(eye(size(A)) - A - B*Kc)*B); % reference

%% Initials
ref = waypoints(:,1);
% draw circles around waypoints with threshhold radius
a = linspace(0,2*pi,720);
ax = threshold*cos(a);
ay = threshold*sin(a);
% break on diverge
max_wp = [max(waypoints(1,:)) + os;max(waypoints(2,:)) + os];
min_wp = [min(waypoints(1,:)) - os;min(waypoints(2,:)) - os];

x0 = [0 0 0 0]';
X(:,1) = x0;
% initialize sizes
Y(:,1) = C*x0; % meas
D = Y; % detections for plotting
U = Y; % control signal
T = Y; % reference signals for plotting

% initial state uncertainty for filters
P0 = 0.01*eye(4);
p_ekf = P0;
p_ukf = P0;
% store:
P_ekf(:,:,1) = P0; P_ukf = P_ekf;

m0 = grnd(x0,P0);
m_ekf = m0;
m_ukf = m0;

% UKF preliminaries:
N = size(m_ukf,1);
alpha = 1; beta = 0; kappa = 0;
lambda = alpha^2 * (N + kappa) - N;
WM = zeros(2*N + 1, 1);
WC = zeros(2*N + 1, 1);
for n = 1:2*N+1 % weights calculation
    if n == 1
        wm = lambda / (N + lambda);
        wc = lambda / (N + lambda) + (1 - alpha^2 + beta);
    else
        wm = 1 / (2 * (N + lambda));
        wc = wm;
    end
    WM(n) = wm;
    WC(n) = wc;
end
% sigmapoints matrix
SX0 = [zeros(size(m_ukf)) Lq -Lq]; 

%% simulate data

for e_ = 0:1:2
    estimator = e_;
%for p_ = 0:0.05:0.5
%    ploss = p_;
%for b_ = 16:-1:0
%    bitcount = b_;

    
% add quantizer noise to measurement noise
if bitcount > 0
    fprintf('Using %d-bit quantizer\n',bitcount);
    P_bits = eye(2) * (2^-bitcount /sqrt(12))^2;
    R  = r^2 * eye(2) + P_bits;
    %Lr = chol(R);
end
figure(estimator + 1);clf;
k = 1;
diverged = 0;
wp = 1;

p_ekf = P0;
p_ukf = P0;

m_ekf = m0;
m_ukf = m0;

% clear for multi run:
clear X X_ekf X_ukf U Y D P T Ds Dl
X(:,1) = x0;
X_ekf(:,1) = m0;
X_ukf(:,1) = m0;


while 1
    % true state propagation
    if k > 1
        % True state:
        X(:,k) = A*X(:,k-1) + B*U(:,k-1) + grnd(zeros(4,1),Lq);
        
        % EKF predict:
        m_ekf = A*m_ekf + B*U(:,k-1);
        p_ekf = A*P_ekf(:,:,k-1)*A' + Q;
    end
    % measurement
    Y(:,k) = [h(X(1,k),X(2,k),S{1}.pos); h(X(1,k),X(2,k),S{2}.pos)] + grnd(zeros(2,1),Lr);
    if bitcount > 0 % quantize 
        Y(:,k) = [quant_uni(Y(1,k),-4*pi/5,pi/5,bitcount);
                  quant_uni(Y(2,k),-4*pi/5,pi/5,bitcount)];
    end
    
    % Compute crossing of the measurements
    dx1 = cos(Y(1,k)); dy1 = sin(Y(1,k));
    dx2 = cos(Y(2,k)); dy2 = sin(Y(2,k));
    d = [dx1 dx2; dy1 dy2]\[S{2}.pos(1)-S{1}.pos(1);S{2}.pos(2)-S{1}.pos(2)]; %combine distance
    D(:,k) = S{1}.pos + [dx1;dy1]*d(1); % crossing at this point
    
    % measurement over channel
    pl = urnd(ploss); % packet lost?
    P(:,k) = pl;      % store for plotting
    
    % EKF update:
    if pl  % only run update if we actually got a measurement
        H_ekf = H(m_ekf(1),m_ekf(2)); % measurement jacobian
        S_ekf = H_ekf * p_ekf * H_ekf' + R; % innovation covariance
        K_ekf = p_ekf * H_ekf' / S_ekf;     % ekf gain
        m_ekf = m_ekf + K_ekf * (Y(:,k) - [h(m_ekf(1),m_ekf(2),S{1}.pos);h(m_ekf(1),m_ekf(2),S{2}.pos)]);
        p_ekf = p_ekf - K_ekf * S_ekf * K_ekf';
    end
    % store:
    X_ekf(:,k) = m_ekf;
    P_ekf(:,:,k) = p_ekf;
    
    % UKF:
    % form sigmapoints on current mean
    SX = sqrt(N + lambda)*SX0 + repmat(m_ukf,1,size(SX0,2));
    % propagate through dynamic model:
    HX_ukf = zeros(size(SX));
    for n = 1: size(HX_ukf,2)
        if k > 1 % skip first ones input amount
            HX_ukf(:,n) = A*SX(:,n) + B*U(:,k-1);
        else
            HX_ukf(:,n) = A*SX(:,n) + 0;
        end
    end
    
    % predict...
    m_ukf = zeros(size(m_ukf));
    p_ukf = zeros(size(p_ukf));
    for n = 1:size(HX_ukf,2) % mean
        m_ukf = m_ukf + WM(n) * HX_ukf(:,n);
    end
    for n = 1:size(HX_ukf,2) % covar
       p_ukf = p_ukf + WC(n) * (HX_ukf(:,n) - m_ukf) * (HX_ukf(:,n) - m_ukf)';
    end
    
    % UKF update
    if pl % only run update if we actually got a measurement
        % form sigmapoints for measurement
        SY = sqrt(N + lambda)*SX0 + repmat(m_ukf,1,size(SX0,2));
        HY_ukf = zeros(size(ref,1),size(SY,2));
        for n = 1: size(HY_ukf,2) % propagate through meas model
            if  bitcount > 0
                HY_ukf(:,n) = Hsq(SY(:,n),-pi*4/5,pi/5,bitcount);
            else
                HY_ukf(:,n) = Hs(SY(:,n));
            end
        end
        mu_ukf = zeros(size(ref));
        S_ukf  = zeros(size(ref,1),size(ref,1));
        C_ukf  = zeros(size(m_ukf,1),size(ref,1));
        for n = 1:size(SY,2)
            mu_ukf = mu_ukf + WM(n) * HY_ukf(:,n);
        end
        for n = 1:size(SY,2)
            S_ukf = S_ukf + WC(n) * (HY_ukf(:,n) - mu_ukf) * (HY_ukf(:,n) - mu_ukf)';
            C_ukf = C_ukf + WC(n) * (SX(:,n) - m_ukf) * (HY_ukf(:,n) - mu_ukf)';
        end
        S_ukf = S_ukf + R;      % meas noise
        K_ukf = C_ukf / S_ukf;  % Kalman gian
        m_ukf = m_ukf + K_ukf*(Y(:,k) - mu_ukf); % update mean
        p_ukf = p_ukf - K_ukf * S_ukf * K_ukf';  % update covar
    end
    X_ukf(:,k) = m_ukf;
    P_ukf(:,:,k) = p_ukf;
    
    
    % which estimator to follow?
    switch (estimator)
        case 0
            xu = X(:,k);
            ctrl_str = 'None';
        case 1
            xu = m_ekf;
            ctrl_str = 'EKF';
        case 2
            xu = m_ukf;
            ctrl_str = 'UKF';
    end
    
    %%% FILTERS END
    % control calculation
    u = -Kc*xu - Kr*ref;
    % cap control signal:
    %if (u(1) > 10) u(1) = 10; elseif (u(2) > 10) u(2) = 10; end
    %if (u(1) <-10) u(1) =-10; elseif (u(2) <-10) u(2) =-10; end
    U(:,k) = u;
    
    % trajectory storage
    T(:,k) = ref;
    % reference prop through measurement equation
    %href = [h(ref(1),ref(2),S{1}.pos); h(ref(1),ref(2),S{2}.pos)];
    % Choose next waypoint if we're within threshold of current and
    %   suffieciently slowed down
    if (norm(ref - xu(1:2)) < threshold && norm([0;0] - xu(3:4)) < 0.5) || (wp == 1 == k)
        wp = wp + 1;
        if wp > size(waypoints,2)
            break;
        end
        ref = waypoints(:,wp);
    end
    % break on diverge of real or any valid estimate
    if(X(1,k) > max_wp(1) || X(1,k) < min_wp(1) || ...
       X(2,k) > max_wp(2) || X(2,k) < min_wp(2) || ...
      (estimator ~= 2 && ...
      (X_ekf(1,k) > max_wp(1) || X_ekf(1,k) < min_wp(1) || ...
       X_ekf(2,k) > max_wp(2) || X_ekf(2,k) < min_wp(2)))|| ...
      (estimator ~= 1 && ...
      (X_ukf(1,k) > max_wp(1) || X_ukf(1,k) < min_wp(1) || ...
       X_ukf(2,k) > max_wp(2) || X_ukf(2,k) < min_wp(2))))
        diverged = 1;
        fprintf('Diverged at %d\n',k)
        break;
    end
    k = k + 1;
    if k > maxsteps
        break
    end
end
% measurement plotting stuffs
ind = P == 1;   % succesful transmits
Ds = D;         % succesful measurements
Dl = D;         % lost measurements
Ds(:,~ind) = nan;
Dl(:,ind)  = nan;


%% Store data
% save data to subfolder savedata with filename format:
%   'rng seed'_'used estimator'_'sampling time in ms'_'quantizer bitcount'_'packet loss prob in %'.mat (as integers)
if ~diverged && wp >= size(waypoints,2) && 0
    fname = sprintf('savedata/%d_%d_%d_%d_%d.mat',seed,estimator,dt*1000,bitcount,round(ploss*100))

    tt = 1:1:k; %time steps
    save(fname,'tt','X','X_ekf','X_ukf','T','Y','Ds','Dl','U','P_ekf','P_ukf');
end % save if

%end % quantizer
%end % packetloss
%end % estimator
%end % seed changer
%% Plot
%yl = [min(waypoints(2,:)) - os, max(waypoints(2,:)) + os];
yl = [-2, 12];
xl = yl; % limits
figure(estimator + 1);clf;
subplot(2,2,[1 3])
hold on
grid on
clear l;
plot(X(1,:),X(2,:),'--','Color',c.true,'LineWidth',2);
l = ['True','|'];
if estimator ~= 2
    plot(X_ekf(1,:),X_ekf(2,:),'Color',c.ekf);
    l = [l, 'EKF','|'];
end
if estimator ~= 1
    plot(X_ukf(1,:),X_ukf(2,:),'Color',c.ukf);
    l = [l, 'UKF','|'];
end
plot(Ds(1,:),Ds(2,:),'o','Color',c.meas);
l = [l, 'Meas.','|'];
if ploss > 0
    plot(Dl(1,:),Dl(2,:),'x','Color',c.lost);
    l = [l, 'Lost meas.','|'];
end
plot(S{1}.pos(1),S{1}.pos(2),'*','Color',c.true);
l = [l, 'Sensor 1','|'];
plot(S{2}.pos(1),S{2}.pos(2),'x','Color',c.true);
l = [l, 'Sensor 2','|'];
for n = 1:size(waypoints,2)
    plot(ax + waypoints(1,n),...
         ay + waypoints(2,n),...
         '--','Color',c.ref,...
         'LineWidth',2);
    plot(waypoints(1,n),waypoints(2,n),'o','Color',c.ref);
end
l = [l, 'Waypoint(s)'];
legend(split(l,'|'),'Location','best');
%legend('True','EKF','UKF','Meas.','Lost meas.',...
%    'Sensor 1','Sensor 2','Waypoint',...
%    'Location','best');
RMSE_EKF = [rmse(X(1,:),X_ekf(1,:)); 
            rmse(X(2,:),X_ekf(2,:));
            rmse(X(3,:),X_ekf(3,:)); 
            rmse(X(4,:),X_ekf(4,:))];
RMSE_UKF = [rmse(X(1,:),X_ukf(1,:)); 
            rmse(X(2,:),X_ukf(2,:));
            rmse(X(3,:),X_ukf(3,:)); 
            rmse(X(4,:),X_ukf(4,:))];
US(:,estimator + 1) = sum(sum(U*U'*Cr));
KS(:,estimator + 1) = k;
RMSES(:,:,estimator + 1) = [RMSE_EKF, RMSE_UKF];
switch (estimator)
    case 0
        title_str = sprintf('Trajectory using true state for controller input.\nRMSE EKF: e_x = %.3f, e_y %.3f, RMSE UKF: e_x = %.3f, e_y %.3f\nInput cost: %.3f, total timesteps: %d',...
            RMSE_EKF(1), RMSE_EKF(2),RMSE_UKF(1), RMSE_UKF(2),...
            US(1),KS(1));
    case 1
        title_str = sprintf('Trajectory using EKF estimate for controller input.\nRMSE EKF: e_x = %.3f, e_y %.3f\nInput cost: %.3f, total timesteps: %d',...
            RMSE_EKF(1), RMSE_EKF(2),US(2),KS(2));
    case 2
        title_str = sprintf('Trajectory using UKF estimate for controller input.\nRMSE UKF: e_x = %.3f, e_y %.3f\nInput cost: %.3f, total timesteps: %d',...
            RMSE_UKF(1), RMSE_UKF(2),US(3),KS(3));
end
title(title_str);
xlabel('X');ylabel('Y');
xlim(xl);ylim(yl);
daspect([1 1 1]);
axis equal
if 1 %noplot
% x only
subplot(4,2,2)
hold on
clear l;
plot(T(1,:),'-', 'Color',c.ref);
l = ['Ref','|'];
plot(X(1,:),'--','Color',c.true,'LineWidth',2);
l = [l,'True','|'];
if estimator ~= 2
    plot(X_ekf(1,:), 'Color',c.ekf);
    l = [l,'EKF','|'];
end
if estimator ~= 1
    plot(X_ukf(1,:), 'Color',c.ukf);
    l = [l,'UKF','|'];
end
plot(Ds(1,:),'o','Color',c.meas);
l = [l,'Meas.'];
if ploss > 0
    plot(Dl(1,:),'x','Color',c.lost);
    l = [l,'|','Lost meas.'];
end
legend(split(l,'|'),'Location','best');
xlabel('Time');ylabel('X');
ylim(xl);
% y only
subplot(4,2,4)
hold on
clear l;
plot(T(2,:),'-', 'Color',c.ref);
l = ['Ref','|'];
plot(X(2,:),'--','Color',c.true,'LineWidth',2);
l = [l,'True','|'];
if estimator ~= 2
    plot(X_ekf(2,:), 'Color',c.ekf);
    l = [l,'EKF','|'];
end
if estimator ~= 1
    plot(X_ukf(2,:), 'Color',c.ukf);
    l = [l,'UKF','|'];
end
plot(Ds(2,:),'o','Color',c.meas);
l = [l,'Meas.'];
if ploss > 0
    plot(Dl(2,:),'x','Color',c.lost);
    l = [l,'|','Lost meas.'];
end
legend(split(l,'|'),'Location','best');
xlabel('Time');ylabel('Y');
ylim(yl);
% speeds
subplot(4,2,6)
hold on
clear l;
plot(X(3,:),'--','Color',c.true,'LineWidth',2);
l = ['True','|'];
if estimator ~= 2
    plot(X_ekf(3,:),'-','Color',c.ekf);
    l = [l,'EKF','|'];
end
if estimator ~= 1
    plot(X_ukf(3,:),'-','Color',c.ukf);
    l = [l,'UKF','|'];
end
plot(U(1,:),'-','Color',c.ctrl);
l = [l,'Ctrl'];
legend(split(l,'|'),'Location','best');
xlabel('Time');ylabel('X Dot');

subplot(4,2,8)
hold on
clear l;
plot(X(4,:),'--','Color',c.true,'LineWidth',2);
l = ['True','|'];
if estimator ~= 2
    plot(X_ekf(4,:),'-','Color',c.ekf);
    l = [l,'EKF','|'];
end
if estimator ~= 1
    plot(X_ukf(4,:),'-','Color',c.ukf);
    l = [l,'UKF','|'];
end
plot(U(2,:),'-','Color',c.ctrl);
l = [l,'Ctrl'];
legend(split(l,'|'),'Location','best');
xlabel('Time');ylabel('Y Dot');
end


figure(216798);

subplot(3,1,(estimator + 1))
hold on
grid on
clear l;
plot(X(1,:),X(2,:),'--','Color',c.true,'LineWidth',2);
l = ['True','|'];
if estimator ~= 2
    plot(X_ekf(1,:),X_ekf(2,:),'Color',c.ekf);
    l = [l, 'EKF','|'];
end
if estimator ~= 1
    plot(X_ukf(1,:),X_ukf(2,:),'Color',c.ukf);
    l = [l, 'UKF','|'];
end
plot(Ds(1,:),Ds(2,:),'o','Color',c.meas);
l = [l, 'Meas.','|'];
if ploss > 0
    plot(Dl(1,:),Dl(2,:),'x','Color',c.lost);
    l = [l, 'Lost meas.','|'];
end
%plot(S{1}.pos(1),S{1}.pos(2),'*','Color',c.true);
%l = [l, 'Sensor 1','|'];
%plot(S{2}.pos(1),S{2}.pos(2),'x','Color',c.true);
%l = [l, 'Sensor 2','|'];
for n = 1:size(waypoints,2)
    plot(ax + waypoints(1,n),...
         ay + waypoints(2,n),...
         '--','Color',c.ref,...
         'LineWidth',2);
    plot(waypoints(1,n),waypoints(2,n),'o','Color',c.ref);
end
l = [l, 'Waypoint(s)'];
legend(split(l,'|'),'Location','best');
%legend('True','EKF','UKF','Meas.','Lost meas.',...
%    'Sensor 1','Sensor 2','Waypoint',...
%    'Location','best');
RMSE_EKF = [rmse(X(1,:),X_ekf(1,:)); 
            rmse(X(2,:),X_ekf(2,:));
            rmse(X(3,:),X_ekf(3,:)); 
            rmse(X(4,:),X_ekf(4,:))];
RMSE_UKF = [rmse(X(1,:),X_ukf(1,:)); 
            rmse(X(2,:),X_ukf(2,:));
            rmse(X(3,:),X_ukf(3,:)); 
            rmse(X(4,:),X_ukf(4,:))];
US(:,estimator + 1) = sum(sum(U*U'*Cr));
KS(:,estimator + 1) = k;
RMSES(:,:,estimator + 1) = [RMSE_EKF, RMSE_UKF];
switch (estimator)
    case 0
        title_str = sprintf('Trajectory using true state for controller input.\nRMSE EKF: e_x = %.3f, e_y %.3f, RMSE UKF: e_x = %.3f, e_y %.3f\nInput cost: %.3f, total timesteps: %d',...
            RMSE_EKF(1), RMSE_EKF(2),RMSE_UKF(1), RMSE_UKF(2),...
            US(1),KS(1));
    case 1
        title_str = sprintf('Trajectory using EKF estimate for controller input.\nRMSE EKF: e_x = %.3f, e_y %.3f\nInput cost: %.3f, total timesteps: %d',...
            RMSE_EKF(1), RMSE_EKF(2),US(2),KS(2));
    case 2
        title_str = sprintf('Trajectory using UKF estimate for controller input.\nRMSE UKF: e_x = %.3f, e_y %.3f\nInput cost: %.3f, total timesteps: %d',...
            RMSE_UKF(1), RMSE_UKF(2),US(3),KS(3));
end
title(title_str);
xlabel('X');ylabel('Y');
xlim(xl);ylim([-2 8]);
daspect([1 1 1]);
axis equal


%end % noplot
%fprintf('This is with %d packet loss.\n',p_*100);
%pause
end % estimator loop