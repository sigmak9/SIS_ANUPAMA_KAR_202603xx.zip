clear; clc; close all;

%% ==========================================================
%  AE 544 Project PP01
%  Anupama Kar 
%  Singularity & Ambiguity Study
%  3-2-1 Euler, Quaternion, CRP
%  Integrator Comparison: ode45, ode15s, ode4
% ===========================================================

%% ---------------- Simulation Setup ----------------

tspan = [0 60];

% Angular velocity designed to:
% - Drive pitch toward 90 deg (Euler singularity)
% - Produce large rotation (~180 deg) for CRP blow-up
w_fun = @(t) [0.2*sin(0.02*t);
              0.15;
              0.2*cos(0.02*t)];

% Initial Conditions
Euler0 = [0; 80*pi/180; 0];    % near singularity
q0     = [1; 0; 0; 0];         % identity quaternion

%% ---------------- Integrator Options ----------------

opt45 = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.5);
opt15 = odeset('RelTol',1e-6,'AbsTol',1e-8);

%% ==========================================================
%  EULER INTEGRATION
% ===========================================================

[t45_e, y45_e] = ode45(@(t,y) EulerODE(t,y,w_fun), tspan, Euler0, opt45);
[t15_e, y15_e] = ode15s(@(t,y) EulerODE(t,y,w_fun), tspan, Euler0, opt15);
[t4_e,  y4_e]  = ode4(@(t,y) EulerODE(t,y,w_fun), tspan, Euler0, 0.05);

%% ==========================================================
%  QUATERNION INTEGRATION
% ===========================================================

[t45_q, y45_q] = ode45(@(t,y) QuatODE(t,y,w_fun), tspan, q0, opt45);
[t15_q, y15_q] = ode15s(@(t,y) QuatODE(t,y,w_fun), tspan, q0, opt15);
[t4_q,  y4_q]  = ode4(@(t,y) QuatODE(t,y,w_fun), tspan, q0, 0.05);

% Normalize quaternions
y45_q = normalizeQuat(y45_q);
y15_q = normalizeQuat(y15_q);
y4_q  = normalizeQuat(y4_q);

%% ==========================================================
%  CLASSICAL RODRIGUES PARAMETERS (DIRECT KINEMATICS)
% ===========================================================

r0 = [0;0;0];   % corresponds to identity quaternion

[t45_r, y45_r] = ode45(@(t,y) CRPODE(t,y,w_fun), tspan, r0, opt45);
[t15_r, y15_r] = ode15s(@(t,y) CRPODE(t,y,w_fun), tspan, r0, opt15);
[t4_r,  y4_r]  = ode4(@(t,y) CRPODE(t,y,w_fun), tspan, r0, 0.05);

%% ==========================================================
%  PLOTS
% ===========================================================

%% Euler Angles vs Time
figure
plot(t45_e, y45_e(:,1)*180/pi,'r','LineWidth',2); hold on
plot(t45_e, y45_e(:,2)*180/pi,'g','LineWidth',2)
plot(t45_e, y45_e(:,3)*180/pi,'b','LineWidth',2)
xlabel('Time (s)')
yline(90,'k--','LineWidth',1.5)
yline(-90,'k--','LineWidth',1.5)
ylabel('Angle (deg)')
legend('\psi (Yaw)','\theta (Pitch)','\phi (Roll)')
title('3-2-1 Euler Angles vs Time (ode45)')
grid on

%% Euler Angle Comparison + Step Size

dt45 = diff(t45_e);
dt15 = diff(t15_e);
dt4 = diff(t4_e);

% ---------------- Step Size ----------------
figure
% Left axis: step size
yyaxis left
plot(t45_e(2:end), dt45, 'b', 'LineWidth',2); hold on
plot(t15_e(2:end), dt15, 'r--', 'LineWidth',2)
plot(t4_e(2:end), dt4, 'g:', 'LineWidth',2)
ylabel('Step Size (s)')
xlabel('Time (s)')
title('Integrator Step Size vs Time and Pitch (θ)')

% Right axis: pitch in degrees
yyaxis right
plot(t45_e, y45_e(:,2)*180/pi, 'k','LineWidth',1.5)  % pitch
ylabel('\theta (deg)')
yline(90,'k--'); 
yline(-90,'k--');

legend('ode45','ode15s','ode4','Pitch θ','Location','best')
grid on

%% Quaternion Components
figure
plot(t45_q,y45_q,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Quaternion Components')
legend('b0','b1','b2','b3')
title('Quaternion Components (Look for Sign Flip)')
grid on

%% Quaternion Norm
figure
plot(t45_q,vecnorm(y45_q,2,2),'LineWidth',2)
xlabel('Time (s)')
ylabel('||q||')
title('Quaternion Norm (Should Remain 1)')
grid on
ylim([0.5 1.5])
%% Quaternion Ambiguity Detection
dotTest = zeros(length(t45_q)-1,1);
for k = 2:length(t45_q)
    dotTest(k-1) = dot(y45_q(k,:), y45_q(k-1,:));
end

figure
plot(t45_q(2:end),dotTest,'LineWidth',2)
xlabel('Time (s)')
ylabel('q_k \cdot q_{k-1}')
title('Quaternion Ambiguity Detection (Flip if ~ -1)')
grid on

%% Rodrigues Blow-Up (Direct Integration)
figure
plot(t45_r,vecnorm(y45_r,2,2),'LineWidth',2)
xlabel('Time (s)')
ylabel('||r||')
title('Classical Rodrigues Blow-Up Near 180° (Direct Integration)')
grid on

% Right axis: pitch θ
yyaxis right
plot(t45_e, y45_e(:,2)*180/pi,'k','LineWidth',1.2)
yline(90,'k--'); yline(-90,'k--')
ylabel('\theta (deg)')
%% ==========================================================
%  ================= FUNCTIONS ===============================
% ===========================================================

function dydt = EulerODE(t,y,w_fun)
% y = [psi theta phi]'
psi = y(1);
theta = y(2);
phi = y(3);

w = w_fun(t);

s2 = sin(theta);
c2 = cos(theta);
s3 = sin(phi);
c3 = cos(phi);

B = 1/c2*[0 s3 c3;
          0 c2*c3 -c2*s3;
          c2 s2*s3 s2*c3];

dydt = B*w;

end

% ------------------------------------------------------------

function dqdt = QuatODE(t,q,w_fun)
% q = [b0 b1 b2 b3]'
w = w_fun(t);

Omega = [ 0    -w(1) -w(2) -w(3);
          w(1)  0     w(3) -w(2);
          w(2) -w(3)  0     w(1);
          w(3)  w(2) -w(1)  0];

dqdt = 0.5 * Omega * q;

end
% ------------------------------------------------------------

function drdt = CRPODE(t,r,w_fun)

w = w_fun(t);

r1 = r(1);
r2 = r(2);
r3 = r(3);

r_tilde = [  0   -r3   r2;
             r3    0   -r1;
            -r2   r1    0];

drdt = 0.5 * ( eye(3) + r_tilde + r*r' ) * w;

end
% ------------------------------------------------------------

function qn = normalizeQuat(q)
for i = 1:size(q,1)
    qn(i,:) = q(i,:) / norm(q(i,:));
end
end

% ------------------------------------------------------------

function [t,y] = ode4(odefun,tspan,y0,h)

t = (tspan(1):h:tspan(2))';
y = zeros(length(t),length(y0));
y(1,:) = y0';

for i = 1:length(t)-1
    k1 = odefun(t(i),y(i,:)') ;
    k2 = odefun(t(i)+h/2,y(i,:)'+h/2*k1);
    k3 = odefun(t(i)+h/2,y(i,:)'+h/2*k2);
    k4 = odefun(t(i)+h,y(i,:)'+h*k3);
    
    y(i+1,:) = y(i,:) + h/6*(k1'+2*k2'+2*k3'+k4');
end
end

%% ==========================================================
%  3D SPACECRAFT ANIMATION
% ==========================================================
% Using Euler angles (yaw-pitch-roll) from ode45 integration

% Load spacecraft image (billboard style)
img = imread('airplane.jpg'); % replace with your image file

% Flat rectangular surface for spacecraft
[X,Y] = meshgrid(linspace(-1,1,2), linspace(-0.5,0.5,2));
Z = zeros(size(X));

figure('Color','w')
axis equal
axis([-2 2 -2 2 -2 2])
grid on
view(3)
xlabel('X'); ylabel('Y'); zlabel('Z')
hold on

% Optional GIF filename
filename = 'Spacecraft_Attitude.gif';

for k = 1:length(t45_e)
   
    cla
    
    % Current Euler angles
    psi   = y45_e(k,1);  % yaw
    theta = y45_e(k,2);  % pitch
    phi   = y45_e(k,3);  % roll
    
    % 3-2-1 rotation
    Rz = [cos(psi) -sin(psi) 0;
          sin(psi)  cos(psi) 0;
          0         0        1];
      
    Ry = [cos(theta) 0 sin(theta);
          0          1 0;
         -sin(theta) 0 cos(theta)];
     
    Rx = [1 0          0;
          0 cos(phi)  -sin(phi);
          0 sin(phi)   cos(phi)];
    
    R = Rz*Ry*Rx;
    
    % Rotate spacecraft surface
    pts = R * [X(:)'; Y(:)'; Z(:)'];
    Xr = reshape(pts(1,:), size(X));
    Yr = reshape(pts(2,:), size(Y));
    Zr = reshape(pts(3,:), size(Z));
    
    % Plot textured spacecraft
    surface(Xr,Yr,Zr,img,'FaceColor','texturemap','EdgeColor','none')
    
    % Body axes
    quiver3(0,0,0,R(1,1),R(2,1),R(3,1),1,'r','LineWidth',2)
    quiver3(0,0,0,R(1,2),R(2,2),R(3,2),1,'g','LineWidth',2)
    quiver3(0,0,0,R(1,3),R(2,3),R(3,3),1,'b','LineWidth',2)
    
    % Display time and pitch
    text(-1.8,1.8,1.8, sprintf('t = %.2f s', t45_e(k)), 'FontSize',12,'FontWeight','bold')
    text(-1.8,1.8,1.5, sprintf('\\theta = %.2f deg', theta*180/pi), 'FontSize',12,'FontWeight','bold')
    
    axis equal
    axis([-2 2 -2 2 -2 2])
    view(35,20)
    drawnow
    
    % Capture frame for GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if k == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
    
end