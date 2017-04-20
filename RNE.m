function RNE()
clc
clear all
close all

figure;
hold on;
h1 = plot3(0,0,0, 'b');   % The link
h2 = plot3(0,0,0, 'or');  % The joint/end-effector
h3 = plot3(0,0,0, 'xk');  % The center of mass
xlabel('x')
ylabel('y')
zlabel('z')
view([0 0 1]);
axis equal;
xlim([-2 2]);
ylim([-2 2]);
zlim([-0.5 0.5]);
legend('Link', 'Joint', 'Center of Mass')

% -------------------------------------------------------------------------
% Define the DH parameters
N_DOFS = 3;
dh.theta = [0 0 0];
dh.alpha = [0 0 0];
dh.offset = [0 0 0];
dh.d = [0 0 0];
dh.a = [0.5 0.5 0.5];
dh.type = ['r' 'r' 'r'];

% -------------------------------------------------------------------------
% Rigid body paramaters: inertia, mass, and cener of mass
rb.I =  repmat(eye(3), 1, 1, N_DOFS);
rb.m = [1 1 1];
% In standard DH, COM is mesured respect to its end-effector (using local 
% frame). When COM is [0 0 0]', it means the COM is located exactly at the 
% end-effector. Therefore, COM usually has negative values, which means the
% COM is behind the end-effector
rb.r = [-0.25 0 0; -0.25 0 0; -0.25 0 0]'; 

% -------------------------------------------------------------------------
% Relative transformation is respect to the EE of the previous joint
T_rel = repmat(eye(4),1, 1, N_DOFS);
T = repmat(eye(4),1, 1, N_DOFS);

% -------------------------------------------------------------------------
% Simulation time
ts = 0.001;
time_span = 0:ts:1;

% -------------------------------------------------------------------------
% Arbitrary trajectory as the inputs
qc_ = [pi/3*sin(2*pi*1*time_span)' pi/3*sin(2*pi*1*time_span)' pi/3*sin(2*pi*1*time_span)'];
qcdot_ = gradient(qc_', ts)';
qcddot_ = gradient(qcdot_', ts)';

% -------------------------------------------------------------------------
z0 = [0; 0; 1];

% Simulation starts here!
for k = 1 : length(time_span)
    
    % ---------------------------------------------------------------------
    % Kinematic simulation
    qc = qc_(k, :);
    qcdot = qcdot_(k, :);
    qcddot = qcddot_(k, :);
    
    % End-effector position, form the base which is located at r
    EE = zeros(3, N_DOFS+1);
    COM = zeros(3, N_DOFS);
    for i = 1 : 1 : N_DOFS
        T_rel = calc_rel_transformation(dh, qc);
        T = calc_transformation(dh, qc);
        EE(:,i+1) = T(1:3,4, i);
        COM(:,i) = EE(:,i+1) + T(1:3,1:3, i) * rb.r(:,i);        
    end
    
    % Draw the robot
    set(h1, 'XData', EE(1, :), 'YData', EE(2, :),'ZData', EE(3, :));
    set(h2, 'XData', EE(1, :), 'YData', EE(2, :),'ZData', EE(3, :));
    set(h3, 'XData', COM(1, :), 'YData', COM(2, :),'ZData', COM(3, :));
    drawnow;
    
    % ---------------------------------------------------------------------
    % Forward recursion
    for i = 1 : N_DOFS
        R = T_rel(1:3, 1:3, i);
        p = [dh.a(i); dh.d(i)*sin(dh.alpha(i)); dh.d(i)*cos(dh.alpha(i))];
             
        if i > 1
            w(:, i) =  R'*(w(:, i-1) + z0.*qcdot(i));
            wdot(:, i) = R'*(wdot(:, i-1) +  z0.*qcddot(i) + ...
                cross(w(:, i-1), z0.*qcdot(i)));
            vdot(:,i) = R'*vdot(:,i-1) + cross(wdot(:,i), p) + ...
                cross(w(:,i), cross(w(:,i),p));
        else
            w(:, i) =  R'*(z0.*qcdot(i));
            wdot(:, i) = R'*(z0.*qcddot(i));
            vdot(:,i) = cross(wdot(:,i), p) + ...
                cross(w(:,i), cross(w(:,i),p));
        end
    end
    
    % Dynamic simulation
    % Backward recursion
    for i = N_DOFS:-1:1
        p = [dh.a(i); dh.d(i)*sin(dh.alpha(i)); dh.d(i)*cos(dh.alpha(i))];
                
        vcdot = vdot(:,i) + cross(wdot(:,i),rb.r(:,i)) + ...
            cross(w(:,i),cross(w(:,i),rb.r(:,i)));
        
        F = rb.m(i)*vcdot;
        N = rb.I(:,:,i)*wdot(:,i)+cross(w(:,i),rb.I(:,:,i)*w(:,1));
        
        if i < N_DOFS
            R = T_rel(1:3, 1:3, i+1);
            n(:,i) = R*(n(:,i+1) + cross(R'*p, f(:,i+1))) + ...
                cross(rb.r(:,i)+p,F) + N;
            f(:,i) = R*f(:,i+1) + F;
        else
            n(:,i) = cross(rb.r(:,i)+p,F) + N;
            f(:,i) = F;
        end
        
        R = T_rel(1:3, 1:3, i);
        
        if dh.type(i) == 't'        
            Q(i,k) = f(:,i)'*R'*z0;
        elseif dh.type(i) == 'r'        
            Q(i,k) = n(:,i)'*R'*z0;
        end
    end
end

figure
plot(time_span, Q', '-b')
hold on

% -------------------------------------------------------------------------
% Using RVC Toolbox for comparison
L(1) = Revolute('d', 0, 'a', 0.5, 'alpha', 0, ...
    'm', rb.m(1), 'I', rb.I(:,:,1), 'r',  rb.r(:,1));
L(2) = Revolute('d', 0, 'a', 0.5, 'alpha', 0, ...
    'm', rb.m(2), 'I', rb.I(:,:,2), 'r', rb.r(:,2));
L(3) = Revolute('d', 0, 'a', 0.5, 'alpha', 0, ...
    'm', rb.m(3), 'I', rb.I(:,:,3), 'r', rb.r(:,3));

threelink = SerialLink(L, 'name', 'two link');
torque_by_rvc_toolbox = threelink.rne(qc_, qcdot_, qcddot_);
plot(time_span, torque_by_rvc_toolbox, '--r')

end

% -------------------------------------------------------------------------
function T = calc_rel_transformation(dh, qc)
% Transformation current link end-effector respect the previous link
% end-effector

N_DOFS = length(qc);

for i = 1 : 1 : N_DOFS
    if dh.type(i) == 'r'
        dh.theta(i) = qc(i);
    elseif dh.type(i) == 'p'
        dh.d(i) = qc(i);
    end
    
    ct = cos(dh.theta(i) + dh.offset(i));
    st = sin(dh.theta(i) + dh.offset(i));
    ca = cos(dh.alpha(i));
    sa = sin(dh.alpha(i));
    
    % Relative transformation
    T(:,:,i) = [ ct    -st*ca   st*sa     dh.a(i)*ct ; ...
        st    ct*ca    -ct*sa    dh.a(i)*st ; ...
        0     sa       ca        dh.d(i)    ; ...
        0     0        0         1         ];
end

end

% -------------------------------------------------------------------------
function  T = calc_transformation(dh, qc)
% Transformation respect to the global frame

N_DOFS = length(qc);
temp = eye(4);

for i = 1 : 1 : N_DOFS
    if dh.type(i) == 'r'
        dh.theta(i) = qc(i);
    elseif dh.type(i) == 'p'
        dh.d(i) = qc(i);
    end
    
    ct = cos(dh.theta(i) + dh.offset(i));
    st = sin(dh.theta(i) + dh.offset(i));
    ca = cos(dh.alpha(i));
    sa = sin(dh.alpha(i));
    
    temp = temp * [ ct    -st*ca   st*sa     dh.a(i)*ct ; ...
        st    ct*ca    -ct*sa    dh.a(i)*st ; ...
        0     sa       ca        dh.d(i)    ; ...
        0     0        0         1         ];
    T(:,:,i) = temp;
end

end

