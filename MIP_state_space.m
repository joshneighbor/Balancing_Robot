% Detailed write up for understanding MIP dynamics using state space
% Work done for Linear Systems Thoery course MAE 280A, Fall 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear

states = {'t' 'theta' 'theta_dot' 'phi' 'phi_dot'};
inputs = {'u'};
outputs = {'theta'; 'phi'};


A = [-50 0 -.2 0 -.2; 0 0 1 0 0; -3260.5 146.0849 0 0 0; 0 0 0 0 1; 4641.3 -94.5706 0 0 0];
B = [10;0;0;0;0];
C = [0 1 0 0 0; 0 0 0 1 0];
D = [0;0];

display('*')
Eigenvalues_A = eig(A) % Check eigen values
display('Examining the open-loop eigenvalues of A, we can tell that the')
display('system is unstable since we have a eigenvalue at 13. This value')
display('indicates how quickly the Mip will fall over.')

sys = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);
SYS = ss(A,B,C,D); %  convert to state-space model

display('*')
display('Compute the controllability matrix')
display('If it has full rank than is controllable.')
CO = ctrb(A,B); % Compute the controllability matrix
rank_CO = rank(CO) % if it has full rank than is controllable
display('The rank of CO is full thus is controllable.')
display('Unable to compute the contr or obs gramian because A is unstable, ie not Hurwitz.')
display('This means [A,B] is not reachable. Matlab uses Lyapunov equations: A*Wc + Wc*A'' + BB'' = 0')
display('all eigenvalues of A need negative real part, to have a unique')
display('solution for W_c')

%Wc = gram(SYS,'c') % omputes the controllability gramian
%Wo = gram(SYS,'o') % computes its observability gramian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In both cases, the state-space model SYS should be stable.
%     The gramians are computed by solving the Lyapunov equations:
%  
%        A*Wc + Wc*A' + BB' = 0  and   A'*Wo + Wo*A + C'C = 0
%          for continuous-time systems
%                 dx/dt = A x + B u  ,   y = C x + D u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = [.01; .25; .01; .01; .05]; % ['t' 'theta' 'theta_dot' 'phi' 'phi_dot']

display('*')
display('I tried several possible closed-loop pole positions for MiP.')
P1 = [-10.4, -9.5522+5.9047i, -9.5522-5.9047i, -44.1691, -10.5]
display('I kept the same stable poles from A and only adjusted the unstable poles.')
K1 = place(A,B,P1)
display('I tested these and this did not yeild the responce time I was looking for. Several more iterations')
display('of this process took place to get to poles that I liked')
display('I tried adjusting all poles so are close together:')
P = [-30 -30 -30 -30 -30]
K = acker(A,B,P)
display('Looking at input u = -Kx we note that max input = 12V for Mip')
display('These poles P are my favorite based on how quickly the system was brought back to steady state and')
display('based on the phi_dot and theta_dot measurements relating to a voltage is acceptable.')

input = K*x %must be less than 12v since max input is 12v.
if input > 12
    input
    return % our input is 12v so cannot use more.
end

AA = (A - B*K); % constant linear state-variable feedback (LSVF)
SYS_LSVF = ss(AA,B,C,D); %  convert to state-space model
Eigenvalues_AA = eig(AA) % Check eigen values

% t = 0:0.01:5;
% r =0.2*ones(size(t));
% [y,t,x]=lsim(SYS_LSVF,r,t);
% [AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
% title('Step Response with LSVF')

display('*')
display('Stabilization of MiP means in a physical context:')
display('The pole of at +13 reffered to theta being unstable, which means')
display('the Mip wanted to fall over. By moving the pole to the LHP by using')
display('LSVF I am making the physical system use its input from the battery to')
display('turn the motors which will now correct theta and try and force it back')
display('to its zero position, ie standing up right.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stabilization of MiP means in a physical context: I chose these poles so
% that the MiP uses its motors to correct or balance the body as it starts
% to fall over. I had to have an input torque which is large enough to
% counter balance. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observability and choice of sensor

display('*')
display('Mip is observable from y = [theta phi ] transpose')
ob_LSVF = obsv(SYS_LSVF); % creates ss with LSVF new A matrix
display('Create observability matrix of LSVF and check for full rank:')
observability_LSVF = rank(ob_LSVF) % checks to see if oberservable 
display('Since full rank it is observable.')
ob_SYS = obsv(SYS);
observability = rank(ob_SYS);
display('Now will check observabiltiy for single output theta:')
C_theta = [0 1 0 0 0; 0 0 0 0 0]; % single output theta
C_phi = [0 0 0 0 0; 0 0 0 1 0]; % single output phi
SYS_theta = ss(A,B,C_theta,D); % new ss w/ single output theta
SYS_phi = ss(A,B,C_phi,D); % new ss w/ single output phi
ob_SYS_theta = obsv(SYS_theta);
observability_SYS_theta = rank(ob_SYS_theta) % checks rank = 4
display('Rank 4 is not full and thus not observable from single output theta.')
display('Now I will check observabiltiy for single output phi:')
ob_SYS_phi = obsv(SYS_phi);
observability_SYS_phi = rank(ob_SYS_phi) % checks rank =5
display('Rank 5 is full and thus is observable from single output phi.')

display('*')
display('Unobservable subspace for MiP with only theta as a measurement is given by phi = x4')
display('This is shown by the NULL_space_theta:')
NULL_space_theta = null(ob_SYS_theta);
display('Note that nullspace of the observability matrix is the unobersable subspace.')
display('since phi is the unobservable component, by looking at the phi value of the system')
display('we do not gain any information about any of the other states.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We can see that for a single output theta is is not observable since it
% has a rank of 4 and the single output phi is observable since the rank is
% 5
% The unobservable subspace for the Mip with only theta as a measurement is
% given soley as phi = x4. This is shown by the NULL_space_theta. The null
% space of the observability matrix is the unobersable subspace. Which
% equals [0 0 0 -1 0].' 
% since phi is the unobservable component, by looking at the phi value of
% the sytem we do not gain any information about any of the other states. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('*')
display('Determine the observability properties using solely theta_dot in C matrix')
C_theta_dot = [0 0 1 0 0]; % changed C matrix for single output theta_dot
SYS_theta_dot = ss(A,B,C_theta_dot,D(1)); % new ss w/ single output theta_dot
ob_SYS_theta_dot = obsv(SYS_theta_dot);
observability_SYS_theta_dot = rank(ob_SYS_theta_dot) % checks rank = 4
display('Rank 4 is not full and thus not observable from single output theta_dot.')


display('*')
display('Design an observer for MiP')
LP = 2*[-4 -5 -6 -7 -8]; % choosing poles of L
L = place(A',C',LP)' % Placing poles of L
display('Controller comprised of your state feedback law and your observer:')
display('Obs_A =[A -B*K; L*C A-B*K-L*C]')
display('Obs_B = [B; B]')
display('Obs_C = [C -D*K; 0 -K]')
display('Obs_D = [D; I]')

Obs_A = [A -B*K; L*C A-B*K-L*C];
Obs_B = [B; B];
Obs_C = [C -D*K; zeros(size(K)) -K];
Obs_D = [D; 1];

Obs_C = [Obs_C; eye(5) zeros(5); zeros(5) eye(5)]; % adding in to view states
Obs_D = [Obs_D; zeros(10,1)];

states1 = {'t' 'theta' 'theta_dot' 'phi' 'phi_dot' 'e1' 'e2' 'e3' 'e4' 'e5'};
inputs1 = {'r'};
outputs1 = {'theta'; 'phi'};

SYS_Obs = ss(Obs_A, Obs_B, Obs_C, Obs_D); %  convert to state-space model
ob_observer = obsv(SYS_Obs); 

%Test of the response to a significant non-zero initial
%condition. Examining both the control signal, the states and the outputs.
[zz,t] = step(SYS_Obs);

figure
subplot(3,1,1)
plot(t,zz(:,4:8))
title('Step Response with Observer-Based State-Feedback Control')
legend('t', 'theta', 'thetadot', 'phi', 'phidot')
subplot(3,1,2)
plot(t,zz(:,1:2))
legend('theta', 'phi')
subplot(3,1,3)
plot(t,zz(:,3))
legend('voltage')

display('Setting initial conditions:')
display('set theta to 15 degrees which I feel is a reasonable angle for the Mip to start at')
%X_0 = [.01 .25 .01 .01 .05 0 0 0 0 0] % initial conditions
X_0 = [0 15*pi/180 0 0 0 0 0 0 0 0]
[Y,T,X] = initial(SYS_Obs,X_0);

figure
subplot(3,1,1)
plot(T,Y(:,4:8))
title('Setting Initial Conditions with Observer-Based LSVF')
legend('t', 'theta', 'thetadot', 'phi', 'phidot')
subplot(3,1,2)
plot(T,Y(:,1:2))
legend('theta', 'phi')
subplot(3,1,3);
plot(T,Y(:,3));
legend('voltage')

display('Checking our observer tracking')
figure
plot(T,Y(:,4));
hold on
plot(T,Y(:,9));
title('Tracking with Observer-Based LSVF')
legend('tau','tau_hat')
display('We see from the above graph that our observer is tracking the states.')
display('By changing are our L poles I can adjust the tracking.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bonus
display('*bonus')
display('I Repeated the above feedback controller design using the state and single output theta_dot')

A_bon = [-50 0 -.2 -.2; 0 0 1 0; -3260.5 146.0849 0 0; 4641.3 -94.5706 0 0];
B_bon = [10;0;0;0];
C_bon = [0 0 1 0; 0 0 0 0];
D_bon = [0;0];

P_bon = [-20 -21 -20.1 -20.2];
K_bon = place(A_bon,B_bon,P_bon);

LP_bon = 2*[-4 -5 -6 -7]; % choosing poles of L
L_bon = place(A_bon',C_bon',LP_bon)'; % Placing poles of L

SYS_bon = ss(A_bon,B_bon,C_bon,D_bon);
display('Check to see if this is observable')
ob_bon = obsv(SYS_bon); % creates ss with LSVF new A matrix
observability_bon = rank(ob_bon) % checks to see if oberservable
display('rank = 4 so we are observable')

Obs_A_bon = [A_bon -B_bon*K_bon; L_bon*C_bon A_bon-B_bon*K_bon-L_bon*C_bon];
Obs_B_bon = [B_bon; B_bon];
Obs_C_bon = [C_bon -D_bon*K_bon; zeros(size(K_bon)) -K_bon];
Obs_D_bon = [D_bon; 1];
% 
% Obs_C_bon = [Obs_C_bon; eye(4) zeros(4); zeros(4) eye(4)];
% Obs_D_bon = [Obs_D_bon; zeros(8,1)];
% 
SYS_Obs_bon = ss(Obs_A_bon,Obs_B_bon,Obs_C_bon,Obs_D_bon);
% 
% figure
% [zzz,tt] = step(SYS_Obs_bon);
% plot(tt,zzz(:,4:10))
% title('Bonus Step Response with Observer-Based State-Feedback Control')
% legend('t', 'theta', 'thetadot', 'phidot')

[zeroz,pole] = ss2tf(A_bon,B_bon,C_bon,D_bon); % turning statespace into transfer function
display('Plot of root locus of system')

figure
rlocus(zeroz(1,:),pole)
display('By looking at the root locus we see it does not possess the Parity Interlacing Property.')
display('We know it does not possess PIP since it does not have an even number of poles between each pair of zeros.')
display('There is no way to create a stable controller that woulb be give us a K to stabilize the system.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T6
display('* Task 6')
Aa = [0 0 0 0 0 0; 0 -50 0 -.2 0 -.2; 0 0 0 1 0 0; 100 -3260.5 146.0849 0 0 0; 0 0 0 0 0 1; 0 4641.3 -94.5706 0 0 0];
Bb = [0;10;0;0;0;0];
Cc = [0 0 1 0 0 0; 0 0 0 0 1 0];
Dd = [0;0];

SYS_sys = ss(Aa,Bb,Cc,Dd);

display('Show that the system remains observable from theta and phi, and also from theta alone')
Cp_theta = [0 0 1 0 0 0; 0 0 0 0 0 0]; % single output theta
Cp_phi = [0 0 0 0 0 0; 0 0 0 0 1 0]; % single output phi
SYSp_theta = ss(Aa,Bb,Cp_theta,Dd); % new ss w/ single output theta
SYSp_phi = ss(Aa,Bb,Cp_phi,Dd); % new ss w/ single output phi
ob_SYSp_theta = obsv(SYSp_theta);
display('Check to see if observable from single output theta:')
observability_SYSp_theta = rank(ob_SYSp_theta) % checks rank = 5
display('Rank 5 means it is not full rank and thus is not observable from single output theta.')
ob_SYSp_phi = obsv(SYSp_phi);
display('Check to see if observable from single output phi:')
observability_SYS_phi = rank(ob_SYSp_phi) % checks rank = 6
display('Rank 6 means it is full rank and thus is observable from single output phi.')

display('*')
NULL_space_thetap = null(ob_SYSp_theta);
display('Compute the controllability matrix')
display('If it has full rank than is controllable and thus reachable since we are in CT.')
CO_disturbance = ctrb(Aa,Bb); % Compute the controllability matrix
rank_CO_disturbance = rank(CO_disturbance) % if it has full rank than is controllable
display('System is unreachable because we have no control over the disturbance torque, none of the other state effect it.')
display('This means the system is not controllable.')

%System is unreachable because we have no control over the disturbance
%torque, none of the other state effect it. This means the system is not
%controllable

display('*')
display('New 6-dimensional feedback gain')
K = acker(A,B,[-44 -19 -19 -10 -10]);
display('K = [0 K]')
Kp = [0 K]
%Lp = [0 K]';
Lpp = (1/2)*[-40 -50 -49 -48 -47 -46];
Lp = place(Aa',Cc',Lpp)' % Placing poles of L

Obs_Ap = [Aa -Bb*Kp; Lp*Cc Aa-Bb*Kp-Lp*Cc];
Obs_Bp = [Bb; Bb];
Obs_Cp = [Cc -Dd*Kp; zeros(size(Kp)) -Kp];
Obs_Dp = [Dd; 1];

Obs_Cp = [Obs_Cp; eye(6) zeros(6); zeros(6) eye(6)];
Obs_Dp = [Obs_Dp; zeros(12,1)];

SYS_Obs_p = ss(Obs_Ap,Obs_Bp,Obs_Cp,Obs_Dp);


display('Setting initial conditions:')
X_00 = [1 0 0 0 0 0 0 0 0 0 0 0] % initial conditions, giving wind a value
[YY,TT,XX] = initial(SYS_Obs_p,X_00);

figure
subplot(3,1,1)
plot(TT,YY(:,4:9))
title('Disturbance: Setting Initial Conditions with Observer-Based LSVF')
legend('wind', 't', 'theta', 'thetadot', 'phi', 'phidot')
subplot(3,1,2)
plot(TT,YY(:,1:2))
legend('theta', 'phi')
subplot(3,1,3);
plot(TT,YY(:,3));
legend('voltage')

display('Disturbance: Checking our observer tracking')
figure
plot(TT,YY(:,5));
hold on
plot(TT,YY(:,11));
title('Disturbance: Tracking with Observer-Based LSVF')
legend('tau','tau_hat')
display('We see from the above graph that our observer is tracking the states.')
display('By changing are our L poles I can adjust the tracking.')

display('*')
display('We have to atjust the first value of K inorder to counter act the wind torque disturbance.')
display('Another way of solving this problem would be to use a feedforeward control since we know the disturbance.')

display('*')
display('Disturbance graphs below.')

% Kp = [100 K]
% Kp = place(Aa,Bb,[0 -21 -22 -23 -24 -25]);
% Obs_Ap = [Aa -Bb*Kp; Lp*Cc Aa-Bb*Kp-Lp*Cc];
% Obs_Bp = [Bb; Bb];
% Obs_Cp = [Cc -Dd*Kp; zeros(size(Kp)) -Kp];
% Obs_Dp = [Dd; 1];
% 
% Obs_Cp = [Obs_Cp; eye(6) zeros(6); zeros(6) eye(6)];
% Obs_Dp = [Obs_Dp; zeros(12,1)];

% SYS_Obs_p = ss(Obs_Ap,Obs_Bp,Obs_Cp,Obs_Dp);
% 
% [YY,TT,XX] = initial(SYS_Obs_p,X_00);
% 
% figure
% 
% plot(TT,YY(:,5))
% title('Disturbance: Adjusting K to counter disturbance from wind')
% legend('wind', 't', 'theta', 'thetadot', 'phi', 'phidot')

