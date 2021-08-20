%% Design the FIB-C(C)DMA Beamformer with Jacobi-Anger expansion
%-------------------------------------------------------------------------%
%%References
% [1]. G. Huang, J. Benesty and J. Chen, "On the Design of
% Frequency-Invariant Beampatterns With Uniform Circular Microphone
% Arrays," in IEEE/ACM Transactions on Audio, Speech, and Language
% Processing, vol. 25, no. 5, pp. 1140-1153, May 2017, doi:
% 10.1109/TASLP.2017.2689681.
%
% [2]. G. Huang, J. Chen and J. Benesty, "On the Design of Robust Steerable
% Frequency-Invariant Beampatterns with Concentric Circular Microphone
% Arrays," 2018 IEEE International Conference on Acoustics, Speech and
% Signal Processing (ICASSP), 2018, pp. 506-510, doi:
% 10.1109/ICASSP.2018.8461297.
%-------------------------------------------------------------------------%
% Written by Jeff Wang, Institute of Acoustics, University of Chinese
% Academy of Science
% Update Date: 19/Jul/2021
%-------------------------------------------------------------------------%
%%
clear;close all

%% scanning frequency and spatial angle
dF = 10;
fre_all = 10:dF:10e3;
c_sound = 340;                          %sound speech, m/s
k_all = 2*pi*fre_all/c_sound;           %wave number

%% array setting
load('4-4_3-2.2_aryloc.mat','ary_loc')  % M*3, 3-D location of each microphone
N_max = 1;                        %the cut-off order, should smaller than (M_min-1)/2
a = [1/3,2/3];                    %the coefficients of 1-order Hypercardioid
theta_tar = deg2rad(30);          %target direction, rad
M_used = [8, 8];                  %# Mics at each ring
r_circ = [2.2, 3]*1e-2;           %radius of rings
theta_st = deg2rad([0,0]);        %the start azimuth angle of the 1st microphone of each ring, rad

%% parameter setting
J = diag(1./1j.^(-N_max:N_max));
Gamma_mat = diag(exp(1j*(N_max:-1:-N_max)*theta_tar));
b = zeros(2*N_max+1, 1);
b(1:N_max) = 0.5*a(2:end);
b(N_max+1) = a(1);
b(N_max+2:end)  = flipud(b(1:N_max));
b_appro = J'*Gamma_mat'*b;

P_rings = length(M_used);
Phi = cell(1,P_rings);
WNG = zeros(size(fre_all));
theta_obe = deg2rad(-180:180);
Bp = zeros(size(fre_all,2),size(theta_obe,2)); % The target beampattern
Bp_syn = Bp;                                   % The synthesised beampattern

%% processing
for f=1:size(fre_all,2)
    
    k = k_all(f);
    for p=1:P_rings
        Phi{p} = besselj((-N_max:N_max).',k*r_circ(p)).*exp(-1j*(-N_max:N_max).'*(theta_st(p)+[0:M_used(p)-1]*2*pi/M_used(p)));
    end
    Phi_mat = conj(cat(2,Phi{:}));
    H_mat = Phi_mat'/(Phi_mat*Phi_mat')*b_appro;
    WNG(f) = -10*log10(sum(abs(H_mat.*H_mat)));
    Bp(f,:) = 20*log10(abs(exp(1j*theta_obe'*(-N_max:N_max))*Gamma_mat*b));
    d1 = exp(1j*k*r_circ(1)*cos(2*pi/M_used(1)*(0:M_used(1)-1)'-theta_obe));
    d2 = exp(1j*k*r_circ(2)*cos(2*pi/M_used(2)*(0:M_used(2)-1)'-theta_obe));
    d = [d1;d2];
    Bp_syn(f,:) = 20*log10(abs(H_mat'*d));
    
end
Bp(Bp<=-50) = -50;
Bp_syn = Bp_syn - max(Bp_syn,[],2);
Bp_syn(Bp_syn<=-50) = -50;

%% plot the beampattern & WNG
figure
plot(fre_all,WNG)
grid on
xticks([2:2:10]*1e3)
xticklabels([2:2:10])
ylim([-60,20])
yticks([-60:20:20])
ylabel('WNG (dB)','Interpreter','latex')
xlabel('$f$ (Hz)','Interpreter','latex')

figure;
[X,Y] = meshgrid(rad2deg(theta_obe), fre_all);
surf(X,Y,Bp,'linestyle','none') %%%ideal BP
view([0,0,1])
xticks([-180,-120,-60,0,60,120,180])
xlim([-180,180])
yticks([2:2:10]*1e3)
yticklabels([2:2:10])
colormap('jet')
colorbar('Ticks',[-50:10:0],...
    'TickLabels',{'-50dB','-40dB','-30dB','-20dB','-10dB','0dB'})
ylabel('$f/kHz$','Interpreter','latex','FontSize',12)
xlabel('$\theta$','Interpreter','latex','FontSize',12)
figure
surf(X,Y,Bp_syn,'linestyle','none')
view([0,0,1])
xticks([-180,-120,-60,0,60,120,180])
xlim([-180,180])
yticks([2:2:10]*1e3)
yticklabels([2:2:10])
colormap('jet')
colorbar('Ticks',[-50:10:0],...
    'TickLabels',{'-50dB','-40dB','-30dB','-20dB','-10dB','0dB'})
ylabel('$f/kHz$','Interpreter','latex','FontSize',12)
xlabel('$\theta$','Interpreter','latex','FontSize',12)

