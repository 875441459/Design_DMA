%% Design the UCCA FIB Beampattern using the Spherical Harmonic Transform
%References
% [1]. ZHAO, Xudong, et al. On the Design of 3D Steerable Beamformers with
% Uniform Concentric Circular Microphone Arrays. IEEE/ACM Transactions on
% Audio, Speech, and Language Processing, 2021.
%%
clear;close all;
%% Scanning frequency range
fs = 16e3;
% num_fft = 1024;
% dF = fs/num_fft;
dF = 100;
fre_stt = dF;
fre_stp = fs/2;
fre_all = fre_stt:dF:fre_stp;
fre_bins = length(fre_all);
c_sound = 340;
k_all = 2*pi*fre_all/c_sound;           %wave number

%% Array setting
L = [7.5 5 3];
phi_tar = 120; %the azimuth angle
theta_tar = 60;  % the elevation angle
sc_loc = [phi_tar theta_tar];%%%%%%%%%%%%%%%%%%%%May need modifiyd
phi_tar = deg2rad(phi_tar);
theta_tar = deg2rad(theta_tar);
%%%%%%%%%%%%%%%%the array setting may need modified%%%%%%%%%%%%%%%%%%%%
% %UCCA-IV
r_circ = [0 1.5, 2.5]*1e-2;       %radius of rings
M_mics = [1 7 7];

P_rings = length(r_circ);
type_ary = 'circular';
par_ary.M = M_mics;
par_ary.r = r_circ;
load('1_7_7_UCCA.mat')
par_ary.ary_loc = ary_loc;
par_ary.Sv_tar = Sv_tar;

%% BP setting
N_max = 2; %the truncation order % # constraints=(N_max+1)*(N_max+2)/2;
[h_BP,BP_dB,WNG,DF]= Design_FIB_SHBCCMAv1(par_ary,N_max, sc_loc, k_all);


