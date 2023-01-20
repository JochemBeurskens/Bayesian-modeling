clear all;
clc;

t_max=1; %set number of time steps in the process
n=100; %maybe not discritize, set the number of possible locations, the larger the more accurate, but also the more computation time is required
n_pf = 500; %number of particles that are used in the filter
c=0.8; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
k=0.3; %setting the multiplication for the linear functions
lrange=[-20 20]; %the range of l, actual locations are between -10.5 and 10.5
l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
sig_l_s=20;
sig_l_vs=20;
sig_l_as=sig_l_vs; %when assuming that the original location is drawn from one normal distribution, only one sigma_L^s is needed, same for c. They can now be seperate draws though, which for C=1 was not possible
sig_lax=8;
sig_lvx=15;
sig_e_s=3;
sig_e_ax=3;
sig_e_vx=3;
p_c=0.3;
p_h=0.2;
%% now checking that the sum equals 1
% generating a space to check the probilities over, these represent the
% observations x_a and x_v
L_xar=linspace(-50,50,n);
L_xvr=linspace(-50,50,n);
P_LC1r=zeros(length(L_xar),length(L_xvr));
P_LC2r=zeros(length(L_xar),length(L_xvr));
% using the set parameters to see whether the probabilities indeed add up
% to 1 in the whole range. See the plots for one specific L_Xa or L_Xv
for i=1:length(L_xar)
    L_xa=L_xar(i);
    for o=1:length(L_xvr)
        L_xv=L_xvr(o);
        [P_LC1,P_LC2]=Normalised_Probabilities_multimodal_func(k,t_max,n,l,sig_l_s,sig_l_as,sig_l_vs,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,L_xa,L_xv);
 
%         [P_LC1,P_LC2]=Normalised_Probabilities_func(k,t_max,n,l,sig_l_s,sig_l_as,sig_l_vs,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,L_xa,L_xv);
        P_LC1r(i,o)=P_LC1;
        P_LC2r(i,o)=P_LC2;
    end
end
%% plotting
figure(1)
plot(P_LC2r(50,:))
hold on
plot(P_LC1r(50,:))
legend("C=2","C=1")
figure(2)
plot(P_LC2r(:,40))
hold on
plot(P_LC1r(:,40))
legend("C=2","C=1")
figure(3)
surf(L_xar,L_xvr,P_LC1r,'FaceColor',[0,158/255,115/255]);
hold on; 
surf(L_xar,L_xvr,P_LC2r,'FaceColor',[230/255,159/255,0]);
