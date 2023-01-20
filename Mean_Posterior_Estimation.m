clear all;
clc;

t_max=5; %set number of time steps in the process
n=1000; %maybe not discritize, set the number of possible locations, the larger the more accurate, but also the more computation time is required
n_pf = 100; %number of particles that are used in the filter
c=0.8; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
k=0.3; %setting the multiplication for the linear functions
lrange=[-50 50]; %the range of l, actual locations are between -10.5 and 10.5
l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
sig_l_s=3;
sig_l_vs=3;
sig_l_as=sig_l_vs; %when assuming that the original location is drawn from one normal distribution, only one sigma_L^s is needed, same for c. They can now be seperate draws though, which for C=1 was not possible
sig_lax=.3;
sig_lvx=.3;
sig_e_s=3;
sig_e_ax=3;
sig_e_vx=3;
p_c=0.3;
p_h=0.2;
%% running the generative model to get fake data
Causal_structure=randi([0,1],1,1); %this makes it such that the causal structure is decided upon, ie. C=1 or C=2, which decides which generative model should be used
if Causal_structure == 0
    %sampling from the generative model for the case of C=1
    [L_xa,L_xv,m,a_x,v_x,L_av_s,c_av_s,m_c,m_h]=gen_model(k,t_max,n,l,sig_l_s,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h);
elseif Causal_structure==1
    %sampling from the generative model for the case of C=2
    [L_xa,L_xv,m,a_x,v_x,L_av_vs,L_av_as,c_av_s,m_c,m_h]=gen_model(k,t_max,n,l,sig_l_vs,sig_l_as,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h);
end
%% obtaining the normalised probabilities and running the particle filters
[P_LC1,P_LC2]=Normalised_Probabilities_func(k,t_max,n,l,sig_l_s,sig_l_as,sig_l_vs,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,L_xa,L_xv);

[estav_c1] =    particle_filter_L_C1(c,k,t_max,n_pf,l,sig_l_s,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,L_xa,L_xv);
[esta,estv] =   particle_filter_L_C2(c,k,t_max,n_pf,l,sig_l_s,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,L_xa,L_xv);

%% now for each time step the location estimate can be calculated
est_La = P_LC1.*estav_c1(:,1,1)'+P_LC2.*esta';
est_Lv = P_LC1.*estav_c1(:,1,2)'+P_LC2.*estv';
%% plot
figure(1)
plot(L_xa);
hold on
plot(L_xv);
plot(estav_c1(:,:,1));
plot(estav_c1(:,:,2));
% plot(L_av_as*ones(20));
% plot(L_av_vs*ones(20));
plot(esta)
plot(estv)
plot(est_La)
plot(est_Lv)