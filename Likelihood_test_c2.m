%for the likelihood calculation I will be using almost the same code, but
%as the very large arrays of a and v take too much space, I will not use
%them here. I will only calculate the parameters that decide where in these
%there are non-zero values, and then use this to make a frequentist
%estimation of the likelihood of a given draw

%% setting parameters

rng(0); %setting the random seed, so that the draws will have the same result each time, for now this is useful.

k=1.0; %setting the multiplication for the linear functions
t_max=1; %setting the amount of time steps for the source activity
n=5; %setting the amount of possible locations

lrange=[-1 1];
mrange=[-1 1];
crange=[-1 1];

sig_e_s=0.1;
sig_l_s=0.6;
sig_m_s=0.1;%setting noise parameters for the generating locations step. used for both C=1 and 2

sig_e_as=0.1;
sig_e_vs=0.1;
sig_l_as=0.6;
sig_m_as=0.6;
sig_l_vs=0.6;
sig_m_vs=0.6; %setting the additional noise parameters for the case of c=2

l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
m=linspace(mrange(1),mrange(2),n);
c=linspace(crange(1),crange(2),n);

sig_e_ax=0.1;
sig_e_vx=0.1;
sig_la=0.5;
sig_ma=0.5;
sig_lv=0.5;
sig_mv=0.5; %setting the noise parameters for the subject's observation

%% drawing the a_x and v_x for which the likelihoods are to be estimated
[f_a_x,f_v_x,i_l_x,i_m_x,i_la_x,i_ma_x,i_lv_x,i_mv_x]=Likelihood(k,t_max,n,c,l,m,sig_e_s,sig_l_s,sig_m_s,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);
%% Now  going over the entire space state to find estimates of the probabilities for the parameters in the generated a_x and v_x
%note that it is here useful to also count the instances that are the same
%as the values in the a_x and v_x

repeat=10000;
L=3;
M=2;
c_given=1;
[f_a_s_plt_stat,f_v_s_plt_stat,i_l_plt_stat,i_m_plt_stat,i_la_plt_stat,i_ma_plt_stat,i_lv_plt_stat,i_mv_plt_stat]=Likelihood_given_input_c1(k,t_max,n,c,l,m,L,M,c_given,sig_e_s,sig_l_s,sig_m_s,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);

f_v_sc2_store=zeros(repeat,t_max);
f_a_sc2_store=zeros(repeat,t_max);

i_lac2_store=zeros(repeat,t_max);
i_mac2_store=zeros(repeat,t_max);
i_lvc2_store=zeros(repeat,t_max);
i_mvc2_store=zeros(repeat,t_max);
for q=1:repeat
    %C=2, run model and plot
    [f_a_s_plt,f_v_s_plt,i_las_plt,i_mas_plt,i_lvs_plt,i_mvs_plt,i_la_plt,i_ma_plt,i_lv_plt,i_mv_plt]=Likelihood(k,t_max,n,c,l,m,sig_e_as,sig_e_vs,sig_l_as,sig_l_vs,sig_m_as,sig_m_vs,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);
    i_lac2_store(q,:)=i_la_plt;
    i_mac2_store(q,:)=i_ma_plt;
    i_lvc2_store(q,:)=i_lv_plt;
    i_mvc2_store(q,:)=i_mv_plt;
    f_v_sc2_store(q,:)=f_v_s_plt;
    f_a_sc2_store(q,:)=f_a_s_plt;
end
figure(3);
scatter3(i_la_store(:,1),i_ma_store(:,1),f_a_s_store(:,1));
hold on;
scatter3(i_lv_store(:,1),i_mv_store(:,1),f_v_s_store(:,1));
hold off;
%% now that
% [GC,GR] = groupcounts(i_la_store(:,1)); %this counts the amount of times a value shows up in a parameter, this can then be used to plot the distribution of said parameter

%make the a_x and v_x parameters that will be used in the estimation
[f_a_s_plt,f_v_s_plt,i_l_plt,i_m_plt,i_la_plt,i_ma_plt,i_lv_plt,i_mv_plt]=Likelihood(k,t_max,n,c,l,m,sig_e_s,sig_l_s,sig_m_s,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);
