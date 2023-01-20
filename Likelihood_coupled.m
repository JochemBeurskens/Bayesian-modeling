clear all;
clc;
rng(6); %setting the random number generator seed, for reproducible results

%% parameter specification
t_max=1; %set number of time steps in the process
n=20; %maybe not discritize, set the number of possible locations, the larger the more accurate, but also the more computation time is required
n_r=1; %number of stes of a_x or v_x obtained for each set of L_s, M_s and c_s
n_lm=1000; %number of sets of L^x and M^x that are considered for the integrals
c=0.8; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
k=1; %setting the multiplication for the linear functions
lrange=[-20 20]; %the range of l, actual locations are between -10.5 and 10.5
l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
sig_l_s=5;
sig_l_vs=5;
sig_l_as=sig_l_vs; %when assuming that the original location is drawn from one normal distribution, only one sigma_L^s is needed, same for c. They can now be seperate draws though, which for C=1 was not possible
sig_lax=8;
sig_lvx=5;
sig_e_s=1;
sig_e_ax=1;
sig_e_vx=1;
sig_m_s=1;
sig_m_as=1;
sig_m_vs=sig_m_as;
sig_m_ax=1;
sig_m_vx=0.5;
p_c=0.5; %if x param > p_c: m=1, else m=2
p_h=1;
mean_rw=2; %setting the mean of the random walk
%now we only need to obtain 1 set of a^x and v^x, and then evaluate it all
%from the perspective of the different L^s and M^s values in the grid
%% obtaining the a_x, v_x values
[i_la_plt,i_lv_plt,m_a,m_v,a_x,v_x,i_l,c_av_s,m_av_xa,m_av_xv,m_s]=gen_model_full(k,t_max,n,l,sig_l_s,sig_lvx,sig_lax,sig_m_s,sig_m_as,sig_m_vs,sig_m_ax,sig_m_vx,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,n_r,n_lm,mean_rw);
% surf(v_x)
%% now obtaining the measurement prob. distribution for the data above (Causal=1)
like_ax=zeros(n,2,t_max,n_lm);
like_vx=zeros(n,2,t_max,n_lm);

pc=normcdf(p_c,m_s,sig_m_ax);
p_c=[1-pc pc]; %prob for m=1 and m=2, note that if the m_x draw is larger than p_c it is put in m=1, thus 1-normcdf is the prob
for q=1:n_lm
    for o=1:n
        for p=1:2
            for t =1:t_max
                prob_Lvx=normpdf(l(i_l),l(o),sig_lvx); %N(L^x,L^s,sig)
                prob_vx(o,p,t,q)=prob_Lvx'.*normpdf(v_x(o,p,t,q),mean_rw,sig_e_s).*p_c(p);
                %prob_vx(i_lv_plt(t,o),m_v(t,o),t,o,p)=p_c(m_v(t,o))*prob_Lvx(i_lv_plt(t,o))*normpdf(v_x(i_lv_plt(t,o),m_v(t,o),t,o,p),mean_rw,sig_e_s);
    
                prob_Lax=normpdf(l(i_l),l(o),sig_lax); %N(L^x,L^s,sig)
                prob_ax(o,p,t,q)=prob_Lax'.*normpdf(a_x(o,p,t,q),mean_rw,sig_e_s).*p_c(p);
                %prob_ax(i_la_plt(t,o),m_a(t,o),t,o,p)=p_c(m_a(t,o))*prob_Lax(i_la_plt(t,o))*normpdf(a_x(i_la_plt(t,o),m_a(t,o),t,o,p),mean_rw,sig_e_s);
            end
        end 
    end
end

% probax now contains the probability values for 100 different sets of ax,
% thus we need to sum over these, note that the 5th dimension of the array
% will represent the different probability sets. If we sum over this
% dimension we should get a look at the average probability distribution
% for this value. Expect a peak at the x location in l and m, but less of a
% peak at the other locations (note that this off course depends on the
% discrepancy between the noise and the measured). 
%After this also a sum is performed over all the different options of L_x
%and M_x that were drawn, resulting in a general shape over the (L,M)-space

figure(1)
surf(sum(prob_ax(:,:,:,:),4))
xlabel('Stimulus Location'), ylabel('Stimulus Meaning'), zlabel('Likelihood(stimuli|a^{x})')
figure(2)
surf(sum(prob_vx(:,:,:,:),4))
xlabel('Stimulus Location'), ylabel('Stimulus Meaning'), zlabel('Likelihood(stimuli|a^{x})')