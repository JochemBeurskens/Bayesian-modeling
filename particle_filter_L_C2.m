% in this code I used the theory on a particle filter from the paper: 
% "A tutorial on particle filters", by M. Speekenbrink
% in this we track a latent states of a stochastic process, by calculating
% the a sequence of posterior distributions over consecutive observations
% the process that will be estimated is a Gaussian random walk, showing the
% change over time of the envelope and facial movements made in speech
% In the generative model this case represents the audio track for the case
% of C=2
clear all;
clc;
% rng(5); %set random seed
t_max=20; %set number of time steps in the process
n=50; %set the number of observations made for each time step, this becomes the amount of particles for each step

c=1.0; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
k=0.3; %setting the multiplication for the linear functions
lrange=[-5 5]; %the range of l
l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
%parameters for the distributions for C=2
sig_l_vs=0.5;
sig_l_as=0.4;
sig_lvx=0.5;
sig_lax=0.4;
sig_e_s=0.5;
sig_e_ax=0.5;
sig_e_vx=0.5;
p_c=0.3;
p_h=0.2;
%% getting the distribution of the parameters from the counts in the generative model
%sampling from the generative model for the case of C=1
L=zeros(n,n);
for i=1:1 %repeated draws to get different centers, in case we want ot estimate the average distribution for different centers
    [i_la_plt,i_lv_plt,m,a_x,v_x,L_av_vs,L_av_as,c_av_s,m_c,m_h]=gen_model(k,t_max,n,l,sig_l_vs,sig_l_as,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h);

    % setting the space parameters for the grid search
    lrange=l;%linspace(-l,l,n);
    mrange=linspace(1,4,4);
    trange=linspace(1,-10,10);
    counts_La = hist(i_la_plt,lrange);
    counts_Lv = hist(i_lv_plt,lrange);
    L=L+[counts_La'.*counts_Lv]; %the shape of this distribution should be similar to that of the likelihood
end
 
%% the particle filter 
%setting arrays to store the important parameters
L_sa=zeros(t_max,1); %to be estimated
L_xa=zeros(t_max,1); %noisy observations, want one observation for each time step
L_sv=zeros(t_max,1); %to be estimated
L_xv=zeros(t_max,1); %noisy observations, want one observation for each time step
phi_xa=zeros(t_max,n); %the samples from the distribution will be stored in here, need N samples per time step
phi_xv=zeros(t_max,n); %the samples from the distribution will be stored in here, need N samples per time step
wa=zeros(t_max,n); %the weights of the different particles can be stored in here, one weight for each sample
w_storea=zeros(t_max,n); %to store the before reweighting estimates
esta=zeros(t_max,1); %the array to store the estimates in, these will be formed by the sum over the particles, weighted by their respective weight
wv=zeros(t_max,n); %the weights of the different particles can be stored in here, one weight for each sample
w_storev=zeros(t_max,n); %to store the before reweighting estimates
estv=zeros(t_max,1); %the array to store the estimates in, these will be formed by the sum over the particles, weighted by their respective weight
N_eff=zeros(t_max,1); %array to store the effective sample size
U=zeros(t_max,n); %storing a parameter for deciding when to resample

%initialising the parameters
L_sa=L_av_as*ones(1,t_max);
L_sv=L_av_vs*ones(1,t_max); %setting the actual stimulus locations

L_xa(1)=normrnd(L_sa(1),sig_lax); %creating the parameter that relates the noisy percept to the stimulus
phi_xa(1,:)=normrnd(L_sa(1),sig_lax,1,n);
L_xv(1)=normrnd(L_sv(1),sig_lvx); %creating the parameter that relates the noisy percept to the stimulus
phi_xv(1,:)=normrnd(L_sv(1),sig_lvx,1,n);

%now for the calculation of the weights we run into a small issue if we
%keep this continuous, as then the P(y|phi) would be 0 for every y (as the
%y is defined as a infinitly precise number in that case, thus the area
%under the curve would be 0). Thus we need to discritize the normal
%distribution of P(y|phi), so that the probability of a value is defined.
%It can then be calculated through the difference between the point just
%below and just above y in a cumulative distribution.
x_comp=linspace(-100,100,1000000); %making a large range of points, with small steps in between to evaluate the probabilities at

w_ba=(normcdf(L_xa(1),x_comp(dsearchn(x_comp',phi_xa(1,:)')),sig_lax)); %obtain the cdf for lower estimate
w_aba=(normcdf(L_xa(1),x_comp(dsearchn(x_comp',phi_xa(1,:)')+ones(size(phi_xa(1,:),1))),sig_lax)); %obtain cdf for higher estimate
wa(1,:)=(w_ba-w_aba)/(sum(w_ba-w_aba)); %obtain weights and normalise, note for t=0 no multiplication with previous weights

w_bv=(normcdf(L_xv(1),x_comp(dsearchn(x_comp',phi_xv(1,:)')),sig_lvx)); %obtain the cdf for lower estimate
w_abv=(normcdf(L_xv(1),x_comp(dsearchn(x_comp',phi_xv(1,:)')+ones(size(phi_xv(1,:),1))),sig_lvx)); %obtain cdf for higher estimate
wv(1,:)=(w_bv-w_abv)/(sum(w_bv-w_abv)); %obtain weights and normalise, note for t=0 no multiplication with previous weights

w_storea(1,:)=wa(1,:);
w_storev(1,:)=wv(1,:);
%now we can obtain the particle filter estimates
esta(1)=sum(phi_xa(1,:).*wa(1,:));
estv(1)=sum(phi_xv(1,:).*wv(1,:));
for i=2:t_max
    L_xa(i)=normrnd(L_sa(i-1),sig_lax); %creating the parameter that relates the noisy percept to the stimulus
    phi_xa(i,:)=normrnd(phi_xa(i-1,:),sig_lax,1,n);
    L_xv(i)=normrnd(L_sv(i-1),sig_lvx); %creating the parameter that relates the noisy percept to the stimulus
    phi_xv(i,:)=normrnd(phi_xv(i-1,:),sig_lvx,1,n);
    
    %now for the calculation of the weights we run into a small issue if we
    %keep this continuous, as then the P(y|phi) would be 0 for every y (as the
    %y is defined as a infinitly precise number in that case, thus the area
    %under the curve would be 0). Thus we need to discritize the normal
    %distribution of P(y|phi), so that the probability of a value is defined.
    %It can then be calculated through the difference between the point just
    %below and just above y in a cumulative distribution.
    x_comp=linspace(-100,100,1000000); %making a large range of points, with small steps in between to evaluate the probabilities at
    
    w_ba=(normcdf(L_xa(i),x_comp(dsearchn(x_comp',phi_xa(i,:)')),sig_lax)); %obtain the cdf for lower estimate
    w_aba=(normcdf(L_xa(i),x_comp(dsearchn(x_comp',phi_xa(i,:)')+ones(size(phi_xa(i,:),1))),sig_lax)); %obtain cdf for higher estimate
    wa(i,:)=((w_ba-w_aba).*wa(i-1,:))/((sum(w_ba-w_aba)).*wa(i-1,:)); %obtain weights and normalise, note for t=0 no multiplication with previous weights

    w_bv=(normcdf(L_xv(i),x_comp(dsearchn(x_comp',phi_xv(i,:)')),sig_lvx)); %obtain the cdf for lower estimate
    w_abv=(normcdf(L_xv(i),x_comp(dsearchn(x_comp',phi_xv(i,:)')+ones(size(phi_xv(i,:),1))),sig_lvx)); %obtain cdf for higher estimate
    wv(i,:)=((w_bv-w_abv).*wv(i-1,:))/((sum(w_bv-w_abv)).*wv(i-1,:)); %obtain weights and normalise, note for t=0 no multiplication with previous weights
    w_storea(i,:)=wa(i,:);
    w_storev(i,:)=wv(i,:);
    %now we can obtain the particle filter estimates
    esta(i)=sum(phi_xa(i,:).*wa(i,:));
    estv(i)=sum(phi_xv(i,:).*wv(i,:));

    %resampling step, using a so called systematic resampler
    N_effa(i)=1/(sum(wa(i,:).*wa(i,:))); %the effective weight of the samples
    if N_effa(i) <= c*n
        u=unifrnd(0,1,n,1);
        w_sum=cumsum(wa(i,:));
        wa(i,:)=1/n;
        for i_b=1:n 
            for r=2:n
                U(i_b,:)=((i_b-1)/n)+u;
                if w_sum(r-1)<=U(i_b,r) && U(i_b,r)<w_sum(r)
                    phi_xa(i,i_b)=phi_xa(i,r);
                end
            end
        end
    end
    N_eff(i)=1/(sum(wv(i,:).*wv(i,:))); %the effective weight of the samples
    if N_eff(i) <= c*n
        u=unifrnd(0,1,n,1);
        w_sum=cumsum(wv(i,:));
        wv(i,:)=1/n;
        for i_b=1:n 
            for r=2:n
                U(i_b,:)=((i_b-1)/n)+u;
                if w_sum(r-1)<=U(i_b,r) && U(i_b,r)<w_sum(r)
                    phi_xv(i,i_b)=phi_xv(i,r);
                end
            end
        end
    end
end

%% plotting the trajectories
figure(1)
plot(L_sa);
hold on
plot(L_xa);
plot(esta);
% legend("stimulus v","percept v","estimate v");

figure(2);
plot(L_sv);
hold on
plot(L_xv);
plot(estv);
% legend("stimulus v","percept v","estimate v");
figure(1)
area=(500*w_storea(5,:))+0.1;
scatter(5,phi_xa(5,:),area)
hold on
area=(500*w_storea(10,:))+0.1;
scatter(10,phi_xa(10,:),area)
hold on
area=(500*w_storea(15,:))+0.1;
scatter(15,phi_xa(15,:),area)
figure(2)
area=(500*w_storev(5,:))+0.1;
scatter(5,phi_xv(5,:),area)
hold on
area=(500*w_storev(10,:))+0.1;
scatter(10,phi_xv(10,:),area)
hold on
area=(500*w_storev(15,:))+0.1;
scatter(15,phi_xv(15,:),area)