function [] = particle_filter_L_C1(k,t_max,n,l,sig_l_s,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,L_av_s)
% in this code I used the theory on a particle filter from the paper: 
% "A tutorial on particle filters", by M. Speekenbrink
% in this we track a latent states of a stochastic process, by calculating
% the a sequence of posterior distributions over consecutive observations
% the process that will be estimated is a Gaussian random walk, showing the
% change over time of the envelope and facial movements made in speech
% In the generative model this case represents the audio track for the case
% of C=1
clear all;
clc;
% rng(5); %set random seed
% t_max=20; %set number of time steps in the process
% n=50; %set the number of observations made for each time step, this becomes the amount of particles for each step
% 
% c=1.0; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
% k=0.3; %setting the multiplication for the linear functions
% lrange=[-5 5]; %the range of l
% l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
% %parameters for the distributions for C=1
% sig_l_s=0.5;
% sig_lax=0.4;
% sig_lvx=0.5;
% sig_e_s=0.5;
% sig_e_ax=0.5;
% sig_e_vx=0.5;
% p_c=0.3;
% p_h=0.2;
% %% getting the distribution of the parameters from the counts in the generative model
% %sampling from the generative model for the case of C=1
% L=zeros(n,n);
% for i=1:1 %repeated draws to get different centers, in case we want ot estimate the average distribution for different centers
%     [i_la_plt,i_lv_plt,m,a_x,v_x,L_av_s,c_av_s,m_c,m_h]=gen_model(k,t_max,n,l,sig_l_s,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h);
%                                                                    
%     % setting the space parameters for the grid search
%     lrange=l;%linspace(-l,l,n);
%     mrange=linspace(1,4,4);
%     trange=linspace(1,-10,10);
%     counts_La = hist(i_la_plt,lrange);
%     counts_Lv = hist(i_lv_plt,lrange);
%     L=L+[counts_La'.*counts_Lv]; %the shape of this distribution should be similar to that of the likelihood
% end
 
%% the particle filter 
%setting arrays to store the important parameters
L_s=zeros(t_max,1); %to be estimated
L_xa=zeros(t_max,1); %noisy observations, want one observation for each time step
L_xv=zeros(t_max,1); %noisy observations, want one observation for each time step
phi_xav=zeros(t_max,n,2); %the samples from the distribution will be stored in here, need N samples per time step
%phi_xv=zeros(t_max,n); %the samples from the distribution will be stored in here, need N samples per time step
wa=zeros(t_max,n); %the weights of the different particles can be stored in here, one weight for each sample
wv=zeros(t_max,n); %the weights of the different particles can be stored in here, one weight for each sample
wav=zeros(t_max,n); %the weights of the different particles can be stored in here, one weight for each sample
w_storeav=zeros(t_max,n); %to store the before reweighting estimates
estav=zeros(t_max,1,2); %the array to store the estimates in, these will be formed by the sum over the particles, weighted by their respective weight
%wv=zeros(t_max,n); %the weights of the different particles can be stored in here, one weight for each sample
%w_storev=zeros(t_max,n); %to store the before reweighting estimates
%estv=zeros(t_max,1); %the array to store the estimates in, these will be formed by the sum over the particles, weighted by their respective weight
N_eff=zeros(t_max,1); %array to store the effective sample size
U=zeros(t_max,n); %storing a parameter for deciding when to resample

%initialising the parameters
L_s(:,1)=L_av_s*ones(1,t_max); %setting the actual stimulus locations

L_xa(1)=normrnd(L_s(1),sig_lax); %creating the parameter that relates the noisy percept to the stimulus
phi_xav(1,:,1)=normrnd(L_s(1),sig_lax,1,n);
L_xv(1)=normrnd(L_s(1),sig_lvx); %creating the parameter that relates the noisy percept to the stimulus
phi_xav(1,:,2)=normrnd(L_s(1),sig_lvx,1,n);

%now for the calculation of the weights we run into a small issue if we
%keep this continuous, as then the P(y|phi) would be 0 for every y (as the
%y is defined as a infinitly precise number in that case, thus the area
%under the curve would be 0). Thus we need to discritize the normal
%distribution of P(y|phi), so that the probability of a value is defined.
%It can then be calculated through the difference between the point just
%below and just above y in a cumulative distribution.
x_comp=linspace(-100,100,1000000); %making a large range of points, with small steps in between to evaluate the probabilities at

%the weights are now calculated over the distribution obtained after
%intgerating over L^{s}, which in the case of C=1 is no longer something
%that is a product of Gaussians. The value of the distribution can still be
%evaluated
w_ba=(normcdf(L_xa(1),x_comp(dsearchn(x_comp',phi_xav(1,:,1)')),sig_lax)); %obtain the cdf for lower estimate
w_aba=(normcdf(L_xa(1),x_comp(dsearchn(x_comp',phi_xav(1,:,1)')+ones(size(phi_xav(1,:,1),1))),sig_lax)); %obtain cdf for higher estimate
wa(1,:)=(w_ba-w_aba);%/(sum(w_ba-w_aba)); %obtain weights and normalise, note for t=0 no multiplication with previous weights
 
w_bv=(normcdf(L_xv(1),x_comp(dsearchn(x_comp',phi_xav(1,:,2)')),sig_lvx)); %obtain the cdf for lower estimate
w_abv=(normcdf(L_xv(1),x_comp(dsearchn(x_comp',phi_xav(1,:,2)')+ones(size(phi_xav(1,:,2),1))),sig_lvx)); %obtain cdf for higher estimate
wv(1,:)=(w_bv-w_abv);%/(sum(w_bv-w_abv)); %obtain weights and normalise, note for t=0 no multiplication with previous weights

wav(1,:)=wa(1,:).*wv(1,:);
wav(1,:)=wav(1,:)/sum(wav(1,:));

w_storeav(1,:)=wav(1,:);
% w_storeav(1,:,2)=wav(1,:,2);
%now we can obtain the particle filter estimates
estav(1,1,1)=sum(phi_xav(1,:,1).*wav(1,:));
estav(1,1,2)=sum(phi_xav(1,:,2).*wav(1,:));
for i=2:t_max
    L_xa(i)=normrnd(L_s(i-1),sig_lax); %creating the parameter that relates the noisy percept to the stimulus
    phi_xav(i,:,1)=normrnd(phi_xav(i-1,:,1),sig_lax,1,n);
    L_xv(i)=normrnd(L_s(i-1),sig_lvx); %creating the parameter that relates the noisy percept to the stimulus
    phi_xav(i,:,2)=normrnd(phi_xav(i-1,:,2),sig_lvx,1,n);
    
    %now for the calculation of the weights we run into a small issue if we
    %keep this continuous, as then the P(y|phi) would be 0 for every y (as the
    %y is defined as a infinitly precise number in that case, thus the area
    %under the curve would be 0). Thus we need to discritize the normal
    %distribution of P(y|phi), so that the probability of a value is defined.
    %It can then be calculated through the difference between the point just
    %below and just above y in a cumulative distribution.
    x_comp=linspace(-100,100,1000000); %making a large range of points, with small steps in between to evaluate the probabilities at
    
    w_ba=(normcdf(L_xa(i),x_comp(dsearchn(x_comp',phi_xav(i,:,1)')),sig_lax)); %obtain the cdf for lower estimate
    w_aba=(normcdf(L_xa(i),x_comp(dsearchn(x_comp',phi_xav(i,:,1)')+ones(size(phi_xav(i,:,1),1))),sig_lax)); %obtain cdf for higher estimate
    wa(i,:)=((w_ba-w_aba));%.*wa(i-1,:)); %obtain weights and normalise, note for t=0 no multiplication with previous weights
     
    w_bv=(normcdf(L_xv(i),x_comp(dsearchn(x_comp',phi_xav(i,:,2)')),sig_lvx)); %obtain the cdf for lower estimate
    w_abv=(normcdf(L_xv(i),x_comp(dsearchn(x_comp',phi_xav(i,:,2)')+ones(size(phi_xav(i,:,2),1))),sig_lvx)); %obtain cdf for higher estimate
    wv(i,:)=((w_bv-w_abv));%.*wv(i-1,:)); %obtain weights and normalise, note for t=0 no multiplication with previous weights
    
    wav(i,:)=(wa(i,:).*wv(i,:)).*wav(i-1,:); %need to multiply with previous values here and not in the seperate estimates, as it needs to be a combined value
    wav(i,:)=wav(i,:)/sum(wav(i,:));

    w_storeav(i,:)=wav(i,:);
%     w_storeav(i,:,1)=wav(i,:,1);
%     w_storeav(i,:,2)=wav(i,:,2);
    %now we can obtain the particle filter estimates
    estav(i,1,1)=sum(phi_xav(i,:,1).*wav(i,:));
    estav(i,1,2)=sum(phi_xav(i,:,2).*wav(i,:));

    %resampling step, using a so called systematic resampler
    N_effa(i)=1/(sum(wav(i,:).*wav(i,:))); %the effective weight of the samples
    if N_effa(i) <= c*n
        u=unifrnd(0,1,n,1);
        w_sum=cumsum(wav(i,:));
        wav(i,:)=1/n;
        for s=1:2
            for i_b=1:n 
                for r=2:n
                    U(i_b,:)=((i_b-1)/n)+u;
                    if w_sum(r-1)<=U(i_b,r) && U(i_b,r)<w_sum(r)
                        phi_xav(i,i_b,s)=phi_xav(i,r,s);
%                         phi_xav(i,i_b,2)=phi_xav(i,r,2);
                    end
                end
            end
        end
    end
%     N_eff(i)=1/(sum(wv(i,:).*wv(i,:))); %the effective weight of the samples
%     if N_eff(i) <= c*n
%         u=unifrnd(0,1,n,1);
%         w_sum=cumsum(wv(i,:));
%         wv(i,:)=1/n;
%         for i_b=1:n 
%             for r=2:n
%                 U(i_b,:)=((i_b-1)/n)+u;
%                 if w_sum(r-1)<=U(i_b,r) && U(i_b,r)<w_sum(r)
%                     phi_xv(i,i_b)=phi_xv(i,r);
%                 end
%             end
%         end
%     end
end

% %% plotting the trajectories
% figure(1)
% plot(L_s);
% hold on
% plot(L_xa);
% plot(L_xv);
% plot(estav(:,:,1));
% plot(estav(:,:,2));
% % legend("stimulus","percept a","percept v","estimate a","estimate v");
% 
% % legend("stimulus v","percept v","estimate v");
% figure(1)
% area=(500*w_storeav(5,:,1))+0.1;
% scatter(5,phi_xav(5,:,1),area)
% hold on
% area=(500*w_storeav(10,:,1))+0.1;
% scatter(10,phi_xav(10,:,1),area)
% hold on
% area=(500*w_storeav(15,:,1))+0.1;
% scatter(15,phi_xav(15,:,1),area)
% figure(2)
% area=(500*w_storev(5,:))+0.1;
% scatter(5,phi_xv(5,:),area)
% hold on
% area=(500*w_storev(10,:))+0.1;
% scatter(10,phi_xv(10,:),area)
% hold on
% area=(500*w_storev(15,:))+0.1;
% scatter(15,phi_xv(15,:),area)
end