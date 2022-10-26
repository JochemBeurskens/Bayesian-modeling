% in this code I used the theory on a particle filter from the paper: 
% "A tutorial on particle filters", by M. Speekenbrink
% in this we track a latent states of a stochastic process, by calculating
% the a sequence of posterior distributions over consecutive observations
% the process that will be estimated is a Gaussian random walk, showing the
% change over time of the envelope and facial movements made in speech
% In the generative model this case represents the audio track for the case
% of C=2

% rng(5); %set random seed
% setting parameters
k=0.3; %setting the multiplication for the linear functions
t_max=100; %setting the amount of time steps for the source activity
n=100; %number of steps in the given range
l=linspace(-5,5,n); %setting the range for the locations
%setting arrays to store the important parameters
sig_l_s=0.1;
sig_lax=0.1;
sig_lvx=1;
sig_e_s=0.5;
sig_e_ax=0.5;
sig_e_vx=0.5;
p_c=0.7;
p_h=0.8;

c=1.0; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
%% the particle filter 
%firstly the random walk, and obtain the distribution for the posterior of L
[dist_La,dist_Lv]=grid_search_post_L_C2(k,t_max,n,l,sig_l_s,sig_lax,sig_lvx,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h);

a_s=zeros(t_max,1); %to be estimated
a_x=zeros(t_max,1); %noisy observations, want one observation for each time step
phi_x=zeros(t_max,n); %the samples from the distribution will be stored in here, need N samples per time step
w=zeros(t_max,n); %the weights of the different particles can be stored in here, one weight for each sample
est=zeros(t_max,1); %the array to store the estimates in, these will be formed by the sum over the particles, weighted by their respective weight
N_eff=zeros(t_max,1); %array to store the effective sample size
U=zeros(t_max,n); %storing a parameter for deciding when to resample
sig_e_s=2; %setting the variance of the random walk
sig_e_ax=1; %setting the noise in the observations

%initialising the parameters
f_a_x=a_s(1)+normrnd(0,sig_e_ax); %creating the parameter that relates the noisy percept to the stimulus
a_x(1)=f_a_x;
phi_x(1,:)=normrnd(0,sig_e_s,1,n);
%now for the calculation of the weights we run into a small issue if we
%keep this continuous, as then the P(y|phi) would be 0 for every y (as the
%y is defined as a infinitly precise number in that case, thus the area
%under the curve would be 0). Thus we need to discritize the normal
%distribution of P(y|phi), so that the probability of a value is defined.
%It can then be calculated through the difference between the point just
%below and just above y in a cumulative distribution.
x_comp=linspace(-100,100,1000000); %making a large range of points, with small steps in between to evaluate the probabilities at
w_b=(normcdf(a_x(1),x_comp(dsearchn(x_comp',phi_x(1,:)')),sig_e_ax)); %obtain the cdf for lower estimate
w_ab=(normcdf(a_x(1),x_comp(dsearchn(x_comp',phi_x(1,:)')+ones(size(phi_x(1,:),1))),sig_e_ax)); %obtain cdf for higher estimate
w(1,:)=(w_b-w_ab)/(sum(w_b-w_ab)); %obtain weights and normalise, note for t=0 no multiplication with previous weights
c_av_s=[normrnd(0,sig_e_s)]; %drawing c_av^s from the normal distribution
for i=2:t_max
    %now we have to update the samples from the distribution, and their
    %respective weights for the next time step
    phi_x(i,:)=normrnd(phi_x(i-1),sig_e_s,1,n); %before was normrnd(phi_x(i-1,...) but then get totally different behavior. Thus this is slightly different from the paper, as there they seem to have made a small mistake!!!!
    w_b=(normcdf(a_x(i),x_comp(dsearchn(x_comp',phi_x(i,:)')),sig_e_ax)); %obtain the cdf for lower estimate
    w_ab=(normcdf(a_x(i),x_comp(dsearchn(x_comp',phi_x(i,:)')+ones(size(phi_x(i,:),1))),sig_e_ax)); %obtain cdf for higher estimate
    w(i,:)=((w_b-w_ab).*w(i-1,:))/(sum(((w_b-w_ab).*w(i-1,:)))); %obtain weights and normalise, note that now they have to be multiplied with the previous weigths
    
    %now we can obtain the particle filter estimates
    est(i)=sum(phi_x(i,:).*w(i,:));

    %resampling step, using a so called systematic resampler
    N_eff(i)=1/(sum(w(i,:).*w(i,:))); %the effective weight of the samples
    if N_eff(i) <= c*n
        u=unifrnd(0,1,n,1);
        w_sum=cumsum(w(i,:));
        w(i,:)=1/n;
        for i_b=1:n 
            for r=2:n
                U(i_b,:)=((i_b-1)/n)+u;
                if w_sum(r-1)<=U(i_b,r) && U(i_b,r)<w_sum(r)
                    phi_x(i,i_b)=phi_x(i,r);
                end
            end
        end
    end
end

%% plotting the trajectories
plot(a_s);
hold on
plot(a_x);
% plot(phi_x);
plot(est);
legend("stimulus","percept","estimate");