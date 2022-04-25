%% Setting the model parameters
% rng(0); %setting the random seed, so that the draws will have the same result each time, for now this is useful.

k=1.0; %setting the multiplication for the linear functions
t_max=3; %setting the amount of time steps for the source activity
n=11; %setting the amount of possible locations

P_C=0.5;% setting the probability of selecting C=1 or C=2
P_av_given_LMf=zeros(n,n,n,t_max); %make an array to store the probabilities in, this is an n (range of L) by n (range of M) by n (range of c or f) by t array (so that the different times can be stored as well, but designed s.t. the first entry contains the first prior)
P_a_given_LMf_c2=zeros(n,n,n,t_max);%same but for the C=2 and the auditory percept
P_v_given_LMf_c2=zeros(n,n,n,t_max);%same but for the C=2 and the visual percept

repeat_l=n;
repeat_m=n;
repeat=1;
Percept_L_av=zeros(n,n,repeat); %make arrays to store the percepts for all the repeated trials, note that there will be n,n,repeat entries (n draws from the space of l, for each draw of l 5 draws from m follow, then 10000 times repeated draws for those positions in m,l space)
Percept_M_av=zeros(n,n,repeat);

lrange=[-5 5];
mrange=[-5 5];
crange=[-0.5 0.5];

l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
m=linspace(mrange(1),mrange(2),n);
c=linspace(crange(1),crange(2),n);
lsteps=linspace(-5.5,5.5,2*n); %setting the range of bins for locations and meanings
msteps=linspace(-5.5,5.5,2*n);

lsteps_plt=linspace(-5.5,5.5,2*n-1); %setting the range of bins for locations and meanings
msteps_plt=linspace(-5.5,5.5,2*n-1);

L_av_xa=zeros(2*n-1,2*n-1,t_max);
L_av_xv=zeros(2*n-1,2*n-1,t_max);
M_av_xa=zeros(2*n-1,2*n-1,t_max);
M_av_xv=zeros(2*n-1,2*n-1,t_max);

P_av_a=zeros(2*n-1,2*n-1,t_max);
P_av_v=zeros(2*n-1,2*n-1,t_max);
%parameters for C=1
sig_e_s=0.1;
sig_l_s=0.6;
sig_m_s=0.3;%setting noise parameters for the generating locations step. used for both C=1 and 2

sig_e_ax=0.1;
sig_e_vx=0.1;
sig_la=0.5;
sig_ma=0.5;
sig_lv=0.25;
sig_mv=0.25; %setting the noise parameters for the subject's observation

%parameters for C=2
sig_e_s_c2=0.1;
sig_l_s_c2=0.6;
sig_m_s_c2=0.3;%setting noise parameters for the generating locations step. used for both C=1 and 2

sig_e_as_c2=0.1;
sig_e_vs_c2=0.1;
sig_l_as_c2=0.5;
sig_m_as_c2=0.5;
sig_l_vs_c2=0.25;
sig_m_vs_c2=0.25; %setting the additional noise parameters for the case of c=2

sig_e_ax_c2=0.1;
sig_e_vx_c2=0.1;
sig_la_c2=0.75;
sig_ma_c2=0.5;
sig_lv_c2=0.25;
sig_mv_c2=1.0; %setting the noise parameters for the subject's observation

%% obtain the prior distributions
for ls=1:length(lsteps)-1
    L(ls)=(normcdf(lsteps(ls+1),0,sig_l_s)-normcdf(lsteps(ls),0,sig_l_s));
    M(ls)=(normcdf(lsteps(ls+1),0,sig_m_s)-normcdf(lsteps(ls),0,sig_m_s));
end
%% Performing the grid search
store_percept_av=zeros(n,n,t_max);
store_percept_a=zeros(n,n,t_max);
store_percept_v=zeros(n,n,t_max);
for t=1:t_max
    store_a=zeros(n,n);
    store_v=zeros(n,n);
    store_a_c2=zeros(n,n);
    store_v_c2=zeros(n,n);
   
    li=1;
    %C=1
    for li=1:repeat_l
        mi=1;
        for mi=1:repeat_m    
            f_v_s_store=zeros(repeat,t_max);
            f_a_s_store=zeros(repeat,t_max);
            L_av_s=l(li);
            M_av_s=m(mi);
            %C=1, run model and plot
            for ls =1:length(lsteps)-1
                for ms =1:length(msteps)-1
                    L_av_xa(ls,ms,t)=L_av_xa(ls,ms,t)+(normcdf(lsteps(ls+1),L_av_s,sig_la)-normcdf(lsteps(ls),L_av_s,sig_la));
                    L_av_xv(ls,ms,t)=L_av_xa(ls,ms,t)+(normcdf(lsteps(ls+1),L_av_s,sig_lv)-normcdf(lsteps(ls),L_av_s,sig_lv));
                    M_av_xa(ls,ms,t)=M_av_xa(ls,ms,t)+(normcdf(msteps(ms+1),M_av_s,sig_ma)-normcdf(msteps(ms),M_av_s,sig_ma));
                    M_av_xv(ls,ms,t)=M_av_xa(ls,ms,t)+(normcdf(msteps(ms+1),M_av_s,sig_mv)-normcdf(msteps(ms),M_av_s,sig_mv));
                end
            end
            mi=mi+1;
        end
        li=li+1;
    end
    if t==1
        P_av_a(:,:,t)=P_C*(L_av_xa(:,:,t)*+M_av_xa(:,:,t));
        P_av_v(:,:,t)=P_C*(L_av_xv(:,:,t)+M_av_xv(:,:,t));
    else 
        P_av_a(:,:,t)=P_av_v(:,:,t-1)*P_C*(L_av_xa(:,:,t)+M_av_xa(:,:,t));
        P_av_v(:,:,t)=P_av_v(:,:,t-1)*P_C*(L_av_xv(:,:,t)+M_av_xv(:,:,t));
    end
end

%% Plotting
P_av_a=L.*P_av_a.*M';%multiplication with the prior distributions
P_av_v=L.*P_av_v.*M';
[x,y]=meshgrid(lsteps_plt,msteps_plt);
%first the results from the update equation
%C=1: first the audio-visual
figure(1)
contour(x,y,P_av_a(:,:,1),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(2)
contour(x,y,P_av_a(:,:,2),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(3)
contour(x,y,P_av_a(:,:,3),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
