%% Setting the model parameters
clear all;
% rng(0); %setting the random seed, so that the draws will have the same result each time, for now this is useful.

k=1.0; %setting the multiplication for the linear functions
t_max=3; %setting the amount of time steps for the source activity
n=21; %setting the amount of possible locations

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
lsteps=linspace(-5,5,n); %setting the range of bins for locations and meanings
msteps=linspace(-5,5,n);
%when I used normcdf I needed to distinguish between the lsteps and the
%ones used for plotting. Using normpdf it isn't autmatically normalised
%anymore but there is no longer a need for these to be seperate.

lsteps_plt=linspace(-5,5,n); %setting the range of bins for locations and meanings
msteps_plt=linspace(-5,5,n);

L_av=zeros(n,n,t_max);
L_av_xv=zeros(n,n,t_max);
M_av=zeros(n,n,t_max);
M_av_xv=zeros(n,n,t_max);
C_L=zeros(n,n); %a parameter to fill with the manual check
C_M=zeros(n,n);
I_L=zeros(n,n); %a parameter to fill with the manual check
I_M=zeros(n,n);

L_av_xa_c2=zeros(n,n,t_max);
L_av_xv_c2=zeros(n,n,t_max);
M_av_xa_c2=zeros(n,n,t_max);
M_av_xv_c2=zeros(n,n,t_max);

P_av_a=zeros(n,n,t_max);
P_av_v=zeros(n,n,t_max);

P_av_a_c2=zeros(n,n,t_max);
P_av_v_c2=zeros(n,n,t_max);
%parameters for C=1
sig_e_s=0.1;
sig_l_s=0.6;
sig_m_s=0.3;%setting noise parameters for the generating locations step. used for both C=1 and 2

sig_e_ax=0.1;
sig_e_vx=0.1;
sig_la=0.5;
sig_ma=0.5;
sig_lv=0.5;
sig_mv=0.5; %setting the noise parameters for the subject's observation

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

levels=linspace(0,1,15);
%% obtain the prior distributions
% for ls=1:length(lsteps)-1
%     L(ls)=(normcdf(lsteps(ls+1),0,sig_l_s)-normcdf(lsteps(ls),0,sig_l_s));%*(normcdf(lsteps(ls+1),L_av_s,sig_la)-normcdf(lsteps(ls),L_av_s,sig_la));
%     M(ls)=(normcdf(lsteps(ls+1),0,sig_m_s)-normcdf(lsteps(ls),0,sig_m_s));%*(normcdf(msteps(ls+1),M_av_s,sig_ma)-normcdf(msteps(ls),M_av_s,sig_ma));
%     L_a_c2(ls)=(normcdf(lsteps(ls+1),-2,sig_l_as_c2)-normcdf(lsteps(ls),-2,sig_l_as_c2));
%     M_a_c2(ls)=(normcdf(lsteps(ls+1),-2,sig_m_as_c2)-normcdf(lsteps(ls),-2,sig_m_as_c2));
%     L_v_c2(ls)=(normcdf(lsteps(ls+1),2,sig_l_vs_c2)-normcdf(lsteps(ls),2,sig_l_vs_c2));
%     M_v_c2(ls)=(normcdf(lsteps(ls+1),2,sig_m_vs_c2)-normcdf(lsteps(ls),2,sig_m_vs_c2));
% end
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
    for li=1:repeat_l %looping over all location grid points, think of these as the points vx and ax
        mi=1;
        for mi=1:repeat_m %looping over all semantic meaning grid points
            L_av_x=l(li); 
            M_av_x=m(mi); %placing the centre of the stimuli at each of the grid points
            %C=1, run model and plot
            for ls =1:length(lsteps) %now for each centre of the stimuli the range of observations is looked at
                for ms =1:length(msteps)
                    % Making the P(a|L_s,M_s)=
                    % sum(N(L_xa;mu_La,sig_La)*N(L_xv;mu_Lv,sig_Lv)*prior(L_s))
                    L_av(ms,ls,t)=L_av(ms,ls,t)+normpdf(L_av_x,lsteps(ms),sig_la)*normpdf(L_av_x,lsteps(ls),sig_lv)*normpdf(L_av_x,0,sig_l_s);
                    % same but for M
                    M_av(ms,ls,t)=M_av(ms,ls,t)+normpdf(M_av_x,msteps(ms),sig_ma)*normpdf(M_av_x,msteps(ls),sig_mv)*normpdf(M_av_x,0,sig_m_s);

                end
            end
        end
    end
    if t==1
        %Multiplying the results for L and M with each other to get the
        %probabilities on a grid.
        P_av(:,:,t)=P_C*((L_av(:,:,t).*M_av(:,:,t)));
        
    else 
        P_av(:,:,t)=P_av(:,:,t-1).*(P_C*(L_av(:,:,t).*M_av(:,:,t)));
        
    end

    %C=2
    for li_a=1:repeat_l %looping over all grid points
        for mi_a=1:repeat_m %looping over all meaning grid points
            %li_v=1;
            for li_v=1:repeat_l %looping over all grid points
                %mi_v=1;
                for mi_v=1:repeat_m
                    L_av_xv=l(li_v); 
                    M_av_xv=m(mi_v);
                    L_av_xa=l(li_a); 
                    M_av_xa=m(mi_a); %placing the centre of the noisy observation at each of the grid points
                    for ls =1:length(lsteps) %now for each centre of the stimuli the range of observations is looked at
                        for ms =1:length(msteps)
                            L_av_xa_c2(ms,ls,t)=L_av_xa_c2(ms,ls,t)+normpdf(L_av_xa,lsteps(ms),sig_la_c2)*normpdf(L_av_xa,0,sig_l_as_c2);
                            L_av_xv_c2(ms,ls,t)=L_av_xv_c2(ms,ls,t)+normpdf(L_av_xv,lsteps(ls),sig_lv_c2)*normpdf(L_av_xv,0,sig_l_vs_c2);
                            M_av_xa_c2(ms,ls,t)=M_av_xa_c2(ms,ls,t)+normpdf(M_av_xa,msteps(ms),sig_ma_c2)*normpdf(M_av_xa,0,sig_m_as_c2);
                            M_av_xv_c2(ms,ls,t)=M_av_xv_c2(ms,ls,t)+normpdf(M_av_xv,msteps(ls),sig_mv_c2)*normpdf(M_av_xv,0,sig_m_vs_c2);
                        end
                    end
                end
            end
        end
    end
    if t==1
        P_av_a_c2(:,:,t)=P_C*((L_av_xa_c2(:,:,t).*M_av_xa_c2(:,:,t)));
        P_av_v_c2(:,:,t)=P_C*((L_av_xv_c2(:,:,t).*M_av_xv_c2(:,:,t)));
    
    else 
        P_av_a_c2(:,:,t)=P_av_a_c2(:,:,t-1).*(P_C*((L_av_xa_c2(:,:,t)+M_av_xa_c2(:,:,t))) );
        P_av_v_c2(:,:,t)=P_av_v_c2(:,:,t-1).*(P_C*((L_av_xv_c2(:,:,t)+M_av_xv_c2(:,:,t))) );
    
    end
end

%% Normalisation and plotting
z=P_av_a_c2(:,:,1).*P_av_v_c2(:,:,1);
for o=1:t_max
%     P_av(:,:,o)=P_av(:,:,o)/sum(sum(P_av(:,:,o)));
    z(:,:,o)=z(:,:,o)/sum(sum(z(:,:,o)));

%     P_av_v_c2(:,:,o)=P_av_v_c2(:,:,o)/sum(sum(P_av_v_c2(:,:,o)));
%     P_av_a_c2(:,:,o)=P_av_a_c2(:,:,o)/sum(sum(P_av_a_c2(:,:,o)));
end
%%
[x,y]=meshgrid(lsteps,msteps);
%the results from the update equation
%C=1: first the audio-visual

figure(1)
contour(x,y,z,5,'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
% figure(2)
% contour(x,y,P_av(:,:,2),'ShowText','on')
% hold on;
% xlabel("Location")
% ylabel("Semantic meaning")
% figure(3)
% contour(x,y,P_av(:,:,3),'ShowText','on')
% hold on;
% xlabel("Location")
% ylabel("Semantic meaning")