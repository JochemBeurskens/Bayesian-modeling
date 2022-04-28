%% Setting the model parameters
% rng(0); %setting the random seed, so that the draws will have the same result each time, for now this is useful.

k=1.0; %setting the multiplication for the linear functions
t_max=10; %setting the amount of time steps for the source activity
n=11; %setting the amount of possible locations

P_C=0.5;% setting the probability of selecting C=1 or C=2
P_av_given_LMf=zeros(n,n,n,t_max); %make an array to store the probabilities in, this is an n (range of L) by n (range of M) by n (range of c or f) by t array (so that the different times can be stored as well, but designed s.t. the first entry contains the first prior)
P_a_given_LMf_c2=zeros(n,n,n,t_max);%same but for the C=2 and the auditory percept
P_v_given_LMf_c2=zeros(n,n,n,t_max);%same but for the C=2 and the visual percept

repeat_l=20;
repeat_m=20;
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

L_av_xa_c2=zeros(2*n-1,2*n-1,t_max);
L_av_xv_c2=zeros(2*n-1,2*n-1,t_max);
M_av_xa_c2=zeros(2*n-1,2*n-1,t_max);
M_av_xv_c2=zeros(2*n-1,2*n-1,t_max);

P_av_a=zeros(2*n-1,2*n-1,t_max);
P_av_v=zeros(2*n-1,2*n-1,t_max);

P_av_a_c2=zeros(2*n-1,2*n-1,t_max);
P_av_v_c2=zeros(2*n-1,2*n-1,t_max);
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
    L_a_c2(ls)=(normcdf(lsteps(ls+1),-2,sig_l_as_c2)-normcdf(lsteps(ls),-2,sig_l_as_c2));
    M_a_c2(ls)=(normcdf(lsteps(ls+1),-2,sig_m_as_c2)-normcdf(lsteps(ls),-2,sig_m_as_c2));
    L_v_c2(ls)=(normcdf(lsteps(ls+1),2,sig_l_vs_c2)-normcdf(lsteps(ls),2,sig_l_vs_c2));
    M_v_c2(ls)=(normcdf(lsteps(ls+1),2,sig_m_vs_c2)-normcdf(lsteps(ls),2,sig_m_vs_c2));
end
%In this case we no longer need to multiply with the prior distribution
%explicitly, as algorithm 1 of A tutorial on particle filters by Maarten Speekenbrink
%shows that obtaining integral(f(y)*p(y)dy) can be obtained by sampling
%from p(y) and evaluating f(y) at those locations.y
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
    for li=1:repeat_l %looping over all location grid points
        %mi=1;
        for mi=1:repeat_m %looping over all semantic meaning grid points
            M_av_s=normrnd(0.0,sig_m_s); 
            L_av_s=normrnd(0.0,sig_l_s); %drawing the M and L from the probability distributions for these parameters, this should simulate the effect of multiplication with the P(L_{av}^{s}) etc.
            %C=1, run model and plot
            for ls =1:length(lsteps)-1 %now for each centre of the stimuli the range of observations is looked at
                for ms =1:length(msteps)-1
                    L_av_xa(ls,ms,t)=L_av_xa(ls,ms,t)+(normcdf(lsteps(ls+1),L_av_s,sig_la)-normcdf(lsteps(ls),L_av_s,sig_la));
                    L_av_xv(ls,ms,t)=L_av_xa(ls,ms,t)+(normcdf(lsteps(ls+1),L_av_s,sig_lv)-normcdf(lsteps(ls),L_av_s,sig_lv));
                    M_av_xa(ls,ms,t)=M_av_xa(ls,ms,t)+(normcdf(msteps(ms+1),M_av_s,sig_ma)-normcdf(msteps(ms),M_av_s,sig_ma));
                    M_av_xv(ls,ms,t)=M_av_xa(ls,ms,t)+(normcdf(msteps(ms+1),M_av_s,sig_mv)-normcdf(msteps(ms),M_av_s,sig_mv));
                    %the above calcualtes the spread of the observations
                    %given each of the locations for the centres in the
                    %grid. By summing the probabilities distributions of
                    %the observations, for all possible stimulus locations
                    %the integral is performed in a numerical fashion. Note
                    %that the multiplication with the prior probabilities
                    %for the stimulus locations has been postponed to line
                    %115, so that the code is less clutered. But as
                    %associativity holds for sums of real numbers, this
                    %shouldn't influence the result (i.e. C*a+C*b=C(a+b)).
                end
            end
            %mi=mi+1;
        end
        %li=li+1;
    end
    if t==1
        P_av_a(:,:,t)=P_C*((L_av_xa(:,:,t)+M_av_xa(:,:,t)));
        P_av_v(:,:,t)=P_C*((L_av_xv(:,:,t)+M_av_xv(:,:,t)));
    else 
        P_av_a(:,:,t)=P_av_a(:,:,t-1).*(P_C*((L_av_xa(:,:,t)+M_av_xa(:,:,t))));
        P_av_v(:,:,t)=P_av_v(:,:,t-1).*(P_C*((L_av_xv(:,:,t)+M_av_xv(:,:,t))));
        %here the probability distributions for the locations and the
        %semantic meanings are added, and then multiplied with the prior
        %probabilities for the stimulus locations. The update step is also
        %performed, thus multiplying the probabilities for the current time
        %step with those for the previous time step. At this point all that
        %rests is normalisation but further the integration has been
        %performed.
    end
%     %C=2
%     %in this case each possible visual stimulus location has to be matched
%     %to a possible auditory stimulus location, creating a double for loop
%     %over the stimulus locations, but the inner for loops for estimating
%     %the probability at a observation location can remain as it doesn't
%     %matter where the stimulus is located due to the observation
%     %probability being investigated at each and every possible grid point
%     %anayways.
%     for li_a=1:repeat_l %looping over all location grid points
%         %mi_a=1;
%         for mi_a=1:repeat_m %looping over all semantic meaning grid points
%             %li_v=1;
%             for li_v=1:repeat_l %looping over all location grid points
%                 %mi_v=1;
%                 for mi_v=1:repeat_m
%                     M_a=normrnd(0.25,sig_m_as_c2); 
%                     M_v=normrnd(-1.25,sig_m_vs_c2);
%                     L_a=normrnd(-0.25,sig_l_as_c2);
%                     L_v=normrnd(1.25,sig_l_vs_c2); %drawing the M and L from the probability distributions for these parameters, this should simulate the effect of multiplication with the P(L_{av}^{s}) etc.
%                     %C=1, run model and plot
%                     for ls =1:length(lsteps)-1 %now for each centre of the stimuli the range of observations is looked at
%                         for ms =1:length(msteps)-1
%                             L_av_xa_c2(ls,ms,t)=L_av_xa_c2(ls,ms,t)+(normcdf(lsteps(ls+1),L_av_sa,sig_la_c2)-normcdf(lsteps(ls),L_av_sa,sig_la_c2));
%                             L_av_xv_c2(ls,ms,t)=L_av_xa_c2(ls,ms,t)+(normcdf(lsteps(ls+1),L_av_sv,sig_lv_c2)-normcdf(lsteps(ls),L_av_sv,sig_lv_c2));
%                             M_av_xa_c2(ls,ms,t)=M_av_xa_c2(ls,ms,t)+(normcdf(msteps(ms+1),M_av_sa,sig_ma_c2)-normcdf(msteps(ms),M_av_sa,sig_ma_c2));
%                             M_av_xv_c2(ls,ms,t)=M_av_xa_c2(ls,ms,t)+(normcdf(msteps(ms+1),M_av_sv,sig_mv_c2)-normcdf(msteps(ms),M_av_sv,sig_mv_c2));
%                             %the above calcualtes the spread of the observations
%                             %given each of the locations for the centres in the
%                             %grid. By summing the probabilities distributions of
%                             %the observations, for all possible stimulus locations
%                             %the integral is performed in a numerical fashion. Note
%                             %that the multiplication with the prior probabilities
%                             %for the stimulus locations has been postponed to line
%                             %115, so that the code is less clutered. But as
%                             %associativity holds for sums of real numbers, this
%                             %shouldn't influence the result (i.e. C*a+C*b=C(a+b)).
%                         end
%                     end
%                     %mi_v=mi_v+1;
%                 end
%                 %li_v=li_v+1;
%             end
%             %mi_a=mi_a+1;
%         end
%         %li_a=li_a+1;
%     end
%     if t==1
%         P_av_a_c2(:,:,t)=P_C*((L_av_xa_c2(:,:,t)*+M_av_xa_c2(:,:,t)));
%         P_av_v_c2(:,:,t)=P_C*((L_av_xv_c2(:,:,t)+M_av_xv_c2(:,:,t)));
%     else 
%         P_av_a_c2(:,:,t)=P_av_a_c2(:,:,t-1)*(P_C*((L_av_xa_c2(:,:,t)+M_av_xa_c2(:,:,t))));
%         P_av_v_c2(:,:,t)=P_av_v_c2(:,:,t-1)*(P_C*((L_av_xv_c2(:,:,t)+M_av_xv_c2(:,:,t))));
%         %here the probability distributions for the locations and the
%         %semantic meanings are added, and then multiplied with the prior
%         %probabilities for the stimulus locations. The update step is also
%         %performed, thus multiplying the probabilities for the current time
%         %step with those for the previous time step. At this point all that
%         %rests is normalisation but further the integration has been
%         %performed.
%     end
end

%% Normalisation and plotting
for o=1:t_max
    P_av_v(:,:,o)=P_av_v(:,:,o)/sum(sum(P_av_v(:,:,o)));
    P_av_a(:,:,o)=P_av_a(:,:,o)/sum(sum(P_av_a(:,:,o))); %normalising the probability distributions
    P_av_v_c2(:,:,o)=P_av_v_c2(:,:,o)/sum(sum(P_av_v_c2(:,:,o)));
    P_av_a_c2(:,:,o)=P_av_a_c2(:,:,o)/sum(sum(P_av_a_c2(:,:,o)));
end
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
contour(x,y,P_av_a(:,:,10),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
