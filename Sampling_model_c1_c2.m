%for the likelihood calculation I will be using almost the same code, but
%as the very large arrays of a and v take too much space, I will not use
%them here. I will only calculate the parameters that decide where in these
%there are non-zero values, and then use this to make a frequentist
%estimation of the likelihood of a given draw

%% setting parameters

% rng(0); %setting the random seed, so that the draws will have the same result each time, for now this is useful.

k=1.0; %setting the multiplication for the linear functions
t_max=3; %setting the amount of time steps for the source activity
n=10; %setting the amount of possible locations

P_C=0.5;% setting the probability of selecting C=1 or C=2
P_av_given_LMf=zeros(n,n,n,t_max+1); %make an array to store the probabilities in, this is an n (range of L) by n (range of M) by n (range of c or f) by t array (so that the different times can be stored as well, but designed s.t. the first entry contains the first prior)
P_a_given_LMf_c2=zeros(n,n,n,t_max);%same but for the C=2 and the auditory percept
P_v_given_LMf_c2=zeros(n,n,n,t_max);%same but for the C=2 and the visual percept

repeat_l=5;
repeat_m=5;
repeat=1000;
Percept_L_av=zeros(n,n,repeat); %make arrays to store the percepts for all the repeated trials, note that there will be n,n,repeat entries (n draws from the space of l, for each draw of l 5 draws from m follow, then 10000 times repeated draws for those positions in m,l space)
Percept_M_av=zeros(n,n,repeat);

lrange=[-5 5];
mrange=[-5 5];
crange=[-0.5 0.5];

l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
m=linspace(mrange(1),mrange(2),n);
c=linspace(crange(1),crange(2),n);
%parameters for C=1
sig_e_s=0.1;
sig_l_s=0.6;
sig_m_s=0.1;%setting noise parameters for the generating locations step. used for both C=1 and 2

sig_e_ax=0.1;
sig_e_vx=0.1;
sig_la=0.5;
sig_ma=0.5;
sig_lv=0.5;
sig_mv=0.5; %setting the noise parameters for the subject's observation

%parameters for C=2
sig_e_s_c2=0.1;
sig_l_s_c2=0.6;
sig_m_s_c2=0.1;%setting noise parameters for the generating locations step. used for both C=1 and 2

sig_e_as_c2=0.1;
sig_e_vs_c2=0.1;
sig_l_as_c2=0.6;
sig_m_as_c2=0.6;
sig_l_vs_c2=0.6;
sig_m_vs_c2=0.6; %setting the additional noise parameters for the case of c=2

sig_e_ax_c2=0.1;
sig_e_vx_c2=0.1;
sig_la_c2=0.5;
sig_ma_c2=0.5;
sig_lv_c2=0.5;
sig_mv_c2=0.5; %setting the noise parameters for the subject's observation
%% drawing the a_x and v_x for which the likelihoods are to be estimated
% [f_a_x,f_v_x,i_l_x,i_m_x,i_la_x,i_ma_x,i_lv_x,i_mv_x]=Likelihood(k,t_max,n,c,l,m,sig_e_s,sig_l_s,sig_m_s,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);
%% Now  going over the entire space state to find estimates of the probabilities for the parameters in the generated a_x and v_x
%note that it is here useful to also count the instances that are the same
%as the values in the a_x and v_x
store_percept_av=zeros(n,n,t_max);
store_percept_a=zeros(n,n,t_max);
store_percept_v=zeros(n,n,t_max);
for t=1:t_max
    store_a=zeros(n,n,n);
    store_v=zeros(n,n,n);
    store_a_c2=zeros(n,n,n);
    store_v_c2=zeros(n,n,n);
    %Drawing random numbers from the distributions for L and M, thus in
    %principle taking the probability distributions of these variables into
    %account. This will then be used as the locations for the input variables
    %in the Likelihood_given_input_c1 formula, and thus will be used to
    %generate samples. The amount of times a sample will hit a specific
    %location will then be the fraction of times a specific value is expected,
    %and can then be used as a pseudo-probability (i.e. non-normalised prob.).
    %Not need to loop over c (or f), as this is location independent so can be
    %obtained once and then get the probabilities for specific values.
    % li=1;
    for li=1:repeat_l
    %     mi=1;
        for mi=1:repeat_m
            M=normrnd(2,sig_m_s); %drawing the M and L from the probability distributions for these parameters, this should simulate the effect of multiplication with the P(L_{av}^{s}) etc.
            L=normrnd(0,sig_l_s);
            M_c2=normrnd(2,sig_m_as_c2); %drawing the M and L from the probability distributions for these parameters, this should simulate the effect of multiplication with the P(L_{av}^{s}) etc.
            M_v=normrnd(-2,sig_m_vs_c2);
            L_c2=normrnd(-3,sig_l_as_c2);
            L_v=normrnd(-3,sig_l_vs_c2);
            
            f_v_s_store=zeros(repeat,t_max);
            f_a_s_store=zeros(repeat,t_max);
            for q=1:repeat
                %C=1, run model and plot
                [i_fax,i_fvx,i_l_plt_stat,i_m_plt_stat,i_la_plt_stat,i_ma_plt_stat,i_lv_plt_stat,i_mv_plt_stat,i_percept_L_av,i_percept_M_av]=Sample_model(k,1,n,c,l,m,L,M,sig_e_s,sig_l_s,sig_m_s,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);
%                 [i_fax,i_fvx,i_las_plt,i_mas_plt,i_lvs_plt,i_mvs_plt,i_la_plt,i_ma_plt,i_lv_plt,i_mv_plt]=Likelihood_given_input(k,t_max,n,c,l,m,L,M,sig_e_as,sig_e_vs,sig_l_as,sig_l_vs,sig_m_as,sig_m_vs,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv,L_v,M_v);
                store_a(i_la_plt_stat,i_ma_plt_stat,i_fax)=store_a(i_la_plt_stat,i_ma_plt_stat,i_fax)+1;
                store_v(i_lv_plt_stat,i_mv_plt_stat,i_fvx)=store_v(i_lv_plt_stat,i_mv_plt_stat,i_fvx)+1;
                store_percept_av(i_percept_L_av,i_percept_M_av,t)=store_percept_av(i_percept_L_av,i_percept_M_av,t)+1; %note that in this 3D structure the percept can be stored simultaneously for the L and M and t dimensions
                %repeating the above for C=2
                [i_fax,i_fvx,i_las_plt,i_mas_plt,i_lvs_plt,i_mvs_plt,i_la_plt,i_ma_plt,i_lv_plt,i_mv_plt,i_percept_L_a,i_percept_M_a,i_percept_L_v,i_percept_M_v]=Sample_model(k,1,n,c,l,m,L_c2,M_c2,sig_e_as_c2,sig_e_vs_c2,sig_l_as_c2,sig_l_vs_c2,sig_m_as_c2,sig_m_vs_c2,sig_e_ax_c2,sig_e_vx_c2,sig_la_c2,sig_ma_c2,sig_lv_c2,sig_mv_c2,L_v,M_v);
                store_a_c2(i_la_plt,i_ma_plt,i_fax)=store_a_c2(i_la_plt,i_ma_plt,i_fax)+1;
                store_v_c2(i_lv_plt,i_mv_plt,i_fvx)=store_v_c2(i_lv_plt,i_mv_plt,i_fvx)+1;
                store_percept_a(i_percept_L_a,i_percept_M_a,t)=store_percept_a(i_percept_L_a,i_percept_M_a,t)+1;
                store_percept_v(i_percept_L_v,i_percept_M_v,t)=store_percept_v(i_percept_L_v,i_percept_M_v,t)+1;
            end
            %the storage containers above will contain the distribution of the
            %visual and auditory probabilities, multiplied with the probability
            %for L_{av}^{s} etc. (i.e. the integration terms).
            mi=mi+1;
        end
        li=li+1;
    end
    if t==1
        P_av_given_LMf(:,:,:,t)=(P_C*((store_a+store_v)));
        P_av_given_LMf(:,:,:,t)=P_av_given_LMf(:,:,:,t);
        
        P_a_given_LMf_c2(:,:,:,t)=(P_C*((store_a_c2)));
        P_a_given_LMf_c2(:,:,:,t)=P_a_given_LMf_c2(:,:,:,t);
        
        P_v_given_LMf_c2(:,:,:,t)=(P_C*((store_v_c2)));
        P_v_given_LMf_c2(:,:,:,t)=P_v_given_LMf_c2(:,:,:,t);
        
        %normalising
        normalisation_const=sum(sum(sum(P_av_given_LMf(:,:,:,t)))) + sum(sum(sum(P_a_given_LMf_c2(:,:,:,t)))) + sum(sum(sum(P_v_given_LMf_c2(:,:,:,t))));
    	P_av_given_LMf(:,:,:,t)=P_av_given_LMf(:,:,:,t)/normalisation_const;
    	P_a_given_LMf_c2(:,:,:,t)=P_a_given_LMf_c2(:,:,:,t)/normalisation_const;
    	P_v_given_LMf_c2(:,:,:,t)=P_v_given_LMf_c2(:,:,:,t)/normalisation_const;
    else
        P_av_given_LMf(:,:,:,t)=P_av_given_LMf(:,:,:,t-1).*(P_C*((store_a+store_v)));
        P_av_given_LMf(:,:,:,t)=P_av_given_LMf(:,:,:,t);
        
        P_a_given_LMf_c2(:,:,:,t)=P_a_given_LMf_c2(:,:,:,t-1).*(P_C*((store_a_c2)));
        P_a_given_LMf_c2(:,:,:,t)=P_a_given_LMf_c2(:,:,:,t);
        
        P_v_given_LMf_c2(:,:,:,t)=P_v_given_LMf_c2(:,:,:,t-1).*(P_C*((store_v_c2)));
        P_v_given_LMf_c2(:,:,:,t)=P_v_given_LMf_c2(:,:,:,t);
        
        %normalising
        normalisation_const=sum(sum(sum(P_av_given_LMf(:,:,:,t)))) + sum(sum(sum(P_a_given_LMf_c2(:,:,:,t)))) + sum(sum(sum(P_v_given_LMf_c2(:,:,:,t))));   
    	P_av_given_LMf(:,:,:,t)=P_av_given_LMf(:,:,:,t)/normalisation_const;
    	P_a_given_LMf_c2(:,:,:,t)=P_a_given_LMf_c2(:,:,:,t)/normalisation_const;
    	P_v_given_LMf_c2(:,:,:,t)=P_v_given_LMf_c2(:,:,:,t)/normalisation_const;
    end
end
%% now plotting
[x,y]=meshgrid(l,m);
%first the results from the update equation
%C=1: first the audio-visual
figure(1)
contour(x,y,P_av_given_LMf(:,:,5,1),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(2)
contour(x,y,P_av_given_LMf(:,:,5,2),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(3)
contour(x,y,P_av_given_LMf(:,:,5,3),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
%C=2: first the auditory
figure(4)
contour(x,y,P_a_given_LMf_c2(:,:,5,1),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(5)
contour(x,y,P_a_given_LMf_c2(:,:,5,2),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(6)
contour(x,y,P_a_given_LMf_c2(:,:,5,3),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
%C=2: then the visual
figure(7)
contour(x,y,P_v_given_LMf_c2(:,:,5,1),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(8)
contour(x,y,P_v_given_LMf_c2(:,:,5,2),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(9)
contour(x,y,P_v_given_LMf_c2(:,:,5,3),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")

%% now plotting the percepts
[x,y]=meshgrid(l,m);
%C=1: first the audio-visual
figure(1)
contour(x,y,store_percept_av(:,:,1),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(2)
contour(x,y,store_percept_av(:,:,2),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(3)
contour(x,y,store_percept_av(:,:,3),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
%C=2: first the auditory
figure(4)
contour(x,y,store_percept_a(:,:,1),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(5)
contour(x,y,store_percept_a(:,:,2),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(6)
contour(x,y,store_percept_a(:,:,3),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
%C=2: then the visual
figure(7)
contour(x,y,store_percept_v(:,:,1),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(8)
contour(x,y,store_percept_v(:,:,2),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")
figure(9)
contour(x,y,store_percept_v(:,:,3),'ShowText','on')
hold on;
xlabel("Location")
ylabel("Semantic meaning")