%% Setting the model parameters
% clear all;
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

L_av_xa=zeros(n,n,t_max);
L_av_xv=zeros(n,n,t_max);
M_av_xa=zeros(n,n,t_max);
M_av_xv=zeros(n,n,t_max);
C_L=zeros(n,n); %a parameter to fill with the manual check
C_M=zeros(n,n);
I_L=zeros(n,n); %a parameter to fill with the manual check
I_M=zeros(n,n);

L_av_xa_c2=zeros(n,n,t_max);
L_av_xv_c2=zeros(n,n,t_max);
M_av_xa_c2=zeros(n,n,t_max);
M_av_xv_c2=zeros(n,n,t_max);

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
%% check with manual solution for the result of the integral
%% First C=1
% sig_int_L=sqrt(sig_la^2*sig_lv^2+sig_la^2*sig_l_s^2+sig_l_s^2*sig_lv^2);
% sig_int_M=sqrt(sig_ma^2*sig_mv^2+sig_ma^2*sig_m_s^2+sig_m_s^2*sig_mv^2);
% forefactor_L=1/(2*pi*sig_int_L);
% forefactor_M=1/(2*pi*sig_int_M);
% 
% for iv =1:length(lsteps) %now for each centre of the stimuli the range of observations is looked at
%     for ia =1:length(lsteps)
%         I_L_exponent=-0.5*((((lsteps(iv)-lsteps(ia))^2*sig_l_s^2+(lsteps(iv)^2)*sig_la^2+(lsteps(ia)^2)*sig_lv^2))/sig_int_L^2);
%         I_M_exponent=-0.5*((((msteps(iv)-msteps(ia))^2*sig_m_s^2+(msteps(iv)^2)*sig_ma^2+(msteps(ia)^2)*sig_mv^2))/sig_int_M^2);
%         I_L(ia,iv)=I_L(ia,iv)+forefactor_L*exp(I_L_exponent);
%         I_M(ia,iv)=I_M(ia,iv)+forefactor_M*exp(I_M_exponent);
%     end
% end
% [x,y]=meshgrid(lsteps,msteps);
% 
% I=I_L.*I_M;
% I=I/sum(sum(I(:,:)));
% figure(4);
% contour(x,y,I,5,'ShowText','on')
%% Then C=2
sig_int_L=sqrt( (sig_la_c2^2+sig_l_as_c2^2)*(sig_l_vs_c2^2*sig_lv_c2^2) );
sig_int_M=sqrt( (sig_ma_c2^2+sig_m_as_c2^2)*(sig_m_vs_c2^2*sig_mv_c2^2) );
forefactor_L=1/(2*pi*sig_int_L);
forefactor_M=1/(2*pi*sig_int_M);

for iv =1:length(lsteps) %now for each centre of the stimuli the range of observations is looked at
    for ia =1:length(lsteps)
        I_L_exponent=-0.5*( (( lsteps(ia)^2  )/(sig_la_c2^2+sig_l_as_c2^2)) + (( lsteps(iv)^2  )/(sig_lv_c2^2+sig_l_vs_c2^2)));
        I_M_exponent=-0.5*( (( msteps(ia)^2  )/(sig_ma_c2^2+sig_m_as_c2^2)) + (( msteps(iv)^2  )/(sig_mv_c2^2+sig_m_vs_c2^2)));
        I_L(ia,iv)=I_L(ia,iv)+forefactor_L*exp(I_L_exponent);
        I_M(ia,iv)=I_M(ia,iv)+forefactor_M*exp(I_M_exponent);
    end
end
[x,y]=meshgrid(lsteps,msteps);

I=I_L.*I_M;
I=I/sum(sum(I(:,:)));
figure(5);
contour(x,y,I,5,'ShowText','on')