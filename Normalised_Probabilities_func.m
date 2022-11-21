t_max=20; %set number of time steps in the process
n=1000; %set the number of possible locations, the larger the more accurate, but also the more computation time is required

c=1.0; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
k=0.3; %setting the multiplication for the linear functions
lrange=[-20 20]; %the range of l, actual locations are between -10.5 and 10.5
l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings

Causal_structure=randi([0,1],1,1); %this makes it such that the causal structure is decided upon, ie. C=1 or C=2, which decides which generative model should be used
if Causal_structure == 0
    % getting the distribution of the parameters from the counts in the generative model
    %parameters for the distributions for C=1
    sig_l_s=2;
    sig_lax=1;
    sig_lvx=1;
    sig_e_s=1;
    sig_e_ax=1;
    sig_e_vx=1;
    p_c=0.3;
    p_h=0.2;
    %sampling from the generative model for the case of C=1
    L=zeros(n,n);
    [i_la_plt,i_lv_plt,m,a_x,v_x,L_av_s,c_av_s,m_c,m_h]=gen_model(k,t_max,n,l,sig_l_s,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h);
elseif Causal_structure==1    
    % getting the distribution of the parameters from the counts in the generative model
    %parameters for the distributions for C=2
    sig_l_vs=2;
    sig_l_as=sig_l_vs; %when assuming that the original location is drawn from one normal distribution, only one sigma_L^s is needed, same for c. They can now be seperate draws though, which for C=1 was not possible
    sig_lax=1;
    sig_lvx=1;
    sig_e_s=1;
    sig_e_ax=1;
    sig_e_vx=1;
    p_c=0.3;
    p_h=0.2;
    %sampling from the generative model for the case of C=2
    L=zeros(n,n);
    [i_la_plt,i_lv_plt,m,a_x,v_x,L_av_vs,L_av_as,c_av_s,m_c,m_h]=gen_model(k,t_max,n,l,sig_l_vs,sig_l_as,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h);
end
%% setting the parameters used in the probability distributions
L_xa=i_la_plt(1);
L_xv=i_lv_plt(1);
L_s=0; %the mean of the prior
%% Now the C=1 pseudo-probability for the L parameters
forefac1    = 2*pi*sqrt(sig_lax^2*sig_l_s^2+sig_lvx^2*sig_l_s^2+sig_lvx^2*sig_lax^2);
exponent1   = exp(-0.5* ( ((L_xa-L_xv)^2+(L_xa-L_s)^2+(L_xv-L_s)^2)/forefac1 ) );
P_LC1       = exponent1/forefac1;

%% Now the C=2 pseudo-probability for the L parameters
forefac2    = 2*pi*sqrt( (sig_lax^2+sig_l_s^2)*(sig_lvx^2+sig_l_s^2) );
exponent2   = exp(-0.5* ( ((L_xa-L_s)^2/((sig_lax^2+sig_l_s^2)))+((L_xv-L_s)^2)/((sig_lvx^2*sig_l_s^2))  ) );
P_LC2       = exponent2/forefac2;

%% normalisation
norm_fac    = P_LC2+P_LC1;
P_LC1       = P_LC1/norm_fac;
P_LC2       = P_LC2/norm_fac;

disp(P_LC1)
disp(P_LC2)