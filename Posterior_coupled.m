clear all;
clc;
rng(6); %setting the random number generator seed, for reproducible results

%% parameter specification
t_max=1; %set number of time steps in the process
n=20; %maybe not discritize, set the number of possible locations, the larger the more accurate, but also the more computation time is required
n_r=1; %number of stes of a_x or v_x obtained for each set of L_s, M_s and c_s
n_lm=600; %number of sets of L^x and M^x that are considered for the integrals
c=0.8; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
k=1; %setting the multiplication for the linear functions
lrange=[-20 20]; %the range of l, actual locations are between -10.5 and 10.5
l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
m=[1,2];
% sig_l_s=3;
% sig_l_vs=5;
% sig_l_as=sig_l_vs; %when assuming that the original location is drawn from one normal distribution, only one sigma_L^s is needed, same for c. They can now be seperate draws though, which for C=1 was not possible
% sig_lax=5;
% sig_lvx=5;
% sig_e_s=2;
% sig_e_as=2;
% sig_e_vs=2;
% sig_e_ax=2;
% sig_e_vx=2;
% sig_m_s=1;
% sig_m_as=1;
% sig_m_vs=sig_m_as;
% sig_m_ax=1;
% sig_m_vx=1;
% p_ci=4; %if x param > p_c: m=1, else m=2

sig_l_s=5;
sig_l_vs=5;
sig_l_as=sig_l_vs; %when assuming that the original location is drawn from one normal distribution, only one sigma_L^s is needed, same for c. They can now be seperate draws though, which for C=1 was not possible
sig_lax=8;
sig_lvx=5;
sig_e_s=1;
sig_e_as=2;
sig_e_vs=2;
sig_e_ax=1;
sig_e_vx=1;
sig_m_s=1;
sig_m_as=1;
sig_m_vs=sig_m_as;
sig_m_ax=1;
sig_m_vx=0.5;
p_ci=.5;
p_h=0.5;
mean_rw=0; %setting the mean of the random walk
%now we only need to obtain 1 set of a^x and v^x, and then evaluate it all
%from the perspective of the different L^s and M^s values in the grid
p_causal=0.5;
%% obtaining the a_x, v_x values
[i_la_plt,i_lv_plt,m_a,m_v,a_x,v_x,i_l,c_av_s,m_av_xa,m_av_xv,m_s]=gen_model_full(k,t_max,n,l,sig_l_s,sig_lvx,sig_lax,sig_m_s,sig_m_as,sig_m_vs,sig_m_ax,sig_m_vx,sig_e_s,sig_e_ax,sig_e_vx,p_ci,p_h,n_r,n_lm,mean_rw);
% surf(v_x)
%% now obtaining the measurement prob. distribution for the data above (Causal=1)
n_ls=50;
n_ms=50;
n_cs=50;
i_l=[i_l,normrnd(0,sig_l_s,1,(n_ls-1))];
for il=1:n_ls
    [~,d_l]=dsearchn(i_l(il),l'); %getting the distances from the observations to the possible locations
    i_l(il)=find(d_l==min(d_l(:)));
end
m_s=[m_s,normrnd(0.5,sig_m_s,1,(n_ms-1))];
for im=1:n_ms
    [~,d_m]=dsearchn(m_s(im),m'); %getting the distances from the observations to the possible locations
    m_s(im)=find(d_m==min(d_m(:)));
end
c_av_s=[c_av_s,normrnd(mean_rw,sig_e_s,1,(n_cs-1))];
like_ax=zeros(n,2,t_max,n_lm);
like_vx=zeros(n,2,t_max,n_lm);
posterior_c1=zeros(n,2,t_max,n_lm,n_ls);
pc=normcdf(p_ci,m_s,sig_m_ax);
p_c=[1-pc pc]; %prob for m=1 and m=2, note that if the m_x draw is larger than p_c it is put in m=1, thus 1-normcdf is the prob
for x=1:n_ls
prob_Ls=normpdf(l(i_l(x)),0,sig_l_s);
prob_Ms=normpdf(m_s(x),0.5,sig_m_s);
prob_cs=normpdf(c_av_s(x),mean_rw,sig_e_s);
for q=1:n_lm
    for o=1:n
        for p=1:2
            for t =1:t_max
                prob_Lvx=normpdf(l(i_l(x)),l(o),sig_lvx); %N(L^x,L^s,sig)
                prob_vx(o,p,t,q)=prob_Lvx'.*normpdf(v_x(o,p,t,q),mean_rw,sig_e_s).*p_c(p);
                %prob_vx(i_lv_plt(t,o),m_v(t,o),t,o,p)=p_c(m_v(t,o))*prob_Lvx(i_lv_plt(t,o))*normpdf(v_x(i_lv_plt(t,o),m_v(t,o),t,o,p),mean_rw,sig_e_s);
    
                prob_Lax=normpdf(l(i_l(x)),l(o),sig_lax); %N(L^x,L^s,sig)
                prob_ax(o,p,t,q)=prob_Lax'.*normpdf(a_x(o,p,t,q),mean_rw,sig_e_s).*p_c(p);
                %prob_ax(i_la_plt(t,o),m_a(t,o),t,o,p)=p_c(m_a(t,o))*prob_Lax(i_la_plt(t,o))*normpdf(a_x(i_la_plt(t,o),m_a(t,o),t,o,p),mean_rw,sig_e_s);
            end
        end 
    end
end
posterior_c1(:,:,:,:,x)=p_causal*prob_Ls.*prob_Ms.*prob_cs.*prob_vx.*prob_ax;
end

%% now calculating the probabilities for the same ax and vx but for C=2
i_la=[normrnd(0,sig_l_as,1,(n_ls))];
i_lv=[normrnd(0,sig_l_vs,1,(n_ls))];
for il=1:n_ls
    [~,d_l]=dsearchn(i_la(il),l'); %getting the distances from the observations to the possible locations
    i_la(il)=find(d_l==min(d_l(:)));
    [~,d_l]=dsearchn(i_lv(il),l'); %getting the distances from the observations to the possible locations
    i_lv(il)=find(d_l==min(d_l(:)));
end
m_sa=[normrnd(0.5,sig_m_as,1,(n_ms))];
m_sv=[normrnd(0.5,sig_m_as,1,(n_ms))];
for im=1:n_ms
    [~,d_m]=dsearchn(m_sa(im),m'); %getting the distances from the observations to the possible locations
    m_sa(im)=find(d_m==min(d_m(:)));
    [~,d_m]=dsearchn(m_sv(im),m'); %getting the distances from the observations to the possible locations
    m_sv(im)=find(d_m==min(d_m(:)));
end
c_a_s=[normrnd(mean_rw,sig_e_s,1,(n_cs))];
c_v_s=[normrnd(mean_rw,sig_e_s,1,(n_cs))];
posterior_c2=zeros(n,2,t_max,n_lm,n_ls,n_ls);
pca=normcdf(p_ci,m_sa,sig_m_ax);
p_ca=[1-pca pca]; %prob for m=1 and m=2, note that if the m_x draw is larger than p_c it is put in m=1, thus 1-normcdf is the prob
pcv=normcdf(p_ci,m_sv,sig_m_vx);
p_cv=[1-pcv pcv];
for x=1:n_ls
for y=1:n_ls
prob_Las=normpdf(l(i_la(x)),0,sig_l_as);
prob_Mas=normpdf(m_sa(x),0.7,sig_m_as);
prob_cas=normpdf(c_a_s(x),mean_rw,sig_e_as);
prob_Lvs=normpdf(l(i_lv(y)),5,sig_l_vs);
prob_Mvs=normpdf(m_sv(y),0.7,sig_m_vs);
prob_cvs=normpdf(c_v_s(y),mean_rw,sig_e_vs);
for q=1:n_lm
    for o=1:n
        for p=1:2
            for t =1:t_max
                prob_Lvx=normpdf(l(i_lv(x)),l(o),sig_lvx); %N(L^x,L^s,sig)
                prob_vx(o,p,t,q)=prob_Lvx'.*normpdf(v_x(o,p,t,q),mean_rw,sig_e_vs).*p_cv(p);
                %prob_vx(i_lv_plt(t,o),m_v(t,o),t,o,p)=p_c(m_v(t,o))*prob_Lvx(i_lv_plt(t,o))*normpdf(v_x(i_lv_plt(t,o),m_v(t,o),t,o,p),mean_rw,sig_e_s);
    
                prob_Lax=normpdf(l(i_la(y)),l(o),sig_lax); %N(L^x,L^s,sig)
                prob_ax(o,p,t,q)=prob_Lax'.*normpdf(a_x(o,p,t,q),mean_rw,sig_e_as).*p_ca(p);
                %prob_ax(i_la_plt(t,o),m_a(t,o),t,o,p)=p_c(m_a(t,o))*prob_Lax(i_la_plt(t,o))*normpdf(a_x(i_la_plt(t,o),m_a(t,o),t,o,p),mean_rw,sig_e_s);
            end
        end 
    end
end
posterior_c2(:,:,:,:,x,y)=(1-p_causal)*prob_Las.*prob_Lvs.*prob_Mas.*prob_Mvs.*prob_cas.*prob_cvs.*prob_vx.*prob_ax;
end
end
post_c1=sum(sum(posterior_c1(:,:,:,:,:),5),4);
post_c2=sum(sum(sum(posterior_c2(:,:,:,:,:,:),6),5),4);
post_c1n=post_c1./(post_c1+post_c2);
post_c2n=post_c2./(post_c1+post_c2);
figure(1)
surf(post_c1n)
figure(2)
surf(post_c2n)
% figure(1)
% surf(sum(sum(posterior_c2(:,:,:,:),5),4))
% xlabel('Stimulus Location'), ylabel('Stimulus Meaning'), zlabel('Likelihood(stimuli|a^{x})')
% figure(2)
% surf(sum(sum(posterior_c1(:,:,:,:),5),4))
% xlabel('Stimulus Location'), ylabel('Stimulus Meaning'), zlabel('Likelihood(stimuli|a^{x})')