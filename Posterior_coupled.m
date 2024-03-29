clear all;
clc;
% rng(3); %setting the random number generator seed, for reproducible results

%% parameter specification
t_max=1; %set number of time steps in the process
n=10; %maybe not discritize, set the number of possible locations, the larger the more accurate, but also the more computation time is required
m_n=10;
n_r=1; %number of stes of a_x or v_x obtained for each set of L_s, M_s and c_s
n_lm=1; %number of sets of L^x and M^x that are considered for the integrals
c=0.8; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
k=1; %setting the multiplication for the linear functions
l_r=40;
lrange=[-l_r/2 l_r/2]; %the range of l, actual locations are between -10.5 and 10.5
l=linspace(lrange(1),lrange(2),n);%lrange;%linspace(1,n,n); %setting the range of locations and meanings
m_r=2;
mrange=[-m_r m_r]; %the range of l, actual locations are between -10.5 and 10.5
m=linspace(mrange(1),mrange(2),m_n); %setting the range of locations and meanings

sig_l_s=1;
sig_l_vs=1;
sig_l_as=sig_l_vs; %when assuming that the original location is drawn from one normal distribution, only one sigma_L^s is needed, same for c. They can now be seperate draws though, which for C=1 was not possible
sig_lax=10;
sig_lvx=1;
sig_e_s=1;
sig_e_ax=1;
sig_e_vx=1;
sig_m_s=1;
sig_m_as=1;
sig_m_vs=sig_m_as;
sig_m_ax=1;
sig_m_vx=1;
p_c=0.5; %if x param > p_c: m=1, else m=2
p_h=1;

n_c_s=10;
diff_mean=3;
mean_rw=10; %setting the mean of the random walk
crange=[mean_rw-diff_mean mean_rw+diff_mean]; %the range of l, actual locations are between -10.5 and 10.5
if n_c_s == 1
    c_l=mean_rw;
else
    c_l=linspace(crange(1),crange(2),n_c_s); %setting the range of locations and meanings
end
%% obtaining the a_x, v_x values
[i_la_plt,i_lv_plt,m_a,m_v,a_x,v_x,i_l,c_av_s,m_av_xa,m_av_xv,m_s]=gen_model_full(k,t_max,l,m,sig_l_s,sig_lvx,sig_lax,sig_m_s,sig_m_as,sig_m_vs,sig_m_ax,sig_m_vx,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,n_r,n_lm,mean_rw);
% a_x=[   -1.9161   -0.1085    0.5103   -1.6292
%     1.0938    0.0892    0.4049    0.9844
%     0.9981   -1.4256    1.8917   -0.9620
%     0.8695   -1.1801    0.1592    1.4307
%     0.4572    0.2975    0.2789    0.7485
%     1.4196   -0.7699    0.3905   -1.6745
%     1.4214   -0.1926   -0.2686    1.5724
%     1.2844    0.8921    0.9381    0.6117
%     0.1610   -0.5772    0.6329    0.0262
%     0.7651   -3.3331   -1.9541    0.2779
%    -0.6937    0.5940   -0.7639    0.6489
%     0.7623   -2.1806   -1.5235    0.4395
%    -0.1942   -0.6217   -0.3277   -0.3331
%     0.1585    0.0913   -0.4925    1.1969
%    -2.8892   -0.6512   -1.5526    1.3257
%    -0.8330    0.4425    1.2924    0.1804
%    -0.0976    0.6613   -0.2145    0.6539
%     0.2258   -0.3354   -0.0225   -1.1352
%     0.7087   -0.1652    1.4996    0.0943
%     1.5971    3.4140   -0.9782   -2.2866
%     0.1959   -0.0465   -0.6613    0.9510
%     0.7182    0.8696    0.4592    1.5024
%    -0.5838   -0.9141   -0.2967    0.5192
%     0.4739   -1.6043   -0.7773   -1.0583
%     0.0582   11.8533    0.5037    1.3847
%     1.1070    0.7473   -0.8766    0.5830
%     0.0844   -1.0890   -1.9064   -0.6763
%    -2.4617   -2.1572   -0.5018   -0.5674
%     0.6605   -0.5471   -0.0163   -1.3161
%     0.3431    1.3117   -0.8473    1.0041
%     0.1683   -0.2427   -1.5496    0.8436
%     1.8458   -0.0254    0.3465   -0.5012
%     0.9162   -0.9471    0.3484    0.2725
%    -0.3857   -0.8793   -1.3805    0.0452
%    -0.1856   -2.0659   -0.9295    1.2524
%    -1.0100   -1.0801   -0.5971   -1.0798
%     0.2539   -0.3938   -1.5369   -0.6221
%    -1.7619    0.6442   -1.4600   -0.8555
%     0.0402   -0.7911    0.5240    0.3886
%     1.2653    2.0725   -0.6216    1.2173
%     0.3737    1.4913   -0.5762   -1.2599
%    -0.1230   -0.5632    0.6065    0.3262
%    -1.0139   -0.9937    1.4928    0.1216
%    -1.4826   -1.1763   -2.0269    0.1340
%    -0.7559   -2.2200    0.5662   -0.6066
%     1.0752    0.3568    0.8656   -0.7601
%    -0.1220    0.2016   -1.0297    0.6997
%    -1.9932   -0.2445   -1.9596   -0.6715
%     0.1679    1.2502    0.7514    0.4160
%     0.3652    2.1694    0.1697   -1.9544];
% figure(1)
% % surf
% imagesc(a_x(:,1))
% xlabel('Meaning'), ylabel('Location')%, zlabel('a^{x}')
% figure(2)
% surf(v_x)
% xlabel('Meaning'), ylabel('Location'), zlabel('v^{x}')
%% we look at one point only first
prob_ax=zeros(n,m_n,n_c_s);%m_n);
prob_vx=zeros(n,m_n,n_c_s);%m_n); %we have the dimensions l,m,t,L^{x},M^{x}
prob_avx=zeros(n,m_n,n_c_s);%m_n); %we have the dimensions l,m,t,L^{x},M^{x}

prob_ax_store=zeros(n,m_n);%m_n);
prob_vx_store=zeros(n,m_n);%m_n); %we have the dimensions l,m,t,L^{x},M^{x}

for l_s_i=1:n    
    L_s = l(l_s_i);
    for m_s_i=1:m_n
        M_s = m(m_s_i);    
        prob_ax_t=ones(n,m_n,1);
        prob_vx_t=ones(n,m_n,1);    
%         prob_avx_t=ones(n,m_n,1);    
        for c_s_i=1:n_c_s    
             mean_val=c_l(c_s_i);%c_s, a_x is observed so cannot loop over this, also should not as only integrate/take into account the c_s
%     prob_c_s=abs(normcdf(c_l-((2*diff_mean)/n_c_s),mean_val,sig_e_s)-normcdf(c_l+((2*diff_mean)/n_c_s),mean_val,sig_e_s));   
            for l_x_i=1:n
                L_x=l(l_x_i);
                prob_Lax=abs(normcdf(L_x-(l_r/n),L_s,sig_lax)-normcdf(L_x+(l_r/n),L_s,sig_lax));%
                prob_Lvx=abs(normcdf(L_x-(l_r/n),L_s,sig_lvx)-normcdf(L_x+(l_r/n),L_s,sig_lvx));%
                for m_x_i=1:m_n
                    M_x=m(m_x_i);
                    prob_Max=abs(normcdf(M_x-(m_r/m_n),M_s,sig_m_ax)-normcdf(M_x+(m_r/m_n),M_s,sig_m_ax));%      
                    prob_Mvx=abs(normcdf(M_x-(m_r/m_n),M_s,sig_m_vx)-normcdf(M_x+(m_r/m_n),M_s,sig_m_vx));%                                    
                    for l_i=1:n
                        for m_i=1:m_n
                            if l_x_i==l_i && m_x_i==m_i                                
                                prob_ax_t(l_i,m_i)=(normpdf(a_x(l_i,m_i),mean_val,sig_e_s));
                                prob_vx_t(l_i,m_i)=(normpdf(v_x(l_i,m_i),mean_val,sig_e_s)); 
                            else
                                prob_ax_t(l_i,m_i)=(normpdf(a_x(l_i,m_i),0,sig_e_s));   
                                prob_vx_t(l_i,m_i)=(normpdf(v_x(l_i,m_i),0,sig_e_s));       
                            end
                        end
%                         prob_avx_t(l_i,m_i)=prob_ax_t(l_ia,m_ia)*prob_vx_t(l_iv,m_iv);                                    
                    end                
                    prob_m(m_x_i)=prob_Max;
                    prob_mv(m_x_i)=prob_Mvx;
                    prob_l(l_x_i)=prob_Lax;
                    prob_lv(l_x_i)=prob_Lvx;   
                    product_p_ax(l_x_i,m_x_i,c_s_i)=prob_Max*prob_Lax*(prod(prod(prob_ax_t(:,:),2),1));
                    product_p_vx(l_x_i,m_x_i,c_s_i)=prob_Mvx*prob_Lvx*(prod(prod(prob_vx_t(:,:),2),1));
                    %the above is thus calculated for each combination of
                    %L^{x},M^{x},c^{s}
                end
            end 
            prob_ax(l_s_i,m_s_i,c_s_i)=sum(sum(product_p_ax(:,:,c_s_i),2),1);    %prob_c_s(c_s_i)*
            prob_vx(l_s_i,m_s_i,c_s_i)=sum(sum(product_p_vx(:,:,c_s_i),2),1);      %prob_c_s(c_s_i)*  
            prob_avx(l_s_i,m_s_i,c_s_i)=((prob_ax(l_s_i,m_s_i,c_s_i).*prob_vx(l_s_i,m_s_i,c_s_i)));
            %the above is calculated for each combination of L^{s},M^{s},c^{s}   
            prob_c_s=normpdf(mean_val,mean_rw,sig_e_s);
            prob_l_s=normpdf(L_s,0,sig_l_s);
            prob_m_s=normpdf(M_s,0,sig_m_s);
            post_ax(l_s_i,m_s_i,c_s_i)=prob_c_s*prob_l_s*prob_m_s*prob_ax(l_s_i,m_s_i,c_s_i);
            post_avx(l_s_i,m_s_i,c_s_i)=prob_c_s*prob_l_s*prob_m_s*prob_avx(l_s_i,m_s_i,c_s_i);
        end
    end
end
prob_ax_plot=sum(prob_ax(:,:,:),3);
% prob_ax_ploti=sum(prob_ax(:,:,:),2);    
% prob_ax_plotii=sum(prob_ax(:,:,:),1);    
% 
% for i=1:n
%     ploti(i,:)=prob_ax_ploti(:,i)
%     plotii(i,:)=prob_ax_plotii(:,:,i)
% end
prob_vx_plot=sum(prob_vx(:,:,:),3);
prob_avx_plot=sum(prob_avx(:,:,:),3);

%normalising the posterior
post_avx=post_ax/(sum(sum(sum(post_avx))));
post_avx_plot=sum(post_avx,3);
% post_ax_ploti=sum(post_ax(:,:,:),2);    
% 
% for i=1:n
%     ploti(i,:)=post_ax_ploti(:,i);
% end
if m_n==1 || n==1
    figure; plot(prob_ax_plot);%
    figure; plot(a_x);
    figure; plot(prob_vx_plot);%
    figure; plot(v_x);
    figure; plot(prob_avx_plot);%
else
    figure; imagesc(l,m,prob_ax_plot');xlabel("l"); ylabel("m");%
%     figure; imagesc(l,c_l,ploti);xlabel("l"); ylabel("c^{s}");%
%     figure; imagesc(m,c_l,plotii);xlabel("m"); ylabel("c^{s}");%
    figure; imagesc(l,m,a_x'); xlabel("l"); ylabel("m"); %zlabel("c^{s}");
    figure; imagesc(l,m,prob_vx_plot');xlabel("l"); ylabel("m");%
    figure; imagesc(l,m,v_x'); xlabel("l"); ylabel("m");
    figure; imagesc(l,m,prob_avx_plot'); xlabel("l"); ylabel("m");
    figure; imagesc(l,m,post_avx_plot');xlabel("l"); ylabel("m");%
end