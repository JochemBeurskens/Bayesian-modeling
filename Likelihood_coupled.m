clear all;
clc;
% rng(3); %setting the random number generator seed, for reproducible results

%% parameter specification
t_max=1; %set number of time steps in the process
n=100; %maybe not discritize, set the number of possible locations, the larger the more accurate, but also the more computation time is required
m_n=1; 
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
sig_lax=4;
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
mean_rw=0; %setting the mean of the random walk
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
prob_ax=zeros(n,m_n);%m_n);
prob_vx=zeros(n,m_n);%m_n); %we have the dimensions l,m,t,L^{x},M^{x}
%this will be :,:,t,p,
%v_x and a_x have dimensions n,2
% pc=normcdf(p_c,m_s,sig_m_ax);
% p_c=[1-pc pc]; %prob for m=1 and m=2, note that if the m_x draw is larger than p_c it is put in m=1, thus 1-normcdf is the prob
% a_x=a_x(1,1,1);
% prob_ax(3,4)=1

%when observation is generated from 0
% mean_rw=0;
% a_x_val=[-0.8842 1.0951 -0.3416 0.7105 0.5367 2.1711 -1.3945 -0.1956 -0.0428 -0.0646];%a_x(:,1);

%when observation is generated from 10
% mean_rw=10;
% a_x_val=[-0.8842 10.0951 -0.3416 0.7105 0.5367 2.1711 -1.3945 -0.1956 -0.0428 -0.0646];%a_x(:,1);
% L_x=3;
% L_s=2;
% figure
% for l_s_i=1:n
%     L_s=l(l_s_i);
% %     disp(L_s)
%     prob_ax_t=ones(n,1);
%     prob_Lax=ones(n,1);
%     for l_x_i=1:n
%         L_x=l(l_x_i);
%         prob_Lax(L_x)=normpdf(L_x,L_s,sig_lax);%abs(normcdf(L_x-0.5,L_s,sig_lax)-normcdf(L_x+0.5,L_s,sig_lax));%
%         for l_i=1:n
%             if L_x==l_i 
% %                 disp(L_s)
% %                 disp("and")                
% %                 disp(L_x);
% %                 disp("and")
% %                 disp(normpdf(a_x_val(l_i),mean_rw,sig_e_s));
%                 prob_ax_t(l_i)=(normpdf(a_x_val(l_i),mean_rw,sig_e_s));
%             else
%                 prob_ax_t(l_i)=(normpdf(a_x_val(l_i),0,sig_e_s));       
%             end
%         end
% %         scatter(L_x,prob_Lax);
% %         scatter(L_x,(prob_ax_t(l_x_i))); hold on;
%         disp(prob_ax_t)
% %         prob_ax_t=prob_Lax*prod(prob_ax_t);
%         prob_ax(L_s)=sum(prob_Lax*prod(prob_ax_t));
%     end
% end
% m_n=1
% a_x=a_x(:,2)
for l_s_i=1:n

    L_s = l(l_s_i);
    for m_s_i=1:m_n
        M_s = m(m_s_i);    
        prob_ax_t=ones(n,m_n,1);

    %     disp(L_s)
    %     prob_Lax=ones(n,1);
        for l_x_i=1:n
            L_x=l(l_x_i);
            prob_Lax=abs(normcdf(L_x-(l_r/n),L_s,sig_lax)-normcdf(L_x+(l_r/n),L_s,sig_lax));%
            for m_x_i=1:m_n
                M_x=m(m_x_i);
                prob_Max=abs(normcdf(M_x-(m_r/m_n),M_s,sig_m_ax)-normcdf(M_x+(m_r/m_n),M_s,sig_m_ax));%                
                    for l_i=1:n
                        for m_i=1:m_n
                            if l_x_i==l_i && m_x_i==m_i                                
                %                 disp(L_s)
                %                 disp("and")                
                %                 disp(L_x);
                %                 disp("and")
                %                 disp(normpdf(a_x_val(l_i),mean_rw,sig_e_s));
                                prob_ax_t(l_i,m_i)=(normpdf(a_x(l_i,m_i),mean_rw,sig_e_s));
                            else
                                prob_ax_t(l_i,m_i)=(normpdf(a_x(l_i,m_i),0,sig_e_s));                                
                            end
                        end
                    end
        %         scatter(L_x,prob_Lax);
        %         scatter(L_x,(prob_ax_t(l_x_i))); hold on;
%                 disp(prob_ax_t)
        %         prob_ax_t=prob_Lax*prod(prob_ax_t);
%                 product(:)=(prod(prob_ax_t(:,:),2));
            prob_m(m_x_i)=prob_Max;
            prob_l(l_x_i)=prob_Lax;
            product_p_ax(l_x_i,m_x_i)=prob_Max*prob_Lax*(prod(prod(prob_ax_t(:,:),2)));
            end
%             product_ax(l_x_i,m_x_i)=prod(prod(prob_ax_t(:,:),2));
        end
%         prob_ax_s=sum(product_s);
        prob_ax(l_s_i,m_s_i)=sum(sum(product_p_ax,2),1);    
    end
end

% hold off;
if m_n==1 || n==1
    figure; plot(prob_ax);%
    figure; plot(a_x);
else
    figure; imagesc(prob_ax);%
    figure; imagesc(a_x);
end