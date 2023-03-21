clear all;
clc;
rng(8); %setting the random number generator seed, for reproducible results

%% parameter specification
t_max=1; %set number of time steps in the process
n=10; %maybe not discritize, set the number of possible locations, the larger the more accurate, but also the more computation time is required
m_n=2;
n_r=1; %number of stes of a_x or v_x obtained for each set of L_s, M_s and c_s
n_lm=1; %number of sets of L^x and M^x that are considered for the integrals
c=0.8; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
k=1; %setting the multiplication for the linear functions
% lrange=[1 40]; %the range of l, actual locations are between -10.5 and 10.5
l=linspace(1,n,n); %setting the range of locations and meanings
mrange=[1 2]; %the range of l, actual locations are between -10.5 and 10.5
m=linspace(1,m_n,m_n); %setting the range of locations and meanings

sig_l_s=5;
sig_l_vs=5;
sig_l_as=sig_l_vs; %when assuming that the original location is drawn from one normal distribution, only one sigma_L^s is needed, same for c. They can now be seperate draws though, which for C=1 was not possible
sig_lax=8;
sig_lvx=5;
sig_e_s=1;
sig_e_ax=1;
sig_e_vx=1;
sig_m_s=2;
sig_m_as=2;
sig_m_vs=sig_m_as;
sig_m_ax=2;
sig_m_vx=2;
p_c=0.5; %if x param > p_c: m=1, else m=2
p_h=1;
mean_rw=10; %setting the mean of the random walk
%% obtaining the a_x, v_x values
[i_la_plt,i_lv_plt,m_a,m_v,a_x,v_x,i_l,c_av_s,m_av_xa,m_av_xv,m_s]=gen_model_full(k,t_max,l,m,sig_l_s,sig_lvx,sig_lax,sig_m_s,sig_m_as,sig_m_vs,sig_m_ax,sig_m_vx,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,n_r,n_lm,mean_rw);

figure(1)
% surf
imagesc(a_x(:,1))
xlabel('Meaning'), ylabel('Location')%, zlabel('a^{x}')
figure(2)
surf(v_x)
xlabel('Meaning'), ylabel('Location'), zlabel('v^{x}')
%% we look at one point only first
prob_ax=ones(n,1);%m_n);
prob_ax_t=ones(n,1);
prob_vx=zeros(n,1);%m_n); %we have the dimensions l,m,t,L^{x},M^{x}
%this will be :,:,t,p,
%v_x and a_x have dimensions n,2
pc=normcdf(p_c,m_s,sig_m_ax);
% p_c=[1-pc pc]; %prob for m=1 and m=2, note that if the m_x draw is larger than p_c it is put in m=1, thus 1-normcdf is the prob
% a_x=a_x(1,1,1);
% prob_ax(3,4)=1
a_x_val=a_x(:,2);
% L_x=3;
% L_s=2;
figure
for l_s_i=1:n
    L_s=l(l_s_i);
    for l_x_i=1:n
        L_x=l(l_x_i);
        prob_Lax=abs(normcdf(L_x-0.5,L_s,sig_lax)-normcdf(L_x+0.5,L_s,sig_lax));%normpdf(L_x,L_s,sig_lax);%
        for l_i=1:n
            if L_x==l_i                      
        %         prob_Lax=abs(normcdf(L_x-0.5,L_s,sig_lax)-normcdf(L_x+0.5,L_s,sig_lax));%normpdf(L_x,L_s,sig_lax);%
                prob_ax(l_x_i)=prob_ax(l_x_i)*(normpdf(a_x_val(l_i),mean_rw,sig_e_s));
            else
        %         prob_Lax=(1-abs(normcdf(L_x-0.5,L_s,sig_lax)-normcdf(L_x+0.5,L_s,sig_lax)));%normpdf(L_x,L_s,sig_lax);%
                prob_ax(l_x_i)=prob_ax(l_x_i)*(normpdf(a_x_val(l_i),0,sig_e_s));       
            end
        end
        prob_ax_t(l_x_i)=prob_Lax*prob_ax_t(l_x_i);
        plot(prob_ax_t); hold on;
    end
    prob_ax(l_s_i)=sum(prob_ax_t);
end
hold off;

figure; plot(prob_ax)
% for l_i=1:n %loop over l
%     L_x=l(l_i);
%     for o=1:n %the loop over L_s
%         L_s=l(o);
%         if L_x==l_i                      
%             prob_Lax=abs(normcdf(L_x-0.5,L_s,sig_lax)-normcdf(L_x+0.5,L_s,sig_lax));%normpdf(L_x,L_s,sig_lax);%
%             prob_ax(:,o)=prob_ax(l_i,o)+(prob_Lax*normpdf(a_x_val,mean_rw,sig_e_s));
%         else
%             prob_Lax=(1-abs(normcdf(L_x-0.5,L_s,sig_lax)-normcdf(L_x+0.5,L_s,sig_lax)));%normpdf(L_x,L_s,sig_lax);%
%             prob_ax(:,o)=prob_ax(l_i,o)+(prob_Lax*normpdf(a_x_val,0,sig_e_s));       
%         end
%     end
% end

% L_s=3;
% L_x=1.4;
% [~,d_l_x]=dsearchn(L_x,l');
% l_i=find(d_l_x==min(d_l_x(:)));
% 
% prob_Lax=abs(normcdf(L_x-0.5,L_s,sig_lax)-normcdf(L_x+0.5,L_s,sig_lax));%normpdf(L_x,L_s,sig_lax);%
% prob_ax(l_i,1)=prob_ax(l_i,1)+(prob_Lax*normpdf(a_x_val(l_i),mean_rw,sig_e_s));
% prob_Lax=(1-abs(normcdf(L_x-0.5,L_s,sig_lax)-normcdf(L_x+0.5,L_s,sig_lax)));%normpdf(L_x,L_s,sig_lax);%
% prob_ax(l_i,1)=prob_ax(l_i,1)+(prob_Lax*normpdf(a_x_val(l_i),0,sig_e_s));       

% figure;
% surf(prob_ax);xlabel('L^{s}'), ylabel('l'); %this plots the row on x and the column on y
%% we look at one array
% prob_ax=zeros(n,n,m_n,m_n);%m_n);
% prob_vx=zeros(n,n,m_n,m_n);%m_n); %we have the dimensions l,m,t,L^{x},M^{x}
% %this will be :,:,t,p,
% %v_x and a_x have dimensions n,2
% pc=normcdf(p_c,m_s,sig_m_ax);
% % p_c=[1-pc pc]; %prob for m=1 and m=2, note that if the m_x draw is larger than p_c it is put in m=1, thus 1-normcdf is the prob
% % a_x=a_x(1,1,1);
% for o=1:n %the loop over L_s
%     L_s=o;
% %     for o_m=1:m_n %the loop over M_s
% %         M_s=o_m;
%         for l_i=1:n %loop over l
% %             val_l=l(l_i);
%             for m_i=1:m_n %loop over m
% %                 val_m=m(m_i);
%                 a_x_val=a_x(l_i,m_i);
%                 for p=1:n %the loop over L_x
%                     L_x=p;
%                     for p_m=1:m_n %the loop over M_x
%                         M_x=p_m;
%                         
%                             if L_x==l_i %&& M_x==m_i                        
% %                                 prob_Max=(normcdf(M_x-0.5,M_s,sig_m_ax)-normcdf(M_x+0.5,M_s,sig_m_ax));%normpdf(M_x,M_s,sig_m_ax);
%                                 prob_Lax=abs(normcdf(L_x-0.5,L_s,sig_lax)-normcdf(L_x+0.5,L_s,sig_lax));%normpdf(L_x,L_s,sig_lax);%
% %                                 prob_ax(L_s,M_s)=prob_ax(L_s,M_s)+(prob_Max*prob_Lax*normpdf(a_x_val,mean_rw,sig_e_s));
%                                 prob_ax(l_i,L_s,m_i,M_x)=prob_ax(l_i,L_s,m_i,M_x)+(prob_Lax*normpdf(a_x_val,mean_rw,sig_e_s));
%                             else
% %                                 prob_Max=(normcdf(M_x-0.5,M_s,sig_m_ax)-normcdf(M_x+0.5,M_s,sig_m_ax));%normpdf(M_x,M_s,sig_m_ax);
%                                 prob_Lax=abs(normcdf(L_x-0.5,L_s,sig_lax)-normcdf(L_x+0.5,L_s,sig_lax));%normpdf(L_x,L_s,sig_lax);%
% %                                 prob_ax(L_s,M_s)=prob_ax(L_s,M_s)+(prob_Max*prob_Lax*normpdf(a_x_val,0,sig_e_s));
%                                 prob_ax(l_i,L_s,m_i,M_x)=prob_ax(l_i,L_s,m_i,M_x)+(prob_Lax*normpdf(a_x_val,0,sig_e_s));
%                             end
%         %                     if p==i_lv_plt
%         %                         if p_m==m_v
%         %                             prob_Mvx=normpdf(m(p_m),m(o_m),sig_m_vx);
%         %                             prob_Lvx=normpdf(l(p),l(o),sig_lvx);
%         %                             prob_vx(i_ll,i_mm,t,p,p_m,o,o_m)=prob_Mvx.*prob_Lvx'.*normpdf(v_x(i_ll,i_mm,t),mean_rw,sig_e_s); %v_x and prob_Lvx are build such that they already containt the array for l,m
%         %                         else 
%         %                             prob_Mvx=normpdf(m(p_m),m(o_m),sig_m_vx);
%         %                             prob_Lvx=normpdf(l(p),l(o),sig_lvx);
%         %                             prob_vx(i_ll,i_mm,t,p,p_m,o,o_m)=prob_Mvx.*prob_Lvx'.*normpdf(v_x(i_ll,i_mm,t),0,sig_e_s);
%         %                         end
%         %                     else 
%         %                         prob_Mvx=normpdf(m(p_m),m(o_m),sig_m_vx);
%         %                         prob_Lvx=normpdf(l(p),l(o),sig_lvx);
%         %                         prob_vx(i_ll,i_mm,t,p,p_m,o,o_m)=prob_Mvx.*prob_Lvx'.*normpdf(v_x(i_ll,i_mm,t),0,sig_e_s);
%                             end
%                         
%                     end
%                 end
%             end
%         end
% %     end
% % end 
% figure(4)
% surf(prob_ax)

% %% now obtaining the measurement prob. distribution for the data above (Causal=1)
% prob_ax=zeros(n,m_n,t_max,n,m_n,n,m_n);
% prob_vx=zeros(n,m_n,t_max,n,m_n,n,m_n); %we have the dimensions l,m,t,L^{x},M^{x},L^{s},M^{s}
% %this will be :,:,t,p,p_m,o,o_m
% %v_x and a_x have dimensions n,2
% pc=normcdf(p_c,m_s,sig_m_ax);
% % p_c=[1-pc pc]; %prob for m=1 and m=2, note that if the m_x draw is larger than p_c it is put in m=1, thus 1-normcdf is the prob
% for o=1:n %the loop over L_s
%     for o_m=m_n %the loop over M_s
%         for p=1:n %the loop over L_x
%             for p_m=1:m_n %the loop over M_x
%                 for i_ll=1:length(l) %the loop over l
%                     for i_mm=1:length(m) %the loop over m
%                         for t =1:t_max
%                             if p==i_la_plt
%                                 if p_m==m_a
%                                     prob_Max=normpdf(m(p_m),m(o_m),sig_m_ax);
%                                     prob_Lax=normpdf(l(p),l(o),sig_lax);
%                                     prob_ax(i_ll,i_mm,t,p,p_m,o,o_m)=prob_Max.*prob_Lax'.*normpdf(a_x(i_ll,i_mm,t),mean_rw,sig_e_s);
%                                 else 
%                                     prob_Max=normpdf(m(p_m),m(o_m),sig_m_ax);
%                                     prob_Lax=normpdf(l(p),l(o),sig_lax);
%                                     prob_ax(i_ll,i_mm,t,p,p_m,o,o_m)=prob_Max.*prob_Lax'.*normpdf(a_x(i_ll,i_mm,t),0,sig_e_s);
%                                 end
%                             else 
%                                 prob_Max=normpdf(m(p_m),m(o_m),sig_m_ax);
%                                 prob_Lax=normpdf(l(p),l(o),sig_lax);
%                                 prob_ax(i_ll,i_mm,t,p,p_m,o,o_m)=prob_Max.*prob_Lax'.*normpdf(a_x(i_ll,i_mm,t),0,sig_e_s);
%                             end
%                             if p==i_lv_plt
%                                 if p_m==m_v
%                                     prob_Mvx=normpdf(m(p_m),m(o_m),sig_m_vx);
%                                     prob_Lvx=normpdf(l(p),l(o),sig_lvx);
%                                     prob_vx(i_ll,i_mm,t,p,p_m,o,o_m)=prob_Mvx.*prob_Lvx'.*normpdf(v_x(i_ll,i_mm,t),mean_rw,sig_e_s); %v_x and prob_Lvx are build such that they already containt the array for l,m
%                                 else 
%                                     prob_Mvx=normpdf(m(p_m),m(o_m),sig_m_vx);
%                                     prob_Lvx=normpdf(l(p),l(o),sig_lvx);
%                                     prob_vx(i_ll,i_mm,t,p,p_m,o,o_m)=prob_Mvx.*prob_Lvx'.*normpdf(v_x(i_ll,i_mm,t),0,sig_e_s);
%                                 end
%                             else 
%                                 prob_Mvx=normpdf(m(p_m),m(o_m),sig_m_vx);
%                                 prob_Lvx=normpdf(l(p),l(o),sig_lvx);
%                                 prob_vx(i_ll,i_mm,t,p,p_m,o,o_m)=prob_Mvx.*prob_Lvx'.*normpdf(v_x(i_ll,i_mm,t),0,sig_e_s);
%                             end
%                         end
%                     end
%                 end
%             end
%         end 
%     end
% end
% figure(3)
% surf(sum(sum(sum(sum(prob_ax(:,:,:,:,:,:,:),7),6),5),4))
% % figure(1)
% surf(sum(sum(prob_ax(:,:,:,:,:),5),4))
% xlabel('Meaning'), ylabel('Location'), zlabel('sum of P(a^{x})')
% figure(2)
% surf(sum(sum(prob_vx(:,:,:,:,:),5),4))
% xlabel('Meaning'), ylabel('Location'), zlabel('sum of P(v^{x})')