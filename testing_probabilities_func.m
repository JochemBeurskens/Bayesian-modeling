clear all;
clc;

t_max=1; %set number of time steps in the process
n=500; %maybe not discritize, set the number of possible locations, the larger the more accurate, but also the more computation time is required
n_pf = 500; %number of particles that are used in the filter
c=0.8; %setting a parameter for the resampling, as this is done when N_eff<=c*n, setting c=1 means at every step, setting c=0 means you never resample
k=0.3; %setting the multiplication for the linear functions
lrange=[-20 20]; %the range of l, actual locations are between -10.5 and 10.5
l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
sig_l_s=sqrt(10);
sig_l_vs=sqrt(10);
sig_l_as=sig_l_vs; %when assuming that the original location is drawn from one normal distribution, only one sigma_L^s is needed, same for c. They can now be seperate draws though, which for C=1 was not possible
sig_lax=sqrt(100);
sig_lvx=sqrt(10);
sig_e_s=sqrt(10);
sig_e_ax=sqrt(10);
sig_e_vx=sqrt(10);
p_c=0.5;
p_h=0.2;
%% now checking that the sum equals 1
% generating a space to check the probilities over, these represent the
% observations x_a and x_v
L_xar=linspace(-20,20,n);
L_xvr=linspace(-20,20,n);
lrange=[-20 20];
P_LC1r=zeros(length(L_xar),length(L_xvr));
P_LC2r=zeros(length(L_xar),length(L_xvr));
% using the set parameters to see whether the probabilities indeed add up
% to 1 in the whole range. See the plots for one specific L_Xa or L_Xv
for i=1:length(L_xar)
    L_xa=L_xar(i);
    for o=1:length(L_xvr)
        L_xv=L_xvr(o);
%         [P_LC1,P_LC2]=Normalised_Probabilities_multimodal_func(k,t_max,n,l,sig_l_s,sig_l_as,sig_l_vs,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,L_xa,L_xv);
 
        [P_LC1,P_LC2]=Normalised_Probabilities_func(sig_l_s,sig_lvx,sig_lax,L_xv,L_xa);
        P_LC1r(i,o)=P_LC1;
        P_LC2r(i,o)=P_LC2;
    end
end
%% plotting
lrange=[-20 20];
figure(1)
plot(P_LC2r(250,:))
hold on
plot(P_LC1r(250,:))
legend("C=2","C=1")
figure(2)
plot(P_LC2r(:,250))
hold on
plot(P_LC1r(:,250))
legend("C=2","C=1")
figure(3)
% surf(L_xar,L_xvr,P_LC1r,'FaceColor','b')%[0,158/255,115/255]);
[c,h] = contourm(L_xar,L_xvr,P_LC1r);
hold on
% contourm(L_xar,L_xvr,P_LC2r)
plot(lrange,lrange)
axis square
ylim(lrange); 
xlim(lrange); 
xticks([-50 -40 -30 -20 -10 0 10 20 30 40 50])
yticks([-50 -40 -30 -20 -10 0 10 20 30 40 50])
clegendm(c,h,4)
% hold on; 
% contour(L_xar,L_xvr,P_LC2r,'FaceColor','g')%[230/255,159/255,0]);
xlabel('L_{a}^{x}')
ylabel('L_{v}^{x}')
zlabel('Probability')
% legend('P(C=1)','P(C=2)')
%% now obtaining a repsonse distribution for the case where the meaning is human human, and the location is -3.5 degrees
% n=10000; %number of draws to be used
% resp_loc_predicted=zeros(n,1);
% l_vals=[-10.5,-3.5,3.5,10.5];
% for i=1:n
%     L_xa=normrnd(0,sig_l_s);
%     L_xv=normrnd(0,sig_l_s);
%     [P_LC1,P_LC2]=Normalised_Probabilities_multimodal_func(k,t_max,1,l,sig_l_s,sig_l_as,sig_l_vs,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,L_xa,L_xv);
%     resp_loc_predicted_LC1=((((L_xa/(sig_lax^2))+L_xv/(sig_lvx^2)))/((1/sig_lax^2)+(1/sig_lvx^2)+(1/sig_l_s^2)))*P_LC1;
%     resp_loc_predicted_LC2=(((L_xa/(sig_lax^2))))/((1/sig_lax^2)+(1/sig_l_as^2))*P_LC2;
%     resp_loc_predicted(i)=resp_loc_predicted_LC1+resp_loc_predicted_LC2;
%     [~,d_l_vals]=dsearchn(resp_loc_predicted(i),l_vals');
%     i_l_vals=find(d_l_vals==min(d_l_vals(:)));
%     resp_loc_predicted(i)=i_l_vals;
% end
% c_1=length(find(resp_loc_predicted==1))/n;
% c_2=length(find(resp_loc_predicted==2))/n;
% c_3=length(find(resp_loc_predicted==3))/n;
% c_4=length(find(resp_loc_predicted==4))/n;
% bar([c_1 c_2 c_3 c_4])
% ylabel('Fraction of responses')
% xticklabels(["-10.5^{\circ}" "-3.5^{\circ}" "3.5^{\circ}" "10.5^{\circ}"])