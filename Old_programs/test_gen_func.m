k=0.3; %setting the multiplication for the linear functions
t_max=10; %setting the amount of time steps for the source activity
n=10; %setting the amount of possible locations

lrange=[-5 5];
mrange=[-1 1];

sig_e_s=0.1;
sig_l_s=0.6;
sig_m_s=0.1;%setting noise parameters

sig_e_as=0.1;
sig_e_vs=0.1;
sig_l_as=0.6;
sig_m_as=0.1;
sig_l_vs=0.6;
sig_m_vs=0.1;

l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
m=linspace(mrange(1),mrange(2),n);

sig_e_ax=0.1;
sig_e_vx=0.1;
sig_la=0.5;
sig_ma=0.5;
sig_lv=0.5;
sig_mv=0.5;

[a_s,v_s,c_av_s,i_l_plt,i_m_plt,i_la_plt,i_ma_plt,i_lv_plt,i_mv_plt,a_x,v_x]=func_genmod(k,t_max,n,l,m,sig_e_s,sig_l_s,sig_m_s,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);
% func_genmod(r,k,t_max,n,l,m,sig_e_as,sig_e_vs,sig_l_as,sig_l_vs,sig_m_as,sig_m_vs,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv)
figure(1);
scatter(i_l_plt,i_m_plt); 
hold on;
scatter(i_la_plt,i_ma_plt);
scatter(i_lv_plt,i_mv_plt); %scatter plotting the indices to compare the observations to the actual meaning/location
xlabel('location');
ylabel('semantic meaning');
axis([lrange mrange])
legend('actual','auditory percept','visual percept')