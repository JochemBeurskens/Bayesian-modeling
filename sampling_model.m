% in this code samples from the generative model are obtained, each of the
% different types of data (loc,meaning,content) is sampled separately for
% clarity.
% rng(5); % set random seed

%% setting a random seed and some parameters
k=0.3; %setting the multiplication for the linear functions
t_max=100; %setting the amount of time steps for the source activity
n=20; %number of steps in the given range
lrange=[-5 5]; %the range of l
l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings

%% location of the stimuli for the case of C=1
i_la_plt=zeros(t_max,1);
i_lv_plt=zeros(t_max,1);

sig_l_s=1;
sig_lvx=1;
sig_lax=4;
l_i=[0 0 0 0]; %the location about which the original signal is centered
%generating the signal, denoted by s:   
L_av_s=normrnd(l_i(1),sig_l_s); %drawing a location and meaning from a prior, this is now a Gaussian
[~,d_l]=dsearchn(L_av_s,l'); %getting the distances from the observations to the possible locations
i_l=find(d_l==min(d_l(:))); %obtaining the index of the closest matching location
for i=1:t_max
    %making the observations, denoted by x:
    L_av_xa=normrnd(L_av_s,sig_lax);
    L_av_xv=normrnd(L_av_s,sig_lvx);
   
    [~,d_l_xa]=dsearchn(L_av_xa,l');
    [~,d_l_xv]=dsearchn(L_av_xv,l'); %getting the distances from the observations to the possible locations
    
    i_la=find(d_l_xa==min(d_l_xa(:)));
    i_lv=find(d_l_xv==min(d_l_xv(:))); %obtaining the index of the closest matching location
    
    i_la_plt(i)=l(i_la); %these values represent the result of the delta function from the formulas
    i_lv_plt(i)=l(i_lv); %these values represent the result of the delta function from the formulas
end
i_l=[l(i_l) l(i_l) l(i_l) l(i_l)];
%% the random walk for the content for the case of C=1
a_s=zeros(t_max,1);
v_s=zeros(t_max,1);
a_x=zeros(t_max,1);
v_x=zeros(t_max,1);
sig_e_s=0.5;
sig_e_ax=0.5;
sig_e_vx=0.5;
c_av_s=[normrnd(0,sig_e_s)]; %drawing c_av^s from the normal distribution
for i=1:t_max
    e_s=normrnd(0,sig_e_s);
    c_av_s=[c_av_s c_av_s(i)+e_s]; %storing and adding noise
    f_a_s=c_av_s(i);
    f_v_s=k*c_av_s(i); %applying the linear mapping functions
    a_s(i)=1*1*f_a_s;  
    v_s(i)=1*1*f_v_s; %generating the result, and storing it in the auditory/visual observation lists
    f_a_x=c_av_s(i)+normrnd(0,sig_e_ax);
    f_v_x=k*c_av_s(i)+normrnd(0,sig_e_vx); %using the same linear functions as before, but now noise is added
    a_x(i)=1*1*f_a_x;
    v_x(i)=1*1*f_v_x; %storing the auditory/visual amplituded in their respective meaing/location bins
end

%% meaning of the stimuli for the case of C=1
p_c=0.7;
p_h=0.8;
m=zeros(t_max,1);
m_i=[0 1 2 3];
for i=1:t_max
    m_c=binornd(1,p_c); %getting t_max draws from the binomial distribution determining whether the stimuli should be congruent
    m_h=binornd(1,p_h); %getting t_max draws from the binomial distribution determining whether the stimuli should be human
    if m_c==0 && m_h==0
        m(i)=0;
    elseif m_c==0 && m_h==1
        m(i)=1;
    elseif m_c==1 && m_h==0
        m(i)=2;
    elseif m_c==1 && m_h==1
        m(i)=3;
    end
end

%% plotting for the case of C=1
lr=[-12 12];
mr=[-1 4];
figure(1);
scatter(m,i_la_plt);
hold on;
scatter(m,i_lv_plt,'v'); %scatter plotting the indices to compare the observations to the actual meaning/location
scatter(m_i,i_l,'.'); 
xlabel('semantic meaning');
ylabel('location');
axis([mr lr])
legend('auditory percept', 'visual percept','actual location')

figure(2);
plot(c_av_s);
hold on;
plot(a_s);
plot(v_s);
plot(a_x);
plot(v_x);

%% location of the stimuli for the case of C=2
i_la_plt=zeros(t_max,1);
i_lv_plt=zeros(t_max,1);

sig_l_vs=1;
sig_l_as=4;
sig_lvx=1;
sig_lax=4;
l_i=[0 0 0 0]; %the location about which the original signal is centered
%generating the signal, denoted by s:   
L_av_vs=normrnd(l_i(1),sig_l_vs); %drawing a location and meaning from a prior, this is now a Gaussian
[~,d_lv]=dsearchn(L_av_vs,l'); %getting the distances from the observations to the possible locations
i_lvs=find(d_lv==min(d_lv(:))); %obtaining the index of the closest matching location

L_av_as=normrnd(l_i(1),sig_l_as); %drawing a location and meaning from a prior, this is now a Gaussian
[~,d_la]=dsearchn(L_av_as,l'); %getting the distances from the observations to the possible locations
i_las=find(d_la==min(d_la(:))); %obtaining the index of the closest matching location
for i=1:t_max
    %making the observations, denoted by x:
    L_av_xa=normrnd(L_av_as,sig_lax);
    L_av_xv=normrnd(L_av_vs,sig_lvx);
   
    [~,d_l_xa]=dsearchn(L_av_xa,l');
    [~,d_l_xv]=dsearchn(L_av_xv,l'); %getting the distances from the observations to the possible locations
    
    i_la=find(d_l_xa==min(d_l_xa(:)));
    i_lv=find(d_l_xv==min(d_l_xv(:))); %obtaining the index of the closest matching location
    
    i_la_plt(i)=l(i_la); %these values represent the result of the delta function from the formulas
    i_lv_plt(i)=l(i_lv); %these values represent the result of the delta function from the formulas
end
i_lvs=[l(i_lvs) l(i_lvs) l(i_lvs) l(i_lvs)];
i_las=[l(i_las) l(i_las) l(i_las) l(i_las)];

%% the random walk for the content for the case of C=2
a_s=zeros(t_max,1);
v_s=zeros(t_max,1);
a_x=zeros(t_max,1);
v_x=zeros(t_max,1);
sig_e_s=0.5;
sig_e_ax=0.5;
sig_e_vx=0.5;
c_av_s=[normrnd(0,sig_e_s)]; %drawing c_av^s from the normal distribution
for i=1:t_max
    e_s=normrnd(0,sig_e_s);
    c_av_s=[c_av_s c_av_s(i)+e_s]; %storing and adding noise
    f_a_s=c_av_s(i);
    if i >= t_max/2.5
        f_v_s=k*c_av_s(i); %applying the linear mapping functions
    else
        f_v_s=normrnd(0,0.05);
    end
    a_s(i)=1*1*f_a_s;  
    v_s(i)=1*1*f_v_s; %generating the result, and storing it in the auditory/visual observation lists
    f_a_x=f_a_s+normrnd(0,sig_e_ax);
    f_v_x=f_v_s+normrnd(0,sig_e_vx); %using the same linear functions as before, but now noise is added
    a_x(i)=1*1*f_a_x;
    v_x(i)=1*1*f_v_x; %storing the auditory/visual amplituded in their respective meaing/location bins
end

%% meaning of the stimuli for the case of C=2
p_c=0.3;
p_h=0.2;
m=zeros(t_max,1);
m_i=[0 1 2 3];
for i=1:t_max
    m_c=binornd(1,p_c); %getting t_max draws from the binomial distribution determining whether the stimuli should be congruent
    m_h=binornd(1,p_h); %getting t_max draws from the binomial distribution determining whether the stimuli should be human
    if m_c==0 && m_h==0
        m(i)=0;
    elseif m_c==0 && m_h==1
        m(i)=1;
    elseif m_c==1 && m_h==0
        m(i)=2;
    elseif m_c==1 && m_h==1
        m(i)=3;
    end
end

%% plotting for the case of C=2
lr=[-12 12];
mr=[-1 4];
figure(3);
scatter(m,i_la_plt);
hold on;
scatter(m,i_lv_plt,'v'); %scatter plotting the indices to compare the observations to the actual meaning/location
scatter(m_i,i_l,'.'); 
xlabel('semantic meaning');
ylabel('location');
axis([mr lr])
legend('auditory percept', 'visual percept','actual location')

figure(4);
plot(c_av_s);
hold on;
plot(a_s);
plot(v_s);
plot(a_x);
plot(v_x);