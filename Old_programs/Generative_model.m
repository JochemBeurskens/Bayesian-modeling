%% C=1
%the code below lets a person estimate a cue for three consecutive
%timesteps. For all of these the location of the target remains the same.

r=3; %how many times should the values of the meaning or location shift
i_l_plt=zeros(r,1);
i_m_plt=zeros(r,1);
i_la_plt=zeros(r,1);
i_ma_plt=zeros(r,1);
i_lv_plt=zeros(r,1);
i_mv_plt=zeros(r,1); %creating arrays to store the meaning/location value for later plotting

for y=1:r
    %generating the signal, denoted by s:
    k=0.3; %setting the multiplication for the linear functions
    t_max=10; %setting the amount of time steps for the source activity
    n=10; %setting the amount of possible locations
    
    l_low=-5;
    l_high=5;
    m_low=-2;
    m_high=2;

    sig_e_s=0.1;
    sig_l_s=0.6;
    sig_m_s=0.1;%setting noise parameters
    
    l=linspace(l_low,l_high,n); %setting the range of locations and meanings
    m=linspace(m_low,m_high,n); %creating two lists to better specify the location and meaning (note that this is made discrete here)
    
    L_av_s=normrnd(0,sig_l_s);
    M_av_s=normrnd(0,sig_m_s); %drawing a location and meaning from a prior, this is now a Gaussian, might be better to make this a discrete stepped function?
    [~,d_l]=dsearchn(L_av_s,l');
    [~,d_m]=dsearchn(M_av_s,m'); 
    i_l=find(d_l==min(d_l(:))); 
    i_m=find(d_m==min(d_m(:))); %finding the location and meaning closest to the ones drawn from the prior, these are then assumed as the meant location/meaning
    
    i_l_plt(y)=l(i_l);
    i_m_plt(y)=m(i_m); %storing the meaing/location for plotting

    a_s=zeros(n,t_max+1,n);
    v_s=zeros(n,t_max+1,n); %creating the auditory and visual observation lists, now still empty/zero
    
    c_av_s=[normrnd(0,sig_e_s)]; %drawing c_av^s from the normal distribution
    for i=1:t_max
        e_s=normrnd(0,sig_e_s);
        c_av_s=[c_av_s c_av_s(i)+e_s]; %storing and adding noise
    end
    f_a_s=c_av_s;
    f_v_s=k*c_av_s; %applying the linear mapping functions
    
    a_s(i_m,:,i_l)=1*1*f_a_s;
    v_s(i_m,:,i_l)=1*1*f_v_s; %generating the result, and storing it in the auditory/visual observation lists
    
    %making the observations, denoted by x:
    sig_e_ax=0.1;
    sig_e_vx=0.1;
    sig_la=0.5;
    sig_ma=0.5;
    sig_lv=0.5;
    sig_mv=0.5; %setting some parameters
    
    L_av_xa=normrnd(L_av_s,sig_la);
    L_av_xv=normrnd(L_av_s,sig_lv);
    M_av_xa=normrnd(M_av_s,sig_ma);
    M_av_xv=normrnd(M_av_s,sig_mv); %generating the meaning and location in the observation, based on a gaussian centered on the actual meaning/location and the above specified noise levels
    
    [~,d_l_xa]=dsearchn(L_av_xa,l');
    [~,d_m_xa]=dsearchn(M_av_xa,m');
    [~,d_l_xv]=dsearchn(L_av_xv,l');
    [~,d_m_xv]=dsearchn(M_av_xv,m'); %getting the distances from the observations to the possible locations/meanings, l and m
    
    i_la=find(d_l_xa==min(d_l_xa(:)));
    i_ma=find(d_m_xa==min(d_m_xa(:)));
    i_lv=find(d_l_xv==min(d_l_xv(:)));
    i_mv=find(d_m_xv==min(d_m_xv(:))); %obtaining the index of the closest matching location/meaning
    
    i_la_plt(y)=l(i_la); %these values represent the result of the delta function from the formulas
    i_ma_plt(y)=m(i_ma);
    i_lv_plt(y)=l(i_lv); %these values represent the result of the delta function from the formulas
    i_mv_plt(y)=m(i_mv); %storing the values for plotting

    a_x=zeros(n,t_max+1,n);
    v_x=zeros(n,t_max+1,n); %creating the auditory and visual amplitude arrays, to store in. Note these are still zeros (might better be Nan's or empty as 0 is also a legitimate location)
    
    f_a_x=c_av_s+sig_e_ax;
    f_v_x=k*c_av_s+sig_e_vx; %using the same linear functions as before, but now noise is added
    
    a_x(i_ma,:,i_la)=1*1*f_a_x;
    v_x(i_mv,:,i_lv)=1*1*f_v_x; %storing the auditory/visual amplituded in their respective meaing/location bins
end
figure(1);
scatter(i_l_plt,i_m_plt); 
hold on;
scatter(i_la_plt,i_ma_plt);
scatter(i_lv_plt,i_mv_plt); %scatter plotting the indices to compare the observations to the actual meaning/location
xlabel('location');
ylabel('semantic meaning');
axis([l_low l_high m_low m_high])
legend('actual','auditory percept','visual percept')
%% C=2

r=3; %how many times should the values of the meaning or location shift
i_lvs_plt=zeros(r,1);
i_mvs_plt=zeros(r,1);
i_las_plt=zeros(r,1);
i_mas_plt=zeros(r,1);
i_la_plt=zeros(r,1);
i_ma_plt=zeros(r,1);
i_lv_plt=zeros(r,1);
i_mv_plt=zeros(r,1); %creating arrays to store the meaning/location value for later plotting

for y=1:r
    k=0.3;
    t_max=10;
    n=10;
    
    l_low=-5;
    l_high=5;
    m_low=-2;
    m_high=2;

    sig_e_as=0.1;
    sig_e_vs=0.1;
    sig_l_as=0.6;
    sig_m_as=0.1;
    sig_l_vs=0.6;
    sig_m_vs=0.1;

    l=linspace(l_low,l_high,n);
    m=linspace(m_low,m_high,n);
    
    L_a_s=normrnd(0,sig_l_as);
    M_a_s=normrnd(0,sig_m_as);
    L_v_s=normrnd(0,sig_l_vs);
    M_v_s=normrnd(0,sig_m_vs); %generating the meaning and location in the observation, based on a gaussian centered on the actual meaning/location and the above specified noise levels

    [~,d_la]=dsearchn(L_a_s,l');
    [~,d_ma]=dsearchn(M_a_s,m');
    i_las=find(d_la==min(d_la(:)));
    i_mas=find(d_ma==min(d_ma(:)));
    [~,d_lv]=dsearchn(L_v_s,l');
    [~,d_mv]=dsearchn(M_v_s,m');
    i_lvs=find(d_lv==min(d_lv(:)));
    i_mvs=find(d_mv==min(d_mv(:))); %obtaining the index of the closest matching location/meaning
    
    i_las_plt(y)=l(i_las);
    i_mas_plt(y)=m(i_mas);
    i_lvs_plt(y)=l(i_lvs);
    i_mvs_plt(y)=m(i_mvs); %storing the meaing/location for plotting

    a_s=zeros(n,t_max+1,n);
    v_s=zeros(n,t_max+1,n); %creating the auditory and visual observation lists, now still empty/zero
    
    c_a_s=[normrnd(0,sig_e_as)];
    c_v_s=[normrnd(0,sig_e_vs)]; %drawing c_av^s from the normal distribution, for both the visual as auditory parts
    for i=1:t_max
        e_as=normrnd(0,sig_e_as);
        c_a_s=[c_a_s c_a_s(i)+e_as];
        e_vs=normrnd(0,sig_e_vs);
        c_v_s=[c_v_s c_v_s(i)+e_vs]; %storing and adding noise
    end
    f_a_s=c_a_s;
    f_v_s=k*c_v_s; %applying the linear mapping functions
    
    a_s(i_mas,:,i_las)=1*1*f_a_s;
    v_s(i_mvs,:,i_lvs)=1*1*f_v_s; %generating the result, and storing it in the auditory/visual observation lists
    
    %now the observation, denoted by x:
    sig_e_ax=0.1;
    sig_e_vx=0.1;
    sig_la=0.5;
    sig_ma=0.5;
    sig_lv=0.5;
    sig_mv=0.5; %setting some parameters
    
    L_av_xa=normrnd(L_a_s,sig_la);
    L_av_xv=normrnd(L_v_s,sig_lv);
    M_av_xa=normrnd(M_a_s,sig_ma);
    M_av_xv=normrnd(M_v_s,sig_mv); %generating the meaning and location in the observation, based on a gaussian centered on the actual meaning/location and the above specified noise levels

    [~,d_l_xa]=dsearchn(L_av_xa,l');
    [~,d_m_xa]=dsearchn(M_av_xa,m');
    [~,d_l_xv]=dsearchn(L_av_xv,l');
    [~,d_m_xv]=dsearchn(M_av_xv,m'); %getting the distances from the observations to the possible locations/meanings, l and m
    i_lax=find(d_l_xa==min(d_l_xa(:)));
    i_max=find(d_m_xa==min(d_m_xa(:)));
    i_lvx=find(d_l_xv==min(d_l_xv(:)));
    i_mvx=find(d_m_xv==min(d_m_xv(:))); %obtaining the index of the closest matching location/meaning
    
    i_la_plt(y)=l(i_lax); %these values represent the result of the delta function from the formulas
    i_ma_plt(y)=m(i_max);
    i_lv_plt(y)=l(i_lvx); %these values represent the result of the delta function from the formulas
    i_mv_plt(y)=m(i_mvx); %storing the values for plotting

    a_x=zeros(n,t_max+1,n);
    v_x=zeros(n,t_max+1,n); %creating the auditory and visual amplitude arrays, to store in. Note these are still zeros (might better be Nan's or empty as 0 is also a legitimate location)
    
    f_a_x=c_a_s+sig_e_ax;
    f_v_x=k*c_v_s+sig_e_vx; %using the same linear functions as before, but now noise is added
    
    a_x(i_max,:,i_lax)=1*1*f_a_x;
    v_x(i_mvx,:,i_lvx)=1*1*f_v_x; %storing the auditory/visual amplituded in their respective meaing/location bins
end
figure(1);
scatter(i_las_plt,i_mas_plt); 
hold on;
scatter(i_lvs_plt,i_mvs_plt); 
scatter(i_la_plt,i_ma_plt);
scatter(i_lv_plt,i_mv_plt); %scatter plotting the indices to compare the observations to the actual meaning/location
xlabel('location');
ylabel('semantic meaning');
axis([l_low l_high m_low m_high])
legend('actual','auditory percept','visual percept')