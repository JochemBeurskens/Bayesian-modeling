function [a_s,v_s,varargout] = func_genmod(k,t_max,n,l,m,varargin)
if nargin==14
    sig_e_s=varargin{1};
    sig_l_s=varargin{2};
    sig_m_s=varargin{3};
    sig_e_ax=varargin{4};
    sig_e_vx=varargin{5};
    sig_la=varargin{6};
    sig_ma=varargin{7};
    sig_lv=varargin{8};
    sig_mv=varargin{9};
    % C=1
    %the code below lets a person estimate a cue for three consecutive
    %timesteps. For all of these the location of the target remains the same.
    
    i_l_plt=zeros(t_max,1);
    i_m_plt=zeros(t_max,1);
    i_la_plt=zeros(t_max,1);
    i_ma_plt=zeros(t_max,1);
    i_lv_plt=zeros(t_max,1);
    i_mv_plt=zeros(t_max,1); %creating arrays to store the meaning/location value for later plotting

    %generating the signal, denoted by s:   
    L_av_s=normrnd(0,sig_l_s);
    M_av_s=normrnd(0,sig_m_s); %drawing a location and meaning from a prior, this is now a Gaussian, might be better to make this a discrete stepped function?
    [~,d_l]=dsearchn(L_av_s,l');
    [~,d_m]=dsearchn(M_av_s,m'); 
    i_l=find(d_l==min(d_l(:))); 
    i_m=find(d_m==min(d_m(:))); %finding the location and meaning closest to the ones drawn from the prior, these are then assumed as the meant location/meaning

    a_s=zeros(n,t_max,n);
    v_s=zeros(n,t_max,n); %creating the auditory and visual observation lists, now still empty/zero
    
    a_x=zeros(n,t_max,n);
    v_x=zeros(n,t_max,n); %creating the auditory and visual amplitude arrays, to store in. Note these are still zeros (might better be Nan's or empty as 0 is also a legitimate location)
    
    c_av_s=[normrnd(0,sig_e_s)]; %drawing c_av^s from the normal distribution

    for i=1:t_max
        i_l_plt(i)=l(i_l);
        i_m_plt(i)=m(i_m); %storing the meaing/location for plotting
        e_s=normrnd(0,sig_e_s);
        c_av_s=[c_av_s c_av_s(i)+e_s]; %storing and adding noise
        f_a_s=c_av_s(i);
        f_v_s=k*c_av_s(i); %applying the linear mapping functions
        
        a_s(i_m,i,i_l)=1*1*f_a_s;
        v_s(i_m,i,i_l)=1*1*f_v_s; %generating the result, and storing it in the auditory/visual observation lists
        
        %making the observations, denoted by x:
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
        
        i_la_plt(i)=l(i_la); %these values represent the result of the delta function from the formulas
        i_ma_plt(i)=m(i_ma);
        i_lv_plt(i)=l(i_lv); %these values represent the result of the delta function from the formulas
        i_mv_plt(i)=m(i_mv); %storing the values for plotting
       
        f_a_x=c_av_s(i)+sig_e_ax;
        f_v_x=k*c_av_s(i)+sig_e_vx; %using the same linear functions as before, but now noise is added
        
        a_x(i_ma,i,i_la)=1*1*f_a_x;
        v_x(i_mv,i,i_lv)=1*1*f_v_x; %storing the auditory/visual amplituded in their respective meaing/location bins
    end
    varargout{1}=c_av_s;
    varargout{2}=i_l_plt;
    varargout{3}=i_m_plt;

    varargout{4}=i_la_plt;
    varargout{5}=i_ma_plt;
    varargout{6}=i_lv_plt;
    varargout{7}=i_mv_plt;

    varargout{8}=a_x;
    varargout{9}=v_x;
end
if nargin==17
    sig_e_as=varargin{1};
    sig_e_vs=varargin{2};
    sig_l_as=varargin{3};
    sig_l_vs=varargin{4};
    sig_m_as=varargin{5};
    sig_m_vs=varargin{6};
    sig_e_ax=varargin{7};
    sig_e_vx=varargin{8};
    sig_la=varargin{9};
    sig_ma=varargin{10};
    sig_lv=varargin{11};
    sig_mv=varargin{12};
    %C=2
    i_lvs_plt=zeros(t_max,1);
    i_mvs_plt=zeros(t_max,1);
    i_las_plt=zeros(t_max,1);
    i_mas_plt=zeros(t_max,1);
    i_la_plt=zeros(t_max,1);
    i_ma_plt=zeros(t_max,1);
    i_lv_plt=zeros(t_max,1);
    i_mv_plt=zeros(t_max,1); %creating arrays to store the meaning/location value for later plotting
    
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
    
    

    a_s=zeros(n,t_max,n);
    v_s=zeros(n,t_max,n); %creating the auditory and visual observation lists, now still empty/zero

    a_x=zeros(n,t_max,n);
    v_x=zeros(n,t_max,n); %creating the auditory and visual amplitude arrays, to store in. Note these are still zeros (might better be Nan's or empty as 0 is also a legitimate location)
        
    c_a_s=[normrnd(0,sig_e_as)];
    c_v_s=[normrnd(0,sig_e_vs)]; %drawing c_av^s from the normal distribution, for both the visual as auditory parts
    
    for i=1:t_max
        i_las_plt(i)=l(i_las);
        i_mas_plt(i)=m(i_mas);
        i_lvs_plt(i)=l(i_lvs);
        i_mvs_plt(i)=m(i_mvs); %storing the meaning/location for plotting

        e_as=normrnd(0,sig_e_as);
        c_a_s=[c_a_s c_a_s(i)+e_as];
        e_vs=normrnd(0,sig_e_vs);
        c_v_s=[c_v_s c_v_s(i)+e_vs]; %storing and adding noise

        f_a_s=c_a_s(i);
        f_v_s=k*c_v_s(i); %applying the linear mapping functions
        
        a_s(i_mas,i,i_las)=1*1*f_a_s;
        v_s(i_mvs,i,i_lvs)=1*1*f_v_s; %generating the result, and storing it in the auditory/visual observation lists
    
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
        
        i_la_plt(i)=l(i_lax); %these values represent the result of the delta function from the formulas
        i_ma_plt(i)=m(i_max);
        i_lv_plt(i)=l(i_lvx); %these values represent the result of the delta function from the formulas
        i_mv_plt(i)=m(i_mvx);

        f_a_x=c_a_s(i)+sig_e_ax;
        f_v_x=k*c_v_s(i)+sig_e_vx; %using the same linear functions as before, but now noise is added
        
        a_x(i_max,i,i_lax)=1*1*f_a_x;
        v_x(i_mvx,i,i_lvx)=1*1*f_v_x; %storing the auditory/visual amplituded in their respective meaing/location bins
    end
    
    %now the observation, denoted by x:

    varargout{1}=c_a_s;
    varargout{2}=c_v_s;


    varargout{3}=i_las_plt;
    varargout{4}=i_mas_plt;
    varargout{5}=i_lvs_plt;
    varargout{6}=i_mvs_plt;
    
    varargout{7}=i_la_plt;
    varargout{8}=i_ma_plt;
    varargout{9}=i_lv_plt;
    varargout{10}=i_mv_plt;

    varargout{11}=a_x;
    varargout{12}=v_x;
end
end