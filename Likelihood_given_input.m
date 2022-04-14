function [varargout] = Likelihood_given_input(k,t_max,n,c,l,m,L,M,varargin)
if nargin==17
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
    percept_L_denominator= (1/sig_la)+(1/sig_lv)+(1/sig_l_s);
    percept_M_denominator= (1/sig_ma)+(1/sig_mv)+(1/sig_m_s);

    i_l_plt=zeros(t_max,1);
    i_m_plt=zeros(t_max,1);
    i_la_plt=zeros(t_max,1);
    i_ma_plt=zeros(t_max,1);
    i_lv_plt=zeros(t_max,1);
    i_mv_plt=zeros(t_max,1); %creating arrays to store the meaning/location value for later plotting
    f_a_s_plt=zeros(t_max,1);
    f_v_s_plt=zeros(t_max,1);
    
    %generating the signal, denoted by s:   
    L_av_s=L;%normrnd(0,sig_l_s);
    M_av_s=M;%normrnd(0,sig_m_s); %drawing a location and meaning from a prior, this is now a Gaussian, might be better to make this a discrete stepped function?
    [~,d_l]=dsearchn(L_av_s,l');
    [~,d_m]=dsearchn(M_av_s,m'); 
    i_l=find(d_l==min(d_l(:))); 
    i_m=find(d_m==min(d_m(:))); %finding the location and meaning closest to the ones drawn from the prior, these are then assumed as the meant location/meaning
    
    c_av_s=[0 + normrnd(0,sig_e_s)]; %drawing c_av^s from the normal distribution
    
    for i=1:t_max
        i_l_plt(i)=l(i_l);
        i_m_plt(i)=m(i_m); %storing the meaing/location for plotting
        e_s=normrnd(0,sig_e_s);
        c_av_s=[c_av_s c_av_s(i)+e_s];
        [~,d_c]=dsearchn(c_av_s(i),c');
        c_binned=find(d_c==min(d_c(:)));
        f_a_s=c(c_binned);
        f_v_s=k*c(c_binned); %applying the linear mapping functions
        f_a_s_plt(i)=f_a_s;
        f_v_s_plt(i)=f_v_s;
    
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
        
        i_la_plt(i)=(i_la); %these values represent the result of the delta function from the formulas
        i_ma_plt(i)=(i_ma);
        i_lv_plt(i)=(i_lv); %these values represent the result of the delta function from the formulas
        i_mv_plt(i)=(i_mv); %storing the values for plotting
       
        f_a_x=c_av_s(i)+sig_e_ax;
        f_v_x=k*c_av_s(i)+sig_e_vx; %using the same linear functions as before, but now noise is added
        [~,d_f_ax]=dsearchn(f_a_x,c');
        [~,d_f_vx]=dsearchn(f_v_x,(k*c)');
        i_fax=find(d_f_ax==min(d_f_ax(:)));
        i_fvx=find(d_f_vx==min(d_f_vx(:)));

        percept_L_av=( (L_av_xa/sig_la) + (L_av_xv/sig_lv) + (L_av_s/sig_l_s) )/percept_L_denominator;
        percept_M_av=( (M_av_xa/sig_ma) + (M_av_xv/sig_mv) + (M_av_s/sig_m_s) )/percept_M_denominator;
        [~,d_percept_L_av]=dsearchn(percept_L_av,l');
        [~,d_percept_M_av]=dsearchn(percept_M_av,m');
        i_percept_L_av=find(d_percept_L_av==min(d_percept_L_av(:)));
        i_percept_M_av=find(d_percept_M_av==min(d_percept_M_av(:)));
    end
    
    varargout{1}=i_fax;
    varargout{2}=i_fvx;
    varargout{3}=i_l_plt;
    varargout{4}=i_m_plt;
    
    varargout{5}=i_la_plt;
    varargout{6}=i_ma_plt;
    varargout{7}=i_lv_plt;
    varargout{8}=i_mv_plt;
    varargout{9}=i_percept_L_av;
    varargout{10}=i_percept_M_av;

end
if nargin==22
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
    L_v=varargin{13};
    M_v=varargin{14};
    %C=2
    percept_La_denominator= (1/sig_la)+(1/sig_l_as);
    percept_Ma_denominator= (1/sig_ma)+(1/sig_m_as);
    percept_Lv_denominator= (1/sig_lv)+(1/sig_l_vs);
    percept_Mv_denominator= (1/sig_mv)+(1/sig_m_vs);

    i_lvs_plt=zeros(t_max,1);
    i_mvs_plt=zeros(t_max,1);
    i_las_plt=zeros(t_max,1);
    i_mas_plt=zeros(t_max,1);
    i_la_plt=zeros(t_max,1);
    i_ma_plt=zeros(t_max,1);
    i_lv_plt=zeros(t_max,1);
    i_mv_plt=zeros(t_max,1); %creating arrays to store the meaning/location value for later plotting
    
    L_a_s=L;%normrnd(0,sig_l_as);
    M_a_s=M;%normrnd(0,sig_m_as);
    L_v_s=L_v;%normrnd(0,sig_l_vs);
    M_v_s=M_v;%normrnd(0,sig_m_vs); %generating the meaning and location in the observation, based on a gaussian centered on the actual meaning/location and the above specified noise levels
    
    [~,d_la]=dsearchn(L_a_s,l');
    [~,d_ma]=dsearchn(M_a_s,m');
    i_las=find(d_la==min(d_la(:)));
    i_mas=find(d_ma==min(d_ma(:)));
    [~,d_lv]=dsearchn(L_v_s,l');
    [~,d_mv]=dsearchn(M_v_s,m');
    i_lvs=find(d_lv==min(d_lv(:)));
    i_mvs=find(d_mv==min(d_mv(:))); %obtaining the index of the closest matching location/meaning
     
    c_a_s=[0+normrnd(0,sig_e_as)];
    c_v_s=[0+normrnd(0,sig_e_vs)]; %drawing c_av^s from the normal distribution, for both the visual as auditory parts
    
    for i=1:t_max
        i_las_plt(i)=l(i_las);
        i_mas_plt(i)=m(i_mas);
        i_lvs_plt(i)=l(i_lvs);
        i_mvs_plt(i)=m(i_mvs); %storing the meaning/location for plotting

        e_as=normrnd(0,sig_e_as);
        c_a_s=[c_a_s c_a_s(i)+e_as];
        [~,d_c_a]=dsearchn(c_a_s(i),c');
        c_binned_a=find(d_c_a==min(d_c_a(:)));
        e_vs=normrnd(0,sig_e_vs);
        c_v_s=[c_v_s c_v_s(i)+e_vs]; %storing and adding noise
        [~,d_c_v]=dsearchn(c_v_s(i),c');
        c_binned_v=find(d_c_v==min(d_c_v(:)));

        f_a_s=c(c_binned_a);
        f_v_s=k*c(c_binned_v); %applying the linear mapping functions
        f_a_s_plt(i)=f_a_s;
        f_v_s_plt(i)=f_v_s;    

        L_a_xa=normrnd(L_a_s,sig_la);
        L_v_xv=normrnd(L_v_s,sig_lv);
        M_a_xa=normrnd(M_a_s,sig_ma);
        M_v_xv=normrnd(M_v_s,sig_mv); %generating the meaning and location in the observation, based on a gaussian centered on the actual meaning/location and the above specified noise levels
    
        [~,d_l_xa]=dsearchn(L_a_xa,l');
        [~,d_m_xa]=dsearchn(M_a_xa,m');
        [~,d_l_xv]=dsearchn(L_v_xv,l');
        [~,d_m_xv]=dsearchn(M_v_xv,m'); %getting the distances from the observations to the possible locations/meanings, l and m
        i_lax=find(d_l_xa==min(d_l_xa(:)));
        i_max=find(d_m_xa==min(d_m_xa(:)));
        i_lvx=find(d_l_xv==min(d_l_xv(:)));
        i_mvx=find(d_m_xv==min(d_m_xv(:))); %obtaining the index of the closest matching location/meaning
        
        i_la_plt(i)=(i_lax); %these values represent the result of the delta function from the formulas
        i_ma_plt(i)=(i_max);
        i_lv_plt(i)=(i_lvx); %these values represent the result of the delta function from the formulas
        i_mv_plt(i)=(i_mvx);

        f_a_x=c(c_binned_a)+sig_e_ax;
        f_v_x=k*c(c_binned_v)+sig_e_vx; %using the same linear functions as before, but now noise is added

        [~,d_f_ax]=dsearchn(f_a_x,c');
        [~,d_f_vx]=dsearchn(f_v_x,(k*c)');
        i_fax=find(d_f_ax==min(d_f_ax(:)));
        i_fvx=find(d_f_vx==min(d_f_vx(:)));

        percept_L_a=( (L_a_xa/sig_la) + (L_a_s/sig_l_as) )/percept_La_denominator;
        percept_M_a=( (M_a_xa/sig_ma) + (M_a_s/sig_m_as) )/percept_Ma_denominator;
        [~,d_percept_L_a]=dsearchn(percept_L_a,l');
        [~,d_percept_M_a]=dsearchn(percept_M_a,m');
        i_percept_L_a=find(d_percept_L_a==min(d_percept_L_a(:)));
        i_percept_M_a=find(d_percept_M_a==min(d_percept_M_a(:)));

        percept_L_v=( (L_v_xv/sig_lv) + (L_v_s/sig_l_vs) )/percept_Lv_denominator;
        percept_M_v=( (M_v_xv/sig_mv) + (M_v_s/sig_m_vs) )/percept_Mv_denominator;
        [~,d_percept_L_v]=dsearchn(percept_L_v,l');
        [~,d_percept_M_v]=dsearchn(percept_M_v,m');
        i_percept_L_v=find(d_percept_L_v==min(d_percept_L_v(:)));
        i_percept_M_v=find(d_percept_M_v==min(d_percept_M_v(:)));
    end
    
    %now the observation, denoted by x:

    varargout{1}=i_fax;
    varargout{2}=i_fvx;

    varargout{3}=i_las_plt;
    varargout{4}=i_mas_plt;
    varargout{5}=i_lvs_plt;
    varargout{6}=i_mvs_plt;
    
    varargout{7}=i_la_plt;
    varargout{8}=i_ma_plt;
    varargout{9}=i_lv_plt;
    varargout{10}=i_mv_plt;
    varargout{11}=i_percept_L_a;
    varargout{12}=i_percept_M_a;
    varargout{13}=i_percept_L_v;
    varargout{14}=i_percept_M_v;
end
end