function[varargout]=gen_model_full(k,t_max,l,m,varargin)
    % in this code samples from the generative model are obtained, each of the
    % different types of data (loc,meaning,content) is sampled separately for
    % clarity.
    % rng(5); % set random seed
if nargin==20
    sig_l_s=varargin{1};
    sig_lvx=varargin{2};
    sig_lax=varargin{3};
    sig_m_s=varargin{4};
    sig_m_as=varargin{5};
    sig_m_vs=varargin{6};
    sig_m_ax=varargin{7};
    sig_m_vx=varargin{8};
    sig_e_s=varargin{9};
    sig_e_ax=varargin{10};
    sig_e_vx=varargin{11};
    p_c=varargin{12};
    p_h=varargin{13};
    n_r=varargin{14};
    n_lm=varargin{15};
    mean_rw=varargin{16};
    % C=1
    %the code below lets a person estimate a cue for three consecutive
    %timesteps. For all of these the location of the target remains the same.   
    % location of the stimuli for the case of C=1
    i_la_plt=zeros(t_max,n_lm);
    i_lv_plt=zeros(t_max,n_lm);
    l_i=[0 0 0 0]; %the location about which the original signal is centered
    %generating the signal, denoted by s:   
    L_av_s=normrnd(l_i(1),sig_l_s); %drawing a location and meaning from a prior, this is now a Gaussian
    [~,d_l]=dsearchn(L_av_s,l'); %getting the distances from the observations to the possible locations
    i_l=find(d_l==min(d_l(:))); %obtaining the index of the closest matching location
    
    %the M_s variables
    m_s=normrnd(0.5,sig_m_s);
    m_a=zeros(t_max,n_lm);
    m_v=zeros(t_max,n_lm);
    %the c_s variables
    for i=1:t_max
        e_s=normrnd(mean_rw,sig_e_s);
        if i==1
            c_av_s=[normrnd(mean_rw,sig_e_s)]; %start point of random walk
        elseif i>1
            c_av_s=[c_av_s c_av_s(i-1)+e_s]; %storing and adding noise
        end
    end
    a_s=zeros(length(l),length(m),t_max,n_lm,n_r);
    v_s=zeros(length(l),length(m),t_max,n_lm,n_r);
    a_x=zeros(length(l),length(m),t_max,n_lm,n_r);
    v_x=zeros(length(l),length(m),t_max,n_lm,n_r);
    
    for o=1:n_lm
        for i=1:t_max
            %making the observations, denoted by x:
            L_av_xa=normrnd(L_av_s,sig_lax);
            L_av_xv=normrnd(L_av_s,sig_lvx);
           
            [~,d_l_xa]=dsearchn(L_av_xa,l');
            [~,d_l_xv]=dsearchn(L_av_xv,l'); %getting the distances from the observations to the possible locations
            
            i_la=find(d_l_xa==min(d_l_xa(:)));
            i_lv=find(d_l_xv==min(d_l_xv(:))); %obtaining the index of the closest matching location
            
            i_la_plt(i,o)=(i_la); %these values represent the result of the delta function from the formulas
            i_lv_plt(i,o)=(i_lv); %these values represent the result of the delta function from the formulas
        end
    %         i_l=[l(i_l) l(i_l) l(i_l) l(i_l)];
        
        % meaning of the stimuli for the case of C=1
    
        for i=1:t_max
            %making the observations, denoted by x:
            m_av_xa=normrnd(m_s,sig_m_ax);
            m_av_xv=normrnd(m_s,sig_m_vx);
           
            [~,d_m_xa]=dsearchn(m_av_xa,m');
            [~,d_m_xv]=dsearchn(m_av_xv,m'); %getting the distances from the observations to the possible locations
            
            i_ma=find(d_m_xa==min(d_m_xa(:)));
            i_mv=find(d_m_xv==min(d_m_xv(:))); %obtaining the index of the closest matching location
            
            m_a(i,o)=(i_ma); %these values represent the result of the delta function from the formulas
            m_v(i,o)=(i_mv); %these values represent the result of the delta function from the formulas
            %making the observations, denoted by x:
%             m_av_xa=normrnd(m_s,sig_m_ax);
%             m_av_xv=normrnd(m_s,sig_m_vx);
%             if m_av_xa >= p_c
%                 m_a(i,o)=1;
%             elseif m_av_xa < p_c
%                 m_a(i,o)=2;
%             end
%             if m_av_xv >= p_c
%                 m_v(i,o)=1;
%             elseif m_av_xv < p_c
%                 m_v(i,o)=2;
%             end
        end
    
        % the random walk for the content for the case of C=1
        for p=1:n_r   
        %     c_av_s=[normrnd(5,sig_e_s)]; %drawing c_av^s from the normal distribution, startpoint of random walk
            for i=1:t_max
                f_a_s=c_av_s(i);
                f_v_s=k*c_av_s(i); %applying the linear mapping functions
        %     end
                a_s(i_la_plt(i,o),m_a(i,o),i)=1*1*f_a_s;  
                v_s(i_lv_plt(i,o),m_v(i,o),i)=1*1*f_v_s; %generating the result, and storing it in the auditory/visual observation lists
                f_a_x=f_a_s+normrnd(0,sig_e_ax);
                f_v_x=f_v_s+normrnd(0,sig_e_vx); %using the same linear functions as before, but now noise is added
                a_x(:,:,i,o,p)=normrnd(0,sig_e_ax,length(l),length(m));
                v_x(:,:,i,o,p)=normrnd(0,sig_e_vx,length(l),length(m));
                a_x(i_la_plt(i,o),m_a(i,o),i,o,p)=1*1*f_a_x;
                v_x(i_lv_plt(i,o),m_v(i,o),i,o,p)=1*1*f_v_x; %storing the auditory/visual amplituded in their respective meaing/location bins
            end
        end
    end
    varargout{1}=i_la_plt;
    varargout{2}=i_lv_plt;
    varargout{3}=m_a;
    varargout{4}=m_v;
    varargout{5}=a_x;
    varargout{6}=v_x;
    varargout{7}=i_l;
    varargout{8}=c_av_s;
    varargout{9}=m_av_xa;
    varargout{10}=m_av_xv;   
    varargout{11}=m_s;
end
if nargin == 22
    sig_l_vs=varargin{1};
    sig_l_as=varargin{2};
    sig_lvx=varargin{3};
    sig_lax=varargin{4};
    sig_m_s=varargin{5};
    sig_m_as=varargin{6};
    sig_m_vs=varargin{7};
    sig_m_ax=varargin{8};
    sig_m_vx=varargin{9};
    sig_e_as=varargin{10};
    sig_e_vs=varargin{11};
    sig_e_ax=varargin{12};
    sig_e_vx=varargin{13};
    p_c=varargin{14};
    p_h=varargin{15};
    n_r=varargin{16};
    n_lm=varargin{17};
    mean_rw=varargin{18};
    % C=2
    %the code below lets a person estimate a cue for three consecutive
    %timesteps. For all of these the location of the target remains the same.   
    % location of the stimuli for the case of C=1
    i_la_plt=zeros(t_max,n_lm);
    i_lv_plt=zeros(t_max,n_lm);
    l_i=[0 0 0 0]; %the location about which the original signal is centered
    %generating the signal, denoted by s:   
    L_a_s=normrnd(l_i(1),sig_l_as); %drawing a location and meaning from a prior, this is now a Gaussian
    [~,d_l]=dsearchn(L_a_s,l'); %getting the distances from the observations to the possible locations
    i_la=find(d_l==min(d_l(:))); %obtaining the index of the closest matching location

    L_v_s=normrnd(l_i(1),sig_l_vs); %drawing a location and meaning from a prior, this is now a Gaussian
    [~,d_l]=dsearchn(L_v_s,l'); %getting the distances from the observations to the possible locations
    i_lv=find(d_l==min(d_l(:))); %obtaining the index of the closest matching location
    %the M_s variables
    m_sa=normrnd(0.5,sig_m_as);
    m_sv=normrnd(0.5,sig_m_vs);

    m_a=zeros(t_max,n_lm);
    m_v=zeros(t_max,n_lm);
    %the c_s variables
    for i=1:t_max
        e_as=normrnd(mean_rw,sig_e_as);
        if i==1
            c_a_s=[normrnd(mean_rw,sig_e_as)]; %start point of random walk
        elseif i>1
            c_a_s=[c_a_s c_a_s(i-1)+e_as]; %storing and adding noise
        end
        e_vs=normrnd(mean_rw,sig_e_vs);
        if i==1
            c_v_s=[normrnd(mean_rw,sig_e_vs)]; %start point of random walk
        elseif i>1
            c_v_s=[c_v_s c_v_s(i-1)+e_vs]; %storing and adding noise
        end
    end
    a_s=zeros(length(l),2,t_max,n_lm,n_r);
    v_s=zeros(length(l),2,t_max,n_lm,n_r);
    a_x=zeros(length(l),2,t_max,n_lm,n_r);
    v_x=zeros(length(l),2,t_max,n_lm,n_r);
    
    for o=1:n_lm
        for i=1:t_max
            %making the observations, denoted by x:
            L_a_x=normrnd(L_a_s,sig_lax);
            L_v_x=normrnd(L_v_s,sig_lvx);
           
            [~,d_l_xa]=dsearchn(L_a_x,l');
            [~,d_l_xv]=dsearchn(L_v_x,l'); %getting the distances from the observations to the possible locations
            
            i_la=find(d_l_xa==min(d_l_xa(:)));
            i_lv=find(d_l_xv==min(d_l_xv(:))); %obtaining the index of the closest matching location
            
            i_la_plt(i,o)=(i_la); %these values represent the result of the delta function from the formulas
            i_lv_plt(i,o)=(i_lv); %these values represent the result of the delta function from the formulas
        end
    %         i_l=[l(i_l) l(i_l) l(i_l) l(i_l)];
        
        % meaning of the stimuli for the case of C=1
    
        for i=1:t_max
            m_a_x=normrnd(m_sa,sig_m_ax);
            m_v_x=normrnd(m_sv,sig_m_vx);
           
            [~,d_m_xa]=dsearchn(m_a_x,m');
            [~,d_m_xv]=dsearchn(m_v_x,m'); %getting the distances from the observations to the possible locations
            
            i_ma=find(d_m_xa==min(d_m_xa(:)));
            i_mv=find(d_m_xv==min(d_m_xv(:))); %obtaining the index of the closest matching location
            
            m_a(i,o)=(i_ma); %these values represent the result of the delta function from the formulas
            m_v(i,o)=(i_mv); %these values represent the result of the delta function from the formulas

%             %making the observations, denoted by x:
%             m_a_x=normrnd(m_s,sig_m_ax);
%             m_v_x=normrnd(m_s,sig_m_vx);
%             if m_a_x >= p_c
%                 m_a(i,o)=1;
%             elseif m_a_x < p_c
%                 m_a(i,o)=2;
%             end
%             if m_v_x >= p_c
%                 m_v(i,o)=1;
%             elseif m_v_x < p_c
%                 m_v(i,o)=2;
%             end
        end

    
        % the random walk for the content for the case of C=1
        for p=1:n_r   
        %     c_av_s=[normrnd(5,sig_e_s)]; %drawing c_av^s from the normal distribution, startpoint of random walk
            for i=1:t_max
                f_a_s=c_a_s(i);
                f_v_s=k*c_v_s(i); %applying the linear mapping functions
        %     end
                a_s(i_la_plt(i,o),m_a(i,o),i)=1*1*f_a_s;  
                v_s(i_lv_plt(i,o),m_v(i,o),i)=1*1*f_v_s; %generating the result, and storing it in the auditory/visual observation lists
                f_a_x=f_a_s+normrnd(0,sig_e_ax);
                f_v_x=f_v_s+normrnd(0,sig_e_vx); %using the same linear functions as before, but now noise is added
                a_x(:,:,i,o,p)=normrnd(0,sig_e_ax,length(l),2);
                v_x(:,:,i,o,p)=normrnd(0,sig_e_vx,length(l),2);
                a_x(i_la_plt(i,o),m_a(i,o),i,o,p)=1*1*f_a_x;
                v_x(i_lv_plt(i,o),m_v(i,o),i,o,p)=1*1*f_v_x; %storing the auditory/visual amplituded in their respective meaing/location bins
            end
        end
    end

    varargout{1}=i_la_plt;
    varargout{2}=i_lv_plt;
    varargout{3}=m_a;
    varargout{4}=m_v;
    varargout{5}=a_x;
    varargout{6}=v_x;
    varargout{7}=i_la;
    varargout{8}=i_lv;
    varargout{9}=c_a_s;    
    varargout{10}=c_v_s;
    varargout{11}=m_a_x;
    varargout{12}=m_v_x;   
    varargout{13}=m_sa;
    varargout{14}=m_sv;
end
end