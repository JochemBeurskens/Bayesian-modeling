function[varargout]=gen_model_full(k,t_max,n,l,varargin)
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
    a_s=zeros(length(l),2,t_max,n_lm,n_r);
    v_s=zeros(length(l),2,t_max,n_lm,n_r);
    a_x=zeros(length(l),2,t_max,n_lm,n_r);
    v_x=zeros(length(l),2,t_max,n_lm,n_r);
    
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
            if m_av_xa >= p_c
                m_a(i,o)=1;
            elseif m_av_xa < p_c
                m_a(i,o)=2;
            end
            if m_av_xv >= p_c
                m_v(i,o)=1;
            elseif m_av_xv < p_c
                m_v(i,o)=2;
            end
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
    varargout{7}=i_l;
    varargout{8}=c_av_s;
    varargout{9}=m_av_xa;
    varargout{10}=m_av_xv;   
    varargout{11}=m_s;
end
if nargin == 6
    sig_l_vs=varargin{1};
    sig_l_as=varargin{2};
    sig_lvx=varargin{3};
    sig_lax=varargin{4};
    sig_e_s=varargin{5};
    sig_e_ax=varargin{6};
    sig_e_vx=varargin{7};
    p_c=varargin{8};
    p_h=varargin{9};
    % location of the stimuli for the case of C=2
    i_la_plt=zeros(t_max,1);
    i_lv_plt=zeros(t_max,1);
    
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
        
        i_la_plt(i)=i_la; %these values represent the result of the delta function from the formulas
        i_lv_plt(i)=i_lv; %these values represent the result of the delta function from the formulas
    end
    i_lvs=[l(i_lvs) l(i_lvs) l(i_lvs) l(i_lvs)];
    i_las=[l(i_las) l(i_las) l(i_las) l(i_las)];
    
    % meaning of the stimuli for the case of C=2
    m_a=zeros(t_max,1);
    m_v=zeros(t_max,1);
    m_i=[1 2 3 4];
    for i=1:t_max
        m_c=binornd(1,p_c); %getting t_max draws from the binomial distribution determining whether the stimuli should be congruent
        m_h=binornd(1,p_h); %getting t_max draws from the binomial distribution determining whether the stimuli should be human
        if m_c==0 && m_h==0
            m_a(i)=1;
            m_v(i)=1;
        elseif m_c==0 && m_h==1
            m_a(i)=2;
            m_v(i)=1;
        elseif m_c==1 && m_h==0
            m_a(i)=1;
            m_v(i)=2;
        elseif m_c==1 && m_h==1
            m_a(i)=2;
            m_v(i)=2;
        end
    end

    % the random walk for the content for the case of C=2
    a_s=zeros(length(l),length(m_i),t_max);
    v_s=zeros(length(l),length(m_i),t_max);
    a_x=zeros(length(l),length(m_i),t_max);
    v_x=zeros(length(l),length(m_i),t_max);
    
%     c_av_s=[normrnd(5,sig_e_s)]; %drawing c_av^s from the normal distribution
    for i=1:t_max
        e_s=normrnd(10,sig_e_s);
        if i==1
            c_av_s=[normrnd(10,sig_e_s)]; %start point of random walk
        elseif i>1
            c_av_s=[c_av_s c_av_s(i-1)+e_s]; %storing and adding noise
        end
        c_av_s=[c_av_s c_av_s(i)+e_s]; %storing and adding noise
        f_a_s=c_av_s(i);
%         if i >= t_max/2.5
        f_v_s=k*c_av_s(i); %applying the linear mapping functions
%         else
%             f_v_s=normrnd(0,0.05);
%         end
%     end
        a_s(i_la_plt(i),m_a(i),i)=1*1*f_a_s;  
        v_s(i_lv_plt(i),m_v(i),i)=1*1*f_v_s; %generating the result, and storing it in the auditory/visual observation lists
        f_a_x=f_a_s+normrnd(0,sig_e_ax);
        f_v_x=f_v_s+normrnd(0,sig_e_vx); %using the same linear functions as before, but now noise is added
        a_x(:,:,i)=normrnd(0,sig_e_ax,length(l),length(m_a));
        v_x(:,:,i)=normrnd(0,sig_e_vx,length(l),length(m_v));
        a_x(i_la_plt(i),m_a(i),i)=1*1*f_a_x;
        v_x(i_lv_plt(i),m_v(i),i)=1*1*f_v_x; %storing the auditory/visual amplituded in their respective meaing/location bins
    end

    varargout{1}=i_la_plt;
    varargout{2}=i_lv_plt;
    varargout{3}=m_a;
    varargout{4}=m_v;
    varargout{5}=a_x;
    varargout{6}=v_x;
    varargout{7}=L_av_vs;
    varargout{8}=L_av_as;
    varargout{9}=c_av_s;
    varargout{10}=m_c;
    varargout{11}=m_h;
end
end