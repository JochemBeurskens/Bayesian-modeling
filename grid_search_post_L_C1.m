function[varargout]=grid_search_post_L_C1(k,t_max,n,l,varargin)    
    % in this program a grid search will be performed in order to get the
    % liklihood
    
    %parameters for C=1
    sig_l_s=varargin{1};
    sig_lax=varargin{2};
    sig_lvx=varargin{3};
    sig_e_s=varargin{4};
    sig_e_ax=varargin{5};
    sig_e_vx=varargin{6};
    p_c=varargin{7};
    p_h=varargin{8};
    % %parameters for C=2
    % sig_l_vs=1;
    % sig_l_as=4;
    % sig_lvx=1;
    % sig_lax=4;
    % sig_e_s=0.5;
    % sig_e_ax=0.5;
    % sig_e_vx=0.5;
    % p_c=0.3;
    % p_h=0.2;
    % sampling from the generative model for the case of C=1
    L=zeros(n,n);
    for i=1:1 %repeated draws to get different centers, in case we want ot estimate the average distribution for different centers
        [i_la_plt,i_lv_plt,m,a_x,v_x,L_av_s,c_av_s,m_c,m_h]=gen_model(k,t_max,n,l,sig_l_s,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h);
    
        % setting the space parameters for the grid search
        lrange=linspace(-l,l,n);
        mrange=linspace(1,4,4);
        trange=linspace(1,-10,10);
        counts_La = hist(i_la_plt,lrange);
        counts_Lv = hist(i_lv_plt,lrange);
        L=L+[counts_La'.*counts_Lv]; %the shape of this distribution should be similar to that of the likelihood
    end
    % the grid search of the posterior for the case of C=1
    % in this case the equation used represents P(C|L^{x})
    dist_L=zeros(length(lrange),length(lrange));
    f_LI=(sig_lax.^2.*sig_l_s.^2+sig_lvx.^2.*sig_l_s.^2+sig_lax.^2.*sig_lvx.^2);
    f_L=1/(2*pi*sqrt(f_LI));
    f_CI=(sig_e_ax.^2.*sig_e_s.^2+sig_e_vx.^2.*sig_e_s.^2+sig_e_ax.^2.*sig_e_vx.^2);
    f_C=1/(2*pi*sqrt(f_CI));
    for la_step=1:length(lrange)
        for lv_step=1:length(lrange)
            dist_L(la_step,lv_step)=exp(-(1/(2*f_LI))*((lrange(la_step)-lrange(lv_step)).^2*sig_l_s.^2+(lrange(lv_step)).^2*sig_lax.^2+(lrange(la_step)).^2*sig_lvx.^2));
        end
    end
    dist_L=f_L*dist_L; 
    
    % % now plotting the distributions
    % figure(1)
    % % surf(lrange,lrange,dist_L); %make a 3dimensional plot of the probability distribution
    % contour(lrange,lrange,dist_L);
    % % note that the same stuff happens for c, only now it changes over time.
    % % But the shape of the equations is very similar, and the time dynamics are
    % % taken into account in the equations.
    varargout{1}=dist_L;

end