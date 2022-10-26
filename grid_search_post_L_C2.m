% function[varargout]=grid_search_post_L_C2(k,t_max,n,l,varargin)    
    % in this program a grid search will be performed in order to get the
    % liklihood
    
    %parameters for C=2
    sig_l_vs=0.5;
    sig_l_as=0.4;
    sig_lvx=0.5;
    sig_lax=0.4;
    sig_e_s=0.5;
    sig_e_ax=0.5;
    sig_e_vx=0.5;
    p_c=0.3;
    p_h=0.2;
    %sampling from the generative model for the case of C=1
    L=zeros(n,n);
    for i=1:1 %repeated draws to get different centers, in case we want ot estimate the average distribution for different centers
        [i_la_plt,i_lv_plt,m,a_x,v_x,L_av_vs,L_av_as,c_av_s,m_c,m_h]=gen_model(k,t_max,n,l,sig_l_vs,sig_l_as,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h);
    
        % setting the space parameters for the grid search
        lrange=l;%linspace(-l,l,n);
        mrange=linspace(1,4,4);
        trange=linspace(1,-10,10);
        counts_La = hist(i_la_plt,lrange);
        counts_Lv = hist(i_lv_plt,lrange);
        L=L+[counts_La'.*counts_Lv]; %the shape of this distribution should be similar to that of the likelihood
    end
    % the grid search of the posterior for the case of C=1
    % in this case the equation used represents P(C|L^{x})
    dist_La=zeros(length(lrange),1);    
    dist_Lv=zeros(length(lrange),1);

    f_LIa=(sig_lax.^2+sig_l_as.^2);
    f_La=1/(2*pi*sqrt(f_LIa));
    f_LIv=(sig_lvx.^2+sig_l_vs.^2);
    f_Lv=1/(2*pi*sqrt(f_LIv));

    for la_step=1:length(lrange)
            dist_La(la_step)=exp(-(1/(2*f_LIa)) *(lrange(la_step)).^2);
    end
    for lv_step=1:length(lrange)
            dist_Lv(lv_step)=exp(-(1/(2*f_LIv)) *(lrange(lv_step)).^2);
    end
    dist_La=f_La*dist_La;
    dist_Lv=f_Lv*dist_Lv;
    
    % now plotting the distributions
    figure(1)
%     surf(lrange,lrange,L); %make a 3dimensional plot of the probability distribution
    plot(dist_La)
    hold on
    plot(dist_Lv)

%     varargout{1}=dist_La;
%     varargout{2}=dist_Lv;
% end