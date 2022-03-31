%for the likelihood calculation I will be using almost the same code, but
%as the very large arrays of a and v take too much space, I will not use
%them here. I will only calculate the parameters that decide where in these
%there are non-zero values, and then use this to make a frequentist
%estimation of the likelihood of a given draw

%% setting parameters

rng(0); %setting the random seed, so that the draws will have the same result each time, for now this is useful.

k=1.0; %setting the multiplication for the linear functions
t_max=1; %setting the amount of time steps for the source activity
n=5; %setting the amount of possible locations

lrange=[-5 5];
mrange=[-5 5];
crange=[-0.5 0.5];

sig_e_s=0.1;
sig_l_s=0.6;
sig_m_s=0.1;%setting noise parameters for the generating locations step. used for both C=1 and 2

l=linspace(lrange(1),lrange(2),n); %setting the range of locations and meanings
m=linspace(mrange(1),mrange(2),n);
c=linspace(crange(1),crange(2),n);

sig_e_ax=0.1;
sig_e_vx=0.1;
sig_la=0.5;
sig_ma=0.5;
sig_lv=0.5;
sig_mv=0.5; %setting the noise parameters for the subject's observation

%% drawing the a_x and v_x for which the likelihoods are to be estimated
[f_a_x,f_v_x,i_l_x,i_m_x,i_la_x,i_ma_x,i_lv_x,i_mv_x]=Likelihood(k,t_max,n,c,l,m,sig_e_s,sig_l_s,sig_m_s,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);
%% Now  going over the entire space state to find estimates of the probabilities for the parameters in the generated a_x and v_x
%note that it is here useful to also count the instances that are the same
%as the values in the a_x and v_x
frac_la=zeros(n,n,n);
frac_ma=zeros(n,n,n);
frac_lv=zeros(n,n,n);
frac_mv=zeros(n,n,n);
%looping over the l and then m space, as these depend on the location, do
%not need to loop over c (or f), as this is location independent so can be
%obtained once and then get the probabilities for specific values.
li=1;
for lr=l
    mi=1;
    for mr=m
        repeat=10000;
        L=lr;
        M=mr;
        
        f_v_s_store=zeros(repeat,t_max);
        f_a_s_store=zeros(repeat,t_max);
        
        i_la_store=zeros(repeat,t_max);
        i_ma_store=zeros(repeat,t_max);
        i_lv_store=zeros(repeat,t_max);
        i_mv_store=zeros(repeat,t_max);
        for q=1:repeat
            %C=1, run model and plot
            [f_a_s_plt_stat,f_v_s_plt_stat,i_l_plt_stat,i_m_plt_stat,i_la_plt_stat,i_ma_plt_stat,i_lv_plt_stat,i_mv_plt_stat]=Likelihood_given_input_c1(k,t_max,n,c,l,m,L,M,sig_e_s,sig_l_s,sig_m_s,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);
        %     [f_a_s_plt,f_v_s_plt,i_l_plt,i_m_plt,i_la_plt,i_ma_plt,i_lv_plt,i_mv_plt]=Likelihood(k,t_max,n,c,l,m,sig_e_s,sig_l_s,sig_m_s,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);
            i_la_store(q,:)=i_la_plt_stat;
            i_ma_store(q,:)=i_ma_plt_stat;
            i_lv_store(q,:)=i_lv_plt_stat;
            i_mv_store(q,:)=i_mv_plt_stat;
            f_v_s_store(q,:)=f_v_s_plt_stat;
            f_a_s_store(q,:)=f_a_s_plt_stat;
        end
        [C_la,R_la] = groupcounts(i_la_store(:,1)); 
        [~,d_la]=dsearchn(R_la,l');
        i_la=find(d_la==min(d_la(:))); 
        frac_la(li,mi,i_la)=C_la;
        
        [C_ma,R_ma] = groupcounts(i_ma_store(:,1)); 
        [~,d_ma]=dsearchn(R_ma,m');
        i_ma=find(d_ma==min(d_ma(:))); 
        frac_ma(li,mi,i_ma)=C_ma;

        [C_lv,R_lv] = groupcounts(i_lv_store(:,1)); 
        [~,d_lv]=dsearchn(R_lv,l');
        i_lv=find(d_lv==min(d_lv(:))); 
        frac_lv(li,mi,i_lv)=C_lv;

        [C_mv,R_mv] = groupcounts(i_mv_store(:,1)); 
        [~,d_mv]=dsearchn(R_mv,m');
        i_mv=find(d_mv==min(d_mv(:))); 
        frac_mv(li,mi,i_mv)=C_mv;
        %so the above grid means: when we put the actual stimulus at
        %locations given by the l(li), m(mi) then we get a distribution for
        %the auditory/visual stimulus (after taking noise into account)
        %which has a similar distribution as the (li,mi,:) parameter
        %belonging to that stimulus (so for la it would be the location
        %distribution for a specific match of l(li) and m(mi))

        mi=mi+1;
    end
    li=li+1;
end

%% now that we have obtained what fractions of the totcal amount will fall where, when the L^s and M^s are moved about, we can count the probability for a specific a^x and v^x

%make the a_x and v_x parameters that will be used in the estimation
[f_a_s_plt,f_v_s_plt,i_l_plt,i_m_plt,i_la_plt,i_ma_plt,i_lv_plt,i_mv_plt]=Likelihood(k,t_max,n,c,l,m,sig_e_s,sig_l_s,sig_m_s,sig_e_ax,sig_e_vx,sig_la,sig_ma,sig_lv,sig_mv);
