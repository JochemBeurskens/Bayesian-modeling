function[varargout] = Normalised_Probabilities_multimodal_func(k,t_max,n,l,sig_l_s,sig_l_as,sig_l_vs,sig_lvx,sig_lax,sig_e_s,sig_e_ax,sig_e_vx,p_c,p_h,i_la_plt,i_lv_plt)
%% setting the parameters used in the probability distributions
L_xa=i_la_plt;
L_xv=i_lv_plt;
L_s1=-10.5; %the means of the priors
L_s2=-3.5;
L_s3=3.5;
L_s4=10.5;
p_p=1/4;
%% Now the C=1 pseudo-normalised probability for the L parameters
forefac1    = 2*pi*sqrt((sig_lax^2)*(sig_l_s^2)+(sig_lvx^2)*(sig_l_s^2)+(sig_lvx^2)*(sig_lax^2));
exponent11(1)   = exp(-0.5* ( ( ((L_xa(1)-L_xv(1))^2)*sig_l_s+ ((L_xa(1)-L_s1)^2)*sig_lvx+ ((L_xv(1)-L_s1)^2)*sig_lax)/forefac1 ) );
exponent12(1)   = exp(-0.5* ( ( ((L_xa(1)-L_xv(1))^2)*sig_l_s+ ((L_xa(1)-L_s2)^2)*sig_lvx+ ((L_xv(1)-L_s2)^2)*sig_lax)/forefac1 ) );
exponent13(1)   = exp(-0.5* ( ( ((L_xa(1)-L_xv(1))^2)*sig_l_s+ ((L_xa(1)-L_s3)^2)*sig_lvx+ ((L_xv(1)-L_s3)^2)*sig_lax)/forefac1 ) );
exponent14(1)   = exp(-0.5* ( ( ((L_xa(1)-L_xv(1))^2)*sig_l_s+ ((L_xa(1)-L_s4)^2)*sig_lvx+ ((L_xv(1)-L_s4)^2)*sig_lax)/forefac1 ) );
for i=2:length(L_xa)
    exponent11(i)   = exponent11(i-1) *( exp(-0.5* ( ( ((L_xa(i)-L_xv(i))^2)*sig_l_s^2+ ((L_xa(i)-L_s1)^2)*sig_lvx^2+ ((L_xv(i)-L_s1)^2)*sig_lax^2)/forefac1 ) ) )/forefac1;
    exponent12(i)   = exponent12(i-1) *( exp(-0.5* ( ( ((L_xa(i)-L_xv(i))^2)*sig_l_s^2+ ((L_xa(i)-L_s2)^2)*sig_lvx^2+ ((L_xv(i)-L_s2)^2)*sig_lax^2)/forefac1 ) ) )/forefac1;
    exponent13(i)   = exponent13(i-1) *( exp(-0.5* ( ( ((L_xa(i)-L_xv(i))^2)*sig_l_s^2+ ((L_xa(i)-L_s3)^2)*sig_lvx^2+ ((L_xv(i)-L_s3)^2)*sig_lax^2)/forefac1 ) ) )/forefac1;
    exponent14(i)   = exponent14(i-1) *( exp(-0.5* ( ( ((L_xa(i)-L_xv(i))^2)*sig_l_s^2+ ((L_xa(i)-L_s4)^2)*sig_lvx^2+ ((L_xv(i)-L_s4)^2)*sig_lax^2)/forefac1 ) ) )/forefac1;
end
P_LC1       = p_p* exponent11+p_p* exponent12+p_p* exponent13+p_p* exponent14;

%% Now the C=2 pseudo-normalised probability for the L parameters
forefac2    = 2*pi*sqrt( ((sig_lax^2)+(sig_l_s^2))*((sig_lvx^2)+(sig_l_s^2)) );
exponent21(1)   = exp(-0.5* ( ( ((L_xa(1)-L_s1)^2)/((sig_lax^2)*(sig_l_s^2)) )+(((L_xv(1)-L_s1)^2)/((sig_lvx^2)*(sig_l_s^2))  ) ) );
exponent22(1)   = exp(-0.5* ( ( ((L_xa(1)-L_s2)^2)/((sig_lax^2)*(sig_l_s^2)) )+(((L_xv(1)-L_s2)^2)/((sig_lvx^2)*(sig_l_s^2))  ) ) );
exponent23(1)   = exp(-0.5* ( ( ((L_xa(1)-L_s3)^2)/((sig_lax^2)*(sig_l_s^2)) )+(((L_xv(1)-L_s3)^2)/((sig_lvx^2)*(sig_l_s^2))  ) ) );
exponent24(1)   = exp(-0.5* ( ( ((L_xa(1)-L_s4)^2)/((sig_lax^2)*(sig_l_s^2)) )+(((L_xv(1)-L_s4)^2)/((sig_lvx^2)*(sig_l_s^2))  ) ) );
for i=2:length(L_xv)
    exponent21(i)   = exponent21(i-1) *( exp(-0.5* ( ( ((L_xa(i)-L_s1)^2)/((sig_lax^2)*(sig_l_s^2)) )+(((L_xv(i)-L_s1)^2)/((sig_lvx^2)*(sig_l_s^2))  ) ) ) )/forefac2;
    exponent21(i)   = exponent22(i-1) *( exp(-0.5* ( ( ((L_xa(i)-L_s2)^2)/((sig_lax^2)*(sig_l_s^2)) )+(((L_xv(i)-L_s2)^2)/((sig_lvx^2)*(sig_l_s^2))  ) ) ) )/forefac2;
    exponent21(i)   = exponent23(i-1) *( exp(-0.5* ( ( ((L_xa(i)-L_s3)^2)/((sig_lax^2)*(sig_l_s^2)) )+(((L_xv(i)-L_s3)^2)/((sig_lvx^2)*(sig_l_s^2))  ) ) ) )/forefac2;
    exponent21(i)   = exponent24(i-1) *( exp(-0.5* ( ( ((L_xa(i)-L_s4)^2)/((sig_lax^2)*(sig_l_s^2)) )+(((L_xv(i)-L_s4)^2)/((sig_lvx^2)*(sig_l_s^2))  ) ) ) )/forefac2;
end
P_LC2       = p_p*exponent21+p_p*exponent22+p_p*exponent23+p_p*exponent24;

%% normalisation, this is done by dividing by the sum of the two, as P(xa,xv)=P(xa,xv|C=1)P(C=1)+P(xa,xv|C=2)P(C=2)
norm_fac    = P_LC2+P_LC1;
P_LC1       = P_LC1./norm_fac;
P_LC2       = P_LC2./norm_fac;

%% output 
varargout{1}=P_LC1;
varargout{2}=P_LC2;
end