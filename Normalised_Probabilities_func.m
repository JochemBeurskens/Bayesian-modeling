function[varargout] = Normalised_Probabilities_func(sig_l_s,sig_lvx,sig_lax,i_lv_plt,i_la_plt)
%% setting the parameters used in the probability distributions
L_xa=i_la_plt;
L_xv=i_lv_plt;
L_s=0; %the mean of the prior
%% Now the C=1 pseudo-normalised probability for the L parameters
forefac1    = 2*pi*sqrt(((sig_lax^2)*(sig_l_s^2))+((sig_lvx^2)*(sig_l_s^2))+((sig_lvx^2)*(sig_lax^2)));
forefac11    = (((sig_lax^2)*(sig_l_s^2))+((sig_lvx^2)*(sig_l_s^2))+((sig_lvx^2)*(sig_lax^2)));

% exponent1(1)   = exp(-0.5* ( ( ((L_xv(1)-L_xa(1))^2)*sig_l_s+ ((L_xa(1)-L_s)^2)*sig_lvx+ ((L_xv(1)-L_s)^2)*sig_lax)/forefac1 ) );
for i=1:length(L_xa)
    exponent1(i)   = ( exp(-0.5* ( ( ((L_xv(i)-L_xa(i))^2)*sig_l_s^2+ ((L_xa(i)-L_s)^2)*sig_lvx^2+ ((L_xv(i)-L_s)^2)*sig_lax^2))/forefac11  ) )/forefac1;%exponent1(i-1) *
end
P_LC1       = exponent1;

%% Now the C=2 pseudo-normalised probability for the L parameters
forefac2    = 2*pi*sqrt( ((sig_lax^2)+(sig_l_s^2))*((sig_lvx^2)+(sig_l_s^2)) );
% exponent2(1)   = exp(-0.5* ( ( ((L_xa(1)-L_s)^2)/((sig_lax^2)+(sig_l_s^2)) )+(((L_xv(1)-L_s)^2)/((sig_lvx^2)+(sig_l_s^2))  ) ) );
for i=1:length(L_xv)
    exponent2(i)   = ( exp(-0.5* ( ( ((L_xa(i)-L_s)^2)/((sig_lax^2)+(sig_l_s^2)) )+(((L_xv(i)-L_s)^2)/((sig_lvx^2)+(sig_l_s^2))  ) ) ) )/forefac2;%exponent2(i-1) *
end
P_LC2       = exponent2;

%% normalisation, this is done by dividing by the sum of the two, as P(xa,xv)=P(xa,xv|C=1)P(C=1)+P(xa,xv|C=2)P(C=2)
norm_fac    = P_LC2+P_LC1;
P_LC1       = P_LC1./norm_fac;
P_LC2       = P_LC2./norm_fac;

%% output 
varargout{1}=P_LC1;
varargout{2}=P_LC2;
end