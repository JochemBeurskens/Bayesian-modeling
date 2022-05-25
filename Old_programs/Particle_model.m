%% Algorithm 1 MC for expected value of P(y)=N(2,5), 10000 draws:
y=[]; %Make list, y, to store the random draws
N=10000;
mu=2;
sig=5;
for i=1:N
    yi=normrnd(mu,sig);
    y=[y,yi];
end
exv=(1/N)*sum(y);
%% Algorithm 2 importance sampling for the expected value, performed on problem of example 1.1:
mu=0.4; %mean of normal distribution
sig=0.1; %std of normal distribution
tau=0.5; %rate of exponential
% firstly we estimate the probability that the value is larger than 3 by
% using alogorithm 1:
y=[]; %Make list, y, to store the random draws
N=2000;
for i=1:N
    yi=normrnd(mu,sig);
    y=[y,yi];
end
e=exprnd(tau,1,N); %get different result from author, here the mean is 1/rate, so should be 2, then get different distribution. But if I use mean of 0.5 I get the author's result
y=y+e;
prob_larger=sum(y>3.0)/N;
%now using importance sampling with an importance distr also exponential:
imp=exprnd(tau,1,N)+3.0; %sampled from the given importance distribution
p_of_imp=exppdf(imp,tau);
weight=p_of_imp./imp;
exv_IS=(1/N)*sum(weight.*imp); %get 0.0055 only when adding 2.6 instead of 3 to the imp function, but do get a lot less variation
% checking my result:
x = 3;
n=1-normcdf(x,mu,sig);
ee=1-expcdf(x,tau);
est=n+ee; %now have a cumulative distr for the exponential gaussian
cest=sum(est); %my estimate is indeed very correct