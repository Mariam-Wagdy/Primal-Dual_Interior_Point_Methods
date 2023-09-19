function [x0, s0, lambda0]=Initial(A,b,c)

[m,n]=size(A);
At=transpose(A);

x_star=At*((A*At)\b);
lambda_star=((A*At)\A)*c;
s_star=c-At*lambda_star;

%% ensure nonnegative values%%
sigma_x= max(-3/2*min(x_star), 0);
sigma_s= max(-3/2*min(s_star), 0);

en1=ones(n,1);
x_hat=x_star+sigma_x*en1;
s_hat=s_star+sigma_s*en1;

%% ensure not too close to zero and not too dissimilar%%
e1n=ones(1,n);
sigma_x_hat=0.5*(transpose(x_hat)*s_hat)/(e1n*s_hat);
sigma_s_hat=0.5*(transpose(x_hat)*s_hat)/(e1n*x_hat);

%% initial point%%
x0=x_hat+sigma_x_hat*en1;
s0=s_hat+sigma_s_hat*en1;
lambda0=lambda_star;
end
