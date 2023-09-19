function [primal,dual,xstar]=CenPa(A,b,c,alpha,si,a)
[m,n] = size(A);
At=transpose(A);
en1 = ones(n,1);

In=eye(n);
z11=zeros(n);
z22=zeros(m);
z23=zeros(m,n);
z32=zeros(n,m);
zn1=zeros(n,1);
zm1=zeros(m,1);

%delta_step=linsolve(LHS,RHS);
%delta_step=inverse(LHS)*RHS;
[x,s,y] =Initial(A,b,c);
k=1;
condition=1;
while condition
    step= [x(:,k);y(:,k);s(:,k)];
    
    p(k)=transpose(c)*x(:,k);
    d(k)=transpose(b)*y(:,k);
    
    mu(k) = transpose(x(:,k))*s(:,k)/n;
    condition= mu(k)> 0.0001;
    
    S = diag(s(:,k));
    X = diag(x(:,k));
    
    rc=At*y(:,k)+s(:,k)-c;
    rb=A*x(:,k)-b;
    
    rXS=x(:,k).*s(:,k)-si*mu(k)*en1;
    
    LHS = [z11 transpose(A) In;A z22 z23; S z32 X];
    RHS = [-rc;-rb;-rXS];
    
    dLHS = decomposition(LHS);
    delta_step(:,k) = dLHS\RHS;
    if a==0 %% fixed, user defined alpa and sigma
        next_step(:,k) = step+alpha*delta_step(:,k);
        
        x(:,k+1) = next_step(1:n,k);
        y(:,k+1) = next_step(n+1:n+m,k);
        s(:,k+1) = next_step(n+m+1:2*n+m,k);
    else %% adaptive alpa and sigma
       
        alpha=1; %% begin with maximum step size
        next_step(:,k) = step+alpha*delta_step(:,k);
        x(:,k+1) = next_step(1:n,k);
        y(:,k+1) = next_step(n+1:n+m,k);
        s(:,k+1) = next_step(n+m+1:2*n+m,k);
        while any(x(:,k).*s(:,k)<0) %% ensure positivity condition holds
            e=min(abs(x(:,k)))/1000; %% difference in step size is 0.001 of the smallest variable in x
            alpha=alpha-e;
            next_step(:,k) = step+alpha*delta_step(:,k);
            x(:,k+1) = next_step(1:n,k);
            y(:,k+1) = next_step(n+1:n+m,k);
            s(:,k+1) = next_step(n+m+1:2*n+m,k);
        end
        
    end
    k=k+1;
end

primal=p(k-1);
dual=d(k-1);
xstar=x(:,k-1);
figure
plot([1:k-1],p,'-',[1:k-1],d,'--');
legend('cTx primal function result','bTy dual function result')
xlabel('Iterations')
ylabel('f(x)')

figure
plot(x(1,:), x(2,:),'b--o');
hold on
for j=1:m
    x1=[0, b(j)/A(j,1)];
    x2=[b(j)/A(j,2),0];
    plot(x1,x2);
end
hold off
end
%%% fixed
%%% adaptive a & s