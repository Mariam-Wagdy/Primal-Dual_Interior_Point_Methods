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
        rXS=x(:,k).*s(:,k)-si*mu(k)*en1;
        
        LHS = [z11 transpose(A) In;A z22 z23; S z32 X];
        RHS = [-rc;-rb;-rXS];
        delta_step(:,k) = dLHS\RHS;
        
        delta_x(:,k)=delta_step(1:n,k);
        delta_y(:,k)=delta_step(n+1:n+m,k);
        delta_s(:,k)=delta_step(n+m+1:2*n+m,k);
        
        RHS_aff = [-rc;-rb;-x(:,k).*s(:,k)];
        delta_step_aff(:,k) = dLHS\RHS_aff;
        
        delta_x_aff(:,k)=delta_step_aff(1:n,k);
        delta_y_aff(:,k)=delta_step_aff(n+1:n+m,k);
        delta_s_aff(:,k)=delta_step_aff(n+m+1:2*n+m,k);
        
       
        ratio_x= -min(min(0,x(:,k)./delta_x_aff(:,k)));
        ratio_s= -min(min(0,s(:,k)./delta_s_aff(:,k)));
        a_pri_aff=min(1,ratio_x);
        a_dual_aff=min(1,ratio_s);
   
        mu_aff=transpose(x(:,k)+a_pri_aff*delta_x_aff(:,k)).*(s(:,k)+a_dual_aff*delta_s_aff(:,k))/n;
        si=(mu_aff/mu(k))^3;
        
      
        x(:,k+1) = x(:,k)+a_pri_aff*delta_x(:,k);
        y(:,k+1) = y(:,k)+a_dual_aff*delta_y(:,k);
        s(:,k+1) = s(:,k)+a_dual_aff*delta_s(:,k);
        
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
title('Duality Test')

figure
plot(x(1,:), x(2,:),'b--o');
hold on
for j=1:m
    x1=[0, b(j)/A(j,1)];
    x2=[b(j)/A(j,2),0];
    plot(x1,x2);
end
hold off
title('Feasible Region and steps taken')

xx=x(1,1:k-1).*s(1,1:k-1);
yy=x(2,1:k-1).*s(2,1:k-1);
figure
plot(xx,yy,'b--o');
xlabel('x1s1')
ylabel('x2s2')
title('Central Path')

figure
plot(1:k-1,mu,'b--o');
xlabel('iterations')
ylabel('mu')
title('Duality Measure')

end
%%% fixed
%%% adaptive a & s