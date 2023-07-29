function [t,f]=simulation_SIQR_twocontrol(A,lambda,beta,delta,mu,r,epsilon,d,psi,S0,I0,Q0,R0,T,A1,A2)
N=1000;
t=linspace(0,T,N+1);
h=T/N;
h2=h/2;

%Initial value of state
S=zeros(1,N+1);
I=zeros(1,N+1);
Q=zeros(1,N+1);
R=zeros(1,N+1);

%Initial value of state, co-state, and control
stateS=zeros(1,N+1);
stateI=zeros(1,N+1);
stateQ=zeros(1,N+1);
stateR=zeros(1,N+1);
co_stateS=zeros(1,N+1);
co_stateI=zeros(1,N+1);
co_stateQ=zeros(1,N+1);
co_stateR=zeros(1,N+1);
u=zeros(1,N+1);

% Two control
u1=zeros(1,N+1);
u2=zeros(1,N+1);
two_stateS=zeros(1,N+1);
two_stateI=zeros(1,N+1);
two_stateQ=zeros(1,N+1);
two_stateR=zeros(1,N+1);
two_co_stateS=zeros(1,N+1);
two_co_stateI=zeros(1,N+1);
two_co_stateQ=zeros(1,N+1);
two_co_stateR=zeros(1,N+1);


%State without control
S(1)=S0;
I(1)=I0;
Q(1)=Q0;
R(1)=R0;

%State with control
stateS(1)=S0;
stateI(1)=I0;
stateQ(1)=Q0;
stateR(1)=R0;

%State with two control
two_stateS(1)=S0;
two_stateI(1)=I0;
two_stateQ(1)=Q0;
two_stateR(1)=R0;

iteration=0;
for i=1:N
    iteration=iteration+1;
    oldu=u;
    oldu1=u1;
    oldu2=u2;
    for i=1:N
        valueS=lambda-beta*delta*S(i)*I(i)-mu*S(i);
        valueI=beta*delta*S(i)*I(i)-(r+epsilon+mu+d)*I(i);
        valueQ=epsilon*I(i)-(psi+d+mu)*Q(i);
        valueR=r*I(i)+psi*Q(i)-mu*R(i);
        %------------------------------------------------------------------
        value2S=lambda-beta*delta*(S(i)+h2*valueS)*(I(i)+h2*valueI)-mu*(S(i)+h2*valueS);
        value2I=beta*delta*(S(i)+h2*valueS)*(I(i)+h2*valueI)-(r+epsilon+mu+d)*(I(i)+h2*valueI);
        value2Q=epsilon*(I(i)+h2*valueI)-(psi+d+mu)*(Q(i)+h2*valueQ);
        value2R=r*(I(i)+h2*valueI)+psi*(Q(i)+h2*valueQ)-mu*(R(i)+h2*valueR);
        %------------------------------------------------------------------
        value3S=lambda-beta*delta*(S(i)+h2*value2S)*(I(i)+h2*value2I)-mu*(S(i)+h2*value2S);
        value3I=beta*delta*(S(i)+h2*value2S)*(I(i)+h2*value2I)-(r+epsilon+mu+d)*(I(i)+h2*value2I);
        value3Q=epsilon*(I(i)+h2*value2I)-(psi+d+mu)*(Q(i)+h2*value2Q);
        value3R=r*(I(i)+h2*value2I)+psi*(Q(i)+h2*value2Q)-mu*(R(i)+h2*value2R);
        %------------------------------------------------------------------
        value4S=lambda-beta*delta*(S(i)+h*value3S)*(I(i)+h*value3I)-mu*(S(i)+h*value3S);
        value4I=beta*delta*(S(i)+h*value3S)*(I(i)+h*value3I)-(r+epsilon+mu+d)*(I(i)+h*value3I);
        value4Q=epsilon*(I(i)+h*value3I)-(psi+d+mu)*(Q(i)+h*value3Q);
        value4R=r*(I(i)+h*value3I)+psi*(Q(i)+h*value3Q)-mu*(R(i)+h*value3R);
        %------------------------------------------------------------------
        S(i+1)=S(i)+(h/6)*(valueS+2*value2S+2*value3S+value4S);
        I(i+1)=I(i)+(h/6)*(valueI+2*value2I+2*value3I+value4I);
        Q(i+1)=Q(i)+(h/6)*(valueQ+2*value2Q+2*value3Q+value4Q);
        R(i+1)=R(i)+(h/6)*(valueR+2*value2R+2*value3R+value4R);
    end
    %% state one control
    for i=1:N
        value_stateS=lambda-beta*delta*stateS(i)*stateI(i)*(1-u(i))-mu*stateS(i);
        value_stateI=beta*delta*stateS(i)*stateI(i)*(1-u(i))-(r+epsilon+mu+d)*stateI(i);
        value_stateQ=epsilon*stateI(i)-(psi+d+mu)*stateQ(i);
        value_stateR=r*stateI(i)+psi*stateQ(i)-mu*stateR(i);
        %------------------------------------------------------------------
        value_state2S=lambda-beta*delta*(stateS(i)+h2*value_stateS)*(stateI(i)+h2*value_stateI)*(1-0.5*(u(i)+u(i+1)))-mu*(stateS(i)+h2*value_stateS);
        value_state2I=(1-0.5*(u(i)+u(i+1)))*beta*delta*(stateS(i)+h2*value_stateS)*(stateI(i)+h2*value_stateI)-(r+epsilon+mu+d)*(stateI(i)+h2*value_stateI);
        value_state2Q=epsilon*(stateI(i)+h2*value_stateI)-(psi+d+mu)*(stateQ(i)+h2*value_stateQ);
        value_state2R=r*(stateI(i)+h2*value_stateI)+psi*(stateQ(i)+h2*value_stateQ)-mu*(stateR(i)+h2*value_stateR);
        %------------------------------------------------------------------
        value_state3S=lambda-(1-0.5*(u(i)+u(i+1)))*beta*delta*(stateS(i)+h2*value_state2S)*(stateI(i)+h2*value_state2I)-mu*(stateS(i)+h2*value_state2S);
        value_state3I=(1-0.5*(u(i)+u(i+1)))*beta*delta*(stateS(i)+h2*value_state2S)*(stateI(i)+h2*value_state2I)-(r+epsilon+mu+d)*(stateI(i)+h2*value_state2I);
        value_state3Q=epsilon*(stateI(i)+h2*value_state2I)-(psi+d+mu)*(stateQ(i)+h2*value_state2Q);
        value_state3R=r*(stateI(i)+h2*value_state2I)+psi*(stateQ(i)+h2*value_state2Q)-mu*(stateR(i)+h2*value_state2R);
        %------------------------------------------------------------------
        value_state4S=lambda-(1-u(i))*beta*delta*(stateS(i)+h*value_state3S)*(stateI(i)+h*value_state3I)-mu*(stateS(i)+h*value_state3S);
        value_state4I=(1-u(i))*beta*delta*(stateS(i)+h*value_state3S)*(stateI(i)+h*value_state3I)-(r+epsilon+mu+d)*(stateI(i)+h*value_state3I);
        value_state4Q=epsilon*(stateI(i)+h*value_state3I)-(psi+d+mu)*(stateQ(i)+h*value_state3Q);
        value_state4R=r*(stateI(i)+h*value_state3I)+psi*(stateQ(i)+h*value_state3Q)-mu*(stateR(i)+h*value_state3R);
        %------------------------------------------------------------------
        stateS(i+1)=stateS(i)+(h/6)*(value_stateS+2*value_state2S+2*value_state3S+value_state4S);
        stateI(i+1)=stateI(i)+(h/6)*(value_stateI+2*value_state2I+2*value_state3I+value_state4I);
        stateQ(i+1)=stateQ(i)+(h/6)*(value_stateQ+2*value_state2Q+2*value_state3Q+value_state4Q);
        stateR(i+1)=stateR(i)+(h/6)*(value_stateR+2*value_state2R+2*value_state3R+value_state4R);
    end
    %% co-state one control
    for i=1:N
        j=N+2-i;
        value_stateS=mu*co_stateS(j)+(1-u(j))*beta*delta*stateI(j)*(co_stateS(j)-co_stateI(j));
        value_stateI= -1+(1-u(j))*beta*delta*stateS(j)*(co_stateS(j)-co_stateI(j))...
            +(r+epsilon+mu+d)*co_stateI(j)-epsilon*co_stateQ(j)-r*co_stateR(j);
        value_stateQ=co_stateQ(j)*(psi+d+mu)-co_stateR(j)*psi;
        value_stateR=co_stateR(j)*mu;
        %------------------------------------------------------------------
        value_state2S=mu*(co_stateS(j)-h2*value_stateS)+(1-0.5*(u(j)+u(j-1)))*beta*delta*0.5*(stateI(j)+stateI(j-1))*((co_stateS(j)-h2*value_stateS)-(co_stateI(j)-h2*value_stateI));
        value_state2I= -1+(1-0.5*(u(j)+u(j-1)))*beta*delta*0.5*(stateI(j)+stateI(j-1))*((co_stateS(j)-h2*value_stateS)-(co_stateI(j)-h2*value_stateI))...
            +(r+epsilon+mu+d)*(co_stateI(j)-h2*value_stateI)-epsilon*(co_stateQ(j)-h2*value_stateQ)-r*(co_stateR(j)-h2*value_stateR);
        value_state2Q=(co_stateQ(j)-h2*value_stateQ)*(psi+d+mu)-(co_stateR(j)-h2*value_stateR)*psi;
        value_state2R=(co_stateR(j)-h2*value_stateR)*mu;
        %------------------------------------------------------------------
        value_state3S=mu*(co_stateS(j)-h2*value_state2S)+(1-0.5*(u(j)+u(j-1)))*beta*delta*0.5*(stateI(j)+stateI(j-1))*((co_stateS(j)-h2*value_state2S)-(co_stateI(j)-h2*value_state2I));
        value_state3I= -1+(1-0.5*(u(j)+u(j-1)))*beta*delta*0.5*(stateI(j)+stateI(j-1))*((co_stateS(j)-h2*value_state2S)-(co_stateI(j)-h2*value_state2I))...
            +(r+epsilon+mu+d)*(co_stateI(j)-h2*value_state2I)-epsilon*(co_stateQ(j)-h2*value_state2Q)-r*(co_stateR(j)-h2*value_state2R);
        value_state3Q=(co_stateQ(j)-h2*value_state2Q)*(psi+d+mu)-(co_stateR(j)-h2*value_state2R)*psi;
        value_state3R=(co_stateR(j)-h2*value_state2R)*mu;
        %------------------------------------------------------------------
        value_state4S=mu*(co_stateS(j)-h*value_state3S)+(1-u(j))*beta*delta*stateI(j-1)*((co_stateS(j)-h*value_state3S)-(co_stateI(j)-h*value_state3I));
        value_state4I= -1+(1-u(j))*beta*delta*stateI(j-1)*((co_stateS(j)-h*value_state3S)-(co_stateI(j)-h*value_state3I))...
            +(r+epsilon+mu+d)*(co_stateI(j)-h*value_state3I)-epsilon*(co_stateQ(j)-h*value_state3Q)-r*(co_stateR(j)-h*value_state3R);
        value_state4Q=(co_stateQ(j)-h*value_state3Q)*(psi+d+mu)-(co_stateR(j)-h*value_state3R)*psi;
        value_state4R=(co_stateR(j)-h*value_state3R)*mu;
        %------------------------------------------------------------------
        co_stateS(j-1)=co_stateS(j)-(h/6)*(value_stateS+2*value_state2S...
            +2*value_state3S+value_state4S);
        co_stateI(j-1)=co_stateI(j)-(h/6)*(value_stateI+2*value_state2I...
            +2*value_state3I+value_state4I);
        co_stateQ(j-1)=co_stateQ(j)-(h/6)*(value_stateQ+2*value_state2Q...
            +2*value_state3Q+value_state4Q);
        co_stateR(j-1)=co_stateR(j)-(h/6)*(value_stateR+2*value_state2R...
            +2*value_state3R+value_state4R);
        %------------------------------------------------------------------
        gh(j)=beta*delta*stateS(j)*stateI(j)...
            *(-co_stateS(j)+co_stateI(j));
        u(j)=min(1,max(0,gh(j)/A));
    end
    %% state two control
    for i=1:N
        two_value_stateS=(1-u1(i))*lambda-beta*delta*two_stateS(i)*two_stateI(i)*(1-u2(i))-mu*two_stateS(i);
        two_value_stateI=beta*delta*two_stateS(i)*two_stateI(i)*(1-u2(i))-(r+epsilon+mu+d)*two_stateI(i);
        two_value_stateQ=epsilon*two_stateI(i)-(psi+d+mu)*two_stateQ(i);
        two_value_stateR=u1(i)*lambda+r*two_stateI(i)+psi*two_stateQ(i)-mu*two_stateR(i);
        %------------------------------------------------------------------
        two_value_state2S=(1-0.5*(u1(i)+u1(i+1)))*lambda-beta*delta*(two_stateS(i)+h2*two_value_stateS)*(two_stateI(i)+h2*two_value_stateI)*(1-0.5*(u2(i)+u2(i+1)))-mu*(two_stateS(i)+h2*two_value_stateS);
        two_value_state2I=(1-0.5*(u2(i)+u2(i+1)))*beta*delta*(two_stateS(i)+h2*two_value_stateS)*(two_stateI(i)+h2*two_value_stateI)-(r+epsilon+mu+d)*(two_stateI(i)+h2*two_value_stateI);
        two_value_state2Q=epsilon*(two_stateI(i)+h2*two_value_stateI)-(psi+d+mu)*(two_stateQ(i)+h2*two_value_stateQ);
        two_value_state2R=0.5*(u1(i)+u1(i+1))*lambda+r*(two_stateI(i)+h2*two_value_stateI)+psi*(two_stateQ(i)+h2*two_value_stateQ)-mu*(two_stateR(i)+h2*two_value_stateR);
        %------------------------------------------------------------------
        two_value_state3S=(1-0.5*(u1(i)+u1(i+1)))*lambda-(1-0.5*(u2(i)+u2(i+1)))*beta*delta*(two_stateS(i)+h2*two_value_state2S)*(two_stateI(i)+h2*two_value_state2I)-mu*(two_stateS(i)+h2*two_value_state2S);
        two_value_state3I=(1-0.5*(u2(i)+u2(i+1)))*beta*delta*(two_stateS(i)+h2*two_value_state2S)*(two_stateI(i)+h2*two_value_state2I)-(r+epsilon+mu+d)*(two_stateI(i)+h2*two_value_state2I);
        two_value_state3Q=epsilon*(two_stateI(i)+h2*two_value_state2I)-(psi+d+mu)*(two_stateQ(i)+h2*two_value_state2Q);
        two_value_state3R=0.5*(u1(i)+u1(i+1))*lambda+r*(two_stateI(i)+h2*two_value_state2I)+psi*(two_stateQ(i)+h2*two_value_state2Q)-mu*(two_stateR(i)+h2*two_value_state2R);
        %------------------------------------------------------------------
        two_value_state4S=(1-u1(i))*lambda-(1-u2(i))*beta*delta*(two_stateS(i)+h*two_value_state3S)*(two_stateI(i)+h*two_value_state3I)-mu*(two_stateS(i)+h*two_value_state3S);
        two_value_state4I=(1-u2(i))*beta*delta*(two_stateS(i)+h*two_value_state3S)*(two_stateI(i)+h*two_value_state3I)-(r+epsilon+mu+d)*(two_stateI(i)+h*two_value_state3I);
        two_value_state4Q=epsilon*(two_stateI(i)+h*two_value_state3I)-(psi+d+mu)*(two_stateQ(i)+h*two_value_state3Q);
        two_value_state4R=u1(i)*lambda+r*(two_stateI(i)+h*two_value_state3I)+psi*(two_stateQ(i)+h*two_value_state3Q)-mu*(two_stateR(i)+h*two_value_state3R);
        %------------------------------------------------------------------
        two_stateS(i+1)=two_stateS(i)+(h/6)*(two_value_stateS+2*two_value_state2S+2*two_value_state3S+two_value_state4S);
        two_stateI(i+1)=two_stateI(i)+(h/6)*(two_value_stateI+2*two_value_state2I+2*two_value_state3I+two_value_state4I);
        two_stateQ(i+1)=two_stateQ(i)+(h/6)*(two_value_stateQ+2*two_value_state2Q+2*two_value_state3Q+two_value_state4Q);
        two_stateR(i+1)=two_stateR(i)+(h/6)*(two_value_stateR+2*two_value_state2R+2*two_value_state3R+two_value_state4R);
    end
    %% co-state two control
    for i=1:N
        j=N+2-i;
        two_value_stateS=mu*two_co_stateS(j)+(1-u2(j))*beta*delta*two_stateI(j)*(two_co_stateS(j)-two_co_stateI(j));
        two_value_stateI= -1+(1-u2(j))*beta*delta*two_stateS(j)*(two_co_stateS(j)-two_co_stateI(j))...
            +(r+epsilon+mu+d)*two_co_stateI(j)-epsilon*two_co_stateQ(j)-r*two_co_stateR(j);
        two_value_stateQ=two_co_stateQ(j)*(psi+d+mu)-two_co_stateR(j)*psi;
        two_value_stateR=two_co_stateR(j)*mu;
        %------------------------------------------------------------------
        two_value_state2S=mu*(two_co_stateS(j)-h2*two_value_stateS)+(1-0.5*(u2(j)+u2(j-1)))*beta*delta*0.5*(two_stateI(j)+two_stateI(j-1))*((two_co_stateS(j)-h2*two_value_stateS)-(two_co_stateI(j)-h2*two_value_stateI));
        two_value_state2I= -1+(1-0.5*(u2(j)+u2(j-1)))*beta*delta*0.5*(two_stateS(j)+two_stateS(j-1))*((two_co_stateS(j)-h2*two_value_stateS)-(two_co_stateI(j)-h2*two_value_stateI))...
            +(r+epsilon+mu+d)*(two_co_stateI(j)-h2*two_value_stateI)-epsilon*(two_co_stateQ(j)-h2*two_value_stateQ)-r*(two_co_stateR(j)-h2*two_value_stateR);
        two_value_state2Q=(two_co_stateQ(j)-h2*two_value_stateQ)*(psi+d+mu)-(two_co_stateR(j)-h2*two_value_stateR)*psi;
        two_value_state2R=(two_co_stateR(j)-h2*two_value_stateR)*mu;
        %------------------------------------------------------------------
        two_value_state3S=mu*(two_co_stateS(j)-h2*two_value_state2S)+(1-0.5*(u2(j)+u2(j-1)))*beta*delta*0.5*(two_stateI(j)+two_stateI(j-1))*((two_co_stateS(j)-h2*two_value_state2S)-(two_co_stateI(j)-h2*two_value_state2I));
        two_value_state3I= -1+(1-0.5*(u2(j)+u2(j-1)))*beta*delta*0.5*(two_stateS(j)+two_stateS(j-1))*((two_co_stateS(j)-h2*two_value_state2S)-(two_co_stateI(j)-h2*two_value_state2I))...
            +(r+epsilon+mu+d)*(two_co_stateI(j)-h2*two_value_state2I)-epsilon*(two_co_stateQ(j)-h2*two_value_state2Q)-r*(two_co_stateR(j)-h2*two_value_state2R);
        two_value_state3Q=(two_co_stateQ(j)-h2*two_value_state2Q)*(psi+d+mu)-(two_co_stateR(j)-h2*two_value_state2R)*psi;
        two_value_state3R=(two_co_stateR(j)-h2*two_value_state2R)*mu;
        %------------------------------------------------------------------
        two_value_state4S=mu*(two_co_stateS(j)-h*two_value_state3S)+(1-u2(j))*beta*delta*two_stateI(j-1)*((two_co_stateS(j)-h*two_value_state3S)-(two_co_stateI(j)-h*two_value_state3I));
        two_value_state4I= -1+(1-u2(j))*beta*delta*two_stateS(j-1)*((two_co_stateS(j)-h*two_value_state3S)-(two_co_stateI(j)-h*two_value_state3I))...
            +(r+epsilon+mu+d)*(two_co_stateI(j)-h*two_value_state3I)-epsilon*(two_co_stateQ(j)-h*two_value_state3Q)-r*(two_co_stateR(j)-h*two_value_state3R);
        two_value_state4Q=(two_co_stateQ(j)-h*two_value_state3Q)*(psi+d+mu)-(two_co_stateR(j)-h*two_value_state3R)*psi;
        two_value_state4R=(two_co_stateR(j)-h*two_value_state3R)*mu;
        %------------------------------------------------------------------
        two_co_stateS(j-1)=two_co_stateS(j)-(h/6)*(two_value_stateS+2*two_value_state2S...
            +2*two_value_state3S+two_value_state4S);
        two_co_stateI(j-1)=two_co_stateI(j)-(h/6)*(two_value_stateI+2*two_value_state2I...
            +2*two_value_state3I+two_value_state4I);
        two_co_stateQ(j-1)=two_co_stateQ(j)-(h/6)*(two_value_stateQ+2*two_value_state2Q...
            +2*two_value_state3Q+two_value_state4Q);
        two_co_stateR(j-1)=two_co_stateR(j)-(h/6)*(two_value_stateR+2*two_value_state2R...
            +2*two_value_state3R+two_value_state4R);
        %------------------------------------------------------------------
        gh1(j)=lambda*(-two_co_stateR(j)+two_co_stateS(j));
        gh2(j)=beta*delta*two_stateS(j)*two_stateI(j)...
            *(-two_co_stateS(j)+two_co_stateI(j));
        u1(j)=min(1,max(0,gh1(j)/A1));
        u2(j)=min(1,max(0,gh2(j)/A2));
    end
    temp=(beta*delta*stateS.*stateI.*(-co_stateS+co_stateI))/A;
    ul=min(1,max(0,temp));
    u=0.5*(ul+oldu);
    
    temp1=(lambda.*(-co_stateR+co_stateS))/A1;
    u11=min(1,max(0,temp1));
    u1=0.5*(u11+oldu1);
    
    temp2=(beta*delta*stateS.*stateI.*(-co_stateS+co_stateI))/A2;
    ul2=min(1,max(0,temp2));
    u2=0.5*(ul2+oldu2);
end
f(1,:)=S;
f(2,:)=I;
f(3,:)=Q;
f(4,:)=R;
f(5,:)=stateS;
f(6,:)=stateI;
f(7,:)=stateQ;
f(8,:)=stateR;
f(9,:)=two_stateS;
f(10,:)=two_stateI;
f(11,:)=two_stateQ;
f(12,:)=two_stateR;
f(13,:)=u;
f(14,:)=u1;
f(15,:)=u2;
end