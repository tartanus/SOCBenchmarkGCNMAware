function [J,OV,SetTime,yout]=testBenchSOCNMConverGS(x)
theta=x;
s=tf('s');
ts=0.1;
t=0:ts:300;     %time vector
kp=theta(1);
ki=theta(2);
L=1;
P=(1/(s+1))*exp(-L*s);
C=kp+ki/s;
f=1/300;          %repetition frequency
lambda=1/0.5e3;      %regularization coefficient
r=sign(sin(2*pi*f*t));  %reference signal
T1=feedback(P*C,1);
yout=lsim(T1,r,t);
OV=OVCalc(yout,t,1);
 try
    SetTime=settlingTimeSOSYS(yout(1:round(end/2)),5,t(1:round(end/2))',1e-3);
 catch
     SetTime=300;
 end
 J=0.1*(1/(2))*sum( (r-yout').^2  )+ 5*OV+SetTime;% + lambda*((theta(1))^2+(theta(2))^2);
end
    






