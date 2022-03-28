%perform optimization with the modified nelder mead for the convergence
%analysis
%parameters initialization
clear; clc; close all;
kpo=1;  %10  0.06
kio=0.5; %3    0.05
%optimization gains
maxKp=1; minKp=0.01;
maxKi=1; minKi=0.001;
%Globalized Constrained Nelder Mead GCNM parameters
optimPeriod=300;        %optimization execution period (s)
resetThreshold=0.1;     %probabilistic restart threshold
refSignalPeriod=300;    %Reference signal period
refAmplitude=1;        %reference amplitude
%constraints limits
OVLim=2;             %5% overshoot
TSLim=20;               %Settling time limit (s)
L=0.1;                  %delay value L=0.1,1,10;
N=100;                  %derivative filter N (Zero for PI controller)
timeVec=0:0.01:30000;
x0_1=[1,1];
constraints=1;
% JG_1=10000;
ts=0.1;
t=0:ts:300;  
restartsCount=0;
%% GCNM SOC with NM convergent variant 
for repetitions=1:10
%     figure()
    restartsCount=0;
    constraints=1;
    while constraints==1  
        globalSearch=1;
        %global search stage
        while globalSearch==1
            x0Aux=[maxKp*rand(), maxKi*rand()];
            [J,OV,SetTime,yout]=testBenchSOCNMConverGS(x0Aux);
            [J_1,OV_1,SetTime_1,yout]=testBenchSOCNMConverGS(x0_1);

            if J<J_1
                x0=x0Aux;
                globalSearch=0;
                restartsCount=restartsCount+1;
            end
            if (OV>OVLim-1 && OV<OVLim) && (SetTime>TSLim-1  && SetTime<TSLim)%OV<OVLim && SetTime<TSLim
                constraints=0;           %stopping condition for all the algorithm (improbable)   
            end
            %update random search variables
            J_1=J;
            x0_1=x0Aux;
            x0Hist=x0;
        end
        
        if constraints==0
            break
        end
        
        %execute local search with convergent NM
        [x,fval,xHist,fHist,exitflag,output]= NM4_s3d6k1n8(@testBenchSOCNMConver,x0);
        %evaluate system performance after local search
        [J,OV,SetTime,yout]=testBenchSOCNMConverGS(x);



        if (OV>0 && OV<OVLim) && (SetTime>TSLim-0.5*TSLim  && SetTime<TSLim)
            constraints=0;           %stopping condition for all the algorithm (improbable)   
        end

               restartsCount-1
                [J,OV,SetTime]
                xHistPRS(:,restartsCount)=x;
    end
    
    figure()
    subplot(2,1,1)
    plot(t,yout)
    hold on
    subplot(2,1,2)
    plot(xHist(1,:))
    hold on
    plot(xHist(2,:))
    legend('kp','ki')
    
    
    restartsCountHistRep(repetitions)=restartsCount-1;
    SetTimeHistRep(repetitions)=   SetTime;
    OVHistRep(repetitions)=      OV;
    JHistRep(repetitions)=    J;
    xHistRep(:,repetitions)=x;
    NMIterHistRep(repetitions)=output.iterations;
    NMfuncCountHistRep(repetitions)=output.funcCount;
    
    
end
%%




% f=testBenchSOCNMConver([5,5])


% 
% options = optimoptions('patternsearch','UseCompletePoll',true);
% [x fval] = patternsearch(@my_funTest,[1 1],[],[],[],[],[],[],...
%     [],options);
% 
% options = optimoptions('ga');
% [x fval] = gamultiobj(@my_funTest,2,[],[],[],[],[],[],options);