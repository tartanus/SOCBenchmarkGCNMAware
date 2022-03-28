function OV=OVCalc(y,time,refAmp)
    OV=abs((max(abs(y))-refAmp)*100);
 end