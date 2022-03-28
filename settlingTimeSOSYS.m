%% settling time measure function
%y: system output
%w: sliding window size
function ts=settlingTimeSOSYS(y,w,time,epsilon)

    for i=(length(y)):-1:w
        
        if abs(std(y(i-w:i)))>epsilon
            ts=time(i);
            break
        else
            ts=time(end);
        end
        
        ts;
    end




end

