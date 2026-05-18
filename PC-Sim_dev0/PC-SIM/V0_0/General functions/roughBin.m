%EventsIn = rand(10,1) * 100
%TimesIn = [0:10:90];
%WW = 10;

%output = roughBin2(EventsIn,TimesIn,WW)

function RoughBinnedEvents = roughBin(SIGNALS,TIMES,WindowWidth)
sigsOut = zeros(size(SIGNALS(:),1),1);
s = 1;

tempsig = 0;
curTimeWindow = 0;

t = 1;
while t <= size(TIMES(:),1)
    if TIMES(t) - curTimeWindow < WindowWidth
        tempsig = tempsig + SIGNALS(t);
        t = t+1;
    else
        sigsOut(s) = tempsig;
        s = s+1;
        tempsig = 0;
        curTimeWindow = TIMES(t);
    end
end

RoughBinnedEvents = sigsOut(1:s-1);

end