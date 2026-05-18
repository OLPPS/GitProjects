
%A = rand(100000,1) * 120;
%bns = [0:10:100,inf]
%B = SpecAnalysis2(A,1,bns);

%figure(1)
%bar(bns(1:end-1),B,'histc')

function [Cnts] = SpecAnalysis(PixEvents, PixCIEs, bins)

bins = [bins, inf];

Sigs = PixEvents(:).*PixCIEs(:);

Cnts = histcounts(Sigs,bins);

end