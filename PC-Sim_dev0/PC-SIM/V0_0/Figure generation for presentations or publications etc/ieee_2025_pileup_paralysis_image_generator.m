counts(1:8) = 0;
jj = 1;
for N = [10000,100000,1000000,10000000,20000000,50000000,70000000,100000000]
N
Events = rand(N,1);
Sorted = sort(Events);

for ii = 1:N-1
if Sorted(ii+1) - Sorted(ii) > 10^(-7.5)
    counts(jj) = counts(jj)+1;
end
end
jj = jj+1;
end
%%
figure(1)
hold off
scatter([10000,100000,1000000,10000000,20000000,50000000,70000000,100000000],counts(1:8))
hold on
%plot([10000,100000,1000000,10000000,20000000,50000000,70000000,100000000],[10000,100000,1000000,10000000,20000000,50000000,70000000,100000000])
hold off
%%
y = counts(1:8);
x = [10000,100000,1000000,10000000,20000000,50000000,70000000,100000000];
% Get coefficients of a line fit through the data.
coefficients = polyfit(x, y, 4);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% Plot everything.
plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
hold on; % Set hold on so the next plot does not blow away the one we just drew.
plot(x, y, 'b.', 'MarkerSize', 15); % Plot training data.
hold off

xlabel('x-ray flux (photons mm^-^2 s^-^1')
ylabel('Counts recorded')