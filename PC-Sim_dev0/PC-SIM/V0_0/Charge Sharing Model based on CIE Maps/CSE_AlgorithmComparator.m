clear
close all
clc


data = dlmread('C:\Users\oliep\OneDrive\Work\PDTF 3 - 2023-2027 (EPSRC SPCS Project)\SPCS codes\edgePixel');
[result_edge, Y_grid_edge, Z_grid_edge, Value_grid_edge, Y_unique_edge, Z_unique_edge, Y_increments_edge, Z_increments_edge, Points_assessed] = dataExtractor(data,0,0,2);

data = dlmread('C:\Users\oliep\OneDrive\Work\PDTF 3 - 2023-2027 (EPSRC SPCS Project)\SPCS codes\cornerPixel');
[result_corner, Y_grid_corner, Z_grid_corner, Value_grid_corner, Y_unique_corner, Z_unique_corner, Y_increments_corner, Z_increments_corner, Points_assessed] = dataExtractor(data,0,0,3);

data = dlmread('C:\Users\oliep\OneDrive\Work\PDTF 3 - 2023-2027 (EPSRC SPCS Project)\SPCS codes\centrePixel');
[result_centre, Y_grid_centre, Z_grid_centre, Value_grid_centre, Y_unique, Z_unique, Y_increments, Z_increments, Points_assessed] = dataExtractor(data,0,0,1);


%% Newer attempt based on tracking all pixels

Pix(1:9,:,:) = zeros(9,size(Value_grid_centre,1),size(Value_grid_centre,2));

Pix(5,:,:) = Value_grid_centre; 

Pix(8,:,:) = Value_grid_edge;
Pix(2,:,:) = rot90(Value_grid_edge,2);
Pix(4,:,:) = rot90(Value_grid_edge,3);
Pix(6,:,:) = rot90(Value_grid_edge,1);

Pix(7,:,:) = Value_grid_corner;
Pix(9,:,:) = rot90(Value_grid_corner,1);
Pix(3,:,:) = rot90(Value_grid_corner,2);
Pix(1,:,:) = rot90(Value_grid_corner,3);



% Create a 3x3 grid of subplots with surf plots
figure(10);

% Loop through the 3x3 grid
for i = 1:3
    for j = 1:3
        % Calculate subplot index
        subplotIndex = (i - 1) * 3 + j;
    switch subplotIndex
        case 1
            P = 7;
        case 2
            P = 8;
        case 3
            P = 9;
        case 4
            P = 4;
        case 5
            P = 5;
        case 6
            P = 6;
        case 7
            P = 1;
        case 8
            P = 2;
        case 9
            P = 3;
    end

        % Create subplot
        subplot(3, 3, subplotIndex);

        % Plot the surface
        mesh(Y_unique, Z_unique, squeeze(Pix(P,:,:)));
        view(0,90)
        
        % Add labels and title
        xlabel('Z');
        ylabel('Y');
        zlabel('CIE');
        title(sprintf('Pix (%d)', P));
    end
end

%%
% Approach1: only apply all pix correction outside of anode
Approach1 = squeeze(sum(max(Pix(:,:,:),0)));
Approach1(5:end-4,5:end-4) = squeeze(Pix(5,5:end-4,5:end-4));


%%
figure(11)
Y_increments = [1:10]
Y_increments(11:20) = [11:10:101]
Y_increments(21:31) = [111:121]
for row = Y_increments
plot(Z_unique(1:10),squeeze(Pix(5,row,1:10)),'-k','DisplayName','Centre')
hold on

plot(Z_unique(1:10),squeeze(Pix(1,row,1:10)),'.-b','DisplayName','CornerSW')
plot(Z_unique(1:10),squeeze(Pix(7,row,1:10)),'b','DisplayName','CornerNW')
plot(Z_unique(1:10),squeeze(Pix(9,row,1:10)),'--b','DisplayName','CornerNE')
plot(Z_unique(1:10),squeeze(Pix(3,row,1:10)),':b','DisplayName','CornerSE')

plot(Z_unique(1:10),squeeze(Pix(2,row,1:10)),'-r','DisplayName','EdgeS')
plot(Z_unique(1:10),squeeze(Pix(4,row,1:10)),'-.r','DisplayName','EdgeW')
plot(Z_unique(1:10),squeeze(Pix(6,row,1:10)),'--r','DisplayName','EdgeE')
plot(Z_unique(1:10),squeeze(Pix(8,row,1:10)),':r','DisplayName','EdgeN')

AllPix = sum(max(Pix(:,row,1:10),0),1);
plot(Z_unique(1:10),squeeze(AllPix),'-g','DisplayName','AllPix')

hold off
legend('Location','south','NumColumns',4)
%title(row)
ylim([-0.4 1])

pause(0.1)
end
%%
figure(12)
for col = 1:Z_increments
plot(Y_unique(1:10),squeeze(Pix(5,1:10,col)),'-k','DisplayName','Centre')
hold on

plot(Y_unique(1:10),squeeze(Pix(1,1:10,col)),'.-b','DisplayName','CornerSW')
plot(Y_unique(1:10),squeeze(Pix(7,1:10,col)),'b','DisplayName','CornerNW')
plot(Y_unique(1:10),squeeze(Pix(9,1:10,col)),'--b','DisplayName','CornerNE')
plot(Y_unique(1:10),squeeze(Pix(3,1:10,col)),':b','DisplayName','CornerSE')

plot(Y_unique(1:10),squeeze(Pix(2,1:10,col)),'-r','DisplayName','EdgeS')
plot(Y_unique(1:10),squeeze(Pix(4,1:10,col)),'-.r','DisplayName','EdgeW')
plot(Y_unique(1:10),squeeze(Pix(6,1:10,col)),'--r','DisplayName','EdgeE')
plot(Y_unique(1:10),squeeze(Pix(8,1:10,col)),':r','DisplayName','EdgeN')

AllPix = sum(max(Pix(:,1:10,col),0),1);
plot(Y_unique(1:10),squeeze(AllPix),'-g','DisplayName','AllPix')

plot(Y_unique(1:10),Approach1(1:10,col),'--r','DisplayName','Approach1')

hold off
legend('Location','southeast')
title(col)
ylim([-0.25 1.25])
pause(1)
end



%% Older running (but working incorrectly) code
%{
figure(5)
maxSignal = Value_grid_edge(round(Y_increments/2), round(Z_increments/2));
for ii = 1:Y_increments
CentreValues = Value_grid_centre(1:10,ii);
EdgeValues = Value_grid_edge(end:-1:end-9,ii);
CornerValues = Value_grid_corner(end:-1:end-9,ii);
CplusE = CentreValues + EdgeValues;
AllPix = CplusE + CornerValues;

plot(Z_unique(1:10),CentreValues(:),'-ob')
hold on
plot(Z_unique(1:10),EdgeValues(:),'-or')
plot(Z_unique(1:10),CornerValues(:),'-og')
plot(Z_unique(1:10),CplusE(:),'-oc')
plot(Z_unique(1:10),AllPix(:),'-ok')
legend('Centre', 'Edge','Corner', 'CplusE', 'AllPix')
title(ii)
hold off
waitforbuttonpress()
%pause(1)
end


%%
figure(6)
maxSignal = Value_grid_edge(round(Y_increments/2), round(Z_increments/2));
for ii = 1:Y_increments
CentreValues = Value_grid_centre(1:10,ii);
EdgeValues = Value_grid_edge(end:-1:end-9,ii);
CornerValues = Value_grid_corner(end:-1:end-9,ii);
if ii < 5
CplusE = CentreValues + 2*EdgeValues;
else
CplusE = CentreValues + EdgeValues;
end
AllPix = CplusE + CornerValues;

plot(Z_unique(1:10),CentreValues(:),'-ob')
hold on
plot(Z_unique(1:10),EdgeValues(:),'-or')
plot(Z_unique(1:10),CornerValues(:),'-og')
plot(Z_unique(1:10),CplusE(:),'-oc')
plot(Z_unique(1:10),AllPix(:),'-ok')
legend('Centre', 'Edge','Corner', 'CplusE', 'AllPix')
title(ii)
hold off
waitforbuttonpress()
%pause(1)
end

%%

clc
maxSignal = Value_grid_edge(round(Y_increments/2), round(Z_increments/2));
step = 8;

CSEApproach_OutsideAnodeFullCoGI = Value_grid_centre;

for yy = 1:Y_increments
    for zz = 1:Z_increments
if yy < step+1 || pos1 > Y_increments-step
if zz < Z_increments/2
    CSEApproach_OutsideAnodeFullCoGI(zz,yy) = CSEApproach_OutsideAnodeFullCoGI(zz,yy) + Value_grid_edge(yy,end+1-zz);
else
    CSEApproach_OutsideAnodeFullCoGI(zz,yy) = CSEApproach_OutsideAnodeFullCoGI(zz,yy) + Value_grid_edge(yy,zz);
end
end

if zz < step+1 || zz > Z_increments-step 
    if zz < Z_increments/2
    CSEApproach_OutsideAnodeFullCoGI(zz,yy) = CSEApproach_OutsideAnodeFullCoGI(zz,yy) + Value_grid_edge(yy, end+1-zz);
    else
        CSEApproach_OutsideAnodeFullCoGI(zz,yy) = CSEApproach_OutsideAnodeFullCoGI(zz,yy) + Value_grid_edge(yy,zz);
    end
end




if (yy < step+1 || yy > Y_increments-step) && (zz < step+1 || zz > Z_increments-step)
    if zz < Z_increments/2
    CSEApproach_OutsideAnodeFullCoGI(zz,yy) = CSEApproach_OutsideAnodeFullCoGI(zz,yy) + Value_grid_corner(end+1-zz,yy);
    else
    CSEApproach_OutsideAnodeFullCoGI(zz,yy) = CSEApproach_OutsideAnodeFullCoGI(zz,yy) + Value_grid_corner(zz,yy);
    end    
end

    end
end






for yy = 1:Y_increments
    for zz = 1:Z_increments

        CSEApproach_AllThresholded(zz,yy) = Value_grid_centre(zz,yy);
        %code for first half
        if zz < Z_increments/2
            if Value_grid_edge(end+1-zz,yy) > 0
                numberOfEdges = 0;
                if yy < 5 || yy > Y_increments-4
                    numberOfEdges = numberOfEdges +1;
                end
                if zz < 5 || zz > Z_increments-4
                    numberOfEdges = numberOfEdges +1;
                end

                CSEApproach_AllThresholded(zz,yy) = CSEApproach_AllThresholded(zz,yy) + Value_grid_edge(end+1-zz,yy)*numberOfEdges;
            end

            if Value_grid_corner(end+1-zz,yy) > 0
                CSEApproach_AllThresholded(zz,yy) = CSEApproach_AllThresholded(zz,yy) + Value_grid_corner(end+1-zz,yy);
            end

            %code for second half (firt reflected)
        else

            if Value_grid_edge(zz,yy) > 0
                numberOfEdges = 0;
                if yy < step+1 || yy > Y_increments-step
                    numberOfEdges = numberOfEdges +1;
                end
                if zz < step+1 || zz > Z_increments-step
                    numberOfEdges = numberOfEdges +1;
                end

                CSEApproach_AllThresholded(zz,yy) = CSEApproach_AllThresholded(zz,yy) + Value_grid_edge(zz,yy)*numberOfEdges;
            end

            if Value_grid_corner(zz,yy) > 0
                CSEApproach_AllThresholded(zz,yy) = CSEApproach_AllThresholded(zz,yy) + Value_grid_corner(zz,yy);
            end

        end

    end
end

%{
for yy = 1:Y_increments
    for zz = 1:Z_increments
if yy < step+1 || 

        CSEApproach_AllThresholded(zz,yy) = Value_grid_centre(zz,yy);

    end
end
%}


figure(8)
surf(Z_grid_centre,Y_grid_centre,CSEApproach_AllThresholded)

figure(9)
surf(Z_grid_centre,Y_grid_centre,CSEApproach_OutsideAnodeFullCoGI)

figure(7)

for ii = 1:Y_increments
CentreValues = Value_grid_centre(ii,1:10);
EdgeValues = Value_grid_edge(end:-1:end-9,ii);
CornerValues = Value_grid_corner(end:-1:end-9,ii);
CplusE = CentreValues + EdgeValues;
AllPix = CplusE + CornerValues;
if yy < 5 || yy > Y_increments - 4
AllPix_OutsideAnode = CplusE + CornerValues*2;
else
AllPix_OutsideAnode = CplusE + CornerValues;
end

CSEAValues_AT = CSEApproach_AllThresholded(1:10,ii);
CSEAValues_OAFC = CSEApproach_OutsideAnodeFullCoGI(1:10,ii);

%{
START WITH BASE
ADD EDGES IF LESS THAN 3 IN X OR Y
ADD CORNER IF LESS THAN 3 IN X AND Y
%}

plot(Z_unique(1:10),CentreValues(:),'-ob','DisplayName','Centre')
hold on
plot(Z_unique(1:10),EdgeValues(:),'-or','DisplayName','Edge')
plot(Z_unique(1:10),CornerValues(:),'-og','DisplayName','Corner')
%plot(Z_unique(1:10),CplusE(:),'-om','DisplayName','CplusE')
%plot(Z_unique(1:10),AllPix(:),':oc','DisplayName','AllPixRaw')
%plot(Z_unique(1:10),AllPix_OutsideAnode(:),':og','DisplayName','AllPix_OA')
plot(Z_unique(1:10),CSEAValues_AT(:),'-.ok','DisplayName','CoGI_Thresh')
plot(Z_unique(1:10),CSEAValues_OAFC(:),'-.or','DisplayName','CoGI_FullOA')

legend show
title(ii)
hold off
if ii < 10 || ii > Y_increments - 9
    waitforbuttonpress()
end
pause(0.025)
end




%{

CSEApproach(zz,yy) = Value_grid_centre(zz,yy);
if yy < step+1
    CSEApproach(zz,yy) = CSEApproach(zz,yy) + Value_grid_edge(end+1-zz,yy);
end


if yy > Y_increments-step
    CSEApproach(zz,yy) = CSEApproach(zz,yy) + Value_grid_edge(end+1-zz,yy);
end

if zz < step+1 
    CSEApproach(zz,yy) = CSEApproach(zz,yy) + Value_grid_edge(end+1-zz,yy);
end

if zz > Z_increments-step 
    CSEApproach(zz,yy) = CSEApproach(zz,yy) + Value_grid_edge(zz,yy);
end



if (yy < 5 || yy > Y_increments-4) && (zz < 5 || zz > Z_increments-4)
    CSEApproach(zz,yy) = CSEApproach(zz,yy) + Value_grid_corner(end+1-zz,yy);
end

%}



%{
CSEApproach(zz,yy) = Value_grid_centre(zz,yy);
if yy < step+1 || yy > Y_increments-step
    CSEApproach(zz,yy) = CSEApproach(zz,yy) + Value_grid_edge(end+1-zz,yy);
end

if zz < step+1 || zz > Z_increments-step
    CSEApproach(zz,yy) = CSEApproach(zz,yy) + Value_grid_edge(end+1-zz,yy);
end

if (yy < step+1 || yy > Y_increments-step) && (zz < step+1 || zz > Z_increments-step)
    CSEApproach(zz,yy) = CSEApproach(zz,yy) + Value_grid_corner(end+1-zz,yy);
end
%}


%}

%% TEST AREA
figure(10)
surf(Value_grid_corner(:,1:round(end/10)))


%%



function [result, Y_grid, Z_grid, Value_grid, Y_unique, Z_unique, Y_increments, Z_increments, Points_assessed] = dataExtractor(RawData,desiredX,desiredT, ImNum)

%
% Extract the X, Y, Z coordinates
X = RawData(:, 1);
Y = RawData(:, 2);
Z = RawData(:, 3);

% Find the range of X, Y, and Z
xRange = [min(X), max(X)];
yRange = [min(Y), max(Y)];
zRange = [min(Z), max(Z)];

fprintf('Range of X: [%f, %f]\n', xRange(1), xRange(2));
fprintf('Range of Y: [%f, %f]\n', yRange(1), yRange(2));
fprintf('Range of Z: [%f, %f]\n', zRange(1), zRange(2));

% Calculate the column index for the specified time t (4th column onward)
valueColumn = 4 + desiredT;

% Extract the slice where X equals the desired value
slice = RawData(abs(X - desiredX) < 1e-6, :); % Tolerance for floating-point comparison

% Extract the Y, Z, and value data for the given time step
Y_slice = slice(:, 2);
Z_slice = slice(:, 3);
values = slice(:, valueColumn);

result = [Y_slice, Z_slice, values];
%disp(result);



clc
Y = result(:, 1); % First column for Y
Z = result(:, 2); % Second column for Z
Value = result(:, 3); % Third column for the values

Y_unique = unique(Y);
Z_unique = unique(Z);
[Y_grid, Z_grid] = meshgrid(Z_unique, Y_unique);

Value_grid = reshape(Value, length(Y_unique), length(Z_unique));

figure(ImNum)
surf(Y_grid,Z_grid,Value_grid)
view(0,90)

Y_increments = size(Y_unique,1);
Z_increments = size(Z_unique,1);
Points_assessed = size(result,1);

end