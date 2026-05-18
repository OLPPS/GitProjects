%%
clear
clc

xpos = 2;
ypos = 4;
zpos = 3;
energy = 5;
time = 1;

UniqueID = 'WO_Sm_RC_Meta_10e6_0_10e6';
inputfilename = strcat(UniqueID, '.hits.txt');
fid = fopen(inputfilename, 'r');
fgetl(fid);  % Skip the first line
data = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter', ',', 'CollectOutput', true);  % Adjust format as needed
fclose(fid);

EventsData = zeros(size(data{1},1),1);

EventsData(:,1) = data{1}(:,xpos);
EventsData(:,2) = data{1}(:,ypos);
EventsData(:,3) = data{1}(:,zpos);
EventsData(:,4) = data{1}(:,time);
EventsData(:,5) = data{1}(:,energy);

outputfilename = strcat(UniqueID, '__FormattedCorrectlyAndCompressed');
save(outputfilename,'-v7.2');

%%

range(data{1}(:,4))
