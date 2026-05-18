clc
for dir = 1:3
    switch dir
        case 1
            cd 'E:\CIE Maps From COMSOL\Central';
        case 2
            cd 'E:\CIE Maps From COMSOL\Edge';
        case 3
            cd 'E:\CIE Maps From COMSOL\Corner';
    end
    for ii = 100:50:600
        for tt = 1:5
            clearvars -except ii tt dir
            switch tt
                case 1
                    thick = '1'
                case 2
                    thick = '1pt5'
                case 3
                    thick = '2'
                case 4
                    thick = '2pt5'
                case 5
                    thick = '3'
            end
            filename = strcat('thickness_',thick,'mm_',num2str(ii),'MicronPitch_Central')
            filename2 = strcat(filename,'.txt');
            %%%check if exists!!!
            fid = fopen(filename2, 'r');
            C = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' ,'HeaderLines',9);
            tline = fgets(fid);
            fclose(fid);
            MAPDATA = zeros(size(C{1},1),39);
            for jj = 1:39
            MAPDATA(:,jj) = C{jj};
            end

            [SelectedSignalStructured3DMap, SelectedSignalTimeIndices3DMap,Structured4DMap_AllTimes,xVals, yVals, zVals] = GenerateCIEMap(MAPDATA,'Max');

            outname = strcat(filename,'_Max_saved')
            save(outname,'SelectedSignalStructured3DMap','SelectedSignalTimeIndices3DMap','xVals', 'yVals', 'zVals','-v7.3')
        end
    end
end