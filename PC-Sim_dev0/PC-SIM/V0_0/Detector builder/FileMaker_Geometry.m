%%Uncomment the below code and edit it if you want to generate a geometry
%%file from within this script

            PixelWidth = 0.170; %mm
            PixelPitch = 0.2; %mm
            PixelDepth = 1.5; %mm
            xOffset = 0; %mm x translation of detector from origin in x
            yOffset = 0; %mm y translation of detector from origin in y
            zOffset = 50; %mm % z translation of detector from origin in z
            MagXEdgeSpace = 72.5/2.0; %mm HALF OF Detector length in x (Detector is pixels area used, assuming symmetric about centre of detector)
            MagYEdgeSpace = 0.6/2.0; %mm %HALF OF Detector height in y (Detector is pixels area used, assuming symmetric about centre of detector)
            PixelsInXDirection = floor(72.5/PixelPitch); %max pixels in x direction
            PixelsInYDirection = 3; %max pixels in y direction

GenerateGeometryFile(PixelWidth, PixelPitch, PixelDepth, xOffset, yOffset, zOffset, MagXEdgeSpace, MagYEdgeSpace, PixelsInXDirection, PixelsInYDirection,'rev')

%}


function GenerateGeometryFile(PixelWidth, PixelPitch, PixelDepth, xOffset, yOffset, zOffset, MagXEdgeSpace, MagYEdgeSpace, PixelsInXDirection, PixelsInYDirection,SourceDetectorOrientation)
            GEOMETRIES = zeros(10,1);
            GEOMETRIES(1) = PixelWidth;
            GEOMETRIES(2) = PixelPitch;
            GEOMETRIES(3) = PixelDepth;
            GEOMETRIES(4) = xOffset;
            GEOMETRIES(5) = yOffset;
            GEOMETRIES(6) = zOffset;
            GEOMETRIES(7) = MagXEdgeSpace; 
            GEOMETRIES(8) = MagYEdgeSpace;
            GEOMETRIES(9) = PixelsInXDirection; 
            GEOMETRIES(10) = PixelsInYDirection;

            outputfilename = strcat("GEOMETRIESFile_W",num2str(GEOMETRIES(1)),'_P',num2str(GEOMETRIES(2)),'_D',num2str(GEOMETRIES(3)),'_Xo',num2str(GEOMETRIES(4)),'_Yo',num2str(GEOMETRIES(5)),'_Zo',num2str(GEOMETRIES(6)),'_MX',num2str(GEOMETRIES(7)),'_MY',num2str(GEOMETRIES(8)),'_Xpix',num2str(GEOMETRIES(9)),'_Ypix',num2str(GEOMETRIES(10)),'_',SourceDetectorOrientation)
            outputfilename = strrep(outputfilename, '.', 'pt')
            save(outputfilename,"GEOMETRIES","SourceDetectorOrientation");
end