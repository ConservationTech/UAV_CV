%% Sam Kelly
% Duke University Marine Lab
% Unmanned Autonomous Systems Lab

% Description: Using UAS technology to find adapted QR Codes from a 
%              collection of drone (or any pictures). The underlying
%              purpose is of this code is for use in wildlife conservation
%              (particularly sea turtles)

% clear all
warning('off','all')
%% Parameters
    
    HomeFolder  = 'D:\OneDrive - Duke University\Research\QR Drone Tags\Images\UAV\DJI_Okareka_10m';
    SaveLoc     = 'D:\OneDrive - Duke University\Research\QR Drone Tags\Images\Analyzed Images';
    PhotoLog    = sprintf('%s\\photolog.txt', HomeFolder);
    TrialName   = '3_DJI_Okareka_10m';
    addpath(genpath(HomeFolder))
    addpath(genpath(SaveLoc))
    addpath(genpath('D:\OneDrive - Duke University\Research\QR Drone Tags\MATLAB\v1.4'))
    
    StartPhoto  = 1;              % Photo analysis starts at...
    Aircraft    = 1;              % 0 = Photo Log 1 = Geotagged Photos 2 = No Geo Tag Available 
    AllKML      = 1;              % 1 = All KML Files, 2 = Unique KML Files 
    TagMode     = 1;              % 1 = BOW Tags, 2 = WOB Tags
    Bradley     = 1;              % 1 = Extensive, 0 = Simple
    TagSize     = 950;            % Minimum number of pixels of tags
    ColorMode   = 1;              % 0 = BW, 1 = Gray, 2 = Color
    Darkness    = .3;            % 0-3 Lower Makes the Image Darker
    Hemisphere  = 0;              % 1 = Northern Hemisphere, 0 = Southern
    % codelist    = [1 9 33 14 19];
    
%% Set Directories
    PictureDir  = dir(sprintf('%s\\*.JPG', HomeFolder));
    PictureDirC = struct2cell(PictureDir);              % Convert to Cells
    mkdir(SaveLoc,TrialName)
    mkdir(sprintf('%s\\%s',SaveLoc,TrialName), 'KML Tags')
    mkdir(sprintf('%s\\%s',SaveLoc,TrialName), 'Pictures')
    mkdir(sprintf('%s\\%s',SaveLoc,TrialName), 'Tag Log')
    addpath(genpath(SaveLoc))


%% Parallel Loop to identify QR Codes
parfor j = StartPhoto:length(PictureDir)
    
    % Load Pictures
     PicLoad = imread([PictureDir(j).name]);    % load pictures
     PicInfo = imfinfo([PictureDir(j).name]);   % Load Picture Info
    
     
    % Locate QR Code
    if TagMode == 2
        TagInfo = locateCodes(PicLoad, 'colMode', ColorMode, 'threshMode', Bradley, ...
            'bradleyThreshold', Darkness,'sizeThresh', TagSize);
     
    else 
        TagInfo = locateCodes(PicLoad, 'colMode', ColorMode, 'threshMode', Bradley, ...
            'bradleyThreshold', Darkness,'sizeThresh', TagSize);
        
    end
    
    print(sprintf('%s\\%s\\Pictures\\QR-%s',SaveLoc,TrialName, ...
            PictureDir(j).name),'-djpeg');
         
        TagInfoC = struct2cell(TagInfo);  % Convert Output to Cells
        K        = length(TagInfo);       % Necessary for loops
    
    if Aircraft == 2
        if isempty(TagInfoC)  % end loop if empty (no pics/KML tags)

            fprintf('%s : No tags found \n',PictureDir(j).name)
        
        else                  % continue with outputs otherwise
        % Create Log
            Time    = cell(1, K);
            Time(:) = {PictureDir(j).date};

            PhotoNum    = cell(1, K);
            PhotoNum(:) = {PictureDir(j).name};

            TagEntry  = TagInfoC(8,:)';
            SizeCheck = TagInfoC(1,:)';
            Orientation = TagInfoC(7,:)';

            fprintf('%s : %d Tags found \n',PictureDir(j).name,K)  % Tag Data

             %% Create Table Log
             FinalTable = cell2table([TagEntry, PhotoNum', Time', Orientation, SizeCheck]); %Orientation, 
             FinalTable.Properties.VariableNames{'Var1'} = 'Tag_Number';
             FinalTable.Properties.VariableNames{'Var2'} = 'Photo';
             FinalTable.Properties.VariableNames{'Var3'} = 'Time';
             FinalTable.Properties.VariableNames{'Var4'} = 'Orientation';
             FinalTable.Properties.VariableNames{'Var5'} = 'Size';
             % Save Table
             filename = sprintf('%s\\%s\\Tag Log\\%d.csv',SaveLoc,TrialName,j);
             writetable(FinalTable, filename, 'Delimiter', ',');
        end

    else
        if isempty(TagInfoC)  % end loop if empty (no pics/KML tags)

            fprintf('%s : No tags found \n',PictureDir(j).name)
        
        else                  % continue with outputs otherwise
            % Load Picture Data

            if Aircraft == 1   % GeoTagged Images
                GPSData = PicInfo.GPSInfo;
                if Hemisphere == 1
                    Lat     = dms2degrees(GPSData.GPSLatitude);
                    Long    = -dms2degrees(GPSData.GPSLongitude);
                else
                    Lat     = -dms2degrees(GPSData.GPSLatitude);
                    Long    = dms2degrees(GPSData.GPSLongitude);
                end
                
                Height  = GPSData.GPSAltitude;
                % Convert Numbers to Cells
                LatC    = num2cell(ones(K,1).*Lat);
                LongC   = num2cell(ones(K,1).*Long);
                HeightC = num2cell(ones(K,1).*Height);

            else   % Photo Log Images
                GPSData = importdata(PhotoLog,'\t',2);
                Lat     = GPSData.data(j,1);
                Long    = GPSData.data(j,2);
                Height  = GPSData.data(j,3);
                % Convert Numbers to Cells
                LatC    = num2cell((ones(K,1).*GPSData.data(j,1)));
                LongC   = num2cell((ones(K,1).*GPSData.data(j,2)));
                HeightC = num2cell((ones(K,1).*GPSData.data(j,3)));
            end
            
            Time    = cell(1, K);
            Time(:) = {PictureDir(j).date};

            PhotoNum    = cell(1, K);
            PhotoNum(:) = {PictureDir(j).name};

            TagEntry  = TagInfoC(8,:)';
            SizeCheck = TagInfoC(1,:)';
            Orientation = TagInfoC(7,:)';
            
            fprintf('%s : %d Tags found (Height: %.1f meters) \n',PictureDir(j).name,K,Height)  % Tag Data

            %% KML Tag Entry
            for z = 1:K
                description = sprintf('Image: %s /n Time: %d',...
                                      PictureDir(j).name,cell2mat(PictureDirC(2,j)'));
                TagName = sprintf('%d',cell2mat(TagEntry(z)));
    %             if AllKML == 1
                    KMLName = sprintf('%s\\%s\\KML Tags\\%s(%d-%d)',...
                                      SaveLoc,TrialName,TagName,j,z);  
    %             else
    %                 KMLName = sprintf('%s\\%s\\KML Tags\\%s',...
    %                                   SaveLoc,TrialName,TagName); 
    %             end
                kmlwritepoint(KMLName,Lat,Long,'Name', TagName, ...
                              'Description',description,'IconScale', 2);
            end
             %% Create Table Log
             FinalTable = cell2table([TagEntry, LatC, LongC, PhotoNum', Time', HeightC, Orientation, SizeCheck]); 
             FinalTable.Properties.VariableNames{'Var1'} = 'Tag_Number';
             FinalTable.Properties.VariableNames{'Var2'} = 'Latitude';
             FinalTable.Properties.VariableNames{'Var3'} = 'Longitude';
             FinalTable.Properties.VariableNames{'Var4'} = 'Photo';
             FinalTable.Properties.VariableNames{'Var5'} = 'Time';
             FinalTable.Properties.VariableNames{'Var6'} = 'Height';
             FinalTable.Properties.VariableNames{'Var7'} = 'Orientation';
             FinalTable.Properties.VariableNames{'Var8'} = 'Size';
             % Save Table
             filename = sprintf('%s\\%s\\Tag Log\\%d.csv',SaveLoc,TrialName,j);
             writetable(FinalTable, filename, 'Delimiter', ',');
        end   % end of tag, no-tag if statement 
    end
end   % end of picture loop


%% Table Analysis
    addpath(genpath(HomeFolder))
    LogDir  = dir(sprintf('%s\\%s\\Tag Log\\*.csv',SaveLoc,TrialName)); 
if Aircraft == 2
        if size(LogDir,1) ~= 0
         % Get list of files
             LogM = csv2cell(LogDir(1).name, 'fromfile');  % First file
             for n = 2:numel(LogDir)
               new = csv2cell(LogDir(n).name,'fromfile');  % Read the nth file
               new(1,:) = [];
               LogM = [LogM;new];  
             end
        % Overview
            fprintf('**ANALYSIS COMPLETE**\n')
            LogA = LogM;
            LogA(1,:) = [];
            Total_Tags = size(LogA,1)
            % Unique ID
            [~, ind] = unique(LogA(:, 1), 'rows');
            CleanLog = LogA(ind,:);
            UniqueTagsFound = CleanLog(:,1)
            UTags_Q  = size(CleanLog,1);
            fprintf('%d Unique Tag(s) Identified \n%d Total Tag(s) Found\n \n', UTags_Q, Total_Tags)
        % Save
             % All tags
             Log = cell2table(LogM);
             Log(1,:) = [];
             Log.Properties.VariableNames{'LogM1'} = 'TagNumber';
             Log.Properties.VariableNames{'LogM2'} = 'Photo';
             Log.Properties.VariableNames{'LogM3'} = 'Time';
             Log.Properties.VariableNames{'LogM4'} = 'Orientation';
             Log.Properties.VariableNames{'LogM5'} = 'Size';
             filename = sprintf('%s\\%s\\%s - All Tags.csv',SaveLoc,TrialName,TrialName);
             writetable(Log, filename, 'Delimiter', ',');
             % Unique Tags
             U  = [LogM(1,:);CleanLog];
             ULog = cell2table(U);
             ULog(1,:) = [];
             ULog.Properties.VariableNames{'U1'} = 'TagNumber';
             ULog.Properties.VariableNames{'U2'} = 'Photo';
             ULog.Properties.VariableNames{'U3'} = 'Time';
             ULog.Properties.VariableNames{'U4'} = 'Orientation';
             ULog.Properties.VariableNames{'U5'} = 'Size';
             filename = sprintf('%s\\%s\\%s - Unique Tags.csv',SaveLoc,TrialName,TrialName);
             writetable(ULog, filename, 'Delimiter', ',');
        else
            fprintf('\n \n **ANALYSIS COMPLETE** \n No tags found during mission\n')
        end
    
else
    if size(LogDir,1) ~= 0
         % Get list of files
         LogM = csv2cell(LogDir(1).name, 'fromfile');  % First file
         for n = 2:numel(LogDir)
           new = csv2cell(LogDir(n).name,'fromfile');  % Read the nth file
           new(1,:) = [];
           LogM = [LogM;new];  
         end
    % Overview
        fprintf('**ANALYSIS COMPLETE**\n')
        LogA = LogM;
        LogA(1,:) = [];
        Total_Tags = size(LogA,1)
        % Unique ID
        [~, ind] = unique(LogA(:, 1), 'rows');
        CleanLog = LogA(ind,:);
        UniqueTagsFound = CleanLog(:,1)
        UTags_Q  = size(CleanLog,1);
        fprintf('%d Unique Tag(s) Identified \n%d Total Tag(s) Found\n \n', UTags_Q, Total_Tags)
    % Save
         % All tags
         Log = cell2table(LogM);
         Log(1,:) = [];
         Log.Properties.VariableNames{'LogM1'} = 'TagNumber';
         Log.Properties.VariableNames{'LogM2'} = 'Latitude';
         Log.Properties.VariableNames{'LogM3'} = 'Longitude';
         Log.Properties.VariableNames{'LogM4'} = 'Photo';
         Log.Properties.VariableNames{'LogM5'} = 'Time';
         Log.Properties.VariableNames{'LogM6'} = 'Height';
         Log.Properties.VariableNames{'LogM7'} = 'Orientation';
         Log.Properties.VariableNames{'LogM8'} = 'Size';
         filename = sprintf('%s\\%s\\%s - All Tags.csv',SaveLoc,TrialName,TrialName);
         writetable(Log, filename, 'Delimiter', ',');
         % Unique Tags
         U  = [LogM(1,:);CleanLog];
         ULog = cell2table(U);
         ULog(1,:) = [];
         ULog.Properties.VariableNames{'U1'} = 'TagNumber';
         ULog.Properties.VariableNames{'U2'} = 'Latitude';
         ULog.Properties.VariableNames{'U3'} = 'Longitude';
         ULog.Properties.VariableNames{'U4'} = 'Photo';
         ULog.Properties.VariableNames{'U5'} = 'Time';
         ULog.Properties.VariableNames{'U6'} = 'Height';
         ULog.Properties.VariableNames{'U7'} = 'Orientation';
         ULog.Properties.VariableNames{'U8'} = 'Size';
         filename = sprintf('%s\\%s\\%s - Unique Tags.csv',SaveLoc,TrialName,TrialName);
         writetable(ULog, filename, 'Delimiter', ',');
    else
        fprintf('\n \n**ANALYSIS COMPLETE** \n No tags found during mission\n')
    end
end

%end
