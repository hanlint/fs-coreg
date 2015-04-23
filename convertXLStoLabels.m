% function convertXLStoLabels(subject,fileName,debug)
%
% SUBJECT - subject ID
% FILENAME = filename of .xls file within subject directory
% DEBUG - wether to save folders as all_debug to avoid overwriting data.
%         this parameter is optional.
%
%
% e.g. convertXLStoLabels('m00105','m00105.xls');
%
%
%
% this function will take in an XLS file and create labels

% the output should look something like this...
%
% #!ascii label  , from subject m00105 vox2ras=TkReg
% 96
% 101316  -42.091  28.162  24.914 1.0000000000
% 104113  -40.798  38.324  23.707 1.0000000000
% 105604  -32.661  43.547  23.376 1.0000000000
% 107322  -22.234  50.440  22.542 1.0000000000
%
% *** CHANGELOG ***
%
% 09/08/2014 HT: modified to print the RAS coordinates, not the volume
% index into the label file. The parcelation.mat file now also saves
% the labels dataset
%
% 10/09/2014 HT: added ability to print to FCSV files for display in Slicer
% 10/10/2014 HT: function now creats a .py file to load the files into
% Slicer
%
% 10/13/2014 HT: modified the handling of depth electrodes. Depth
% electrodes not require surface coordinates as well. They are in the
% SurfaceRAS space. This allows the electrodes to be plotted in Slicer as
% well.
%
% 10/25/2014 HT: Added handling of depth information. Depths are now
% assigned a region label from the aparc2009s+aseg.mgz file. Electrodes
% that are 'missing' from the CT or just unable to find can also be given a
% 'Missing' attribute.
%
% NOTE: region_codes 0 = unknown
%                    1-75 = surface parcelation
%                    76+ = subcortical parcelation

function labels = convertXLStoLabels(subject,fileName,debug)
addpath(genpath('/Applications/freesurfer/'))

if(nargin <3)
    debug = 0;
end

%subject = 'm00105';
%fileName = 'm00105.xls';
%debug = 1;


basePath = '/KLAB/coregistration';
subjectPath = [basePath '/' subject];
filePath = [basePath '/' subject '/' fileName];

% load subcortical map
x = load(sprintf('%s/scripts/freesurfer_subcortical_map.mat',basePath));
freesurfer_subcortical_map = x.freesurfer_subcortical_map;

x = load(sprintf('%s/scripts/freesurfer_surface2009_map.mat',basePath));
freesurfer_surface2009_map = x.freesurfer_surface2009_map;

warning off;
fprintf('*******************\n');

% confirm that subject folder and file exists
if(~exist(subjectPath,'dir'))
    error(sprintf('Error: subject folder %s does not exist!',subjectPath));
end
if(~exist(filePath,'file'))
    error(sprintf('Error: file %s does not exist!', filePath));
end

% is excel file readable?
%if(isempty(xlsfinfo(filePath)))
%    error('File not a valid Excel file! Make sure it is saved as a XLS file in Excel 97-2004 format');
%end

% attempt to read the data
[~,~,data] = xlsread(filePath);
[nrows ncols] = size(data);

% find patient ID
idx = find(strcmp(data,'Patient:'),1,'first');

if(isempty(idx))
    error('Patient name not found, please make sure the Patient: field exists');
end

patient = data{idx+size(data,1)};
if(isnan(patient))
    error('Patient name not found, please make sure that data field is filled out');
else
    fprintf('Found Subject Name: %s\n',patient);
end

% note: remove trailing whitespace
% note: skip empty rows


% find header rows and idx
% note: must correct for trailing whitespaces!!!
[row0 header_name] = find(strcmp(data,'Name'),1,'first');
[row1 header_channel] = find(strcmp(data,'Channel'),1,'first');
[row2 header_hemi] = find(strcmp(data,'Hemisphere'),1,'first');
[row3 header_surfaceIdx] = find(strcmp(data,'Surface Index'),1,'first');
[row4 header_volumeIdx] = find(strcmp(data,'Volume Index'),1,'first');
[row5 header_regionIdx] = find(strcmp(data,'Region'),1,'first');

isRegionOverride = ~isempty(header_regionIdx);

if(isempty(header_channel) || isempty(header_hemi) || isempty(header_surfaceIdx) || isempty(header_name))
    error('One of the following headers is missing: Channel, Hemisphere, or Surface Index');
elseif(row1 ~= row2 || row2 ~= row3 || row1 ~= row0)
    error('Rows of headers do not agree!');
else
    fprintf('Header row: %d\n',row1);
    fprintf('Name column: %d\n',header_name);
    fprintf('Channel column: %d\n',header_channel);
    fprintf('Hemisphere column: %d\n', header_hemi);
    fprintf('Volume Index column: %d\n',header_volumeIdx);
    fprintf('Surface Index column: %d\n', header_surfaceIdx);
end

fprintf('*******************\n');

% process data into dataset file

row_data_start = row1+1; % data starts at the row after the header
numLabels = nrows-row_data_start+1;

fprintf('# of Labels: %d\n',numLabels);

for i = 1:numLabels
    r = (i-1)+row_data_start; % row
    fprintf('Processing row %d\n',r);
    labels.name{i} = data{r,header_name};
    labels.channel(i) = data{r,header_channel};
    labels.hemi{i} = data{r,header_hemi};
    
    % depending on the label type, read the data differently
    labelType = strtrim(data{r,header_surfaceIdx-1});
    
    labels.label_type{i,1} = labelType;

    switch((labelType))
        
        case {'Vertex','Floating','Depth'}
            
            
            % read the volume index
            str = data{r,header_volumeIdx};
            
            if(isnan(str))
                error(sprintf('row %d does not have any volume index entry', r));
            else
                try
                    C = textscan(str, '[%f, %f, %f]');
                    labels.volume_idx(i,:) = double([C{1} C{2} C{3}]);
                catch
                    error(sprintf('row %d, volume index not formatted properly!'));
                end
            end
            
            
            % read the surface index
            str = data{r,header_surfaceIdx};

            try
                if(strcmpi(labelType,'Vertex'))
                    pat = '%d  [%.2f, %.2f, %.2f]';
                    C = textscan(str, pat);
                    labels.vertex(i) = C{1};
                    labels.coords(i,:) = double([C{2} C{3} C{4}]);
                else
                    pat = '[%.2f, %.2f, %.2f]';
                    C = textscan(str, pat);
                    labels.vertex(i) = -1;
                    labels.coords(i,:) = double([C{1} C{2} C{3}]);
                end
            catch
                error(sprintf('row %d, surface index not formatted properly!',r));
                disp(lasterr);
            end
 
        case 'Missing'
            labels.vertex(i) = -1;
            labels.coords(i,:) = [NaN NaN NaN];
            labels.volume_idx(i,:) = [NaN NaN NaN];          
        case 'SurfaceRAS'
            error('row %d says SurfaceRAS, please change to Depth, Floating, or Missing!',r);
        otherwise
            error('row %d, label type must be one of: Vertex, Floating, Depth, or Missing! Capitalization matters.',r);    
    end
    
    if(isRegionOverride && ~isnan(data{r,header_regionIdx}))
        labels.regionOverride(i,1) = data{r,header_regionIdx};  
    else
        labels.regionOverride(i,1) = NaN;
    end
end

labels.name = labels.name';
labels.channel = labels.channel';
labels.hemi = labels.hemi';
labels.vertex = labels.vertex';

labels = dataset(labels);
labels.hemi = nominal(labels.hemi);

% check that hemi were entered correctly
hemigroup = unique(labels.hemi);
for i = 1:length(hemigroup)
    if(hemigroup(i) ~= 'lh' && hemigroup(i) ~= 'rh')
        error('hemisphere field only have lh or rh entries');
    end
end

% check that the channel numbers are unique
if(length(unique(labels.channel)) ~= length(labels.channel))
    error('channel numbers have repeats');
end

% make label files
header = sprintf('#!ascii label  , from subject %s vox2ras=TkReg',subject);
fprintf('*******************\n');
fprintf('Making label files...\n');
fprintf('Header: %s\n',header);
fprintf('Using folders:\n');

% use left/right folders, or just use all folder?
if(length(hemigroup)== 2)
    use_hemi_folders = 1;
elseif(length(hemigroup) == 1)
    use_hemi_folders = 0;
else
    error('hemisphere field only have lh or rh entries');
end

if(use_hemi_folders && debug)
    folders.lh = [subjectPath '/label/left_debug'];
    folders.rh = [subjectPath '/label/right_debug'];
elseif(use_hemi_folders && ~debug)
    folders.lh = [subjectPath '/label/left'];
    folders.rh = [subjectPath '/label/right'];
elseif(~use_hemi_folders && debug)
    folders.all = [subjectPath '/label/all_debug'];
elseif(~use_hemi_folders && ~debug)
    folders.all = [subjectPath '/label/all'];
end

disp(folders)

% TO DO: clear all data from the label folders before doing anything!
delete(sprintf('%s/label/all.label',subjectPath));
delete(sprintf('%s/label/left.label',subjectPath));
delete(sprintf('%s/label/right.label',subjectPath));
delete(sprintf('%s/label/all_surf.label',subjectPath));
delete(sprintf('%s/label/left_surf.label',subjectPath));
delete(sprintf('%s/label/right_surf.label',subjectPath));

% make directories
if(use_hemi_folders)
    mkdir(folders.lh); mkdir(folders.rh);
else
    mkdir(folders.all);
end

% make individual label files!
for i = 1:length(labels)
    label = labels(i,:);
    
    % only print if vertex, depth, or floating
    if(ismember(label.label_type,{'Vertex','Depth','Floating'}))
        
        
        if(use_hemi_folders && label.hemi == 'lh')
            labelPath = sprintf('%s/%d.label',folders.lh,label.channel);
        elseif(use_hemi_folders && label.hemi == 'rh');
            labelPath = sprintf('%s/%d.label',folders.rh,label.channel);
        elseif(~use_hemi_folders)
            labelPath = sprintf('%s/%d.label',folders.all,label.channel);
        end
        
        fprintf('%d.label | %s | %s\n',label.channel,label.name{1},labelPath);
        
        [fid, message] = fopen(labelPath,'w');
        
        if(fid == -1)
            error(message);
        end
        
        fprintf(fid,'%s\n',header);
        fprintf(fid,'1\n');
        fprintf(fid,'%d  %.3f  %.3f  %.3f %.10f\n',label.vertex, label.coords(1), label.coords(2), label.coords(3),label.channel);
        fclose(fid);
    end
end

% make all labels
if(use_hemi_folders)
    fprintf('Writing %s\n',sprintf('%s.label',folders.lh));
    writeFreesurferLabelFile(sprintf('%s.label',folders.lh),labels(~strcmp(labels.label_type,'Missing') & labels.hemi== 'lh',:));
    
    fprintf('Writing %s\n',sprintf('%s.label',folders.rh));
    writeFreesurferLabelFile(sprintf('%s.label',folders.rh),labels(~strcmp(labels.label_type,'Missing') & labels.hemi== 'rh',:));
    
else
    fprintf('Writing %s\n',sprintf('%s.label',folders.all));
    writeFreesurferLabelFile(sprintf('%s.label',folders.all),labels);
end


% make surf labels

if(use_hemi_folders)
    %if(~isempty(labels(strcmp(labels.label_type,'Vertex') & labels.hemi== 'lh',:)))
        fprintf('Writing %s\n',sprintf('%s_surf.label',folders.lh));
        writeFreesurferLabelFile(sprintf('%s_surf.label',folders.lh),labels(strcmp(labels.label_type,'Vertex') & labels.hemi== 'lh',:));
    %end
    
    %if(~isempty(labels(strcmp(labels.label_type,'Vertex') & labels.hemi== 'rh',:)))
        fprintf('Writing %s\n',sprintf('%s_surf.label',folders.rh));
        writeFreesurferLabelFile(sprintf('%s_surf.label',folders.rh),labels(strcmp(labels.label_type,'Vertex') & labels.hemi== 'rh',:));
    %end
    
elseif(use_hemi_folders)
    fprintf('Writing %s\n',sprintf('%s_surf.label',folders.all));
    writeFreesurferLabelFile(sprintf('%s_surf.label',folders.all),labels(strcmp(labels.label_type,'Vertex') ,:));
end


% assign each electrode a parcelation

fprintf('******\n');
% does the annotation file exist?
annotPath_lh = sprintf('%s/label/lh.aparc.a2009s.annot',subjectPath);
annotPath_rh = sprintf('%s/label/rh.aparc.a2009s.annot',subjectPath);
annotAsegPath = sprintf('%s/mri/aparc.a2009s+aseg.mgz',subjectPath);

if(~exist(annotPath_lh) || ~exist(annotPath_rh))
    fprintf('Files not found: %s\n',annotPath_lh);
    fprintf('Files not found: %s\n',annotPath_rh);
    error('annotation files do not exist - did you run autorecon3?');
else
    [lh.vtx, lh.lbl, lh.colortable] = read_annotation(annotPath_lh);
    [rh.vtx, rh.lbl, rh.colortable] = read_annotation(annotPath_rh);
    aseg = MRIread(annotAsegPath);

    labels.region_code = NaN(length(labels),1);
    labels.region_name = cell(length(labels),1);
    
    for i = 1:length(labels)
        if(labels.hemi(i) == 'rh')
            vtx = rh.vtx;
            lbl = rh.lbl;
            ctable = rh.colortable;
        else
            vtx = lh.vtx;
            lbl = lh.lbl;
            ctable = lh.colortable;
        end
        
        idx = find(labels.vertex(i) == vtx);
        
        switch(labels.label_type{i})
            case 'Missing'
            	labels.region_code(i) = 0;
                labels.region_name{i} = 'Unknown';
            case 'Vertex'
                if(isempty(idx))
                    error('Channel %d vertex %d not found!',i,labels.vertex(i));
                else
                    idx = find(lbl(idx) == ctable.table(:,5));
                    labels.region_code(i) = idx-1;
                    labels.region_name{i} = ctable.struct_names{idx};
                end
                
            case {'Floating','Depth'}
                labels.region_code(i) = getASEGregion([labels.coords(i,:) 1],aseg);
                
                % was mapped to a 2009s parcelation
                if(labels.region_code(i) >= 11100)
                     labels.region_code(i) = mod(labels.region_code(i),100);
                     idx = find(freesurfer_surface2009_map.region_code == labels.region_code(i));
                     if(isempty(idx))
                        error('Surface code: %d not found!',labels.region_code(i));
                    else
                        labels.region_name{i} = freesurfer_surface2009_map.region_name{idx};
                     end  
                else % was mapped to a subcortical segmentation (aseg.mgz)
                    idx = find(freesurfer_subcortical_map.region_code == labels.region_code(i));
                    if(isempty(idx))
                        error('ASEG code: %d not found!',labels.region_code(i));
                    else
                        labels.region_name{i} = freesurfer_subcortical_map.region_name{idx};
                    end
                    if(labels.region_code(i) ~= 0) 
                        labels.region_code(i) = labels.region_code(i) + 75;
                    end
                    
                end
        end
        
        if(~isnan(labels.regionOverride(i)))
            labels.region_code(i) = labels.regionOverride(i);
            
            idx = find(freesurfer_surface2009_map.region_code == labels.regionOverride(i));
            labels.region_name{i} = freesurfer_surface2009_map.region_name{idx};
        end
        

    end
    
end



% Use regexp to get a list of grid names and channel names
for i = 1:length(labels)
    name = labels.name{i};
    tok1 = regexp(name,'(.+)-(\d+)','tokens');
    tok2 = regexp(name,'(\S+) (\d+)','tokens');
    
    isFormat1 = ~(isempty(tok1) || length(tok1{1}) ~= 2);
    isFormat2 = ~(isempty(tok2) || length(tok2{1}) ~= 2);
    
    if(isFormat1)
        tok = tok1;
    elseif(isFormat2)
        tok = tok2;
    else
        error(sprintf('Error: channel name %s not parsable, must be in the right format (e.g. LFT-1 or 4x5G 01)',name));
    end
    
    labels.group_name{i} = tok{1}{1};
    labels.group_channel(i) = str2num(tok{1}{2});
    
end


% make parcelation file
parcelPath = sprintf('%s/label/parcelation.mat',subjectPath);

% sort data
labels = sortrows(labels,'channel');
parcel.channels = labels.channel;
parcel.orig_channels = labels.channel;
parcel.hemisphere = 1+(labels.hemi == 'rh'); % 1 = left hemi, 2 = right hemi
parcel.volumeind = labels.volume_idx;
parcel.vertex = labels.vertex;
parcel.channel_name = labels.name;
parcel.is_depth = labels.vertex == -1;
parcel.labels = labels;
parcel.regions = labels.region_name;
parcel.region_codes = labels.region_code;

fprintf('Writing parcelation file to: %s\n',parcelPath);
save(parcelPath,'-struct','parcel');

% write FCSV file for display in Slicer 3D
fcsvPath = sprintf('%s/label/labels.fcsv',subjectPath);
fid = fopen(fcsvPath,'w');

% print header information
fprintf(fid,'# Markups fiducial file version = 4.3\n');
fprintf(fid,'# CoordinateSystem = 0\n');
fprintf(fid,'# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');

% print for each label
for i = 1:length(labels)
    if(ismember(labels.label_type{i},{'Depth','Vertex','Floating'})) % not a depth or floating or missing
        fprintf(fid,'FiducialNode_%d,%.2f,%.2f,%.2f,0,0,0,1,1,1,1,,%s,,\n',i,labels.coords(i,1),labels.coords(i,2),labels.coords(i,3),labels.name{i});
    end
end

fclose(fid);

fprintf('Writing FCSV file to: %s\n',fcsvPath);

% now, write FCSV files for each grid/strip/depth for display in Slicer 3D
% (with different colors!)
uniqueGroupNames = unique(labels.group_name);

% color cycle: 'blue', 'red',
% 'green','purple','teal','orange','yellow','pink','light green','dark
% teal',
colorCycle = [0 0 255;
    255 0 0;
    50 150 50;
    175 0 255;
    0 255 255;
    255 160 0;
    255 255 0;
    255 0 255;
    0  255 0;
    50 128 128;]/255;


% make label files for each group name and write header information
fid = NaN(length(uniqueGroupNames),1);
for i = 1:length(uniqueGroupNames)
    
    cidx = mod(i,size(colorCycle,1))+1;
    fidPath{i} = sprintf('%s/label/labels_%s.fcsv',subjectPath,uniqueGroupNames{i});
    
    fid(i) = fopen(fidPath{i},'w');
    fprintf(fid(i),'# Markups fiducial file version = 4.3\n');
    fprintf(fid(i),'# CoordinateSystem = 0\n');
    fprintf(fid(i),'# selectedColor = %.2f,%.2f,%.2f\n',colorCycle(cidx,1),colorCycle(cidx,2),colorCycle(cidx,3));
    fprintf(fid(i),'# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');
    
end

% go through labels and print the appropriate file
ngroup = zeros(length(fid),1);
for i = 1:length(labels)
    if(ismember(labels.label_type{i},{'Depth','Vertex','Floating'}))
        group_idx = strcmp(labels.group_name{i},uniqueGroupNames); % get the group it belongs to
        ngroup(group_idx) = ngroup(group_idx)+1;
        
        fprintf(fid(group_idx),'FiducialNode_%d,%.2f,%.2f,%.2f,0,0,0,1,1,1,1,,%s,,\n',i,labels.coords(i,1),labels.coords(i,2),labels.coords(i,3),labels.name{i});
        
    end
    
end

% close all the fids
for i = 1:length(fid)
    fprintf('I wrote to %s, entered for group %s: %d labels\n',fidPath{i},uniqueGroupNames{i},ngroup(i));
    fclose(fid(i));
end

% now make the .py file for slicer
pySlicerFileName = sprintf('%s/load_slicer_%s.py',subjectPath,strrep(fileName,'.xls',''));
fid = fopen(pySlicerFileName,'w');

fprintf('Writing Python script to load slicer modules to %s\n',pySlicerFileName);

% first, adjust the background colors
fprintf(fid,'viewNode = slicer.app.layoutManager().threeDWidget(0).threeDView().mrmlViewNode()\n');
fprintf(fid,'viewNode.SetBackgroundColor(0,0,0)\n');
fprintf(fid,'viewNode.SetBackgroundColor2(0,0,0)\n');
fprintf(fid,'viewNode.SetAxisLabelsVisible(0)\n');
fprintf(fid,'viewNode.SetBoxVisible(0)\n');
fprintf(fid,'\n');

% now position the camera
if(sum(labels.hemi == 'lh') > sum(labels.hemi == 'rh'))
    camPosition = [-601,-4.0,6.0];
else
    camPosition = [599,-4.0,6.0];
end

fprintf(fid,'nodes = slicer.util.getNode("vtkMRMLCameraNode*")\n');
fprintf(fid,'cam = nodes.GetCamera()\n');
fprintf(fid,'cam.SetPosition(%.2f,%.2f,%.2f)\n',camPosition(1),camPosition(2),camPosition(3));

% load the surf files
if(any(labels.hemi == 'lh'))
    fprintf(fid,'slicer.util.loadModel("%s/surf/lh.pial")\n',subjectPath);
end
if(any(labels.hemi == 'rh'))
    fprintf(fid,'slicer.util.loadModel("%s/surf/rh.pial")\n',subjectPath);
end

% load the fiducial files
for i = 1:length(fidPath)
    if(ngroup(i) > 0)
        fprintf(fid,'slicer.util.loadMarkupsFiducialList("%s")\n',fidPath{i});
    end
end

fprintf(fid,'scene = slicer.mrmlScene\n');

% set each fiducial to its appropriate color
for i = 1:length(fidPath)
    if(ngroup(i) > 0)
        cidx = mod(i,size(colorCycle,1))+1;
        
        fprintf(fid,'fiducial = scene.GetFirstNodeByName("labels_%s")\n',uniqueGroupNames{i});
        fprintf(fid,'fiducialDisplay = fiducial.GetMarkupsDisplayNode()\n');
        fprintf(fid,'fiducialDisplay.SetSelectedColor(%.2f,%.2f,%.2f)\n',colorCycle(cidx,1),colorCycle(cidx,2),colorCycle(cidx,3));
        
        % stupid thing to refresh the GUI display .. better way needed
        fprintf(fid,'fiducial.SetNthFiducialVisibility(1,0)\n');
        fprintf(fid,'fiducial.SetNthFiducialVisibility(1,1)\n');
        fprintf(fid,'\n');
    end
end

fclose(fid);


% done!
fprintf('convertXLStoLabels completed without error! congratulations!\n');
end






















