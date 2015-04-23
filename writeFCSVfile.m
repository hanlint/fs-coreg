function writeFCSVfile(fcsvPath,coords,name)

fid = fopen(fcsvPath,'w');

% print header information
fprintf(fid,'# Markups fiducial file version = 4.3\n');
fprintf(fid,'# CoordinateSystem = 0\n');
fprintf(fid,'# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');

% print for each label
for i = 1:length(name)
    fprintf(fid,'FiducialNode_%d,%.2f,%.2f,%.2f,0,0,0,1,1,1,1,,%s,,\n',i,coords(i,1),coords(i,2),coords(i,3),name{i});
end

fclose(fid);

fprintf('Writing FCSV file to: %s\n',fcsvPath);

end