function writeFreesurferLabelFile(filePath,vertex,coords,val)
header = sprintf('#!ascii label  , from subject  vox2ras=TkReg');

if(nargin == 2);
    coords = vertex.coords;
    val = vertex.channel;
    vertex = vertex.vertex;
end


[fid, message] = fopen(filePath,'w');

if(fid == -1)
    error(message);
end

numLabels = length(vertex);

fprintf(fid,'%s\n',header);
fprintf(fid,'%d\n',numLabels);

for k = 1:length(vertex)
    fprintf(fid,'%d  %.3f  %.3f  %.3f %.10f\n',vertex(k), coords(k,1), coords(k,2), coords(k,3),val(k));
end

fclose(fid);


end