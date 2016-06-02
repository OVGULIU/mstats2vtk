function mstats2vtk(ms, varargin)
    % *mstats2vtk* writes microstructure and 
    % its statistics to a vtk file for visualization in ParaView
    %
    % % Syntax
    % * mstats2vtk(ms) - writes microstructure only to a vtk file in the
    %   current folder with default file name "mstats.vtk"
    % * mstats2vtk(ms,'fname',vtkFilePath) - writes microstructure only to a vtk file
    %   with full path specified in vtkFilePath
    % * mstats2vtk(ms,'tps',tps) - writes microstructure and two-point
    %   statistics 
    % * mstats2vtk(ms,'basis',basis) writes microstructure and PC bases
    % * mstats2vtk(ms,'tps',tps,'basis',basis) - writes microstructure, two-point
    %   statistics, and PC bases
    %
    % % Input
    % * ms - array containing states (e.g. phase IDs) with [xDim,yDim,zDim]
    %   dimensions
    % * tps - array containing spatial correlations with 
    %   - [xDim,yDim,zDim] dimensions for the case of only one correlation to export, or
    %   - [numCorr,xDim,yDim,zDim] dimensions with numCorr being the number of correlations
    %   (e.g. numCorr=2 for 2 autocorrelations, etc.) for the case of multiple correlations
    % * basis - array containing PC bases with [xDim*yDim*zDim,numPCs] dimensions,
    %   where numPCs is the number of principal components
    
    % Author: Marat I. Latypov
    % while at Georgia Tech Lorraine
    % E-mail: marat.latypov@georgiatech-metz.fr
    % Website: latmarat.net
	% June 2016
    
    %% Digest user input
    basis = get_option(varargin,'basis',[],'double');
    tps = get_option(varargin,'tps',[],'double');
    vtkFileName = get_option(varargin,'fname','mstats.vtk','char');
    
    %% Preparations
    
    % open the file
    vtkFile = fopen(vtkFileName,'wt');
        
    % domain size
    dsize = size(ms);
    % dimensions of the domain
    dims = ones(1,3);
    dims(1:numel(dsize)) = dsize;
    
    %% Write header and mesh
    head = ['# vtk DataFile Version 2.0\n'...
            'Data set from mstats2vtk\n' ...
            'ASCII\n\n'];
    
    fprintf(vtkFile,head);
    fprintf(vtkFile,'DATASET RECTILINEAR_GRID\n');
    fprintf(vtkFile,'DIMENSIONS %d %d %d',dims(1)+1,dims(2)+1,dims(3)+1);
    
    % formats for writing 10 numbers per line
    fmtF = [repmat('%.6f ',1,9) '%.6f \n'];
    fmtD = [repmat('%d ',1,9) '%d \n'];
    fmtE = [repmat('%e ',1,9) '%e \n'];
    
    % write coordinates
    dimName = 'XYZ';
    for ii = 1:numel(dims)
        fprintf(vtkFile,'\n%s_COORDINATES %d float\n',dimName(ii), dims(ii)+1);
        nodes = -0.5:dims(ii)-0.5;
        fprintf(vtkFile,fmtF,nodes);
    end
    % write data size
    fprintf(vtkFile,'\nCELL_DATA %d\n',numel(ms));
    
    %% Write the data
    % Write states
    fprintf('\nWriting states...\n')
    writeScalars(vtkFile,'states','int',fmtD,reshape(ms,[numel(ms),1]))
%     fprintf(vtkFile,'SCALARS states int 1\n');
%     fprintf(vtkFile,'LOOKUP_TABLE default\n');
%     fprintf(vtkFile,fmtD,reshape(ms,[numel(ms),1]));
    
    % Write spatial correlations
    % check if stats size matches microstructure domain
    if (numel(tps)/size(tps,1)) == prod(dims) 
        fprintf('Writing two-point statistics...\n');
        for ii = 1:size(tps,1)
            dname = sprintf('f_%03d',ii);
            writeScalars(vtkFile,dname,'float',fmtE,reshape(tps,[numel(tps),1]))
%             fprintf(vtkFile,'\nSCALARS f_%03d float 1\n',ii);
%             fprintf(vtkFile,'LOOKUP_TABLE default\n');
%             fprintf(vtkFile,fmtE,reshape(tps,[numel(tps),1]));
        end
    % case of only one correlation
    elseif isequal(size(tps),size(ms))
        fprintf('Writing two-point statistics...\n');
        writeScalars(vtkFile,'f','float',fmtE,reshape(tps,[numel(tps),1]))
    elseif (numel(tps)/size(tps,1)) ~= prod(dims) && ~isempty(tps)
        fprintf('Microstructure domain and stats sizes do NOT match! Stats data skipped\n');
    end
    
    % Write PC bases
    % check if basis size matches microstructure domain
    if size(basis,1) == prod(dims)
        fprintf('Writing PC bases...\n');
        for ii = 1:size(basis,2)
            dname = sprintf('basis_%02d',ii);
            writeScalars(vtkFile,dname,'float',fmtE,basis(:,ii))
%             fprintf(vtkFile,'\nSCALARS basis_%02d float 1\n',ii);
%             fprintf(vtkFile,'LOOKUP_TABLE default\n');
%             fprintf(vtkFile,fmtE,basis(:,ii));
        end        
    elseif size(basis,1) ~= prod(dims) && ~isempty(basis)
        fprintf('Microstructure domain and bases sizes do NOT match! Basis data skipped\n'); 
    end

    %% Closing
    fclose(vtkFile);
    fprintf('\n*** Microstructure data are written to "%s" ***\n',vtkFileName);
end

function writeScalars(vtkFile,dname,dtype,fmt,data)
    % helper function for writing scalar data to vtk files
    
    fprintf(vtkFile,'\nSCALARS %s %s 1\n',dname,dtype);
    fprintf(vtkFile,'LOOKUP_TABLE default\n');
    fprintf(vtkFile,fmt,data);
end 