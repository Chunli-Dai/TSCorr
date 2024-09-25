function [xo,yo,id] = getxyq(rang0)
%Creat a file xyq1.dat, which selecting all points from Higlist2023.gmt that are within the tile.
% Output: x, y coordinates in the given projection, point id index (line number).

    % Replace 'path_to_file.gmt' with the path to your .gmt file
    filePath ='/blue/chunlidai/apps/horizontalv/code/Higlist2023.gmt';
    constant
    % Initialize output arrays
    x = [];
    y = [];
    id = [];

    % Open the .gmt file
    fileID = fopen(filePath, 'r');
    
    % Check if the file is opened successfully
    if fileID == -1
        error('File cannot be opened: %s', filePath);
    end

if exist(filePath,'file')
fidq1 = fopen(filePath);
nq1 = linecount(fidq1);
fidq1 = fopen(filePath);
lonlat=fscanf(fidq1, '%f', [2, nq1])';
end

[x,y]=latlon2xy(lonlat(:,2),lonlat(:,1),projgdal);

M=x>=rang0(1)&x<=rang0(2)&y>=rang0(3)&y<=rang0(4);

xo=x(M);
yo=y(M);
id=find(M==1);
output=[xo(:),yo(:),id(:)];

fprintf(['\n Total of ',num2str(length(xo)),' Points within selected tile Higlist2023. \n'])

    % Close the file
    fclose(fileID);

    if 0
    %Does not work; when files are first compiled, it only read the old data.

% Output: x, y coordinates in the given projection, point id index (line number).
    filenameq1='xyq1.dat'; %

    %add lines to the top of existing files
    if exist(filenameq1,'file')
       fprintf(['\n getxyq.m reading existing file, ',filenameq1,'.\n']);
       fidq1 = fopen(filenameq1);
       nq1 = linecount(fidq1);
       fprintf(['\n ',filenameq1,' has ',num2str(nq1),' lines of data. \n']);
       fidq1 = fopen(filenameq1);
       xyq1=fscanf(fidq1, '%f', [3, nq1])';
       [nj,ni]=size(xyq1);
       if nj~=nq1
	       fprintf(['\n The format of ',filenameq1,' is wrong, should be n by 3.\n']);
	       fprintf(['\n Writing ',num2str(length(xo)),' points and replacing file, ',filenameq1,'.\n']);

       else
	     fprintf(['\n Adding ',num2str(length(xo)),' more points to file, ',filenameq1,'.\n']);
	     output=[output;xyq1];
       end
    end

    % Write the matrix to a text file
    % replace the old file;
    fidq1 = fopen(filenameq1,'w');
    writematrix(output, filenameq1, 'Delimiter', '\t');
    end % if 0

    return
end
