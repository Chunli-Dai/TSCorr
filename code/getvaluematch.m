
function values = getvaluematch(targetString, instr)
    % Given input string (char), instr, find the two values after the match
    % e.g., Upper Left  (  321381.000, 4165044.000) ( 36d58'34.29"E, 37d36'54.59"N)

    % Find the line containing the target string
    matchStart = strfind(instr, targetString);

    if ~isempty(matchStart)
        % Extract the substring starting from the targetString
	matchedLine = instr(matchStart:matchStart+70) ;

	matchedLine=strrep(matchedLine,targetString,'');

        % Define the format to read two floating-point numbers
        %formatSpec = '%*s  (%f, %f)';
        formatSpec = '%f, %f)';

        % Use sscanf to read the numbers
        values = sscanf(matchedLine, formatSpec);
    else
        disp(['Target string "', targetString, '" not found in the input data.']);
        values = [];
    end

    return
end

