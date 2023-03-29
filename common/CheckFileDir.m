function CheckFileDir(fileName)

    check = exist([cd filesep, fileName], 'dir');   % Returns 7 if folder exists and MATLAB can access it
    
    if check == 0
        listing = dir(fileName);
        if isempty(listing)
           mkdir(fileName)
        end
    end
end