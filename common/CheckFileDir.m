function CheckFileDir(fileName)

    check = exist([cd filesep, fileName], 'dir');
    if check == 0
           mkdir(fileName)
    end
end