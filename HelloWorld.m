function HelloWorld()
    disp('Hello World')
    x = 1;
    CheckFileDir('Test')
    save(['Test' filesep, 'test'], 'x')
end