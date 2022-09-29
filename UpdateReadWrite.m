function [read, write] = UpdateReadWrite(read, write, overlap, delay, i)
    read = read - overlap + (delay(max(1, i - 1)) - delay(i));
    write = write - overlap;
end