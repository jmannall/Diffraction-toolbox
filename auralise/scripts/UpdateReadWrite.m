function [read, write] = UpdateReadWrite(read, write, overlap, delay, k)
    read = read - overlap + (delay(max(1, k - 1)) - delay(k));
    write = write - overlap;
end