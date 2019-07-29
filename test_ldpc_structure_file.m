clear;

fileID = fopen('test.dat','w','n');

fwrite(fileID, 2,'uint16');

fclose(fileID);

fileID = fopen('test.dat','r','n');

fwrite(fileID, 2, 'uint16');

