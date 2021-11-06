clear;
data = importdata('output.dat');
Nt = size(data,1);
Nb = (size(data,2) -1)/6;

figure(1);
hold on;
plot(data(1:Nt,1), data(1:Nt,2));

figure(2);
hold on;
for i=1:Nb,
    plot(data(1:Nt,2 + (i-1)*6), data(1:Nt,3 + (i-1)*6));
end;