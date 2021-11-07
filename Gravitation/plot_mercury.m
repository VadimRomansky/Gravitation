clear;
data = importdata('output.dat');
mercury = importdata('mercury.dat');
Nt = size(data,1);
Nb = (size(data,2) -1)/6;
Np = size(mercury,1);

figure(1);
hold on;
plot(data(1:Nt,8), data(1:Nt,9));

figure(2);
hold on;
plot(mercury(1:Np,2), mercury(1:Np,3));

figure(3);
hold on;
plot(mercury(1:Np,1), mercury(1:Np,5));

delta = (mercury(Np,5) - mercury(1,5))*3600;