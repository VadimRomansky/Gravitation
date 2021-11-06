clear;
data = importdata('general.dat');

Nt = size(data,1);

figure(1);
plot(1:Nt, data(1:Nt,1),'Color','red','LineWidth',2);
title('E');

figure(2);
hold on;
plot(1:Nt, data(1:Nt,2),'Color','red','LineWidth',2);
plot(1:Nt, data(1:Nt,3),'Color','green','LineWidth',2);
plot(1:Nt, data(1:Nt,4),'Color','blue','LineWidth',2);
title('p');

figure(3);
hold on;
plot(1:Nt, data(1:Nt,5),'Color','red','LineWidth',2);
plot(1:Nt, data(1:Nt,6),'Color','green','LineWidth',2);
plot(1:Nt, data(1:Nt,7),'Color','blue','LineWidth',2);
title('rc');

figure(4);
hold on;
plot(1:Nt, data(1:Nt,8),'Color','red','LineWidth',2);
plot(1:Nt, data(1:Nt,9),'Color','green','LineWidth',2);
plot(1:Nt, data(1:Nt,10),'Color','blue','LineWidth',2);
title('L');