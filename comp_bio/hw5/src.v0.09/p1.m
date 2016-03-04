clear all; close all

density = [0.845 .7 .6 .5 .4 .3 .2];
d845 = load('coll_time_845.dat');
d7 = load('coll_time_7.dat');
d6 = load('coll_time_6.dat');
d5 = load('coll_time_5.dat');
d4 = load('coll_time_4.dat');
d3 = load('coll_time_3.dat');
d2 = load('coll_time_2.dat');
d1 = load('coll_time_1.dat');

charTimes = [];

[h, x] = hist(d845(end/2:end,2),100);
gidx = find(h==0);
gidx = gidx(1)-1;
figure; plot(x(1:gidx(1)), log(h(1:gidx(1))) );
fit845 = polyfit(x(1:gidx(1)), log(h(1:gidx(1))),1);
charTimes = [charTimes -1/fit845(1)];

[h, x] = hist(d7(end/2:end,2),100);
gidx = find(h==0);
gidx = gidx(1)-1;
figure; plot(x(1:gidx(1)), log(h(1:gidx(1))) );
fit7 = polyfit(x(1:gidx(1)), log(h(1:gidx(1))),1);
charTimes = [charTimes -1/fit7(1)];

[h, x] = hist(d6(end/2:end,2),100);
gidx = find(h==0);
gidx = gidx(1)-1;
figure; plot(x(1:gidx(1)), log(h(1:gidx(1))) );
fit6 = polyfit(x(1:gidx(1)), log(h(1:gidx(1))),1);
charTimes = [charTimes -1/fit6(1)];

[h, x] = hist(d5(end/2:end,2),100);
gidx = find(h==0);
gidx = gidx(1)-1;
figure; plot(x(1:gidx(1)), log(h(1:gidx(1))) );
fit5 = polyfit(x(1:gidx(1)), log(h(1:gidx(1))),1);
charTimes = [charTimes -1/fit5(1)];

[h, x] = hist(d4(end/2:end,2),100);
gidx = find(h==0);
gidx = gidx(1)-1;
figure; plot(x(1:gidx(1)), log(h(1:gidx(1))) );
fit4 = polyfit(x(1:gidx(1)), log(h(1:gidx(1))),1);
charTimes = [charTimes -1/fit4(1)];

[h, x] = hist(d3(end/2:end,2),100);
gidx = find(h==0);
gidx = gidx(1)-1;
figure; plot(x(1:gidx(1)), log(h(1:gidx(1))) );
fit3 = polyfit(x(1:gidx(1)), log(h(1:gidx(1))),1);
charTimes = [charTimes -1/fit3(1)];

[h, x] = hist(d2(end/2:end,2),100);
gidx = find(h==0);
gidx = gidx(1)-1;
figure; plot(x(1:gidx(1)), log(h(1:gidx(1))) );
fit2 = polyfit(x(1:gidx(1)), log(h(1:gidx(1))),1);
charTimes = [charTimes -1/fit2(1)];

figure; plot(density, charTimes, 'x-')