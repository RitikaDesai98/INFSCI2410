%Question 1 & 2
s1=[1 -1 1 -1 1 -1 1 -1];
s2=[1 1 1 1 -1 -1 -1 -1];
s3=[1 1 -1 -1 1 1 -1 -1];
T1 = [0 2 2];
T2 = [-2 0 2];
T3 = [1 1 1];

T = [T1',T2',T3'];

s = [s1',s2',s3'];

sn = normc(s)

ort = sn'*sn

x = T*sn';
r = x*sn

sa = (sn(:,1)+sn(:,2))/2

r = x*sa

ta = (T(:,1)+T(:,2))/2

%Question 3
randomS = rand(100,3)-.5;
sn = normc(randomS);
sn'* sn

x = T*sn';
r = x*sn

randomS = rand(10000,3)-.5;
sn = normc(randomS);
sn'* sn

x = T*sn';
r = x*sn