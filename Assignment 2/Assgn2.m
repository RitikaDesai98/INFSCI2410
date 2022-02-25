%Question - 1
rng('default')
w = rand(1,3);
bias = 1;

for k=1:5
    perception(0,0,0) % False AND False
    perception(0,1,0) % False AND True
    perception(1,0,0) % True AND False
    perception(1,1,1) % True AND True
end

a = 1;
b = 1;
op = a*w(1,1) + b*w(1,2) + bias*w(1,3);
if op > 0
   op = 1;
else
   op = 0;
end
disp(op)

contour(a,b,op)

function perception(ip1, ip2, op)

w = rand(1,3);
bias = 1;
lr = 1; %learning rate

op1 = (ip1*w(1,1)) + (ip1*w(1,2)) + (bias*w(1,3));
    if op1 > 0  %activation function
        op1 = 1;
    else
        op1 = 0;
    end
error = op - op1;
w(1,1) = w(1,1) + (error * ip1 * lr);
w(1,2) = w(1,2) + (error * ip2 * lr);
w(1,3) = w(1,3) + (error * bias * lr);
end


%Question - 2

[a b]=meshgrid(-3:.1:3)
for j=1:61
for k=1:61
z(j,k)=rbu([a(j,k),b(j,k)],[1,-1]);
end
end
figure
contour(a,b,z)
hold on;


% Question - 3

[x1 y1]=meshgrid(-3:.1:3);
for a=1:61
for b=1:61
z1(a,b)=rbu([x1(a,b),y1(a,b)],[1,-1]);
end
end
figure
contour(x1,y1,z1)
hold on

zans = z+z1;
figure
contour(zans)


function rr = rbu( ss,proto )
% Radial Basis Unit
% exponential of negative squared distance bewteen ss and proto
rr = exp(-sum((ss-proto).*(ss-proto))) ;
end

% Question - 5

rng('default')
w = rand(1,7);

for k=1:5
    perception(0,0,1,1,0,1,1)

end

x1 = 0;
x2 = 0;
x3 = 1;
x4 = 1;
x5 = 0;
x6 = 1;

op = x1*w(1,1) + x2*w(1,2) + x3*w(1,3) + x4*w(1,4) + x5*w(1,5) + x6*w(1,6) + bias*w(1,7);
if op > 0
   op = 1;
else
   op = 0;
end
disp(op)


function perception(ip1, ip2, ip3, ip4, ip5, ip6, op)

w = rand(1,7);
bias = 1;
lr = 1; %learning rate

op1 = (ip1*w(1,1)) + (ip2*w(1,2)) + (ip3*w(1,3)) + (ip4*w(1,4)) + (ip5*w(1,5)) + (ip6*w(1,6)) + (bias*w(1,7)) ;
    if op1 > 0  %activation function
        op1 = 1;
    else
        op1 = 0;
    end
error = op - op1;
w(1,1) = w(1,1) + (error * ip1 * lr);
w(1,2) = w(1,2) + (error * ip2 * lr);
w(1,3) = w(1,3) + (error * ip3 * lr);
w(1,4) = w(1,4) + (error * ip4 * lr);
w(1,5) = w(1,5) + (error * ip5 * lr);
w(1,6) = w(1,6) + (error * ip6 * lr);
w(1,7) = w(1,7) + (error * bias * lr);
end

% Question - 6

rng('default')
w = rand(1,7);

for k=1:5
    perception(0,0,1,1,0,1,1)
    perception(1,0,1,1,0,1,1)

end

x1 = 0;
x2 = 0;
x3 = 1;
x4 = 1;
x5 = 0;
x6 = 1;

op = x1*w(1,1) + x2*w(1,2) + x3*w(1,3) + x4*w(1,4) + x5*w(1,5) + x6*w(1,6) + bias*w(1,7);
if op > 0
   op = 1;
else
   op = 0;
end
disp(op)


function perception(ip1, ip2, ip3, ip4, ip5, ip6, op)

w = rand(1,7);
bias = 1;
lr = 1; %learning rate

op1 = (ip1*w(1,1)) + (ip2*w(1,2)) + (ip3*w(1,3)) + (ip4*w(1,4)) + (ip5*w(1,5)) + (ip6*w(1,6)) + (bias*w(1,7)) ;
    if op1 > 0  %activation function
        op1 = 1;
    else
        op1 = 0;
    end
error = op - op1;
w(1,1) = w(1,1) + (error * ip1 * lr);
w(1,2) = w(1,2) + (error * ip2 * lr);
w(1,3) = w(1,3) + (error * ip3 * lr);
w(1,4) = w(1,4) + (error * ip4 * lr);
w(1,5) = w(1,5) + (error * ip5 * lr);
w(1,6) = w(1,6) + (error * ip6 * lr);
w(1,7) = w(1,7) + (error * bias * lr);
end

