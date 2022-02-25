M = [[-1 2];[-1 -2]];
rng = [-3:0.1:3];

%Question 3
xdata = [0 0.1 1];
ydata = [0 1 1];
scatter(xdata, ydata, 200, 'filled')
hold on;
xvals=linspace(0,5,51);
p2=polyfit(xdata,ydata,2)
yvals=polyval(p2,xvals);
plot(xvals,yvals) 

p1=polyfit(xdata,ydata,1) 
yvals=polyval(p1,xvals);
plot(xvals,yvals)

% Question 1

dynamo4x(rng, 0.02, 1000, M, [0 0]')
hold on;

% Question 2

dynamo4x(rng, 0.02, 1000, -M, [0 0]')
hold on;

function ptseq = dynamo4x( rang,step,niter,cmtx,cvct )
figure
xinit=[rang rang rang(1)*ones(size(rang)) rang(end)*ones(size(rang))];
yinit=[rang(1)*ones(size(rang)) rang(end)*ones(size(rang)) rang rang];
[xcr,ycr] = meshgrid(linspace(rang(1),rang(end),100));
nullx=cmtx(1,1)*xcr+cmtx(1,2)*ycr+cvct(1,1);
nully=cmtx(2,1)*xcr+cmtx(2,2)*ycr+cvct(2,1);
% size(xinit)
for k=1:max(size(xinit))
    ptseq=[xinit(k);yinit(k)] ;
    curr=[xinit(k);yinit(k)] ;
    for i=1:niter
        curr=onestep(curr,cmtx,cvct,step) ;
        ptseq=[ptseq,curr] ;
    end
    plot(ptseq(1,:),ptseq(2,:),'blue')
    axis equal
    hold on ;
end
axis([rang(1) rang(end) rang(1) rang(end)])
contour(xcr,ycr,nullx,[0 0],'red','LineWidth',2) ;
contour(xcr,ycr,nully,[0 0],'red','LineWidth',2) ;
qplot(cmtx,cvct,[rang(1):.5:rang(end)])
end

function newx = onestep( oldx,cmat,cvec,dt)
newx=oldx+(cmat*oldx+cvec)*dt ;
end

function qv = qplot( cmat,cvec,rang )
%This function generates a quiver plot for the linear system
%specified by cmat. rang is the range (for both x and y)
[xx yy] = meshgrid(rang) ;
qv=[] ;
for xp=rang
    tcol=[] ;
    for yp=rang
        temp=.1*(cmat*[xp yp]'+cvec) ;
        tcol=[tcol,temp] ;
    end
    qv=cat(3,qv,tcol) ;
end
qx=reshape(qv(1,:,:),size(xx,1),size(xx,2)) ; 
qy=reshape(qv(2,:,:),size(yy,1),size(yy,2)) ;
quiver(xx,yy,qx,qy,'black') ;
axis equal
axis([[rang(1) rang(end)],[rang(1) rang(end)]]) ;
end