%training set
cxort.smat=rand(20,2)
cxort.tmat=((cxort.smat(:,1)-.5).*(cxort.smat(:,2)-.5))>0
cxort.tmat=double(cxort.tmat)
colort=[cxort.tmat,zeros(20,1),1-cxort.tmat];
figure
scatter(cxort.smat(:,1),cxort.smat(:,2),200,colort,'filled')
%validation set
cxorv.smat=rand(20,2)
cxorv.tmat=((cxorv.smat(:,1)-.5).*(cxorv.smat(:,2)-.5))>0
cxorv.tmat=double(cxorv.tmat)
colorv=[cxorv.tmat,zeros(20,1),1-cxorv.tmat];
hold on
scatter(cxorv.smat(:,1),cxorv.smat(:,2),50,colorv,'filled')
%test set
cxors.smat=rand(20,2)
cxors.tmat=((cxors.smat(:,1)-.5).*(cxors.smat(:,2)-.5))>0
cxors.tmat=double(cxors.tmat)
colors=[cxors.tmat,zeros(20,1),1-cxors.tmat];
scatter(cxors.smat(:,1),cxors.smat(:,2),200,colors)

cxort.smat=rand(20,2)
cxort.tmat=((cxort.smat(:,1)-.5).*(cxort.smat(:,2)-.5))>0
cxort.tmat=double(cxort.tmat)
cxorv.smat=rand(20,2)
cxorv.tmat=((cxorv.smat(:,1)-.5).*(cxorv.smat(:,2)-.5))>0
cxorv.tmat=double(cxorv.tmat)
cxors.smat=rand(20,2)
cxors.tmat=((cxors.smat(:,1)-.5).*(cxors.smat(:,2)-.5))>0
cxors.tmat=double(cxors.tmat)
scatter(cxort.smat(:,1),cxort.smat(:,2),50+200*cxort.tmat,'b','filled')
hold on
scatter(cxorv.smat(:,1),cxorv.smat(:,2),50+200*cxorv.tmat,'r','filled')
scatter(cxors.smat(:,1),cxors.smat(:,2),50+200*cxors.tmat,'g','filled')

inet=initnet3(2,2,1,2,2)

[bnet,fnet, errs]=bp3trvats(inet,cxort,cxorv,cxors,20000,.2,0,.04);
errgr=errs([1:50:10001],:);
figure
plot(errgr(:,1),errgr(:,2),'b')
hold on
plot(errgr(:,1),errgr(:,3),'r')
plot(errgr(:,1),errgr(:,4),'g')

cancer=[19 21 21 21 21 23 23 23 23 23 23 23 23 23 23 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25, 25 25 25 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 29 29 29 29 29 29 29 29 29 29 29 29 29 31 31 31 31 31 31 31 31 31 31 31 31 31 33 33 33 33 35 35 35 35 35 35 35 35 35 35 37 37 37 37 37 37 39];
no_cancer=[17 17 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 29 29 29 29 29 29 29 29 29 29 31 31 31 33 33 33 33 33 33 33 33 35 35 37 37 37 39];
clas=[ones(1,100),zeros(1,200)];
[x,y,t,auc]=perfcurve(clas,[cancer,no_cancer],1);
figure
plot(x,y)
hold on 
scatter(x,y,200,'filled')

% Compute the number of true positives, true negatives, false positives and false negatives at each BMI level.
x
y

train.smat=2*rand(30,2)-1
dtrain = sqrt(diag(train.smat*train.smat'));
train.tmat=(dtrain<.85).*(dtrain>.45)
ktrain=[train.tmat,zeros(30,1),1-train.tmat];
figure
scatter(train.smat(:,1),train.smat(:,2),300,ktrain,'filled')
test.smat=2*rand(20,2)-1
dtest = sqrt(diag(test.smat*test.smat'));
test.tmat=(dtest<.85).*(dtest>.45)
ktest=[test.tmat,zeros(20,1),1-test.tmat];
hold on
scatter(test.smat(:,1),test.smat(:,2),100,ktest,'filled')

trainet=initnet3(2,4,1,2,2)
act0=forw3(trainet,train)
nf0=bp3(trainet,train,20000,.02,0)

testnet=initnet3(2,4,1,2,2)
act1=forw3(testnet,test)
nf1=bp3(testnet,test,20000,.02,0)
act10k=forw3(nf,test)

function [bestnet,finalnet,errmat]=bp3trvats(net0,trlist,valist,tslist,niter,eta,nlev,eps)
netk=net0;
maxiter=niter;
ase=2*eps ;
valbest=10;
i=0 ;
errmat=[] ;
bestnet=netk ;
while ((ase>eps)&&(i<maxiter))
    netk=cyc3(netk,trlist,eta,nlev) ;
    activ=forw3(netk,trlist) ;
    diff = trlist.tmat-activ.out;
    sse=sum(sum(diff.*diff));
    trerr=sse/prod(size(diff)) ;
    activ=forw3(netk,valist) ;
    diff = valist.tmat-activ.out;
    sse=sum(sum(diff.*diff));
    valerr=sse/prod(size(diff)) ;
    if (valerr<valbest) 
        ibest=i ;
        valbest=valerr;
        bestnet=netk;
    end
    activ=forw3(netk,tslist) ;
    diff = tslist.tmat-activ.out;
    sse=sum(sum(diff.*diff));
    tserr=sse/prod(size(diff)) ;
    errmat=[errmat;[i,trerr,valerr,tserr]];
    i=i+1 ;
end
finalnet=netk;
finalnet.iter=i;
finalnet.err=ase;
end


function netact=forw3(netwk,pats)
netact.stim=pats.smat;
netact.hid=layersigpn(netact.stim,netwk.wih,netwk.hbias) ;
netact.out=layersig01(netact.hid,netwk.whout,netwk.obias) ;
end

function netstruct=initnet3(n1,n2,n3,uamp,vamp)
netstruct.wih=uamp*(rand(n2,n1)-0.5) ;
netstruct.hbias=uamp*(rand(1,n2)-0.5) ;
netstruct.whout=vamp*(rand(n3,n2)-0.5) ;
netstruct.obias=vamp*(rand(1,n3)-0.5);
end

function lout=layersigpn(x,w,b)
lsz=size(x) ;
lll=b'*ones(lsz(1),1)'+w*x' ;
lout=sigpn(lll)';
end

function squash=sigpn(x)
squash=(1-exp(-x))./(1+exp(-x)) ;
end

function lout=layersig01(x,w,b)
lsz=size(x) ;
lll=b'*ones(lsz(1),1)'+w*x' ;
lout=sig01(lll)';
end

function squash=sig01(x)
squash=1./(1+exp(-x)) ;
end

function finalnet=bp3(net0,patstr,niter,eta,nlev)
netk=net0;
for i=10000:niter
    netk=cyc3(netk,patstr,eta,nlev) ;
end
finalnet=netk;
end

function newstruct=cyc3(nstruct,pstruct,dt,noi)
newstruct=nstruct;
szs=size(pstruct.smat);
patk=ceil(szs(1)*rand());
activity=forw1p3(nstruct,pstruct,patk,noi);
odelt=(pstruct.tmat(patk,:)-activity.out);
hdelt=0.5*(nstruct.whout'*odelt').*(1+activity.hid').*(1-activity.hid');
newstruct.whout=newstruct.whout+dt*odelt'*activity.hid;
newstruct.obias=newstruct.obias+dt*odelt ;
newstruct.wih=newstruct.wih+dt*hdelt*activity.stim;
newstruct.hbias=newstruct.hbias+dt*hdelt' ;
end

function netact=forw1p3(netwk,pats,patno,nois)
netact.stim=pats.smat(patno,:);
netact.hid=nlayersigpn(netact.stim,netwk.wih,netwk.hbias,nois) ;
netact.out=layersig01(netact.hid,netwk.whout,netwk.obias) ;
end

function lout=nlayersigpn(x,w,b,noise)
lsz=size(x) ;
% lll=b'*ones(lsz(1),1)'+w*x'
bb=repmat(b',1,lsz(1)) ;
lll=bb+w*x' ;
lout=sigpn(lll)'+noise*(rand(size(lll'))-0.5) ;
end