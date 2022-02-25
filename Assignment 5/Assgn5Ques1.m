xor.smat = [[0 0];[0 1];[1 0];[1 1]]
xor.tmat=[0 1 1 0]'
xornet0=initnet3(2,1,1,2,2)
act0=forw3(xornet0,xor)
nf=bp3(xornet0,xor,10000,.02,0)
act10k=forw3(nf,xor)
act10k.out
act0.stim
act0.hid
act0.out

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
for i=1:niter
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
bb=repmat(b',1,lsz(1)) ;
lll=bb+w*x' ;
lout=sigpn(lll)'+noise*(rand(size(lll'))-0.5) ;
end