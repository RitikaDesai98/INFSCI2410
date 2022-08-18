reber.prob = [0 1 0 0 0 0 0; 0 0 0.5 0.5 0 0 0; 0 0 0.5 0 0.5 0 0; 0 0 0 0.5 0 0.5 0; 0 0 0 0.5 0 0 0.5; 0 0 0.5 0 0 0 0.5; 1 0 0 0 0 0 0]
reber.ind = [0 1 0 0 0 0 0; 0 0 2 3 0 0 0; 0 0 4 0 5 0 0; 0 0 0 4 0 5 0; 0 0 0 6 0 0 2; 0 0 6 0 0 0 3; 7 0 0 0 0 0 0]
reber.labels = 'BPQXRJE'

reber5=makestringlist(reber,4)

n0=initnet3srn(5,3,5,2,2,4312)

reber.labels(reber5.list)

stringproc(n0,[1 2 4 4 5 6 4 4 4 5 3 7],reber)

function strings = makestringlist(tgram,nstrings)
jj=1; % initial state
strings.list=[];
strings.states=[];
nstates=size(tgram.prob,2);
for ii=1:nstrings
    strings.ind(ii)=jj ; %index into superstring
    seq=[]; %initialize one string
    st=1;
    stlist=[] ; %initial state list
    while (st<nstates)
        rr=rand();
        cumu=0; i=0;
        while (cumu<rr)
            i=i+1;
            cumu=cumu+tgram.prob(st,i);
        end
        letter=tgram.ind(st,i) ;
        seq=[seq letter];
        stlist=[stlist st];
        st=i;
    end
    seq=[seq nstates] ; % append end character to seq
    strings.list=[strings.list seq];
    jj=jj+size(seq,2);
    strings.states=[strings.states stlist];
end
end


function netstruct=initnet3srn(n1,n2,n3,uamp,vamp,rs)
rng(rs);
netstruct.wih=uamp*(rand(n2,n1)-0.5) ;
netstruct.hh=uamp*(rand(n2,n2)-0.5) ;
netstruct.hbias=uamp*(rand(1,n2)-0.5) ;
netstruct.whout=vamp*(rand(n3,n2)-0.5) ;
netstruct.obias=vamp*(rand(1,n3)-0.5);
netstruct.context=zeros(1,n2);
end

function finalnet=bp3srn(net0,strlist,niter,eta,nlev)
netk=net0;
for i=1:niter
    ts=selectstring(strlist) ; % choses a new string from the training set
    netk.context=zeros(1,size(netk.wih,1)); % rests the context for a new string
    for j=1:size(ts,2)-1 % this loop trains a single string
        netk=cyc3srn(netk,ts(j),ts(j+1),eta,nlev) ;
    end
end
finalnet=netk;
end

function [sout,hlist,slist] = stringproc(netwk,strg,gramm)
hlist=[] ;
ctxinp=zeros(1,size(netwk.wih,1));
slist=[] ;
lets=gramm.labels(strg) 
%STRINGS!!!
s1=[] ;
sout=[];
 
for j=1:size(strg,2)-1
    hhh=hidlayersrn(strg(j),ctxinp,netwk.wih,netwk.hh,netwk.hbias,0.0);
    ou=layersig01(hhh,netwk.whout,netwk.obias);
    hlist=[hlist;hhh] ;
    s1=[s1,lets(j)];
    scell=cellstr(s1) ;
    slist=[slist;scell] ;
    sout=[sout, sprintf('%c %c',gramm.labels(strg(j)), gramm.labels(strg(j+1)))] ;
    for kk=1:size(netwk.whout,1)
        sout=[sout sprintf('%6.3f',ou(kk))];
    end
    sout=[sout, sprintf('\n')];
    
    ctxinp=hhh;
end
end
