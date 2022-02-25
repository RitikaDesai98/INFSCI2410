rng('default')
locations=rand(6,2); 
scatter(locations(:,1),locations(:,2),400,'r','filled')
dmat=squareform(pdist(locations)) 

tspmat=hopfieldwts(6,10,dmat);
iacn(.05*rand(6,6),0.2*rand(6,6),tspmat,.05,100000)


rseq=[6 3 5 1 2 4];
hold on
plot(locations(rseq,1),locations(rseq,2),'k')




function hmtx=hopfieldwts(nc,mag,distances)
hmtx=zeros(nc,nc,nc,nc) ;
mdist=max(max(distances)) ;
nsd=10.0 ;
revdist=10*(mdist-distances+1).*(mdist-distances+1)/(mdist*mdist) ;
for j=1:nc
hmtx(j,:,j,:)=mag*(eye(nc)-1) ;
hmtx(:,j,:,j)=mag*(eye(nc)-1) ;
end
for j=2:nc 
for k=1:(j-1)
for m1=1:nc
m2=1+mod(m1,nc) ;
hmtx(j,m1,k,m2) = revdist(j,k) ;
hmtx(k,m1,j,m2) = revdist(j,k) ;
hmtx(j,m2,k,m1) = revdist(j,k) ;
hmtx(k,m2,j,m1) = revdist(j,k) ;
end
end
end
end


function newact=iaciter(extinp,oldact,cmat,lr)
del=mul4d2d(cmat,oldact) + extinp ;
ldel = 2./exp(-del) - 1.0 ; 
newact = oldact + lr*(ldel>0).*(1-oldact) - lr*(ldel<0).*oldact ;
end


function finalact=iacn(extin,initact,conmat,dt,niter)
finalact=initact;
for k=1:niter
finalact=iaciter(extin,finalact,conmat,dt);
end
end


function outm=mul4d2d(m4d,m2d)
newdim=prod(size(m2d));
m4d2=reshape(m4d,newdim,newdim) ;
m2d1=reshape(m2d,newdim,1);
penout=m4d2*m2d1 ;
outm=reshape(penout,size(m2d)) ;
end
