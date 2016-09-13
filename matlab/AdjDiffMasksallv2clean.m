function [averadpixsh,maxradpixsh,avecolpixsh,SigM] =AdjDiffMasksallv2clean(DWI,ndir)

[Nx,Ny,numimgs]=size(DWI);
nb0=numimgs-ndir;

nrow=9;
ncol=8;
N=Nx;
mask=DWI*0;
lenfillchk=7000*(N/128)^2;

% find mask
lenfillb0(1:nb0)=0;
for i=1:ndir+nb0

    DWItmp=DWI(:,:,i);
    DWItmpsm=medfilt2(DWItmp,[11 11]);
    DWItmpsm2=medfilt2(DWItmpsm,[11 11]);
    DWItmpsm=DWItmpsm2;

    e=edge(DWItmpsm,'canny'); %finds lots of edges, works for ALL 3 sites

    e(1:15,:)=0;
    e(:,1:15)=0;
    e(N-15:N,:)=0;
    e(:,N-15:N)=0;
    radc=30;
    circ=makecirc(N,N/2,N/2,radc);
    e(find(circ))=0; % remove some central edges from 'canny'

    ef2=imfill(e,[N/2,N/2]); % flood from central points (to overcome closed edges within phantom)
    lenfill=length(find(ef2));
    if i<nb0+1
        lenfillb0(i)=lenfill;
        AVElenfill=mean(lenfillb0(find(lenfillb0)));
    end

    se=strel('disk',3,0);
    if lenfill>lenfillchk

        % all got filled so border is not closed
        e2=imdilate(e,se);
        locations=find(circ);
        ef2=imfill(e2,locations);
        ef2=imfill(ef2,'holes');
        ef2=imerode(ef2,se); % erode the first imdilate around phantom
        lenfill=length(find(ef2));
        lenfillb0(i)=lenfill;
        AVElenfill=mean(lenfillb0(find(lenfillb0)));

        if lenfill>lenfillchk
            e2=imdilate(e,se);% try to dilate x2
            e2=imdilate(e2,se);
            ef2=imfill(e2,locations);
            ef2=imerode(ef2,se); % erode x2
            ef2=imerode(ef2,se);
            lenfill=length(find(ef2));
            if i==1
                lenfillb0(i)=lenfill;
                AVElenfill=lenfill;
            end
        end
    end

    if lenfill<AVElenfill*0.97
        clear locations circ
        radc=radc+2;
        circ=makecirc(N,N/2,N/2,radc);
        locations=find(circ);
        ef3=imfill(ef2,locations);
        lenfill=length(find(ef3));
        lenfill=length(find(ef3));

	    if lenfill>lenfillchk
            ef3d=imdilate(ef2,se);% try to dilate x2
            ef3d=imdilate(ef3d,se);
            ef3f=imfill(ef3d,locations);
            ef3f=imerode(ef3f,se); % erode x2
            ef3f=imerode(ef3f,se);
            lenfill=length(find(ef3f))
            ef3=ef3f;
        end

        while lenfill<AVElenfill*0.97
            clear locations circ
            radc=radc+2;
            circ=makecirc(N,N/2,N/2,radc);
            locations=find(circ);
            ef3f=imfill(ef3,locations);
            lenfill=length(find(ef3f));

            if lenfill>lenfillchk

                ef3d=imdilate(ef3,se);% try to dilate x2
                ef3d=imdilate(ef3d,se);
                ef3f=imfill(ef3d,locations);
                ef3f=imerode(ef3f,se); % erode x2
                ef3f=imerode(ef3f,se);
                lenfill=length(find(ef3f));
            end
            ef3=ef3f;
        end
        ef2=ef3;
    end

    ef=ef2;
    ef1=imerode(ef,se); % erode outer edges
    ef2=imdilate(ef1,se); % get back the phantom edge
    ef2=imfill(ef2,'holes'); % fill any leftover holes in phantom
    lenfill=length(find(ef2));

    if lenfill>lenfillchk
        'ERROR: edge did not close!'
        pause
    end

    if i<nb0+1, lenfillb0(i)=lenfill;
        AVElenfill=mean(lenfill(find(lenfill)));
    end
    mask(:,:,i)=ef2(:,:);
end

SigM=mask;

for i=1:numimgs
    DY(:,:,i)=abs(SigM(:,:,i)-SigM(:,:,1));
    clear tmp
    tmp(:,:)=DY(:,:,i);
    tmp2=tmp;
    tmp2(find(tmp==2))=0;
    DY(:,:,i)=tmp2(:,:);
end

clear tmp*

% radial mask thickness search

numrad=[];
whichrad=[];
nntotmax=[];
angstep=2;
numradlines=length([1:angstep:359]);
Mask1=mask(:,:,1);

[rmask,cmask]=find(Mask1);
rmaskmin=min(rmask);
rmaskmax=max(rmask);
radx=(rmaskmax-rmaskmin)/2;
xlen=radx+floor(5*Nx/128);
xctr=rmaskmin+radx;

cmaskmax=max(cmask);
cmaskmin=min(cmask);
rady=(cmaskmax-cmaskmin)/2;
ylen=rady+(5*Nx/128);
yctr=cmaskmin+rady;

PhantomCentre=[xctr,yctr];
PhantomB0Dist=rady/radx;

for i=2:numimgs
    tmp(:,:)=DY(:,:,i);

    [nrad,wrad,nm]=checkradpixsh(tmp,angstep,xlen);

    numrad(i)=nrad;
    whichrad{i}=wrad(:);
    nntotmax{i}=nm(:);

    clear nrad wrad nm

end

averadpixsh(2:numimgs)=0;
maxradpixsh(2:numimgs)=0;
for i=2:numimgs
    % find the max pixshift value for each grad dir and the radial line with max pixel shift
    clear tmp
    tmp=nntotmax{i};
    [maxpixsh,maxradnj]=max(tmp);
    maxnumpixshift(i)=maxpixsh;
    whrad=whichrad{i};
    maxradj(i)=whrad(maxradnj);
    if maxradj(i)>0
         maxang(i)=(maxradj(i)-1)*angstep;
    else
         maxang(i)=0;
    end
    maxradpixsh(i)=max(tmp);
end

for i=2:numimgs
    tmp=nntotmax{i};
    whrad=whichrad{i};
    averadpixsh(i)=mean(tmp);
end
clear tmp nntot*

% find adj mask pixels for thickness along y
for i=2:numimgs
    tmp(:,:)=DY(:,:,i);
    for cj=1:Ny
        n(cj)=0;
        ntot=[];
        for rj=1:Nx
            if tmp(rj,cj)==1
                 n(cj)=n(cj)+1;
            elseif n(cj)>0
                ntot=[ntot n(cj)];
                n(cj)=0;
            end
        end
        nntot{cj}=ntot;
    end

    ncj=0;
    for cj=1:Ny
        if nntot{cj}
            ncj=ncj+1;
            col(i,ncj)=cj;
            nntotmax(i,ncj)=max(nntot{cj});
        end
    end
    numcol(i)=ncj;
end

avecolpixsh(1:numimgs)=0;

for i=2:numimgs
    tmp(3:numcol(i)-2)=nntotmax(i,3:numcol(i)-2);
    if isempty(tmp(3:numcol(i)-2))
        avecolpixsh(i)=0;
        maxcolpixsh(i)=0;
    else
        avecolpixsh(i)=mean(tmp(3:numcol(i)-2));
        maxcolpixsh(i)=max(tmp(3:numcol(i)-2));
    end
end

