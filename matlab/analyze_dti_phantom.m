%%  analyze_dti_phantom(dwi, fa, bval, output, nyqopt)
%%
%%  'dwi':    4D diffusion weighted image
%%  'fa':     FA map from DTIfit
%%  'bval':   B value files from dcm2nii
%%  'output': folder to store outputs
%%  'nyqopt':  (1, 0) 1 to measure nyquist ghost (use on not accelerated data).

function analyze_dti_phantom(dwi, fa, bval, output, nyqopt)

%% Part 1: load data
bval = dlmread(bval);
dwi = load_nifti(dwi);
fa = load_nifti(fa);

% calculate number of directions / b0 volumes
ndir = length(find(bval>0));
nb0 = length(find(bval==0));

% create output directory
if~isdir(output)
    mkdir(output)
end

% initalize outputs
outname = strcat(output, '/main_stats.csv');
outname01 = strcat(output, '/SNR_each-b0.csv');
outname03 = strcat(output, '/SNR_each-DWI.csv');
outname04 = strcat(output, '/Nyq_each-b0.csv');
outname06 = strcat(output, '/PXShift.csv');

fid = fopen(outname, 'w');
fid01 = fopen(outname01, 'w');
fid03 = fopen(outname03, 'w');
fid04 = fopen(outname04, 'w');
fid06 = fopen(outname06, 'w');

% print headers
fprintf(fid01, '%s,%s,%s,%s,%s,%s\n', ...
               'b=0#', 'stdNoise', 'aveNoise','aveSignal', 'stdSignal', 'SNR(=aveS/stdN)');
fprintf(fid03, '%s,%s,%s,%s,%s,%s\n', ...
               'GradDir#', 'stdNoise', 'aveNoise', 'aveSignal','stdSignal','SNR(=aveS/stdN)');
fprintf(fid04, '%s,%s,%s,%s\n', ...
               'b#', 'aveNyq', 'avenoise','Nyqratio');
fprintf(fid06, '%s,%s,%s,%s\n', ...
               'GradDir#','Ave.radial pixsh', 'Max.radial pixsh', 'Ave.col pixsh');

% load DWI, averaging over 3 central slices. dim1,2=AXIAL, dim3=all directions)
dims = size(dwi.vol);
central_slice = ceil(dims(3)/2);
DWI = mean(dwi.vol(:,:, central_slice-1:central_slice+1, :), 3);
[Nx, Ny, numimgs] = size(DWI);

% load FA, taking central slice only
FA = fa.vol(:,:,central_slice);

clear fa dwi

%% Part 2: getting masks and pixel shifts
[averadpsh,maxradpsh,avecolpsh,SigM] = AdjDiffMasksallv2clean(DWI,ndir);

figure(5); set(gcf,'Visible', 'off');
    plot(averadpsh,'b*-')
    title(['Pixel shifts'])
    hold on
    plot(maxradpsh,'b*--') % gives info about a large pixel shift which may be lost in average
    plot(avecolpsh,'r*-') % don't plot max pix shift for col because misleadingly large at vertical boundaries
    legend('ave radial pixsh','max radial pixsh','ave column pixsh')

aveaverad=mean(averadpsh);
avemaxrad=mean(maxradpsh);
aveavecol=mean(avecolpsh);

for i=1:numimgs-nb0
    j=nb0+i;
    fprintf(fid06,'%s,%6.2f,%d,%6.2f\n',num2str(i,'%02d'),averadpsh(j),maxradpsh(j),avecolpsh(j));
end

%% Part 3: measure noise from difference image of b0 (phantom rad=35/128 px)
phantrad=35*(Nx/128);

Diff(1:Nx,1:Ny,1:nb0)=0;
Diff(1:Nx,1:Ny,1)=double(DWI(1:Nx,1:Ny,2)-DWI(1:Nx,1:Ny,1));
nd2TOT=[];

for b=2:nb0

    Diff(1:Nx,1:Ny,b)=double(DWI(1:Nx,1:Ny,b)-DWI(1:Nx,1:Ny,1));
    N2=floor(Nx/2);

    noisemskctr=makecirc(Nx,N2,N2,phantrad);
    Diff2d(1:Nx,1:Ny)=Diff(1:Nx,1:Ny,b);
    nd2{b-1}=Diff2d(find(noisemskctr));
    std2ALL(b)=std(nd2{b-1});
    ave2ALL(b)=mean(nd2{b-1});

    % check that noise is not Rician
    noiseratio=ave2ALL(b)/std2ALL(b);
    [nd2hist,xhist2]=hist(nd2{b-1},20);

    if b==2
        XX=Diff(1:Nx,1:Ny,b);
        XX(find(noisemskctr))=max(max(Diff(1:Nx,1:Ny,b)))*.8;
        figure(100); set(gcf,'Visible', 'off');
        subplot(2,2,1)
            imagesc(Diff(1:Nx,1:Ny,b))
            axis image; axis off
            title(['T2w image b=0 Img#',num2str(b)])
        subplot(2,2,2)
            imagesc(XX)
            axis image; axis off
        subplot(2,2,3)
            plot(nd2{b-1})
            title(['ave(noise)=',num2str(ave2ALL(b),'%5.2f'),' std(noise)=',num2str(std2ALL(b),'%5.2f'),' noiseratio=',num2str(noiseratio,'%5.2f')])
        subplot(2,2,4)
            plot(xhist2,nd2hist)
            title('noise histogram')
    end
    nd2TOT=cat(2,nd2TOT, nd2{b-1});
end

% use same noise (from img#1-img#2) for std(noise) calc for images #1 & #2
std2ALL(1)=std2ALL(2);
ave2ALL(1)=ave2ALL(2);

[nd2hist,xhist2]=hist(nd2TOT,20);
clear ave2 std2 nd2
nd2=nd2TOT(:);
ave2b0=mean(nd2);
std2b0=std(nd2);
noiseratio=ave2b0/std2b0;

figure(150); set(gcf,'Visible', 'off');
    subplot(1,2,1)
        plot(nd2TOT)
        title(['ave(noise)=',num2str(ave2b0,'%5.2f'),' std(noise)=',num2str(std2b0,'%5.2f'),' noiseratio=',num2str(noiseratio,'%5.2f')])
    subplot(1,2,2)
        plot(xhist2,nd2hist)
        title('TOTAL noise histogram')

clear Diff nd2*

ngrad=numimgs-nb0;
Diff(1:Nx,1:Ny,1:ngrad)=0;
Diff(1:Nx,1:Ny,1)=double(DWI(1:Nx,1:Ny,2+nb0)-DWI(1:Nx,1:Ny,1+nb0));
nd2TOT=[];

% loop through gradients to calculate
for b=2:ngrad
    % difference always w.r.t. first DWI
    Diff(1:Nx,1:Ny,b)=double(DWI(1:Nx,1:Ny,b+nb0)-DWI(1:Nx,1:Ny,1+nb0));

    N2=floor(Nx/2);
    noisemskctr=makecirc(Nx,N2,N2,phantrad);
    Diff2d(1:Nx,1:Ny)=Diff(1:Nx,1:Ny,b);
    nd2{b-1}=Diff2d(find(noisemskctr));
    std2ALL(nb0+b)=std(nd2{b-1});
    ave2ALL(nb0+b)=mean(nd2{b-1});

    % check that noise is not Rician
    noiseratio=ave2ALL(nb0+b)/std2ALL(nb0+b);
    [nd2hist,xhist2]=hist(nd2{b-1},20);

    % plot for the first gradient
    if b==2
        % calculate XX, whatever that is
        figure(200); set(gcf,'Visible', 'off');
        XX=Diff(1:Nx,1:Ny,b);
        XX(find(noisemskctr==0))=0;
        subplot(2,2,1)
            imagesc(Diff(1:Nx,1:Ny,b))
            axis image; axis off
            title(['DiffImg grad dir b>0 Img#',num2str(b)])
        subplot(2,2,2)
            imagesc(XX)
            axis image; axis off
        subplot(2,2,3)
            plot(nd2{b-1})
            title(['ave(noise)=',num2str(ave2ALL(nb0+b),'%5.2f'),' std(noise)=',num2str(std2ALL(nb0+b),'%5.2f'),' noiseratio=',num2str(noiseratio,'%5.2f')])
        subplot(2,2,4)
            plot(xhist2,nd2hist)
            title('noise histogram')
    end
    nd2TOT=cat(2,nd2TOT, nd2{b-1});
end

% use same noise (from dwi#1-dwi#2) for std(noise) calc for images
std2ALL(nb0+1)=std2ALL(nb0+2);
ave2ALL(nb0+1)=ave2ALL(nb0+2);

[nd2hist,xhist2]=hist(nd2TOT,20);
clear ave2 std2 nd2
nd2=nd2TOT(:);
ave2=mean(nd2);
std2=std(nd2);
noiseratio=ave2/std2;

figure(300); set(gcf,'Visible', 'off');
    subplot(1,2,1);
        plot(nd2TOT)
        title(['b>0 ave(noise)=',num2str(ave2,'%5.2f'),' std(noise)=',num2str(std2,'%5.2f'),' noiseratio=',num2str(noiseratio,'%5.2f')])
    subplot(1,2,2)
        plot(xhist2,nd2hist)
        title('TOTAL noise histogram')

% Getting SNR measurements
Stot=[];
for b=1:numimgs
    I(1:Nx,1:Ny)=DWI(:,:,b);
    Smat=I(find(noisemskctr));
    S=Smat(:)';
    aveS(b)=mean(S);
    stdS(b)=std(S);

    SNR(b)= aveS(b)/std2ALL(b);
    II=I*0;
    II(find(noisemskctr))=I(find(noisemskctr));
    Stot=cat(1,Stot,II(:)');
end

SNRb0TOT=mean(SNR(1:nb0));
STD_SNRb0TOT=std(SNR(1:nb0));
SNRTOT=mean(SNR(nb0+1:numimgs));
STD_SNRTOT=std(SNR(nb0+1:numimgs));

figure(400); set(gcf,'Visible', 'off');
    plot(SNR,'bo-'); hold on; plot([nb0+1:numimgs],SNR(nb0+1:numimgs),'b*-')
    title(['SNR b=0 images:',num2str(SNRb0TOT,'%05.2f'),'(',num2str(STD_SNRb0TOT),')  ||    SNR DWIs:',num2str(SNRTOT,'%05.2f'),'(',num2str(STD_SNRTOT),')'])

AVESNRb0=mean(SNR(1:nb0));
STDSNRb0=std(SNR(1:nb0));
CV_SNRb0=STDSNRb0*100/AVESNRb0;
AVESNRdwi=mean(SNR(1+nb0:numimgs));
STDSNRdwi=std(SNR(1+nb0:numimgs));
CV_SNRdwi=STDSNRdwi*100/AVESNRdwi;

for b=1:nb0
    fprintf(fid01,'%d,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f\n', b, std2ALL(b), ave2ALL(b),aveS(b), stdS(b),SNR(b));
end

for b=nb0+1:numimgs
    fprintf(fid03,'%d,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f\n', b-nb0, std2ALL(b), ave2ALL(b), aveS(b),stdS(b),SNR(b));
end

% Plotting results for SNR for DWIs
figure(500); set(gcf,'Visible', 'off');
    subplot(3,1,1)
        plot([1:numimgs-nb0],aveS(nb0+1:numimgs),'b*-')
        xlabel('Gradient Direction'); axis([1 numimgs-nb0 0 max(aveS(nb0+1:numimgs))*1.2])
        title('Ave. Signal')
    subplot(3,1,2)
        plot([1:numimgs-nb0], std2ALL(nb0+1:numimgs),'b*-')
        xlabel('Gradient Direction'); axis([1 numimgs-nb0 0 max(std2ALL(nb0+1:numimgs))*1.2])
        title('Std noise')
    subplot(3,1,3)
        plot([1:numimgs-nb0],SNR(nb0+1:numimgs),'b*-')
        xlabel('Gradient Direction'); axis([1 numimgs-nb0 0 max(SNR(nb0+1:numimgs))*1.2])
        title('SNR')

% Part 4: getting nyquist ghost info
if nyqopt == 1
    disp('MSG: Nyquist Ghost Measurement: ON')
    for b=1:nb0
        I_B0(1:Nx,1:Ny)=DWI(:,:,b);
        Smask=SigM(:,:,1);
        [ndNyq,XX2]=noisedistNyqSmask(I_B0,Smask);
        aveNyq(b)=mean(ndNyq);

        [ndnoise,XX3]=noisedist_DTISmask(I_B0,Smask);
        avenoise(b)=mean(ndnoise);

        XX=(XX2+XX3)/2;
        maxv=max(max(XX));

        ratNyq(b) = aveNyq(b)./avenoise(b);
        fprintf(fid04,'%d,%6.4f,%6.4f,%6.4f\n', b, aveNyq(b), avenoise(b), ratNyq(b));
    end

    STDNyqrat=std(ratNyq);
    AVENyqrat=mean(ratNyq);
    CV_Nyqrat=STDNyqrat*100/AVENyqrat;


    % Plotting results for SNR and Nyquist Ratio (b=0)
    figure(450); set(gcf,'Visible', 'off');
        subplot(4,1,1)
            plot([1:nb0],aveS(1:nb0),'b*-');
            xlabel('b=0 Image# '); axis([1 nb0 0 max(aveS(1:nb0))*1.2])
            title('Ave. Signal')
        subplot(4,1,2)
            plot([1:nb0], std2ALL(1:nb0),'b*-');
            xlabel('b=0 Image# '); axis([1 nb0 0 max(std2ALL(1:nb0))*1.2])
            title('Std noise')
        subplot(4,1,3)
            plot([1:nb0],SNR(1:nb0),'b*-');
            xlabel('b=0 Image#'); axis([1 nb0 0 max(SNR(1:nb0))*1.2])
            title('SNR')
        subplot(4,1,4)
            plot([1:nb0],ratNyq(1:nb0),'b*-');
            xlabel('b=0 Image#'); axis([1 nb0 0 max(ratNyq(1:nb0))*1.2])
            title(['Ratio of Nyquist ghost: AVE(STD) =',num2str(AVENyqrat,'%03.3f'),'(',num2str(STDNyqrat,'%03.3f'),')'])
else
    STDNyqrat = NaN
    AVENyqrat = NaN
    CV_Nyqrat = NaN
end


%% Part 5: plot FA
fav=FA(find(FA));
aveFA=mean(fav);
stdFA=std(fav);

figure(14); set(gcf,'Visible', 'off');
    subplot(2,2,1)
        imagesc(FA, [0 0.1]);
        set(gca,'DataAspectRatio',[1 1 1]); colorbar;
        title('FA map');
    subplot(2,2,2)
        plot(fav);
        title('FA values');
    subplot(2,2,3)
        [n,x]=hist(fav,30); plot(x,n,'k-');
        title(['FA mean(std)=', num2str(aveFA,'%5.3f'), '(',num2str(stdFA,'%5.3f'),')']);

fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
             'stdSNR_b0','aveSNR_b0','CoV_SNR_b0(%)','stdSNR_dwi','aveSNR_dwi','CoV_SNR_dwi(%)','AVE Ave.radpixsh','AVE Max.radpixsh','AVE Ave.colpixsh','AVE FA','STD FA','STDNyqratio','AVENyqratio','CoV_Nyqratio');
fprintf(fid, '%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.2f,%6.2f,%6.2f,%5.3f,%5.3f,%6.4f,%6.4f,%6.4f\n', ...
             STDSNRb0,AVESNRb0,CV_SNRb0,STDSNRdwi,AVESNRdwi,CV_SNRdwi,aveaverad,avemaxrad,aveavecol,aveFA,stdFA,STDNyqrat,AVENyqrat,CV_Nyqrat);

%% Part 6: print figures
fig5name=strcat(output, '/PXShift_ave.jpg'); print('-f5','-djpeg',fig5name);

% SNR
fig100name = strcat(output, '/b0_diff-roi-noise-hist.jpg');  print('-f100','-djpeg',fig100name);
fig150name = strcat(output, '/b0_noise-hist.jpg');           print('-f150','-djpeg',fig150name);
fig200name = strcat(output, '/DWI_diff-roi-noise-hist.jpg'); print('-f200','-djpeg',fig200name);
fig300name = strcat(output, '/DWI_noise-hist.jpg');          print('-f300','-djpeg',fig300name);
fig400name = strcat(output, '/SNR_avg-std.jpg');             print('-f400','-djpeg',fig400name);
fig500name = strcat(output, '/SNR_individual.jpg');          print('-f500','-djpeg',fig500name);
if nyqopt == 1
    fig450name = strcat(output, '/SNR_Nyq_eachb0.jpg'); print('-f450','-djpeg',fig450name);
end

% FA
fig14name=strcat(output, '/FAvalues.jpg'); print('-f14','-djpeg',fig14name);

% le fin
fclose all;
exit;
end
