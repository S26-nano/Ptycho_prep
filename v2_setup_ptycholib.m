%ptydir = '/CNMshare/savedata/2016R1/20160202/Analysis/ptycholib/src/';
%ptydir = '/CNMshare/savedata/pythonscripts/ptycholib/ptycholib/src/';
%Scan 9 setup - transmission pattern 20180705
%%{
scandir = '/CNMshare/savedata/2018R2/20180705/';
scannum = 18;
spiral = 0;
numpts = [41 41];
imagechan = 19;
numiters =20;
samdetdist = 0.8; %meters
stepsize = [10 10]; %nm
%probeguess = 'pg_foc.txt';
%probeguess = 'probe_10kev_30df_film.txt';
probeguess = '10_4_kev_foc.txt';
pxsz = 55; %um
imgsperpt = 5;
imgctr = [128 128];
imgsize = 256;
imstart = 121452; % use this to avoid errors with ccdnum readback

createh5 =1;
runptycho = 0;
showptycho = 0;
%}



%-----%

% Create h5 images from experimental data
h5dir = [scandir 'Images/' num2str(scannum,'%d') '/h5/']

if(createh5)
    unix(['mkdir ' h5dir]);
    imagenums = loadmda([scandir 'mda/26idbSOFT_' num2str(scannum,'%4.4d') '.mda'],imagechan,0,0);
    if(spiral)
        h=waitbar(0,'Converting tiff images to h5');
        for ii=1:numpts
            waitbar(ii/numpts);
            %CHECK scannum+1 next line
            filename = [scandir 'Images/' num2str(scannum) '/scan_' num2str(scannum) ...
                '_img_' num2str(imagenums(ii,2),'%6.6d') '.tif'];
            %img1 = resample_pixirad(filename,0);
            img1 = double(imread(filename));
            sz1 = size(img1);
            h5create([h5dir 'image_' num2str((ii-1),'%06d') '.h5'],'/entry/instrument/detector/data',[sz1(1) sz1(2)]);
            h5write([h5dir 'image_' num2str((ii-1),'%06d') '.h5'],'/entry/instrument/detector/data',img1);
        end
    else
        h=waitbar(0,'Converting tiff images to h5');
        ind1 = 0;%imstart = 15231;
        for ii=1:size(imagenums,1)
            waitbar(ii/size(imagenums,1));
            for jj=1:size(imagenums,2)
                clear img2;
                for kk= 1:imgsperpt
                    filename = [scandir 'Images/' num2str(scannum) '/scan_' num2str(scannum) ...
                        '_img_' num2str(imstart,'%6.6d') '.tif'];
                    %filename = [scandir 'Images/' num2str(scannum) '/scan_' num2str(scannum) ...
                    %    '_img_' num2str(imagenums(ii,jj,1),'%6.6d') '.tif'];
                    
                    if(kk==1)
                        img2 = double(imread(filename));
                    else
                        img2 = img2 + double(imread(filename));
                    end
                    imstart=imstart+1;
                end
                % Hot pixel removal here
                %%{
                img3 = zeros(size(img2,1),size(img2,2),4);
                img3(:,:,1) = circshift(img2,[0 1]);
                img3(:,:,2) = circshift(img2,[1 0]);
                img3(:,:,3) = circshift(img2,[0 -1]);
                img3(:,:,4) = circshift(img2,[-1 0]);
                img4 = median(img3,3);
                ccdmask = abs(img2-img4)>50;   %CHANGE THRESHOLD HERE
                img2 = img2.*(1-ccdmask)+img4.*ccdmask;
                %}
                % Blank out pixels here TLX TLY BRX BRY
                %%{
				masks=[200 370 225 414;227 409 272 430;216 453 253 483;496 331 514 352;324 209 348 234;...
					    256 256 260 260];
                for kk = 1:size(masks,1)
                    img2(masks(kk,2):masks(kk,4),masks(kk,1):masks(kk,3)) = 0;
                end
                %}
                % Center image here
                %%{
                %img2 = flipud(fliplr(circshift(img2,[258-313,258-311]))); %for scan 7 transmission
                %img2 = fliplr(circshift(img2,[258-313,258-311])); %for scan 9 transmission
                img2 = fliplr(circshift(img2,[258-377,258-384])); %scan 9 20180705
                %}
                % Take 512x512 subset
                %img1 = img2(3:514,3:514);
                % Take 256x256 subset
                img1 = img2(131:387,131:387);
                sz1 = size(img1);
                h5create([h5dir 'image_' num2str((ind1),'%06d') '.h5'],'/entry/instrument/detector/data',[sz1(1) sz1(2)]);
                h5write([h5dir 'image_' num2str((ind1),'%06d') '.h5'],'/entry/instrument/detector/data',img1);
                ind1=ind1+1;
            end
        end
    end
    close(h)
end

%-----%
