%Inclass 14
%GB comments
100 Could run a loop on the threshold parameters and define the “best” mask
2a 100 
2b 100

%Work with the image stemcells_dapi.tif in this folder

% (1) Make a binary mask by thresholding as best you can
% 
img=imread('stemcells_dapi.tif');
img_mask = img > 300; %works better than 250
figure (2); imshow(img_mask, []);
cprop=regionprops(img_mask,'PixelIdxList','Area','Centroid');
areas=[cprop.Area];

average=mean(areas);
stdev=std(areas);
ids=find(areas>average/4);
isize=size(img);
nmask=false(isize);
lids=length(ids);

for ii=1:lids
nmask(cprop(ids(ii)).PixelIdxList)=true;
end
fmask=imdilate(nmask,strel('disk',2));
figure(3); imshow(fmask,[]);

% (2) Try to separate touching objects using watershed. Use two different
% ways to define the basins. (A) With erosion of the mask (B) with a distance transform. Which works better in this case?

img3=fmask;
CC=bwconncomp(img3);
stats=regionprops(CC,'Area');
area3=[stats.Area];
fusedCandidates=area3>mean(area3)+std(area3);
sublist=CC.PixelIdxList(fusedCandidates);
sublist=cat(1,sublist{:});
fusedMask=false(size(img3));
fusedMask(sublist)=1;
%figure (5); imshow(fusedMask, []);

%This is with distance transformation

distt=bwdist(fusedMask);
%imshow(fmask,[]);
D=bwdist(~fusedMask);
%imshow(D,[], 'InitialMagnification','fit');
D=-D;
D(~fusedMask)=-inf;
%imshow(D,[], 'InitialMagnification','fit');
L=watershed(D);
rgb=label2rgb(L,'jet',[.5 .5 .5 ]);
dmask=L>1 | (img3-fusedMask);
figure (4); imshow(dmask, 'InitialMagnification','fit');

% with erosion mask

s=round(0.75*sqrt(mean(area3))/pi); %had to decrease the erosion
nucmin=imerode(fusedMask,strel('disk',s));
%figure (6); imshow(nucmin,[]);

outside= ~imdilate(fusedMask,strel('disk',1));
%figure(7); imshow(outside,[]);

basin=imcomplement(bwdist(outside));
basin=imimposemin(basin, nucmin | outside);
pcolor (basin); shading flat;

L3=watershed(basin);
%figure (8); imshow(L,[]); colormap('jet');
newnmask=L3>1 | (img3-fusedMask);
figure(9); imshow(newnmask,[]);

%Miguel Angel: The Erosion mask seems to be cleaner, although a lot of big
%cells don't get watershed, this has to be to the relationship between the
%area and what I denomated "big". I used the mean area + stdev, but perhaps
%there are fused cells closely lower to that value and that's why they
%don't get picked up. For the distance transformation it seems to have
%picked up more big cells, however I encountered the problem discussed in
%class where we will get a lot of oversegmentation and that will need to
%raise the minima or tell which minima to look for.
