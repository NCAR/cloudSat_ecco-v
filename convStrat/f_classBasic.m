function classBasic=f_classBasic(conv,stratMixed,mixedConv,melt)
% Find basic classification and merge suitable data
classBasic=nan(size(conv));

%% Handle mixed
% Make mixed+conv mask
maskMixedOrig=conv>=stratMixed;

% Check if all stratiform
if max(max(maskMixedOrig))==0;
    classBasic(~isnan(conv))=1;
    return
end

% Remove areas that are small
maskMixed=bwareaopen(maskMixedOrig,5);

% Remove data from within too small areas so the convective mask won't pick
% them up
conv(maskMixed==0 & maskMixedOrig==1)=0;

% Check if most of region is below melting layer and below large
% stratiform, i.e., check if it is in rain

mixedAreas=bwconncomp(maskMixed);

for ii=1:mixedAreas.NumObjects
    pixInds=mixedAreas.PixelIdxList{ii};
    meltArea=melt(pixInds);
    belowFrac=sum(meltArea<20)./length(pixInds);
    if belowFrac>0.8
        thisMat=zeros(size(conv));
        thisMat(pixInds)=1;
        [r c]=ind2sub(size(conv),pixInds);
        ucols=unique(c);
        convCols=conv(:,ucols);
        meltCols=melt(:,ucols);
        thisCols=thisMat(:,ucols);
        % We turn them upside down because that is what the original code
        % does
        convCols=flipud(convCols);
        meltCols=flipud(meltCols);
        thisCols=flipud(thisCols);

        checkCols=nan(size(meltCols));
        for jj=1:length(ucols);
            meltCol=meltCols(:,jj);
            checkCol=convCols(:,jj);
            meltCol(1:min(find(~isnan(meltCol))))=10;
            firstInd=min(find(meltCol>=20));
            checkCol(1:firstInd)=1;
            lastInd=min(find(isnan(checkCol)))-1;
            if isempty(lastInd)
                lastInd=length(meltCol);
            end
            % Remove not connected
            lastMixed=max(find(thisCols(:,jj)==1));
            testConv=convCols(:,jj);
            testConv(1:lastMixed)=1;
            firstNan=min(find(isnan(testConv)));
            if lastInd>firstNan
                lastInd=max([firstInd,firstNan]);
            end
            checkCols(firstInd:lastInd,jj)=convCols(firstInd:lastInd,jj);            
        end
        stratPerc=length(find(checkCols<stratMixed))/length(find(~isnan(checkCols)));

        medThick=median(sum(~isnan(checkCols),1)); % Median thickness of above melting layer area

        if stratPerc>0.8 & medThick>20
            maskMixed(pixInds)=0;
            conv(pixInds)=0;
        end
    end
end

% Set up check
horLarge=imdilate(maskMixed, strel('line', 50,0));%100

% Enlarge mixedConv
mixedLarge1=imdilate(maskMixed, strel('disk', 15)); %25
mixedLarge=imclose(mixedLarge1,strel('disk', 30)); %50
mixedLarge(isnan(conv))=0;
mixedLarge=imfill(mixedLarge,'holes');
mixedLarge=imerode(mixedLarge,strel('disk', 3));

% Make sure we don't enlarge into unconnected areas
mixedRays=find(any(mixedLarge==1,1));
for ii=1:length(mixedRays)
    mixedCol=mixedLarge(:,mixedRays(ii));
    mixedHorCol=horLarge(:,mixedRays(ii));
    rayPieces=bwconncomp(mixedCol);
    if rayPieces.NumObjects>1
        for jj=1:rayPieces.NumObjects
            if ~any(mixedHorCol(rayPieces.PixelIdxList{jj})==1)
                mixedCol(rayPieces.PixelIdxList{jj})=0;
            end
        end
        mixedLarge(:,mixedRays(ii))=mixedCol;
    end
end

mixedLarge=imdilate(mixedLarge,strel('disk', 3));

%% Handle conv
% Make conv mask
maskConv=conv>=mixedConv;

% Remove areas that are small
%maskConv=bwareaopen(maskConv,100);

% Set up check
horLarge2=imdilate(maskConv, strel('line', 50,0));%100

% Enlarge conv
convLarge1=imdilate(maskConv, strel('disk', 5)); %15
convLarge=imclose(convLarge1,strel('disk', 10));%50
convLarge(isnan(conv))=0;
convLarge=imfill(convLarge,'holes');
convLarge=imerode(convLarge,strel('disk', 3));

% Make sure we don't enlarge into unconnected areas
convRays=find(any(convLarge==1,1));
for ii=1:length(convRays)
    convCol=convLarge(:,convRays(ii));
    convHorCol=horLarge2(:,convRays(ii));
    rayPieces2=bwconncomp(convCol);
    if rayPieces2.NumObjects>1
        for jj=1:rayPieces2.NumObjects
            if ~any(convHorCol(rayPieces2.PixelIdxList{jj})==1)
                convCol(rayPieces2.PixelIdxList{jj})=0;
            end
        end
        convLarge(:,convRays(ii))=convCol;
    end
end

convLarge=imdilate(convLarge,strel('disk', 3));

classBasic(convLarge)=3; % Convective
classBasic(isnan(classBasic) & mixedLarge)=2; % Mixed
classBasic(isnan(classBasic))=1; % Stratiform
classBasic(isnan(conv))=nan;

end

