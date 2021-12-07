function classSub=f_classSub(classIn,asl,topo,melt,temp)

% 14 strat low
% 16 strat mid
% 18 strat high
% 25 mixed
% 30 conv
% 32 conv elevated
% 34 conv shallow
% 36 conv mid
% 38 conv deep

classSub=nan(size(classIn));

% Loop through convective areas and check if they are near the surface
convMask=zeros(size(classIn));
convMask(classIn==3)=1;

convAreas=bwconncomp(convMask);

% Calculate distance between asl and topo
distAslTopo=asl-topo;

for ii=1:convAreas.NumObjects
    % Check if next to aircraft

    [row1 col1]=ind2sub(size(classIn),convAreas.PixelIdxList{ii});
    planePix=sum(row1==18); % Number of next to plane pixels

    % Check if near surface
    aslArea=distAslTopo(convAreas.PixelIdxList{ii});
    nearSurfPix=sum(aslArea<500);
    % Shallow, mid, or deep
    meltMax=max(melt(convAreas.PixelIdxList{ii}));
    if meltMax<20 % Below melting layer: shallow
        % Near plane check
        if planePix>10
            classSub(convAreas.PixelIdxList{ii})=30;
        else
            classSub(convAreas.PixelIdxList{ii})=34;
        end
    else
        minTemp=min(temp(convAreas.PixelIdxList{ii}));
        if minTemp>=-25 % Below divergence level: mid
            classSub(convAreas.PixelIdxList{ii})=36;
        else % Deep
            classSub(convAreas.PixelIdxList{ii})=38;

        end
    end
end

% Stratiform
% Low
classSub(classIn==1 & melt<20)=14;
% Mid
classSub(classIn==1 & melt>=20 & temp>=-25)=16;
% High
classSub(classIn==1 & melt>=20 & temp<-25)=18;

% Mixed
classSub(classIn==2)=25;
end

