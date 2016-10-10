function I = DispDictionary(D,ImageType,IsOriginalData)
%Display the dictionary elements as a image
%Version 1: 10/21/2009
%Version 2: 10/28/2009
%Written by Mingyuan Zhou, Duke ECE, mz1@ee.duke.edu
if nargin<2
    ImageType = 'gray';
end
if nargin<3
    IsOriginalData = false;
end
[P,K] = size(D);
RowNum = fix(sqrt(K));
ColNum = ceil(K/RowNum);
if strcmp(ImageType,'gray')==1
    PatchSize = ceil(sqrt(P));
    D = [D;zeros(PatchSize^2-P,K)];
elseif strcmp(ImageType,'rgb')==1
    PatchSize = ceil(sqrt(P/3));
else
    disp('error');
    return;
end
I = zeros(RowNum*PatchSize + RowNum + 1, ColNum*PatchSize + ColNum + 1, 3);
I(:,:,3)=1;
for row = 1:RowNum
    for col = 1:ColNum
        k = (row-1)*ColNum + col;
        if k<=K
            if ~IsOriginalData
                D(:,k) = D(:,k) - min(D(:,k));
                D(:,k) = D(:,k)/max(D(:,k));
            end
            if strcmp(ImageType,'gray')==1
                I( (row-1)*(PatchSize+1)+2 : row*(PatchSize+1) , (col-1)*(PatchSize+1)+2 : col*(PatchSize+1), 1) = reshape(D(:,k),PatchSize,PatchSize);
                I( (row-1)*(PatchSize+1)+2 : row*(PatchSize+1) , (col-1)*(PatchSize+1)+2 : col*(PatchSize+1), 2) = reshape(D(:,k),PatchSize,PatchSize);
                I( (row-1)*(PatchSize+1)+2 : row*(PatchSize+1) , (col-1)*(PatchSize+1)+2 : col*(PatchSize+1), 3) = reshape(D(:,k),PatchSize,PatchSize);
            else
                I( (row-1)*(PatchSize+1)+2 : row*(PatchSize+1) , (col-1)*(PatchSize+1)+2 : col*(PatchSize+1), 1) = reshape(D(1:PatchSize^2,k),PatchSize,PatchSize);
                I( (row-1)*(PatchSize+1)+2 : row*(PatchSize+1) , (col-1)*(PatchSize+1)+2 : col*(PatchSize+1), 2) = reshape(D(PatchSize^2+(1:PatchSize^2),k),PatchSize,PatchSize);
                I( (row-1)*(PatchSize+1)+2 : row*(PatchSize+1) , (col-1)*(PatchSize+1)+2 : col*(PatchSize+1), 3) = reshape(D(2*PatchSize^2+(1:PatchSize^2),k),PatchSize,PatchSize);
            end
        else
            I( (row-1)*(PatchSize+1)+2 : row*(PatchSize+1) , (col-1)*(PatchSize+1)+2 : col*(PatchSize+1), 1) = 0;
            I( (row-1)*(PatchSize+1)+2 : row*(PatchSize+1) , (col-1)*(PatchSize+1)+2 : col*(PatchSize+1), 2) = 0;
            I( (row-1)*(PatchSize+1)+2 : row*(PatchSize+1) , (col-1)*(PatchSize+1)+2 : col*(PatchSize+1), 3) = 0;
        end
    end
end
% imshow(I,[]);
imagesc(I(:,:,1))