function [IOut,output] = denoiseImageDCT(Image,sigma,K,varargin)
%==========================================================================
%Overcomplete discrete cosine dictionary DCT
% ===================================================================
% input parameters : Image - noisy image (grayscale)
%                    sigma - Gaussian random noise window
%                    K - indicates the number of atoms in the dictionary 256.
%                    varargin - list of variable length input parameters.
% Adjustable parameters : blockSize - the size of the block on which the algorithm works.
%                         errorFactor - The allowed error factor to be represented. 
%                         maxBlocksToConsider - the maximum number of blocks that can be processed.
%                         slidingFactor - The sliding distance between processed blocks. 
%                         waitBarOn - Shows the progress of the algorithm.
% output parameters :IOut - Reconstructs the image.
%                    output - String structure:
%                    D - dictionary for denoising
% =========================================================================
%% parameter setting
Reduce_DC = 1;
[NN1,NN2] = size(Image);
C = 1.15;
waitBarOn = 1;
maxBlocksToConsider = 260000;
slidingDis = 1;
bb = 8;

%% Preprocessing
for argI = 1:2:length(varargin)
    if (strcmp(varargin{argI}, 'slidingFactor'))
        slidingDis = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'errorFactor'))
        C = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'maxBlocksToConsider'))
        maxBlocksToConsider = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'blockSize'))
        bb = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'waitBarOn'))
        waitBarOn = varargin{argI+1};
    end
end
errT = C*sigma;
%% DCT Framework Initial Dictionary
Pn=ceil(sqrt(K));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); 
    end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);
while (prod(floor((size(Image)-bb)/slidingDis)+1)>maxBlocksToConsider)
    slidingDis = slidingDis+1;
end

%% Image Processing
[blocks,idx] = my_im2col(Image,[bb,bb],slidingDis);

if (waitBarOn)
    counterForWaitBar = size(blocks,2);
    h = waitbar(0,'Denoising In Process ...');
end

for jj = 1:10000:size(blocks,2)
    if (waitBarOn)
        waitbar(jj/counterForWaitBar);
    end
    jumpSize = min(jj+10000-1,size(blocks,2));
    if (Reduce_DC)
        vecOfMeans = mean(blocks(:,jj:jumpSize));
        blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);
    end

%% Optimized OMP algorithm
tic
Coefs = OMPerr(DCT,blocks(:,jj:jumpSize),errT);
    if (Reduce_DC)
        blocks(:,jj:jumpSize)= DCT*Coefs + ones(size(blocks,1),1) * vecOfMeans;
    else
        blocks(:,jj:jumpSize)= DCT*Coefs ;
    end
end
toc
time = toc

count = 1;
Weight= zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(size(Image)-bb+1,idx);
% save date1 rows;
% save date cols;

toc
time1=toc

for i  = 1:length(cols)
    col = cols(i); 
    row = rows(i);
    block =reshape(blocks(:,count),[bb,bb]);
    IMout(row:row+bb-1,col:col+bb-1)=IMout(row:row+bb-1,col:col+bb-1)+block;
    Weight(row:row+bb-1,col:col+bb-1)=Weight(row:row+bb-1,col:col+bb-1)+ones(bb);
    count = count+1;
end;

toc

time2=toc
if (waitBarOn)
    close(h);
end

%% Output
IOut = (Image+0.034*sigma*IMout)./(1+0.034*sigma*Weight);
output.D = DCT;

