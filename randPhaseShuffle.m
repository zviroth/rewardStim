% randPhaseShuffle
% 
% Creates a phase-randomized version of an image. Image can be 2D (height x
% width) or 3D (height x width x rgb). 
%
%   usage: [shuf, ind] = randPhaseShuffle(orig);
%   by: zvi roth
%   date: 10/22/2018
%   purpose: creates a phase-randomized version of an image
%
%   NOTE: the Fourier amplitude of shuf will not be the same as the
%   spectrum of orig because of rounding error. This stems from the 
%   intermediate result, ifft2(ifftshift(F_shuffle)), including tiny
%   complex values. The Fourier amplitude of that intermediate result
%   is identical to the original spectrum.
%

function [shuf, ind] = randPhaseShuffle(orig)
shuf = uint8(zeros(size(orig)));
    if size(orig,3)>1 %shuffle phases randomly for first channel
        [shuf(:,:,1), ind] = randPhaseShuffleBW(orig(:,:,1));
        for ch=2:size(orig,3) %use indices from first channel for other channels
            [shuf(:,:,ch), ind] = randPhaseShuffleBW(orig(:,:,ch),ind);
        end
    else%single channel
        [shuf, ind] = randPhaseShuffleBW(orig);
    end
end

function [shuf, ind] = randPhaseShuffleBW(orig,ind)
    F = fft2(orig);
    F = fftshift(F);
    Fph = angle(F);
    Famp = abs(F);
    if exist('ind','var')%random shuffle
        Fph_shuffle = Fph(ind);
    else %we already have the indices to shuffle by
        [Fph_shuffle, ind] = phaseShuffle(Fph);
    end
    F_shuffle = Famp.*exp(1i*Fph_shuffle);
    shuf=abs(ifft2(ifftshift(F_shuffle)));
end
   

function [randPhase, ind] = phaseShuffle(orig)
    sz = size(orig);
    randPhase = zeros(sz);
    ind = zeros(sz);
    trueInd = 1:numel(orig);
    trueInd = reshape(trueInd,sz(1),sz(2));%the original indices
    if mod(sz(1),2)==0 % even rows 
        if mod(sz(2),2)==0 % even columns and rows
            [randPhase(1,1), ind(1,1)] = phaseShuffle(orig(1,1));%first element
            [randPhase(2:end,1), tempInd] = phaseShuffle(orig(2:end,1));%first column, excluding first element
            colInd = trueInd(2:end,1);
            ind(2:end,1) = colInd(tempInd);%tempInd are relative indices, only within colInd. Here we convert them to the absolute original indices
            [randPhase(1,2:end), tempInd] = phaseShuffle(orig(1,2:end));%first row, excluding first element
            rowInd = trueInd(1,2:end);
            ind(1,2:end) = rowInd(tempInd);
            [randPhase(2:end,2:end), tempInd] = phaseShuffle(orig(2:end,2:end));%everything else
            restInd = trueInd(2:end,2:end);
            ind(2:end,2:end) = restInd(tempInd);
            
        else %only even rows, not columns
            [randPhase(1,:), tempInd] = phaseShuffle(orig(1,:));%first row
            rowInd = trueInd(1,:);
            ind(1,:) = rowInd(tempInd);
            [randPhase(2:end,:), tempInd] = phaseShuffle(orig(2:end,:));%everything excluding first row
            restInd = trueInd(2:end,:);
            ind(2:end,:) = restInd(tempInd);
            
        end
        return
    elseif mod(sz(2),2)==0 % even columns but not rows
        [randPhase(:,1), tempInd] = phaseShuffle(orig(:,1));%first column
        colInd = trueInd(:,1);
        ind(:,1) = colInd(tempInd);
        [randPhase(:,2:end), tempInd] = phaseShuffle(orig(:,2:end));% everything else
        restInd = trueInd(:,2:end);
        ind(:,2:end) = restInd(tempInd);
        return
    end
    
    %orig must now have odd rows and odd columns
    lfull = numel(orig);
    lhalf = floor(lfull/2);
    ind(1:lhalf) = randperm(lhalf);%shuffled indiced for left half
    ind(lhalf+1) = lhalf+1;%central element does not move
    ind(lhalf+2:lfull) = lfull + 1 - ind(lhalf:-1:1);%shuffled indices for right half
    randPhase = orig(ind);
    return
end

