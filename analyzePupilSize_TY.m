% dataFolder = '/Users/rothzn/rewardData/';
dataFolder = '/Users/rothzn/data/reward/';
% datafiles = {'181012_stim07','181012_stim07','181012_stim08'};
samplerate=500;

% runs = length(datafiles);
prefix = '190104_stim';%HanLoo
runs = [14:19]';

prefix = '181012_stim';%Tenzin
runs = [6:8]';

prefix = '190108_stim';%Tenzin
runs = [85:92]';

% prefix = '181018_stim';%me, tired
% runs=[3:5]';

runsString = num2str(runs,'%02d');
rows=8;
cols=1;
xmin=0;
xmax=9000;
ymin=800;
ymax=1700;
% for r=1:length(runs)
%     stimfile = [dataFolder prefix runsString(r,:) '.mat'];
%     e{r} = getTaskEyeTraces(stimfile, 'removeBlink=3');
%     s{r} = load(stimfile);
%     x{r} = getTaskParameters(stimfile);
% end
%%
for r=1:length(runs)
    RT(r,:) = x{r}{1}.reactionTime;
    nullTrials(r,:) = x{r}{2}.randVars.nullTrial;
    correctness(r,:) = s{r}.task{1}{1}.correctness;
    correctResponse(r,:) = s{r}.task{1}{1}.correctResponse;
    response(r,:) = s{r}.task{1}{1}.response;
    runSize(r,:) = size(e{r}.eye.pupil);
    meanPupil = nanmean(e{r}.eye.pupil)';
    figure(1)
    subplot(rows,cols,r)
    plot(meanPupil)
    axis([xmin xmax ymin ymax]);
    hold all
    stimTime = s{r}.fixStimulus.stimTime * samplerate;
    interTime = s{r}.fixStimulus.interTime * samplerate;
    responseTime = s{r}.fixStimulus.responseTime * samplerate;
    allTimes = [stimTime interTime stimTime interTime responseTime];
    plotRwdStimSegments(meanPupil, allTimes)
    title(['run ' num2str(r)]);
%     %%
%     figure
%     plot(e.eye.pupil')
%     hold on
%     plot(meanPupil,'color','k','linewidth',3);
%     %%
%     mean1 = nanmean(e.eye.pupil(s.task{1}{1}.correctResponse==1,:));
%     std1 = std(e.eye.pupil(s.task{1}{1}.correctResponse==1,:),0,1,'omitnan');
%     sem1 = std1./sqrt(size(e.eye.pupil(s.task{1}{1}.correctResponse==1,:),1));
%     mean2 = nanmean(e.eye.pupil(s.task{1}{1}.correctResponse==2,:));
%     std2 = std(e.eye.pupil(s.task{1}{1}.correctResponse==2,:),0,1,'omitnan');
%     ste2 = std2./sqrt(size(e.eye.pupil(s.task{1}{1}.correctResponse==2,:),1));
%     figure(2)
%     plot(mean1);
%     hold all
%     plot(mean2)
end
%%
figure(2)
rows=1;
cols=1;
for rwd=1%:2
    numRuns(rwd) = length(1:length(runs));
    numTrials(rwd) = sum(runSize(1:end,1));
    trialLength(rwd) = max(runSize(1:end,2));
    rwdPupil{rwd} = NaN(numTrials(rwd),trialLength(rwd));
    trialCounter=0;
    for r=rwd:2:length(runs)
        rwdPupil{rwd}(trialCounter+1:trialCounter+runSize(r,1),1:runSize(r,2)) = e{r}.eye.pupil;
    end
    subplot(rows,cols,rwd)
    meanPupil = nanmean(rwdPupil{rwd})';
    plot(meanPupil)
    axis([xmin xmax ymin ymax]);
    hold all
%     stimTime = s{1}.fixStimulus.stimTime * samplerate;
%     interTime = s{1}.fixStimulus.interTime * samplerate;
%     responseTime = s{1}.fixStimulus.responseTime * samplerate;
    allTimes = [stimTime interTime stimTime interTime responseTime];
    plotRwdStimSegments(meanPupil, allTimes)
    
    stdRwd = std(rwdPupil{rwd},0,1,'omitnan');
    stdRwd(isnan(stdRwd)) = ones;
    
    meanPupil(isnan(meanPupil)) = zeros;
    
    dsErrorsurface(1:length(stdRwd), meanPupil, stdRwd, [0.8 0.8 0.8], 0.2);
    if rwd==1
       title('high'); 
    else
        title('low');
    end
    temp = correctness(1:end,:)';
    correctnessRwd(rwd,:) = temp(:)';
    temp = correctResponse(1:end,:)';
    correctRespRwd(rwd,:) = temp(:)';
    temp = response(1:end,:)';
    responseRwd(rwd,:) = temp(:)';
    temp = RT(1:end,:)';
    rtRwd(rwd,:) = temp(:)';
    temp=nullTrials(1:end,:)';
    nullRwd(rwd,:) = temp(:)';
end

%%
figure(3)
rows=1;
cols=2;
for rwd=1%:2
    rwdTitle = 'high - ';
    if rwd==2
        rwdTitle = 'low - ';
    end
    for c =  1:2
        if c==1
            subplotTitle = [rwdTitle 'incorrect'];
        else
            subplotTitle = [rwdTitle 'correct'];
        end
        subplot(rows,cols,(rwd-1)*2+c)
        rwdCorrectPupil = rwdPupil{rwd}(correctnessRwd(rwd,:)==(c*2-3),:);
        meanPupil = nanmean(rwdCorrectPupil);
        plot(meanPupil)
        axis([xmin xmax ymin ymax]);
        hold all
        plotRwdStimSegments(meanPupil, allTimes)
        
        stdRwd = std(rwdCorrectPupil,0,1,'omitnan');
        stdRwd(isnan(stdRwd)) = ones;
        meanPupil(isnan(meanPupil)) = zeros;
        
        dsErrorsurface(1:length(stdRwd), meanPupil, stdRwd, [0.8 0.8 0.8], 0.2);
        title(subplotTitle);
        axis([xmin 2000 ymin ymax]);
   end
end


%%
figure(4)
rows=1;
cols=2;
for rwd=1%:2
    rwdTitle = 'high - ';
    if rwd==2
        rwdTitle = 'low - ';
    end
    for c =  1:2
        if c==1
            subplotTitle = [rwdTitle 'first'];
        else
            subplotTitle = [rwdTitle 'second'];
        end
        subplot(rows,cols,(rwd-1)*2+c)
        correctRespPupil = rwdPupil{rwd}(correctRespRwd(rwd,:)==c,:);
        meanPupil = nanmean(correctRespPupil);
        plot(meanPupil)
        axis([xmin xmax ymin ymax]);
        hold all
        plotRwdStimSegments(meanPupil, allTimes)
        
        stdRwd = std(correctRespPupil,0,1,'omitnan');
        stdRwd(isnan(stdRwd)) = ones;
        meanPupil(isnan(meanPupil)) = zeros;
        
        dsErrorsurface(1:length(stdRwd), meanPupil, stdRwd, [0.8 0.8 0.8], 0.2);
        title(subplotTitle);
        axis([xmin 2000 ymin ymax]);
   end
end

%%
figure(5)
rows=1;
cols=2;
for rwd=1%:2
    rwdTitle = 'high - ';
    if rwd==2
        rwdTitle = 'low - ';
    end
    for c =  1:2
        if c==1
            subplotTitle = [rwdTitle 'stim'];
        else
            subplotTitle = [rwdTitle 'null'];
        end
        subplot(rows,cols,(rwd-1)*2+c)
        
        nullPupil = rwdPupil{rwd}(nullRwd(rwd,:)==c-1,:);
        meanPupil = nanmean(nullPupil);
        plot(meanPupil)
        axis([xmin xmax ymin ymax]);
        hold all
        plotRwdStimSegments(meanPupil, allTimes)
        
        stdRwd = std(nullPupil,0,1,'omitnan');
        stdRwd(isnan(stdRwd)) = ones;
        meanPupil(isnan(meanPupil)) = zeros;
        
        dsErrorsurface(1:length(stdRwd), meanPupil, stdRwd, [0.8 0.8 0.8], 0.2);
        title(subplotTitle);
        axis([xmin 2000 ymin ymax]);
   end
end


%%
figure(6)
for r=1:length(runs)
   plot( e{r}.eye.pupil')
   hold all
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plotRwdStimSegments(meanPupil, allTimes)

itime = 1;
startT=1;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4);

end