dataFolder = '/Users/rothzn/rwdRapidData/';
% datafiles = {'181012_stim07','181012_stim07','181012_stim08'};
samplerate=500;

baselineNorm = 0;
baselineTime1 = 2100;
baselineTime2 = 100;



% AleFis 2019-01-22
% big phasic, no tonic
rwdLevels=2;
prefix = '190122_stim';
runs=[20 22:26]';
ymin=500;
ymax=1200;

% % SteGie 2019-01-23
% % unclear phasic response
% rwdLevels=2;
% prefix = '190123_stim';
% runs=[3:8]';
% ymin=600;
% ymax=1800;

% % BraRoa 2019-01-24
% % Nice tonic difference, slight (?) phasic response
% rwdLevels=2;
% prefix = '190124_stim';
% runs=[15:22]';
% ymin=600;
% ymax=1100;

 
% % CelMar 2019-01-25
% % nice phasic
% rwdLevels=2;
% prefix = '190125_stim';
% runs=[2:7]';
% ymin=1000;
% ymax=1700;

% % OluAbi 2019-01-25
% % BAD pupil data
% rwdLevels=2;
% prefix = '190125_stim';
% runs=[11:18]';
% ymin=600;
% ymax=1700;

% % MohMat 2019-01-25
% % zero tonic, but mixed phasic, especially post response
% rwdLevels=2;
% prefix = '190125_stim';
% runs=[20:25]';
% ymin=500;
% ymax=1500;

% % SopRam 2019-01-25
% % clear phasic difference, tonic might be opposite
% rwdLevels=2;
% prefix = '190125_stim';
% runs=[27:34]';
% ymin=800;
% ymax=1700;

% % HanLoo 2019-01-28
% % noisy. possible phasic response, no tonic
% rwdLevels=2;
% prefix = '190128_stim';
% runs=[2:9]';
% ymin=1000;
% ymax=1900;

% % MitPat 2019-01-28
% % very noisy. if anything, opposite effects
% rwdLevels=2;
% prefix = '190128_stim';
% runs=[12:15]';
% ymin=1000;
% ymax=2800;

% % HilUna 2019-01-29
% % slight phasic, no tonic
% rwdLevels=2;
% prefix = '190129_stim';
% runs=[11:18]';
% yminSub=500;
% ymaxSub=1500;

% % KieNew 2019-01-29
% % nothing
% rwdLevels=2;
% prefix = '190129_stim';
% runs=[20:27]';
% yminSub=500;
% ymaxSub=1500;

% % EyaKes 2019-01-30
% % clear tonic and very weak phasic
% rwdLevels=2;
% prefix = '190130_stim';
% runs=[11:16]';

% % VirOlu 2019-01-31
% % bad subject
% rwdLevels=2;
% prefix = '190131_stim';
% runs=[2:9]';
% ymin=1000;
% ymax=2200;

% % CarCuo 2019-01-31
% % nice phasic response
% rwdLevels=2;
% prefix = '190131_stim';
% runs=[13:18]';
% ymin=1000;
% ymax=1800;

% % TucCos 2019-02-01
% % tiny opposite tonic and phasic effect
% rwdLevels=2;
% prefix = '190201_stim';
% runs=[69:72]';
% ymin=1000;
% ymax=1800;

% % FabLla 2019-02-04
% % slight tonic and phasic
% rwdLevels=2;
% prefix = '190204_stim';
% runs=[12:17]';
% ymin=1100;
% ymax=2400;

% % AleElc 2019-02-04
% % tiny phasic
% rwdLevels=2;
% prefix = '190204_stim';
% runs=[19:24]';
% ymin=500;
% ymax=1400;

% % CarMay 2019-02-06
% % nice tonic and phasic
% rwdLevels=2;
% prefix = '190206_stim';
% runs=[22:29]';
% ymin=800;
% ymax=1800;

% % BriMar 2019-02-07
% % opposite phasic
% rwdLevels=2;
% prefix = '190207_stim';
% runs=[24:29]';
% ymin=800;
% ymax=1800;

% % SteGie 2019-02-08
% % no response
% rwdLevels=2;
% prefix = '190208_stim';
% runs=[2:7]';
% ymin=800;
% ymax=1800;

% % JesSti 2019-02-08
% % Nice tonic and phasic
% rwdLevels=2;
% prefix = '190208_stim';
% runs=[10:15]';
% ymin=800;
% ymax=1800;

% % MarNem 2019-02-12
% % no tonic, maybe tiny phasic
% rwdLevels=2;
% prefix = '190212_stim';
% runs=[17:22]';
% ymin=1400;
% ymax=2200;

% % RoxNem 2019-02-12
% % tiny tonic effect
% rwdLevels=2;
% prefix = '190212_stim';
% runs=[10:13]';
% ymin=1300;
% ymax=2300;

% 
% % ClaDry 2019-02-13
% % tiny but clear phasic
% rwdLevels=2;
% prefix = '190213_stim';
% runs=[2:9]';
% ymin=900;
% ymax=1200;

% % EliAli 2019-02-15
% % small tonic, no phasic
% rwdLevels=2;
% prefix = '190215_stim';
% runs=[8:15]';
% ymin=1100;
% ymax=1900;

% % GilBin 2019-02-15
% % nothing (maybe slight phasic, and slight opposite tonic?)
% rwdLevels=2;
% prefix = '190215_stim';
% runs=[17:24]';
% ymin=600;
% ymax=1100;

% % LieCha 2019-02-22
% % small tonic, tiny phasic
% rwdLevels=2;
% prefix = '190222_stim';
% runs=[2:7]';
% ymin=1300;
% ymax=1900;

% % ChaFea 2019-02-27
% % tiny phasic, perhaps opposite tonic
% rwdLevels=2;
% prefix = '190227_stim';
% runs=[13:18]';
% ymin=1100;
% ymax=1800;

% % SamOnu 2019-02-27
% % no phasic or tonic. low reward noisy.
% rwdLevels=2;
% prefix = '190227_stim';
% runs=[4:6 8:10]';
% ymin=900;
% ymax=1500;

% % GreGre 2019-03-01
% % no phasic, maybe opposite tonic
% rwdLevels=2;
% prefix = '190301_stim';
% runs=[4:9]';
% ymin=900;
% ymax=1500;

% % CarDun 2019-03-01
% % small phasic, opposite tonic
% rwdLevels=2;
% prefix = '190306_stim';
% runs=[4:11]';
% ymin=700;
% ymax=1550;

% AriMck 2019-03-21
% tonic effect, no phasic
rwdLevels=2;
prefix = '190321_stim';
runs=[31:38]';
ymin=700;
ymax=1550;

% LirTha 2019-03-21
% clear tonic, slight phasic
rwdLevels=2;
prefix = '190322_stim';
runs=[40:47]';
ymin=700;
ymax=1550;

% ChaHar 2019-03-25
% opposite phasic!
rwdLevels=2;
prefix = '190325_stim';
runs=[28 30 32:36]';
ymin=700;
ymax=1550;

% DanSch 2019-03-25
% early phasic, no tonic
rwdLevels=2;
prefix = '190325_stim';
runs=[16:20 22:24]';
ymin=700;
ymax=1550;

% ChrUel 2019-03-26
% unclear phasic, no tonic
rwdLevels=2;
prefix = '190326_stim';
runs=[38:47]';
ymin=700;
ymax=1550;

% GabFre 2019-03-27
% nice phasic and tonic
rwdLevels=2;
prefix = '190327_stim';
runs=[49:52 54:57]';
ymin=1100;
ymax=2000;

% SheYu 2019-03-29
% tiny phasic, no tonic
rwdLevels=2;
prefix = '190329_stim';
runs=[17:22]';
ymin=700;
ymax=1550;

% AngFal 2019-04-01
% noisy, unclear phasic
rwdLevels=2;
prefix = '190401_stim';
runs=[31:36]';
ymin=1200;
ymax=1900;

% ChaObe 2019-04-02
% small but clear phasic 
rwdLevels=2;
prefix = '190402_stim';
runs=[36:45]';
ymin=700;
ymax=1500;

% PriDel 2019-04-04
% tiny tonic, small unclear phasic
rwdLevels=2;
prefix = '190404_stim';
runs=[19:24]';
ymin=700;
ymax=1500;

% NatKes 2019-04-04
% 
rwdLevels=2;
prefix = '190405_stim';
runs=[28:33]';
ymin=700;
ymax=1500;


runsString = num2str(runs,'%02d');
rows=ceil(length(runs)/rwdLevels);
cols=rwdLevels;
xmin=0;
xmax=9000;

clear e s x
for r=1:length(runs)
    stimfile = [dataFolder prefix runsString(r,:) '.mat'];
    e{r} = myGetTaskEyeTraces(stimfile, 'removeBlink=3');
    s{r} = load(stimfile);
    x{r} = getTaskParameters(stimfile);
    if isfield(s{r}.myscreen,'stimulus')
        rewardType{r} = s{r}.myscreen.stimulus.rewardType;
    else
        if mod(r,2)==0
            rewardType{r} = 'H';
        else
            rewardType{r} = 'L';
        end
    end
end
%%
close all
clear RT nullTrials correctness correctResponse response runSize rwdRuns rwdLevel



for r=1:length(runs)
    if rewardType{r}=='H'
            rwdLevel(r) = 1;
        else
            rwdLevel(r) = 2;
    end
    
    %     e{r}.eye= e{r}.eyes{1};
    if baselineNorm
        baseline = nanmean(e{r}.eye.pupil(:,end-baselineTime1:end-baselineTime2)')';
        e{r}.eye.pupil(2:end,:) = e{r}.eye.pupil(2:end,:) - baseline(1:end-1);
    end
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
    %     stimTime = s{r}.fixStimulus.stimTime * samplerate;
    cueTime = s{r}.fixStimulus.cueTime * samplerate;
    interTime = s{r}.fixStimulus.interTime * samplerate;
    responseTime = s{r}.fixStimulus.responseTime * samplerate;
    %     allTimes = [stimTime interTime stimTime interTime responseTime];
    allTimes = [cueTime interTime cueTime responseTime];
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
xmax=4000;
figure(2)
rows=rwdLevels+1;
cols=2;
clear correctnessRwd correctRespRwd responseRwd rtRwd nullRwd
for rwd=1:rwdLevels
    rwdRuns = find(rwdLevel(:)==rwd);
    numRuns(rwd) = length(rwdRuns);
    numTrials(rwd) = sum(runSize(rwdRuns,1));
    trialLength(rwd) = max(runSize(rwdRuns,2));

    rwdPupil{rwd} = NaN(numTrials(rwd),trialLength(rwd));
    trialCounter=0;
    for r=rwdRuns(1):rwdRuns(end)
        rwdPupil{rwd}(trialCounter+1:trialCounter+runSize(r,1),1:runSize(r,2)) = e{r}.eye.pupil;
        trialCounter = trialCounter + runSize(r,1);
    end
    subplot(rows,cols,rwd)
    meanPupil = nanmean(rwdPupil{rwd})';
    plot(meanPupil)
    axis([xmin xmax ymin ymax]);
    hold all
    %     stimTime = s{1}.fixStimulus.stimTime * samplerate;
    %     interTime = s{1}.fixStimulus.interTime * samplerate;
    %     responseTime = s{1}.fixStimulus.responseTime * samplerate;
%     allTimes = [stimTime interTime stimTime interTime responseTime];
    
    
    stdRwd = std(rwdPupil{rwd},0,1,'omitnan');
    stdRwd(isnan(stdRwd)) = ones;
    
    meanPupil(isnan(meanPupil)) = zeros;
    
    dsErrorsurface(1:length(stdRwd), meanPupil, stdRwd, [0.8 0.8 0.8], 0.2);
    plotRwdStimSegments(meanPupil, allTimes)
    if rwd==1
        title('high');
    else
        title('low');
    end
%     temp = correctness(rwd:rwdLevels:end,:)';
%     correctnessRwd(rwd,:) = temp(:)';
%     temp = correctResponse(rwd:rwdLevels:end,:)';
%     correctRespRwd(rwd,:) = temp(:)';
%     temp = response(rwd:rwdLevels:end,:)';
%     responseRwd(rwd,:) = temp(:)';
%     temp = RT(rwd:rwdLevels:end,:)';
%     rtRwd(rwd,:) = temp(:)';
%     temp=nullTrials(rwd:rwdLevels:end,:)';
%     nullRwd(rwd,:) = temp(:)';
end
subplot(rows,cols,rwdLevels+1)
for rwd=1:rwdLevels
    plot(nanmean(rwdPupil{rwd})');
    hold all
    axis([xmin xmax ymin ymax]);
end

if rwdLevels>1
    legend('high','low');
    subplot(rows,cols,rwdLevels+2)
    rwd1 = nanmean(rwdPupil{1});
    rwd2 = nanmean(rwdPupil{2});
    l = min(length(rwd1),length(rwd2));
    meanDiff = rwd1(1:l) - rwd2(1:l);
%     tempDiff = rwdPupil{1}(:,1:l) - rwdPupil{2}(:,1:l);
%     stdDiff = std(tempDiff,0,1,'omitnan');
%     stdDiff(isnan(stdDiff)) = ones;
%     stdDiff(stdDiff==0) = 0.001;
    plot(meanDiff');
    hold on
    plot(zeros(l,1),'--','Color',[0.5,0.5,0.5]);
    
%     dsErrorsurface(1:l, meanDiff, stdDiff, [0.8 0.8 0.8], 0.2);
    plotRwdStimSegments(meanDiff, allTimes)
    xlim([xmin xmax]);
    title('high minus low');
    
    %plot again, wider view of entire trial length
    subplot(rows,cols,rwdLevels+3)
    for rwd=1:rwdLevels
        plot(nanmean(rwdPupil{rwd})');
        hold all
        ylim([ymin ymax]);
    end
    legend('high','low');
    subplot(rows,cols,rwdLevels+4)
    rwd1 = nanmean(rwdPupil{1});
    rwd2 = nanmean(rwdPupil{2});
    l = min(length(rwd1),length(rwd2));
    meanDiff = rwd1(1:l) - rwd2(1:l);
%     tempDiff = rwdPupil{1}(:,1:l) - rwdPupil{2}(:,1:l);
%     stdDiff = std(tempDiff,0,1,'omitnan');
%     stdDiff(isnan(stdDiff)) = ones;
%     stdDiff(stdDiff==0) = 0.001;
    plot(meanDiff');
    hold on
    plot(zeros(l,1),'--','Color',[0.5,0.5,0.5]);
    
%     dsErrorsurface(1:l, meanDiff, stdDiff, [0.8 0.8 0.8], 0.2);
    plotRwdStimSegments(meanDiff, allTimes)
    title('high minus low');
    
end
set(gcf,'position',[100 100 1000 700])

% %% correct vs. incorrect trials
% xmax=3000;
% figure(3)
% rows=2;
% cols=2;
% for rwd=1:rwdLevels
%     rwdTitle = 'high - ';
%     if rwd==2
%         rwdTitle = 'low - ';
%     end
%     for c =  1:2
%         if c==1
%             subplotTitle = [rwdTitle 'incorrect'];
%         else
%             subplotTitle = [rwdTitle 'correct'];
%         end
%         subplot(rows,cols,(rwd-1)*2+c)
%         rwdCorrectPupil = rwdPupil{rwd}(correctnessRwd(rwd,:)==(c*2-3),:);
%         meanPupil = nanmean(rwdCorrectPupil);
%         plot(meanPupil)
%         axis([xmin xmax ymin ymax]);
%         hold all
%         
% 
%         stdRwd = std(rwdCorrectPupil,0,1,'omitnan');
%         stdRwd(isnan(stdRwd)) = ones;
%         meanPupil(isnan(meanPupil)) = zeros;
% 
%         dsErrorsurface(1:length(stdRwd), meanPupil, stdRwd, [0.8 0.8 0.8], 0.2);
%         plotRwdStimSegments(meanPupil, allTimes)
%         title(subplotTitle);
%         axis([xmin xmax ymin ymax]);
% %         if c==1 & rwd==2
% %            keyboard
% %         end
%    end
% end
% 

% %% first dimming vs. second dimming
% figure(4)
% rows=2;
% cols=2;
% for rwd=1:rwdLevels
%     rwdTitle = 'high - ';
%     if rwd==2
%         rwdTitle = 'low - ';
%     end
%     for c =  1:2
%         if c==1
%             subplotTitle = [rwdTitle 'first'];
%         else
%             subplotTitle = [rwdTitle 'second'];
%         end
%         subplot(rows,cols,(rwd-1)*2+c)
%         correctRespPupil = rwdPupil{rwd}(correctRespRwd(rwd,:)==c,:);
%         meanPupil = nanmean(correctRespPupil);
%         plot(meanPupil)
%         axis([xmin xmax ymin ymax]);
%         hold all
%         
%         
%         stdRwd = std(correctRespPupil,0,1,'omitnan');
%         stdRwd(isnan(stdRwd)) = ones;
%         meanPupil(isnan(meanPupil)) = zeros;
%         
%         dsErrorsurface(1:length(stdRwd), meanPupil, stdRwd, [0.8 0.8 0.8], 0.2);
%         plotRwdStimSegments(meanPupil, allTimes)
%         title(subplotTitle);
%         axis([xmin xmax ymin ymax]);
%     end
% end

% %% stimulated vs. null trials
% figure(5)
% rows=2;
% cols=2;
% for rwd=1:rwdLevels
%     rwdTitle = 'high - ';
%     if rwd==2
%         rwdTitle = 'low - ';
%     end
%     for c =  1:2
%         if c==1
%             subplotTitle = [rwdTitle 'stim'];
%         else
%             subplotTitle = [rwdTitle 'null'];
%         end
%         subplot(rows,cols,(rwd-1)*2+c)
%         
%         nullPupil = rwdPupil{rwd}(nullRwd(rwd,:)==c-1,:);
%         meanPupil = nanmean(nullPupil);
%         plot(meanPupil)
%         axis([xmin xmax ymin ymax]);
%         hold all
%         
%         
%         stdRwd = std(nullPupil,0,1,'omitnan');
%         stdRwd(isnan(stdRwd)) = ones;
%         meanPupil(isnan(meanPupil)) = zeros;
%         
%         dsErrorsurface(1:length(stdRwd), meanPupil, stdRwd, [0.8 0.8 0.8], 0.2);
%         plotRwdStimSegments(meanPupil, allTimes)
%         title(subplotTitle);
%         axis([xmin xmax ymin ymax]);
%     end
% end

%%
% figure(6)
% rows=2;
% cols=2;
% clear responsePupil
% for rwd=1:rwdLevels
%     rtTitle = 'high - RT';
%     if rwd==2
%         rtTitle = 'low - RT';
%     end
% 
%     subplot(rows,cols,rwd)
%     for trial=1:numTrials(rwd)
%         if ~isnan(rtRwd(rwd,trial))
%             startSample = ceil(rtRwd(rwd,trial)*500);
%             responsePupil{rwd}(trial,:) = rwdPupil{rwd}(trial,startSample:6000+startSample);
%         else
%             responsePupil{rwd}(trial,:) = NaN;
%         end
%     end
%     meanPupil = nanmean(responsePupil{rwd});
%     plot(meanPupil)
%     axis([xmin xmax ymin ymax]);
%     hold all
%     
% 
%     stdRwd = std(responsePupil{rwd},0,1,'omitnan');
%     stdRwd(isnan(stdRwd)) = ones;
%     meanPupil(isnan(meanPupil)) = zeros;
% 
%     dsErrorsurface(1:length(stdRwd), meanPupil, stdRwd, [0.8 0.8 0.8], 0.2);
%     plotRwdStimSegments(meanPupil, allTimes)
%     title(rtTitle);
%     axis([xmin xmax ymin ymax]);
% 
% end

%%
figure(7)
for r=1:length(runs)
    plot(e{r}.eye.pupil')
    hold all
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plotRwdStimSegments(meanPupil, allTimes)

itime = 1;
startT=1;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0 1 1]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0.1 0.7 0.6]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[0 1 1]);
itime = itime+1;
startT=endT;
endT = startT+allTimes(itime);
plot(startT:endT-1,meanPupil(startT:endT-1),'lineWidth',4,'color',[1 1 0.0]);

end