%
%        $Id: rwdStim.m,v 1.2 2016/01/22 15:30:22 eli Exp $
%      usage: rwdStim('rewardType=''H''','runNum=1','useStaircase=1','currBal=30','numTrials=15','displayName=''rm315''', 'fixThresh=0.2');
%         by: zvi roth
%       date: 10/02/2018
%    purpose: measure arousal effect on stimulus responses


% getStimvolFromVarname('orientation', s.myscreen, s.task{2})
% e = getTaskParameters(s.myscreen,s.task);
% mrPrintSurf

function [] = rwdStim(varargin)

% check arguments
if ~any(nargin == [0:10])
    help otopySzSF
    return
end

% % evaluate the input arguments
getArgs(varargin, [], 'verbose=0');

% set default parameters
% if ieNotDefined('direction'),direction = -1;end

if ieNotDefined('displayName'), displayName = 'rm315'; end
if ieNotDefined('waitForBacktick')
    if strcmp(displayName, 'rm315') || strcmp(displayName, 'laptop')
        waitForBacktick = 0;
    else
        waitForBacktick = 1;
    end
end
if ieNotDefined('useStaircase'), useStaircase = 1; end
if ieNotDefined('threshStair1'), threshStair1 = 0; end
if ieNotDefined('threshStair2'), threshStair2 = 0.3; end
interTime = 0.7;
stimTime = 0.3;
responseTime=1;
stimLen = 2*interTime + 2*stimTime + responseTime;
if ieNotDefined('stimLen'),stimLen = 2*interTime + 2*stimTime + responseTime;end %should be equal to fixation trial length
%also, should be a multiple of frameLen
if ieNotDefined('trialLen'),trialLen = 18;end %in seconds
if ieNotDefined('frameLen'),frameLen = 0.25; end
if ieNotDefined('innerEdge'),innerEdge = 1.3; end
if ieNotDefined('outerEdge'),outerEdge = 25; end
if ieNotDefined('rewardType'), rewardType = 'H'; end
% if ieNotDefined('rewardValue'), rewardValue = 1.0; end %max reward per trial. replaces incr
if ieNotDefined('runNum'), runNum = 1; end
if ieNotDefined('probRwd'), probRwd = 0; end
if ieNotDefined('currBal'), currBal = 30; end
% if ieNotDefined('fixThresh'), fixThresh = 0.2; end
if ieNotDefined('numTRs'), numTRs = 204; end
TR=1.5;
if ieNotDefined('numTrials'), numTrials = ceil(TR*numTRs/trialLen); end

incrRwdL = -0.01;%reward decreases on every low reward run
incrRwdH = 0.05;%reward increases on every high reward run
initRwdL = 0.09;%reward for first low reward run
initRwdH = 1.0;%reward for first high reward run

if rewardType == 'H'
    rewardValue = initRwdH + incrRwdH * runNum/2;
elseif rewardType == 'L'
    rewardValue = initRwdL + incrRwdL * runNum/2;
    rewardValue = min(rewardValue,0.01);%don't want to reach zero
end

global stimulus;

% store parameters in stimulus variable
% update stimulus every 250 ms
stimulus.frameLen = frameLen;
% stimulus is on for 1.5 seconds
stimulus.stimLen = stimLen;

% trial is 15 seconds
stimulus.trialLen = trialLen;
% inner and outer edges
stimulus.inner = innerEdge;
stimulus.outer = outerEdge;


%reward parameters
stimulus.currBal = currBal;
stimulus.rewardType = rewardType;
stimulus.rewardValue = rewardValue;
stimulus.probRwd = probRwd;
stimulus.runNum = runNum;
% global incr;%max reward per trial
% if stimulus.rewardType == 'H'
%     incr = 1.1;
% elseif stimulus.rewardType == 'L'
%     incr = 0.01;
% end


% initalize the screen
myscreen.background = 'gray';
myscreen.autoCloseScreen = 0;
myscreen.allowpause = 1;
myscreen.saveData = 1;
myscreen.displayName = displayName;
%myscreen.displayName = '3tb';
% myscreen.displayName = '7t';
% myscreen.displayName = 'rm315';
% myscreen.displayName = 'laptop';

myscreen = initScreen(myscreen);

global fixStimulus
fixStimulus.useStaircase = useStaircase;
fixStimulus.diskSize = 0;
fixStimulus.fixWidth = 0.75;
fixStimulus.fixLineWidth = 3;
fixStimulus.stairStepSize = 0.05;
fixStimulus.threshStair1 = threshStair1;
fixStimulus.threshStair2 = threshStair2;
% fixStimulus.threshold = fixThresh;
fixStimulus.responseTime = responseTime;
fixStimulus.stimTime = stimTime;
fixStimulus.interTime = interTime;
fixStimulus.trialTime = trialLen;%same trial length for both tasks
fixStimulus.waitForBacktick = waitForBacktick;
fixStimulus.fixWidth = 1;
fixStimulus.fixLineWidth = 3;

[task{1} myscreen] = myFixStairInitTask(myscreen);
task{1}{1}.numTrials = numTrials;
% task{1}{1}.nTrials = numTrials;
task{1}{1}.response = zeros(numTrials,1);
task{1}{1}.correctResponse = zeros(numTrials,1);
task{1}{1}.correctness = zeros(numTrials,1);
%both tasks should wait for backtick to begin
task{1}{1}.waitForBacktick = waitForBacktick;
task{2}{1}.waitForBacktick = waitForBacktick;

%each segment is a different phase. Entire trial is a single orientation,
%contrast,and spatial frequency.
seglen = [stimulus.frameLen * ones(1,stimulus.stimLen/stimulus.frameLen) stimulus.trialLen-stimulus.stimLen];
seglen(end) = stimulus.trialLen-stimulus.stimLen - 0.5;%shorten blank segment length to be sure we finish before the trigger
task{2}{1}.seglen = seglen;
task{2}{1}.synchToVol = zeros(size(seglen));
if waitForBacktick
    task{2}{1}.synchToVol(end) = 1;
end

orientations = linspace(0, 180, 7);
orientations = orientations(1:end-1);

% stimulus properties, block randomized
task{2}{1}.randVars.block.contrast = logspace(-0.7,0,5);
% task{2}{1}.randVars.block.spatFreq = [0.15 0.3 0.6 1.2 2.4 4.8];
task{2}{1}.randVars.block.spatFreq = logspace(-0.8,0.5,5);
% task{2}{1}.randVars.block.orientation = 1:length(orientations);
task{2}{1}.randVars.block.nullTrial = [0 0 0 1]; %determines proportion of null trials

task{2}{1}.random = 0;
task{2}{1}.numTrials = numTrials;
% task{2}{1}.nTrials = numTrials;
task{2}{1}.collectEyeData = true;

% initialize the task
for phaseNum = 1:length(task{2})
    [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% do our initialization which creates the gratings
stimulus = myInitStimulus(stimulus,myscreen,task);
stimulus.orientations = orientations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % run the eye calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);


%initial screen
mglClearScreen;
totalRwd = task{1}{1}.numTrials * stimulus.rewardValue; %incr;
mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
%     if stimulus.rewardType == 'H';
%         text = sprintf('This is a high-reward run');
text = sprintf('You will gain or lose up to');
mglTextDraw(text,[0 2]);
text = sprintf('$%0.2f',totalRwd);
mglTextSet('Helvetica',80,[1 1 1],0,0,0,0,0,0,0);
mglTextDraw(text,[0 0]);
text = sprintf('in this run based on performance');
mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
mglTextDraw(text,[0 -2]);
mglFlush;

mglWaitSecs(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mglClearScreen(); mglFlush; mglClearScreen();
mglFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,[0 1 1], [0 0]) % default fixation cross
mglFlush
mglFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,[0 1 1], [0 0]) % default fixation cross
if ~waitForBacktick
    mglWaitSecs(3);
end

phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
    % update the fixation task
    [task{1} myscreen] = updateTask(task{1},myscreen,1);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end




% print out command for next run

%calculate reward for this run
if probRwd
    randP = randperm(task{1}{1}.numTrials);
    rwd = task{1}{1}.correctness(randP(1)) * task{1}{1}.numTrials * stimulus.rewardValue;%incr;
else
    rwd = sum(task{1}{1}.correctness) * stimulus.rewardValue;
end
stimulus.currBal = stimulus.currBal + rwd;
disp(sprintf('\n% --------------------------------------------- %\n'));

%calculate reward for next run
newRunNum = stimulus.runNum+1;
if stimulus.rewardType == 'H'
    newRewardType = 'L';
    %     newRewardValue = initRwdL + incrRwdL * newRunNum/2;
elseif stimulus.rewardType == 'L'
    newRewardType = 'H';
    %     newRewardValue = initRwdH + incrRwdH * newRunNum/2;
end

% disp(['rwdStim(''rewardType=''''' newRewardType ...
%     ''''''',''runNum=' num2str(newRunNum) ...
%     ''',''useStaircase=' num2str(fixStimulus.useStaircase) ...
%     ''',''currBal=' num2str(stimulus.currBal) ...
%     ''',''numTrials=' num2str(task{1}{1}.numTrials) ...
%     ''',''displayName=''''' myscreen.displayName ...
%     ''''''', ''fixThresh=' num2str(fixStimulus.threshold) ''');']);
disp(['rwdStim(''rewardType=''''' newRewardType ...
    ''''''',''runNum=' num2str(newRunNum) ...
    ''',''useStaircase=0' ...
    ''',''currBal=' num2str(stimulus.currBal) ...
    ''',''numTrials=' num2str(task{1}{1}.numTrials) ...
    ''',''displayName=''''' myscreen.displayName ...
    ''''''', ''fixThresh1=' num2str(fixStimulus.threshStair1) ...
    ''', ''fixThresh2=' num2str(fixStimulus.threshStair2) ''');']);


disp(sprintf('\n% --------------------------------------------- %\n'));




mglClearScreen;
% mglTextSet('Helvetica',50,[0 0.5 1 1],0,0,0,0,0,0,0);
mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
% text = sprintf('You got %0.2f%% correct', stimulus.percentCorrect*100);
% mglTextDraw(text,[0 3]);
text = sprintf('You gained');
mglTextDraw(text,[0 3]);
text = sprintf('$%0.2f', rwd);
mglTextSet('Helvetica',70,[1 1 1],0,0,0,0,0,0,0);
mglTextDraw(text,[0 1.5]);
text = sprintf('in this run');
mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
mglTextDraw(text,[0 0]);
text = sprintf('Current Balance: $%0.2f',stimulus.currBal);
mglTextDraw(text,[0 -4]);
mglFlush;
mglWaitSecs(3);


% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task, myscreen)
global stimulus;



% if stimulus.Counter >= 40
%             task.numTrials = task.trialnum-1;
%             return
% end



% randomize the current phase of the stimulus
newPhase = ceil(rand(1)*stimulus.numPhases);
while stimulus.phaseNum == newPhase
    newPhase = ceil(rand(1)*stimulus.numPhases);
end
stimulus.phaseNum = newPhase;
% randomize the current orientation of the stimulus
numOrientations = length(stimulus.orientations);
newOri = ceil(rand(1)*numOrientations);
while stimulus.oriNum == newOri
    newOri = ceil(rand(1)*numOrientations);
end
stimulus.oriNum = newOri;

% make the grating
% grating = mglMakeGrating(stimulus.width/stimulus.sFac, stimulus.height/stimulus.sFac, task.thistrial.spatFreq*stimulus.sFac, ...
%     task.thistrial.orientation, stimulus.phases(newPhase), stimulus.pixRes, stimulus.pixRes);
grating = mglMakeGrating(stimulus.width/stimulus.sFac, stimulus.height/stimulus.sFac, task.thistrial.spatFreq*stimulus.sFac, ...
    stimulus.orientations(newOri), stimulus.phases(newPhase), stimulus.pixRes, stimulus.pixRes);
grating = grating*task.thistrial.contrast;

% scale to range of display
grating = 255*(grating+1)/2;

% make it rgba
grating = uint8(permute(repmat(grating, [1 1 4]), [3 1 2]));
grating(4,:,:) = 256;

% update the texture
mglBindTexture(stimulus.tex, grating);

if task.numTrials == task.trialnum && task.thistrial.thisseg==length(task.seglen)
    disp('last trial');
end

% screen = mglFrameGrab;
% load('screens.mat','screens');
% screens{end+1} = screen(:,:,1)';
% save('screens.mat','screens');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = screenUpdateCallback(task, myscreen)

global stimulus;

% clear the screen
mglClearScreen;

if any(task.thistrial.thisseg == stimulus.stimulusSegments) && (~task.thistrial.nullTrial)
    % draw the texture
    mglBltTexture(stimulus.tex, [0 0 stimulus.height stimulus.height], 0, 0, 0);
    mglBltTexture(stimulus.innerMaskTex, [0 0 stimulus.height stimulus.height], 0, 0, 0);
    mglBltTexture(stimulus.outerMaskTex, [0 0 stimulus.height stimulus.height], 0, 0, 0);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)
% keep an array that lists which of the segments we
% are presenting the stimulus in.
% stimulus.stimulusSegments = 1:stimulus.stimLen/stimulus.frameLen;
stimulus.stimulusSegments = 1:(length(task{2}{1}.seglen) - 1);

stimulus.pixRes = min(myscreen.screenHeight/myscreen.imageHeight, myscreen.screenWidth/myscreen.imageWidth);

% scale factor
stimulus.sFac = 1;

% spatial frequency
stimulus.sf = 1.4;

% which phases we will have
stimulus.numPhases = 16;
stimulus.phases = 0:(360-0)/stimulus.numPhases:360;

% size of stimulus
stimulus.height = 0.5*floor(myscreen.imageHeight/0.5)+2;
stimulus.width = stimulus.height;

% size of annulus
% stimulus.outer = stimulus.height;
stimulus.outTransition = 0;

% stimulus.inner = 0.75;
stimulus.inTransition = 0;

% chose a sin or square
stimulus.square = 0;

% initial phase number
stimulus.phaseNum = 1;
% initial orientation number
stimulus.oriNum = 1;

% make a grating just to get the size
tmpGrating = mglMakeGrating(stimulus.width, stimulus.height, stimulus.sf, 0, stimulus.phases(1), stimulus.pixRes, stimulus.pixRes);
sz = size(tmpGrating,2);


% create mask for fixation and edge
out = stimulus.outer/stimulus.width;
in = stimulus.inner/stimulus.width;
twOut = stimulus.outTransition/stimulus.width;
twIn = stimulus.inTransition/stimulus.width;
outerMask = mkDisc(sz,(out*sz)/2,[(sz+1)/2 (sz+1)/2],twOut*sz,[1 0]);
innerMask = mkDisc(sz,(in*sz)/2,[(sz+1)/2 (sz+1)/2],twIn*sz,[0 1]);

% rescale mask to max out at 1
outerMask = outerMask/max(outerMask(:));
innerMask = innerMask/max(innerMask(:));

outerMask(:,:,4) = (-1*(outerMask*255))+255;
innerMask(:,:,4) = (-1*(innerMask*255))+255;

outerMask(:,:,1:3) = 128;
innerMask(:,:,1:3) = 128;

innerMask = uint8(permute(innerMask, [3 1 2]));
outerMask = uint8(permute(outerMask, [3 1 2]));

stimulus.innerMaskTex = mglCreateTexture(innerMask);
stimulus.outerMaskTex = mglCreateTexture(outerMask);





% make a grating again, but now scale it
tmpGrating = mglMakeGrating(stimulus.width/stimulus.sFac, stimulus.height/stimulus.sFac, stimulus.sf, 0, stimulus.phases(1), stimulus.pixRes, stimulus.pixRes);
r = uint8(permute(repmat(tmpGrating, [1 1 4]), [3 1 2]));
stimulus.tex = mglCreateTexture(r,[],1);



% fixStairInitTask.m
%
%        $Id$
%      usage: [fixTask myscreen] = fixStairInitTask(myscreen)
%         by: justin gardner
%       date: 09/07/06
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: Implements a fixation task. In this task, the fixation cross
%             starts out cyan and then darkens twice. The fixation cross then
%             turns yellow to indicate the response interval, and the subject
%             is required to press 1 or 2 to indicate in which interval the cross
%             appeared darker. The cross will then turn green or red to indicate
%             correct or incorrect responses. The dimness of the target is
%             controlled by a 2 down 1 up staircase to keep task difficulty
%             the same.
%
%             See testExperiment.m for how this is used in a task. If you want
%             to change parameters, before you call fixStairInitTask, set
%             appropriate fields of the global variable fixStimulus. e.g.:
%
%             global fixStimulus
%             fixStimulus.interTime = 1;
%
%             See the code, for a list of all parameters that can be changed.
%
function [task, myscreen] = myFixStairInitTask(myscreen)

% check arguments
if ~any(nargin == [1])
    help fixDispStairInitTask
    return
end

% create the stimulus for the experiment, use defaults if they are
% not already set
global fixStimulus;
myscreen = initStimulus('fixStimulus',myscreen);
% if ~isfield(fixStimulus,'threshold'); fixStimulus.threshold = 0.5; end
if ~isfield(fixStimulus,'pedestal'); fixStimulus.pedestal = 0.4; end
if ~isfield(fixStimulus,'stairUp'); fixStimulus.stairUp = 1; end
if ~isfield(fixStimulus,'stairDown'); fixStimulus.stairDown = 2; end
if ~isfield(fixStimulus,'stairStepSize'); fixStimulus.stairStepSize = 0.05; end
if ~isfield(fixStimulus,'stairUseLevitt'); fixStimulus.stairUseLevitt = 0; end
if ~isfield(fixStimulus,'stimColor'); fixStimulus.stimColor = [0 1 1]; end
if ~isfield(fixStimulus,'responseColor'); fixStimulus.responseColor = [1 1 0]; end
if ~isfield(fixStimulus,'interColor'); fixStimulus.interColor = [0 1 1]; end
if ~isfield(fixStimulus,'correctColor'); fixStimulus.correctColor = [0 0.8 0]; end
if ~isfield(fixStimulus,'incorrectColor'); fixStimulus.incorrectColor = [0.8 0 0]; end
if ~isfield(fixStimulus,'responseTime'); fixStimulus.responseTime = 1; end
if ~isfield(fixStimulus,'stimTime'); fixStimulus.stimTime = 0.4; end
if ~isfield(fixStimulus,'interTime'); fixStimulus.interTime = 0.8; end
if ~isfield(fixStimulus,'diskSize'); fixStimulus.diskSize = 1; end
if ~isfield(fixStimulus,'pos'); fixStimulus.pos = [0 0]; end
if ~isfield(fixStimulus,'fixWidth'); fixStimulus.fixWidth = 1; end
if ~isfield(fixStimulus,'fixLineWidth'); fixStimulus.fixLineWidth = 3; end
if ~isfield(fixStimulus,'trainingMode'); fixStimulus.trainingMode = 0;end
if ~isfield(fixStimulus,'verbose'); fixStimulus.verbose = 1;end

% for trainingMode set text
if fixStimulus.trainingMode
    mglTextSet('Helvetica',64,[1 1 1],0,0,0,0,0,0,0);
end

% create a fixation task
% task{1}.seglen = [fixStimulus.interTime fixStimulus.stimTime fixStimulus.interTime fixStimulus.stimTime fixStimulus.interTime fixStimulus.responseTime...
%     fixStimulus.trialTime-(3*fixStimulus.interTime + 2*fixStimulus.stimTime + fixStimulus.responseTime) - 0.5];
task{1}.seglen = [fixStimulus.stimTime fixStimulus.interTime fixStimulus.stimTime fixStimulus.interTime fixStimulus.responseTime...
    fixStimulus.trialTime-(2*fixStimulus.interTime + 2*fixStimulus.stimTime + fixStimulus.responseTime) - 0.5];
% task{1}.getResponse = [0 0 0 0 0 1];
task{1}.getResponse = [0 0 0 0 1 0];
task{1}.synchToVol = zeros(size(task{1}.seglen));
% task{1}.synchToVol(end) = 1;
if fixStimulus.waitForBacktick
    task{1}.synchToVol(end) = 1;
end


[task{1}, myscreen] = addTraces(task{1}, myscreen, 'segment', 'phase', 'response');

% % init the staircase
% fixStimulus.staircase = upDownStaircase(fixStimulus.stairUp,fixStimulus.stairDown,fixStimulus.threshold,fixStimulus.stairStepSize,fixStimulus.stairUseLevitt);
% fixStimulus.staircase.minThreshold = 0;
% fixStimulus.staircase.maxThreshold = 1;

% init a 2 down 1 up staircase
% if initStair == 1 | ieNotDefined('stimulus.stair{1}.threshold') | ieNotDefined('stimulus.stair{2}.threshold')
%     disp(sprintf('\nATTENTION: Initalizing new staircase: threshold1 = %0.2f, threshold2 = %0.2f\n', threshStair1, threshStair2));
fixStimulus.stair{1} = upDownStaircase(1,2,fixStimulus.threshStair1,fixStimulus.stairStepSize,1);
fixStimulus.stair{1}.minThreshold = 0;
fixStimulus.stair{2} = upDownStaircase(1,2,fixStimulus.threshStair2,fixStimulus.stairStepSize,1);
fixStimulus.stair{2}.minThreshold = 0;
% else
%     disp(sprintf('\nATTENTION: Continuing old staircase: threshold1 = %0.2f, threshold2 = %0.2f\n', stimulus.stair{1}.threshold, stimulus.stair{2}.threshold));
% end



% init the task
[task{1}, myscreen] = initTask(task{1},myscreen,@fixStartSegmentCallback,@fixDrawStimulusCallback,@fixTrialResponseCallback,@fixTrialStartCallback);

[task{1}, myscreen] = addTraces(task{1}, myscreen, 'fixStair');

% keep the correct and incorrect counts
task{1}.correct = 0;
task{1}.incorrect = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = fixTrialStartCallback(task, myscreen)

global fixStimulus;

% choose stimulus interval
task.thistrial.sigInterval = 1+(rand > 0.5);

% which staircase will be used this trial
fixStimulus.whichStair = round(rand)+1;

if fixStimulus.verbose
    disp(['trial = ' num2str(task.trialnum) ' stair = ' num2str(fixStimulus.whichStair) ' threshold = ' num2str(fixStimulus.stair{fixStimulus.whichStair}.threshold) ' sig = ' num2str(task.thistrial.sigInterval)]);
    %   disp(sprintf('sigint = %i threshold = %0.2f',task.thistrial.sigInterval,fixStimulus.threshold));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = fixStartSegmentCallback(task, myscreen)

global fixStimulus;

if isfield(fixStimulus,'displayText')
    % delete the texture
    if ~isempty(fixStimulus.displayText)
        mglDeleteTexture(fixStimulus.displayText);
    end
end
fixStimulus.displayText = [];

% choose what color the fixation cross will be
whichInterval = find(task.thistrial.thisseg == [1 3]);%[2 4]);

% if this is the signal interval
if ~isempty(whichInterval)
    if task.thistrial.sigInterval == whichInterval
        %     fixStimulus.thisStrength = (1-(fixStimulus.pedestal+fixStimulus.threshold));
        fixStimulus.thisStrength = (1-(fixStimulus.pedestal+fixStimulus.stair{fixStimulus.whichStair}.threshold));
    else
        fixStimulus.thisStrength = (1-fixStimulus.pedestal);
    end
    fixStimulus.thisColor = fixStimulus.stimColor*fixStimulus.thisStrength;
    % write out what the strength is
    myscreen = writeTrace(fixStimulus.thisStrength,task.fixStairTrace,myscreen);
    % if this is the response interval
elseif task.thistrial.thisseg == 5%6
    fixStimulus.thisColor = fixStimulus.responseColor;
    % if this is the inter stimulus interval
else
    fixStimulus.thisColor = fixStimulus.interColor;
    myscreen = writeTrace(0,task.fixStairTrace,myscreen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called every frame udpate to draw the fixation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = fixDrawStimulusCallback(task, myscreen)

global fixStimulus;

if ~isempty(fixStimulus.displayText)
    mglBltTexture(fixStimulus.displayText,fixStimulus.displayTextLoc);
end
mglGluDisk(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize*[1 1],myscreen.background,60);

mglFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,fixStimulus.thisColor,fixStimulus.pos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called when subject responds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = fixTrialResponseCallback(task, myscreen)

global fixStimulus;

% get correct or incorrect
response = find(task.thistrial.buttonState) == task.thistrial.sigInterval;
response = response(1);
%save response and correct response
task.response(task.trialnum) = find(task.thistrial.buttonState);
task.correctResponse(task.trialnum) = task.thistrial.sigInterval;
task.correctness(task.trialnum) = 2*response(1)-1;%1 or -1

if response
    % set to correct fixation color
    fixStimulus.thisColor = fixStimulus.correctColor;
    % set trace to 2 to indicate correct response
    myscreen = writeTrace(2,task.fixStairTrace,myscreen);
    % and update correct count
    task.correct = task.correct+1;
    if fixStimulus.verbose
        disp('YES');
    end
else
    % set to incorrect fixation color
    fixStimulus.thisColor = fixStimulus.incorrectColor;
    % set trace to -2 to indicate incorrect response
    myscreen = writeTrace(-2,task.fixStairTrace,myscreen);
    % and update incorrect count
    task.incorrect = task.incorrect+1;
    if fixStimulus.verbose
        disp('NO');
    end
end

% update staircase
if fixStimulus.useStaircase
    fixStimulus.stair{fixStimulus.whichStair} = upDownStaircase(fixStimulus.stair{fixStimulus.whichStair}, 1);
    %     fixStimulus.staircase = upDownStaircase(fixStimulus.staircase,response);
    %     fixStimulus.threshold = fixStimulus.staircase.threshold;
end