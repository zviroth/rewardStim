%
%        $Id: rwdStim.m,v 1.2 2016/01/22 15:30:22 eli Exp $
%      usage: rwdStim('useStaircase',1);
%         by: eli merriam
%       date: 01/27/07
%    purpose: oriented grating stimulus times a radial or angular modulator


% getStimvolFromVarname('orientation', s.myscreen, s.task{2})
% e = getTaskParameters(s.myscreen,s.task);
% mrPrintSurf

function myscreen = rwdStim(varargin)

% check arguments
if ~any(nargin == [0:10])
    help otopySzSF
    return
end

% % evaluate the input arguments
getArgs(varargin, [], 'verbose=0');

% set default parameters
% if ieNotDefined('direction'),direction = -1;end

if ieNotDefined('waitForBacktick'), waitForBacktick = 0; end
if ieNotDefined('useStaircase'), useStaircase = 0; end
interTime = 0.7;
stimTime = 0.3;
responseTime=1;
stimLen = 2*interTime + 2*stimTime + responseTime;
if ieNotDefined('stimLen'),stimLen = 1.5;end %should be equal to fixation trial length: 3*fixStimulus.interTime + 2*fixStimulus.stimTime + fixStimulus.responseTime
if ieNotDefined('trialLen'),trialLen = 9;end %in seconds
if ieNotDefined('frameLen'),frameLen = 0.25; end
if ieNotDefined('innerEdge'),innerEdge = 1.2; end
if ieNotDefined('outerEdge'),outerEdge = 12; end
if ieNotDefined('rewardVal'), rewardVal = 'H'; end
if ieNotDefined('probRwd'), probRwd = 0; end
if ieNotDefined('currBal'), currBal = 30; end
if ieNotDefined('fixThresh'), fixThresh = 0.2; end
TR=1.5;
numTRs = 15;%168;
if ieNotDefined('numTrials'), numTrials = ceil(TR*numTRs/trialLen); end

% check to see whether screen is still open
global stimulus;

% store parameters in stimulus variable
% update stimulus every 250 ms
stimulus.frameLen = frameLen;
% stimulus is on for 30 seconds
stimulus.stimLen = stimLen;
% % clockwise or counterclockwise
% stimulus.directon = direction;
% trial is 15 seconds
stimulus.trialLen = trialLen;
% inner and outer edges
stimulus.inner = innerEdge;
stimulus.outer = outerEdge;


%reward parameters
stimulus.currBal = currBal;
stimulus.rewardVal = rewardVal;
stimulus.probRwd = probRwd;
global incr;
if stimulus.rewardVal == 'H'
    incr = 1.1;
elseif stimulus.rewardVal == 'L'
    incr = 0.01;
end


% initalize the screen
myscreen.background = 'gray';
myscreen.autoCloseScreen = 0;
myscreen.allowpause = 1;
myscreen.saveData = 1;
%myscreen.displayName = '3tb';
% myscreen.displayName = '7t';
myscreen.displayName = 'rm315';
% myscreen.displayName = 'laptop';

myscreen = initScreen(myscreen);

global fixStimulus
fixStimulus.useStaircase = useStaircase;
fixStimulus.diskSize = 0;
fixStimulus.fixWidth = 0.75;
fixStimulus.fixLineWidth = 3;
fixStimulus.stairStepSize = 0.05;
fixStimulus.threshold = fixThresh;
fixStimulus.responseTime = responseTime;
fixStimulus.stimTime = stimTime;
fixStimulus.interTime = interTime;
fixStimulus.trialTime = trialLen;%same trial length for both tasks
fixStimulus.waitForBacktick = waitForBacktick
fixStimulus.fixWidth = 1;
fixStimulus.fixLineWidth = 3;

[task{1} myscreen] = myFixStairInitTask(myscreen);
task{1}{1}.numTrials = numTrials; 
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
% seglen = stimulus.frameLen * ones(1,stimulus.stimLen/stimulus.frameLen);
% seglen(end) = 0.1;
task{2}{1}.seglen = seglen;
task{2}{1}.synchToVol = zeros(size(seglen));
if waitForBacktick
    task{2}{1}.synchToVol(end) = 1;
end

orientations = linspace(0, 180, 17);
orientations = orientations(1:end-1);
stimulus.orientations = orientations;

% stimulus properties, block randomized
task{2}{1}.randVars.block.contrast = logspace(-0.7,0,5);
task{2}{1}.randVars.block.spatFreq = [0.15 0.3 0.6 1.2 2.4 4.8];
task{2}{1}.randVars.block.orientation = 1:length(orientations);
task{2}{1}.randVars.block.nullTrial = [0 0 0 1]; %determines proportion of null trials

task{2}{1}.random = 0;
task{2}{1}.numTrials = numTrials;
task{2}{1}.collectEyeData = true;

% initialize the task
for phaseNum = 1:length(task{2})
    [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% do our initialization which creates the gratings
stimulus = myInitStimulus(stimulus,myscreen,task);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % run the eye calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);


%initial screen
mglClearScreen;
totalRwd = task{1}{1}.numTrials * incr;
mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
%     if stimulus.rewardVal == 'H';
%         text = sprintf('This is a high-reward run');
text = sprintf('You will gain or lose up to');
        mglTextDraw(text,[0 2]);
        text = sprintf('$%0.2f',totalRwd);
        mglTextSet('Helvetica',80,[1 1 1],0,0,0,0,0,0,0);
        mglTextDraw(text,[0 0]);
        text = sprintf('in this run based on performance');
        mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
        mglTextDraw(text,[0 -2]);
%         text = sprintf('You will gain or lose up to /n $%0.2f /n in this run based on performance',totalRwd);
%         mglTextDraw(text,[0 -3]);
%     elseif stimulus.rewardVal == 'L';
%         text = sprintf('This is a low-reward run');
%         mglTextDraw(text,[0 3]);
%         text = sprintf('You will gain or lose up to $%0.2f in this run based on performance',totalRwd);
%         mglTextDraw(text,[0 -3]);
%     end
    mglFlush;


mglWaitSecs(3);

% mglFixationCross(0.7, 2, [0 1 1], [0 0]) % default fixation cross

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

% stimulus.percentCorrect = nansum(stimulus.trialCorrectness)/length(stimulus.trialCorrectness); disp('Percent Correct:'); disp(stimulus.percentCorrect*100);
% stimulus.avgRT = nanmean(stimulus.trialRT); disp('Average RT:'); disp(stimulus.avgRT)
% stimulus.currBal = stimulus.currBal + reward - punish; disp('Current balance:'); disp(stimulus.currBal);
% 
if probRwd
    randP = randperm(task{1}{1}.numTrials);
   rwd = task{1}{1}.correctness(randP(1)) * task{1}{1}.numTrials * incr;
else
%     rwd = (task{1}{1}.correct - task{1}{1}.incorrect)*incr;
    rwd = sum(task{1}{1}.correctness) * incr;
end
stimulus.currBal = stimulus.currBal + rwd;
disp(sprintf('\n% --------------------------------------------- %\n'));
if stimulus.rewardVal == 'H'
  disp(sprintf('arousalRwd(''rewardVal=''''L'''''', ''currBal=%0.2f'', ''fixThresh=%0.2f'');', stimulus.currBal, fixStimulus.threshold));
elseif stimulus.rewardVal == 'L'
  disp(sprintf('arousalRwd(''rewardVal=''''H'''''', ''currBal=%0.2f'', ''fixThresh=%0.2f'');', stimulus.currBal, fixStimulus.threshold));
end
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
text = sprintf('in this run', rwd); 
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
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

% randomize the current phase of the stimulus
newPhase = ceil(rand(1)*stimulus.numPhases);
while stimulus.phaseNum == newPhase;
    newPhase = ceil(rand(1)*stimulus.numPhases);
end
stimulus.phaseNum = newPhase;

% make the grating
grating = mglMakeGrating(stimulus.width/stimulus.sFac, stimulus.height/stimulus.sFac, task.thistrial.spatFreq*stimulus.sFac, ...
    task.thistrial.orientation, stimulus.phases(newPhase), stimulus.pixRes, stimulus.pixRes);
grating = grating*task.thistrial.contrast;
% make it a square wave
if stimulus.square
    grating = sign(grating);
end

% scale to range of display
grating = 255*(grating+1)/2;

% make it rgba
grating = uint8(permute(repmat(grating, [1 1 4]), [3 1 2]));
grating(4,:,:) = 256;

% update the texture
mglBindTexture(stimulus.tex, grating);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

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

% make a grating just to get the size
tmpGrating = mglMakeGrating(stimulus.width, stimulus.height, stimulus.sf, 0, stimulus.phases(1), stimulus.pixRes, stimulus.pixRes);
sz = size(tmpGrating,2);

% outer = linspace(8.5, 12, 8);
% inner = linspace(1, 4.5, 8);

% outer = linspace(7, 10.5, 8);
% inner = linspace(0.5, 4, 8);


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
function [task myscreen] = myFixStairInitTask(myscreen)

% check arguments
if ~any(nargin == [1])
  help fixDispStairInitTask
  return
end

% create the stimulus for the experiment, use defaults if they are
% not already set
global fixStimulus;
myscreen = initStimulus('fixStimulus',myscreen);
if ~isfield(fixStimulus,'threshold') fixStimulus.threshold = 0.5; end
if ~isfield(fixStimulus,'pedestal') fixStimulus.pedestal = 0.4; end
if ~isfield(fixStimulus,'stairUp') fixStimulus.stairUp = 1; end
if ~isfield(fixStimulus,'stairDown') fixStimulus.stairDown = 2; end
if ~isfield(fixStimulus,'stairStepSize') fixStimulus.stairStepSize = 0.05; end
if ~isfield(fixStimulus,'stairUseLevitt') fixStimulus.stairUseLevitt = 0; end
if ~isfield(fixStimulus,'stimColor') fixStimulus.stimColor = [0 1 1]; end
if ~isfield(fixStimulus,'responseColor') fixStimulus.responseColor = [1 1 0]; end
if ~isfield(fixStimulus,'interColor') fixStimulus.interColor = [0 1 1]; end
if ~isfield(fixStimulus,'correctColor') fixStimulus.correctColor = [0 0.8 0]; end
if ~isfield(fixStimulus,'incorrectColor') fixStimulus.incorrectColor = [0.8 0 0]; end
if ~isfield(fixStimulus,'responseTime') fixStimulus.responseTime = 1; end
if ~isfield(fixStimulus,'stimTime') fixStimulus.stimTime = 0.4; end
if ~isfield(fixStimulus,'interTime') fixStimulus.interTime = 0.8; end
if ~isfield(fixStimulus,'diskSize') fixStimulus.diskSize = 1; end
if ~isfield(fixStimulus,'pos') fixStimulus.pos = [0 0]; end
if ~isfield(fixStimulus,'fixWidth') fixStimulus.fixWidth = 1; end
if ~isfield(fixStimulus,'fixLineWidth') fixStimulus.fixLineWidth = 3; end
if ~isfield(fixStimulus,'trainingMode') fixStimulus.trainingMode = 0;end
if ~isfield(fixStimulus,'verbose') fixStimulus.verbose = 1;end

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
    
 
[task{1} myscreen] = addTraces(task{1}, myscreen, 'segment', 'phase', 'response');

% init the staircase
fixStimulus.staircase = upDownStaircase(fixStimulus.stairUp,fixStimulus.stairDown,fixStimulus.threshold,fixStimulus.stairStepSize,fixStimulus.stairUseLevitt);
fixStimulus.staircase.minThreshold = 0;
fixStimulus.staircase.maxThreshold = 1;

% init the task
[task{1} myscreen] = initTask(task{1},myscreen,@fixStartSegmentCallback,@fixDrawStimulusCallback,@fixTrialResponseCallback,@fixTrialStartCallback);

[task{1} myscreen] = addTraces(task{1}, myscreen, 'fixStair');

% keep the correct and incorrect counts
task{1}.correct = 0;
task{1}.incorrect = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixTrialStartCallback(task, myscreen)

global fixStimulus;
% choose stimulus interval
task.thistrial.sigInterval = 1+(rand > 0.5);
if fixStimulus.verbose
  disp(sprintf('sigint = %i threshold = %0.2f',task.thistrial.sigInterval,fixStimulus.threshold));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixStartSegmentCallback(task, myscreen)

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
    fixStimulus.thisStrength = (1-(fixStimulus.pedestal+fixStimulus.threshold));
  else
    fixStimulus.thisStrength = (1-fixStimulus.pedestal);
  end
  fixStimulus.thisColor = fixStimulus.stimColor*fixStimulus.thisStrength;
  % write out what the strength is
  myscreen = writeTrace(fixStimulus.thisStrength,task.fixStairTrace,myscreen);
  % training mode text
  if fixStimulus.trainingMode
    if whichInterval == 1
      fixStimulus.displayText = mglText('Interval 1');
      fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
    else
      fixStimulus.displayText = mglText('Interval 2');
      fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
    end
  end
% if this is the response interval
elseif task.thistrial.thisseg == 5%6
  fixStimulus.thisColor = fixStimulus.responseColor;
  if fixStimulus.trainingMode
    fixStimulus.displayText = mglText('Which interval was darker (1 or 2)?');
    fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
  end
% if this is the inter stimulus interval
else
  % training mode, clear sceen here
  if fixStimulus.trainingMode,mglClearScreen;end
  fixStimulus.thisColor = fixStimulus.interColor;
  myscreen = writeTrace(0,task.fixStairTrace,myscreen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called every frame udpate to draw the fixation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixDrawStimulusCallback(task, myscreen)

global fixStimulus;

if fixStimulus.trainingMode,mglClearScreen;end

if ~isempty(fixStimulus.displayText)
  mglBltTexture(fixStimulus.displayText,fixStimulus.displayTextLoc);
end
mglGluDisk(fixStimulus.pos(1),fixStimulus.pos(2),fixStimulus.diskSize*[1 1],myscreen.background,60);

mglFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,fixStimulus.thisColor,fixStimulus.pos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called when subject responds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = fixTrialResponseCallback(task, myscreen)

global fixStimulus;

% get correct or incorrect
response = find(task.thistrial.buttonState) == task.thistrial.sigInterval;
response = response(1);
%save response and correct response
task.response(task.trialnum) = find(task.thistrial.buttonState);
task.correctResponse(task.trialnum) = task.thistrial.sigInterval;
task.correctness(task.trialnum) = 2*response(1)-1;%1 or -1

if response
  % for training mode
  if fixStimulus.trainingMode
    mglDeleteTexture(fixStimulus.displayText);
    fixStimulus.displayText = mglText('Correct!');
    fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
  end
  % set to correct fixation color
  fixStimulus.thisColor = fixStimulus.correctColor;
  % set trace to 2 to indicate correct response
  myscreen = writeTrace(2,task.fixStairTrace,myscreen);
  % and update correct count
  task.correct = task.correct+1;

else
  if fixStimulus.trainingMode
    mglDeleteTexture(fixStimulus.displayText);
    fixStimulus.displayText = mglText('Incorrect.');
    fixStimulus.displayTextLoc = [fixStimulus.pos(1) fixStimulus.pos(2)+fixStimulus.diskSize*2];
  end
  % set to incorrect fixation color
  fixStimulus.thisColor = fixStimulus.incorrectColor;
  % set trace to -2 to indicate incorrect response
  myscreen = writeTrace(-2,task.fixStairTrace,myscreen);
  % and update incorrect count
  task.incorrect = task.incorrect+1;
end

% update staircase
if fixStimulus.useStaircase
    fixStimulus.staircase = upDownStaircase(fixStimulus.staircase,response);
    fixStimulus.threshold = fixStimulus.staircase.threshold;
end