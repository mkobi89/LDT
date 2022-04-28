function [trl, event] = ft_trialfun_example1(cfg)

% FT_TRIALFUN_EXAMPLE1 is an example trial function. It searches for events
% of type "trigger" and specifically for a trigger with value 7, followed
% by a trigger with value 64.
% 
% You would use this function as follows
%   cfg           = [];   
%   cfg.dataset   = string, containing filename or directory
%   cfg.trialfun  = 'ft_trialfun_example1';
%   cfg           = definetrial(cfg);
%   data          = preprocessing(cfg);
%
% You can use this example trial function as template for your own
% conditial trial definitions.
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING
load(cfg.dataset)

% read the header information and the events from the data
hdr   = data.hdr;
event = data.event;

% cfg.dataformat = 'matlab';
% cfg.headerformat = 'matlab';

% search for "trigger" events
value  = [event(find(strcmp('trigger', {event.type}))).value]';
sample = [event(find(strcmp('trigger', {event.type}))).sample]';

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
posttrig =  round(cfg.trialdef.poststim * hdr.Fs);

% look for the combination of a trigger "7" followed by a trigger "64" 
% for each trigger except the last one
trl = [];
cond = 0;
for j = 1:(length(value)-1)
  
  trg1 = value(j);
  trg2 = value(j+1);
  
  if trg1==150 && trg2==160 || trg1==151 && trg2==161 || trg1==152 && trg2==162 ...
     || trg1==150 && trg2==165 || trg1==151 && trg2==166 || trg1==152 && trg2==167 ...   
     || trg1==250 && trg2==260 || trg1==251 && trg2==261 || trg1==252 && trg2==262 ...
     || trg1==250 && trg2==265 || trg1==251 && trg2==266 || trg1==252 && trg2==267 ...
     || trg1==301 && trg2==311 || trg1==302 && trg2==312 || trg1==303 && trg2==313 || trg1==304 && trg2==314 ...
     || trg1==301 && trg2==321 || trg1==302 && trg2==322 || trg1==303 && trg2==323 || trg1==304 && trg2==324 ...
     || trg1==305 && trg2==315 || trg1==306 && trg2==316 || trg1==307 && trg2==317 || trg1==308 && trg2==318 || trg1==309 && trg2==319 ...
     || trg1==305 && trg2==325 || trg1==306 && trg2==326 || trg1==307 && trg2==327 || trg1==308 && trg2==328 || trg1==309 && trg2==329

     
    trlbegin = sample(j) + pretrig;       
    trlend   = sample(j) + posttrig;       
    offset   = pretrig;
    
    rt(j)       = ((sample(j+1)-sample(j))/hdr.Fs)*1000;
     if trg2 == 160  || trg2==161 || trg2==162 || trg2 == 260  || trg2==261 || trg2==262
         correct = 1;
     elseif trg2 >= 311 && trg2 <= 319 
         correct = 1;
     else
         correct = 0;
     end
%     allcond(j) = trg1(j);
    
    newtrl   = [trlbegin trlend offset rt(j) trg1 trg2 correct];
    trl      = [trl; newtrl];
  end
    
end

