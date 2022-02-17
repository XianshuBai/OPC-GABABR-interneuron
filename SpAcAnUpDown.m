function SpAcAnUpDown
    upordown = questdlg('Signal goes up or down?','Signal direction', ...
        'Up','Down','Up');
    upordown = strcmpi(upordown,'Up');
    %% Main Interface.
    GuiSize=FigSize;
    S.fh = figure('units','pixels','position',GuiSize,'menubar','none',...
        'name','SpAcAn','numbertitle','off','resize','off',...
        'Color',[0.95 0.95 0.95],'HandleVisibility','off');
    
    %% Load file.
    S.pbLoadFile = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[12 GuiSize(4)-35 120 30],'string','Load','fontWeight','bold');
    S.pbSetPara = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[12 GuiSize(4)-70 120 25],'string','Set Parameters','fontWeight','bold');
    S.pbSaveAll = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[12 GuiSize(4)-95 120 25],'string','Save All','fontWeight','bold');
    %% next trace.
    S.pbNext = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[12 GuiSize(4)-140 120 25],'string','Next','fontWeight','bold');
    S.pbPre = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[12 GuiSize(4)-170 120 25],'string','Previous','fontWeight','bold');
    S.pbBegin = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[12 GuiSize(4)-200 30 25],'string','<=','fontWeight','bold');
    S.pbGoTo = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[45 GuiSize(4)-200 50 25],'string','Go To','fontWeight','bold');
    S.pbEnd = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[102 GuiSize(4)-200 30 25],'string','=>','fontWeight','bold');
    %% analysis
    S.pbAnalyze = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[12 GuiSize(4)-245 120 30],'string','Analyze','fontWeight','bold');
    %% Combine all avg
    S.pbAvg = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[12 GuiSize(4)-285 120 25],'string','Combine Avg','fontWeight','bold');
    S.pbAvgSmooth = uicontrol('Parent',S.fh,'style','push','unit','pix',...
        'position',[12 GuiSize(4)-310 120 25],'string','Combine AvgSmooth','fontWeight','bold');
    %%
    set(S.pbLoadFile,'callback',{@pbLoadFile,S});
    set(S.pbPre,'callback',{@pbPlotNext,S});
    set(S.pbNext,'callback',{@pbPlotNext,S});
    set(S.pbBegin,'callback',{@pbPlotNext,S});
    set(S.pbGoTo,'callback',{@pbPlotNext,S});
    set(S.pbEnd,'callback',{@pbPlotNext,S});
    set(S.pbAnalyze,'callback',{@pbAnalyze,S})
    set(S.pbAvg,'callback',{@combineAvg,S});
    set(S.pbAvgSmooth,'callback',{@combineAvg,S});
    set(S.pbSetPara,'callback',{@setPara,S});
    set(S.pbSaveAll,'callback',{@SaveFile,S});
    setappdata(S.fh,'upordown',upordown);
end


function Para=detectionParameters
    leading=0.005; % second
    ScaleN=[0.0002,0.001,0.02];  % 1 ms, 0.5ms, 2ms.
    scale=1e12;
    smoothWindow=0.005; % second
    threshold=4.5;
    lengthAvgAll=20;
    
    Para.leading=leading;
    Para.ScaleN=ScaleN;
    Para.scale=scale;
    Para.smoothWindow=smoothWindow;
    Para.threshold=threshold;
    Para.lengthAvgAll=lengthAvgAll;
end

function setPara(varargin)
    S = varargin{3}; % varargin = [3.0015]    []    [1x1 struct];
    S_fig = S.fh;
    if isappdata(S_fig,'Para')
        Para=getappdata(S_fig,'Para');
    else
        Para=detectionParameters;
    end
    upordown = inputdlg({'Leading time (s):','ScaleN (s)','scale',...
        'smooth window (s)','threshold','length for Average (xleading):'},...
        'Parameters',...
        1,...
        {num2str(Para.leading),num2str(Para.ScaleN),num2str(Para.scale),...
        num2str(Para.smoothWindow),num2str(Para.threshold),num2str(Para.lengthAvgAll)});
    Para.leading=str2double(upordown{1});
    Para.ScaleN=str2num(upordown{2}); %#ok<ST2NM>
    Para.scale=str2double(upordown{3});
    Para.smoothWindow=str2double(upordown{4});
    Para.threshold=str2double(upordown{5});
    Para.lengthAvgAll=str2double(upordown{6});
    setappdata(S_fig,'Para',Para);
end

function SaveFile(varargin)
    S = varargin{3}; % varargin = [3.0015]    []    [1x1 struct];
    S_fig = S.fh;
    file=getappdata(S_fig,'file'); %#ok<NASGU>
    data=getappdata(S_fig,'data'); %#ok<NASGU>
    Para=getappdata(S_fig,'Para'); %#ok<NASGU>
    [filename,path]=uiputfile('*.mat');
    save(fullfile(path,filename),'file','data','Para');
end

function pbLoadFile(varargin)
    S = varargin{3}; % varargin = [3.0015]    []    [1x1 struct];
    S_fig = S.fh;
    
    [data,file]=HEKARead();
    if ~isempty(data) && ~isempty(file)
        setappdata(S_fig,'file',file);
        setappdata(S_fig,'data',data);
        setappdata(S_fig,'No',1);
        Para=detectionParameters;
        setappdata(S_fig,'Para',Para);
    end
end

function pbAnalyze(varargin)
    S = varargin{3};
    S_fig = S.fh;
    currNo=getappdata(S_fig,'No');
    data=getappdata(S_fig,'data');
    Para=getappdata(S_fig,'Para');
    No=labeldata(data);
    
    Istruct=data{currNo,5};
    variableName=['Tr_',num2str(No(currNo,1)),'_',num2str(No(currNo,2)),...
        '_',num2str(No(currNo,3)),'_',num2str(No(currNo,4)),'_',...
        num2str(No(currNo,5))];
    I=Istruct.TrRawData;
    TrXInterval=Istruct.TrXInterval;
    upordown=getappdata(S_fig,'upordown');
    
    [results,Ismooth,avg,avgsmooth]=SpAcAnAnalysis(TrXInterval,I,variableName,currNo,upordown,Para);
    Istruct.variableName=variableName;
    Istruct.results=results;
    Istruct.Ismooth=Ismooth;
    Istruct.Iavg=avg;
    Istruct.Iavgsmooth=avgsmooth;
    data{currNo,5}=Istruct;
    setappdata(S_fig,'data',data);
end

function combineAvg(varargin)
    S = varargin{3};
    S_fig=S.fh;
    currButton=varargin{1};
    
    data=getappdata(S_fig,'data');
    currNo=0;
    all=[];
    if currButton==S.pbAvg
        while currNo<size(data,1)
            currNo=currNo+1;
            if ~isempty(data{currNo,5})
                Istruct=data{currNo,5};
                if isfield(Istruct,'Iavg')
                    all=cat(2,all,Istruct.Iavg);
                    disp([Istruct.variableName,'_Avg']);
                end
            end
        end
        if numel(all)>0
            assignin('base','CombinedAvg',all);
            evalin('base','open(''CombinedAvg'')');
        end
    elseif currButton==S.pbAvgSmooth
        while currNo<size(data,1)
            currNo=currNo+1;
            if ~isempty(data{currNo,5})
                Istruct=data{currNo,5};
                if isfield(Istruct,'Iavgsmooth')
                    all=cat(2,all,Istruct.Iavgsmooth);
                    disp([Istruct.variableName,'_AvgSmooth']);
                end
            end
        end
        if numel(all)>0
            assignin('base','CombinedAvgSmooth',all);
            evalin('base','open(''CombinedAvgSmooth'')');
        end
    end
end


function pbPlotNext(varargin)
    S = varargin{3};
    currButton=varargin{1};
    S_fig = S.fh;
    currNo=getappdata(S_fig,'No');
    data=getappdata(S_fig,'data');
    No=labeldata(data);
    if currButton==S.pbNext
        while currNo<size(data,1)
            currNo=currNo+1;
            if ~isempty(data{currNo,5})
                break;
            end
        end
        if currNo<=size(data,1)
            setappdata(S_fig,'No',currNo);
        else
            warning('Aready at the End of data file!!!');
            return;
        end
    elseif currButton==S.pbPre
        while (currNo>1)
            currNo=currNo-1;
            if ~isempty(data{currNo,5})
                break;
            end
        end
        if currNo>1
            setappdata(S_fig,'No',currNo);
        else
            warning('Aready at the Beginning of data file!!!');
            return;
        end
    elseif currButton==S.pbBegin
        currNo=1;
        while (currNo<size(data,1)) && isempty(data{currNo,5})
            currNo=currNo+1;
        end
        setappdata(S_fig,'No',currNo);
    elseif currButton==S.pbEnd
        currNo=size(data,1);
        while currNo>1 && isempty(data{currNo,5})
            currNo=currNo-1;
        end
        setappdata(S_fig,'No',currNo);
    elseif currButton==S.pbGoTo
        x=inputdlg({'Root','Group No.:','Series No.:','Sweep No.:','Trace No.:'},'Goto...',1,{'1','1','1','1','1'});
        Group=str2double(x{2});
        Series=str2double(x{3});
        Sweep=str2double(x{4});
        Trace=str2double(x{5});
        bw=(No(:,2)==Group) & (No(:,3)==Series) & (No(:,4)==Sweep) & (No(:,5)==Trace);
        currNo=1:size(No,1);
        currNo=currNo(bw);
        if numel(currNo)>1
            currNo=currNo(1);
        elseif numel(currNo)==0
            return;
        end
        setappdata(S_fig,'No',currNo);
        % return;
    end

    variableName=['Tr_',num2str(No(currNo,1)),'_',num2str(No(currNo,2)),...
        '_',num2str(No(currNo,3)),'_',num2str(No(currNo,4)),'_',...
        num2str(No(currNo,5))];
    
    
    Istruct=data{currNo,5};
    I=Istruct.TrRawData;
    TrXInterval=Istruct.TrXInterval;
    t=(1:numel(I))*TrXInterval;
    plot(t,I);
    axis tight;
    xlabel('Time (s)')
    ylabel('Current Density (A)')
    title(['Row No.: ',num2str(currNo),'  ',variableName],'Interpreter','none')
end

function Pos=FigSize
    Screen=get(0,'ScreenSize');
    FigHeight=320;
    FigWidth=150;
    Pos=[Screen(3)-FigWidth-300 round((Screen(4)-FigHeight-50)) FigWidth FigHeight];
end
function No=labeldata(data)
    No=zeros(size(data));
    for j=1:size(data,2)
        m=0;
        for k=1:size(data,1)
            if ~isempty(data{k,j})
                m=m+1;
                if (k>1) && (j>1) && (~isempty(data{k-1,j-1}))
                    m=1;
                end
            end
            No(k,j)=m;
        end
    end
    
    for k=1:size(data,1)
        for j=1:4
            if ~isempty(data{k,j})
                for m=j:5
                    No(k,m)=NaN;
                end
            end
        end
    end
end




% =====================================================================%
% =====================================================================%
% ===========   Original Analysis is here==============================%
% =====================================================================%
% =====================================================================%

function [results,Ismooth,avg,avgsmooth]=SpAcAnAnalysis(tStep,I,variableName,currNo,upordown,Para)
    % spontaneous activity detection of mFPC.
    % Usage: SpAcAn(Trace_1_17_1_1)
    % Changing the threshold value in line 12 can change the detection limits.
    
    %% Parameters.
    % leading=0.005; % second
    % ScaleN=[0.0005,0.001,0.02];  % 1 ms, 0.5ms, 2ms.
    % scale=1e12;
    % smoothWindow=0.005; % second
    % threshold=4.5;
    
    % [leading,ScaleN,scale,smoothWindow,threshold,lengthAvgAll]=detectionParameters;
    leading=Para.leading;
    ScaleN=Para.ScaleN;
    scale=Para.scale;
    smoothWindow=Para.smoothWindow;
    threshold=Para.threshold;
    lengthAvgAll=Para.lengthAvgAll;
    
    %% Scale current and make it positive.
    if upordown
        Ismooth=I;
    else
        Ismooth=-I;
    end
    Im=median(Ismooth);
    Ismooth=Ismooth-Im;
    Ismooth=Ismooth*scale;
    %% parameters in points.
    % tStep=median(diff(t));
    t=(1:numel(I))*tStep;
    smoothWindow=round(smoothWindow/tStep);
    ScaleN=round(ScaleN/tStep);
    leadingPnt=round(leading/tStep);
    %% Smooth data
    % Ismooth=smoothdata(I,'lowess',11);
    Ismooth=smoothdata(Ismooth,'sgolay',smoothWindow);
    
    
    %% Detect spikes
    bw=PeakThresholding(Ismooth,'ScaleN',ScaleN,'Threshold',threshold);
    bw=labelBW(bw);
    
    %% Get location and amplitude
    tlen=numel(bw);
    x=(1:tlen)';
    bwMax=max(bw);
    results=zeros(bwMax,5);
    for k=1:bwMax
        currbw=bw==k;
        tStart=min(x(currbw))-leadingPnt;
        if tStart<1
            tStart=1;
        end
        
        xmax=max(x(currbw));
        currBase=Ismooth(tStart:xmax);
        currBasePos=find(currBase==min(currBase));
        if numel(currBasePos)>1
            currBasePos=round(mean(currBasePos));
        end
        currBase=median(currBase); % currBase(currBasePos);
        currBasePos=currBasePos+tStart-1;
        
        % Find next peak.
        while ((xmax+1)<=tlen) && (~(bw(xmax+1)))
            xmax=xmax+1;
        end
        tStop=xmax;
        
        currPeak=Ismooth(tStart:tStop);
        currPeakPos=find(currPeak==max(currPeak));
        if numel(currPeakPos)>1
            currPeakPos=round(mean(currPeakPos));
        end
        currPeak=currPeak(currPeakPos);
        currPeakPos=currPeakPos+tStart-1;
        
        results(k,1)=t(currBasePos);
        results(k,2)=t(currPeakPos);
        results(k,3)=currBase;
        results(k,4)=currPeak;
        results(k,5)=currPeak-currBase;
    end
    
    if upordown
        results(:,3)=results(:,3)/scale+Im;
        results(:,4)=results(:,4)/scale+Im;
        results(:,5)=results(:,4)-results(:,3);
    else
        results(:,3)=results(:,3)/scale+Im;
        results(:,3)=-results(:,3);
        results(:,4)=results(:,4)/scale+Im;
        results(:,4)=-results(:,4);
        results(:,5)=results(:,4)-results(:,3);
    end
    
    assignin('base',variableName,I);
    assignin('base',[variableName,'_smooth'],Ismooth);
    assignin('base',[variableName,'_results'],results);
    
    Ismooth=Ismooth/scale+Im;
    if ~upordown
        Ismooth=-Ismooth;
    end
    [avg,avgsmooth]=averageFit(Ismooth,I,results,leadingPnt,tStep,lengthAvgAll);
    
    assignin('base',[variableName,'_Avg'],avg);
    assignin('base',[variableName,'_AvgSmooth'],avgsmooth);
    
    subplot(1,6,1:4)
    plot(t,I,'color',[0.8,0.8,0.8]);
    hold all
    plot(t,Ismooth,'color',[0.8,0,0.0]);
    plot(results(:,1),results(:,3),'o','color',[0,0,1]);
    plot(results(:,2),results(:,4),'*','color',[0,0.5,1]);
    hold off
    axis tight
    title(['Row No.: ',num2str(currNo),'  ',variableName],'Interpreter','none')
    xlabel('Time (s)')
    ylabel('Current Density (A)')
    
    subplot(1,6,5)
    if upordown
        hist(results(:,5),linspace(0,max(results(:,5)),30));
    else
        hist(results(:,5),linspace(min(results(:,5)),0,30));
    end
    axis tight
    xlabel('Current (A)')
    ylabel('Number of Spikes')
    
    subplot(1,6,6)
    t=1:numel(avg);
    t=t*tStep*1000;
    plot(t,avg);hold all
    plot(t,avgsmooth); hold off;
    axis tight
    xlabel('Time (ms)')
    ylabel('Current Density (A)')
    subplot(1,6,1:4)
    
    % set(gcf,'position',[134,440,2334,830]);
end


function [avg,avgsmooth]=averageFit(Ismooth,I,results,leading,tStep,lengthAvgAll)
    % lengthAll=10;
    avg=zeros(leading*lengthAvgAll,1);
    avgsmooth=zeros(leading*lengthAvgAll,1);
    
    m=0;
    for k=1:size(results,1)
        x1=round(results(k,2)/tStep)-leading*3;
        x2=x1+leading*lengthAvgAll-1;
        if x1<1 || x2>numel(I)
            continue;
        end
        avg=avg+I(x1:x2);
        avgsmooth=avgsmooth+Ismooth(x1:x2);
        m=m+1;
    end
    avg=avg/m;
    avgsmooth=avgsmooth/m;
end


function bw=PeakThresholding(I,varargin)
    % Thresholding sparks.
    %   Syntax:
    %       bw=SparkThresholding('ScaleN',[2,10],'Threshold',3,'mask',cellmask);
    %   Input:
    %       ScaleN:     2-10 points in the third dimension.
    %       Threshold:  the threshold to sperate sparks from their background.
    %   Output:
    %       bw: the binayr image stack marking positive pixels.
    
    % Inputs.
    p=inputParser;
    p.addParameter('ScaleN',[2,1,6],@(x)numel(x)==3 && min(x)>=1 && x(1)<=x(3));
    p.addParameter('Threshold',3.5,@(x) isscalar(x)||size(x,1)==size(I,1)&&size(x,2)==size(I,2));
    parse(p, varargin{:});
    p=p.Results;
    clear('varargin');
    
    bw=zeros(size(I));
    
    for k=p.ScaleN(1):p.ScaleN(2):p.ScaleN(3)
        h=[ones(1,k) -ones(1,k)];   h=reshape(h,[numel(h),1]);
        cA=padarray(I,k*2-1,'replicate','both');
        
        cA=convn(cA,h,'same');
        cA=cA(k*2:k*2-1+numel(I));
        
        cA=cA-median(cA);
        sigma=median(abs(cA))/0.6745;
        bw=bw+double(cA>sigma.*p.Threshold);
    end
    bw=bw>=1;
end




function bwLabel=labelBW(bw)
    bwLabel=zeros(size(bw));
    if bw(1)
        m=1;
    else
        m=0;
    end
    bwLabel(1)=m;
    for k=2:numel(bw)
        if (bw(k-1) && bw(k))
            bwLabel(k)=m;
        elseif (~bw(k-1) && bw(k))
            m=m+1;
            bwLabel(k)=m;
        end
    end
end






















function varargout=HEKARead(string)
    % HEKARead imports HEKA PatchMaster and ChartMaster*.DAT files.
    %
    % Example:
    % tree=HEKARead(FILENAME);
    % tree=HEKARead();
    %
    % FILENAME is the path and name of the HEKA DAT file to import.
    %
    % The kcl file generated will be placed in TARGETPATH if supplied. If not,
    % the file will be created in the directory taken from FILENAME.
    %
    % HEKARead has been tested with Windows generated. DAT files on Windows,
    % Linux and Mac OS10.4.
    %
    % Both bundled and unbundled data files are supported. If your files are
    % unbundled,they must all be in the same folder.
    %
    %
    % Notes:
    % Timestamps from the data file are rounded to the nearest nanonsecond for
    % sigTOOL.
    % Waveform data are scaled to SI units of A or V in HEKA files. For
    % sigTOOL,they are scaled to pA,pV,nA,nV... etc as appropriate given
    % the data range.
    %
    % The HEKA DAT format is versatile and not all combinations of settings may
    % have been anticipated here. If you encounter problems importing files
    % please report the bug and send a sample DAT file using Help->Bug Report
    % in the sigTOOL GUI
    %
    % Details of the HEKA file format are available from
    %       ftp://server.hekahome.de/pub/FileFormat/Patchmasterv9/
    %
    %--------------------------------------------------------------------------
    % Author: Malcolm Lidierth 12/09
    % Copyright ? The Author & King's College London 2009-
    %--------------------------------------------------------------------------
    %
    % Revisions
    % 17.04.10  TrXUnit: see within
    % 28.11.11  TrXUnit: see within
    % 15.08.12  Updated to support interleaved channels and PatchMaster 2.60
    %               files dated 24-Jan-2011 onwards.
    
    
    %% Check DataFile.
    varargout(1)={[]};
    varargout(2)={[]};
    switch (nargin)
        case 0                                            % If there is no input.
            [filename,pathname] = uigetfile({'*.dat','HEKA PatchMaster File (*.dat)'},...
                'Please select a HEKA file (*.dat):');
        case 1                                            % If there is one string input.
            if exist(string,'file') ==2                   % If the string input is a file name.
                [pathname,filename,~] = fileparts(string);
            elseif exist(string,'dir')==7                 % If the string input is a folder name.
                [filename,pathname] = uigetfile({'*.dat','HEKA PatchMaster File (*.dat)'},...
                    'Please select a HEKA file (*.dat):',string);
            else
                [filename,pathname] = uigetfile({'*.dat','HEKA PatchMaster File (*.dat)'},...
                    'Please select a HEKA file (*.dat):');
            end
    end
    clear('string','ext')
    
    
    if filename==0
        return;
    else
        [~,filename,~]=fileparts(filename);
    end
    if isempty(pathname); pathname=pwd; end
    datafile=fullfile(pathname,[filename,'.dat']);
    
    
    %% Open file and get bundle header. Assume little-endian to begin with
    endian='ieee-le';
    fh=fopen(datafile,'r',endian);
    [bundle,littleendianflag,isBundled]=getBundleHeader(fh);
    
    
    % Big endian so repeat process
    if ~isempty(littleendianflag) && littleendianflag==false
        fclose(fh);
        endian='ieee-be';
        fh=fopen(datafile,'r',endian);
        bundle=getBundleHeader(fh);
    end
    clear('littleendianflag','thisfile')
    
    % Find the pulse data
    if isBundled
        ext={bundle.oBundleItems(1:12).oExtension};
        idx=strcmp('.pul',ext);  %15.08.2012 - change from strmatch
        start=bundle.oBundleItems(idx).oStart;
    else    % Or open pulse file if not bundled
        fclose(fh);
        start=0;
        fh=fopen(fullfile(pathname,[filename,'.pul']),'r',endian);
    end
    
    fseek(fh,start,'bof');
    Magic=fread(fh,4,'uint8=>char')';
    if ~strcmpi('eerT',Magic); fclose(fh); disp('Pulse file wrong.');return;end
    
    Levels=fread(fh,1,'int32=>int32');
    Sizes=fread(fh,double(Levels),'int32=>int32');
    
    
    tree=getTree(fh,Sizes,ftell(fh));   % Get the tree form the pulse file
    
    tree{1,1}.RoBasicInfo.FilePath=pathname;
    tree{1,1}.RoBasicInfo.FileName=filename;
    tree{1,1}.RoBasicInfo.isLittleEndian=endian;
    tree{1,1}.RoBasicInfo.isBundled=isBundled;
    tree{1,1}.RoFileHeader=bundle;
    tree{1,1}.RawTree=tree;
    tree{1,1}=orderfields(tree{1,1});
    clear('Levels','Magic','Sizes','ans','filename','pathname','idx','start')
    
    %% Set offset for data
    if ~isBundled
        fclose(fh);
        fh=fopen(datafile,'r',endian);
    end
    clear('endian','ext')
    
    
    % Count channels.
    tree_chan=false(size(tree));
    for k=1:size(tree,1)
        for j=1:size(tree,2)
            if ~isempty(tree{k,j}); tree_chan(k,j)=true; end
        end
    end
    clear('k','j')
    
    % read every trace.
    for k=1:size(tree,1)
        if tree_chan(k,5)
            [format_raw,nbytes]=LocalFormatToString(tree{k,5}.TrDataFormat);
            readfmt=[format_raw '=>double'];
            fseek(fh,tree{k,5}.TrData,'bof');
            if tree{k,5}.TrInterleaveSize==0
                data(1:tree{k,5}.TrDataPoints)=fread(fh,double(tree{k,5}.TrDataPoints),readfmt);
                data=data*tree{k,5}.TrDataScaler;
                tree{k,5}.TrRawData=data';
            else
                offset=1;
                nelements= double(tree{k,5}.TrInterleaveSize/nbytes);
                data=zeros(tree{k,5}.TrDataPoints,1);
                for nread=1:floor(tree{k,5}.TrDataPoints/double(tree{k,5}.TrInterleaveSize/nbytes))
                    [data(offset:offset+nelements-1,1),N]=fread(fh,nelements,readfmt);
                    if (N<nelements)
                        disp('End of file reached unexpectedly');
                    end
                    offset=offset+nelements;
                    fseek(fh,double(tree{k,5}.TrInterleaveSkip-tree{k,5}.TrInterleaveSize),'cof');
                end
                data=data*tree{k,5}.TrDataScaler;
                tree{k,5}.TrRawData=data';
            end
            tree{k,5}.TrRecordingMode=patchType(tree{k,5}.TrRecordingMode);
            clear('data','format_raw','nbytes','readfmt')
        end
    end
    
    %% Close the file and return data.
    fclose(fh);
    
    varargout(1)={tree};
    varargout(2)={datafile};
    
end






%--------------------------------------------------------------------------
function [h,littleendianflag,isBundled]=getBundleHeader(fh)
    %--------------------------------------------------------------------------
    % Get the bundle header from a HEKA .dat file
    fseek(fh,0,'bof');
    h.oSignature=deblank(fread(fh,8,'uint8=>char')');
    switch h.oSignature
        case 'DATA'
            % Old format: nothing to do
            h.oVersion=[];
            h.oTime=[];
            h.oItems=[];
            h.oIsLittleEndian=[];
            h.oBundleItems(1:12)=[];
            h.BundleHeaderSize=0;
            isBundled=false;
        case {'DAT1' 'DAT2'}
            % Newer format
            h.oVersion=fread(fh,32,'uint8=>char')';
            h.oTime=fread(fh,1,'double');
            h.oItems=fread(fh,1,'int32=>int32');
            h.oIsLittleEndian=fread(fh,1,'uint8=>logical');
            h.BundleHeaderSize=256;
            switch h.oSignature
                case 'DAT1'
                    h.oBundleItems=[];
                    isBundled=false;
                case 'DAT2'
                    fseek(fh,64,'bof');
                    for k=1:12
                        h.oBundleItems(k).oStart=fread(fh,1,'int32=>int32');
                        h.oBundleItems(k).oLength=fread(fh,1,'int32=>int32');
                        h.oBundleItems(k).oExtension=deblank(fread(fh,8,'uint8=>char')');
                        h.oBundleItems(k).BundleItemSize=16;
                    end
                    isBundled=true;
            end
        otherwise
            error('This legacy file format is not supported');
    end
    littleendianflag=h.oIsLittleEndian;
    return
end

%--------------------------------------------------------------------------
function [Tree,Counter]=getTree(fh,Sizes,Position)
    %--------------------------------------------------------------------------
    % Main entry point for loading tree
    [Tree,Counter]=getTreeReentrant(fh,{},Sizes,0,Position,0);
    return
end

%--------------------------------------------------------------------------
function [Tree,Position,Counter]=getTreeReentrant(fh,Tree,Sizes,Level,Position,Counter)
    %--------------------------------------------------------------------------
    % Recursive routine called from LoadTree
    [Tree,Position,Counter,nchild]=getOneLevel(fh,Tree,Sizes,Level,Position,Counter);
    for k=1:double(nchild)
        [Tree,Position,Counter]=getTreeReentrant(fh,Tree,Sizes,Level+1,Position,Counter);
    end
    return
end

%--------------------------------------------------------------------------
function [Tree,Position,Counter,nchild]=getOneLevel(fh,Tree,Sizes,Level,Position,Counter)
    %--------------------------------------------------------------------------
    % Gets one record of the tree and the number of children
    [s, Counter]=getOneRecord(fh,Level,Counter);
    Tree{Counter,Level+1}=s;
    Position=Position+Sizes(Level+1);
    fseek(fh,Position,'bof');
    nchild=fread(fh,1,'int32=>int32');
    Position=ftell(fh);
    return
end

%--------------------------------------------------------------------------
function [rec, Counter]=getOneRecord(fh,Level,Counter)
    %--------------------------------------------------------------------------
    % Gets one record
    Counter=Counter+1;
    switch Level
        case 0
            rec=getRoot(fh);
        case 1
            rec=getGroupInHEKARead(fh);
        case 2
            rec=getSeries(fh);
        case 3
            rec=getSweep(fh);
        case 4
            rec=getTrace(fh);
        otherwise
            error('Unexpected Level');
    end
    return
end

% The functions below return data as defined by the HEKA PatchMaster
% specification

%--------------------------------------------------------------------------
function p=getRoot(fh)
    %--------------------------------------------------------------------------
    p.RoVersion=fread(fh,1,'int32=>int32');
    p.RoMark=fread(fh,1,'int32=>int32');%                         =   4; (* INT32 *)
    p.RoVersionName=deblank(fread(fh,32,'uint8=>char')');%        =   8; (* String32Type *)
    p.RoAuxFileName=deblank(fread(fh,80,'uint8=>char')');%        =  40; (* String80Type *)
    p.RoRootText=deblank(fread(fh,400,'uint8=>char')');% (* String400Type *)
    p.RoStartTime=fread(fh,1,'double=>double') ;%        = 520; (* LONGREAL *)
    p.RoStartTimeMATLAB=time2date(p.RoStartTime);
    p.RoMaxSamples=fread(fh,1,'int32=>int32'); %        = 528; (* INT32 *)
    p.RoCRC=fread(fh,1,'int32=>int32'); %                = 532; (* CARD32 *)
    p.RoFeatures=fread(fh,1,'int16=>int16'); %           = 536; (* SET16 *)
    p.RoFiller1=fread(fh,1,'int16=>int16');%         = 538; (* INT16 *)
    p.RoFiller2=fread(fh,1,'int32=>int32');%         = 540; (* INT32 *)
    p.RootRecSize= 544;
    p=orderfields(p);
    return
end

%--------------------------------------------------------------------------
function g=getGroupInHEKARead(fh)
    %--------------------------------------------------------------------------
    % Group
    g.GrMark=fread(fh,1,'int32=>int32');%               =   0; (* INT32 *)
    g.GrLabel=deblank(fread(fh,32,'uint8=>char')');%               =   4; (* String32Size *)
    g.GrText=deblank(fread(fh,80,'uint8=>char')');%                =  36; (* String80Size *)
    g.GrExperimentNumber=fread(fh,1,'int32=>int32');%   = 116; (* INT32 *)
    g.GrGroupCount=fread(fh,1,'int32=>int32');%         = 120; (* INT32 *)
    g.GrCRC=fread(fh,1,'int32=>int32');%                = 124; (* CARD32 *)
    g.GroupRecSize=128;%     (* = 16 * 8 *)
    g=orderfields(g);
    return
end

%--------------------------------------------------------------------------
function s=getSeries(fh)
    %--------------------------------------------------------------------------
    s.SeMark=fread(fh,1,'int32=>int32');%               =   0; (* INT32 *)
    s.SeLabel=deblank(fread(fh,32,'uint8=>char')');%              =   4; (* String32Type *)
    s.SeComment=deblank(fread(fh,80,'uint8=>char')');%            =  36; (* String80Type *)
    s.SeSeriesCount=fread(fh,1,'int32=>int32');%        = 116; (* INT32 *)
    s.SeNumbersw=fread(fh,1,'int32=>int32');%       = 120; (* INT32 *)
    s.SeAmplStateOffset=fread(fh,1,'int32=>int32');%    = 124; (* INT32 *)
    s.SeAmplStateSeries=fread(fh,1,'int32=>int32');%    = 128; (* INT32 *)
    s.SeSeriesType=fread(fh,1,'uint8=>uint8');%         = 132; (* BYTE *)
    
    % Added 15.08.2012
    s.SeUseXStart=logical(fread(fh,1,'uint8=>uint8'));%         = 133; (* BYTE *)
    
    s.SeFiller2=fread(fh,1,'uint8=>uint8');%         = 134; (* BYTE *)
    s.SeFiller3=fread(fh,1,'uint8=>uint8');%         = 135; (* BYTE *)
    s.SeTime=fread(fh,1,'double=>double') ;%               = 136; (* LONGREAL *)
    s.SeTimeMATLAB=time2date(s.SeTime);
    s.SePageWidth=fread(fh,1,'double=>double') ;%          = 144; (* LONGREAL *)
    for k=1:4
        s.SeSwUserParamDescr(k).Name=deblank(fread(fh,32,'uint8=>char')');%
        s.SeSwUserParamDescr(k).Unit=deblank(fread(fh,8,'uint8=>char')');%
    end
    s.SeFiller4=fread(fh,32,'uint8=>uint8');%         = 312; (* 32 BYTE *)
    s.SeSeUserParams=fread(fh,4,'double=>double');%       = 344; (* ARRAY[0..3] OF LONGREAL *)
    s.SeLockInParams=getSeLockInParams(fh);%       = 376; (* SeLockInSize = 96,see "Pulsed.de" *)
    s.SeAmplifierState=getAmplifierState(fh);%     = 472; (* AmplifierStateSize = 400 *)
    s.SeUsername=deblank(fread(fh,80,'uint8=>char')');%           = 872; (* String80Type *)
    for k=1:4
        s.SeSeUserParamDescr(k).Name=deblank(fread(fh,32,'uint8=>char')');% (* ARRAY[0..3] OF UserParamDescrType = 4*40 *)
        s.SeSeUserParamDescr(k).Unit=deblank(fread(fh,8,'uint8=>char')');%
    end
    s.SeFiller5=fread(fh,1,'int32=>int32');%         = 1112; (* INT32 *)
    s.SeCRC=fread(fh,1,'int32=>int32');%                = 1116; (* CARD32 *)
    
    % Added 15.08.2012
    s.SeSeUserParams2=fread(fh,4,'double=>double');
    for k=1:4
        s.SeSeUserParamDescr2(k).Name=deblank(fread(fh,32,'uint8=>char')');%
        s.SeSeUserParamDescr2(k).Unit=deblank(fread(fh,8,'uint8=>char')');%
    end
    s.SeScanParams=fread(fh,96,'uint8=>uint8');
    s.SeriesRecSize=1408;%      (* = 176 * 8 *)
    s=orderfields(s);
    return
end

%--------------------------------------------------------------------------
function sw=getSweep(fh)
    %--------------------------------------------------------------------------
    sw.SwMark=fread(fh,1,'int32=>int32');%               =   0; (* INT32 *)
    sw.SwLabel=deblank(fread(fh,32,'uint8=>char')');%              =   4; (* String32Type *)
    sw.SwAuxDataFileOffset=fread(fh,1,'int32=>int32');%  =  36; (* INT32 *)
    sw.SwStimCount=fread(fh,1,'int32=>int32');%          =  40; (* INT32 *)
    sw.SwSweepCount=fread(fh,1,'int32=>int32');%         =  44; (* INT32 *)
    sw.SwTime=fread(fh,1,'double=>double');%               =  48; (* LONGREAL *)
    sw.SwTimeMATLAB=time2date(sw.SwTime);% Also add in MATLAB datenum format
    sw.SwTimer=fread(fh,1,'double=>double');%              =  56; (* LONGREAL *)
    sw.SwSwUserParams=fread(fh,4,'double=>double');%       =  64; (* ARRAY[0..3] OF LONGREAL *)
    sw.SwTemperature=fread(fh,1,'double=>double');%        =  96; (* LONGREAL *)
    sw.SwOldIntSol=fread(fh,1,'int32=>int32');%          = 104; (* INT32 *)
    sw.SwOldExtSol=fread(fh,1,'int32=>int32');%          = 108; (* INT32 *)
    sw.SwDigitalIn=fread(fh,1,'int16=>int16');%          = 112; (* SET16 *)
    sw.SwSweepKind=fread(fh,1,'int16=>int16');%          = 114; (* SET16 *)
    sw.SwFiller1=fread(fh,1,'int32=>int32');%         = 116; (* INT32 *)
    sw.SwMarkers=fread(fh,4,'double=>double');%            = 120; (* ARRAY[0..3] OF LONGREAL *)
    sw.SwFiller2=fread(fh,1,'int32=>int32');%         = 152; (* INT32 *)
    sw.SwCRC=fread(fh,1,'int32=>int32');%                = 156; (* CARD32 *)
    sw.SweepRecSize         = 160;%      (* = 20 * 8 *)
    sw=orderfields(sw);
    return
end

%--------------------------------------------------------------------------
function tr=getTrace(fh)
    %--------------------------------------------------------------------------
    tr.TrMark=fread(fh,1,'int32=>int32');%               =   0; (* INT32 *)
    tr.TrLabel=deblank(fread(fh,32,'uint8=>char')');%              =   4; (* String32Type *)
    tr.TrTraceCount=fread(fh,1,'int32=>int32');%         =  36; (* INT32 *)
    tr.TrData=fread(fh,1,'int32=>int32');%               =  40; (* INT32 *)
    tr.TrDataPoints=fread(fh,1,'int32=>int32');%         =  44; (* INT32 *)
    tr.TrInternalSolution=fread(fh,1,'int32=>int32');%   =  48; (* INT32 *)
    tr.TrAverageCount=fread(fh,1,'int32=>int32');%       =  52; (* INT32 *)
    tr.TrLeakCount=fread(fh,1,'int32=>int32');%          =  56; (* INT32 *)
    tr.TrLeakTraces=fread(fh,1,'int32=>int32');%         =  60; (* INT32 *)
    tr.TrDataKind=fread(fh,1,'uint16=>uint16');%           =  64; (* SET16 *) NB Stored unsigned
    tr.TrFiller1=fread(fh,1,'int16=>int16');%         =  66; (* SET16 *)
    tr.TrRecordingMode=fread(fh,1,'uint8=>uint8');%      =  68; (* BYTE *)
    tr.TrAmplIndex=fread(fh,1,'uint8=>uint8');%          =  69; (* CHAR *)
    tr.TrDataFormat=fread(fh,1,'uint8=>uint8');%         =  70; (* BYTE *)
    tr.TrDataAbscissa=fread(fh,1,'uint8=>uint8');%       =  71; (* BYTE *)
    tr.TrDataScaler=fread(fh,1,'double=>double');%         =  72; (* LONGREAL *)
    tr.TrTimeOffset=fread(fh,1,'double=>double');%         =  80; (* LONGREAL *)
    tr.TrZeroData=fread(fh,1,'double=>double');%           =  88; (* LONGREAL *)
    tr.TrYUnit=deblank(fread(fh,8,'uint8=>char')');%              =  96; (* String8Type *)
    tr.TrXInterval=fread(fh,1,'double=>double');%          = 104; (* LONGREAL *)
    tr.TrXStart=fread(fh,1,'double=>double');%             = 112; (* LONGREAL *)
    % 17.04.10 TrXUnit bytes may include some trailing characters after NULL
    % byte
    tr.TrXUnit=deblank(fread(fh,8,'uint8=>char')');%              = 120; (* String8Type *)
    tr.TrYRange=fread(fh,1,'double=>double');%             = 128; (* LONGREAL *)
    tr.TrYOffset=fread(fh,1,'double=>double');%            = 136; (* LONGREAL *)
    tr.TrBandwidth=fread(fh,1,'double=>double');%          = 144; (* LONGREAL *)
    tr.TrPipetteResistance=fread(fh,1,'double=>double');%  = 152; (* LONGREAL *)
    tr.TrCellPotential=fread(fh,1,'double=>double');%      = 160; (* LONGREAL *)
    tr.TrSealResistance=fread(fh,1,'double=>double');%     = 168; (* LONGREAL *)
    tr.TrCSlow=fread(fh,1,'double=>double');%              = 176; (* LONGREAL *)
    tr.TrGSeries=fread(fh,1,'double=>double');%            = 184; (* LONGREAL *)
    tr.TrRsValue=fread(fh,1,'double=>double');%            = 192; (* LONGREAL *)
    tr.TrGLeak=fread(fh,1,'double=>double');%              = 200; (* LONGREAL *)
    tr.TrMConductance=fread(fh,1,'double=>double');%       = 208; (* LONGREAL *)
    tr.TrLinkDAChannel=fread(fh,1,'int32=>int32');%      = 216; (* INT32 *)
    tr.TrValidYrange=fread(fh,1,'uint8=>logical');%        = 220; (* BOOLEAN *)
    tr.TrAdcMode=fread(fh,1,'uint8=>uint8');%            = 221; (* CHAR *)
    tr.TrAdcChannel=fread(fh,1,'int16=>int16');%         = 222; (* INT16 *)
    tr.TrYmin=fread(fh,1,'double=>double');%               = 224; (* LONGREAL *)
    tr.TrYmax=fread(fh,1,'double=>double');%               = 232; (* LONGREAL *)
    tr.TrSourceChannel=fread(fh,1,'int32=>int32');%      = 240; (* INT32 *)
    tr.TrExternalSolution=fread(fh,1,'int32=>int32');%   = 244; (* INT32 *)
    tr.TrCM=fread(fh,1,'double=>double');%                 = 248; (* LONGREAL *)
    tr.TrGM=fread(fh,1,'double=>double');%                 = 256; (* LONGREAL *)
    tr.TrPhase=fread(fh,1,'double=>double');%              = 264; (* LONGREAL *)
    tr.TrDataCRC=fread(fh,1,'int32=>int32');%            = 272; (* CARD32 *)
    tr.TrCRC=fread(fh,1,'int32=>int32');%                = 276; (* CARD32 *)
    tr.TrGS=fread(fh,1,'double=>double');%                 = 280; (* LONGREAL *)
    tr.TrSelfChannel=fread(fh,1,'int32=>int32');%        = 288; (* INT32 *)
    
    % Added 15.08.2012
    tr.TrInterleaveSize=fread(fh,1,'int32=>int32');%        = 292; (* INT32 *)
    tr.TrInterleaveSkip=fread(fh,1,'int32=>int32');%        = 296; (* INT32 *)
    tr.TrImageIndex=fread(fh,1,'int32=>int32');%        = 300; (* INT32 *)
    tr.TrMarkers=fread(fh,10,'double=>double');%        = 304; (* ARRAY[0..9] OF LONGREAL *)
    tr.TrSECM_X=fread(fh,1,'double=>double');%        = 384; (* LONGREAL *)
    tr.TrSECM_Y=fread(fh,1,'double=>double');%        = 392; (* LONGREAL *)
    tr.TrSECM_Z=fread(fh,1,'double=>double');%        = 400; (* LONGREAL *)
    tr.TraceRecSize=408;
    
    tr=orderfields(tr);
    return
end

%--------------------------------------------------------------------------
function L=getSeLockInParams(fh)
    %--------------------------------------------------------------------------
    offset=ftell(fh);
    L.loExtCalPhase=fread(fh,1,'double=>double') ;%        =   0; (* LONGREAL *)
    L.loExtCalAtten=fread(fh,1,'double=>double') ;%        =   8; (* LONGREAL *)
    L.loPLPhase=fread(fh,1,'double=>double') ;%            =  16; (* LONGREAL *)
    L.loPLPhaseY1=fread(fh,1,'double=>double') ;%          =  24; (* LONGREAL *)
    L.loPLPhaseY2=fread(fh,1,'double=>double') ;%          =  32; (* LONGREAL *)
    L.loUsedPhaseShift=fread(fh,1,'double=>double') ;%     =  40; (* LONGREAL *)
    L.loUsedAttenuation=fread(fh,1,'double=>double');%    =  48; (* LONGREAL *)
    skip=fread(fh,1,'double=>double'); %#ok<NASGU>
    L.loExtCalValid=fread(fh,1,'uint8=>logical') ;%        =  64; (* BOOLEAN *)
    L.loPLPhaseValid=fread(fh,1,'uint8=>logical') ;%       =  65; (* BOOLEAN *)
    L.loLockInMode=fread(fh,1,'uint8=>uint8') ;%         =  66; (* BYTE *)
    L.loCalMode=fread(fh,1,'uint8=>uint8') ;%            =  67; (* BYTE *)
    L.LockInParamsSize=96;
    fseek(fh,offset+L.LockInParamsSize,'bof');
    return
end

%--------------------------------------------------------------------------
function A=getAmplifierState(fh)
    %--------------------------------------------------------------------------
    offset=ftell(fh);
    A.E9StateVersion=fread(fh,1,'double=>double');%       =   0; (* 8 = SizeStateVersion *)
    A.E9RealCurrentGain=fread(fh,1,'double=>double');%    =   8; (* LONGREAL *)
    A.E9RealF2Bandwidth=fread(fh,1,'double=>double');%    =  16; (* LONGREAL *)
    A.E9F2Frequency=fread(fh,1,'double=>double');%        =  24; (* LONGREAL *)
    A.E9RsValue=fread(fh,1,'double=>double');%            =  32; (* LONGREAL *)
    A.E9RsFraction=fread(fh,1,'double=>double');%         =  40; (* LONGREAL *)
    A.E9GLeak=fread(fh,1,'double=>double');%              =  48; (* LONGREAL *)
    A.E9CFastAmp1=fread(fh,1,'double=>double');%          =  56; (* LONGREAL *)
    A.E9CFastAmp2=fread(fh,1,'double=>double');%          =  64; (* LONGREAL *)
    A.E9CFastTau=fread(fh,1,'double=>double');%           =  72; (* LONGREAL *)
    A.E9CSlow=fread(fh,1,'double=>double');%              =  80; (* LONGREAL *)
    A.E9GSeries=fread(fh,1,'double=>double');%            =  88; (* LONGREAL *)
    A.E9StimDacScale=fread(fh,1,'double=>double');%       =  96; (* LONGREAL *)
    A.E9CCStimScale=fread(fh,1,'double=>double');%        = 104; (* LONGREAL *)
    A.E9VHold=fread(fh,1,'double=>double');%              = 112; (* LONGREAL *)
    A.E9LastVHold=fread(fh,1,'double=>double');%          = 120; (* LONGREAL *)
    A.E9VpOffset=fread(fh,1,'double=>double');%           = 128; (* LONGREAL *)
    A.E9VLiquidJunction=fread(fh,1,'double=>double');%    = 136; (* LONGREAL *)
    A.E9CCIHold=fread(fh,1,'double=>double');%            = 144; (* LONGREAL *)
    A.E9CSlowStimVolts=fread(fh,1,'double=>double');%     = 152; (* LONGREAL *)
    A.E9CCtr.TrackVHold=fread(fh,1,'double=>double');%       = 160; (* LONGREAL *)
    A.E9TimeoutLength=fread(fh,1,'double=>double');%      = 168; (* LONGREAL *)
    A.E9SearchDelay=fread(fh,1,'double=>double');%        = 176; (* LONGREAL *)
    A.E9MConductance=fread(fh,1,'double=>double');%       = 184; (* LONGREAL *)
    A.E9MCapacitance=fread(fh,1,'double=>double');%       = 192; (* LONGREAL *)
    A.E9SerialNumber=fread(fh,1,'double=>double');%       = 200; (* 8 = SizeSerialNumber *)
    A.E9E9Boards=fread(fh,1,'int16=>int16');%           = 208; (* INT16 *)
    A.E9CSlowCycles=fread(fh,1,'int16=>int16');%        = 210; (* INT16 *)
    A.E9IMonAdc=fread(fh,1,'int16=>int16');%            = 212; (* INT16 *)
    A.E9VMonAdc=fread(fh,1,'int16=>int16');%            = 214; (* INT16 *)
    A.E9MuxAdc=fread(fh,1,'int16=>int16');%             = 216; (* INT16 *)
    A.E9TstDac=fread(fh,1,'int16=>int16');%             = 218; (* INT16 *)
    A.E9StimDac=fread(fh,1,'int16=>int16');%            = 220; (* INT16 *)
    A.E9StimDacOffset=fread(fh,1,'int16=>int16');%      = 222; (* INT16 *)
    A.E9MaxDigitalBit=fread(fh,1,'int16=>int16');%      = 224; (* INT16 *)
    A.E9SpareInt1=fread(fh,1,'int16=>int16');%       = 226; (* INT16 *)
    A.E9SpareInt2=fread(fh,1,'int16=>int16');%       = 228; (* INT16 *)
    A.E9SpareInt3=fread(fh,1,'int16=>int16');%       = 230; (* INT16 *)
    
    A.E9AmplKind=fread(fh,1,'uint8=>uint8');%           = 232; (* BYTE *)
    A.E9IsEpc9N=fread(fh,1,'uint8=>uint8');%            = 233; (* BYTE *)
    A.E9ADBoard=fread(fh,1,'uint8=>uint8');%            = 234; (* BYTE *)
    A.E9BoardVersion=fread(fh,1,'uint8=>uint8');%       = 235; (* BYTE *)
    A.E9ActiveE9Board=fread(fh,1,'uint8=>uint8');%      = 236; (* BYTE *)
    A.E9Mode=fread(fh,1,'uint8=>uint8');%               = 237; (* BYTE *)
    A.E9Range=fread(fh,1,'uint8=>uint8');%              = 238; (* BYTE *)
    A.E9F2Response=fread(fh,1,'uint8=>uint8');%         = 239; (* BYTE *)
    
    A.E9RsOn=fread(fh,1,'uint8=>uint8');%               = 240; (* BYTE *)
    A.E9CSlowRange=fread(fh,1,'uint8=>uint8');%         = 241; (* BYTE *)
    A.E9CCRange=fread(fh,1,'uint8=>uint8');%            = 242; (* BYTE *)
    A.E9CCGain=fread(fh,1,'uint8=>uint8');%             = 243; (* BYTE *)
    A.E9CSlowToTstDac=fread(fh,1,'uint8=>uint8');%      = 244; (* BYTE *)
    A.E9StimPath=fread(fh,1,'uint8=>uint8');%           = 245; (* BYTE *)
    A.E9CCtr.TrackTau=fread(fh,1,'uint8=>uint8');%         = 246; (* BYTE *)
    A.E9WasClipping=fread(fh,1,'uint8=>uint8');%        = 247; (* BYTE *)
    
    A.E9RepetitiveCSlow=fread(fh,1,'uint8=>uint8');%    = 248; (* BYTE *)
    A.E9LastCSlowRange=fread(fh,1,'uint8=>uint8');%     = 249; (* BYTE *)
    A.E9Locked=fread(fh,1,'uint8=>uint8');%             = 250; (* BYTE *)
    A.E9CanCCFast=fread(fh,1,'uint8=>uint8');%          = 251; (* BYTE *)
    A.E9CanLowCCRange=fread(fh,1,'uint8=>uint8');%      = 252; (* BYTE *)
    A.E9CanHighCCRange=fread(fh,1,'uint8=>uint8');%     = 253; (* BYTE *)
    A.E9CanCCtr.Tracking=fread(fh,1,'uint8=>uint8');%      = 254; (* BYTE *)
    A.E9HasVmonPath=fread(fh,1,'uint8=>uint8');%        = 255; (* BYTE *)
    
    A.E9HasNewCCMode=fread(fh,1,'uint8=>uint8');%       = 256; (* BYTE *)
    A.E9Selector=fread(fh,1,'uint8=>char');%           = 257; (* CHAR *)
    A.E9HoldInverted=fread(fh,1,'uint8=>uint8');%       = 258; (* BYTE *)
    A.E9AutoCFast=fread(fh,1,'uint8=>uint8');%          = 259; (* BYTE *)
    A.E9AutoCSlow=fread(fh,1,'uint8=>uint8');%          = 260; (* BYTE *)
    A.E9HasVmonX100=fread(fh,1,'uint8=>uint8');%        = 261; (* BYTE *)
    A.E9TestDacOn=fread(fh,1,'uint8=>uint8');%          = 262; (* BYTE *)
    A.E9QMuxAdcOn=fread(fh,1,'uint8=>uint8');%          = 263; (* BYTE *)
    
    A.E9RealImon1Bandwidth=fread(fh,1,'double=>double');% = 264; (* LONGREAL *)
    A.E9StimScale=fread(fh,1,'double=>double');%          = 272; (* LONGREAL *)
    
    A.E9Gain=fread(fh,1,'uint8=>uint8');%               = 280; (* BYTE *)
    A.E9Filter1=fread(fh,1,'uint8=>uint8');%            = 281; (* BYTE *)
    A.E9StimFilterOn=fread(fh,1,'uint8=>uint8');%       = 282; (* BYTE *)
    A.E9RsSlow=fread(fh,1,'uint8=>uint8');%             = 283; (* BYTE *)
    A.E9Old1=fread(fh,1,'uint8=>uint8');%            = 284; (* BYTE *)
    A.E9CCCFastOn=fread(fh,1,'uint8=>uint8');%          = 285; (* BYTE *)
    A.E9CCFastSpeed=fread(fh,1,'uint8=>uint8');%        = 286; (* BYTE *)
    A.E9F2Source=fread(fh,1,'uint8=>uint8');%           = 287; (* BYTE *)
    
    A.E9TestRange=fread(fh,1,'uint8=>uint8');%          = 288; (* BYTE *)
    A.E9TestDacPath=fread(fh,1,'uint8=>uint8');%        = 289; (* BYTE *)
    A.E9MuxChannel=fread(fh,1,'uint8=>uint8');%         = 290; (* BYTE *)
    A.E9MuxGain64=fread(fh,1,'uint8=>uint8');%          = 291; (* BYTE *)
    A.E9VmonX100=fread(fh,1,'uint8=>uint8');%           = 292; (* BYTE *)
    A.E9IsQuadro=fread(fh,1,'uint8=>uint8');%           = 293; (* BYTE *)
    A.E9SpareBool4=fread(fh,1,'uint8=>uint8');%      = 294; (* BYTE *)
    A.E9SpareBool5=fread(fh,1,'uint8=>uint8');%      = 295; (* BYTE *)
    
    A.E9StimFilterHz=fread(fh,1,'double=>double');%       = 296; (* LONGREAL *)
    A.E9RsTau=fread(fh,1,'double=>double');%              = 304; (* LONGREAL *)
    A.E9FilterOffsetDac=fread(fh,1,'int16=>int16');%    = 312; (* INT16 *)
    A.E9ReferenceDac=fread(fh,1,'int16=>int16');%       = 314; (* INT16 *)
    A.E9SpareInt6=fread(fh,1,'int16=>int16');%       = 316; (* INT16 *)
    A.E9SpareInt7=fread(fh,1,'int16=>int16');%       = 318; (* INT16 *)
    A.E9Spares1=320;
    
    A.E9CalibDate=fread(fh,2,'double=>double');%          = 344; (* 16 = SizeCalibDate *)
    A.E9SelHold=fread(fh,1,'double=>double');%            = 360; (* LONGREAL *)
    A.AmplifierStateSize   = 400;
    fseek(fh,offset+A.AmplifierStateSize,'bof');
    return
end






%--------------------------------------------------------------------------
function [fmt,nbytes]=LocalFormatToString(n)
    %--------------------------------------------------------------------------
    switch n
        case 0
            fmt='int16';
            nbytes=2;
        case 1
            fmt='int32';
            nbytes=4;
        case 2
            fmt='single';
            nbytes=4;
        case 3
            fmt='double';
            nbytes=8;
    end
    return
end


%--------------------------------------------------------------------------
function str=time2date(t)
    %--------------------------------------------------------------------------
    t=t-1580970496;
    if t<0
        t=t+4294967296;
    end
    t=t+9561652096;
    str=datestr(t/(24*60*60)+datenum(1601,1,1));
    return
end
%--------------------------------------------------------------------------

function str=patchType(n)
    switch n
        case 0
            str='Inside-out';
        case 1
            str='Cell-attached';
        case 2
            str='Outside-out';
        case 3
            str='Whole=cell';
        case 4
            str='Current-lamp';
        case 5
            str='Voltage-clamp';
        otherwise
            str='External/Unknown';
    end
    return
end

