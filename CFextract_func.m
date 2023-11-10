function [CF,thr,Q10,msp] = CFextract_func(curvedata, curvesettings)
%% 
% input:
%   curvedata: a curvedata struct, containing the spike times
%   curvesettings: a curvesettings struct, containing the recording metadata
% output:
%   CF: the estimated characteristic frequency in Hz, based on the tuning 
%       curve, that is visually checked, i.e. the frequency with the lowest
%       threshold in the tuning curve
%   thr: the threshold at CF
%   Q10: the Q10dB value, calculated by the CF/bandwidth of the tuning
%       curve at 10 dB above threshold
%   msp: the spontaneous rate in spikes/s, based on the firing rate during 
%       the silent trials

% By: Go Ashida and Amarins Heeringa

%%
[tlen,StartTime,EndTime,nwSpikCnt,nwSpikRat]=getSpikeCounts(curvedata,curvesettings);

% loop variables
x1 = curvedata.depvars_sort(:,1,1);
x2 = curvedata.depvars_sort(:,1,2);

% calculate mean and std of spike rates 
% y = curvedata.spike_counts * 1000 / tlen;
y = nwSpikRat;
m1 = mean(y,2);  % average wrt 2nd index
s1 = std(y,1,2); % std wrt 2nd index 

% data for non-spont trials
j = (curvedata.isspont(:,1)==0); % find out non-spont trials

% data for spont trials
jsp = (curvedata.isspont(:,1)==1); % find out spont trials
msp = m1(jsp);
ssp = s1(jsp);

% only procede if there are two looping variables, frequency and level
if curvesettings.stimcache.nloopvars == 2 

    % get dependent variables (without spont)
    v1 = sort(unique(x1(j)));
    v2 = sort(unique(x2(j)));
    
    % matrix to store mean rates
    m = zeros(length(v1),length(v2));
    s = zeros(length(v1),length(v2));
    for i2 = 1:length(v2)
        m(:,i2) = m1( x2==v2(i2) );
        s(:,i2) = s1( x2==v2(i2) );
    end
    
    % Tuning curve: a response is mean sfr + 1.2x std sfr
    resp=msp+1.2*ssp;
    
    done=NaN;
    while isnan(done)
    tc=zeros(1,length(v1)); % pre-allocate the tuning curve (tc)
    for i=1:length(v1)
        xi=find(m(i,:)>resp & m(i,:)>15);
        if ~isempty(xi) && length(xi)>1
            xn=find(diff(xi)==1,1,'first');
            if ~isempty(xn)
            tc(1,i)=v2(xi(xn));
            else
                tc(1,i)=NaN;
            end
        elseif length(xi)==1
            tc(1,i)=v2(xi);
        else
            tc(1,i)=NaN;
        end
    end
    [thr,CF,Q10,indx]=getQ10dB(tc,v1);
    CFs=v1(tc==thr);
    
    % Rate-level function
    if length(CFs)>1
        RLF=zeros(length(v2),length(CFs));
        for i=1:length(CFs)
            RLF(:,i)=m(v1==CFs(i),:);
        end
        RLF=mean(RLF,2);
    else
        RLF=m(v1==CF,:);
    end
    
    % Plot it
    figure;
    subplot(2,1,1)
    surf(v1,v2,m')
    view(2)
    xlabel('Frequency (Hz)')
    ylabel('Level (dB SPL)')
    h=colorbar;
    ylabel(h,'Rate (spikes/s)')
    xlim([min(v1) max(v1)]);
    ylim([min(v2) max(v2)]);
    title('Receptive field')
    
    h1 = subplot(2,2,3);
    plot(v1,tc,'k')
    xlim([min(v1) max(v1)]);
    ylim([min(v2) max(v2)]);
    hold on
    line([CF CF],ylim,'Color','m','LineStyle','--')
    if ~isnan(Q10)
    line(v1(indx),[thr+10 thr+10],'Color','r','LineStyle','--')
    end
    xlabel('Frequency (Hz)')
    ylabel('Level (dB SPL)')
    title('Tuning curve')
    text(CF+100,max(v2)-10,sprintf('CF = %d',CF))
    text(CF+100,max(v2)-20,sprintf('Q10dB = %0.2f',Q10))
    
    h2 = subplot(2,2,4);
    if ~isempty(RLF)
        plot(v2,RLF,'k')
        hold on
        errorbar(min(v2)-10,msp,ssp,'or','MarkerEdgeColor','k','MarkerFaceColor','r')
        line([thr thr],ylim,'Color','b','LineStyle','--')
        xlabel('Intensity (dB SPL)')
        ylabel('Spike rate (sp/s)')
        title('RLF at CF')
        text(thr+10,msp,sprintf('Thr = %d',thr))
    end
    
    ok=input('does this look ok? (y/n) ','s');
    if strcmp(ok,'y')
        done=1;
    else
        fprintf('The response threshold was %0.2f\n',resp);
        resp=input('What is the new response threshold: \n');
        close gcf
    end
    end
    
    % Check if the CF was correctly determined, e.g. when the tail
    % frequency was taken as CF, you can correct this here.
    done2=NaN;
    while isnan(done2)
        ok2=input('does the cf look ok? (y/n) ','s');
        if strcmp(ok2,'y')
            done2=1;
        else
            disp('These are the frequencies used: ')
            disp(v1)
            CF=input('What is the new cf: \n');
            [thr,Q10,indx]=getQ10dB_withCF(tc,v1,CF);
            cla(h1)
            h1 = subplot(2,2,3);
            plot(v1,tc,'k')
            xlim([min(v1) max(v1)]);
            ylim([min(v2) max(v2)]);
            hold on
            line([CF CF],ylim,'Color','m','LineStyle','--')
            if ~isnan(Q10)
            line(v1(indx),[thr+10 thr+10],'Color','r','LineStyle','--')
            end
            xlabel('Frequency (Hz)')
            ylabel('Level (dB SPL)')
            title('Tuning curve')
            text(CF+100,max(v2)-10,sprintf('CF = %d',CF))
            text(CF+100,max(v2)-20,sprintf('Q10dB = %0.2f',Q10))
            
            RLF=m(v1==CF,:);
            cla(h2)
            h2 = subplot(2,2,4);
            if ~isempty(RLF)
                plot(v2,RLF,'k')
                hold on
                errorbar(min(v2)-10,msp,ssp,'or','MarkerEdgeColor','k','MarkerFaceColor','r')
                line([thr thr],ylim,'Color','b','LineStyle','--')
                xlabel('Level (dB SPL)')
                ylabel('Rate (spikes/s)')
                title('RLF at CF')
                text(thr+10,msp,sprintf('Thr = %d',thr))
            end
        end

    end
    
% if there are not two looping variables
else 
    disp('Not a valid CF file')
    CF = NaN; thr = NaN; Q10 = NaN; msp = NaN;
end
end