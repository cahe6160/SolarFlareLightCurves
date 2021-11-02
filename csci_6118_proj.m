%Cole Tamburri, 2020, University of Colorado at Boulder APS Department
%
%Performs impulsivity calculations for solar flares in the 304A line as
%recorded by SDO/EVE MEGS-B Level 2 data.  The input data from this
%instrument are downloaded using "downloadsdoevefit.m" as .fit files, and 
%then processed for use by calculating solar quiet values with fft304.m
%(very simple, only one line, which filters the data over 5 hours in order
%to remove flare signatures).  Flares are identified according to the
%ribbondb database (http://solarmuri.ssl.berkeley.edu/~kazachenko/RibbonDB/)
%and the script can be changed to study any of the 2048 flares within this
%database.  

%Start time, peak irradiance, and end time values for the flares are
%calculated in order to obtain values for peak irradiance and full width at
%half height in time.  In addition, models are fit to the rise and decay
%periods of each flare in order to improve the temporal resolution of the
%light curve and obtain a more functional value for FWHH in time (without
%this, the particularly high-impulsivity flares would be subject to error
%as a result of the large differences in irradiance value between
%consecutive time steps, leading to error in FWHH

%An effort is also made to filter out and correct for low SNR flares.

%clear workspace
clear;
starti=100000; %initialize to something outside range
%load the file sdoeve304, which containts the data for the 304A line as
%recorded by SDO/EVE for the extent of the MEGS-B experiment, between 2010
%and 2014
filename='/Users/owner/Desktop/csci_proj_arr.mat';
load(filename);
mkdir('/Users/owner/Desktop/fake/');
ribbondb_info = readtable('/Users/owner/Desktop/ribbondb_v1.0.csv');
%initialize arrays
starttimes=NaN(2048,1);
endtimes=NaN(2048,1);
maxtimes=NaN(2048,1);
curly_Is = zeros(2048,1);
decay_gof=cell(2048,1);
impulse_gof=cell(2048,1);
decay_chi = cell(2048,1);
decay_pval = cell(2048,1);
decay_chistats = cell(2048,1);
rsqurs_exp2_four2 = NaN(2048,2);
rsqurs_3_four2 = NaN(2048,2);
decay_chistats_byhand = cell(2048,1);
errs = cell(2048,1);
chisquareds = cell(2049,2);
chisquareds{1}={'Chi-Squared, sdo/eve error','Reduced Chi-Squared, sdo/eve error','Reduced Chi-Squared, verification, sdo/eve error','Reduced Chi-Squared, errors from sq'};
significance = 0.05;
endj=NaN;
tst=0;
rise_chi = cell(2048,1);
rise_pval = cell(2048,1);
rise_chistats = cell(2048,1);
rise_chistats_byhand = cell(2048,1);


% %initialize possible functions
% modelfun_atan = @(b,x)b(1) + b(2)*atan(b(3)*x(:,1));
% modelfun_cuberoot = @(b,x) b(1)+b(2)*(x(:,1)).^(1/3);
% %0.5*sqrt(pi)*a*c*exp{[(d*(b-t))+(c*c*d*d)/4]}*{erf(z)-erf[z-(x/c)]};
% %modelfun_gausspwr = @(b,x) b(1)*exp(-((x-b(2))/b(3))^2)*
% convolution = '0.5*sqrt(pi)*a*c*exp((d*(b-x)+(c*c*d*d)/4))*(erf(z)-erf(z-(x/c)))';

%initialize figures
f=figure(1)
clf
set(gcf,'Position',[100 100 1500 1500])
%initialize certain indices

%use i array to decide which flares in ribbondb are to be used


for i=1:2049

    i
    

    sqerrstd=[];
    sqerrstd2=[];
    sqerrstd3=[];
    %"windowstart2" is the beginning of the window as defined by ribbondb
    %define our window to be one hour before "windowstart2" to 5 hours after
    %that point (essentially just expand window, having downloaded all
    %SDO/EVE data
    wst=windowstart2(i);
    wstn=windowstart2(i+1);
    %wst - (1.5/24 works for flares 1-25)
    wstb=wst-(2/24);
    wend=wst+(3.5/24);
    
    %find the indices of times in vectime (from SDO/EVE) corresponding to
    %this window
    eventi=find(vectime>wstb & vectime<wend);
    eventi_halfsolarrot=find(vectime>wst-(12/24) & vectime<wend+(12/24));
    
    %assume flares do not overlap
    if wend>wstn
        eventi=find(vectime>wstb & vectime<wstn);
        eventi_halfsolarrot=find(vectime>wst-(12/24) & vectime<(wstn+(12/24)));
    end
    %extract irradiance/time/sq/difference between irradiance and sq
    %values for this event in our defined window
    
    irrev=vec304(eventi);
    
    inan=find(isnan(irrev));
        timeev=vectime(eventi);

    
   

    sqev=filt304(eventi);
    
    errev=vecerr(eventi);
    
    
    

    passage = vec304(eventi_halfsolarrot);
    
    passaget = vectime(eventi_halfsolarrot);
    %standard deviations of both solar quiet and irradiance values,
    %ignoring NaN
    sqstd=nanstd(sqev);
    irrstd=nanstd(irrev);
    %irrstd=nanstd(passage);
    irrmean=nanmean(passage);
    diff=irrev-sqev;
 
    for p=2:length(irrev)-1
        if irrev(p)>irrev(p-1)+2*irrstd & irrev(p)>irrev(p+1)+2*irrstd 
            irrev(p)=NaN;
        end
    end
 
    for p=2:length(irrev)-1
        if irrev(p)<irrev(p-1)-irrstd  
            irrev(p)=NaN;
        end
    end    

          
            
    %N.B.: use solar quiet standard deviation to avoid picking up later flares
    %which might be included in the same window (use irrsq, and later,
    %larger flares might skew the std so that the first, less strong flare
    %might be ignored by the calculation)
    
    %find flares with low SNR - NEEDS WORK
%     diffmaxi=find(diff==max(diff));
%     avgdiffma=mean(diff(diffmaxi-5:diffmaxi+5));
%     if avgdiffma<1.5*abs(min(diff))
%         flag=1
%         i
%     end
    
    nans=isnan(irrev);
    
        
    %now begin logic to find the start time of the flare
    %(1) ideal start time scenario - find 40 points which lie consecutively
    %above the solar quiet by 6 times the solar quiet standard deviation;
    %if this is sufficiently far into the window, subtract 80 from that
    %index to find the start time. If this is not sufficiently far into the
    %window, simply subtract 39 and the start time will be roughly the
    %beginning of the window.
    %(2) and (3) perform the same process, but lower the criteria to 30
    %points consecutively above 6*sqstd in (2) and 20 points consecutively
    %above 6*sqstd.  Always subtract an extra 40 points when moving back to
    %find the calculated start time
    %
    n=0;
    m=0;
    lastj=0;
    tst=0;
    for j=1:length(irrev)
        if diff(j)>(2*irrstd)
            n=n+1;
            if n == 1;
                lastj=j;
            end
            if n==40 && j>n+40 
                lastj=j;
                tst=timeev(j-80);  
                starti=j-80;
                1
                break
            end
            if n==40 && j<n+40 
                lastj=j;
                tst=timeev(j-39); 
                starti=j-39;
                2
                break  
            end
          
        end 
    end
   
    %(2)
    n=0;
    m=0;
    if tst==0 || starti >1500
        for j=1:length(irrev)
            if diff(j)>(2*irrstd)
                n=n+1;
                if n==30 && j>n+30
                    tst=timeev(j-60);
                    starti=j-60;
                    3
                    break  
                end
                if n==30 && j<n+30
                    tst=timeev(j-29); 
                    starti=j-29;
                    4
                    break
                end
            end
        end
    end
    
    %(3)
    n=0;
    m=0;
    if tst==0 || starti >1500
        for j=1:length(irrev)
            if diff(j)>(2*irrstd)
                n=n+1;
                if n==20 && j>n+20
                    tst=timeev(j-40); 
                    starti=j-40;
                    5
                    break
                end
                if n==20 && j<n+20
                    tst=timeev(j-19);  
                    starti=j-19;
                    6
                    break
                end
            end
        end
    end
    
    %(4)
    
    if tst==0 || starti >1500
        for j=1:length(irrev)
            if diff(j)>(2*irrstd)
                n=n+1;
                if n==10 && j>n+10
                    tst=timeev(j-20); 
                    starti=j-20;
                    7
                    break
                end
                if n==10 && j<n+10
                    tst=timeev(j-9);  
                    starti=j-9;
                    8
                    break
                end
            end
        end
    end
    


    %at this point, j will be the index at which the criteria (whether
    %in 1, 2, or 3) was satisfied.  Carry this along: 
    startj = j;
    
    %(5)-(7) account for if the found start time was already counted as an 
    %event. If it has been, it will have been stored in "starttimes" below already.
    %Then start from startj (the point at which the criteria was met for
    %the first start time) and move forwards.  For (4) through (6), simply
    %repeat the steps above, but start after the end time of the last
    %flare (calculated below as endj).  (5) and (6) lower
    %"number-of-points" criteria in the same way as (2) and (3)
    if starti > length(timeev) || endj > length(timeev)
        'No classification of flare!'
        
    
    %if the time of start is actually below the first element in
    %timeev (the window of time being studied), the start time will lie
    %outside.  Skip these events for now, still plotting solar quiet and 
    %irradiance, but not start or end times - but will need to be fixed later.
    elseif tst<timeev(1) 
        
        'Start too early!'
        
    else
        tst=timeev(starti);
        %Now the fun bit.  Find irradiance peak for well-behaved data,
        %first determining the next time that a light curve returns below
        %the solar quiet for a certain number of consecutive points (50, in
        %this case)

        %N.B. this first piece of logic only works if the start time 
        %was found using the first condition! (and doesn't work great if 
        %the increase in light curve is so fast that you're quickly on the decay phase)
        m=0;
        for j=(starti+40):length(diff)
            if diff(j)<0
                m=m+1;
                if m==50
                    endj=j;
                    
                    break
                end
            end
        end
        if isnan(endj)
            disp('endj is NaN')
        else
        
        %find the tend, the end time for the event, by identifying the
        %position within timeev (the window of time for the event being
        %studied) at which endj lies
        tend=timeev(endj);
        
        %now the window of the flare is consolidated into "window,"
        %consisting of the irradiance values between startj and endj
        window=irrev(starti:endj);
        windowt=timeev(starti:endj);
        sqwindow = sqev(starti:endj);
        sq=sqev(1:starti);
                
        %find the average irradiance values within a 15-point spread around
        %each point (greater than 7 positions from the start and less than
        %7 positions from the end, of course) - this will give a better
        %idea of the actual average at the peak
        avgmeans=NaN(length(window)-7,1);
        for k=8:length(window)-7
            avgmeans(k-7)=nanmean(window(k-7:k+7));
        end
        
        %find the maximum in this average of the window
        maxind1=find(avgmeans==max(avgmeans));
        
        %find where the maximum lies in the original window, calculated at
        %the beginning of studying the flare
        maxind=starti+maxind1;
        
        %find the time of the maximum
        maxt=timeev(maxind);
        end
    end
    
    %find if the start time of this event is before the end time of the
    %last or the flare is really after this calculated window
    if i>1 
        if tst < endtimes(i-1) %|| max(irrev(endj:end))-max(window) > 2*irrstd
            
            %reset tst
            tst=0;
            
            %(5)
            for j=endj:length(irrev)
                if diff(j)>(2*irrstd)
                    n=n+1;
                    if n==40 && j>n+40
                        tst=timeev(j-80);
                        starti=j-80;
                        break
                    end
                    if n==40 && j<n+40
                        tst=timeev(j-39); 
                        starti=j-39; 
                        break
                    end
                else
                    n=0;
                end
            end
            
            %(6)
            if tst==0
                for j=endj:length(irrev)
                    if diff(j)>(2*irrstd)
                        n=n+1;
                        if n==30 && j>n+30
                            tst=timeev(j-60);
                            starti=j-70;
                            break
                        end
                        if n==30 && j<n+30
                            tst=timeev(j-29);
                            starti=j-29;
                            break
                        end
                    else
                        n=0;
                    end
                end
            end

            %(7)
            if tst==0
                for j=endj:length(irrev)
                    if diff(j)>(2*irrstd)
                        n=n+1;
                        if n==20 && j>n+20
                            tst=timeev(j-40)
                            starti=j-40;
                            break
                        end
                        if n==20 && j<n+20
                            tst=timeev(j-19); 
                            starti=j-19;
                            break
                        end
                    else
                        n=0;
                    end
                end
            end
            
            %(8)
            if tst==0
                for j=endj:length(irrev)
                    if diff(j)>(2*irrstd)
                        n=n+1;
                        if n==10 && j>n+10
                            tst=timeev(j-20)
                            starti=j-20;; 
                            break
                        end
                        if n==10 && j<n+10
                            tst=timeev(j-9); 
                            starti=j-9;
                            break
                        end
                    else
                        n=0;
                    end
                end
            end
            if tst==0
                starti=100000;
            end
        end
    end
    
    %at this point, j will be the index at which the criteria (whether
    %in 1, 2, or 3) was satisfied.  Carry this along: 
    startj = j;
    
    %(5)-(7) account for if the found start time was already counted as an 
    %event. If it has been, it will have been stored in "starttimes" below already.
    %Then start from startj (the point at which the criteria was met for
    %the first start time) and move forwards.  For (4) through (6), simply
    %repeat the steps above, but start after the end time of the last
    %flare (calculated below as endj).  (5) and (6) lower
    %"number-of-points" criteria in the same way as (2) and (3)

    
    
    %if the time of start is actually below the first element in
    %timeev (the window of time being studied), the start time will lie
    %outside.  Skip these events for now, still plotting solar quiet and 
    %irradiance, but not start or end times - but will need to be fixed later.
    if starti > length(timeev) || endj > length(timeev) 
        'No classification of flare!'
        
    
    
    elseif tst<timeev(1) 
        
        'Start too early!'
        
        
    else
        tst=timeev(starti);
        %Now the fun bit.  Find irradiance peak for well-behaved data,
        %first determining the next time that a light curve returns below
        %the solar quiet for a certain number of consecutive points (50, in
        %this case)

        %N.B. this first piece of logic only works if the start time 
        %was found using the first condition! (and doesn't work great if 
        %the increase in light curve is so fast that you're quickly on the decay phase)
        m=0;
        for j=(starti+40):length(diff)
            if diff(j)<0.5*irrstd
                m=m+1;
                if m==50
                    endj=j;
                    i
                    
                    break
                end
            end
        end
        if isnan(endj)
            
            disp('endj is NaN')
        else
        
        %find the tend, the end time for the event, by identifying the
        %position within timeev (the window of time for the event being
        %studied) at which endj lies
        tend=timeev(endj);
        
        %now the window of the flare is consolidated into "window,"
        %consisting of the irradiance values between startj and endj
        window=irrev(starti:endj);
        windowt=timeev(starti:endj);
        sqwindow = sqev(starti:endj);
        sq=sqev(1:starti);
                
        %find the average irradiance values within a 15-point spread around
        %each point (greater than 7 positions from the start and less than
        %7 positions from the end, of course) - this will give a better
        %idea of the actual average at the peak
        avgmeans=NaN(length(window)-7,1);
        for k=8:length(window)-7
            avgmeans(k-7)=nanmean(window(k-7:k+7));
        end
        
        %find the maximum in this average of the window
        maxind1=find(avgmeans==max(avgmeans));
        
        %find where the maximum lies in the original window, calculated at
        %the beginning of studying the flare
        maxind=starti+maxind1;
        
        %find the time of the maximum
        maxt=timeev(maxind);
        if length(maxt)>1
            maxt=maxt(1)
        end
        end
    end
    

    
    %You're done with finding critical times for that flare (well, assuming
    %it was nicely behaved enough for the code)! 
    
    %Store the times in arrays initialized above
    if ~isnan(endj) && starti<length(timeev) && endj <length(timeev) && starti<endj
    starttimes(i)=tst;
    endtimes(i)=tend;
    maxtimes(i)=maxt;
    

%%%%%END ORIGINAL MODEL-BUILDING%%%%%

    %calculate the impulsivity value! Now the window is between the start
    %and end - repeat finding the avgmeans/average values with a 15 point
    %window - this can be changed to a smaller or larger window as you
    %wish.
    window=irrev(starti:endj);
    windowt=timeev(starti:endj);
    sqwind=sqev(starti:endj);
    errwind=errev(starti:endj);
    sqpre=irrev(1:starti);
    sqerr = nanstd(sqpre);
    
    windowthr=windowt(1):((windowt(2)-windowt(1))/8):windowt(length(windowt));
    
    avgmeans=NaN(length(window)-7,1);
    sqav = mean(sqev);
    for k=8:length(window)-7
        avgmeans(k-7)=nanmean(window(k-7:k+7));
    end
    models = NaN(length(windowt),1);

    %identify where the maximum values actually lie...
    maxI_304 = max(window);
    
    %...and find the index corresponding to this in the window...
    %(1) option 1 is to use the raw maximum of the dataset
    maxind1=find(window==maxI_304);
    maxind=starti+maxind1;
    %(2) option 2 is to use the average maximum over 15 point window
%     maxind1=find(avgmeans==max(avgmeans));
%     maxind=startj+8+maxind1;
    %...but also identify where that is in the original sqev, irrev, and
    %timeev arrays.
    
    sq_at_max = sqev(maxind);
    peakt304 = timeev(maxind);
    
    %now the irradiance at half height (relative to solar quiet) is just
    %the average of the solar quiet value at the max and the max irradiance
    %itself - this needs some work, consider other options.
    irr_hh = (sqav+maxI_304)/2;
    
    %if this half height actually has an irradiance value in it, find the
    %times corresponding to irradiance points which lie roughly around this
    %irradiance, and take the mean time in this window to find the low and
    %high values of time at half height
    if isempty(irr_hh)==0
        %base impulsivity calculation only on the spline - see original
        %script for other options
        for k=2:length(window)
            if isnan(window(k))
                window(k)=window(k-1);
            end
        end
        for k=length(window)-1:-1:1
            if isnan(window(k))
                window(k)=window(k+1)
            end
        end
        spl = interp1(windowt,window,windowthr);
        splsq = interp1(windowt,sqwind,windowthr);
        splerr = interp1(windowt,errwind,windowthr);
        
        
        smspl = smooth(spl,5);
        vsmspl = smooth(spl,20);        
        
        %define error of each data point as the residual between a very
        %smoothed curve with the only slightly smoothed curve
        sqspl = smooth(splsq,5);
        errorbars = abs(vsmspl-spl');
        errm = nanmean(errorbars);
       	for o=1:length(errorbars)
            if errorbars(o) == 0
                errorbars(o) = errm;
            end
        end
        resid=smspl-sqspl;
        lenres=1:length(resid);
         %identify where the maximum values actually lie...
         
            maxI_304 = max(smspl);
            maxI_3042 = max(window);

            %...and find the index corresponding to this in the window...
            %use the raw maximum of the dataset - see other options in
            %original script
            maxind1=find(smspl==maxI_304);
            maxind2=find(window==maxI_3042);
            maxind=starti+maxind2;

        if maxind1 >1
           
            %...but also identify where that is in the original sqev, irrev, and
            %timeev arrays.
            
            sq_at_max = sqev(maxind);
            peakt304 = timeev(maxind);

            %now the irradiance at half height (relative to solar quiet) is just
            %the average of the solar quiet value at the max and the max irradiance
            %itself - this needs some work, consider other options.
            irr_hh = (sqav+maxI_304)/2;

            [d,ix] = min(abs(smspl(1:maxind1)-irr_hh));
            t_hh_low = windowthr(ix);

            [d,ix] = min(abs(smspl(maxind1:length(smspl))-irr_hh));
            t_hh_high = windowthr(maxind1+ix);

            if length(t_hh_low) > 1
                t_hh_low = t_hh_low(1);
            end

            if length(t_hh_high) > 1
                t_hh_high = t_hh_high(1);
            end        
            %distance between these two will be the FWHH
            t_hh = t_hh_high-t_hh_low;

            curly_I = maxI_304/t_hh;

            %store the impulsivity in the appropriate array
            curly_Is(i) = curly_I;

            %redo modeling process with the high-res window to compare to
            %spline


            risehr_t = windowthr(1:maxind1)';
            decayhr_t = windowthr((maxind1+1):length(spl))';
            risehr = spl(1:maxind1)';
            decayhr = spl((maxind1+1):length(spl))';
            riseerr=splerr(1:maxind1)';
            decayerr=splerr((maxind1+1):length(splerr));
            risehr_mirror=NaN(2*length(risehr),1);
            risehr_mirrort=NaN(2*length(risehr_t),1);
            dt = risehr_t(2) - risehr_t(1);

          %fit is poorly conditioned if datetime values used - switched to
            %length
            %[rise_fit_obj_hr,gof_3_hr] = fit(risehr,risehr_t,'poly3');
            if length(risehr_t) > 3 && length(decayhr_t) > 3
                [rise_fit_obj_exp1_hr,gof_exp1_hr] = fit((1:length(risehr_t))',risehr,'exp1');
                [rise_fit_obj_exp2_hr,gof_exp2_hr] = fit((1:length(risehr_t))',risehr,'exp2');

                n_coeff(1) = length(coeffvalues(rise_fit_obj_exp1_hr));
                n_coeff(2) = length(coeffvalues(rise_fit_obj_exp2_hr));


                rise_fit_exp2_hr = feval(rise_fit_obj_exp2_hr,(1:length(risehr_t))');
                %rise_fit3_hr = feval(rise_fit_obj_hr,risehr);
                rise_fit_exp1_hr = feval(rise_fit_obj_exp1_hr,(1:length(risehr_t))');

                [decay_fit_obj_pwr1_hr,gofdecpwr1_hr] = fit(decayhr_t,decayhr,'power1');
                [decay_fit_obj_pwr2_hr,gofdecpwr2_hr] = fit(decayhr_t,decayhr,'power2');
                [decay_fit_obj_four4_hr,gofdec_four4_hr] = fit(decayhr_t,decayhr,'Fourier4');
                [decay_fit_obj_four3_hr,gofdec_four3_hr] = fit(decayhr_t,decayhr,'Fourier3');
                [decay_fit_obj_four2_hr,gofdec_four2_hr] = fit(decayhr_t,decayhr,'Fourier2');
                [decay_fit_obj_four1_hr,gofdec_four1_hr] = fit(decayhr_t,decayhr,'Fourier1');

                n_coeff(3) = length(coeffvalues(decay_fit_obj_pwr1_hr));
                n_coeff(4) = length(coeffvalues(decay_fit_obj_pwr2_hr));
                n_coeff(5) = length(coeffvalues(decay_fit_obj_four1_hr));
                n_coeff(6) = length(coeffvalues(decay_fit_obj_four2_hr));
                n_coeff(7) = length(coeffvalues(decay_fit_obj_four3_hr));
                n_coeff(8) = length(coeffvalues(decay_fit_obj_four4_hr));

                decay_fitpwr1= feval(decay_fit_obj_pwr1_hr,decayhr_t);
                decay_fitpwr2= feval(decay_fit_obj_pwr2_hr,decayhr_t);
                decay_fitfour1= feval(decay_fit_obj_four1_hr,decayhr_t);
                decay_fitfour2= feval(decay_fit_obj_four2_hr,decayhr_t);
                decay_fitfour3= feval(decay_fit_obj_four3_hr,decayhr_t);
                decay_fitfour4= feval(decay_fit_obj_four4_hr,decayhr_t);

                fits{i,1}=rise_fit_obj_exp1_hr;
                fits{i,2}=rise_fit_obj_exp2_hr;
                fits{i,3}=decay_fit_obj_four1_hr;
                fits{i,4}=decay_fit_obj_four2_hr;
                fits{i,5}=decay_fit_obj_four3_hr;
                fits{i,6}=decay_fit_obj_four3_hr;

                spl_riset = windowthr(1:maxind1)';
                spl_rise = spl(1:maxind1)';
                spl_decay = spl((maxind1)+1:length(spl))';
                err_rise = errorbars(1:maxind1);
                err_decay = errorbars(maxind1+1:length(spl));
                riser2arr_exp2 = corrcoef(spl_rise,rise_fit_exp2_hr);
                riser2arr_exp1 = corrcoef(spl_rise,rise_fit_exp1_hr);
                %rise_3_r2arr = corrcoef(spl_riset,rise_fit3_hr);
                decay_four2_r2arr = corrcoef(spl_decay,decay_fitfour2);

                rise_exp2_r2 = riser2arr_exp2(1,2);
                rise_exp1_r2 = riser2arr_exp1(1,2);
                %rise_3_r2 = rise_3_r2arr(1,2);
                decay_four2_r2 = decay_four2_r2arr(1,2);

                adjrsq=[gof_exp1_hr.adjrsquare,gof_exp2_hr.adjrsquare,gofdecpwr1_hr.adjrsquare,gofdecpwr2_hr.adjrsquare,gofdec_four1_hr.adjrsquare,...
                    gofdec_four2_hr.adjrsquare,gofdec_four3_hr.adjrsquare,gofdec_four4_hr.adjrsquare];

                %possible error performance
                errstr={'mae','mse','rmse','mare','msre','rmsre','mape','mspe','rmspe'};
                errstrwr2={'mae','mse','rmse','mare','msre','rmsre','mape','mspe','rmspe','r2','adjr2'};
                decaymods=[decay_fitpwr1,decay_fitpwr2,...
                    decay_fitfour1,decay_fitfour2,decay_fitfour3,decay_fitfour4];
                %%error performance metrics
                errtab = NaN(8,length(errstr)+2);
                chisqs = NaN(8,4);

                %vetting by comparison of models - flag of 0 if we have
                %issues
                if ~isempty(fits{i,4}) && length(decay_fitfour2) > 50
                    if decay_fitfour2(52)>decay_fitfour2(51) && decay_fitfour2(53)>decay_fitfour2(52) && decay_fitfour2(54)>decay_fitfour2(53) && decay_fitfour2(55)>decay_fitfour2(54) && decay_fitfour2(2)>decay_fitfour2(1) && decay_fitfour2(3)>decay_fitfour2(2) && decay_fitfour2(4)>decay_fitfour2(3) && decay_fitfour2(5)>decay_fitfour2(4)
                        flags{i} = ribbondb_info{i,2};
                    end
                end

                if ~isempty(fits{i,4})
                    if decay_fitfour2(1)-min(rise_fit_exp2_hr) < 0.6*(max(rise_fit_exp2_hr)-min(rise_fit_exp2_hr))
                        flags{i} = ribbondb_info{i,2};
                    end
                end

                if ~isempty(fits{i,4})
                    if 0.3*(max(decay_fitfour2)-min(decay_fitfour2)) > (max(rise_fit_exp2_hr)-min(rise_fit_exp2_hr))
                        flags{i} = ribbondb_info{i,2};
                    end
                end

                for u=1:8 %all decay models and all rise models
                    for m=1:length(errstr)
                        if u == 1
                            met = errperf(spl_rise,rise_fit_exp1_hr,errstr{m});
                        elseif u == 2
                            met = errperf(spl_rise,rise_fit_exp2_hr,errstr{m});
                        elseif u > 2
                            met = errperf(spl_decay,decaymods(:,u-2),errstr{m});
                        end
                        errtab(u,m) = met;
                    end
                    m=10;

                    if u == 1
                        cc = corrcoef(spl_rise,rise_fit_exp1_hr');   
                    elseif u == 2
                        cc = corrcoef(spl_rise,rise_fit_exp2_hr');
                    elseif u > 2
                        cc = corrcoef(spl_decay,decaymods(:,u-2));
                    end
                    errtab(u,m) = cc(1,2);                    

                    if u == 1
                        resid = spl_rise-rise_fit_exp1_hr;
                        stand_dev = nanstd(resid);
                        chisq = 0;
                        %option 3 using errorbars from sdo/eve file - see other
                        %options in original file
                        for t=1:length(resid)
                            chi=(resid(t)/riseerr(t))^2;
                            chisq=chisq+chi;
                        end
                        df = length(spl_riset)-n_coeff(u);
                        redchisq1 = chisq/df;
                        %option 4 using errorbars from sdo/eve file but using the
                        %equation from formula
                        redchisq2=chi_squared(spl_rise,rise_fit_exp1_hr,n_coeff(u),riseerr);
                    elseif u == 2
                        resid = spl_rise - rise_fit_exp2_hr;
                        stand_dev = nanstd(resid);
                        chisq=0;
                        %option 3 using errorbars from file
                        for t=1:length(resid)
                            chi=(resid(t)/riseerr(t))^2;
                            chisq=chisq+chi;
                        end
                         %option 4 using errorbars from sdo/eve file but using the
                        %equation from formula
                        redchisq2=chi_squared(spl_rise,rise_fit_exp2_hr,n_coeff(u),riseerr);
                        df = length(spl_rise)-n_coeff(u);
                        redchisq1 = chisq/df;
                    elseif u > 2

                        resid = spl_decay - decaymods(:,u-2);
                            stand_dev = nanstd(resid);
                        for h=1:length(spl_decay)
                            sqerrstd3(h) = sqerr;
                        end
                        chisq=0;
                        for t=1:length(resid)
                            chi=(resid(t)/decayerr(t))^2;
                            chisq=chisq+chi;
                        end
                        df = length(spl_decay)-n_coeff(u-2);
                        redchisq1 = chisq/df;
                         %option 4 using errorbars from sdo/eve file but using the
                        %equation from formula
                        redchisq2=chi_squared(spl_decay,decaymods(:,u-2),n_coeff(u),decayerr');
                    end
                    chisqs(u,1) = chisq;
                    chisqs(u,2) = redchisq1;
                    chisqs(u,3) = redchisq2;
                end
                errtab(:,11) = adjrsq;
                errtable = array2table(errtab,'VariableNames',errstrwr2,'RowNames',{'RiseExp1','RiseExp2','DecayPwr1','DecayPwr2','DecayFour1','DecayFour2','DecayFour3','DecayFour4'});

                errs{i}=errtable;
                chisquareds{i+1} = chisqs;

                modresid_rise_exp2 = spl_rise - rise_fit_exp2_hr;
                modresid_rise_exp1 = spl_rise - rise_fit_exp1_hr;
                modresid_decay = spl_decay - decay_fitfour2;
                mkdir('/Users/owner/Desktop/modeleval');

                rsqurs_exp2_four2(i,1) = rise_exp2_r2;
                rsqurs_exp2_four2(i,2) = decay_four2_r2;
                rsqurs_3_four2(i,1) = rise_exp2_r2;
                rsqurs_3_four2(i,2) = decay_four2_r2;
            else
                flags{i}='1' %if the fit was too short on either side
            end
        end
            
    end
    
    flid = char(ribbondb_info{i,2});
    fltitle = flid(1:13);
    

    %plotting - indices leave something to be desired
    if maxind(1) > 1 && strcmp(flags{i},'1') == 0
            
            figin = ceil(i/25);
            figfl = floor(i/26);
            figure(1)
            a=subplot(5,5,i-25*(figin-1))
            scatter(windowthr,smspl'/(1e-3),'r','.');
            hold on
            scatter(risehr_t,rise_fit_exp2_hr/(1e-3),'.');
            scatter(decayhr_t,decay_fitfour2/(1e-3),'.');
            title([fltitle ' [' num2str(curly_Is(i)) ']'],'Interpreter','None');
          
            flags{i} ='0'; %%% WHY THIS LINE?
    

    elseif maxind(1) ~= 1 && strcmp(flags{i},'1') == 0
            
            figin = ceil(i/25);
            figfl = floor(i/26);
            figure(1)
            subplot(5,5,i-25*(figin-1))
            cla(f)
            scatter(timeev,irrev/(1e-3),'r','.');

            title(fltitle,'Interpreter', 'none');
     
    
    elseif strcmp(flags{i},'1') == 1
        
        figin = ceil(i/25);
        figfl = floor(i/26);
        f=figure(1)
        subplot(5,5,i-25*(figin-1))
        cla(f)
        scatter(timeev,zeros(length(timeev),1),'r','.');
  
         title(fltitle,'Interpreter', 'none');
        
   
    end
    
    else
        
        figin = ceil(i/25);
        figfl = floor(i/26);
        f=figure(1)
        subplot(5,5,i-25*(figin-1))
        cla(f)
        scatter(timeev,zeros(length(timeev),1),'r','.');
  
         title(fltitle,'Interpreter', 'none');
        
    
        
    end
    if i-25*(figin-1) == 1
        subplot(5,5,i-25*(figin-1));
        xlabel('Datetime');
        ylabel('Flux (mW/m^2)');
    end
    %if rem(i,length(i))==0
        sgtitle('Format: YYYYMMDD_(starttime) [impulsivity]','Interpreter','none');
    %end
    if rem(i,25) == 0 || i == length(windowstart2)
        saveas(f,['/Users/owner/Desktop/flaresummary_ver/flares',num2str(i-24),'to',num2str(i),'.png']);
        f=figure(1);
        clf
        set(gcf,'Position',[100 100 1500 1500])
    end
    
    
end

%%save(filename,'fits','-append');

%  save(filename,'chisquareds','-append');
%  save(filename,'errs','-append');
%  save(filename,'rsqurs_exp2_four2','-append');
%  save(filename,'flags','-append');
%  %save(filename,'rsqurs_3_four2','-append');
% %save(filename,'rise_chi','rise_pval','rise_chistats','-append');
%     
% %once all flares have been dealt with, consider saving all the parameters
% %you've calculated to an array and put in the same file which was loaded
% %before
% save(filename2,'starttimes','endtimes','maxtimes','curly_Is','windowstart2');

    