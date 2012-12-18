classdef  magfd < handle
%  MAGFD [Last Updated 08 April 2011]
%  Function to compute Earths magnetic field
%  and components: X,Y,Z,T for a given latitude
%  and longitude, date and altitude.
%  Uses MATLAB MAT files sh1900.mat to sh2010.mat in 5 yr
%  intervals. ONLY uses harmonic expansion to order 10 at present.
%
%  Usage: out=magfd(DATE,ITYPE,ALT,COLAT,ELONG);
%
%  DATE = date of survey (decimal years)
%  ITYPE=1 for geodetic coordinates (usual case)
%  ITYPE=2 for geocentric coordinates
%  ALT = (for ITYPE=1) altitude of survey relative to sealevel (km +ve up)
%  ALT = (for ITYPE=2) radial distance from center of earth in km
%  COLAT=90-latitude (decimal degrees)
%  ELONG=longitude of survey (decimal degrees)
%
%  Output array out contains components X,Y,Z,T in nanoteslas
%   X north component
%   Y east component     
%   Z vertical component +ve down
%   T total field magnitude
%
%  ref: IAGA, Division V, Working Group VMOD, 
%   The 10th generation International Geomagnetic 
%   Reference Field, Geophys. J. Int, 161, 561-565, 2005.
%
% Maurice A. Tivey March 1997
% Mod Dec 1999 (add igrf2000 and y2k compliance
% Mod Nov 2000 (use up to degree 10 sh coefficients)
% Mod Apr 2005 added 2005 coeffs
% Mod Sep 2006 some clean up and info added
% Mod Jan 2010 added 2010 coefficients
% http://deeptow.whoi.edu/matlab.html
% Copyright: Maurice A. Tivey, 2005
%  Woods Hole Oceanographic Institution
    
    properties
        datYear
        agh
        dgh
        last
        NMAX=13; % Max number of harmonic degrees
    end
    %J=magfd(DATE,ITYPE,ALT,COLAT,ELONG,a)
    properties (Access = private)
        DGRF=1000:5:2010;
        igrfyear=2010;
        igrffile='sh2010';
    end
    methods
        function mc=magfd(date)
            if nargin<1
                date=mc.igrfyear;
            end
            mc.readCoeff(date)
        end
        function readCoeff(mc,DATE)
            DATE=abs(DATE);
            % Determine year for base DGRF to use.
             if DATE < mc.igrfyear,
              error('No data for past years');
             else
              S=load(mc.igrffile);   % load in igrf data file
              mc.datYear=mc.igrfyear;
             end
            % combine spherical harmonic coefficients from first 8 degrees 
            % with degrees 9 thru 13 
             mc.agh=[S.agh,S.agh41];
             mc.dgh=[S.dgh,S.dgh41];
        end
        function checkDate(mc,d)
            %TODO: at some point this should load the appropreate data file
            if d < mc.igrfyear,
              error('No data for past years');
            end
        end
        function Fixed=lned(mc,DATE,ALT,COLAT,ELONG)
            
            T     = DATE - mc.datYear;
           
            mc.checkDate(DATE);
            
            
            
            D2R   = pi/180;
            R     = ALT;
            SLAT  = cos(COLAT*D2R);
            CLAT  = sin(COLAT*D2R);
            X     = 0.0;
            Y     = 0.0;
            Z     = 0.0;
            L     = 1;
            M     = 1;
            N     = 0;
            RE    = 6371.2; % Earth's mean radius
            % if geocentric coordinates desired then only need to define the following
            RATIO = RE/R;
           
            %number of loop iterations
            NPQ=(mc.NMAX*(mc.NMAX+3))/2;
            
            %preallocate for speed
            P=zeros(1,NPQ);
            Q=zeros(1,NPQ);
            SL=zeros(1,mc.NMAX);
            CL=zeros(1,mc.NMAX);
            
            
            CL(1) = cos(ELONG*D2R);
            SL(1) = sin(ELONG*D2R);
            %
            %     COMPUTATION OF SCHMIDT QUASI-NORMAL COEFFICIENTS  P AND X(=Q)
            %
              P(1)  = 2.0*SLAT;
              P(2)  = 2.0*CLAT;
              P(3)  = 4.5*SLAT*SLAT - 1.5;
              P(4)  = sqrt(27)*CLAT*SLAT;
              Q(1)  = -CLAT;
              Q(2)  =  SLAT;
              Q(3)  = -3.0*CLAT*SLAT;
              Q(4)  = sqrt(3)*(SLAT*SLAT - CLAT*CLAT);

      
            for K=1:NPQ,
                if N < M 
                    M     = 0;
                    N     = N + 1;
                    RR    = RATIO^(N + 2);
                    FN    = N;
                end
                FM    = M;
                if K >= 5 %8,5,5
                    if (M-N) == 0 %,7,6,7
                        ONE   = sqrt(1.0 - 0.5/FM);
                        J     = K - N - 1;
                        P(K)  = (1.0 + 1.0/FM)*ONE*CLAT*P(J);
                        Q(K)  = ONE*(CLAT*Q(J) + SLAT/FM*P(J));
                        SL(M) = SL(M-1)*CL(1) + CL(M-1)*SL(1);
                        CL(M) = CL(M-1)*CL(1) - SL(M-1)*SL(1);
                    else
                        ONE   = sqrt(FN*FN - FM*FM);
                        TWO   = sqrt((FN - 1.0)^2 - FM*FM)/ONE;
                        THREE = (2.0*FN - 1.0)/ONE;
                        I     = K - N;
                        J     = K - 2*N + 1;
                        P(K)  = (FN + 1.0)*(THREE*SLAT/FN*P(I) - TWO/(FN - 1.0)*P(J));
                        Q(K)  = THREE*(SLAT*Q(I) - CLAT/FN*P(I)) - TWO*Q(J);
                    end
                    %
                    %     SYNTHESIS OF X, Y AND Z IN GEOCENTRIC COORDINATES
                    %
                end
                ONE   = (mc.agh(L) + mc.dgh(L)*T)*RR;

                if M == 0 %10,9,10
                    X     = X + ONE*Q(K);
                    Z     = Z - ONE*P(K);
                    L     = L + 1;
                else
                    TWO   = (mc.agh(L+1) + mc.dgh(L+1)*T)*RR;
                    THREE = ONE*CL(M) + TWO*SL(M);
                    X     = X + THREE*Q(K);
                    Z     = Z - THREE*P(K);
                    if CLAT > 0 %12,12,11
                        Y = Y+(ONE*SL(M)-TWO*CL(M))*FM*P(K)/((FN + 1.0)*CLAT);
                    else
                        Y = Y + (ONE*SL(M) - TWO*CL(M))*Q(K)*SLAT;
                    end
                    L     = L + 2;
                end
                M     = M + 1;
            end
            
            X     = X * 10^-9;
            Y     = Y * 10^-9;
            Z     = Z * 10^-9;
            Fixed = [X,Y,Z]';
            mc.last=Fixed;
        end
        function B=body(mc,DATE,ALT,COLAT,ELONG,a)
            %default to identity matrix transform
            if nargin<6
                a=eye(3);
            end
            % Calulate latitude and longitude in radians
            latr  = (90-COLAT)*pi/180;
            longr = ELONG*pi/180;

            % Calculate rotation matrix for LNED Equation (2) from Thesis
            T(1,1) = -sin(latr)*cos(longr); T(1,2) = -sin(longr); T(1,3) = -cos(latr)*cos(longr);
            T(2,1) = -sin(latr)*sin(longr); T(2,2) =  cos(latr);  T(2,3) = -cos(latr)*sin(longr);
            T(3,1) =  cos(latr);            T(3,2) = 0;           T(3,3) = -sin(latr);

            
            B=a*T*mc.lned(DATE,ALT,COLAT,ELONG);
        end
    end
end


