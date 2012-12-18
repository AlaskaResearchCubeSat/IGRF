function IGRFtst(com,baud)
    if(nargin<3)
        baud=9600;
    end
    if(nargin<2)
        com='COM23';
    end
    
    %load magnetic field data
    mc=magfd();
    
    %mag model for reduced degree model
    mcl=magfd();
    mcl.NMAX=3;
    
    year=2013.5;
    
    Re = 6.378137*10^6;                 % Radius of the Earth (m)[IUGG value of equatorial radius]
    h  = 600000;                        % Orbit Altitude (m)
    R  =  (Re + h)/1000;                

    
    lat=0:1:360;
    lon=-90:1:90;
    %preallocate angle
    err_ang=zeros(length(lon),length(lat));
    
    try
        %open serial port
        %ser=serial(com,'BaudRate',baud);
        %set timeout to 15s
        %set(ser,'Timeout',15);
        %open port
        %fopen(ser);
    
        for k=1:length(lat)
            for kk=1:length(lon)
                %field_pc=mc.body(year,R,90-lat(k),lon(kk));
                field_pc=mc.lned(year,R,90-lat(k),lon(kk));
%                 line='';
%                 while ~strncmp(line,'unknown',length('unknown'))   
%                     %format command
%                     cmd=sprintf('shval3 %0.2f %0.4f %f %f',year,R,lat(k)*pi/180,lon(kk)*pi/180)
%                     fprintf(ser,'%s\n',cmd);
%                     %capture echo
%                     fgetl(ser);
%                     %get result
%                     line=fgetl(ser);
%                     if strncmp(line,'unknown',length('unknown'))
%                         line
%                         continue;
%                     end
%                     %wait for command to complete
%                     if ~waitReady(ser,5)
%                         error('MSP430 not responding');
%                     end
%                     %extract values
%                     field_msp=sscanf(line,'%f %f %f');
%                     if(length(field_msp)~=3)
%                         line
%                         continue;
%                     end
                    field_msp=mcl.lned(year,R,90-lat(k),lon(kk));
                    err_ang(kk,k)=abs(180/pi*acos(dot(field_pc,field_msp)/(norm(field_pc)*norm(field_msp))));
%                 end
            end
        end
    catch err
        if exist('ser','var')
            if strcmp(ser.Status,'open')
                fclose(ser);
            end
            delete(ser);
        end
        rethrow(err);
    end
    %close serial port
    %fclose(ser);
    %delete serial object
    %delete(ser);
    %plot error and stuff
    pcolor(lat,lon,err_ang);
    colormap(hsv(512));
    shading('interp');
    %show colorbar for plot
    colorbar('eastOutside')
    hold on;
    %load topomap data
    load('topo.mat','topo','topomap1');
    
    %plot continent outlines
    contour(0:359,-89:90,topo,[0 0],'k');
    axis('square');
    axis('equal');
    axis([0 360 -90 90]);
    hold off;
    fprintf('Max Error = %fdeg\nMin Error = %fdeg\n',max(max(err_ang)),min(min(err_ang)));
end

function [success]=waitReady(sobj,timeout,output)
    if nargin<3
        output=false;
    end
    if nargin<2
        timeout=5;
    end
    msg=0;
    count=0;
    while msg(end)~='>'
        len=sobj.BytesAvailable;
        if len==0
            if count*3>=timeout
                success=false;
                return
            end
            pause(3);
            count=count+1;
            continue;
        end
        [msg,~,~]=fread(sobj,len);
        if output
            fprintf('%s\n',char(msg'));
        end
    end
    success=true;
end