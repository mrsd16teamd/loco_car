% /*
%  * Copyright (c) 2008, Maxim Likhachev
%  * All rights reserved.
%  * 
%  * Redistribution and use in source and binary forms, with or without
%  * modification, are permitted provided that the following conditions are met:
%  * 
%  *     * Redistributions of source code must retain the above copyright
%  *       notice, this list of conditions and the following disclaimer.
%  *     * Redistributions in binary form must reproduce the above copyright
%  *       notice, this list of conditions and the following disclaimer in the
%  *       documentation and/or other materials provided with the distribution.
%  *     * Neither the name of the University of Pennsylvania nor the names of its
%  *       contributors may be used to endorse or promote products derived from
%  *       this software without specific prior written permission.
%  * 
%  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%  * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  * POSSIBILITY OF SUCH DAMAGE.
%  */

%{
function y = radians(deg)
   y = deg * pi / 180;
end;
%}

function[] = genmprim_unicycleplussidewaysplusbackturn(outfilename)

%
%generates motion primitives and saves them into file
%
%written by Maxim Likhachev
%---------------------------------------------------
%

%defines

UNICYCLE_MPRIM_32DEGS = 1;

ANGLE_THRESHOLD_DEG = 2;

if UNICYCLE_MPRIM_32DEGS == 1
    resolution = 0.05;
    numberofangles = 32; %preferably a power of 2, definitely multiple of 8
    numberofprimsperangle = 7;

    %multipliers (multiplier is used as costmult*cost)
    forwardcostmult = 1;
    backwardcostmult = 10;
    forwardandturncostmult = 2;
    sidestepcostmult = 20;
    turninplacecostmult = 20;
    
    %note, what is shown x,y,theta *changes* (that is, dx,dy,dtheta and not absolute numbers)
    
    %0 degreees
    basemprimendpts0_c = zeros(numberofprimsperangle, 4); %x,y,theta,costmult 
    %angles are positive counterclockwise
    %0 theta change
    basemprimendpts0_c(1,:) = [1 0 0 forwardcostmult];
    basemprimendpts0_c(2,:) = [8 0 0 forwardcostmult];
    basemprimendpts0_c(3,:) = [-1 0 0 backwardcostmult];    
    %1/32 theta change
    basemprimendpts0_c(4,:) = [8 1 1 forwardandturncostmult];
    basemprimendpts0_c(5,:) = [8 -1 -1 forwardandturncostmult];
    %1/32 theta change going backward
    basemprimendpts0_c(6,:) = [-8 -1 1 backwardcostmult];
    basemprimendpts0_c(7,:) = [-8 1 -1 backwardcostmult];
 %{ 
    %11.25 degrees
    basemprimendpts11p25_c = zeros(numberofprimsperangle, 4); %x,y,theta,costmult (multiplier is used as costmult*cost)
    %angles are positive counterclockwise
    %0 theta change     
    basemprimendpts11p25_c(1,:) = [4 1 0 forwardcostmult];
    basemprimendpts11p25_c(2,:) = [8 2 0 forwardcostmult];    
    basemprimendpts11p25_c(3,:) = [-4 -1 0 backwardcostmult];     
    %1/32 theta change
    basemprimendpts11p25_c(4,:) = [8 3 1 forwardandturncostmult];
    basemprimendpts11p25_c(5,:) = [9 1 -1 forwardandturncostmult];    
    %1/32 theta change going back
    basemprimendpts11p25_c(6,:) = [-8 -3 1 backwardcostmult];
    basemprimendpts11p25_c(7,:) = [-9 -1 -1 backwardcostmult];    

    %22.5 degrees
    basemprimendpts22p5_c = zeros(numberofprimsperangle, 4); %x,y,theta,costmult (multiplier is used as costmult*cost)
    %angles are positive counterclockwise
    %0 theta change 
    basemprimendpts22p5_c(1,:) = [2 1 0 forwardcostmult];
    basemprimendpts22p5_c(2,:) = [6 3 0 forwardcostmult];
    basemprimendpts22p5_c(3,:) = [-2 -1 0 backwardcostmult];    
    %1/32 theta change
    basemprimendpts22p5_c(4,:) = [5 4 1 forwardandturncostmult];
    basemprimendpts22p5_c(5,:) = [7 2 -1 forwardandturncostmult];    
    %1/32 theta change going back
    basemprimendpts22p5_c(6,:) = [-5 -4 1 backwardcostmult];
    basemprimendpts22p5_c(7,:) = [-7 -2 -1 backwardcostmult];    

    %33.75 degrees
    basemprimendpts33p75_c = zeros(numberofprimsperangle, 4); %x,y,theta,costmult (multiplier is used as costmult*cost)
    %angles are positive counterclockwise
    %0 theta change 
    basemprimendpts33p75_c(1,:) = [3 2 0 forwardcostmult];
    basemprimendpts33p75_c(2,:) = [6 4 0 forwardcostmult];
    basemprimendpts33p75_c(3,:) = [-3 -2 0 backwardcostmult];    
    %1/32 theta change
    basemprimendpts33p75_c(4,:) = [5 5 1 forwardandturncostmult];
    basemprimendpts33p75_c(5,:) = [7 3 -1 forwardandturncostmult];    
    %1/32 theta change going back
    basemprimendpts33p75_c(6,:) = [-5 -5 1 backwardcostmult];
    basemprimendpts33p75_c(7,:) = [-7 -3 -1 backwardcostmult];    

    %45 degrees
    basemprimendpts45_c = zeros(numberofprimsperangle, 4); %x,y,theta,costmult (multiplier is used as costmult*cost)
    %angles are positive counterclockwise
    %0 theta change 
    basemprimendpts45_c(1,:) = [1 1 0 forwardcostmult];
    basemprimendpts45_c(2,:) = [6 6 0 forwardcostmult];
    basemprimendpts45_c(3,:) = [-1 -1 0 backwardcostmult];
    %1/32 theta change
    basemprimendpts45_c(4,:) = [5 7 1 forwardandturncostmult];
    basemprimendpts45_c(5,:) = [7 5 -1 forwardandturncostmult];
    %1/32 theta change going back
    basemprimendpts45_c(6,:) = [-5 -7 1 backwardcostmult];
    basemprimendpts45_c(7,:) = [-7 -5 -1 backwardcostmult];
%}
   else
    fprintf(1, 'ERROR: undefined mprims type\n');
    return;    
end;
    
    
fout = fopen(outfilename, 'w');


%write the header
fprintf(fout, 'resolution_m: %f\n', resolution);
fprintf(fout, 'numberofangles: %d\n', numberofangles);
fprintf(fout, 'totalnumberofprimitives: %d\n', numberofprimsperangle*numberofangles);

%iterate over angles
for angleind = 1:numberofangles
    
    figure(1);
    hold off;

    text(0, 0, int2str(angleind));
    
    %iterate over primitives    
    for primind = 1:numberofprimsperangle
        fprintf(fout, 'primID: %d\n', primind-1);
        fprintf(fout, 'startangle_c: %d\n', angleind-1);

        %current angle
        currentangle = (angleind-1)*2*pi/numberofangles;
        currentangle_36000int = round((angleind-1)*36000/numberofangles);
        
        fprintf(1, 'currentangle: %f currentangle_3600int: %f\n', currentangle, currentangle_36000int);
        basemprimendpts_c = basemprimendpts0_c(primind,:);    
        angle = currentangle;

%{
        %compute which template to use
        if (rem(currentangle_36000int, 9000) == 0)
            basemprimendpts_c = basemprimendpts0_c(primind,:);    
            angle = currentangle;
            fprintf(1, '90\n');
        elseif (rem(currentangle_36000int, 4500) == 0)
            basemprimendpts_c = basemprimendpts45_c(primind,:);
            angle = currentangle - 45*pi/180;
            fprintf(1, '45\n');
        elseif (rem(currentangle_36000int-7875, 9000) == 0)
            basemprimendpts_c = basemprimendpts33p75_c(primind,:);
            basemprimendpts_c(1) = basemprimendpts33p75_c(primind, 2); %reverse x and y
            basemprimendpts_c(2) = basemprimendpts33p75_c(primind, 1);
            basemprimendpts_c(3) = -basemprimendpts33p75_c(primind, 3); %reverse the angle as well
            angle = currentangle - 78.75*pi/180;
            fprintf(1, '78p75\n');
        elseif (rem(currentangle_36000int-6750, 9000) == 0)
            basemprimendpts_c = basemprimendpts22p5_c(primind,:);
            basemprimendpts_c(1) = basemprimendpts22p5_c(primind, 2); %reverse x and y
            basemprimendpts_c(2) = basemprimendpts22p5_c(primind, 1);
            basemprimendpts_c(3) = -basemprimendpts22p5_c(primind, 3); %reverse the angle as well
            %fprintf(1, '%d %d %d onto %d %d %d\n', basemprimendpts22p5_c(1), basemprimendpts22p5_c(2), basemprimendpts22p5_c(3), ...
            %    basemprimendpts_c(1), basemprimendpts_c(2), basemprimendpts_c(3));
            angle = currentangle - 67.5*pi/180;
            fprintf(1, '67p5\n');            
        elseif (rem(currentangle_36000int-5625, 9000) == 0)
            basemprimendpts_c = basemprimendpts11p25_c(primind,:);
            basemprimendpts_c(1) = basemprimendpts11p25_c(primind, 2); %reverse x and y
            basemprimendpts_c(2) = basemprimendpts11p25_c(primind, 1);
            basemprimendpts_c(3) = -basemprimendpts11p25_c(primind, 3); %reverse the angle as well
            fprintf(1, 'basemprimendpts11p25_c(primind, 2): %d\n', basemprimendpts11p25_c(primind, 2));
            angle = currentangle - 56.25*pi/180;
            fprintf(1, '56p25\n');
        elseif (rem(currentangle_36000int-3375, 9000) == 0)
            basemprimendpts_c = basemprimendpts33p75_c(primind,:);
            angle = currentangle - 33.75*pi/180;
            fprintf(1, '33p75\n');
        elseif (rem(currentangle_36000int-2250, 9000) == 0)
            basemprimendpts_c = basemprimendpts22p5_c(primind,:);
            angle = currentangle - 22.5*pi/180;
            fprintf(1, '22p5\n');
        elseif (rem(currentangle_36000int-1125, 9000) == 0)
            basemprimendpts_c = basemprimendpts11p25_c(primind,:);
            angle = currentangle - 11.25*pi/180;
            fprintf(1, '11p25\n');
        else
            fprintf(1, 'ERROR: invalid angular resolution. angle = %d\n', currentangle_36000int);
            return;
        end;
%}
        
        %now figure out what action will be        
        baseendpose_c = basemprimendpts_c(1:3);
        additionalactioncostmult = basemprimendpts_c(4);        
        endx_c = baseendpose_c(1)*cos(angle) - baseendpose_c(2)*sin(angle);        
        endy_c = baseendpose_c(1)*sin(angle) + baseendpose_c(2)*cos(angle);
        endtheta_c = rem(angleind - 1 + baseendpose_c(3), numberofangles);
        fprintf(1, 'round(sqrt(%f * %f + %f * %f)) = %f\n', endx_c, endx_c, endy_c, endy_c, round(sqrt(endx_c * endx_c + endy_c * endy_c)));

        if round(sqrt(endx_c * endx_c + endy_c * endy_c)) == 1
%{
           if endy_c < endx_c
              multiplier = round(endx_c / endy_c);
           else
              multiplier = round(endy_c / endx_c);
           end;
           endy_c = round(endy_c * multiplier);
           endx_c = round(endx_c * multiplier);
           
           fprintf(1, 'is one!\n');
%}
           if endx_c == 0
              diff_angle = 90 * pi / 180;
              if endy_c < 0
                 diff_angle = -90 * pi / 180;
              end;
           else
              fprintf(1,'calc diff angle\n');
              diff_angle = atan(endy_c / endx_c);
           end;
          
           % continue until we found a valid angle
           while abs(diff_angle - angle) > (ANGLE_THRESHOLD_DEG * pi / 180)
             endx_c = cos(diff_angle);
             endy_c = sin(diff_angle);
              
           end;

        end;

        fprintf(1, 'unrounded endx_c: %f, unrounded endy_c: %f\n', endx_c, endy_c);
        endx_c = round(endx_c);
        endy_c = round(endy_c);
        endpose_c = [endx_c endy_c endtheta_c];
        
        fprintf(1, 'rotation angle=%f\n', angle*180/pi);
        fprintf(1, 'endx_c: %d, endy_c: %d, endtheta_c: %d\n', endx_c, endy_c, endtheta_c);
        
        if baseendpose_c(2) == 0 & baseendpose_c(3) == 0
            %fprintf(1, 'endpose=%d %d %d\n', endpose_c(1), endpose_c(2), endpose_c(3));
        end;
        
        %generate intermediate poses (remember they are w.r.t 0,0 (and not
        %centers of the cells)
        numofsamples = 10;
        intermcells_m = zeros(numofsamples,3);
        if UNICYCLE_MPRIM_32DEGS == 1
            startpt = [0 0 currentangle];
            endpt = [endpose_c(1)*resolution endpose_c(2)*resolution ...
                rem(angleind - 1 + baseendpose_c(3), numberofangles)*2*pi/numberofangles];
            intermcells_m = zeros(numofsamples,3);
            if ((endx_c == 0 & endy_c == 0) | baseendpose_c(3) == 0) %turn in place or move forward            
                for iind = 1:numofsamples
                    intermcells_m(iind,:) = [startpt(1) + (endpt(1) - startpt(1))*(iind-1)/(numofsamples-1) ...
                                            startpt(2) + (endpt(2) - startpt(2))*(iind-1)/(numofsamples-1) ...
                                            0];
                    rotation_angle = (baseendpose_c(3) ) * (2*pi/numberofangles);
                    intermcells_m(iind,3) = rem(startpt(3) + (rotation_angle)*(iind-1)/(numofsamples-1), 2*pi);
                end;            
            else %unicycle-based move forward or backward
                R = [cos(startpt(3)) sin(endpt(3)) - sin(startpt(3));
                    sin(startpt(3)) -(cos(endpt(3)) - cos(startpt(3)))];
                S = pinv(R)*[endpt(1) - startpt(1); endpt(2) - startpt(2)];
                l = S(1); 
                tvoverrv = S(2);
                rv = (baseendpose_c(3)*2*pi/numberofangles + l/tvoverrv);
                tv = tvoverrv*rv;
                         
                if ((l < 0 & tv > 0) | (l > 0 & tv < 0))
                    fprintf(1, 'WARNING: l = %d < 0 -> bad action start/end points\n', l);
                    l = 0;
                end;
                %compute rv
                %rv = baseendpose_c(3)*2*pi/numberofangles;
                %compute tv
                %tvx = (endpt(1) - startpt(1))*rv/(sin(endpt(3)) - sin(startpt(3)))
                %tvy = -(endpt(2) - startpt(2))*rv/(cos(endpt(3)) - cos(startpt(3)))
                %tv = (tvx + tvy)/2.0;              
                %generate samples
                for iind = 1:numofsamples                                        
                    dt = (iind-1)/(numofsamples-1);
                                        
                    %dtheta = rv*dt + startpt(3);
                    %intermcells_m(iind,:) = [startpt(1) + tv/rv*(sin(dtheta) - sin(startpt(3))) ...
                    %                        startpt(2) - tv/rv*(cos(dtheta) - cos(startpt(3))) ...
                    %                        dtheta];
                    
                    if(abs(dt*tv) < abs(l))
                        intermcells_m(iind,:) = [startpt(1) + dt*tv*cos(startpt(3)) ...
                                                 startpt(2) + dt*tv*sin(startpt(3)) ...
                                                 startpt(3)];
                    else
                        dtheta = rv*(dt - l/tv) + startpt(3);
                        intermcells_m(iind,:) = [startpt(1) + l*cos(startpt(3)) + tvoverrv*(sin(dtheta) - sin(startpt(3))) ...
                                                 startpt(2) + l*sin(startpt(3)) - tvoverrv*(cos(dtheta) - cos(startpt(3))) ...
                                                 dtheta];
                    end;
                end; 
                %correct
                errorxy = [endpt(1) - intermcells_m(numofsamples,1) ... 
                           endpt(2) - intermcells_m(numofsamples,2)];
                fprintf(1, 'l=%f errx=%f erry=%f\n', l, errorxy(1), errorxy(2));
                interpfactor = [0:1/(numofsamples-1):1];
                intermcells_m(:,1) = intermcells_m(:,1) + errorxy(1)*interpfactor';
                intermcells_m(:,2) = intermcells_m(:,2) + errorxy(2)*interpfactor';
            end;                                        
        end;
    
        %write out
        fprintf(fout, 'endpose_c: %d %d %d\n', endpose_c(1), endpose_c(2), endpose_c(3));
        fprintf(fout, 'additionalactioncostmult: %d\n', additionalactioncostmult);
        fprintf(fout, 'intermediateposes: %d\n', size(intermcells_m,1));
        for interind = 1:size(intermcells_m, 1)
            fprintf(fout, '%.4f %.4f %.4f\n', intermcells_m(interind,1), intermcells_m(interind,2), intermcells_m(interind,3));
        end;
        
        plot(intermcells_m(:,1), intermcells_m(:,2));
        axis([-0.3 0.3 -0.3 0.3]);
        text(intermcells_m(numofsamples,1), intermcells_m(numofsamples,2), int2str(endpose_c(3)));
        hold on;
        
    end;
    grid;
    fprintf(1,'\n');
    pause;
end;
        
fclose('all');
