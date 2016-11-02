% this script is designed to be run *after* setup.m
% and prepares the workspace for a line-by-line walk through of 
% sanim_XY_vehicle_viz.m
%
% Marc Compere, comperem@gmail.com
% created : 28 Dec 2015
% modified: 11 Jan 2016

anim_frame_name_str=char(datetime('now','Format','yyyy-MM-dd_HH_mm_ss'));

% for veh_object2.m
object_type=1;
obj_length=vehicle_length;
obj_width=vehicle_width;

plotAxisLimits = plotAxisLimits;

% for mdlInitializeSizes()
Config.axisMode=plotAxisLimitMode;
Config.ax=plotAxisLimits;
Config.Ts=animation_update_interval;
Config.enable_CG_trace=enable_CG_trace;
Config.enable_rearAxle_trace=enable_rearAxle_trace;
Config.save_anim_frames=save_anim_frames;
Config.anim_frame_name_str=anim_frame_name_str;
Config.L=vehicle_length;
Config.W=vehicle_width;
%N=1;


% for mdlUpdate()

% since N=1, just put in scalars
X   = 2; % (m) object positions in global XY frame
Y   = 6; % (m) object positions in global XY frame
yaw = pi/10; % (rad) object yaw orientations about global Z axis
delta = yaw/2;


% init:
%[sys,x0,str,ts] = sanim_XY_vehicle_viz(0,0,0,0,Config)


% mdlUpdate()
%[sys,x,str,ts] = sanim_XY_vehicle_viz(0,x0,[X_ic, Y_ic, yaw_ic, 0],2,Config)









