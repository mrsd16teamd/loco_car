function out=veh_object2(object_type,obj_length, obj_width)
%
% This file was developed from a copy of wheel_object.m.
%
% Generate vertices and faces for a top-view of a simple 2D vehicle patch
% object or tire.
%
% Assumed CG is at object's geometric center.
%
% Marc Compere
% comperem@gmail.com
% created : 22 August 2002
% modified: 28 Dec 2015

L = obj_length; % (m) L for vehicle or tire (depending on object_type)
W = obj_width; % (m) width of vehicle or tire


% --------------------------------------------------------------
% -------------------- vehicle object --------------------------
% --------------------------------------------------------------
% define 2D top-down view of vehicle body, starting from the nose, centered
% at vehicle's geometric center, g
% pt# 1       2     3           4      5     6        7       8     9  10
x_veh = (1/2)*[ 1,  0.95,   0.85,      -0.9,    -1,   -1,    -0.9,   0.85, 0.95,  1];
y_veh = (1/2)*[ 0,   0.6,      1,         1,   0.8, -0.8,      -1,     -1, -0.6,  0];

%x_veh = x_veh + 0.5; % (m) translate all vertices forward such that the
%                     %     origin (x,y)=(0,0) is at the rear axle center


% make the vehicel chassis a little longer than the wheelbase
veh_obj.vertices = 1.3*[L*x_veh;
                        W*y_veh];

% define the object faces
veh_obj.faces = 1:length(x_veh);



% --------------------------------------------------------------
% ---------------------  tire object ---------------------------
% --------------------------------------------------------------
% define 2D top-down view of vehicle body, starting from the nose, centered
% at geometric center of tire
% pt#        1       2     3         4      5     6         7      8     9        10     11     12
x_tire = (1/2)*[ 1   ,  0.98,  0.95,    -0.95, -0.98,   -1,      -1, -0.98, -0.95,    0.95,  0.98,     1];
y_tire = (1/2)*[ 0.60,  0.90,     1,        1,  0.90, 0.60,    -0.60, -0.90,    -1,      -1, -0.90, -0.60];

L_tire = obj_length;
W_tire = obj_width;
tire_obj.vertices = [L_tire*x_tire;
                     W_tire*y_tire];

% define the object faces
tire_obj.faces = 1:length(x_tire);






% Assign the function output
if object_type==1, % output vehicle vertices and faces
   out.vertices = veh_obj.vertices;
   out.faces    = veh_obj.faces;
elseif object_type==2, % output tire vertices and face
   out.vertices = [ tire_obj.vertices ];
   out.faces    = [ tire_obj.faces ];
end


% draw?
draw=0;
if draw==1,
   clf
   veh_handle = patch('Vertices',out.vertices','Faces',out.faces,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0 0 1],'FaceAlpha',0.1);
   xlabel('X-axis')
   ylabel('Y-axis')
   zlabel('Z-axis')
   axis equal
   grid on
end











