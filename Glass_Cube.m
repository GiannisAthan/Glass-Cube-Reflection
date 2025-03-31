% ATHIOANNIS 1444
% n1: Refractive index of outside material
% n2: Refractive index of cube 
% dim: dimentions of cube in cm
% y: coordinates of arrays impact point on cube in cm
%ang_inc: angle of incidence in degrees (θin)
function Glass_Cube = Glass_Cube(n1, n2 ,dim , y , ang_inc)
%setting default values
% n2=1.5; % Refractive index of glass cube
% n1 = 1.0003; % Refractive index of air

while (ang_inc <=1) && (ang_inc >=89)
    ang_inc = input( 'angle of incident must be between 1-89 degrees\n' );
end

while dim/2 < 0.01
    dim = input( 'The dimencions must be inside the optics field (above 0.01 cm)\n' );
end

while (dim/2 > -y) && (dim/2 < y)
    y = input( 'y coordinate of arrays impact point must be on the cube\n' );
end

% create a figure with black background
figure;
set(gcf, 'Color' , [0 0 0]); % sets figure background color to RGB=[0 0 0] value (black)

% define the vertices of the square
vertices = [-dim/2 -dim/2; -dim/2 dim/2; dim/2 dim/2; dim/2 -dim/2];
% define the face of the square
faces = [1 2 3 4];

% create the square patch object
h = patch('Faces', faces, 'Vertices', vertices);
% set the square fill color to gray 
set(h, 'FaceColor', [0.5 0.5 0.5]);
hold on;
axis off;

%creating the mirror on the side of the cube
plot( [dim/2 dim/2], [dim/2 -dim/2] , 'Color' , [1 1 1], 'Linewidth' , 10 ); % plots a thick white line

u = [cos(ang_inc*pi/180), sin(ang_inc*pi/180)]; % unit vector pointing in the direction of the line
% Calculate the starting point coordinates
d = dim/2; % distance from ending point to starting point (1 cm)
x0 = -dim/2 - d * u(1);
y0 = y - d * u(2);

%creating the incidence array
plot( [x0 -dim/2],[y0 y] , 'Color' , [1 1 1], 'Linewidth' , 2 ); % plots a thick white line

%checking critical angle for total reflection
crit_angle = asin(n2 / n1);
if (n1 > n2) && (crit_angle < ang_inc)
    disp('We have total reflection!');
    u1 = [cos(crit_angle*pi/180), sin(crit_angle*pi/180)]; % unit vector pointing in the direction of the line
    % Calculate the starting point coordinates
    d = dim/2; % distance from ending point to starting point (1 cm)
    x1 = -dim/2 - d * u1(1);
    y1 = y - d * u1(2);

    %creating reflaction array
    plot( [x1 -dim/2],[y1 y] , 'Color' , [1 0 0], 'Linewidth' , 1 ); % plots red line
end

%checking Brewster angle for total     
Br_ang = atan(n2/n1);
if Br_ang == ang_inc
    disp('We have total refraction!')
        u1 = [cos(Br_ang*pi/180), sin(Br_ang*pi/180)]; % unit vector pointing in the direction of the line
    % Calculate the starting point coordinates
    d = dim/2; % distance from ending point to starting point (1 cm)
    x1 = -dim/2 - d * u1(1);
    y1 = y - d * u1(2);

    %creating reflaction array
    plot( [x1 -dim/2],[y1 y] , 'Color' , [1 1 0], 'Linewidth' , 1 ); % plots a yellow line
end

% Snell's Law at the n1->n2 interface:
refr_angle =asin((n1/n2) * sin(ang_inc*pi/180)); %[rad]

u2 = [cos(refr_angle*pi/180), sin(refr_angle*pi/180)]; % unit vector pointing in the direction of the line
% Calculate the starting point coordinates
d = dim; % distance from ending point to starting point (1 cm)
x2 = -dim/2 + d * u2(1);
y2 = y + d * u2(2);

%creating reflaction array
plot( [-dim/2 x2],[y y2] , 'Color' , [0 1 0], 'Linewidth' , 1 ); % plots a green line

%distance of array changed to dim to check if the array 'hits' the mirror 
if x2 >= dim/2 - 0.2
      disp('We have total reflection!');
    u1 = [cos(crit_angle*pi/180), sin(crit_angle*pi/180)]; % unit vector pointing in the direction of the line
    % Calculate the starting point coordinates
    d = dim/2; % distance from ending point to starting point (1 cm)
    x1 = -dim/2 - d * u1(1);
    y1 = y - d * u1(2);

    %creating reflaction array
    plot( [x1 -dim/2],[y1 y] , 'Color' , [1 0 0], 'Linewidth' , 1 ); % plots red line
end
end