%ATHIOANNHS 1444
function Glass_Cube_basic = Glass_Cube_basic(ang_inc)
%setting default values
n2=1.5; % Refractive index of glass cube
n1 = 1.0003; % Refractive index of air

disp('------------------------------')

while (ang_inc <0) || (ang_inc >89)
    ang_inc = input( 'angle of incident must be between 0-89 degrees\n' );
end

% create a figure with black background
figure;
set(gcf, 'Color' , [0 0 0]); % sets figure background color to RGB=[0 0 0] value (black)

% define the vertices of the square
vertices = [-1 -1; -1 1; 1 1; 1 -1];
% define the face of the square
faces = [1 2 3 4];

% create the square patch object
h = patch('Faces', faces, 'Vertices', vertices);
% set the square fill color to gray 
set(h, 'FaceColor', [0.5 0.5 0.5]);
hold on;
axis off;

% Calculate the position for the text
text_x = -1;
text_y = 1.2;
% add the ang_inc value to the plot
text(text_x, text_y, sprintf('Incident angle: %.2f degrees', ang_inc), ...
     'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', 12);

%creating the mirror on the side of the cube
plot( [1 1], [1 -1] , 'Color' , [1 1 1], 'Linewidth' , 10 ); % plots a thick white line

u = [cos(ang_inc*pi/180), sin(ang_inc*pi/180)]; % unit vector pointing in the direction of the line
% Calculate the starting point coordinates
x0 = -1 -2*u(1);
y0 =  - 2*u(2);

%creating the incidence array
plot( [x0 -1],[y0 0] , 'Color' , [1 1 1], 'Linewidth' , 2 ); % plots a thick white line

% Snell's Law at the n1->n2 interface:
refr_angle1 = asin((n1/n2) * sin(ang_inc*pi/180)) *180/pi; %[deg]

u2 = [cos(refr_angle1*pi/180), sin(refr_angle1*pi/180)]; % unit vector pointing in the direction of the line
% Calculate the starting point coordinates
x2 = -1 + 2*u2(1);
y2 = 2*u2(2);

disp('First array enters the Cube!')

%stopping rays so they touch the other surfaces of the cube
%1st ray hits the mirror
if (abs(x2 - 1) < 2.2) && (y2<1)
    
    x2=1;
    % Snell's Law at the n2->n1 interface:
    refr_angle2 = asin((n2/n1) * sin(refr_angle1*pi/180)) *180/pi; %[deg]

    %checking critical angle for total reflection
    crit_angle = asin(n2 / n1);
    if crit_angle <refr_angle2
        
        disp('We have total reflection!(second array hits mirror)');
        u3 = [cos(refr_angle2*pi/180), sin(refr_angle2*pi/180)]; % unit vector pointing in the direction of the line
        % Calculate the starting point coordinates
        x3 = x2 - 2*u3(1);
        y3 = y2 +  2* u3(2);
    end
    
    %stopping rays so they touch the other surfaces of the cube
    %checking if 3rd ray hits the left side of the cube
    if (x3 - 1 < 0.1) && (y3<1)
        
        disp('3rd ray hits 1st face of cube')
        x3 = -1;
        % Snell's Law at the n2->n1 interface:
        refr_angle3 = asin((n2/n1) * sin(refr_angle1*pi/180)) *180/pi; %[deg]

        %checking critical angle for total reflection
        crit_angle = asin(n2 / n1);
        if crit_angle < refr_angle3
            
            u4 = [cos(refr_angle2*pi/180), sin(refr_angle2*pi/180)]; % unit vector pointing in the direction of the line
            % Calculate the starting point coordinates
            x4 = x3 + 2*u4(1);
            y4 = y3 +  2* u4(2);
        end
        
        %checking if 4th ray hits left side of cube
        if (abs(x4 )) > 1
        
            disp('4th array hits 1st side of cube)');
            x4 = -1;
            y4 = 2*u4(1) + x4;
        end
    
    %checking if 4th ray hits upper side of the cube
        if y4 > 1
        
            disp('4th array hits upper side of cube)');
            y4 = 1;
            x4 = 2*u4(2) - y4;
        end
    
        plot( [x3 x4],[y3 y4] , 'Color' , [1 1 0], 'Linewidth' , 1 ); % plots a yellow line
    end
    
    %checking if 3rd ray hits upper side of cube
    if y3>1
        disp('3rd ray hits hits upper side of cube')
        y3 = 1;
            % Snell's Law at the n2->n1 interface:
        refr_angle3 = asin((n2/n1) * sin(refr_angle1*pi/180)) *180/pi; %[deg]

        %checking critical angle for total reflection
        crit_angle = asin(n2 / n1);
        if crit_angle < refr_angle3
            u4 = [cos(refr_angle3*pi/180), sin(refr_angle3*pi/180)]; % unit vector pointing in the direction of the line
            % Calculate the starting point coordinates
            x4 = x3 - 2*u4(1);
            y4 = y3 -  2* u4(2);
        end
        
        %checking if 4th ray hits the left side of the cube
        if x4 < -1
            
            disp('4th ray hits the left side of the cube')
            x4 = -1;
            y4 = 2* u4(1) +x4;
        end
        plot( [x3 x4],[y3 y4] , 'Color' , [1 1 0], 'Linewidth' , 1 ); % plots a yellow line
    end
    plot( [x2 x3],[y2 y3] , 'Color' , [1 0 0], 'Linewidth' , 1 ); % plots a red line       
    axis equal;
end

%2nd ray hits upper side of cube
if (abs(x2 - 1) > 0.1) && (y2>1)

    disp('second array hits the upper side of cube');
    y2=1;
    %Snell's Law at the n2->n1 interface:
    refr_angle2 = asin((n2/n1) * sin(refr_angle1*pi/180)) *180/pi; %[deg]

    %checking critical angle for total reflection
    crit_angle = asin(n2 / n1);
    if crit_angle <refr_angle2
        u3 = [cos(refr_angle1*pi/180), sin(refr_angle1*pi/180)]; % unit vector pointing in the direction of the line
        % Calculate the starting point coordinates
        x3 = x2 + 2*u3(1);
        y3 = y2 -  2* u3(2);
    end
    
    %stopping rays so they touch the other surfaces of the cube
    %checking if 3rd ray hits mirror
    if (x3>1) && (y3<1)
        
        disp('3rd ray hits mirror')
        x3=1;
        y3 = 2* u3(1) - x3;
        
        %Snell's Law at the n2->n1 interface:    
        refr_angle3 = asin((n2/n1) * sin(refr_angle2*pi/180)) *180/pi; %[deg]

        %checking critical angle for total reflection
        crit_angle = asin(n2 / n1);
        if crit_angle <refr_angle3
            u4 = [cos(refr_angle2*pi/180), sin(refr_angle2*pi/180)]; % unit vector pointing in the direction of the line
            % Calculate the starting point coordinates
            x4 = x3 - 2*u4(1);
            y4 = y3 -  2* u4(2);
        end
        
        %checking if 4th ray hits lower side of cube
        if y4< -1
            disp('4th ray hits lower side of cube')
            y4 = -1;
            x4 = -2*u4(2) - y4;
        end
        plot( [x3 x4],[y3 y4] , 'Color' , [1 1 0], 'Linewidth' , 1 ); % plots a yellow line       

    end
%     
%     %checking if 2nd ray hits the lower side of the cube
%     if abs(y3 - 1)>0.1
%          
%         disp('3rd ray hits the lower side of the cube')
%         y3= -1;
%         x3 = 2* u3(2) + y3;
%     end
    plot( [x2 x3],[y2 y3] , 'Color' , [1 0 0], 'Linewidth' , 1 ); % plots a red line       
end

%creating reflaction array
plot( [-1 x2],[0 y2] , 'Color' , [0 1 0], 'Linewidth' , 1 ); % plots a green line
axis equal

%checking where the 4th part of the ray hits and creating the part of refraction
if x4 == -1
    refr_angle4 = asin((n2/n1) * sin(refr_angle3*pi/180)) *180/pi; %[deg]
    u5 = [cos(refr_angle4*pi/180), sin(refr_angle4*pi/180)]; % unit vector pointing in the direction of the line
    x5 = x4 - 1*u5(1);
    y5 = y4 -  1* u5(2);
    
elseif y4 == -1 
    refr_angle4 = asin((n2/n1) * sin(refr_angle3*pi/180)) *180/pi; %[deg]
    u5 = [cos(refr_angle4*pi/180), sin(refr_angle4*pi/180)] % unit vector pointing in the direction of the line
    x5 = x4 - 1*u5(1);
    y5 = y4 -  1* u5(2);
    
elseif y4 == 1 
    refr_angle4 = asin((n2/n1) * sin(refr_angle3*pi/180)) *180/pi; %[deg]
    u5 = [cos(refr_angle4*pi/180), sin(refr_angle4*pi/180)]; % unit vector pointing in the direction of the line
    x5 = x4 + 1*u5(1);
    y5 = y4 +  1* u5(2);
end

%plotting the refraction part of the ray 
plot( [x4 x5],[y4 y5] , 'Color' , 'c', 'Linewidth' , 2 ); % plots a cyan  line 

disp('ray exits the cube!')
disp('-------------------------')
end