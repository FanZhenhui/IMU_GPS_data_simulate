%% -- statellite orbit just for show --
clear;clc;close all;

%% -- init data
DtoR=2*pi/360;
resG_countX = 1; 
tSum_stae = 12*2600;
vStae = 2 * pi * 26560*1000/tSum_stae;
dTheta = 360/12/3600/100*DtoR;
T = (745)*100*50;
matG1_G24_fs = zeros(T,3*24);
matV1_V24_fs = zeros(T,3*24);
axis_arrow_x = [0,0,0];
axis_arrow_y=[100000000,0,0];
for time=1:1500:T
    resG_countY = 1;
    rEarth = 26560000;
    eEarth = 0.00;
    E = [0:0.1:2*pi+0.05];
    x = rEarth*(cos(E)-eEarth);
    y = rEarth*sqrt((1-eEarth^2))*sin(E);
    z = 0*E;
    axisA1=[32.8 92.8 152.8 212.8 272.8 332.8 ];
    for k=1:6
        angA=axisA1(k)*DtoR;
        angB=55*DtoR;
        angC=pi/100;
        transR3=[cos(angA) -sin(angA) 0;
                 sin(angA)  cos(angA) 0;
                 0          0         1;];
        transR1=[1     0           0;
                 0     cos(angB)  -sin(angB);
                 0     sin(angB)  cos(angB);];
        transR2=[cos(angC) -sin(angC) 0;
                 sin(angC) cos(angC) 0; 
                 0  0  1;];
        L1=length(E);
        transR312=transR3*transR1*transR2;
        Ans=transR312*[x;y; z;];
        x1=Ans(1,:);
        y1=Ans(2,:);
        z1=Ans(3,:);
        plot3(x1,y1,z1,'k-');
        hold on;
        grid on;
    end
    Ctable = [10 50 160 260 ;
              80 180 220 320 ;
              10 130 250 340;
              50 150 170 300;
              100 210 310 340 ;
              140 150 240 350;];
    Wx = ones(1,1);Wy = ones(1,1);Wz = ones(1,1);
    for k = 1:6
        n = k+1;
        angA = axisA1(k)*DtoR; 
        angB = 55*DtoR;
        for m = 1:4
           angC = Ctable(k,m)*DtoR+dTheta*time;
           v    = [vStae*cos(angC);vStae*sin(angC);0]; 
           x    = rEarth*(cos(angC)-eEarth);
           y    = rEarth*sqrt((1-eEarth^2))*sin( angC);
           z    = 0*angC;
           transR3 = [cos(angA)  -sin(angA)  0;
                      sin(angA)  cos(angA) 0;
                      0        0     1;];
           transR1 = [1         0    0;
                      0       cos(angB)  -sin(angB);
                      0       sin(angB) cos(angB);];
           transR2 = [cos(angC) -sin(angC) 0;
                      sin(angC) cos(angC) 0; 
                      0          0  1;];
           L1 = length(E);
           transR312 = transR3*transR1*transR2;
           Ans = transR312*[x;y;z;];
           vAns = transR312*v;
           Wx = [Wx Ans(1,:)];
           Wy = [Wy Ans(2,:)];
           Wz = [Wz Ans(3,:)];
           x1 = Ans(1,:);
           y1 = Ans(2,:);
           z1 = Ans(3,:);
           DrawSatellite(x1,y1,z1);
           matG1_G24_fs(resG_countX,(resG_countY-1)*3 + 1) = x1;
           matG1_G24_fs(resG_countX,(resG_countY-1)*3 + 2) = y1;
           matG1_G24_fs(resG_countX,(resG_countY-1)*3 + 3) = z1;
           matV1_V24_fs(resG_countX,(resG_countY-1)*3 + 1) = vAns(1);
           matV1_V24_fs(resG_countX,(resG_countY-1)*3 + 2) = vAns(2);
           matV1_V24_fs(resG_countX,(resG_countY-1)*3 + 3) = vAns(3);
           resG_countY = resG_countY+1;
        end
    end
resG_countX = resG_countX+1;
hold off;
drawnow;
end
axis equal;
axis off;

function DrawSatellite(movex,movey,movez)
    plot3(movex,movey,movez,'b.','MarkerSize',20)
end