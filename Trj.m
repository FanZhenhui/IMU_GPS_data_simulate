function [att,Vel,Pos,Wibb,fb]=Trj(att,Vel,Pos,Wt,At,TS)
glvs;

%% 姿态更新
Cnt1=[cos(att(3)),cos(att(1))*sin(att(3)),-sin(att(1))*sin(att(3));
    -sin(att(3)),cos(att(1))*cos(att(3)),-sin(att(1))*cos(att(3));
    0,sin(att(1)),cos(att(1))
    ];  %不考虑横滚角的姿态矩阵
        %轨迹坐标系变换到导航坐标系，轨迹坐标系不考虑横滚角，y轴指向轨迹向前，x轴水平向右，z轴与xy轴成右手定则

% attm=att+Wt*TS/2;
attm(1)=att(1)+Wt(1)*TS/2;     %attm中储存这个TS周期一半时的姿态
attm(2)=att(2)+Wt(2)*TS/2;
attm(3)=att(3)-Wt(3)*TS/2;      %航向角方向和逆时针相反，航向角顺时针转动为正


Cntm=[cos(attm(3)),cos(attm(1))*sin(attm(3)),-sin(attm(1))*sin(attm(3));
    -sin(attm(3)),cos(attm(1))*cos(attm(3)),-sin(attm(1))*cos(attm(3));
    0,sin(attm(1)),cos(attm(1))
    ];  %轨迹坐标系变换到导航坐标系，轨迹坐标系不考虑横滚角，y轴指向轨迹向前，x轴水平向右，z轴与xy轴成右手定则

% att2=att+Wt*TS;
att2(1)=att(1)+Wt(1)*TS;      %att2中储存这个TS周期结束时的姿态
att2(2)=att(2)+Wt(2)*TS;
att2(3)=att(3)-Wt(3)*TS;      %航向角方向和逆时针相反，航向角顺时针转动为正

if att2(3)<0
   att2(3)=att2(3)+2*pi;
end
if att2(3)>2*pi
   att2(3)=att2(3)-2*pi;
end  
    
att=att2;
Cnt2=[cos(att2(3)),cos(att2(1))*sin(att2(3)),-sin(att2(1))*sin(att2(3));
    -sin(att2(3)),cos(att2(1))*cos(att2(3)),-sin(att2(1))*cos(att2(3));
    0,sin(att2(1)),cos(att2(1))
    ];  %轨迹坐标系变换到导航坐标系，轨迹坐标系不考虑横滚角，y轴指向轨迹向前，x轴水平向右，z轴与xy轴成右手定则
Cbn=[cos(att2(2))*cos(att2(3))+sin(att2(2))*sin(att2(3))*sin(att2(1)),-cos(att2(2))*sin(att2(3))+sin(att2(2))*cos(att2(3))*sin(att2(1)),-sin(att2(2))*cos(att2(1));
    sin(att2(3))*cos(att2(1)),cos(att2(3))*cos(att2(1)),sin(att2(1));
    sin(att2(2))*cos(att2(3))-cos(att2(2))*sin(att2(3))*sin(att2(1)),- sin(att2(2))*sin(att2(3))-cos(att2(2))*cos(att2(3))*sin(att2(1)),cos(att2(2))*cos(att2(1))
    ];  %导航坐标系变换到载体坐标系的变换矩阵v
Cba=[cos(att2(2)),0,-sin(att2(2))*cos(att2(1));
    0,1,sin(att2(1));
    sin(att2(2)),0,cos(att2(2))*cos(att2(1))
    ];  %姿态角速率变换到载体坐标系下的变换矩阵
        %Cba*Wt=Wbnb，即载体相对于导航系的角速度在载体系下的投影

%% 速度更新 四阶龙格库塔法
K1=Cnt1*At';
%Velm1=Vel+K1'*TS/2;
K2=Cntm*At';
%Velm2=Vel+K2'*TS/2;
K3=Cntm*At';
%Vel2=Vel+K3'*TS;
K4=Cnt2*At';
K=1/6*(K1+2*K2+2*K3+K4)';
Vel1=Vel;
Vel=Vel1+K*TS;
Velm=Vel1+K*TS/2;

%% 位置更新 四阶龙格库塔法
%经纬度及高度变化率（即速度）为[VN/RM、VE/(RN*cosL)、VU]，可以近似的认为RN=RM=Re
K1=[Vel1(2)/glv.Re;Vel1(1)/(glv.Re*cos(Pos(1)));Vel1(3)];
posm1=Pos+K1'*TS/2;
K2=[Velm(2)/glv.Re;Velm(1)/(glv.Re*cos(posm1(1)));Velm(3)];
posm2=Pos+K2'*TS/2;
K3=[Velm(2)/glv.Re;Velm(1)/(glv.Re*cos(posm2(1)));Velm(3)];
pos2=Pos+K3'*TS;
K4=[Vel(2)/glv.Re;Vel(1)/(glv.Re*cos(pos2(1)));Vel(3)];
%Pos1=Pos;
K=1/6*(K1+2*K2+2*K3+K4)';
Pos=Pos+K*TS;

%% 产生惯性器件信息
gn=[0,0,-glv.G];
Wnie=[0,glv.Wie*cos(Pos(1)),glv.Wie*sin(Pos(1))];
Wnen=[-Vel(2)/glv.Re,Vel(1)/glv.Re,Vel(1)*tan(Pos(1))/glv.Re];
Wibb=(Cbn*(Wnie+Wnen)'+Cba*Wt')';
fb=(Cbn*((Cnt2*At')'+cross(2*Wnie+Wnen,Vel)-gn)')'; %比力方程
Wibb=Wibb+0.01*pi/180/3600+0.001*pi/180/3600*randn(1,3);  
% Wibb(1)=Wibb(1)+10*pi/180/3600+0.001*pi/180/3600*randn(1,1); 
% Wibb(2)=Wibb(2)+0.01*pi/180/3600+0.001*pi/180/3600*randn(1,1); 
% Wibb(3)=Wibb(3)+0.01*pi/180/3600+0.001*pi/180/3600*randn(1,1); 
fb=fb+1e-4*glv.G+1e-5*glv.G*randn(1,3);
% Wibb=(Cbn*(Wnie+Wnen)'+Cba*Wt')';
% fb=(Cbn*((Cnt2*At')'+cross(2*Wnie+Wnen,Vel)-gn)')';
