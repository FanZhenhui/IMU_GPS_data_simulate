function [att,Vel,Pos,Wibb,fb]=Trj(att,Vel,Pos,Wt,At,TS)
glvs;

%% ��̬����
Cnt1=[cos(att(3)),cos(att(1))*sin(att(3)),-sin(att(1))*sin(att(3));
    -sin(att(3)),cos(att(1))*cos(att(3)),-sin(att(1))*cos(att(3));
    0,sin(att(1)),cos(att(1))
    ];  %�����Ǻ���ǵ���̬����
        %�켣����ϵ�任����������ϵ���켣����ϵ�����Ǻ���ǣ�y��ָ��켣��ǰ��x��ˮƽ���ң�z����xy������ֶ���

% attm=att+Wt*TS/2;
attm(1)=att(1)+Wt(1)*TS/2;     %attm�д������TS����һ��ʱ����̬
attm(2)=att(2)+Wt(2)*TS/2;
attm(3)=att(3)-Wt(3)*TS/2;      %����Ƿ������ʱ���෴�������˳ʱ��ת��Ϊ��


Cntm=[cos(attm(3)),cos(attm(1))*sin(attm(3)),-sin(attm(1))*sin(attm(3));
    -sin(attm(3)),cos(attm(1))*cos(attm(3)),-sin(attm(1))*cos(attm(3));
    0,sin(attm(1)),cos(attm(1))
    ];  %�켣����ϵ�任����������ϵ���켣����ϵ�����Ǻ���ǣ�y��ָ��켣��ǰ��x��ˮƽ���ң�z����xy������ֶ���

% att2=att+Wt*TS;
att2(1)=att(1)+Wt(1)*TS;      %att2�д������TS���ڽ���ʱ����̬
att2(2)=att(2)+Wt(2)*TS;
att2(3)=att(3)-Wt(3)*TS;      %����Ƿ������ʱ���෴�������˳ʱ��ת��Ϊ��

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
    ];  %�켣����ϵ�任����������ϵ���켣����ϵ�����Ǻ���ǣ�y��ָ��켣��ǰ��x��ˮƽ���ң�z����xy������ֶ���
Cbn=[cos(att2(2))*cos(att2(3))+sin(att2(2))*sin(att2(3))*sin(att2(1)),-cos(att2(2))*sin(att2(3))+sin(att2(2))*cos(att2(3))*sin(att2(1)),-sin(att2(2))*cos(att2(1));
    sin(att2(3))*cos(att2(1)),cos(att2(3))*cos(att2(1)),sin(att2(1));
    sin(att2(2))*cos(att2(3))-cos(att2(2))*sin(att2(3))*sin(att2(1)),- sin(att2(2))*sin(att2(3))-cos(att2(2))*cos(att2(3))*sin(att2(1)),cos(att2(2))*cos(att2(1))
    ];  %��������ϵ�任����������ϵ�ı任����v
Cba=[cos(att2(2)),0,-sin(att2(2))*cos(att2(1));
    0,1,sin(att2(1));
    sin(att2(2)),0,cos(att2(2))*cos(att2(1))
    ];  %��̬�����ʱ任����������ϵ�µı任����
        %Cba*Wt=Wbnb������������ڵ���ϵ�Ľ��ٶ�������ϵ�µ�ͶӰ

%% �ٶȸ��� �Ľ����������
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

%% λ�ø��� �Ľ����������
%��γ�ȼ��߶ȱ仯�ʣ����ٶȣ�Ϊ[VN/RM��VE/(RN*cosL)��VU]�����Խ��Ƶ���ΪRN=RM=Re
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

%% ��������������Ϣ
gn=[0,0,-glv.G];
Wnie=[0,glv.Wie*cos(Pos(1)),glv.Wie*sin(Pos(1))];
Wnen=[-Vel(2)/glv.Re,Vel(1)/glv.Re,Vel(1)*tan(Pos(1))/glv.Re];
Wibb=(Cbn*(Wnie+Wnen)'+Cba*Wt')';
fb=(Cbn*((Cnt2*At')'+cross(2*Wnie+Wnen,Vel)-gn)')'; %��������
Wibb=Wibb+0.01*pi/180/3600+0.001*pi/180/3600*randn(1,3);  
% Wibb(1)=Wibb(1)+10*pi/180/3600+0.001*pi/180/3600*randn(1,1); 
% Wibb(2)=Wibb(2)+0.01*pi/180/3600+0.001*pi/180/3600*randn(1,1); 
% Wibb(3)=Wibb(3)+0.01*pi/180/3600+0.001*pi/180/3600*randn(1,1); 
fb=fb+1e-4*glv.G+1e-5*glv.G*randn(1,3);
% Wibb=(Cbn*(Wnie+Wnen)'+Cba*Wt')';
% fb=(Cbn*((Cnt2*At')'+cross(2*Wnie+Wnen,Vel)-gn)')';
