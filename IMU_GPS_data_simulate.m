%% Step.1 GPS 位置与速率
clear;clc;close all;
DtoR=2*pi/360;
resG_countX = 1; 
tSum_stae = 12*2600;
vStae = 2 * pi * 26560*1000/tSum_stae;
dTheta = 360/12/3600/100*DtoR;
T = (745)*100;
matG1_G24 = zeros(T,3*24);
matV1_V24 = zeros(T,3*24);
hwait = waitbar(0,'IMU-GPS-data-simulate'); 
for time=0:T
    resG_countY = 1;
    a=26560000;
    e=0.02;
    E=[0:0.1:2*pi];
    x=a*(cos(E)-e);
    y=a*sqrt((1-e^2))*sin(E);
    z=0*E;
    A1=[32.8 92.8 152.8 212.8 272.8 332.8 ];
    for k=1:6
        A=A1(k)*DtoR;
        B=55*DtoR;
        C=pi/100;
        R3=[cos(A) -sin(A) 0;
        sin(A)  cos(A) 0;
        0        0     1;];
    R1=[1         0    0;
        0       cos(B)  -sin(B);
        0 sin(B) cos(B);];
    R2=[cos(C) -sin(C) 0;
        sin(C) cos(C) 0; 
        0  0  1;];
    L1=length(E);
    R312=R3*R1*R2;
    Ans=R312*[x;y; z;];
    x1=Ans(1,:);
    y1=Ans(2,:);
    z1=Ans(3,:);
end
Ctable=[10 50 160 260 ;
        80 180 220 320 ;
        10 130 250 340;
        50 150 170 300;
        100 210 310 340 ;
        140 150 240 350;];
    
Wx=ones(1,1);Wy=ones(1,1);Wz=ones(1,1);
   
for k=1:6
     n=k+1;
     A=A1(k)*DtoR; 
     B=55*DtoR;
     for m=1:4
           C=Ctable(k,m)*DtoR+dTheta*time;
           v = [vStae*cos(C);vStae*sin(C);0]; 
           x=a*(cos(C)-e);
           y=a*sqrt((1-e^2))*sin( C);
           z=0*C;
           R3=[cos(A)  -sin(A)  0;
                sin(A)  cos(A) 0;
                 0        0     1;];
           R1=[1         0    0;
               0       cos(B)  -sin(B);
               0       sin(B) cos(B);];
           R2=[cos(C) -sin(C) 0;
               sin(C) cos(C) 0; 
               0          0  1;];
      L1=length(E);
      R312=R3*R1*R2;
      Ans=R312*[x;y;z;];
      vAns = R312*v;
      Wx=[Wx Ans(1,:)];
      Wy=[Wy Ans(2,:)];
      Wz=[Wz Ans(3,:)];
      x1=Ans(1,:);
      y1=Ans(2,:);
      z1=Ans(3,:);
      matG1_G24(resG_countX,(resG_countY-1)*3 + 1) = x1;
      matG1_G24(resG_countX,(resG_countY-1)*3 + 2) = y1;
      matG1_G24(resG_countX,(resG_countY-1)*3 + 3) = z1;
      matV1_V24(resG_countX,(resG_countY-1)*3 + 1) = vAns(1);
      matV1_V24(resG_countX,(resG_countY-1)*3 + 2) = vAns(2);
      matV1_V24(resG_countX,(resG_countY-1)*3 + 3) = vAns(3);
      resG_countY = resG_countY+1;
     end
  end
resG_countX = resG_countX+1;
end
save('matG1_G24.mat','matG1_G24');
save('matV1_V24.mat','matV1_V24');
disp('--- Step.1 卫星位置速度生成 ---');
%% Step.2 IMU 比力，角速度，导航位置速度姿态仿真值及真值
glvs;             
TS=0.01;          

initPos=[10*glv.D2R,110*glv.D2R,20];                              
vb=[0,0,0];                                                       
inittheta=0*glv.D2R;  initgama=0*glv.D2R;  initpsi=0.0001*glv.D2R;      
initVel=(Trans_att2attm([inittheta initgama initpsi])*vb')';                   

att=[inittheta,initgama,initpsi];
Vel=initVel;
Pos=initPos;

W1=90*glv.D2R/50; 
V1=10;            
A1=V1*W1;        

W2=30*glv.D2R/10; 
V2=10;  
A2=V2*W2;

W3=40*glv.D2R/10;

w_1 = 90*glv.D2R/5;
v_1 = 10;
A_1 = v_1 * w_1;
Para=[% 持续时间，姿态角速率Wx Wy Wz，载体加速度在载体系下的投影Ax Ay Az； 
    300, 0,0,0, 0,0,0;
    10,  0,0,0, 1,0,0;
    120, 0,0,0, 0,0,0;
    5,   0,0,w_1, 0,A_1,0;
    60, 0,0,0, 0,0,0;
    5, 0,0,-w_1, -A_1,0,0;
    120, 0,0,0, 0,0,0;
    5, 0,0,w_1, 0,-A_1,0;
    120, 0,0,0, 0,0,0;
        ];
% IMU and GPS raw data simulate   
profile(1,1:21)=[att*glv.R2D,initVel,initPos,zeros(1,12)];
kk=1;
len=sum(Para(:,1))/TS;

countPx = 1;
for k=1:size(Para,1)
    Wt=Para(k,2:4);
    At=Para(k,5:7);
    for j=1:Para(k,1)/TS
    [att,Vel,Pos,Wibb,fb]=Trj(att,Vel,Pos,Wt,At,TS);
    [x_u,y_u,z_u] = TransN2Ecef(Pos(1),Pos(2),Pos(3));
    matPosEcef(countPx,:) = [x_u,y_u,z_u];
    for countP = 1:24
        matPI(countPx,countP) = sqrt((x_u - matG1_G24(countPx,countP*3-2))^2 + (y_u - matG1_G24(countPx,countP*3-1))^2+(z_u - matG1_G24(countPx,countP*3-0))^2); 
    end
    countPx = countPx+1;
    Pos_GPS(1)=Pos(1)+2.5/glv.Re*randn(1,1); 
    Pos_GPS(2)=Pos(2)+2.5/glv.Re*randn(1,1); 
    Pos_GPS(3)=Pos(3)+5*randn(1,1);          
    profile(kk,1:21)=[att*glv.R2D,Vel,Pos,Wibb,fb,vb,Pos_GPS];
    kk=kk+1;
    waitbar(kk/len);
    end
end

Profile(:,1:6)=profile(:,10:15);            
Profile(:,7:15)=profile(:,1:9); 
Profile(:,16:21)=profile(:,16:21); 
imu=Profile(1:len,:);                     
%save('trace.dat','imu','-ascii','-double'); % IMU数据保存成dat格式文档，双精度
% save('20210316_745s.txt','imu','-ascii','-double');  % IMU数据保存成txt格式文档，双精度

% plot result
% t=[1:length(profile)]*TS;
% figure name '姿态角'
% subplot(311);plot(t,profile(:,1));xlabel('t/s');ylabel('pitch/°');title('俯仰角');grid on;
% subplot(312);plot(t,profile(:,2));xlabel('t/s');ylabel('roll/°');title('横滚角');grid on;
% subplot(313);plot(t,profile(:,3));xlabel('t/s');ylabel('heading/°');title('航向角');grid on;
% 
% figure name '速度'
% subplot(311);plot(t,profile(:,4));xlabel('t/s');ylabel('VE/m/s');title('东向速度');grid on;
% subplot(312);plot(t,profile(:,5));xlabel('t/s');ylabel('VNm/s');title('北向速度');grid on;
% subplot(313);plot(t,profile(:,6));xlabel('t/s');ylabel('VU/m/s');title('天向速度');grid on;
% 
% figure name '位置'
% subplot(311);plot(t,profile(:,7)*glv.R2D);xlabel('t/s');ylabel('\itL\rm/°');title('纬度');grid on;
% subplot(312);plot(t,profile(:,8)*glv.R2D);xlabel('t/s');ylabel('\it\lambda\rm/°');title('经度');grid on;
% subplot(313);plot(t,profile(:,9));xlabel('t/s');ylabel('\ith/m');title('高度');grid on;
% 
% figure name 'GPS位置'
% subplot(311);plot(t,profile(:,19)*glv.R2D);xlabel('t/s');ylabel('\itL\rm/°');title('GPS纬度');grid on;
% subplot(312);plot(t,profile(:,20)*glv.R2D);xlabel('t/s');ylabel('\it\lambda\rm/°');title('GPS经度');grid on;
% subplot(313);plot(t,profile(:,21));xlabel('t/s');ylabel('\ith/m');title('GPS高度');grid on;
% 
% figure name 'GPS经纬度'
% plot((profile(:,20)*glv.R2D),profile(:,19)*glv.R2D);xlabel('经度');ylabel('维度\rm/°');title('GPS纬度');grid on;
save('IMUdata.mat','imu');
save('matPI.mat','matPI');        
disp('--- Step.2 惯导数据、载体伪距生成 ---');       
%% Step.3 卫星伪距伪距率、误差值
% PIs, tu, tru
[lenX_,lenY] = size(matG1_G24);
matTru = zeros(lenX_-1,1); 
matTu  = zeros(lenX_-1,1);
countP = 0;
Trj_Tru_a = 100;
matPIs   = zeros(lenX_-1,24);

for i = 1:lenX_-1
    %-- step.1: cal param
    true_PosIMU = [imu(i,13),imu(i,14),imu(i,15)];
    Rm=glv.Re*(1-2*glv.e+3*glv.e*sin(true_PosIMU(1))^2);
    Rn=glv.Re*(1+glv.e*sin(true_PosIMU(1))^2);
    e = glv.e;
    sL  = sin(true_PosIMU(1)); cL = cos(true_PosIMU(1));  sL2 = sL^2;   
    RMh = Rm + true_PosIMU(3); RNh = Rn + true_PosIMU(3);
    RN = Rn;
    sJ  = sin(true_PosIMU(2)); cJ = cos(true_PosIMU(2)); H = true_PosIMU(3);
    CneT = [-sL*cJ,cL*cJ,-sJ;-sL*sJ,cL*sJ,cJ;cL,sL,0]; 
    Cne  = [CneT(3,3),CneT(3,1),CneT(3,2);
            CneT(1,3),CneT(1,1),CneT(1,2);
            CneT(2,3),CneT(2,1),CneT(2,2)];
        
    %-- step.2: ENU->ecef
    true_VnIMU = [imu(i,10);imu(i,11);imu(i,12)];
    XI_star = Cne*true_VnIMU;
    
    %-- step.3: cal weijulv
    for j = 1:24
        [xI(1),xI(2),xI(3)] = TransN2Ecef(true_PosIMU(1),true_PosIMU(2),true_PosIMU(3));
        rI = sqrt((xI(1)-matG1_G24(i,j*3-2))^2+(xI(2)-matG1_G24(i,j*3-1))^2+(xI(3)-matG1_G24(i,j*3-0))^2);
        E  = [xI(1)-matG1_G24(i,j*3-2),xI(2)-matG1_G24(i,j*3-1),xI(3)-matG1_G24(i,j*3-0)]/rI;
        matPIs(i,j) = E(1)*(XI_star(1)-matV1_V24(i,j*3-2))+E(2)*(XI_star(2)-matV1_V24(i,j*3-1))+E(3)*(XI_star(3)-matV1_V24(i,j*3-0));
    end
end

for countX = 1:(lenX_-1)/100
    Trj_Tru_sumA = 0;
    for countY = 1:100
          matTru((countX-1)*100+countY,1) = Trj_Tru_a;% + 100*(1/(countY/100+1))-100;
          matTu ((countX-1)*100+countY,1) = Trj_Tru_a*((countX-1)*100+countY)*0.01;% + 100*log(countY/100+1);
    end
end
save('matPIs.mat','matPIs');
save('matTru.mat','matTru');  
save('matTu.mat','matTu');   
disp('--- Step.3 载体伪距率、误差生成 ---');        
