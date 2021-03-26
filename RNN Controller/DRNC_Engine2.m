%将DRNC用于发动机转子转速的控制
%更改输出
%Author:Taylon
%Email:taylonupup@163.com
%Create: Dec 27, 2020
%------------------------------------------------------------
clearvars;clc;

ALR=readfis('LR_engine2');

A=[-2.01 2.06;-0.0044 -2.92];
B=[0.84 1.2;0.53 0.28];
C=[1 0;-0.86 0.51];
D=[0 0;-0.21 0.56];
engine=ss(A,B,C,D);
Ts=0.04;
engineD=c2d(engine,Ts);
%step(engine,'-',engineD,'--')

%初始化控制对象
%ts=0.001;   %采样步长
L=1000;     %仿真时长
r=zeros(L,1);
yr=zeros(L,1);
y=zeros(L,1);
u=zeros(L,1);
ec=zeros(L,1);
em=zeros(L,1);
yu=zeros(L,1);
ym=zeros(L,1);
alpha=0;  %动量系数
etaC=0.2;   %学习率
etaI=0.4;
%maxS=2;     %灵敏度最大值

%初始化DRNC(控制器)
Nc=6;pc=3;qc=1;
wc_in1=zeros(Nc,L);wc_in1(:,1)=0.1*ones(Nc,1);
wc_in2=zeros(Nc,L);wc_in2(:,1)=0.1*ones(Nc,1);
wc_in3=zeros(Nc,L);wc_in3(:,1)=0.1*ones(Nc,1);
wc_r=zeros(Nc,L);wc_r(:,1)=0.1*ones(Nc,1);
wc_out=zeros(Nc,L);wc_out(:,1)=0.1*ones(Nc,1);
Sc=zeros(Nc,L);diffSc=zeros(Nc,L);
Hc=zeros(Nc,L);
bhc=zeros(Nc,L);
boc=0.1*ones(1,L);
%权值的分配:[p*N N N q*N]，一共分为6个部分
gradWc=0.01*ones((pc+qc+2)*Nc,L);
deltaWc=0.01*ones((pc+qc+2)*Nc,L);
gradBoc=0.01*ones(1,L);
deltaBoc=0.2*ones(1,L);
etac=0.08*ones(6,L);
p1c=Nc;p2c=2*Nc;
p3c=3*Nc;p4c=4*Nc;p5c=5*Nc;
Pc=zeros(Nc,L);
Q1c=zeros(Nc,L);Q2c=zeros(Nc,L);Q3c=zeros(Nc,L);


%初始化DRNI(识别器)
Ni=4;pi=2;qi=1;
wi_in1=zeros(Ni,L);wi_in1(:,1)=0.1*ones(Ni,1);
wi_in2=zeros(Ni,L);wi_in2(:,1)=0.1*ones(Ni,1);
wi_r=zeros(Ni,L);wi_r(:,1)=0.1*ones(Ni,1);
wi_out=zeros(Ni,L);wi_out(:,1)=0.1*ones(Ni,1);
Si=0.1*ones(Ni,L);diffSi=0.1*ones(Ni,L);
Hi=zeros(Ni,L);
bi=0.1*ones(Ni,L);
%权值的分配:[p*N N N q*N]，一共分为5个部分
gradWi=zeros((pi+qi+2)*Ni,L);
deltaWi=zeros((pi+qi+2)*Ni,L);
etai=zeros(5,L);
p1i=Ni;p2i=2*Ni;
p3i=3*Ni;p4i=4*Ni;
Pi=zeros(Ni,L);
Q1i=zeros(Ni,L);Q2i=zeros(Ni,L);

PITr=[ones(L/4,1);zeros(L/4,1);ones(L/4,1);zeros(L/4,1)];     %目标转速
nL=zeros(L,1);      %实际低压转速
nH=zeros(L,1);      %实际高压转速
PIT=zeros(L,1);      %实际涡轮落压比
wf=zeros(L,1);      %燃油量

%仿真开始
for k=2:1:L
    %控制对象更新
    r(k)=PITr(k);
    yr(k)=PITr(k);
    
    xk=engineD.A*[nL(k-1);nH(k-1)]+engineD.B*[wf(k);0];
    nL(k)=xk(1);nH(k)=xk(2);
    yk=engineD.C*xk+engineD.D*[wf(k);0];
    nL(k)=yk(1);PIT(k)=yk(2);
    
    y(k)=PIT(k);u(k)=wf(k);
    
    %DRNC更新
    Hc(:,k)=wc_in1(:,k-1)*r(k)+wc_in2(:,k-1)*u(k-1)+wc_in3(:,k-1)*y(k-1)+wc_r(:,k-1).*Sc(:,k-1)+bhc(:,k-1);
    Sc(:,k)=Acfun(Hc(:,k));
    u(k)=(wc_out(:,k-1))'*Sc(:,k)+boc(k);
    wf(k+1)=u(k);
    
    %DRNI更新
    Hi(:,k)=wi_in1(:,k-1)*u(k)+wi_in2(:,k-1)*y(k-1)+wi_r(:,k-1).*Si(:,k-1)+bi(:,k-1);
    Si(:,k)=Acfun(Hi(:,k));
    ym(k)=(wi_out(:,k-1))'*Si(:,k);
    
    ec(k)=yr(k)-y(k);
    em(k)=y(k)-ym(k);
    %反馈阶段
    %求取灵敏度yu(k)
    yu(k)=wi_out(:,k-1)'*((1-Si(:,k).^2).*wi_in1(:,k-1));

    
    %计算DRNC权值梯度
    gradWc(p5c+1:end,k)=ec(k)*yu(k)*Sc(:,k);     %输出权值梯度
    diffSc(:,k)=diffAcfun(Hc(:,k));
    gradWc(p4c+1:p5c,k)=ec(k)*yu(k)*wc_out(:,k-1).*diffSc(:,k);   %隐含层偏置梯度
    Pc(:,k)=diffSc(:,k).*(Sc(:,k-1)+wc_r(:,k-1).*Pc(:,k-1));
    gradWc(p3c+1:p4c,k)=ec(k)*yu(k)*wc_out(:,k-1).*Pc(:,k);       %循环权值梯度
    Q3c(:,k)=diffSc(:,k).*(y(k-1)+wc_r(:,k-1).*Q3c(:,k-1));
    gradWc(p2c+1:p3c,k)=ec(k)*yu(k)*wc_out(:,k-1).*Q3c(:,k);       %第三个输入权值梯度
    Q2c(:,k)=diffSc(:,k).*(u(k-1)+wc_r(:,k-1).*Q2c(:,k-1));
    gradWc(p1c+1:p2c,k)=ec(k)*yu(k)*wc_out(:,k-1).*Q2c(:,k);       %第二个输入权值梯度
    Q1c(:,k)=diffSc(:,k).*(r(k)+wc_r(:,k-1).*Q1c(:,k-1)); 
    gradWc(1:p1c,k)=ec(k)*yu(k)*wc_out(:,k-1).*Q1c(:,k);       %第一个输入权值梯度
    gradBoc(k)=ec(k)*yu(k);
    %计算DRNC自适应学习率
    etac(6,k)=evalfis([abs(ec(k)) max(abs(gradWc(p5c+1:end,k)))],ALR);
    etac(5,k)=evalfis([abs(ec(k)) max(abs(gradWc(p4c+1:p5c,k)))],ALR);
    etac(4,k)=evalfis([abs(ec(k)) max(abs(gradWc(p3c+1:p4c,k)))],ALR);
    etac(3,k)=evalfis([abs(ec(k)) max(abs(gradWc(p2c+1:p3c,k)))],ALR);
    etac(2,k)=evalfis([abs(ec(k)) max(abs(gradWc(p1c+1:p2c,k)))],ALR);
    etac(1,k)=evalfis([abs(ec(k)) max(abs(gradWc(1:p1c,k)))],ALR);
    %etac(:,k)=etaC*(2*ec(k))^2;

    %计算权值变化量
    deltaWc(p5c+1:end,k)=etac(6,k)*gradWc(p5c+1:end,k);
    deltaWc(p4c+1:p5c,k)=etac(5,k)*gradWc(p4c+1:p5c,k);
    deltaWc(p3c+1:p4c,k)=etac(4,k)*gradWc(p3c+1:p4c,k);
    deltaWc(p2c+1:p3c,k)=etac(3,k)*gradWc(p2c+1:p3c,k);
    deltaWc(p1c+1:p2c,k)=etac(2,k)*gradWc(p1c+1:p2c,k);
    deltaWc(1:p1c,k)=etac(1,k)*gradWc(1:p1c,k);
    deltaBoc(k)=etac(6,k)*gradBoc(k);
    
    %权值更新
    wc_out(:,k)=wc_out(:,k-1)+deltaWc(p5c+1:end,k);
    bhc(:,k)=bhc(:,k-1)+deltaWc(p4c+1:p5c,k);
    wc_r(:,k)=wc_r(:,k-1)+deltaWc(p3c+1:p4c,k);
    wc_in3(:,k)=wc_in3(:,k-1)+deltaWc(p2c+1:p3c,k);
    wc_in2(:,k)=wc_in2(:,k-1)+deltaWc(p1c+1:p2c,k);
    wc_in1(:,k)=wc_in1(:,k-1)+deltaWc(1:p1c,k);
    boc(k)=boc(k-1)+deltaBoc(k);
    
    %计算DRNI权值梯度
    gradWi(p4i+1:end,k)=em(k)*Si(:,k);     %输出权值梯度
    diffSi(:,k)=diffAcfun(Hi(:,k));
    gradWi(p3i+1:p4i,k)=em(k)*wi_out(:,k-1).*diffSi(:,k);   %隐含层偏置梯度
    Pi(:,k)=diffSi(:,k).*(Si(:,k-1)+wi_r(:,k-1).*Pi(:,k-1));
    gradWi(p2i+1:p3i,k)=em(k)*wi_out(:,k-1).*Pi(:,k);       %循环权值梯度
    Q2i(:,k)=diffSi(:,k).*(y(k-1)+wi_r(:,k-1).*Q2i(:,k-1));
    gradWi(p1i+1:p2i,k)=em(k)*wi_out(:,k-1).*Q2i(:,k);       %第二个输入权值梯度
    Q1i(:,k)=diffSi(:,k).*(u(k)+wi_r(:,k-1).*Q1i(:,k-1)); 
    gradWi(1:p1i,k)=em(k)*wi_out(:,k-1).*Q1i(:,k);       %第一个输入权值梯度
    %计算DRNI自适应学习率
%     etai(5,k)=1./max(gradWi(p4i+1:end,k).^2);
%     etai(4,k)=1./max(gradWi(p3i+1:p4i,k).^2);
%     etai(3,k)=1./max(gradWi(p2i+1:p3i,k).^2);
%     etai(2,k)=1./max(gradWi(p1i+1:p2i,k).^2);
%     etai(1,k)=1./max(gradWi(1:p1i,k).^2);
    etai(:,k)=etaI;
    %计算权值变化量
    deltaWi(p4i+1:end,k)=etai(5,k)*gradWi(p4i+1:end,k)+alpha*deltaWi(p4i+1:end,k-1);
    deltaWi(p3i+1:p4i,k)=etai(4,k)*gradWi(p3i+1:p4i,k)+alpha*deltaWi(p3i+1:p4i,k-1);
    deltaWi(p2i+1:p3i,k)=etai(3,k)*gradWi(p2i+1:p3i,k)+alpha*deltaWi(p2i+1:p3i,k-1);
    deltaWi(p1i+1:p2i,k)=etai(2,k)*gradWi(p1i+1:p2i,k)+alpha*deltaWi(p1i+1:p2i,k-1);
    deltaWi(1:p1i,k)=etai(1,k)*gradWi(1:p1i,k)+alpha*deltaWi(1:p1i,k-1);
    
    %权值更新
    wi_out(:,k)=wi_out(:,k-1)+deltaWi(p4i+1:end,k);
    bi(:,k)=bi(:,k-1)+deltaWi(p3i+1:p4i,k);
    wi_r(:,k)=wi_r(:,k-1)+deltaWi(p2i+1:p3i,k);
    wi_in2(:,k)=wi_in2(:,k-1)+deltaWi(p1i+1:p2i,k);
    wi_in1(:,k)=wi_in1(:,k-1)+deltaWi(1:p1i,k);
end

%展示
kk=1:L;
P1=figure(1);
set(P1,'position',[200 400 800 600]);
plot(kk,y,'b--',kk,yr,'k-',kk,ym,'r:');
legend('实际输出','目标输出','识别输出');

P2=figure(2);
set(P2,'position',[200 0 800 600]);
plot(kk,ec,'b--',kk,em,'r:');
legend('控制误差','识别误差');

P3=figure(3);
set(P3,'position',[900 400 800 600]);
plot(yu);hold on;
yy=diff(y);uu=diff(u);
yyuu=yy./uu;
plot(yyuu);
legend('辨识器估计灵敏度','差分所得灵敏度')
title('控制对象灵敏度');

P4=figure(4);
set(P4,'position',[900 0 800 600]);
plot(u);
title('输入量');

P5=figure(5);
set(P5,'position',[1400 400 800 600]);
plot(etac(6,:));
title('学习率');
%数据和参数处理
clear p*


function y=Acfun(x)
    y=2./(1+exp(-x))-1;
end

function y=diffAcfun(x)
%福大那篇论文中激励函数的导数
    y=2*exp(-x)./(1+exp(-x)).^2;
end