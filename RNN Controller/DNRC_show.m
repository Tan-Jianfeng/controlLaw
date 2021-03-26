%画图，展示控制器效果
%%
%输出
Ts=0.04;
L=600;
kk=1:L;
tt=kk*Ts;
P1=figure(1);
plot(tt,y(kk),'b--*',tt,yr(kk),'k-',tt,ym(kk),'r:o');
legend('实际输出','目标输出','识别输出');
xlabel('{\it{t}}/s');ylabel('{\it{n}}_L')
ylim([-0.1 1.2])
set(P1,'position',[200 400 800 600]);
set(gca,'fontsize',10.5);
%%
P2=figure(2);
set(P2,'position',[200 0 800 600]);
plot(tt,ec(kk),'b--',tt,em(kk),'r:','linewidth',2);
legend('控制误差','识别误差');
xlabel('{\it{t}}/s');ylabel('误差')
set(P2,'position',[200 400 800 600]);
set(gca,'fontsize',10.5);
%%
P3=figure(3);
set(P3,'position',[900 0 800 600]);
plot(tt,u(kk),'k-p');
xlabel('{\it{t}}/s');ylabel('{\it{w}}_f/kg/s')
set(P3,'position',[200 400 800 600]);
set(gca,'fontsize',10.5);
%%
P4=figure(4);
set(P4,'position',[1400 400 800 600]);
plot(tt,etac(6,kk),'k-h');
xlabel('{\it{t}}/s');ylabel('\eta')
set(P4,'position',[200 400 800 600]);
set(gca,'fontsize',10.5);
%%
%画模糊推理模型
LR=readfis('LR_engine');
figure(5)
plotmf(LR,'input',1)
xlabel('e_c(k)');ylabel('隶属度');
set(gca,'fontsize',10.5);

figure(6)
plotmf(LR,'input',2)
xlabel('\Delta W(k)');ylabel('隶属度');
set(gca,'fontsize',10.5);

figure(7)
gensurf(LR);
xlabel('e_c(k)');ylabel('\Delta W(k)');zlabel('\eta (k)');
set(gca,'fontsize',10.5);