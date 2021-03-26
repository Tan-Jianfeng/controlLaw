%利用PID来控制低压转速回路
clearvars;clc;
A=[-2.01 2.06;-0.0044 -2.92];
B=[0.84 1.2;0.53 0.29];
C=[1 0;-0.86 5.08];
D=[0 0;-0.21 0.56];
[b,a]=ss2tf(A,B,C,D,1);
clearvars -regexp [A-Z]
engine=tf(b(2),a);

stepplot(engine)
C = pidtune(engine,'PI');
pidTuner(engine,C)