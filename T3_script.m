clear all;
clc;
p=6;
Rs=0.294; % Resistência do enrolamento do estator
Lls=1.39e-3; %Indutância de dispersão do estator
Lm=41e-3; %Indutância de magnetização
Ls=Lm+Lls;%Indutância equivalente do estator
Rr=0.156; %Resistência do enrolamento do rotor
Llr=0.740e-3; %Indutância de dispersão do rotor
Lr=Lm+Llr; %Indutância equivalente do rotor
JM=0.4; %Momento de inércia do motor
ys=0+j*0; %Fluxo concatenado do estator
yr=0+j*0; %fluxo concatenado do rotor
wm=0; %Velocidade inicial.
Q=0; 
VLL=220;
Vs=VLL/sqrt(3);
Vm=Vs*sqrt(2);
f=60;
T = 1/f;
w=2*pi*f;
passo=1e-5;
tfinal=2;
c=0;
TL = 0;   
integral_Iefa = 0;
integral_Iefb = 0;
integral_Iefc = 0;
k=sqrt(2/3);
for t=0:passo:tfinal
    c=c+1;

    %Definição das tensões no referencial abc
    vas=Vm*cos(w*t); 
    vbs=Vm*cos(w*t-2*pi/3);
    vcs=Vm*cos(w*t+2*pi/3);
    
    %Declaração das constantes de deslocamento angular.
    alfa=exp(j*2*pi/3);
    alfa2=exp(j*4*pi/3);

    %Declaração do vetor resultante de tensões no estator.
    vs=k*(vas+alfa*vbs+alfa2*vcs);
    %Declaração das componentes no referencial alfa-beta
    valfas=real(vs);
    vbetas=imag(vs);
    
    %Vetor de fluxo concatenado, fluxo estator-fluxo rotor.
    Y=[ys;yr];
    L=[Lls+Lm Lm;Lm Llr+Lm];
    %Vetor de correntes no estator.
    I=inv(L)*Y;
    %Vetor das correntes resultantes no estator
    is=I(1);
    ialfas=real(is);
    ibetas=imag(is);
    ir=I(2);
    
    %Torque do motor 
    T=(1/(k^2))*(2/3)*(p/2)*imag(yr*ir');
    %Velocidade angular elétrica
    w0=(p/2)*wm;
    
    %Rotação do motor;
    nm=wm*60/(2*pi);
    %Corrente no referêncial abc;
    ias=(1/k)*(2/3)*real(is);
    ibs=(1/k)*(2/3)*real(alfa2*is);
    ics=(1/k)*(2/3)*real(alfa*is);

    %Potência elétrica abc
    pele=vas*ias+vbs*ibs+vcs*ics;
    

    %Potência elétrica no alfa-beta
    palfabeta=valfas*ialfas+vbetas*ibetas;
    
    %palphabeta 

    pelealfabeta2 = 1/k^2*2/3*real(vs*is);
    %Definição da carga JM = JL;
    JL=0.4;
    
    %Momento de inércial total
    J=JM+JL;
    

    %Definição da carga aplicada 

    
    if t > 0.05
        TL = 7460/2/(2*pi*1164/60);
    else
        TL = 0;
    end


    %------ Aplicação do método de euler-------
    %Derivada do vetor fluxo concatenado no estator. 
    pys=vs-Rs*is;
    
    %Derivada do fluxo concatenado no rotor.
    pyr=j*w0*yr-Rr*ir;
    

    %Derivada da velocidade angular (aceleração angular)
    pwm=(T-TL)/J;
    
    %Velocidade 
    pQ=wm;
    
    % x(n) = x(n-1) + dx/dt*deltaT 
    ys=ys+pys*passo;
    yr=yr+pyr*passo;
    wm=wm+pwm*passo;
    Q=Q+pQ*passo;

   
     %Plotagem dos vetores.
    v_t(c)=c;
    v_T(c)=T;
    v_wm(c)=wm;
    v_TL(c)=TL;
    v_Q(c)=Q;
    v_ias(c)=ias;
    v_ibs(c)=ibs;
    v_ics(c)=ics;
    v_ialfas(c)=ialfas;
    v_ibetas(c)=ibetas;
    v_valfas(c)=valfas;
    v_vbetas(c)=vbetas;
    v_pele(c)=pele;
    v_palfabeta(c)=palfabeta;

    v_pelealfabeta2(c)=pelealfabeta2;

end

 

figure(1)
plot(v_t,v_wm*60/(2*pi)),grid
ylabel('n_m (rpm)'),xlabel('t (s)')


figure(2)
plot(v_t,v_T,v_t,v_TL),grid
ylabel('T e T_L (N.m)'),xlabel('t (s)')
legend('T','T_L')

figure(3)
plot(v_t,v_ias,v_t,v_ibs,v_t,v_ics),grid
xlabel('t (s)'),ylabel('i_a, i_b, i_c (A)')
legend('i_a','i_b','i_c')

% plotar=0;
% if plotar==1
    figure(4)
    plot(v_t,v_ialfas,v_t,v_ibetas),grid
    ylabel('I'),xlabel('t(s)')
    legend('I alfaS','IbetaS');

   figure(5)
   plot(v_t,v_pele)
   ylabel('Pot.Elétrica'),xlabel('t(s)'),grid
    legend('P_ele');

   figure(6)
   plot(v_t,v_palfabeta),grid
   ylabel('Pot.alphaBeta'),xlabel('t(s)')
   legend('Potência AlphaBeta');


   figure(7)
   plot(v_t,v_pelealfabeta2),grid
   ylabel('Pot.alphaBeta2'),xlabel('t(s)')
   legend('Potência AlphaBeta2');
 

  