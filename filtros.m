fc = 2*pi/3;                    %frequencia de corte
[b,a] = butter(6, fc, 's');     %lowpass Butterworth filter analogic.
fButter = tf(b,a);              %funcao de transferencia do filtro analógico

%a partir desse filtro, projetar um filtro digital
fprintf('Método da transformação bilinear');
                                %método da transformação bilinear
                                %it is a standard method of mapping the s or analog plane into the z or digital plane
fs = 1;                         %sampling frequency in hertz
[numd,dend] = bilinear(b, a, fs);
fButterDig = tf(numd, dend, fs); %funcao de transferencia do filtro digital na forma direta
[z,p,k] = tf2zp(numd, dend);    %polos, zeros e ganhos

%implementar para forma direta e cascata
%Direta
% Já foi gerada na forma direta na linha anterior

%Cascata
%sys = zpk(Z,p,k,Ts) creates a discrete-time zero-pole-gain model with sample time Ts (in seconds). 
%Set Ts = -1 or Ts = [] to leave the sample time unspecified. The input arguments Z, P, K are as in 
%the continuous-time case.
sys_6 = zpk(z(6),p(6),1);
sys_5 = zpk(z(5),p(5),1);
sys_4 = zpk(z(4),p(4),1);
sys_3 = zpk(z(3),p(3),1);
sys_2 = zpk(z(2),p(2),1);
sys_1 = zpk(z(1),p(1),1);

fButterDig_Cast = sys_6*sys_5*sys_4*sys_3*sys_2*sys_1*k;

[numCas, denCas, Ts] = tfdata(fButterDig_Cast, 'v');

[h,w] = freqz(numd, dend);
%gráfico do filtro modulo da resposta em frequencia forma direta
semilogy(w,abs(h));
set(gcf,'NumberTitle','off');
set(gcf,'Name','Modulo da resposta em frequencia precisão infinita - forma direta/cascata');
title('Precisão infinita - Forma Direta/Cascata')
xlabel('Frequência(\pi rad)')
ylabel('Modulo da resposta em frequência(dB)')
    
hold on
[hC,wC] = freqz(numCas, denCas);
%gráfico do filtro modulo da resposta em frequencia forma cascata
semilogy(wC,abs(hC));
grid on



%gráfico do filtro digital butterworth forma direta e cascata
dirCas = fvtool(numd,dend,numCas,denCas);
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','Forma direta e Forma em cascata');
legend(dirCas, 'Forma direta','Forma em cascata')

%Representação em ponto flutuante
[newnumd, newdend] = tfdata(fButterDig, 'v');           %'v' para gerar cada valor em uma célula
[zn,pn,kn] = tf2zp(numCas, denCas);                     %polos, zeros e ganhos

%direta
for i=4:-1:1                                            %de 4 casas decimais até 1
    newnumd = round(newnumd,i);                         %parametro = valor e numero de casas decimais
    newdend = round(newdend,i);
    
    [h1,w1] = freqz(newnumd, newdend);
    %gráfico do filtro modulo da resposta em frequencia forma direta
    semilogy(w1,abs(h1));
    set(gcf,'NumberTitle','off'); 
    set(gcf,'Name',['Modulo da resposta em frequencia precisão finita, ' num2str(i) ' casas decimais - forma direta/cascata']);
    title('Precisão finita - Forma Direta/Cascata')
    xlabel('Frequência(\pi rad)')
    ylabel('Modulo da resposta em frequência(dB)')
    hold on

%end

%cascata
%for i=4:-1:1
    
    for j=1:6
    z2(j) = round(zn(j),i);                             %parametro = valor e numero de casas decimais
    p2(j) = round(pn(j),i);   
    end
    
    sys_61 = zpk(z2(6),p2(6),1,-1);
    sys_51 = zpk(z2(5),p2(5),1,-1);
    sys_41 = zpk(z2(4),p2(4),1,-1);
    sys_31 = zpk(z2(3),p2(3),1,-1);
    sys_21 = zpk(z2(2),p2(2),1,-1);
    sys_11 = zpk(z2(1),p2(1),1,-1);

    fButterDig_Castnew = sys_61*sys_51*sys_41*sys_31*sys_21*sys_11*kn;
    [numCasn, denCasn] = tfdata(fButterDig_Castnew, 'v');


    [hC1,wC1] = freqz(numCasn, denCasn);
    %gráfico do filtro modulo da resposta em frequencia forma direta
    semilogy(wC1,abs(hC1));
    grid on
    
    newdirCas = fvtool(newnumd,newdend,numCasn,denCasn);
    set(gcf,'NumberTitle','off'); 
    set(gcf,'Name',['Forma Direta/Cascata: ' num2str(i) ' casas decimais']);    
    legend(newdirCas, 'Forma direta','Forma em cascata')
    
end

%Transformação em frequência na representação da forma direta tipo Z^-1 = -z^-1

[numPB,denPB] = tfdata(fButterDig,'v');
newnumPA = numPB;
newdenPA = denPB;
newnumPA(2) = -numPB(2);
newdenPA(2) = -denPB(2);
newnumPA(4) = -numPB(4);
newdenPA(4) = -denPB(4);
newnumPA(6) = -numPB(6);
newdenPA(6) = -denPB(6);

%fprintf('Função de transferência do filtro transformado em passa-alta')
%filtroFDPA = tf (newnumPA,newdenPA,1);

[hR1,wR1] = freqz(newnumPA, newdenPA);
%gráfico do filtro modulo da resposta em frequencia forma direta
%transformação em frequencia
semilogy(wR1,abs(hR1));
set(gcf,'NumberTitle','off'); 
set(gcf,'Name', 'Modulo da resposta em frequencia - Tranformação em passa alta forma direta');
title('Tranformação em passa alta forma direta')
xlabel('Frequência(\pi rad)')
ylabel('Modulo da resposta em frequência(dB)')

%Transformação em frequência na representação da forma direta tipo Z^-1 = -z^-2

[numPB,denPB] = tfdata (fButterDig,'v');
newnumPF = numPB;
newdenPF = denPB;
newnumPF(2) = 0;
newdenPF(2) = 0;
newnumPF(4) = 0;
newdenPF(4) = 0;
newnumPF(6) = 0;
newdenPF(6) = 0;
newnumPF(3) = -numPB(2);
newdenPF(3) = -denPB(2);
newnumPF(5) = numPB(3);
newdenPF(5) = denPB(3);
newnumPF(7) = -numPB(4);
newdenPF(7) = -denPB(4);
newnumPF(9) = numPB(5);
newdenPF(9) = denPB(5);
newnumPF(11) = -numPB(6);
newdenPF(11) = -denPB(6);
newnumPF(13) = numPB(7);
newdenPF(13) = denPB(7);

%fprintf('Função de transferência do filtro transformado passa-faixa');
%filtroFDPF = tf (newnumPF,newdenPF,1);

%hold on
[hR2,wR2] = freqz(newnumPF, newdenPF);
%gráfico do filtro modulo da resposta em frequencia forma direta
%transformação em frequencia
%%semilogy(wR2,abs(hR2));
%%set(gcf,'NumberTitle','off'); 
%%set(gcf,'Name', 'Modulo da resposta em frequencia - Tranformação em passa faixa forma direta');
%%title('Tranformação em passa faixa forma direta')
%%xlabel('Frequência(\pi rad)')
%%ylabel('Modulo da resposta em frequência(dB)')
grid on

fvtool(newnumPF,newdenPF);
fvtool(newnumPA,newdenPA);

%referencias
%https://www.passeidireto.com/arquivo/2444345/aula-09---projeto-de-filtros-iir
%https://www.mathworks.com/help/signal/ref/bilinear.html