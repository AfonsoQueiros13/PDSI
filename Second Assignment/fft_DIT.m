% AUTOR: Prof. Anï¿½bal Joï¿½o de Sousa Ferreira, FEUP/INESC Tec
% DATA : 27/04/2019
% Este programa ï¿½ parte integrante dos trabalhos prï¿½ticos
% da UC de Processamento Digital de Sinal do
% 3ï¿½ ano do MIEEC (FEUP)
%
% ï¿½ Anï¿½bal Ferreira
%
% algoritmo FFT "Decimation in Time", radix 2
%
% NOTA: este cï¿½digo ï¿½ meramente ilustrativo de uma soluï¿½ï¿½o
%       acadï¿½mica de implementaï¿½ï¿½o da FFT-DIT.
%
% Hï¿½, entre outras, duas formas simples de testar o algoritmo da FFT.
%
% Uma consiste em fornecer ï¿½ entrada da FFT, um vector de dados que
% corresponda a uma sinusï¿½ide (ou cossinusï¿½ide) cuja frequï¿½ncia ï¿½
% exatamente amostrada pela DFT. Neste caso, no intervalo zero-PI,
% existirï¿½ uma ï¿½nica linha espectral diferente de zero.
%
% Uma outra consiste em apresentar uma rampa ï¿½ entrada da FFT,
% normalizada entre zero e a unidade. Neste caso, a parte real
% de todas as linhas espectrais, exceto a primeira, valerï¿½ exactamente -1/2.
%
% Por ï¿½ltimo, pode-se naturalmente confirmar atravï¿½s da funï¿½ï¿½o Matlab fft()
% se se obtï¿½m os mesmos resultados ...


% INICIALIZAï¿½ï¿½O
% Nï¿½mero de estï¿½gios. Este parï¿½metro define todos os outros na implementaï¿½ï¿½o da FFT.
est=4;
N=2^est

% vector de teste para o caso de uma rampa de magnitude normalizada.
%re=1/N*([0:N-1]);
%im=zeros(1,N);
%Z = [re im];
%Y = fft(Z);


% vector de teste para uma sinusï¿½ide pura ï¿½ entrada da FFT.
% (comentar se se preferir o sinal de rampa)
%bin_teste=2  % 0 <= bin_teste <= N/2 = k
re=sin(bin_teste*2*pi/N*[0:N-1]); im=zeros(1,N);

xnovo=re+j*im;
Z = [re im];
Y = fft(Z);

%
% preparaï¿½ï¿½o do vector que define as trocas do "bit reversal"
% (nï¿½o ï¿½ necessï¿½rio analisar em detalhe, nem preciso responder ï¿½s 3 perguntas)
%
nchg=0; N2=N/2; flag=zeros(1,N);
for i=1:N-2,                         % excluem-se os casos de i==0 e i==N-1, porquï¿½ ?
   if(flag(i)==0)
      dv=N2; ic=1; final=0; pos=i;
      for k=1:est,                   % que faz este ciclo "for" ?
         if (pos/dv >= 1)
            pos=pos-dv;
            final=final+ic;
         end
         dv=dv/2; ic=ic*2;
      end
      if (i~=final)                  % que se pretende evitar aqui ?
         rev(1+nchg)=i;
         rev(N2+nchg)=final;
         flag(final)=1;
         nchg=nchg+1;
      end
   end   
end

% este ï¿½ o nï¿½mero de trocas a fazer
nchg

% este ï¿½ o vector que detalha as trocas a fazer
rev(1:nchg)        % pontos de um lado
rev(N2:N2+nchg-1)  % pontos do outro lado ("bit reversed")

%
% "bit reversal"
%
for i=1:nchg,
   tmp=re(1+rev(i)); re(1+rev(i))=re(1+rev(N2+i-1)); re(1+rev(N2+i-1))=tmp;
   tmp=im(1+rev(i)); im(1+rev(i))=im(1+rev(N2+i-1)); im(1+rev(N2+i-1))=tmp;
end

%re   % antes da FFT
%im   % antes da FFT

grup=N; butf=1; phi=2*pi/N;
for i=0:est-1,                            % avanï¿½a nos estï¿½gios
   grup=grup/2;
   for j=0:grup-1,                        % avanï¿½a nos grupos
      for k=0:butf-1,                     % avanï¿½a nas borboletas
         ptr1=2*j*butf+k;
         ptr2=ptr1+butf;
         arg=phi*grup*k;
         C=cos(arg); S=sin(arg);
         retmp=C*re(1+ptr2)+S*im(1+ptr2);
         imtmp=C*im(1+ptr2)-S*re(1+ptr2);
         re(1+ptr2)=re(1+ptr1)-retmp;
         im(1+ptr2)=im(1+ptr1)-imtmp;
         re(1+ptr1)=re(1+ptr1)+retmp;
         im(1+ptr1)=im(1+ptr1)+imtmp;
      end
   end
   butf=butf*2;
end

%re  % depois da FFT
%im  % depois da FFT

y=(re.^2+im.^2)/N^2;


plot([0:N-1],y);
title('Densidade Espectral');
xlabel('Linha Espectral');
ylabel('Quadrado da Magnitude');
figure(2);
plot(Y);
title('Densidade Espectral usando a função fft()');
xlabel('Linha Espectral');
ylabel('Quadrado da Magnitude');

