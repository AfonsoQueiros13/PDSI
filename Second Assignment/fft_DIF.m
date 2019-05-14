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
bin_teste=2  % 0 <= bin_teste <= N/2 = k
re=sin(bin_teste*2*pi/N*[0:N-1]); im=zeros(1,N);
xnovo=re+j*im;
%Z = [re im];
%Y = fft(Z);

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

%re   % antes da FFT
%im   % antes da FFT

S=log2(N);                                                                  % computing the number of conversion stages
Half=N/2;                                                                   % half the length of the array
for stage=1:S;                                                              % stages of transformation
    for index=0:(N/(2^(stage-1))):(N-1);                                    % series of "butterflies" for each stage
        for n=0:(Half-1);                                                   % creating "butterfly" and saving the results
            pos=n+index+1;                                                  % index of the data sample
            pow=(2^(stage-1))*n;                                            % part of power of the complex multiplier
            w=exp((-1i)*(2*pi)*pow/N);                                      % complex multiplier
            a=re(pos)+re(pos+Half);                                           % 1-st part of the "butterfly" creating operation
            b=(re(pos)-re(pos+Half)).*w;                                      % 2-nd part of the "butterfly" creating operation
            re(pos)=a;                                                       % saving computation of the 1-st part
            re(pos+Half)=b;                                                  % saving computation of the 2-nd part
            c=im(pos)+im(pos+Half);                                           % 1-st part of the "butterfly" creating operation
            d=(im(pos)-im(pos+Half)).*w;                                      % 2-nd part of the "butterfly" creating operation
            im(pos)=c;                                                       % saving computation of the 1-st part
            im(pos+Half)=d;                                                  % saving computation of the 2-nd part   
        end;
    end;
    Half=Half/2;                                                                % computing the next "Half" value
end;
%re  % depois da FFT
%im  % depois da FFT
y=(re.^2+im.^2)/N^2;
z=bitrevorder(abs(y));     
plot([0:N-1],z);
title('Densidade Espectral');
xlabel('Linha Espectral');
ylabel('Quadrado da Magnitude');
%figure(2);
%plot(Y);
%title('Densidade Espectral usando a função fft()');
%xlabel('Linha Espectral');
%ylabel('Quadrado da Magnitude');





                                        % entry point to the function of FFT decimation in frequency


