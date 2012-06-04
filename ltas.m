function ltas
  x=wavread('/home/guido/Escritorio/1-sobrio.wav');
  x2=wavread('/home/guido/Escritorio/1-ebrio.wav');
  
  %sobrio
  L=10;
  N=length(x);
  TamVent = round(N/L);
  
  h=hanning(TamVent);
  Transf=zeros(L,N);
  ltas=zeros(1,N);
  for i=1:L
      xaux = x(((i-1)*(TamVent) +1):i*TamVent);
      xaux = [zeros(1,(i-1)*(TamVent)) xaux zeros(1,N-(i*TamVent))];
      haux = [zeros(1,(i-1)*(TamVent)) h'        zeros(1,N-(i*TamVent))];
      Transf(i,:)=((fft(xaux(1:N).*haux(1:N))).^2 )./TamVent;
  end
  
  
  for j=1:N
      ltas(j)=sum(Transf(1:L,j))/L;
      ltas(j)=10*log(ltas(j));    
  end
  
  figure(1)
  subplot(2,1,1)
  stem(abs(ltas(1:N)));
  subplot(2,1,2)
  stem(abs(fft(x)));

 
 %ebrio
  N2=length(x2);
  TamVent2 = round(N2/L);
  
  h2=hanning(TamVent2);
  Transf2=zeros(L,N2);
  ltas2=zeros(1,N2);
  for i=1:L
      xaux2 = x2(((i-1)*(TamVent2) +1):i*TamVent2);
      xaux2 = [zeros(1,(i-1)*(TamVent2)) xaux2 zeros(1,N2-(i*TamVent2))];
      haux2 = [zeros(1,(i-1)*(TamVent2)) h2'        zeros(1,N2-(i*TamVent2))];
      Transf2(i,:)=((fft(xaux2(1:N2).*haux2(1:N2))).^2 )./TamVent2;
  end
  
  
  for j=1:N2
      ltas2(j)=sum(Transf2(1:L,j))/L;
      ltas2(j)=10*log(ltas2(j));    
  end
 
  figure(2)
  subplot(2,1,1)
  stem(abs(ltas2(1:N2)));
  subplot(2,1,2)
  stem(abs(fft(x2)));