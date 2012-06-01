function LTAS
  
  x=wavread('/home/guido/Escritorio/guido.wav',10000);
   
  
  //ejemplo con senoidales conocidas
//  t=[0:1/10000:1-1/10000];
//  x=5*sin(2*%pi*50*t) + 10*sin(2*%pi*120*t) + 2*sin(2*%pi*30*t);

  L=100; //cantidad de ventanas
  N=length(x);
  //ventana de hamming
  //definimos el tamanio de la ventana
   
  TamVent = round(N/L);
  h=window('hm',TamVent);
  
  Transf=zeros(L,N);
  ltas=zeros(1,N);
  
  for i=1:L
      xaux = x(((i-1)*(TamVent) +1):i*TamVent);//tomo valores de x en la ventana
      xaux = [zeros(1,(i-1)*(TamVent)) xaux zeros(1,N-(i*TamVent))];
      haux = [zeros(1,(i-1)*(TamVent)) h zeros(1,N-(i*TamVent))];
      Transf(i,:)=((fft(xaux.*haux)).^2 )./TamVent;
  end
  
  //Transf es i corresponde a la ventana y j corresponde a los valores de la trasnformada
  
  for j=1:N
      ltas(j)=sum(Transf(1:L,j))/L;
      ltas(j)=10*log(ltas(j));    
  end
  figure(1)
  subplot(2,1,1)
  plot2d3(abs(ltas(1:N)));
  subplot(2,1,2)
  //figure(2)
  plot2d3(abs(fft(x)));
endfunction