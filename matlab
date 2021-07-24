if isempty(varargin)
        
    % simulacion de parametros
    
    par.a = 'exp'; % simulation name suffix: 'exp' experimental
    par.runId = 0; % simulation ID (used to reproduce results)
    par.R = 100; % numero de antenas de recepcion
    par.T = 100; % numero de antenas de transmision
  
par.mod = 'BPSK'; % tipo de modulacion 'BPSK','QPSK','16QAM','64QAM'
    par.nombre = ['ERR_' num2str(par.R) 'x' num2str(par.T) '_' par.mod '_' par.a] ;  % nombre de la simulacion (usado para guardar resultados)/ num2str( Convertir números en matriz de caracteres)
    par.muestras = 10000; % numero de muestras El enfoque de Monte Carlo implica la simulación repetida de muestras dentro de las funciones de densidad deprobabilidad de los datos de entrada
    par.SNR = 0:10:60; % rango SNR (dB) eje x
    par.detector ={'ZF-BLAST','MMSE-BLAST'}; % detector usado
     par.alg.maxiter = 5;
     else
      
    disp('usar los ajustes y parámetros de simulación personalizados...')    
    par = varargin{1}; % el único argumento es la estructura de paridad
  end
  rng(par.runId); % creacion de datos aleatorios
 
  % creacion de la constelacion de mapas de Gray
   switch (par.mod)
  case 'BPSK', 
      par.symbols = [ -1 1 ];
   end
 
par.Es = mean(abs(par.symbols).^2);% media de los simbolos
 par.Q = log2(length(par.symbols)); % numero de bits por simbolos
  par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');
  
  time_elapsed = 0;% tiempo de simulacion
  
  %Inicializar las matrices de resultados (detector x SNR)
  
  res.BER = zeros(length(par.detector),length(par.SNR)); % tasa de error de bit
  
  %generar un flujo de bits aleatorios (antena x bit x muestras)
  bits = randi([0 1],par.T,par.Q,par.muestras);
  
   tic
  for t=1:par.muestras
  
    % generacion de simbolos de tx
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
  
   % generar id matriz de canal Gaussiano y vector de ruido
   
    n = sqrt(0.5)*(randn(par.R,1)+1i*randn(par.R,1));
    H = sqrt(0.5)*(randn(par.R,par.T)+1i*randn(par.R,par.T));
    
     x = H*s; % transmision de un canal sin ruido
     
      % SNR bucle
    for k=1:length(par.SNR)
      
      N0 = par.T*par.Es*10^(-par.SNR(k)/10); %calcular la variación del ruido (el promedio de SNR por antena receptora es: SNR=MT*Es/N0)
      y = x+sqrt(N0)*n; %transmitir datos a través de un canal ruidoso
    
      % bucle de algoritmo      
      for d=1:length(par.detector)
          
          switch (par.detector{d}) % seleciona algoritmo
        case 'ZF-BLAST', % zero-forcing detection
            [idxhat,bithat] = ZF(par,H,y);
        case 'MMSE-BLAST', % MMSE detector
            [idxhat,bithat] = MMSE(par,H,y,N0);
        otherwise,
            error('par.detector no definido')           
          end
          
          %% metricas de error
          err = (idx~=idxhat);
            
        res.BER(d,k) = res.BER(d,k) + sum(sum(bits(:,:,t)~=bithat))/(par.T*par.Q);
         end % algorithm loop
                 
    end % SNR loop    
    end % trials loop
 
  %% normalizacion de resultados
  res.BER = res.BER/par.muestras;
  %% guardar resultados
  save([ par.nombre '_' num2str(par.runId) ],'par','res'); 
  %% mostrar los resultados
  ;
  figure(1)
  for d=1:length(par.detector)
    if d==1
    
    semilogy(par.SNR,res.BER(d,:),'LineWidth',2)
      hold on
    else
 
     semilogy(par.SNR,res.BER(d,:),'LineWidth',2)
    end
  end
  hold off
  grid on
  title('MODULACION BPSK')
  xlabel('SNR[dB]','FontSize',10)
  ylabel('(BER)','FontSize',10)
  axis([min(par.SNR) max(par.SNR) 1e-3 1])
  legend(par.detector,'FontSize',10)
  set(gca,'FontSize',10)
  
end
 
 %% zero-forcing (ZF-BLAST) detector
function [idxhat,bithat] = ZF(par,H,y)
  xhat = H\y; 
  
  % xhat =(H'*H)\eye(H'*y);
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.T,1)*par.symbols).^2,[],2);% disminucion de numero de señales de interferencia
  bithat = par.bits(idxhat,:);
end 
        
 %% MMSE detector (MMSE-BLAST)
function [idxhat,bithat] = MMSE(par,H,y,N0)
  xhat = (H'*H+(N0/par.Es)*eye(par.T))\(H'*y);    
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.T,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);  
end
