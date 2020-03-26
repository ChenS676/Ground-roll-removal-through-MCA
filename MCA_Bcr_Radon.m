function parts=MCA_Bcr_Radon(signal,wave1,wave2,Cweight,thdtype,itermax,expdecrease,stop,sigma,display,dt)
% MCA_Bcr: Morphological Component Analysis of a 1D signal (a vector) using highly redundant dictionaries.
%	   The optimization bp is solved using a modified version of the BCR algorithm.
%	   MCA_bcr solves the following optimization problem
%		(part_i,for all i) = argmin Sum_i || \Phi^+_i part_i ||_p + lambda * ||img - Sum_i part_i||_2^2
%			Sparsity is promoted by the lp norm. Ideally the l_0 norm but this is relaxed to the l_1 norm.	
%				p = 1 (l_1 norm: Soft thesholding as a solution).
%				p = 0 (l_0 norm: difficult but approximated with a Hard thresholding).
%	   Each component is supposed to be sparsely described in its corresponding dictionary \Phi.
%  Usage:
%    part=MCA_Bcr(signal,dict,pars1,pars2,pars3,itermax,gamma,comptv,expdecrease,stop,mask)
%  Inputs:
%    signal     	1D signal vector, n = 2^J, if n ~= 2^J, use zero-padding and crop at the end.
%    dict		Names of dictionaries for each part (see directory Dictionary for available transforms)
%    parsi		Parameters of dictionary (built using the MakeList function)
%    itermax		Nb of relaxation iterations
%    expdecrease	Exponential/Linear decrease of the regularization parameter
%    stop	 	Stop criterion, the algorithm stops when lambda <= stop*sigma (typically k=3), sigma is the noise WGN std
%    sigma		Value of noise std. If not provided, it will be estimated (default).
%    display		Display algorithm progress. 
%  Outputs:
%    parti 		Estimated ith semantic component (1D signal).
%    options		Structure containing all options and parameters of the algorithm, including default choices.
%
%  Description
%    The dictionaries and their parameters can be built using the MakeList function.
%    A demo GUI (MCADemo) can be called to guide the user in these steps. 
%  See Also
%    FastLA, FastLS, MCADemo


% initializations. Put the signal as a column vector.
numberofdicts = 2;              % no. of signal components
sigpad=signal;
[nt,nx]=size(sigpad);
n=nt*nx;
part	      = zeros(n,numberofdicts);
CQ1           = Cweight.CQ1;    % l2 norm normalization constants
CQ2           = Cweight.CQ2;    

% dictionaries - wavelet transform
p1=wave1.p1;
x1=wave1.x1;

% dictionaries - chirplet transform
p2=wave2.p2;
x2=wave2.x2;

np=length(p1);
stopcriterion = stop*sigma;
t=(0:nt-1).*dt;
w1coef=cgnr_radon_hyper(sigpad,dt,x1,p1);
figure
pcolor(p1,t,w1coef),shading interp;set(gca,'XAxisLocation','top');axis ij;
set(gcf, 'Renderer', 'ZBuffer');
xlabel('扫描斜率');ylabel('时间/s');
w2coef=cgnr_radon_linear(sigpad,dt,x2,p2);
figure
pcolor(p2,t,w2coef),shading interp;set(gca,'XAxisLocation','top');axis ij;
set(gcf, 'Renderer', 'ZBuffer');
xlabel('扫描斜率');ylabel('时间/s');
for i=1:(np)
     c(i)=max(abs(w1coef(:,i)));        
 end
maxW1=max(c)
for i=1:(np)
     d(i)=max(abs(w2coef(:,i)));            
 end
maxW2=max(d)
deltamax = min([maxW1 maxW2])


% calculate the starting thd, which is the minimum of maximal coefficients 
% of the signal in each dictionary.
delta=deltamax ;
if expdecrease	
   lambda=(deltamax/stopcriterion)^(1/(1-itermax)); % Exponential decrease.
else	 	
   lambda=(deltamax-stopcriterion)/(itermax-1);	 % Slope of the linear decrease. 
end
	
if display
   % Create and return a handle on the waitbar.
   h = waitbar(0,'MCA in progress: Please wait...');
   nbpr=ceil(sqrt(numberofdicts+2));
end

part = zeros(n,2); 
% start the modified Block Relaxation Algorithm.
for iter=0:itermax-1
    tic;
    iter
   % calculate the residual signal.
   summ=reshape(sum(part,2),nt,nx);
   residual1=sigpad-summ;
        
   % cycle over dictionaries
   % Update Parta assuming other parts fixed.
   % Solve for Parta the marginal penalized minimization problem (Hard thesholding, l_1 -> Soft).
  residual=reshape(residual1,n,1) ;
   Ra1    = part(:,1)+residual;
   Ra=reshape(Ra1,nt,nx); 
   w1coef=cgnr_radon_hyper(Ra,dt,x1,p1);
   w1=w1coef;
    for i=1:(np)
   w1coef(:,i) = eval([thdtype 'Thresh(w1coef(:,i),CQ1*delta)']); 
    end    
   a= invfwd_tx_sstackn_hyper(w1coef,dt,p1,x1);
   part(:,1)=reshape(a,n,1);
   Ra2    = part(:,2)+residual;
   Ra=reshape(Ra2,nt,nx); 
   w2coef=cgnr_radon_linear(Ra,dt,x2,p2);
   w2=w2coef;
    for i=1:(np)
   w2coef(:,i) = eval([thdtype 'Thresh(w2coef(:,i),CQ2*delta)']); 
    end
    b=invfwd_tx_sstackn_linear(w2coef,dt,p2,x2);
    part(:,2)=reshape(b,n,1);
   % Update the regularization parameter delta.
   if expdecrease	delta=delta*lambda; % Exponential decrease.	
   else		delta=delta-lambda; % Linear decrease.
   end
	% Displays the progress.
	if display,
	   waitbar((iter+1)/itermax,h);	    
	   figure(display);
	   subplot(nbpr,nbpr,1);plot(signal);axis tight;drawnow;
	   subplot(nbpr,nbpr,2);plot(sum(part(1:N,:),2));title('\Sigma_i Part_i');axis tight;drawnow;
	   for np=1:numberofdicts
	      subplot(nbpr,nbpr,np+2);plot(part(1:N,np));title(sprintf('Part_%d',np));axis tight;drawnow;
	   end
    end
    toc;
end
	
if display,
   % Closes the waitbar window
   close(h);
end

% Crop data to original size.
parts = part(1:n,:);
