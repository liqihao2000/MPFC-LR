function pde = ex19_2_MPFCdata(para)

if nargin == 0
    epsilon = 1;
    M = 1;
    alpha = 1;    
    beta_bar    = 1;
    beta  = 1;
    sigma = 0;    
else
    if ~isstruct(para)
        exit('we need a struct data');
    end
    if ~isfield(para,'epsilon') || isempty(para.epsilon)
        epsilon = 1;
    else
        epsilon = para.epsilon;
    end
    if ~isfield(para,'M') || isempty(para.M)
        M = 1;
    else
        M = para.M;
    end    
    if ~isfield(para,'alpha') || isempty(para.alpha)
        alpha = 1;
    else
        alpha = para.alpha;
    end
    if ~isfield(para,'beta') || isempty(para.beta)
        beta = 1;
    else
        beta = para.beta;
    end  
    if ~isfield(para,'beta_bar') || isempty(para.beta_bar)
        beta_bar = 1;
    else
        beta_bar = para.beta_bar;
    end     
    if ~isfield(para,'sigma') || isempty(para.sigma)
        sigma = 0;
    else
        sigma = para.sigma;
    end     
    if ~isfield(para,'name') || isempty(para.name)
        name = 'ex19_2_MPFCdata';
    else
        name = para.name;
    end     
end

pde = struct('epsilon',epsilon, ...
             'M',M, ...
             'alpha',alpha, ...
             'beta_bar',beta_bar, ...
             'beta',beta, ...
             'sigma',sigma, ... 
             'init',@init, ...
             'name',name);

    function z = init(x,y,z) 
%         %         
%         n   = size(x(:),1);
%         N   = round(n^(1/3)); 
%         z = -0.35 +0.01*randi([-1 1],N,N,N);        
%         save('ex19_2_MPFCdata_Init_u0.mat','z');
        
        load('ex19_2_MPFCdata_Init_u0.mat');
    end
end