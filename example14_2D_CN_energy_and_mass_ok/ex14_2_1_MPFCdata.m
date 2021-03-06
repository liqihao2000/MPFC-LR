function pde = ex14_2_1_MPFCdata(para)

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
end

pde = struct('epsilon',epsilon, ...
             'M',M, ...
             'alpha',alpha, ...
             'beta_bar',beta_bar, ...
             'beta',beta, ...
             'sigma',sigma, ...         
             'init',@init, ...
             'name','ex14_2_1_MPFCdata');

    function z = init(x,y) 
%         z = sin(2*pi*x/64).*cos(2*pi*y/64);
        z =   0.07 - 0.02*cos(2*pi*(x-12)/32).*sin(2*pi*(y-1)/32) ...
            + 0.02*cos(pi*(x+10)/32).^2.*cos(pi*(y+3)/32).^2 ...
            -0.01*sin(4*pi*x/32).^2.*sin(4*pi*(y-6)/32).^2;     
    end


end