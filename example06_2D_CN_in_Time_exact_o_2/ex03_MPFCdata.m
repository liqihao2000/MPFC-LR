function pde = ex03_MPFCdata(para)

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
        sigma = 1;
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
             'exact',@exact, ...
             'exactpsi',@exactpsi, ...
             'init',@init, ...
             'rhs',@rhs, ...
             'name','ex03_MPFCdata');

    function z = init(x,y) 
        z = sin(2*pi*x/64).*cos(2*pi*y/64);
    end

    function z = exact(x,y,t)  
        z = sin(2*pi*x/64).*cos(2*pi*y/64).*cos(t);
    end

    function z = exactpsi(x,y,t)  
        z = -sin(2*pi*x/64).*cos(2*pi*y/64).*sin(t);
    end

    function z = rhs(x,y,t)  
        z= -(cos((y.*pi)./32).*sin((x.*pi)./32).*(134217728.*beta_bar.*cos(t) + 134217728.*beta.*sin(t) - 262144.*M.*pi.^2.*cos(t) + 1024.*M.*pi.^4.*cos(t) - M.*pi.^6.*cos(t) + 786432.*M.*pi.^2.*cos(t).^3 - 134217728.*M.*sigma.*cos(t) - 786432.*M.*pi.^2.*cos((x.*pi)./32).^2.*cos(t).^3 - 1572864.*M.*pi.^2.*cos((y.*pi)./32).^2.*cos(t).^3 + 262144.*M.*epsilon.*pi.^2.*cos(t) + 2359296.*M.*pi.^2.*cos((x.*pi)./32).^2.*cos((y.*pi)./32).^2.*cos(t).^3))./134217728;
    end

end