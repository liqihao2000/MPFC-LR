function pde = ex21_2_MPFCdata(para)

if nargin == 0
    epsilon = 1;
    M = 1;
    alpha = 1;    
    beta_bar = 1;
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
        name = 'ex21_MPFCdata';
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

    function z = init(x,y)
        eta = 0;
        s   = 0;
        
        d = 1/20*(max(x(:)) - min(x(:)));
        kx = sqrt(3/4*(1+eta^2*s^2)/(1+(3/2)*eta^2*s^2));
        ky = sqrt(3/4*(1+eta)^2/(1+(3/2)*eta^2*s^2));
%         phi_bar = 0.285;
        phi_bar = 0.14;
        x0 = (max(x(:)) + min(x(:)))/2;
        y0 = (max(y(:)) + min(y(:)))/2;        
        w  = (1-((x-x0).^2+(y-y0).^2)/d^2).^2;
        w = w.*((x-x0).^2+(y-y0).^2<=d^2);
        R = (1/2)*eta^2*s^2/(1+(3/2)*eta^2*s^2);
        phis = cos(kx*x).*cos(ky*(y-s*x)/sqrt(3))+1/2*cos(2*ky*(y-s*x)/sqrt(3));
        A_cof = -4/5*(phi_bar+sqrt(15*(epsilon+R)-36*phi_bar.^2)/3);
        z = phi_bar+w.*(A_cof*phis);
    end

end