function pde = ex18_2_MPFCdata(para)

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
        name = 'ex18_2_MPFCdata';
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
        phi_0 = 0.12;
        A_cof = 0.15;
        a = 0.66;
        b = 0.38;
        c = 0.46;
        
        theta = pi/6;
        x2 = x;
        y2 =  cos(theta)*y - sin(theta)*z;
        z2 =  sin(theta)*y + cos(theta)*z;
        
        theta = -pi/6;
        x3 = x;
        y3 =  cos(theta)*y - sin(theta)*z;
        z3 =  sin(theta)*y + cos(theta)*z;       

        r2 = cos(a*x2).*cos(b*y2)+cos(a*x2).*cos(c.*z2)+cos(b.*y2).*cos(c.*z2)-0.5*cos(b.*y2);
        r3 = cos(a*x3).*cos(b*y3)+cos(a*x3).*cos(c.*z3)+cos(b.*y3).*cos(c.*z3)-0.5*cos(b.*y3);

        d = 15;
        x0 =  45; y0 =  45; z0 = 64;
        w2 = (1-((x-x0).^2+(y-y0).^2+(z-z0).^2)/d^2).^2;
        w2 = w2.*((x-x0).^2+(y-y0).^2+(z-z0).^2<=d^2);
        
        x0 = 75; y0 = 75; z0 = 64;
        w3 = (1-((x-x0).^2+(y-y0).^2+(z-z0).^2)/d^2).^2;
        w3 = w3.*((x-x0).^2+(y-y0).^2+(z-z0).^2<=d^2);

        z=phi_0 + w2.*(A_cof*r2) + w3.*(A_cof*r3);
%         save('ex18_1_2_MPFCdata_Init_u0.mat','z');
%         load('ex18_1_2_MPFCdata_Init_u0.mat');
    end
end