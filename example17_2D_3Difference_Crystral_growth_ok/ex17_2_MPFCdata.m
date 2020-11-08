function pde = ex17_2_MPFCdata(para)

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
        name = 'ex17_2_MPFCdata';
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
        phi_0 = 0.14;
        d = 20;
        
        x0 = 350; y0=400;
        w1  = (1-((x-x0).^2+(y-y0).^2)/d^2).^2;
        w1 = w1.*((x-x0).^2+(y-y0).^2<=d^2);
        
        x0 = 200; y0=200;
        w2  = (1-((x-x0).^2+(y-y0).^2)/d^2).^2;
        w2 = w2.*((x-x0).^2+(y-y0).^2<=d^2);
        
        x0 = 600; y0=300;
        w3  = (1-((x-x0).^2+(y-y0).^2)/d^2).^2;
        w3 = w3.*((x-x0).^2+(y-y0).^2<=d^2);
        
        theta = pi/4;
        x2 =  cos(theta)*x + sin(theta)*y;
        y2 = -sin(theta)*x + cos(theta)*y;
        
        theta = -pi/4;
        x3 =  cos(theta)*x + sin(theta)*y;
        y3 = -sin(theta)*x + cos(theta)*y;
        
%         z1 = C*(cos(q*y/sqrt(3)).*cos(q*x) - 0.5*cos(2*q*y/sqrt(3)) );
%         z2 = C*(cos(q*y2/sqrt(3)).*cos(q*x2) - 0.5*cos(2*q*y2/sqrt(3)) );
%         z3 = C*(cos(q*y3/sqrt(3)).*cos(q*x3) - 0.5*cos(2*q*y3/sqrt(3)) );
        
        eta = 0;
        s   = 0;
        
        kx = sqrt(3/4*(1+eta^2*s^2)/(1+(3/2)*eta^2*s^2));
        ky = sqrt(3/4*(1+eta)^2/(1+(3/2)*eta^2*s^2));
        phi_bar = 0.285;       
        R = (1/2)*eta^2*s^2/(1+(3/2)*eta^2*s^2);        
        A_cof = -4/5*(phi_bar+sqrt(15*(epsilon+R)-36*phi_bar.^2)/3);        
        z1 = cos(kx*x).*cos(ky*(y-s*x)/sqrt(3))+1/2*cos(2*ky*(y-s*x)/sqrt(3));
        z2 = cos(kx*x2).*cos(ky*(y2-s*x2)/sqrt(3))+1/2*cos(2*ky*(y2-s*x2)/sqrt(3));
        z3 = cos(kx*x3).*cos(ky*(y3-s*x3)/sqrt(3))+1/2*cos(2*ky*(y3-s*x3)/sqrt(3));  
        
        z=phi_0 + w1.*(A_cof*z1) + w2.*(A_cof*z2)+ w3.*(A_cof*z3);
%         save('ex17_1_MPFCdata_Init_u0.mat','z');
%         load('ex17_1_MPFCdata_Init_u0.mat');
    end
end