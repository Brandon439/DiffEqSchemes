function u = conjugate_gradient(A,f,tol)
%
%   Example:  
%   x = conjugate_gradient(A,f,tol) 
%   
MAXITS = length(f);

u = 0*f;
%inserted what r is
r = f - A*u;
p = r;
for k = 1:MAXITS
    w = A*p;
    alpha = (r'*r)/(p'*w);
    %inserted line iterating u
    u = u + alpha*p;
    rnew = r - alpha*w;
    if( norm(rnew) < tol )
        fprintf('Converged! its= %7.0f, tol=%10.3e\n', [k tol]);
        return;
    end
    beta = (rnew'*rnew)/(r'*r);
    p = rnew + beta*p;
    r = rnew;
end
fprintf('Caution: CG went to max iterations without converging!\n');
fprintf('MAXITS = %7.0f, tol =%10.3e\n', [MAXITS tol]);

end
