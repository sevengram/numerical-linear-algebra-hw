function [ kappa ] = estimate_cond_no( A )
    [n,~] = size(A);
    norm_A = norm(A,1);
    norm_inv_A = NaN;
    for i=1:5
        w = rand(n,1);
        norm_inv_A = max(norm_inv_A, norm(A\w,1) / norm(w,1));
    end
    kappa = norm_A * norm_inv_A;
end
