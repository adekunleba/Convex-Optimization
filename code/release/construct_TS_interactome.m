function [ TS_interactome ] = construct_TS_interactome( interactome, expression_profile, varargin )
    params = inputParser;    
    params.addParameter('method', 'ActPro_L1', @(x) ischar(x)); 
    params.addParameter('expression_threshold', 0.75, @(x) isscalar(x)); 
    params.addParameter('rw', 1e-1, @(x) isscalar(x)); % from Magger, O., Waldman, Y. Y., Ruppin, B., & Sharan, R. (2012). Enhancing the prioritization of disease-causing genes through tissue specific protein interaction networks. PLoS Computational Biology, 8(9), e1002690. http://doi.org/10.1371/journal.pcbi.1002690    
    params.addParameter('alpha', 0.5, @(x) isscalar(x) & x >= 0 & x <= 1);
    params.addParameter('lambda', 1, @(x) isscalar(x) & x >= 0);

    params.parse(varargin{:});
    par = params.Results;
    if(par.alpha ~= 0.5) 
        par.lambda = (1 - par.alpha) / par.alpha;
    end
    
    n = size(interactome, 1);
    m = nnz(interactome) / 2;
    
    TS_interactome = interactome;
    expressed_genes = par.expression_threshold <= expression_profile;
    
    switch(par.method)
        case('NR') % induced subgraph of expressed genes
            TS_interactome(expressed_genes, expressed_genes) = interactome(expressed_genes, expressed_genes);
            TS_interactome(~expressed_genes, :) = 0;
            TS_interactome(:, ~expressed_genes) = 0;
            
        case('ERW') % induced subgraph of expressed genes
            penalty = ones(size(interactome, 1), 1);
            penalty(~expressed_genes, 1) =  par.rw;
            TS_interactome = bsxfun(@times, interactome, penalty); % Penalize columns (destinations)
            TS_interactome = bsxfun(@times, TS_interactome, penalty'); % Penalize rows (sources)
            
        case('adaptive_ERW')
            x = expression_profile;
            smoothed_activity = n*x ./ sum(x);
                                   
            TS_interactome = bsxfun(@times, interactome, smoothed_activity'); % Penalize rows (sources)            
            TS_interactome = bsxfun(@times, TS_interactome, smoothed_activity); % Penalize columns (destinations)             
            

        case('PenPro')
            n = size(interactome, 1);
            W = interactome;
            D = diag(sum(W));
            L = D - W;
            e = double(par.expression_threshold <= expression_profile); % ==> Binarize expression values first
            lambda = (n*par.alpha) / (m*(1 - par.alpha));
            smoothed_penalty = (speye(n) + lambda*L)\e; 
            
            TS_interactome = bsxfun(@times, interactome, smoothed_penalty'); % Penalize rows (sources)            
            TS_interactome = bsxfun(@times, TS_interactome, smoothed_penalty); % Penalize columns (destinations)            
            
        case('ActPro')
            n = size(interactome, 1);
            W = interactome;
            D = diag(sum(W));
            L = D - W;
            e = expression_profile;
            lambda = (n*par.alpha) / (m*(1 - par.alpha));
            smoothed_activity = (speye(n) + lambda*L)\e; 
            
            TS_interactome = bsxfun(@times, interactome, smoothed_activity'); % Penalize rows (sources)            
            TS_interactome = bsxfun(@times, TS_interactome, smoothed_activity); % Penalize columns (destinations)            
            
        case('ActPro_L1')
            W = interactome;
            
            [ii, jj, vv] = find(tril(W));
            B = sparse([1:m, 1:m]', [ii; jj], [vv; -vv], m, n);
            
            expression_profile(expression_profile < 0) = 0;
            z = n*expression_profile ./ norm(expression_profile, 1);
            cvx_par.budget = n;
            
            tic
            cvx_begin quiet
                cvx_precision best

                variable x(n, 1)
                variable y(m, 1)

                minimize( n*sum_square(y) + par.lambda*m*norm(x-z, 1) );
                subject to
                    y == B*x;
                    0 <= x <= cvx_par.budget;
                    sum(x) == cvx_par.budget;
            cvx_end                          
            toc;

            x(x < 0) = 0;
            smoothed_activity = n*x ./ sum(x);
                                   
            TS_interactome = bsxfun(@times, interactome, smoothed_activity'); % Penalize rows (sources)            
            TS_interactome = bsxfun(@times, TS_interactome, smoothed_activity); % Penalize columns (destinations)            


        case('ActPro_L1_unweighted')
            W = interactome;           
            [ii, jj] = find(tril(W));
            
            B = sparse([1:m, 1:m]', [ii; jj], [ones(m,1); -ones(m,1)], m, n);
            
            expression_profile(expression_profile < 0) = 0;
            z = n*expression_profile ./ norm(expression_profile, 1);
            cvx_par.budget = n;
            
            tic
            cvx_begin quiet
                cvx_precision best

                variable x(n, 1)
                variable y(m, 1)

                minimize( n*sum_square(y) + par.lambda*m*norm(x-z, 1) );
                subject to
                    y == B*x;
                    0 <= x <= cvx_par.budget;
                    sum(x) == cvx_par.budget;
            cvx_end                          
            toc;
            
            Delta = z - x;
            x(x < 0) = 0;
            smoothed_activity = n*x ./ sum(x);
                                   
            TS_interactome = bsxfun(@times, interactome, smoothed_activity'); % Penalize rows (sources)            
            TS_interactome = bsxfun(@times, TS_interactome, smoothed_activity); % Penalize columns (destinations)            

        case('PenPro_L1')
            W = interactome;                      
            [ii, jj, vv] = find(tril(W));            
            B = sparse([1:m, 1:m]', [ii; jj], [vv; -vv], m, n);            

            expression_profile(expression_profile < 0) = 0;
            e = double(par.expression_threshold <= expression_profile); % ==> Binarize expression values first            
            z = n*e ./ norm(e, 1);
            
            cvx_par.budget = n;
            
            tic
            cvx_begin quiet
                cvx_precision best

                variable x(n, 1)
                variable y(m, 1)

                minimize( n*sum_square(y) + par.lambda*m*norm(x-z, 1) );
                subject to
                    y == B*x;
                    0 <= x <= cvx_par.budget;
                    sum(x) == cvx_par.budget;
            cvx_end                          
            toc;
            
            Delta = z - x;
            x(x < 0) = 0;
            smoothed_activity = n*x ./ sum(x);
                                   
            TS_interactome = bsxfun(@times, interactome, smoothed_activity'); % Penalize rows (sources)            
            TS_interactome = bsxfun(@times, TS_interactome, smoothed_activity); % Penalize columns (destinations)            

            
        case('ActPro_L1_POGS')
            W = interactome;
            
            [ii, jj, vv] = find(tril(W));
            B = sparse([1:m, 1:m]', [ii; jj], [vv; -vv], m, n);
            
            expression_profile(expression_profile < 0) = 0;
            z = n*expression_profile ./ norm(expression_profile, 1);
            cvx_par.budget = n;
            

            A = [B; ones(1, n); speye(n)];
            f.h = [kSquare(m); kIndEq0; kAbs(n)];
            f.b = [zeros(m, 1); cvx_par.budget; z];
            f.c = [n*ones(m, 1); 1; par.lambda*m*ones(n, 1)];
            g.h = kIndGe0;


            pogs_params.adaptive_rho = true;
            pogs_params.quiet = true;
            tic    
            [smoothed_activity, ~, ~, ~, ~, status] = pogs(A, f, g, pogs_params);            
            toc
            if status > 0
              fprintf('non-zero exit code %d\n', status);
            end            
            
            TS_interactome = bsxfun(@times, interactome, smoothed_activity'); % Penalize rows (sources)            
            TS_interactome = bsxfun(@times, TS_interactome, smoothed_activity); % Penalize columns (destinations)            

        otherwise
            fprintf('method %s is not recognized. Retunring original interactome', par.method);
    end
end

