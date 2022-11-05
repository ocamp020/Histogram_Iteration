%%% VFI Toolbox 

classdef VFI_Toolbox

    methods (Static)

        
        % -------------------------------------------------------------------------------------------
        % -------------------------------------------------------------------------------------------
        

        % Rouwenhorst (1995) - from Karen Kopecky (2010)

        %[zgrid, P] = rouwenhorst(rho, sigma_eps, n)
        %
        % rho is the 1st order autocorrelation
        % sigma_eps is the standard deviation of the error term
        % n is the number of points in the discrete approximation
        %
        function [zgrid, P] = Kopecky_rouwenhorst(rho,sigma_eps,n)
        
            mu_eps = 0;
            
            q = (rho+1)/2;
            nu = ((n-1)/(1-rho^2))^(1/2) * sigma_eps;
            
            P = [q 1-q;1-q q];
            
            
            for i=2:n-1
               P = q*[P zeros(i,1);zeros(1,i+1)] + (1-q)*[zeros(i,1) P;zeros(1,i+1)] + ...
                   (1-q)*[zeros(1,i+1); P zeros(i,1)] + q*[zeros(1,i+1); zeros(i,1) P];
               P(2:i,:) = P(2:i,:)/2;
            end
            
            zgrid = linspace(mu_eps/(1-rho)-nu,mu_eps/(1-rho)+nu,n);

        end

        

        % -------------------------------------------------------------------------------------------
        % -------------------------------------------------------------------------------------------


        function [N,grid,Pi,PDF,CDF] = Rouwenhorst95(rho, sigma, N)

            p = (1+rho)/2;
            q = p;
            psi = sigma*sqrt((N-1)/(1-rho^2));
            s = (1-q)/(2-(p+q));

	        % Check if N>1
		    if N==1
                    z=0;
                    Pi = 1;
                    PDF_z = 1;
                    
                    % Return 
                    grid = 0;
                    Pi = 1;
                    PDF = 1;
                    CDF = 1;
		    
		    else
			    P = [p 1-p;1-p p];
			    for i=2:N-1
				    P = q*[P zeros(i,1);zeros(1,i+1)] + (1-q)*[zeros(i,1) 	P;zeros(1,i+1)] + ...
                        (1-q)*[zeros(1,i+1); P zeros(i,1)] + q*[zeros(1,i+1); zeros(i,1) P];
        	        P(2:i,:) = P(2:i,:)/2;
                end
		    
		    % Distribution
                    PDF_z = binopdf(0:N-1,N-1,1-s);
                    CDF_z = cumsum(PDF_z);
                    % Create z grid
                    z = linspace(-psi,psi,N);
                    % Return    
                    grid = z;
                    Pi = P;
                    PDF = PDF_z;
                    CDF = CDF_z;
            end
        end

        % -------------------------------------------------------------------------------------------
        % -------------------------------------------------------------------------------------------
        
        function [N,grid,Pi,PDF,CDF] = Tauchen86(rho, sigma, N,Omega)
            % Check if N>1
            if N==1
                z=0;
                Pi=1;
                PDF=1;
                grid = z;
                CDF=1;
            else
                % Create z grid
                z = linspace(-Omega*sigma/((1-rho^2)^(1/2)),Omega*sigma/((1-rho^2)^(1/2)),N);
                % Define intermediate step length
                h = (z(2)-z(1))/2;
                % Define auxiliary matrices
                z_0 = repmat(transpose(z),1,N); % Matrix of today's z each row is a value, columns are equal
                z_1 = repmat(z,N,1); % Matrix of tomorrow's z each column is a value, rows are equal
                % Define intervals
                z_lim = zeros(N,N,2); % First matrix is lower bounds. Second matrix is upper bounds.
                z_lim(:,1,1)    = - Inf;
                z_lim(:,2:end  ,1) = (z_1(:,2:end  ) - rho.*z_0(:,2:end  ) - h )./sigma;
                z_lim(:,1:end-1,2) = (z_1(:,1:end-1) - rho.*z_0(:,1:end-1) + h)./sigma;
                z_lim(:,end,2)  = Inf;
                % Fill in transition matrix
                Pi_z = normcdf(z_lim(:,:,2)) - normcdf(z_lim(:,:,1));
                Pi_z = Pi_z./sum(Pi_z,2);
                Pi_z;
                %Get stationary distribution of markov chain
                [V,D,W] = eig(transpose(Pi_z));
                PDF_z = real(V(:,1));
                PDF_z = PDF_z./sum(PDF_z);
                CDF_z = cumsum(PDF_z);
                grid = z;
                Pi = Pi_z;
                PDF = PDF_z;
                CDF = CDF_z;
            end
        end

        % -------------------------------------------------------------------------------------------
        % -------------------------------------------------------------------------------------------

        % Make Grid 
        % Make_Grid(n,θ_x,x_min,x_max,scale_type="Poly")
        % n: number of points
        % θ_x: Curvature parameter
        % x_min: minimum value for grid
        % x_max: maximum value for grid
        % scale_type = "Poly" for Polynomial Scaling
        % scale_type = "Exp" for Exponential Scaling
        function x_grid = Make_Grid(n,theta_x, x_min, x_max, scale_type)
            if theta_x~=1
                grid = linspace(0,1,n);
                if scale_type=="Poly"
                    x_grid = x_min + (x_max - x_min)*grid.^(theta_x);
                elseif scale_type=="Exp"
                    x_grid = x_min + (x_max - x_min)*((exp(theta_x*grid)-1)/(exp(theta_x)-1));
                end
            else
                x_grid = linspace(x_min, x_max, n);
            end
        end

        % -------------------------------------------------------------------------------------------
        % -------------------------------------------------------------------------------------------
        
        % Grid Inv
        % Grid_Inv(x,n,θ_x,x_min,x_max,scale_type="Poly")
        % x: point of the grid
        % n: number of points in the grid
        % θ_x: Curvature parameter
        % x_min: minimum value for grid
        % x_max: maximum value for grid
        % scale_type = "Poly" for Polynomial Scaling
        % scale_type = "Exp" for Exponential Scaling

        function x_ind = Grid_Inv(x, n, theta_x, x_min, x_max, scale_type)
            % Check corner solution
            if x<x_min
                x_ind = 1
            elseif x>x_max
                x_ind = n
            else
                if theta_x~=1
                    if scale_type=="Poly"
                        x_ind = int32(floor(((x-x_min)/(x_max-x_min)).^(1/theta_x) * (n-1)+1));
                    elseif scale_type=="Exp"
                        x_ind = int32(floor((1/theta_x).*log.(((x-x_min)*(exp(theta_x)-1)/(x_max-x_min))+1)*(n-1)+1));
                    else
                        error("Scale_type must be either Poly or Exp")
                    end
                else
                    x_ind = int32(floor(((x-x_min)/(x_max-x_min))*(n-1)+1))
                end
            end
        end


        % -------------------------------------------------------------------------------------------
        % -------------------------------------------------------------------------------------------
        
      

        
    end % Ends Static Methods

    


end % Ends Class VFI_Toolbox

