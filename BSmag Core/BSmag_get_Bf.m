function [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_Bf(BSmag,X,Y,Z)
%---------------------------------------------------
%  NAME:      BSmag_get_Bf.m
%  WHAT:      Calculates B at field points.
%  REQUIRED:  BSmag Toolbox 20150407
%  Original AUTHOR:    20150407, L. Queval (loic.queval@gmail.com)
%  COPYRIGHT: 2015, Loic Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause).
%  MODIFIED 2019 from BSmag_get_B to fix offset-by-half error, by P. Picot (ppicot@uwo.ca) 
%  USE:
%    [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_Bf(BSmag,X,Y,Z)
%
%  INPUTS:
%    BSmag      = BSmag data structure
%    X          = Field points x-coordinate vector or matrix
%    Y          = Field points y-coordinate vector or matrix
%    Z          = Field points z-coordinate vector or matrix
%
%  OUTPUTS:
%    BSmag      = BSmag data structure (no update)
%    X          = Field points x-coordinate vector or matrix
%    Y          = Field points y-coordinate vector or matrix
%    Z          = Field points z-coordinate vector or matrix
%    BX         = Field points B x-component vector or matrix
%    BY         = Field points B y-component vector or matrix
%    BZ         = Field points B z-component vector or matrix
%----------------------------------------------------


mu0 = 4*pi*1e-7; % vacuum permeability [N/A^2]
           
BX = zeros(size(X,1),size(X,2),size(X,3));
BY = zeros(size(X,1),size(X,2),size(X,3));
BZ = zeros(size(X,1),size(X,2),size(X,3));

for nF = 1:BSmag.Nfilament % Loop on each filament

    Gamma = BSmag.filament(nF).Gamma;

%    dGamma repurposed as oversample param here, in a way to not break (most) existing code
    if ( (BSmag.filament(nF).dGamma > 2)  && (BSmag.filament(nF).dGamma < 10) )
        os = floor(BSmag.filament(nF).dGamma);
        fprintf('Note: dGamma=%d overrides default(=3) oversampling of Gamma.', os);  
    else
        os = 3;  % oversampling factor
    end
    
    I = BSmag.filament(nF).I;

    % Discretization of Gamma
 
    N = size(Gamma,1); % Number of points defining Gamma
    gamlen = N * os - 1; % number of points defining oversampled Gamma
    
    x_P = interp1(Gamma(:,1), linspace(1,N,gamlen), 'spline'); % discretization of Gamma for x component
    y_P = interp1(Gamma(:,2), linspace(1,N,gamlen), 'spline'); % discretization of Gamma for y component
    z_P = interp1(Gamma(:,3), linspace(1,N,gamlen), 'spline'); % discretization of Gamma for z component
    
    % delta-el for each point straddles the point.  Principal difference from parent code.
    for k=2:gamlen-1
        dx_P(k) = ( x_P(k+1) - x_P(k-1) )/2;
        dy_P(k) = ( y_P(k+1) - y_P(k-1) )/2;
        dz_P(k) = ( z_P(k+1) - z_P(k-1) )/2;
    end
        
    % end points are special cases:
    %if ( (norm(x_P(1)-x_P(gamlen)) + norm(y_P(1)-y_P(gamlen)) + norm(z_P(1)-z_P(gamlen)) ) < 1e-12 )
    if ( norm(Gamma(1) - Gamma(N)) < 1e-12 )
       % disp('Closed loop case with overlap found. Averaging at overlap.');
        % if the two endpoints are coincident and duplicate (a closed loop), then compute d*_P 
        % from the adjacent points on both sides, and assign half the current 
        % to each point by halving d*_P (i.e., quarter of the 2-point span).
        dx_P(1) = (x_P(2)-x_P(gamlen-1))/4;
        dy_P(1) = (y_P(2)-y_P(gamlen-1))/4;
        dz_P(1) = (z_P(2)-z_P(gamlen-1))/4;

        dx_P(gamlen) = dx_P(1);
        dy_P(gamlen) = dy_P(1);
        dz_P(gamlen) = dz_P(1);


    else
       % disp('Open loop case found.  Extrapolating past ends.');
        % else, if the the two endpoints are not coincident, then we can't 
        % extrapolate to the next (unknown) point, so just assign
        % the d* from the last two known points:
        dx_P(1) = x_P(2)-x_P(1);
        dy_P(1) = y_P(2)-y_P(1);
        dz_P(1) = z_P(2)-z_P(1);

        dx_P(gamlen) = x_P(gamlen) - x_P(gamlen-1);
        dy_P(gamlen) = y_P(gamlen) - y_P(gamlen-1);
        dz_P(gamlen) = z_P(gamlen) - z_P(gamlen-1);
    end
    
    % Add contribution of each source point P on each field point M (where we want to calculate the field)
    for m = 1:size(X,1)
        for n = 1:size(X,2)
            for p = 1:size(X,3)

            % M is the field point
            x_M = X(m,n,p);
            y_M = Y(m,n,p);
            z_M = Z(m,n,p);

            % Loop on each discretized segment of Gamma PkPk+1
            for k = 1:gamlen
                PkM3 = (sqrt((x_M-x_P(k))^2 + (y_M-y_P(k))^2 + (z_M-z_P(k))^2))^3;
                DBx(k) = (dy_P(k)*(z_M-z_P(k)) - dz_P(k)*(y_M-y_P(k)))/PkM3;
                DBy(k) = (dz_P(k)*(x_M-x_P(k)) - dx_P(k)*(z_M-z_P(k)))/PkM3;
                DBz(k) = (dx_P(k)*(y_M-y_P(k)) - dy_P(k)*(x_M-x_P(k)))/PkM3;
            end
            % Sum
            BX(m,n,p) = BX(m,n,p) + mu0*I/4/pi*sum(DBx);
            BY(m,n,p) = BY(m,n,p) + mu0*I/4/pi*sum(DBy);
            BZ(m,n,p) = BZ(m,n,p) + mu0*I/4/pi*sum(DBz);

            end
        end
    end
end
