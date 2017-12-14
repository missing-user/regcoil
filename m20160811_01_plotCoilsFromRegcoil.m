% Averaged
%regcoilFilenames = {'/Users/elizabethpaul/Documents/Research/Fall_2017/20171121_comparing_W7X_offset_and_actual/W7X_averaged/regcoil_out.w7x.nc'};
%nescinFilenames = {'/Users/elizabethpaul/Documents/Research/Fall_2017/20171121_comparing_W7X_offset_and_actual/W7X_averaged/nescin.averaged'};
%ilambda = 9;
% Offset
% regcoilFilenames = {'/Users/elizabethpaul/Documents/Research/Fall_2017/20171121_comparing_W7X_offset_and_actual/W7X_offset/regcoil_out.w7x.nc'};
% nescinFilenames = {'/Users/elizabethpaul/Documents/Research/Fall_2017/20171121_comparing_W7X_offset_and_actual/W7X_offset/nescin.w7x_winding_surface_from_Drevlak_5013_116'};
% ilambda = 10;
% Actual
regcoilFilenames = {'/Users/elizabethpaul/Documents/Research/Fall_2017/20171121_comparing_W7X_offset_and_actual/W7X_init/regcoil_out.w7x.nc'};
nescinFilenames = {'/Users/elizabethpaul/Documents/Research/Fall_2017/20171121_comparing_W7X_offset_and_actual/W7X_init/nescin.w7x_winding_surface_from_Drevlak_0'};
ilambda = 8;

thetaShift = 0;

coilsPerHalfPeriod=5;
numHalfPeriodsToPlot=1;

coil_thickness = 0.06;

colors = [1,0,0;
    1,0.7,0;
    0,0.8,0;
    0,0,1;
    1,0,1];

ntheta=150;
nzeta=160;
figure(4)
clf

for whichFile = 1:1
    % Read regcoil_out file:
    filename = regcoilFilenames{whichFile};
    fprintf(['Reading ',filename,'\n'])

    nfp = double(ncread(filename,'nfp'));
    chi2_B = ncread(filename,'chi2_B');
    chi2_K = ncread(filename,'chi2_K');
    fprintf('chi2_B: %g,  chi2_K: %g\n',chi2_B(end),chi2_K(end))
    net_poloidal_current_Amperes = ncread(filename,'net_poloidal_current_Amperes');
    theta = ncread(filename,'theta_coil');
    theta = circshift(theta,thetaShift);
    for itheta=1:(length(theta)-1)
       if (theta(itheta) > theta(itheta+1))
           theta(itheta+1) = theta(itheta+1)+2*pi;
       end
    end

    nzeta = double(ncread(filename,'nzeta_coil'));
    nzetal=nzeta*nfp;
    zetal = linspace(0,2*pi,nzetal+1);
    zetal(end)=[];
    [zetal_2D, theta_2D] = meshgrid(zetal,theta);
    potential0 = ncread(filename,'current_potential');
    potential1 = potential0(:,:,ilambda);
    potential1 = circshift(potential1,thetaShift,1);
    %size(potential0)
    %fprintf('Here comes current potential:\n')
    %potential0(:,:,ilambda)
    potential = kron(ones(1,nfp),potential1) + kron(((1:nfp)-1)*net_poloidal_current_Amperes/nfp,ones(numel(theta),nzeta));
    potential = potential / net_poloidal_current_Amperes * nfp;
    %fprintf('Here comes final potential:\n')
    %potential
    %fprintf('Here comes size(potential):\n')
    %size(potential)
    %return
    fprintf('min/max of potential1: %g / %g\n',min(min(potential1)), max(max(potential1)))
    fprintf('min/max of potential:  %g / %g\n',min(min(potential)), max(max(potential)))

    % Read surface from nescin file:
    filename = nescinFilenames{whichFile};
    fprintf(['Reading ',filename,'\n'])
    fid = fopen(filename,'r');
    %{
    line = fgetl(fid);  % Grid spatial dimensions
    line = fgetl(fid);  % nu,nv,...
    line = fgetl(fid);  %  256 64 256 64 64 10 F
    line = fgetl(fid);  % <blank>
    line = fgetl(fid);  % Fourier Dimensions
    line = fgetl(fid);  % mf, nf, md, nd
    line = fgetl(fid);  % 4 4 24 20
    data = sscanf(line,'%d %d %d %d');
    mf=data(1);
    nf=data(2);
    line = fgetl(fid); % <blank>
    line = fgetl(fid); % Plasma information from VMEC
    line = fgetl(fid); % np, iota_edge, phip_edge, curpol
    line = fgetl(fid); % <values>
    line
    data = sscanf(line,'%d %g %g %g');
    nfp=data(1);
    iota_edge=data(2);
    phip_edge = data(3);
    curpol = data(4);
    fprintf('Read nfp=%d\n',nfp)
    line = fgetl(fid); % <blank>
    line = fgetl(fid); % Current controls
    line = fgetl(fid); % cut, cup, ibex
    line = fgetl(fid); % <values>
    data = sscanf(line,'%g %g %g');
    cut=data(1);
    cup=data(2);
    %}
    %{
    nzetal=nzeta*nfp;
    theta = linspace(0,2*pi,ntheta+1);
    theta(end)=[];
    zetal = linspace(0,2*pi,nzetal+1);
    zetal(end)=[];
    [zetal_2D, theta_2D] = meshgrid(zetal,theta);
    %}

    %search_string = '----- Coil Surface';
    search_string = '------ Current Surface'
    while true
        line = fgetl(fid);
        if strncmp(line,search_string,numel(search_string))
            break
        end
    end
    line = fgetl(fid); %Number of fourier modes in table
    line = fgetl(fid);
    mnmax_nescin = sscanf(line,'%d');
    fprintf('  Reading %d modes from nescin file %s\n',mnmax_nescin,filename)
    line = fgetl(fid); %Table of fourier coefficients
    line = fgetl(fid); %m,n,crc2,czs2,crs2,czc2
    xm_nescin = zeros(mnmax_nescin,1);
    xn_nescin = zeros(mnmax_nescin,1);
    rmnc_nescin = zeros(mnmax_nescin,1);
    zmns_nescin = zeros(mnmax_nescin,1);
    for i=1:mnmax_nescin
        line = fgetl(fid);
        data = sscanf(line,'%d %d %g %g %g %g %g %g');
        xm_nescin(i) = data(1);
        xn_nescin(i) = data(2);
        rmnc_nescin(i) = data(3);
        zmns_nescin(i) = data(4);
    end

    %{
    % Build current potential
    v = zetal_2D*nfp/(2*pi);
    u = theta_2D/(2*pi);
    It = cut;
    Ip = cup*nfp;
    %potential = -Ip*v/nfp - It*u; %The signs are wrong here!!
    potential = Ip*v/nfp + It*u;
    search_string = '---- Phi(m,n) for least squares ---';
    while true
        line = fgetl(fid);
        if strncmp(line,search_string,numel(search_string))
            break
        end
    end
    for m = 0:mf
        for n = (-nf):nf
            line = fgetl(fid);
            data = sscanf(line,'%d %d %g');
            potential = potential + data(3)*sin(2*pi*(m*u+n*v));
        end
    end
    %}
    fclose(fid);
    % Done reading nescin file.
    fprintf(['Done reading ',filename,'\n'])
    
    
    figure(whichFile*10+1)
    clf
    
    % W7-X has 50 modular coils
    %contours = linspace(0,2*pi,51);
    %contours = linspace(-nfp,0,51); %This line works if the sign for the secular part of the current potential is wrong.
    %contours = linspace(0,nfp,51);
    contours = linspace(0,nfp,1+coilsPerHalfPeriod*2*nfp);
    contours(end)= [];
    dc = contours(2)-contours(1);
    contours = contours + 0.5*dc;
    
    contourf(zetal_2D,theta_2D,potential,contours)
    %contourf(zetal_2D,theta_2D,potential,25)
    hold on
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title('Current potential')
    set(gcf,'Position',[12         374        1248         313])
    %return
    
    contours_theta = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    contours_zeta = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    contours_x = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    contours_y = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    contours_z = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    contours_dxdtheta = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    contours_dydtheta = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    contours_dzdtheta = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    contours_dxdzeta = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    contours_dydzeta = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    contours_dzdzeta = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    coils_x = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    coils_y = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    coils_z = cell(coilsPerHalfPeriod*numHalfPeriodsToPlot,1);
    for j=1:coilsPerHalfPeriod*numHalfPeriodsToPlot
        this_contour = contours(j+2*coilsPerHalfPeriod);
        C = contourc(zetal,theta,potential,[this_contour,this_contour]);
        N = C(2,1);
        if N ~= size(C,2)-1
            fprintf('It appears there are multiple disconnected contours. This program presently cannot handle this.\n')
            N
            size(C)
        end
        this_zeta = C(1,2:end)';
        this_theta = C(2,2:end)';
        contours_zeta{j} = [this_zeta; this_zeta(1)];
        contours_theta{j}  = [this_theta; this_theta(1)];
        plot(contours_zeta{j},contours_theta{j},'r','LineWidth',2)
        contours_x{j} = zeros(size(contours_theta{j}));
        contours_y{j} = zeros(size(contours_theta{j}));
        contours_z{j} = zeros(size(contours_theta{j}));
        contours_dxdtheta{j} = zeros(size(contours_theta{j}));
        contours_dydtheta{j} = zeros(size(contours_theta{j}));
        contours_dzdtheta{j} = zeros(size(contours_theta{j}));
        contours_dxdzeta{j} = zeros(size(contours_theta{j}));
        contours_dydzeta{j} = zeros(size(contours_theta{j}));
        contours_dzdzeta{j} = zeros(size(contours_theta{j}));
    end
    
    x = zeros(size(theta_2D));
    y = zeros(size(theta_2D));
    z = zeros(size(theta_2D));
    dxdtheta = zeros(size(theta_2D));
    dydtheta = zeros(size(theta_2D));
    dzdtheta = zeros(size(theta_2D));
    dxdzeta = zeros(size(theta_2D));
    dydzeta = zeros(size(theta_2D));
    dzdzeta = zeros(size(theta_2D));
    
    for i = 1:mnmax_nescin
        angle = xm_nescin(i)*theta_2D + xn_nescin(i)*zetal_2D*nfp;
        %angle2 = zetal_2D + pi; % Rotate everything in the x-y plane
        angle2 = zetal_2D; 
        
        x = x + rmnc_nescin(i)*cos(angle).*cos(angle2);
        y = y + rmnc_nescin(i)*cos(angle).*sin(angle2);
        z = z + zmns_nescin(i)*sin(angle);
        
        %{
    dxdtheta = dxdtheta - xm_nescin(i)*rmnc_nescin(i)*sin(angle).*cos(angle2);
    dydtheta = dydtheta - xm_nescin(i)*rmnc_nescin(i)*sin(angle).*sin(angle2);
    dzdtheta = dzdtheta + xm_nescin(i)*zmns_nescin(i)*cos(angle);
    
    dxdzeta = dxdzeta - nfp*xn_nescin(i)*rmnc_nescin(i)*sin(angle).*cos(angle2) ...
        - rmnc_nescin(i)*cos(angle).*sin(angle2);
    dydzeta = dydzeta - nfp*xn_nescin(i)*rmnc_nescin(i)*sin(angle).*sin(angle2) ...
        + rmnc_nescin(i)*cos(angle).*cos(angle2);
    dzdzeta = dzdzeta + nfp*xn_nescin(i)*zmns_nescin(i)*cos(angle);
        %}
        
        for j=1:coilsPerHalfPeriod*numHalfPeriodsToPlot
            angle = xm_nescin(i)*contours_theta{j} + xn_nescin(i)*contours_zeta{j}*nfp;
            angle2 = contours_zeta{j};
            
            contours_x{j} = contours_x{j} + rmnc_nescin(i)*cos(angle).*cos(angle2);
            contours_y{j} = contours_y{j} + rmnc_nescin(i)*cos(angle).*sin(angle2);
            contours_z{j} = contours_z{j} + zmns_nescin(i)*sin(angle);
            
            contours_dxdtheta{j} = contours_dxdtheta{j} - xm_nescin(i)*rmnc_nescin(i)*sin(angle).*cos(angle2);
            contours_dydtheta{j} = contours_dydtheta{j} - xm_nescin(i)*rmnc_nescin(i)*sin(angle).*sin(angle2);
            contours_dzdtheta{j} = contours_dzdtheta{j} + xm_nescin(i)*zmns_nescin(i)*cos(angle);
            
            contours_dxdzeta{j} = contours_dxdzeta{j} - nfp*xn_nescin(i)*rmnc_nescin(i)*sin(angle).*cos(angle2) ...
                - rmnc_nescin(i)*cos(angle).*sin(angle2);
            contours_dydzeta{j} = contours_dydzeta{j} - nfp*xn_nescin(i)*rmnc_nescin(i)*sin(angle).*sin(angle2) ...
                + rmnc_nescin(i)*cos(angle).*cos(angle2);
            contours_dzdzeta{j} = contours_dzdzeta{j} + nfp*xn_nescin(i)*zmns_nescin(i)*cos(angle);
            
        end
    end
    
    for j=1:coilsPerHalfPeriod*numHalfPeriodsToPlot
        % Compute normal direction:
        Nx = contours_dydzeta{j} .* contours_dzdtheta{j} - contours_dzdzeta{j} .* contours_dydtheta{j};
        Ny = contours_dzdzeta{j} .* contours_dxdtheta{j} - contours_dxdzeta{j} .* contours_dzdtheta{j};
        Nz = contours_dxdzeta{j} .* contours_dydtheta{j} - contours_dydzeta{j} .* contours_dxdtheta{j};
        norm_normal = sqrt(Nx.*Nx + Ny.*Ny + Nz.*Nz);
        Nx = Nx ./ norm_normal;
        Ny = Ny ./ norm_normal;
        Nz = Nz ./ norm_normal;
        
        % Compute tangent direction:
        indices = (1:numel(contours_x{j}))';
        next_index = circshift(indices,[-1,0]);
        prev_index = circshift(indices,[1,0]);
        Tx = contours_x{j}(next_index) - contours_x{j}(prev_index);
        Ty = contours_y{j}(next_index) - contours_y{j}(prev_index);
        Tz = contours_z{j}(next_index) - contours_z{j}(prev_index);
        norm_tangent = sqrt(Tx.*Tx + Ty.*Ty + Tz.*Tz);
        Tx = Tx ./ norm_tangent;
        Ty = Ty ./ norm_tangent;
        Tz = Tz ./ norm_tangent;
        
        % Compute binormal:
        Bx = Ty .* Nz - Tz .* Ny;
        By = Tz .* Nx - Tx .* Nz;
        Bz = Tx .* Ny - Ty .* Nx;
        
        coils_x{j} = [...
            contours_x{j} + coil_thickness*(Nx+Bx), ...
            contours_x{j} + coil_thickness*(Nx-Bx), ...
            contours_x{j} + coil_thickness*(-Nx-Bx), ...
            contours_x{j} + coil_thickness*(-Nx+Bx), ...
            contours_x{j} + coil_thickness*(Nx+Bx)];
        
        coils_y{j} = [...
            contours_y{j} + coil_thickness*(Ny+By), ...
            contours_y{j} + coil_thickness*(Ny-By), ...
            contours_y{j} + coil_thickness*(-Ny-By), ...
            contours_y{j} + coil_thickness*(-Ny+By), ...
            contours_y{j} + coil_thickness*(Ny+By)];
        
        coils_z{j} = [...
            contours_z{j} + coil_thickness*(Nz+Bz), ...
            contours_z{j} + coil_thickness*(Nz-Bz), ...
            contours_z{j} + coil_thickness*(-Nz-Bz), ...
            contours_z{j} + coil_thickness*(-Nz+Bz), ...
            contours_z{j} + coil_thickness*(Nz+Bz)];
    end
    
    figure(whichFile*10+2)
    surf(x,y,z)
    daspect([1,1,1])
    axis vis3d
    hold on
    
    for j=1:coilsPerHalfPeriod*numHalfPeriodsToPlot
        plot3(contours_x{j},contours_y{j},contours_z{j},'r','LineWidth',2)
    end
    
    light
    lighting gouraud
    %zoom(1.6)
    %campos([  574.9370 -457.0244  424.3304])
    %camva(1.0271)
    axis off
    
    
    figure(whichFile*10+3)
    clf
    for j=1:coilsPerHalfPeriod*numHalfPeriodsToPlot
        plot3(contours_x{j},contours_y{j},contours_z{j},'.-r','LineWidth',2,'MarkerSize',15)
        hold on
    end
    daspect([1,1,1])
    axis vis3d
    axis off
    
    
    figure(4)
    if whichFile==1
        clf
        edgeColor='k';
        offset=0;
    else
        edgeColor=':k';
        alpha=0.5;
        colors = alpha*colors + (1-alpha)*ones(5,3);
        offset=0;
    end
    set(gcf,'Color','w')
    ambientStrength = 0.5;
    diffuseStrength = 1;
    for j=1:coilsPerHalfPeriod*numHalfPeriodsToPlot
        nextColor = mod(j-1,size(colors,1))+1;
        surf(-coils_x{j}+offset,-coils_y{j},coils_z{j},'EdgeColor','none','FaceColor',colors(nextColor,:),'AmbientStrength',ambientStrength,'DiffuseStrength',diffuseStrength)
        hold on
        for k=1:4
            plot3(-coils_x{j}(:,k)+offset,-coils_y{j}(:,k),coils_z{j}(:,k),edgeColor,'LineWidth',1.3)
        end
    end
end

daspect([1,1,1])
axis vis3d
axis off
%light
%campos([-1.4267   18.9647   29.7379])
campos([-2.3196   24.0172  -17.4164])
camlight
%lighting gouraud
camva(6)

%annotation(gcf,'textbox',[0.19 0.85 0.5 0.03],'String','(d)','FontSize',18,'FitBoxToText','off','LineStyle','none');
%annotation(gcf,'textbox',[0.1 0.85 0.9 0.03],'String','REGCOIL','FontSize',18,'FitBoxToText','off','LineStyle','none','HorizontalAlignment','center');

%figureFilename = ['/Users/mattland/Box Sync/work16/REGCOIL paper/',mfilename];
%export_fig(figureFilename,'-png','-r200')
