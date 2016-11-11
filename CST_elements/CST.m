%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
% Elastic Constant Strain Triangular Elements                            %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace
clc
clear all

nod1= csvread('Nodes_1.csv');
nod2= csvread('Nodes_2.csv');
nod3= csvread('Nodes_3.csv');
nod4= csvread('Nodes_4.csv');
nod5= csvread('Nodes_5.csv');

Elm1=csvread('Elements_1.csv');
Elm2=csvread('Elements_2.csv');
Elm3=csvread('Elements_3.csv');
Elm4=csvread('Elements_4.csv');
Elm5=csvread('Elements_5.csv');

for no=1:5
    
    % Read nodes and coords
    if no==1
        Nodes = nod1;
    end
    if no==2
        Nodes = nod2;
    end
    if no==3
        Nodes = nod3;
    end
    if no==4
        Nodes = nod4;
    end
    if no==5
        Nodes = nod5;
    end
    [N,l] = size(Nodes);
    
    % Read element material id, thickness and nodal connectivity
    if no==1
        Elems = Elm1;
    end
    if no==2
        Elems = Elm2;
    end
    if no==3
        Elems = Elm3;
    end
    if no==4
        Elems = Elm4;
    end
    if no==5
        Elems = Elm5;
    end
    
    [E,l] = size(Elems);
    j_dbc=1;
    j_nbc=1;
    
    % Read material info
    Mats = load('Materials.txt');
    [M,l] = size(Mats);
    
    %Determine Derichlet BC
    for (i=1:N)
        if (Nodes(i,2)==0)
            DBC(j_dbc,1)=Nodes(i,1);
            DBC(j_dbc,2)=1;
            DBC(j_dbc,3)=0;
            j_dbc=j_dbc+1;
        end
        if (Nodes(i,3)==0)
            DBC(j_dbc,1)=Nodes(i,1);
            DBC(j_dbc,2)=2;
            DBC(j_dbc,3)=0;
            j_dbc=j_dbc+1;
        end
    end
    [P,l] = size(DBC);
    % Determine Neumann BC
    for (i=1:N)
        if (Nodes(i,2)==3)
            right(j_nbc,1)=Nodes(i,1);
            right(j_nbc,2)=1;
            right(j_nbc,3)=0;
            j_nbc=j_nbc+1;
        end
    end
    j_nbc=1;
    
    
    comb=combnk(right(:,1),2);
    
    for i=1:E
        for j=1:size(comb(:,1))
            if comb(j,:) == Elems(i,4:5) | comb(j,:) == Elems(i,5:6) | comb(j,:) == Elems(i,[4 6])| comb(j,:) == Elems(i,[6 4]) | comb(j,:) == Elems(i,[5 4]) | comb(j,:) == Elems(i,[6 5])
                NBC(j_nbc,1)=Elems(i,1);
                NBC(j_nbc,2:3)=comb(j,:);
                NBC(j_nbc,4)=1;
                NBC(j_nbc,5)=1;
                j_nbc=j_nbc+1;
            end
        end
    end
    [Q,l] = size(NBC);
    
    % Determining the hole nodes
    i_hol=1;
    for i=1:N
        if (Nodes(i,2)<=0.25 && Nodes(i,3)<=0.25)
            hole(i_hol)=Nodes(i,1);
            i_hol=i_hol+1;
        end
    end
    
    % Determining the hole elements
    i_hol=1;
    for i=1:E
        if (hole(i_hol)==Elems(i,4)||hole(i_hol)==Elems(i,5)||hole(i_hol)==Elems(i,6))
            hol_el(i_hol)=Elems(i,1);
            i_hol=i_hol+1;
        end
    end
    hol_el=unique(hol_el);
    
    % Identify out-of-plane conditions
    %   ipstrn = 1    Plane strain
    %   ipstrn = 2    Plane stress
    ipstrn = 2;
    
    % Determine total number of degrees-of-freedom
    udof = 2;     % Degrees-of-freedom per node
    NDOF = N*udof;
    
    % Initialize global matrix and vectors
    K = zeros(NDOF,NDOF);   % Stiffness matrix
    U = zeros(NDOF,1);      % Displacement vector
    F = zeros(NDOF,1);      % Force vector
    
    % Set penalty for displacement constraints
    Klarge = 10^10;
    
    % Loop over CST element
    for e = 1:E
        
        % Establish element connectivity and coordinates
        Nnums = Elems(e,4:6);
        xy = Nodes(Nnums(:),2:3);
        
        % Extract element thickness for plane stress
        h = Elems(e,3);
        
        % Extract element elastic Young's modulus and Poisson's ratio
        Y = Mats(Elems(e,2),2);
        nu = Mats(Elems(e,2),3);
        
        % Construct element stiffness matrix
        [Ke] = CST_El_Stiff(ipstrn,xy,h,Y,nu);
        
        % Assemble element stiffness matrix into global stiffness matrix
        ig = udof*(Nnums(:)-1);
        for ni = 1:3
            i0 = udof*(ni-1);
            for nj = 1:3
                j0 = udof*(nj-1);
                for i = 1:udof
                    for j = 1:udof
                        K(ig(ni)+i,ig(nj)+j) = K(ig(ni)+i,ig(nj)+j) + Ke(i0+i,j0+j);
                    end
                end
            end
        end
    end
    % K
    
    % Construct global force vector
    for q = 1:Q
        
        % Determine loaded edge
        e = NBC(q,1);
        in1 = NBC(q,2);
        in2 = NBC(q,3);
        idof = NBC(q,4);
        tval = NBC(q,5);
        h = Elems(e,3);
        
        % Establish edge length
        xlen2 = (Nodes(in2,2)-Nodes(in1,2))^2;
        ylen2 = (Nodes(in2,3)-Nodes(in1,3))^2;
        elen = sqrt(xlen2+ylen2);
        fval = tval*elen*h/2;
        iloc1 = udof*(in1-1)+idof;
        iloc2 = udof*(in2-1)+idof;
        F(iloc1) = F(iloc1) + fval;
        F(iloc2) = F(iloc2) + fval;
        F;
        
    end
    
    for p = 1:P
        inode = DBC(p,1);
        idof = DBC(p,2);
        idiag = udof*(inode-1) + idof;
        K(idiag,idiag) = Klarge;
        F(idiag) = Klarge*DBC(p,3);
    end
    % K
    % F
    
    % Solve system to determine displacements
    U = K\F;
    
    % Recover internal element displacement, strains and stresses
    Disp = zeros(E,6);
    Eps = zeros(E,3);
    Sig = zeros(E,3);
    
    for e = 1:E
        
        % Establish element connectivity and coordinates
        Nnums = Elems(e,4:6);
        xy = Nodes(Nnums(:),2:3);
        
        % Extract element thickness for plane stress
        h = Elems(e,3);
        
        % Extract element elastic Young's modulus and Poisson's ratio
        Y = Mats(Elems(e,2),2);
        nu = Mats(Elems(e,2),3);
        
        % Extract element nodal displacements
        inode1 = Nnums(1);
        inode2 = Nnums(2);
        inode3 = Nnums(3);
        Disp(e,1) = U(udof*(inode1-1)+1);
        Disp(e,2) = U(udof*inode1);
        Disp(e,3) = U(udof*(inode2-1)+1);
        Disp(e,4) = U(udof*inode2);
        Disp(e,5) = U(udof*(inode3-1)+1);
        Disp(e,6) = U(udof*inode3);
        
        u = Disp(e,:)';
        [eps,sig] = CST_El_Str(ipstrn,xy,u,h,Y,nu);
        
        % Store element strains
        Eps(e,:) = eps;
        
        % Store element stresses
        Sig(e,:) = sig;
        
    end
    
    % Computing Strain concentration factor
    
    sig_nom= 1/0.75;
    
    for i=1:size(hol_el)
        sig_max=max(Sig(hol_el(i),:,:));
    end
    SCF(no)=mean(sig_max)/sig_nom;
    
    
    PE(no,1)=0.5*U'*K*U;
    
    PE(no,2)=3/N;
    if (no==1)
        str=sprintf('Original plate vs deformed plate using constant strain triangle elements for %d elements',E);
        figure;
        Plot_deformation;
        title(str);
        xlabel('\leftarrow 2L \rightarrow');
        ylabel('\leftarrow 2H \rightarrow');
        max(U)
    end
%     E
    clearvars -except nod1 nod2 nod3 nod4 nod5 Elm1 Elm2 Elm3 Elm4 Elm5 PE SCF;
end

for i=1:5
    if (PE(i,2)==min(PE(:,2)))
        PE_ex=PE(i,1);
    end
end
PE(:,1)=abs(PE(:,1)-PE_ex)/abs(PE_ex);

% figure;
% plot(log(PE(:,2)),log(PE(:,1)), '-o');
% title('Error in Energy norm');
% xlabel('$log(h)$','Interpreter','latex');
% ylabel('$log(\frac{|U_{FE}-U_{EX}|}{|U_{EX}|})$','Interpreter','latex');
% axis square;
%
% Disp
% Eps
% Sig
