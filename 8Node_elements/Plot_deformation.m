%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     PLOTTING THE DEFORMATION                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  PLotting the initial position nodes

plot(Nodes(:,2),Nodes(:,3),'ob'); axis equal; axis tight; hold on;
plot(-Nodes(:,2),-Nodes(:,3),'ob'); axis equal; axis tight; hold on;
plot(-Nodes(:,2),Nodes(:,3),'ob'); axis equal; axis tight; hold on;
plot(Nodes(:,2),-Nodes(:,3),'ob'); axis equal; axis tight; hold on;

% Finding the final position of the nodes
j=1;
for i=1:2:size(U)
    n_disp(j,1)=U(i);
    n_disp(j,2)=U(i+1);
    j=j+1;
end
n_final(:,1)=Nodes(:,2)+n_disp(:,1);
n_final(:,2)=Nodes(:,3)+n_disp(:,2);


% Plotting the final positions of the nodes

plot(n_final(:,1),n_final(:,2),'*r'); axis equal; axis tight; hold on;
plot(-n_final(:,1),n_final(:,2),'*r'); axis equal; axis tight; hold on;
plot(n_final(:,1),-n_final(:,2),'*r'); axis equal; axis tight; hold on;
plot(-n_final(:,1),-n_final(:,2),'*r'); axis equal; axis tight; hold on;
% set(gca,'color','black')