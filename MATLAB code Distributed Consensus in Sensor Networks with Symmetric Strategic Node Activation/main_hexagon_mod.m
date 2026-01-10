%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MATLAB code for the work "Distributed Consensus in Wireless Sensor
% Networks with Strategic Node Activation"
% 
% Sameed Ahmed, Leonardo Badia, Giorgio Gnecco, Daniela Selvi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

clear all
close all
clc

set(0,'defaultAxesFontSize',14);
set(0,'defaultLegendFontSize',14);
set(0,'defaultTextFontSize',16);

tic

%n_rows=3;
%n_columns=3;
%nodes=[1:n_rows*n_columns];

%N_nodes=size(nodes,2);

%W_full=zeros(N_nodes,N_nodes);

%for i1=1:n_rows
%    for j1=1:n_columns
%        for i2=1:n_rows
%            for j2=1:n_columns
%                if (i1==i2 && abs(j1-j2)==1) || (j1==j2 && abs(i1-i2)==1) 
%                    W_full((i1-1)*n_columns+j1,(i2-1)*n_columns+j2)=1;
%                end
%            end
%        end
%    end
%end

N_nodes=19;

W_full=zeros(N_nodes,N_nodes);

list{1}=[2 3 4 5 6 7];
list{2}=[1 7 8 9 10 3];
list{3}=[1 2 10 11 12 4];
list{4}=[1 3 12 13 14 5];
list{5}=[1 4 14 15 16 6];
list{6}=[1 5 16 17 18 7];
list{7}=[1 6 18 19 8 2];
list{8}=[2 7 19 9];
list{9}=[2 8 10];
list{10}=[2 9 11 3];
list{11}=[3 10 12];
list{12}=[3 11 13 4];
list{13}=[4 12 14];
list{14}=[4 13 15 5];
list{15}=[5 14 16];
list{16}=[5 15 17 6];
list{17}=[6 16 18];
list{18}=[6 17 19 7];
list{19}=[7 18 8];

for i=1:size(list,2)
    for j=1:size(list{i},2)
        W_full(i,list{i}(j))=1;
    end
end

figure(1)
G=digraph(W_full','omitselfloops');
%g=plot(G,'Layout','subspace')
g=plot(G,'NodeFontSize',14)

%g.XData = [1 2 3 1 2 3 1 2 3];
%g.YData = [-1 -1 -1 0 0 0 1 1 1];

%G.Edges.Weight(1:5)=4;
%g.LineWidth = G.Edges.Weight;

%highlight(g,1,'MarkerSize',8,'NodeColor','red')
%highlight(g,2,'MarkerSize',8,'NodeColor','red')
%highlight(g,1,2,'EdgeColor','red')
%highlight(g,1,4,'EdgeColor','red')
%highlight(g,2,1,'EdgeColor','red')
%highlight(g,2,3,'EdgeColor','red')
%highlight(g,2,5,'EdgeColor','red')

%camroll(-5)
%title('network topology and example of activation','interpreter','latex')

set(gcf, 'Position', get(0, 'Screensize'));
%saveas(gcf,'figure1_hexagon_mod.pdf');
%print('figure1_hexagon_mod','-dpdf','-bestfit');

%there exists an arborescence, with node 1 as root

OutDegree_full=zeros(N_nodes,N_nodes);
for k=1:N_nodes
    OutDegree_full(k,k)=sum(W_full(:,k));
end

InDegree_full=zeros(N_nodes,N_nodes);
for k=1:N_nodes
    InDegree_full(k,k)=sum(W_full(k,:));
end

L_full=InDegree_full-W_full;

eig(L_full);

W_full_normalized_active=W_full;
for k=1:N_nodes
    if OutDegree_full(k,k)>0
        W_full_normalized_active(:,k)=W_full(:,k)/OutDegree_full(k,k);
    end
end

for h=1:N_nodes
    for k=1:N_nodes
        if (h==11 || k==11)
            W_full_normalized_active(h,k)= W_full_normalized_active(h,k)/2;
        end
    end
end

F=0.1;
W_full_normalized_inactive=F*W_full_normalized_active;

%%%%

W=W_full_normalized_inactive;
InDegree=zeros(N_nodes);
L=zeros(N_nodes);

value=zeros(1,2^N_nodes);
N_active_nodes=zeros(1,2^N_nodes);

old_code=zeros(1,N_nodes);

for k=0:2^N_nodes-1
    code=dec2gc_mod(k,N_nodes);
    N_active_nodes(k+1)=sum(code);
    l=find(code~=old_code);
    old_code=code;
    if code(l)
        W(:,l)=W_full_normalized_active(:,l);
    else
        W(:,l)=W_full_normalized_inactive(:,l);
    end
    for l=1:N_nodes
        InDegree(l,l)=sum(W(l,:));
    end

    L=InDegree-W;

    temp=sort(eig(L),'ComparisonMethod','real');
    value(k+1)=-1/real(temp(2));
end

for i=0:N_nodes
    optimal_value(i+1)=max(value(find(N_active_nodes==i)));
end

tolerance=10^(-10);
optimal_configuration=dec2gc_mod(find(value>=max(optimal_value)-tolerance)-1,N_nodes)

%%% plots

figure(2)
cc=hsv(N_nodes+3);

hold on
c_vector=[0:0.01:30];
%c_vector=[0];
%c_vector=[30];
legend_string=[];
for i=0:N_nodes
    y_1_2(i+1,:)=-optimal_value(i+1)+c_vector.*i;
    y_3(i+1,:)=-optimal_value(i+1)+c_vector*(i>0);
    plot(c_vector,y_1_2(i+1,:),'color',cc(i+1,:),'Linewidth',1);
    %legend_string=strvcat(legend_string,sprintf('no. active nodes=%.1d',i));
    legend_string=strvcat(legend_string,sprintf('$n_\\textrm{active}$=%.1d',i));
end
for k=1:size(y_1_2,2)
    best_y_1_2(k)=min(y_1_2(:,k));
    best_y_3(k)=min(y_3(:,k));
end
plot(c_vector,best_y_1_2,'--black','Linewidth',2);
legend_string=strvcat(legend_string,sprintf('lower envelope',i));

legend(legend_string,'FontSize',6,'Location','Northwest','interpreter','latex');

xlabel('unit activation cost $k$','interpreter','latex','FontSize',16)
%ylabel('minimum penalty 2 (convergence time + total activation cost over the nodes)','interpreter','latex')
ylabel('minimum of $c_\textrm{global}$','interpreter','latex','FontSize',16)

%title('penalty values of optimal configurations for given numbers of active nodes','interpreter','latex')

grid on

set(gcf, 'Position', get(0, 'Screensize'));
%saveas(gcf,'figure2_hexagon_mod.pdf');
%print('figure2_hexagon_mod','-dpdf','-bestfit');

%%% Nash

tic

InDegree_active=zeros(N_nodes,N_nodes);
InDegree_inactive=zeros(N_nodes,N_nodes);
W_active=zeros(N_nodes,N_nodes);
W_inactive=zeros(N_nodes,N_nodes);
L_active=zeros(N_nodes,N_nodes);
L_inactive=zeros(N_nodes,N_nodes);

W=W_full_normalized_inactive;

old_code=zeros(1,N_nodes);

for k=0:2^N_nodes-1
    code=dec2gc_mod(k,N_nodes);
    %N_active_nodes(k+1)=sum(code);
    l=find(code~=old_code);
    old_code=code;
    if code(l)
        W(:,l)=W_full_normalized_active(:,l);
    else
        W(:,l)=W_full_normalized_inactive(:,l);
    end
    for i=1:N_nodes
        W_active=W;
        W_inactive=W;
        W_active(:,i)=W_full_normalized_active(:,i);
        W_inactive(:,i)=W_full_normalized_inactive(:,i);
        for l=1:N_nodes
            InDegree_active(l,l)=sum(W_active(l,:));
            InDegree_inactive(l,l)=sum(W_inactive(l,:));
        end

        L_active=InDegree_active-W_active;
        L_inactive=InDegree_inactive-W_inactive;

        temp_active=sort(eig(L_active),'ComparisonMethod','real');
        temp_inactive=sort(eig(L_inactive),'ComparisonMethod','real');

        value_active(i,k+1)=-1/real(temp_active(2));
        value_inactive(i,k+1)=-1/real(temp_inactive(2));

    end
end

%%%%

difference=value_active-value_inactive;

for h=1:size(c_vector,2)
    c=c_vector(h);
    flag_matrix_1=difference>c;
    flag_matrix_indifference=difference==c;
    flag_matrix_0=difference<c;

    Nash_code{h}=[];
    indices{h}=[];
    for k=0:2^N_nodes-1
        code=dec2gc_mod(k,N_nodes);

        Nash_flag=1;
        for l=1:N_nodes
            if code(l) && flag_matrix_0(l,k+1)==1
                Nash_flag=0;
            elseif ~code(l) && flag_matrix_1(l,k+1)==1
                Nash_flag=0;
            end
        end
        if Nash_flag==1
            Nash_code{h}=[Nash_code{h};code];
            indices{h}=[indices{h};k];
        end
        objective_1_Nash{h}=-value(indices{h}+1);
        objective_2_Nash{h}=objective_1_Nash{h}+c*sum(Nash_code{h},2)';
        objective_3_Nash{h}=objective_1_Nash{h}+(c*(sum(Nash_code{h},2)'>0));
    end

%display(c,'cost if active')

display(Nash_code{h},'pure Nash equilibria')

%display(value(find(N_active_nodes==N_nodes)),'second smallest eigenvalue for full activity')

%display(optimal_value,'maximum second smallest eigenvalue for variable activity (from 0 to all active nodes)')

%display(value(indices+1),'second smallest eigenvalue at pure Nash equilibria')

end

%%% symmetric Nash equilibrium

p_train=[0:1/N_nodes:1];
%p_train=rand(1,N_nodes+1);

expected_utility_inactive=zeros(size(p_train,2),N_nodes);
expected_utility_active=zeros(size(p_train,2),N_nodes);

for p_index=1:size(p_train,2) %keep N_nodes+1 values for p, equally spaced (with the idea of fitting a polynomial of degree n_nodes)
    for i=1:N_nodes %here one supposes that one of the N_nodes transmitters decides to be inactive, and is assigned randomly
        % to one of the N_nodes

        symmetric_Nash_old_code=zeros(1,N_nodes-1);

        expected_utility_inactive(p_index,i)=0;
        W=W_full_normalized_inactive;
        W(:,i)=W_full_normalized_inactive(:,i);
        for k=0:2^(N_nodes-1)-1
            symmetric_Nash_code=dec2gc_mod(k,N_nodes-1);%the first Gray codeword is made only by zeros;
            N_active_nodes_symmetric_Nash(k+1)=sum(symmetric_Nash_code);
            probability_configuration(k+1)=p_train(p_index)^N_active_nodes_symmetric_Nash(k+1)*(1-p_train(p_index))^(N_nodes-1-N_active_nodes_symmetric_Nash(k+1));
            l=find(symmetric_Nash_code~=symmetric_Nash_old_code);
            symmetric_Nash_old_code=symmetric_Nash_code;
            if l<i
                if symmetric_Nash_code(l)
                    W(:,l)=W_full_normalized_active(:,l);
                else
                    W(:,l)=W_full_normalized_inactive(:,l);
                end
            elseif l>=i
                if symmetric_Nash_code(l)
                    W(:,l+1)=W_full_normalized_active(:,l+1);
                else
                    W(:,l+1)=W_full_normalized_inactive(:,l+1);
                end
            end

            for l=1:N_nodes
                InDegree(l,l)=sum(W(l,:));
            end

            L=InDegree-W;

            temp=sort(eig(L),'ComparisonMethod','real');
            value_symmetric_Nash_inactive(i,k+1)=-1/real(temp(2));

            expected_utility_inactive(p_index,i)=expected_utility_inactive(p_index,i)+probability_configuration(k+1)*value_symmetric_Nash_inactive(i,k+1);
        end
          expected_utility_inactive(p_index,i)=expected_utility_inactive(p_index,i)/N_nodes;
    end

    for i=1:N_nodes %here one supposes that one of the N_nodes transmitters decides to be active, and is assigned randomly
        % to one of the N_nodes

        symmetric_Nash_old_code=zeros(1,N_nodes-1);

        expected_utility_active(p_index,i)=0;
        W=W_full_normalized_inactive;
        W(:,i)=W_full_normalized_active(:,i);
        for k=0:2^(N_nodes-1)-1
            symmetric_Nash_code=dec2gc_mod(k,N_nodes-1);%the first Gray codeword is made only by zeros;
            N_active_nodes_symmetric_Nash(k+1)=sum(symmetric_Nash_code);
            probability_configuration(k+1)=p_train(p_index)^N_active_nodes_symmetric_Nash(k+1)*(1-p_train(p_index))^(N_nodes-1-N_active_nodes_symmetric_Nash(k+1));
            l=find(symmetric_Nash_code~=symmetric_Nash_old_code);
            symmetric_Nash_old_code=symmetric_Nash_code;
            if l<i
                if symmetric_Nash_code(l)
                    W(:,l)=W_full_normalized_active(:,l);
                else
                    W(:,l)=W_full_normalized_inactive(:,l);
                end
            elseif l>=i
                if symmetric_Nash_code(l)
                    W(:,l+1)=W_full_normalized_active(:,l+1);
                else
                    W(:,l+1)=W_full_normalized_inactive(:,l+1);
                end
            end

            for l=1:N_nodes
                InDegree(l,l)=sum(W(l,:));
            end

            L=InDegree-W;

            temp=sort(eig(L),'ComparisonMethod','real');
            value_symmetric_Nash_active(i,k+1)=-1/real(temp(2));

            expected_utility_active(p_index,i)=expected_utility_active(p_index,i)+probability_configuration(k+1)*value_symmetric_Nash_active(i,k+1);
        end
        expected_utility_active(p_index,i)=expected_utility_active(p_index,i)/N_nodes;
    end

    additional_average_on_the_positions_expected_utility_inactive(p_index)=mean(expected_utility_inactive(p_index,:));
    additional_average_on_the_positions_expected_utility_active(p_index)=mean(expected_utility_active(p_index,:));

end

marginal_utility=additional_average_on_the_positions_expected_utility_active-additional_average_on_the_positions_expected_utility_inactive;

polynomial_coefficients=polyfit(p_train,marginal_utility,N_nodes);

%%% plots

figure(3)

p_test=[0:0.0001:1];
for p_index=1:size(p_test,2)
    fitted_polynomial(p_index)=polyval(polynomial_coefficients,p_test(p_index));
end

hold on

plot(p_test,fitted_polynomial,'b','Linewidth',1);

c=0.5;
plot(p_test,c*ones(size(p_test)),'r','Linewidth',1);

%legend('marginal expected utility of activation','marginal cost of activation','FontSize',14,'interpreter','latex');
legend(['expected utility of one-node activation' newline 'for zero unit activation cost ($\beta(p)-\alpha(p)$)'],'unit activation cost ($k$)','FontSize',14,'Location','South','interpreter','latex');

xlabel('common probability of activation $p$','interpreter','latex','FontSize',16)
ylabel('$\beta-\alpha$, $k$','interpreter','latex','FontSize',16)

%title('symmetric mixed strategy Nash equilibria for a given unit activation cost $k$','interpreter','latex')

grid on

set(gcf, 'Position', get(0, 'Screensize'));
%saveas(gcf,'figure3_hexagon_mod.pdf');
%print('figure3_hexagon_mod','-dpdf','-bestfit');

%%% optimal symmetric mixed strategy

p_train=[0:1/N_nodes:1];
%p_train=rand(1,N_nodes+1);

expected_utility=zeros(size(p_train,2),1);

for p_index=1:size(p_train,2)
    optimal_symmetric_mixed_strategy_old_code=zeros(1,N_nodes);
    expected_utility(p_index)=0;
    W=W_full_normalized_inactive;
    for k=0:2^N_nodes-1
        optimal_symmetric_mixed_strategy_code=dec2gc_mod(k,N_nodes);%the first Gray codeword is made only by zeros;
        N_active_nodes_optimal_symmetric_mixed_strategy(k+1)=sum(optimal_symmetric_mixed_strategy_code);
        probability_configuration(k+1)=p_train(p_index)^N_active_nodes_optimal_symmetric_mixed_strategy(k+1)*(1-p_train(p_index))^(N_nodes-N_active_nodes_optimal_symmetric_mixed_strategy(k+1));
        l=find(optimal_symmetric_mixed_strategy_code~=optimal_symmetric_mixed_strategy_old_code);
        optimal_symmetric_mixed_strategy_old_code=optimal_symmetric_mixed_strategy_code;
        if optimal_symmetric_mixed_strategy_code(l)
            W(:,l)=W_full_normalized_active(:,l);
        else
            W(:,l)=W_full_normalized_inactive(:,l);
        end

        for l=1:N_nodes
            InDegree(l,l)=sum(W(l,:));
        end

        L=InDegree-W;

        temp=sort(eig(L),'ComparisonMethod','real');
        value_optimal_symmetric_mixed_strategy(k+1)=-1/real(temp(2));

        expected_utility(p_index,1)=expected_utility(p_index)+probability_configuration(k+1)*value_optimal_symmetric_mixed_strategy(k+1);
    end
end

polynomial_coefficients_2=polyfit(p_train,expected_utility',N_nodes);

%%% plots

figure(4)

%p_test=[0:0.0001:1];
for p_index=1:size(p_test,2)
    fitted_polynomial_2(p_index)=polyval(polynomial_coefficients_2,p_test(p_index));
end

hold on

plot(p_test,-fitted_polynomial_2,'b','Linewidth',1);

c=4;
total_cost_test=c*p_test*N_nodes;

plot(p_test,total_cost_test,'r','Linewidth',1);

plot(p_test,total_cost_test-fitted_polynomial_2,'g','Linewidth',1);

legend('expected global cost for zero unit activation cost $\delta(p)$','expected total cost of activation $k n p$','sum $\delta(p) + k n p$','Location','Northeast','FontSize',14,'interpreter','latex');

ylim([0,200]);

xlabel('common probability of activation $p$','interpreter','latex','FontSize',16)
ylabel('$\delta$, $k n p$, $\delta+k n p$','interpreter','latex','FontSize',16)

%title('optimal symmetric mixed strategies for a given unit cost of activation','interpreter','latex')

grid on

set(gcf, 'Position', get(0, 'Screensize'));
%saveas(gcf,'figure4_hexagon_mod.pdf');
%print('figure4_hexagon_mod','-dpdf','-bestfit');

%%% symmetric Nash equilibria obtained by varying c

derivative=diff(fitted_polynomial);

counter=1;
if derivative(1)<0
    sign_old_derivative=-1;
else sign_old_derivative=1;
end
intervals{counter}=[1];
intervals_p{counter}=[0];

for k=1:size(fitted_polynomial,2)-1
    if derivative(k)*sign_old_derivative<0
        intervals{counter}=[intervals{counter},k];
        intervals_p{counter}=[intervals_p{counter},p_test(k)];
        counter=counter+1;
        sign_old_derivative=-sign_old_derivative;
        intervals{counter}=k;
        intervals_p{counter}=p_test(k);
    end
end
intervals{counter}=[intervals{counter},size(fitted_polynomial,2)];
intervals_p{counter}=[intervals_p{counter},1];

%%%%

for k=1:size(c_vector,2)
    c=c_vector(k);
    symmetric_Nash{k}=[];
    if fitted_polynomial(1)<=c
        symmetric_Nash{k}=[symmetric_Nash{k},0];
    end
    for l=1:size(intervals,2)
        if (fitted_polynomial(intervals{l}(1))-c)*(fitted_polynomial(intervals{l}(2))-c)<=0
            [temp_min,index]=min(abs(fitted_polynomial(intervals{l}(1):intervals{l}(2))-c));
            if index ~=1
                symmetric_Nash{k}=[symmetric_Nash{k},p_test(index+intervals{l}(1)-1)];
            end
        end
        if l==size(intervals,2) && fitted_polynomial(end)>=c
            symmetric_Nash{k}=[symmetric_Nash{k},1];
        end
    end
    objective_1_symmetric_Nash{k}=-polyval(polynomial_coefficients_2,symmetric_Nash{k});
    objective_2_symmetric_Nash{k}=objective_1_symmetric_Nash{k}+c*symmetric_Nash{k}*N_nodes;
    objective_3_symmetric_Nash{k}=objective_1_symmetric_Nash{k}+c*(1-(1-symmetric_Nash{k}).^N_nodes);
end

%%% symmetric optimal mixed strategy obtained by varying c

derivative_2=diff(fitted_polynomial_2);

counter=1;
if derivative_2(1)<0
    sign_old_derivative_2=-1;
else sign_old_derivative_2=1;
end
intervals_2{counter}=[1];
intervals_2_p{counter}=[0];

for k=1:size(fitted_polynomial_2,2)-1
    if derivative_2(k)*sign_old_derivative_2<0
        intervals_2{counter}=[intervals_2{counter},k];
        intervals_2_p{counter}=[intervals_2_p{counter},p_test(k)];
        counter=counter+1;
        sign_old_derivative_2=-sign_old_derivative_2;
        intervals_2{counter}=k;
        intervals_2_p{counter}=p_test(k);
    end
end
intervals_2{counter}=[intervals_2{counter},size(fitted_polynomial_2,2)];
intervals_2_p{counter}=[intervals_2_p{counter},1];

%%%%

for k=1:size(c_vector,2)
    c=c_vector(k);
    total_cost_test_1=0;
    objective_function_1=-fitted_polynomial_2+total_cost_test_1;
    [temp_min_1,temp_min_index_1]=min(objective_function_1);
    symmetric_optimal_mixed_strategy_1{k}=temp_min_index_1;
    symmetric_optimal_mixed_strategy_p_1{k}=p_test(temp_min_index_1);
    y_objective_1_optimal_symmetric_mixed_strategy(k)=temp_min_1;

    total_cost_test_2=c*p_test*N_nodes;
    objective_function_2=-fitted_polynomial_2+total_cost_test_2;
    [temp_min_2,temp_min_index_2]=min(objective_function_2);
    symmetric_optimal_mixed_strategy_2{k}=temp_min_index_2;
    symmetric_optimal_mixed_strategy_p_2{k}=p_test(temp_min_index_2);
    y_objective_2_optimal_symmetric_mixed_strategy(k)=temp_min_2;

    total_cost_test_3=c*(1-(1-p_test).^N_nodes);
    objective_function_3=-fitted_polynomial_2+total_cost_test_3;
    [temp_min_3,temp_min_index_3]=min(objective_function_3);
    symmetric_optimal_mixed_strategy_3{k}=temp_min_index_3;
    symmetric_optimal_mixed_strategy_p_3{k}=p_test(temp_min_index_3);
    y_objective_3_optimal_symmetric_mixed_strategy(k)=temp_min_3;
end

%%% objective 1

cc=hsv(12);

best_y_objective_1=best_y_1_2(1)*ones(size(best_y_1_2));

for h=1:size(objective_1_Nash,2)
    worst_objective_1_Nash(h)=max(objective_1_Nash{h});
    best_objective_1_Nash(h)=min(objective_1_Nash{h});
end

for h=1:size(objective_1_symmetric_Nash,2)
    worst_objective_1_symmetric_Nash(h)=max(objective_1_symmetric_Nash{h});
    best_objective_1_symmetric_Nash(h)=min(objective_1_symmetric_Nash{h});
end

figure(5)
hold on

tolerance=2*10^(-1);

plot(c_vector,best_objective_1_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(3,:));
plot(c_vector,worst_objective_1_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(5,:));
plot(c_vector,best_y_objective_1,'x--','Linewidth',1,'Markersize',1,'color',cc(1,:));
plot(c_vector,best_objective_1_symmetric_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(7,:));
plot(c_vector,worst_objective_1_symmetric_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(9,:));
plot(c_vector,y_objective_1_optimal_symmetric_mixed_strategy+tolerance,'x--','Linewidth',1,'Markersize',1,'color',cc(11,:));

legend_string={'best NE in pure strategies for LC','worst NE in pure strategies for LC','optimal $n$-tuple of pure strategies for GC'''' ','best symmetric NE in mixed strategies for F-LC','worst symmetric NE in mixed strategies for F-LC','optimal symmetric mixed strategy for F-GC'''' '};
legend(legend_string,'FontSize',14,'Location','Northeast','interpreter','latex');

ylim([0,200]);

xlabel('unit activation cost $k$','interpreter','latex','FontSize',16)
%ylabel('penalty 1 (convergence time)','interpreter','latex')
ylabel('(expected) $c_\textrm{convergence}$','interpreter','latex','FontSize',16)

%title('comparison of penalty values of strategies (case of penalty 1)','interpreter','latex')

grid on

set(gcf, 'Position', get(0, 'Screensize'));
%saveas(gcf,'figure5_hexagon_mod.pdf');
%print('figure5_hexagon_mod','-dpdf','-bestfit');

%%% objective 2

best_y_objective_2=best_y_1_2;

for h=1:size(objective_2_Nash,2)
    worst_objective_2_Nash(h)=max(objective_2_Nash{h});
    best_objective_2_Nash(h)=min(objective_2_Nash{h});
end

for h=1:size(objective_2_symmetric_Nash,2)
    worst_objective_2_symmetric_Nash(h)=max(objective_2_symmetric_Nash{h});
    best_objective_2_symmetric_Nash(h)=min(objective_2_symmetric_Nash{h});
end

figure(6)
hold on

plot(c_vector,best_objective_2_Nash+tolerance,'x--','Linewidth',1,'Markersize',1,'color',cc(3,:));
plot(c_vector,worst_objective_2_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(5,:));
plot(c_vector,best_y_objective_2,'x--','Linewidth',1,'Markersize',1,'color',cc(1,:));
plot(c_vector,best_objective_2_symmetric_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(7,:));
plot(c_vector,worst_objective_2_symmetric_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(9,:));
plot(c_vector,y_objective_2_optimal_symmetric_mixed_strategy,'x--','Linewidth',1,'Markersize',1,'color',cc(11,:));

legend_string={'best NE in pure strategies for LC','worst NE in pure strategies for LC','optimal $n$-tuple of pure strategies for GC','best symmetric NE in mixed strategies for F-LC','worst symmetric NE in mixed strategies for F-LC','optimal symmetric mixed strategy for F-GC'};
legend(legend_string,'FontSize',14,'Location','Northwest','interpreter','latex');

ylim([0,550]);

xlabel('unit activation cost $k$','interpreter','latex','FontSize',16)
%ylabel('penalty 2 (convergence time + total activation cost over the players)','interpreter','latex')
ylabel('(expected) $c_\textrm{global}$','interpreter','latex','FontSize',16)

%title('comparison of penalty values of strategies (case of penalty 2)','interpreter','latex')

grid on

set(gcf, 'Position', get(0, 'Screensize'));
%saveas(gcf,'figure6_hexagon_mod.pdf');
%print('figure6_hexagon_mod','-dpdf','-bestfit');

%%% objective 3

best_y_objective_3=best_y_3;

for h=1:size(objective_3_Nash,2)
    worst_objective_3_Nash(h)=max(objective_3_Nash{h});
    best_objective_3_Nash(h)=min(objective_3_Nash{h});
end

for h=1:size(objective_3_symmetric_Nash,2)
    worst_objective_3_symmetric_Nash(h)=max(objective_3_symmetric_Nash{h});
    best_objective_3_symmetric_Nash(h)=min(objective_3_symmetric_Nash{h});
end

figure(7)
hold on

plot(c_vector,best_objective_3_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(3,:));
plot(c_vector,worst_objective_3_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(5,:));
plot(c_vector,best_y_objective_3,'x--','Linewidth',1,'Markersize',1,'color',cc(1,:));
plot(c_vector,best_objective_3_symmetric_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(7,:));
plot(c_vector,worst_objective_3_symmetric_Nash,'x--','Linewidth',1,'Markersize',1,'color',cc(9,:));
plot(c_vector,y_objective_3_optimal_symmetric_mixed_strategy+tolerance,'x--','Linewidth',1,'Markersize',1,'color',cc(11,:));

legend_string={'best NE in pure strategies for LC','worst NE in pure strategies for LC','optimal $n$-tuple of pure strategies for GC'' ','best symmetric NE in mixed strategies for F-LC','worst symmetric NE in mixed strategies for F-LC','optimal symmetric mixed strategy for F-GC'' '};
legend(legend_string,'FontSize',14,'Location','Northwest','interpreter','latex');

ylim([0,250]);

xlabel('unit activation cost $k$','interpreter','latex','FontSize',16)
%ylabel('penalty 3 (convergence time + maximum activation cost over the players)','interpreter','latex')
ylabel('(expected) $c_\textrm{individual,max}$','interpreter','latex','FontSize',16)

%title('comparison of penalty values of strategies (case of penalty 3)','interpreter','latex')

grid on

set(gcf, 'Position', get(0, 'Screensize'));
%saveas(gcf,'figure7_hexagon_mod.pdf');
%print('figure7_hexagon_mod','-dpdf','-bestfit');

toc

save workspace_hexagon_mod
%load workspace_hexagon_mod