
d = 3;
K = 6;
T = 200000; % time in days
h = 10;
N = T/h;
q = zeros(N,K,d); % iteration, planet, coordinate
p = zeros(N,K,d); % iteration, planet, coordinate
m = zeros(K,1);
G = 2.95912208286*10^(-4);

% Sun
m(1) = 1.00000597682;
p(1,1,:) = 0;
q(1,1,:) = 0;


% Jupiter
m(2) = 0.000954786104043;
p(1,2,:) = [0.00565429, -0.00412490, -0.00190589].*m(2); % Jupiter
q(1,2,:) = [-3.5023653, -3.8169847, -1.5507963];


% Saturn
m(3) = 0.000285583733151;
p(1,3,:) = [0.00168318, 0.00483525, 0.001922462].*m(3); % Saturn
q(1,3,:) = [9.0755314, -3.0458353, -1.6483708];

% Uranus
m(4) = 0.0000437273164546;
p(1,4,:) = [0.00354178, 0.00137102, 0.00055029]*m(4);
q(1,4,:) = [8.3101420, -16.2901086, -7.2521278];

% Neptune
m(5) = 0.0000517759138449;
p(1,5,:) = [0.00288930, 0.00114527, 0.00039677]*m(5);
q(1,5,:) = [11.4707666, -25.7294829, -10.8169456];

% Pluto
m(6) = 1/(1.3 * 10^8);
p(1,6,:) = [0.00276725, -0.00170702, -0.00136504]*m(6);
q(1,6,:) = [-15.5387357, -25.2225594, -3.1902382];



for n = 2:N
    for k = 1:K
        term = 0;
        q_re = reshape(q(n-1,:,:),[K,d]);
        for i = 1:K
            if i == k
                continue
            end
            term1 = (m(k)*m(i)/(norm(q_re(k,:)-q_re(i,:))^3)) * (q_re(k,:)-q_re(i,:));
            term = term + term1;
        end
        term = G.*reshape(term,[1,1,d]);
        p(n,k,:) = p(n-1,k,:) - h*term;
        %p(n,k,:) = p(n-1,k,:) - h*terms(G,K,k,d,m,n,q);
        q(n,k,:) = q(n-1,k,:) + h*p(n-1,k,:)./m(k);
    end
end

sun = [q(:,1,1), q(:,1,2), q(:,1,3)];
jupiter = [q(:,2,1), q(:,2,2), q(:,2,3)];
saturn = [q(:,3,1), q(:,3,2), q(:,3,3)];
uranus = [q(:,4,1), q(:,4,2), q(:,4,3)];
neptune = [q(:,5,1), q(:,5,2), q(:,5,3)];
pluto = [q(:,6,1), q(:,6,2), q(:,6,3)];

figure(1)
plot3(sun(:,1), sun(:,2), sun(:,3),'o')
hold on
plot3(jupiter(:,1), jupiter(:,2z), jupiter(:,3))
hold on
plot3(saturn(:,1), saturn(:,2), saturn(:,3))
hold on
plot3(uranus(:,1), uranus(:,2), uranus(:,3))
hold on
plot3(neptune(:,1), neptune(:,2), neptune(:,3))
hold on 
plot3(pluto(:,1), pluto(:,2), pluto(:,3))

function term = terms(G,K,k,d,m,n,q)
    term = 0;
    q_re = reshape(q(n-1,:,:),[K,d]);
    for i = 1:K
        if i == k
            continue
        end
        term = term + (m(i)/(norm(q_re(k,:)-q_re(i,:))^3)) * (q_re(k,:)-q_re(i,:));
    end
    term = G.*reshape(term,[1,1,d]);
end
