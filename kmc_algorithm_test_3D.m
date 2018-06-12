clearvars
close all

field = ones(12,12,12); %field for possible events

vacancies = randi([1,1728],10,1);

field(vacancies) = 0;
time = 0;
time_vec = zeros(1001,1);
pos_vec = zeros(size(vacancies,1),3,1001);

for i = 1:size(vacancies,1)
    [pos_vec(i,2,1), pos_vec(i,1,1),pos_vec(i,3,1)] = ind2sub([12 12 12],vacancies(i));
    plot3(pos_vec(i,1,1),pos_vec(i,2,1),pos_vec(i,3,1),'.','MarkerSize',14)
    hold on
end
axis([1 12 1 12 1 12])
grid on
grid minor

    rate_up = 3;
    rate_down = 3;
    rate_left = 3;
    rate_right = 3;
    rate_forward = 3;
    rate_back = 3;

for j = 1:10000


    rates_vec = zeros(6*size(vacancies,1),1); %4n-3 = up, 4n-2 = down, 4n-1 = left, 4n = right for nth particle

    rates_vec(~mod(1:size(rates_vec,1),6)) = rate_right;
    rates_vec(~(mod(1:size(rates_vec,1),6)-5)) = rate_left;
    rates_vec(~(mod(1:size(rates_vec,1),6)-4)) = rate_down;
    rates_vec(~(mod(1:size(rates_vec,1),6)-3)) = rate_up;
    rates_vec(~(mod(1:size(rates_vec,1),6)-2)) = rate_forward;
    rates_vec(~(mod(1:size(rates_vec,1),6)-1)) = rate_back;
    
    tot_rate = sum(rates_vec);

    partitions = [0; cumsum(rates_vec)]./tot_rate;
    
    num = rand;
    event = find(sort([partitions; num]) == num)-1;

    %action sequence
    particle = ceil(event/6);
    action = mod(event,6);

    [p_y,p_x,p_z] = ind2sub([12,12,12], vacancies(particle));

    field(p_y,p_x,p_z) = 1;

    switch action
        case 3 %up
            if p_y-1>0
                field(p_y-1,p_x,p_z) = 0;
            else
                field(p_y,p_x,p_z) = 0;
            end
        case 4 %down
            if p_y+1<size(field,1)+1
                field(p_y+1,p_x,p_z) = 0;
            else
                field(p_y,p_x,p_z) = 0;
            end
        case 5 %left
            if p_x-1>0
                field(p_y,p_x-1,p_z) = 0;
            else
                field(p_y,p_x,p_z) = 0;
            end
        case 0 %right
            if p_x+1<size(field,2)+1
                field(p_y,p_x+1,p_z) = 0;
            else
                field(p_y,p_x,p_z) = 0;
            end
        case 1 %back
            if p_z-1>0
                field(p_y,p_x,p_z-1) = 0;
            else
                field(p_y,p_x,p_z) = 0;
            end
        case 2 %forward
            if p_z+1<size(field,3)+1
                field(p_y,p_x,p_z+1) = 0;
            else
                field(p_y,p_x,p_z) = 0;
            end
    end

    d_time = (1/tot_rate)*log(1/rand);
    time = time+d_time;
    time_vec(j+1) = time_vec(j)+d_time;
    vacancies = find(field(:) == 0);
    
    [pos_vec(1:size(vacancies,1),1,j+1), pos_vec(1:size(vacancies,1),2,j+1), pos_vec(1:size(vacancies,1),3,j+1) ] =...
        ind2sub([12,12,12],vacancies);

    clf

    plot3(pos_vec(:,2,j+1),pos_vec(:,1,j+1),pos_vec(:,3,j+1),'.','MarkerSize',14)
    
    axis([1 12 1 12 1 12])
    grid on
    grid minor
    drawnow
end
