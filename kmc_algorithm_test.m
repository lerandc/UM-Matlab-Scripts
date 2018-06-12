clearvars
close all

field = ones(50); %field for possible events

vacancies = randi([1,2500],100,1);

field(vacancies) = 0;
time = 0;
time_vec = zeros(1001,1);
pos_vec = zeros(size(vacancies,1),2,1001);

for i = 1:size(vacancies,1)
    [pos_vec(i,2,1), pos_vec(i,1,1)] = ind2sub([50 50],vacancies(i));
    plot(pos_vec(i,1,1),-1*pos_vec(i,2,1),'b.','MarkerSize',14)
    hold on
end
axis([0 50 -50 0])

    rate_up = 3;
    rate_down = 3;
    rate_left = 3;
    rate_right = 3;

for j = 1:1000


    rates_vec = zeros(4*size(vacancies,1),1); %4n-3 = up, 4n-2 = down, 4n-1 = left, 4n = right for nth particle

    rates_vec(~mod(1:size(rates_vec,1),4)) = rate_right;
    rates_vec(~(mod(1:size(rates_vec,1),4)-3)) = rate_left;
    rates_vec(~(mod(1:size(rates_vec,1),4)-2)) = rate_down;
    rates_vec(~(mod(1:size(rates_vec,1),4)-1)) = rate_up;
    
    tot_rate = sum(rates_vec);

    partitions = [0; cumsum(rates_vec)]./tot_rate;
    
    num = rand;
    event = find(sort([partitions; num]) == num)-1;

    %action sequence
    particle = ceil(event/4);
    action = mod(event,4);

    [p_y,p_x] = ind2sub([50 50], vacancies(particle));

    field(p_y,p_x) = 1;

    switch action
        case 1 %up
            if p_y-1>0
                field(p_y-1,p_x) = 0;
            else
                field(p_y,p_x) = 0;
            end
        case 2 %down
            if p_y+1<size(field,1)+1
                field(p_y+1,p_x) = 0;
            else
                field(p_y,p_x) = 0;
            end
        case 3 %left
            if p_x-1>0
                field(p_y,p_x-1) = 0;
            else
                field(p_y,p_x) = 0;
            end
        case 0 %right
            if p_x+1<size(field,2)+1
                field(p_y,p_x+1) = 0;
            else
                field(p_y,p_x) = 0;
            end
    end

    d_time = (1/tot_rate)*log(1/rand);
    time = time+d_time;
    time_vec(j+1) = time_vec(j)+d_time;
    vacancies = find(field(:) == 0);
    
    [pos_vec(1:size(vacancies,1),1,j+1), pos_vec(1:size(vacancies,1),2,j+1)] = ind2sub([50 50],vacancies);

    clf

    plot(pos_vec(:,2,j+1),-1*pos_vec(:,1,j+1),'b.','MarkerSize',14)
    
    axis([0 50 -50 0])
    drawnow
end
