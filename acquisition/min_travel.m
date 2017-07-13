

function out = min_travel(arr, n_iter, temperature)
    %simulated annealing
    cooling_rate=0.0001;
    cost = 0;
    for i=1:length(arr)-1
        %cost is the rotation time spent
        cost = cost+max(dist(arr(i+1,1),arr(i,1)),dist(arr(i+1,2),arr(i,2)));
    end
    disp(['Starting cost ', num2str(cost)])
    for n=1:n_iter
        temperature = temperature*exp(-cooling_rate*(n-1));
        %swap random rows and check cost
        rows=int16((length(arr)-1)*rand(1,2)+1);
        arr_new = arr;
        tmp=arr_new(rows(1),:);
        arr_new(rows(1),:)=arr_new(rows(2),:);
        arr_new(rows(2),:)=tmp;
            
        cost_new = 0;
        for i=1:length(arr_new)-1
            cost_new = cost_new+max(dist(arr_new(i+1,1),arr_new(i,1)),dist(arr_new(i+1,2),arr_new(i,2)));
        end
        diff = double(cost_new-cost);
        if diff<0
            arr = arr_new;
            cost=cost_new;
        elseif exp(-diff/temperature)>rand()
            arr = arr_new;
            cost=cost_new;
        end
    end
    disp(['End cost ', num2str(cost_new)])
    out = arr;
end

function d = dist(p1,p2)
    p1=int16(p1);
    p2=int16(p2);
    d1=abs(p1-p2);
    d2=abs(p1+360-p2);
    d3=abs(p2-p1-360);  
    d=min([d1,d2,d3]);
end