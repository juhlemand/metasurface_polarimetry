%meas_points = int32(360*rand(100,2));
%disp('Annealing travel path...')
%meas_points = min_travel2(meas_points, 300000, 10);

function out = min_travel(arr, n_iter, temperature)
    %simulated annealing
    cooling_rate=0.00002;
    
    %cost is the rotation time spent
    c = cost(arr);
    disp(['Starting cost ', num2str(c)])
    for n=1:n_iter
        temperature = temperature*exp(-cooling_rate*(n-1));
        %swap random rows and check cost
        rows=uint32((length(arr)-1)*rand(1,2)+1);
        arr_new = arr;
        tmp=arr_new(rows(1),:);
        arr_new(rows(1),:)=arr_new(rows(2),:);
        arr_new(rows(2),:)=tmp;
            
        c_new = cost(arr_new);
        delta = double(c_new-c);
        if delta<0
            arr = arr_new;
            c=c_new;
        elseif exp(-delta/temperature)>rand()
            arr = arr_new;
            c=c_new;
        end
    end
    disp(['End cost ', num2str(c_new)])
    out = arr;
end

function c = cost(arr)
    %distance on circular coordinates
    d1=abs(diff(arr));
    d2=abs(diff(arr)+360);
    d3=abs(diff(arr)-360);  
    c=min(min(d1,d2),d3);
    c=sum(max(c,[],2));
end