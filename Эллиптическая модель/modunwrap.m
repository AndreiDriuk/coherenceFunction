%% function eliminates phase jump exceeding an arbitrary value
function Res = modunwrap(A, tol, dim)

Res = A;

    if dim == 1 
        diff = A(:, 2:end)-A(:, 1:end-1);
        log = abs(diff)>=tol;
        [row, col] = find(log);
        
        for i = 1:length(row)
                Res(row(i),:) =[Res(row(i), 1:col(i)), Res(row(i),col(i)+1:end )-pi*sign(Res(row(i), col(i)+1)-Res(row(i), col(i)))];
        end
    else
        diff = A(2:end,:)-A(1:end-1,:);
        log = abs(diff)>=tol;
        [row, col] = find(log);
        
        for i = 1:length(col)
                Res(:,col(i)) =[Res(1:row(i), col(i)); Res(row(i)+1:end, col(i) )-pi*sign(Res(row(i)+1, col(i))-Res(row(i), col(i)))];
        end
        
    end
    
    

end