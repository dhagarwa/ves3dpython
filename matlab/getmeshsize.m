function [h_min, h_max] = getmeshsize(p)
%p is patch 
    mat_x = reshape(p.x, [p.Nu, p.Nv]);
    mat_y = reshape(p.y, [p.Nu, p.Nv]);
    mat_z = reshape(p.z, [p.Nu, p.Nv]);
    h_min = zeros(p.Nu, p.Nv);
    h_max = zeros(p.Nu, p.Nv);
    for row = 1:p.Nu
        for col=1:p.Nv

            minh = 100000;
            maxh = -100000;
            %upper left
            if (row-1) > 0 && (col-1) > 0
                dist = sqrt((mat_x(row,col) - mat_x(row-1, col-1))^2 + (mat_y(row,col) - mat_y(row-1, col-1))^2 + (mat_z(row,col) - mat_z(row-1, col-1))^2);
                if dist < minh
                    minh = dist;
                end
                if dist > maxh
                    maxh = dist;
                end

            end
            
            %upper
            if (row-1) > 0 && (col) > 0
                dist = sqrt((mat_x(row,col) - mat_x(row-1, col))^2 + (mat_y(row,col) - mat_y(row-1, col))^2 + (mat_z(row,col) - mat_z(row-1, col))^2);
                if dist < minh
                    minh = dist;
                end
                if dist > maxh
                    maxh = dist;
                end

            end
            
            %upper right
            if (row-1) > 0 && (col+1) <= p.Nv
                dist = sqrt((mat_x(row,col) - mat_x(row-1, col+1))^2 + (mat_y(row,col) - mat_y(row-1, col+1))^2 + (mat_z(row,col) - mat_z(row-1, col+1))^2);
                if dist < minh
                    minh = dist;
                end
                if dist > maxh
                    maxh = dist;
                end

            end
            
            %left
            if (row) > 0 && (col-1) > 0
                dist = sqrt((mat_x(row,col) - mat_x(row, col-1))^2 + (mat_y(row,col) - mat_y(row, col-1))^2 + (mat_z(row,col) - mat_z(row, col-1))^2);
                if dist < minh
                    minh = dist;
                end
                if dist > maxh
                    maxh = dist;
                end

            end
            
            %right
            if (row) > 0 && (col+1) <= p.Nv
                dist = sqrt((mat_x(row,col) - mat_x(row, col+1))^2 + (mat_y(row,col) - mat_y(row, col+1))^2 + (mat_z(row,col) - mat_z(row, col+1))^2);
                if dist < minh
                    minh = dist;
                end
                if dist > maxh
                    maxh = dist;
                end

            end
            
            %lower left
            if (row+1) <= p.Nu && (col-1) > 0
                dist = sqrt((mat_x(row,col) - mat_x(row+1, col-1))^2 + (mat_y(row,col) - mat_y(row+1, col-1))^2 + (mat_z(row,col) - mat_z(row+1, col-1))^2);
                if dist < minh
                    minh = dist;
                end
                if dist > maxh
                    maxh = dist;
                end

            end

            %lower
            if (row+1) <= p.Nu && (col) > 0
                dist = sqrt((mat_x(row,col) - mat_x(row+1, col))^2 + (mat_y(row,col) - mat_y(row+1, col))^2 + (mat_z(row,col) - mat_z(row+1, col))^2);
                if dist < minh
                    minh = dist;
                end
                if dist > maxh
                    maxh = dist;
                end

            end

            %lower right
            if (row+1) <= p.Nu && (col+1) <= p.Nv
                dist = sqrt((mat_x(row,col) - mat_x(row+1, col+1))^2 + (mat_y(row,col) - mat_y(row+1, col+1))^2 + (mat_z(row,col) - mat_z(row+1, col+1))^2);
                if dist < minh
                    minh = dist;
                end
                if dist > maxh
                    maxh = dist;
                end

            end  
            
            h_min(row, col) = minh;
            h_max(row, col) = maxh;
        end

    end



end