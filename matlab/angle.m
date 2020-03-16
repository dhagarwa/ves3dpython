%Function to find the angle between two vectors a, b
function val = angle(a, b)
    val = atan2(norm(cross(a,b)), dot(a,b));
        
end
