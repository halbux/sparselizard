function matrix = getallpossiblevectorswithsumequaln(numberofterms, n)

% Every row in 'matrix' is a possible length 'numberofterms' vector of integers whose sum is n.
%
% This function is written in a recursive way!


if n == 0
    matrix = zeros(1, numberofterms);
end

if numberofterms == 0
    matrix = [];
else

    % numberofterms = 1 will be the stopping criterion for the recursive call:
    if numberofterms > 1

        currentindex = 1;
        for i = (0:n)

            allpossiblevectorsfornminusi = getallpossiblevectorswithsumequaln(numberofterms - 1, n - i);

            matrix(currentindex:currentindex + size(allpossiblevectorsfornminusi, 1) - 1, 1) = i;
            matrix(currentindex:currentindex + size(allpossiblevectorsfornminusi, 1) - 1, 2:1 + size(allpossiblevectorsfornminusi, 2)) = allpossiblevectorsfornminusi;    

            currentindex = currentindex + size(allpossiblevectorsfornminusi, 1);

        end

    else

        matrix = n;

    end
    
end

end