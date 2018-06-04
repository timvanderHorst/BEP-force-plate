function printOptimisationResults(R0, R1, V0, V1)
fprintf('       Start guess                       Result\n');
fprintf('      -------------                     --------\n');

for rr = 1 : 3
    fprintf('%6.3f  %6.3f  %6.3f',R0(:,rr));
    if(rr == 2)
        fprintf('    -->    ');
    else
        fprintf('           ');
    end
    fprintf( '%6.3f  %6.3f  %6.3f\n',R1(:,rr));
end
fprintf( '\n');
for rr = 1 : 3
    fprintf('               %7.3f',V0(rr));
    if(rr == 2)
        fprintf('    -->    ');
    else
        fprintf('           ');
    end
    fprintf( '%7.3f\n',V1(rr));
end
end