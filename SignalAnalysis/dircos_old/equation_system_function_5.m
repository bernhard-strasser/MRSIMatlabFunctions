function equation_system = equation_system_function_5(dircos)

global z

equation_system = [dircos(1)^2 + dircos(4)^2 + z(1)^2 - 1; ...
                   dircos(2)^2 + dircos(5)^2 + z(2)^2 - 1; ...
                   dircos(3)^2 + dircos(6)^2 + z(3)^2 - 1; ...
                   dircos(1)^2 + dircos(2)^2 + dircos(3)^2 - 1; ...
                   dircos(4)^2 + dircos(5)^2 + dircos(6)^2 - 1; ...
                   dircos(1)*dircos(2) + dircos(4)*dircos(5) + z(1)*z(2); ...
                   dircos(1)*dircos(3) + dircos(4)*dircos(6) + z(1)*z(3); ...
                   dircos(2)*dircos(3) + dircos(5)*dircos(6) + z(2)*z(3); ...
                   dircos(1)*dircos(4) + dircos(2)*dircos(5) + dircos(3)*dircos(6); ...
                   dircos(1)*z(1) + dircos(2)*z(2) + dircos(3)*z(3); ...
                   dircos(4)*z(1) + dircos(5)*z(2) + dircos(6)*z(3); ...
                   dircos(2)*dircos(6) - dircos(3)*dircos(5) - z(1); ...
                   dircos(3)*dircos(4) - dircos(1)*dircos(6) - z(2); ...
                   dircos(1)*dircos(5) - dircos(2)*dircos(4) - z(3)];