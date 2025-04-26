function field = solve(coeff, field, num_of_iters, sweep_dirs, relax_factor)
    sz=size(field);
    sweep_dirs_len=numel(sweep_dirs);

    % The system of equation is solved iteratively num_of_iters number of times
    for s=1:num_of_iters
        % TDMA algorithm is used to calculate the new temperatures for each of the columns seperately
        switch sweep_dirs(mod(s-1,sweep_dirs_len)+1)
            case 1
                for i=1:sz(2)
                    % Calculating the right hand side of the equation
                    rhs=coeff(:,i,6);

                    % Making sure that the index is in-bounds 
                    if i~=1
                        rhs=rhs+coeff(:,i,3).*field(:,i-1);
                    end
                    if i~=sz(2)
                        rhs=rhs+coeff(:,i,4).*field(:,i+1);
                    end

                    if relax_factor==1
                        field(:,i)=tdma(-coeff(:,i,1),coeff(:,i,5),-coeff(:,i,2),rhs);
                    else
                        field(:,i)=field(:,i)*(1-relax_factor)+relax_factor*tdma(-coeff(:,i,1),coeff(:,i,5),-coeff(:,i,2),rhs);
                    end
                end
            case 2
                for i=sz(2):-1:1
                    % Calculating the right hand side of the equation
                    rhs=coeff(:,i,6);

                    % Making sure that the index is in-bounds 
                    if i~=1
                        rhs=rhs+coeff(:,i,3).*field(:,i-1);
                    end
                    if i~=sz(2)
                        rhs=rhs+coeff(:,i,4).*field(:,i+1);
                    end

                    if relax_factor==1
                        field(:,i)=tdma(-coeff(:,i,1),coeff(:,i,5),-coeff(:,i,2),rhs);
                    else
                        field(:,i)=field(:,i)*(1-relax_factor)+relax_factor*tdma(-coeff(:,i,1),coeff(:,i,5),-coeff(:,i,2),rhs);
                    end
                end
            case 3
                for i=1:sz(1)
                    % Calculating the right hand side of the equation
                    rhs=coeff(i,:,6);

                    % Making sure that the index is in-bounds 
                    if i~=1
                        rhs=rhs+coeff(i,:,1).*field(i-1,:);
                    end
                    if i~=sz(1)
                        rhs=rhs+coeff(i,:,2).*field(i+1,:);
                    end

                    if relax_factor==1
                        field(i,:)=tdma(-coeff(i,:,3),coeff(i,:,5),-coeff(i,:,4),rhs);
                    else
                        field(i,:)=field(i,:)*(1-relax_factor)+relax_factor*tdma(-coeff(i,:,3),coeff(i,:,5),-coeff(i,:,4),rhs);
                    end
                end
            case 4
                for i=1:sz(1)
                    % Calculating the right hand side of the equation
                    rhs=coeff(i,:,6);

                    % Making sure that the index is in-bounds 
                    if i~=1
                        rhs=rhs+coeff(i,:,1).*field(i-1,:);
                    end
                    if i~=sz(1)
                        rhs=rhs+coeff(i,:,2).*field(i+1,:);
                    end
                    
                    if relax_factor==1
                        field(i,:)=tdma(-coeff(i,:,3),coeff(i,:,5),-coeff(i,:,4),rhs);
                    else
                        field(i,:)=field(i,:)*(1-relax_factor)+relax_factor*tdma(-coeff(i,:,3),coeff(i,:,5),-coeff(i,:,4),rhs);
                    end
                end
        end
    end
end