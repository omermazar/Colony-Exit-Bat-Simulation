function  [ X0, Y0, Teta0, T0] =  InitPosition(RandomTime, Xmax, Ymax, xyResolution, Terrain, AllParams, BatOrPrey, ObjectNum)
% Function [] = InitPosition()
% will return the starting posinion and angle s of the prey or bat

%%% Nov2021

TestMode=  AllParams.SimParams.TestMode;

isIn = 1; %flag, inside structure or not

sz = size(Terrain);
sz=sz(1); %How many objects inside the strucure, first is always the walls
   
switch  TestMode
    case 'swarm'
        X0 = round(AllParams.SimParams.swarm_initial_xbat_max*rand(1) ./ xyResolution);
        Y0 = round(AllParams.SimParams.swarm_initial_ybat_std*randn(1) ./ xyResolution);
        Teta0 = AllParams.SimParams.swarm_initial_teta_mean + AllParams.SimParams.swarm_initial_teta_std*randn(1) ;
    
%     case 'caveExit'
%         %% New Jun22
%         if any(strcmp(AllParams.TerrainParams.Type, {'Room-U', 'Room-L'}))
%             Ymin = 0;
%             maxX = Xmax-0.5; minX = Xmax-2;
%             maxY = Ymin+1.5; minY = Ymin+1.5;
%             X0 = (minX + (maxX- minX)*rand(1,1)) ./ xyResolution;
%             Y0 = (minY + (maxY- minY)*rand(1,1)) ./ xyResolution;
%             Teta0 = 5*pi/6 + pi/3*rand(1);
% 
%         end % if any()
        %%
    otherwise
        %% New Jun22
        if contains(AllParams.TerrainParams.Type, {'Room-U', 'Room-L'}) &&  strcmp(TestMode, 'caveExit') && strcmp(BatOrPrey, 'Bat')
%             any(strcmp(AllParams.TerrainParams.Type, {'Room-U', 'Room-L'})) &&  strcmp(TestMode, 'caveExit') && strcmp(BatOrPrey, 'Bat')
            %%% random at starting zone
            Ymin = 0;
            % New June 2024 - Inition Bats in caves are at least 25 cm off walls
            % maxX = Xmax-0.5; minX = Xmax-2;
            % maxY = Ymin+0.5; minY = Ymin+2.5;
            maxX = Xmax-0.75; minX = Xmax-2.25;
            maxY = Ymin+0.75; minY = Ymin+2.75;
            X0 = (minX + (maxX- minX)*rand(1,1)) ./ xyResolution;
            Y0 = (minY + (maxY- minY)*rand(1,1)) ./ xyResolution;
            Teta0 = 5*pi/6 + pi/3*rand(1);
        
        % Simple Cave
        elseif contains(AllParams.TerrainParams.Type, {'Simple Cave'}) &&  strcmp(TestMode, 'caveExit') && strcmp(BatOrPrey, 'Bat')
            Ymin = 0;
            maxX = 5; minX = 5.5;
            maxY = 3; minY = 3.5;
            X0 = (minX + (maxX- minX)*rand(1,1)) ./ xyResolution;
            Y0 = (minY + (maxY- minY)*rand(1,1)) ./ xyResolution;
            Teta0 = 2*pi*rand(1);
         
        % Simple 'Free Space'
        elseif contains(AllParams.TerrainParams.Type, {'Free Space'})
            X0 = round((0.1*Xmax + 0.8*Xmax*rand(1)) ./ xyResolution);
            Y0 = round((0.1*Ymax + 0.8*Ymax*rand(1))./ xyResolution);
            Teta0 = -pi + 2*pi*rand(1);
       
        else %%% Random at any point in the room
            while isIn == 1

                X0 = round((0.1*Xmax + 0.8*Xmax*rand(1)) ./ xyResolution);
                Y0 = round((0.1*Ymax + 0.8*Ymax*rand(1))./ xyResolution);
                Teta0 = -pi + 2*pi*rand(1);
                
                % check if start point is in empty terrain
                if sz>2 %If not only walls present
                    isIn = 0;
                    for i=3:+2:sz-1 %for every object loop
                        if inpolygon(X0, Y0, Terrain(i,:), Terrain(i+1,:)) == 1
                            isIn=1;
                            % check also for minimum distance
                        elseif p_poly_dist(X0, Y0, Terrain(i,:), Terrain(i+1,:)) < AllParams.BatFlightParams.MinDistanceAllowed/xyResolution *2
                            isIn=1;
                        end
                    end
                end
            end
        end % end % if any()
end % switch

T0 = round(RandomTime*rand(1));