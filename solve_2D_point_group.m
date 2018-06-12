% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luis Rangel DaCosta
% Extra Credit Project 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ ME:
% This script needs to be run in version 2016b of Matlab (or newer) as it
% utilizes several instances of implicit expansion.
%
% The program will request the name of the file to be taken 
% in as input through the standard input interface. Script outputs the 
% result to the command window/interface. 
% If point group is unable to be solved (I could not implement solving
% glide symmetries on time), prints out the result 'unknown'.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all
format long
prompt = 'Please enter the name of the file with the structure parameters: ';
name = input(prompt,'s');
param = load(name);

%creation of super lattice to enable examination of local symmetries and
%geometries
matrix = zeros(100,3);

t1 = param(1,1).*[1 0];
gamma = pi*param(1,3)/180;

if gamma == pi/2
   t2 = param(1,2).*[0 1];
elseif gamma > pi/2
   theta = gamma - pi/2; 
   t2 = [-1 -1/tan(theta)]; t2 = param(1,2).*t2./norm(t2); 
else
   t2 = [1 -tan(gamma)]; t2 = param(1,2).*t2./norm(t2); 
end

count = 1;
for iter = 2:size(param,1)
    for j = -5:1:5
        for k = -5:1:5
            matrix(count,1) = param(iter,1);
            matrix(count,2:3) = j.*t1+k.*t2 + param(iter,2).*t1+param(iter,3).*t2;
            count = count + 1;
        end
    end
end

%finds atom closest to origin of superlattice to define as current orgin;
o_x = find(matrix(:,2) == 0);
o_id = o_x(find(matrix(o_x,3) == 0));

dist_mat = sqrt(matrix(:,2).^2+matrix(:,3).^2);
if isempty(o_id)
    [~,o_id] =min(dist_mat);
    dist_mat = sqrt((matrix(:,2)-matrix(o_id,2)).^2+(matrix(:,3)-matrix(o_id,3)).^2);
end
[~,nearest_ids] = sort(dist_mat);

%adding 20 empty space midpoints as possible lattice points
midpoints = get_midpoints([0 0],matrix(nearest_ids(1:20),2:3));

matrix = round([matrix; midpoints].*10000)./10000;

dist_mat = sqrt(matrix(:,2).^2+matrix(:,3).^2);
%grabs list of abitrary number of nearest neighbors to point chosen as
%original origin
[nearest_dists,nearest_ids] = sort(dist_mat);

highest_rot_order = 1;
%holders for 2, 4, 3 fold rotational symmetry indentifiers
rot_id_matrix = zeros(40,3);
rotations = [pi pi/2 2*pi/3];

%solves for rotational symmetries of all chosen nearest neighbors and atom
%midpoints
for iter = 1:40
    local_origin = nearest_ids(iter);
    local_ids = return_nearest_unique_neighbor(matrix,local_origin);
    for j = 1:3
        for k = 1:length(local_ids)
            check_reflect = rotate_coord(rotations(j),matrix(local_ids(k),2:3),matrix(local_origin,2:3));
            check1 = find( abs(matrix(1:end-20,2)-check_reflect(1)) < 1e-2);
            check2 = find( abs(matrix(check1,3)-check_reflect(2)) < 1e-2);
            if ~isempty(check2)
                if matrix(check1(check2),1) == matrix(local_ids(k),1)
                    rot_id_matrix(iter,j) = rot_id_matrix(iter,j)+1/length(local_ids);
                end
            end
        end
    end
end

rot_id_matrix = round(rot_id_matrix.*1000)./1000;
rot_id_matrix = ~(rot_id_matrix-1);
%labels for 2, 3, 4, 6 fold 
rot_id_mat_2 = [rot_id_matrix(:,1) rot_id_matrix(:,3) rot_id_matrix(:,2) rot_id_matrix(:,1).*rot_id_matrix(:,3)]; 

%chooses random origin point of set of points with highest rotational order
%if lattice points without atoms are in set with points with items, points
%with atoms are preffered

list = nearest_ids(1:40);

if ~isempty(find(rot_id_mat_2(:,4)))
    valid_origins = list(find(rot_id_mat_2(:,4)));
    option = 6;
elseif ~isempty(find(rot_id_mat_2(:,3)))
    valid_origins = list(find(rot_id_mat_2(:,3)));
    option = 4;
elseif ~isempty(find(rot_id_mat_2(:,2)))
    valid_origins = list(find(rot_id_mat_2(:,2)));
    option = 3;
elseif ~isempty(find(rot_id_mat_2(:,1)))
    valid_origins = list(find(rot_id_mat_2(:,1)));
    option = 2;
else %no rotational symmetries
    option = 1;
    valid_origins = o_id;
    origin_id = o_id;
end

valid_types = matrix(valid_origins,1);
check = ~~(valid_types);
    
if any(check) && (length(valid_origins)>2)
   valid_with_atom = valid_origins(check);
   valid_no_atom = valid_origins(~check);
   dist = sqrt(matrix(valid_with_atom,2).^2+matrix(valid_with_atom,3).^2);
   [~, order] = sort(dist);
   origin_id = valid_with_atom(order(1));
elseif (length(valid_origins)>2) %all non-atom origins by defintion close to origin, eachother
   origin_id = valid_origins(randi(length(valid_origins),1));
end

%origin is now established
%now will switch through the possible plane lattices according to
%rotiational symmetry and search for reflection, glide lines; if searching
%within no reflection lattices, will establish translational vectors first
%then 

%close all
flag = [];
switch option
    case 6 %4 fold rotational symmetry
       %midpoint of possible origins
       new_origin = valid_origins(valid_origins ~= origin_id);
       new_origin = new_origin(randi(length(new_origin)));
       mirror = matrix(new_origin,2:3)-matrix(origin_id,2:3);
       
       check = is_mirror(mirror,origin_id,matrix);
       
       if check
           point_group = 'p6mm';
       else
           point_group = 'p6';
       end
    case 4 %4 fold rotational symmetry
      dist = sqrt((matrix(valid_origins,2)-matrix(origin_id,2)).^2+(matrix(valid_origins,3)-matrix(origin_id,3)).^2);
      [dist, order] = sort(dist);
      new_origin_1 = valid_origins(order(2));
      [~,next_closest] = unique(dist);
      new_origin_2 = valid_origins(order(next_closest(3)));

      mirror1 = matrix(new_origin_1,2:3)-matrix(origin_id,2:3);
      check1 = is_mirror(mirror1,origin_id,matrix);
      if check1
          point_group = 'p4mm';
      else
          mirror2 = matrix(new_origin_2,2:3)-matrix(origin_id,3:3);
          check2 = is_mirror(mirror2,origin_id,matrix);
          if check2
              point_group = 'p4gm';
          else
              point_group = 'p4';
          end
      end
      
    case 3 %need to check lines between empty space origins and atom origins
      new_origin_no_atom_1 = valid_no_atom(1);
      dist = [0 0 0]';
      for iter = 1:length(valid_no_atom)
          dist(iter) = norm(matrix(new_origin_no_atom_1,2:3)-matrix(valid_no_atom(iter),2:3));
      end
      [~,order] = sort(dist);
      new_origin_no_atom_2 = valid_no_atom(order(2));
      
      mirror1 = matrix(new_origin_no_atom_1,2:3) - matrix(origin_id,2:3);
      
      check1 = is_mirror(mirror1,origin_id,matrix);
      if check1
          point_group = 'p3m1';
      else
          mirror2 = matrix(new_origin_no_atom_1,2:3) - matrix(new_origin_no_atom_2,2:3);
          check2 = is_mirror(mirror2,origin_id,matrix);
          if check2
              point_group = 'p31m';
          else
              point_group = 'p3';
          end
      end
    case 2 %2-fold rotational symmetry
      dist = sqrt((matrix(valid_with_atom,2)-matrix(origin_id,2)).^2+(matrix(valid_with_atom,3)-matrix(origin_id,3)).^2);
      [~, order] = sort(dist);
      new_atom_origin_1 = valid_with_atom(order(2));
      [~,next_closest] = unique(dist);
      iter = 3;
      
      while(matrix(origin_id,1) ~= matrix(valid_with_atom(order(next_closest(iter))),1))
        iter = iter+1;
        new_atom_origin_2 = valid_with_atom(order(next_closest(iter)));
      end
      
      dist = sqrt((matrix(valid_no_atom,2)-matrix(origin_id,2)).^2+(matrix(valid_no_atom,3)-matrix(origin_id,3)).^2);
      [~, order] = sort(dist);
      new_vacant_origin = valid_no_atom(order(2));
      
      mirror1 = matrix(origin_id,2:3) - matrix(new_atom_origin_1,2:3);
      mirror2 = matrix(origin_id,2:3) - matrix(new_atom_origin_2,2:3);
      
      check1 = is_mirror(mirror1,origin_id,matrix);
      check2 = is_mirror(mirror2,origin_id,matrix);
      check3 = is_mirror(mirror2,new_vacant_origin,matrix);
      if check1
          if check2
              check_list = zeros(length(valid_origins),1);
              for i = 1:length(check_list)
                  check_list(i) = is_mirror(mirror1,valid_origins(i),matrix);
              end
              
              if all(check_list)
                 point_group = 'p2mm'; 
              else
                 point_group = 'c2mm';
              end
          else
              point_group = 'p2mg';
          end
      else
        point_group = 'unknown';
      end
    case 1
      point_group = 'unknown';
end

%prints result
fprintf(strcat('The point group is: \n',point_group,'.\n'))

function [near_ids] = return_nearest_unique_neighbor(matrix, current_id)
    current_coord = matrix(current_id,2:3);
    matrix(:,2:3) = matrix(:,2:3)-current_coord;
    dist_mat = sqrt(matrix(1:end-20,2).^2+matrix(1:end-20,3).^2);
    [~,nearest_ids] = sort(dist_mat);
    near_ids = nearest_ids(2:13);
end

function [new] = rotate_coord(theta, coord, origin)
    Rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    coord = coord(:)-origin(:);
    new = Rot*coord+origin(:);
    new = round(new.*1000)./1000; %truncate decimals
end

function [midpoints] = get_midpoints(origin, nearest_neighbors)
    midpoints = zeros(length(nearest_neighbors(2:end,1)),3);
    for i = 2:length(nearest_neighbors)
        midpoints(i-1,2:3) = (origin+nearest_neighbors(i,1:2))./2;
    end
end

function [flag] = is_mirror(mirror,origin_id,matrix)
       if any(~mirror)
           points = return_nearest_unique_neighbor(matrix, origin_id);
           points_mat = matrix(points,:);
           points_coords = (points_mat(:,2:3)-matrix(origin_id,2:3))';
           new_points_coords = (round(100.*(points_coords))./100)';

           new_mat_points_coords = (round(100.*(((matrix(1:end-20,2:3)-matrix(origin_id,2:3))')))./100)';
           
           if ~mirror(1)
               pos = new_points_coords(new_points_coords(:,1) > 0,:);
               neg = new_points_coords(new_points_coords(:,1) < 0,:);
               pos_id = points(new_points_coords(:,1) > 0);
               neg_id = points(new_points_coords(:,1) < 0);

               flag = [];
               for i = 1:length(pos)
                   check = find(new_mat_points_coords(:,1) == -pos(i,1));
                   if ~isempty(check)
                       check2 = find(new_mat_points_coords(check,2) == pos(i,2));
                       if isempty(check2)
                           flag = 0;
                           return;
                       elseif (matrix(check(check2),1) ~= matrix(pos_id(i),1))
                           flag = 0;
                           return;
                       end       
                   else
                       flag = 0;
                       return;
                   end
               end

               if isempty(flag)
                   for i = 1:length(neg)
                       check = find(new_mat_points_coords(:,1) == neg(i,1));
                       if ~isempty(check)
                           check2 = find(new_mat_points_coords(check,2) == neg(i,2));
                           if isempty(check2)
                               flag = 0;
                               return;
                       elseif (matrix(check(check2),1) ~= matrix(neg_id(i),1))
                           flag = 0;
                           return;                               
                           end       
                       else
                           flag = 0;
                           return;
                       end
                   end
               end 

               if isempty(flag)
                   flag = 1;
                   return;
               end
           else
               pos = new_points_coords(new_points_coords(:,2) > 0,:);
               neg = new_points_coords(new_points_coords(:,2) < 0,:);
               pos_id = points(new_points_coords(:,2) > 0);
               neg_id = points(new_points_coords(:,2) < 0);
               
               flag = [];
               for i = 1:length(pos)
                   check = find(new_mat_points_coords(:,1) == pos(i,1));
                   if ~isempty(check)
                       check2 = find(new_mat_points_coords(check,2) == -pos(i,2));
                       if isempty(check2)
                           flag = 0;
                           return;
                       elseif (matrix(check(check2),1) ~= matrix(pos_id(i),1))
                           flag = 0;
                           return;                           
                       end       
                   else
                       flag = 0;
                       return;
                   end
               end

               if isempty(flag)
                   for i = 1:length(neg)
                       check = find(new_mat_points_coords(:,1) == neg(i,1));
                       if ~isempty(check)
                           check2 = find(new_mat_points_coords(check,2) == -neg(i,2));
                           if isempty(check2)
                               flag = 0;
                               return;
                       elseif (matrix(check(check2),1) ~= matrix(neg_id(i),1))
                           flag = 0;
                           return;                               
                           end       
                       else
                           flag = 0;
                           return;
                       end
                   end
               end 

               if isempty(flag)
                   flag = 1;
               end    
           end
       else
           theta = atan(mirror(1)/mirror(2));
           rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
           vec = rot*mirror(:); vec = round(1000*vec./(norm(vec)))./1000;
           points = return_nearest_unique_neighbor(matrix, origin_id);
           points_mat = matrix(points,:);
           points_coords = (points_mat(:,2:3)-matrix(origin_id,2:3))';
           new_points_coords = (round(100.*(rot*points_coords))./100)';
           
           new_mat_points_coords = (round(100.*(rot*((matrix(1:end-20,2:3)-matrix(origin_id,2:3))')))./100)';
           
           pos = new_points_coords(new_points_coords(:,1) > 0,:);
           neg = new_points_coords(new_points_coords(:,1) < 0,:);
           pos_id = points(new_points_coords(:,1) > 0);
           neg_id = points(new_points_coords(:,1) < 0);
           
           flag = [];
           for i = 1:length(pos)
               check = find(new_mat_points_coords(:,1) == -pos(i,1));
               if ~isempty(check)
                   check2 = find(new_mat_points_coords(check,2) == pos(i,2));
                   if isempty(check2)
                       flag = 0;
                       return;
                   elseif (matrix(check(check2),1) ~= matrix(pos_id(i),1))
                       flag = 0;
                       return;
                   end       
               else
                   flag = 0;
                   return;
               end
           end
           
           if isempty(flag)
               for i = 1:length(neg)
                   check = find(new_mat_points_coords(:,1) == -neg(i,1));
                   if ~isempty(check)
                       check2 = find(new_mat_points_coords(check,2) == neg(i,2));
                       if isempty(check2)
                           flag = 0;
                           return;
                       elseif (matrix(check(check2),1) ~= matrix(neg_id(i),1))
                           flag = 0;
                           return;       
                       end       
                   else
                       flag = 0;
                       return;
                   end
               end
           end 
           
           if isempty(flag)
               flag = 1;
           end
       end
end
