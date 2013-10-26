function [Ev,Ep,Ov,Op,V,P] = readVTK(vtkfile)
%for VTK files starting with 0000. e.g.
%VTKfile.0000.vtk,VTKfile.0001.vtk...
%NOT
%VTKfile.0001.vtk,VTKfile.0002.vtk....

for ii = 1:10
      
       
  %fid = fopen(vtkfile,'r','b');
  fid = fopen(sprintf('%s%02d.vtk',vtkfile,ii-1),'r','b');

  if fid == -1
    break;
  end

  fgetl(fid); % # vtk Datafile Version 3.6
  fgetl(fid); % comments
  fgetl(fid); % BINARY
  fgetl(fid); % DATASET UNSTRUCTURED_GRID
  s = fgetl(fid); % POINTS N double

  n_pts = sscanf(s, '%*s %d %*s');

  % V = fread(fid,prod(n_pts*3),'double');

  while (~feof(fid))
    s = fgetl(fid);
    found = strfind(s,'VECTORS u double');
    if (found)
      V(:,ii) = fread(fid,n_pts*3,'double');
      %V = reshape(V,3,n_pts)';
      break;
    end
  end

  while (~feof(fid))
    s = fgetl(fid);
    found = strfind(s,'SCALARS pressure double');
    if (found)
      s = fgetl(fid); % LOOKUP_TABLE default
      P(:,ii) = fread(fid,n_pts,'double');
      %V = reshape(V,3,n_pts)';
      break;
    end
  end

  fclose(fid);

end


%% compute errors in pressure and velocity
if (size(P,2) < 3)
  return ;
end


ncols = size(P,2);

Rtable = [1       ,  0    , 0  , 0    , 0    ;
          2       , -1    , 0  , 0    , 0    ;
          8/3     , -2    , 1/3, 0    , 0    ;
          64/21   , -8/3  , 2/3, -1/21, 0    ;
          1024/315, -64/21, 8/9, -2/21, 1/315
         ];


pref = zeros(size(P(:,end)));
vref = zeros(size(V(:,end)));

Npts = 3;

for k=1:min(ncols,Npts)
  pref = pref + Rtable(min(ncols,Npts),k)*P(:,end-k+1);
  vref = vref + Rtable(min(ncols,Npts),k)*V(:,end-k+1);
end


%if (ncols == 1)
%  pref = P(:,end);
%  vref = V(:,end);
%elseif (ncols == 2)
%  pref = 2*P(:,end) - P(:,end-1);
%  vref = 2*V(:,end) - V(:,end-1);
%elseif (ncols >= 3)
%  pref = 8*P(:,end)/3 - 2*P(:,end-1) + P(:,end-2)/3;
%  vref = 8*V(:,end)/3 - 2*V(:,end-1) + V(:,end-2)/3;
%elseif (ncols == 4)
%  pref = 64*P(:,end)/21 - 8*P(:,end-1)/3 + 2*P(:,end-2)/3 - P(:,end-3)/21;
%  vref = 64*V(:,end)/21 - 8*V(:,end-1)/3 + 2*V(:,end-2)/3 - V(:,end-3)/21;
%end


%pref = P(:,end);
%vref = V(:,end);

if (0)
  for k=1:size(P,2)
  
    Ep(k) = norm(P(:,k) - pref,Inf);
    Ev(k) = norm(V(:,k) - vref,Inf);

  end

  Ep=Ep';
  Ev=Ev';

  for k=1:length(Ep)-1

    Op(k,1) = log2( Ep(k)/Ep(k+1) );
    Ov(k,1) = log2( Ev(k)/Ev(k+1) );

  end
else

  for k=1:size(P,2)-1
    Ep(k) = norm(P(:,k+1) - P(:,k),Inf);
    Ev(k) = norm(V(:,k+1) - V(:,k),Inf);
  end

  for k=1:length(Ep)-1

    Op(k,1) = log2( Ep(k)/Ep(k+1) );
    Ov(k,1) = log2( Ev(k)/Ev(k+1) );

  end

end













