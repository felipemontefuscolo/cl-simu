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

for k=1:size(P,2)-1

  Ep(k) = norm(P(:,k) - P(:,end),'inf');
  Ev(k) = norm(V(:,k) - V(:,end),'inf');

end

Ep=Ep';
Ev=Ev';

for k=1:length(Ep)-1  

  Op(k,1) = log2( Ep(k)/Ep(k+1) );
  Ov(k,1) = log2( Ev(k)/Ev(k+1) );

end
  

















