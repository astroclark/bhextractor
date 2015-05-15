function [U,Values] = pc_evectors2(A,numvec,norm_flag,sort_flag)
%PC_EVECTORS Get the top numvecs eigenvectors of the covariance matrix
%            of A, using Turk and Pentland's trick for numrows >> numcols
%            Returns the eigenvectors as the colums of Vectors and a
%            vector of ALL the eigenvectors in Values.

% Matthew Dailey 2000
% Adjusted by ISH 2007

  % Check arguments

  if nargin ~= 4
    error('usage: pc_evectors2(A,numvecs,norm_flag,sort_flag)');
  end;
  
  nexamp = size(A,2);

  % Now compute the eigenvectors of the covariance matrix, using
  % a little trick from Turk and Pentland 1991

  % Compute the "average" vector
  % mean(A) gives you a row vector containing the mean of each column of A

  fprintf(1,'Computing average vector and vector differences from avg...\n');
  Psi = mean(A')';

  % Compute difference with average for each vector

  %for i = 1:nexamp
  %  A(:,i) = A(:,i) - Psi;
  %end;

  % Get the patternwise (nexamp x nexamp) covariance matrix

  fprintf(1,'Calculating L=A''A\n');
  L = A'*A;

  % Get the eigenvectors (columns of Vectors) and eigenvalues (diag of Values)

  fprintf(1,'Calculating eigenvectors of L...\n');
  [V,tmp] = eig(L); %%% NOTE: if you use [U,S,V]=svd(A), then the SVD output
                    %%% vector, V, is equal to fliplr(V) here!!!
  
  % Get the eigenvalues out of the diagonal matrix and
  % normalize them so the evalues are specifically for cov(A'), not A*A'.
  
  S = tmp.^(1/2);
  Values = diag(tmp);
  %Values = Values/nexamp;
  
  % Sort the vectors/values according to size of eigenvalue
  % NOTE: sorting screws up the order of the original matrix, so you cannot
  % have A = U*S*V' after sorting!!
  
  if sort_flag
      
    fprintf(1,'Sorting evectors/values...\n');
    [Values,ind] = sort(Values,'descend');
    V = V(:,ind);   
    S = flipud(S(:,ind));
  else
    V = fliplr(V);          %flip matrix to look like SVD output
    S = fliplr(flipud(S));
  end
  
  % Convert the eigenvectors of A'*A into eigenvectors of A*A'

  fprintf(1,'Computing eigenvectors of the real covariance matrix..\n');
  U = A*V;
  
  % normalises the eigenvalues so they are of the covariance matrix and not
  % L
  if norm_flag
     for i = 1:nexamp
      U(:,i) = U(:,i)/norm(U(:,i));
     end
  end
  
  % Normalize Vectors to unit length, kill vectors corr. to tiny evalues

%   num_good = 0;
%   for i = 1:nexamp
%     Vectors(:,i) = Vectors(:,i)/norm(Vectors(:,i));
%     if Values(i) < 0.00001
%       % Set the vector to the 0 vector; set the value to 0.
%       Values(i) = 0;
%       Vectors(:,i) = zeros(size(Vectors,1),1);
%     else
%       num_good = num_good + 1;
%     end;
%   end;
%   if (numvecs > num_good)
%     fprintf(1,'Warning: numvecs is %d; only %d exist.\n',numvecs,num_good);
%     numvecs = num_good;
%   end;


U = U(:,1:numvec);

