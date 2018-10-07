function M = tensor(varargin)

M = 1;
for j = 1:nargin
  if iscell(varargin{j})
    for k = 1:varargin{j}{2}
      M = kron(M,varargin{j}{1});
    end
  else
    M = kron(M,varargin{j});
  end
end
