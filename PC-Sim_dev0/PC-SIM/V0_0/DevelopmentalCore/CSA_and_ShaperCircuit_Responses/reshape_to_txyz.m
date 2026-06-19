function Q_txyz = reshape_to_txyz(Q)

% Input: (x,y,z,t)
% Output: (t,x,y,z)

Q_txyz = permute(Q, [4 1 2 3]);

end
