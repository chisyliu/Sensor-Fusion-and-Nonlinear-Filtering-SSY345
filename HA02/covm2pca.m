function [pc,eigval] = covm2pca(covm)
    % calculate eigenvector of the covariance matrix
    [eigVec, eigval] = eig(covm);
    % sort eigenvectors by descendent order of eigenvalues
    % the eigenvector with the highest eigenvalue is the major axis
    [eigval,I] = sort(diag(eigval),'descend');
    pc = eigVec(:,I);
end

