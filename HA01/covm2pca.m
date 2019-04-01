function [pc,eigval] = covm2pca(covm)
    [eigVec, eigval] = eig(covm);
    [eigval,I] = sort(diag(eigval),'descend');
    pc = eigVec(:,I);
end

