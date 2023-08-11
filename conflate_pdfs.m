function conflated = conflate_pdfs(x, pdfs)
    conflated = ones(size(x)) / (x(end) - x(1));
    for i = 1:length(pdfs)
        pdf = pdfs{i} / trapz(x, pdfs{i});
        conflated = conflated .* pdf;
    end
    conflated = conflated / trapz(x, conflated);
end

